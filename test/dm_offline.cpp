#define BOOST_TEST_MODULE dm_offline

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "dm_offline_evaluator.h"

using namespace dmAsyncAsterisk;
namespace bdata = boost::unit_test::data; 

constexpr int TEST_DATA_MAX_VAL = 1000;
constexpr int SECURITY_PARAM = 128;

struct GlobalFixture {
  GlobalFixture() {
    ZZ_p::init(conv<ZZ>("17816577890427308801"));
  }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);


BOOST_AUTO_TEST_SUITE(dm_offline_evaluator)

BOOST_AUTO_TEST_CASE(random_share) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 5;
  
  std::vector<std::future<RepShare<Field>>> parties;  
  parties.reserve(nP+1);
  for (int i = 0; i <= nP; i++) {    
    parties.push_back(std::async(std::launch::async, [&, i]() { 
      ZZ_p_ctx.restore();
      RepShare<Field> shares;
      std::vector<Field> keySh;

      RandGenPool vrgen(i, nP);
      auto tpool = std::make_shared<ThreadPool>(4);
      std::vector<std::unique_ptr<OTProvider>> ot;
      
      auto network1 = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      auto network2 = std::make_shared<NetIOMP>(i, nP+1, 10000+2000, nullptr, true);
      if (i == 0) {
        for (size_t j = 1; j <= nP; j++) {
          ot.emplace_back(std::make_unique<OTProvider>(i, j, network2->getRecvChannel(j)));
        }
      }
      else {
        ot.emplace_back(std::make_unique<OTProvider>(i, 0, network2->getSendChannel(0)));
      }

      if(i == 0) {
        keySh.resize(nP);
        for(size_t j = 1; j <= nP; j++) {
          randomizeZZp(vrgen.pi(j), keySh[j-1], sizeof(Field));
        }
      }
      else {
        keySh.resize(2);
        randomizeZZp(vrgen.p0(), keySh[0], sizeof(Field));
        randomizeZZp(vrgen.all_minus_0(), keySh[1], sizeof(Field));
      }      

      Field temp;
      OfflineEvaluator::randomShare(nP, i, vrgen, shares, keySh, temp, false);

      return shares;
    }));
    
  }
  int i = 0;
  std::vector<RepShare<Field>> v_shares(nP+1);
  for (auto& p : parties) { 
    v_shares[i] = p.get();     
    i++;
  }
  for (size_t i=0; i < nP; i++) {
    BOOST_TEST(v_shares[0].getValues()[i] == v_shares[i+1].getValues()[0]);
  }
  Field temp_val = v_shares[1].getValues()[1];
  for (size_t i=2; i <= nP; i++) {
    BOOST_TEST(v_shares[i].getValues()[1] == temp_val);
  }
}

BOOST_AUTO_TEST_CASE(depth_2_circuit) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 5;
  Circuit<Field> circ;
  std::vector<wire_t> input_wires;
  std::unordered_map<wire_t, int> input_pid_map;
  
  for (size_t i = 0; i < 4; ++i) {
    auto winp = circ.newInputWire();
    input_wires.push_back(winp);
    input_pid_map[winp] = 1;
  }
  
  auto w_aab =
      circ.addGate(GateType::kAdd, input_wires[0], input_wires[1]);
  auto w_cmd =
      circ.addGate(GateType::kMul, input_wires[2], input_wires[3]);
  auto w_cons = circ.addConstOpGate(GateType::kConstAdd, w_aab, Field(2));
  auto w_cons_m = circ.addConstOpGate(GateType::kConstMul, w_cmd, Field(2));
  auto w_mout = circ.addGate(GateType::kMul, w_aab, w_cons);
  auto w_aout = circ.addGate(GateType::kAdd, w_aab, w_cmd);
      // auto w_cons = circ.addConstOpGate(GateType::kConstAdd, w_aout, 2);
  circ.setAsOutput(w_mout);
  circ.setAsOutput(w_aout);
  auto level_circ = circ.orderGatesByLevel();

  std::vector<std::future<PreprocCircuit<Field>>> parties;
  parties.reserve(nP+1);
  for (int i = 0; i <= nP; ++i) {
    parties.push_back(std::async(std::launch::async, [&, i, input_pid_map]() {
      ZZ_p_ctx.restore();
      auto network1 = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      auto network2 = std::make_shared<NetIOMP>(i, nP+1, 10000+2000, nullptr, true);
      RandGenPool vrgen(i, nP);
      OfflineEvaluator eval(nP, i, SECURITY_PARAM, std::move(network1), std::move(network2), level_circ, nP+1);
      return eval.run(input_pid_map);
    }));
  }

  std::vector<PreprocCircuit<Field>> v_preproc;
  v_preproc.reserve(parties.size());
  for (auto& f : parties) {
    v_preproc.push_back(f.get());
  }

  BOOST_TEST(v_preproc[0].gates.size() == level_circ.num_gates);
  const auto& preproc_0 = v_preproc[0];
  std::vector<Field> temp_val, temp_tag;
  for(int j = 0; j < 4; j++) {
    temp_val.push_back(v_preproc[1].gates[j]->mask.getValues()[1]);
    temp_tag.push_back(v_preproc[1].gates[j]->mask.getTags()[1]);
  }
  for (int i = 1; i <= nP; ++i) {
    BOOST_TEST(v_preproc[i].gates.size() == level_circ.num_gates);
    const auto& preproc_i = v_preproc[i];
    for(int j = 0; j < 4; j++) {
      auto tpmask = preproc_0.gates[j]->mask;
      auto mask_i = preproc_i.gates[j]->mask;
      BOOST_TEST(mask_i.getValues()[0] == tpmask.getValues()[i-1]);
      BOOST_TEST(mask_i.getTags()[0] == tpmask.getTags()[i-1]);
      BOOST_TEST(mask_i.getValues()[1] == temp_val[j]);
      BOOST_TEST(mask_i.getTags()[1] == temp_tag[j]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
