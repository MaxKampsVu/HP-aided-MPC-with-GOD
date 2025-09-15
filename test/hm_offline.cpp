#define BOOST_TEST_MODULE hm_offline

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "hm_offline_evaluator.h"

using namespace hmAsyncAsterisk;
namespace bdata = boost::unit_test::data; 

constexpr int TEST_DATA_MAX_VAL = 1000;
constexpr int SECURITY_PARAM = 128;

struct GlobalFixture {
  GlobalFixture() {
    ZZ_p::init(conv<ZZ>("17816577890427308801"));
  }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);


BOOST_AUTO_TEST_SUITE(hm_offline_evaluator)

BOOST_AUTO_TEST_CASE(random_share) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 5;
  
  std::vector<std::future<AuthShamirShare<Field>>> parties;
  parties.reserve(nP+1);  
  TPShamirShare<Field> TPShamirShares;  
  for (int i = 0; i <= nP; i++) {    
    parties.push_back(std::async(std::launch::async, [&, i]() { 
      ZZ_p_ctx.restore();
      AuthShamirShare<Field> shares;

      int th = (nP-1)/2;
      std::vector<std::vector<Field>> rand_sh(nP-th);
      std::vector<size_t> idx(nP-th, 0);
      RandGenPool vrgen(i, nP);

      std::vector<Field> keySh(nP + 1);
      Field key = Field(0);
      
      auto network = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);

      if(i == 0)  {
        key = 0;
        keySh[0] = 0;
        std::vector<Field> keyShs;
        randomizeZZp(vrgen.self(), key, sizeof(Field));
        keyShs.push_back(key);
        for(int j = 1; j <= th; j++) {
          randomizeZZp(vrgen.pi(j), keySh[j], sizeof(Field));
          keyShs.push_back(keySh[j]);
        }            

        std::vector<Field> evalPoints(th + 1);
        std::iota(evalPoints.begin(), evalPoints.end(), 0);
        fieldPoly poly_key = reconstructPolynomial(evalPoints, keyShs);
        for(int j = th+1; j <= nP; j++) {
          keySh[j] = eval(poly_key, Field(j));
          network->send(j, &keySh[j], sizeof(Field));
        }
      }
      else if (i <= th) {
        randomizeZZp(vrgen.p0(), key, sizeof(Field));
      }
      else {
        network->recv(0, &key, sizeof(Field));
      }

      
      if(i <= th) {
      OfflineEvaluator::randomShare(nP, i, vrgen, *network, 
                                 shares, TPShamirShares, key, keySh, rand_sh, idx);
        if(i == 0) {
          for(int j = th+1; j <= nP; j++) {
            size_t rand_sh_num = rand_sh[j-th-1].size();
            network->send(j, &rand_sh_num, sizeof(size_t));
            network->send(j, rand_sh[j-th-1].data(), sizeof(Field) * rand_sh_num);
          }       
        }
      }
      else {
        size_t rand_sh_num;
        network->recv(0, &rand_sh_num, sizeof(size_t));       
        rand_sh[i-th-1].resize(rand_sh_num);
        network->recv(0, rand_sh[i-th-1].data(), sizeof(Field) * rand_sh_num);
        OfflineEvaluator::randomShare(nP, i, vrgen, *network, 
                                  shares, TPShamirShares, key, keySh, rand_sh, idx);
      }
      return shares;
    }));
    
  }
  int i = 0;
  for (auto& p : parties) { 
    auto res = p.get();
    BOOST_TEST(res.valueAt() == TPShamirShares.commonValueWithParty(i));
    BOOST_TEST(res.tagAt() == TPShamirShares.commonTagWithParty(i));     
    i++;
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
      auto network = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      RandGenPool vrgen(i, nP);
      OfflineEvaluator eval(nP, i, SECURITY_PARAM, std::move(network), level_circ, nP+1);
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
  
  for (int i = 1; i <= nP; ++i) {
    BOOST_TEST(v_preproc[i].gates.size() == level_circ.num_gates);
    const auto& preproc_i = v_preproc[i];
    for(int j = 0; j < 4; j++) {
      auto tpmask = preproc_0.gates[j]->tpmask;
      auto mask_i = preproc_i.gates[j]->mask;
      BOOST_TEST(mask_i.valueAt() == tpmask.commonValueWithParty(i));
      BOOST_TEST(mask_i.tagAt() == tpmask.commonTagWithParty(i));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
