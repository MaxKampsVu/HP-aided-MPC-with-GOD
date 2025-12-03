#define BOOST_TEST_MODULE dm_online

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "dmgod_offline_evaluator.h"
#include "dmgod_online_evaluator.h"

using namespace dmAsyncAsteriskGOD;
namespace bdata = boost::unit_test::data;

constexpr int TEST_DATA_MAX_VAL = 1000;
constexpr int SECURITY_PARAM = 128;

struct GlobalFixture {
  GlobalFixture() {
    ZZ_p::init(conv<ZZ>("17816577890427308801"));
  }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_SUITE(dm_online_evaluator)

BOOST_AUTO_TEST_CASE(add) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 4;
  auto seed_block = emp::makeBlock(0, 200);
  emp::PRG prg(&seed_block);
  std::mt19937 gen(rand());
  std::uniform_int_distribution<uint> distrib(0, TEST_DATA_MAX_VAL);
  Circuit<Field> circ;
  std::vector<wire_t> input_wires;
  std::unordered_map<wire_t, int> input_pid_map;
  std::unordered_map<wire_t, Field> inputs;

  for (size_t i = 0; i < 2; ++i) {
    auto winp = circ.newInputWire();
    input_wires.push_back(winp);
    input_pid_map[winp] = 1;    
    inputs[winp] = Field(distrib(gen));
  }
  auto w_amb = circ.addGate(GateType::kAdd, input_wires[0], input_wires[1]);  
  circ.setAsOutput(w_amb);  
  auto level_circ = circ.orderGatesByLevel();

  auto exp_output = circ.evaluate(inputs);

  std::vector<std::future<std::vector<Field>>> parties;
  parties.reserve(nP+1);
  for (int i = 0; i <= nP; ++i) {
      parties.push_back(std::async(std::launch::async, [&, i, input_pid_map, inputs]() {
      ZZ_p_ctx.restore();
      auto network1 = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      auto network2 = std::make_shared<NetIOMP>(i, nP+1, 10000+100, nullptr, true);

      PreprocCircuit<Field> preproc;
      {
        OfflineEvaluator eval(nP, i, SECURITY_PARAM, network1, network2, level_circ, nP+1);
        preproc = eval.run(input_pid_map);
      }
      
      std::vector<Field> res;
      {
        OnlineEvaluator online_eval(nP, i, SECURITY_PARAM, std::move(network1), std::move(preproc), level_circ, nP+1);      
        res = online_eval.evaluateCircuit(inputs);
      }
      
      return res;      
    }));
  } 
  int i = 0;
  for (auto& p : parties) {
    if(i != 0) {
      auto output = p.get();
      BOOST_TEST(exp_output == output);
    }
    i++;
  }
}

BOOST_AUTO_TEST_CASE(mult) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 5;
  auto seed_block = emp::makeBlock(0, 200);
  emp::PRG prg(&seed_block);
  std::mt19937 gen(rand());
  std::uniform_int_distribution<uint> distrib(0, TEST_DATA_MAX_VAL);
  Circuit<Field> circ;
  std::vector<wire_t> input_wires;
  std::unordered_map<wire_t, int> input_pid_map;
  std::unordered_map<wire_t, Field> inputs;

  for (size_t i = 0; i < 2; ++i) {
    auto winp = circ.newInputWire();
    input_wires.push_back(winp);
    input_pid_map[winp] = 1;    
    inputs[winp] = Field(distrib(gen));
  }
  auto w_amb = circ.addGate(GateType::kMul, input_wires[0], input_wires[1]);  
  circ.setAsOutput(w_amb);  
  auto level_circ = circ.orderGatesByLevel();

  auto exp_output = circ.evaluate(inputs);

  std::vector<std::future<std::vector<Field>>> parties;
  parties.reserve(nP+1);
  for (int i = 0; i <= nP; ++i) {
      parties.push_back(std::async(std::launch::async, [&, i, input_pid_map, inputs]() {
      ZZ_p_ctx.restore();
      auto network1 = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      auto network2 = std::make_shared<NetIOMP>(i, nP+1, 10000+100, nullptr, true);

      PreprocCircuit<Field> preproc;
      {
        OfflineEvaluator eval(nP, i, SECURITY_PARAM, network1, network2, level_circ, nP+1);
        preproc = eval.run(input_pid_map);
      }
      
      std::vector<Field> res;
      {
        OnlineEvaluator online_eval(nP, i, SECURITY_PARAM, std::move(network1), std::move(preproc), level_circ, nP+1);      
        res = online_eval.evaluateCircuit(inputs);
      }
      
      return res;      
    }));
  } 
  int i = 0;
  for (auto& p : parties) {
    if(i != 0) {
      auto output = p.get();
      BOOST_TEST(exp_output == output);
    }
    i++;
  }
}

BOOST_AUTO_TEST_CASE(dotp2) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 5;
  auto seed_block = emp::makeBlock(0, 200);
  emp::PRG prg(&seed_block);
  std::mt19937 gen(rand());
  std::uniform_int_distribution<uint> distrib(0, TEST_DATA_MAX_VAL);
  Circuit<Field> circ;
  std::vector<wire_t> input_wires;
  std::unordered_map<wire_t, int> input_pid_map;
  std::unordered_map<wire_t, Field> inputs;


  // input wires for vectors a(2) b(2)
  for (size_t i = 0; i < 4; ++i) {
    auto winp = circ.newInputWire();
    input_wires.push_back(winp);
    input_pid_map[winp] = 1;    
    inputs[winp] = Field(i+1);
  }
  auto w_amb = circ.addGate(GateType::kDotprod, {input_wires[0], input_wires[1]}, {input_wires[2], input_wires[3]});  
  circ.setAsOutput(w_amb);  
  auto level_circ = circ.orderGatesByLevel();

  auto exp_output = circ.evaluate(inputs);

  std::vector<std::future<std::vector<Field>>> parties;
  parties.reserve(nP+1);
  for (int i = 0; i <= nP; ++i) {
      parties.push_back(std::async(std::launch::async, [&, i, input_pid_map, inputs]() {
      ZZ_p_ctx.restore();
      auto network1 = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      auto network2 = std::make_shared<NetIOMP>(i, nP+1, 10000+100, nullptr, true);

      PreprocCircuit<Field> preproc;
      {
        OfflineEvaluator eval(nP, i, SECURITY_PARAM, network1, network2, level_circ, nP+1);
        preproc = eval.run(input_pid_map);
      }
      
      std::vector<Field> res;
      {
        OnlineEvaluator online_eval(nP, i, SECURITY_PARAM, std::move(network1), std::move(preproc), level_circ, nP+1);      
        res = online_eval.evaluateCircuit(inputs);
      }
      
      return res;      
    }));
  } 
  int i = 0;
  for (auto& p : parties) {
    if(i != 0) {
      auto output = p.get();
      std::cout << "size: " << output.size() << std::endl;
      if(output.size() > 0) {
        std::cout << "output: " << output[0] << std::endl;
        std::cout << "exp output: " << exp_output[0] << std::endl;
      } 
      BOOST_TEST(exp_output == output);
    }
    i++;
  }
}

BOOST_AUTO_TEST_CASE(depth_2_circuit) {
  ZZ_pContext ZZ_p_ctx;
  ZZ_p_ctx.save();
  int nP = 10;
  auto seed_block = emp::makeBlock(0, 200);
  emp::PRG prg(&seed_block);
  std::mt19937 gen(rand());
  std::uniform_int_distribution<uint> distrib(0, TEST_DATA_MAX_VAL);
  Circuit<Field> circ;
  std::vector<wire_t> input_wires;
  std::unordered_map<wire_t, int> input_pid_map;
  std::unordered_map<wire_t, Field> inputs;

  for (size_t i = 0; i < 4; ++i) {
    auto winp = circ.newInputWire();
    input_wires.push_back(winp);
    input_pid_map[winp] = 1;    
    inputs[winp] = Field(distrib(gen));
  }
  auto w_aab = circ.addGate(GateType::kAdd, input_wires[0], input_wires[1]);
  auto w_cmd = circ.addGate(GateType::kMul, input_wires[2], input_wires[3]);
  auto w_mout = circ.addGate(GateType::kMul, w_aab, w_cmd);
  auto w_aout = circ.addGate(GateType::kAdd, w_aab, w_cmd);
  circ.setAsOutput(w_cmd);
  circ.setAsOutput(w_mout);
  circ.setAsOutput(w_aout);
  auto level_circ = circ.orderGatesByLevel();

  auto exp_output = circ.evaluate(inputs);
  std::vector<std::future<std::vector<Field>>> parties;
  parties.reserve(nP+1);
  for (int i = 0; i <= nP; ++i) {
      parties.push_back(std::async(std::launch::async, [&, i, input_pid_map, inputs]() {
      ZZ_p_ctx.restore();
      auto network1 = std::make_shared<NetIOMP>(i, nP+1, 10000, nullptr, true);
      auto network2 = std::make_shared<NetIOMP>(i, nP+1, 10000+2000, nullptr, true);

      PreprocCircuit<Field> preproc;
      {
        OfflineEvaluator eval(nP, i, SECURITY_PARAM, network1, network2, level_circ, nP+1);
        preproc = eval.run(input_pid_map);
      }
      
      std::vector<Field> res;
      {
        OnlineEvaluator online_eval(nP, i, SECURITY_PARAM, std::move(network1), std::move(preproc), level_circ, nP+1);      
        res = online_eval.evaluateCircuit(inputs);
      }
           
      return res;      
    }));
  }
  int i = 0;
  for (auto& p : parties) {
    if(i != 0) {
      auto output = p.get();
      if(i > 0) {
        BOOST_TEST(exp_output == output);
      }
      else {
        for (int j=0; j<2; j++) { // because for the last output, HP doesn't have the result.
          BOOST_TEST(exp_output[j] == output[j]);
        }
      }
    }
    i++;
  }
}

BOOST_AUTO_TEST_SUITE_END()
