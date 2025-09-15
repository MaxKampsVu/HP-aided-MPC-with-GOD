#define BOOST_TEST_MODULE utils

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "circuit.h"

using namespace utils;
namespace bdata = boost::unit_test::data;

constexpr int TEST_DATA_MAX_VAL = 1000;

BOOST_AUTO_TEST_SUITE(circuit)

BOOST_DATA_TEST_CASE(no_op_circuit,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  circ.setAsOutput(wa);

  auto output = circ.evaluate({{wa, input}});

  BOOST_TEST(output[0] == input);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 1);
  BOOST_TEST(level_circ.count[GateType::kInp] == 1);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
}

BOOST_DATA_TEST_CASE(add_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wsum = circ.addGate(GateType::kAdd, wa, wb);
  circ.setAsOutput(wsum);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});

  BOOST_TEST(output[0] == input_a + input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 3);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 1);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
}

BOOST_DATA_TEST_CASE(sub_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wdiff = circ.addGate(GateType::kSub, wa, wb);
  circ.setAsOutput(wdiff);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});

  BOOST_TEST(output[0] == input_a - input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 3);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 1);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
}

BOOST_DATA_TEST_CASE(mul_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wprod = circ.addGate(GateType::kMul, wa, wb);
  circ.setAsOutput(wprod);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});

  BOOST_TEST(output[0] == input_a * input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 2);
  BOOST_TEST(level_circ.num_gates == 3);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 1);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
}

BOOST_DATA_TEST_CASE(const_add_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wsum = circ.addConstOpGate(GateType::kConstAdd, wa, input_b);
  circ.setAsOutput(wsum);

  auto output = circ.evaluate({{wa, input_a}});

  BOOST_TEST(output[0] == input_a + input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 2);
  BOOST_TEST(level_circ.count[GateType::kInp] == 1);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 1);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
}

BOOST_DATA_TEST_CASE(const_mul_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wprod = circ.addConstOpGate(GateType::kConstMul, wa, input_b);
  circ.setAsOutput(wprod);

  auto output = circ.evaluate({{wa, input_a}});

  BOOST_TEST(output[0] == input_a * input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 2);
  BOOST_TEST(level_circ.count[GateType::kInp] == 1);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 1);
}

BOOST_DATA_TEST_CASE(depth_2_circuit,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, input_c, input_d, idx) {
  std::vector<int> inputs = {input_a, input_b, input_c, input_d};

  Circuit<int> circ;
  std::vector<wire_t> input_wires;
  for (size_t i = 0; i < inputs.size(); ++i) {
    input_wires.push_back(circ.newInputWire());
  }

  auto w_aab = circ.addGate(GateType::kAdd, input_wires[0], input_wires[1]);
  auto v_aab = inputs[0] + inputs[1];

  auto w_cmd = circ.addGate(GateType::kMul, input_wires[2], input_wires[3]);
  auto v_cmd = inputs[2] * inputs[3];

  auto w_mout = circ.addGate(GateType::kMul, w_aab, w_cmd);
  auto v_mout = v_aab * v_cmd;

  auto w_aout = circ.addGate(GateType::kAdd, w_aab, w_cmd);
  auto v_aout = v_aab + v_cmd;

  circ.setAsOutput(w_mout);
  circ.setAsOutput(w_aout);

  std::unordered_map<wire_t, int> input_map;
  for (size_t i = 0; i < inputs.size(); ++i) {
    input_map[input_wires[i]] = inputs[i];
  }
  auto outputs = circ.evaluate(input_map);

  BOOST_TEST(outputs[0] == v_mout);
  BOOST_TEST(outputs[1] == v_aout);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 3);
  BOOST_TEST(level_circ.num_gates == 8);
  BOOST_TEST(level_circ.count[GateType::kInp] == 4);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 2);
  BOOST_TEST(level_circ.count[GateType::kMul] == 2);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
}

BOOST_AUTO_TEST_SUITE_END()
