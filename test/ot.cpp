#define BOOST_TEST_MODULE ot_provider

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "ot_provider.h"
#include "netmp.h"

using namespace asyncAsterisk;
using namespace io;
namespace bdata = boost::unit_test::data;

constexpr int TEST_DATA_MAX_VAL = 1000;
constexpr int SECURITY_PARAM = 128;

struct GlobalFixture {
    GlobalFixture() {
        ZZ_p::init(conv<ZZ>("17816577890427308801"));
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_SUITE(ot_provider)

BOOST_AUTO_TEST_CASE(OPE) {
    ZZ_pContext ZZ_p_ctx;
    ZZ_p_ctx.save();
    std::mt19937 gen(rand());
    std::uniform_int_distribution<uint> distrib(0, TEST_DATA_MAX_VAL);
    int numOPE = 4;

    std::vector<std::vector<Field>> inputToOPE(2);
    std::vector<Field> exp_output;
    std::vector<Field> output(numOPE, Field(0));
    for (size_t i = 0; i < numOPE; ++i) {
        inputToOPE[0].push_back(Field(distrib(gen)));
        inputToOPE[1].push_back(Field(distrib(gen)));
        exp_output.push_back(inputToOPE[0][i] * inputToOPE[1][i]);
    }    

    std::vector<std::future<std::vector<Field>>> parties;
    parties.reserve(2);

    for (int i = 0; i < 2; ++i) {
        parties.push_back(std::async(std::launch::async, [&, i]() {
        ZZ_p_ctx.restore();
        auto network = std::make_shared<NetIOMP>(i, 2, 10000+100, nullptr, true);

        std::unique_ptr<OTProvider> otProvider;
        if (i ==0) {
            otProvider = std::make_unique<OTProvider>(i, 1-i, network->getRecvChannel(1-i));
        }
        else {
            otProvider = std::make_unique<OTProvider>(i, 1-i, network->getSendChannel(1-i));
        } 
        
        if (i == 0) {
            auto res = otProvider->multiplyRecv(inputToOPE[0]); 
            return res;
        } else {
            PRG prg(&zero_block, 200);
            auto res = otProvider->multiplySend(inputToOPE[1], prg); 
            return res;            
        }       
      }));
    }
    for (auto& p : parties) {
      auto shares = p.get();
      for (size_t i = 0; i < numOPE; ++i) {
        output[i] += shares[i];
      }
    }
    BOOST_TEST(exp_output == output);
}

BOOST_AUTO_TEST_SUITE_END()