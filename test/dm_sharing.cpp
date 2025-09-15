#define BOOST_TEST_MODULE dm_sharing

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "types.h"
#include "dm_sharing.h"

using namespace dmAsyncAsterisk;
namespace bdata = boost::unit_test::data;

struct GlobalFixture {
  GlobalFixture() {
    ZZ_p::init(conv<ZZ>("17816577890427308801"));
  }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

constexpr int TEST_DATA_MAX_VAL = 1000;
constexpr int NUM_SAMPLES = 1;
std::random_device rd;
std::mt19937 engine(rand());
std::uniform_int_distribution<uint64_t> distrib;
// Field MAC_key = distrib(engine);

std::vector<RepShare<Field>> generateAuthRepShares(Field secret, size_t nP) {
  Field MAC_key = Field(5);
  std::random_device rd;
  std::mt19937 engine(rd());
  std::uniform_int_distribution<uint64_t> distrib; 
  
  Field tag = secret * MAC_key;

  std::vector<Field> key_shares(nP+1);
  std::vector<Field> values(nP+1);
  std::vector<Field> tags(nP+1);  

  Field sum1 = Field(0), sum2 = Field(0), sum3 = Field(0);
  int i;
  for (i = 0; i < nP; ++i) {
    key_shares[i] = Field(distrib(engine));
    sum1 += key_shares[i];
    values[i] = Field(distrib(engine));
    sum2 += values[i];
    tags[i] = Field(distrib(engine));
    sum3 += tags[i];
  }
  key_shares[i] = MAC_key - sum1;
  values[i] = secret - sum2;
  tags[i] = tag - sum3;

  std::vector<RepShare<Field>> ARS;
  for(size_t i = 1; i <= nP ; i++) {
    std::vector<Field> keyVec, valVec, tagVec;
    keyVec.push_back(key_shares[i]);
    keyVec.push_back(key_shares[0]);
    valVec.push_back(values[i]);
    valVec.push_back(values[0]);
    tagVec.push_back(tags[i]);
    tagVec.push_back(tags[0]);
    RepShare<Field> temp(keyVec, valVec, tagVec);
    ARS.push_back(temp);
  }

  return ARS;
}

Field reconstructAuthRepShares(std::vector<RepShare<Field>>& v_ars, size_t nP) {
  Field secret = v_ars[0].getValues()[1], key = v_ars[0].getKeySh()[1], tag = v_ars[0].getTags()[1];
  for (size_t i=0; i<nP; i++) {
    secret += v_ars[i].getValues()[0];
    key += v_ars[i].getKeySh()[0];
    tag += v_ars[i].getTags()[0];
  }

  if(secret * key == tag) { return secret; }
  else {
    std::cout<< "Incorrect sharing !!!" << std::endl;
    return Field(0);
  }
}

BOOST_AUTO_TEST_SUITE(authenticated_Replicated_sharing)

BOOST_DATA_TEST_CASE(reconstruction, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), secret_val, idx) {
  size_t nP = 4;
  Field secret = Field(secret_val);

  auto v_aas = generateAuthRepShares(secret, nP);
  
  auto recon_value = reconstructAuthRepShares(v_aas, nP);
    
  BOOST_TEST(recon_value == secret);
}

BOOST_DATA_TEST_CASE(share_arithmetic, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), vala, valb, idx) {
  size_t nP = 4;
  Field a = Field(vala);
  Field b = Field(valb);
  auto v_ars_a = generateAuthRepShares(a, nP);
  auto v_ars_b = generateAuthRepShares(b, nP);

  std::vector<RepShare<Field>> v_ars_c(nP);

  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_ars_c[i] = v_ars_a[i] + v_ars_b[i];
  }

  auto sum = reconstructAuthRepShares(v_ars_c, nP);

  BOOST_TEST(sum == a + b);
  
  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_ars_c[i] = v_ars_a[i] - v_ars_b[i];
  }

  auto difference = reconstructAuthRepShares(v_ars_c, nP);
  BOOST_TEST(difference == a - b);  
}

BOOST_DATA_TEST_CASE(share_const_arithmetic, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), secret_val, const_val, idx) {
  size_t nP = 6;
  Field secret = Field(secret_val);
  Field constant = Field(const_val);
  auto v_ars = generateAuthRepShares(secret, nP);

  std::vector<RepShare<Field>> v_ars_res(nP);
  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_ars_res[i] = v_ars[i].add(constant, i+1);
  }

  auto sum = reconstructAuthRepShares(v_ars_res, nP);
  BOOST_TEST(sum == secret + constant);
}

BOOST_DATA_TEST_CASE(const_addition, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), secret_val, const_val, idx) {
  size_t nP = 6;
  Field secret = Field(secret_val);
  Field constant = Field(const_val);
  auto v_ars = generateAuthRepShares(secret, nP);

  std::vector<RepShare<Field>> v_ars_res(nP);
  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_ars_res[i] = v_ars[i] * constant;
  }

  auto product = reconstructAuthRepShares(v_ars_res, nP);
  BOOST_TEST(product == secret * constant);  
}


BOOST_AUTO_TEST_SUITE_END()
