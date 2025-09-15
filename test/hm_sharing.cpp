#define BOOST_TEST_MODULE hm_sharing

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include "types.h"
#include "hm_sharing.h"

using namespace hmAsyncAsterisk;
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

fieldPoly randomPolynomial(long degree, const Field& constant) {
  fieldPoly poly;
  poly.SetLength(degree + 1); // Set the length of the polynomial

  std::random_device rd;
  std::mt19937 engine(rd());
  ZZ p = conv<ZZ>("17816577890427308801"); // Example prime number
  std::uniform_int_distribution<uint64_t> distrib(0, conv<uint64_t>(p) - 1);

  poly[0] = constant;
  for (long i = 1; i <= degree; ++i) {
      poly[i] = Field(distrib(engine));
  }

  return poly;
}

std::vector<AuthShamirShare<Field>> generateAuthShamirShares(Field secret, size_t nP) {
  Field MAC_key = Field(5);
  std::random_device rd;
  std::mt19937 engine(rd());
  std::uniform_int_distribution<uint64_t> distrib; 
  
  Field tag = secret * MAC_key;

  long degree = (nP-1)/2;
  fieldPoly poly_key = randomPolynomial(degree, MAC_key);
  fieldPoly poly_val = randomPolynomial(degree, secret);
  fieldPoly poly_tag = randomPolynomial(degree, tag);

  std::vector<Field> key_shares(nP);
  std::vector<Field> values(nP);
  std::vector<Field> tags(nP);  

  for (int i = 0; i < nP; ++i) {
    key_shares[i] = eval(poly_key, Field(i+1));
    values[i] = eval(poly_val, Field(i+1));
    tags[i] = eval(poly_tag, Field(i+1));
  }

  std::vector<AuthShamirShare<Field>> AAS;
  for(size_t i = 0; i < nP ; i++) {
    AuthShamirShare<Field> temp(key_shares[i], values[i], tags[i]);
    AAS.push_back(temp);
  }

  return AAS;
}

Field reconstructAuthShamirShares(const std::vector<AuthShamirShare<Field>>& v_aas, size_t nP) {
  long degree = (nP-1)/2;
  std::vector<Field> evalPoints, key_vec, val_vec, tag_vec;
  for(size_t i = 1; i <= degree+1; i++) {
    evalPoints.push_back(Field(i));
    key_vec.push_back(v_aas[i-1].keySh());
    val_vec.push_back(v_aas[i-1].valueAt());
    tag_vec.push_back(v_aas[i-1].tagAt());
  }
  fieldPoly poly_key = reconstructPolynomial(evalPoints, key_vec);
  fieldPoly poly_val = reconstructPolynomial(evalPoints, val_vec);
  fieldPoly poly_tag = reconstructPolynomial(evalPoints, tag_vec);
  
  Field key = coeff(poly_key, 0);
  Field secret = coeff(poly_val, 0);
  Field tag = coeff(poly_tag, 0);

  if(secret * key == tag) { return secret; }
  else {
    std::cout<< "Incorrect sharing !!!" << std::endl;
    return Field(0);
  }
}

BOOST_AUTO_TEST_SUITE(authenticated_Shamir_sharing)

BOOST_DATA_TEST_CASE(reconstruction, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), secret_val, idx) {
  size_t nP = 4;
  Field secret = Field(secret_val);

  auto v_aas = generateAuthShamirShares(secret, nP);
  
  auto recon_value = reconstructAuthShamirShares(v_aas, nP);
    
  BOOST_TEST(recon_value == secret);
}

BOOST_DATA_TEST_CASE(share_arithmetic, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), vala, valb, idx) {
  size_t nP = 4;
  Field a = Field(vala);
  Field b = Field(valb);
  auto v_aas_a = generateAuthShamirShares(a, nP);
  auto v_aas_b = generateAuthShamirShares(b, nP);

  std::vector<AuthShamirShare<Field>> v_aas_c(nP);

  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_aas_c[i] = v_aas_a[i] + v_aas_b[i];
  }

  auto sum = reconstructAuthShamirShares(v_aas_c, nP);

  BOOST_TEST(sum == a + b);
  
  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_aas_c[i] = v_aas_a[i] - v_aas_b[i];
  }

  auto difference = reconstructAuthShamirShares(v_aas_c, nP);
  BOOST_TEST(difference == a - b);  
}

BOOST_DATA_TEST_CASE(share_const_arithmetic, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), secret_val, const_val, idx) {
  size_t nP = 6;
  Field secret = Field(secret_val);
  Field constant = Field(const_val);
  auto v_aas = generateAuthShamirShares(secret, nP);

  std::vector<AuthShamirShare<Field>> v_aas_res(nP);
  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_aas_res[i] = v_aas[i].add(constant);
  }

  auto sum = reconstructAuthShamirShares(v_aas_res, nP);
  BOOST_TEST(sum == secret + constant);
}

BOOST_DATA_TEST_CASE(const_addition, bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::xrange(NUM_SAMPLES), secret_val, const_val, idx) {
  size_t nP = 6;
  Field secret = Field(secret_val);
  Field constant = Field(const_val);
  auto v_aas = generateAuthShamirShares(secret, nP);

  std::vector<AuthShamirShare<Field>> v_aas_res(nP);
  for (size_t i = 0; i < nP; ++i) {
    // This implicitly checks compound assignment operators too.
    v_aas_res[i] = v_aas[i] * constant;
  }

  auto product = reconstructAuthShamirShares(v_aas_res, nP);
  BOOST_TEST(product == secret * constant);  
}


BOOST_AUTO_TEST_SUITE_END()
