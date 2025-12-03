#include "helpers.h"

namespace utils {
    void randomizeZZp(PRG& prg, Field& val, int nbytes) {
        uint64_t var;
        prg.random_data(&var, nbytes);
        val = Field(var);
    }

    /*fieldPoly reconstructPolynomial(const std::vector<Field>& x, const std::vector<Field>& y) {
        int t = x.size() - 1;
        fieldPoly result;
        for (int i = 0; i <= t; i++) {
            fieldPoly term(1); // Initialize term as 1
            Field denom = Field(1);    
            for (int j = 0; j <= t; j++) {
                if (i != j) {
                    term *= (fieldPoly(1, 1) - x[j]); // term *= (x - x[j])
                    denom *= (x[i] - x[j]); // denom *= (x[i] - x[j])
                }
            }
            result += term * (y[i] / denom); // Add the term to the result
        }
        return result;
    }*/

    fieldPoly reconstructPolynomial(const std::vector<Field>& x, const std::vector<Field>& y) {
        Vec<Field> x_vec, y_vec;
        x_vec.SetLength(x.size());
        y_vec.SetLength(y.size());
        for (size_t i = 0; i < x.size(); ++i) {
            x_vec[i] = x[i];
            y_vec[i] = y[i];
        }
        // Use NTL's interpolate function
        fieldPoly result;
        interpolate(result, x_vec, y_vec);
        return result;
    }

    fieldPoly randomPolynomial(PRG& prg, long degree, const Field& constant) {
        fieldPoly poly;
        poly.SetLength(degree + 1);
        poly[0] = constant;
        for (long i = 1; i <= degree; ++i) {
            randomizeZZp(prg, poly[i], sizeof(Field));
        }
        return poly;
    }

uint64_t to_uint64(Field a) {
    ZZ z = a.LoopHole();  
    return conv<uint64_t>(z); 
}

fieldDig hashFields(const std::vector<Field>& fields) {
    Hash hash;

    std::vector<uint64_t> input_vec(fields.size());
    for (size_t i = 0; i < fields.size(); ++i)
        input_vec[i] = to_uint64(fields[i]);

    const char* data_ptr = reinterpret_cast<const char*>(input_vec.data());
    size_t byte_len = input_vec.size() * sizeof(uint64_t);
    hash.put(data_ptr, static_cast<int>(byte_len)); 
    
    uint8_t dig_bytes[Hash::DIGEST_SIZE]; 
    hash.digest(dig_bytes); 
    
    std::vector<Field> result(4);
    for (size_t i = 0; i < 4; ++i) {
        uint64_t val;
        std::memcpy(&val, dig_bytes + i * sizeof(uint64_t), sizeof(uint64_t));
        result[i] = Field(val);
    }

    return result;
}

void appendFieldDig(std::vector<Field>& vec, fieldDig& dig) {
        if (vec.size() < dig.size()) {
        throw std::runtime_error("vector too small to overwrite last 4 elements");
    }

    std::copy(dig.begin(), dig.end(), vec.end() - 4);
}

// -------------- Asterisk methods ----------------

int pidFromOffset(int id, int offset) {
  int pid = (id + offset) % 4;
  if (pid < 0) {
    pid += 4;
  }
  return pid;
}

int offsetFromPid(int id, int pid) {
  if (id < pid) {
    return pid - id;
  }

  return 4 + pid - id;
}

size_t upperTriangularToArray(size_t i, size_t j) {
  // (i, j) co-ordinate in upper triangular matrix (without diagonal) to array
  // index in column major order.
  auto mn = std::min(i, j);
  auto mx = std::max(i, j);
  auto idx = (mx * (mx - 1)) / 2 + mn;
  return idx;
}

std::vector<uint64_t> packBool(const bool* data, size_t len) {
  std::vector<uint64_t> res;
  for (size_t i = 0; i < len;) {
    uint64_t temp = 0;
    for (size_t j = 0; j < 64 && i < len; ++j, ++i) {
      if (data[i]) {
        temp |= (0x1ULL << j);
      }
    }
    res.push_back(temp);
  }

  return res;
}

void unpackBool(const std::vector<uint64_t>& packed, bool* data, size_t len) {
  for (size_t i = 0, count = 0; i < len; count++) {
    uint64_t temp = packed[count];
    for (int j = 63; j >= 0 && i < len; ++i, --j) {
      data[i] = (temp & 0x1) == 0x1;
      temp >>= 1;
    }
  }
}

// void randomizeZZpE(emp::PRG& prg, NTL::ZZ_pE& val) {
//   std::vector<Ring> coeff(NTL::ZZ_pE::degree());
//   prg.random_data(coeff.data(), sizeof(Ring) * coeff.size());

//   NTL::ZZ_pX temp;
//   temp.SetLength(NTL::ZZ_pE::degree());

//   for (size_t i = 0; i < coeff.size(); ++i) {
//     temp[i] = coeff[i];
//   }

//   NTL::conv(val, temp);
// }

// void randomizeZZpE(emp::PRG& prg, NTL::ZZ_pE& val, Ring rval) {
//   std::vector<Ring> coeff(NTL::ZZ_pE::degree() - 1);
//   prg.random_data(coeff.data(), sizeof(Ring) * coeff.size());

//   NTL::ZZ_pX temp;
//   temp.SetLength(NTL::ZZ_pE::degree());

//   temp[0] = rval;
//   for (size_t i = 1; i < coeff.size(); ++i) {
//     temp[i] = coeff[i];
//   }

//   NTL::conv(val, temp);
// }

// void receiveZZpE(emp::NetIO* ios, NTL::ZZ_pE* data, size_t length) {
//   auto degree = NTL::ZZ_pE::degree();
//   // Assumes that every co-efficient of ZZ_pE is same range as Ring.
//   std::vector<uint8_t> serialized(sizeof(Ring));

//   NTL::ZZ_pX poly;
//   poly.SetLength(degree);
//   for (size_t i = 0; i < length; ++i) {
//     for (size_t d = 0; d < degree; ++d) {
//       ios->recv_data(serialized.data(), serialized.size());
//       auto coeff = NTL::conv<NTL::ZZ_p>(
//           NTL::ZZFromBytes(serialized.data(), serialized.size()));
//       poly[d] = coeff;
//     }
//     NTL::conv(data[i], poly);
//   }
// }

// void sendZZpE(emp::NetIO* ios, const NTL::ZZ_pE* data, size_t length) {
//   auto degree = NTL::ZZ_pE::degree();
//   // Assumes that every co-efficient of ZZ_pE is same range as Ring.
//   std::vector<uint8_t> serialized(sizeof(Ring));

//   for (size_t i = 0; i < length; ++i) {
//     const auto& poly = NTL::rep(data[i]);
//     for (size_t d = 0; d < degree; ++d) {
//       const auto& coeff = NTL::rep(NTL::coeff(poly, d));
//       NTL::BytesFromZZ(serialized.data(), coeff, serialized.size());
//       ios->send_data(serialized.data(), serialized.size());
//     }
//   }
//   ios->flush();
// }

} // namespace utils