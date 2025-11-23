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

} // namespace utils