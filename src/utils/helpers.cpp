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
} // namespace utils