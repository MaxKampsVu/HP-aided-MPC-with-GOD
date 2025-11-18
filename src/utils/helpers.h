#pragma once

#include <emp-tool/emp-tool.h>
#include "types.h"

using namespace emp;
using namespace NTL;

namespace utils {
    void randomizeZZp(PRG& prg, Field& val, int nbytes);
    fieldPoly reconstructPolynomial(const std::vector<Field>& x, const std::vector<Field>& y);
    fieldPoly randomPolynomial(PRG& prg, long degree, const Field& constant);
    std::vector<Field> hashFields(const std::vector<Field>& fields);
} // namespace utils