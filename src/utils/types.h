#pragma once

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

using namespace NTL;

namespace utils {
    using Field = ZZ_p;
    using fieldPoly = ZZ_pX;
    constexpr uint64_t FIELDSIZE = 8; // bytes
} // namespace utils