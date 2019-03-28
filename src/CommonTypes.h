#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <array>

static constexpr int X_DIM = 0;
static constexpr int Y_DIM = 1;
static constexpr int Z_DIM = 2;
static constexpr int DIM_  = 3;

using Real  = double;
using Real3 = std::array<Real, DIM_>;

using Box    = std::array<Real3, DIM_>;
using Matrix = Box;

using Rvec = double[DIM_];

using Int3  = std::array<int, DIM_>;

using Range  = std::array<Real, 2>;
using Range3 = std::array<Range, DIM_>;

#endif /* COMMON_TYPES_H */
