#ifndef __PHYSICALTYPES_HPP__
#define __PHYSICALTYPES_HPP__

#include "UnitCheckedType.hpp"

// For clarity
using dblscalar = STScalar<double>;
using dbl3vec = STVector<double, 3>;
using dbl3tens = STTensor<double, 3>;

// No physical units
using UCScalar = UnitCheckedType<SF{0,1}, SF{0,1}, SF{0,1}, dblscalar>;
using UCVector = UnitCheckedType<SF{0,1}, SF{0,1}, SF{0,1}, dbl3vec>;

// Single units, scalar and (where sensible), 3-vector

using Length = UnitCheckedType<SF{1,1}, SF{0,1}, SF{0,1}, dblscalar>;
using Mass = UnitCheckedType<SF{0,1}, SF{1,1}, SF{0,1}, dblscalar>;
using Time = UnitCheckedType<SF{0,1}, SF{0,1}, SF{1,1}, dblscalar>;
using Position = UnitCheckedType<SF{1,1}, SF{0,1}, SF{0,1}, dbl3vec>;

// Common physical quantities for dynamics
using Speed = UnitCheckedType<SF{1,1}, SF{0,1}, SF{-1,1}, dblscalar>;
using Velocity = UnitCheckedType<SF{1,1}, SF{0,1}, SF{-1,1}, dbl3vec>;
using Momentum = UnitCheckedType<SF{1,1}, SF{1,1}, SF{-1,1}, dbl3vec>;

using Acceleration = UnitCheckedType<SF{1,1}, SF{0,1}, SF{-2,1}, dbl3vec>;
using Force = UnitCheckedType<SF{1,1}, SF{1,1}, SF{-2,1}, dbl3vec>;

using Energy = UnitCheckedType<SF{2,1}, SF{1,1}, SF{-2,1}, dblscalar>;

#endif