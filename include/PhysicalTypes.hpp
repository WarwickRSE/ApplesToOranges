#ifndef __PHYSICALTYPES_HPP__
#define __PHYSICALTYPES_HPP__

#include "UnitCheckedType.hpp"

// For clarity
using dblscalar = STScalar<double>;
using dbl3vec = STVector<double, 3>;
using dbl4vec = STVector<double, 4>;
using dbl3tens = STTensor<double, 3>;
using dbl4tens = STTensor<double, 4>;

const dblscalar sc_ident = dblscalar::identity();
const dbl3vec vec3_ident = dbl3vec::identity();
const dbl3tens tens3_ident = dbl3tens::identity();

// No physical units
using UCScalar = UnitCheckedType<0, 0, 0, dblscalar>;
using UCVector = UnitCheckedType<0, 0, 0, dbl3vec>;
using UCTensor = UnitCheckedType<0, 0, 0, dbl3tens>;

// Single units, scalar and (where sensible), 3-vector

using Length = UnitCheckedType<1, 0, 0, dblscalar>;
using Mass = UnitCheckedType<0, 1, 0, dblscalar>;
using Time = UnitCheckedType<0, 0, 1, dblscalar>;
using Frequency = UnitCheckedType<0, 0, -1, dblscalar>;
using Position = UnitCheckedType<1, 0, 0, dbl3vec>;

// Common physical quantities for dynamics
using Speed = UnitCheckedType<1, 0, -1, dblscalar>;
using Velocity = UnitCheckedType<1, 0, -1, dbl3vec>;
using Momentum = UnitCheckedType<1, 1, -1, dbl3vec>;

using Acceleration = UnitCheckedType<1, 0, -2, dbl3vec>;
using Force = UnitCheckedType<1, 1, -2, dbl3vec>;
using Pressure = UnitCheckedType<-1, 1, -2, dblscalar>;

using Energy = UnitCheckedType<2, 1, -2, dblscalar>;
using Power = UnitCheckedType<2, 1, -3, dblscalar>;

// Common Tensor types
using Stress = UnitCheckedType<-1, 1, -2, dblscalar>;
using StressVector = UnitCheckedType<-1, 1, -2, dbl3vec>;
using StressTensor = UnitCheckedType<-1, 1, -2, dbl3tens>;


// Other basic units
using Temperature = UnitCheckedTypeThermodynamic<0, 0, 0, 1, dblscalar>;
using Current =    UnitCheckedTypeElectrodynamic<0, 0, 0, 0, 1, dblscalar>;
using AmountOfSubstance =    UnitCheckedTypeFull<0, 0, 0, 0, 0, 1, 0, dblscalar>;
using LuminousIntensity =    UnitCheckedTypeFull<0, 0, 0, 0, 0, 0, 1, dblscalar>;

// Some more common derived units - most of the rest of the International System of Units 22 derived types
using Charge =     UnitCheckedTypeElectrodynamic<0, 0, 1, 0, 1, dblscalar>;
using Voltage =    UnitCheckedTypeElectrodynamic<2, 1, -3, 0, -1, dblscalar>;
using Capacitance =UnitCheckedTypeElectrodynamic<-2, -1, 4, 0, 2, dblscalar>;
using Resistance = UnitCheckedTypeElectrodynamic<2, 1, -3, 0, -2, dblscalar>;
using Inductance = UnitCheckedTypeElectrodynamic<2, 1, -2, 0, -2, dblscalar>;
// Tesla:
using MagneticFluxDensity = UnitCheckedTypeElectrodynamic<0, 1, -2, 0, -1, dblscalar>;

using Illuminance = UnitCheckedTypeFull<-2, 0, 0, 0, 0, 0, 1, dblscalar>;

//For convenience, a function to set nice unit strings for commonest SI compounds
inline void registerSIUnitStrings(){
    UnitChecked::registerUnits<Force>("N");
    UnitChecked::registerUnits<Pressure>("Pa");
    UnitChecked::registerUnits<Power>("W");
    UnitChecked::registerUnits<Energy>("J");
    UnitChecked::registerUnits<Charge>("C");
    UnitChecked::registerUnits<Voltage>("V");
    UnitChecked::registerUnits<Capacitance>("F");
    UnitChecked::registerUnits<Resistance>("Ohm");
}

#endif
