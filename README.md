# Physical Unit Checking for C++

Adds an awareness of Physical (SI) units to C++ codes using custom data
structures.
Units checking is strictly enforced at compile time - a Length cannnot be added to
a Time as this is not physically meaningful. 
To achieve this, every physical quantity has a different type, and the inter-operation of these is prescribed - for example a Length divided by a Time becomes a Velocity. 
Two Velocities can be added together, or Velocity can be divided by Time again to yield Acceleration etc. 
Using C++20 facilities this can enforce units checking at compile time for arbitrary rational exponents. See below for details.

Thus, _any physical equation which compiles is shown to be dimensionally correct_. See the test code for a host of examples.

The UnitCheckedType class enforces unit checking on all operations and the inner storage class holds numeric data. 
Numeric operations and some functions such as dot product are avalable. 
In general, any class can be used as the inner data-storage class, as long as it exposes sufficient function - see below for details. 

Notes:
* Some operations may be missing 
* Only scalar, vector and tensor storage is implemented
    * These may be incomplete


Created June 2023 by Warwick RSE

Version 1.0 (Only Length Mass and Time) July 2023

Version 2.0 (All 7 SI base units) February 2024
## How to use this code
You need the 5 header files to be included/accessible in your build, and can then include either PhysicalTypes.hpp to get common Physical units and access to UnitCheckedType to construct your own, or include UnitCheckedType.hpp directly if you do not want the extra defines.
You will need at least C++17 support - C++20 for full power (see section below about fractional exponent support).

The PhysicalTypes.hpp header provides common Physical named units, such as Position (a double precision 3-element vector with units of m) or Force (a double precision 3-element vector with units of kgms<sup>-2</sup>)
After including the header you can also instantiate a type of arbitrary units. Syntactic sugar is given for "dynamical" (Length, Mass, Time), "thermodynamical" (Length, Mass, Time, Temperature) and "electrodynamical" (Length, Mass, Time, Temperature, Current) types to avoid having to provide all 7 units directly. 
So for example we can create a UnitCheckedType\<L, M, T, Storage\> where
L is the exponent of length (m), M of mass (kg) and T of time (s). Storage is the data storage type, which can be one of the built ins, scalar (dblscalar), 3-component or 4-component vector (dbl3vec or dbl4vec) or 3x3 or 4x4 tensor (dbl3tens or dbl4tens) all double precision; or can be ay other suitable storage type (see below). 
For more general quantities, the complete order (along with the names used in the code) is
* L Length (m)
* M Mass (kg)
* T Time (s)
* K Temperature (K)
* A Current (A)
* MO Amount (mole)
* CD Luminous intensity (cd)

A UnitCheckedType (name for backwards compatibility) or UnitCheckedTypeDynamic exposes the first 3, UnitCheckedTypeThermodynamic the first 4, and UnitCheckedTypeElectrodynamic the first 5. UnitCheckedTypeFull exposes all 7 of these.

### Fractional exponents
C++20 introduces support which allows us to cleanly implement fractional exponents, such as m<sup>1/2</sup>. While we believe this could be done pre-20, we do not support it - instead we roll over to Integer exponents only. To enable fractional exponent support with C++20 supporting compilers, set the define -DUSE\_FRACTIONAL\_POWERS at the compile step.

Fractional powers mostly arise due to square-rooting etc, but can be instantiated directly if desired by creating an SF (a simple fraction type) type with 2 parameters, the numerator and the denominator, such as SF{1,2} for 1/2, and creating a UnitChecked type with this as the relevant parameter. See example code for an example.

Note: Fraction equality is strict - "1/2" is not the same entity as "2/4". As long as you never instantiate an unsimplified fraction, they will never arise from arithmetic, so this should not matter. 

### Example code
The provided src/test.cpp and Makefile builds an example code showing how to use the code and demonstrating many of the features. Build with `make` and run with `./Test`. A _failing_ example build is also given, which shows compile time errors such as trying to add incompatible units. Build this with `make clean && make FAIL=true` - note it _will not compile_.



### Storage types
A storage type is a class handling data storage, which must implement any function that you want the resulting UnitChecked type to expose to a final user. Any other function can be omitted - they will not be available to use, but code that does not use them will compile and run fine. _HOWEVER_ it may give better/ more instructive compile errors to implement invalid functions and use something like a static\_assert or =delete to force non-compilation. 

Functions that the UnitChecked type will wrap if available:
* Any of the following constructors
    * Default constructor
    * One element constructor (all values set)
    * Initializer list constructor handling up to 3 layers of nesting
* Copy, move etc constructors
* 'get' function taking as many arguments as wanted, returning a reference to value, and const variant returning copy (Note: this is exposed as unsafeGet, as it effectively casts away units)
* Optional [] operator if 1D access makes sense and is desired
* << stream operator
* Standard arithmetic operators for self binary ops
* Comparison ops as desired
* Optional arithmetic operators for scalar ops (e.g. multiply with double or int)
* Heterogenous ops for any desired interactions (such as scalar-vector or vector-tensor multiplication which we have implemented)

## But why the name?
Because comparing a distance to a time is like comparing apples to oranges, and adding one of each together doesn't get you two, it gets you an apple and an orange. My PhD supervisor loved this phrase, and it's the first thing I think of whenever I am checking units make sense or are compatible.

## Which operations make sense?
In theory, the unit-checked type could wrap arbitrary functions on the underlying data. However, a lot of operations can only be applied to a unit-less value - for example, sin, exp etc can only be given a dimensionless argument. For all of these, one should combine unit-ed quantities until a dimensionless value is reached. Then, either extract the value with .get if there is an overload for the underlying storage type, or use a cast to a suitable numeric type if there is one.
