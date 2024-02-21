# Physical Unit Checking for C++

Adds an awareness of Physical (SI) units to C++ codes using custom data
structures. Units checking is strictly enforced at compile time - a Length cannnot be added to
a Time as this is not physically meaningful.

The UnitCheckedType class enforces unit checking (currently for 3 units, length,
time and mass) and the inner storage class holds numeric data.

This is a WORK IN PROGRESS but can be used in real codes with some caveats:
In particular:
* Only 3 units are implemented (Length, Time, Mass)
* Some constructors may be missing
* Some operators may be missing
* Error checking is minimal
* Only scalar, vector and tensor storage is implemented
    * These may be incomplete
In short, basic physical equations can definitely be coded and checked with the current state, but there may be bugs or things may fail to compile which should be OK.

However, using C++20 facilities this can enforce units checking at compile time
for arbirary rational exponents, and for any data-storage inner class which exposes
sufficient functionality.

Created June 2023 by Warwick RSE
Version 1.0 July 2023
## How to use this code
You need the 5 header files to be included/accessible in your build, and can then include either PhysicalTypes.hpp to get common Physical units and access to UnitCheckedType to construct your own, or include UnitCheckedType.hpp directly if you do not want the extra defines.
You will need at least C++17 support - C++20 for full power (see section below about fractional exponent support).

The PhysicalTypes.hpp header provides common Physical named units, such as Position (a double precision 3-element vector with units of m) or Force (a double precision 3-element vector with units of kgms<sup>-2</sup>)
After including the header you can also instantiate a type of arbitrary units; create a UnitCheckedType\<L, M, T, Storage\>. L is the exponent of length (m), M of mass (kg) and T of time (s). Storage is the data storage type, which can be one of the built ins, scalar (dblscalar), 3-component or 4-component vector (dbl3vec or dbl4vec) or 3x3 or 4x4 tensor (dbl3tens or dbl4tens) all double precision; or can be ay other suitable storage type (see below). 

### Fractional exponents
C++20 introduces support which allows us to cleanly implement fractional exponents, such as m<sup>1/2</sup>. While we believe this could be done pre-20, we do not support it - instead we roll over to Integer exponents only. To enable fractional exponent support with C++20 supporting compilers, set the define -DUSE\_FRACTIONAL\_POWERS at the compile step. 

Fractional powers mostly arise due to square-rooting etc, but can be instantiated directly if desired by creating an SF type with 2 parameter, the numerator and the denominator, such as SF{1,2} for 1/2, and creating a UnitChecked type with this as the relevant parameter. See example code for an example.

Note: Fraction equality is strict - 1/2 is not the same entity as 2/4. As long as you never instantiate an unsimplified fraction, they will never arise from arithmetic, so this should not matter. 

### Example code
The provided src/test.cpp and Makefile builds an example code showing how to use the code and demonstrating some of the features.


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
