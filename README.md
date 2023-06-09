# Physical Unit Checking for C++

Adds an awareness of Physical (SI) units to C++ codes using custom data
structures. Units checking is strictly enforced - a Length cannnot be added to
a Time as this is not physically meaningful.

The UnitCheckedType class enforces unit checking (currently for 3 units, length,
time and mass) and the inner storage class holds numeric data.

This is a WORK IN PROGRESS and is not ready for use - currently it
is at proof-of-concept stage.
In particular:
* Only 3 units are implemented
* Many constructors are missing
* Nearly all operators are missing
* Error checking is minimal
* Only scalar, vector and tensor storage is available
** These are incomplete - missing most constructors, operators etc

However, using C++20 facilities this can enforce units checking at compile time
for arbirary rational exponents, and for any data-storage inner class which exposes
sufficient functionality.

Created June 2023 by Warwick RSE
## How to use this code
EITHER copy the 5 header files from the include directory into an appropriate place, or copy and unpack the tarball and add to your build. 

The PhysicalTypes.hpp header provides common Physical named units, such as Position (a double precision 3-element vector with units of m) or Force (a double precision 3-element vector with units of kgms^-2)
After including the header you can also instantiate a type of arbitrary units; create a UnitCheckedType\<L, M, T, Storage\> where L is the exponent of length (m), M of mass (kg) and T of time (s). Storage is the data storage type, which can be one of the built ins, scalar (dblscalar), 3-component or 4-component vector (dbl3vec or dbl4vec) or 3x3 or 4x4 tensor (dbl3tens or dbl4tens) all double precision; or can be ay other suitable storage type (see below). 

The code supports either Integral exponents for units (m^2 but not m^(1/2)) for C++17 or newer. Fractional exponent support requires C++20 and is enabled by setting -DUSE\_FRACTIONAL\_POWERS at the compile step. 

The provided src/test.cpp and Makefile builds an example code showing how to use the code and demonstrating some of the features.


### Storage types
A storage type is a class handling data storage, which must implement any function that you want the resulting UnitChecked type to expose to a final user. Any other function can be omitted - they will not be available to use, but code that does not use them will compile and run fine. _HOWEVER_ it may give better/ more instructive compile errors to implement invalid functions and use something like a static\_assert or =delete to force non-compilation. 

Functions that the UnitChecked type will wrap if available:
* Any of the following constructors
  ** Default constructor
  ** One element constructor (all values set)
  ** Initializer list constructor handling any depth of nesting required
* Copy, move etc constructors
* 'get' function taking as many arguments as wanted, returning a reference to value, and const variant returning copy
* Optional [] operator if 1D access makes sense and is desired
* << stream operator
* Standard arithmetic operators for self binary ops
* Comparison ops as desired
* Optional arithmetic operators for scalar ops (e.g. multiply with double or int)
* Heterogenous ops for any desired interactions (such as scalar-vector or vector-tensor multiplication which we have implemented)


