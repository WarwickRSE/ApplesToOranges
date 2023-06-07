## Physical Unit Checking for C++

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

