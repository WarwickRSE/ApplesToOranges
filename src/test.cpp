
#include <iostream>
#include "PhysicalTypes.hpp"

int main(){

  // Create three related physical quantities

  // A time
  Time t{0.1};
  std::cout<<"Defined time t= "<<t<<t.units()<<std::endl;
  
  // A position
  Position x{1.0, 2.0, 3.0};  
  std::cout<<"Defined position x= "<<x<<x.units()<<std::endl;

  // A velocity
  Velocity v{1.0, 2.0, 3.0};
  std::cout<<"Defined velocity v= "<<v<<v.units()<<std::endl;

  // Update position using x = x_0 + v * t
  x = x + v * t;
  std::cout<<"Position at t ="<<t<<t.units()<<" is "<<x<<x.units()<<std::endl;

  Length l{1.0};
  // Demo of comparison operators
  std::cout<<"l > 0? "<< (l > Length{0.0})<<std::endl;
  std::cout<<"x == x? "<< (x == x)<<std::endl;

  std::cout<< "x > l? "<< (x.magnitude() > l.magnitude())<<std::endl;
  std::cout<< "x > l? "<< (x > l)<<std::endl;

  // Try gridded operations and show how to initialise
  UCScalar Arr[3][3];
  Arr[0][0] = UCScalar{1.0, 2.0, 3.0};
  Arr[1][1] = UCScalar{1.0, 2.0, 3.0};
  Arr[2][2] = UCScalar{1.0, 2.0, 3.0};

  auto Arr2 = Arr;
  auto Arr3 = Arr2;
  std::cout<<"Arr3[0][0] = "<<Arr3[0][0]<<std::endl;

#ifdef USE_FRACTIONAL_POWERS

  // Creating your own custom physical type to use elsewhere
  using UnsquareM = UnitCheckedType<SF{1,2}, 0, 0, dblscalar>;
  UnsquareM lu{0.3};
  std::cout<<" lu = "<<lu<<lu.units()<<std::endl;
  std::cout<<" l+ lu^2 = "<<l + lu*lu<<l.units()<<std::endl;

#endif

  return 0;
};