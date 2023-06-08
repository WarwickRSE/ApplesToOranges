
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


  // Try gridded operations and show how to initialise
  UCScalar Arr[3][3];
  Arr[0][0] = UCScalar{1.0, 2.0, 3.0};
  Arr[1][1] = UCScalar{1.0, 2.0, 3.0};
  Arr[2][2] = UCScalar{1.0, 2.0, 3.0};

  auto Arr2 = Arr;
  auto Arr3 = Arr2;

  return 0;
};