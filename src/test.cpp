
#include <iostream>
#include "UnitCheckedType.hpp"

int main(){

  // Create three related physical quantities

  // A time
  UnitCheckedType<SF{0,1}, SF{0,1}, SF{1,1}, dblscalar> t{0.1};
  
  std::cout<<"Defined time t= "<<t<<t.units()<<std::endl;
  // A position
  UnitCheckedType<SF{1,1}, SF{0,1}, SF{0,1}, dbl3vec> x{1.0, 2.0, 3.0};  
  std::cout<<"Defined position x= "<<x<<x.units()<<std::endl;

  // A velocity
  UnitCheckedType<SF{1,1}, SF{0,1}, SF{-1,1}, dbl3vec> v{1.0, 2.0, 3.0};
  std::cout<<"Defined velocity v= "<<v<<v.units()<<std::endl;

  // Update position using x = x_0 + v * t
  auto x_new = x + v * t;
  std::cout<<"Position at t ="<<t<<t.units()<<" is "<<x_new<<x_new.units()<<std::endl;



  return 0;
};