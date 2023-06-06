
#include <iostream>
#include "UnitCheckedType.hpp"

int main(){

  SimpleFrac a{1,2}, b{1,2};
  std::cout<< (a+b).num << "/" << (a+b).denom << std::endl;
  std::cout<< simplifiable(a+b)<<std::endl;
  std::cout<< (a==b)<<std::endl;


  ABC c{1,2,3}, d{2,3,4};

  auto e = c*d;

  std::cout<< "a is "<<c<<c.units()<<std::endl;
  std::cout<< "b is "<<d<<d.units()<<std::endl;
  
  std::cout<<"a*b is "<<e<<e.units()<<std::endl;

  return 0;
};