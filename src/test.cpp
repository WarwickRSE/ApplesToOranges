
#include <iostream>
#include "UnitCheckedType.hpp"

int main(){

  SimpleFrac a{1,2}, b{1,2};
  std::cout<< (a+b).num << "/" << (a+b).denom << std::endl;
  std::cout<< simplifiable(a+b)<<std::endl;
  std::cout<< (a==b)<<std::endl;


  UnitCheckedType<SF{1,1},SF{1,1},SF{1,1},1> c{1};
  UnitCheckedType<SF{1,2},SF{1,1},SF{1,1},1> d{1};

  auto e = c*d;

  std::cout<< "a is "<<c.units()<<std::endl;
  std::cout<< "b is "<<d.units()<<std::endl;
  
  std::cout<<"a*b is "<<e.units()<<std::endl;

  return 0;
};