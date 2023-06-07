
#include <iostream>
#include "UnitCheckedType.hpp"

int main(){

  SimpleFrac a{1,2}, b{1,2};
  std::cout<< (a+b).num << "/" << (a+b).denom << std::endl;
  std::cout<< simplifiable(a+b)<<std::endl;
  std::cout<< (a==b)<<std::endl;


  ABC c{1.0,2.0,3.0}, d{2.0,3.0,4.0};

  auto e = c*d;

  std::cout<<"thrid val of first is "<<c[2]<<std::endl;

  std::cout<< "first is "<<c<<c.units()<<std::endl;
  std::cout<< "second is "<<d<<d.units()<<std::endl;
  
  std::cout<<"first*second is "<<e<<e.units()<<std::endl;

  ABCs f{7.0};
  auto g = f*d;
  
  std::cout<< "thrid is "<<f<<f.units()<<std::endl;
  std::cout<<"thrid*first is "<<g<<g.units()<<std::endl;


  ABCt h{11.0};
  std::cout<<"testing gets "<<c.get(0)<<" "<<f.get()<<" "<<h.get(0, 0)<<std::endl;

  ABCt jj{{1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}};
  std::cout<<"testing gets "<<jj.get(0, 0)<<" "<<jj.get(1, 1)<<" "<<jj.get(2, 2)<<std::endl;
  return 0;
};