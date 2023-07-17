
#include <iostream>
#include <vector>
#include "PhysicalTypes.hpp"

int main(){

#ifdef USE_FRACTIONAL_POWERS 
#ifdef DEBUG
  // Quick test of fractions
  SimpleFrac aa{1,2}, bb{2,4}, cc{1,3};
  std::cout<<"aa = "<<aa<<std::endl;
  std::cout<<"bb = "<<bb<<std::endl;
  std::cout<<"cc = "<<cc<<std::endl;
  std::cout<<"aa == bb? "<<(aa==bb)<<std::endl;
  std::cout<<"aa equal bb? "<<is_equal(aa, bb)<<std::endl;

  std::cout<<"aa < cc? "<<is_less(aa,cc)<<std::endl;
  std::cout<<"aa > cc? "<<is_greater(aa,cc)<<std::endl;
#endif
#endif
  // Create three related physical quantities

  // A time
  Time t{0.1};
  std::cout<<"Defined time t= "<<t<<t.units()<<std::endl;

  std::vector<Time> t_steps;
  t_steps.push_back(t);

  Time t2(t);

  // A position
  Position x{1.0, 2.0, 3.0};  
  std::cout<<"Defined position x= "<<x<<x.units()<<std::endl;

  // A velocity
  Velocity v{1.0, 2.0, 3.0};
  std::cout<<"Defined velocity v= "<<v<<v.units()<<std::endl;

  // Update position using x = x_0 + v * t
  x = x + v * t;
  std::cout<<"Position at t ="<<t<<t.units()<<" is "<<x<<x.units()<<std::endl;

  // Dot product
  std::cout<<"x.v="<<x.dot(v)<<std::endl;

  // Numeric multiply
  std::cout<<"2x= "<<2.0*x<<x.units()<<"=="<<x*2.0<<x.units()<<std::endl;
  std::cout<<"x/2= "<<x/2.0<<x.units()<<std::endl;

  Length l{1.0};
  // Demo of comparison operators
  std::cout<<"l > 0? "<< (l > Length{0.0})<<std::endl;
  std::cout<<"x == x? "<< (x == x)<<std::endl;

  std::cout<<"x.magnitude()= "<<x.magnitude()<<std::endl;
  std::cout<<"l.magnitude()= "<<l.magnitude()<<std::endl;

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

#ifdef DEBUG
  UnitCheckedType<0, 0, 0, STDummy> dummy{1.0};
#endif

  return 0;
};