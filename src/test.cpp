
#include <iostream>
#include <vector>
#include <chrono>
#include "PhysicalTypes.hpp"

void run_timer_test();
void run_addition_timer_test();

int main(){

#ifdef USE_FRACTIONAL_POWERS
#ifdef DEBUG

  std::cout<<"___________________________________________________________"<<std::endl;
  // Quick test of fractions
  SimpleFrac aa{1,2}, bb{2,4}, cc{1,3};
  std::cout<<"aa = "<<aa<<std::endl;
  std::cout<<"bb = "<<bb<<std::endl;
  std::cout<<"cc = "<<cc<<std::endl;
  std::cout<<"aa == bb? "<<(aa==bb)<<std::endl;
  std::cout<<"aa equal bb? "<<is_equal(aa, bb)<<std::endl;

  std::cout<<"aa < cc? "<<is_less(aa,cc)<<std::endl;
  std::cout<<"aa > cc? "<<is_greater(aa,cc)<<std::endl;

  std::cout<<"___________________________________________________________"<<std::endl;
#endif
#endif
  // Create three related physical quantities

  // A time
  Time t{0.1};
  std::cout<<"Defined time t= "<<t<<t.units()<<std::endl;

  auto OneOverT = 1.0/t;
  std::cout<<"Inverse of t is "<<OneOverT<<OneOverT.units()<<std::endl;

  auto OneOverPos = 1.0/Position{1.0, 1.0, 1.0};
  std::cout<<"Inverse of pos is "<<OneOverPos<<OneOverPos.units()<<std::endl;

  std::vector<Time> t_steps;
  t_steps.push_back(t);
  std::cout<<"A vector of type time, containing "<<t_steps.size()<<" elements"<<std::endl;

  Time t2(t);

  std::cout<<"Two times are the same type: "<<t.isSameType(t2)<<std::endl;
  std::cout<<"t and 1/t are not: "<<t.isSameType(OneOverT)<<std::endl;

  // A position
  Position x{1.0, 2.0, 3.0};
  std::cout<<"Defined position x= "<<x<<x.units()<<std::endl;

  // A velocity
  Velocity v{1.0, 2.0, 3.0};
  std::cout<<"Defined velocity v= "<<v<<v.units()<<std::endl;

  // Update position using x = x_0 + v * t
  x = x + v * t;
  std::cout<<"Position at t ="<<t<<t.units()<<" is "<<x<<x.units()<<std::endl;

#ifdef FAIL_DEMO
  // Invalid - trying to add a position to a velocity
  auto bad_add = x + v;
#endif

  // Construct a dimensionless value (still unit checked, but strictly all units 0)
  auto dimless = 2.0 * t/t;

  // Powers
  // NOTE: since exponentiation changes the units, we have to supply the exponent at compile time to get the compile time unit checking. Hence, this is a templated function. As a convenience, there is also a runtime version for DIMENSIONLESS types only
  auto tsq = pow<2>(t);
  std::cout<<"t^2= "<<tsq<<tsq.units()<<std::endl;
  auto tsqsqrt = sqrt(tsq);
  std::cout<<"sqrt(t^2)= "<<tsqsqrt<<tsqsqrt.units()<<std::endl;

#ifdef USE_FRACTIONAL_POWERS
  auto tsqrt2 = pow<SF{1,2}>(tsq);
  std::cout<<"         = "<<tsqrt2<<tsqrt2.units()<<std::endl;
#endif
  // Valid only with fractional powers enabled, as units are s^(1/2)
#if defined USE_FRACTIONAL_POWERS || defined FAIL_DEMO
  auto tsqrt = sqrt(t);
  std::cout<<"sqrt(t)= "<<tsqrt<<tsqrt.units()<<std::endl;
  auto tcbrt = cbrt(t);
  std::cout<<"cbrt(t)= "<<tcbrt<<tcbrt.units()<<std::endl;
#endif


#ifdef FAIL_DEMO
  // Not valid - t has units:
  auto tsq2 = pow(t, 2);
#endif
  // Valid - dimensionless type supplied
  std::cout<<"t/t is dimensionless, so can be raised to a runtime power: 2.0^2 = "<< pow(dimless, 2)<<" and 2.0^2.0 = "<<pow(dimless, 2.0)<<std::endl;



  // Dot product
  std::cout<<"x.v="<<x.dot(v)<<std::endl;

  Position x2{4.0, 5.0, 6.0};
  // Cross product
  std::cout<<"x cross x2="<<x.cross(x2)<<std::endl;

  //Normalize
  x2.normalize();
  std::cout << "x2 normalised is " << x2 << std::endl;
  std::cout << "magnitude of x2 normalised is " << x2.magnitude() << std::endl;

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

  //Using casting in comparison
  UCScalar l2{1.0};
  std::cout<<"Comparison of Scalar (no units) to double: "<<std::endl;
  std::cout<<"Construct comparable type, val > UCScalar{1.0}"<<(l2 > UCScalar{1.0})<<std::endl;
  std::cout<<"Cast UCScalar to double, static_cast<double>(val) > 1.0 "<<(static_cast<double>(l2) > 1.0)<<std::endl;

  // Using casting to construct
  double tmp2{l2};
  std::cout<<"Casting to double: "<< typeid(tmp2).name()<<std::endl;
#if !defined NO_NARROWING_CONVERSIONS || defined FAIL_DEMO
  float tmp1{l2};  //Disallowed by -DNO_NARROWING_CONVERSIONS
  std::cout<<"Narrowing to float: "<< typeid(tmp1).name()<<std::endl;
#endif

  // Using special functions on dimensionless values
  // For scalars, cast to, or store to, a double
  std::cout<<"sin(l/l_0): "<<sin(static_cast<double>(l/Length{1.0}))<<std::endl;
  // OR use specialisation (if provided) on the underlying value
  std::cout<<"sin(l/l_0): "<<mysinfunction((l/Length{1.0}).data())<<std::endl;
  // Storage type can provide operations on higher rank types, either elementwise or whatever
  std::cout<<"sin(x/x_0): "<<mysinfunction((x/Length{1.0}).data())<<std::endl;

  // Accessor functions
  std::cout<<"t's value, ignoring units, is "<<t.unsafeGet()<<" and can be added to a simple double like this:"<<1.0 + t.unsafeGet()<< std::endl;
  auto scal = UCScalar{7.5};
  std::cout<<"Type without units supports simple get like this: "<<scal.get()<<std::endl;


  // Try gridded operations and show how to initialise
  UCScalar Arr[3][3];
  Arr[0][0] = UCScalar{1.0, 2.0, 3.0};
  Arr[1][1] = UCScalar{1.0, 2.0, 3.0};
  Arr[2][2] = UCScalar{1.0, 2.0, 3.0};

  auto Arr2 = Arr;

#ifdef USE_FRACTIONAL_POWERS

  // Creating your own custom physical type to use elsewhere
  using UnsquareM = UnitCheckedType<SF{1,2}, 0, 0, dblscalar>;
  UnsquareM lu{0.3};
  std::cout<<"A quantity with fractional exponent units: lu = "<<lu<<lu.units()<<std::endl;
  std::cout<<" l+ lu^2 = "<<l + lu*lu<<l.units()<<std::endl;

#endif

#ifdef DEBUG
  // Debug check making sure that we can instantiate a type without any methods, and thus do not force any methods to be defined
  UnitCheckedType<0, 0, 0, STDummy> dummy{1.0};
#endif

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Performance tests"<<std::endl;
  UCScalar Arr3[100];
  auto sz = sizeof(Arr3)/sizeof(double);
  if(sz == 100){
    std::cout<<"There is no storage overhead, a 100 element array of UCScalar is the size of "<<sz<<" doubles"<<std::endl;
  }else{
    std::cout<<"There is storage overhead, a 100 element array of UCScalar is the size of "<<sz<<" doubles"<<std::endl;
  }

#ifdef RUN_TIMER_TEST
  run_timer_test();
  run_addition_timer_test();
#endif

  return 0;
};


class timer{
    private:
        std::chrono::time_point<std::chrono::steady_clock> start, stop;
    public:
        void start_timer(){
            start = std::chrono::steady_clock::now();
        }
        void stop_timer(){
            stop = std::chrono::steady_clock::now();
        }
        long time_taken(){
            return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        }
        void print_time(){
            std::cout<<"Time taken: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()<<" ms"<<std::endl;
        }
};

void run_timer_test(){

  // Stack arrays first, keep to about 1000 elements in 1D

  std::cout<<"Checking stack array performance -------------------"<<std::endl;
  const long sz = 1000;
  UCScalar Arr[sz];
  long iters = 10000000;
  timer t;

  for(int j = 0; j < sz; j++){
    Arr[j] = UCScalar{1.0};
  }

  t.start_timer();
  for (long i = 0; i < iters; i++){
    for(int j = 0; j < sz; j++){
      Arr[j] += UCScalar{1.0};
    }
  }
  t.stop_timer();
  t.print_time();
  auto first_time = t.time_taken();
  // Use Arr to avoid optimiser removing entirely
  std::cout<<Arr[10]<<std::endl;

  double Arr2[sz];
  for(int j = 0; j < sz; j++){
    Arr2[j] = 1.0;
  }
  t.start_timer();
  for (long i = 0; i < iters; i++){
    for(int j = 0; j < sz; j++){
      Arr2[j] += 1.0;
    }
  }
  t.stop_timer();
  t.print_time();
  auto second_time = t.time_taken();
  std::cout<<Arr2[10]<<std::endl;

  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"         Time overhead for unit checked use: "<<(((double)first_time / (double)second_time)-1)*100.0<<" %"<<std::endl;

  std::cout<<"Checking heap array performance -------------------"<<std::endl;
  const long hsz = 1000000;
  UCScalar * hArr = new UCScalar[hsz];
  long hiters = 10000;

  t.start_timer();
  for (long i = 0; i < hiters; i++){
    for(int j = 0; j < hsz; j++){
      hArr[j] += UCScalar{1.0};
    }
  }
  t.stop_timer();
  t.print_time();
  first_time = t.time_taken();
  // Use Arr to avoid optimiser removing entirely
  std::cout<<hArr[10]<<std::endl;

  double* hArr2 = new double[hsz];
  t.start_timer();
  for (long i = 0; i < hiters; i++){
    for(int j = 0; j < hsz; j++){
      hArr2[j] += 1.0;
    }
  }
  t.stop_timer();
  t.print_time();
  second_time = t.time_taken();
  std::cout<<hArr2[10]<<std::endl;

  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"         Time overhead for unit checked use: "<<(((double)first_time / (double)second_time)-1)*100.0<<" %"<<std::endl;

};

void run_addition_timer_test(){
  // Compare the time for += and A = A + b

  std::cout<<"Checking relative performance of + and += -------------------"<<std::endl;
  const long sz = 1000;
  UCScalar Arr[sz];
  long iters = 10000000;
  timer t;

  for(int j = 0; j < sz; j++){
    Arr[j] = UCScalar{1.0};
  }

  t.start_timer();
  for (long i = 0; i < iters; i++){
    for(int j = 0; j < sz; j++){
      Arr[j] += UCScalar{1.0};
    }
  }
  t.stop_timer();
  t.print_time();
  auto first_time = t.time_taken();
  // Use Arr to avoid optimiser removing entirely
  std::cout<<Arr[10]<<std::endl;

  t.start_timer();
  for (long i = 0; i < iters; i++){
    for(int j = 0; j < sz; j++){
      Arr[j] = Arr[j] + UCScalar{1.0};
    }
  }
  t.stop_timer();
  t.print_time();
  auto second_time = t.time_taken();
  // Use Arr to avoid optimiser removing entirely
  std::cout<<Arr[10]<<std::endl;

  std::cout<<"         Time benefit from best optimization: "<<(((double)first_time / (double)second_time)-1)*100.0<<" %"<<std::endl;
};
