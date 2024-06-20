
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include "PhysicalTypes.hpp"

void debug_checks();

void basic_demo();
void initialisation_and_access_demo();
void more_initialisation_demo();
void exponents_and_roots();
void products_and_functions();
void comparison_demo();

void constexpr_checks();
void last_bits();
void io_checks();
void internal_checks();
void performance_tests();

void run_timer_test();
void run_addition_timer_test();

int main(){

  debug_checks();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Basic demo of Unit Checking"<<std::endl;
  basic_demo();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Ways to initialise Unit Checked variables"<<std::endl;
  initialisation_and_access_demo();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"MORE ways to initialise Unit Checked variables"<<std::endl;
  more_initialisation_demo();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Exponentiation and roots"<<std::endl;
  exponents_and_roots();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Products and special functions"<<std::endl;
  products_and_functions();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Comparison operators"<<std::endl;
  comparison_demo();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Constexpr checks and EVEN MORE ways to initialise"<<std::endl;
  constexpr_checks();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"A few final bits"<<std::endl;
  last_bits();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"IO checks"<<std::endl;
  io_checks();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Some internal stuff"<<std::endl;
  internal_checks();

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Performance tests"<<std::endl;
  performance_tests();

  return 0;

};

void debug_checks(){
  /** Useful checks of fractions and storage types for debugging
  */

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

#ifdef DEBUG
  // Initialization checks etc on Storage types
  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Storage type checks"<<std::endl;

  STScalar s1{1.0};
  std::cout<<"s1 = "<<s1<<std::endl;
  STScalar s2{s1};
  std::cout<<"s2 = "<<s2<<std::endl;

  STVector<double, 3> v1{1.0, 2.0, 3.0};
  std::cout<<"v1 = "<<v1<<std::endl;
  STVector v2{v1};
  std::cout<<"v2 = "<<v2<<std::endl;
  STVector<double, 3> v3{s1, s1, s1};
  std::cout<<"v3 = "<<v3<<std::endl;

  STTensor<double, 3> test_t1{{1.0, 2.0, 3.0},{4.0, 5.0, 6.0},{7.0, 8.0, 9.0}};
  std::cout<<"test_t1 = "<<test_t1<<std::endl;
  STTensor<double, 3> test_t11{1.0, 2.0, 3.0, 4.0, 5.0, 6.0,  7.0, 8.0, 9.0};
  std::cout<<"test_t11 = "<<test_t11<<std::endl;
  STTensor test_t2{test_t1};
  std::cout<<"test_t2 = "<<test_t2<<std::endl;
  STTensor test_t3{v1, v1, v1};
  std::cout<<"test_t3 = "<<test_t3<<std::endl;
  std::cout<<"___________________________________________________________"<<std::endl;

  //Scalar
  STScalar ss = s1 * s1;
  STScalar ssd = ss / s1;
  //Vector
  STVector vv = v1 * v1;
  STVector vvd = vv / v1;
  //Tensor
  STTensor tt = test_t1 * test_t1;
  STTensor ttd = tt / test_t1;
  //Scalar-vector
  STVector m1 = ssd * vvd;
  STVector m2 = m1 * s1;
  STVector m3 = s1 / m2;
  STVector m4 = m3 / s1;
  //Scalar-tensor - only multiply implemented
  STTensor m5 = s1 * test_t1;
  STTensor m6 = m5 * s1;
  //Vector-tensor - only multiply implemented
  STTensor m9 = m4 * m6;
  STTensor m10 = m9 * v1;
  std::cout<<"Squashing some unused variable flags: "<<m10<<" "<<ttd<<std::endl;

#endif

#ifdef DEBUG
  // Debug check making sure that we can instantiate a type without any methods, and thus do not force any methods to be defined
  UnitCheckedType<0, 0, 0, STDummy> dummy;
#endif


};

void basic_demo(){
  /**
   * @brief A basic demo of the unit checking system
   *
   * Initialises some physical quantities, performs some operations, and demonstrates the use of the unit checking system. Mostly shows how x = x + v * t can be performed, and how the units are checked.
   */

  // Create three related physical quantities

  // A time
  Time t{0.1};
  std::cout<<"Defined time t= "<<t<<t.units()<<std::endl;

  std::vector<Time> t_steps;
  t_steps.push_back(t);
  t_steps.push_back(Time{0.2});
  std::cout<<"Defined a vector of type time, containing "<<t_steps.size()<<" elements"<<std::endl;

  // A position
  Position x{1.0, 2.0, 3.0};
  std::cout<<"Defined position x= "<<x<<x.units()<<std::endl;

  // A velocity
  Velocity v{1.0, 2.0, 3.0};
  std::cout<<"Defined velocity v= "<<v<<v.units()<<std::endl;

  Time t2 = t;
  std::cout<<"Two times are the same type?: "<<t.isSameType(t2)<<std::endl;
  std::cout<<"Position and Velocity are not?: "<<x.isSameType(v)<<std::endl;

  // Update position using x = x_0 + v * t
  std::cout<<"Using SUVAT equations, update position using x = x_0 + v * t"<<std::endl;
  x = x + v * t;
  std::cout<<"Position at t ="<<t<<t.units()<<" is "<<x<<x.units()<<std::endl;

  // Update position using x += v * t
  x += v * (t + t2);
  std::cout<<"Position at t ="<<t*2+t2<<t.units()<<" is "<<x<<x.units()<<std::endl;

  // Subtraction test
  auto x_sub = x - x/2.0;
  std::cout<<"For x="<<x<<", x - x/2 = "<<x_sub<<std::endl;


#ifdef FAIL_DEMO
  // Invalid - trying to add a position to a velocity
  auto bad_add = x + v;
#endif

  std::cout<<"______________________________"<<std::endl;

  std::cout<<"One can also multiply by numeric types, and perform a one-over operation"<<std::endl;
  // Numeric multiply
  std::cout<<"2x= "<<2.0*x<<x.units()<<"=="<<x*2.0<<x.units()<<std::endl;
  std::cout<<"x/2= "<<x/2.0<<x.units()<<std::endl;

  auto OneOverT = 1.0/t;
  std::cout<<"Inverse of t is "<<OneOverT<<OneOverT.units()<<std::endl;

  auto OneOverPos = 1.0/Position{1.0, 1.0, 1.0};
  std::cout<<"Inverse of pos is "<<OneOverPos<<OneOverPos.units()<<std::endl;

  auto t_over_t = t/t;
  std::cout<<"t/t is dimensionless: "<<t_over_t<<t_over_t.units()<<std::endl;

  // All SI base units included:
  Energy E{1.0};
  Charge Q{1.0};
  Voltage V = E/Q;
  std::cout<<"As of version 2.0, all 7 SI base units are available: Energy, Charge, Voltage, Time, Length, Mass, Temperature\n";
  std::cout<<"For example Voltage = Energy/Charge = "<<V<<V.units()<<std::endl;

};

void initialisation_and_access_demo(){
  // Construct a vector from a scalar
  std::cout<<"Values can be constructed or uniform initialised ({}-idiom) from underlying data type (here double) \n";
  std::cout<<" Or they can be set equal to other values of correct units\n";
  Length l = Length{1.0};
  Length l2{l};
  std::cout<<"Vectors can be set from numbers, scalars of the same units, or other vectors\n";
  Position xl{l, l, l};
  Position x2{xl};
  Position x3{0.0, 0.0, 0.0};
  std::cout<<"Initialising from the wrong units is not allowed, nor is setting _equal_ to a number, such as x = 0.1;"<<std::endl;
#ifdef FAIL_DEMO
  // Initialising from wrong numeric underlying type not allowed
  Position xlll{0, 0, 0};
  // Initialising from other units not allowed
  Velocity vl{l, Length{0.0}, Length{0.0}};
  // Scalar from vector makes no sense...
  Length l1{xl};
  // Can't set equal to scalar as escapes unit checking again
  Length l3 = 0.2;
#endif

  // Accessor functions
  std::cout<<"It is not recommended, but units can be ignored like this:\n";
  std::cout<<"l's value, forcibly ignoring units, is "<<l.unsafeGet()<<" and can be added to a simple double like this:"<<1.0 + l.unsafeGet()<< std::endl;
  auto scal = UCScalar{7.5};
  std::cout<<"A type without units safely allows getting of the value: "<<scal.get()<<std::endl;

};

void exponents_and_roots(){

  Time t{0.1};

  // Powers
  // NOTE: since exponentiation changes the units, we have to supply the exponent at compile time to get the compile time unit checking. Hence, this is a templated function. As a convenience, there is also a runtime version for DIMENSIONLESS types only
  std::cout<<"Unit checked types with units can be raised to powers, but the power must be known at compile time as this is when the unit checking occurs \n";
  auto tsq = pow<2>(t);
  std::cout<<"t^2= "<<tsq<<tsq.units()<<std::endl;
  auto tsqsqrt = sqrt(tsq);
  std::cout<<"sqrt(t^2)= "<<tsqsqrt<<tsqsqrt.units()<<std::endl;

#ifdef USE_FRACTIONAL_POWERS
  // Here we need FRACTIONAL_POWERS to be able to specify the power as a fraction
  // See below for nth roots using only integers which work as long as the resulting units are integer only
  auto tsqrt2 = pow<SF{1,2}>(tsq);
  std::cout<<"         = "<<tsqrt2<<tsqrt2.units()<<std::endl;
#endif
std::cout<<"Note that non-integer units can only occur if the precompiler define USE_FRACTIONAL_POWERS is set"<<std::endl;
  // Valid only with fractional powers enabled, as units are s^(1/N) (N=2, 3, 5 shown)
#if defined USE_FRACTIONAL_POWERS || defined FAIL_DEMO
  std::cout<<"With the ability to use fractional exponents in the units, we can do any sorts of roots\n";
  auto tsqrt = sqrt(t);
  std::cout<<"sqrt(t)= "<<tsqrt<<tsqrt.units()<<std::endl;
  auto tcbrt = cbrt(t);
  std::cout<<"cbrt(t)= "<<tcbrt<<tcbrt.units()<<std::endl;
  auto fifthroot_t = nthrt<5>(t);
  std::cout<<"5th root of t= "<<fifthroot_t<<fifthroot_t.units()<<std::endl;
#endif
  // But this is always valid as we make sure to have integral units
  std::cout<<"Otherwise we can only do ones where the units are come out to integer powers\n";
  auto fifthpow = pow<5>(t);
  auto fifthroot = nthrt<5>(fifthpow);
  std::cout<<"5th root of t^5= "<<fifthroot<<fifthroot.units()<<std::endl;


#ifdef FAIL_DEMO
  // Not valid - t has units:
  auto tsq2 = pow(t, 2);
#endif
  // Valid - dimensionless type supplied
  auto dimless = t/t;
  std::cout<<"t/t is dimensionless, so can be raised to a runtime power: 2.0^2 = "<< pow(dimless, 2)<<" and 2.0^2.0 = "<<pow(dimless, 2.0)<<std::endl;

};

void products_and_functions(){

  Length l{0.1};
  Position x{1.0, 2.0, 3.0};
  Velocity v{1.0, 1.0, 1.0};

  std::cout<<"A variety of linear algrbra type operations are defined: \n";
  // Dot product
  auto dt = x.dot(v);
  std::cout<<"x.v="<<dt<<dt.units()<<std::endl;

  // Cross product
  auto crs = x.cross(v);
  std::cout<<"x cross v="<<crs<<crs.units()<<std::endl;

  // Outer product
  auto out = x.outer(v);
  std::cout<<"x outer v="<<out<<out.units()<<std::endl;


  //Normalize, two ways
  auto x2 = x.normalized();
  x.normalize();
  std::cout <<"x normalised is " << x << std::endl;
  std::cout <<"              or " << x2 << std::endl;
  std::cout << "magnitude of x normalised is " << x.magnitude() << std::endl;

  std::cout<<"Using special functions can be done two ways, depending on what is defined by the stored data types\n";
  // Using special functions on dimensionless values
  // For scalars, cast to, or store to, a double
  std::cout<<"sin(l/l_0): "<<sin(static_cast<double>(l/Length{1.0}))<<std::endl;
  // OR use specialisation (if provided) on the underlying value
  std::cout<<"sin(l/l_0): "<<mysinfunction((l/Length{1.0}).data())<<std::endl;
  // Storage type can provide operations on higher rank types, either elementwise or whatever
  std::cout<<"sin(x/x_0): "<<mysinfunction((x/Length{1.0}).data())<<std::endl;

};

void more_initialisation_demo(){

  std::cout<<"For convenience, there are some initialisation shortcuts\n";
  // Identities
  Length il=makeIdentity<Length>();
  std::cout<<"Identity Scalar :"<<il<<std::endl;
  Position ip=makeIdentity<Position>();
  std::cout<<"Identity Vector :"<< ip<<std::endl;
  UCTensor it=makeIdentity<UCTensor>();
  std::cout<<"Identity Tensor :\n"<<it<<std::endl;
  UCTensor it2{tens3_ident};
  std::cout<<"Identity Tensor another way :\n"<<it2<<std::endl;

  // Using casting to construct
  UCScalar l{1.0};
  double tmp2{l};
  std::cout<<"Casting to double from type with no units: "<< typeid(tmp2).name()<<std::endl;
#if !defined NO_NARROWING_CONVERSIONS || defined FAIL_DEMO
  float tmp1{il};  //Disallowed by -DNO_NARROWING_CONVERSIONS
  std::cout<<"Narrowing to float: "<< typeid(tmp1).name()<<std::endl;
#endif

};

void comparison_demo(){
  Length l{1.0};
  Position x{1.0, 2.0, 3.0};
  // Demo of comparison operators
  std::cout<<"We can compare values of the same units\n";

  std::cout<<"l > 0? "<< (l > Length{0.0})<<std::endl;
  std::cout<<"x == x? "<< (x == x)<<std::endl;
  std::cout<<"x != x? "<< (x != x)<<std::endl;

  std::cout<<"x.magnitude()= "<<x.magnitude()<<std::endl;
  std::cout<<"l.magnitude()= "<<l.magnitude()<<std::endl;

  std::cout<< "x > l? (using magnitude) "<< (x.magnitude() > l.magnitude())<<std::endl;
  std::cout<< "x > l? (using scalar-vector comparator (magnitude behind the scenes))"<< (x > l)<<std::endl;

  // The next block just checks resolution of all the comparisons for completeness
  std::cout<<"Doing all the Scalar-Scalar combos: "<<(l < l)<<" "<<(l<=l)<<" "<<(l==l)<<" "<<(l!=l)<<" "<<(l>=l)<<" " <<(l > l)<<std::endl;
  std::cout<<"Doing all the Vector-Vector combos: "<<(x < x)<<" "<<(x<=x)<<" "<<(x==x)<<" "<<(x!=x)<<" "<<(x>=x)<<" " <<(x > x)<<std::endl;
  std::cout<<"Doing all the Scalar-vector combos: "<<(x < l)<<" "<<(x<=l)<<" "<<(x==l)<<" "<<(x!=l)<<" "<<(x>=l)<<" " <<(x > l)<<std::endl;
  std::cout<<"And doing them the other way around: "<<(l < x)<<" "<<(l<=x)<<" "<<(l==x)<<" "<<(l!=x)<<" "<<(l>=x)<<" " <<(l > x)<<std::endl;

  //Tensors only have element-wise equality defined
  UCTensor t{{1.0, 2.0, 3.0},{4.0, 5.0, 6.0},{7.0, 8.0, 9.0}};
  std::cout<<"Doing Tensor-Tensor equality: "<<(t == t)<<" "<<(t!=t)<<std::endl;


  //Using casting in comparison
  std::cout<<"A dimensionless value can be compared to a double in two ways: \n";
  UCScalar l2{1.0};
  std::cout<<"Construct comparable type, val > UCScalar{1.0} ?"<<(l2 > UCScalar{1.0})<<std::endl;
  std::cout<<"Cast UCScalar to double, static_cast<double>(val) > 1.0 ? "<<(static_cast<double>(l2) > 1.0)<<std::endl;

};

void constexpr_checks(){
  // Checking that we can apply constexpr for compile-time values
  constexpr Length l3{1.0};
  std::cout<<"Values can be used in constexpr contexts: constexpr l3 = "<<l3<<l3.units()<<std::endl;
  constexpr Position p{0.1, 0.2}, p2{0.1}, p3(0.1);
  std::cout<<"Missing values are zero initialised: "<<p<<std::endl;
  std::cout<<"But a single value performs Scalar initialisation: "<<p2<<std::endl;
  std::cout<<"Whichever way you do it: "<<p3<<std::endl;
  constexpr Position p4{l3, l3, l3}, p5{l3};
  std::cout<<"And you can initialise a vector from 1, or n, Scalars: "<<p4<<" "<<p5<<std::endl;
  // Show shortcut Scalar initialisation on a Tensor
  constexpr UCTensor TT{1.0};
  std::cout<<"Which is very useful on a Tensor TT{1.0} -> \n"<<TT<<std::endl;

  // More initialisation tests on tensors
  constexpr UCTensor TT2{{1.0, 2.0, 3.0},{4.0, 5.0, 6.0},{7.0, 8.0, 9.0}};
  constexpr UCTensor TT3{1.0, 2.0, 3.0,  4.0, 5.0, 6.0,  7.0, 8.0, 9.0}; // Don't have to provide the inner braces IF the number of elements is correct (fills as many as given)
  std::cout<<"Initialisation from nested lists: \n"<<TT2<<std::endl;
  std::cout<<"Initialisation from flat list (length > 1): \n"<<TT3<<std::endl;
  std::cout<<"Transpose of TT3: \n"<<TT3.transpose()<<std::endl;


  Stress s{1.0}; // A Scalar
  StressVector sv{s,s,s}; // A 3-vector constructed from a scalar
  StressTensor st{{s,s,s},{s,s,s},{s,s,s}}; // Tensor from scalars
  StressTensor st2{sv, sv, sv}; // Tensor from vectors


  // Try gridded operations and show how to initialise
  UCScalar Arr[3][3];
  Arr[0][0] = UCScalar{1.0};
  Arr[1][1] = UCScalar{1.0};
  Arr[2][2] = UCScalar{1.0};
  std::cout<<"Arr[0][0] = "<<Arr[0][0]<<std::endl;
  std::cout<<"Arr[0][1] = "<<Arr[0][1]<<std::endl;
  std::cout<<"Arr[1][1] = "<<Arr[1][1]<<std::endl;

  auto Arr2 = Arr;
  std::cout<<"Arr2[0][0] = "<<Arr2[0][0]<<std::endl;
  std::cout<<"Arr2[0][1] = "<<Arr2[0][1]<<std::endl;
  std::cout<<"Arr2[1][1] = "<<Arr2[1][1]<<std::endl;

};

void last_bits(){

  Length l{0.1};

#ifdef USE_FRACTIONAL_POWERS

  // Creating your own custom physical type to use elsewhere
  using UnsquareMeters = UnitCheckedType<SF{1,2}, 0, 0, dblscalar>;
  UnsquareMeters lu{0.3};
  std::cout<<"A quantity with fractional exponent units: lu = "<<lu<<lu.units()<<std::endl;
  std::cout<<" l+ lu^2 = "<<l + lu*lu<<l.units()<<std::endl;

#endif

  //If we want to have different names for our units, we can manipulate the _strings used for display_
  UnitChecked::unitNames[0] = "mic";
  std::cout<<"Units can be customised: Length will now show as "<<l<<l.units()<<std::endl;
};

void io_checks(){

  // Setup a stream, push some entities into it, and then extract them again
  std::stringstream ss;

  Time t{12.34};
  ss<<t;
  Time t2;
  ss>>t2;
  std::cout<<"Wrote "<<t<<" and read "<<t2<<std::endl;
  if(t != t2){
    throw std::runtime_error("Time read from stream is not the same as written");
  }

  ss.clear(); // Because we hit the end of the stream before
  Position x{1.2, 3.4, 5.6};
  ss<<x;
  Position x2;
  ss>>x2;
  std::cout<<"Wrote "<<x<<" and read "<<x2<<std::endl;
  if(x != x2){
    throw std::runtime_error("Position read from stream is not the same as written");
  }

  ss.clear();
  UCTensor t3{{1.0, 2.0, 3.0},{4.0, 5.0, 6.0},{7.0, 8.0, 9.0}};
  ss<<t3;
  UCTensor t4;
  ss>>t4;
  std::cout<<"Wrote "<<t3<<" and read "<<t4<<std::endl;
  if(t3 != t4){
    throw std::runtime_error("Tensor read from stream is not the same as written");
  }

  // Now do multiple
  ss.clear();
  ss<<t<<x<<t<<x;
  Time t5, t6;
  Position x5, x6;
  ss>>t5>>x5>>t6>>x6;
  std::cout<<"Wrote "<<t<<" "<<x<<" "<<t<<" "<<x<<" and read "<<t5<<" "<<x5<<" "<<t6<<" "<<x6<<std::endl;
  if(t5 != t || x5 != x || t6 != t || x6 != x){
    throw std::runtime_error("Multiple reads from stream failed");
  }

};

void internal_checks(){
  // Checks on stuff user shouldn't need to know about but that should be done somewhere

  // Mostly to do with constness, constexpr, and references (lval, rval etc)

  std::cout<<"Checking if Time is Unitchecked "<<is_unitchecked_v<Time><<std::endl;
  std::cout<<"Checking if double is Unitchecked "<<is_unitchecked_v<double><<std::endl;
  std::cout<<"Checking if Time is numeric (storage type) "<<is_unitchecked_numeric_v<Time><<std::endl;

  UnitCheckedType<0, 0, 0, std::string> uct;
  std::cout<<"Checking a non-numeric storage type "<<is_unitchecked_numeric_v<decltype(uct)><<std::endl;


  //Verifying all of the initialiser list options
  Length l{0.0}; // From underlying numeric type
  Length l2{l}; // From self-type
  Position x{l, l, l}; // From smaller storage, same units
  dbl3vec d{1.0, 1.0, 1.0};
  Position x2{d}; // From storage type, no units
  constexpr UCVector v{1.0}; // Convenience - from scalar
  static_assert(v.get(1) == 1.0);
  UCTensor t{v, v, v};
  dbl3tens tt{tens3_ident};
  UCTensor t2{tt};
  UCScalar s{1.0};
  UCTensor t3{{s, s, s}, {s, s, s}, {s, s, s}};
 
#ifdef FAIL
  Position x2{l}; // No - one length is not OK
  //UCTensor t33{s, s, s}; // Works, so that we can allow single nest with correct total number
  UCTensor t4{v}; // Ditto - too small
#endif
}


void performance_tests(){

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
