
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include <cassert>
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
  std::cout<<"Checking fractions"<<std::endl;
  SimpleFrac aa{1,2}, bb{2,4}, cc{1,3};
  assert(divides(aa, 2));
  assert(SimpleFrac{3} == (SimpleFrac{3,1}));
  SimpleFrac zero{0,0}; // A bad number, but make sure we can construct it
  assert(zero + zero == zero); //Catches that final branch
  std::cout<<"aa = "<<aa<<std::endl;
  std::cout<<"bb = "<<bb<<std::endl;
  std::cout<<"cc = "<<cc<<std::endl;
  std::cout<<"cc as a double "<<(double)cc<<std::endl;
  std::cout<<"aa == bb? "<<(aa==bb)<<std::endl;
  assert(aa != bb);
  std::cout<<"aa equal bb? "<<is_equal(aa, bb)<<std::endl;
  assert(is_equal(aa, bb));
  std::cout<<"bb is simplifiable? "<<simplifiable(bb)<<std::endl;
  assert(simplifiable(bb) && !simplifiable(aa));
  SimpleFrac dd{3,2};
  assert(is_greater(dd, 1) && is_less(dd, 2));
  assert(is_equal(dd*2, 3));
  std::cout<<dd<<" + "<<aa<<" = "<<dd+aa<<std::endl; // Output for integer
  std::cout<<"aa < cc? "<<is_less(aa,cc)<<std::endl;
  std::cout<<"aa > cc? "<<is_greater(aa,cc)<<std::endl;
  assert(is_greater(aa, cc));
  std::cout<<"aa + bb - cc = "<<aa+bb-cc<<std::endl;
  assert(aa+bb-cc == (SimpleFrac{2,3}));
  assert(aa + (-cc) == aa - cc);
  std::cout<<"aa*bb/cc = "<<aa*bb/cc<<std::endl;
  assert(aa*bb/cc == (SimpleFrac{3,4}));

  std::cout<<"___________________________________________________________"<<std::endl;
#ifdef FAIL_DEMO
  // Check that we can't construct a type with fraction with a zero denominator
  UnitCheckedType<SF{1,0}, 0, 0, STDummy> dd1;
  // Check that we can't construct a fraction in non-reduced terms
  UnitCheckedType<0, SF{2,4}, 0, STDummy> dd2;
  //Check that we can't have negative denominator
  UnitCheckedType<0, 0, SF{1,-1}, STDummy> dd3;
#endif
#endif
#endif

#ifdef DEBUG
  // Initialization checks etc on Storage types
  std::cout<<"Storage type checks"<<std::endl;

  STScalar s1{1.5};
  STScalar s2{s1};
  assert(s1 == s2);
  assert(s1 != -s2);
  STScalar<double> s3;
  s3 = s1;

  STVector<double, 3> v1{1.0, 2.0, 3.0};
  STVector v2{v1};
  assert(v1 == v2);
  STVector<double, 3> v3{s1, s1, s1};
  assert(v3[0] == s1 && v3[1] == s1 && v3[2] == s1);
  STVector<double, 3> v4;
  v4 = v1;

  STTensor<double, 3> t1{{1.0, 2.0, 3.0},{4.0, 5.0, 6.0},{7.0, 8.0, 9.0}};
  STTensor<double, 3> t11{1.0, 2.0, 3.0, 4.0, 5.0, 6.0,  7.0, 8.0, 9.0};
  STTensor t2{t1};
  assert(t1 == t11);
  assert(t1 == t2);
  STTensor t3{v1, v1, v1};
  assert(t3.get(0,0) == v1[0] && t3.get(1,0) == v1[0] && t3.get(2,0) == v1[0]);
  assert(t3.get(0,1) == v1[1] && t3.get(1,1) == v1[1] && t3.get(2,1) == v1[1]);
  assert(t3.get(0,2) == v1[2] && t3.get(1,2) == v1[2] && t3.get(2,2) == v1[2]);
  STTensor<double, 3> t4;
  t4 = t1;

  // Quick plus and minus - these get touched later once wrapped but are useful to check here too
  s3 = s1 + s2;
  assert(s3 - s2 == s1);
  s3 += s2;
  s3 -= s1;
  assert(s3 == 2.0 * s2);
  v4 = v1 + v2;
  assert(v4 - v2 == v1);
  v4 += v1;
  v4 -= v2;
  assert(v4 == 2.0 * v2);
  t4 = t1 + t2;
  assert(t4 - t2 == t1);
  t4 += t2;
  t4 -= t1;
  assert(t4 == 2.0 * t2);

  //Scalar
  STScalar ss = s1 * -s1;
  STScalar ssd = ss / s1;
  assert(ssd == -s1);
  //Vector
  STVector vv = v1 * -v1;
  STVector vvd = vv / v1;
  assert(vvd == -v1);
  //Tensor
  STTensor tt = t1 * -t1;
  STTensor ttd = tt / t1;
  assert(ttd == -t1);
  //Scalar-vector
  STVector m1 = ssd * vvd;
  STVector m2 = m1 * s1;
  STVector m3 = m2 / s1;
  STVector m4 = s1 / m3;
  assert(m3 == m1);
  assert(m4[0] == s1.get()/m3[0] && m4[1] == s1.get()/m3[1] && m4[2] == s1.get()/m3[2]);
  //Scalar-tensor - only multiply implemented
  STTensor m5 = s1 * t1;
  STTensor m6 = m5 * s1;
  assert(m6.get(0,0) == s1.get()*s1.get()*t1.get(0,0) && m6.get(1,2) == s1.get()*s1.get()*t1.get(1,2));
  // m6.magnitude() == s1.get()*s1.get()*test_t1.magnitude());
  //Vector-tensor - only multiply implemented
  STTensor m9 = v1 * m6;
  STTensor m10 = m9 * v1;
  assert(m10.get(0,0) == m6.get(0,0)*v1[0]*v1[0] && m10.get(1,2) == m6.get(1,2)*v1[1]*v1[1]);

  //Casting checks. Use a float since generally exclude narrowing conversions
  STScalar<float> sf = 2.3;
  double d = (double)sf;
  float f = (float)sf;
  assert(d == f);

  //Some other functions for Scalars
  STScalar sabs = s1.magnitude();
  assert(sabs == s1);
  auto rt_s = s1.sqrt();
  auto pow_s = rt_s.pow(2.0);
  assert((pow_s - sabs) < 1e-15);//Can only be approx...
  rt_s = s1.cbrt();
  pow_s = rt_s.pow(3.0);
  assert((pow_s - sabs) < 1e-15);//Ditto

  //And for Vectors
  v1[1] = v1[0];// Write access
  assert(v1[0] == v1[1]);
  v1.normalize();
  assert(v1.magnitude() == 1.0); // Normalise and magnitude
  //Make sure its representible exactly and a nice value
  v1 = {2.0, 3.0, 4.0};
  v2 = {3.0, 4.0, 5.0};
  auto v_dot = v1.dot(v2);
  assert(v_dot == 38.0);
  auto v_cross = v1.cross(v2);
  assert(v_cross[0] == -1.0 && v_cross[1] == 2.0 && v_cross[2] == -1.0);
  auto v_outer = v1.outer(v2);
  assert(v_outer.get(0,0) == 6.0 && v_outer.get(1,2) == 15.0);

  //And for Tensors
  t1.get(1,0) = t1.get(0, 0);// Write access
  assert(t1.get(1,0) == t1.get(0, 0));
  auto t1_t = t1.transpose();
  assert(t1_t != t1);
  assert(t1_t.get(0,0) ==t1.get(0,0) && t1_t.get(0,1) == t1.get(1,0) && t1_t.get(2,2) == t1.get(2,2) && t1_t.get(2,1) == t1.get(1,2));

  //Rest of the comparators on Scalars
  s2 = s1/2.0;
  assert( s2 != s1 && s2 < s1 && !(s2 >= s1) && s2 <= s1 && !(s2 > s1));
  // And vectors
  v2 = v1/2.0;
  assert(v2 != v1 && v2 < v1 && !(v2 >= v1) && v2 <= v1 && !(v2 > v1));
  // And cross compare
  s1 = v1.magnitude() * 2.0;
  assert(s1 != v1 && s1 > v1 && !(s1 <= v1) && s1 >= v1 && !(s1 < v1));

  //Streaming
  std::stringstream str_s;
  s1 = 1.5;
  str_s<<s1;
  STScalar<double> sst;
  str_s>>sst;
  assert(s1 == sst);

  str_s.clear();
  str_s<<v1;
  STVector<double, 3> vst;
  str_s>>vst;
  assert(v1 == vst);

  str_s.clear();
  str_s<<t1;
  STTensor<double, 3> tst;
  str_s>>tst;
  assert(t1 == tst);

#endif

#ifdef DEBUG
  // Debug check making sure that we can instantiate a type without any methods, and thus do not force any methods to be defined
  const UnitCheckedType<0, 0, 0, STDummy> dummy;
  std::cout<<&dummy<<std::endl; //Use dummy to make sure it must exist

  std::cout<<"Done"<<std::endl;
  std::cout<<"___________________________________________________________"<<std::endl;

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
  assert(t.isSameType(t2) && !x.isSameType(v));

  // Update position using x = x_0 + v * t
  std::cout<<"Using SUVAT equations, update position using x = x_0 + v * t"<<std::endl;
  x = x + v * t;
  std::cout<<"Position at t ="<<t<<t.units()<<" is "<<x<<x.units()<<std::endl;

  // Update position using x += v * t
  t += t;
  x += v * (t + t2);
  std::cout<<"Position at t ="<<(t+t2)<<t2.units()<<" is "<<x<<x.units()<<std::endl;
  std::cout<<(Position{1.0, 2.0, 3.0} + 4.0*v*Time{0.1})<<std::endl;
  assert((x - (Position{1.0, 2.0, 3.0} + 4.0*v*Time{0.1})).magnitude() < Length{1e-15});

  // Subtraction test
  auto x_sub = x - x/2.0;
  std::cout<<"For x="<<x<<", x - x/2 = "<<x_sub<<std::endl;

  //Quick results check on operators
  assert( x+x == x*2);
  assert(x_sub == x/2.0);
  t -= Time{0.1};
  assert(t == Time{0.1});
  t /= 2.0;
  assert(t == Time{0.05});
  auto y = Position{1.0, 2.0, 3.0};
  y -= Position{0.5, 1.0, 1.5};
  assert(y == (Position{0.5, 1.0, 1.5}));
  y /= 5.0;
  assert(y == (Position{0.1, 0.2, 0.3}));

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
  Length l = Length{1.7};
  Length l2{l};
  std::cout<<"Vectors can be set from numbers, scalars of the same units, or other vectors\n";
  Position xl{l, l, l};
  Position x2{xl};
  Position x3{0.0, 0.0, 0.0};
  assert(x3+x3 == x3);//Using the value
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

  //Sometimes we want to access an element of a longer type RETAINING units
  // Currently this WILL produce a copy, so use with caution in performance critical code
  auto l0 = xl.getElement(0);
  std::cout<<"Accessing an element of a vector, retaining units: l0 = "<<l0<<l0.units()<<std::endl;
  std::cout<<"This is useful sometimes for comparisons etc, avoiding the unsafe Idiom below. E.g. compare l and xl[0]: "<<(l == l0)<<std::endl;


  // Accessor functions
  std::cout<<"It is not recommended, but units can be ignored like this:\n";
  std::cout<<"l's value, forcibly ignoring units, is "<<l.unsafeGet()<<" and can be added to a simple double like this:"<<1.0 + l.unsafeGet()<< std::endl;
  auto scal = UCScalar{7.5};
  std::cout<<"A type without units safely allows getting of the value: "<<scal.get()<<std::endl;

#ifdef FAIL_DEMO
  // Can't access the value of a type with units
  auto val = l.get();
#endif


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
  assert((tsqsqrt-t).magnitude() < Time{1e-15});

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
  std::cout<<"Otherwise we can only do ones where the units come out to integer powers\n";
  auto fifthpow = pow<5>(t);
  auto fifthroot = nthrt<5>(fifthpow);
  std::cout<<"5th root of t^5= "<<fifthroot<<fifthroot.units()<<std::endl;
  assert((fifthroot - t).magnitude() < Time{1e-15});


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

  std::cout<<"A variety of linear algebra type operations are defined: \n";
  // Dot product
  auto dt = x.dot(v);
  std::cout<<"x.v="<<dt<<dt.units()<<std::endl;

  // Cross product
  auto crs = x.cross(v);
  std::cout<<"x cross v="<<crs<<crs.units()<<std::endl;

  // Outer product
  auto out = x.outer(v);
  std::cout<<"x outer v="<<out<<out.units()<<std::endl;

  //All of the above multiply the units, hence:
  assert(crs.isSameUnits(out) && crs.isSameUnits(dt));
  //Note numeric results have been checked by StorageType already


  //Normalize, two ways
  auto x2 = x.normalized();
  x.normalize();
  std::cout <<"x normalised is " << x << std::endl;
  std::cout <<"              or " << x2 << std::endl;
  std::cout << "magnitude of x normalised is " << x.magnitude() << std::endl;
  assert(x.isSameType(x2));

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
  float tmp1{tmp2};  //Disallowed by -DNO_NARROWING_CONVERSIONS: Error in Clang, Warning in GCC
  std::cout<<"Narrowing to float: "<< typeid(tmp1).name()<<std::endl;
#endif

};

void comparison_demo(){
  Length l{1.0};
  Position x{1.0, 2.0, 3.0};
  // Demo of comparison operators
  std::cout<<"We can compare values of the same units\n";

  std::cout<<"l > 0? "<< (l > Length{0.0})<<std::endl;
  assert(l > Length{0.0});
  std::cout<<"x == x? "<< (x == x)<<std::endl;
  std::cout<<"x != x? "<< (x != x)<<std::endl;
  assert(x == x && !(x != x));

  std::cout<<"x.magnitude()= "<<x.magnitude()<<std::endl;
  std::cout<<"l.magnitude()= "<<l.magnitude()<<std::endl;

  std::cout<< "x > l? (using magnitude) "<< (x.magnitude() > l.magnitude())<<std::endl;
  std::cout<< "x > l? (using scalar-vector comparator (magnitude behind the scenes))"<< (x > l)<<std::endl;

  // The next block just checks resolution of all the comparisons for completeness
  assert(!(l<l) && (l<=l) && (l==l) && !(l!=l) && (l>=l) && !(l>l));
  assert(!(x<x) && (x<=x) && (x==x) && !(x!=x) && (x>=x) && !(x>x));
  assert(!(x<l) && !(x<=l) && !(x==l) && (x!=l) && (x>=l) && (x>l));
  assert((l<x) && (l<=x) && !(l==x) && (l!=x) && !(l>=x) && !(l>x));

  //Tensors only have element-wise equality defined
  UCTensor t{{1.0, 2.0, 3.0},{4.0, 5.0, 6.0},{7.0, 8.0, 9.0}};
  std::cout<<"Doing Tensor-Tensor equality: "<<(t == t)<<" "<<(t!=t)<<std::endl;
  assert(t == t && !(t != t));


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

  Force f{1.0};
  Voltage v{0.1};
  //If we want to have different names for our units, we can manipulate the _strings used for display_
  //Fundamental types can be changed, and will print this way for ALL compounds
  std::cout<<"Length is Fundamental (checking 2 ways)? "<<l.isFundamentalType()<<" or "<<Length::isFundamentalType()<<std::endl;
  std::cout<<"Force is not? "<<f.isFundamentalType()<<" or "<<Force::isFundamentalType()<<std::endl;
  std::cout<<"Voltage is also not? "<<v.isFundamentalType()<<" or "<<Voltage::isFundamentalType()<<std::endl;

  UnitChecked::registerUnits<Length>("mic");
  std::cout<<"Units can be customised: Length will now show as "<<l<<l.units()<<" and in compounds such as Force as "<<f.units()<<std::endl;
  //Or with the convenience function like this
  UnitChecked::registerUnitsForTypeOf(l, "mm");
  std::cout<<"And now Length will show as "<<l<<l.units()<<" and Force as "<<f.units()<<std::endl;

  // Or for compound types can provide a fully custom string per TYPE (not instance) used from now on
  // NOTE: this is for this type and ONLY this type, not any further compounded types
  UnitChecked::registerUnits<Force>("N");
  std::cout<<"Force will now show as "<<f<<f.units()<<std::endl;
  //Or use the "helper"
  UnitChecked::registerUnitsForTypeOf(f, "kN ");
  Force f2{2.0};
  std::cout<<"Similarly ALL forces will now show as "<<f2<<f2.units()<<std::endl;
  auto fl = f*l;
  std::cout<<"But compounded types will still show as "<<fl.units()<<std::endl;

  //Using a string variable rather than a literal
  std::string units = "lbs";
  UnitChecked::registerUnits<Force>(units);
  std::cout<<"Force will now show as "<<f<<f.units()<<std::endl;

  UnitChecked::registerUnits<Pressure>("kPa");
  Pressure p{1.0};
  std::cout<<"And now Pressure will show as "<<p<<p.units()<<std::endl;


  //Using the helper for SI common set
  registerSIUnitStrings();
  std::cout<<"And now we have some nice _printing_ for common SI units such as a Voltage like "<<v<<v.units()<<" and Pressure like "<<p<<p.units()<<std::endl;
  std::cout<<"But remember this wont affect further compounds, such as Voltage/Length which is still "<<(v/l).units()<<std::endl;
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
  // Finishes up the code coverage too

  std::cout<<"Checking if Time is Unitchecked "<<is_unitchecked_v<Time><<std::endl;
  static_assert(is_unitchecked_v<Time>);
  std::cout<<"Checking if double is Unitchecked "<<is_unitchecked_v<double><<std::endl;
  static_assert(!is_unitchecked_v<double>);
  std::cout<<"Checking if Time is numeric (storage type) "<<is_unitchecked_numeric_v<Time><<std::endl;
  static_assert(is_unitchecked_numeric_v<Time>);

  class blank{}; // A type which is not numeric
  UnitCheckedType<0, 0, 0, STScalar<blank> > uct;
  std::cout<<"Checking a non-numeric storage type "<<is_unitchecked_numeric_v<decltype(uct)><<std::endl;
  static_assert(!is_unitchecked_numeric_v<decltype(uct)>);

  //Verifying all of the initialiser list options compile
  Length l{1.5}; // From underlying numeric type
  Length l2{l}; // From self-type
  assert(l == l2);
  Position x{l, l, l}; // From smaller storage, same units
  dbl3vec d{1.5, 1.5, 1.5};
  Position x2{d}; // From storage type, no units
  assert(x == x2);
  constexpr UCVector v{1.0}; // Convenience - from scalar
  static_assert(v.get(1) == 1.0); // Checks constexpr-ness - else can't static_assert
  UCTensor t{v, v, v};
  dbl3tens tt{tens3_ident};
  UCTensor t2{tt};
  UCScalar s{1.0};
  UCTensor t3{{s, s, s}, {s, s, s}, {s, s, s}};
  assert(t == t3);

  //Checking undersized lists - check values and zeros
  UCVector shrt{s, s};
  for(int i = 0; i < 3; i++){
    if(i < 2){
      assert(shrt.get(i) == s.get());
    }else{
      assert(shrt.get(i) == 0.0);
    }
  }
  UCTensor t4{shrt, shrt, shrt}, t5{shrt, shrt};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      if(j < 2){
        assert(t4.get(i, j) == s.get());
      }else{
        assert(t4.get(i, j) == 0.0);
      }
      if(i < 2 && j < 2){
        assert(t5.get(i, j) == s.get());
      }else{
        assert(t5.get(i, j) == 0.0);
      }
    }
  }
  UCVector shrt2{1.0}; // Single element sets all
  for(int i = 0; i < 3; i++){
    assert(shrt2.get(i) == s.get());
  }
  UCTensor t4_2{{1.0, 1.0}, {1.0, 1.0}};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      if(i < 2 && j < 2){
        assert(t4_2.get(i, j) == 1.0);
      }else{
        assert(t4_2.get(i, j) == 0.0);
      }
    }
  }


  // This line is just for coverage purposes, we already checked this works
  // but it seems we have to forcibly use it to trigger the coverage checker
  std::cout<<l.hasNoUnits()<<std::endl;
  // These cover some more branches of the units construction/print
  UnitCheckedTypeFull<1, 1, 1, 1, 1, 1, 1, dblscalar> allOne;
  UnitCheckedTypeFull<2, 2, 2, 2, 2, 2, 2, dblscalar> allTwo;
  UnitCheckedTypeFull<-1, -1, -1, -1, -1, -1, -1, dblscalar> allMinus;
  std::cout<<allOne.units()<<" "<<allTwo.units()<<" "<<allMinus.units()<<std::endl;

  //Verify that hasNoUnits and isSameUnits catch "all" the cases (first mismatch)
  UnitCheckedTypeFull<0, 1, 1, 1, 1, 1, 1, dblscalar> U1;
  UnitCheckedTypeFull<0, 0, 1, 1, 1, 1, 1, dblscalar> U2;
  UnitCheckedTypeFull<0, 0, 0, 1, 1, 1, 1, dblscalar> U3;
  UnitCheckedTypeFull<0, 0, 0, 0, 1, 1, 1, dblscalar> U4;
  UnitCheckedTypeFull<0, 0, 0, 0, 0, 1, 1, dblscalar> U5;
  UnitCheckedTypeFull<0, 0, 0, 0, 0, 0, 1, dblscalar> U6;
  assert(!U1.hasNoUnits() && !U2.hasNoUnits() && !U3.hasNoUnits() && !U4.hasNoUnits() && !U5.hasNoUnits() && !U6.hasNoUnits());
  assert(!allOne.isSameUnits(U1) && !U1.isSameUnits(U2) && !U2.isSameUnits(U3) && !U3.isSameUnits(U4) && !U4.isSameUnits(U5));
  std::cout<<U1.hasNoUnits()<<" "<<U2.hasNoUnits()<<" "<<U3.hasNoUnits()<<" "<<U4.hasNoUnits()<<" "<<U5.hasNoUnits()<<" "<<U6.hasNoUnits()<<std::endl;
  std::cout<<U1.isSameUnits(U2)<<" "<<U2.isSameUnits(U3)<<" "<<U3.isSameUnits(U4)<<" "<<U4.isSameUnits(U5)<<" "<<U5.isSameUnits(U6)<<std::endl;

  //Ditto for isFundamentalType
  assert(!allOne.isFundamentalType() && !U1.isFundamentalType() && !U2.isFundamentalType() && !U3.isFundamentalType() && !U4.isFundamentalType() && !U5.isFundamentalType() && U6.isFundamentalType());
  // Cases with some exponents not unity, esp cases where they sum to zero
  UnitCheckedTypeFull<1, -1, 1, -1, 1, -1, 0, dblscalar> U7;
  UnitCheckedTypeFull<0, 2, 2, 2, -2, -2, -2, dblscalar> U8;
  UnitCheckedTypeFull<1, 2, 3, 4, 5, 6, 7, dblscalar> U9;
  UnitCheckedTypeFull<1, 2, 0, 3, -3, -2, -1, dblscalar> U10;
  UnitCheckedTypeFull<0, 0, 0, 2, 0, 0, 0, dblscalar> U11;
  UnitCheckedTypeFull<0, -1, 0, 0, 0, 0, 0, dblscalar> U12;
  assert(!U7.isFundamentalType() && !U8.isFundamentalType() && !U9.isFundamentalType() && !U10.isFundamentalType() && !U11.isFundamentalType() && !U12.isFundamentalType());
  // And check we get the right index from whichFundamentalType
  assert(l.whichFundamentalType() == 1 && (U1/U2).whichFundamentalType() == 0 && (U2/U3).whichFundamentalType() == 2 && (U3/U4).whichFundamentalType() == 3 && (U4/U5).whichFundamentalType() == 4 && (U5/U6).whichFundamentalType() == 5 && U6.whichFundamentalType() == 6);


  //Verify default init behaviour
  UCScalar def{};
  assert(def.get() == 0.0);

  //Fully Verify getElement
  auto el = x.getElement(0);
  assert(x.isSameUnits(el) && !x.isSameType(el));
  assert(el.unsafeGet() == x.unsafeGet(0));
  //Tensor getting vector \todo Is this the useful way round?
  auto el2 = t.getElement(0);
  assert(t.isSameUnits(el2) && !t.isSameType(el2));
  assert(el2.unsafeGet(0) == t.unsafeGet(0, 0));
  //Tensor getting Scalar
  auto el3 = t.getElement(1, 2);
  assert(t.isSameUnits(el3) && !t.isSameType(el3));
  assert(el3.unsafeGet() == t.unsafeGet(1, 2));

  // Trying out a reference getter
  auto el4 = x.getElementRef(1), el5 = x.getElementRef(2);
  el4=Length{2.0};
  el5 = el4+Length{1.0};
  el4 += Length{0.2};
  std::cout<<"Now we have set x using a reference "<<x<<std::endl;
  assert(x == (Position{1.5, 2.2, 3.0}));
  el4 -= Length{0.2};
  el5 = el5 - Length{1.0};
  std::cout<<"And we've changed x some more "<<x<<std::endl;
  assert(x == (Position{1.5, 2.0, 2.0}));

  //Testing some other things
  //Summing multiple references
  auto el_sum = el4 + el5;
  assert(el_sum == Length{4.0});
  //Subtraction
  auto el_sub = el4 - el5;
  assert(el_sub == Length{0.0});

  using LengthSq = UnitCheckedType<2, 0, 0, dblscalar>;
  using LengthInv = UnitCheckedType<-1, 0, 0, dblscalar>;

  auto el_mult = el4 * el5; //Multiplying references
  auto el_mult2 = el * el4; //Multiplying reference and value
  auto el_mult3 = el4 * el; //Multiplying value and reference
  assert(el_mult == LengthSq{4.0} && el_mult2 == LengthSq{3.0} && el_mult3 == LengthSq{3.0});
  auto el_div = el4 / el5; //Dividing references
  auto el_div2 = el_mult / el4; //Dividing reference and value
  auto el_div3 = el4 / el_mult; //Dividing value and reference
  assert(el_div == UCScalar{1.0} && el_div2 == Length{2.0} && el_div3 == LengthInv{0.5});

  //Cursory check with tensors
  auto el6 = t.getElementRef(1, 2);
  el6 = UCScalar{1.7};
  std::cout<<"And now we have modified t using a reference-to-element\n "<<t<<std::endl;
  assert(t.get(1, 2) == 1.7);

  //Get const ref - this works fine and we can access it
  const auto y=x;
  auto el8 = y.getElementRef(0);
  assert(el8.isConstRef());
  std::cout<<"Accesing an element via const-ref obj "<<el8<<" ( "<<el8.isConstRef()<<" )"<<std::endl;
  assert(el8 == Length{1.5});

  // Getting a const ref to a temporary
  auto el7b = (x + x).getElementRef(0);
  assert(el7b.isConstRef());
  // NOTE: does not mean we can do anything useful with el7b like this, not even read it, but is the same behaviour as regular types

  //Turn ref back into value (a copy)
  Length el4b = el4;
  assert(el4b == Length{2.0});
  Length el9 = el8;
  assert(el9 == Length{1.5});

#ifdef FAIL_DEMO

  el8 = Length{1.0};// el8 is constref, can't write

  //Can't do this - no non-const reference to a temporary (rvalue)
  // Note if we use auto we get a const reference, which is correct - see earlier
  UnitCheckedType<1, 0, 0, STScalarRef<double, false> > el7 = (x + x).getElementRef(0);

  el7b = Length{1.0}; // el7b is constref, can't write

  Length el10 = el7b; //Capturing reference is OK, using the value is not
  assert(el10 == Length{3.0}); // Might "work", depending on optimiser, but is UNDEFINED behaviour

  //Disallow reference access to r-values with get too
  double & el7a = (v+v).get(0);
#endif

  std::cout<<"Internal checks OK\n";
}


void performance_tests(){

  std::cout<<"___________________________________________________________"<<std::endl;
  std::cout<<"Performance tests"<<std::endl;
  Length Arr1[100];
  auto sz = sizeof(Arr1)/sizeof(double);
  if(sz == 100){
    std::cout<<"Case 1: There is no storage overhead, a 100 element array of Lengths is the size of "<<sz<<" doubles"<<std::endl;
  }else{
    std::cout<<"Case 1: There is storage overhead, a 100 element array of Lengths is the size of "<<sz<<" doubles"<<std::endl;
  }
  Position Arr2[100];
  sz = sizeof(Arr2)/(sizeof(double));
  if(sz == 300){
    std::cout<<"Case 2: There is no storage overhead, a 100 element array of Positions is the size of "<<sz<<" doubles"<<std::endl;
  }else{
    std::cout<<"Case 2: There is storage overhead, a 100 element array of Positions is the size of "<<sz<<" doubles"<<std::endl;
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
