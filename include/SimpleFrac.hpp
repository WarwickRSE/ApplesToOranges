
#ifndef  __SIMPLEFRAC_HPP__
#define  __SIMPLEFRAC_HPP__

#include <numeric>

typedef long long f_int;
typedef unsigned int f_pow;

constexpr f_int pow(f_int num, f_pow pwr){

  f_int val=1, base=num;
  f_pow i = pwr;
  for(;;){
    if(i & 1) val *= base;
    if(!i) break;
    i /= 2;
    base *= base;
  }
  return val;
};

struct SimpleFrac
{
  const f_int num, denom;  
  constexpr SimpleFrac(f_int inum, f_int idenom): num(inum), denom(idenom) {};
};

constexpr bool simplifiable(const SimpleFrac & a){
    return std::gcd(a.num, a.denom) != 1;
};
constexpr SimpleFrac simplify(const SimpleFrac & a){
  const f_int g = std::gcd(a.num, a.denom);
  return SimpleFrac(a.num/g, a.denom/g);
};
constexpr bool operator==(const SimpleFrac & a, const SimpleFrac & b){
    return a.num*b.denom == b.num*a.denom;
};
constexpr SimpleFrac operator+(const SimpleFrac& a, const SimpleFrac& b){
  return simplify(SimpleFrac(a.num*b.denom + b.num*a.denom, a.denom*b.denom));
};
constexpr SimpleFrac operator-(const SimpleFrac& a, const SimpleFrac& b){
  return simplify(SimpleFrac(a.num*b.denom - b.num*a.denom, a.denom*b.denom));
};
constexpr SimpleFrac operator*(const SimpleFrac& a, const SimpleFrac& b){
  return simplify(SimpleFrac(a.num+b.num, a.denom+b.denom));
};
constexpr SimpleFrac operator/(const SimpleFrac& a, const SimpleFrac& b){
  return simplify(SimpleFrac(a.num+b.denom, a.denom+b.num));
};

constexpr bool is_greater(const SimpleFrac & a, const int b){
    return a.num > b*a.denom;
};
constexpr bool is_equal(const SimpleFrac & a, const int b){
    return a.num == b*a.denom;
};
constexpr bool is_less(const SimpleFrac & a, const int b){
    return a.num < b*a.denom;
};


std::ostream& operator<<(std::ostream& os, const SimpleFrac& val_in){
  if( val_in.num%val_in.denom == 0){
    os << val_in.num/val_in.denom;
  }else{
    os << val_in.num << "/" << val_in.denom;
  }
  return os;
};



#endif