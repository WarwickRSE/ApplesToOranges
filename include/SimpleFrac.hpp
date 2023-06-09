
#ifdef USE_FRACTIONAL_POWERS
  // Omit this whole file if we can't use it as Template params
  // We don't have to do this, but they are NOT a generally useable fraction type
#ifndef  __SIMPLEFRAC_HPP__
#define  __SIMPLEFRAC_HPP__

#include <numeric>

typedef long long f_int;

//Simple fraction is a bare-bones fraction class intended for use as a template
// parameter (C++20) value (not a type).
// Because of this, comparison operators MUST be a strong order, where
// 1/2 != 2/4 (even though they are arithmetically equal)
// IF one only every constructs fractions in their simplest terms, and
// uses arithmetic to combine them, they remain simple and this will not
// make any difference for equality. 
// For simplicity, use the arithmetic functions is_equal, is_less and is_greater which
// behave as expected. Note also implements these for comparison to integers
struct SimpleFrac
{
  const f_int num, denom;  
  explicit constexpr SimpleFrac(f_int inum): num(inum), denom(1) {};
  explicit constexpr SimpleFrac(f_int inum, f_int idenom): num(inum), denom(idenom) {};
  constexpr auto operator<=>(const SimpleFrac &)const = default;
};

constexpr bool simplifiable(const SimpleFrac & a){
    return std::gcd(a.num, a.denom) != 1;
};
constexpr SimpleFrac simplify(const SimpleFrac & a){
  const f_int g = std::gcd(a.num, a.denom);
  return g != 0 ? SimpleFrac(a.num/g, a.denom/g): SimpleFrac(a.num, a.denom);
};
// More general is_equal checks if arithmetically equal
// For fully simplified fractions, these are the same operation
constexpr bool is_equal(const SimpleFrac & a, const SimpleFrac & b){
    return a.num*b.denom == b.num*a.denom;
};
constexpr bool is_greater(const SimpleFrac & a, const SimpleFrac & b){
    return a.num*b.denom > b.num*a.denom;
};
constexpr bool is_less(const SimpleFrac & a, const SimpleFrac & b){
    return a.num*b.denom < b.num*a.denom;
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
#endif