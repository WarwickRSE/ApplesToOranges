#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>

#include <helper.hpp>
#include <SimpleFrac.hpp>
#include <StorageTypes.hpp>

#ifdef USE_FRACTIONAL_POWERS
 using SF = SimpleFrac;
#else
 using SF = int;
 inline constexpr bool is_equal(int a, int b){return a==b;}
 inline constexpr bool is_less(int a, int b){return a<b;}
 inline constexpr bool is_greater(int a, int b){return a>b;}
#endif


template <SF L, SF M, SF T, typename ST>
class UnitCheckedType{

    // Verify anything necessary about fractions
#ifdef USE_FRACTIONAL_POWERS
 static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0, "Fraction cannot have zero denom");
#endif

    // Verify that ST is a valid storage type
    // \TODO

  public:
    ST val; //\TODO Should we try to protect val from direct modification?

    UnitCheckedType():val(0){};

    // Constructors are very permissive - expect ST to restrict to valid values if necessary
    template <typename Tl>
    explicit UnitCheckedType(Tl x):val(x){}


    UnitCheckedType(const UnitCheckedType& src):val(src.val){}
    UnitCheckedType& operator=(const UnitCheckedType& src){
        val=src.val;
        return *this;
    }

    // Only arithmetic types for this
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > > 
    UnitCheckedType& operator=(const num& src){
        val=src;
        return *this;
    }// Scalar equals for arithmetic types


    // Initializer list type deduction requires that we have
    // a constructor for each layer of nesting we want to accept
    // After 2-3 layers the syntax is so ugly and error prone that
    // there it little point going further
    template <typename Tl>
    UnitCheckedType(std::initializer_list<Tl> l):val(l){}
    template <typename Tl>
    UnitCheckedType(std::initializer_list<std::initializer_list<Tl> > l):val(l){}
    template <typename Tl>
    UnitCheckedType(std::initializer_list<std::initializer_list<std::initializer_list<Tl> > > l):val(l){}

    // Accessors
    // unsafeGet - no unit checks, just return bare value. Use with CAUTION
    template <typename... Args>
    auto& unsafeGet(Args ... args_in){
        return val.get(args_in...);
    }
    template <typename... Args>
    auto unsafeGet(Args ... args_in)const{
        return val.get(args_in...);
    }

    std::string to_string()const{
      std::stringstream ss;
      ss<<val;
      return ss.str();
    };

    std::string units()const{
      std::stringstream ss;
      if constexpr(is_greater(M, 0) && ! is_equal(M, 1)) ss<<"kg^(" << M << ")";
      if constexpr(is_greater(L, 0) && ! is_equal(L, 1)) ss<<"m^(" << L << ")";
      if constexpr(is_greater(T, 0) && ! is_equal(T, 1)) ss<<"s^(" << T << ")";
      if constexpr(is_equal(M, 1)) ss<<"kg";
      if constexpr(is_equal(L, 1)) ss<<"m";
      if constexpr(is_equal(T, 1)) ss<<"s";
      if constexpr(is_less(M, 0)) ss<<"kg^(" << M << ")";
      if constexpr(is_less(L, 0)) ss<<"m^(" << L << ")";
      if constexpr(is_less(T, 0)) ss<<"s^(" << T << ")";
      return ss.str();
    };


    // Operators
    // Unary minus
    UnitCheckedType operator-()const{
        UnitCheckedType tval;
        tval.val = -val;
        return tval;
    }

    // Addition/subtraction for same type and dims only
    UnitCheckedType operator+=(const UnitCheckedType & other){
        val += other.val;
        return *this;
    }
    friend UnitCheckedType operator+(UnitCheckedType lhs, const UnitCheckedType & other){
        return lhs+=other;
    }
    UnitCheckedType operator-=(const UnitCheckedType & other){
        val -= other.val;
        return *this;
    }
    friend UnitCheckedType operator-(UnitCheckedType lhs, const UnitCheckedType & other){
        return lhs-=other;
    }

    template<typename Ts>
    using ReturnTypeMultiply = decltype(operator*(std::declval<ST>(), std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, typename STi>
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeMultiply<STi> > operator*(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeMultiply<STi> > tval;
        tval.val = this->val*other.val;
        return tval;
    }
    template<typename Ts>
    using ReturnTypeDivide = decltype(operator/(std::declval<ST>(), std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, typename STi>
        UnitCheckedType<L-Li, M-Mi, T-Ti, ReturnTypeDivide<STi> > operator/(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L-Li, M-Mi, T-Ti, ReturnTypeDivide<STi> > tval;
        tval.val = this->val/other.val;
        return tval;
    }

    // Numeric multiply/divide
    template<typename num>
    UnitCheckedType operator*=(const num & other){
        val *= other;
        return *this;
    }
    // Only arithmetic types for this
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > > 
    friend UnitCheckedType operator*(const num & lhs, UnitCheckedType rhs){
        static_assert(std::is_arithmetic_v<num>, "Can only multiply by arithmetic types");
        return rhs*=lhs;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > > 
    friend UnitCheckedType operator*(UnitCheckedType lhs, const num & rhs){
        static_assert(std::is_arithmetic_v<num>, "Can only multiply by arithmetic types");
        return lhs*=rhs;
    }
    template<typename num>
    UnitCheckedType operator/=(const num & other){
        val /= other;
        return *this;
    }
    // Only arithmetic types for this
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > > 
    friend UnitCheckedType operator/(UnitCheckedType lhs, const num & rhs){
        return lhs/=rhs;
    }
    //Divide arithmetic type by unitted type - remember units!
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
    friend UnitCheckedType<-L,-M,-T, ST> operator/(const num & lhs, UnitCheckedType rhs){
        UnitCheckedType<-L,-M,-T,ST> tval;
        tval.val = lhs/rhs.val;
        return tval;
    }

    // Other products - dot, cross, etc \TODO implement
    // \TODO currently requires ST implement this - should not force this
    template<typename Ts>
    using ReturnTypeDot = decltype((std::declval<ST>()).dot(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, typename STi>
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeDot<STi> > dot(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeDot<STi> > tval;
        tval.val = this->val.dot(other.val);
        return tval;
    }    
    // Magnitude function
    using ReturnTypeMagnitude = decltype((std::declval<ST>()).magnitude());
    UnitCheckedType<L, M, T, ReturnTypeMagnitude > magnitude()const{
        UnitCheckedType<L, M, T, ReturnTypeMagnitude > tval;
        tval.val = val.magnitude(); // \TODO should use plain = here and throughout operators?
        return tval;
    }

    // Comparison operators
    // Implement these for matching units only, but allow different
    // storage types providing they implement a comparison
    // Equality
    template<typename STi>
    friend bool operator==(const UnitCheckedType<L, M, T, ST> & first, const UnitCheckedType<L, M, T, STi> & other){
        return first.val == other.val;
    }
    template<typename STi>
    friend bool operator!=(const UnitCheckedType<L, M, T, ST> & first, const UnitCheckedType<L, M, T, STi> & other){
        return first.val != other.val;
    }
    // Comparisons
    template<typename STi>
    friend bool operator<(const UnitCheckedType<L, M, T, ST> & first, const UnitCheckedType<L, M, T, STi> & other){
        return first.val < other.val;
    }
    template<typename STi>
    friend bool operator<=(const UnitCheckedType<L, M, T, ST> & first, const UnitCheckedType<L, M, T, STi> & other){
        return first.val <= other.val;
    }
    template<typename STi>
    friend bool operator>(const UnitCheckedType<L, M, T, ST> & first, const UnitCheckedType<L, M, T, STi> & other){
        return first.val > other.val;
    }
    template<typename STi>
    friend bool operator>=(const UnitCheckedType<L, M, T, ST> & first, const UnitCheckedType<L, M, T, STi> & other){
        return first.val >= other.val;
    }


};

template <SF L, SF M, SF T, typename ST>
std::ostream& operator<<(std::ostream& os, const UnitCheckedType<L, M, T, ST>& val_in){
  os << val_in.to_string();
  return os;
}

#endif