#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>

#include <helper.hpp>
#include <SimpleFrac.hpp>
#include <StorageTypes.hpp>

using SF = SimpleFrac;

template <SF L, SF M, SF T, typename ST>
class UnitCheckedType{

    // Verify anything necessary about fractions
    static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0, "Fraction cannot have zero denom");

    // Verify that ST is a valid storage type
    // \TODO

  public:
    ST val; //\TODO Should we try to protect val from direct modification?

    UnitCheckedType():val(0){};

    // Constructors are very permissive - expect ST to restrict to valid values if necessary
    template <typename Tl>
    explicit UnitCheckedType(Tl x):val(x){}

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
    auto& operator[](size_t i){
        return val[i];
    }
    auto operator[](size_t i)const{
        return val[i];
    }

    template <typename... Args>
    auto& get(Args ... args_in){
        return val.get(args_in...);
    }
    template <typename... Args>
    auto get(Args ... args_in)const{
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

    template<typename... Ts>
    using WrapTypeMultiply = decltype(operator*(std::declval<Ts>()...))(Ts...);
    template<SF Li, SF Mi, SF Ti, typename STi>
        UnitCheckedType<L+Li, M+Mi, T+Ti, typename std::invoke_result_t<WrapTypeMultiply<ST, STi>, ST, STi> > operator*(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L+Li, M+Mi, T+Ti, typename std::invoke_result_t<WrapTypeMultiply<ST, STi>, ST, STi> > tval;
        tval.val = this->val*other.val;
        return tval;
    }
    template<typename... Ts>
    using WrapTypeDivide = decltype(operator/(std::declval<Ts>()...))(Ts...);
    template<SF Li, SF Mi, SF Ti, typename STi>
        UnitCheckedType<L+Li, M+Mi, T+Ti, typename std::invoke_result_t<WrapTypeDivide<ST, STi>, ST, STi> > operator/(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L+Li, M+Mi, T+Ti, typename std::invoke_result_t<WrapTypeMultiply<ST, STi>, ST, STi> > tval;
        tval.val = this->val/other.val;
        return tval;
    }

    // Other products - dot, cross, etc \TODO implement

    // Magnitude/ norm

    auto magnitude()const{
        return val.magnitude();
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