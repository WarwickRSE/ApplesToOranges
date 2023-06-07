#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>

#include <helper.hpp>
#include <SimpleFrac.hpp>
#include <StorageTypes.hpp>

using SF = SimpleFrac;
using dblscalar = STScalar<double>;
using dbl3vec = STVector<double, 3>;
using dbl3tens = STTensor<double, 3>;

template <SF L, SF M, SF T, typename ST>
class UnitCheckedType{

    // Verify anything necessary about fractions
    static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0, "Fraction cannot have zero denom");

    // Verify that ST is a valid storage type
    // \TODO

    SF ids[3] = {L,M,T};

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

  template<typename... Ts>
  using WrapTypeMultiply = decltype(operator*(std::declval<Ts>()...))(Ts...);
  template<SF Li, SF Mi, SF Ti, typename STi>
    UnitCheckedType<L+Li, M+Mi, T+Ti, typename std::invoke_result_t<WrapTypeMultiply<ST, STi>, ST, STi> > operator*(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
      UnitCheckedType<L+Li, M+Mi, T+Ti, typename std::invoke_result_t<WrapTypeMultiply<ST, STi>, ST, STi> > tval;
      tval.val = this->val*other.val;
      return tval;
    }
};

template <SF L, SF M, SF T, typename ST>
std::ostream& operator<<(std::ostream& os, const UnitCheckedType<L, M, T, ST>& val_in){
  os << val_in.to_string();
  return os;
}


//\TODO any way to make this more readable?
using UCDouble     = UnitCheckedType<SF{0,1}, SF{0,1}, SF{0,1}, dblscalar>;
using UCDouble3Vec = UnitCheckedType<SF{0,1}, SF{0,1}, SF{0,1}, dbl3vec>;


using ABCs          = UnitCheckedType<SF{1,1}, SF{1,1}, SF{1,1}, dblscalar>;
using ABC          = UnitCheckedType<SF{1,1}, SF{1,1}, SF{1,1}, dbl3vec>;
using ABCt          = UnitCheckedType<SF{1,1}, SF{1,1}, SF{1,1}, dbl3tens>;

#endif