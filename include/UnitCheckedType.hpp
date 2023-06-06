#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>

#include <helper.hpp>
#include <SimpleFrac.hpp>
#include <StorageTypes.hpp>

using SF = SimpleFrac;
using dblscalar = STScalar<double>;
using dbl3vec = ST3Vector<double>;
//using dbl3vec = ST3VectorDouble;

template <SF L, SF M, SF T, class ST>
class UnitCheckedType{

    // Verify anything necessary about fractions
    static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0, "Fraction cannot have zero denom");

    // Verify that ST is a valid storage type
    // \TODO

    // Extract its core type
    typedef typename extract_value_type<ST>::value_type core_type;

    ST val;

    SF ids[3] = {L,M,T};

  public:
    UnitCheckedType():val(0){};

    template<typename Td>
    explicit UnitCheckedType(Td x):val(x){}

    template<typename Td>
    UnitCheckedType(Td x, Td y, Td z):val(x,y,z){}

    // Accessors
    core_type& operator[](size_t i){
        return val[i];
      //return i<val.size()? val[i]: val[0];
    }
    core_type operator[](size_t i)const{
        return val[i];
      //return i<val.size()?val[i]:0;
    }

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

  
  template<SF Li, SF Mi, SF Ti>
    UnitCheckedType<L+Li, M+Mi, T+Ti, ST> operator*(const UnitCheckedType<Li, Mi, Ti, ST> &other)const{
      UnitCheckedType<L+Li, M+Mi, T+Ti, ST> tval;
      for(size_t i=0;i<val.size();++i){
        tval[i]=val[i]*other[i];
      }
      return tval;
    }

};

//\TODO any way to make this more readable?
using UCDouble     = UnitCheckedType<SF{0,1}, SF{0,1}, SF{0,1}, dblscalar>;
using UCDouble3Vec = UnitCheckedType<SF{0,1}, SF{0,1}, SF{0,1}, dbl3vec>;

using ABC          = UnitCheckedType<SF{1,1}, SF{1,1}, SF{1,1}, dbl3vec>;

#endif