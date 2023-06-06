#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>

#include <SimpleFrac.hpp>

using SF = SimpleFrac;

template <SF L, SF M, SF T, int dim>
class UnitCheckedType{
  static_assert(dim == 1 || dim == 3, "Checking only supports scalar or 3D vector types"); // Not necessary but do this for now
  private:

    static const size_t dims=dim;
    double val[dims];
    SF ids[3] = {L,M,T};


  public:
    UnitCheckedType(){
      for(size_t i=0;i<dim;++i){
        val[i]=0;
      }
    }
    explicit UnitCheckedType(double x){
      static_assert(dims == 1, "Cannot construct vector from single value"); 
      val[0]=x;
    }
    UnitCheckedType(double x, double y, double z){
      static_assert(dims == 3, "Scalar type cannot be constructed from 3 elements"); 
      val[0]=x;val[1]=y;val[2]=z;
    }

    // Accessors
    double& operator[](size_t i){
      return i<dims? val[i]: val[0];
    }
    double operator[](size_t i)const{
      return i<dims?val[i]:0;
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
    UnitCheckedType<L+Li, M+Mi, T+Ti, dim> operator*(const UnitCheckedType<Li, Mi, Ti, dim> &other)const{
      UnitCheckedType<L+Li, M+Mi, T+Ti, dim> tval;
      for(size_t i=0;i<dim;++i){
        tval[i]=this->val[i]*other[i];
      }
      return tval;
    }

};


#endif