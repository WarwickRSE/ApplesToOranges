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

// \TODO Add any other operators to complete set
// \TODO Any other functions we want?

template <SF L, SF M, SF T, typename ST>
class UnitCheckedType{

    // Friend all other unit-checked-types for heterogeneous functions
    template<SF,SF,SF,typename> friend class UnitCheckedType;

    // Verify anything necessary about fractions
#ifdef USE_FRACTIONAL_POWERS
 static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0, "Fraction cannot have zero denom");
#endif

  private:
    ST val;

    static constexpr bool hasNoUnits(){return is_equal(L, 0) && is_equal(M, 0) && is_equal(T, 0);}

  public:
    UnitCheckedType():val(){};

    // Constructors are very permissive - expect ST to restrict to valid values if necessary
    template <typename Tl>
    explicit UnitCheckedType(Tl x):val(x){}


    UnitCheckedType(const UnitCheckedType& src):val(src.val){}
    UnitCheckedType& operator=(const UnitCheckedType& src){
        val=src.val;
        return *this;
    }

    // Scalar initialisation for single value of arithmetic types only
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > > 
    UnitCheckedType& operator=(const num& src){
        val=src;
        return *this;
    }


    // Initializer list type deduction requires that we have
    // a constructor for each layer of nesting we want to accept
    // After 2-3 layers the syntax to initialise this way is so ugly
    // and error prone that there it little point going further
    template <typename Tl>
    UnitCheckedType(std::initializer_list<Tl> l):val(l){}
    template <typename Tl>
    UnitCheckedType(std::initializer_list<std::initializer_list<Tl> > l):val(l){}
    template <typename Tl>
    UnitCheckedType(std::initializer_list<std::initializer_list<std::initializer_list<Tl> > > l):val(l){}

    // Accessors
    // Value extraction for all units 0 only
    template <typename... Args>
    auto& get(Args ... args_in){
      static_assert(hasNoUnits());
      return val.get(args_in...);
    }

    template <typename... Args>
    auto get(Args ... args_in)const{
      static_assert(hasNoUnits());
      return val.get(args_in...);
    }
 
    auto& data(){
      static_assert(hasNoUnits());
      return val;
    }
    auto data()const{
      static_assert(hasNoUnits());
      return val;
    }

    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
    explicit operator num() const{
      static_assert(hasNoUnits());
      return static_cast<num>(val);
    }

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

    template<SF Li, SF Mi, SF Ti, typename STi>
    bool isSameUnits(const UnitCheckedType<Li, Mi, Ti, STi> & other)const{
        return is_equal(L, Li) && is_equal(M, Mi) && is_equal(T, Ti);
    }
    template<SF Li, SF Mi, SF Ti, typename STi>
    bool isSameRank(const UnitCheckedType<Li, Mi, Ti, STi> & other)const{
        return std::is_same_v<ST, STi>;
    }
    template<SF Li, SF Mi, SF Ti, typename STi>
    bool isSameType(const UnitCheckedType<Li, Mi, Ti, STi> & other)const{
        return isSameUnits(other) && isSameRank(other);
    }

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

    // Allow only arithmetic types and exclude char
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    UnitCheckedType operator*=(const num & other){
        val *= other;
        return *this;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedType operator*(const num & lhs, UnitCheckedType rhs){
        return rhs*=lhs;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedType operator*(UnitCheckedType lhs, const num & rhs){
        return lhs*=rhs;
    }

    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    UnitCheckedType operator/=(const num & other){
        val /= other;
        return *this;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedType operator/(UnitCheckedType lhs, const num & rhs){
        return lhs/=rhs;
    }
    //Divide arithmetic type by unitted type
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedType<-L,-M,-T, ST> operator/(const num & lhs, UnitCheckedType rhs){
        UnitCheckedType<-L,-M,-T,ST> tval;
        tval.val = lhs/rhs.val;
        return tval;
    }

    // Exponentiation
    template<SF exp>
    friend UnitCheckedType<L*exp, M*exp, T*exp, ST> pow(const UnitCheckedType & base){
        UnitCheckedType<L*exp, M*exp, T*exp, ST> tval;
        if constexpr(exp.denom == 1){
            tval.val = base.val.pow((long)exp.num);
        }else{
            tval.val = base.val.pow((double)exp.num/exp.denom);
        }
        return tval;
    }

    friend UnitCheckedType<L/2, M/2, T/2, ST> sqrt(const UnitCheckedType & base){
        UnitCheckedType<L/2, M/2, T/2, ST> tval;
        tval.val = base.val.sqrt();
        return tval;
    }

    // Other products - dot, cross, etc \TODO implement others
    template<typename Ts, typename Q=ST>
    using ReturnTypeDot = decltype((std::declval<Q>()).dot(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, typename STi, typename Q=ST>
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeDot<Q,STi> > dot(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeDot<Q,STi> > tval;
        tval.val = this->val.dot(other.val);
        return tval;
    }    
    // Magnitude function
    template <typename Q=ST>
    using ReturnTypeMagnitude = decltype((std::declval<Q>()).magnitude());
    template <typename Q=ST>
    UnitCheckedType<L, M, T, ReturnTypeMagnitude<Q> > magnitude()const{
        UnitCheckedType<L, M, T, ReturnTypeMagnitude<Q> > tval;
        tval.val = val.magnitude();
        return tval;
    }
    // Cross product
    template<typename Ts, typename Q=ST>
    using ReturnTypeCross = decltype((std::declval<Q>()).cross(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, typename STi, typename Q=ST>
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeCross<Q,STi> > cross(const UnitCheckedType<Li, Mi, Ti, STi> &other)const{
        UnitCheckedType<L+Li, M+Mi, T+Ti, ReturnTypeCross<Q,STi> > tval;
        tval.val = this->val.cross(other.val);
        return tval;
    }

    //Normalize function
    template <typename Q=ST>
    void normalize(){
        val.normalize();
    }

    // Comparison operators
    // Implement these for matching units only, but allow different
    // storage types providing they implement a comparison
    // Equality
    template<typename STi>
    bool operator==(const UnitCheckedType<L, M, T, STi> & other)const{
        return val == other.val;
    }
    template<typename STi>
    bool operator!=(const UnitCheckedType<L, M, T, STi> & other)const{
        return val != other.val;
    }
    // Comparisons
    template<typename STi>
    bool operator<(const UnitCheckedType<L, M, T, STi> & other)const{
        return val < other.val;
    }
    template<typename STi>
    bool operator<=(const UnitCheckedType<L, M, T, STi> & other)const{
        return val <= other.val;
    }
    template<typename STi>
    bool operator>(const UnitCheckedType<L, M, T, STi> & other)const{
        return val > other.val;
    }
    template<typename STi>
    bool operator>=(const UnitCheckedType<L, M, T, STi> & other)const{
        return val >= other.val;
    }


};

template <SF L, SF M, SF T, typename ST>
std::ostream& operator<<(std::ostream& os, const UnitCheckedType<L, M, T, ST>& val_in){
  os << val_in.to_string();
  return os;
}

#endif
