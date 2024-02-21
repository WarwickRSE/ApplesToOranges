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

/* Units are in order with short names:
L Length (m)
M Mass (kg)
T Time (s)
K Temperature (K)
A Current (A)
MO Amount (mole)
CD Luminous intensity (cd)
*/
template <SF L, SF M, SF T, SF K, SF A, SF MO, SF CD, typename ST>
class UnitCheckedTypeFull{
    // Friend all other unit-checked-types for heterogeneous functions
    template<SF,SF,SF,SF,SF,SF,SF, typename> friend class UnitCheckedTypeFull;

    // Verify anything necessary about fractions
#ifdef USE_FRACTIONAL_POWERS
 static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0 && K.denom !=0 && A.denom !=0 && MO.denom !=0 && CD.denom !=0, "Fraction cannot have zero denom");
#endif

  private:
    ST val;

    static constexpr bool hasNoUnits(){return is_equal(L, 0) && is_equal(M, 0) && is_equal(T, 0) && is_equal(K, 0) && is_equal(A, 0) && is_equal(MO, 0) && is_equal(CD, 0);}

  public:
    UnitCheckedTypeFull():val(){};

    // Constructors are very permissive - expect ST to restrict to valid values if necessary
    template <typename Tl>
    explicit UnitCheckedTypeFull(Tl x):val(x){}


    UnitCheckedTypeFull(const UnitCheckedTypeFull& src):val(src.val){}
    UnitCheckedTypeFull& operator=(const UnitCheckedTypeFull& src){
        val=src.val;
        return *this;
    }

    // Scalar initialisation for single value of arithmetic types only
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > > 
    UnitCheckedTypeFull& operator=(const num& src){
        val=src;
        return *this;
    }


    // Initializer list type deduction requires that we have
    // a constructor for each layer of nesting we want to accept
    // After 2-3 layers the syntax to initialise this way is so ugly
    // and error prone that there it little point going further
    template <typename Tl>
    UnitCheckedTypeFull(std::initializer_list<Tl> l):val(l){}
    template <typename Tl>
    UnitCheckedTypeFull(std::initializer_list<std::initializer_list<Tl> > l):val(l){}
    template <typename Tl>
    UnitCheckedTypeFull(std::initializer_list<std::initializer_list<std::initializer_list<Tl> > > l):val(l){}

    // Accessors
    // Value extraction for all units 0 only
    template <typename... Args>
    auto& get(Args ... args_in){
      static_assert(hasNoUnits());
      if constexpr(sizeof...(args_in) == 0){
        return val;
      }else{
        return val.get(args_in...);
      }
    }

    template <typename... Args>
    auto get(Args ... args_in)const{
      static_assert(hasNoUnits());
      if constexpr(sizeof...(args_in) == 0){
        return val;
      }else{
        return val.get(args_in...);
      }
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
      if constexpr(is_greater(K, 0) && ! is_equal(K, 1)) ss<<"K^(" << K << ")";
      if constexpr(is_greater(A, 0) && ! is_equal(A, 1)) ss<<"A^(" << A << ")";
      if constexpr(is_greater(MO, 0) && ! is_equal(MO, 1)) ss<<"mol^(" << MO << ")";
      if constexpr(is_greater(CD, 0) && ! is_equal(CD, 1)) ss<<"cd^(" << CD << ")";
      if constexpr(is_equal(M, 1)) ss<<"kg";
      if constexpr(is_equal(L, 1)) ss<<"m";
      if constexpr(is_equal(T, 1)) ss<<"s";
      if constexpr(is_equal(K, 1)) ss<<"K";
      if constexpr(is_equal(A, 1)) ss<<"A";
      if constexpr(is_equal(MO, 1)) ss<<"mol";
      if constexpr(is_equal(CD, 1)) ss<<"cd";
      if constexpr(is_less(M, 0)) ss<<"kg^(" << M << ")";
      if constexpr(is_less(L, 0)) ss<<"m^(" << L << ")";
      if constexpr(is_less(T, 0)) ss<<"s^(" << T << ")";
      if constexpr(is_less(K, 0)) ss<<"K^(" << K << ")";
      if constexpr(is_less(A, 0)) ss<<"A^(" << A << ")";
      if constexpr(is_less(MO, 0)) ss<<"mol^(" << MO << ")";
      if constexpr(is_less(CD, 0)) ss<<"cd^(" << CD << ")";
      return ss.str();
    };


    // Operators
    // Unary minus
    UnitCheckedTypeFull operator-()const{
        UnitCheckedTypeFull tval;
        tval.val = -val;
        return tval;
    }

    // Addition/subtraction for same type and dims only
    UnitCheckedTypeFull operator+=(const UnitCheckedTypeFull & other){
        val += other.val;
        return *this;
    }
    friend UnitCheckedTypeFull operator+(UnitCheckedTypeFull lhs, const UnitCheckedTypeFull & other){
        return lhs+=other;
    }
    UnitCheckedTypeFull operator-=(const UnitCheckedTypeFull & other){
        val -= other.val;
        return *this;
    }
    friend UnitCheckedTypeFull operator-(UnitCheckedTypeFull lhs, const UnitCheckedTypeFull & other){
        return lhs-=other;
    }

    template<typename Ts>
    using ReturnTypeMultiply = decltype(operator*(std::declval<ST>(), std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeMultiply<STi> > operator*(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeMultiply<STi> > tval;
        tval.val = this->val*other.val;
        return tval;
    }
    template<typename Ts>
    using ReturnTypeDivide = decltype(operator/(std::declval<ST>(), std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
        UnitCheckedTypeFull<L-Li, M-Mi, T-Ti, K-Ki, A-Ai, MO-MOi, CD-CDi, ReturnTypeDivide<STi> > operator/(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L-Li, M-Mi, T-Ti, K-Ki, A-Ai, MO-MOi, CD-CDi, ReturnTypeDivide<STi> > tval;
        tval.val = this->val/other.val;
        return tval;
    }

    // Numeric multiply/divide

    // Allow only arithmetic types and exclude char
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    UnitCheckedTypeFull operator*=(const num & other){
        val *= other;
        return *this;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedTypeFull operator*(const num & lhs, UnitCheckedTypeFull rhs){
        return rhs*=lhs;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedTypeFull operator*(UnitCheckedTypeFull lhs, const num & rhs){
        return lhs*=rhs;
    }

    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    UnitCheckedTypeFull operator/=(const num & other){
        val /= other;
        return *this;
    }
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedTypeFull operator/(UnitCheckedTypeFull lhs, const num & rhs){
        return lhs/=rhs;
    }
    //Divide arithmetic type by unitted type
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> && !std::is_same_v<num,char> > >
    friend UnitCheckedTypeFull<-L,-M,-T, -K, -A, -MO, -CD,ST> operator/(const num & lhs, UnitCheckedTypeFull rhs){
        UnitCheckedTypeFull<-L,-M,-T,-K,-A,-MO,-CD,ST> tval;
        tval.val = lhs/rhs.val;
        return tval;
    }

    // Exponentiation
    template<SF exp>
    friend UnitCheckedTypeFull<L*exp, M*exp, T*exp, K*exp, A*exp, MO*exp, CD*exp, ST> pow(const UnitCheckedTypeFull & base){
        UnitCheckedTypeFull<L*exp, M*exp, T*exp, K*exp, A*exp, MO*exp, CD*exp,ST> tval;
        if constexpr(exp.denom == 1){
            tval.val = base.val.pow((long)exp.num);
        }else{
            tval.val = base.val.pow((double)exp.num/exp.denom);
        }
        return tval;
    }

    friend UnitCheckedTypeFull<L/2, M/2, T/2, K/2, A/2, MO/2, CD/2, ST> sqrt(const UnitCheckedTypeFull & base){
        UnitCheckedTypeFull<L/2, M/2, T/2, K/2, A/2, MO/2, CD/2,ST> tval;
        tval.val = base.val.sqrt();
        return tval;
    }

    // Other products - dot, cross, etc \TODO implement others
    template<typename Ts, typename Q=ST>
    using ReturnTypeDot = decltype((std::declval<Q>()).dot(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeDot<Q,STi> > dot(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeDot<Q,STi> > tval;
        tval.val = this->val.dot(other.val);
        return tval;
    }    
    // Magnitude function
    template <typename Q=ST>
    using ReturnTypeMagnitude = decltype((std::declval<Q>()).magnitude());
    template <typename Q=ST>
    UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ReturnTypeMagnitude<Q> > magnitude()const{
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ReturnTypeMagnitude<Q> > tval;
        tval.val = val.magnitude();
        return tval;
    }
    // Cross product
    template<typename Ts, typename Q=ST>
    using ReturnTypeCross = decltype((std::declval<Q>()).cross(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi,ReturnTypeCross<Q,STi> > cross(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi,ReturnTypeCross<Q,STi> > tval;
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
    bool operator==(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val == other.val;
    }
    template<typename STi>
    bool operator!=(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val != other.val;
    }
    // Comparisons
    template<typename STi>
    bool operator<(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val < other.val;
    }
    template<typename STi>
    bool operator<=(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val <= other.val;
    }
    template<typename STi>
    bool operator>(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val > other.val;
    }
    template<typename STi>
    bool operator>=(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val >= other.val;
    }


};

template <SF L, SF M, SF T, SF K, SF A, SF MO, SF CD, typename ST>
std::ostream& operator<<(std::ostream& os, const UnitCheckedTypeFull<L, M, T, K, A, MO, CD,ST>& val_in){
  os << val_in.to_string();
  return os;
}

// Note - all checked types interoperate
// Default is "dynamics only" - i.e. Length Mass and Time
template<SF L, SF M, SF T, typename ST>
using UnitCheckedTypeDynamic = UnitCheckedTypeFull<L, M, T, 0, 0, 0, 0, ST>;
template<SF L, SF M, SF T, typename ST>
using UnitCheckedType = UnitCheckedTypeDynamic<L, M, T, ST>;

// Dynamics plus temperature
template<SF L, SF M, SF T, SF K, typename ST>
using UnitCheckedTypeThermodynamic = UnitCheckedTypeFull<L, M, T, K, 0, 0, 0, ST>;

// Electrodynamics - Length, Mass, Time, Temperature, Charge
template<SF L, SF M, SF T, SF K, SF A, typename ST>
using UnitCheckedTypeElectrodynamic = UnitCheckedTypeFull<L, M, T, K, A, 0, 0, ST>;


#endif
