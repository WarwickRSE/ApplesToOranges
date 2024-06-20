#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>
#include <numeric>
#include <initializer_list>

#include "helper.hpp"
#include "SimpleFrac.hpp"
#include "StorageTypes.hpp"

#ifdef USE_FRACTIONAL_POWERS
 using SF = SimpleFrac;
 inline constexpr bool divides(SF a, int b){return true;}// Stub
#else
 using SF = int;
 inline constexpr bool is_equal(int a, int b){return a==b;}
 inline constexpr bool is_less(int a, int b){return a<b;}
 inline constexpr bool is_greater(int a, int b){return a>b;}
 inline constexpr bool divides(int a, int divisor){return divisor!=0 ? std::gcd(a, divisor)==divisor : false;}
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

The following is used for _printing_ the units.
Prefixes (e.g. milli, kilo etc) are not included in this code. Mixing prefixes is risky and can lead to numerical unexpectedness, so we do not implement it.
But if you use, for example, microns for length, then use microns everywhere, in all compounds, and your results will correctly come out in microns too, so you may wish to have microns as the label.
Or modify this array for localisation.
If you do change this, probably do it once, at the start, and then leave it alone.
*/
/*Unit names*/
namespace UnitChecked{
    inline std::string unitNames[7] = {"m", "kg", "s", "K", "A", "mol", "cd"};
};
template <SF L, SF M, SF T, SF K, SF A, SF MO, SF CD, typename ST>
class UnitCheckedTypeFull{
    // Friend all other unit-checked-types for heterogeneous functions
    template<SF,SF,SF,SF,SF,SF,SF, typename> friend class UnitCheckedTypeFull;

    // Verify anything necessary about fractions
#ifdef USE_FRACTIONAL_POWERS
 static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0 && K.denom !=0 && A.denom !=0 && MO.denom !=0 && CD.denom !=0, "Fraction cannot have zero denom");
#endif

#ifndef NO_SUPPORT_TEMPLATE_FRIENDSHIP
  private:
#else
  public:
#endif
    ST val;

  private:
    static constexpr bool hasNoUnits(){return is_equal(L, 0) && is_equal(M, 0) && is_equal(T, 0) && is_equal(K, 0) && is_equal(A, 0) && is_equal(MO, 0) && is_equal(CD, 0);}

  public:
    constexpr UnitCheckedTypeFull():val(){};

    // Constructors are very permissive - expect ST to restrict to valid values if necessary
    template <typename Tl>
    constexpr explicit UnitCheckedTypeFull(Tl x):val(x){}


    constexpr UnitCheckedTypeFull(const UnitCheckedTypeFull& src):val(src.val){}
    constexpr UnitCheckedTypeFull(const ST& src):val(src){}
    constexpr UnitCheckedTypeFull& operator=(const UnitCheckedTypeFull& src){
        val=src.val;
        return *this;
    }

    // Initializer list type deduction requires that we have
    // a constructor for each layer of nesting we want to accept
    // After 2-3 layers the syntax to initialise this way is so ugly
    // and error prone that there it little point going further
    // Exclude self as initialiser list type so that copy constructor is used instead
    // NOTE: Valid types for list are either base data type of underlying storage
    // or the storage itself (although only one element will be read)
    // OR a unitchecked type of SAME units but any storage type (deferring to storage type for validity checks)
    typedef typename extract_value_type<ST>::value_type ST_t;

    // Use one-element constructor above for init from single ST value, so exclude that specific case here
    template<typename Tl, typename std::enable_if_t<std::is_same_v<ST_t, Tl> && !std::is_same_v<ST, Tl>, bool> =false >
    constexpr UnitCheckedTypeFull(std::initializer_list<Tl> l):val(l){}
    template<typename Tl, typename std::enable_if_t<std::is_same_v<ST_t, Tl>, bool> =false >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<Tl> > l):val(l){}
    template<typename Tl, typename std::enable_if_t<std::is_same_v<ST_t, Tl>, bool> =false >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<std::initializer_list<Tl> > > l):val(l){}

    // Get a copy of the value, with units gone
    constexpr ST stripUnits()const{
        return val;
    }

    // Overload to give nicer error message
    // In single nest case we exclude both self (use copy contructor) and type ST (use single element constructor). In neither case is more than one element meaningful
    template<typename Tl, typename std::enable_if_t<!std::is_same_v<ST_t, Tl> && !std::is_same_v<ST, Tl> && !std::is_same_v<UnitCheckedTypeFull, Tl>, bool> =true >
    constexpr UnitCheckedTypeFull(std::initializer_list<Tl> l){
        // Condition always false BUT this is only known after substitution of Tl
        static_assert(std::is_same_v<ST_t, Tl>, "Initialiser list type must match storage data-type or units");
    }
    template<typename Tl, typename std::enable_if_t<!std::is_same_v<ST_t, Tl>, bool> =true >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<Tl> > l){
        static_assert(std::is_same_v<ST_t, Tl>, "Initialiser list type must match storage data-type or units");
    }
    template<typename Tl, typename std::enable_if_t<!std::is_same_v<ST_t, Tl>, bool> =true >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<std::initializer_list<Tl> > > l){
        static_assert(std::is_same_v<ST_t, Tl>, "Initialiser list type must match storage data-type or units");
    }

    template <typename STi, typename=std::enable_if_t<!std::is_same_v<ST, STi> > >
    constexpr UnitCheckedTypeFull(std::initializer_list< UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> > l):val(l){}
    template <typename STi, typename=std::enable_if_t<!std::is_same_v<ST, STi> > >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list< UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> > > l):val(l){}
    template <typename STi, typename=std::enable_if_t<!std::is_same_v<ST, STi> > >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<std::initializer_list<UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> > > > l):val(l){}

    // Accessors

    // Value setter
    void set(const ST & val_in){
      val = val_in;
    }

    // Value extraction for all units 0 only
    template <typename... Args>
    auto& get(Args ... args_in){
      static_assert(hasNoUnits());
      return val.get(args_in...);
    }

    template <typename... Args>
    constexpr auto get(Args ... args_in)const{
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
      using namespace UnitChecked;
      std::stringstream ss;
      if constexpr(is_greater(M, 0) && ! is_equal(M, 1)) ss<<unitNames[1]<<"^(" << M << ")";
      if constexpr(is_greater(L, 0) && ! is_equal(L, 1)) ss<<unitNames[0]<<"^(" << L << ")";
      if constexpr(is_greater(T, 0) && ! is_equal(T, 1)) ss<<unitNames[2]<<"^(" << T << ")";
      if constexpr(is_greater(K, 0) && ! is_equal(K, 1)) ss<<unitNames[3]<<"^(" << K << ")";
      if constexpr(is_greater(A, 0) && ! is_equal(A, 1)) ss<<unitNames[4]<<"^(" << A << ")";
      if constexpr(is_greater(MO, 0) && ! is_equal(MO, 1)) ss<<unitNames[5]<<"^(" << MO << ")";
      if constexpr(is_greater(CD, 0) && ! is_equal(CD, 1)) ss<<unitNames[6]<<"^(" << CD << ")";
      if constexpr(is_equal(M, 1)) ss<<unitNames[1]<<" ";
      if constexpr(is_equal(L, 1)) ss<<unitNames[0]<<" ";
      if constexpr(is_equal(T, 1)) ss<<unitNames[2]<<" ";
      if constexpr(is_equal(K, 1)) ss<<unitNames[3]<<" ";
      if constexpr(is_equal(A, 1)) ss<<unitNames[4]<<" ";
      if constexpr(is_equal(MO, 1)) ss<<unitNames[5]<<" ";
      if constexpr(is_equal(CD, 1)) ss<<unitNames[6]<<" ";
      if constexpr(is_less(M, 0)) ss<<unitNames[1]<<"^(" << M << ")";
      if constexpr(is_less(L, 0)) ss<<unitNames[0]<<"^(" << L << ")";
      if constexpr(is_less(T, 0)) ss<<unitNames[2]<<"^(" << T << ")";
      if constexpr(is_less(K, 0)) ss<<unitNames[3]<<"^(" << K << ")";
      if constexpr(is_less(A, 0)) ss<<unitNames[4]<<"^(" << A << ")";
      if constexpr(is_less(MO, 0)) ss<<unitNames[5]<<"^(" << MO << ")";
      if constexpr(is_less(CD, 0)) ss<<unitNames[6]<<"^(" << CD << ")";
      return ss.str();
    };

    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    bool isSameUnits(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> & other)const{
        return is_equal(L, Li) && is_equal(M, Mi) && is_equal(T, Ti) && is_equal(K, Ki) && is_equal(A, Ai) && is_equal(MO, MOi) && is_equal(CD, CDi);
    }
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    bool isSameRank(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> & other)const{
        return std::is_same_v<ST, STi>;
    }
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    bool isSameType(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> & other)const{
        return isSameUnits(other) && isSameRank(other);
    }

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
    auto operator*(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeMultiply<STi> > tval;
        tval.val = this->val*other.val;
        return tval;
    }
    template<typename Ts>
    using ReturnTypeDivide = decltype(operator/(std::declval<ST>(), std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    auto operator/(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
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
    friend auto operator/(const num & lhs, UnitCheckedTypeFull rhs){
        UnitCheckedTypeFull<-L,-M,-T,-K,-A,-MO,-CD,ST> tval;
        tval.val = lhs/rhs.val;
        return tval;
    }

    // Exponentiation (note that pre-20 C++ does not support non-integral types so we rely on our defining SF as int in that case)
    template<SF exp>
    friend auto pow(const UnitCheckedTypeFull & base){
        UnitCheckedTypeFull<L*exp, M*exp, T*exp, K*exp, A*exp, MO*exp, CD*exp,ST> tval;
        tval.val = base.val.pow((double)exp);
        return tval;
    }

    friend auto sqrt(const UnitCheckedTypeFull & base){
        static_assert(divides(L, 2) && divides(M,2) && divides(T,2) && divides(K,2) && divides(A,2) && divides(MO,2) && divides(CD,2));
        UnitCheckedTypeFull<L/2, M/2, T/2, K/2, A/2, MO/2, CD/2,ST> tval;
        tval.val = base.val.sqrt();
        return tval;
    }
    friend auto cbrt(const UnitCheckedTypeFull & base){
        static_assert(divides(L, 3) && divides(M,3) && divides(T,3) && divides(K,3) && divides(A,3) && divides(MO,3) && divides(CD,3));
        UnitCheckedTypeFull<L/3, M/3, T/3, K/3, A/3, MO/3, CD/3, ST> tval;
        tval.val = base.val.cbrt();
        return tval;
    }
    template <long exp>
    friend auto nthrt(const UnitCheckedTypeFull & base){
        static_assert(divides(L, exp) && divides(M,exp) && divides(T,exp) && divides(K,exp) && divides(A,exp) && divides(MO,exp) && divides(CD,exp));
        UnitCheckedTypeFull<L/exp, M/exp, T/exp, K/exp, A/exp, MO/exp, CD/exp, ST> tval;
        tval.val = base.val.pow(1.0/exp);
        return tval;
    }
    // Exponentiation for unitless types - runtime allowable
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
    friend auto pow(const UnitCheckedTypeFull & base, const num & exp){
        static_assert(hasNoUnits());
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ST> tval;
        tval.val = base.val.pow(exp);
        return tval;
    }

    // Other products and linear algebra functions - magnitude, normalize, transpose, dot, cross, etc
    // Again, storage types are expected to define these IFF applicable
    template<typename Ts, typename Q=ST>
    using ReturnTypeDot = decltype((std::declval<Q>()).dot(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
    auto dot(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeDot<Q,STi> > tval;
        tval.val = this->val.dot(other.val);
        return tval;
    }    
    // Magnitude function
    template <typename Q=ST>
    using ReturnTypeMagnitude = decltype((std::declval<Q>()).magnitude());
    template <typename Q=ST>
    auto magnitude()const{
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ReturnTypeMagnitude<Q> > tval;
        tval.val = val.magnitude();
        return tval;
    }
    // Cross product
    template<typename Ts, typename Q=ST>
    using ReturnTypeCross = decltype((std::declval<Q>()).cross(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
    auto cross(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi,ReturnTypeCross<Q,STi> > tval;
        tval.val = this->val.cross(other.val);
        return tval;
    }

    //Outer product
    template<typename Ts, typename Q=ST>
    using ReturnTypeOuter = decltype((std::declval<Q>()).outer(std::declval<Ts>()));
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
        auto outer(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi,ReturnTypeOuter<Q,STi> > tval;
        tval.val = this->val.outer(other.val);
        return tval;
    }

    //Normalize function
    template <typename Q=ST>
    void normalize(){
        val.normalize();
    }
    template <typename Q=ST>
    UnitCheckedTypeFull normalized()const{
        UnitCheckedTypeFull tval;
        tval.val = val;
        tval.val.normalize();
        return tval;
    }

    //Transpose
    template <typename Q=ST>
    auto transpose()const{
        UnitCheckedTypeFull<M, L, T, K, A, MO, CD, Q> tval;
        tval.val = val.transpose();
        return tval;
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


  friend std::ostream& operator<<(std::ostream& os, const UnitCheckedTypeFull& val_in){
    os << val_in.to_string();
    return os;
  }

  friend std::istream& operator>>(std::istream& is, UnitCheckedTypeFull& val_in){
    is >> val_in.val;
    return is;
  }
};

template <typename U>
using WrappedType = decltype(std::declval<U>().stripUnits());
template <typename U>
U makeIdentity(){
  U tval;
  tval.set(WrappedType<U>::identity());
  return tval;
}

// Template magic to identify is_unitchecked_type
// Primary template that matches anything and defaults to false
template<typename T>
struct is_unitchecked{
    static constexpr bool value = false;
    static constexpr bool numeric = std::is_arithmetic_v<T>;
};

// Specialization for myclass that sets the value to true
template <SF L, SF M, SF T, SF K, SF A, SF MO, SF CD, typename ST>
struct is_unitchecked<UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ST> >{
    static constexpr bool value = true;
    typedef typename extract_value_type<ST>::value_type core_type;
    static constexpr bool numeric = std::is_arithmetic_v<core_type>;
};

// Helper variable templates
template<typename T>
inline constexpr bool is_unitchecked_v = is_unitchecked<T>::value;

template<typename T>
inline constexpr bool is_unitchecked_numeric_v = is_unitchecked<T>::numeric;

// Needed for template friend function name resolution pre c++20
// Just need compiler to find template function with right name, any type
template <typename T> void pow();
template <typename T> void nthrt();

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
