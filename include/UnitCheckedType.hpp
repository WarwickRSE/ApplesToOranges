#ifndef UNITCHECKEDTYPE_HPP
#define UNITCHECKEDTYPE_HPP

#include <sstream>
#include <numeric>
#include <initializer_list>
#include <map>
#include <typeindex>

#include "helper.hpp"
#include "SimpleFrac.hpp"
#include "StorageTypes.hpp"

#ifdef USE_FRACTIONAL_POWERS
 using SF = SimpleFrac;
 inline constexpr bool divides(SF a, int b){return true;}///<Stub - fractions are divisible by any integer except 0, caller should not pass 0!
#else
 using SF = int;///<Fallback typedef if USE_FRACTIONAL_POWERS not defined
 inline constexpr bool is_equal(int a, int b){return a==b;}///< Integer equality helper for fallback when USE_FRACTIONAL_POWERS not defined
 inline constexpr bool is_less(int a, int b){return a<b;}///< Integer less-than helper for fallback when USE_FRACTIONAL_POWERS not defined
 inline constexpr bool is_greater(int a, int b){return a>b;}///< Integer greater-than helper for fallback when USE_FRACTIONAL_POWERS not defined
 inline constexpr bool divides(int a, int divisor){return divisor!=0 ? std::gcd(a, divisor)==divisor : false;}///< Integer divisibility helper for fallback when USE_FRACTIONAL_POWERS not defined
#endif

/// \todo Add any other operators to complete set
/// \todo Any other functions we want?

/** @file
 * Physical Units checked types, wrapping some sort of data storage type
 *
 * SI Units are based on expressing the units for any physical quantity as a product of powers of the fundamental units. The fundamental units and their short names are:
 *- M Mass (kg)
 *- L Length (m)
 *- T Time (s)
 *- K Temperature (K)
 *- A Current (A)
 *- MO Amount (mole)
 *- CD Luminous intensity (cd)
*
* NOTE: For clarity, we use "Unitless" to mean "a UnitCheckedTypeFull" for which all exponents are zero, and use phrases such as "numeric type" or "ordinary type" to mean types which are not UnitCheckedTypeFull at all. We also use "Fundamental" to mean one of the SI base units above (exactly one of these has an exponent of exactly one, all others are zero), and "Compound" to describe any other combination or power of these.
*
* NOTE 2: Prefixes (e.g. milli, kilo etc) are not included in this code. Mixing prefixes is risky and can lead to numerical unexpectedness, so we do not implement it.
But if you use, for example, microns for length, and then use microns everywhere, in all compounds, then your results will correctly come out in microns too. See unitNames for how you can change the displayed unit names for this and other purposes.
*/
/*
The following is used for _printing_ the units.
As noted, if you use, for example, microns for length, and then use microns everywhere, in all compounds, then your results will correctly come out in microns too, so you may wish to have microns as the label. You can also modify the names array for localisation.
If you do change this, probably do it once, at the start of a program, and then leave it alone.
*/
/// Extra stuff for display
namespace UnitChecked{
    /// Strings for printing Fundamental names
    inline std::string unitNames[7] = {"kg", "m", "s", "K", "A", "mol", "cd"};
    /// Strings for printing Compound names
    inline std::map<std::type_index, std::string> customUnitStrings;
    template<typename T>
    inline constexpr void registerUnits(std::string_view newUnits){
        /// Register a new unit string for a type
        /** Stores a Units string. If T is a Fundamental unit, this overrides the text in all cases. If T is a Compound unit, this overrides the text only for that exact type, and takes precedence over the Fundamental unit text.
         * @param newUnits The new unit string (likely a string literal)
        */
        if constexpr(T::isFundamentalType()){
          unitNames[T::whichFundamentalType()] = newUnits;
        }else{
          customUnitStrings[std::type_index(typeid(T))] = newUnits;
        }
    }
    template<typename T>
    inline constexpr void registerUnitsForTypeOf(const T & theType, std::string_view newUnits){
        /// Register a new unit string for a type
        /**
         * Uses an instance rather than requiring user specify the template parameter but NOTE this applies to ALL instances, not just the passed one.
         * @param theType An instance of the type to be registered
         * @param newUnits The new unit string (likely a string literal)
        */
        registerUnits<T>(newUnits);
    }
};
/** @brief The basic Unit Checked wrapper type
 *
 * This is the wrapper layer which takes care of the units checking in operations such as equality, addition etc. This type is very permissive and as far as possible deferrs to the underlying storage type for validity of operations.
 *
 * @tparam L Exponent of Length dimension
 * @tparam M Exponent of Mass dimension
 * @tparam T Exponent of Time dimension
 * @tparam K Exponent of Temperature dimension
 * @tparam A Exponent of Current dimension
 * @tparam MO Exponent of AmountOfSubstance dimension
 * @tparam CD Exponent of LuminousIntensity dimension
 * @tparam ST Underlying storage type
 */
template <SF L, SF M, SF T, SF K, SF A, SF MO, SF CD, typename ST>
class UnitCheckedTypeFull{
    /// Friend all other unit-checked-types for heterogeneous functions (NOTE: this seems to be the correct way to do this, but it is not supported by some compilers, so we have a NO_SUPPORT_TEMPLATE_FRIENDSHIP define to allow for this.)
    template<SF,SF,SF,SF,SF,SF,SF, typename> friend class UnitCheckedTypeFull;

    // Verify anything necessary about fractions
#ifdef USE_FRACTIONAL_POWERS
  static_assert(L.denom !=0 && M.denom !=0 && T.denom !=0 && K.denom !=0 && A.denom !=0 && MO.denom !=0 && CD.denom !=0, "Fraction cannot have zero denom");
  static_assert(L.denom > 0 && M.denom > 0 && T.denom > 0 && K.denom > 0 && A.denom > 0 && MO.denom > 0 && CD.denom > 0, "Denominator must be positive");
  static_assert(!simplifiable(L) && !simplifiable(M) && !simplifiable(T) && !simplifiable(K) && !simplifiable(A) && !simplifiable(MO) && !simplifiable(CD),
  "Fraction must be defined in simplest form");
#endif

#ifndef NO_SUPPORT_TEMPLATE_FRIENDSHIP
  private:
#else
  public:
#endif
    ST val;///<Stored value(s)

  public:
/**
 * @name Constructors
 *
 * UnitCheckedTypeFull implements generic constructor wrappers, for the obvious set. Note that these are intended to be permissive, and the underlying storage type is expected to enforce any necessary restrictions.
 *@{
 */
    constexpr UnitCheckedTypeFull():val(){}///<Default constructor

    template <typename Tl>
    constexpr explicit UnitCheckedTypeFull(Tl x):val(x){}///<Prospective single element constructor

    constexpr UnitCheckedTypeFull(const UnitCheckedTypeFull& src):val(src.val){}///<Copy constructor
    ///Copy constructor fake to make a nicer error message if we try and remove the const-ness from a ref
    //NOTE: actually setting val like this is _WRONG_ but with the assert this code will never compile, so we do it to squelch a confusing warning
    template <typename STm, typename=std::enable_if_t<!std::is_same_v<ST, STm> > >
    constexpr UnitCheckedTypeFull(const UnitCheckedTypeFull<L,M,T,K,A,MO,CD, STm> & src):val(src.val){
      if constexpr(STm::is_const_v && !ST::is_const_v){
        static_assert(!STm::is_const_v, "Error: Trying to remove const from a reference!");
      }
    }

    constexpr UnitCheckedTypeFull(const ST& src):val(src){}///<Constructor from storage type value
    constexpr UnitCheckedTypeFull& operator=(const UnitCheckedTypeFull& src){
        val=src.val;
        return *this;
    }///<Copy assignment
    template<typename U, typename=std::enable_if<std::is_convertible_v<U, ST> > >
    const UnitCheckedTypeFull& operator=(const UnitCheckedTypeFull<L,M,T,K,A,MO,CD, U>& src){
        val=src.val;
        return *this;
    }///<Assignment from a related type (different but convertible storage type)

///@}
    // Initializer list type deduction requires that we have
    // a constructor for each layer of nesting we want to accept
    // After 2-3 layers the syntax to initialise this way is so ugly
    // and error prone that there it little point going further
    // Exclude self as initialiser list type so that copy constructor is used instead
    // NOTE: Valid types for list are either base data type of underlying storage
    // or the storage itself (although only one element will be read)
    // OR a unitchecked type of SAME units but any storage type (deferring to storage type for validity checks)
    /// Helper typedef for the (probably numeric) type contained in the storage type
    typedef typename extract_value_type<ST>::value_type ST_t;

    // Use one-element constructor above for init from single ST value, so exclude that specific case here
/**@{
 @name Initialiser list constructors
 *
 * These constructors allow initialisation of the UnitCheckedTypeFull from initialiser lists. We allow up to 3 levels of nesting, but feel the syntax is so ugly and error-prone that there is little point going further - whether or not any particular level and/or type is available, and how it handles various sizes of list, depends on the storage type, except that if the source type is UnitCheckedTypeFull, then units are required to match.
 */
    template<typename Tl, typename std::enable_if_t<std::is_same_v<ST_t, Tl> && !std::is_same_v<ST, Tl>, bool> =false >
    constexpr UnitCheckedTypeFull(std::initializer_list<Tl> l):val(l){}
    template<typename Tl, typename std::enable_if_t<std::is_same_v<ST_t, Tl>, bool> =false >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<Tl> > l):val(l){}
    template<typename Tl, typename std::enable_if_t<std::is_same_v<ST_t, Tl>, bool> =false >
    constexpr UnitCheckedTypeFull(std::initializer_list<std::initializer_list<std::initializer_list<Tl> > > l):val(l){}

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

///@}

/**
 * @name Accessor functions
 *
 * These cover all the possible ways to access the stored values, for reading and writing. We divide "getting" into two versions - a simple "get" for types without units, and an "unsafeGet" which explicitly abandons the unit checking. The operators and operations provided by UnitCheckedTypeFull should cover most use-cases without needing to access the value directly. The main reason to need to do the latter is to pass something to a function which expects a number - if this is a special function (e.g. sin) then the argument ought to be dimensionally unit-less for physical correctness, so the solution is to construct this argument and use "get". If a needed function just doesn't understand UnitCheckedTypeFull (external library code etc), the only option is to use "unsafeGet" and then construct the proper Unitted type from the return.
 * @{
 */

    /// Value setter
    void set(const ST & val_in){
      val = val_in;
    }

    /// Value access for Unitless type
    template <typename... Args>
    auto& get(Args ... args_in)&{
      static_assert(hasNoUnits());
      return val.get(args_in...);
    }
    /// Value access for Unitless type (const ref version)
    template <typename... Args>
    constexpr const auto & get(Args ... args_in)const&{
      static_assert(hasNoUnits());
      return val.get(args_in...);
    }
    /// Value access for Unitless type (const version)
    template <typename... Args>
    constexpr auto get(Args ... args_in)const&&{
      static_assert(hasNoUnits());
      return val.get(args_in...);
    }
    /// Value extraction for Unitless type
    auto& data(){
      static_assert(hasNoUnits());
      return val;
    }
    /// Value extraction for Unitless type (const version)
    auto data()const{
      static_assert(hasNoUnits());
      return val;
    }
    /// Casting to any arithmetic type, for Unitless type
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
    explicit operator num() const{
      static_assert(hasNoUnits());
      return static_cast<num>(val);
    }

    /// Value access preserving units but making a COPY of the value
    template <typename... Args>
    auto getElement(Args ... args_in)const{
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ReturnTypeElement<ST, Args...> > tval;
        tval.set(val.getElement(args_in...));
        return tval;
    }
    /// Reference value access, preserving units. NOTE: only valid for l-values ("real" variables, not temporaries)
    template <typename... Args>
    auto getElementRef(Args ... args_in)&{
        using STm = ReturnTypeElementRef<ST, Args...>;
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STm> tval(val.getElementRef(args_in...));
        return tval;
    }
    template <typename... Args>
    auto getElementRef(Args ... args_in)const&{
        using STm = typename ST_add_const<ReturnTypeElementRef<ST, Args...> >::modified_type;
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STm> tval(val.getElementRef(args_in...));
        return tval;
    }

    /// "Unsafe" value access - ignores units
    template <typename... Args>
    auto& unsafeGet(Args ... args_in){
        return val.get(args_in...);
    }
    /// "Unsafe" value access - ignores units (const version)
    template <typename... Args>
    auto unsafeGet(Args ... args_in)const{
        return val.get(args_in...);
    }

    /// Get a copy of the value, with units gone (Not intended for general use, but needed internally and cannot be private)
    constexpr ST stripUnits()const{
        return val;
    }
///@}

/** @name Functions for Display
 *@{
*/
    std::string to_string()const{
      std::stringstream ss;
      ss<<val;
      return ss.str();
    }///<Stringify using stream operator

    static std::string units(){
      const std::string theUnits = make_unit_str();
      if(UnitChecked::customUnitStrings.count(std::type_index(typeid(UnitCheckedTypeFull)))>0)
        return UnitChecked::customUnitStrings[std::type_index(typeid(UnitCheckedTypeFull))];
      return theUnits;
    }

    static std::string make_unit_str(){
      using namespace UnitChecked;
      std::stringstream ss;
      // Positive exponents first, then negative, ommitted if zero
      if constexpr(is_greater(M, 0)){
        ss<<unitNames[0];
        if constexpr(is_equal(M, 1)) ss<<" ";
        else ss<<"^(" << M << ")";
      }
      if constexpr(is_greater(L, 0)){
        ss<<unitNames[1];
        if constexpr(is_equal(L, 1)) ss<<" ";
        else ss<<"^(" << L << ")";
      }
      if constexpr(is_greater(T, 0)){
        ss<<unitNames[2];
        if constexpr(is_equal(T, 1)) ss<<" ";
        else ss<<"^(" << T << ")";
      }
      if constexpr(is_greater(K, 0)){
        ss<<unitNames[3];
        if constexpr(is_equal(K, 1)) ss<<" ";
        else ss<<"^(" << K << ")";
      }
      if constexpr(is_greater(A, 0)){
        ss<<unitNames[4];
        if constexpr(is_equal(A, 1)) ss<<" ";
        else ss<<"^(" << A << ")";
      }
      if constexpr(is_greater(MO, 0)){
        ss<<unitNames[5];
        if constexpr(is_equal(MO, 1)) ss<<" ";
        else ss<<"^(" << MO << ")";
      }
      if constexpr(is_greater(CD, 0)){
        ss<<unitNames[6];
        if constexpr(is_equal(CD, 1)) ss<<" ";
        else ss<<"^(" << CD << ")";
      }
      if constexpr(is_less(M, 0)) ss<<unitNames[0]<<"^(" << M << ")";
      if constexpr(is_less(L, 0)) ss<<unitNames[1]<<"^(" << L << ")";
      if constexpr(is_less(T, 0)) ss<<unitNames[2]<<"^(" << T << ")";
      if constexpr(is_less(K, 0)) ss<<unitNames[3]<<"^(" << K << ")";
      if constexpr(is_less(A, 0)) ss<<unitNames[4]<<"^(" << A << ")";
      if constexpr(is_less(MO, 0)) ss<<unitNames[5]<<"^(" << MO << ")";
      if constexpr(is_less(CD, 0)) ss<<unitNames[6]<<"^(" << CD << ")";
      return ss.str();
    }///<Produce string for units. Skips any unit with zero exponent, and uses the UnitChecked namespace for unit names
///@}

/** @name Type comparison helpers
 *@{
*/
    static constexpr bool hasNoUnits(){return is_equal(L, 0) && is_equal(M, 0) && is_equal(T, 0) && is_equal(K, 0) && is_equal(A, 0) && is_equal(MO, 0) && is_equal(CD, 0);}///<Check if this type has all Unit exponents zero

    static constexpr bool isFundamentalType(){
        // First check only one is non-zero
        constexpr bool tmp = ((int)(!is_equal(L,0)) + (int)(!is_equal(M,0)) + (int)(!is_equal(T,0)) + (int)(!is_equal(K,0)) + (int)(!is_equal(A,0)) + (int)(!is_equal(MO,0)) + (int)(!is_equal(CD,0)))!=1;
        if constexpr(tmp) return false;
        //Now check that that one is 1
        else return (L+M+T+K+A+MO+CD)==1;
    }///<Check if this type is a fundamental unit type (only one exponent is 1, all others are zero)

    static constexpr size_t whichFundamentalType(){
        //I am sorry, but the compiler parameter order and the usual SI units order do not match and it's too late to fix. These are the order of the strings in unitNames
        return is_equal(L,1)? 1 : is_equal(M,1)? 0 : is_equal(T,1)? 2 : is_equal(K,1)? 3 : is_equal(A,1)? 4 : is_equal(MO,1)? 5 : 6;}///<Get the index of the fundamental unit type (0-6), junk if not fundamental

    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    constexpr bool isSameUnits(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> & other)const{
        return is_equal(L, Li) && is_equal(M, Mi) && is_equal(T, Ti) && is_equal(K, Ki) && is_equal(A, Ai) && is_equal(MO, MOi) && is_equal(CD, CDi);
    }///<Compare units with another type
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    constexpr bool isSameRank(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> & other)const{
        return std::is_same_v<ST, STi>;
    }///<Compare storage types with another type
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    constexpr bool isSameType(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> & other)const{
        return isSameUnits(other) && isSameRank(other);
    }///<Compare Units and StorageType
///@}

/** @name Return typedefs
 *
 * Typedefs used to get the return types of operations implemented by the Storage Type ST, so that we can construct the UnitCheckedFull entity with the correct type to return. NOTE: the template parameter Q is essential so that we do not require any of these operations to be present unless they are used.  \todo Add Q to mult and dvide for consistency
 * @{
 */
    template <typename Q, typename...Args>
    using ReturnTypeElement = decltype(std::declval<Q>().getElement(std::declval<Args>()...));
    template<typename Q, typename...Args>
    using ReturnTypeElementRef = decltype(std::declval<Q>().getElementRef(std::declval<Args>()...));
    template<typename Q=ST>
    using ReturnTypeStripRef = decltype(STStripReference(std::declval<Q>()));
    template<typename Ts>
    using ReturnTypeMultiply = decltype(operator*(std::declval<ST>(), std::declval<Ts>()));
    template<typename Ts>
    using ReturnTypeDivide = decltype(operator/(std::declval<ST>(), std::declval<Ts>()));
    template <typename Q=ST>
    using ReturnTypeMagnitude = decltype((std::declval<Q>()).magnitude());
    template<typename Ts, typename Q=ST>
    using ReturnTypeDot = decltype((std::declval<Q>()).dot(std::declval<Ts>()));
    template<typename Ts, typename Q=ST>
    using ReturnTypeCross = decltype((std::declval<Q>()).cross(std::declval<Ts>()));
    template<typename Ts, typename Q=ST>
    using ReturnTypeOuter = decltype((std::declval<Q>()).outer(std::declval<Ts>()));

///@}

/** @name Arithmetic operators
 *
 * These functions implement basic maths between two UnitCheckedTypeFull types. See also the arithmetic overloads for simple numeric types. Note that only Dimensionally-valid cases are allowed - so for addition and subtraction we must match units, and for multiplication and divison we will produce a value with combined units (note this means *= and /= are invalid unless the other type has no units). \todo Implement *= and /= for other with no units
@{
*/
    /// Unary minus
    UnitCheckedTypeFull operator-()const{
        UnitCheckedTypeFull tval;
        tval.val = -val;
        return tval;
    }

    /*Lets use Addition as the model describing what is going on here.
    For +=, we just modify in place and return the reference, as always
    But to avoid implementing all the operators on reference types, we avoid the $= operators and use the binary ones, which will allow the casting to convert away from Reference types. 
    Since ST should really be "numerical" we entirely expect += and  + to do the same thing, so we only lose a little efficiency by introducing the additional equals operation (I am not sure how far the optimiser will fix this) and we save a lot of code duplication.
    For binary addition, we have to pass both arguments by const ref (can no longer use the nice pass-by-value-to-maximise-copy-elision) because a copy wrong in the reference case. Instead we pass by const ref and construct the return type by stripping off any reference-hood.
    We also use is_convertible_v to restrict to storage types which can convert to ours, although this is mostly to try and get nicer errors, because the storage type wont implement addition in that case.
    */
    template<typename U, typename=std::enable_if<std::is_convertible_v<U, ST> > >
    UnitCheckedTypeFull operator+=(const UnitCheckedTypeFull<L,M,T,K,A,MO,CD,U> & other){
        val = val + other.val;
        return *this;
    }
    template<typename U, typename=std::enable_if<std::is_convertible_v<U, ST> > >
    friend UnitCheckedTypeFull<L,M,T,K,A,MO,CD, ReturnTypeStripRef<> > operator+(const UnitCheckedTypeFull & lhs, const UnitCheckedTypeFull<L,M,T,K,A,MO,CD,U> & other){
        //Case where StripRef has an effect - implies it's a ref...
        if constexpr(!std::is_same_v<ST, ReturnTypeStripRef<> >){
          UnitCheckedTypeFull<L,M,T,K,A,MO,CD, ReturnTypeStripRef<> > tval;
          tval.val = STStripReference(lhs.val) + other.val;
          return tval;
        }else{
          UnitCheckedTypeFull tval;
          tval.val = lhs.val + other.val;
          return tval;
        }
    }
    template<typename U, typename=std::enable_if<std::is_convertible_v<U, ST> > >
    UnitCheckedTypeFull operator-=(const UnitCheckedTypeFull<L,M,T,K,A,MO,CD,U> & other){
        val = val - other.val;
        return *this;
    }
    template<typename U, typename=std::enable_if<std::is_convertible_v<U, ST> > >
    friend UnitCheckedTypeFull<L,M,T,K,A,MO,CD, ReturnTypeStripRef<> > operator-(const UnitCheckedTypeFull & lhs, const UnitCheckedTypeFull<L,M,T,K,A,MO,CD,U> & other){
        if constexpr(!std::is_same_v<ST, ReturnTypeStripRef<> >){
          UnitCheckedTypeFull<L,M,T,K,A,MO,CD, ReturnTypeStripRef<> > tval;
          tval.val = STStripReference(lhs.val) - other.val;
          return tval;
        }else{
          UnitCheckedTypeFull tval;
          tval.val = lhs.val - other.val;
          return tval;
        }
    }
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    auto operator*(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
      if constexpr(!std::is_same_v<ST, ReturnTypeStripRef<> >){
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeMultiply<ReturnTypeStripRef<STi> > > tval;
        tval.val = STStripReference(val)*other.val;
        return tval;
      }else{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeMultiply<STi> > tval;
        tval.val = val*other.val;
        return tval;
      }
    }
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi>
    auto operator/(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
      if constexpr(!std::is_same_v<ST, ReturnTypeStripRef<> >){
        UnitCheckedTypeFull<L-Li, M-Mi, T-Ti, K-Ki, A-Ai, MO-MOi, CD-CDi, ReturnTypeMultiply<ReturnTypeStripRef<STi> > > tval;
        tval.val = STStripReference(val)/other.val;
        return tval;
      }else{
        UnitCheckedTypeFull<L-Li, M-Mi, T-Ti, K-Ki, A-Ai, MO-MOi, CD-CDi, ReturnTypeDivide<STi> > tval;
        tval.val = val/other.val;
        return tval;
      }
    }
///@}

/** @name Arithmetic operators with numeric types
 *
 * These operators allow multiplication and division of the UnitCheckedTypeFull by an arithmetic type. The arithmetic type must be numeric, and char is excluded.
 * @{
*/
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
///@}

/** @name Powers and Roots
 *
 * Unit checking occurs at compile time, so we cannot have a generic power function for a type with units. We provide the most common (sqrt and cbrt), a genral version, and a runtime version which can only accept types with no units.
 * As always, the results depend on what the underlying StorageType supplies.
 * @{
 */
    // Exponentiation (note that pre-20 C++ does not support non-integral types in function templates as well as class templates, so we rely on our having defined SF as int in that case)
    /// Exponentiation for arbitrary exponent
    template<SF exp>
    friend auto pow(const UnitCheckedTypeFull & base){
        UnitCheckedTypeFull<L*exp, M*exp, T*exp, K*exp, A*exp, MO*exp, CD*exp,ST> tval;
        tval.val = base.val.pow((double)exp);
        return tval;
    }
    /// Square root
    friend auto sqrt(const UnitCheckedTypeFull & base){
        static_assert(divides(L, 2) && divides(M,2) && divides(T,2) && divides(K,2) && divides(A,2) && divides(MO,2) && divides(CD,2));
        UnitCheckedTypeFull<L/2, M/2, T/2, K/2, A/2, MO/2, CD/2,ST> tval;
        tval.val = base.val.sqrt();
        return tval;
    }
    ///Cube root
    friend auto cbrt(const UnitCheckedTypeFull & base){
        static_assert(divides(L, 3) && divides(M,3) && divides(T,3) && divides(K,3) && divides(A,3) && divides(MO,3) && divides(CD,3));
        UnitCheckedTypeFull<L/3, M/3, T/3, K/3, A/3, MO/3, CD/3, ST> tval;
        tval.val = base.val.cbrt();
        return tval;
    }
    ///Nth root
    template <long exp>
    friend auto nthrt(const UnitCheckedTypeFull & base){
        static_assert(divides(L, exp) && divides(M,exp) && divides(T,exp) && divides(K,exp) && divides(A,exp) && divides(MO,exp) && divides(CD,exp));
        UnitCheckedTypeFull<L/exp, M/exp, T/exp, K/exp, A/exp, MO/exp, CD/exp, ST> tval;
        tval.val = base.val.pow(1.0/exp);
        return tval;
    }
    /// Exponentiation for unitless types - runtime allowable
    template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
    friend auto pow(const UnitCheckedTypeFull & base, const num & exp){
        static_assert(hasNoUnits());
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ST> tval;
        tval.val = base.val.pow(exp);
        return tval;
    }
///@}

/** @name Linear Algebra
 *
 * These functions supply the basic Linear Algebra operations such as magnitude, transpose and various products. Units are handled here, but the actual operations are, as always, deferred to the underlying StorageType, and are available ONLY if this considers that they make sense. If extending the set, see the Typedefs for how we ensure that we only require the operations calling code uses.
 * @{
 */
    /// Dot product (Unitted like Multiplication)
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
    auto dot(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi, ReturnTypeDot<Q,STi> > tval;
        tval.val = this->val.dot(other.val);
        return tval;
    }    
    /// Cross product (Unitted like Multiplication)
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
    auto cross(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi,ReturnTypeCross<Q,STi> > tval;
        tval.val = this->val.cross(other.val);
        return tval;
    }

    /// Outer product (Unitted like Multiplication)
    template<SF Li, SF Mi, SF Ti, SF Ki, SF Ai, SF MOi, SF CDi, typename STi, typename Q=ST>
        auto outer(const UnitCheckedTypeFull<Li, Mi, Ti, Ki, Ai, MOi, CDi, STi> &other)const{
        UnitCheckedTypeFull<L+Li, M+Mi, T+Ti, K+Ki, A+Ai, MO+MOi, CD+CDi,ReturnTypeOuter<Q,STi> > tval;
        tval.val = this->val.outer(other.val);
        return tval;
    }

    /// Magnitude function (Same units as input)
    template <typename Q=ST>
    auto magnitude()const{
        UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ReturnTypeMagnitude<Q> > tval;
        tval.val = val.magnitude();
        return tval;
    }

    ///In-place normalize function
    template <typename Q=ST>
    void normalize(){
        val.normalize();
    }
    ///Normalized function - returns normalized copy
    template <typename Q=ST>
    UnitCheckedTypeFull normalized()const{
        UnitCheckedTypeFull tval;
        tval.val = val;
        tval.val.normalize();
        return tval;
    }

    ///Transpose function
    template <typename Q=ST>
    auto transpose()const{
        UnitCheckedTypeFull<M, L, T, K, A, MO, CD, Q> tval;
        tval.val = val.transpose();
        return tval;
    }
///@}

/**
 * @name Comparison operators
 *
 * Comparisons are implemented only for matching units (Apples To Oranges, the namesake of this code), but as always, defer to the StorageType for the actual comparison, and allow anything which is implemented.
 @{
 */
    template<typename STi>
    bool operator==(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val == other.val;
    }
    template<typename STi>
    bool operator!=(const UnitCheckedTypeFull<L, M, T, K, A, MO, CD, STi> & other)const{
        return val != other.val;
    }
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
///@}

  friend std::ostream& operator<<(std::ostream& os, const UnitCheckedTypeFull& val_in){
    os << val_in.to_string();
    return os;
  }

  friend std::istream& operator>>(std::istream& is, UnitCheckedTypeFull& val_in){
    is >> val_in.val;
    return is;
  }
};

///Helper typedef to get the underlying storage type of a UnitCheckedTypeFull
template <typename U>
using WrappedType = decltype(std::declval<U>().stripUnits());
/**
 * @brief Function to make an identity entity with specified units
 *
 * @tparam U - The type of the entity to be created. Note this cannot be inferred, so must be given explicitly
 * @return The Identity Entity with given Units for given StorageType
 */
template <typename U>
U makeIdentity(){
  U tval;
  tval.set(WrappedType<U>::identity());
  return tval;
}

/** @brief Template magic to identify is_unitchecked_type at compile time
 *
 * This primary template matches anything and defaults to false for the ::value field and defers to is_arithmetic for ::numeric field
*/
template<typename T>
struct is_unitchecked{
    static constexpr bool value = false;
    static constexpr bool numeric = std::is_arithmetic_v<T>;
};

/** @brief Template magic to identify is_unitchecked_type at compile time
 *
 * This specialisation only matches UnitCheckedTypeFull, giving true for the ::value field and defers to is_arithmetic applied to the underlying numeric value for the ::numeric field
*/
template <SF L, SF M, SF T, SF K, SF A, SF MO, SF CD, typename ST>
struct is_unitchecked<UnitCheckedTypeFull<L, M, T, K, A, MO, CD, ST> >{
    static constexpr bool value = true;
    typedef typename extract_value_type<ST>::value_type core_type;
    static constexpr bool numeric = std::is_arithmetic_v<core_type>;
};

/// Shortcut for is_unitchecked value field
template<typename T>
inline constexpr bool is_unitchecked_v = is_unitchecked<T>::value;
/// Shortcut for is_unitchecked numeric field
template<typename T>
inline constexpr bool is_unitchecked_numeric_v = is_unitchecked<T>::numeric;

// Pre C++20, need compiler to find these as template functions with right name, although any type will do
///Stub - needed for template friend function name resolution pre c++20
template <typename T> void pow();
///Stub - needed for template friend function name resolution pre c++20
template <typename T> void nthrt();

// Note - all checked types interoperate, these typedefs are just shortcuts
///"Dynamics" shortcut specialisation - Length, Mass, Time
template<SF L, SF M, SF T, typename ST>
using UnitCheckedTypeDynamic = UnitCheckedTypeFull<L, M, T, 0, 0, 0, 0, ST>;
/// Default shortcut specialisation - matches Dynamics
template<SF L, SF M, SF T, typename ST>
using UnitCheckedType = UnitCheckedTypeDynamic<L, M, T, ST>;

/// Dynamics plus temperature shortcut specialisation
template<SF L, SF M, SF T, SF K, typename ST>
using UnitCheckedTypeThermodynamic = UnitCheckedTypeFull<L, M, T, K, 0, 0, 0, ST>;

/// Electrodynamics shortcut specialisation - Length, Mass, Time, Temperature, Charge
template<SF L, SF M, SF T, SF K, SF A, typename ST>
using UnitCheckedTypeElectrodynamic = UnitCheckedTypeFull<L, M, T, K, A, 0, 0, ST>;


#endif
