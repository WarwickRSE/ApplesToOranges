#ifndef __STORAGE_TYPES_HPP__
#define __STORAGE_TYPES_HPP__

#include <algorithm>
#include <cstring>
#include <cmath>

/** @file
 * Basic data types and Linear Algebra type operations on them */

// Valid Storage Types probably implement (UnitChecking will attempt to use if a user does)
// As many of the following as make sense:
  // Default constructor
  // One element constructor (all values set)
  // Initializer list constructor handling any depth of nesting required
// Copy, move etc constructors
// 'get' function taking as many arguments as wanted, returning a reference to value, and const variant returning copy
// Optional [] operator if 1D access makes sense and is desired
// << stream operator
// Standard arithmetic operators for self binary ops
// Comparison ops as desired
// Optional arithmetic operators for scalar ops
// Heterogenous ops for any desired interactions (such as scalar-vector or vector-tensor)

#ifdef DEBUG
/** @brief A dummy class which implements nothing
     *
     * Used for development to make sure we do not make demands on storage types, by instantiating a UnitChecked type using this, i.e. `UnitCheckedType<0, 0, 0, STDummy> dummy;`
*/
class STDummy{
  public:
    STDummy()=default;
    int get()const{return 0;}
};
#endif

/** @brief A Scalar type
   *
   * I.e. a single value
   * @tparam T The underlying numeric type of the value
*/
template <typename T>
class STScalar{
    private:
        T val{};///<The value

    public:
        ///Default constructor
        constexpr STScalar(){}
        ///Single element constructor
        constexpr STScalar(T val_in):val(val_in){}
        ///Initializer list constructor (assumes one element long, or zero if empty)
        constexpr STScalar(std::initializer_list<T> l){
          l.size()>0? val=*(l.begin()):val=0;
        }
        ///Copy constructor
        constexpr STScalar(const STScalar &a) = default;
        STScalar operator=(const STScalar &a){
          val=a.val;
          return *this;
        }///<Copy assignment

        constexpr T& get(){
            return val;
        }///<Get the value (by reference)
        constexpr T get()const{
            return val;
        }///<Get the value (by copy)

        constexpr auto getElement()const{
            return *this;
        }///< Get element (as a valid StorageType)

        constexpr explicit operator T() const{return val;}///<Cast to stored type

        template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
        explicit operator num() const{
#ifdef NO_NARROWING_CONVERSIONS
          return num{val};
# else
          return static_cast<num>(val);
 #endif
        }///<Cast to any arithmetic type (if possible)

        static STScalar<T> identity(){
          return STScalar<T>{1};
        }///<Identity entity

        STScalar<T> magnitude()const{
          return std::abs(val);
        }///<Magnitude of the value

        STScalar operator-()const{
          return -val;
        }///<Unary minus
        STScalar operator +=(const STScalar & other){
          val+=other.val;
          return *this;
        }
        friend STScalar operator+(const STScalar & a, const STScalar & b){
          return a.val+b.val;
        }
        STScalar operator-=(const STScalar & other){
          val-=other.val;
          return *this;
        }
        friend STScalar operator-(const STScalar & a, const STScalar & b){
          return a.val-b.val;
        }
        STScalar operator *=(const STScalar & other){
            val *= other.val;
            return *this;
        }
        friend STScalar operator*(STScalar lhs, const STScalar & other){
            return lhs*=other;
        }
        STScalar operator/=(const STScalar & other){
            val /= other.val;
            return *this;
        }
        friend STScalar operator/(STScalar lhs, const STScalar & other){
            return lhs/=other;
        }

        STScalar sqrt()const{
          return std::sqrt(val);
        }///<Square root
        STScalar cbrt()const{
          return std::cbrt(val);
        }///<Cube root
        STScalar pow(double p)const{
          return std::pow(val, p);
        }///<Power for real exponent

        //Comparisons
        friend bool operator==(const STScalar & first, const STScalar & other){
          return first.val==other.val;
        }
        friend bool operator!=(const STScalar & first, const STScalar & other){
          return first.val!=other.val;
        }
        friend bool operator<(const STScalar & first, const STScalar & other){
          return first.val<other.val;
        }
        friend bool operator>(const STScalar & first, const STScalar & other){
          return first.val>other.val;
        }
        friend bool operator<=(const STScalar & first, const STScalar & other){
          return first.val<=other.val;
        }
        friend bool operator>=(const STScalar & first, const STScalar & other){
          return first.val>=other.val;
        }
};
///Stream output for Scalars
template <typename T>
std::ostream& operator<<(std::ostream& os, const STScalar<T>& val_in){
  os <<val_in.get();
  return os;
};
///Stream input for Scalars
template <typename T>
std::istream& operator>>(std::istream& is, STScalar<T>& val_in){
  is >>val_in.get();
  return is;
};

/// A special scalar reference type
/* This is needed for the LHS of reference assignments etc, but mostly
silently converts to a scalar in other operations
*/
template <typename T>
class STScalarRef{
  public:
    T &val;
    STScalarRef(T &val_in):val(val_in){}
    STScalarRef(const STScalarRef &a) = default;

    STScalarRef operator=(const STScalarRef &a){
      val=a.val;
      return *this;
    }
    template<typename U>
    STScalarRef operator=(const U &a){
      val=a.get();
      return *this;
    }
    operator STScalar<T>()const{
      return STScalar<T>(val);
    }
};

///Generic function to "strip reference" from a type - no-op for most Storage types
template<typename ST>
ST STStripReference(const ST& a){
  return a;
}
///Specialisation for STScalarRef decaying to STScalar
template<typename T>
STScalar<T> STStripReference(const STScalarRef<T> &a){
  return STScalar<T>(a);
}


//Forward declare for use in outer product
template <typename T, int dim>
class STTensor;

/** @brief A Vector (in linear algebra sense) type
   *
   * Fixed length 1-D array of values
   * @tparam T The underlying numeric type of the values
   * @tparam dim The number of elements in the vector
*/
template <typename T, int dim>
class STVector{
    private:
        T val[dim]{};///< The values

        T normSq()const{
          /// \internal Square-and-add all elements
          T sum=0;
          for(size_t i = 0; i<dim; i++){
            sum+=val[i]*val[i];
          }
          return sum;
        }
    public:
        ///Default constructor
        constexpr STVector(){};
        ///Single element constructor (all values set)
        constexpr STVector(T val_in){for(size_t i = 0; i<dim; i++){val[i]=val_in;}};
        ///Copy constructor
        constexpr STVector(const STVector &a) = default;
        ///Initializer list constructor from numeric type
        constexpr STVector(std::initializer_list<T> l){
          /** Either one element long (all values set same), or for any other length up to dim values are set, any not supplied are zero, any excess are ignored
          */
            const size_t ct=std::min((int)l.size(), dim);
            if(ct == 1){
              for(size_t i = 0; i < dim; i++){val[i]=*(l.begin());}
            }else{
              for(size_t i = 0; i<ct; i++){val[i]=*(l.begin()+i);}
              for(size_t i = ct; i<dim; i++){val[i]=0;} // Zero any excess values
            }
        }
        /// Initializer list constructor from STScalars of same numeric type
        constexpr STVector(std::initializer_list< STScalar<T> > l){
           /** \sa STVector(std::initializer_list<T> l) for behaviour of various lengths of list*/
            const size_t ct=std::min((int)l.size(), dim);
            if(ct == 1){
              for(size_t i = 0; i < dim; i++){val[i]=(l.begin())->get();}
            }else{
              for(size_t i = 0; i<ct; i++){val[i]=(l.begin()+i)->get();}
              for(size_t i = ct; i<dim; i++){val[i]=0;} // Zero any excess values
            }
        }

        // Allow initialisation from wrapped STScalars, where wrapper is assumed to have a stripUnits method to get this back
        /// Typedef to get storage type from arbitrary UnitChecked type
        template <typename U>
        using WrappedType = decltype(std::declval<U>().stripUnits());
        /// Initializer list constructor from doubly nested list of wrapped STScalars
        /** This is needed to allow constructing a Vector of UnitChecked types from a list of Scalar UnitChecked types. Unit checking is done BEFORE this point is reached
         * \sa STVector(std::initializer_list< STScalar<T> > l) for behaviour of various lengths of list
         * \param l The list of lists of wrapped STScalars
         * \tparam U Type of list elements, which must be a Wrapped STScalar<T> for same T as this Vector. That is, when we apply stripUnits to U we must get STScalar<T>
        */
        template <typename U, typename=std::enable_if_t<std::is_same_v<WrappedType<U>, STScalar<T> > > >
        constexpr STVector(std::initializer_list<U> l){
            const size_t ct=std::min((int)l.size(), dim);
            if(ct == 1){
              for(size_t i = 0; i < dim; i++){val[i]=(l.begin())->stripUnits().get();}
            }else{
              for(size_t i = 0; i<ct; i++){val[i]=(l.begin()+i)->stripUnits().get();}
              for(size_t i = ct; i<dim; i++){val[i]=0;} // Zero any excess values
            }
        }

        ///Copy assignment
        STVector operator=(const STVector &a){
          for(size_t i = 0; i<dim; i++){
            val[i]=a[i];
          }
          return *this;
        }

        /// Standard [] operator by reference
        constexpr T& operator[](size_t i){
            return val[i];
        }
        /// Standard [] operator by copy
        constexpr T operator[](size_t i)const{
            return val[i];
        }
        /// Get by reference
        constexpr T& get(size_t i){
            return val[i];
        }
        /// Get by copy
        constexpr T get(size_t i)const{
            return val[i];
        }

        /// Get element (as a valid StorageType)
        constexpr auto getElement(size_t i)const{
            return STScalar<T>(val[i]);
        }
        /// Get reference to element (as a valid reference StorageType)
        constexpr STScalarRef<T> getElementRef(size_t i){
            return STScalarRef<T>(val[i]);
        }

        /// Identity entity
        static STVector<T,dim> identity(){
          return STVector<T, dim>{1, 1, 1};
        }

        /// Magnitude of the vector (square root of sum of squares)
        STScalar<T> magnitude()const{
          return STScalar<T>{std::sqrt(normSq())};
        }

        /// Unary minus
        STVector operator-()const{
          STVector out;
          for(size_t i = 0; i<dim; i++){
            out[i]=-val[i];
          }
          return out;
        }
        STVector operator +=(const STVector & other){
          for(size_t i = 0; i<dim; i++){
            val[i]+=other[i];
          }
          return *this;
        }
        friend STVector operator+(STVector lhs, const STVector & other){
          return lhs+=other;
        }
        STVector operator -=(const STVector & other){
          for(size_t i = 0; i<dim; i++){
            val[i]-=other[i];
          }
          return *this;
        }
        friend STVector operator-(STVector lhs, const STVector & other){
          return lhs-=other;
        }

        STVector operator *=(const STVector & other){
            for(size_t i = 0; i<dim; i++){
                val[i] *= other.val[i];
            }
            return *this;
        }
        friend STVector operator*(STVector lhs, const STVector & other){
            return lhs*=other;
        }
        STVector operator/=(const STVector & other){
            for(size_t i = 0; i<dim; i++){
                val[i] /= other.val[i];
            }
            return *this;
        }
        friend STVector operator/(STVector lhs, const STVector & other){
            return lhs/=other;
        }

        /// Linear algebra - standard dot product A dot B
        STScalar<T> dot(const STVector & other)const{
          STScalar<T> sum=0;
          for(size_t i = 0; i<dim; i++){
            sum+=val[i]*other.val[i];
          }
          return sum;
        }
        /// Linear algebra - standard cross product A cross B, well-defined only for 3-vectors
        STVector<T,dim> cross(const STVector & other) const{
          static_assert(dim ==3); // Well defined only for dim 3 vector
          STVector<T,dim> out;
          for(size_t i = 0; i<dim; i++){
            out[i] = val[(i+1)%dim]*other.val[(i+2)%dim] - val[(i+2)%dim]*other.val[(i+1)%dim];
          }
          return out;
        }
        /// Linear algebra - standard outer product A outer B
        STTensor<T, dim> outer(const STVector & other) const{
          STTensor<T, dim> out;
          for(size_t i = 0; i<dim; i++){
            for(size_t j = 0; j<dim; j++){
              out[i*dim+j] = val[i]*other.val[j];
            }
          }
          return out;
        }
        /// In-place normalization
        void normalize(){
          const STScalar<T> mag = magnitude();
          for(size_t i = 0; i<dim; i++){
            val[i]/=mag.get();
          }
        }

        // Comparisons - in terms of ordering of the norm only
        /// Equality element wise
        friend bool operator==(const STVector & first, const STVector & other){
          for(size_t i = 0; i<dim; i++){
            if(first.val[i]!=other.val[i]) return false;
          }
          return true;
        }
        /// Less-than using norm
        friend bool operator<(const STVector & first, const STVector & other){
          return first.normSq()<other.normSq();
        }

        /// Equality with scalars - compares magnitudes
        friend bool operator==(const STVector & a, const STScalar<T> & b){
          return a.normSq()==b.magnitude()*b.magnitude();
        }
        friend bool operator==(const STScalar<T> & a, const STVector & b){
          return b==a;
        }
        /// Less-than with scalars - compares magnitudes
        friend bool operator<(const STVector & a, const STScalar<T> & b){
          return a.normSq()<b.magnitude()*b.magnitude();
        }
        friend bool operator<(const STScalar<T> & a, const STVector & b){
          return a.magnitude()*a.magnitude() < b.normSq();
        }

        // Implement the rest in terms of < for easy tweaking/ expansion
        template<typename T1, typename T2, typename=std::enable_if_t<std::is_same_v<T1, STVector> || std::is_same_v<T2, STVector>, int> >
        friend bool operator!=(const T1 & first, const T2 & other){
          return !(first == other);
        }
        template<typename T1, typename T2, typename=std::enable_if_t<std::is_same_v<T1, STVector> || std::is_same_v<T2, STVector>, int> >
        friend bool operator>(const T1 & first, const T2 & other){
          return other < first;
        }
        template<typename T1, typename T2, typename=std::enable_if_t<std::is_same_v<T1, STVector> || std::is_same_v<T2, STVector>, int> >
        friend bool operator<=(const T1 & first, const T2 & other){
          return !(first > other);
        }
        template<typename T1, typename T2, typename=std::enable_if_t<std::is_same_v<T1, STVector> || std::is_same_v<T2, STVector>, int> >
        friend bool operator>=(const T1 & first, const T2 & other){
          return !(first < other);
        }

};

///Stream output for Vectors
/** Outputs the vector as a comma-separated list in brackets, e.g. (0.0, 1.0, 2.0)
*/
template <typename T, int dim>
std::ostream& operator<<(std::ostream& os, const STVector<T, dim>& val_in){
  os << "(";
  for(size_t i = 0; i<dim; i++){
    os << val_in[i];
    if(i<dim-1) os << ',';
  }
  os<<")";
  return os;
};
///Stream input for Vectors
/** Reads a vector from stream. Expects the format output by operator<< so e.g. (1.0, 2.0, 3.0) */
template <typename T, int dim>
std::istream& operator>>(std::istream& is, STVector<T, dim>& val_in){
  for(size_t i = 0; i<dim; i++){
    if(!is.good()) return is;
    char seps;
    is>>seps; // Read and toss the bracket and commas
    if(seps != '(' && seps != ',') return is;
    is >>val_in.get(i);
  }
  char seps_last;
  is>>seps_last; // Read and toss the closing ')'
  return is;
};

/**
 * @brief A Tensor type
 *
 * An (ideally small) square matrix of fixed size
 *
 * @tparam T The underlying data type
 * @tparam dim The size (dim x dim) of the matrix
 */
template <typename T, int dim>
class STTensor{
    private:
        T val[dim*dim]{};///<The values

    public:
        constexpr STTensor(){};
        constexpr STTensor(T val_in){for(size_t i = 0; i<dim*dim; i++){val[i]=val_in;}};
        constexpr STTensor(const STTensor &a) = default;
        constexpr STTensor operator=(const STTensor &a){
          for(size_t i = 0; i<dim*dim; i++){
            val[i]=a[i];
          }
          return *this;
        }

        /// Initialisation from a list, singly or doubly nested, of values of type T
        /** For a single list of length 1, we set all values. For a single list of length >1 we set as many elements as are given in 1-D ordering. For a nested list, we set as many values as are given in as many rows as given. All non-set elements are explicitly zeroed.
        */
        template<typename Tl, typename std::enable_if_t<std::is_same_v<Tl, T> || std::is_same_v<Tl, std::initializer_list<T> >, int > =0 >
        constexpr STTensor(std::initializer_list<Tl> l){
            const size_t ct=std::min((int)l.size(), dim);
            if constexpr (std::is_same_v<Tl, std::initializer_list<T>>){
                for(size_t i = 0; i<ct; i++){
                    auto l2 = *(l.begin()+i);
                    const size_t ct2=std::min((int)l2.size(), dim);
                    for(size_t j=0; j<ct2; j++){
                      // Assign as many elements as are given
                      val[i*dim+j]= *(l2.begin()+j);
                    }
                    for(size_t j = ct2; j<dim; j++){
                        val[i*dim+j]=0;
                    }
                }
                for(size_t i = ct; i<dim; i++){
                    for(size_t j = 0; j<dim; j++){
                        val[i*dim+j]=0;
                    }
                }
            }else if constexpr(std::is_same_v<Tl, T>){
                // Single layer, single element
                if(l.size() == 1){
                  for(size_t i = 0; i<dim*dim; i++){val[i]=*(l.begin());}
                }else{
                  // Single layer - as many as given in 1-D order
                  const size_t ct=std::min((int)l.size(), dim*dim);
                  for(size_t i = 0; i<dim*dim; i++){val[i]= *(l.begin()+i);}
                  for(size_t i = ct; i<dim*dim; i++){val[i]=0;}
                }
            }
        }

        /// Initialisation from a list of STScalars
        /** \sa STTensor(std::initializer_list<Tl> l) for details of how elements are mapped */
        template<typename Tl, typename std::enable_if_t<std::is_same_v<Tl, STScalar<T> > || std::is_same_v<Tl, std::initializer_list<STScalar<T> > >, int > =0 >
        constexpr STTensor(std::initializer_list<Tl> l){
            const size_t ct=std::min((int)l.size(), dim);
            if constexpr (std::is_same_v<Tl, std::initializer_list<STScalar<T>>>){
                for(size_t i = 0; i<ct; i++){
                    auto l2 = *(l.begin()+i);
                    const size_t ct2=std::min((int)l2.size(), dim);
                    for(size_t j=0; j<ct2; j++){
                      // Assign as many elements as are given
                      val[i*dim+j]= *(l2.begin()+j)->get();
                    }
                    for(size_t j = ct2; j<dim; j++){
                        val[i*dim+j]=0;
                    }
                }
                for(size_t i = ct; i<dim; i++){
                    for(size_t j = 0; j<dim; j++){
                        val[i*dim+j]=0;
                    }
                }
            }else{
                // Single layer, single element
                if(l.size() == 1){
                  for(size_t i = 0; i<dim*dim; i++){val[i]=(l.begin())->get();}
                }else{
                  // Single layer - as many as given in 1-D order
                  const size_t ct=std::min((int)l.size(), dim*dim);
                  for(size_t i = 0; i<ct; i++){val[i]= *(l.begin()+i)->get();}
                  for(size_t i = ct; i<dim*dim; i++){val[i]=0;}
                }
            }
        }

        /// Initialisation from a list of STVectors
        /** In this case we allow a list of vectors of the correct size (dim), and set as many rows as are given, zeroing any missing ones.
        */
        constexpr STTensor(std::initializer_list<STVector<T, dim> > l){
            const size_t ct=std::min((int)l.size(), dim);
            for(size_t i = 0; i<ct; i++){
              for(size_t j = 0; j<dim; j++){
                val[i*dim+j]=(l.begin()+i)->get(j);
              }
            }
            for(size_t i = ct; i<dim; i++){
              for(size_t j = 0; j<dim; j++){
                val[i*dim+j] = 0;
              }
            }
        }

        /// Typedef to get storage type from arbitrary UnitChecked type
        template <typename U>
        using WrappedType = decltype(std::declval<U>().stripUnits());
        /// Initializer list constructor from doubly nested list of wrapped STScalars
        /** This is needed to allow constructing a Tensor of UnitChecked types from a list of Scalar UnitChecked types. Unit checking is done BEFORE this point is reached
         * \sa STTensor(std::initializer_list< STScalar<T> > l) for behaviour of various lengths of list
         * \param l The list of lists of wrapped STScalars
         * \tparam U Type of list elements, which must be a Wrapped STScalar<T> for same T as this Tensor. That is, when we apply stripUnits to U we must get STScalar<T>
        */
        template <typename U, typename std::enable_if_t<std::is_same_v<WrappedType<U>, STScalar<T> >, int > =0 >
        constexpr STTensor(std::initializer_list<std::initializer_list<U> > l){
            const size_t ct=std::min((int)l.size(), dim);
            for(size_t i = 0; i<ct; i++){
              const size_t ct2=std::min((int)(l.begin()+i)->size(), dim);
              for(size_t j = 0; j<ct2; j++){
                val[i*dim+j]=(l.begin()+i)->begin()[j].stripUnits().get();
              }
              for(size_t j = ct2; j<dim; j++){
                val[i*dim+j]=0;
              }
            }
            for(size_t i = ct; i<dim; i++){
              for(size_t j = 0; j<dim; j++){
                val[i*dim+j] = 0;
              }
            }
        }

        /// Initializer list constructor from a list of wrapped STVectors
        /** This is needed to allow constructing a Tensor of UnitChecked types from a list of Vector UnitChecked types. Unit checking is done BEFORE this point is reached
         * \sa STTensor(std::initializer_list< STVector<T, dim> > l) for behaviour of various lengths of list
         * \param l The list of wrapped STVectors
         * \tparam U Type of list elements, which must be a Wrapped STVector<T, dim> for same T and dim as this Tensor. That is, when we apply stripUnits to U we must get STVector<T, dim>
        */
        template <typename U, typename std::enable_if_t<std::is_same_v<WrappedType<U>, STVector<T, dim> >, int > =0 >
        constexpr STTensor(std::initializer_list<U> l){
            const size_t ct=std::min((int)l.size(), dim);
            for(size_t i = 0; i<ct; i++){
              auto tmp = (l.begin()+i)->stripUnits();
              for(size_t j = 0; j<dim; j++){
                val[i*dim+j]=tmp.get(j);
              }
            }
            for(size_t i = ct; i<dim; i++){
              for(size_t j = 0; j<dim; j++){
                val[i*dim+j] = 0;
              }
            }
        }
        /// 1-D access (0<=i <dim*dim) by reference
        constexpr T& operator[](size_t i){
            return val[i];
        }
        /// 1-D access (0<=i <dim*dim) by copy
        constexpr T operator[](size_t i)const{
            return val[i];
        }
        /// 2-D access (0<=i,j<dim) by reference
        constexpr T& get(size_t i, size_t j){
            return val[i*dim+j];
        }
        /// 2-D access (0<=i,j<dim) by copy
        constexpr T get(size_t i, size_t j)const{
            return val[i*dim+j];
        }

        /// Get element (as a valid StorageType) - both indices -> scalar
        constexpr STScalar<T> getElement(size_t i, size_t j)const{
            return STScalar<T>{val[i*dim+j]};
        }
        /// Get element (as a valid StorageType) - single index -> vector \todo Is this the right way round (column vs row)
        constexpr STVector<T, dim> getElement(size_t i)const{
            STVector<T, dim> out;
            for(size_t j=0; j<dim; j++){
              out[j]=val[i*dim+j];
            }
            return out;
        }
        ///Get reference to an element
        constexpr STScalarRef<T> getElementRef(size_t i, size_t j){
            return STScalarRef<T>(val[i*dim+j]);
        }

        /// Identity entity (i.e. diagonal ones)
        static STTensor<T,dim> identity(){
          return STTensor<T, dim>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        }

        /// Linear algebra transpose (i<-> j)
        constexpr STTensor transpose()const{
          STTensor out;
          for(size_t i = 0; i<dim; i++){
            for(size_t j = 0; j<dim; j++){
              out[i*dim+j]=val[j*dim+i];
            }
          }
          return out;
        }

        //// Equality element wise. No ordering for tensors
        friend bool operator==(const STTensor & first, const STTensor & other){
          for(size_t i = 0; i<dim*dim; i++){
            if(first.val[i]!=other.val[i]) return false;
          }
          return true;
        }
        /// Inequality element wise
        friend bool operator!=(const STTensor & first, const STTensor & other){
          return ! (first == other);
        }

        /// Unary minus
        STTensor operator-()const{
          STTensor out;
          for(size_t i = 0; i<dim*dim; i++){
            out[i]=-val[i];
          }
          return out;
        }
        STTensor operator +=(const STTensor & other){
          for(size_t i = 0; i<dim*dim; i++){
            val[i]+=other[i];
          }
          return *this;
        }
        friend STTensor operator+(STTensor lhs, const STTensor & other){
          return lhs+=other;
        }
        STTensor operator -=(const STTensor & other){
          for(size_t i = 0; i<dim*dim; i++){
            val[i]-=other[i];
          }
          return *this;
        }
        friend STTensor operator-(STTensor lhs, const STTensor & other){
          return lhs-=other;
        }

        STTensor operator *=(const STTensor & other){
            for(size_t i = 0; i<dim*dim; i++){
                val[i] *= other.val[i];
            }
            return *this;
        }
        friend STTensor operator*(STTensor lhs, const STTensor & other){
            return lhs*=other;
        }
        STTensor operator/=(const STTensor & other){
            for(size_t i = 0; i<dim*dim; i++){
                val[i] /= other.val[i];
            }
            return *this;
        }
        friend STTensor operator/(STTensor lhs, const STTensor & other){
            return lhs/=other;
        }
};
/// Stream output for Tensors
/** Outputs a bracketed, comma separated, and new-line broken array, e.g.
 * (1.0, 2.0)
 * (3.0, 4.0)
*/
template <typename T, int dim>
std::ostream& operator<<(std::ostream& os, const STTensor<T, dim>& val_in){
  for(size_t i = 0; i<dim; i++){
    os << "(";
    for(size_t j = 0; j<dim; j++){
      os << val_in[i*dim + j];
      if(j<dim-1) os << ',';
    }
    if(i < dim-1) os << ")\n";
  }
  os<<")";
  return os;
};

///Stream input for Tensors
/** Expects the format output by operator<< so e.g.
 * (1.0, 2.0)
 * (3.0, 4.0)
*/
template <typename T, int dim>
std::istream& operator>>(std::istream& is, STTensor<T, dim>& val_in){
  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      if(!is.good()) return is;
      char seps;
      is>>seps; // Read and toss the bracket and commas
      if(seps != '(' && seps != ',') return is;
      is >>val_in.get(i, j);
    }
    char line_seps;
    is>>line_seps; // Read and toss the ')\n' part
  }
  return is;
};



///Element-wise multiplication of a Scalar and a Vector
template <typename T, int dim>
STVector<T, dim> operator*(const STScalar<T> &a, const STVector<T, dim> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a.get()*b.get(i);
  }
  return out;
}
///Element-wise multiplication of a Vector and a Scalar
template <typename T, int dim>
STVector<T, dim> operator*(const STVector<T, dim> &a, const STScalar<T> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a.get(i)*b.get();
  }
  return out;
}
///Element-wise division of a Scalar and a Vector
/** Implemented mostly for completeness, assumes we want the vector result where each element is the scalar divided by the corresponding element of the vector
*/
template <typename T, int dim>
STVector<T, dim> operator/(const STScalar<T> &a, const STVector<T, dim> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a.get()/b.get(i);
  }
  return out;
}
///Element-wise division of a Vector and a Scalar
template <typename T, int dim>
STVector<T, dim> operator/(const STVector<T, dim> &a, const STScalar<T> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a.get(i)/b.get();
  }
  return out;
}

///Element-wise multiplication of a Scalar and a Tensor
template <typename T, int dim>
STTensor<T, dim> operator*(const STScalar<T> &a, const STTensor<T, dim> &b){
  STTensor<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      out.get(i,j) = a.get()*b.get(i,j);
    }
  }
  return out;
}
///Element-wise multiplication of a Tensor and a Scalar
template <typename T, int dim>
STTensor<T, dim> operator*(const STTensor<T, dim> &a, const STScalar<T> &b){
  STTensor<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      out.get(i,j) = a.get(i,j)*b.get();
    }
  }
  return out;
}

///Vector-Tensor multiply, Vector on left
/** This is the operation where we treat the vector as a row vector and multiply it by the tensor on the right, i.e. vT. Note that we don't check/care if vector _is_ a row vector, we just treat it as such
 * \todo Consider adding co-contra variance to the vector type to make this more explicit. If so, implement transpose too
*/
template <typename T, int dim>
STTensor<T, dim> operator*(const STVector<T, dim> &a, const STTensor<T, dim> &b){
  STTensor<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      out[j*dim + i] += a[j]*b[j*dim+i]; // Is this the right way round???
    }
  }
  return out;
}
///Vector-Tensor multiply, Vector on right
/** This is the operation where we treat the vector as a column vector and multiply it by the tensor on the left, i.e. Tv. Note that we don't check/care if vector _is_ a column vector, we just treat it as such
*/
template <typename T, int dim>
STTensor<T, dim> operator*(const STTensor<T, dim> &a, const STVector<T, dim> &b){
  STTensor<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      out[j*dim + i] += a[j*dim+i]*b[j]; // Is this the right way round???
    }
  }
  return out;
}

// Special functions
// Provide any function taking dimensionless argument(s) like this
// EXAMPLES ONLY
///Example of a Sin function on Scalars
template<typename T>
T mysinfunction(const STScalar<T> &a){
  return std::sin(a.get());
}
///Example of a Sin function on Vectors
template<typename T>
STVector<T, 3> mysinfunction(const STVector<T, 3> &a){
  STVector<T, 3> out;
  for(size_t i = 0; i<3; i++){
    out[i] = std::sin(a.get(i));
  }
  return out;
}



#endif
