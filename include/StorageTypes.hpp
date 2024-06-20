#ifndef __STORAGE_TYPES_HPP__
#define __STORAGE_TYPES_HPP__

#include <algorithm>
#include <cstring>
#include <cmath>

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
// Implements nothing - for development purposes
class STDummy{
  public:
    STDummy()=default;
};
#endif

template <typename T>
class STScalar{
    private:
        T val{};

    public:
        constexpr STScalar(){};
        constexpr STScalar(T val_in):val(val_in){};
        constexpr STScalar(std::initializer_list<T> l){l.size()>0? val=*(l.begin()):val=0;};

        constexpr STScalar(const STScalar &a) = default;
        STScalar operator=(const STScalar &a){
          val=a.val;
          return *this;
        }
        constexpr T& operator[](size_t i){
            return val;
        }// Include these for consistency
        constexpr T operator[](size_t i)const{
            return val;
        }
        constexpr T& get(){
            return val;
        }
        constexpr T get()const{
            return val;
        }

        constexpr explicit operator T() const{return val;}

        template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
        explicit operator num() const{
#ifdef NO_NARROWING_CONVERSIONS
          return num{val};
# else
          return static_cast<num>(val);
 #endif
        }

        static STScalar<T> identity(){
          return STScalar<T>{1};
        }

        STScalar<T> magnitude()const{
          return std::abs(val);
        }

        STScalar operator-()const{
          return -val;
        }
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
        }
        STScalar cbrt()const{
          return std::cbrt(val);
        }
        STScalar pow(double p)const{
          return std::pow(val, p);
        }

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
template <typename T>
std::ostream& operator<<(std::ostream& os, const STScalar<T>& val_in){
  os <<val_in[0];
  return os;
};
template <typename T>
std::istream& operator>>(std::istream& is, STScalar<T>& val_in){
  is >>val_in.get();
  return is;
};

//Forward declare for use in outer product
template <typename T, int dim>
class STTensor;

template <typename T, int dim>
class STVector{
    private:
        T val[dim]{};

        T normSq()const{
          T sum=0;
          for(size_t i = 0; i<dim; i++){
            sum+=val[i]*val[i];
          }
          return sum;
        }
    public:
        constexpr STVector(){};
        constexpr STVector(T val_in){for(size_t i = 0; i<dim; i++){val[i]=val_in;}};
        constexpr STVector(const STVector &a) = default;
        constexpr STVector(std::initializer_list<T> l){
            const size_t ct=std::min((int)l.size(), dim);
            if(ct == 1){
              for(size_t i = 0; i < dim; i++){val[i]=*(l.begin());}
            }else{
              for(size_t i = 0; i<ct; i++){val[i]=*(l.begin()+i);}
              for(size_t i = ct; i<dim; i++){val[i]=0;} // Zero any excess values
            }
        }
        // From STScalars
        constexpr STVector(std::initializer_list< STScalar<T> > l){
            const size_t ct=std::min((int)l.size(), dim);
            if(ct == 1){
              for(size_t i = 0; i < dim; i++){val[i]=(l.begin())->get();}
            }else{
              for(size_t i = 0; i<ct; i++){val[i]=(l.begin()+i)->get();}
              for(size_t i = ct; i<dim; i++){val[i]=0;} // Zero any excess values
            }
        }

        // Allow initialisation from wrapped STScalars, where wrapper is assumed to have a stripUnits method to get this back
        template <typename U>
        using WrappedType = decltype(std::declval<U>().stripUnits());
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

        STVector operator=(const STVector &a){
          for(size_t i = 0; i<dim; i++){
            val[i]=a[i];
          }
          return *this;
        }

        constexpr T& operator[](size_t i){
            return val[i];
        }
        constexpr T operator[](size_t i)const{
            return val[i];
        }
        constexpr T& get(int i){
            return val[i];
        }
        constexpr T get(int i)const{
            return val[i];
        }

        static STVector<T,dim> identity(){
          return STVector<T, dim>{1, 1, 1};
        }

        STScalar<T> magnitude()const{
          return STScalar<T>{std::sqrt(normSq())};
        }

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

        STScalar<T> dot(const STVector & other)const{
          STScalar<T> sum=0;
          for(size_t i = 0; i<dim; i++){
            sum+=val[i]*other.val[i];
          }
          return sum;
        }

        STVector<T,dim> cross(const STVector & other) const{
          static_assert(dim ==3); // Well defined only for dim 3 vector
          STVector<T,dim> out;
          for(size_t i = 0; i<dim; i++){
            out[i] = val[(i+1)%dim]*other.val[(i+2)%dim] - val[(i+2)%dim]*other.val[(i+1)%dim];
          }
          return out;
        }

        STTensor<T, dim> outer(const STVector & other) const{
          STTensor<T, dim> out;
          for(size_t i = 0; i<dim; i++){
            for(size_t j = 0; j<dim; j++){
              out[i*dim+j] = val[i]*other.val[j];
            }
          }
          return out;
        }

        void normalize(){
          const STScalar<T> mag = magnitude();
          for(size_t i = 0; i<dim; i++){
            val[i]/=mag[0];
          }
        }

        // Comparisons - in terms of ordering of the norm only
        friend bool operator==(const STVector & first, const STVector & other){
          for(size_t i = 0; i<dim; i++){
            if(first.val[i]!=other.val[i]) return false;
          }
          return true;
        }
        friend bool operator!=(const STVector & first, const STVector & other){
          return ! (first == other);
        }
        friend bool operator<(const STVector & first, const STVector & other){
          return first.normSq()<other.normSq();
        }

        // Compare with Scalars
        friend bool operator==(const STVector & a, const STScalar<T> & b){
          return a.normSq()==b.magnitude()*b.magnitude();
        }
        friend bool operator==(const STScalar<T> & a, const STVector & b){
          return b==a;
        }
        friend bool operator<(const STVector & a, const STScalar<T> & b){
          return a.normSq()<b.magnitude()*b.magnitude();
        }
        friend bool operator<(const STScalar<T> & a, const STVector & b){
          return a.magnitude()*a.magnitude() < b.normSq();
        }

        // Implement the rest in terms of < for easy tweaking/ expansion
        template<typename T2>
        friend bool operator>(const STVector & first, const T2 & other){
          return other < first;
        }
        template<typename T2>
        friend bool operator<=(const STVector & first, const T2 & other){
          return !(first > other);
        }
        template<typename T2>
        friend bool operator>=(const STVector & first, const T2 & other){
          return !(first < other);
        }

};

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

template <typename T, int dim>
class STTensor{
    private:
        T val[dim*dim]{};

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

        // Allow initialisation from single element, or from doubly nested list, only of type T
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

        // From STScalars, ditto
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

        // And from STVectors, single nest only
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

        // Allow initialisation from (doubly nested list of) wrapped STScalars, where wrapper is assumed to have a stripUnits method to get this back
        template <typename U>
        using WrappedType = decltype(std::declval<U>().stripUnits());
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

        // Ditto but from STVectors (single nested list)
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

        constexpr T& operator[](size_t i){
            return val[i];
        }
        constexpr T operator[](size_t i)const{
            return val[i];
        }
        constexpr T& get(size_t i, size_t j){
            return val[i*dim+j];
        }
        constexpr T get(size_t i, size_t j)const{
            return val[i*dim+j];
        }

        static STTensor<T,dim> identity(){
          return STTensor<T, dim>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        }

        constexpr STTensor transpose()const{
          STTensor out;
          for(size_t i = 0; i<dim; i++){
            for(size_t j = 0; j<dim; j++){
              out[i*dim+j]=val[j*dim+i];
            }
          }
          return out;
        }

        // Comparision, but no ordering
        friend bool operator==(const STTensor & first, const STTensor & other){
          for(size_t i = 0; i<dim*dim; i++){
            if(first.val[i]!=other.val[i]) return false;
          }
          return true;
        }
        friend bool operator!=(const STTensor & first, const STTensor & other){
          return ! (first == other);
        }

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



// Scalar-vector multiply, both ways round
template <typename T, int dim>
STVector<T, dim> operator*(const STScalar<T> &a, const STVector<T, dim> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a[0]*b[i];
  }
  return out;
}
template <typename T, int dim>
STVector<T, dim> operator*(const STVector<T, dim> &a, const STScalar<T> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a[i]*b[0];
  }
  return out;
}
// Scalar-vector divide, both ways round
template <typename T, int dim>
STVector<T, dim> operator/(const STScalar<T> &a, const STVector<T, dim> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a[0]/b[i];
  }
  return out;
}
template <typename T, int dim>
STVector<T, dim> operator/(const STVector<T, dim> &a, const STScalar<T> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a[i]/b[0];
  }
  return out;
}

// Scalar-Tensor multiply, both ways round
template <typename T, int dim>
STTensor<T, dim> operator*(const STScalar<T> &a, const STTensor<T, dim> &b){
  STTensor<T, dim> out;
  for(size_t i = 0; i<dim*dim; i++){
    out[i] = a[0]*b[i];
  }
  return out;
}
template <typename T, int dim>
STTensor<T, dim> operator*(const STTensor<T, dim> &a, const STScalar<T> &b){
  STTensor<T, dim> out;
  for(size_t i = 0; i<dim*dim; i++){
    out[i] = a[i]*b[0];
  }
  return out;
}

// Vector-Tensor multiply, both ways round
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
template<typename T>
T mysinfunction(const STScalar<T> &a){
  return std::sin(a[0]);
}
template<typename T>
STVector<T, 3> mysinfunction(const STVector<T, 3> &a){
  STVector<T, 3> out;
  for(size_t i = 0; i<3; i++){
    out[i] = std::sin(a[i]);
  }
  return out;
}



#endif
