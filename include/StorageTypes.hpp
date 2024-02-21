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
    double val;
  public:
    STDummy()=default;
    STDummy(double a){val=a;}
    STDummy(std::initializer_list<double> l){l.size()>0? val=*(l.begin()):val=0;};

};
#endif

template <typename T>
class STScalar{
    private:
        T val{};

    public:
        STScalar(){};
        STScalar(T val_in):val(val_in){};
        STScalar(std::initializer_list<T> l){l.size()>0? val=*(l.begin()):val=0;};

        STScalar(const STScalar &a) = default;
        STScalar operator=(const STScalar &a){
          val=a.val;
          return *this;
        }
        T& operator[](size_t i){
            return val;
        }// Include these for consistency
        T operator[](size_t i)const{
            return val;
        }
        T& get(){
            return val;
        }
        T get()const{
            return val;
        }

        explicit operator T() const{return val;}

        template<typename num, typename=std::enable_if_t<std::is_arithmetic_v<num> > >
        explicit operator num() const{
#ifdef NO_NARROWING_CONVERSIONS
          return num{val};
# else
          return static_cast<num>(val);
 #endif
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
        STScalar pow(long p)const{
          return std::pow(val, p);
        }
        STScalar pow(double p)const{
          return std::pow(val, p);
        }

        STScalar dot(const STScalar & other)const{
          return val*other.val;
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
        STVector(){};
        STVector(T val_in){for(size_t i = 0; i<dim; i++){val[i]=val_in;}};
        STVector(std::initializer_list<T> l){
            memset(val, 0, sizeof(T)*dim);
            const size_t ct=std::min((int)l.size(), dim);
            for(size_t i = 0; i<ct; i++){val[i]=*(l.begin()+i);}
        }
        //Allow for initialisation from any other type which offers a get() method, such as an STScalar, but potentially wrapped in Units etc
        template <typename U>
        STVector(std::initializer_list<U> l){
            memset(val, 0, sizeof(T)*dim);
            const size_t ct=std::min((int)l.size(), dim);
            for(size_t i = 0; i<ct; i++){val[i]=(l.begin()+i)->unsafeGet();}
        }
        STVector(const STVector &a) = default;
        STVector operator=(const STVector &a){
          for(size_t i = 0; i<dim; i++){
            val[i]=a[i];
          }
          return *this;
        }

        T& operator[](size_t i){
            return val[i];
        }
        T operator[](size_t i)const{
            return val[i];
        }
        T& get(int i){
            return val[i];
        }
        T get(int i)const{
            return val[i];
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
          for(size_t i = 0; i<dim; i++){
            if(first.val[i]!=other.val[i]) return true;
          }
          return false;
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
class STTensor{
    private:
        T val[dim*dim]{};

    public:
        STTensor(){};
        STTensor(T val_in){for(size_t i = 0; i<dim*dim; i++){val[i]=val_in;}};
        STTensor(const STTensor &a) = default;
        STTensor operator=(const STTensor &a){
          for(size_t i = 0; i<dim*dim; i++){
            val[i]=a[i];
          }
          return *this;
        }

        // Allow initialisation from single element, or from doubly nested list
        template<typename Tl>
        STTensor(std::initializer_list<Tl> l){
            // Force everything to zero
            memset(val, 0, sizeof(T)*dim*dim);
            const size_t ct=std::min((int)l.size(), dim);
            if constexpr (std::is_same_v<Tl, std::initializer_list<T>>){
                for(size_t i = 0; i<ct; i++){
                    auto l2 = *(l.begin()+i);
                    const size_t ct2=std::min((int)l2.size(), dim);
                    for(size_t j=0; j<ct2; j++){
                      // Assign as many elements as are given
                      val[i*dim+j]= *(l2.begin()+j);
                    }
                }
            }else if constexpr(std::is_same_v<Tl, T>){
                // Single layer - assume single element
                for(size_t i = 0; i<dim*dim; i++){val[i]=*(l.begin());}
            }
        }
        T& operator[](size_t i){
            return val[i];
        }
        T operator[](size_t i)const{
            return val[i];
        }
        T& get(size_t i, size_t j){
            return val[i*dim+j];
        }
        T get(size_t i, size_t j)const{
            return val[i*dim+j];
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
  os << "(";
  for(size_t i = 0; i<dim*dim; i++){
    os << val_in[i];
    if(i<dim-1) os << ',';
  }
  os<<")";
  return os;
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

#endif
