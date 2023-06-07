#ifndef __STORAGE_TYPES_HPP__
#define __STORAGE_TYPES_HPP__

#include <algorithm>
#include <cstring>

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

 

template <typename T>
class STScalar{
    private:
        T val;

    public:
        STScalar(T val_in):val(val_in){};
        STScalar(std::initializer_list<T> l){l.size()>0? val=*(l.begin()):val=0;};

        STScalar(const STScalar &a) = default;
        STScalar operator=(const STScalar &a){
          val=a.val;
          return *this;
        }
        T& operator[](size_t i){
            return val;
        }// Dumb but required
        T operator[](size_t i)const{
            return val;
        }
        T& get(){
            return val;
        }
        T get()const{
            return val;
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
        T val[dim];

    public:
        STVector(){};
        STVector(T val_in){for(size_t i = 0; i<dim; i++){val[i]=val_in;}};
        STVector(std::initializer_list<T> l){
            memset(val, 0, sizeof(T)*dim);
            const size_t ct=std::min((int)l.size(), dim);
            for(size_t i = 0; i<ct; i++){val[i]=*(l.begin()+i);}
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
        
        STVector operator +=(const STVector & other){
          for(size_t i = 0; i<dim; i++){
            val[i]+=other[i];
          }
          return *this;
        }
        friend STVector operator+(STVector lhs, const STVector & other){
          return lhs+=other;
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
        T val[dim*dim];

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

// Scalar value multiply
template <typename T>
STScalar<T> operator*(const STScalar<T> &a, const STScalar<T> &b){
  return STScalar<T>{a[0]*b[0]};
}

// Element wise vector multiply - same length only
template <typename T, int dim>
STVector<T, dim> operator*(const STVector<T, dim> &a, const STVector<T, dim> &b){
  STVector<T, dim> out;
  for(size_t i = 0; i<dim; i++){
    out[i] = a[i]*b[i];
  }
  return out;
}

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


#endif
