#ifndef __STORAGE_TYPES_HPP__
#define __STORAGE_TYPES_HPP__

// To be valid storage types for the UnitChecked class, these must implement the following:
//

template <typename T>
class STScalar{
    private:
        T val;

    public:
        STScalar(T val_in):val(val_in){};
        STScalar(std::initializer_list<T> l):val(l){};

        int size()const{
            return 1;
        }
        T& operator[](size_t i){
            return val;
        }// Dumb but required
        T operator[](size_t i)const{
            return val;
        }
};
template <typename T>
std::ostream& operator<<(std::ostream& os, const STScalar<T>& val_in){
  os <<val_in[0];
  return os;
};

template <typename T>
class ST3Vector{
    private:
        T val[3];

    public:
        ST3Vector(){};
        ST3Vector(T val_in):val{val_in,val_in,val_in}{};
        ST3Vector(std::initializer_list<T> l){l.size() ==3? (val[0]=*(l.begin()),val[1]=*(l.begin()+1),val[2]=*(l.begin()+2)) : (val[0]=0,val[1]=0,val[2]=0);};

        int size()const{
            return 3;
        }

        T& operator[](size_t i){
            return val[i];
        }
        T operator[](size_t i)const{
            return val[i];
        }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const ST3Vector<T>& val_in){
  os << "("<<val_in[0]<<','<<val_in[1]<<','<<val_in[2]<<")";
  return os;
};



#endif