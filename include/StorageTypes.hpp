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
class ST3Vector{
    private:
        T val[3];

    public:
        ST3Vector(){};
        ST3Vector(T val_in):val{val_in,val_in,val_in}{};
        ST3Vector(T x, T y, T z):val{x,y,z}{};
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

/*class ST3VectorDouble{
    private:
      double val[3];

    public:
        constexpr ST3VectorDouble(double x, double y, double z){val[0]=x;val[1]=y;val[2]=z;};
        constexpr int size()const{
            return 3;
        }

        double& operator[](size_t i){
            return val[i];
        }
        double operator[](size_t i)const{
            return val[i];
        }
};
*/

#endif