#ifndef HELPER_H
#define HELPER_H

// Template magic from https://stackoverflow.com/questions/11056714/c-type-traits-to-extract-template-parameter-class to extract T from X<T>
// Template sequence to extract the underlying value and integer parameter from a templated class. This is overloaded for the patterns in StorageTypes

/// Generic template for any class T: will match any T and extract T as value_type and dim of 0
template<typename T>
struct extract_value_type
{
    typedef T value_type;
    static const int dim = 0;
};
/// Overload for classes of the form X<T> - will extract T as value_type and dim of 0
template<template<typename> class X, typename T>
struct extract_value_type<X<T>>
{
    typedef T value_type;
    static const int dim = 0;
};
/// Overload for classes of the form X<T, int> - will extract T as value_type and second param as dim
template<template<typename,int> class X, typename T, int dim_>
struct extract_value_type<X<T, dim_>>
{
    typedef T value_type;
    static const int dim = dim_;
};

#endif