#ifndef HELPER_H
#define HELPER_H

// Template magic from https://stackoverflow.com/questions/11056714/c-type-traits-to-extract-template-parameter-class to extract T from X<T>
template<typename T>
struct extract_value_type //lets call it extract_value_type
{
    typedef T value_type;
    static const int dim = 0;
};
// Overload for classes with just a type
template<template<typename> class X, typename T>
struct extract_value_type<X<T>>   //specialization
{
    typedef T value_type;
    static const int dim = 0;
};
// Overload for classes with a type and an int (rank)
template<template<typename,int> class X, typename T, int dim_>
struct extract_value_type<X<T, dim_>>   //specialization
{
    typedef T value_type;
    static const int dim = dim_;
};

#endif