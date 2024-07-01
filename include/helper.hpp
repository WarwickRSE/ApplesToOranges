#ifndef HELPER_H
#define HELPER_H

// Template magic from https://stackoverflow.com/questions/11056714/c-type-traits-to-extract-template-parameter-class to extract T from X<T>
// Several variations made here

///Specific helpers for template extraction and modification of classes of the form we use later
namespace STUtils{

/// Default to extract values a class is templated on - will match any T and extract T as value_type and dim of 0
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

///Default for re-templating on another type
template<typename T, typename T2>
struct modify_template_type
{
    typedef T2 modified_type;
};
///Overload to re-template X<ST> to X<STm>
template<template<typename> class X, typename ST ,typename STm>
struct modify_template_type<X<ST>, STm>
{
    typedef X<STm> modified_type;
};
///Default for adding const to a type (no-op)
template<typename T>
struct add_const
{
    typedef T modified_type;
};
///Overload to add const via X<ST, cc> -> X<ST, true>
template<template<typename, bool> class X, typename ST, bool cc>
struct add_const<X<ST, cc>>
{
    typedef X<ST, true> modified_type; //Bool param always true
};
};

#endif