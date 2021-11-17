//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/TypeSelectors.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
// Reference:   A. Alexandrescu. 2001. Modern C++ Design: Generic Programming and Design
//              Patterns Applied, 1st edition. Addison-Wesley Professional, USA.
//--------------------------------------------------------------------------------------

#ifndef TYPESELECTORS_H
#define TYPESELECTORS_H

namespace cagd
{
    (*@\Green{// based on the logical value of the compile time constant select\_first\_type, the class selects one of the}@*)
    (*@\Green{// types U and V}@*)
    template <bool select_first_type, typename U, typename V>
    class TypeSelector
    {
    public:
        typedef U Result;
    };

    (*@\Green{// specializations for all possible cases when the value of select\_first\_type is false}@*)
    template <typename U, typename V>
    class TypeSelector<false, U, V>
    {
    public:
        typedef V Result;
    };

    (*@\Green{// logical template function that decides whether the given typename T denotes a primitive type or not}@*)
    template <typename T>
    inline bool primitive(T)
    {
        return false; (*@\Green{// with the exception of the specializations listed below, the function always returns false}@*)
    }

    (*@\Green{// specializations that belong to primitive types always return true}@*)
    template <>
    inline bool primitive(int)
    {
        return true;
    }

    template <>
    inline bool primitive(unsigned int)
    {
        return true;
    }

    template <>
    inline bool primitive(short int)
    {
        return true;
    }

    template <>
    inline bool primitive(unsigned short int)
    {
        return true;
    }

    template <>
    inline bool primitive(long int)
    {
        return true;
    }

    template <>
    inline bool primitive(unsigned long int)
    {
        return true;
    }

    template <>
    inline bool primitive(long long int)
    {
        return true;
    }

    template <>
    inline bool primitive(unsigned long long int)
    {
        return true;
    }

    template <>
    inline bool primitive(signed char)
    {
        return true;
    }

    template <>
    inline bool primitive(char)
    {
        return true;
    }

    template <>
    inline bool primitive(unsigned char)
    {
        return true;
    }

    template <>
    inline bool primitive(float)
    {
        return true;
    }

    template <>
    inline bool primitive(double)
    {
        return true;
    }

    template <>
    inline bool primitive(long double)
    {
        return true;
    }

    template <>
    inline bool primitive(wchar_t)
    {
        return true;
    }

    (*@\Green{// if we have omitted some primitive types, then this file should be appended by following the given examples above\ldots}@*)
}

#endif (*@\Green{// TYPESELECTORS\_H}@*)
