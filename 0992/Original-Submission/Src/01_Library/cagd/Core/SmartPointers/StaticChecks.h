//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/StaticChecks.h
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

#ifndef STATICCHECKS_H
#define STATICCHECKS_H

namespace cagd
{
    //(*@\Green{// The non-specialized variant of the class CompileTimeError (i.e., which belongs to the compile time constant}@*)
    //(*@\Green{// false) has only a declaration. Compared to this, the specialization that corresponds to the compile time}@*)
    //(*@\Green{// constant true also provides a definition. Therefore, in case of false values any instantiation will generate a}@*)
    //(*@\Green{// compile time error.}@*)
    template <bool> class CompileTimeError;
    template <>     class CompileTimeError<true>{};

    //(*@\Green{// The macro STATIC\_CHECK will generate a compile time error whenever the given expression evaluates to false.}@*)
    //(*@\Green{// In such cases the given message will appear in the list of compile time errors.}@*)
    #define STATIC_CHECK(expression, message){\
        CompileTimeError<((expression) != false)> ERROR_##message;}
}

#endif // STATICCHECKS_H
