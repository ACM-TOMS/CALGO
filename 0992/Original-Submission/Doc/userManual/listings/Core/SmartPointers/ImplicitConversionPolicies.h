//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/ImplicitConversionPolicies.h
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

#ifndef IMPLICITCONVERSIONPOLICIES_H
#define IMPLICITCONVERSIONPOLICIES_H

namespace cagd
{
    (*@\Green{// Using implicit conversion policies, one can allow or disallow to convert smart pointers to raw pointers of}@*)
    (*@\Green{// stored type.}@*)
    class ImplicitConversionPolicy
    {
    public:
        class Allowed
        {
        public:
            enum {ENABLED = true};
        };

        class Disallowed
        {
        public:
            enum {ENABLED = false};
        };
    };
}

#endif (*@\Green{// IMPLICITCONVERSIONPOLICIES\_H}@*)
