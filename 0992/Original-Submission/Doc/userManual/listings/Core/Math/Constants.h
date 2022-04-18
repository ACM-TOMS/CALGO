//----------------------------------------------------------------------------------
// File:        Core/Math/Constants.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "PascalTriangles.h"
#include <limits>

namespace cagd
{
    const double         PI             = 3.1415926535897932384626433832795;
    const double         TWO_PI         = 2.0 * PI;
    const double         HALF_PI        = PI / 2.0;
    const double         DEG_TO_RADIAN  = PI / 180.0;
    const double         EPS            = 1.0e-9;
    const double         MACHINE_EPS    = std::numeric_limits<double>::epsilon();
    const double         TINY           = std::numeric_limits<double>::min();

    const int            MAX_DIFF_ORDER = 500; (*@\Green{// maximal differentiation order}@*)
    const PascalTriangle BC(MAX_DIFF_ORDER);   (*@\Green{// Pascal triangle of binomial coefficients}@*)

    namespace variable
    {
        enum Type{U = 0, V = 1};
    }
}

#endif (*@\Green{// CONSTANTS\_H}@*)
