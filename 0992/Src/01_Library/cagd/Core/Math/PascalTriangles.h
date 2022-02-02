//----------------------------------------------------------------------------------
// File:        Core/Math/PascalTriangles.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef PASCALTRIANGLES_H
#define PASCALTRIANGLES_H

#include "Matrices.h"

namespace cagd
{
    class PascalTriangle: public TriangularMatrix<double>
    {
    public:
        //(*@\Green{// default/special constructor}@*)
        PascalTriangle(int order = 0);

        //(*@\Green{// redeclared virtual resizing method}@*)
        bool resizeRows(int row_count);
    };
}

#endif //(*@\Green{// PASCALTRIANGLES\_H}@*)
