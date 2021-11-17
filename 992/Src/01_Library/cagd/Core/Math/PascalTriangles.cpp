//----------------------------------------------------------------------------------
// File:        Core/Math/PascalTriangles.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "PascalTriangles.h"

namespace cagd
{
    //(*@\Green{// default/special constructor}@*)
    PascalTriangle::PascalTriangle(int order): TriangularMatrix<double>(order + 1)
    {
        _data[0][0] = 1.0;

        for (int k = 1; k <= order; k++)
        {
            _data[k][0] = _data[k][k] = 1.0;

            #pragma omp parallel for
            for (int l = 1; l <= k / 2; l++)
            {
                _data[k][l]     = _data[k - 1][l - 1] + _data[k - 1][l];
                _data[k][k - l] = _data[k][l];
            }
        }
    }

    //(*@\Green{// redefined virtual resizing method}@*)
    bool PascalTriangle::resizeRows(int row_count)
    {
        int old_row_count = _row_count;

        if (!TriangularMatrix::resizeRows(row_count))
        {
            return false;
        }

        for (int k = old_row_count; k < row_count; k++)
        {
            _data[k][0] = _data[k][k] = 1.0;

            #pragma omp parallel for
            for (int l = 1; l <= k / 2; l++)
            {
                _data[k][l]     = _data[k - 1][l - 1] + _data[k - 1][l];
                _data[k][k - l] = _data[k][l];
            }
        }

        return true;
    }
}
