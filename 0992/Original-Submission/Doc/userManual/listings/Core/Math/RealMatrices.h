//----------------------------------------------------------------------------------
// File:        Core/Math/RealMatrices.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef REALMATRICES_H
#define REALMATRICES_H

#include "Matrices.h"

namespace cagd
{
    (*@\Green{// General (rectangular) real matrices.}@*)
    class RealMatrix: public Matrix<double>
    {
    public:
        (*@\Green{// default/special constructor}@*)
        explicit RealMatrix(int row_count, int column_count);

        (*@\Green{// matrix addition}@*)
        const RealMatrix operator +(const RealMatrix &rhs) const;

        (*@\Green{// matrix subtraction}@*)
        const RealMatrix operator -(const RealMatrix &rhs) const;

        (*@\Green{// matrix multiplication}@*)
        const RealMatrix operator *(const RealMatrix &rhs) const;

        (*@\Green{// multiplication by scalar from right}@*)
        const RealMatrix operator *(const double &rhs) const;

        (*@\Green{// multiplication by scalar from left}@*)
        friend const RealMatrix operator *(const double &lhs, const RealMatrix &rhs);

        (*@\Green{// multiplication by real column vector from right}@*)
        const RealMatrix operator *(const ColumnMatrix<double> &rhs) const;

        (*@\Green{// multiplication by real row vector from left}@*)
        friend const RealMatrix operator *(const RowMatrix<double> &lhs,
                                           const RealMatrix &rhs);

        (*@\Green{// multiplication of a real column matrix by a real row matrix from right}@*)
        friend const RealMatrix operator *(const ColumnMatrix<double> &lhs,
                                           const RowMatrix<double> &rhs);

        (*@\Green{// multiplication of a real row matrix by a real matrix from right}@*)
        friend RowMatrix<double>& operator *=(RowMatrix<double> &lhs,
                                              const RealMatrix &rhs);

        (*@\Green{// division by scalar from right}@*)
        const RealMatrix operator /(const double &rhs) const;

        (*@\Green{// add to *this}@*)
        RealMatrix& operator +=(const RealMatrix &rhs);

        (*@\Green{// subtract from *this}@*)
        RealMatrix& operator -= (const RealMatrix &T);

        (*@\Green{// multiplicate *this by a real matrix}@*)
        RealMatrix& operator *= (const RealMatrix &rhs);

        (*@\Green{// multiplicate *this by a constant}@*)
        RealMatrix& operator *= (const double &rhs);

        (*@\Green{// multiplicate *this by a real column vector from right}@*)
        RealMatrix& operator *=(const ColumnMatrix<double>& rhs);

        (*@\Green{// divide *this by a constant}@*)
        RealMatrix& operator /= (const double &rhs);

        (*@\Green{// sets each matrix element to zero}@*)
        void loadNullMatrix();

        (*@\Green{// loads a matrix with ones on the first main diagonal and zeros elsewhere}@*)
        void loadIdentityMatrix();

        (*@\Green{// returns the transpose of the stored matrix}@*)
        RealMatrix transpose() const;

        (*@\Green{// decides whether the matrix is a square one}@*)
        bool isSquare() const;

        (*@\Green{// decides whether the matrix consists of a single row}@*)
        bool isRowMatrix() const;

        (*@\Green{// decides whether the matrix consists of a single column}@*)
        bool isColumnMatrix() const;

        (*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        RealMatrix* clone() const;
    };

    (*@\Green{// Auxiliary real row and column matrix related arithmetical binary operators.}@*)

    double operator *(const RowMatrix<double> &lhs, const ColumnMatrix<double> &rhs);

    const RowMatrix<double> operator +(
        const RowMatrix<double> &lhs, const RowMatrix<double> &rhs);
    const RowMatrix<double> operator -(
        const RowMatrix<double> &lhs, const RowMatrix<double> &rhs);
    const RowMatrix<double> operator *(const RowMatrix<double> &lhs, const double &rhs);
    const RowMatrix<double> operator *(const double &lhs, const RowMatrix<double> &rhs);

    RowMatrix<double>& operator +=(RowMatrix<double> &lhs, const RowMatrix<double> &rhs);
    RowMatrix<double>& operator -=(RowMatrix<double> &lhs, const RowMatrix<double> &rhs);

    RowMatrix<double>& operator *=(RowMatrix<double> &lhs, const double &rhs);
    RowMatrix<double>& operator /=(RowMatrix<double> &lhs, const double &rhs);

    const ColumnMatrix<double> operator +(
        const ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs);
    const ColumnMatrix<double> operator -(
        const ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs);
    const ColumnMatrix<double> operator *(
        const ColumnMatrix<double> &lhs, const double &rhs);
    const ColumnMatrix<double> operator *(
        const double &lhs, const ColumnMatrix<double> &rhs);

    ColumnMatrix<double>& operator +=(
        ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs);
    ColumnMatrix<double>& operator -=(
        ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs);

    ColumnMatrix<double>& operator *=(ColumnMatrix<double> &lhs, const double &rhs);
    ColumnMatrix<double>& operator /=(ColumnMatrix<double> &lhs, const double &rhs);
}

#endif (*@\Green{// REALMATRICES\_H}@*)
