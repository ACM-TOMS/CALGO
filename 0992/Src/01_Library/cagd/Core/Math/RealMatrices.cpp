//----------------------------------------------------------------------------------
// File:        Core/Math/RealMatrices.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "RealMatrices.h"

#include "../Exceptions.h"

#include <new>
using namespace std;

namespace cagd
{
    //(*@\Green{// Implementation of general (rectangular) real matrices.}@*)

    //(*@\Green{// default/special constructor}@*)
    RealMatrix::RealMatrix(int row_count, int column_count):
        Matrix<double>(row_count, column_count)
    {
    }

    //(*@\Green{// matrix addition}@*)
    const RealMatrix RealMatrix::operator +(const RealMatrix &rhs) const
    {
        if (_row_count != rhs._row_count || _column_count != rhs._column_count)
        {
            throw Exception("RealMatrix::operator +(const RealMatrix &rhs) : "
                            "incompatible matrix dimensions!");
        }

        RealMatrix result(*this);

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            result._data[i] += rhs._data[i];
        }

        return result;
    }

    //(*@\Green{// matrix subtraction}@*)
    const RealMatrix RealMatrix::operator -(const RealMatrix &rhs) const
    {
        if (_row_count != rhs._row_count || _column_count != rhs._column_count)
        {
            throw Exception("RealMatrix::operator -(const RealMatrix &rhs) : "
                            "incompatible matrix dimensions!");
        }

        RealMatrix result(*this);

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            result._data[i] -= rhs._data[i];
        }

        return result;
    }

    //(*@\Green{// matrix multiplication}@*)
    const RealMatrix RealMatrix::operator *(const RealMatrix &rhs) const
    {
        if (_column_count != rhs._row_count)
        {
            throw Exception("RealMatrix::operator *(const RealMatrix &rhs) : "
                            "incompatible inner matrix dimensions!");
        }

        RealMatrix result(_row_count, rhs._column_count);

        #pragma omp parallel for
        for (int r_c = 0; r_c < _row_count * rhs._column_count; r_c++)
        {
            int r             = r_c / rhs._column_count;
            int c             = r_c % rhs._column_count;
            int offset        = r * _column_count;
            int result_offset = r * rhs._column_count + c;

            for (int i = 0; i < _column_count; i++)
            {
                result._data[result_offset] += _data[offset + i] *
                                               rhs._data[i * rhs._column_count + c];
            }
        }

        return result;
    }

    //(*@\Green{// multiplication by scalar from right}@*)
    const RealMatrix RealMatrix::operator *(const double &rhs) const
    {
        RealMatrix result(*this);

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            result._data[i] *= rhs;
        }

        return result;
    }

    //(*@\Green{// multiplication by scalar from left}@*)
    const RealMatrix operator *(const double &lhs, const RealMatrix &rhs)
    {
        RealMatrix result(rhs);

        #pragma omp parallel for
        for (int i = 0; i < result._row_count * result._column_count; i++)
        {
            result._data[i] *= lhs;
        }

        return result;
    }

    //(*@\Green{// multiplication by real column vector from right}@*)
    const RealMatrix RealMatrix::operator *(const ColumnMatrix<double> &rhs) const
    {
        if (_column_count != rhs.rowCount())
        {
            throw Exception("RealMatrix::operator *(const ColumnMatrix<double> &rhs): "
                            "incompatible inner matrix dimensions!");
        }

        RealMatrix result(_row_count, 1);

        #pragma omp parallel for
        for (int r = 0; r < _row_count; r++)
        {
            for (int c = 0, offset = r * _column_count; c < _column_count; c++)
            {
                result._data[r] += _data[offset + c] * rhs[c];
            }
        }

        return result;
    }

    //(*@\Green{// multiplication by real row vector from left}@*)
    const RealMatrix operator *(const RowMatrix<double> &lhs, const RealMatrix &rhs)
    {
        if (lhs.columnCount() != rhs._row_count)
        {
            throw Exception("operator *(const RowMatrix<double> &lhs, "
                            "const RealMatrix &rhs): "
                            "incompatible inner matrix dimensions!");
        }

        RealMatrix result(1, rhs._column_count);

        #pragma omp parallel for
        for (int c = 0; c < rhs._column_count; c++)
        {
            for (int r = 0; r < rhs._row_count; r++)
            {
                result._data[c] += lhs[r] * rhs._data[r * rhs._column_count + c];
            }
        }

        return result;
    }

    //(*@\Green{// multiplication of a real column matrix by a real row matrix from right}@*)
    const RealMatrix operator *(
        const ColumnMatrix<double> &lhs, const RowMatrix<double> &rhs)
    {
        RealMatrix result(lhs.rowCount(), rhs.columnCount());

        #pragma omp parallel for
        for (int r = 0; r < result._row_count; r++)
        {
            result.setRow(r, rhs);

            double multiplier = lhs[r];

            for (int c = 0, offset = r * result._column_count; c < result._column_count;
                 c++)
            {
                result._data[offset + c] *= multiplier;
            }
        }

        return result;
    }

    //(*@\Green{// multiplication of a real row matrix by a real matrix from right}@*)
    RowMatrix<double>& operator *=(RowMatrix<double> &lhs, const RealMatrix &rhs)
    {
        if (lhs.columnCount() != rhs._row_count)
        {
            throw Exception("operator *=(RowMatrix<double> &lhs, "
                            "const RealMatrix &rhs): "
                            "incompatible inner matrix dimensions!");
        }

        std::vector<double> result(rhs._column_count);

        #pragma omp parallel for
        for (int i = 0; i < rhs._column_count; i++)
        {
            for (int c = 0; c < lhs.columnCount(); c++)
            {
                result[i] += lhs[c] * rhs._data[c * rhs._column_count + i];
            }
        }

        lhs.resizeColumns(rhs._column_count);

        #pragma omp parallel for
        for (int i = 0; i < lhs.columnCount(); i++)
        {
            lhs[i] = result[i];
        }

        return lhs;
    }

    //(*@\Green{// division by scalar from right}@*)
    const RealMatrix RealMatrix::operator /(const double &rhs) const
    {
        if (!rhs)
        {
            throw Exception("RealMatrix::operator /(const double &rhs): "
                            "division by zero!");
        }

        RealMatrix result(*this);

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            result._data[i] /= rhs;
        }

        return result;
    }

    //(*@\Green{// add to *this}@*)
    RealMatrix& RealMatrix::operator +=(const RealMatrix &rhs)
    {
        if (_row_count != rhs._row_count || _column_count != rhs._column_count)
        {
            throw Exception("RealMatrix::operator +=(const RealMatrix &rhs) : "
                            "incompatible matrix dimensions!");
        }

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            _data[i] += rhs._data[i];
        }

        return *this;
    }

    //(*@\Green{// subtract from *this}@*)
    RealMatrix& RealMatrix::operator -=(const RealMatrix &rhs)
    {
        if (_row_count != rhs._row_count || _column_count != rhs._column_count)
        {
            throw Exception("RealMatrix::operator -=(const RealMatrix &rhs) : "
                            "incompatible matrix dimensions!");
        }

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            _data[i] -= rhs._data[i];
        }

        return *this;
    }

    //(*@\Green{// multiplicate *this by a real matrix}@*)
    RealMatrix& RealMatrix::operator *=(const RealMatrix &rhs)
    {
        if (_column_count != rhs._row_count)
        {
            throw Exception("RealMatrix::operator *=(const RealMatrix &rhs) : "
                            "incompatible inner matrix dimensions!");
        }

        RealMatrix result(_row_count, rhs._column_count);

        #pragma omp parallel for
        for (int r_c = 0; r_c < _row_count * rhs._column_count; r_c++)
        {
            int r             = r_c / rhs._column_count;
            int c             = r_c % rhs._column_count;
            int offset        = r * _column_count;
            int result_offset = r * rhs._column_count + c;

            for (int i = 0; i < _column_count; i++)
            {
                result._data[result_offset] += _data[offset + i] *
                                               rhs._data[i * rhs._column_count + c];
            }
        }

        _data.swap(result._data);

        _column_count = result._column_count;

        return *this;
    }

    //(*@\Green{// multiplicate *this by a constant}@*)
    RealMatrix& RealMatrix::operator *=(const double &rhs)
    {
        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            _data[i] *= rhs;
        }

        return *this;
    }

    //(*@\Green{// multiplicate *this by a real column vector from right}@*)
    RealMatrix& RealMatrix::operator *=(const ColumnMatrix<double>& rhs)
    {
        if (_column_count != rhs.rowCount())
        {
            throw Exception("RealMatrix::operator *=(const ColumnMatrix<double> &rhs): "
                            "incompatible inner matrix dimensions!");
        }

        vector<double> result(_row_count);

        #pragma omp parallel for
        for (int r = 0; r < _row_count; r++)
        {
            for (int c = 0, offset = r * _column_count; c < _column_count; c++)
            {
                result[r] += _data[offset + c] * rhs[c];
            }
        }

        _column_count = 1;
        _data.swap(result);

        return *this;
    }

    //(*@\Green{// divide *this by a constant}@*)
    RealMatrix& RealMatrix::operator /=(const double &rhs)
    {
        if (!rhs)
        {
            throw Exception("RealMatrix::operator /=(const double &rhs) : "
                            "division by zero!");
        }

        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            _data[i] /= rhs;
        }

        return *this;
    }

    //(*@\Green{// sets each matrix element to zero}@*)
    void RealMatrix::loadNullMatrix()
    {
        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            _data[i] = 0.0;
        }
    }

    //(*@\Green{// loads a matrix with ones on the first main diagonal and zeros elsewhere}@*)
    void RealMatrix::loadIdentityMatrix()
    {
        #pragma omp parallel for
        for (int i = 0; i < _row_count * _column_count; i++)
        {
            _data[i] = 0.0;
        }

        int step = _column_count + 1;

        #pragma omp parallel for
        for (int offset = 0; offset < min(_row_count, _column_count) * _column_count;
             offset += step)
        {
            _data[offset] = 1.0;
        }
    }

    //(*@\Green{// returns the transpose of the stored matrix}@*)
    RealMatrix RealMatrix::transpose() const
    {
        RealMatrix result(_column_count, _row_count);

        #pragma omp parallel for
        for (int r_c = 0; r_c < _row_count * _column_count; r_c++)
        {
            int r = r_c / _column_count;
            int c = r_c % _column_count;

            result._data[c * _row_count + r] = _data[r_c];
        }

        return result;
    }

    //(*@\Green{// decides whether the matrix is a square one}@*)
    bool RealMatrix::isSquare() const
    {
        return (_row_count == _column_count);
    }

    //(*@\Green{// decides whether the matrix consists of a single row}@*)
    bool RealMatrix::isRowMatrix() const
    {
        return (_row_count == 1);
    }

    //(*@\Green{// decides whether the matrix consists of a single column}@*)
    bool RealMatrix::isColumnMatrix() const
    {
        return (_column_count == 1);
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    RealMatrix* RealMatrix::clone() const
    {
        return new (nothrow) RealMatrix(*this);
    }

    //(*@\Green{// Implementation of auxiliary real row and column matrix related arithmetical binary operators.}@*)

    double operator *(const RowMatrix<double> &lhs, const ColumnMatrix<double> &rhs)
    {
        if (lhs.columnCount() != rhs.rowCount())
        {
            throw Exception("operator *(const RowMatrix<double> &lhs, "
                            "const ColumnMatrix<double> &rhs) : "
                            "incompatible inner matrix dimensions!");
        }

        double result = 0.0;

        #pragma omp parallel for reduction(+:result)
        for (int i = 0; i < lhs.columnCount(); i++)
        {
            result += lhs[i] * rhs[i];
        }

        return result;
    }

    const RowMatrix<double> operator +(
        const RowMatrix<double> &lhs, const RowMatrix<double> &rhs)
    {
        if (lhs.columnCount() != rhs.columnCount())
        {
            throw Exception("operator +(const RowMatrix<double> &lhs, "
                            "const RowMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        RowMatrix<double> result(lhs);

        #pragma omp parallel for
        for (int i = 0; i < result.columnCount(); i++)
        {
            result[i] += rhs[i];
        }

        return result;
    }

    const RowMatrix<double> operator -(
        const RowMatrix<double> &lhs, const RowMatrix<double> &rhs)
    {
        if (lhs.columnCount() != rhs.columnCount())
        {
            throw Exception("operator -(const RowMatrix<double> &lhs, "
                            "const RowMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        RowMatrix<double> result(lhs);

        #pragma omp parallel for
        for (int i = 0; i < result.columnCount(); i++)
        {
            result[i] -= rhs[i];
        }

        return result;
    }

    const RowMatrix<double> operator *(const RowMatrix<double> &lhs, const double &rhs)
    {
        RowMatrix<double> result(lhs);

        #pragma omp parallel for
        for (int i = 0; i < result.columnCount(); i++)
        {
            result[i] *= rhs;
        }

        return result;
    }

    const RowMatrix<double> operator *(const double &lhs, const RowMatrix<double> &rhs)
    {
        RowMatrix<double> result(rhs);

        #pragma omp parallel for
        for (int i = 0; i < result.columnCount(); i++)
        {
            result[i] *= lhs;
        }

        return result;
    }

    RowMatrix<double>& operator +=(RowMatrix<double> &lhs, const RowMatrix<double> &rhs)
    {
        if (lhs.columnCount() != rhs.columnCount())
        {
            throw Exception("operator +=(RowMatrix<double> &lhs, "
                            "const RowMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        #pragma omp parallel for
        for (int i = 0; i < lhs.columnCount(); i++)
        {
            lhs[i] += rhs[i];
        }

        return lhs;
    }

    RowMatrix<double>& operator -=(RowMatrix<double> &lhs, const RowMatrix<double> &rhs)
    {
        if (lhs.columnCount() != rhs.columnCount())
        {
            throw Exception("operator -=(RowMatrix<double> &lhs, "
                            "const RowMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        #pragma omp parallel for
        for (int i = 0; i < lhs.columnCount(); i++)
        {
            lhs[i] -= rhs[i];
        }

        return lhs;
    }

    RowMatrix<double>& operator *=(RowMatrix<double> &lhs, const double &rhs)
    {
        #pragma omp parallel for
        for (int i = 0; i < lhs.columnCount(); i++)
        {
            lhs[i] *= rhs;
        }

        return lhs;
    }

    RowMatrix<double>& operator /=(RowMatrix<double> &lhs, const double &rhs)
    {
        #pragma omp parallel for
        for (int i = 0; i < lhs.columnCount(); i++)
        {
            lhs[i] /= rhs;
        }

        return lhs;
    }

    const ColumnMatrix<double> operator +(
        const ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs)
    {
        if (lhs.rowCount() != rhs.rowCount())
        {
            throw Exception("operator +(const ColumnMatrix<double> &lhs, "
                            "const ColumnMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        ColumnMatrix<double> result(lhs);

        #pragma omp parallel for
        for (int i = 0; i < result.rowCount(); i++)
        {
            result[i] += rhs[i];
        }

        return result;
    }

    const ColumnMatrix<double> operator -(
        const ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs)
    {
        if (lhs.rowCount() != rhs.rowCount())
        {
            throw Exception("operator -(const ColumnMatrix<double> &lhs, "
                            "const ColumnMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        ColumnMatrix<double> result(lhs);

        #pragma omp parallel for
        for (int i = 0; i < result.rowCount(); i++)
        {
            result[i] -= rhs[i];
        }

        return result;
    }

    const ColumnMatrix<double> operator *(
        const ColumnMatrix<double> &lhs, const double &rhs)
    {
        ColumnMatrix<double> result(lhs);

        #pragma omp parallel for
        for (int i = 0; i < result.rowCount(); i++)
        {
            result[i] *= rhs;
        }

        return result;
    }

    const ColumnMatrix<double> operator *(
        const double &lhs, const ColumnMatrix<double> &rhs)
    {
        ColumnMatrix<double> result(rhs);

        #pragma omp parallel for
        for (int i = 0; i < result.rowCount(); i++)
        {
            result[i] *= lhs;
        }

        return result;
    }

    ColumnMatrix<double>& operator +=(
        ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs)
    {
        if (lhs.rowCount() != rhs.rowCount())
        {
            throw Exception("operator +=(ColumnMatrix<double> &lhs, "
                            "const ColumnMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        #pragma omp parallel for
        for (int i = 0; i < lhs.rowCount(); i++)
        {
            lhs[i] += rhs[i];
        }

        return lhs;
    }

    ColumnMatrix<double>& operator -=(
        ColumnMatrix<double> &lhs, const ColumnMatrix<double> &rhs)
    {
        if (lhs.rowCount() != rhs.rowCount())
        {
            throw Exception("operator -=(ColumnMatrix<double> &lhs, "
                            "const ColumnMatrix<double> &rhs) : "
                            "incompatible matrix dimensions!");
        }

        #pragma omp parallel for
        for (int i = 0; i < lhs.rowCount(); i++)
        {
            lhs[i] -= rhs[i];
        }

        return lhs;
    }

    ColumnMatrix<double>& operator *=(ColumnMatrix<double> &lhs, const double &rhs)
    {
        #pragma omp parallel for
        for (int i = 0; i < lhs.rowCount(); i++)
        {
            lhs[i] *= rhs;
        }

        return lhs;
    }

    ColumnMatrix<double>& operator /=(ColumnMatrix<double> &lhs, const double &rhs)
    {
        #pragma omp parallel for
        for (int i = 0; i < lhs.rowCount(); i++)
        {
            lhs[i] /= rhs;
        }

        return lhs;
    }
}
