//----------------------------------------------------------------------------------
// File:        Core/Math/Matrices.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef MATRICES_H
#define MATRICES_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace cagd
{
    //(*@\Green{// forward declaration of the template class Matrix}@*)
    template <typename T>
    class Matrix;

    //(*@\Green{// forward declaration of the template class RowMatrix}@*)
    template <typename T>
    class RowMatrix;

    //(*@\Green{// forward declaration of the template class ColumnMatrix}@*)
    template <typename T>
    class ColumnMatrix;

    //(*@\Green{// forward declaration of the template class TriangularMatrix}@*)
    template <typename T>
    class TriangularMatrix;

    //(*@\Green{// forward declarations of overloaded and templated input/output from/to stream operators}@*)
    template <typename T>
    std::ostream& operator << (std::ostream& lhs, const Matrix<T>& rhs);

    template <typename T>
    std::istream& operator >>(std::istream& lhs, Matrix<T>& rhs);

    template <typename T>
    std::ostream& operator << (std::ostream& lhs, const TriangularMatrix<T>& rhs);

    template <typename T>
    std::istream& operator >>(std::istream& lhs, TriangularMatrix<T>& rhs);

    //(*@\Green{// Interface of the template class Matrix:}@*)
    template <typename T>
    class Matrix
    {
        friend std::ostream& cagd::operator << <T>(std::ostream&, const Matrix<T>& rhs);
        friend std::istream& cagd::operator >> <T>(std::istream&, Matrix<T>& rhs);

    protected:
        int            _row_count;
        int            _column_count;
        std::vector<T> _data;

    public:
        //(*@\Green{// default/special constructor}@*)
        explicit Matrix(int row_count = 1, int column_count = 1);

        //(*@\Green{// get element by non-constant reference}@*)
        T& operator ()(int row, int column);

        //(*@\Green{// get element by constant reference}@*)
        const T& operator ()(int row, int column) const;

        //(*@\Green{// get dimensions}@*)
        int rowCount() const;
        int columnCount() const;

        //(*@\Green{// set dimensions}@*)
        virtual bool resizeRows(int row_count);
        virtual bool resizeColumns(int column_count);

        //(*@\Green{// update}@*)
        bool setRow(int index, const RowMatrix<T>& row);
        bool setColumn(int index, const ColumnMatrix<T>& column);

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual Matrix* clone() const;

        //(*@\Green{// destructor}@*)
        virtual ~Matrix();
    };

    //(*@\Green{// Interface of the template class RowMatrix:}@*)
    template <typename T>
    class RowMatrix: public Matrix<T>
    {
    public:
        //(*@\Green{// default/special constructor}@*)
        explicit RowMatrix(int column_count = 1);

        //(*@\Green{// get element by non-constant reference}@*)
        T& operator ()(int column);
        T& operator [](int column);

        //(*@\Green{// get element by constant reference}@*)
        const T& operator ()(int column) const;
        const T& operator [](int column) const;

        //(*@\Green{// a row matrix consists of a single row}@*)
        bool resizeRows(int row_count);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        RowMatrix* clone() const;
    };

    //(*@\Green{// Interface of the template class ColumnMatrix:}@*)
    template <typename T>
    class ColumnMatrix: public Matrix<T>
    {
    public:
        //(*@\Green{// default/special constructor}@*)
        explicit ColumnMatrix(int row_count = 1);

        //(*@\Green{// get element by non-constant reference}@*)
        T& operator ()(int row);
        T& operator [](int row);

        //(*@\Green{// get element by constant reference}@*)
        const T& operator ()(int row) const;
        const T& operator [](int row) const;

        //(*@\Green{// a column matrix consists of a single column}@*)
        bool resizeColumns(int column_count);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        ColumnMatrix* clone() const;
    };

    //(*@\Green{// Interface of the template class TriangularMatrix:}@*)
    template <typename T>
    class TriangularMatrix
    {
        friend std::ostream& cagd::operator << <T>(std::ostream&,
                                                   const TriangularMatrix<T>& rhs);
        friend std::istream& cagd::operator >> <T>(std::istream&,
                                                   TriangularMatrix<T>& rhs);

    protected:
        int                           _row_count;
        std::vector< std::vector<T> > _data;

    public:
        //(*@\Green{// default/special constructor}@*)
        explicit TriangularMatrix(int row_count = 1);

        //(*@\Green{// get element by non-constant reference}@*)
        T& operator ()(int row, int column);

        //(*@\Green{// get element by constant reference}@*)
        const T& operator ()(int row, int column) const;

        //(*@\Green{// get dimension}@*)
        int rowCount() const;

        //(*@\Green{// set dimension}@*)
        virtual bool resizeRows(int row_count);

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual TriangularMatrix* clone() const;

        //(*@\Green{// destructor}@*)
        virtual ~TriangularMatrix();
    };

    //(*@\Green{// Implementation of the template class Matrix:}@*)

    //(*@\Green{// default/special constructor}@*)
    template <typename T>
    Matrix<T>::Matrix(int row_count, int column_count):
        _row_count(row_count < 0 || column_count <= 0 ? 0 : row_count),
        _column_count(column_count < 0 || !_row_count ? 0 : column_count),
        _data(_row_count * _column_count)
    {
        assert("The row count of a matrix should be non-negative!" && row_count >= 0);
        assert("The column count of a matrix should be non-negative!" &&
               column_count >= 0);
    }

    //(*@\Green{// get element by non-constant reference}@*)
    template <typename T>
    inline T& Matrix<T>::operator ()(int row, int column)
    {
        assert("The given row index is out of bounds!" && (row >= 0 && row < _row_count));
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < _column_count));
        return _data[row * _column_count + column];
    }

    //(*@\Green{// get element by constant reference}@*)
    template <typename T>
    inline const T& Matrix<T>::operator ()(int row, int column) const
    {
        assert("The given row index is out of bounds!" && (row >= 0 && row < _row_count));
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < _column_count));
        return _data[row * _column_count + column];
    }

    //(*@\Green{// get dimensions}@*)
    template <typename T>
    inline int Matrix<T>::rowCount() const
    {
        return _row_count;
    }

    template <typename T>
    inline int Matrix<T>::columnCount() const
    {
        return _column_count;
    }

    //(*@\Green{// set dimensions}@*)
    template <typename T>
    bool Matrix<T>::resizeRows(int row_count)
    {
        assert("The row count of a matrix should be non-negative!" && row_count >= 0);

        if (row_count < 0)
        {
            return false;
        }

        if (_row_count != row_count)
        {
            _data.resize(row_count * _column_count);
            _row_count = row_count;

            if (_row_count == 0)
            {
                _column_count = 0;
            }
        }

        return false;
    }

    template <typename T>
    bool Matrix<T>::resizeColumns(int column_count)
    {
        assert("The column count of a matrix should be non-negative!" &&
               column_count >= 0);

        if (column_count < 0)
        {
            return false;
        }

        if (_column_count != column_count)
        {
            std::vector<T> new_data(_row_count * column_count);

            int preserved_column_count = std::min(column_count, _column_count);

            #pragma omp parallel for
            for (int r = 0; r < _row_count; r++)
            {
                for (int c = 0,
                     old_offset = r * _column_count,
                     new_offset = r * column_count;
                     c < preserved_column_count;
                     c++)
                {
                    int old_index       = old_offset + c;
                    int new_index       = new_offset + c;
                    new_data[new_index] = _data[old_index];
                }
            }

            _data.swap(new_data);

            _column_count = column_count;

            if (_column_count == 0)
            {
                _row_count = 0;
            }
        }

        return true;
    }

    //(*@\Green{// update}@*)
    template <class T>
    bool Matrix<T>::setRow(int row, const RowMatrix<T>& entire_row)
    {
        assert("The given column index is out of bounds!" &&
               (row >= 0 && row < _row_count));
        assert("The column counts of the underlying row matrices should be equal!" &&
                _column_count == entire_row.columnCount());

        if (row < 0 || row >= _row_count || _column_count != entire_row.columnCount())
        {
            return false;
        }

        int offset = row * _column_count;

        #pragma omp parallel for
        for (int c = 0; c < _column_count; c++)
        {
            _data[offset + c] = entire_row._data[c];
        }

        return true;
    }

    template <class T>
    bool Matrix<T>::setColumn(int column, const ColumnMatrix<T>& entire_column)
    {
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < _column_count));
        assert("The row counts of the underlying column matrices should be equal!" &&
               _row_count == entire_column.rowCount());

        if (column < 0 || column >= _column_count ||
            _row_count != entire_column.rowCount())
        {
            return false;
        }

        #pragma omp parallel for
        for (int r = 0; r < _row_count; r++)
        {
            _data[r * _column_count + column] = entire_column._data[r];
        }

        return true;
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    template <typename T>
    Matrix<T>* Matrix<T>::clone() const
    {
        return new (std::nothrow) Matrix<T>(*this);
    }

    //(*@\Green{// destructor}@*)
    template <typename T>
    Matrix<T>::~Matrix()
    {
        _row_count = _column_count = 0;
        _data.clear();
    }

    //(*@\Green{// Implementation of the template class RowMatrix:}@*)

    //(*@\Green{// default/special constructor}@*)
    template <typename T>
    RowMatrix<T>::RowMatrix(int column_count): Matrix<T>(1, column_count)
    {
    }

    //(*@\Green{// get element by non-constant reference}@*)
    template <typename T>
    inline T& RowMatrix<T>::operator ()(int column)
    {
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < this->_column_count));
        return this->_data[column];
    }

    template <typename T>
    inline T& RowMatrix<T>::operator [](int column)
    {
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < this->_column_count));
        return this->_data[column];
    }

    //(*@\Green{// get element by constant reference}@*)
    template <typename T>
    inline const T& RowMatrix<T>::operator ()(int column) const
    {
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < this->_column_count));
        return this->_data[column];
    }

    template <typename T>
    inline const T& RowMatrix<T>::operator [](int column) const
    {
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < this->_column_count));
        return this->_data[column];
    }

    //(*@\Green{// a row matrix consists of a single row}@*)
    template <typename T>
    bool RowMatrix<T>::resizeRows(int row_count)
    {
        return (row_count == 1);
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    template <typename T>
    RowMatrix<T>* RowMatrix<T>::clone() const
    {
        return new (std::nothrow) RowMatrix<T>(*this);
    }

    //(*@\Green{// Implementation of the template class ColumnMatrix:}@*)

    //(*@\Green{// default/special constructor}@*)
    template <typename T>
    ColumnMatrix<T>::ColumnMatrix(int row_count): Matrix<T>(row_count, 1)
    {
    }

    //(*@\Green{// get element by non-constant reference}@*)
    template <typename T>
    inline T& ColumnMatrix<T>::operator ()(int row)
    {
        assert("The given row index is out of bounds!" &&
               (row >= 0 && row < this->_row_count));
        return this->_data[row];
    }

    template <typename T>
    inline T& ColumnMatrix<T>::operator [](int row)
    {
        assert("The given row index is out of bounds!" &&
               (row >= 0 && row < this->_row_count));
        return this->_data[row];
    }

    //(*@\Green{// get element by constant reference}@*)
    template <typename T>
    inline const T& ColumnMatrix<T>::operator ()(int row) const
    {
        assert("The given row index is out of bounds!" &&
               (row >= 0 && row < this->_row_count));
        return this->_data[row];
    }

    template <typename T>
    inline const T& ColumnMatrix<T>::operator [](int row) const
    {
        assert("The given row index is out of bounds!" &&
               (row >= 0 && row < this->_row_count));
        return this->_data[row];
    }

    //(*@\Green{// a column matrix consists of a single column}@*)
    template <typename T>
    bool ColumnMatrix<T>::resizeColumns(int column_count)
    {
        return (column_count == 1);
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    template <typename T>
    ColumnMatrix<T>* ColumnMatrix<T>::clone() const
    {
        return new (std::nothrow) ColumnMatrix<T>(*this);
    }

    //(*@\Green{// Implementation of the template class TriangularMatrix:}@*)

    //(*@\Green{// default/special constructor}@*)
    template <typename T>
    TriangularMatrix<T>::TriangularMatrix(int row_count):
        _row_count(row_count < 0 ? 0 : row_count),
        _data(_row_count)
    {
        assert("The row count of a triangular matrix should be a non-negative "
               "integer!" && row_count >= 0);

        for (int r = 0; r < _row_count; r++)
        {
            _data[r].resize(r + 1);
        }
    }

    //(*@\Green{// get element by non-constant reference}@*)
    template <typename T>
    inline T& TriangularMatrix<T>::operator ()(int row, int column)
    {
        assert("The given row index is out of bounds!" && (row >= 0 && row < _row_count));
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column <= row));
        return _data[row][column];
    }

    //(*@\Green{// get element by constant reference}@*)
    template <typename T>
    inline const T& TriangularMatrix<T>::operator ()(int row, int column) const
    {
        assert("The given row index is out of bounds!" && (row >= 0 && row < _row_count));
        assert("The given column index is out of bounds!" &&
               (column >= 0 && column <= row));
        return _data[row][column];
    }

    //(*@\Green{// get dimension}@*)
    template <typename T>
    inline int TriangularMatrix<T>::rowCount() const
    {
        return _row_count;
    }

    //(*@\Green{// set dimension}@*)
    template <typename T>
    bool TriangularMatrix<T>::resizeRows(int row_count)
    {
        assert("The row count of a triangular matrix should be a non-negative "
               "integer!" && row_count >= 0);

        if (row_count < 0)
        {
            return false;
        }

        if (_row_count != row_count)
        {

            _data.resize(row_count);

            for (int r = _row_count; r < row_count; r++)
            {
                _data[r].resize(r + 1);
            }

            _row_count = row_count;
        }

        return true;
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    template <typename T>
    TriangularMatrix<T>* TriangularMatrix<T>::clone() const
    {
        return new (std::nothrow) TriangularMatrix<T>(*this);
    }

    //(*@\Green{// destructor}@*)
    template <typename T>
    TriangularMatrix<T>::~TriangularMatrix()
    {
        _row_count = 0;
        _data.clear();
    }

    //(*@\Green{// Definitions of overloaded and templated input/output from/to stream operators:}@*)

    //(*@\Green{// output to stream}@*)
    template <typename T>
    std::ostream& operator <<(std::ostream& lhs, const Matrix<T>& rhs)
    {
        lhs << rhs._row_count << " " << rhs._column_count << std::endl;

        for (int r = 0; r < rhs._row_count; r++)
        {
            for (int c = 0, offset = r * rhs._column_count; c < rhs._column_count; c++)
            {
                lhs << rhs._data[offset + c] << " ";
            }

            lhs << std::endl;
        }

        return lhs;
    }

    //(*@\Green{// input from stream}@*)
    template <typename T>
    std::istream& operator >>(std::istream& lhs, Matrix<T>& rhs)
    {
        lhs >> rhs._row_count >> rhs._column_count;

        int size = rhs._row_count * rhs._column_count;

        rhs._data.resize(size);

        for (int i = 0; i < size; i++)
        {
            lhs >> rhs._data[i];
        }

        return lhs;
    }

    //(*@\Green{// output to stream}@*)
    template <typename T>
    std::ostream& operator <<(std::ostream& lhs, const TriangularMatrix<T>& rhs)
    {
        lhs << rhs._row_count << std::endl;

        for (int r = 0; r < rhs._row_count; r++)
        {
            for (int c = 0; c <= r; c++)
            {
                lhs << rhs._data[r][c] << " ";
            }
            lhs << std::endl;
        }

        return lhs;
    }

    //(*@\Green{// input from stream}@*)
    template <typename T>
    std::istream& operator >>(std::istream& lhs, TriangularMatrix<T>& rhs)
    {
        lhs >> rhs._row_count;

        rhs._data.resize(rhs._row_count);
        for (int r = 0; r < rhs._row_count; r++)
        {
            rhs._data[r].resize(r + 1);

            for (int c = 0; c <= r; r++)
            {
                lhs >> rhs._data[r][c];
            }
        }

        return lhs;
    }
}

#endif //(*@\Green{// MATRICES\_H}@*)
