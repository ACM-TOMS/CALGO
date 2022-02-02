//---------------------------------------------------------------------------------------
// File:        Core/Math/RealMatrixDecompositions.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
// Reference:   W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery. 2007.
//              Numerical Recipes 3rd Edition: The art of Scientific Computing.
//              Cambridge University Press, New York.
//---------------------------------------------------------------------------------------

#ifndef REALMATRIXDECOMPOSITIONS_H
#define REALMATRIXDECOMPOSITIONS_H

#include "../Exceptions.h"
#include "Matrices.h"
#include "RealMatrices.h"
#include <cmath>
#include <limits>

namespace cagd
{
    //(*@\Green{// Generates and stores the pivoted Doolittle-type $PLU$ decomposition of a given regular square matrix.}@*)
    class PLUDecomposition
    {
    private:
        bool           _decomposition_is_done;
        double         _determinant;
        RowMatrix<int> _P;
        RealMatrix     _LU;

    public:
        //(*@\Green{// Tries to determine the pivoted Doolittle-type $LU$ decomposition of the given real}@*)
        //(*@\Green{// square matrix $M$. If $M$ is singular then its $LU$ factorization cannot be determined.}@*)
        //(*@\Green{// If the decomposition is unsuccessful, the boolean variable \_decomposition\_is\_done}@*)
        //(*@\Green{// will be set to false, and if the input variable generate\_exceptions is set to true,}@*)
        //(*@\Green{// then the method will also generate an exception with a meaningful reason.}@*)
        PLUDecomposition(const RealMatrix &M, bool generate_exceptions = false);

        //(*@\Green{// Returns either the successful or the unsuccessful state of the decomposition.}@*)
        bool   isCorrect() const;

        //(*@\Green{// If the decomposition was successful the method will return the product of the diagonal entries of}@*)
        //(*@\Green{// the upper triangular matrix $U$ (i.e., the determinant of $M$), otherwise it will return zero.}@*)
        double determinant() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        PLUDecomposition* clone() const;

        //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $M \cdot X = B$, where $M$ is a regular}@*)
        //(*@\Green{// square matrix, while $B$ can be either a column or rectangular matrix of multiple columns with the same}@*)
        //(*@\Green{// row dimension as $M$. In case of success, the solution(s) will be stored in $X$.}@*)
        //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
        //(*@\Green{// the transpose of $B$ should be equal to the row dimension of $M$.}@*)
        //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
        //(*@\Green{// mathematical and boolean operators.}@*)
        template <typename T>
        bool solveLinearSystem(const Matrix<T> &B, Matrix<T> &X,
                               bool represent_solutions_as_columns = true);
    };

    //(*@\Green{// Generates and stores the factors $L$ and $U$ of the unpivoted Doolittle-type $LU$ decomposition of a}@*)
    //(*@\Green{// given regular matrix.}@*)
    class FactorizedUnpivotedLUDecomposition
    {
    private:
        bool       _decomposition_is_done;
        double     _determinant;
        RealMatrix _L, _U;

    public:
        //(*@\Green{// The constructor tries to determine the unpivoted Doolittle-type $LU$ decomposition of the given real}@*)
        //(*@\Green{// square matrix $M$. If $M$ is singular then its $LU$ factorization cannot be determined.}@*)
        //(*@\Green{// If the decomposition is unsuccessful, the boolean variable \_decomposition\_is\_done will be}@*)
        //(*@\Green{// set to false, and if the input variable generate\_exceptions is set to true, then the method will}@*)
        //(*@\Green{// also generate an exception with a meaningful reason.}@*)
        FactorizedUnpivotedLUDecomposition(const RealMatrix &M,
                                           bool generate_exceptions = false);

        //(*@\Green{// Returns either the successful or the unsuccessful state of the decomposition.}@*)
        bool isCorrect() const;

        //(*@\Green{// If the decomposition was successful the method will return the product of the diagonal entries of}@*)
        //(*@\Green{// the upper triangular matrix $U$ (i.e., the determinant of $M$), otherwise it will return zero.}@*)
        double determinant() const;

        //(*@\Green{// Returns a constant reference to the lower triangular matrix $L$.}@*)
        const RealMatrix& L() const;

        //(*@\Green{// Returns a constant reference to the upper triangular matrix $U$.}@*)
        const RealMatrix& U() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        FactorizedUnpivotedLUDecomposition* clone() const;


        //(*@\Green{// If $M$ is a regular real square matrix (i.e., the decomposition $M = L \cdot U$ can be constructed),}@*)
        //(*@\Green{// then the next two methods can be used to solve the system $M \cdot X = B$ of linear equations in}@*)
        //(*@\Green{// two steps: at first one has to solve the system $L \cdot Y = B$, then one should solve the system}@*)
        //(*@\Green{// $U \cdot X = Y$.}@*)

        //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $L \cdot Y = B$, where $B$ can}@*)
        //(*@\Green{// be either a column or rectangular matrix of multiple columns with the same row dimension as $L$.}@*)
        //(*@\Green{// In case of success, the solution(s) will be stored in $Y$.}@*)
        //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
        //(*@\Green{// the transpose of B should be equal to the row dimension of $L$.}@*)
        //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
        //(*@\Green{// mathematical and boolean operators.}@*)
        template <typename T>
        bool solveLLinearSystem(const Matrix<T> &B, Matrix<T> &Y,
                                bool represent_solutions_as_columns = true);

        //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $U \cdot X = Y$, where $Y$ can}@*)
        //(*@\Green{// be either a column or rectangular matrix of multiple columns with the same row dimension as $U$.}@*)
        //(*@\Green{// In case of success, the solution(s) will be stored in $X$.}@*)
        //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
        //(*@\Green{// the transpose of Y should be equal to the row dimension of $U$.}@*)
        //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
        //(*@\Green{// mathematical and boolean operators.}@*)
        template <typename T>
        bool solveULinearSystem(const Matrix<T> &Y, Matrix<T> &X,
                                bool represent_solutions_as_columns = true);
    };
    
    //(*@\Green{// Generates and stores the singular value decomposition of a given real (not necessarily square) matrix $M$.}@*)
    class SVDecomposition
    {
    private:
        bool              _decomposition_is_done;
        double            _product_of_singular_values;
        RealMatrix        _U;
        RowMatrix<double> _S;
        RealMatrix        _V;

    public:
        //(*@\Green{// Tries to determine the singular value decomposition of the given real matrix $M$.}@*)
        //(*@\Green{// If the decomposition is unsuccessful, the boolean variable \_decomposition\_is\_done}@*)
        //(*@\Green{// will be set to false, and if the input variable generate\_exceptions is set to true,}@*)
        //(*@\Green{// then the method will also generate an exception with a meaningful reason.}@*)
        SVDecomposition(const RealMatrix &M, bool generate_exceptions = false);

        //(*@\Green{// Returns either the successful or the unsuccessful state of the decomposition.}@*)
        bool   isCorrect() const;

        //(*@\Green{// Returns the ratio of the largest and smalles singular values.}@*)
        double conditionNumber() const;

        //(*@\Green{// Returns the ratio of the smallest and largest singular values.}@*)
        double reciprocalConditionNumber() const;

        //(*@\Green{// Returns the product of the obtained singular values. If $M$ is a regular real square matrix,}@*)
        //(*@\Green{// then this product coincides with the determinant of $M$.}@*)
        double productOfSingularValues() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        SVDecomposition* clone() const;
        
        //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $M \cdot X = B$, where $M$ is not}@*)
        //(*@\Green{// necessarily a real square matrix, while $B$ can be either a column or rectangular matrix of multiple}@*)
        //(*@\Green{// columns with the same row dimension as $M$. In case of success, the solution(s) will be stored in $X$.}@*)
        //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
        //(*@\Green{// the transpose of $B$ should be equal to the row dimension of $M$.}@*)
        //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
        //(*@\Green{// mathematical and boolean operators.}@*)
        template <typename T>
        bool solveLinearSystem(const Matrix<T> &B, Matrix<T> &X,
                               bool represent_solutions_as_columns = true);
    };

    //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $M \cdot X = B$, where $M$ is a regular}@*)
    //(*@\Green{// square matrix, while $B$ can be either a column or rectangular matrix of multiple columns with the same row}@*)
    //(*@\Green{// dimension as $M$. In case of success, the solution(s) will be stored in $X$.}@*)
    //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
    //(*@\Green{// the transpose of $B$ should be equal to the row dimension of $M$.}@*)
    //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
    //(*@\Green{// mathematical and boolean operators.}@*)
    template <typename T>
    bool PLUDecomposition::solveLinearSystem(const Matrix<T> &B, Matrix<T> &X,
                                             bool represent_solutions_as_columns)
    {
        if (!_decomposition_is_done)
        {
            return false;
        }

        if (represent_solutions_as_columns)
        {
            int dimension = _LU.columnCount();

            if (B.rowCount() != dimension)
            {
                return false;
            }

            X = B;

            #pragma omp parallel for
            for (int k = 0; k < B.columnCount(); k++)
            {
                int ii = 0;

                for (int i = 0; i < dimension; i++)
                {
                    int ip = _P[i];
                    T sum    = X(ip, k);
                    X(ip, k) = X(i, k);

                    if (ii != 0)
                    {
                        for (int j = ii - 1; j < i; j++)
                        {
                            sum -= _LU(i, j) * X(j, k);
                        }
                    }
                    else
                    {
                        if (sum != 0.0)
                        {
                            ii = i + 1;
                        }
                    }

                    X(i, k) = sum;
                }

                for (int i = dimension - 1; i >= 0; i--)
                {
                    T sum = X(i, k);

                    for (int j = i + 1; j < dimension; j++)
                    {
                        sum -= _LU(i, j) * X(j, k);
                    }

                    X(i, k) = sum /= _LU(i, i);
                }
            }
        }
        else
        {
            int dimension = _LU.rowCount();
            if (B.columnCount() != dimension)
            {
                return false;
            }

            X = B;

            #pragma omp parallel for
            for (int k = 0; k < B.rowCount(); k++)
            {
                int ii = 0;

                for (int i = 0; i < dimension; i++)
                {
                    int ip = _P[i];
                    T sum    = X(k, ip);
                    X(k, ip) = X(k, i);

                    if (ii != 0)
                    {
                        for (int j = ii - 1; j < i; j++)
                        {
                            sum -= _LU(i, j) * X(k, j);
                        }
                    }
                    else
                    {
                        if (sum != 0.0)
                        {
                            ii = i + 1;
                        }
                    }

                    X(k, i) = sum;
                }

                for (int i = dimension - 1; i >= 0; i--)
                {
                    T sum = X(k, i);

                    for (int j = i + 1; j < dimension; j++)
                    {
                        sum -= _LU(i, j) * X(k, j);
                    }

                    X(k, i) = sum /= _LU(i, i);
                }
            }
        }

        return true;
    }

    //(*@\Green{// If $M$ is a regular real square matrix (i.e., the decomposition $M = L \cdot U$ can be constructed),}@*)
    //(*@\Green{// then the next two methods can be used to solve the system $M \cdot X = B$ of linear equations in}@*)
    //(*@\Green{// two steps: at first one has to solve the system $L \cdot Y = B$, then one should solve the system}@*)
    //(*@\Green{// $U \cdot X = Y$.}@*)

    //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $L \cdot Y = B$, where $B$ can}@*)
    //(*@\Green{// be either a column or rectangular matrix of multiple columns with the same row dimension as $L$.}@*)
    //(*@\Green{// In case of success, the solution(s) will be stored in $Y$.}@*)
    //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
    //(*@\Green{// the transpose of B should be equal to the row dimension of $L$.}@*)
    //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
    //(*@\Green{// mathematical and boolean operators.}@*)
    template <typename T>
    bool FactorizedUnpivotedLUDecomposition::solveLLinearSystem(
            const Matrix<T> &B, Matrix<T> &Y, bool represent_solutions_as_columns)
    {
        if (!_decomposition_is_done)
        {
            return false;
        }

        int dimension = _L.columnCount();

        RowMatrix<int> P(dimension);

        for (int i = 0; i < dimension; i++)
        {
            P[i] = i;
        }

        if (represent_solutions_as_columns)
        {
            if (B.rowCount() != dimension)
            {
                return false;
            }

            Y = B;

            #pragma omp parallel for
            for (int k = 0; k < B.columnCount(); k++)
            {
                int ii = 0;

                for (int i = 0; i < dimension; i++)
                {
                    int ip = P[i];
                    T sum    = Y(ip, k);
                    Y(ip, k) = Y(i, k);

                    if (ii != 0)
                    {
                        for (int j = ii - 1; j < i; j++)
                        {
                            sum -= _L(i, j) * Y(j, k);
                        }
                    }
                    else
                    {
                        if (sum != 0.0)
                        {
                            ii = i + 1;
                        }
                    }

                    Y(i, k) = sum;
                }

                for (int i = dimension - 1; i >= 0; i--)
                {
                    T sum = Y(i, k);

                    for (int j = i + 1; j < dimension; j++)
                    {
                        sum -= _L(i, j) * Y(j, k);
                    }

                    Y(i, k) = sum /= _L(i, i);
                }
            }
        }
        else
        {
            if (B.columnCount() != dimension)
            {
                return false;
            }

            Y = B;

            #pragma omp parallel for
            for (int k = 0; k < B.rowCount(); k++)
            {
                int ii = 0;

                for (int i = 0; i < dimension; i++)
                {
                    int ip = P[i];
                    T sum    = Y(k, ip);
                    Y(k, ip) = Y(k, i);

                    if (ii != 0)
                    {
                        for (int j = ii - 1; j < i; j++)
                        {
                            sum -= _L(i, j) * Y(k, j);
                        }
                    }
                    else
                    {
                        if (sum != 0.0)
                        {
                            ii = i + 1;
                        }
                    }

                    Y(k, i) = sum;
                }

                for (int i = dimension - 1; i >= 0; i--)
                {
                    T sum = Y(k, i);

                    for (int j = i + 1; j < dimension; j++)
                    {
                        sum -= _L(i, j) * Y(k, j);
                    }

                    Y(k, i) = sum /= _L(i, i);
                }
            }
        }

        return true;
    }

    //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $U \cdot X = Y$, where $Y$ can}@*)
    //(*@\Green{// be either a column or rectangular matrix of multiple columns with the same row dimension as $U$.}@*)
    //(*@\Green{// In case of success, the solution(s) will be stored in $X$.}@*)
    //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
    //(*@\Green{// the transpose of Y should be equal to the row dimension of $U$.}@*)
    //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
    //(*@\Green{// mathematical and boolean operators.}@*)
    template <typename T>
    bool FactorizedUnpivotedLUDecomposition::solveULinearSystem(
            const Matrix<T> &Y, Matrix<T> &X, bool represent_solutions_as_columns)
    {
        if (!_decomposition_is_done)
        {
            return false;
        }

        int dimension = _U.columnCount();

        RowMatrix<int> P(dimension);

        for (int i = 0; i < dimension; i++)
        {
            P[i] = i;
        }

        if (represent_solutions_as_columns)
        {
            if (Y.rowCount() != dimension)
            {
                return false;
            }

            X = Y;

            #pragma omp parallel for
            for (int k = 0; k < Y.columnCount(); k++)
            {
                int ii = 0;

                for (int i = 0; i < dimension; i++)
                {
                    int ip = P[i];
                    T sum    = X(ip, k);
                    X(ip, k) = X(i, k);

                    if (ii != 0)
                    {
                        for (int j = ii - 1; j < i; j++)
                        {
                            sum -= _U(i, j) * X(j, k);
                        }
                    }
                    else
                    {
                        if (sum != 0.0)
                        {
                            ii = i + 1;
                        }
                    }

                    X(i, k) = sum;
                }

                for (int i = dimension - 1; i >= 0; i--)
                {
                    T sum = X(i, k);

                    for (int j = i + 1; j < dimension; j++)
                    {
                        sum -= _U(i, j) * X(j, k);
                    }

                    X(i, k) = sum /= _U(i, i);
                }
            }
        }
        else
        {
            if (Y.columnCount() != dimension)
            {
                return false;
            }

            X = Y;

            #pragma omp parallel for
            for (int k = 0; k < Y.rowCount(); k++)
            {
                int ii = 0;

                for (int i = 0; i < dimension; i++)
                {
                    int ip = P[i];
                    T sum    = X(k, ip);
                    X(k, ip) = X(k, i);

                    if (ii != 0)
                    {
                        for (int j = ii - 1; j < i; j++)
                        {
                            sum -= _U(i, j) * X(k, j);
                        }
                    }
                    else
                    {
                        if (sum != 0.0)
                        {
                            ii = i + 1;
                        }
                    }

                    X(k, i) = sum;
                }

                for (int i = dimension - 1; i >= 0; i--)
                {
                    T sum = X(k, i);

                    for (int j = i + 1; j < dimension; j++)
                    {
                        sum -= _U(i, j) * X(k, j);
                    }

                    X(k, i) = sum /= _U(i, i);
                }
            }
        }

        return true;
    }

    //(*@\Green{// Using multi-threading, tries to solve systems of linear equations of type $M \cdot X = B$, where $M$ is not necessarily}@*)
    //(*@\Green{// a real square matrix, while $B$ can be either a column or rectangular matrix of multiple columns with the same}@*)
    //(*@\Green{// row dimension as $M$. In case of success, the solution(s) will be stored in $X$.}@*)
    //(*@\Green{// If the input boolean parameter represent\_solutions\_as\_columns is false, then the column dimension of}@*)
    //(*@\Green{// the transpose of $B$ should be equal to the row dimension of $M$.}@*)
    //(*@\Green{// Note that, $T$ can be either double or Cartesian3, or any other custom type which has similar overloaded}@*)
    //(*@\Green{// mathematical and boolean operators.}@*)
    template <typename T>
    bool SVDecomposition::solveLinearSystem(
            const Matrix<T> &B, Matrix<T> &X, bool represent_solutions_as_columns)
    {
        if (!_decomposition_is_done)
        {
            X.resizeRows(0);
            return false;
        }

        int m = _U.rowCount();
        int n = _U.columnCount();

        double epsilon = std::numeric_limits<double>::epsilon();
        double threshold = 0.5 * sqrt(m + n + 1.0) * _S[0] * epsilon;

        if (represent_solutions_as_columns)
        {
            if (B.rowCount() != m)
            {
                throw Exception("SVDecomposition::solveLinearSystem : bad sizes!");
            }

            X.resizeRows(n);
            X.resizeColumns(B.columnCount());

            RowMatrix<double> tmp(n);

            #pragma omp parallel for
            for (int c = 0; c < B.columnCount(); c++)
            {
                for (int j = 0; j < n; j++)
                {
                    double s = 0.0;

                    if (_S[j] > threshold)
                    {
                        for (int i = 0; i < m; i++)
                        {
                            s += _U(i, j) * B(i, c);
                        }

                        s /= _S[j];
                    }

                    tmp[j] = s;
                }

                for (int j = 0; j < n; j++)
                {
                    double s = 0.0;

                    for (int jj = 0; jj < n; jj++)
                    {
                        s += _V(j, jj) * tmp[jj];
                    }

                    X(j, c) = s;
                }
            }
        }
        else
        {
            if (B.columnCount() != m)
            {
                throw Exception("SVDecomposition::solveLinearSystem : bad sizes!");
            }

            X.resizeRows(B.rowCount());
            X.resizeColumns(n);

            RowMatrix<double> tmp(n);

            #pragma omp parallel for
            for (int r = 0; r < B.rowCount(); r++)
            {
                for (int j = 0; j < n; j++)
                {
                    double s = 0.0;

                    if (_S[j] > threshold)
                    {
                        for (int i = 0; i < m; i++)
                        {
                            s += _U(i, j) * B(r, i);
                        }

                        s /= _S[j];
                    }

                    tmp[j] = s;
                }

                for (int j = 0; j < n; j++)
                {
                    double s = 0.0;

                    for (int jj = 0; jj < n; jj++)
                    {
                        s += _V(j, jj) * tmp[jj];
                    }

                    X(r, j) = s;
                }
            }
        }

        return true;
    }
}

#endif //(*@\Green{// REALMATRIXDECOMPOSITIONS\_H}@*)
