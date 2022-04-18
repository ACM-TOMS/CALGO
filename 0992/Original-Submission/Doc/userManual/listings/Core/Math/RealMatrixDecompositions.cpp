//---------------------------------------------------------------------------------------
// File:        Core/Math/RealMatrixDecompositions.cpp
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

#include "Constants.h"
#include "RealMatrixDecompositions.h"

#include "../Exceptions.h"

#include <algorithm>
#include <cmath>

using namespace std;

namespace cagd
{
    (*@\Green{// Tries to determine the pivoted Doolittle-type $LU$ decomposition of the given real square matrix $M$.}@*)
    (*@\Green{// If $M$ is singular then its $LU$ factorization cannot be determined.}@*)
    (*@\Green{// If the decomposition is unsuccessful, the boolean variable \_decomposition\_is\_done will be set to false, and if}@*)
    (*@\Green{// the input variable generate\_exceptions is set to true, then the method will also generate an exception with a}@*)
    (*@\Green{// meaningful reason.}@*)
    PLUDecomposition::PLUDecomposition(const RealMatrix &M,
                                       bool generate_exceptions):
        _decomposition_is_done(false),
        _determinant(0.0),
        _P(M.rowCount()),
        _LU(M)
    {
        if (!_LU.isSquare() && generate_exceptions)
        {
            throw Exception("Only square real matrices can have PLU decompositions!");
        }

        int dimension = _LU.rowCount();

        if (dimension <= 1 && generate_exceptions)
        {
            throw Exception("PLUDecomposition : the dimension of the given real "
                            "square matrix should be at least 2 x 2!");

            return;
        }

        for (int i = 0; i < dimension; i++)
        {
            _P[i] = i;
        }

        vector<double> implicit_scaling_of_each_row(dimension);

        double row_interchanges = 1.0;

        (*@\Green{// loop over the rows to get the implicit scaling information}@*)
        for (int r = 0; r < dimension; r++)
        {
            double big = 0.0;
            for (int c = 0; c < dimension; c++)
            {
                double temp = abs(_LU(r, c));
                if (temp > big)
                {
                    big = temp;
                }
            }

            if (big == 0.0)
            {
                _decomposition_is_done = false;
                _determinant           = 0.0;

                if (generate_exceptions)
                {
                    throw Exception("PLUDecomposition : the given square matrix is "
                                    "singular!");
                }

                return;
            }

            implicit_scaling_of_each_row[r] = 1.0 / big;
        }

        (*@\Green{// search for the largest pivot element}@*)
        for (int k = 0; k < dimension; k++)
        {
            int i_max = k;
            double big = 0.0;
            for (int i = k; i < dimension; i++)
            {
                double temp = implicit_scaling_of_each_row[i] * abs(_LU(i, k));

                if (temp > big)
                {
                    big   = temp;
                    i_max = i;
                }
            }

            (*@\Green{// do we need to interchange rows?}@*)
            if (k != i_max)
            {
                for (int j = 0; j < dimension; j++)
                {
                    double temp   = _LU(i_max, j);
                    _LU(i_max, j) = _LU(k, j);
                    _LU(k, j)     = temp;
                }

                (*@\Green{// change the parity of row\_interchanges}@*)
                row_interchanges = -row_interchanges;

                (*@\Green{// also interchange the scale factor}@*)
                implicit_scaling_of_each_row[i_max] = implicit_scaling_of_each_row[k];
            }

            _P[k] = i_max;

            if (_LU(k, k) == 0.0)
            {
                _LU(k, k) = TINY;
            }

            for (int i = k + 1; i < dimension; i++)
            {
                (*@\Green{// divide by pivot element}@*)
                double temp = _LU(i, k) /= _LU(k, k);

                (*@\Green{// reduce remaining submatrix}@*)
                for (int j = k + 1; j < dimension; j++)
                {
                    _LU(i, j) -= temp * _LU(k, j);
                }
            }
        }

        _determinant = row_interchanges;
        for (int i = 0; i < dimension; i++)
        {
            _determinant *= _LU(i, i);
        }

        _decomposition_is_done = true;
    }

    (*@\Green{// Returns either the successful or the unsuccessful state of the decomposition.}@*)
    bool PLUDecomposition::isCorrect() const
    {
        return _decomposition_is_done;
    }

    (*@\Green{// If the decomposition was successful the method will return the product of the diagonal entries of}@*)
    (*@\Green{// the upper triangular matrix $U$ (i.e., the determinant of $M$), otherwise it will return zero.}@*)
    double PLUDecomposition::determinant() const
    {
        if (_decomposition_is_done)
        {
            return _determinant;
        }

        return 0.0;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    PLUDecomposition* PLUDecomposition::clone() const
    {
        return new (nothrow) PLUDecomposition(*this);
    }

    (*@\Green{// The next constructor tries to determine the unpivoted Doolittle-type $LU$ decomposition of the given real}@*)
    (*@\Green{// square matrix $M$. If $M$ is singular then its $LU$ factorization cannot be determined.}@*)
    (*@\Green{// If the decomposition is unsuccessful, the boolean variable \_decomposition\_is\_done will be set to false,}@*)
    (*@\Green{// and if the input variable generate\_exceptions is set to true, then the method will}@*)
    (*@\Green{// also generate an exception with a meaningful reason.}@*)
    FactorizedUnpivotedLUDecomposition::FactorizedUnpivotedLUDecomposition(
            const RealMatrix &M, bool generate_exceptions):
        _decomposition_is_done(false),
        _determinant(0.0),
        _L(M.rowCount(), M.columnCount()), _U(M.rowCount(), M.columnCount())
    {
        if (!M.isSquare())
        {
            throw Exception("Only square real matrices can have unpivoted "
                            "LU decompositions!");
        }

        int dimension = M.rowCount();

        _L.loadIdentityMatrix();
        _U.loadNullMatrix();

        for (int j = 0; j < dimension; j++)
        {
            _U(0, j) = M(0, j);
        }

        for (int i = 1; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                for (int k = 0; k <= i - 1; k++)
                {
                    double s1 = 0.0;
                    if (k == 0)
                    {
                        s1 = 0.0;
                    }
                    else
                    {
                        for (int p = 0; p <= k - 1; p++)
                        {
                            s1 += _L(i, p) * _U(p, k);
                        }
                    }

                    if (!_U(k, k))
                    {
                        _L.resizeRows(0);
                        _U.resizeRows(0);

                        _decomposition_is_done = false;
                        _determinant = 0.0;

                        if (generate_exceptions)
                        {
                            throw Exception("FactorizedUnpivotedLUDecomposition : the "
                                            "given square matrix is singular!");
                        }

                        return;
                    }

                    _L(i, k) = (M(i, k) - s1) / _U(k, k);
                }

                for (int k=i; k < dimension; k++)
                {
                    double s2 = 0.0;

                    for (int p = 0; p <= i - 1; p++)
                    {
                        s2 += _L(i, p) * _U(p, k);
                    }

                    _U(i, k) = M(i, k) - s2;
                }
            }
        }

        _determinant = 1.0;
        for (int i = 0; i < dimension; i++)
        {
            _determinant *= _U(i, i);
        }

        _decomposition_is_done = true;
    }

    (*@\Green{// Returns either the successful or the unsuccessful state of the decomposition.}@*)
    bool FactorizedUnpivotedLUDecomposition::isCorrect() const
    {
        return _decomposition_is_done;
    }

    (*@\Green{// If the decomposition was successful the method will return the product of the diagonal entries of}@*)
    (*@\Green{// the upper triangular matrix $U$ (i.e., the determinant of $M$), otherwise it will return zero.}@*)
    double FactorizedUnpivotedLUDecomposition::determinant() const
    {
        if (_decomposition_is_done)
        {
            return _determinant;
        }

        return 0.0;
    }

    (*@\Green{// Returns a constant reference to the lower triangular matrix $L$.}@*)
    const RealMatrix& FactorizedUnpivotedLUDecomposition::L() const
    {
        return _L;
    }

    (*@\Green{// Returns a constant reference to the upper triangular matrix $U$.}@*)
    const RealMatrix& FactorizedUnpivotedLUDecomposition::U() const
    {
        return _U;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    FactorizedUnpivotedLUDecomposition* FactorizedUnpivotedLUDecomposition::clone() const
    {
        return new (nothrow) FactorizedUnpivotedLUDecomposition(*this);
    }

    (*@\Green{// Some auxiliar functions for singular value decomposition.}@*)
    template <typename T>
    inline T sign(const T &a, const T &b)
    {
        return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
    }

    double EuclideanNorm(const double &a, const double &b)
    {
        double abs_a  = abs(a), abs_b = abs(b),
               factor = abs_b / abs_a, rfactor = abs_a / abs_b;

        return (abs_a > abs_b ? abs_a * sqrt(1.0 + factor * factor) :
                                (abs_b == 0.0 ?
                                 0.0 : abs_b * sqrt(1.0 + rfactor * rfactor)));
    }

    (*@\Green{// Tries to determine the singular value decomposition of the given real matrix $M$.}@*)
    (*@\Green{// If the decomposition is unsuccessful, the boolean variable \_decomposition\_is\_done will be set to false, and if}@*)
    (*@\Green{// the input variable generate\_exceptions is set to true,}@*)
    (*@\Green{// then the method will also generate an exception with a meaningful reason.}@*)
    SVDecomposition::SVDecomposition(const RealMatrix &M, bool generate_exceptions):
        _decomposition_is_done(false),
        _product_of_singular_values(0.0),
        _U(M),
        _S(M.columnCount()),
        _V(M.columnCount(), M.columnCount())
    {
        int m = _U.rowCount();
        int n = _U.columnCount();

        (*@\Green{// decompose}@*)
        {
            bool              flag;
            int               i, iterations, j, jj, k, l, nm;
            double            M_norm = 0.0, c, f, g = 0.0, h, s, scale = 0.0, x, y, z;
            RowMatrix<double> rv1(n);

            for (i = 0; i < n; i++)
            {
                l = i + 2;
                rv1[i] = scale * g;
                g = s = scale = 0.0;

                if (i < m)
                {
                    for (k = i; k < m; k++)
                    {
                        scale += abs(_U(k, i));
                    }

                    if (scale != 0.0)
                    {
                        for (k = i; k < m; k++)
                        {
                            _U(k, i) /= scale;
                            s += _U(k, i) * _U(k, i);
                        }

                        f = _U(i, i);
                        g = -sign(sqrt(s), f);
                        h = f * g - s;

                        _U(i, i) = f - g;

                        for (j = l - 1; j < n; j++)
                        {
                            s = 0.0;
                            for (k = i; k < m; k++)
                            {
                                s += _U(k, i) *_U(k, j);
                            }

                            f = s / h;
                            for (k = i; k < m; k++)
                            {
                                _U(k, j) += f * _U(k, i);
                            }
                        }

                        for (k = i; k < m; k++)
                        {
                            _U(k, i) *= scale;
                        }
                    }
                }

                _S[i] = scale * g;
                g = s = scale = 0.0;

                if (i+1 <= m && i+1 != n)
                {
                    for (k = l - 1; k < n; k++)
                    {
                        scale += abs(_U(i, k));
                    }

                    if (scale != 0.0)
                    {
                        for (k = l - 1; k < n; k++)
                        {
                            _U(i, k) /= scale;
                            s += _U(i, k) * _U(i, k);
                        }

                        f = _U(i, l - 1);
                        g = -sign(sqrt(s),f);
                        h = f * g - s;

                        _U(i, l - 1) = f - g;

                        for (k = l - 1; k < n; k++)
                        {
                            rv1[k] = _U(i, k) / h;
                        }

                        for (j = l - 1; j < m; j++)
                        {
                            s = 0.0;
                            for (k = l - 1; k < n; k++)
                            {
                                s += _U(j, k) * _U(i, k);
                            }

                            for (k = l - 1; k < n; k++)
                            {
                                _U(j, k) += s * rv1[k];
                            }
                        }

                        for (k = l - 1; k < n; k++)
                        {
                            _U(i, k) *= scale;
                        }
                    }
                }

                M_norm = max(M_norm, (abs(_S[i]) + abs(rv1[i])));
            }

            for (i = n - 1; i >= 0; i--)
            {
                if (i < n - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                        {
                            _V(j, i) = (_U(i, j) / _U(i, l)) / g;
                        }

                        for (j = l; j < n; j++)
                        {
                            s = 0.0;
                            for (k = l; k < n; k++)
                            {
                                s += _U(i, k) * _V(k, j);
                            }

                            for (k = l; k < n; k++)
                            {
                                _V(k, j) += s * _V(k, i);
                            }
                        }
                    }

                    for (j = l; j < n; j++)
                    {
                        _V(i, j) = _V(j, i) = 0.0;
                    }
                }

                _V(i, i) = 1.0;
                g = rv1[i];
                l = i;
            }

            for (i = min(m, n) - 1; i >= 0; i--)
            {
                l = i+1;
                g = _S[i];
                for (j = l; j < n; j++)
                {
                    _U(i, j) = 0.0;
                }

                if (g != 0.0)
                {
                    g = 1.0 / g;
                    for (j = l; j < n; j++)
                    {
                        s = 0.0;
                        for (k = l; k < m; k++)
                        {
                            s += _U(k, i) * _U(k, j);
                        }

                        f = (s / _U(i, i)) * g;
                        for (k = i; k < m; k++)
                        {
                            _U(k, j) += f * _U(k, i);
                        }
                    }

                    for (j = i; j < m; j++)
                    {
                        _U(j, i) *= g;
                    }
                }
                else
                {
                    for (j = i; j < m; j++)
                    {
                        _U(j, i) = 0.0;
                    }
                }

                ++_U(i, i);
            }

            for (k = n - 1; k >= 0; k--)
            {
                for (iterations = 0; iterations < 30; iterations++)
                {
                    flag = true;
                    for (l = k; l >= 0; l--)
                    {
                        nm = l - 1;

                        if (l == 0 || abs(rv1[l]) <= MACHINE_EPS * M_norm)
                        {
                            flag = false;
                            break;
                        }

                        if (abs(_S[nm]) <= MACHINE_EPS * M_norm)
                        {
                            break;
                        }
                    }

                    if (flag)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i < k + 1; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];

                            if (abs(f) <= MACHINE_EPS * M_norm)
                            {
                                break;
                            }

                            g     = _S[i];
                            h     = EuclideanNorm(f, g);
                            _S[i] = h;
                            h     = 1.0 / h;
                            c     = g * h;
                            s     = -f * h;

                            for (j = 0; j < m; j++)
                            {
                                y         = _U(j, nm);
                                z         = _U(j, i);
                                _U(j, nm) = y * c + z * s;
                                _U(j, i)  = z * c - y * s;
                            }
                        }
                    }

                    z = _S[k];

                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            _S[k] = -z;
                            for (j = 0; j < n; j++)
                            {
                                _V(j, k) = -_V(j, k);
                            }
                        }
                        break;
                    }

                    if (iterations == 29)
                    {
                        _decomposition_is_done = false;

                        if (generate_exceptions)
                        {
                            throw Exception("No convergence in 30 singular value "
                                            "decomposition iterations!");
                        }

                        return;
                    }

                    x  = _S[l];
                    nm = k - 1;
                    y  = _S[nm];
                    g  = rv1[nm];
                    h  = rv1[k];
                    f  = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g  = EuclideanNorm(f, 1.0);
                    f  = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
                    c  = s = 1.0;

                    for (j = l; j <= nm; j++)
                    {
                        i      =  j + 1;
                        g      =  rv1[i];
                        y      =  _S[i];
                        h      =  s * g;
                        g      =  c * g;
                        z      =  EuclideanNorm(f, h);
                        rv1[j] =  z;
                        c      =  f / z;
                        s      =  h / z;
                        f      =  x * c + g * s;
                        g      =  g * c - x * s;
                        h      =  y * s;
                        y      *= c;

                        for (jj = 0; jj < n; jj++)
                        {
                            x        = _V(jj, j);
                            z        = _V(jj, i);
                            _V(jj, j) = x * c + z * s;
                            _V(jj, i) = z * c - x * s;
                        }

                        z     = EuclideanNorm(f, h);
                        _S[j] = z;

                        if (z)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }

                        f = c * g + s * y;
                        x = c * y - s * g;

                        for (jj = 0; jj < m; jj++)
                        {
                            y         = _U(jj, j);
                            z         = _U(jj, i);
                            _U(jj, j) = y * c + z * s;
                            _U(jj, i) = z * c - y * s;
                        }
                    }

                    rv1[l] = 0.0;
                    rv1[k] = f;
                    _S[k]  = x;
                }
            }
        }

        (*@\Green{// reorder}@*)
        {
            int               i, j, k, s, inc = 1;
            double            sw;
            RowMatrix<double> su(m), sv(n);

            do
            {
                inc *= 3;
                inc++;
            }
            while (inc <= n);

            do
            {
                inc /= 3;
                for (i = inc; i < n; i++)
                {
                    sw = _S[i];
                    for (k = 0; k < m; k++)
                    {
                        su[k] = _U(k, i);
                    }

                    for (k = 0; k < n; k++)
                    {
                        sv[k] = _V(k, i);
                    }

                    j = i;

                    while (_S[j - inc] < sw)
                    {
                        _S[j] = _S[j - inc];

                        for (k=0;k<m;k++)
                        {
                            _U(k, j) = _U(k, j - inc);
                        }

                        for (k = 0; k < n; k++)
                        {
                            _V(k, j) = _V(k, j - inc);
                        }

                        j -= inc;

                        if (j < inc)
                        {
                            break;
                        }
                    }

                    _S[j] = sw;

                    for (k = 0; k < m; k++)
                    {
                        _U(k, j) = su[k];
                    }

                    for (k = 0; k < n; k++)
                    {
                        _V(k, j) = sv[k];
                    }
                }
            }
            while (inc > 1);

            for (k = 0; k < n; k++)
            {
                s=0;

                for (i = 0; i < m; i++)
                {
                    if (_U(i, k) < 0.0)
                    {
                        s++;
                    }
                }

                for (j = 0; j < n; j++)
                {
                    if (_V(j, k) < 0.0)
                    {
                        s++;
                    }
                }

                if (s > (m + n) / 2)
                {
                    for (i = 0; i < m; i++)
                    {
                        _U(i, k) = -_U(i, k);
                    }

                    for (j = 0; j < n; j++)
                    {
                        _V(j, k) = -_V(j, k);
                    }
                }
            }
        }

        (*@\Green{// product of singular values (i.e., the determinant of M provided that it is a square matrix)}@*)
        _product_of_singular_values = 1.0;
        for (int i = 0; i < _S.columnCount(); i++)
        {
            _product_of_singular_values *= _S[i];
        }

        _decomposition_is_done = true;
    }

    (*@\Green{// Returns either the successful or the unsuccessful state of the decomposition.}@*)
    bool SVDecomposition::isCorrect() const
    {
        return _decomposition_is_done;
    }

    (*@\Green{// Returns the ratio of the largest and smalles singular values.}@*)
    double SVDecomposition::conditionNumber() const
    {
        int last = _S.columnCount() - 1;
        return (_S[0] <= 0.0 || _S[last] <= 0.0) ? 0.0 : _S[0] / _S[last];
    }

    (*@\Green{// Returns the ratio of the smallest and largest singular values.}@*)
    double SVDecomposition::reciprocalConditionNumber() const
    {
        int last = _S.columnCount() - 1;
        return (_S[0] <= 0.0 || _S[last] <= 0.0) ? 0.0 : _S[last] / _S[0];
    }

    (*@\Green{// Returns the product of the obtained singular values. If $M$ is a regular real square matrix,}@*)
    (*@\Green{// then this product coincides with the determinant of $M$.}@*)
    double SVDecomposition::productOfSingularValues() const
    {
        if (_decomposition_is_done)
        {
            return _product_of_singular_values;
        }

        return 0.0;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    SVDecomposition* SVDecomposition::clone() const
    {
        return new (nothrow) SVDecomposition(*this);
    }
}
