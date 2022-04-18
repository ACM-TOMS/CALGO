//----------------------------------------------------------------------------------
// File:        EC/BSurfaces3.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "BSurfaces3.h"

#include "../Core/Exceptions.h"
#include "BCurves3.h"

#include <algorithm>
#include <cmath>
#include <new>

using namespace std;

namespace cagd
{
    //(*@\Green{// special constructor}@*)
    OrdinarySurfaceCoefficients::OrdinarySurfaceCoefficients(
            GLint u_dimension, GLint v_dimension, const std::vector<GLint> &sigma)
    {
        if (sigma.size() != 3)
        {
            throw Exception("The vector sigma should contain 3 elements!");
        }

        for (GLint k = 0; k < 3; k++)
        {
            if (sigma[k] <= 0)
            {
                throw Exception("Each element of the vector should be positive "
                                "integers!");
            }
        }

        if (!u_dimension || !v_dimension)
        {
            throw Exception("Dimensions of underlying EC spaces should be at least 1 in "
                            "all directions!");
        }

        _lambda.resize(3);

        for (GLint k = 0; k < 3; k++)
        {
            for (GLint zeta = 0; zeta < sigma[k]; zeta++)
            {
                _lambda[k].push_back(pair< RowMatrix<GLdouble>, RowMatrix<GLdouble> >(
                    RowMatrix<GLdouble>(u_dimension), RowMatrix<GLdouble>(v_dimension)));
            }
        }
    }

    //(*@\Green{// based on the given variable type, returns the constant reference of either $\lambda_{n_0,i_0}^{k,\zeta}$ or $\lambda_{n_1,i_1}^{k,\zeta}$, where the}@*)
    //(*@\Green{// input variable index denotes either $i_0$ or $i_1$}@*)
    const GLdouble& OrdinarySurfaceCoefficients::operator()(
            GLint k, GLint zeta, variable::Type type, GLint index) const
    {
        assert("The given coordinate function index is out of bounds!" &&
               (k >= 0 && k < (GLint)_lambda.size()));

        assert("The given summation term index is out of bounds!" &&
               (zeta >= 0 && zeta < (GLint)_lambda[k].size()));

        assert("The given coefficient index is out of bounds!" &&
               (index >= 0 && (type == variable:: U ?
                               index < _lambda[k][zeta].first.columnCount() :
                               index < _lambda[k][zeta].second.columnCount())));

        return (type == variable::U ? _lambda[k][zeta].first[index] :
                                      _lambda[k][zeta].second[index]);
    }

    //(*@\Green{// based on the given variable type, returns non-constant references to either $\lambda_{n_0,i_0}^{k,\zeta}$ or $\lambda_{n_1,i_1}^{k,\zeta}$, where the}@*)
    //(*@\Green{// input variable index denotes either $i_0$ or $i_1$}@*)
    GLdouble& OrdinarySurfaceCoefficients::operator()(
            GLint k, GLint zeta, variable::Type type, GLint index)
    {
        assert("The given coordinate function index is out of bounds!" &&
               (k >= 0 && k < (GLint)_lambda.size()));

        assert("The given summation term index is out of bounds!" &&
               (zeta >= 0 && zeta < (GLint)_lambda[k].size()));

        assert("The given coefficient index is out of bounds!" &&
               (index >= 0 && (type == variable:: U ?
                               index < _lambda[k][zeta].first.columnCount() :
                               index < _lambda[k][zeta].second.columnCount())));

        return (type == variable::U ? _lambda[k][zeta].first[index] :
                                      _lambda[k][zeta].second[index]);
    }

    //(*@\Green{// returns the number of multiplied linear combination pairs appearing in the coordinate function $s^{k}$}@*)
    GLint OrdinarySurfaceCoefficients::sigma(GLint k) const
    {
        assert("The given coordinate function index is out of bounds!" &&
               (k >= 0 && k < (GLint)_lambda.size()));

        return (GLint)_lambda[k].size();
    }

    //(*@\Green{// based on the given variable type, returns the value either of $n_0+1$ or $n_1+1$}@*)
    GLint OrdinarySurfaceCoefficients::dimension(variable::Type type) const
    {
        return (type == variable::U ? _lambda[0][0].first.columnCount() :
                                      _lambda[0][0].second.columnCount());
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    OrdinarySurfaceCoefficients* OrdinarySurfaceCoefficients::clone() const
    {
        return new (nothrow) OrdinarySurfaceCoefficients(*this);
    }

    //(*@\Green{// special constructor}@*)
    BSurface3::BSurface3(const ECSpace &u_space, const ECSpace &v_space):
        TensorProductSurface3(
                u_space.alpha(), u_space.beta(),
                v_space.alpha(), v_space.beta(),
                u_space.dimension(), v_space.dimension())
    {
        _S[variable::U] = u_space;
        _S[variable::V] = v_space;

        for (GLint i = variable::U; i <= variable::V; i++)
        {
            _T[i] = SP<RealMatrix>::Default(
                    _S[i].basisTransformationFromNBToOrdinary());

            if (!_T[i])
            {
                throw Exception("One of the transformation matrices that maps the "
                                "normalized B-basis of an EC space to its ordinary "
                                "basis could not be created!");
            }
        }
    }

    //(*@\Green{// Inherited virtual method that has to be redeclared and redefined, since if one alters the endpoints of}@*)
    //(*@\Green{// the definition domain $\left[\alpha_r, \beta_r\right]$ both bases $\mathcal{F}_{n_r}^{\alpha_r,\beta_r}$ and $\mathcal{B}_{n_r}^{\alpha_r,\beta_r}$ of $\mathbb{S}_{n_r}^{\alpha_r,\beta_r}$ have to be updated. Note that,}@*)
    //(*@\Green{// the new interval length should also satisfy the condition $0 < \beta_r - \alpha_r < \ell^{\prime}\left(\mathbb{S}_{n_r}^{\alpha_r,\beta_r}\right)$. The corresponding}@*)
    //(*@\Green{// basis transformation matrix $[t_{i_r,j_r}^{n_r}]_{i_r=0,\,j_r=0}^{n_r,\,n_r}$ that maps $\mathcal{B}_{n_r}^{\alpha_r, \beta_r}$ to $\mathcal{F}_{n_r}^{\alpha_r, \beta_r}$ will also be updated.}@*)
    GLboolean BSurface3::setInterval(
            variable::Type type,
            GLdouble alpha, GLdouble beta,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits)
    {
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);

        if (expected_correct_significant_digits < 0 ||
            !TensorProductSurface3::setInterval(type, alpha, beta) ||
            !_S[type].setDefinitionDomain(alpha, beta,
                                          check_for_ill_conditioned_matrices,
                                          expected_correct_significant_digits))
        {
            return GL_FALSE;
        }

        _T[type] = SP<RealMatrix>::Default(
                _S[type].basisTransformationFromNBToOrdinary());

        return (_T[type] ? GL_TRUE : GL_FALSE);
    }

    //(*@\Green{// The next tree member functions are inherited pure virtual methods that have to be redeclared and defined.}@*)

    //(*@\Green{// Evaluates all normalized B-basis function values $\left\{b_{n_r,i_r}\left(u_r\right)\right\}_{i_r = 0}^{n_r}$, where $u_r \in \left[\alpha_r, \beta_r\right]$ is specified by the}@*)
    //(*@\Green{// variable 'parameter\_value' and, based on the variable 'type', $r$ is either $0$, or $1$. In case of success, the obtained}@*)
    //(*@\Green{// values will be stored in the real row matrix referenced by the variable 'values'.}@*)
    GLboolean BSurface3::blendingFunctionValues(
            variable::Type type, GLdouble parameter_value,
            RowMatrix<GLdouble> &values) const
    {
        if (parameter_value < _S[type].alpha() || parameter_value > _S[type].beta())
        {
            values.resizeColumns(0);
            return GL_FALSE;
        }

        GLint dimension = _S[type].dimension();

        values.resizeColumns(dimension);

        #pragma omp parallel for
        for (GLint i = 0; i < dimension; i++)
        {
            values[i] = _S[type](ECSpace::B_BASIS, i, 0, parameter_value);
        }

        return GL_TRUE;
    }

    //(*@\Green{// Calculates the partial derivatives $\left[\frac{\partial^{r}}{\partial u_0^{r-j} \partial u_1^{j}}\mathbf{s}_{n_0,n_1}\left(u_0, u_1\right)\right]_{j=0}^{r}$ for all $r=0,\ldots,\rho$, where $\rho\in\mathbb{N}$ denotes}@*)
    //(*@\Green{// the maximum order of partial derivatives. Parameter values $u_0\in\left[\alpha_0,\beta_0\right]$ and $u_1\in\left[\alpha_1,\beta_1\right]$ correspond to}@*)
    //(*@\Green{// the input variables 'u' and 'v', respectively, while, in case of success, the obtained partial derivatives will}@*)
    //(*@\Green{// be stored in a triangular matrix referenced by the variable 'd'.}@*)
    GLboolean BSurface3::calculateAllPartialDerivatives(
            GLint maximum_order_of_partial_derivatives,
            GLdouble u, GLdouble v, PartialDerivatives& d) const
    {
        assert("The given maximum order of partial derivatives should be "
               "non-negative!" && maximum_order_of_partial_derivatives >= 0);

        if (maximum_order_of_partial_derivatives < 0 ||
            u < _S[variable::U].alpha() || u > _S[variable::U].beta() ||
            v < _S[variable::V].alpha() || v > _S[variable::V].beta())
        {
            d.resizeRows(0);
            return GL_FALSE;
        }

        d.resizeRows(maximum_order_of_partial_derivatives + 1);
        d.loadNullVectors();

        GLint u_dimension = _S[variable::U].dimension();
        GLint v_dimension = _S[variable::V].dimension();

        #pragma omp parallel for
        for (GLint order = 0; order <= maximum_order_of_partial_derivatives; order++)
        {
            #pragma omp parallel for
            for (GLint j_u = order; j_u >= 0; j_u--)
            {
                GLint j_v = order - j_u;

                Cartesian3 &derivative = d(order, j_v);

                #pragma omp parallel for
                for (GLint i_u = 0; i_u < u_dimension; i_u++)
                {
                    Cartesian3 aux;

                    for (GLint i_v = 0; i_v < v_dimension; i_v++)
                    {
                        aux += _data(i_u, i_v) *
                               _S[variable::V](ECSpace::B_BASIS, i_v, j_v, v);
                    }

                    derivative += aux * _S[variable::U](ECSpace::B_BASIS, i_u, j_u, u);
                }
            }
        }

        return GL_TRUE;
    }

    //(*@\Green{// Computes the partial derivatives $\left[\frac{\partial^{j}}{\partial u_r^{j}}\mathbf{s}_{n_0,n_1}\left(u_0, u_1\right)\right]_{j=0}^{\rho}$, where $\rho\in\mathbb{N}$ denotes the maximum order of}@*)
    //(*@\Green{// partial derivatives. Parameter values $u_0\in\left[\alpha_0,\beta_0\right]$ and $u_1\in\left[\alpha_1,\beta_1\right]$ correspond to the input variables 'u'}@*)
    //(*@\Green{// and 'v', respectively. The index $r$ is determined by the variable 'direction'. In case of success, the obtained}@*)
    //(*@\Green{// partial derivatives will be stored in a column matrix referenced by the variable 'd'.}@*)
    GLboolean BSurface3::calculateDirectionalDerivatives(
            variable::Type direction, GLint maximum_order_of_directional_derivatives,
            GLdouble u, GLdouble v, DirectionalDerivatives &d) const
    {
        assert("The given maximum order of directional derivatives should be "
               "non-negative!" && maximum_order_of_directional_derivatives >= 0);

        if (maximum_order_of_directional_derivatives < 0)
        {
            d.resizeRows(0);
            return GL_FALSE;
        }

        if (direction == variable::U)
        {
            if (u < _S[variable::U].alpha() || u > _S[variable::U].beta())
            {
                d.resizeRows(0);
                return GL_FALSE;
            }

            d.resizeRows(maximum_order_of_directional_derivatives + 1);
            d.loadNullVectors();

            GLint u_dimension = _S[variable::U].dimension();
            GLint v_dimension = _S[variable::V].dimension();

            #pragma omp parallel for
            for (GLint order = 0; order <= maximum_order_of_directional_derivatives;
                 order++)
            {
                Cartesian3 &u_derivative = d[order];

                #pragma omp parallel for
                for (GLint i = 0; i < u_dimension; i++)
                {
                    Cartesian3 aux;

                    for (GLint j = 0; j < v_dimension; j++)
                    {
                        aux += _data(i, j) * _S[variable::V](ECSpace::B_BASIS, j, 0, v);
                    }

                    u_derivative += aux * _S[variable::U](ECSpace::B_BASIS, i, order, u);
                }
            }

            return GL_TRUE;
        }
        else
        {
            if (v < _S[variable::V].alpha() || v > _S[variable::U].beta())
            {
                d.resizeRows(0);
                return GL_FALSE;
            }

            d.resizeRows(maximum_order_of_directional_derivatives + 1);
            d.loadNullVectors();

            GLint u_dimension = _S[variable::U].dimension();
            GLint v_dimension = _S[variable::V].dimension();

            #pragma omp parallel for
            for (GLint order = 0; order <= maximum_order_of_directional_derivatives;
                 order++)
            {
                Cartesian3 &v_derivative = d[order];

                #pragma omp parallel for
                for (GLint i = 0; i < u_dimension; i++)
                {
                    Cartesian3 aux;

                    for (GLint j = 0; j < v_dimension; j++)
                    {
                        aux += _data(i, j) *
                               _S[variable::V](ECSpace::B_BASIS, j, order, v);
                    }

                    v_derivative += aux * _S[variable::U](ECSpace::B_BASIS, i, 0, u);
                }
            }

            return GL_TRUE;
        }
    }

    //(*@\Green{// Overloaded function operator that evaluates the derivative $b_{n_r,i_r}^{\left(j\right)}\left(u_r\right)$, where $u_r \in \left[\alpha_r, \beta_r\right]$, $i_r=0,1,\ldots,n_r$,}@*)
    //(*@\Green{// and $j \in \mathbb{N}$.}@*)
    GLdouble BSurface3::operator ()(
            variable::Type type, GLint function_index,
            GLint differentiation_order, GLdouble parameter_value) const
    {
        assert("The given B-basis function index is out of bounds!" &&
               (function_index >= 0 && function_index < _S[type].dimension()));
        assert("The given differentiation order should be non-negative!" &&
               differentiation_order >= 0);

        return _S[type](ECSpace::B_BASIS, function_index,
                        differentiation_order, parameter_value);
    }

    //(*@\Green{// Overloaded member functions for general order elevation. In order to generate a higher dimensional EC}@*)
    //(*@\Green{// space $\mathbb{S}_{n_r+q}^{\alpha_r,\beta_r}$ that fulfills the condition $\mathbb{S}_{n_r}^{\alpha_r,\beta_r}\subset\mathbb{S}_{n_r+q}^{\alpha_r,\beta_r}$ where $q\geq 1$ and $0<\beta_r-\alpha_r < \min\left\{\ell^{\prime}\left(\mathbb{S}_{n_r}^{\alpha_r,\beta_r}\right),\right.$}@*)
    //(*@\Green{// $\left.\ell^{\prime}\left(\mathbb{S}_{n_r+q}^{\alpha_r,\beta_r}\right)\right\}$ one has either to define a new complex root $z=a\pm\mathbf{i}b$ of multiplicity $m \geq 1$, or to specify}@*)
    //(*@\Green{// an existing root with a multiplicity that is greater than its former value $m'\geq 1$. If $z \in \mathbb{C}\setminus\mathbb{R}$, then}@*)
    //(*@\Green{// $q=2\left(m-m'\right)$, otherwise $q=m-m'$, where $m'=0$ whenever $z$ denotes a nonexisting former root.}@*)
    //(*@\Green{// Their implementations are based on the extension of Lemma \mref{lem:general_order_elevation}.}@*)
    BSurface3* BSurface3::performOrderElevation( //(*@\label{BSurface3::performOrderElevation::start}@*)
            variable::Type type,
            GLdouble a, GLdouble b, GLint multiplicity,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits) const
    {
        assert("The given order of multiplicity should be non-negative!" &&
               multiplicity >= 0);
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);


        return performOrderElevation(type, 
                                     CharacteristicPolynomial::Zero(a, b, multiplicity),
                                     check_for_ill_conditioned_matrices,
                                     expected_correct_significant_digits);
    }

    BSurface3* BSurface3::performOrderElevation(
            variable::Type type,
            const CharacteristicPolynomial::Zero zero,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits) const
    {
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);

        if (expected_correct_significant_digits < 0)
        {
            return nullptr;
        }

        BSurface3* result = clone();

        if (!result)
        {
            return nullptr;
        }

        if (type == variable::U)
        {
            try
            {
                if (!result->_S[variable::U].insertZero(
                        zero,
                        true,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits))
                {
                    delete result, result = nullptr;
                    return nullptr;
                }
            }
            catch (Exception &e)
            {
                delete result, result = nullptr;
                throw e;
            }

            GLint u_old_dimension      = _S[variable::U].dimension();
            GLint u_elevated_dimension = result->_S[variable::U].dimension();

            if (u_old_dimension == u_elevated_dimension)
            {
                return result;
            }

            GLdouble u_alpha = _S[variable::U].alpha(), u_beta = _S[variable::U].beta();

            result->_data.resizeRows(u_elevated_dimension);

            if (u_elevated_dimension - u_old_dimension == 1)
            {
                GLint n = u_old_dimension - 1;

                #pragma omp parallel for
                for (GLint column = 0; column < _S[variable::V].dimension(); column++)
                {
                    (*result)(0, column)               = _data(0, column);
                    (*result)(u_old_dimension, column) = _data(n, column);

                    #pragma omp parallel for
                    for (GLint i = 1; i <= n / 2; i++)
                    {
                        GLdouble ratio       = (*this)(variable::U, i, i, u_alpha) /
                                               (*result)(variable::U, i, i, u_alpha);
                        (*result)(i, column) = (1 - ratio) * _data(i - 1, column) +
                                               ratio * _data(i, column);
                    }

                    #pragma omp parallel for
                    for (GLint i = 1; i <= u_old_dimension / 2; i++)
                    {
                        GLdouble ratio                         =
                                (*this)(variable::U, n - i, i, u_beta) /
                                (*result)(variable::U, u_old_dimension - i, i, u_beta);
                        (*result)(u_old_dimension - i, column) =
                                ratio * _data(n - i, column) +
                                (1 - ratio) * _data(u_old_dimension - i, column);
                    }
                }

                return result;
            }

            GLint v_dimension = _S[variable::V].dimension();
            GLdouble v_alpha  = _S[variable::V].alpha(), v_beta = _S[variable::V].beta();

            RowMatrix<GLdouble>     u_knot_vector(u_elevated_dimension);
            ColumnMatrix<GLdouble>  v_knot_vector(v_dimension);

            Matrix<Cartesian3>      sample(u_elevated_dimension, v_dimension);

            GLdouble u_step = (u_beta - u_alpha) / (u_elevated_dimension - 1);
            GLdouble v_step = (v_beta - v_alpha) / (v_dimension - 1);

            GLboolean sampling_aborted = GL_FALSE;

            #pragma omp parallel for
            for (GLint i_j = 0; i_j < u_elevated_dimension * v_dimension; i_j++)
            {
                #pragma omp flush (sampling_aborted)
                if (!sampling_aborted)
                {
                    GLint i = i_j / v_dimension;
                    GLint j = i_j % v_dimension;

                    u_knot_vector[i] = min(u_alpha + i * u_step, u_beta);
                    v_knot_vector[j] = min(v_alpha + j * v_step, v_beta);

                    PartialDerivatives d;
                    if (!calculateAllPartialDerivatives(
                            0, u_knot_vector[i], v_knot_vector[j], d))
                    {
                        sampling_aborted = GL_TRUE;
                        #pragma omp flush (sampling_aborted)
                    }

                    sample(i, j) = d(0, 0);
                }
            }

            if (sampling_aborted ||
                !result->updateDataForInterpolation(
                        u_knot_vector, v_knot_vector, sample))
            {
                delete result, result = nullptr;
            }

            return result;
        }
        else
        {
            try
            {
                if (!result->_S[variable::V].insertZero(
                        zero,
                        true,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits))
                {
                    delete result, result = nullptr;
                    return nullptr;
                }
            }
            catch(Exception &e)
            {
                delete result, result = nullptr;
                throw e;
            }

            GLint v_old_dimension      = _S[variable::V].dimension();
            GLint v_elevated_dimension = result->_S[variable::V].dimension();

            if (v_old_dimension == v_elevated_dimension)
            {
                return result;
            }

            GLdouble v_alpha = _S[variable::V].alpha(), v_beta = _S[variable::V].beta();

            result->_data.resizeColumns(v_elevated_dimension);

            if (v_elevated_dimension - v_old_dimension == 1)
            {
                GLint n = v_old_dimension - 1;

                #pragma omp parallel for
                for (GLint row = 0; row < _S[variable::U].dimension(); row++)
                {
                    (*result)(row, 0)               = _data(row, 0);
                    (*result)(row, v_old_dimension) = _data(row, n);

                    #pragma omp parallel for
                    for (GLint i = 1; i <= n / 2; i++)
                    {
                        GLdouble ratio    = (*this)(variable::V, i, i, v_alpha) /
                                            (*result)(variable::V, i, i, v_alpha);
                        (*result)(row, i) = (1 - ratio) * _data(row, i - 1) +
                                            ratio * _data(row, i);
                    }

                    #pragma omp parallel for
                    for (GLint i = 1; i <= v_old_dimension / 2; i++)
                    {
                        GLdouble ratio                      =
                                (*this)(variable::V, n - i, i, v_beta) /
                                (*result)(variable::V, v_old_dimension - i, i, v_beta);
                        (*result)(row, v_old_dimension - i) =
                                ratio * _data(row, n - i) +
                                (1 - ratio) * _data(row, v_old_dimension - i);
                    }
                }

                return result;
            }

            GLint u_dimension  = _S[variable::U].dimension();
            GLdouble u_alpha  = _S[variable::U].alpha(), u_beta = _S[variable::U].beta();

            RowMatrix<GLdouble>     u_knot_vector(u_dimension);
            ColumnMatrix<GLdouble>  v_knot_vector(v_elevated_dimension);

            Matrix<Cartesian3>      sample(u_dimension, v_elevated_dimension);

            GLdouble u_step = (u_beta - u_alpha) / (u_dimension - 1);
            GLdouble v_step = (v_beta - v_alpha) / (v_elevated_dimension - 1);

            GLboolean sampling_aborted = GL_FALSE;

            #pragma omp parallel for
            for (GLint ij = 0; ij < u_dimension * v_elevated_dimension; ij++)
            {
                #pragma omp flush (sampling_aborted)
                if (!sampling_aborted)
                {
                    GLint i = ij / v_elevated_dimension;
                    GLint j = ij % v_elevated_dimension;

                    u_knot_vector[i] = min(u_alpha + i * u_step, u_beta);
                    v_knot_vector[j] = min(v_alpha + j * v_step, v_beta);

                    PartialDerivatives d;
                    if (!calculateAllPartialDerivatives(
                            0, u_knot_vector[i], v_knot_vector[j], d))
                    {
                        sampling_aborted = GL_TRUE;
                        #pragma omp flush (sampling_aborted)
                    }
                    else
                    {
                        sample(i, j) = d(0, 0);
                    }
                }
            }

            if (sampling_aborted ||
                !result->updateDataForInterpolation(
                        u_knot_vector, v_knot_vector, sample))
            {
                delete result, result = nullptr;
            }

            return result;
        }
    } //(*@\label{BSurface3::performOrderElevation::end}@*)

    //(*@\Green{// Using the extension of the general B-algorithm formulated in Theorem \mref{thm:general_subdivision}, the method subdivides the}@*)
    //(*@\Green{// given B-surface along the specified direction into two patches of the same order at the parameter value}@*)
    //(*@\Green{// $\gamma \in \left(\alpha_r,\beta_r\right)$. In case of success, the output is a nonzero raw pointer to a $2$-element row matrix that stores}@*)
    //(*@\Green{// $2$ smart pointers (based on the deep copy ownership policy) that point to the left and right patches of the}@*)
    //(*@\Green{// subdivided surface.}@*)
    RowMatrix<SP<BSurface3>::Default>* BSurface3::performSubdivision( //(*@\label{BSurface3::performSubdivision::start}@*)
            variable::Type type,
            GLdouble gamma,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits) const
    {
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);

        if (expected_correct_significant_digits < 0)
        {
            return nullptr;
        }

        if (type == variable::U)
        {
            GLdouble u_alpha = _S[variable::U].alpha(), u_beta = _S[variable::U].beta();

            if (gamma <= u_alpha || gamma >= u_beta)
            {
                return nullptr;
            }

            RowMatrix<SP<BSurface3>::Default>* result =
                    new (nothrow) RowMatrix<SP<BSurface3>::Default>(2);

            if (!result)
            {
                return nullptr;
            }

            ECSpace lambda_space(_S[variable::U]), rho_space(_S[variable::U]);

            try
            {
                if (!lambda_space.setDefinitionDomain(
                        u_alpha, gamma,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits) ||
                    !rho_space.setDefinitionDomain(
                            gamma, u_beta,
                            check_for_ill_conditioned_matrices,
                            expected_correct_significant_digits))
                {
                    delete result, result = nullptr;
                    return nullptr;
                }
            }
            catch (Exception &e)
            {
                delete result, result = nullptr;
                throw e;
            }

            (*result)[0] = SP<BSurface3>::Default(
                    new (nothrow) BSurface3(lambda_space, _S[variable::V]));

            (*result)[1] = SP<BSurface3>::Default(
                    new (nothrow) BSurface3(rho_space, _S[variable::V]));

            if (!(*result)[0] || !(*result)[1])
            {
                delete result, result = nullptr;
                return nullptr;
            }

            SP<BSurface3>::Default &L = (*result)[0];
            SP<BSurface3>::Default &R = (*result)[1];

            GLint u_dimension = _S[variable::U].dimension();
            GLint u_n         = u_dimension - 1;
            GLint half_u_n    = u_n / 2;

            GLint  v_dimension = _S[variable::V].dimension();
            GLint v_n         = v_dimension - 1;

            for (GLint l = 0; l <= v_n; l++)
            {
                ColumnMatrix<Cartesian3> c_alpha(half_u_n + 1),
                                         c_gamma(half_u_n + 1),
                                         c_beta(half_u_n + 1);

                for (GLint j = 0; j <= half_u_n; j++)
                {
                    #pragma omp parallel for
                    for (GLint i = 0; i < u_dimension; i++)
                    {
                        const Cartesian3 &cp = _data(i, l);

                        GLdouble b_i_j_u_alpha = _S[variable::U](ECSpace::B_BASIS,
                                                                 i, j, u_alpha);
                        GLdouble b_i_j_gamma   = _S[variable::U](ECSpace::B_BASIS,
                                                                 i, j, gamma);
                        GLdouble b_i_j_u_beta  = _S[variable::U](ECSpace::B_BASIS,
                                                                 i, j, u_beta);

                        #pragma omp critical (c_derivatives_alpha_gamma_beta)
                        c_alpha[j] += cp * b_i_j_u_alpha;
                        c_gamma[j] += cp * b_i_j_gamma;
                        c_beta [j] += cp * b_i_j_u_beta;
                    }
                }

                L->_data(u_n, l) = c_gamma[0], L->_data(0, l)   = _data(0,   l);
                R->_data(0, l)   = c_gamma[0], R->_data(u_n, l) = _data(u_n, l);

                for (GLint i = 1; i <= half_u_n; i++)
                {
                    Cartesian3 &lambda_u_n_minus_i = L->_data(u_n - i, l);
                    Cartesian3 &rho_i              = R->_data(i, l);

                    lambda_u_n_minus_i = rho_i = c_gamma[i];

                    #pragma omp parallel for
                    for (GLint j = 0; j < i; j++)
                    {
                        Cartesian3 backward = L->_data(u_n - j, l);
                        Cartesian3 forward  = R->_data(j, l);

                        backward *= lambda_space(ECSpace::B_BASIS, u_n - j, i,
                                                 gamma);
                        forward  *= rho_space   (ECSpace::B_BASIS,       j, i,
                                                 gamma);

                        #pragma omp critical (subdivision_points_I)
                        lambda_u_n_minus_i -= backward;
                        rho_i              -= forward;
                    }

                    #pragma omp critical (last_division)
                    lambda_u_n_minus_i /= lambda_space(ECSpace::B_BASIS, u_n - i, i,
                                                       gamma);
                    rho_i              /= rho_space   (ECSpace::B_BASIS,       i, i,
                                                       gamma);
                }

                for (GLint i = 1; i <= (u_n - 1) / 2; i++)
                {
                    Cartesian3 &lambda_i         = L->_data(i, l);
                    Cartesian3 &rho_u_n_minus_i  = R->_data(u_n - i, l);

                    lambda_i         = c_alpha[i];
                    rho_u_n_minus_i  = c_beta [i];

                    #pragma omp parallel for
                    for (GLint j = 0; j < i; j++)
                    {
                        Cartesian3 backward = R->_data(u_n - j, l);
                        Cartesian3 forward  = L->_data(j, l);

                        backward *= rho_space   (ECSpace::B_BASIS, u_n - j, i, u_beta);
                        forward  *= lambda_space(ECSpace::B_BASIS,       j, i, u_alpha);

                        #pragma omp critical (subdivision_points_II)
                        lambda_i        -= forward;
                        rho_u_n_minus_i -= backward;
                    }


                    #pragma omp critical (last_division)
                    lambda_i         /= lambda_space(ECSpace::B_BASIS,       i, i,
                                                     u_alpha);
                    rho_u_n_minus_i  /= rho_space   (ECSpace::B_BASIS, u_n - i, i,
                                                     u_beta);
                }
            }

            return result;
        }
        else
        {
            GLdouble v_alpha = _S[variable::V].alpha(), v_beta = _S[variable::V].beta();

            if (gamma <= v_alpha || gamma >= v_beta)
            {
                return nullptr;
            }

            RowMatrix<SP<BSurface3>::Default>* result =
                    new (nothrow) RowMatrix<SP<BSurface3>::Default>(2);

            if (!result)
            {
                return nullptr;
            }

            ECSpace lambda_space(_S[variable::V]), rho_space(_S[variable::V]);

            try
            {
                if (!lambda_space.setDefinitionDomain(
                        v_alpha, gamma,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits) ||
                    !rho_space.setDefinitionDomain(
                            gamma, v_beta,
                            check_for_ill_conditioned_matrices,
                            expected_correct_significant_digits))
                {
                    delete result, result = nullptr;
                    return nullptr;
                }
            }
            catch (Exception &e)
            {
                delete result, result = nullptr;
                throw e;
            }

            (*result)[0] = SP<BSurface3>::Default(
                    new (nothrow) BSurface3(_S[variable::U], lambda_space));

            (*result)[1] = SP<BSurface3>::Default(
                    new (nothrow) BSurface3(_S[variable::U], rho_space));

            if (!(*result)[0] || !(*result)[1])
            {
                delete result, result = nullptr;
                return nullptr;
            }

            SP<BSurface3>::Default &L = (*result)[0];
            SP<BSurface3>::Default &R = (*result)[1];


            GLint v_dimension = _S[variable::V].dimension();
            GLint v_n         = v_dimension - 1;
            GLint half_v_n    = v_n / 2;

            GLint u_dimension = _S[variable::U].dimension();
            GLint u_n         = u_dimension - 1;

            for (GLint k = 0; k <= u_n; k++)
            {
                ColumnMatrix<Cartesian3> c_alpha(half_v_n + 1),
                                         c_gamma(half_v_n + 1),
                                         c_beta(half_v_n + 1);

                for (GLint j = 0; j <= half_v_n; j++)
                {
                    #pragma omp parallel for
                    for (GLint i = 0; i < v_dimension; i++)
                    {
                        const Cartesian3 &cp   = _data(k, i);
                        GLdouble b_i_j_v_alpha = _S[variable::V](ECSpace::B_BASIS, i, j,
                                                                 v_alpha);
                        GLdouble b_i_j_gamma   = _S[variable::V](ECSpace::B_BASIS, i, j,
                                                                 gamma);
                        GLdouble b_i_j_v_beta  = _S[variable::V](ECSpace::B_BASIS, i, j,
                                                                 v_beta);

                        #pragma omp critical
                        c_alpha[j] += cp * b_i_j_v_alpha;
                        c_gamma[j] += cp * b_i_j_gamma;
                        c_beta [j] += cp * b_i_j_v_beta;
                    }
                }

                L->_data(k, v_n) = c_gamma[0], L->_data(k, 0)   = _data(k, 0);
                R->_data(k, 0)   = c_gamma[0], R->_data(k, v_n) = _data(k, v_n);

                for (GLint i = 1; i <= half_v_n; i++)
                {
                    Cartesian3 &lambda_v_n_minus_i = L->_data(k, v_n - i);
                    Cartesian3 &rho_i              = R->_data(k, i);

                    lambda_v_n_minus_i = rho_i = c_gamma[i];

                    #pragma omp parallel for
                    for (GLint j = 0; j < i; j++)
                    {
                        Cartesian3 backward = L->_data(k, v_n - j);
                        Cartesian3 forward  = R->_data(k, j);

                        backward *= lambda_space(ECSpace::B_BASIS, v_n - j, i, gamma);
                        forward  *= rho_space   (ECSpace::B_BASIS,       j, i, gamma);

                        #pragma omp critical (subdivision_points_I)
                        lambda_v_n_minus_i -= backward;
                        rho_i              -= forward;
                    }

                    #pragma omp critical (last_division)
                    lambda_v_n_minus_i /= lambda_space(ECSpace::B_BASIS, v_n - i, i,
                                                       gamma);
                    rho_i              /= rho_space   (ECSpace::B_BASIS,       i, i,
                                                       gamma);
                }

                for (GLint i = 1; i <= (v_n - 1) / 2; i++)
                {
                    Cartesian3 &lambda_i         = L->_data(k, i);
                    Cartesian3 &rho_v_n_minus_i  = R->_data(k, v_n - i);

                    lambda_i         = c_alpha[i];
                    rho_v_n_minus_i  = c_beta [i];

                    #pragma omp parallel for
                    for (GLint j = 0; j < i; j++)
                    {
                        Cartesian3 backward = R->_data(k, v_n - j);
                        Cartesian3 forward  = L->_data(k, j);

                        backward *= rho_space   (ECSpace::B_BASIS, v_n - j, i, v_beta);
                        forward  *= lambda_space(ECSpace::B_BASIS,       j, i, v_alpha);

                        #pragma omp critical (subdivision_points_II)
                        lambda_i        -= forward;
                        rho_v_n_minus_i -= backward;
                    }

                    #pragma omp critical (last_division)
                    lambda_i         /= lambda_space(ECSpace::B_BASIS,       i, i,
                                                     v_alpha);
                    rho_v_n_minus_i  /= rho_space   (ECSpace::B_BASIS, v_n - i, i,
                                                     v_beta);
                }
            }

            return result;
        }
    } //(*@\label{BSurface3::performSubdivision::end}@*)

    //(*@\Green{// Given an ordinary integral surface of type (\mref{eq:ordinary_integral_surface}) the method updates the control points $[\mathbf{p}_{i_0,i_1}]_{i_0 = 0,\,i_1=0}^{n_0,\,n_1}$}@*)
    //(*@\Green{// of the B-surface (\mref{eq:B-surface}) in order to ensure control point based exact description. The ordinary surface}@*)
    //(*@\Green{// (\mref{eq:ordinary_integral_surface}) is described by an instance of the class OrdinarySurfaceCoefficients, while the control points}@*)
    //(*@\Green{// $[\mathbf{p}_{i_0,i_1}]_{i_0 = 0,\,i_1=0}^{n_0,\,n_1}$ that have to be updated are stored in the inherited data structure}@*)
    //(*@\Green{// Matrix$<$Cartesian3$>$ TensorProductSurface3::\_data. The method is based on formulas of Theorem \mref{thm:integral_surfaces}.}@*)
    GLboolean BSurface3::updateControlPointsForExactDescription(//(*@\label{BSurface3::updateControlPointsForExactDescription::start}@*)
            const OrdinarySurfaceCoefficients &lambda)
    {
        if (!_T[variable::U] || !_T[variable::V] ||
            _S[variable::U].dimension() != lambda.dimension(variable::U) ||
            _S[variable::V].dimension() != lambda.dimension(variable::V))
        {
            return GL_FALSE;
        }

        #pragma omp parallel for
        for (GLint j_0_j_1 = 0;
             j_0_j_1 < _S[variable::U].dimension() * _S[variable::V].dimension();
             j_0_j_1++)
        {
            GLint j_0 = j_0_j_1 / _S[variable::V].dimension();
            GLint j_1 = j_0_j_1 % _S[variable::V].dimension();

            Cartesian3 &p = _data(j_0, j_1);

            p[0] = p[1] = p[2] = 0.0;

            for (GLint k = 0; k < 3; k++)
            {
                for (GLint zeta = 0; zeta < lambda.sigma(k); zeta++)
                {
                    GLdouble u_sum = 0.0;

                    for (GLint i_0 = 0; i_0 < _S[variable::U].dimension(); i_0++)
                    {
                        u_sum += lambda(k, zeta, variable::U, i_0) *
                                 (*_T[variable::U])(i_0, j_0);
                    }

                    GLdouble v_sum = 0.0;

                    for (GLint i_1 = 0; i_1 < _S[variable::V].dimension(); i_1++)
                    {
                        v_sum += lambda(k, zeta, variable::V, i_1) *
                                 (*_T[variable::V])(i_1, j_1);
                    }

                    p[k] += u_sum * v_sum;
                }
            }
        }

        return GL_TRUE;
    }//(*@\label{BSurface3::updateControlPointsForExactDescription::end}@*)

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    BSurface3* BSurface3::clone() const
    {
        return new (std::nothrow) BSurface3(*this);
    }
}
