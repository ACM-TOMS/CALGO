//----------------------------------------------------------------------------------
// File:        EC/BCurves3.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "BCurves3.h"

#include "../Core/Exceptions.h"

#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

namespace cagd
{
    //(*@\Green{// special constructor}@*)
    BCurve3::BCurve3(const ECSpace &S, GLenum data_usage_flag):
            LinearCombination3(S.alpha(), S.beta(), S.dimension(), data_usage_flag),
            _S(S),
            _T(_S.basisTransformationFromNBToOrdinary())
    {
        if (!_T)
        {
            throw Exception("The transformation matrix that maps the normalized "
                            "B-basis of the given EC space to its ordinary basis "
                            "could not be created!");
        }
    }

    //(*@\Green{// Calculates all normalized B-basis functions at the parameter value $u \in \left[\alpha, \beta\right]$.}@*)
    GLboolean BCurve3::blendingFunctionValues(
            GLdouble u, RowMatrix<GLdouble> &values) const
    {
        if (u < _S.alpha() || u > _S.beta())
        {
            values.resizeColumns(0);
            return GL_FALSE;
        }

        GLint dimension = _S.dimension();

        values.resizeColumns(dimension);

        #pragma omp parallel for
        for (GLint i = 0; i < dimension; i++)
        {
            values[i] = _S(ECSpace::B_BASIS, i, 0, u);
        }

        return GL_TRUE;
    }

    //(*@\Green{// Differentiates the given B-curve up to a maximal order at the parameter value $u \in \left[\alpha,\beta\right]$. The zeroth and}@*)
    //(*@\Green{// higher order derivatives will be stored in the rows of the column matrix referenced by the variable d.}@*)
    GLboolean BCurve3::calculateDerivatives(
            GLint maximum_order_of_derivatives, GLdouble u, Derivatives &d) const
    {
        assert("The given maximum order of derivatives should be non-negative!" &&
               maximum_order_of_derivatives >= 0);

        if (maximum_order_of_derivatives < 0 || u < _S.alpha() || u > _S.beta())
        {
            d.resizeRows(0);
            return GL_FALSE;
        }

        GLint dimension = _S.dimension();

        d.resizeRows(maximum_order_of_derivatives + 1);
        d.loadNullVectors();

        //(*@\Green{// implementation of the formula $\mathbf{c}_n^{\left(j\right)}\left(u\right) = \sum_{i=0}^n \mathbf{p}_i b_{n,i}^{\left(j\right)}\left(u\right)$, where the control points $[\mathbf{p}_i]_{i=0}^n$ are stored}@*)
        //(*@\Green{// in the inherited column matrix \_data}@*)
        #pragma omp parallel for
        for (GLint j = 0; j <= maximum_order_of_derivatives; j++)
        {
            for (GLint i = 0; i < dimension; i++)
            {
                d[j] += _data[i] * _S(ECSpace::B_BASIS, i, j, u);
            }
        }

        return GL_TRUE;
    }

    //(*@\Green{// Overloaded function operator that evaluates the derivative $b_{n,i}^{\left(j\right)}\left(u\right)$, where $u \in \left[\alpha, \beta\right], ~i=0,1,\ldots,n,~j\in\mathbb{N}$.}@*)
    GLdouble BCurve3::operator ()(GLint i, GLint j, GLdouble u) const
    {
        assert("The given B-basis function index is out of bounds!" &&
               (i >= 0 && i < _S.dimension()));
        assert("The given differentiation order should be non-negative!" && j >= 0);

        return _S(ECSpace::B_BASIS, i, j, u);
    }

    //(*@\Green{// Overloaded member functions for general order elevation. In order to generate a higher dimensional EC}@*)
    //(*@\Green{// space $\mathbb{S}_{n+q}^{\alpha,\beta}$ that fulfills the condition $\mathbb{S}_n^{\alpha,\beta}\subset\mathbb{S}_{n+q}^{\alpha,\beta}$, where $q\geq 1$ and $0<\beta-\alpha<\min\left\{\ell^{\prime}\left(\mathbb{S}_n^{\alpha,\beta}\right), \ell^{\prime}\left(\mathbb{S}_{n+q}^{\alpha,\beta}\right)\right\}$,}@*)
    //(*@\Green{// one has either to define a new complex root $z=a\pm\mathbf{i}b$ of multiplicity $m \geq 1$, or to specify an existing}@*)
    //(*@\Green{// root with a multiplicity that is greater than its former value $m'\geq 1$. If $z \in \mathbb{C}\setminus\mathbb{R}$, then $q=2\left(m-m'\right)$,}@*)
    //(*@\Green{// otherwise $q=m-m'$, where $m'=0$ whenever $z$ denotes a nonexisting former root.}@*)
    //(*@\Green{// The implementation of the next two methods are based on Lemma \mref{lem:general_order_elevation}.}@*)
    BCurve3* BCurve3::performOrderElevation(
            GLdouble a, GLdouble b, GLint multiplicity,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits) const
    {
        return performOrderElevation(CharacteristicPolynomial::Zero(a, b, multiplicity),
                                     check_for_ill_conditioned_matrices,
                                     expected_correct_significant_digits);
    }

    BCurve3* BCurve3::performOrderElevation(
            const CharacteristicPolynomial::Zero &zero,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits) const
    {
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);

        if (expected_correct_significant_digits < 0)
        {
            return nullptr;
        }

        //(*@\Green{//At first, we try to clone our curve object together with is embedded EC space.}@*)
        BCurve3* result = clone();

        if (!result)
        {
            return nullptr;
        }

        try
        {
            //(*@\Green{// Then, we try to update the factorization of the characteristic polynomial (\mref{eq:characteristic_polynomial}) associated with the}@*)
            //(*@\Green{// differential equation (\mref{eq:differential_equation}). If the factorization is changed, the ordinary basis and the normalized}@*)
            //(*@\Green{// B-basis stored in the clone will be automatically updated.}@*)
            if (!result->_S.insertZero(zero,
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

        //(*@\Green{// In what follows, we compare the dimensions of the original and new EC spaces.}@*)
        GLint old_dimension      = _S.dimension();
        GLint elevated_dimension = result->_S.dimension();

        //(*@\Green{// If there is no difference, we return the pointer that points to the cloned object.}@*)
        if (old_dimension == elevated_dimension)
        {
            return result;
        }

        GLdouble alpha = _S.alpha(), beta = _S.beta();

        result->_data.resizeRows(elevated_dimension);

        //(*@\Green{// If the difference is $1$, we apply the order elevation formulas of Lemma \mref{lem:general_order_elevation}:}@*)
        if (elevated_dimension - old_dimension == 1)//(*@\label{BCurve3::performOrderElevation::start}@*)
        {
            GLint n = old_dimension - 1;

            (*result)[0]             = _data[0];
            (*result)[old_dimension] = _data[n];

            //(*@\Green{// apply formula (\mref{eq:order_elevation_inner_points_first_half}), i.e., $\mathbf{p}_{1,i}=\left(1-\frac{b_{n,i}^{\left(i\right)}\left(\alpha\right)}{b_{n+1,i}^{\left(i\right)}\left(\alpha\right)}\right)\mathbf{p}_{i-1}+\frac{b_{n,i}^{\left(i\right)}\left(\alpha\right)}{b_{n+1,i}^{\left(i\right)}\left(\alpha\right)}\mathbf{p}_i,~i=1,\ldots,\left\lfloor \frac{n}{2} \right\rfloor$;}@*)
            #pragma omp parallel for
            for (GLint i = 1; i <= n / 2; i++)
            {
                GLdouble ratio = (*this)(i, i, alpha) / (*result)(i, i, alpha);
                (*result)[i]   = (1 - ratio) * _data[i - 1] + ratio * _data[i];
            }

            //(*@\Green{// apply formula (\mref{eq:order_elevation_inner_points_second_half}), i.e., $\mathbf{p}_{1,n+1-i} = \frac{b_{n,n-i}^{\left(i\right)}\left(\beta\right)}{b_{n+1,n+1-i}^{\left(i\right)}\left(\beta\right)} \mathbf{p}_{n-i} + \left(1-\frac{b_{n,n-i}^{\left(i\right)}\left(\beta\right)}{b_{n+1,n+1-i}^{\left(i\right)}\left(\beta\right)}\right)\mathbf{p}_{n+1-i},$}@*)
            //(*@\Green{// $i=1,\ldots,\left\lfloor\frac{n+1}{2}\right\rfloor$.}@*)
            #pragma omp parallel for
            for (GLint i = 1; i <= old_dimension / 2; i++)
            {
                GLdouble ratio               = (*this)(n - i, i, beta) /
                                               (*result)(old_dimension - i, i, beta);
                (*result)[old_dimension - i] = ratio * _data[n - i] +
                                               (1 - ratio) * _data[old_dimension - i];
            }

            return result;
        }//(*@\label{BCurve3::performOrderElevation::end}@*)

        //(*@\Green{// If the difference between the dimensions is strictly greater than $1$, we reduce the order elevation to a curve}@*)
        //(*@\Green{// interpolation problem.}@*)
        ColumnMatrix<GLdouble>    knot_vector(elevated_dimension);
        ColumnMatrix<Cartesian3>  sample(elevated_dimension);

        GLboolean sampling_aborted = GL_FALSE;
        GLdouble  step = (beta - alpha) / (elevated_dimension - 1);

        #pragma omp parallel for
        for (GLint i = 0; i < elevated_dimension; i++)
        {
            #pragma omp flush (sampling_aborted)
            if (!sampling_aborted)
            {
                //(*@\Green{// uniform subdivision points of the interval $\left[\alpha,\beta\right]$}@*)
                knot_vector[i] = min(alpha + i * step, beta);

                Derivatives d;
                if (!calculateDerivatives(0, knot_vector[i], d))
                {
                    sampling_aborted = GL_TRUE;
                    #pragma omp flush (sampling_aborted)
                }
                else
                {
                    sample[i] = d[0]; //(*@\Green{// sample point of the original curve}@*)
                }
            }
        }

        //(*@\Green{// If the sampling was aborted or the curve interpolation problem cannot be solved, the clone will be}@*)
        //(*@\Green{// deleted and a null pointer will be returned.}@*)
        if (sampling_aborted || !result->updateDataForInterpolation(knot_vector, sample))
        {
            delete result, result = nullptr;
        }

        return result;
    }

    //(*@\Green{// Using the general B-algorithm presented in Theorem \mref{thm:general_subdivision}, the method subdivides the given B-curve}@*)
    //(*@\Green{// into two arcs of the same order at the parameter value $\gamma \in \left(\alpha,\beta\right)$. In case of success, the output is a}@*)
    //(*@\Green{// nonzero raw pointer to a $2$-element row matrix that stores $2$ smart pointers (based on the deep copy}@*)
    //(*@\Green{// ownership policy) that point to the left and right arcs of the subdivided curve.}@*)
    RowMatrix<SP<BCurve3>::Default>* BCurve3::performSubdivision(//(*@\label{BCurve3::performSubdivision::start}@*)
            GLdouble gamma,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits) const
    {
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);

        GLdouble alpha  = _S.alpha(), beta = _S.beta();

        //(*@\Green{// check whether $\gamma \in \left(\alpha,\beta\right)$}@*)
        if (gamma <= alpha || gamma >= beta || expected_correct_significant_digits < 0)
        {
            return nullptr;
        }

        //(*@\Green{// try to dynamically allocate a 2-element row matrix that stores smart pointers to the left and right arcs}@*)
        RowMatrix<SP<BCurve3>::Default> *result =
                new (nothrow) RowMatrix<SP<BCurve3>::Default>(2);

        if (!result)
        {
            return nullptr;
        }

        //(*@\Green{// Using the same EC space \_S, the left and right arcs of the subdivided curve are defined over the}@*)
        //(*@\Green{// subintervals $\left[\alpha,\gamma\right]$ and $\left[\gamma,\beta\right]$, respectively.}@*)
        (*result)[0] = SP<BCurve3>::Default(new (nothrow) BCurve3(_S, _data_usage_flag));
        (*result)[1] = SP<BCurve3>::Default(new (nothrow) BCurve3(_S, _data_usage_flag));

        if (!(*result)[0] || !(*result)[1])
        {
            delete result, result = nullptr;
            return nullptr;
        }

        //(*@\Green{// In order to simplify the lines below and to ease the comparison to the proposed formulas, we will use}@*)
        //(*@\Green{// references to the left and right arcs.}@*)
        BCurve3 &lambda        = *(*result)[0], &rho          = *(*result)[1];
        BCurve3 &b_alpha_gamma = lambda,        &b_gamma_beta = rho;

        try
        {
            //(*@\Green{// Set the definition domains of the arcs, i.e., try to generate on both sections the corresponding}@*)
            //(*@\Green{// ordinary bases, the normalized B-bases and the transformation matrices as well.}@*)
            if (!lambda.setDefinitionDomain(alpha, gamma,
                                            check_for_ill_conditioned_matrices,
                                            expected_correct_significant_digits) ||
                !rho.setDefinitionDomain(gamma, beta,
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

        GLint dimension = _S.dimension();
        GLint n         = dimension - 1;

        //(*@\Green{// Try to calculate the derivatives $\{\mathbf{c}_n^{\left(i\right)}\left(\alpha\right),\mathbf{c}_n^{\left(i\right)}\left(\gamma\right),\mathbf{c}_n^{\left(i\right)}\left(\beta\right)\}_{i=0}^{\left\lfloor\frac{n}{2}\right\rfloor}$.}@*)
        Derivatives c_alpha, c_gamma, c_beta;

        Derivatives *c_ptr[3] = {&c_alpha, &c_gamma, &c_beta};
        GLdouble    u[3] = {alpha, gamma, beta};

        GLboolean c_alpha_gamma_beta_derivatives_aborted = GL_FALSE;

        #pragma omp parallel for
        for (GLint i = 0; i < 3; i++)
        {
            #pragma omp flush (c_alpha_gamma_beta_derivatives_aborted)
            if (!c_alpha_gamma_beta_derivatives_aborted)
            {
                if (!calculateDerivatives(n / 2, u[i], *c_ptr[i]))
                {
                    c_alpha_gamma_beta_derivatives_aborted = GL_TRUE;
                    #pragma omp flush (c_alpha_gamma_beta_derivatives_aborted)
                }
            }
        }

        if (c_alpha_gamma_beta_derivatives_aborted)
        {
            delete result, result = nullptr;
            return nullptr;
        }

        //(*@\Green{// Set the initial conditions of recursive formulas (\mref{eq:left_subdivision_points_second_half}) and (\mref{eq:right_subdivision_points_first_half}).}@*)
        lambda[n] = c_gamma[0], lambda[0]  = _data[0];
        rho[0]    = c_gamma[0], rho[n]     = _data[n];

        //(*@\Green{// Apply recursive formulas (\mref{eq:left_subdivision_points_second_half}) and (\mref{eq:right_subdivision_points_first_half}), i.e., }@*)
        //(*@\Green{// $\Red{\boldsymbol{\lambda}_{n-i}\left(\gamma\right)} = \frac{1}{b_{n,n-i}^{\left(i\right)}\left(\gamma;\alpha,\gamma\right)}\left(\mathbf{c}_n^{\left(i\right)}\left(\gamma\right)-\sum_{j = 0}^{i-1}\Red{\boldsymbol{\lambda}_{n-j}\left(\gamma\right)}\cdot b_{n,n-j}^{\left(i\right)}\left(\gamma; \alpha, \gamma\right)\right),~i=1,\ldots,\left\lfloor\frac{n}{2}\right\rfloor$,}@*)
        //(*@\Green{// \hspace{0.35cm}$\Blue{\boldsymbol{\varrho}_i\left(\gamma\right)}=\frac{1}{b_{n,i}^{\left(i\right)}\left(\gamma;\gamma,\beta\right)}\left(\mathbf{c}_n^{\left(i\right)}\left(\gamma\right)-\sum_{j = 0}^{i-1}\Blue{\boldsymbol{\varrho}_j\left(\gamma\right)}\cdot b_{n,j}^{\left(i\right)}\left(\gamma;\gamma,\beta\right)\right),~i=1,\ldots,\left\lfloor\frac{n}{2}\right\rfloor$.}@*)
        for (GLint i = 1; i <= n / 2; i++)
        {
            Cartesian3 &lambda_n_minus_i = lambda[n - i];
            Cartesian3 &rho_i            = rho   [i];

            lambda_n_minus_i = rho_i = c_gamma[i];

            #pragma omp parallel for
            for (GLint j = 0; j < i; j++)
            {
                Cartesian3 backward = lambda[n - j] * b_alpha_gamma(n - j, i, gamma);
                Cartesian3 forward  = rho   [j]     * b_gamma_beta (    j, i, gamma);

                #pragma omp critical (inner_points)
                lambda_n_minus_i -= backward;
                rho_i            -= forward;
            }

            lambda_n_minus_i /= b_alpha_gamma(n - i, i, gamma);
            rho_i            /= b_gamma_beta (    i, i, gamma);
        }

        //(*@\Green{// Apply recursive formulas (\mref{eq:left_subdivision_points_first_half}) and (\mref{eq:right_subdivision_points_second_half}), i.e.,}@*)
        //(*@\Green{// \hspace{0.295cm}$\Red{\boldsymbol{\lambda}_i\left(\gamma\right)}=\frac{1}{b_{n,i}^{\left(i\right)}\left(\alpha;\alpha,\gamma\right)}\left(\mathbf{c}_n^{\left(i\right)}\left(\alpha\right)-\sum_{j = 0}^{i-1} \Red{\boldsymbol{\lambda}_j\left(\gamma\right)}\cdot b_{n,j}^{\left(i\right)}\left(\alpha; \alpha, \gamma\right)\right),~i = 1,\ldots,\left\lfloor\frac{n-1}{2}\right\rfloor$,}@*)
        //(*@\Green{// $\Blue{\boldsymbol{\varrho}_{n-i}\left(\gamma\right)}=\frac{1}{b_{n,n-i}^{\left(i\right)}\left(\beta;\gamma,\beta\right)}\left(\mathbf{c}_n^{\left(i\right)}\left(\beta\right)-\sum_{j=0}^{i-1}\Blue{\boldsymbol{\varrho}_{n-j}\left(\gamma\right)}\cdot b_{n,n-j}^{\left(i\right)}\left(\beta;\gamma,\beta\right)\right),~i=1,\ldots,\left\lfloor\frac{n-1}{2}\right\rfloor$.}@*)
        for (GLint i = 1; i <= (n - 1) / 2; i++)
        {
            Cartesian3 &lambda_i         = lambda [i];
            Cartesian3 &rho_n_minus_i    = rho[n - i];

            lambda_i         = c_alpha[i];
            rho_n_minus_i    = c_beta [i];

            #pragma omp parallel for
            for (GLint j = 0; j < i; j++)
            {
                Cartesian3 backward = rho[n - j] * b_gamma_beta (n - j, i, beta);
                Cartesian3 forward  = lambda [j] * b_alpha_gamma(    j, i, alpha);

                #pragma omp critical (inner_points)
                lambda_i      -= forward;
                rho_n_minus_i -= backward;
            }

            lambda_i       /= b_alpha_gamma(    i, i, alpha);
            rho_n_minus_i  /= b_gamma_beta (n - i, i, beta);
        }

        return result;
    }//(*@\label{BCurve3::performSubdivision::end}@*)

    //(*@\Green{// A redeclared and redefined inherited virtual method. If one alters the endpoints of the definition domain $\left[\alpha, \beta\right]$}@*)
    //(*@\Green{// both bases $\mathcal{F}_{n}^{\alpha,\beta}$ and $\mathcal{B}_{n}^{\alpha,\beta}$ of $\mathbb{S}_n^{\alpha,\beta}$ have to be updated. Note that, the new interval length should also satisfy}@*)
    //(*@\Green{// the condition $0<\beta-\alpha<\ell^{\prime}\left(\mathbb{S}_{n}^{\alpha,\beta}\right)$. The transformation matrix $[t_{i,j}^n]_{i=0,\,j=0}^{n,\,n}$ that maps $\mathcal{B}_{n}^{\alpha,\beta}$ to $\mathcal{F}_{n}^{\alpha,\beta}$ is also}@*)
    //(*@\Green{// updated.}@*)
    GLboolean BCurve3::setDefinitionDomain(
            GLdouble alpha, GLdouble beta,
            bool check_for_ill_conditioned_matrices,
            GLint expected_correct_significant_digits)
    {
        assert("The number of expected correct significant digits should be "
               "non-negative!" && expected_correct_significant_digits >= 0);

        if (expected_correct_significant_digits < 0 ||
            !LinearCombination3::setDefinitionDomain(alpha, beta) ||
            !_S.setDefinitionDomain(alpha, beta,
                                    check_for_ill_conditioned_matrices,
                                    expected_correct_significant_digits))
        {
            return GL_FALSE;
        }

        _T = SP<RealMatrix>::Default(_S.basisTransformationFromNBToOrdinary());

        return (_T ? GL_TRUE : GL_FALSE);
    }

    //(*@\Green{// Given an ordinary integral curve $\mathbf{c}\left(  u\right)  =\sum_{i=0}^{n}\boldsymbol{\lambda}_{i}\varphi_{n,i}\left(  u\right)  ,~u\in\left[  \alpha,\beta\right],~0<\beta-\alpha<\ell^{\prime}\left(\mathbb{S}_{n}^{\alpha,\beta}\right),~\boldsymbol{\lambda}_{i}\in\mathbb{R}^{3}$,}@*)
    //(*@\Green{// the method determines the coefficients $[\mathbf{p}_j]_{j = 0}^n$ of the normalized B-basis functions $\left\{b_{n,j}\left(u\right) : u \in \left[\alpha, \beta\right]\right\}_{j=0}^{n}$}@*)
    //(*@\Green{// such that $\mathbf{c}\left(u\right) \equiv \sum_{j = 0}^n \mathbf{p}_j b_{n,j}\left(u\right)$, $\forall u \in \left[\alpha, \beta\right]$. Note that, the control points $[\mathbf{p}_j]_{j = 0}^n$ that have to be}@*)
    //(*@\Green{// updated are stored in the inherited data structure ColumnMatrix$<$Cartesian3$>$ LinearCombination3::\_data.}@*)
    //(*@\Green{// The method's implementation is based on Theorems \mref{thm:efficient_basis_transformation} and \mref{thm:integral_curves}.}@*)
    GLboolean BCurve3::updateControlPointsForExactDescription(//(*@\label{BCurve3::updateControlPointsForExactDescription::start}@*)
            const RowMatrix<Cartesian3> &lambda)
    {
        if (!_T || _S.factorizationOfTheCharacteristicPolynomialChanged())
        {
            return GL_FALSE;
        }

        GLint dimension = _S.dimension();

        if (_T->rowCount() != dimension || lambda.columnCount() != dimension )
        {
            return GL_FALSE;
        }

        //(*@\Green{// Unknown control points $[\mathbf{p}_j]_{j=0}^n$ are calculated by means of \cite[Corollary 2.1, p. 43]{Roth2015b}, i.e.,}@*)
        //(*@\Green{// $\mathbf{p}_j = \sum_{i=0}^n\boldsymbol{\lambda}_i t_{i,j}^n,~j=0,1,\ldots,n.$ (Consider also Theorem \mref{thm:integral_curves}.)}@*)
        #pragma omp parallel for
        for (GLint j = 0; j < dimension; j++)
        {
            Cartesian3 &p_j = _data[j];

            p_j[0] = p_j[1] = p_j[2] = 0.0;

            for (GLint i = 0; i < dimension; i++)
            {
                p_j += lambda[i] * (*_T)(i, j);
            }
        }

        return GL_TRUE;
    }//(*@\label{BCurve3::updateControlPointsForExactDescription::end}@*)

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    BCurve3* BCurve3::clone() const
    {
        return new (nothrow) BCurve3(*this);
    }
}

