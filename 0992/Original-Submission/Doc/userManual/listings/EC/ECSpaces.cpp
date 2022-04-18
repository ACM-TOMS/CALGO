//----------------------------------------------------------------------------------
// File:        EC/ECSpaces.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "ECSpaces.h"

#include "../Core/Exceptions.h"
#include "../Core/Math/Constants.h"
#include "../Core/Math/Matrices.h"
#include "../Core/Math/RealMatrixDecompositions.h"

#include "../Core/Utilities.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <new>

using namespace std;

namespace cagd
{
    (*@\Green{// default/special constructors based on a possible complex root $z = a \pm \mathbf{i}b$ of order $m\geq 1$, where}@*)
    (*@\Green{// $r\in\left\{0,\ldots, m-1\right\}$; if the boolean parameter cosine is true, the generated ordinary basis function}@*)
    (*@\Green{// corresponds to $u^{r} \cdot \left[e^{au}\right] \cdot \cos(bu)$, otherwise to $u^{r} \cdot \left[e^{au}\right] \cdot \sin(bu)$}@*)
    ECSpace::_OrdinaryBasisFunction::_OrdinaryBasisFunction(
            double a, double b, int r, bool cosine):
            _type(cosine ? AET_COSINE : AET_SINE),
            _a(a), _b(b),
            _r(r)
    {
        if (_a != 0.0 && _b == 0.0)
        {
            if (cosine)
            {
                _type = AE_COSINE;
            }
            else
            {
                _type = AE_SINE; (*@\Green{// Note that the enumerator AE\_SINE denotes the constant zero function.}@*)
            }
        }

        if (_a == 0.0 && _b != 0.0)
        {
            if (cosine)
            {
                _type = AT_COSINE;
            }
            else
            {
                _type = AT_SINE;
            }
        }

        if (_a == 0.0 && _b == 0.0)
        {
            if (cosine)
            {
                _type = P_COSINE;
            }
            else
            {
                _type = P_SINE; (*@\Green{// Note that the enumeratorr P\_SINE denotes the constant zero function.}@*)
            }
        }
    }

    (*@\Green{// overloaded function operator that evaluates the $j$th order ($j\geq 0$) derivative of the stored ordinary}@*)
    (*@\Green{// basis function at the parameter value $u \in \left[\alpha, \beta\right]$}@*)
    double ECSpace::_OrdinaryBasisFunction::operator ()(int j, double u) const
    {
        double result = 0.0;

        if (j == 0) (*@\Green{// zeroth order derivatives of ordinary basis functions}@*)
        {
            switch (_type)
            {
            case AET_COSINE:
                result = pow(u, _r) * exp(_a * u) * cos(_b * u);
                break;

            case AET_SINE:
                result = pow(u, _r) * exp(_a * u) * sin(_b * u);
                break;

            case AE_COSINE:
                result = pow(u, _r) * exp(_a * u);
                break;

            case AE_SINE:
                result = 0.0;
                break;

            case AT_COSINE:
                result = pow(u, _r) * cos(_b * u);
                break;

            case AT_SINE:
                result = pow(u, _r) * sin(_b * u);
                break;

            case P_COSINE:
                result = pow(u, _r);
                break;

            case P_SINE:
                result = 0.0;
                break;
            }
        }
        else        (*@\Green{// higher order derivatives of ordinary basis functions}@*)
        {
            switch (_type)
            {
            case AET_COSINE:
                {
                    if (!_r)
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} e^{au} \cos(bu)$}@*)
                        double derivative = 0.0;

                        double exp_a_u    = exp(_a * u);
                        double a_power    = 1.0;

                        double b_power    = pow(_b, j);
                        int    b_exponent = j;

                        double cb         = cos(_b * u);
                        double sb         = sin(_b * u);

                        for (int k = 0; k <= j; k++)
                        {
                            double aux = 1.0;

                            switch (b_exponent % 4)
                            {
                            case 0: aux =  cb; break;
                            case 1: aux = -sb; break;
                            case 2: aux = -cb; break;
                            case 3: aux =  sb; break;
                            }

                            derivative += BC(j, k) * a_power * b_power * aux;

                            a_power    *= _a;
                            b_power    /= _b;
                            b_exponent --;
                        }

                        derivative *= exp_a_u;

                        result = derivative;
                    }
                    else
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} u^r e^{au} \cos(bu)$}@*)
                        double derivative = 0.0;

                        double u_coefficient = 1.0;
                        double u_power       = pow(u, _r);
                        int    u_exponent    = _r;

                        double cb            = cos(_b * u);
                        double sb            = sin(_b * u);

                        for (int l = 0; l <= min(j, _r); l++)
                        {
                            (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^{j-l}}{\mathrm{d}u^{j-l}} e^{au} \cos(bu)$}@*)
                            double exp_au_cos_bu_derivative = 0.0;

                            {
                                double exp_a_u    = exp(_a * u);
                                double a_power    = 1.0;

                                double b_power    = pow(_b, j - l);
                                int    b_exponent = j - l;

                                for (int k = 0; k <= j - l; k++)
                                {
                                    double aux = 1.0;

                                    switch (b_exponent % 4)
                                    {
                                    case 0: aux =  cb; break;
                                    case 1: aux = -sb; break;
                                    case 2: aux = -cb; break;
                                    case 3: aux =  sb; break;
                                    }

                                    exp_au_cos_bu_derivative += BC(j - l, k) * a_power *
                                                                b_power * aux;

                                    a_power    *= _a;
                                    b_power    /= _b;
                                    b_exponent --;
                                }

                                exp_au_cos_bu_derivative *= exp_a_u;
                            }

                            derivative += BC(j, l) * u_coefficient * u_power *
                                          exp_au_cos_bu_derivative;

                            u_coefficient *= u_exponent;
                            u_exponent    --;

                            if (u)
                            {
                                u_power /= u;
                            }
                            else
                            {
                                u_power = u_exponent ? 0.0 : 1.0;
                            }
                        }

                        result = derivative;
                    }
                }
                break;

            case AET_SINE:
                {
                    if (!_r)
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} e^{au} \sin(bu)$}@*)
                        double derivative = 0.0;

                        double exp_a_u    = exp(_a * u);
                        double a_power    = 1.0;

                        double b_power    = pow(_b, j);
                        int    b_exponent = j;

                        double cb         = cos(_b * u);
                        double sb         = sin(_b * u);

                        for (int k = 0; k <= j; k++)
                        {
                            double aux = 1.0;

                            switch (b_exponent % 4)
                            {
                            case 0: aux =  sb; break;
                            case 1: aux =  cb; break;
                            case 2: aux = -sb; break;
                            case 3: aux = -cb; break;
                            }

                            derivative += BC(j, k) * a_power * b_power * aux;

                            a_power    *= _a;
                            b_power    /= _b;
                            b_exponent --;
                        }

                        derivative *= exp_a_u;

                        result = derivative;
                    }
                    else
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} u^r e^{au} \sin(bu)$}@*)
                        double derivative    = 0.0;

                        double u_coefficient = 1.0;
                        double u_power       = pow(u, _r);
                        int    u_exponent    = _r;

                        double cb            = cos(_b * u);
                        double sb            = sin(_b * u);

                        for (int l = 0; l <= min(j, _r); l++)
                        {
                            (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^{j-l}}{\mathrm{d}u^{j-l}} e^{au} \sin(bu)$}@*)
                            double exp_au_sin_bu_derivative = 0.0;

                            {
                                double exp_a_u    = exp(_a * u);
                                double a_power    = 1.0;

                                double b_power    = pow(_b, j - l);
                                int    b_exponent = j - l;

                                for (int k = 0; k <= j - l; k++)
                                {
                                    double aux = 1.0;

                                    switch (b_exponent % 4)
                                    {
                                    case 0: aux =  sb; break;
                                    case 1: aux =  cb; break;
                                    case 2: aux = -sb; break;
                                    case 3: aux = -cb; break;
                                    }

                                    exp_au_sin_bu_derivative += BC(j - l, k) * a_power *
                                                                b_power * aux;

                                    a_power    *= _a;
                                    b_power    /= _b;
                                    b_exponent --;
                                }

                                exp_au_sin_bu_derivative *= exp_a_u;
                            }

                            derivative += BC(j, l) * u_coefficient * u_power *
                                          exp_au_sin_bu_derivative;

                            u_coefficient *= u_exponent;
                            u_exponent    --;

                            if (u)
                            {
                                u_power /= u;
                            }
                            else
                            {
                                u_power = u_exponent ? 0.0 : 1.0;
                            }
                        }

                        result = derivative;
                    }
                }
                break;

            case AE_COSINE:
                {
                    if (!_r)
                    {
                        return pow(_a, j) * exp(_a * u);
                    }
                    else
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} u^{r} e^{au}$}@*)
                        double derivative    = 0.0;

                        double u_coefficient = 1.0;
                        int    u_exponent    = _r;
                        double u_power       = pow(u, _r);

                        double exp_a_u       = exp(_a * u);
                        double a_power       = pow(_a, j);

                        for (int k = 0; k <= min(j, _r); k++)
                        {
                            derivative    += BC(j, k) * u_coefficient * u_power *
                                             a_power;

                            u_coefficient *= u_exponent;
                            u_exponent    --;

                            if (u)
                            {
                                u_power   /= u;
                            }
                            else
                            {
                                u_power = u_exponent ? 0.0 : 1.0;
                            }

                            a_power       /= _a;
                        }

                        derivative *= exp_a_u;

                        result = derivative;
                    }
                }
                break;

            case AE_SINE:
                result = 0.0;
                break;

            case AT_COSINE:
                {
                    if (!_r)
                    {
                        double sb    = (j > 0 ? sin(_b * u) : 0.0);
                        double cb    = cos(_b * u);
                        double power = pow(_b, j);

                        switch (j % 4)
                        {
                        case 0: return  power * cb;
                        case 1: return -power * sb;
                        case 2: return -power * cb;
                        case 3: return  power * sb;
                        }
                    }
                    else
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} u^{r} \cos(bu)$}@*)
                        double derivative    = 0.0;
                        double u_coefficient = 1.0;
                        double u_power       = pow(u, _r);
                        int    u_exponent    = _r;

                        double b_power       = pow(_b, j);
                        int    b_exponent    = j;

                        double cb            = cos(_b * u);
                        double sb            = (j > 0 ? sin(_b * u) : 0.0);

                        for (int k = 0; k <= min(j, _r); k++)
                        {
                            double aux = 1.0;

                            switch (b_exponent % 4)
                            {
                            case 0: aux =  cb; break;
                            case 1: aux = -sb; break;
                            case 2: aux = -cb; break;
                            case 3: aux =  sb; break;
                            }

                            derivative    += BC(j, k) * u_coefficient * u_power *
                                             b_power * aux;

                            u_coefficient *= u_exponent;
                            u_exponent    --;

                            if (u)
                            {
                                u_power   /= u;
                            }
                            else
                            {
                                u_power = u_exponent ? 0.0 : 1.0;
                            }

                            b_power       /= _b;
                            b_exponent    --;
                        }

                        result = derivative;
                    }
                }
                break;

            case AT_SINE:
                {
                    if (!_r)
                    {
                        double cb    = (j > 0 ? cos(_b * u) : 0.0);
                        double sb    = sin(_b * u);
                        double power = pow(_b, j);

                        switch (j % 4)
                        {
                        case 0: return  power * sb;
                        case 1: return  power * cb;
                        case 2: return -power * sb;
                        case 3: return -power * cb;
                        }
                    }
                    else
                    {
                        (*@\Green{// applying the Leibniz rule for the evaluation of the derivative $\frac{\mathrm{d}^j}{\mathrm{d}u^j} u^{r} \cos(bu)$}@*)
                        double derivative    = 0.0;
                        double u_coefficient = 1.0;
                        double u_power       = pow(u, _r);
                        int    u_exponent    = _r;

                        double b_power       = pow(_b, j);
                        int    b_exponent    = j;

                        double cb            = (j > 0 ? cos(_b * u) : 0.0);
                        double sb            = sin(_b * u);

                        for (int k = 0; k <= min(j, _r); k++)
                        {
                            double aux = 1.0;

                            switch (b_exponent % 4)
                            {
                            case 0: aux =  sb; break;
                            case 1: aux =  cb; break;
                            case 2: aux = -sb; break;
                            case 3: aux = -cb; break;
                            }

                            derivative    += BC(j, k) * u_coefficient * u_power *
                                             b_power * aux;

                            u_coefficient *= u_exponent;
                            u_exponent    --;

                            if (u)
                            {
                                u_power   /= u;
                            }
                            else
                            {
                                u_power = u_exponent ? 0.0 : 1.0;
                            }

                            b_power       /= _b;
                            b_exponent    --;
                        }

                        result = derivative;
                    }
                }
                break;

            case P_COSINE:
                {
                    if (j > _r)
                    {
                        return 0.0;
                    }

                    double u_coefficient = 1.0;
                    int    u_exponent    = _r - j;

                    for (int k = _r; k > u_exponent; k--)
                    {
                        u_coefficient *= k;
                    }

                    result = u_coefficient * pow(u, u_exponent);
                }
                break;

            case P_SINE:
                result = 0.0;
                break;
            }
        }

        return result;
    }

    (*@\Green{// getters}@*)
    ECSpace::_OrdinaryBasisFunction::Type ECSpace::_OrdinaryBasisFunction::type() const
    {
        return _type;
    }

    double ECSpace::_OrdinaryBasisFunction::a() const
    {
        return _a;
    }

    double ECSpace::_OrdinaryBasisFunction::b() const
    {
        return _b;
    }

    int ECSpace::_OrdinaryBasisFunction::r() const
    {
        return _r;
    }

    (*@\Green{// Our smart pointers automatically delete the dynamically allocated objects referenced by them when they}@*)
    (*@\Green{// go out of scope. However, it may happen that we need to manually delete the referenced objects and}@*)
    (*@\Green{// to re-initialize the smart pointers to their default null values in order to highlight that some}@*)
    (*@\Green{// mathematical operation could not be performed successfully. For such special cases we}@*)
    (*@\Green{// will use the next private method.}@*)
    void ECSpace::_deleteAllDynamicallyAllocatedObjects()
    {
        _rho                       = SP<RealMatrix>::Default();
        _reversed_v_Wronskian_beta = SP<RealMatrix>::Default();
        _L                         = SP<RealMatrix>::Default();
        _U                         = SP<RealMatrix>::Default();
        _lambda                    = SP< ColumnMatrix<double> >::Default();
        _mu                        = SP<RealMatrix>::Default();
    }

    (*@\Green{// default/special constructor}@*)
    ECSpace::ECSpace(
            double alpha, double beta,
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits):
            _alpha(alpha), _beta(beta)
    {
        if (_alpha >= _beta)
        {
            throw Exception("The starting point of the definition domain should be "
                            "less than its ending point!");
        }

        (*@\Green{// we have to ensure that $z = 0$ is at least a first order root of the characteristic polynomial,}@*)
        (*@\Green{// otherwise the underlying EC vector space of functions would not contain the constants and}@*)
        (*@\Green{// consequently would not possess normalizable bases}@*)
        _polynomial.insertZero(0.0, 0.0, 1);

        updateBothBases(check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits);
    }

    (*@\Green{// Inserts an $m$th order complex root of the form $z = a \pm \mathbf{i}b$ into the factorization of the characteristic polynomial.}@*)
    bool ECSpace::insertZero(
            double a, double b, int m,
            bool update_both_bases,
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits)
    {
        _polynomial.insertZero(a, b, m);

        if (update_both_bases)
        {
            return updateBothBases(check_for_ill_conditioned_matrices,
                                   expected_correct_significant_digits);
        }

        return true;
    }

    bool ECSpace::insertZero(
            const CharacteristicPolynomial::Zero &zero,
            bool update_both_bases,
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits)
    {
        _polynomial.insertZero(zero);

        if (update_both_bases)
        {
            return updateBothBases(check_for_ill_conditioned_matrices,
                                   expected_correct_significant_digits);
        }

        return true;
    }

    (*@\Green{// If exists, deletes a complex root of the form  $z = a \pm \mathbf{i}b$ independently of its order. Note that the}@*)
    (*@\Green{// root $z = 0$ will never be deleted since its existence its critical.}@*)
    bool ECSpace::deleteZero(
            double a, double b,
            bool update_both_bases,
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits)
    {
        if (a == 0.0 && b == 0.0)
        {
            return false;
        }

        if (!_polynomial.deleteZero(a, b))
        {
            return false;
        }
        else
        {
            if (update_both_bases)
            {
                return updateBothBases(check_for_ill_conditioned_matrices,
                                       expected_correct_significant_digits);
            }

            return true;
        }
    }

    bool ECSpace::deleteZero(
            const CharacteristicPolynomial::Zero &zero,
            bool update_both_bases,
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits)
    {
        if (zero.real == 0.0 && zero.imaginary == 0.0)
        {
            return false;
        }

        if (!_polynomial.deleteZero(zero))
        {
            return false;
        }
        else
        {
            if (update_both_bases)
            {
                return updateBothBases(check_for_ill_conditioned_matrices,
                                       expected_correct_significant_digits);
            }

            return true;
        }
    }

    (*@\Green{// Generates/updates all ordinary basis functions $\left\{\varphi_{n,i}\left(u\right) : u \in \left[\alpha, \beta\right]\right\}_{i=0}^{n}$ and also calculates all}@*)
    (*@\Green{// information necessary for the evaluation and differentiation of all normalized B-basis functions}@*)
    (*@\Green{// $\left\{b_{n,i}\left(u\right) : u \in \left[\alpha, \beta\right]\right\}_{i=0}^{n}$.}@*)
    (*@\Green{// The method's implementation is based on the construction process ((\mref{eq:forward_Wronskian})--(\mref{eq:boundary_conditions}), (\mref{eq:LU_factorization_of_reversed_system})--(\mref{eq:construction})).}@*)
    bool ECSpace::updateBothBases((*@\label{ECSpace::updateBothOrdinaryAndNNBBases::start}@*)
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits)
    {
        if (!_polynomial.factorizationChanged())
        {
            return true;
        }

        _is_reflection_invariant = _polynomial.isEvenOrOddFunction();

        _phi.clear();

        for (int z = 0; z < (int)_polynomial._factorization.size(); z++)
        {
            CharacteristicPolynomial::Zero &zero = _polynomial._factorization[z];
            if (zero.order >= 1)
            {
                for (int r = 0; r <= zero.order - 1; r++)
                {
                    _phi.push_back(_OrdinaryBasisFunction(
                            zero.real, zero.imaginary, r, true));

                    if (zero.imaginary != 0.0)
                    {
                        _phi.push_back(_OrdinaryBasisFunction(
                                zero.real, zero.imaginary, r, false));
                    }
                }
            }
        }

        int dimension = (int)_phi.size();
        int n         = dimension - 1;

        if (!dimension)
        {
            _deleteAllDynamicallyAllocatedObjects();
            return false;
        }

        _rho    = SP<RealMatrix>::Default(
                new (nothrow) RealMatrix(dimension, dimension));
        _reversed_v_Wronskian_beta = SP<RealMatrix>::Default(
                new (nothrow) RealMatrix(dimension, dimension));
        _L      = SP<RealMatrix>::Default(
                new (nothrow) RealMatrix(dimension, dimension));
        _U      = SP<RealMatrix>::Default(
                new (nothrow) RealMatrix(dimension, dimension));
        _lambda = SP< ColumnMatrix<double> >::Default(
                new (nothrow) ColumnMatrix<double>(dimension));
        _mu     = SP<RealMatrix>::Default(
                new (nothrow) RealMatrix(
                        dimension,
                        n / (_is_reflection_invariant ? 2 : 1) +1));

        if (!_rho || !_reversed_v_Wronskian_beta || !_L || !_U || !_lambda || !_mu)
        {
            _deleteAllDynamicallyAllocatedObjects();
            return false;
        }

        if (dimension == 1)
        {
            (*_rho)(0, 0) = (*_reversed_v_Wronskian_beta)(0, 0) = (*_L)(0, 0)
                          = (*_U)(0, 0) = (*_lambda)[0] = (*_mu)(0, 0) = 1.0;

            if (_polynomial.factorizationChanged())
            {
                _polynomial.flipFactorizationChangedState();
            }

            return true;
        }

        RealMatrix phi_Wronskian_alpha(dimension, dimension);
        RealMatrix phi_Wronskian_beta(n, dimension);

        #pragma omp parallel for
        for (int j_k = 0; j_k < dimension * dimension; j_k++)
        {
            int j = j_k / dimension;
            int k = j_k % dimension;

            phi_Wronskian_alpha(j, k) = _phi[k](j, _alpha);

            if (j < n)
            {
                phi_Wronskian_beta(j, k)  = _phi[k](j, _beta);
            }
        }

        vector<double> condition_numbers(
                check_for_ill_conditioned_matrices ? dimension : 0);

        bool particular_integrals_aborted = false;

        #pragma omp parallel for
        for (int i = 0; i <= n; i++)
        {
            #pragma omp flush (particular_integrals_aborted)
            if (!particular_integrals_aborted)
            {
                RealMatrix Hermite(dimension, dimension);

                RowMatrix<double> initial_conditions(dimension);
                initial_conditions[i] = 1.0;

                for (int j = 0; j <= i; j++)
                {
                    for (int k = 0; k <= n; k++)
                    {
                        Hermite(j, k) = phi_Wronskian_alpha(j, k);
                    }
                }

                for (int j = 0; j <= n - 1 - i; j++)
                {
                    for (int k = 0; k <= n; k++)
                    {
                        Hermite(i + 1 + j, k) = phi_Wronskian_beta(j, k);
                    }
                }

                if (check_for_ill_conditioned_matrices)
                {
                    SVDecomposition SVD(Hermite);

                    if (!SVD.isCorrect())
                    {
                        particular_integrals_aborted = true;
                        #pragma omp flush (particular_integrals_aborted)
                    }
                    else
                    {
                        condition_numbers[i] = SVD.conditionNumber();
                    }
                }

                PLUDecomposition  PLUD(Hermite);
                RowMatrix<double> particular_integral_coefficients;

                if (!PLUD.solveLinearSystem(
                        initial_conditions, particular_integral_coefficients, false))
                {
                    particular_integrals_aborted = true;
                    #pragma omp flush (particular_integrals_aborted)
                }
                else
                {
                    _rho->setRow(i, particular_integral_coefficients);
                }
            }
        }

        if (particular_integrals_aborted)
        {
            _deleteAllDynamicallyAllocatedObjects();
            return false;
        }

        if (check_for_ill_conditioned_matrices)
        {
            double maximal_condition_number = -std::numeric_limits<double>::max();
            for (int i = 0; i < dimension; i++)
            {
                if (maximal_condition_number < condition_numbers[i])
                {
                    maximal_condition_number = condition_numbers[i];
                }
            }

            double estimated_correct_significant_digits =
                    -log10(2.0 * MACHINE_EPS * maximal_condition_number);

            if (estimated_correct_significant_digits <
                expected_correct_significant_digits)
            {
                _deleteAllDynamicallyAllocatedObjects();
                throw Exception("Based on one of the calculated condition numbers (" +
                                toString(maximal_condition_number) +
                                "), the estimated number (" +
                                toString(estimated_correct_significant_digits) +
                                ") of correct significant digits is less than the "
                                "number (" +
                                toString(expected_correct_significant_digits) +
                                ") of the expected ones! (When a condition number is "
                                "large, the corresponding system is considered ill-"
                                "conditioned and the solution may not be accurate.)\n"
                                 "Try one of the followings:\n"
                                 "  1) lower the number of expected correct significant "
                                 "digits,\n"
                                 "  2) decrease the dimension of the underlying EC "
                                 "space,\n"
                                 "  3) change the endpoints of the definition domain,\n"
                                 "  4) run your code without testing for ill-"
                                 "conditioned matrices and hope for the best (the "
                                 "standard condition number may lead to an overly "
                                 "pessimistic estimate for the overall error).");
            }
        }

        #pragma omp parallel for
        for (int j = 0; j <= n; j++)
        {
            for (int i = n; i >= 0; i--)
            {
                for (int k = 0; k <= n; k++)
                {
                    (*_reversed_v_Wronskian_beta)(j, n - i) += (*_rho)(i, k) *
                                                               _phi[k](j, _beta);
                }
            }
        }

        if (check_for_ill_conditioned_matrices)
        {
            SVDecomposition SVD(*_reversed_v_Wronskian_beta);

            if (!SVD.isCorrect())
            {
                _deleteAllDynamicallyAllocatedObjects();
                return false;
            }

            double condition_number = SVD.conditionNumber();

            double estimated_correct_significant_digits =
                    -log10(2.0 * MACHINE_EPS * condition_number);

            if (estimated_correct_significant_digits <
                expected_correct_significant_digits)
            {
                _deleteAllDynamicallyAllocatedObjects();
                throw Exception("Based on one of the calculated condition numbers (" +
                                toString(condition_number) +
                                "), the estimated number (" +
                                toString(estimated_correct_significant_digits) +
                                ") of correct significant digits is less than the "
                                "number (" +
                                toString(expected_correct_significant_digits) +
                                ") of the expected ones! (When a condition number is "
                                "large, the corresponding system is considered ill-"
                                "conditioned and the solution may not be accurate.)\n"
                                "Try one of the followings:\n"
                                "  1) lower the number of expected correct significant "
                                "digits,\n"
                                "  2) decrease the dimension of the underlying EC "
                                "space,\n"
                                "  3) change the endpoints of the definition domain,\n"
                                "  4) run your code without testing for ill-"
                                "conditioned matrices and hope for the best (the "
                                "standard condition number may lead to an overly "
                                "pessimistic estimate for the overall error).");
            }
        }

        FactorizedUnpivotedLUDecomposition FLUD(*_reversed_v_Wronskian_beta);

        if (!FLUD.isCorrect())
        {
            _deleteAllDynamicallyAllocatedObjects();
            return false;
        }

        ColumnMatrix<double> e0(dimension);
        e0[0] = 1.0;

        Matrix<GLdouble> identity(dimension, n / (_is_reflection_invariant ? 2 : 1) + 1);

        for (int i = 0; i < identity.columnCount(); i++)
        {
            identity(i, i) = 1.0;
        }

        if (!FLUD.solveLLinearSystem(e0, *_lambda) ||
            !FLUD.solveULinearSystem(identity, *_mu))
        {
            _deleteAllDynamicallyAllocatedObjects();
            return false;
        }

        _polynomial.flipFactorizationChangedState();

        return true;
    }(*@\label{ECSpace::updateBothOrdinaryAndNNBBases::end}@*)

    (*@\Green{// Returns the dimension of the underlying EC vector space of functions.}@*)
    int ECSpace::dimension() const
    {
        if (_polynomial.factorizationChanged())
        {
            throw Exception("Before checking the dimension both the ordinary basis and "
                            "the normalized B-basis of the underlying EC space should "
                            "be updated!");
        }

        return (int)_phi.size();
    }

    (*@\Green{// Returns the starting point of the definition domain.}@*)
    double ECSpace::alpha() const
    {
        return _alpha;
    }

    (*@\Green{// Returns the ending point of the definition domain.}@*)
    double ECSpace::beta() const
    {
        return _beta;
    }

    (*@\Green{// Overloaded function operator that evaluates the $j$th order ($j\geq 0$) derivative either of the $i$th ordinary basis }@*)
    (*@\Green{// function or of the normalized B-basis function at the parameter value $u \in \left[\alpha, \beta\right]$.}@*)
    (*@\Green{// The method's implementation is based on Corollary \mref{cor:B_basis_derivatives}.}@*)
    double ECSpace::operator ()(ECSpace::BasisFunctionType type,(*@\label{ECSpace::operator()::start}@*)
                                int i, int j, double u) const
    {
        if (_polynomial.factorizationChanged())
        {
            throw Exception("Before differentiation both the ordinary basis and the "
                            "normalized B-basis of the underlying EC space should be "
                            "updated!");
        }

        int dimension = (int)_phi.size();

        if (!dimension)
        {
            throw Exception("The dimension of the underlying EC space of functions "
                            "should be at least 1!");
        }

        if (i < 0 || i > dimension)
        {
            throw Exception("The function index i should be a non-negative integer that "
                            "is strictly less than the dimension of the underlying EC "
                            "vector space of functions!");
        }

        if (j < 0)
        {
            throw Exception("The differentiation order j should be a non-negative "
                            "integer!");
        }

        if (type == ORDINARY_BASIS)
        {
            return _phi[i](j, u);
        }
        else
        {
            if (!_rho || !_reversed_v_Wronskian_beta || !_L || !_U || !_lambda || !_mu)
            {
                throw Exception("The normalized B-basis of the underlying EC vector "
                                "space of functions was not updated!");
            }

            int n         = dimension - 1;
            int index     = n - i;
            double result = 0.0;

            if (!_is_reflection_invariant || (_is_reflection_invariant && (i > n / 2)))
            {
                #pragma omp parallel for
                for (int r = 0; r <= index; r++)
                {
                    double c_sum = 0.0;
                    for (int c = 0; c <= n; c++)
                    {
                        c_sum += (*_rho)(n - r, c) * _phi[c](j, u);
                    }

                    #pragma omp atomic
                    result += (*_mu)(r, index) * c_sum;
                }

                result *= (*_lambda)[index];
            }
            else
            {
                if (i <= n / 2)
                {
                    u = _alpha + _beta - u;
                    #pragma omp parallel for
                    for (int r = 0; r <= i; r++)
                    {
                        double c_sum = 0.0;
                        for (int c = 0; c <= n; c++)
                        {
                            c_sum += (*_rho)(n - r, c) * _phi[c](j, u);
                        }

                        #pragma omp atomic
                        result += (*_mu)(r, i) * c_sum;
                    }

                    result *= (*_lambda)[i];
                    result *= ((j % 2) ? -1.0 : 1.0);
                }
            }

            return result;
        }
    }(*@\label{ECSpace::operator()::end}@*)

    (*@\Green{// Calculates the entries of the matrix $[t_{i,j}^{n}]_{i=0,j=0}^{n,n}$ that appears in the general basis transformation (\mref{eq:basis_transformation})}@*)
    (*@\Green{// that maps the normalized B-basis $\left\{b_{n,i}\left(u\right) : u \in \left[\alpha, \beta\right]\right\}_{i=0}^{n}$ of the underlying EC vector space of functions}@*)
    (*@\Green{// to its ordinary basis $\left\{\varphi_{n,i}\left(u\right) : u \in \left[\alpha, \beta\right]\right\}_{i=0}^{n}$. The method's implementation is based on Theorem \mref{thm:efficient_basis_transformation}.}@*)
    RealMatrix *ECSpace::basisTransformationFromNBToOrdinary() const (*@\label{ECSpace::basisTransformationFromNNBToOrdinary::start}@*)
    {
        if (_polynomial.factorizationChanged())
        {
            throw Exception("Before the evaluation of the transformation matrix both "
                            "the ordinary basis and the normalized B-basis of the "
                            "underlying EC space should be updated!");
        }

        int dimension = (int)_phi.size();
        int n         = dimension - 1;

        RealMatrix *T = new (nothrow) RealMatrix(dimension, dimension);

        if (!T)
        {
            return nullptr;
        }

        #pragma omp parallel for
        for (int j = 0; j <= n; j++)
        {
            (*T)(0, j) = 1.0;
        }

        #pragma omp parallel for
        for (int i = 1; i <= n; i++)
        {
            (*T)(i, 0) = _phi[i](0, _alpha);
            (*T)(i, n) = _phi[i](0, _beta);
        }

        for (int j = 1; j <= n / 2; j++)
        {
            RowMatrix<double> d_j_b_0_to_j_at_alpha(j + 1);

            #pragma omp parallel for
            for (int k = 0; k <= j; k++)
            {
                d_j_b_0_to_j_at_alpha[k] = operator ()(B_BASIS, k, j, _alpha);
            }

            #pragma omp parallel for
            for (int i = 1; i <= n; i++)
            {
                (*T)(i, j) =  _phi[i](j, _alpha);

                for (int k = 0; k <= j - 1; k++)
                {
                    (*T)(i, j) -= (*T)(i, k) * d_j_b_0_to_j_at_alpha[k];
                }

                (*T)(i, j) /= d_j_b_0_to_j_at_alpha[j];
            }
        }

        for (int j = 1; j <= (n - 1) / 2; j++)
        {
            RowMatrix<double> d_j_b_n_downto_n_minus_j_at_beta(j + 1);

            #pragma omp parallel for
            for (int k = 0; k <= j; k++)
            {
                d_j_b_n_downto_n_minus_j_at_beta[k] = operator ()(B_BASIS,
                                                                  n - k, j, _beta);
            }

            #pragma omp parallel for
            for (int i = 1; i <= n; i++)
            {
                (*T)(i, n - j) =  _phi[i](j, _beta);

                for (int k = 0; k <= j - 1; k++)
                {
                    (*T)(i, n - j) -= (*T)(i, n - k) *
                                      d_j_b_n_downto_n_minus_j_at_beta[k];
                }

                (*T)(i, n - j) /= d_j_b_n_downto_n_minus_j_at_beta[j];
            }
        }

        return T;
    }(*@\label{ECSpace::basisTransformationFromNNBToOrdinary::end}@*)

    (*@\Green{// Determines whether the specified EC vector space is reflection invariant.}@*)
    bool ECSpace::isReflectionInvariant() const
    {
        if (_polynomial.factorizationChanged())
        {
            throw Exception("Before checking the parity of the characteristic "
                            "polynomial both the ordinary basis and the normalized "
                            "B-basis of the underlying EC space should be updated!");
        }

        return _is_reflection_invariant;
    }

    (*@\Green{// If possible, generates the \LaTeX{} expression of the $i$th ordinary basis functions.}@*)
    bool ECSpace::LaTeXExpression(int i, string &expression) const
    {
        if (_polynomial.factorizationChanged() || i >= (int)_phi.size())
        {
            expression.clear();
            return false;
        }

        expression.clear();

        switch (_phi[i].type())
        {
        case _OrdinaryBasisFunction::AET_COSINE:

            if (_phi[i].r())
            {
                if (_phi[i].r() == 1)
                {
                    expression += "u \\cdot ";
                }
                else
                {
                    expression += "u^{";
                    expression += toString(_phi[i].r());
                    expression += "} \\cdot ";
                }
            }

            expression += "e^{";
            if (abs(_phi[i].a()) != 1.0)
            {
                expression += toString(_phi[i].a());
                expression += " \\cdot u} \\cdot \\cos\\left(";
            }
            else
            {
                expression += (_phi[i].a() > 0.0 ? "" : "-");
                expression += "u} \\cdot \\cos\\left(";
            }

            if (abs(_phi[i].b()) != 1.0)
            {
                expression += toString(abs(_phi[i].b()));
                expression += " \\cdot u\\right)";
            }
            else
            {
                expression += "u\\right)";
            }

            break;

        case _OrdinaryBasisFunction::AET_SINE:

            if (_phi[i].b() < 0.0)
            {
                expression += "-";
            }

            if (_phi[i].r())
            {
                if (_phi[i].r() == 1)
                {
                    expression += "u \\cdot ";
                }
                else
                {
                    expression += "u^{";
                    expression += toString(_phi[i].r());
                    expression += "} \\cdot ";
                }
            }

            expression += "e^{";
            if (abs(_phi[i].a()) != 1.0)
            {
                expression += toString(_phi[i].a());
                expression += " \\cdot u} \\cdot \\sin\\left(";
            }
            else
            {
                expression += (_phi[i].a() > 0.0 ? "" : "-");
                expression += "u} \\cdot \\sin\\left(";
            }

            if (abs(_phi[i].b()) != 1.0)
            {
                expression += toString(abs(_phi[i].b()));
                expression += " \\cdot u\\right)";
            }
            else
            {
                expression += "u\\right)";
            }

            break;

        case _OrdinaryBasisFunction::AE_COSINE:

            if (_phi[i].r())
            {
                if (_phi[i].r() == 1)
                {
                    expression += "u \\cdot ";
                }
                else
                {
                    expression += "u^{";
                    expression += toString(_phi[i].r());
                    expression += "} \\cdot ";
                }
            }

            expression += "e^{";
            if (abs(_phi[i].a()) != 1.0)
            {
                expression += toString(_phi[i].a());
                expression += " \\cdot u}";
            }
            else
            {
                expression += (_phi[i].a() > 0.0 ? "" : "-");
                expression += "u}";
            }

            break;

        case _OrdinaryBasisFunction::AE_SINE:
            expression += "0";
            break;

        case _OrdinaryBasisFunction::AT_COSINE:

            if (_phi[i].r())
            {
                if (_phi[i].r() == 1)
                {
                    expression += "u \\cdot ";
                }
                else
                {
                    expression += "u^{";
                    expression += toString(_phi[i].r());
                    expression += "} \\cdot ";
                }
            }

            if (abs(_phi[i].b()) != 1.0)
            {
                expression += "\\cos\\left(";
                expression += toString(abs(_phi[i].b()));
                expression += " \\cdot u\\right)";
            }
            else
            {
                expression += "\\cos\\left(u\\right)";
            }

            break;

        case _OrdinaryBasisFunction::AT_SINE:

            if (_phi[i].b() < 0.0)
            {
                expression += "-";
            }

            if (_phi[i].r())
            {
                if (_phi[i].r() == 1)
                {
                    expression += "u \\cdot ";
                }
                else
                {
                    expression += "u^{";
                    expression += toString(_phi[i].r());
                    expression += "} \\cdot ";
                }
            }

            if (abs(_phi[i].b()) != 1.0)
            {
                expression += "\\sin\\left(";
                expression += toString(abs(_phi[i].b()));
                expression += " \\cdot u\\right)";
            }
            else
            {
                expression += "\\sin\\left(u\\right)";
            }

            break;

        case _OrdinaryBasisFunction::P_COSINE:

            if (_phi[i].r())
            {
                if (_phi[i].r() > 1)
                {
                    expression += "u^{";
                    expression += toString(_phi[i].r());
                    expression += "}";
                }
                else
                {
                    expression += "u";
                }
            }
            else
            {
                expression += "1";
            }

            break;

        case _OrdinaryBasisFunction::P_SINE:

            expression += "0";

            break;
        }

        return true;
    }

    (*@\Green{// Sets the endpoints of the definition domain (in case of modifications the normalized B-basis has to be udated).}@*)
    bool ECSpace::setDefinitionDomain(
            double alpha, double beta,
            bool check_for_ill_conditioned_matrices,
            int expected_correct_significant_digits)
    {
        if ((_alpha != alpha || _beta != beta) && _alpha < beta)
        {
            _alpha = alpha;
            _beta  = beta;
            _polynomial.flipFactorizationChangedState();
        }

        return updateBothBases(check_for_ill_conditioned_matrices,
                               expected_correct_significant_digits);
    }

    (*@\Green{// Returns whether the ordninary basis and the normalized B-basis have to be updated.}@*)
    bool ECSpace::factorizationOfTheCharacteristicPolynomialChanged() const
    {
        return _polynomial.factorizationChanged();
    }

    (*@\Green{// Generates images of all ordinary or normalized B-basis functions.}@*)
    RowMatrix<SP<GenericCurve3>::Default>* ECSpace::generateImagesOfAllBasisFunctions(
            BasisFunctionType type,
            int maximum_order_of_derivatives, int div_point_count) const
    {
        assert("The given maximum order of derivatives should be non-negative!" &&
               maximum_order_of_derivatives >= 0);

        assert("The number of subdivision points should be greater than 1!" &&
               div_point_count > 1);

        if (maximum_order_of_derivatives < 0 || div_point_count <= 1 ||
            _polynomial.factorizationChanged())
        {
            return nullptr;
        }

        int dimension = (int)_phi.size();
        int n         = dimension - 1;

        RowMatrix<SP<GenericCurve3>::Default> *result =
                new (nothrow) RowMatrix<SP<GenericCurve3>::Default>(dimension);

        if (!result)
        {
            return nullptr;
        }

        GLdouble step = (_beta - _alpha) / (div_point_count - 1);

        RowMatrix<GLdouble> u(div_point_count);

        #pragma omp parallel for
        for (int k = 0; k < div_point_count; k++)
        {
            u[k] = min(_alpha + k * step, _beta);
        }

        for (int i = 0; i <= n; i++)
        {
            (*result)[i] = SP<GenericCurve3>::Default(
                    new (nothrow) GenericCurve3(maximum_order_of_derivatives,
                                                div_point_count));

            if (!(*result)[i])
            {
                delete result, result = nullptr;
                return result;
            }

            GenericCurve3 &function_i = *(*result)[i];

            #pragma omp parallel for
            for (int j_k = 0; j_k < (maximum_order_of_derivatives + 1) * div_point_count;
                 j_k++)
            {
                int j = j_k / div_point_count;
                int k = j_k % div_point_count;

                Cartesian3 &derivative = function_i(j, k);

                switch (j) {
                case 0:
                    derivative[0] = u[k];
                    break;
                case 1:
                    derivative[0] = 1.0;
                    break;
                default:
                    derivative[0] = 0.0;
                    break;
                }

                derivative[1] = operator ()(type, i, j, u[k]);
            }
        }

        return result;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    ECSpace* ECSpace::clone() const
    {
        return new ECSpace(*this);
    }

    (*@\Green{// destructor}@*)
    ECSpace::~ECSpace()
    {
        _deleteAllDynamicallyAllocatedObjects();
    }


    (*@\Green{// overloaded output to stream operator}@*)
    ostream& operator <<(ostream &lhs, const ECSpace &rhs)
    {
        return lhs << rhs._polynomial << endl;
    }
}

