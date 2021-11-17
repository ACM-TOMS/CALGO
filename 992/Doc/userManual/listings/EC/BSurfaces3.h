//----------------------------------------------------------------------------------
// File:        EC/BSurfaces3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef BSURFACES3_H
#define BSURFACES3_H

#include "../Core/SmartPointers/SpecializedSmartPointers.h"
#include "../Core/Geometry/Surfaces/TensorProductSurfaces3.h"
#include "../Core/Math/Constants.h"
#include "../Core/Math/Matrices.h"
#include "ECSpaces.h"
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

namespace cagd
{
    (*@\Green{// Let $\mathcal{F}_{n_{r}}^{\alpha_{r},\beta_{r}}=\left\{  \varphi_{n_{r},i_{r}}\left(  u_{r}\right)  :u_{r}\in\left[\alpha_{r},\beta_{r}\right]  \right\}  _{i_{r}=0}^{n_{r}},~\varphi_{n_{r},0} \equiv 1,~0<\beta_r - \alpha_r < \ell^{\prime}\left(\mathbb{S}_{n_r}^{\alpha_r, \beta_r}\right)$ be the ordinary basis}\label{src:BSurfaces3.h:OrdinarySurfaceCoefficients:start}@*)
    (*@\Green{// and $\mathcal{B}_{n_{r}}^{\alpha_{r},\beta_{r}}=\left\{  b_{n_{r},j_{r}}\left(u_{r}\right)  :u_{r}\in\left[  \alpha_{r},\beta_{r}\right]  \right\}_{j_{r}=0}^{n_{r}}$ be the normalized B-basis of some EC vector space $\mathbb{S}_{n_{r}}^{\alpha_{r},\beta_{r}}$ of}@*)
    (*@\Green{// functions and denote by $[  t_{i_{r},j_{r}}^{n_{r}}]  _{i_{r}=0,~j_{r}=0}^{n_{r},~n_{r}}$ the regular square matrix that transforms $\mathcal{B}_{n_{r}}^{\alpha_{r},\beta_{r}}$ to $\mathcal{F}_{n_{r}}^{\alpha_{r},\beta_{r}}$, where}@*)
    (*@\Green{// $r=0,1$. Consider also the ordinary integral surface}@*)
    (*@\Green{// $\mathbf{s}\left(  u_0, u_1\right)  =\left[\begin{array}[c]{ccc}s^{0}\left(  u_0, u_1\right)   & s^{1}\left(  u_0, u_1\right)   & s^{2}\left(  u_0, u_1\right)\end{array}\right]  ^{T}\in\mathbb{R}^{3},~\left(u_0, u_1\right)\in\left[  \alpha_{0},\beta_{0}\right]  \times\left[  \alpha_{1},\beta_{1}\right],$}@*)
    (*@\Green{// where $s^{k}\left(  u_0, u_1\right)  =\sum_{\zeta=0}^{\sigma_{k}-1}\prod_{r=0}^{1}\left(  \sum_{i_{r}=0}^{n_{r}}\lambda_{n_r,i_{r}}^{k,\zeta}\varphi_{n_{r},i_{r}}\left(  u_{r}\right)  \right)  ,~\sigma_{k}\geq 1, ~k=0,1,2$.}@*)

    (*@\Green{// The class OrdinarySurfaceCoefficients stores the pairs $\left\{\left(\left[\lambda_{n_0,i_0}^{k,\zeta}\right]_{i_0=0}^{n_0}, \left[\lambda_{n_1,i_1}^{k,\zeta}\right]_{i_1=0}^{n_1}\right)\right\}_{\zeta = 0,\,k=0}^{\sigma_{k}-1,\,2}$ of row matrices}@*)
    (*@\Green{// that consist of those coefficients which appear in the multiplied linear combinations of the coordinate functions.}@*)
    class OrdinarySurfaceCoefficients
    {
    protected:
        (*@\Green{$// \hspace{1.3cm}k$\hspace{2.1cm}$\zeta$}@*)
        std::vector< std::vector< pair<RowMatrix<GLdouble>,
                                       RowMatrix<GLdouble> > > > _lambda;

    public:
        (*@\Green{// special constructor\hspace{3.05cm}$n_0+1$,\hspace{2.3cm}$n_1+1$}@*)
        OrdinarySurfaceCoefficients(GLint u_dimension, GLint v_dimension,
                                    const std::vector<GLint> &sigma); (*@\Green{// $\left[\sigma_{k}\right]_{k=0}^2$}@*)

        (*@\Green{// based on the given variable type, returns the constant reference of either $\lambda_{n_0,i_0}^{k,\zeta}$ or $\lambda_{n_1,i_1}^{k,\zeta}$, where the}@*)
        (*@\Green{// input variable index denotes either $i_0$ or $i_1$}@*)
        const GLdouble& operator()(GLint k, GLint zeta, variable::Type type,
                                   GLint index) const;

        (*@\Green{// based on the given variable type, returns non-constant references to either $\lambda_{n_0,i_0}^{k,\zeta}$ or $\lambda_{n_1,i_1}^{k,\zeta}$, where the}@*)
        (*@\Green{// input variable index denotes either $i_0$ or $i_1$}@*)
        GLdouble& operator()(GLint k, GLint zeta, variable::Type type, GLint index);

        (*@\Green{// returns the number of multiplied linear combination pairs appearing in the coordinate function $s^{k}$}@*)
        GLint sigma(GLint k) const;

        (*@\Green{// based on the given variable type, returns the value either of $n_0+1$ or $n_1+1$}@*)
        GLint dimension(variable::Type type) const;

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual OrdinarySurfaceCoefficients* clone() const;
    }; (*@\label{src:BSurfaces3.h:OrdinarySurfaceCoefficients:end}@*)

    (*@\Green{// The class BSurface3 represents tensor product of surfaces of type (\mref{eq:B-surface}), i.e.,}@*)
    (*@\Green{// $\mathbf{s}_{n_0,n_1}\left(u_0,u_1\right)=\sum_{i_0=0}^{n_0} \sum_{i_1 = 0}^{n_1}\mathbf{p}_{i_0,i_1}b_{n_0,i_0}\left(u_0;\alpha_0, \beta_0\right)b_{n_1,i_1}\left(u_1;\alpha_1, \beta_1\right),\:\mathbf{p}_{i_0,i_1}=\left[p_{i_0,i_1}^{k}\right]_{k=0}^2\in\mathbb{R}^3$}@*)
    class BSurface3: public TensorProductSurface3
    {
    protected:
        ECSpace                 _S[2]; (*@\Green{// EC spaces $\mathbb{S}_{n_r}^{\alpha_r,\beta_r}$, $r=0,1$}@*)
        SP<RealMatrix>::Default _T[2]; (*@\Green{// transformation matrices $[  t_{i_{r},j_{r}}^{n_{r}}]_{i_{r}=0,~j_{r}=0}^{n_{r},~n_{r}}$, $r=0,1$}@*)

    public:
        (*@\Green{// special constructor}@*)
        BSurface3(const ECSpace &u_space, const ECSpace &v_space);

        (*@\Green{// Inherited virtual method that has to be redeclared and redefined, since if one alters the endpoints of}@*)
        (*@\Green{// the definition domain $\left[\alpha_r, \beta_r\right]$ both bases $\mathcal{F}_{n_r}^{\alpha_r,\beta_r}$ and $\mathcal{B}_{n_r}^{\alpha_r,\beta_r}$ of $\mathbb{S}_{n_r}^{\alpha_r,\beta_r}$ have to be updated. Note that,}@*)
        (*@\Green{// the new interval length should also satisfy the condition $0 < \beta_r - \alpha_r < \ell^{\prime}\left(\mathbb{S}_{n_r}^{\alpha_r,\beta_r}\right)$.}@*)
        GLboolean setInterval(
                variable::Type type, GLdouble alpha, GLdouble beta,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5);
        
        (*@\Green{// inherited pure virtual methods that have to be redeclared and defined}@*)
        GLboolean blendingFunctionValues(
                variable::Type type, GLdouble parameter_value,
                RowMatrix<GLdouble> &values) const;

        GLboolean calculateAllPartialDerivatives(
                GLint maximum_order_of_partial_derivatives,
                GLdouble u, GLdouble v, PartialDerivatives& pd) const;

        GLboolean calculateDirectionalDerivatives(
                variable::Type direction,
                GLint maximum_order_of_directional_derivatives, GLdouble u, GLdouble v,
                DirectionalDerivatives &d) const;

        (*@\Green{// const and non-const indexing operators of the control points}@*)
        using TensorProductSurface3::operator ();

        (*@\Green{// Overloaded function operator that evaluates the derivative $b_{n_r,i_r}^{\left(j\right)}\left(u_r\right)$, where $u_r \in \left[\alpha_r, \beta_r\right]$, $i_r=0,1,\ldots,n_r$,}@*)
        (*@\Green{// and $j \in \mathbb{N}$.}@*)
        GLdouble operator ()(
                variable::Type type, GLint function_index,
                GLint differentiation_order, GLdouble parameter_value) const;

        (*@\Green{// Overloaded member functions for general order elevation. In order to generate a higher dimensional EC}@*)
        (*@\Green{// space $\mathbb{S}_{n_r+q}^{\alpha_r,\beta_r}$ that fulfills the condition $\mathbb{S}_{n_r}^{\alpha_r,\beta_r}\subset\mathbb{S}_{n_r+q}^{\alpha_r,\beta_r}$ where $q\geq 1$ and $0<\beta_r-\alpha_r < \min\left\{\ell^{\prime}\left(\mathbb{S}_{n_r}^{\alpha_r,\beta_r}\right),\right.$}@*)
        (*@\Green{// $\left.\ell^{\prime}\left(\mathbb{S}_{n_r+q}^{\alpha_r,\beta_r}\right)\right\}$ one has either to define a new complex root $z=a\pm\mathbf{i}b$ of multiplicity $m \geq 1$, or to specify}@*)
        (*@\Green{// an existing root with a multiplicity that is greater than its former value $m'\geq 1$. If $z \in \mathbb{C}\setminus\mathbb{R}$, then}@*)
        (*@\Green{// $q=2\left(m-m'\right)$, otherwise $q=m-m'$, where $m'=0$ whenever $z$ denotes a nonexisting former root.}@*)
        (*@\Green{// Their implementations are based on the extension of Lemma \mref{lem:general_order_elevation}.}@*)
        BSurface3* performOrderElevation(
                variable::Type type,
                GLdouble a, GLdouble b, GLint multiplicity,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5) const;
        
        BSurface3* performOrderElevation(
                variable::Type type,
                const CharacteristicPolynomial::Zero zero,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5) const;
        
        (*@\Green{// Using the extension of the general B-algorithm formulated in Theorem \mref{thm:general_subdivision}, the method subdivides the}@*)
        (*@\Green{// given B-surface along the specified direction into two patches of the same order at the parameter value}@*)
        (*@\Green{// $\gamma \in \left(\alpha_r,\beta_r\right)$. In case of success, the output is a nonzero raw pointer to a $2$-element row matrix that stores}@*)
        (*@\Green{// $2$ smart pointers (based on the deep copy ownership policy) that point to the left and right patches of the}@*)
        (*@\Green{// subdivided surface.}@*)
        RowMatrix<SP<BSurface3>::Default>* performSubdivision(
                variable::Type type,
                GLdouble gamma,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5) const;

        (*@\Green{// Given an ordinary integral surface of type (\mref{eq:ordinary_integral_surface}) the method updates the control points $[\mathbf{p}_{i_0,i_1}]_{i_0 = 0,\,i_1=0}^{n_0,\,n_1}$}@*)
        (*@\Green{// of the B-surface (\mref{eq:B-surface}) in order to ensure control point based exact description. The ordinary surface}@*)
        (*@\Green{// (\mref{eq:ordinary_integral_surface}) is described by an instance of the class OrdinarySurfaceCoefficients, while the control points}@*)
        (*@\Green{// $[\mathbf{p}_{i_0,i_1}]_{i_0 = 0,\,i_1=0}^{n_0,\,n_1}$ that have to be updated are stored in the inherited data structure}@*)
        (*@\Green{// Matrix$<$Cartesian3$>$ TensorProductSurface3::\_data. The method is based on formulas of Theorem \mref{thm:integral_surfaces}.}@*)
        GLboolean updateControlPointsForExactDescription(
                const OrdinarySurfaceCoefficients &lambda);

        (*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual BSurface3* clone() const;
    };
}

#endif (*@\Green{// BSURFACES3\_H}@*)
