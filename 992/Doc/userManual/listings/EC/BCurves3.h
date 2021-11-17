//----------------------------------------------------------------------------------
// File:        EC/BCurves3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef BCURVES3_H
#define BCURVES3_H

#include "../Core/SmartPointers/SpecializedSmartPointers.h"
#include "../Core/Geometry/Curves/LinearCombinations3.h"
#include "ECSpaces.h"

namespace cagd
{
    (*@\Green{// The class BCurve3 below represents a B-curve of type (\mref{eq:B_curve}). It is derived from the abstract base class}@*)
    (*@\Green{// LinearCombination3 that provides two pure virtual methods which have to be redeclared and defined in}@*)
    (*@\Green{// this derived class, otherwise one cannot instantiate curve objects from it. One of these methods is}@*)
    (*@\Green{//}@*)
    (*@\Green{// \hspace{0.5cm}GLboolean blendingFunctionValues(GLdouble u, RowMatrix$<$GLdouble$>$ \&values) const}@*)
    (*@\Green{//}@*)
    (*@\Green{// that has to evaluate all normalized B-basis functions of the underlying EC space represented by the variable \_S.}@*)
    (*@\Green{// This method is required by the base class LinearCombination3 when it tries to solve user-defined curve}@*)
    (*@\Green{// interpolation problems by calling its method}@*)
    (*@\Green{//}@*)
    (*@\Green{// \hspace{0.5cm}virtual GLboolean updateDataForInterpolation(const ColumnMatrix$<$GLdouble$>$\& knot\_vector,}@*)
    (*@\Green{// \hspace{5.3cm}const ColumnMatrix$<$Cartesian3$>$\& data\_points\_to\_interpolate).}@*)
    (*@\Green{//}@*)
    (*@\Green{// The other pure virtual method that has to be redeclared and implemented is}@*)
    (*@\Green{//}@*)
    (*@\Green{// \hspace{0.5cm}GLboolean calculateDerivatives(GLuint max\_order\_of\_derivatives, GLdouble u, Derivatives \&d) const}@*)
    (*@\Green{//}@*)
    (*@\Green{// that is used by the base class in order to generate the image of the B-curve by invoking its method}@*)
    (*@\Green{//}@*)
    (*@\Green{// \hspace{0.5cm}virtual GenericCurve3* generateImage(GLuint max\_order\_of\_derivatives, GLuint div\_point\_count,}@*)
    (*@\Green{// \hspace{4.5cm}GLenum usage\_flag = GL\_STATIC\_DRAW) const.}@*)
    (*@\Green{//}@*)
    (*@\Green{// Apart from the aforementioned two pure virtual methods, from image generation and solution of curve}@*)
    (*@\Green{// interpolation problems, the base class LinearCombination3 also provides other fully implemented}@*)
    (*@\Green{// functionalities that are responsible for the control point and vertex buffer object management of}@*)
    (*@\Green{// the control polygon:}@*)
    (*@\Green{//}@*)
    (*@\Green{// \hspace{0.5cm}const Cartesian3\&\ \ \ \,operator [](const GLint \&index) const; // get control point by constant reference}@*)
    (*@\Green{// \hspace{0.5cm}Cartesian3\&\ \ \ \ \ \ \ \ \ \ \ \,operator [](const GLint \&index); \ \ \ \ \ \ \ \ // get control point by non-constant reference}@*)
    (*@\Green{//}@*)
    (*@\Green{// \hspace{0.5cm}virtual GLvoid deleteVertexBufferObjectsOfData();}@*)
    (*@\Green{// \hspace{0.5cm}virtual GLboolean renderData(GLenum render\_mode = GL\_LINE\_STRIP) const;}@*)
    (*@\Green{// \hspace{0.5cm}virtual GLboolean updateVertexBufferObjectsOfData(GLenum usage\_flag = GL\_STATIC\_DRAW).}@*)
    class BCurve3: public LinearCombination3
    {
    protected:
        (*@\Green{// represents the EC space $\mathbb{S}_n^{\alpha,\beta}=\left\langle\mathcal{F}_{n}^{\alpha,\beta}:=\{\varphi_{n,i}\left(u\right):u\in\left[\alpha,\beta\right]\}_{i=0}^n\right\rangle = \left\langle\mathcal{B}_{n}^{\alpha,\beta}:=\{b_{n,i}\left(u\right):u\in\left[\alpha,\beta\right]\}_{i=0}^n\right\rangle$,}@*)
        (*@\Green{// where $0<\beta-\alpha<\ell^{\prime}\left(\mathbb{S}_{n}^{\alpha,\beta}\right)$}@*)
        ECSpace _S;

        (*@\Green{// a smart pointer to the transformation matrix $[t_{i,j}^{n}]_{i=0,\,j=0}^{n,\,n}$ that maps the normalized B-basis $\mathcal{B}_{n}^{\alpha,\beta}$ to the}@*)
        (*@\Green{// ordinary basis $\mathcal{F}_{n}^{\alpha,\beta}$}@*)
        SP<RealMatrix>::Default _T;

    public:
        (*@\Green{// special constructor}@*)
        BCurve3(const ECSpace &S, GLenum data_usage_flag = GL_STATIC_DRAW);

        (*@\Green{// inherited pure virtual methods that have to be redeclared and defined}@*)
        GLboolean blendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const;
        GLboolean calculateDerivatives(GLint maximum_order_of_derivatives, GLdouble u,
                                       Derivatives &d) const;

        (*@\Green{// Overloaded function operator that evaluates the derivative $b_{n,i}^{\left(j\right)}\left(u\right)$, where $u \in \left[\alpha, \beta\right]$, $i=0,1,\ldots,n$}@*)
        (*@\Green{// and $j \in \mathbb{N}$.}@*)
        GLdouble operator ()(GLint i, GLint j, GLdouble u) const;

        (*@\Green{// Overloaded member functions for general order elevation. In order to generate a higher dimensional EC}@*)
        (*@\Green{// space $\mathbb{S}_{n+q}^{\alpha,\beta}$ that fulfills the condition $\mathbb{S}_n^{\alpha,\beta}\subset\mathbb{S}_{n+q}^{\alpha,\beta}$, where $q\geq 1$ and $0<\beta-\alpha<\min\left\{\ell^{\prime}\left(\mathbb{S}_n^{\alpha,\beta}\right), \ell^{\prime}\left(\mathbb{S}_{n+q}^{\alpha,\beta}\right)\right\}$,}@*)
        (*@\Green{// one has either to define a new complex root $z=a\pm\mathbf{i}b$ of multiplicity $m \geq 1$, or to specify an existing}@*)
        (*@\Green{// root with a multiplicity that is greater than its former value $m'\geq 1$. If $z \in \mathbb{C}\setminus\mathbb{R}$, then $q=2\left(m-m'\right)$,}@*)
        (*@\Green{// otherwise $q=m-m'$, where $m'=0$ whenever $z$ denotes a nonexisting former root.}@*)
        (*@\Green{// The implementation of the next two methods are based on Lemma \mref{lem:general_order_elevation}.}@*)
        BCurve3* performOrderElevation(
                GLdouble a, GLdouble b, GLint multiplicity,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5) const;

        BCurve3* performOrderElevation(
                const CharacteristicPolynomial::Zero &zero,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5) const;

        (*@\Green{// Using the general B-algorithm presented in Theorem \mref{thm:general_subdivision}, the method subdivides the given B-curve}@*)
        (*@\Green{// into two arcs of the same order at the parameter value $\gamma \in \left(\alpha,\beta\right)$. In case of success, the output is a}@*)
        (*@\Green{// nonzero raw pointer to a $2$-element row matrix that stores $2$ smart pointers (based on the deep copy}@*)
        (*@\Green{// ownership policy) that point to the left and right arcs of the subdivided curve.}@*)
        virtual RowMatrix<SP<BCurve3>::Default>* performSubdivision(
                GLdouble gamma,
                bool check_for_ill_conditioned_matrices = false,
                GLint expected_correct_significant_digits = 5) const;

        (*@\Green{// Inherited virtual method that has to be redeclared and redefined, since if one alters the endpoints of}@*)
        (*@\Green{// the definition domain $\left[\alpha, \beta\right]$ both bases $\mathcal{F}_{n}^{\alpha,\beta}$ and $\mathcal{B}_{n}^{\alpha,\beta}$ of $\mathbb{S}_n^{\alpha,\beta}$ have to be updated. Note that, the new}@*)
        (*@\Green{// interval length should also satisfy the condition $0<\beta-\alpha<\ell^{\prime}\left(\mathbb{S}_{n}^{\alpha,\beta}\right).$}@*)
        GLboolean setDefinitionDomain(GLdouble alpha, GLdouble beta,
                                      bool check_for_ill_conditioned_matrices = false,
                                      GLint expected_correct_significant_digits = 5);

        (*@\Green{// Given an ordinary integral curve $\mathbf{c}\left(  u\right)  =\sum_{i=0}^{n}\boldsymbol{\lambda}_{i}\varphi_{n,i}\left(  u\right)  ,~u\in\left[  \alpha,\beta\right],~0<\beta-\alpha<\ell^{\prime}\left(\mathbb{S}_{n}^{\alpha,\beta}\right),~\boldsymbol{\lambda}_{i}\in\mathbb{R}^{3}$,}@*)
        (*@\Green{// the method determines the coefficients $[\mathbf{p}_j]_{j = 0}^n$ of the normalized B-basis functions $\left\{b_{n,j}\left(u\right) : u \in \left[\alpha, \beta\right]\right\}_{j=0}^{n}$}@*)
        (*@\Green{// such that $\mathbf{c}\left(u\right) \equiv \sum_{j = 0}^n \mathbf{p}_j b_{n,j}\left(u\right)$, $\forall u \in \left[\alpha, \beta\right]$. Note that, the control points $[\mathbf{p}_j]_{j = 0}^n$ that have to be}@*)
        (*@\Green{// updated are stored in the inherited data structure ColumnMatrix$<$Cartesian3$>$ LinearCombination3::\_data.}@*)
        (*@\Green{// The method's implementation is based on Theorems \mref{thm:efficient_basis_transformation} and \mref{thm:integral_curves}.}@*)
        GLboolean updateControlPointsForExactDescription(
                const RowMatrix<Cartesian3> &lambda);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual BCurve3* clone() const;
    };
}

#endif (*@\Green{// BCURVES3\_H}@*)
