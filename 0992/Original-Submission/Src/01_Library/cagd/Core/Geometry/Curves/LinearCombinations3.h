//----------------------------------------------------------------------------------
// File:        Core/Geometry/Curves/LinearCombinations3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef LINEARCOMBINATIONS3_H
#define LINEARCOMBINATIONS3_H

#include "../Coordinates/Cartesians3.h"
#include "GenericCurves3.h"
#include "../../Math/Matrices.h"
#include "../../Shaders/ShaderPrograms.h"

namespace cagd
{
    //(*@\Green{// Represents the linear combination $\mathbf{c}\left(u\right) = \sum_{i=0}^n \mathbf{p}_{i} b_{n,i}\left(u; u_{\min}, u_{\max}\right)$, where the vectors $\left[\mathbf{p}_{i}\right]_{i=0}^{n} \in \mathcal{M}_{1,n+1}\left(\mathbb{R}^3\right)$}@*)
    //(*@\Green{// usually form a control polygon of points (but they may also denote other user-defined information vectors}@*)
    //(*@\Green{// like first and higher order derivatives as in the case of polynomial Hermite curves), while the function system}@*)
    //(*@\Green{// $\left\{b_{n,i}\left(u; u_{\min}, u_{\max}\right) : u \in \left[u_{\min}, u_{\max}\right]\right\}_{i=0}^n$ forms the basis of some not necessarily EC vector space}@*)
    //(*@\Green{// of functions.}@*)
    //(*@\Green{// If vectors $\left[\mathbf{p}_{i}\right]_{i=0}^{n}$ do not describe a control polygon, then the virtual vertex buffer object handling methods}@*)
    //(*@\Green{// deleteVertexBufferObjectsOfData, renderData and updateVertexBufferObjectsOfData should be redeclared}@*)
    //(*@\Green{// and redefined in derived classes.}@*)
    //(*@\Green{// Virtual methods generateImage and updateDataForInterpolation should be redeclared and redefined in}@*)
    //(*@\Green{// derived classes only when the user knows more efficient techniques for these tasks.}@*)
    //(*@\Green{// Pure virtual methods blendingFunctionValues and  calculateDerivatives and should be obligatorily redeclared}@*)
    //(*@\Green{// and defined in derived classes, otherwise the user cannot instantiate objects from them.}@*)
    class LinearCombination3
    {
    public:
        //(*@\Green{// public nested class that represents a curve point and its associated first and higher order derivatives}@*)
        class Derivatives: public ColumnMatrix<Cartesian3>
        {
        public:
            //(*@\Green{// default/special constructor}@*)
            Derivatives(GLint maximum_order_of_derivatives = 2);

            //(*@\Green{// when called, all inherited Cartesian coordinates are set to the origin}@*)
            GLvoid loadNullVectors();

            //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
            Derivatives* clone() const;
        };

    protected:
        GLuint                   _vbo_data;       //(*@\Green{// vertex buffer object of the user-specified}@*)
                                                  //(*@\Green{// data that usually represents a control polygon}@*)
        GLenum                   _data_usage_flag;//(*@\Green{// usage flag of the previous vertex buffer object}@*)
        GLdouble                 _u_min, _u_max;  //(*@\Green{// endpoints of the definition domain}@*)
        ColumnMatrix<Cartesian3> _data;           //(*@\Green{// vectors (usually representing control points)}@*)
                                                  //(*@\Green{// that have to be blended}@*)

    public:
        //(*@\Green{// special constructor}@*)
        LinearCombination3(GLdouble u_min, GLdouble u_max,
                           GLint data_count, GLenum data_usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// copy constructor}@*)
        LinearCombination3(const LinearCombination3& lc);

        //(*@\Green{// assignment operator}@*)
        LinearCombination3& operator =(const LinearCombination3& rhs);

        //(*@\Green{// By default, it is assumed that the vectors stored in ColumnMatrix$<$Cartesian3$>$ \_data}@*)
        //(*@\Green{// represent the vertices of a control polygon. If this is not the case, the next four virtual}@*)
        //(*@\Green{// VBO handling methods can be redeclared and defined in derived classes.}@*)

        //(*@\Green{// Deletes the vertex buffer objects of the control polygon.}@*)
        virtual GLvoid    deleteVertexBufferObjectsOfData();

        //(*@\Green{// Control polygon rendering methods.}@*)

        //(*@\Green{// The next rendering method will return GL\_FALSE if:\label{src:LinearCombination3:renderData:start}}@*)
        //(*@\Green{// \ \ - the given shader program is not active;}@*)
        //(*@\Green{// \ \ - the user-defined position attribute name cannot be found in the list of active attributes}@*)
        //(*@\Green{// \ \ \ \ of the provided shader program, or exists but is of incorrect type;}@*)
        //(*@\Green{// \ \ - the render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP, GL\_LINE\_LOOP}@*)
        //(*@\Green{// \ \ \ \ and GL\_POINTS).}@*)
        virtual GLboolean renderData(const ShaderProgram &program,
                                     GLenum render_mode = GL_LINE_STRIP) const;

        //(*@\Green{// If during rendering one intends to use shader program objects that are not instances of}@*)
        //(*@\Green{// our class ShaderProgram, then one should specify an attribute location associated with}@*)
        //(*@\Green{// positions of type vec3.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
        //(*@\Green{// \ \ - there is no active shader program object;}@*)
        //(*@\Green{// \ \ - the given position location either cannot be found in the list of active attribute locations,}@*)
        //(*@\Green{// \ \ \ \ or exists but is of incorrect type;}@*)
        //(*@\Green{// \ \ - the render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP, GL\_LINE\_LOOP}@*)
        //(*@\Green{// \ \ \ \ and GL\_POINTS).}@*)
        virtual GLboolean renderData(GLenum render_mode = GL_LINE_STRIP,
                                     GLint vec3_position_location = 0) const; //(*@\label{src:LinearCombination3:renderData:end}@*)

        //(*@\Green{// Updates the vertex buffer objects of the control polygon.}@*)
        virtual GLboolean updateVertexBufferObjectsOfData(
                GLenum usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// get data by constant reference}@*)
        const Cartesian3& operator [](const GLint &index) const;

        //(*@\Green{// get data by non-constant reference}@*)
        Cartesian3& operator [](const GLint &index);

        //(*@\Green{// sets/returns the endpoints of the definition domain}@*)
        virtual GLboolean setDefinitionDomain(const GLdouble &u_min,
                                              const GLdouble &u_max);
        GLvoid definitionDomain(GLdouble& u_min, GLdouble& u_max) const;

        //(*@\Green{// get number of data}@*)
        GLint dataCount() const;

        //(*@\Green{// abstract method that should calculate a row matrix consisting of the blending function values}@*)
        //(*@\Green{// $\left\{b_{n,i}(u) : u \in \left[u_{\min}, u_{\max}\right]\right\}_{i=0}^{n}$, where $n+1$ denotes the number of the user-specified data}@*)
        virtual GLboolean blendingFunctionValues(GLdouble u,
                                                 RowMatrix<GLdouble>& values) const = 0;

        //(*@\Green{// abstract method that should calculate the point and its associated higher order order derivatives of the}@*)
        //(*@\Green{// linear combination $\mathbf{c}\left(u\right) = \sum_{i=0}^n \mathbf{p}_i b_{n,i}\left(u\right)$ at the parameter value $u \in \left[u_{\min}, u_{\max}\right]$, where $\left[\mathbf{p}_i\right]_{i=0}^n$}@*)
        //(*@\Green{// denotes the user-specified data}@*)
        virtual GLboolean calculateDerivatives(GLint maximum_order_of_derivatives,
                                               GLdouble u, Derivatives& d) const = 0;

        //(*@\Green{// generates the image of the given linear combination}@*)
        virtual GenericCurve3* generateImage(
                GLint maximum_order_of_derivatives,
                GLint div_point_count, GLenum usage_flag = GL_STATIC_DRAW) const;

        //(*@\Green{// Solves the user-specified curve interpolation problem $\mathbf{c}\left(u_j\right) = \mathbf{d}_j \in \mathbb{R}^3,~j=0,1,\ldots,n$, where $\left[u_j\right]_{j=0}^n$}@*)
        //(*@\Green{// forms a knot vector consisting of strictly increasing subdivision points of the definition domain $[u_{\min}, u_{\max}]$,}@*)
        //(*@\Green{// while $\left[\mathbf{d}_j\right]_{j=0}^n$ denotes the data points that have to be interpolated.}@*)
        //(*@\Green{// In case of success, the solution $\left[\mathbf{p}_i\right]_{i=0}^n$ will be stored in the column matrix \_data.}@*)
        virtual GLboolean updateDataForInterpolation(
                const ColumnMatrix<GLdouble>& knot_vector,
                const ColumnMatrix<Cartesian3>& data_points_to_interpolate);

        //(*@\Green{// destructor}@*)
        virtual ~LinearCombination3();
    };
}

#endif //(*@\Green{// LINEARCOMBINATIONS3\_H}@*)
