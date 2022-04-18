//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/TensorProductSurfaces3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef TENSORPRODUCTSURFACES3_H
#define TENSORPRODUCTSURFACES3_H

#include <GL/glew.h>

#include "../Coordinates/Cartesians3.h"
#include "../../Math/Constants.h"
#include "../../Math/Matrices.h"
#include "../../Geometry/Curves/GenericCurves3.h"
#include "../../SmartPointers/SpecializedSmartPointers.h"
#include "../../Shaders/ShaderPrograms.h"
#include "TriangleMeshes3.h"

#include <iostream>
#include <vector>

namespace cagd
{
    //(*@\Green{// Represents the tensor product surface $\mathbf{s}\left(u,v\right) = \sum_{i=0}^n \sum_{j=0}^m \mathbf{p}_{i,j} b_{n,i}\left(u; u_{\min}, u_{\max}\right) b_{m,j}\left(v; v_{\min}, v_{\max}\right)$,}@*)
    //(*@\Green{// where the vectors $\left[\mathbf{p}_{i,j}\right]_{i=0,\,j=0}^{n,m} \in \mathcal{M}_{n+1,m+1}\left(\mathbb{R}^3\right)$ usually form a control net of points (but they may also}@*)
    //(*@\Green{// denote other user-defined information vectors like first and higher order partial derivatives like in the case}@*)
    //(*@\Green{// of polynomial Hermite patches), while function systems $\left\{b_{n,i}\left(u; u_{\min}, u_{\max}\right) : u \in \left[u_{\min}, u_{\max}\right]\right\}_{i=0}^n$ and}@*)
    //(*@\Green{// $\left\{b_{m,j}\left(v; v_{\min}, v_{\max}\right) : v \in \left[v_{\min}, v_{\max}\right]\right\}_{j=0}^m$ form the bases of some not necessarily EC vector spaces of functions.}@*)
    //(*@\Green{// If vectors $\left[\mathbf{p}_{i,j}\right]_{i=0,\,j=0}^{n,m}$ do not form a control net, then the vertex buffer object handling virtual methods}@*)
    //(*@\Green{// deleteVertexBufferObjectsOfData, renderData and updateVertexBufferObjectsOfData should be redeclared}@*)
    //(*@\Green{// and redefined in derived classes.}@*)
    //(*@\Green{// Virtual methods generateImage and updateDataForInterpolation should be redeclared and redefined in}@*)
    //(*@\Green{// derived classes only when the user knows more efficient techniques for these tasks.}@*)
    //(*@\Green{// Pure virtual methods blendingFunctionValues,  calculateAllPartialDerivatives and calculateDirectional-}@*)
    //(*@\Green{// Derivatives should be obligatorily redeclared and defined in derived classes, otherwise the user cannot}@*)
    //(*@\Green{// instantiate objects from them.}@*)
    class TensorProductSurface3
    {
    public:
        //(*@\Green{// When generating the image of the underlying tensor product surface, the user is able to choose different}\label{src:TensorProductSurface3.h:ImageColorScheme:start}@*)
        //(*@\Green{// color schemes based on the following predefined pointwise zeroth, first and second order energy variations.}@*)
        //(*@\Green{// In each case the obtained color variation is determined by using the utility function coldToHotColorMap}@*)
        //(*@\Green{// presented in Listing \mref{src:Utilities.cpp} that generates a rainbow-like color map based on the minimum and maximum}@*)
        //(*@\Green{// value of the applied energy.}@*)
        //(*@\Green{// If $x(u,v),~y(u,v),~z(u,v),~\left\|\mathbf{n}(u,v)\right\|,~K(u,v)$ and $H(u,v)$ denote the $x$-, $y$-, $z$-coordinates, the length of}@*)
        //(*@\Green{// the normal vector, the Gaussian- and mean curvatures of the surface point $\mathbf{s}(u,v)$, respectively, then the}@*)
        //(*@\Green{// following color schemes reflect the pointwise variation of the energy quantity $\varphi(u,v)$, where:}@*)
        enum ImageColorScheme
        {
            DEFAULT_NULL_FRAGMENT = 0,            //(*@\Green{// $\varphi(u,v) = 0$}@*)
            X_VARIATION_FRAGMENT,                 //(*@\Green{// $\varphi(u,v) = x(u,v)$}@*)
            Y_VARIATION_FRAGMENT,                 //(*@\Green{// $\varphi(u,v) = y(u,v)$}@*)
            Z_VARIATION_FRAGMENT,                 //(*@\Green{// $\varphi(u,v) = z(u,v)$}@*)
            NORMAL_LENGTH_FRAGMENT,               //(*@\Green{// $\varphi(u,v) = \left\|\mathbf{n}(u,v)\right\|$}@*)
            GAUSSIAN_CURVATURE_FRAGMENT,          //(*@\Green{// $\varphi(u,v) = K(u,v) \left\|\mathbf{n}(u,v)\right\|$}@*)
            MEAN_CURVATURE_FRAGMENT,              //(*@\Green{// $\varphi(u,v) = H(u,v) \left\|\mathbf{n}(u,v)\right\|$}@*)
            WILLMORE_ENERGY_FRAGMENT,             //(*@\Green{// $\varphi(u,v) = H^2(u,v) \left\|\mathbf{n}(u,v)\right\|$}@*)
            LOG_WILLMORE_ENERGY_FRAGMENT,         //(*@\Green{// $\varphi(u,v) = \ln\left(1+H^2(u,v)\right)\left\|\mathbf{n}(u,v)\right\|$}@*)
            UMBILIC_DEVIATION_ENERGY_FRAGMENT,    //(*@\Green{// $\varphi(u,v) = 4\left(H^2(u,v) - K(u,v)\right) \left\|\mathbf{n}(u,v)\right\|$}@*)
            LOG_UMBILIC_DEVIATION_ENERGY_FRAGMENT,//(*@\Green{// $\varphi(u,v) = \ln\left(1+4\left(H^2(u,v) - K(u,v)\right)\right) \left\|\mathbf{n}(u,v)\right\|$}@*)
            TOTAL_CURVATURE_ENERGY_FRAGMENT,      //(*@\Green{// $\varphi(u,v) = \left(\frac{3}{2}H^2(u,v) - \frac{1}{2}K(u,v)\right) \left\|\mathbf{n}(u,v)\right\|$}@*)
            LOG_TOTAL_CURVATURE_ENERGY_FRAGMENT   //(*@\Green{// $\varphi(u,v) = \ln\left(1+\frac{3}{2}H^2(u,v) - \frac{1}{2}K(u,v)\right) \left\|\mathbf{n}(u,v)\right\|$}@*)
        }; //(*@\label{src:TensorProductSurface3.h:ImageColorScheme:end}@*)

        //(*@\Green{// a nested public class that represents a triangular matrix that consists of all (mixed) partial derivatives}@*)
        //(*@\Green{// up to a maximum order}@*)
        class PartialDerivatives: public TriangularMatrix<Cartesian3>
        {
        public:
            //(*@\Green{// default/special constructor}@*)
            PartialDerivatives(GLint maximum_order_of_partial_derivatives = 0);
            //(*@\Green{// when called, all inherited Cartesian coordinates are set to the origin}@*)
            GLvoid loadNullVectors();
            //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
            PartialDerivatives* clone() const;
        };

        //(*@\Green{// nested public class that represents a column matrix consisting of either $u$- or $v$-directional partial}@*)
        //(*@\Green{// derivatives up to a maximum order}@*)
        class DirectionalDerivatives: public ColumnMatrix<Cartesian3>
        {
        public:
            //(*@\Green{// default/special constructor}@*)
            DirectionalDerivatives(GLint maximum_order_of_directional_derivatives = 0);
            //(*@\Green{// when called, all inherited Cartesian coordinates are set to the origin}@*)
            GLvoid loadNullVectors();
            //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
            DirectionalDerivatives* clone() const;
        };

    protected:
        //(*@\Green{// endpoints of the $u$- and $v$-intervals of the definition domain}@*)
        std::vector<GLdouble>   _min_value, _max_value;
        //(*@\Green{// boolean flags that indicate whether the surface is closed in a given direction}@*)
        std::vector<GLboolean>  _closed;
        //(*@\Green{// vertex buffer object of the user-specified \_data that usually describes a control net}@*)
        GLuint                  _vbo_data;
        Matrix<Cartesian3>      _data; //(*@\Green{// $\left[\mathbf{p}_{i,j}\right]_{i=0,\,j=0}^{n,m} \in \mathcal{M}_{n+1,m+1}\left(\mathbb{R}^3\right)$}@*)

    public:
        //(*@\Green{// special constructor}@*)
        TensorProductSurface3(
                GLdouble u_min, GLdouble u_max,
                GLdouble v_min, GLdouble v_max,
                GLint row_count = 4, GLint column_count = 4,
                GLboolean u_closed = GL_FALSE, GLboolean v_closed = GL_FALSE);

        //(*@\Green{// copy constructor}@*)
        TensorProductSurface3(const TensorProductSurface3& surface);

        //(*@\Green{// assignment operator}@*)
        TensorProductSurface3& operator =(const TensorProductSurface3& surface);

        //(*@\Green{// sets/queries the definition domain of the surface in a given direction}@*)
        virtual GLboolean setInterval(variable::Type type,
                                      GLdouble min_value, GLdouble max_value);
        GLvoid interval(variable::Type type,
                        GLdouble &min_value, GLdouble &max_value) const;

        //(*@\Green{// get data by constant reference}@*)
        const Cartesian3& operator ()(GLint row, GLint column) const;

        //(*@\Green{// get data by non-constant reference}@*)
        Cartesian3& operator ()(GLint row, GLint column);

        //(*@\Green{// a pure virtual method that should calculate a row matrix that consists of either $u$- or $v$-directional}@*)
        //(*@\Green{// blending function values}@*)
        virtual GLboolean blendingFunctionValues(
                variable::Type type,
                GLdouble parameter_value, RowMatrix<GLdouble>& values) const = 0;

        //(*@\Green{// pure virtual methods that should calculate the point and higher order (mixed) partial derivatives of the}@*)
        //(*@\Green{// given surface}@*)
        virtual GLboolean calculateAllPartialDerivatives(
                GLint maximum_order_of_partial_derivatives,
                GLdouble u, GLdouble v, PartialDerivatives& pd) const = 0;

        virtual GLboolean calculateDirectionalDerivatives(
                variable::Type direction,
                GLint maximum_order_of_directional_derivatives,
                GLdouble u, GLdouble v, DirectionalDerivatives& d) const = 0;

        //(*@\Green{// generates a triangle mesh that approximates the shape of the given tensor product surface}@*)
        virtual TriangleMesh3* generateImage(
                GLint u_div_point_count, GLint v_div_point_count,
                ImageColorScheme color_scheme = DEFAULT_NULL_FRAGMENT,
                GLenum usage_flag = GL_STATIC_DRAW) const;

        //(*@\Green{// generates either $u$- or $v$-directional isoparametric lines of the given tensor product surface}@*)
        virtual RowMatrix<SP<GenericCurve3>::Default>* generateIsoparametricLines(
                variable::Type type, GLint iso_line_count,
                GLint maximum_order_of_derivatives, GLint div_point_count,
                GLenum usage_flag = GL_STATIC_DRAW) const;

        //(*@\Green{// Tries to solve the surface interpolation problem $\mathbf{s}\left(u_i, v_j\right) = \mathbf{d}_{i,j},~i=0,1,\ldots,n,~j=0,1,\ldots,m$,}@*)
        //(*@\Green{// where $\left[u_i\right]_{i=0}^n \subset \left[u_{\min}, u_{\max}\right]$ and $\left[v_j\right]_{j=0}^m \subset \left[v_{\min}, v_{\max}\right]$ are two strictly increasing knot vectors,}@*)
        //(*@\Green{// while $\left[\mathbf{d}_{i,j}\right]_{i=0,\,j=0}^{n,\,m}$ denotes the grid of data points that has to be interpolated.}@*)
        //(*@\Green{// In case of success, the solution will be stored in the member field \_data that corresponds to $\left[\mathbf{p}_{i,j}\right]_{i=0,\,j=0}^{n,m}$.}@*)
        GLboolean updateDataForInterpolation(
                const RowMatrix<GLdouble> &u_knot_vector,
                const ColumnMatrix<GLdouble> &v_knot_vector,
                Matrix<Cartesian3> &data_points_to_interpolate);

        //(*@\Green{// By default, it is assumed that the vectors stored in Matrix$<$Cartesian3$>$ \_data}@*)
        //(*@\Green{// represent the vertices of a control net. If this is not the case, the next virtual}@*)
        //(*@\Green{// VBO handling methods can be redeclared and redefined in derived classes.}@*)

        //(*@\Green{// Deletes the vertex buffer objects of the control net.}@*)
        virtual GLvoid    deleteVertexBufferObjectsOfData();

        //(*@\Green{// Control net rendering methods.}@*)

        //(*@\Green{// The next rendering method will return GL\_FALSE if:\label{src:TensorProductSurface3:render:start}}@*)
        //(*@\Green{// \ \ - the given shader program is not active;}@*)
        //(*@\Green{// \ \ - the user-defined position attribute name cannot be found in the list of active attributes}@*)
        //(*@\Green{// \ \ \ \ of the provided shader program, or exists but is of incorrect type;}@*)
        //(*@\Green{// \ \ - the u\_render\_mode or v\_render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP,}@*)
        //(*@\Green{// \ \ \ \  GL\_LINE\_LOOP and GL\_POINTS).}@*)
        virtual GLboolean renderData(
                const ShaderProgram &program,
                GLenum u_render_mode = GL_LINE_STRIP,
                GLenum v_render_mode = GL_LINE_STRIP) const;

        //(*@\Green{// If during rendering one intends to use shader program objects that are not instances of}@*)
        //(*@\Green{// our class ShaderProgram, then one should specify an attribute location associated with}@*)
        //(*@\Green{// positions of type vec3.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
        //(*@\Green{// \ \ - there is no active shader program object;}@*)
        //(*@\Green{// \ \ - the given position location either cannot be found in the list of active attribute locations,}@*)
        //(*@\Green{// \ \ \ \ or exists but is of incorrect type;}@*)
        //(*@\Green{// \ \ - the u\_render\_mode or v\_render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP,}@*)
        //(*@\Green{// \ \ \ \ GL\_LINE\_LOOP and GL\_POINTS).}@*)
        virtual GLboolean renderData(
                GLenum u_render_mode = GL_LINE_STRIP,
                GLenum v_render_mode = GL_LINE_STRIP,
                GLint vec3_position_location = 0) const; //(*@\label{src:TensorProductSurface3:render:end}@*)

        //(*@\Green{// Updates the vertex buffer objects of the control net.}@*)
        virtual GLboolean updateVertexBufferObjectsOfData(
                GLenum usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// get the row count of the data grid}@*)
        GLint rowCount() const;

        //(*@\Green{// get the column count of the data grid}@*)
        GLint columnCount() const;

        //(*@\Green{// destructor}@*)
        virtual ~TensorProductSurface3();
    };
}

#endif //(*@\Green{// TENSORPRODUCTSURFACES3\_H}@*)
