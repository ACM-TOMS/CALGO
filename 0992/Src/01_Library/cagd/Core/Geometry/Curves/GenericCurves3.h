//----------------------------------------------------------------------------------
// File:        Core/Geometry/Curves/GenericCurves3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef GENERICCURVES3_H
#define GENERICCURVES3_H

#include <GL/glew.h>
#include <iostream>

#include "../Coordinates/Cartesians3.h"
#include "../Coordinates/Colors4.h"
#include "../../Math/Matrices.h"
#include "../../Shaders/ShaderPrograms.h"

namespace cagd
{
    class GenericCurve3
    {
        //(*@\Green{// friend classes that will be discussed later in Listings \mref{src:LinearCombinations3.h}/\mref{src:LinearCombinations3.cpp} and \mref{src:TensorProductSurfaces3.h}/\mref{src:TensorProductSurfaces3.cpp}, respectively}@*)
        friend class LinearCombination3;
        friend class TensorProductSurface3;

        //(*@\Green{// overloaded friend input/output from/to stream operators}@*)
        friend std::ostream& operator <<(std::ostream& lhs, const GenericCurve3& rhs);
        friend std::istream& operator >>(std::istream& lhs, GenericCurve3& rhs);

    protected:
        GLenum              _usage_flag;
        RowMatrix<GLuint>   _vbo_derivative;
        Matrix<Cartesian3>  _derivative;

    public:
        //(*@\Green{// default/special constructor}@*)
        GenericCurve3(GLint maximum_order_of_derivatives = 2,
                      GLint point_count = 0,
                      GLenum usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// copy constructor}@*)
        GenericCurve3(const GenericCurve3& curve);

        //(*@\Green{// assignment operator}@*)
        GenericCurve3& operator =(const GenericCurve3& rhs);

        //(*@\Green{// Deletes the vertex buffer objects of zeroth and higher order derivatives.}@*)
        GLvoid    deleteVertexBufferObjects();

        //(*@\Green{// The next rendering method will return GL\_FALSE if:\label{src:GenericCurve3:renderDerivatives:start}}@*)
        //(*@\Green{// \ - the given shader program is not active;}@*)
        //(*@\Green{// \ - the user-defined position attribute name cannot be found in the list of active attributes}@*)
        //(*@\Green{// \ \ \ of the provided shader program, or exists but is of incorrect type;}@*)
        //(*@\Green{// \ - the specified order of derivatives is greater than or equal to the row count of the}@*)
        //(*@\Green{// \ \ \ Matrix$<$Cartesian3$>$ \_derivative;}@*)
        //(*@\Green{// \ - the render\_mode is incorrect (in case of zeroth order derivatives one should use one of the constants}@*)
        //(*@\Green{// \ \ \ GL\_LINE\_STRIP, GL\_LINE\_LOOP and GL\_POINTS, while in case of higher order derivatives}@*)
        //(*@\Green{// \ \ \ one has to use either GL\_LINES or GL\_POINTS).}@*)
        GLboolean renderDerivatives(const ShaderProgram &program, GLint order,
                                    GLenum render_mode) const;

        //(*@\Green{// If during rendering one intends to use shader program objects that are not instances of}@*)
        //(*@\Green{// our class ShaderProgram, then one should specify an attribute location associated with}@*)
        //(*@\Green{// positions of type vec3.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
        //(*@\Green{// \ - there is no active shader program object;}@*)
        //(*@\Green{// \ - the given position location either cannot be found in the list of active attribute locations,}@*)
        //(*@\Green{// \ \ \ or exists but is of incorrect type;}@*)
        //(*@\Green{// \ - the specified order of derivatives is greater than or equal to the row count of the}@*)
        //(*@\Green{// \ \ \ Matrix$<$Cartesian3$>$ \_derivative;}@*)
        //(*@\Green{// \ - the render\_mode is incorrect (in case of zeroth order derivatives one should use one of the constants}@*)
        //(*@\Green{// \ \ \ GL\_LINE\_STRIP, GL\_LINE\_LOOP and GL\_POINTS, while in case of higher order derivatives}@*)
        //(*@\Green{// \ \ \ one has to use either GL\_LINES or GL\_POINTS).}@*)
        GLboolean renderDerivatives(GLint order, GLenum render_mode,
                                    GLint vec3_position_location = 0) const; //(*@\label{src:GenericCurve3:renderDerivatives:end}@*)

        //(*@\Green{// Updates the vertex buffer objects of the zeroth and higher order derivatives.}@*)
        GLboolean updateVertexBufferObjects(GLenum usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// One can either map or unmap a vertex buffer object associated with a specific order of derivatives.}@*)
        GLfloat*  mapDerivatives(GLint order, GLenum access_mode = GL_READ_ONLY) const;
        GLboolean unmapDerivatives(GLint order) const;

        //(*@\Green{// get derivative by constant reference}@*)
        const Cartesian3& operator ()(GLint order, GLint index) const;

        //(*@\Green{// get derivative by non-constant reference}@*)
        Cartesian3& operator ()(GLint order, GLint index);

        //(*@\Green{// other updating and querying methods}@*)
        GLboolean setDerivative(GLint order, GLint index,
                                GLdouble x, GLdouble y, GLdouble z = 0.0);

        GLboolean setDerivative(GLint order, GLint index, const Cartesian3& d);

        GLboolean derivative(GLint order, GLint index,
                             GLdouble &x, GLdouble &y, GLdouble &z) const;

        GLboolean derivative(GLint order, GLint index, Cartesian3 &d) const;

        GLint     maximumOrderOfDerivatives() const;
        GLint     pointCount() const;
        GLenum    usageFlag() const;

        //(*@\Green{// The next method can be used for generating Matlab code, by using which one can create scalable}@*)
        //(*@\Green{// vector graphic file formats like EPS.}@*)
        //(*@\Green{// Variable file\_name specifies the name of the output file. Make sure that its extension is ``.m".}@*)
        //(*@\Green{// Variable access\_mode can be either std::ios\_base::out, or std::ios\_base::out $|$ std::ios\_base::app, in case of}@*)
        //(*@\Green{// which the given file will be either overwritten, or appended.}@*)
        //(*@\Green{// Variable line\_color specifies the red, green and blue components of the color `Color' property of the}@*)
        //(*@\Green{// Matlab-commands plot and plot3.}@*)
        //(*@\Green{// Variable line\_style can be used to specify formatum strings for Matlab-like line styles, which are used by}@*)
        //(*@\Green{// the Matlab-commands plot and plot3. It accepts one of the formatum strings ``'-'", ``'-.'" and ``':'".}@*)
        //(*@\Green{// Variable line\_width specifies the positive line width which will be used by the Matlab-commands plot}@*)
        //(*@\Green{// and plot3.}@*)
        //(*@\Green{// Variables \{x$|$y$|$z\}\_coordinate\_name store the names of arrays that will contain the $\{x|y|z\}$-coordinates}@*)
        //(*@\Green{// of the curve points. These array names will be passed as input arguments to the plotting commands}@*)
        //(*@\Green{// of Matlab.}@*)
        GLboolean generateMatlabCodeForRendering(
            const std::string &file_name,
            std::ios_base::openmode access_mode = std::ios_base::out | std::ios_base::app,
            const Color4 &line_color = colors::blue,
            const std::string &line_style = "'-'",
            const GLdouble &line_width = 1.0,
            const std::string &x_coordinate_name = "x",
            const std::string &y_coordinate_name = "y",
            const std::string &z_coordinate_name = "z") const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual GenericCurve3* clone() const;

        //(*@\Green{// destructor}@*)
        virtual ~GenericCurve3();
    };
}

#endif // GENERICCURVES3_H
