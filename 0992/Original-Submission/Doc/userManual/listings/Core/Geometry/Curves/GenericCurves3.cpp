//----------------------------------------------------------------------------------
// File:        Core/Geometry/Curves/GenericCurves3.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "GenericCurves3.h"
#include <fstream>

using namespace std;

namespace cagd
{
    (*@\Green{// default/special constructor}@*)
    GenericCurve3::GenericCurve3(
        GLint maximum_order_of_derivatives, GLint point_count, GLenum usage_flag):
        _usage_flag(usage_flag),
        _vbo_derivative(maximum_order_of_derivatives < 0 ?
                        0 : maximum_order_of_derivatives + 1),
        _derivative(maximum_order_of_derivatives < 0 ?
                    0 : maximum_order_of_derivatives + 1,
                    point_count < 0 ? 0 : point_count)
    {
        assert("The maximum order of derivatives should be a non-negative integer!" &&
               maximum_order_of_derivatives >= 0);

        assert("The number of subdivision points should be a non-negative integer!" &&
               point_count >= 0);

        assert("Invalid usage flag!" &&
               (usage_flag == GL_STREAM_DRAW  || usage_flag == GL_STREAM_READ  ||
                usage_flag == GL_STREAM_COPY  ||
                usage_flag == GL_DYNAMIC_DRAW || usage_flag == GL_DYNAMIC_READ ||
                usage_flag == GL_DYNAMIC_COPY ||
                usage_flag == GL_STATIC_DRAW  || usage_flag == GL_STATIC_READ  ||
                usage_flag == GL_STATIC_COPY));
    }

    (*@\Green{// copy constructor}@*)
    GenericCurve3::GenericCurve3(const GenericCurve3& curve):
            _usage_flag(curve._usage_flag),
            _vbo_derivative(curve._vbo_derivative.columnCount()),
            _derivative(curve._derivative)
    {
        GLboolean vbo_update_is_possible = GL_TRUE;
        for (GLint i = 0; i < curve._vbo_derivative.columnCount(); i++)
        {
            vbo_update_is_possible &= curve._vbo_derivative[i];
        }

        if (vbo_update_is_possible)
        {
            updateVertexBufferObjects(_usage_flag);
        }
    }

    (*@\Green{// assignment operator}@*)
    GenericCurve3& GenericCurve3::operator =(const GenericCurve3& rhs)
    {
        if (this != &rhs)
        {
            deleteVertexBufferObjects();

            _usage_flag = rhs._usage_flag;
            _derivative = rhs._derivative;

            GLboolean vbo_update_is_possible = GL_TRUE;
            for (GLint i = 0; i < rhs._vbo_derivative.columnCount(); i++)
            {
                vbo_update_is_possible &= rhs._vbo_derivative[i];
            }

            if (vbo_update_is_possible)
            {
                updateVertexBufferObjects(_usage_flag);
            }
        }
        return *this;
    }

    (*@\Green{// Deletes the vertex buffer objects of zeroth and higher order derivatives.}@*)
    GLvoid GenericCurve3::deleteVertexBufferObjects()
    {
        for (GLint i = 0; i < _vbo_derivative.columnCount(); i++)
        {
            if (_vbo_derivative[i])
            {
                glDeleteBuffers(1, &_vbo_derivative[i]);
                _vbo_derivative[i] = 0;
            }
        }
    }

    (*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
    (*@\Green{// \ - the given shader program is not active;}@*)
    (*@\Green{// \ - the user-defined position attribute name cannot be found in the list of active attributes}@*)
    (*@\Green{// \ \ \ of the provided shader program, or exists but is of incorrect type;}@*)
    (*@\Green{// \ - the specified order of derivatives is greater than or equal to the row count of the}@*)
    (*@\Green{// \ \ \ Matrix$<$Cartesian3$>$ \_derivative;}@*)
    (*@\Green{// \ - the render\_mode is incorrect (in case of zeroth order derivatives one should use one of the constants}@*)
    (*@\Green{// \ \ \ GL\_LINE\_STRIP, GL\_LINE\_LOOP and GL\_POINTS, while in case of higher order derivatives}@*)
    (*@\Green{// \ \ \ one has to use either GL\_LINES or GL\_POINTS).}@*)
    GLboolean GenericCurve3::renderDerivatives(
        const ShaderProgram &program, GLint order, GLenum render_mode) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The vertex buffer object corresponding to the given differentiation "
               "order does not exist!" && _vbo_derivative[order]);

        if (order < 0 || order >= _derivative.rowCount() || !_vbo_derivative[order])
        {
            return GL_FALSE;
        }

        GLint current_program;
        glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);

        if (current_program != (GLint)program._id)
        {
            return GL_FALSE;
        }

        GLint position = program.positionAttributeLocation();

        if (position < 0)
        {
            return GL_FALSE;
        }

        GLsizei point_count = _derivative.columnCount();

        glEnableVertexAttribArray(position);
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_derivative[order]);
        glVertexAttribPointer(position, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

        if (!order)
        {
            if (render_mode != GL_LINE_STRIP &&
                render_mode != GL_LINE_LOOP  &&
                render_mode != GL_POINTS)
            {
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glDisableVertexAttribArray(position);
                return GL_FALSE;
            }

            glDrawArrays(render_mode, 0, point_count);
        }
        else
        {
            if (render_mode != GL_LINES && render_mode != GL_POINTS)
            {
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glDisableVertexAttribArray(position);
                return GL_FALSE;
            }

            glDrawArrays(render_mode, 0, 2 * point_count);
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(position);
        glVertexAttrib3f(position, 0.0f, 0.0f, 0.0f);

        return GL_TRUE;
    }

    (*@\Green{// If during rendering one intends to use shader program objects that are not instances of}@*)
    (*@\Green{// our class ShaderProgram, then one should specify an attribute location associated with}@*)
    (*@\Green{// positions of type vec3.}@*)
    (*@\Green{//}@*)
    (*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
    (*@\Green{// \ - there is no active shader program object;}@*)
    (*@\Green{// \ - the given position location either cannot be found in the list of active attribute locations,}@*)
    (*@\Green{// \ \ \ or exists but is of incorrect type;}@*)
    (*@\Green{// \ - the specified order of derivatives is greater than or equal to the row count of the}@*)
    (*@\Green{// \ \ \ Matrix$<$Cartesian3$>$ \_derivative;}@*)
    (*@\Green{// \ - the render\_mode is incorrect (in case of zeroth order derivatives one should use one of the constants}@*)
    (*@\Green{// \ \ \ GL\_LINE\_STRIP, GL\_LINE\_LOOP and GL\_POINTS, while in case of higher order derivatives}@*)
    (*@\Green{// \ \ \ one has to use either GL\_LINES or GL\_POINTS).}@*)
    GLboolean GenericCurve3::renderDerivatives(
        GLint order, GLenum render_mode, GLint vec3_position_location) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The vertex buffer object corresponding to the given differentiation "
               "order does not exist!" && _vbo_derivative[order]);

        assert("The given position attribute location should be a non-negative "
               "integer!" && vec3_position_location >= 0);

        if (order < 0 || order >= _derivative.rowCount() || !_vbo_derivative[order] ||
            vec3_position_location < 0)
        {
            return GL_FALSE;
        }

        GLint current_program = 0;
        glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);

        if (!current_program)
        {
            return GL_FALSE;
        }

        GLint attribute_count;
        glGetProgramiv(current_program, GL_ACTIVE_ATTRIBUTES, &attribute_count);

        if (!attribute_count)
        {
            return GL_FALSE;
        }

        GLsizei max_attribute_name_length;
        glGetProgramiv(current_program, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH,
                       &max_attribute_name_length);

        GLchar *attribute_name_data = new GLchar[max_attribute_name_length];

        GLboolean given_location_exists_and_is_of_type_vec3 = GL_FALSE;

        for (GLint attribute = 0;
             attribute < attribute_count && !given_location_exists_and_is_of_type_vec3;
             attribute++)
        {
            GLsizei actual_length = 0;
            GLint   array_size    = 0;
            GLenum  type          = 0;

            glGetActiveAttrib(current_program, attribute, max_attribute_name_length,
                              &actual_length, &array_size, &type, attribute_name_data);
            string name(&attribute_name_data[0], actual_length);
            GLint  location = glGetAttribLocation(current_program, name.c_str());

            if (type == GL_FLOAT_VEC3 && vec3_position_location == location)
            {
                given_location_exists_and_is_of_type_vec3 = GL_TRUE;
            }
        }

        delete[] attribute_name_data;

        if (!given_location_exists_and_is_of_type_vec3)
        {
            return GL_FALSE;
        }

        GLint point_count = _derivative.columnCount();

        glEnableVertexAttribArray(vec3_position_location);
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_derivative[order]);
        glVertexAttribPointer(vec3_position_location, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

        if (!order)
        {
            if (render_mode != GL_LINE_STRIP &&
                render_mode != GL_LINE_LOOP  &&
                render_mode != GL_POINTS)
            {
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glDisableVertexAttribArray(vec3_position_location);
                return GL_FALSE;
            }

            glDrawArrays(render_mode, 0, point_count);
        }
        else
        {
            if (render_mode != GL_LINES && render_mode != GL_POINTS)
            {
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glDisableVertexAttribArray(vec3_position_location);
                return GL_FALSE;
            }

            glDrawArrays(render_mode, 0, 2 * point_count);
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(vec3_position_location);
        glVertexAttrib3f(vec3_position_location, 0.0f, 0.0f, 0.0f);

        return GL_TRUE;
    }

    (*@\Green{// Updates the vertex buffer objects of the zeroth and higher order derivatives.}@*)
    GLboolean GenericCurve3::updateVertexBufferObjects(GLenum usage_flag)
    {
        if (usage_flag != GL_STREAM_DRAW  && usage_flag != GL_STREAM_READ  &&
            usage_flag != GL_STREAM_COPY  &&
            usage_flag != GL_DYNAMIC_DRAW && usage_flag != GL_DYNAMIC_READ &&
            usage_flag != GL_DYNAMIC_COPY &&
            usage_flag != GL_STATIC_DRAW  && usage_flag != GL_STATIC_READ  &&
            usage_flag != GL_STATIC_COPY)
        {
            return GL_FALSE;
        }

        deleteVertexBufferObjects();

        _usage_flag = usage_flag;

        for (GLint d = 0; d < _vbo_derivative.columnCount(); d++)
        {
            glGenBuffers(1, &_vbo_derivative[d]);

            if (!_vbo_derivative[d])
            {
                for (GLint i = 0; i < d; i++)
                {
                    glDeleteBuffers(1, &_vbo_derivative[i]);
                    _vbo_derivative[i] = 0;
                }

                return GL_FALSE;
            }
        }

        GLint curve_point_count = _derivative.columnCount();

        GLfloat *coordinate = nullptr;

        (*@\Green{// curve points}@*)
        GLsizeiptr curve_point_byte_size = 3 * curve_point_count * sizeof(GLfloat);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_derivative[0]);
        glBufferData(GL_ARRAY_BUFFER, curve_point_byte_size, nullptr, _usage_flag);

        coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        if (!coordinate)
        {
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            deleteVertexBufferObjects();
            return GL_FALSE;
        }

        for (GLint i = 0; i < curve_point_count; i++)
        {
            for (GLint j = 0; j < 3; j++, coordinate++)
            {
                *coordinate = (GLfloat)_derivative(0, i)[j];
            }
        }

        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
        {
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            deleteVertexBufferObjects();
            return GL_FALSE;
        }

        (*@\Green{// higher order derivatives}@*)
        GLsizeiptr higher_order_derivative_byte_size = 2 * curve_point_byte_size;

        for (GLint d = 1; d < _derivative.rowCount(); d++)
        {
            glBindBuffer(GL_ARRAY_BUFFER, _vbo_derivative[d]);
            glBufferData(GL_ARRAY_BUFFER, higher_order_derivative_byte_size, nullptr,
                         _usage_flag);

            coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

            if (!coordinate)
            {
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                deleteVertexBufferObjects();
                return GL_FALSE;
            }

            for (GLint i = 0; i < curve_point_count; i++, coordinate += 3)
            {
                Cartesian3 sum = _derivative(0, i);
                sum += _derivative(d, i);

                for (GLint j = 0; j < 3; j++, coordinate++)
                {
                    *coordinate       = (GLfloat)_derivative(0, i)[j];
                    *(coordinate + 3) = (GLfloat)sum[j];
                }
            }

            if (!glUnmapBuffer(GL_ARRAY_BUFFER))
            {
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                deleteVertexBufferObjects();
                return GL_FALSE;
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return GL_TRUE;
    }

    (*@\Green{// One can either map or unmap a vertex buffer object associated with a specific order of derivatives.}@*)
    GLfloat* GenericCurve3::mapDerivatives(GLint order, GLenum access_mode) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _vbo_derivative.columnCount()));

        assert("The vertex buffer object corresponding to the given differentiation "
               "order does not exist!" && !_vbo_derivative[order]);

        if (order < 0 || order >= _vbo_derivative.columnCount() ||
            !_vbo_derivative[order])
        {
            return nullptr;
        }

        if (access_mode != GL_READ_ONLY && access_mode != GL_WRITE_ONLY &&
            access_mode != GL_READ_WRITE)
        {
            return nullptr;
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_derivative[order]);

        return (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, access_mode);
    }

    GLboolean GenericCurve3::unmapDerivatives(GLint order) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _vbo_derivative.columnCount()));

        assert("The vertex buffer object corresponding to the given differentiation "
               "order does not exist!" && !_vbo_derivative[order]);

        if (order < 0 || order >= _vbo_derivative.columnCount() ||
            !_vbo_derivative[order])
        {
            return GL_FALSE;
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_derivative[order]);

        return glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    (*@\Green{// get derivative by constant reference}@*)
    const Cartesian3& GenericCurve3::operator ()(GLint order, GLint index) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The given curve point index is out of bounds!" &&
               (index >= 0 && index < _derivative.columnCount()));

        return _derivative(order, index);
    }

    (*@\Green{// get derivative by non-constant reference}@*)
    Cartesian3& GenericCurve3::operator ()(GLint order, GLint index)
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The given curve point index is out of bounds!" &&
               (index >= 0 && index < _derivative.columnCount()));

        return _derivative(order, index);
    }

    (*@\Green{// other updating and querying methods}@*)
    GLboolean GenericCurve3::setDerivative(
        GLint order, GLint index, GLdouble x, GLdouble y, GLdouble z)
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The given curve point index is out of bounds!" &&
               (index >= 0 && index < _derivative.columnCount()));

        if (order < 0 || order >= _derivative.rowCount() ||
            index < 0 || index >= _derivative.columnCount())
        {
            return GL_FALSE;
        }

        Cartesian3 &reference = _derivative(order, index);

        reference[0] = x;
        reference[1] = y;
        reference[2] = z;

        return GL_TRUE;
    }

    GLboolean GenericCurve3::setDerivative(GLint order, GLint index, const Cartesian3& d)
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The given curve point index is out of bounds!" &&
               (index >= 0 && index < _derivative.columnCount()));

        if (order < 0 || order >= _derivative.rowCount() ||
            index < 0 || index >= _derivative.columnCount())
        {
            return GL_FALSE;
        }

        _derivative(order, index) = d;

        return GL_TRUE;
    }

    GLboolean GenericCurve3::derivative(
        GLint order, GLint index, GLdouble& x, GLdouble& y, GLdouble& z) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The given curve point index is out of bounds!" &&
               (index >= 0 && index < _derivative.columnCount()));

        if (order < 0 || order >= _derivative.rowCount() ||
            index < 0 || index >= _derivative.columnCount())
        {
            return GL_FALSE;
        }

        const Cartesian3 &reference = _derivative(order, index);

        x = reference[0];
        y = reference[1];
        z = reference[2];

        return GL_TRUE;
    }

    GLboolean GenericCurve3::derivative(GLint order, GLint index, Cartesian3& d) const
    {
        assert("The given differentiation order is out of bounds!" &&
               (order >= 0 && order < _derivative.rowCount()));

        assert("The given curve point index is out of bounds!" &&
               (index >= 0 && index < _derivative.columnCount()));

        if (order < 0 || order >= _derivative.rowCount() ||
            index < 0 || index >= _derivative.columnCount())
        {
            return GL_FALSE;
        }

        d = _derivative(order, index);

        return GL_TRUE;
    }

    GLint GenericCurve3::maximumOrderOfDerivatives() const
    {
        return _derivative.rowCount() - 1;
    }

    GLint GenericCurve3::pointCount() const
    {
        return _derivative.columnCount();
    }

    GLenum GenericCurve3::usageFlag() const
    {
        return _usage_flag;
    }

    (*@\Green{// The next method can be used for generating Matlab code, by using which one can create scalable}@*)
    (*@\Green{// vector graphic file formats like EPS.}@*)
    (*@\Green{// Variable file\_name specifies the name of the output file. Make sure that its extension is ``.m".}@*)
    (*@\Green{// Variable access\_mode can be either std::ios\_base::out, or std::ios\_base::out $|$ std::ios\_base::app, in case of}@*)
    (*@\Green{// which the given file will be either overwritten, or appended.}@*)
    (*@\Green{// Variable line\_color specifies the red, green and blue components of the color `Color' property of the}@*)
    (*@\Green{// Matlab-commands plot and plot3.}@*)
    (*@\Green{// Variable line\_style can be used to specify formatum strings for Matlab-like line styles, which are used by}@*)
    (*@\Green{// the Matlab-commands plot and plot3. It accepts one of the formatum strings ``'-'", ``'-.'" and ``':'".}@*)
    (*@\Green{// Variable line\_width specifies the positive line width which will be used by the Matlab-commands plot}@*)
    (*@\Green{// and plot3.}@*)
    (*@\Green{// Variables \{x$|$y$|$z\}\_coordinate\_name store the names of arrays that will contain the $\{x|y|z\}$-coordinates}@*)
    (*@\Green{// of the curve points. These array names will be passed as input arguments to the plotting commands}@*)
    (*@\Green{// of Matlab.}@*)
    GLboolean GenericCurve3::generateMatlabCodeForRendering(const string &file_name,
        ios_base::openmode access_mode,
        const Color4 &line_color,
        const std::string &line_style,
        const GLdouble &line_width,
        const string &x_coordinate_name,
        const string &y_coordinate_name,
        const string &z_coordinate_name) const
    {
        if (file_name.empty() || x_coordinate_name.empty() ||
            y_coordinate_name.empty() || z_coordinate_name.empty() ||
            (x_coordinate_name == y_coordinate_name) ||
            (x_coordinate_name == z_coordinate_name) ||
            (y_coordinate_name == z_coordinate_name) ||
            line_width <= 0.0 || _derivative.columnCount() == 0)
        {
            return GL_FALSE;
        }

        if (line_style != "'-'" && line_style != "'-.'" && line_style != "':'")
        {
            return GL_FALSE;
        }

        if ((access_mode != (fstream::out | fstream::app)) &&
            (access_mode != fstream::out))
        {
            return GL_FALSE;
        }

        fstream f(file_name.c_str(), access_mode);

        if (!f || !f.good())
        {
            f.close();
            return GL_FALSE;
        }

        f << endl << "hold all;" << endl << endl;

        f << x_coordinate_name << " = [..." << endl;

        for (GLint i = 0; i < _derivative.columnCount(); i++)
        {
            f << _derivative(0, i)[0] << " ";
        }
        f << "..." << endl << "];" << endl;

        f << endl;

        f << y_coordinate_name << " = [..." << endl;

        for (GLint i = 0; i < _derivative.columnCount(); i++)
        {
            f << _derivative(0, i)[1] << " ";
        }
        f << "..." << endl << "];" << endl;

        f << endl;

        f << z_coordinate_name << " = [..." << endl;

        for (GLint i = 0; i < _derivative.columnCount(); i++)
        {
            f << _derivative(0, i)[2] << " ";
        }
        f << "..." << endl << "];" << endl;

        f << endl << "plot3(" << x_coordinate_name << ", "
                              << y_coordinate_name << ", "
                              << z_coordinate_name << ", "
                              << line_style << ", "
                              << "'Color', ["
                              << line_color[0] << ", "
                              << line_color[1] << ", "
                              << line_color[2] << "], "
                              << "'LineWidth', " << line_width
                              << ");" << endl;


        f << endl << "axis equal;" << endl;

        return GL_TRUE;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    GenericCurve3* GenericCurve3::clone() const
    {
        return new (nothrow) GenericCurve3(*this);
    }

    (*@\Green{// destructor}@*)
    GenericCurve3::~GenericCurve3()
    {
        deleteVertexBufferObjects();
    }

    (*@\Green{// overloaded friend input/output from/to stream operators}@*)
    ostream& operator <<(ostream& lhs, const GenericCurve3& rhs)
    {
        return lhs << rhs._usage_flag << " " << rhs._derivative << endl;
    }

    std::istream& operator >>(std::istream& lhs, GenericCurve3& rhs)
    {
        rhs.deleteVertexBufferObjects();

        return lhs >> rhs._usage_flag >> rhs._derivative;
    }
}
