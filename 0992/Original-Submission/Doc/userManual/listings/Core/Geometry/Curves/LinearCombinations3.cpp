//----------------------------------------------------------------------------------
// File:        Core/Geometry/Curves/LinearCombinations3.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "LinearCombinations3.h"

#include "../../Math/RealMatrices.h"
#include "../../Math/RealMatrixDecompositions.h"

#include <algorithm>

using namespace std;

namespace cagd
{
    (*@\Green{// default/special constructor}@*)
    LinearCombination3::Derivatives::Derivatives(GLint maximum_order_of_derivatives):
        ColumnMatrix<Cartesian3>(maximum_order_of_derivatives + 1)
    {
    }

    (*@\Green{// when called, all inherited Cartesian coordinates are set to the origin}@*)
    GLvoid LinearCombination3::Derivatives::loadNullVectors()
    {
        #pragma omp parallel for
        for (GLint i = 0; i < _row_count; i++)
        {
            Cartesian3 &p_i = _data[i];

            for (GLint j = 0; j < 3; j++)
            {
                p_i[j] = 0.0;
            }
        }
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    LinearCombination3::Derivatives* LinearCombination3::Derivatives::clone() const
    {
        return new (nothrow) Derivatives(*this);
    }

    (*@\Green{// special constructor}@*)
    LinearCombination3::LinearCombination3(GLdouble u_min, GLdouble u_max,
                                           GLint data_count, GLenum data_usage_flag):
        _vbo_data(0),
        _data_usage_flag(data_usage_flag),
        _u_min(u_min), _u_max(u_max),
        _data(data_count)
    {
    }

    (*@\Green{// copy constructor}@*)
    LinearCombination3::LinearCombination3(const LinearCombination3 &lc):
        _vbo_data(0),
        _data_usage_flag(lc._data_usage_flag),
        _u_min(lc._u_min), _u_max(lc._u_max),
        _data(lc._data)
    {
        if (lc._vbo_data)
        {
            updateVertexBufferObjectsOfData(_data_usage_flag);
        }
    }

    (*@\Green{// assignment operator}@*)
    LinearCombination3& LinearCombination3::operator =(const LinearCombination3& rhs)
    {
        if (this != &rhs)
        {
            deleteVertexBufferObjectsOfData();

            _data_usage_flag = rhs._data_usage_flag;
            _u_min = rhs._u_min;
            _u_max = rhs._u_max;
            _data = rhs._data;

            if (rhs._vbo_data)
            {
                updateVertexBufferObjectsOfData(_data_usage_flag);
            }
        }

        return *this;
    }

    (*@\Green{// Deletes the vertex buffer objects of the control polygon.}@*)
    GLvoid LinearCombination3::deleteVertexBufferObjectsOfData()
    {
        if (_vbo_data)
        {
            glDeleteBuffers(1, &_vbo_data);
            _vbo_data = 0;
        }
    }

    (*@\Green{// Control polygon rendering methods.}@*)

    (*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
    (*@\Green{// \ \ - the given shader program is not active;}@*)
    (*@\Green{// \ \ - the user-defined position attribute name cannot be found in the list of active attributes}@*)
    (*@\Green{// \ \ \ \ of the provided shader program, or exists but is of incorrect type;}@*)
    (*@\Green{// \ \ - the render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP, GL\_LINE\_LOOP}@*)
    (*@\Green{// \ \ \ \ and GL\_POINTS).}@*)
    GLboolean LinearCombination3::renderData(const ShaderProgram &program,
                                             GLenum render_mode) const
    {
        GLint current_program;
        glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);

        if (current_program != (GLint)program._id)
        {
            return GL_FALSE;
        }

        GLint position = program.positionAttributeLocation();

        if (!_vbo_data || position < 0)
        {
            return GL_FALSE;
        }

        if (render_mode != GL_LINE_STRIP && render_mode != GL_LINE_LOOP &&
            render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        glEnableVertexAttribArray(position);
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);

        glVertexAttribPointer(position, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glDrawArrays(render_mode, 0, _data.rowCount());

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
    (*@\Green{// \ \ - there is no active shader program object;}@*)
    (*@\Green{// \ \ - the given position location either cannot be found in the list of active attribute locations,}@*)
    (*@\Green{// \ \ \ \ or exists but is of incorrect type;}@*)
    (*@\Green{// \ \ - the render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP, GL\_LINE\_LOOP}@*)
    (*@\Green{// \ \ \ \ and GL\_POINTS).}@*)
    GLboolean LinearCombination3::renderData(GLenum render_mode,
                                             GLint vec3_position_location) const
    {
        if (!_vbo_data || vec3_position_location < 0)
        {
            return GL_FALSE;
        }

        if (render_mode != GL_LINE_STRIP && render_mode != GL_LINE_LOOP &&
            render_mode != GL_POINTS)
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

        glEnableVertexAttribArray(vec3_position_location);
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);

        glVertexAttribPointer(vec3_position_location, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glDrawArrays(render_mode, 0, (GLsizei)_data.rowCount());

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(vec3_position_location);
        glVertexAttrib3f(vec3_position_location, 0.0f, 0.0f, 0.0f);

        return GL_TRUE;
    }

    (*@\Green{// Updates the vertex buffer objects of the control polygon.}@*)
    GLboolean LinearCombination3::updateVertexBufferObjectsOfData(GLenum usage_flag)
    {
        GLint data_count = _data.rowCount();

        if (!data_count)
        {
            return GL_FALSE;
        }

        if (usage_flag != GL_STREAM_DRAW  && usage_flag != GL_STREAM_READ  &&
            usage_flag != GL_STREAM_COPY  &&
            usage_flag != GL_DYNAMIC_DRAW && usage_flag != GL_DYNAMIC_READ &&
            usage_flag != GL_DYNAMIC_COPY &&
            usage_flag != GL_STATIC_DRAW  && usage_flag != GL_STATIC_READ  &&
            usage_flag != GL_STATIC_COPY)
        {
            return GL_FALSE;
        }

        _data_usage_flag = usage_flag;

        deleteVertexBufferObjectsOfData();

        glGenBuffers(1, &_vbo_data);
        if (!_vbo_data)
            return GL_FALSE;

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);
        glBufferData(GL_ARRAY_BUFFER, data_count * 3 * sizeof(GLfloat), nullptr,
                     _data_usage_flag);

        GLfloat *coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        if (!coordinate)
        {
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            deleteVertexBufferObjectsOfData();
            return GL_FALSE;
        }

        for (GLint i = 0; i < data_count; i++)
        {
            Cartesian3 &p_i = _data[i];

            for (GLint j = 0; j < 3; j++, coordinate++)
            {
                *coordinate = (GLfloat)p_i[j];
            }
        }

        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
        {
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            deleteVertexBufferObjectsOfData();
            return GL_FALSE;
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return GL_TRUE;
    }

    (*@\Green{// get data by constant reference}@*)
    const Cartesian3& LinearCombination3::operator [](const GLint &index) const
    {
        assert("The given index is out of bounds!" &&
               (index >= 0 && index < _data.rowCount()));

        return _data[index];
    }

    (*@\Green{// get data by non-constant reference}@*)
    Cartesian3& LinearCombination3::operator [](const GLint &index)
    {
        assert("The given index is out of bounds!" &&
               (index >= 0 && index < _data.rowCount()));

        return _data[index];
    }

    (*@\Green{// sets/returns the endpoints of the definition domain}@*)
    GLboolean LinearCombination3::setDefinitionDomain(
        const GLdouble &u_min, const GLdouble &u_max)
    {
        if (u_min >= u_max)
        {
            return GL_FALSE;
        }

        _u_min = u_min;
        _u_max = u_max;

        return GL_TRUE;
    }

    GLvoid LinearCombination3::definitionDomain(GLdouble& u_min, GLdouble& u_max) const
    {
        u_min = _u_min;
        u_max = _u_max;
    }

    (*@\Green{// get number of data}@*)
    GLint LinearCombination3::dataCount() const
    {
        return _data.rowCount();
    }

    (*@\Green{// generates the image of the given linear combination}@*)
    GenericCurve3* LinearCombination3::generateImage(
        GLint maximum_order_of_derivatives,
        GLint div_point_count, GLenum usage_flag) const
    {                
        GenericCurve3 *result = new (nothrow) GenericCurve3(
                maximum_order_of_derivatives, div_point_count, usage_flag);

        if (!result)
        {
            return nullptr;
        }

        GLdouble step = (_u_max - _u_min) / (div_point_count - 1);

        GLboolean aborted = GL_FALSE;

        #pragma omp parallel for
        for (GLint i = 0; i < div_point_count; i++)
        {
            #pragma omp flush (aborted)
            if (!aborted)
            {
                GLdouble u = min(_u_max, _u_min + i * step);

                Derivatives d;
                if (!calculateDerivatives(maximum_order_of_derivatives, u, d))
                {
                    aborted = GL_TRUE;
                    #pragma omp flush (aborted)
                }

                result->_derivative.setColumn(i, d);
            }
        }

        if (aborted)
        {
            delete result, result = nullptr;
        }

        return result;
    }

    (*@\Green{// Solves the user-specified curve interpolation problem $\mathbf{c}\left(u_j\right) = \mathbf{d}_j \in \mathbb{R}^3,~j=0,1,\ldots,n$, where $\left[u_j\right]_{j=0}^n$}@*)
    (*@\Green{// forms a knot vector consisting of strictly increasing subdivision points of the definition domain $[u_{\min}, u_{\max}]$,}@*)
    (*@\Green{// while $\left[\mathbf{d}_j\right]_{j=0}^n$ denotes the data points that have to be interpolated.}@*)
    (*@\Green{// In case of success, the solution $\left[\mathbf{p}_i\right]_{i=0}^n$ will be stored in the column matrix \_data.}@*)
    GLboolean LinearCombination3::updateDataForInterpolation(
        const ColumnMatrix<GLdouble>& knot_vector,
        const ColumnMatrix<Cartesian3>& data_points_to_interpolate)
    {
        GLint data_count = _data.rowCount();

        if (data_count != knot_vector.rowCount() ||
            data_count != data_points_to_interpolate.rowCount())
        {
            return GL_FALSE;
        }

        GLboolean  collocation_matrix_construction_aborted = GL_FALSE;
        RealMatrix collocation_matrix(data_count, data_count);

        #pragma omp parallel for
        for (GLint r = 0; r < knot_vector.rowCount(); r++)
        {
            #pragma omp flush (collocation_matrix_construction_aborted)
            if (!collocation_matrix_construction_aborted)
            {
                RowMatrix<GLdouble> current_blending_function_values(data_count);
                if (!blendingFunctionValues(knot_vector[r],
                                            current_blending_function_values))
                {
                    collocation_matrix_construction_aborted = GL_TRUE;
                    #pragma omp flush (collocation_matrix_construction_aborted)
                }
                else
                {
                    collocation_matrix.setRow(r, current_blending_function_values);
                }
            }
        }

        if (collocation_matrix_construction_aborted)
        {
            return GL_FALSE;
        }

        PLUDecomposition PLUD(collocation_matrix);

        return PLUD.solveLinearSystem(data_points_to_interpolate, _data);
    }

    (*@\Green{// destructor}@*)
    LinearCombination3::~LinearCombination3()
    {
        deleteVertexBufferObjectsOfData();
    }
}
