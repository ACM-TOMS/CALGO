//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/TensorProductSurfaces3.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "TensorProductSurfaces3.h"

#include "../../Math/RealMatrices.h"
#include "../../Math/RealMatrixDecompositions.h"
#include "../../Utilities.h"

#include <algorithm>
using namespace std;

namespace cagd
{
    (*@\Green{// default/special constructor}@*)
    TensorProductSurface3::PartialDerivatives::PartialDerivatives(
        GLint maximum_order_of_partial_derivatives):
        TriangularMatrix<Cartesian3>(
            maximum_order_of_partial_derivatives < 0 ?
            0 : maximum_order_of_partial_derivatives + 1)
    {
    }

    (*@\Green{// when called, all inherited Cartesian coordinates are set to the origin}@*)
    GLvoid TensorProductSurface3::PartialDerivatives::loadNullVectors()
    {
        #pragma omp parallel for
        for (GLint r = 0; r < _row_count; r++)
        {
            for (GLint c = 0; c <= r; c++)
            {
                Cartesian3 &partial_derivative = _data[r][c];

                for (GLint i = 0; i < 3; i++)
                {
                    partial_derivative[i] = 0.0;
                }
            }
        }
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    TensorProductSurface3::PartialDerivatives*
        TensorProductSurface3::PartialDerivatives::clone() const
    {
        return new (nothrow) TensorProductSurface3::PartialDerivatives(*this);
    }

    (*@\Green{// default/special constructor}@*)
    TensorProductSurface3::DirectionalDerivatives::DirectionalDerivatives(
        GLint maximum_order_of_directional_derivatives):
        ColumnMatrix<Cartesian3>(
            maximum_order_of_directional_derivatives < 0 ?
            0 : maximum_order_of_directional_derivatives + 1)
    {
    }

    (*@\Green{// when called, all inherited Cartesian coordinates are set to the origin}@*)
    GLvoid TensorProductSurface3::DirectionalDerivatives::loadNullVectors()
    {
        #pragma omp parallel for
        for (GLint r = 0; r < _row_count; r++)
        {
            Cartesian3 &directional_derivative = _data[r];

            for (GLint i = 0; i < 3; i++)
            {
                directional_derivative[i] = 0.0;
            }
        }
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    TensorProductSurface3::DirectionalDerivatives*
        TensorProductSurface3::DirectionalDerivatives::clone() const
    {
        return new (nothrow) TensorProductSurface3::DirectionalDerivatives(*this);
    }

    (*@\Green{// special constructor}@*)
    TensorProductSurface3::TensorProductSurface3(GLdouble u_min, GLdouble u_max,
            GLdouble v_min, GLdouble v_max,
            GLint row_count, GLint column_count,
            GLboolean u_closed, GLboolean v_closed):
        _min_value(2), _max_value(2), _closed(2),
        _vbo_data(0),
        _data(row_count, column_count)
    {
        _min_value[variable::U] = u_min,    _max_value[variable::U] = u_max;
        _min_value[variable::V] = v_min,    _max_value[variable::V] = v_max;
        _closed[variable::U]    = u_closed, _closed[variable::V]    = v_closed;
    }

    (*@\Green{// copy constructor}@*)
    TensorProductSurface3::TensorProductSurface3(const TensorProductSurface3 &surface):
        _min_value(surface._min_value),
        _max_value(surface._max_value),
        _closed(surface._closed),
        _vbo_data(0),
        _data(surface._data)
    {
        if (surface._vbo_data)
        {
            updateVertexBufferObjectsOfData();
        }
    }

    (*@\Green{// assignment operator}@*)
    TensorProductSurface3& TensorProductSurface3::operator =(
        const TensorProductSurface3 &rhs)
    {
        if (this != &rhs)
        {
            deleteVertexBufferObjectsOfData();

            _min_value = rhs._min_value;
            _max_value = rhs._max_value;
            _closed    = rhs._closed;

            _data = rhs._data;

            if (rhs._vbo_data)
            {
                updateVertexBufferObjectsOfData();
            }
        }

        return *this;
    }

    (*@\Green{// sets/queries the definition domain of the surface in a given direction}@*)
    GLboolean TensorProductSurface3::setInterval(
        variable::Type type,
        GLdouble min_value, GLdouble max_value)
    {
        _min_value[type] = min_value;
        _max_value[type] = max_value;

        return GL_TRUE;
    }

    void TensorProductSurface3::interval(
        variable::Type type,
        GLdouble &min_value, GLdouble &max_value) const
    {
        min_value = _min_value[type];
        max_value = _max_value[type];
    }

    (*@\Green{// get data by constant reference}@*)
    const Cartesian3& TensorProductSurface3::operator ()(GLint row, GLint column) const
    {
        assert("The given row index is out of bounds!" &&
               (row >= 0 && row < _data.rowCount()));

        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < _data.columnCount()));

        return _data(row, column);
    }

    (*@\Green{// get data by non-constant reference}@*)
    Cartesian3& TensorProductSurface3::operator ()(GLint row, GLint column)
    {
        assert("The given row index is out of bounds!" &&
               (row >= 0 && row < _data.rowCount()));

        assert("The given column index is out of bounds!" &&
               (column >= 0 && column < _data.columnCount()));

        return _data(row, column);
    }

    (*@\Green{// generates a triangle mesh that approximates the shape of the given tensor product surface}@*)
    TriangleMesh3* TensorProductSurface3::generateImage(
        GLint u_div_point_count, GLint v_div_point_count,
        ImageColorScheme color_scheme, GLenum usage_flag) const
    {
        if (u_div_point_count <= 1 || v_div_point_count <= 1)
        {
            return GL_FALSE;
        }

        (*@\Green{// calculating the number of vertices, unit normal vectors and texture coordinates}@*)
        GLint vertex_count = u_div_point_count * v_div_point_count;

        (*@\Green{// calculating the number of triangular faces}@*)
        GLint face_count = 2 * (u_div_point_count - 1) * (v_div_point_count - 1);

        TriangleMesh3 *result = nullptr;
        result = new (nothrow) TriangleMesh3(vertex_count, face_count, usage_flag);

        if (!result)
        {
            return nullptr;
        }

        (*@\Green{// step sizes of a uniform subdivision grid in the definition domain}@*)
        GLdouble du = (_max_value[variable::U] - _min_value[variable::U]) /
                      (u_div_point_count - 1);

        GLdouble dv = (_max_value[variable::V] - _min_value[variable::V]) /
                      (v_div_point_count - 1);

        (*@\Green{// step sizes of a uniform subdivision grid in the unit square}@*)
        GLfloat sdu = 1.0f / (u_div_point_count - 1);
        GLfloat tdv = 1.0f / (v_div_point_count - 1);

        (*@\Green{// values used for color scheme generation}@*)
        GLint maximum_order_of_partial_derivatives = 1;

        switch (color_scheme)
        {
        case DEFAULT_NULL_FRAGMENT:
        case X_VARIATION_FRAGMENT:
        case Y_VARIATION_FRAGMENT:
        case Z_VARIATION_FRAGMENT:
        case NORMAL_LENGTH_FRAGMENT:
            maximum_order_of_partial_derivatives = 1;
            break;

        case GAUSSIAN_CURVATURE_FRAGMENT:
        case MEAN_CURVATURE_FRAGMENT:
        case WILLMORE_ENERGY_FRAGMENT:
        case LOG_WILLMORE_ENERGY_FRAGMENT:
        case UMBILIC_DEVIATION_ENERGY_FRAGMENT:
        case LOG_UMBILIC_DEVIATION_ENERGY_FRAGMENT:
        case TOTAL_CURVATURE_ENERGY_FRAGMENT:
        case LOG_TOTAL_CURVATURE_ENERGY_FRAGMENT:
            maximum_order_of_partial_derivatives = 2;
            break;
        }

        std::vector<GLdouble> fragment(color_scheme > DEFAULT_NULL_FRAGMENT ?
                                       vertex_count : 0);

        (*@\Green{// error flag}@*)
        GLboolean aborted = GL_FALSE;

        #pragma omp parallel for
        for (GLint i_j = 0; i_j < vertex_count; i_j++)
        {
            #pragma omp flush (aborted)
            if (!aborted)
            {
                GLint i = i_j / v_div_point_count;
                GLint j = i_j % v_div_point_count;

                GLdouble u = min(_min_value[variable::U] + i * du,
                                 _max_value[variable::U]);

                GLfloat  s = min(i * sdu, 1.0f);

                GLdouble v = min(_min_value[variable::V] + j * dv,
                                 _max_value[variable::V]);

                GLfloat  t = min(j * tdv, 1.0f);

                GLint vertex_id[4];

                vertex_id[0] = i_j;
                vertex_id[1] = vertex_id[0] + 1;
                vertex_id[2] = vertex_id[1] + v_div_point_count;
                vertex_id[3] = vertex_id[2] - 1;

                (*@\Green{// calculating all needed surface data (i.e., partial derivatives)}@*)
                PartialDerivatives pd;
                if (!calculateAllPartialDerivatives(
                        maximum_order_of_partial_derivatives, u, v, pd))
                {
                    aborted = GL_TRUE;
                    #pragma omp flush (aborted)
                }
                else
                {
                    (*@\Green{// surface point}@*)
                    (*result)._position[vertex_id[0]] = pd(0, 0);

                    (*@\Green{// unit surface normal}@*)
                    Cartesian3 &normal = (*result)._normal[vertex_id[0]];
                    normal  = pd(1, 0);
                    normal ^= pd(1, 1);

                    GLdouble normal_length = normal.length();

                    if (normal_length && normal_length != 1.0)
                    {
                        normal /= normal_length;
                    }

                    (*@\Green{// coefficients of the first fundamental form}@*)
                    GLdouble e_1 = 0.0, f_1 = 0.0, g_1 = 0.0;

                    (*@\Green{// coefficients of the second fundamental form}@*)
                    GLdouble e_2 = 0.0, f_2 = 0.0, g_2 = 0.0;

                    (*@\Green{// Gaussian curvature}@*)
                    GLdouble K   = 0.0;

                    (*@\Green{// mean curvature}@*)
                    GLdouble H   = 0.0;

                    if (color_scheme > NORMAL_LENGTH_FRAGMENT)
                    {
                        (*@\Green{// coefficients of the first fundamental form}@*)
                        e_1 = pd(1, 0) * pd(1, 0);
                        f_1 = pd(1, 0) * pd(1, 1);
                        g_1 = pd(1, 1) * pd(1, 1);

                        (*@\Green{// coefficients of the second fundamental form}@*)
                        e_2 = normal * pd(2, 0);
                        f_2 = normal * pd(2, 1);
                        g_2 = normal * pd(2, 2);
                    }

                    if (color_scheme >= GAUSSIAN_CURVATURE_FRAGMENT)
                    {
                        K = (e_2 * g_2 - f_2 * f_2) / (e_1 * g_1 - f_1 * f_1);
                    }

                    if (color_scheme >= MEAN_CURVATURE_FRAGMENT)
                    {
                        H = (g_1 * e_2 - 2.0 * f_1 * f_2 + e_1 * g_2) /
                            (e_1 * g_1 - f_1 * f_1);
                    }

                    (*@\Green{// texture coordinates}@*)
                    TCoordinate4 &tex = (*result)._tex[vertex_id[0]];
                    tex.s() = s;
                    tex.t() = t;

                    (*@\Green{// connectivity information}@*)
                    if (i < u_div_point_count - 1 && j < v_div_point_count - 1)
                    {
                        GLint face_index = 2 * (i_j - i);

                        TriangularFace &lower_face = (*result)._face[face_index];
                        lower_face[0] = vertex_id[0];
                        lower_face[1] = vertex_id[1];
                        lower_face[2] = vertex_id[2];

                        TriangularFace &upper_face = (*result)._face[face_index + 1];
                        upper_face[0] = vertex_id[0];
                        upper_face[1] = vertex_id[2];
                        upper_face[2] = vertex_id[3];
                    }

                    (*@\Green{// fragment generation based on the selected color scheme}@*)
                    if (color_scheme > DEFAULT_NULL_FRAGMENT)
                    {
                        switch (color_scheme)
                        {
                        case X_VARIATION_FRAGMENT:
                            fragment[vertex_id[0]] = pd(0, 0).x();
                            break;

                        case Y_VARIATION_FRAGMENT:
                            fragment[vertex_id[0]] = pd(0, 0).y();
                            break;

                        case Z_VARIATION_FRAGMENT:
                            fragment[vertex_id[0]] = pd(0, 0).z();
                            break;

                        case NORMAL_LENGTH_FRAGMENT:
                            fragment[vertex_id[0]] = normal_length;
                            break;

                        case GAUSSIAN_CURVATURE_FRAGMENT:
                            fragment[vertex_id[0]] = K * normal_length;
                            break;

                        case MEAN_CURVATURE_FRAGMENT:
                            fragment[vertex_id[0]] = H * normal_length;
                            break;

                        case WILLMORE_ENERGY_FRAGMENT:
                            fragment[vertex_id[0]] = H * H * normal_length;
                            break;

                        case LOG_WILLMORE_ENERGY_FRAGMENT:
                            fragment[vertex_id[0]] = log(1.0 + H * H) * normal_length;
                            break;

                        case UMBILIC_DEVIATION_ENERGY_FRAGMENT:
                            fragment[vertex_id[0]] = 4.0 * (H * H - K) * normal_length;
                            break;

                        case LOG_UMBILIC_DEVIATION_ENERGY_FRAGMENT:
                            fragment[vertex_id[0]] = log(1.0 + 4.0 * (H * H - K)) *
                                                     normal_length;
                            break;

                        case TOTAL_CURVATURE_ENERGY_FRAGMENT:
                            fragment[vertex_id[0]] = (1.5 * H * H - 0.5 * K) *
                                                     normal_length;
                            break;

                        case LOG_TOTAL_CURVATURE_ENERGY_FRAGMENT:
                            fragment[vertex_id[0]] = log(1.0 + 1.5 * H * H - 0.5 * K) *
                                                     normal_length;
                            break;

                        default:
                            break;
                        }
                    }
                }
            }
        }

        if (aborted)
        {
            delete result, result = nullptr;
            return result;
        }

        (*@\Green{// generating energy based pointwise color variation}@*)
        if (color_scheme > DEFAULT_NULL_FRAGMENT)
        {
            GLdouble min_fragment_value =  numeric_limits<GLdouble>::max();
            GLdouble max_fragment_value = -numeric_limits<GLdouble>::max();

            for (GLint i = 0; i < vertex_count; i++)
            {
                if (min_fragment_value > fragment[i])
                {
                    min_fragment_value = fragment[i];
                }

                if (max_fragment_value < fragment[i])
                {
                    max_fragment_value = fragment[i];
                }
            }

            #pragma omp parallel for
            for (GLint i = 0; i < vertex_count; i++)
            {
                (*result)._color[i] = coldToHotColormap((GLfloat)fragment[i],
                                                        (GLfloat)min_fragment_value,
                                                        (GLfloat)max_fragment_value);
            }
        }

        return result;
    }

    (*@\Green{// generates either $u$- or $v$-directional isoparametric lines of the given tensor product surface}@*)
    RowMatrix<SP<GenericCurve3>::Default>*
        TensorProductSurface3::generateIsoparametricLines(
                variable::Type type, GLint iso_line_count,
                GLint maximum_order_of_derivatives, GLint div_point_count,
                GLenum usage_flag) const
    {
        if (iso_line_count < 2 || maximum_order_of_derivatives < 0 ||
            div_point_count < 2)
        {
            return nullptr;
        }

        if (usage_flag != GL_STREAM_DRAW  && usage_flag != GL_STREAM_READ  &&
            usage_flag != GL_STREAM_COPY  &&
            usage_flag != GL_DYNAMIC_DRAW && usage_flag != GL_DYNAMIC_READ &&
            usage_flag != GL_DYNAMIC_COPY &&
            usage_flag != GL_STATIC_DRAW  && usage_flag != GL_STATIC_READ  &&
            usage_flag != GL_STATIC_COPY)
        {
            return nullptr;
        }

        RowMatrix<SP<GenericCurve3>::Default> *result =
            new (nothrow) RowMatrix<SP<GenericCurve3>::Default>(iso_line_count);

        if (!result)
        {
            return nullptr;
        }

        GLdouble iso_line_step = (_max_value[1-type] - _min_value[1-type]) /
                                 (iso_line_count - 1);

        GLdouble div_step = (_max_value[type] - _min_value[type]) /
                            (div_point_count - 1);

        for (GLint l = 0; l < iso_line_count; l++)
        {
            (*result)[l] = SP<GenericCurve3>::Default(new (nothrow) GenericCurve3(
                    maximum_order_of_derivatives, div_point_count, usage_flag));

            if (!(*result)[l])
            {
                delete result, result = nullptr;
                return nullptr;
            }

            GenericCurve3 &iso_line_l = *(*result)[l];

            GLdouble iso_value = min(_min_value[1-type] + l * iso_line_step,
                                     _max_value[1-type]);

            (*@\Green{// error flag}@*)
            GLboolean aborted = GL_FALSE;

            #pragma omp parallel for
            for (GLint k = 0; k < div_point_count; k++)
            {
                #pragma omp flush (aborted)
                if (!aborted)
                {
                    GLdouble div_value = min(_min_value[type] + k * div_step,
                                             _max_value[type]);

                    DirectionalDerivatives d;
                    if (!calculateDirectionalDerivatives(
                            type, maximum_order_of_derivatives,
                            (type == variable::V ? iso_value : div_value),
                            (type == variable::U ? iso_value : div_value), d))
                    {
                        aborted = GL_TRUE;
                        #pragma omp flush (aborted)
                    }
                    else
                    {
                        iso_line_l._derivative.setColumn(k, d);
                    }
                }
            }

            if (aborted)
            {
                delete result, result = nullptr;
                return nullptr;
            }
        }

        return result;
    }

    (*@\Green{// Tries to solve the surface interpolation problem $\mathbf{s}\left(u_i, v_j\right) = \mathbf{d}_{i,j},~i=0,1,\ldots,n,~j=0,1,\ldots,m$,}@*)
    (*@\Green{// where $\left[u_i\right]_{i=0}^n \subset \left[u_{\min}, u_{\max}\right]$ and $\left[v_j\right]_{j=0}^m \subset \left[v_{\min}, v_{\max}\right]$ are two strictly increasing knot vectors,}@*)
    (*@\Green{// while $\left[\mathbf{d}_{i,j}\right]_{i=0,\,j=0}^{n,\,m}$ denotes the grid of data points that has to be interpolated.}@*)
    (*@\Green{// In case of success, the solution will be stored in the member field \_data that corresponds to $\left[\mathbf{p}_{i,j}\right]_{i=0,\,j=0}^{n,m}$.}@*)
    GLboolean TensorProductSurface3::updateDataForInterpolation(
            const RowMatrix<GLdouble> &u_knot_vector,
            const ColumnMatrix<GLdouble> &v_knot_vector,
            Matrix<Cartesian3> &data_points_to_interpolate)
    {
        GLint row_count = _data.rowCount();
        if (!row_count)
        {
            return GL_FALSE;
        }

        GLint column_count = _data.columnCount();
        if (!column_count)
        {
            return GL_FALSE;
        }

        if (u_knot_vector.columnCount() != row_count ||
            v_knot_vector.rowCount() != column_count ||
            data_points_to_interpolate.rowCount() != row_count ||
            data_points_to_interpolate.columnCount() != column_count)
        {
            return GL_FALSE;
        }

        (*@\Green{// 1: calculate the $u$-collocation matrix}@*)
        GLboolean  u_collocation_matrix_construction_aborted = GL_FALSE;
        RealMatrix u_collocation_matrix(row_count, row_count);

        #pragma omp parallel for
        for (GLint i = 0; i < row_count; i++)
        {
            #pragma omp flush (u_collocation_matrix_construction_aborted)
            if (!u_collocation_matrix_construction_aborted)
            {
                RowMatrix<GLdouble> u_blending_values;
                if (!blendingFunctionValues(
                        variable::U, u_knot_vector[i], u_blending_values))
                {
                    u_collocation_matrix_construction_aborted = GL_TRUE;
                    #pragma omp flush (u_collocation_matrix_construction_aborted)
                }
                else
                {
                    u_collocation_matrix.setRow(i, u_blending_values);
                }
            }
        }

        if (u_collocation_matrix_construction_aborted)
        {
            return GL_FALSE;
        }

        (*@\Green{// 2: calculate the $v$-collocation matrix and perform $PLU$-decomposition on it}@*)
        GLboolean  v_collocation_matrix_construction_aborted = GL_FALSE;
        RealMatrix v_collocation_matrix(column_count, column_count);

        #pragma omp parallel for
        for (GLint j = 0; j < column_count; ++j)
        {
            #pragma omp flush (v_collocation_matrix_construction_aborted)
            if (!v_collocation_matrix_construction_aborted)
            {
                RowMatrix<GLdouble> v_blending_values;
                if (!blendingFunctionValues(
                        variable::V, v_knot_vector[j], v_blending_values))
                {
                    v_collocation_matrix_construction_aborted = GL_TRUE;
                    #pragma omp flush (v_collocation_matrix_construction_aborted)
                }
                else
                {
                    v_collocation_matrix.setRow(j, v_blending_values);
                }
            }
        }

        if (v_collocation_matrix_construction_aborted)
        {
            return GL_FALSE;
        }

        (*@\Green{// 3: perfom $PLU$-decomposition on the $u$-collocation matrix and for all fixed $j \in \{0, 1,\ldots, m\}$}@*)
        (*@\Green{// \ \ \ \ determine intermediate vectors $\mathbf{a}_k(v_j) = \sum_{\ell=0}^{m} \mathbf{p}_{k,\ell} b_{m,\ell}(v_j), k = 0, 1,\ldots, n$ such that}@*)
        (*@\Green{// \ \ \ \ $\sum_{k=0}^{n} \mathbf{a}_k(v_j) b_{n,k}(u_i) = \mathbf{d}_{i,j}$ for all $i = 0, 1,\ldots, n$}@*)
        PLUDecomposition   u_PLUD(u_collocation_matrix);
        Matrix<Cartesian3> a(row_count, column_count);
        if (!u_PLUD.solveLinearSystem(data_points_to_interpolate, a))
        {
            return GL_FALSE;
        }

        (*@\Green{// 4: perfom $PLU$-decomposition on the $v$-collocation matrix and for all fixed $i \in \{0, 1,\ldots, n\}$}@*)
        (*@\Green{// \ \ \ \ determine control points $\mathbf{p}_{i,j}, j = 0, 1,..., m$ such that $\sum_{\ell=0}^{m} \mathbf{p}_{i,\ell} b_{m,\ell}(v_j) = \mathbf{a}_i(v_{j})$ for all $j = 0, 1,\ldots, m$}@*)
        PLUDecomposition v_PLUD(v_collocation_matrix);
        if (!v_PLUD.solveLinearSystem(a, _data, GL_FALSE))
        {
            return GL_FALSE;
        }

        return GL_TRUE;
    }

    (*@\Green{// By default, it is assumed that the vectors stored in Matrix$<$Cartesian3$>$ \_data}@*)
    (*@\Green{// represent the vertices of a control net. If this is not the case, the next virtual}@*)
    (*@\Green{// VBO handling methods can be redeclared and redefined in derived classes.}@*)

    (*@\Green{// Deletes the vertex buffer objects of the control net.}@*)
    void TensorProductSurface3::deleteVertexBufferObjectsOfData()
    {
        if (_vbo_data)
        {
            glDeleteBuffers(1, &_vbo_data);
            _vbo_data = 0;
        }
    }

    (*@\Green{// Control net rendering methods.}@*)

    (*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
    (*@\Green{// \ \ - the given shader program is not active;}@*)
    (*@\Green{// \ \ - the user-defined position attribute name cannot be found in the list of active attributes}@*)
    (*@\Green{// \ \ \ \ of the provided shader program, or exists but is of incorrect type;}@*)
    (*@\Green{// \ \ - the u\_render\_mode or v\_render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP,}@*)
    (*@\Green{// \ \ \ \  GL\_LINE\_LOOP and GL\_POINTS).}@*)
    GLboolean TensorProductSurface3::renderData(
        const ShaderProgram &program,
        GLenum u_render_mode, GLenum v_render_mode) const
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

        if (u_render_mode != GL_LINE_STRIP && u_render_mode != GL_LINE_LOOP &&
            u_render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        if (v_render_mode != GL_LINE_STRIP && v_render_mode != GL_LINE_LOOP &&
            v_render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        if (u_render_mode == GL_LINE_STRIP && _closed[variable::U])
        {
            u_render_mode = GL_LINE_LOOP;
        }

        if (v_render_mode == GL_LINE_STRIP && _closed[variable::V])
        {
            v_render_mode = GL_LINE_LOOP;
        }

        if (u_render_mode == GL_LINE_LOOP && !_closed[variable::U])
        {
            u_render_mode = GL_LINE_STRIP;
        }

        if (v_render_mode == GL_LINE_LOOP && !_closed[variable::V])
        {
            v_render_mode = GL_LINE_STRIP;
        }

        GLint row_count = _data.rowCount();

        if (!row_count)
        {
            return GL_FALSE;
        }

        GLint column_count = _data.columnCount();

        if (!column_count)
        {
            return GL_FALSE;
        }

        GLsizei size = row_count * column_count;

        glEnableVertexAttribArray(position);

            glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);

            glVertexAttribPointer(position, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

            if (u_render_mode == GL_POINTS && v_render_mode == GL_POINTS)
            {
                glDrawArrays(u_render_mode, 0, size);
            }
            else
            {
                GLint offset = 0;

                for (GLint i = 0; i < row_count; i++, offset += column_count)
                {
                    glDrawArrays(u_render_mode, offset, column_count);
                }

                for (GLint j = 0; j < column_count; j++, offset += row_count)
                {
                    glDrawArrays(v_render_mode, offset, row_count);
                }
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
    (*@\Green{// \ \ - there is no active shader program object;}@*)
    (*@\Green{// \ \ - the given position location either cannot be found in the list of active attribute locations,}@*)
    (*@\Green{// \ \ \ \ or exists but is of incorrect type;}@*)
    (*@\Green{// \ \ - the u\_render\_mode or v\_render\_mode is incorrect (one should use one of the constants GL\_LINE\_STRIP,}@*)
    (*@\Green{// \ \ \ \ GL\_LINE\_LOOP and GL\_POINTS).}@*)
    GLboolean TensorProductSurface3::renderData(
            GLenum u_render_mode, GLenum v_render_mode,
            GLint vec3_position_location) const
    {
        if (!_vbo_data || vec3_position_location < 0)
        {
            return GL_FALSE;
        }

        if (u_render_mode != GL_LINE_STRIP && u_render_mode != GL_LINE_LOOP &&
            u_render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        if (v_render_mode != GL_LINE_STRIP && v_render_mode != GL_LINE_LOOP &&
            v_render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        if (u_render_mode == GL_LINE_STRIP && _closed[variable::U])
        {
            u_render_mode = GL_LINE_LOOP;
        }

        if (v_render_mode == GL_LINE_STRIP && _closed[variable::V])
        {
            v_render_mode = GL_LINE_LOOP;
        }

        if (u_render_mode == GL_LINE_LOOP && !_closed[variable::U])
        {
            u_render_mode = GL_LINE_STRIP;
        }

        if (v_render_mode == GL_LINE_LOOP && !_closed[variable::V])
        {
            v_render_mode = GL_LINE_STRIP;
        }

        GLint row_count = _data.rowCount();

        if (!row_count)
        {
            return GL_FALSE;
        }

        GLint column_count = _data.columnCount();

        if (!column_count)
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

        GLsizei size = row_count * column_count;

        glEnableVertexAttribArray(vec3_position_location);

            glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);

            glVertexAttribPointer(vec3_position_location, 3,
                                  GL_FLOAT, GL_FALSE, 0, nullptr);

            if (u_render_mode == GL_POINTS && v_render_mode == GL_POINTS)
            {
                glDrawArrays(u_render_mode, 0, size);
            }
            else
            {
                GLint offset = 0;

                for (GLint i = 0; i < row_count; i++, offset += column_count)
                {
                    glDrawArrays(u_render_mode, offset, column_count);
                }

                for (GLint j = 0; j < column_count; j++, offset += row_count)
                {
                    glDrawArrays(v_render_mode, offset, row_count);
                }
            }

            glBindBuffer(GL_ARRAY_BUFFER, 0);

        glDisableVertexAttribArray(vec3_position_location);
        glVertexAttrib3f(vec3_position_location, 0.0f, 0.0f, 0.0f);

        return GL_TRUE;
    }

    (*@\Green{// Updates the vertex buffer objects of the control net.}@*)
    GLboolean TensorProductSurface3::updateVertexBufferObjectsOfData(GLenum usage_flag)
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

        if (!_data.rowCount() || !_data.columnCount())
        {
            return GL_FALSE;
        }

        deleteVertexBufferObjectsOfData();

        glGenBuffers(1, &_vbo_data);

        if (!_vbo_data)
        {
            return GL_FALSE;
        }

        GLsizeiptr vertex_byte_size =
            6 * (GLsizeiptr)(_data.rowCount() * _data.columnCount() * sizeof(GLfloat));

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);
        glBufferData(GL_ARRAY_BUFFER, vertex_byte_size, nullptr, usage_flag);

        GLfloat *vertex_coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,
                                                           GL_WRITE_ONLY);

        for (GLint i = 0; i < _data.rowCount(); i++)
        {
            for (GLint j = 0; j < _data.columnCount(); j++)
            {
                for (GLint c = 0; c < 3; c++, vertex_coordinate++)
                {
                    *vertex_coordinate = (GLfloat)_data(i, j)[c];
                }
            }
        }

        for (GLint j = 0; j < _data.columnCount(); j++)
        {
            for (GLint i = 0; i < _data.rowCount(); i++)
            {
                for (GLint c = 0; c < 3; c++, vertex_coordinate++)
                {
                    *vertex_coordinate = (GLfloat)_data(i, j)[c];
                }
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);

        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
        {
            return GL_FALSE;
        }

        return GL_TRUE;
    }

    (*@\Green{// get the row count of the data grid}@*)
    GLint TensorProductSurface3::rowCount() const
    {
        return _data.rowCount();
    }

    (*@\Green{// get the column count of the data grid}@*)
    GLint TensorProductSurface3::columnCount() const
    {
        return _data.columnCount();
    }

    (*@\Green{// destructor}@*)
    TensorProductSurface3::~TensorProductSurface3()
    {
        deleteVertexBufferObjectsOfData();
    }
}
