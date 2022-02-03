#pragma warning(disable:4503)
#include "YourGLWidget.h"

#include <iostream>
#include <fstream>

#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Geometry/Surfaces/Materials.h>
#include <Core/Math/Constants.h>
#include "../Spaces/SpecializedECSpaces.h"

using namespace std;
using namespace cagd;
using namespace cagd::variable;

(*@\Green{// your default/special constructor}@*)
void YourGLWidget::YourGLWidget(/*...*/)
{
    try
    {
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:light:end} of Listing \mref{src:GLWidgetTestingShaderPrograms.cpp}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// Communicating with our shader programs via uniform variables\ldots}@*)
        _two_sided_lighting.enable();

        if (!_two_sided_lighting.setUniformDirectionalLight("light_source[0]", *_light))
        {
            throw Exception("Two-sided per pixed lighting: could not initialize the "
                            "uniform variable \"light_source[0]\"!");
        }

        if (!_two_sided_lighting.setUniformValue1i("light_source[0].enabled", GL_TRUE))
        {
            throw Exception("Two-sided per pixed lighting: could not initialize the "
                            "uniform variable \"light_source[0].enabled\"!");
        }

        _two_sided_lighting.disable();

        _reflection_lines.enable();

        if (!_reflection_lines.setUniformDirectionalLight("light_source[0]", *_light))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"light_source[0]\"!");
        }

        if (!_reflection_lines.setUniformValue1i("light_source[0].enabled", GL_TRUE))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"light_source[0].enabled\"!");
        }

        if (!_reflection_lines.setUniformValue1f("scale_factor", 3.0f))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"scale_factor\"!");
        }

        if (!_reflection_lines.setUniformValue1f("smoothing", 1.0f))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"smoothing\"!");
        }

        if (!_reflection_lines.setUniformValue1f("shading", 0.01f))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"shading\"!");
        }

        _reflection_lines.disable();

        (*@\Green{// Loading a triangulated unit sphere centered at the origin. Using translation and scaling transformations,}@*)
        (*@\Green{// it will be rendered multiple times at the positions of the control points.}@*)
        if (!_sphere.loadFromOFF("Models/sphere.off"))
        {
            throw Exception("Could not load the model file sphere.off!");
        }

        if (!_sphere.updateVertexBufferObjects())
        {
            throw Exception("Could not update the VBOs of the triangulated sphere!");
        }

        (*@\Green{// Determines the common scaling factor of the unit sphere that has to be rendered at the positions}@*)
        (*@\Green{// of the control points.}@*)
        _control_point_radius = 0.06875;


        (*@\Green{// Determines the scaling transformation of the unit sphere that has to be rendered at the positions}@*)
        (*@\Green{// of the control points.}@*)
        _sphere_S.setScalingFactors(
                _control_point_radius, _control_point_radius, _control_point_radius);

        (*@\Green{// Loading a triangulated right circular cone with unit base radius and appex $(0, 0, 2)$.}@*)
        (*@\Green{// Using translation and scaling transformation, this cone will be rendered multiple times at the endpoints}@*)
        (*@\Green{// of the tangent vectors of the isoparametric lines.}@*)
        if (!_cone.loadFromOFF("Models/cone.off"))
        {
            throw Exception("Could not load model file cone.off!");
        }

        if (!_cone.updateVertexBufferObjects())
        {
            throw Exception("Could not update the VBO of the triangulated cone!");
        }
        _cone_S.setScalingFactors(_control_point_radius / 1.5f,
                                  _control_point_radius / 1.5f,
                                  _control_point_radius / 1.5f);

        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:start}--\mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:end} of Listing \mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// Memory allocations and parameter settings\ldots}@*)
        _dimension.resize(2);
        _surface_div_point_count.resize(2);
        _isoparametric_line_count.resize(2);
        _maximum_order_of_derivatives.resize(2);
        _curve_div_point_count.resize(2);

        _dimension[U]                    =  7;
        _surface_div_point_count[U]      = 50;
        _isoparametric_line_count[U]     =  5;
        _maximum_order_of_derivatives[U] =  1;
        _curve_div_point_count[U]        = 20;

        _dimension[V]                    =  3;
        _surface_div_point_count[V]      = 70;
        _isoparametric_line_count[V]     =  3;
        _maximum_order_of_derivatives[V] =  1;
        _curve_div_point_count[V]        = 10;

        (*@\Green{// select a color scheme\ldots}@*)
        _color_scheme = TensorProductSurface3::DEFAULT_NULL_FRAGMENT;

        (*@\Green{// In order to provide B-representations of the patches of the ordinary integral surface}@*)
        (*@\Green{// $\mathbf{s}\left(u,v\right) =\def\arraystretch{1.3}\left[\begin{array}{c}s^0\left(u, v\right)\\s^1\left(u, v\right)\\s^2\left(u, v\right)\end{array}\right]=\left[\begin{array}{c}\left(1-e^{\omega_0 u}\right) \cos\left(u\right) \left(\frac{5}{4}+\cos\left(v\right)\right)\\\left(e^{\omega_0 u}-1\right) \sin\left(u\right) \left(\frac{5}{4}+\cos\left(v\right)\right)\\7-e^{\omega_1 u} - \sin\left(v\right) + e^{\omega_0 u} \sin\left(v\right)\end{array}\right],\def\arraystretch{1.0}$}@*)
        (*@\Green{// we will define the exponential-trigonometric and pure trigonometric EC spaces}@*)
        (*@\Green{// $\mathbb{ET}_6^{\alpha_0^k,\beta_0^k}=\Big\langle\Big\{\varphi_{6,0}\left(u\right) \equiv 1,\varphi_{6,1}\left(u\right)=\cos\left(u\right),\varphi_{6,2}\left(u\right)=\sin\left(u\right),\varphi_{6,3}\left(u\right)=e^{\omega_0 u},\varphi_{6,4}\left(u\right)=e^{\omega_1 u},$}@*)
        (*@\Green{// \hspace{1.75cm}$\left.\left.\varphi_{6,5}\left(u\right)=e^{\omega_0 u} \cos\left(u\right),\varphi_{6,6}\left(u\right)=e^{\omega_0 u} \sin\left(u\right):u \in \left[\alpha_0^k,\beta_0^k\right]\right\}\right\rangle,~k=0,1,\ldots,4$}@*)
        (*@\Green{// and}@*)
        (*@\Green{// \hspace{0.25cm}$\mathbb{T}_2^{\alpha_1^{\ell}, \beta_1^{\ell}} =\left\langle\left\{\varphi_{2,0}\left(v\right) \equiv 1, \varphi_{2,1}\left(v\right)=\cos\left(v\right),\varphi_{2,2}\left(v\right)=\sin\left(v\right) :v \in \left[\alpha_1^{\ell}, \beta_1^{\ell}\right]\right\}\right\rangle,~\ell = 0,1,\ldots,\left\lfloor 3\tau \right\rfloor,$}@*)
        (*@\Green{// respectively, where $\omega_0 = \frac{1}{6\pi}$, $\omega_1 = \frac{1}{3\pi}$, $\big\{\alpha_0^k\big\}_{k=0}^4 = \left\{\frac{7\pi}{2}+\frac{29\pi}{40}\cdot k\right\}_{k=0}^4$, $\big\{\beta_0^k\big\}_{k=0}^4  = \big\{\alpha_0^k + \frac{29\pi}{40}\big\}_{k=0}^4$,}@*)
        (*@\Green{// $\big\{\alpha_1^{\ell}\big\}_{\ell=0}^{\left\lfloor 3\tau\right\rfloor} = \left\{-\frac{\pi}{3\tau}+\frac{2\pi}{\left\lfloor 3\tau\right\rfloor}\cdot\ell\right\}_{\ell=0}^{\left\lfloor 3\tau\right\rfloor}$, $\big\{\beta_1^{\ell}\big\}_{\ell=0}^{\left\lfloor 3\tau\right\rfloor} = \big\{\alpha_1^{\ell} + \frac{2\pi}{\left\lfloor 3\tau \right\rfloor}\big\}_{\ell=0}^{\left\lfloor 3\tau\right\rfloor}$ and $\tau \geq 1$. Moreover, we have to store}@*)
        (*@\Green{// the function coefficients that appear in the parametric form of the given ordinary integral surface.}@*)
        (*@\Green{// To achieve this, we will instantiate and initialize an object from the class OrdinarySurfaceCoefficients.}@*)
        (*@\Green{// For more details consider the lines \mref{src:BSurfaces3.h:OrdinarySurfaceCoefficients:start}--\mref{src:BSurfaces3.h:OrdinarySurfaceCoefficients:end} of Listing \mref{src:BSurfaces3.h}.}@*)
        std::vector<GLint> sigma(3);
        sigma[0] = 1;
        sigma[1] = 1;
        sigma[2] = 3;

        OrdinarySurfaceCoefficients lambda(_dimension[U], _dimension[V], sigma);

        (*@\Green{// coefficients appearing in the single (i.e., zeroth) seperable product of the coordinate function $s^0\left(u,v\right)$}@*)
        lambda(0, 0, U, 1)    =  1.0;       (*@\Green{// \hspace{0.24cm}$\Red{\mathbf{1}} \cdot \varphi_{6,1}\left(u\right) = \cos\left(u\right)$}@*)
        lambda(0, 0, U, 5)    = -1.0;       (*@\Green{// $\Red{\mathbf{-1}} \cdot \varphi_{6,5}\left(u\right) = -e^{\omega_0 u}\cos\left(u\right)$}@*)
        lambda(0, 0, V, 0)    =  5.0 / 4.0; (*@\Green{// \hspace{0.21cm}$\Red{\mathbf{\frac{5}{4}}} \cdot \varphi_{2,0}\left(v\right) = \frac{5}{4}$}@*)
        lambda(0, 0, V, 1)    =  1.0;       (*@\Green{// \hspace{0.26cm}$\Red{\mathbf{1}} \cdot \varphi_{2,1}\left(v\right) = \cos\left(v\right)$}@*)

        (*@\Green{// coefficients appearing in the single (i.e., zeroth) seperable product of the coordinate function $s^1\left(u,v\right)$}@*)
        lambda(1, 0, U, 2)    = -1.0;       (*@\Green{// \hspace{0.02cm}$\Red{\mathbf{-1}} \cdot \varphi_{6,2}\left(u\right) = -\sin\left(u\right)$}@*)
        lambda(1, 0, U, 6)    =  1.0;       (*@\Green{// \hspace{0.24cm}$\Red{\mathbf{1}} \cdot \varphi_{6,6}\left(u\right) = e^{\omega_0 u}\sin\left(u\right)$}@*)
        lambda(1, 0, V, 0)    =  5.0 / 4.0; (*@\Green{// \hspace{0.21cm}$\Red{\mathbf{\frac{5}{4}}} \cdot \varphi_{2,0}\left(v\right) = \frac{5}{4}$}@*)
        lambda(1, 0, V, 1)    =  1.0;       (*@\Green{// \hspace{0.25cm}$\Red{\mathbf{1}} \cdot \varphi_{2,1}\left(v\right) = \cos\left(v\right)$}@*)

        (*@\Green{// coefficients appearing in the zeroth seperable product of the coordinate function $s^2\left(u,v\right)$}@*)
        lambda(2, 0, U, 0)    =  7.0;       (*@\Green{// \hspace{0.26cm}$\Red{\mathbf{7}} \cdot \varphi_{6,0}\left(u\right) = 7$}@*)
        lambda(2, 0, U, 4)    = -1.0;       (*@\Green{// \hspace{0.03cm}$\Red{\mathbf{-1}} \cdot \varphi_{6,4}\left(u\right) = e^{\omega_1 u}$}@*)
        lambda(2, 0, V, 0)    =  1.0;       (*@\Green{// \hspace{0.23cm}$\Red{\mathbf{1}} \cdot \varphi_{2,0}\left(v\right) = 1$}@*)

        (*@\Green{// coefficients appearing in the first seperable product of the coordinate function $s^2\left(u,v\right)$}@*)
        lambda(2, 1, U, 0)    = -1.0;       (*@\Green{// $\Red{\mathbf{-1}} \cdot \varphi_{6,0}\left(u\right) = -1$}@*)
        lambda(2, 1, V, 2)    =  1.0;       (*@\Green{// \hspace{0.24cm}$\Red{\mathbf{1}} \cdot \varphi_{2,2}\left(v\right) = \sin\left(v\right)$}@*)

        (*@\Green{// coefficients appearing in the second seperable product of the coordinate function $s^2\left(u,v\right)$}@*)
        lambda(2, 2, U, 3)    =  1.0;       (*@\Green{// \hspace{0.24cm}$\Red{\mathbf{1}} \cdot \varphi_{6,3}\left(u\right) = e^{\omega_0 u}$}@*)
        lambda(2, 2, V, 2)    =  1.0;       (*@\Green{// \hspace{0.24cm}$\Red{\mathbf{1}} \cdot \varphi_{2,2}\left(v\right) = \sin\left(v\right)$}@*)

        GLdouble u_min        =  7.0 * PI / 2.0,
                 u_max        = 57.0 * PI / 8.0;
        GLint    row_count    =  5;

        GLdouble tau          =  2.0; (*@\Green{// one can also try other values like 1.0, 1.25, 1.5, 1.75, etc.}@*)

        GLdouble v_min        = -PI / 3.0 / tau,
                 v_max        = TWO_PI + v_min;
        GLint    column_count = floor(3.0 * tau);

        _patches.resizeRows(row_count);
        _patches.resizeColumns(column_count);

        _img_patches.resizeRows(row_count);
        _img_patches.resizeColumns(column_count);

        _u_isoparametric_lines.resizeRows(row_count);
        _u_isoparametric_lines.resizeColumns(column_count);

        _v_isoparametric_lines.resizeRows(row_count);
        _v_isoparametric_lines.resizeColumns(column_count);

        Matrix<SP<RowMatrix< SP<GenericCurve3>::Default> >::Default>
                *matrix_of_isoparametric_lines[2] =
                    {&_u_isoparametric_lines, &_v_isoparametric_lines};

        GLdouble u_step = (u_max - u_min) / row_count;
        GLdouble v_step = (v_max - v_min) / column_count;

        for (GLint k = 0; k < row_count; k++)
        {
            GLdouble u = u_min + k * u_step;

            (*@\Green{// the exponential-trigonometric EC space $\mathbb{ET}_{6}^{\alpha_0^k, \beta_0^k}$}@*)
            SnailUECSpace ET(u, u + u_step,                           (*@\Green{// $\alpha_0^k, ~\beta_0^k$}@*)
                             check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits);

            for (GLint l = 0; l < column_count; l++)
            {
                GLdouble v = v_min + l * v_step;

                (*@\Green{// the first order pure trigonometric EC space $\mathbb{T}_{2}^{\alpha_1^{\ell}, \beta_1^{\ell}}$}@*)
                SnailVECSpace T(v, v + v_step,                        (*@\Green{// $\alpha_1^{\ell}, ~\beta_1^{\ell}$}@*)
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits);

                _patches(k, l) = SP<BSurface3>::Default(new BSurface3(ET, T));

                if (!_patches(k, l))
                {
                    throw Exception("Could not create the B-surface patch in row " +
                                    toString(k) + " and column " + toString(l) + "!");
                }

                BSurface3 &surface = *_patches(k, l);

                if (!surface.updateControlPointsForExactDescription(lambda))
                {
                    throw Exception("Could not perform the B-representation of the "
                                    "ordinary surface patch in row " + toString(k) +
                                    " and column " + toString(l) + "!");
                }

                if (!surface.updateVertexBufferObjectsOfData())
                {
                    throw Exception("Could not update the VBO of the control net of the "
                                    "B-surface patch in row " + toString(k) +
                                    " and column " + toString(l) + "!");
                }

                _img_patches(k, l) = SP<TriangleMesh3>::Default(
                            surface.generateImage(
                                _surface_div_point_count[U], _surface_div_point_count[V],
                                _color_scheme));

                if (!_img_patches(k, l))
                {
                    throw Exception("Could not generate the image of the B-surface in "
                                    "row " + toString(k) + " and column " +
                                    toString(l) + "!");
                }

                if (!_img_patches(k, l)->updateVertexBufferObjects())
                {
                    throw Exception("Could not update the VBOs of the image of the "
                                    "B-surface in row " + toString(k) + " and column " +
                                    toString(l) + "!");
                }

                for (GLint direction = U; direction <= V; direction++)
                {
                    SP< RowMatrix<SP<GenericCurve3>::Default> >::Default
                        &sp_isoparametric_lines =
                            (*matrix_of_isoparametric_lines[direction])(k, l);

                    sp_isoparametric_lines =
                            SP< RowMatrix<SP<GenericCurve3>::Default> >::Default(
                                surface.generateIsoparametricLines(
                                    (variable::Type)direction,
                                    _isoparametric_line_count[direction],
                                    _maximum_order_of_derivatives[direction],
                                    _curve_div_point_count[direction]));

                    if (!sp_isoparametric_lines)
                    {
                        throw Exception("Could not generate the " +
                                        string(direction ? "v" : "u") + "-directional "
                                        "isoparametric lines of the B-surface in row " +
                                        toString(k) + " and column " + toString(l) + "!");
                    }

                    RowMatrix<SP<GenericCurve3>::Default>
                            &isoparametric_lines = *sp_isoparametric_lines;

                    for (GLint i = 0; i < isoparametric_lines.columnCount(); i++)
                    {
                        if (!isoparametric_lines[i]->updateVertexBufferObjects())
                        {
                            string th;

                            switch (i % 4)
                            {
                            case 1:  th = "st"; break;
                            case 2:  th = "nd"; break;
                            case 3:  th = "rd"; break;
                            default: th = "th"; break;
                            }

                            throw Exception("Could not update the VBOs of the " +
                                            toString(i) + th + " " +
                                            string(direction ? "v" : "u") +
                                            "-directional isoparametric line of "
                                            "the B-surface in row " + toString(k) +
                                            " and column " + toString(l) + "!");
                        }
                    }
                }
            }
        }

        _apply_reflection_lines                 = false; (*@\Green{// synchronize it with a check box}@*)
        _show_patches                           = true;  (*@\Green{// synchronize it with a check box}@*)
        _show_control_nets                      = true;  (*@\Green{// synchronize it with a check box}@*)
        _show_u_isoparametric_lines             = false; (*@\Green{// synchronize it with a check box}@*)
        _show_v_isoparametric_lines             = false; (*@\Green{// synchronize it with a check box}@*)
        _show_tangents_of_u_isoparametric_lines = false; (*@\Green{// synchronize it with a check box}@*)
        _show_tangents_of_v_isoparametric_lines = false; (*@\Green{// synchronize it with a check box}@*)
        _transparency                           = 0.0f;  (*@\Green{// synchronize it with a double spin box}@*)

        glClearColor(1.0, 1.0, 1.0, 1.0);          (*@\Green{// set the background color}@*)
        glEnable(GL_DEPTH_TEST);                   (*@\Green{// enable depth testing}@*)
        glEnable(GL_LINE_SMOOTH);                  (*@\Green{// enable the anti-aliasing of line primitives}@*)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_POINT_SMOOTH);                 (*@\Green{// enable the anti-aliasing of point primitives}@*)
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    }
    catch (Exception &e)
    {
        cout << e << endl;
    }
}

GLboolean GLWidget::_updateTransformationMatrices()
{
    (*@\Green{// \ldots}@*)
    (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms:_updateTransformationMatrices:implementation:start}--\mref{src:GLWidgetTestingShaderPrograms:_updateTransformationMatrices:implementation:end} of Listing \mref{src:GLWidgetTestingShaderPrograms.cpp}}@*)
    (*@\Green{// \ldots}@*)
}

(*@\Green{// your rendering method}@*)
void YourGLWidget::render()
{
    (*@\Green{// clear the color and depth buffers}@*)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    (*@\Green{// revert to the original transformation matrices in case of the color shader program}@*)
    _color_shader.enable();
        _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
    _color_shader.disable();

    (*@\Green{// At first, we render non-transparent geometries like control nets, spheres that represent control points,}@*)
    (*@\Green{// isoparametric lines and their tangent vectors:}@*)

    bool show_isoparametric_lines[2] =
        {_show_u_isoparametric_lines, _show_v_isoparametric_lines};

    bool show_tangents_of_isoparametric_lines[2] =
        {_show_tangents_of_u_isoparametric_lines,
         _show_tangents_of_v_isoparametric_lines};

    const Color4 *isoparametric_line_color[2] = {&colors::dark_purple, &colors::blue};

    const Color4 *isoparametric_line_tangent_color[2] = {&colors::purple, &colors::cyan};

    const Matrix<SP<RowMatrix< SP<GenericCurve3>::Default> >::Default>
            *matrix_of_isoparametric_lines[2] =
                {&_u_isoparametric_lines, &_v_isoparametric_lines};

    for (GLint direction = U; direction <= V; direction++)
    {
        (*@\Green{// try to render the isoparametric lines in the selected direction\ldots}@*)
        if (show_isoparametric_lines[direction] &&
            matrix_of_isoparametric_lines[direction])
        {
            const Matrix<SP<RowMatrix< SP<GenericCurve3>::Default> >::Default>
                    &matrix = *matrix_of_isoparametric_lines[direction];

            _color_shader.enable();

                glLineWidth(3.0);

                for (GLint k = 0; k < matrix.rowCount(); k++)
                {
                    for (GLint l = 0; l < matrix.columnCount(); l++)
                    {
                        if (matrix(k, l))
                        {
                            const RowMatrix<SP<GenericCurve3>::Default>
                                    &isoparametric_lines = *matrix(k, l);

                            _color_shader.setUniformColor(
                                        "color", *isoparametric_line_color[(k + l) % 2]);

                            for (GLint i = 0; i < isoparametric_lines.columnCount(); i++)
                            {
                                if (isoparametric_lines[i])
                                {
                                    isoparametric_lines[i]->renderDerivatives(
                                            0, GL_LINE_STRIP);
                                }
                            }
                        }
                    }
                }

                glLineWidth(1.0);

            _color_shader.disable();
        }

        (*@\Green{// try to render the tangent vectors of the isoparametric lines in the selected direction\ldots}@*)
        if (show_tangents_of_isoparametric_lines[direction] &&
            matrix_of_isoparametric_lines[direction])
        {
            const Matrix<SP<RowMatrix< SP<GenericCurve3>::Default> >::Default>
                    &matrix = *matrix_of_isoparametric_lines[direction];

            for (GLint k = 0; k < matrix.rowCount(); k++)
            {
                GLint index_offset = (k == matrix.rowCount() - 1);

                for (GLint l = 0; l < matrix.columnCount(); l++)
                {
                    if (matrix(k, l))
                    {
                        const RowMatrix<SP<GenericCurve3>::Default>
                                &isoparametric_lines = *matrix(k, l);

                        const Color4 &front_color =
                                *isoparametric_line_tangent_color[(k + l) % 2];

                        for (GLint i = 0;
                             i <= direction *
                                  (isoparametric_lines.columnCount() - 2 + index_offset);
                             i++)
                        {
                            if (isoparametric_lines[i])
                            {
                                const GenericCurve3 &curve = *isoparametric_lines[i];

                                _color_shader.enable();
                                    _color_shader.setUniformColor("color", front_color);
                                    curve.renderDerivatives(1, GL_LINES);
                                _color_shader.disable();

                                _renderSpheresAndConesAtEndpointsOfTangentVectors(
                                        curve, front_color, colors::silver);
                            }
                        }
                    }
                }
            }

        }
    }

    (*@\Green{// try to render the control nets of the B-surface patches\ldots}@*)
    if (_show_control_nets)
    {
        for (GLint k = 0; k < _patches.rowCount(); k++)
        {
            for (GLint l = 0; l < _patches.columnCount(); l++)
            {
                if (_patches(k, l))
                {
                    const BSurface3 &surface = *_patches(k, l);

                    if ((k + l) % 2)
                    {
                        _renderSpheresAtControlPoints(
                            surface, colors::ice_blue, colors::silver);
                    }
                    else
                    {
                        _renderSpheresAtControlPoints(
                            surface, colors::purple, colors::silver);
                    }

                    _renderControlNet(surface, colors::gray);
                }
            }
        }
    }

    (*@\Green{// Second, we render possible transparent geometries like the triangle mesh images of all B-surfaces:}@*)
    if (_show_patches)
    {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(30, 1.0);

        (*@\Green{// select the constant reference of either the two-sided lighting or the reflection line generating shader program\ldots}@*)
        const ShaderProgram &surface_shader = _apply_reflection_lines ?
                                              _reflection_lines : _two_sided_lighting;

        (*@\Green{// revert to the original transformation matrices\ldots}@*)
        surface_shader.enable();

            surface_shader.setUniformMatrix4fv("VM",  1, GL_FALSE, _VM.address());
            surface_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
            surface_shader.setUniformMatrix4fv("N",   1, GL_TRUE,  _tN.address());

        surface_shader.disable();

        (*@\Green{// in order to obtain the best coloring effects, color schemes different than the DEFAULT\_NULL\_FRAGMENT}@*)
        (*@\Green{// should be used together with black uniform (color) materials\ldots}@*)
        bool default_color_scheme =
                (_color_scheme == TensorProductSurface3::DEFAULT_NULL_FRAGMENT);

        for (GLint k = 0; k < _img_patches.rowCount(); k++)
        {
            for (GLint l = 0; l < _img_patches.columnCount(); l++)
            {
                if (_img_patches(k, l))
                {
                    const TriangleMesh3 &mesh = *_img_patches(k, l);

                    if ((k + l) % 2)
                    {
                        _renderTransparentMesh(
                            surface_shader, mesh,
                            default_color_scheme ? colors::baby_blue  : colors::black,
                            default_color_scheme ? colors::light_blue : colors::black,
                            _transparency);
                    }
                    else
                    {
                        _renderTransparentMesh(
                            surface_shader, mesh,
                            default_color_scheme ? colors::purple       : colors::black,
                            default_color_scheme ? colors::light_purple : colors::black,
                            _transparency);
                    }
                }
            }
        }

        glDisable(GL_POLYGON_OFFSET_FILL);
    }
}

(*@\Green{// auxiliary private method used for control point rendering}@*)
void YourGLWidget::_renderSpheresAtControlPoints(
        const BSurface3 &surface,
        const Color4 &front_color_material, const Color4 &back_color_material) const
{
    (*@\Green{// \ldots}@*)
    (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingBSurfaceOperations:_renderSpheresAtControlPoints:implementation:start}--\mref{src:GLWidgetTestingBSurfaceOperations:_renderSpheresAtControlPoints:implementation:end} of Listing \mref{src:GLWidgetTestingBSurfaceOperations.cpp}}@*)
    (*@\Green{// \ldots}@*)
}

(*@\Green{// auxiliary private method used for control net rendering}@*)
void YourGLWidget::_renderControlNet(
        const BSurface3 &surface, const Color4 &color, bool use_dashed_line) const
{
    (*@\Green{// \ldots}@*)
    (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingBSurfaceOperations:_renderControlNet:implementation:start}--\mref{src:GLWidgetTestingBSurfaceOperations:_renderControlNet:implementation:end} of Listing \mref{src:GLWidgetTestingBSurfaceOperations.cpp}}@*)
    (*@\Green{// \ldots}@*)
}

(*@\Green{// auxiliary private method used for transparent triangle mesh rendering}@*)
void YourGLWidget::_renderTransparentMesh(
        const ShaderProgram &shader, const TriangleMesh3 &mesh,
        const Color4 &front_color_material, const Color4 &back_color_material,
        GLfloat transparency) const
{
    (*@\Green{// \ldots}@*)
    (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingBSurfaceOperations:_renderTransparentMesh:implementation:start}--\mref{src:GLWidgetTestingBSurfaceOperations:_renderTransparentMesh:implementation:end} of Listing \mref{src:GLWidgetTestingBSurfaceOperations.cpp}}@*)
    (*@\Green{// \ldots}@*)
}

(*@\Green{// auxiliary private method used for tangent vector rendering}@*)
bool GLWidget::_renderSpheresAndConesAtEndpointsOfTangentVectors(
        const GenericCurve3 &curve,
        const Color4 &front_color_material, const Color4 &back_color_material) const
{
    if (curve.maximumOrderOfDerivatives() < 1)
    {
        return false;
    }

    _two_sided_lighting.enable();

        _two_sided_lighting.setUniformColorMaterial(
                "front_material", front_color_material);
        _two_sided_lighting.setUniformColorMaterial(
                "back_material", back_color_material);

        for (int i = 0; i < curve.pointCount(); i++)
        {
            const Cartesian3 &point = curve(0, i);

            Translate        sphere_T(point[0], point[1], point[2]);
            GLTransformation sphere_VM  = _VM * sphere_T * _cone_S;
            GLTransformation sphere_PVM = (*_P) * sphere_VM;

            _two_sided_lighting.setUniformMatrix4fv(
                    "VM",  1, GL_FALSE, sphere_VM.address());
            _two_sided_lighting.setUniformMatrix4fv(
                    "PVM", 1, GL_FALSE, sphere_PVM.address());
            _two_sided_lighting.setUniformMatrix4fv(
                    "N",   1, GL_TRUE,  sphere_VM.inverse().address());

            _sphere.render(_two_sided_lighting);

            const Cartesian3 &tangent = curve(1, i);

            Cartesian3 sum = point; sum += tangent;
            Cartesian3 K = tangent, I((GLdouble)rand() / (GLdouble)RAND_MAX,
                                      (GLdouble)rand() / (GLdouble)RAND_MAX,
                                      (GLdouble)rand() / (GLdouble)RAND_MAX);
            K.normalize();
            I.normalize();

            Cartesian3 J = K ^ I; J.normalize();

            I = J ^ K; I.normalize();

            GLTransformation C;

            C[ 0] =   I[0], C[ 1] =   I[1], C[ 2] =   I[2], C[ 3] = 0.0f;
            C[ 4] =   J[0], C[ 5] =   J[1], C[ 6] =   J[2], C[ 7] = 0.0f;
            C[ 8] =   K[0], C[ 9] =   K[1], C[10] =   K[2], C[11] = 0.0f;
            C[12] = sum[0], C[13] = sum[1], C[14] = sum[2], C[15] = 1.0f;

            GLTransformation cone_VM  = _VM * C * _cone_S;
            GLTransformation cone_PVM = (*_P) * cone_VM;

            _two_sided_lighting.setUniformMatrix4fv(
                    "VM",  1, GL_FALSE, cone_VM.address());
            _two_sided_lighting.setUniformMatrix4fv(
                    "PVM", 1, GL_FALSE, cone_PVM.address());
            _two_sided_lighting.setUniformMatrix4fv(
                    "N",   1, GL_TRUE,  cone_VM.inverse().address());

            _cone.render(_two_sided_lighting);
        }

    _two_sided_lighting.disable();

    return true;
}

(*@\Red{// The following comments specify the parameter settings that have to be applied in the constructor of the}@*)
(*@\Red{// class YourGLWidget in order to obtain the images shown in Figs.\ \mref{fig:control_point_based_exact_description_snail_tau_2}(\textit{a})--(\textit{b}) and \mref{fig:isoparametric_lines_user_manual}(\textit{a})--(c).}@*)

(*@\Green{// Fig.\ \mref{fig:control_point_based_exact_description_snail_tau_2}(\textit{a}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_patches = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_control\_nets = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_v\_isoparametric\_lines = false;}@*)

(*@\Green{// Fig.\ \mref{fig:control_point_based_exact_description_snail_tau_2}(\textit{b}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_patches = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_control\_nets = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_v\_isoparametric\_lines = false;}@*)

(*@\Green{// Fig.\ \mref{fig:isoparametric_lines_user_manual}(\textit{a}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_patches = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_control\_nets = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_isoparametric\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_isoparametric\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_v\_isoparametric\_lines = false;}@*)

(*@\Green{// Fig.\ \mref{fig:isoparametric_lines_user_manual}(\textit{b}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_patches = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_control\_nets = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_isoparametric\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_u\_isoparametric\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_v\_isoparametric\_lines = false;}@*)

(*@\Green{// Fig.\ \mref{fig:isoparametric_lines_user_manual}(\textit{c}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_patches = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_control\_nets = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_isoparametric\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_u\_isoparametric\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_tangents\_of\_v\_isoparametric\_lines = \Red{true};}@*)
