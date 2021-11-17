#if defined(_MSC_VER)
#pragma warning(disable:4503)
#endif

#include "GLWidget.h"

#include <QMessageBox>

#include <iostream>
#include <fstream>

#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Math/Constants.h>
#include <Core/Geometry/Surfaces/Materials.h>
#include "../Spaces/SpecializedECSpaces.h"

using namespace std;
using namespace cagd::variable;

namespace cagd
{
    //--------------------------------
    // special and default constructor
    //--------------------------------
    GLWidget::GLWidget(QWidget *parent, const QGLFormat &format): QGLWidget(format, parent)
    {
    }

    //--------------------------------------------------------------------------------------
    // this virtual function is called once before the first call to paintGL() or resizeGL()
    //--------------------------------------------------------------------------------------
    void GLWidget::initializeGL()
    {
        try
        {
            // try to initialize the OpenGL Extension Wrangler library
            if (glewInit() != GLEW_OK)
            {
                throw Exception("Could not initialize the "
                                "OpenGL Extension Wrangler Library!");
            }

            // test whether your platform is supported from the perspective of the proposed function library
            if (!platformIsSupported())
            {
                throw Exception("The platform is not supported!\n\n"
                                "We assume that the user has: \n"
                                "- a multi-core CPU; and \n"
                                "- a graphics adapter managed by a driver that supports at least OpenGL 3.0.");
            }

            // creating an orthogonal perspective projection matrix
            _aspect = (float)width() / (float)height();
            _left   = _bottom = -10.0f;
            _right  = _top    = +10.0f;
            _near   = -20.0f;
            _far    = +20.0f;

            _P = SP<OrthogonalProjection>::Default(
                        new (nothrow) OrthogonalProjection(
                                _aspect, _left, _right, _bottom, _top, _near, _far));

            if (!_P)
            {
                throw Exception("Could not create the orthogonal projection matrix!");
            }

            // creating a view (or world) transformation matrix
            _eye[0]    = _eye[1]    = 0.0, _eye[2]    = 6.0;
            _center[0] = _center[1] =      _center[2] = 0.0;
            _up[0]     = _up[2]     = 0.0, _up[1]     = 1.0;

            _V = SP<LookAt>::Default(new (nothrow) LookAt(_eye, _center, _up));

            if (!_V)
            {
                throw Exception("Could not create the view/world transformation matrix!");
            }

            // specifying the axes of rotation matrices
            _Rx.setDirection(Cartesian3(1.0, 0.0, 0.0));
            _Ry.setDirection(Cartesian3(0.0, 1.0, 0.0));
            _Rz.setDirection(Cartesian3(0.0, 0.0, 1.0));

            // By default, all rotation angles, translation units and scaling factors are set to 0, 0 and 1, respectively.
            // However, these parameters can easily be modified, e.g.
            _Rx.setAngle((GLfloat)(-61.0 * DEG_TO_RADIAN));
            _Rz.setAngle((GLfloat)(-21.0 * DEG_TO_RADIAN));
            _S.setScalingFactors(1.6f, 1.6f, 1.6f);
            _T.setZDirectionalUnits(2.20f);

            // In what follows, we assume that the vertex and fragment shader files color.vert/frag,
            // two_sided_lighting.vert/frag and reflection_lines.vert/frag are downloaded and placed
            // in a folder named Shaders that is created along side of your executable.

            // If the variable logging_is_enabled is set to true, one obtains logging information at run-time about
            // creating, loading, compiling, attaching and linking different types of shaders.
            // If one is convinced that the applied shaders do not contain bugs, this logging mechanism can be deactivated.
            GLboolean logging_is_enabled = GL_FALSE;

            // By default, logging information appear in the standard console output, but they can be redirected
            // to any output streams as it is shown in the following lines.
            fstream   shader_log;

            if (logging_is_enabled)
            {
                shader_log.open("shader.log", ios_base::out);
            }

            if (!_color_shader.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX, "../../00_Dependencies/Shaders/color.vert",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the vertex shader of the "
                                "color shader program!");
            }

            if (!_color_shader.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT, "../../00_Dependencies/Shaders/color.frag",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the fragment shader of the "
                                "color shader program!");
            }

            if (!_color_shader.linkAttachedShaders(logging_is_enabled, shader_log))
            {
                throw Exception("Could not link the attached shaders of the "
                                "color shader program!");
            }

            if (!_two_sided_lighting.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX,
                        "../../00_Dependencies/Shaders/two_sided_lighting.vert",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the vertex shader of two "
                                "sided lighting shader program!");
            }

            if (!_two_sided_lighting.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT,
                        "../../00_Dependencies/Shaders/two_sided_lighting.frag",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the fragment shader of "
                                "two sided lighting shader program!");
            }

            if (!_two_sided_lighting.linkAttachedShaders(logging_is_enabled, shader_log))
            {
                throw Exception("Could not link the attached shaders of the "
                                "two sided lighting shader program!");
            }

            if (!_reflection_lines.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX,
                        "../../00_Dependencies/Shaders/reflection_lines.vert",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the vertex shader of the "
                                "reflection lines shader program!");
            }

            if (!_reflection_lines.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT,
                        "../../00_Dependencies/Shaders/reflection_lines.frag",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the fragment shader of the "
                                "reflection lines shader program!");
            }

            if (!_reflection_lines.linkAttachedShaders(logging_is_enabled, shader_log))
            {
                throw Exception("Could not link the attached shaders of the "
                                "reflection lines shader program!");
            }

            if (logging_is_enabled)
            {
                shader_log.close();
            }

            // Try to update all required transformation matrices (the method _updateTransformationMatrices() has to
            // be called after every successful transformation related event handling).
            if (!_updateTransformationMatrices())
            {
                throw Exception("Could not update all transformation matrices!");
            }

            // creating a directional light object

            Cartesian3 direction(0.0, 0.0, 1.0);
            Cartesian3 eye = _eye;
            eye.normalize();

            Cartesian3 half_vector = direction;
            half_vector += eye;

            half_vector.normalize();

            _light = SP<DirectionalLight>::Default(
                        new DirectionalLight(
                            Homogeneous3(0.0f, 0.0f, 1.0f, 0.0f), // direction vector of the light
                            Homogeneous3(half_vector),            // homogeneous half vector
                            Color4(0.4f, 0.4f, 0.4f, 1.0f),       // ambient light intensity
                            Color4(0.8f, 0.8f, 0.8f, 1.0f),       // diffuse light intensity
                            Color4(1.0f, 1.0f, 1.0f, 1.0f)));     // specular light intensity

            if (!_light)
            {
                throw Exception("Could not create the directional light object!");
            }

            // examples for communicating with our shader programs via uniform variables

            // set the default color of the shader program
            _color_shader.enable();

            if (!_color_shader.setUniformColor("color", colors::gray))
            {
                throw Exception("Color shader program: could not initialize the "
                                "uniform variable \"color\"!");
            }

            _color_shader.disable();

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

            if (!_sphere.loadFromOFF("../../00_Dependencies/Models/sphere.off"))
            {
                throw Exception("Could not load model file!");
            }

            if (!_sphere.updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBO of the model!");
            }

            _control_point_radius = 0.06875;
            _sphere_S.setScalingFactors(_control_point_radius, _control_point_radius, _control_point_radius);

            if (!_cone.loadFromOFF("../../00_Dependencies/Models/cone.off"))
            {
                throw Exception("Could not load model file!");
            }

            if (!_cone.updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBO of the model!");
            }
            _cone_S.setScalingFactors(_control_point_radius / 1.5f,
                                      _control_point_radius / 1.5f,
                                      _control_point_radius / 1.5f);

            GLboolean check_for_ill_conditioned_matrices  = GL_FALSE;
            GLint     expected_correct_significant_digits = 3;

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

            _dimension[V]                    =   3;
            _surface_div_point_count[V]      = 100;
            _isoparametric_line_count[V]     =   3;
            _maximum_order_of_derivatives[V] =   1;
            _curve_div_point_count[V]        =  13;

            _color_scheme = TensorProductSurface3::DEFAULT_NULL_FRAGMENT;

            std::vector<GLint> sigma(3);
            sigma[0] = 1;
            sigma[1] = 1;
            sigma[2] = 3;

            OrdinarySurfaceCoefficients lambda(_dimension[U], _dimension[V], sigma);

            lambda(0, 0, U, 1)    =  1.0;
            lambda(0, 0, U, 5)    = -1.0;

            lambda(0, 0, V, 0)    =  5.0 / 4.0;
            lambda(0, 0, V, 1)    =  1.0;

            lambda(1, 0, U, 2)    = -1.0;
            lambda(1, 0, U, 6)    =  1.0;

            lambda(1, 0, V, 0)    =  5.0 / 4.0;
            lambda(1, 0, V, 1)    =  1.0;

            lambda(2, 0, U, 0)    =  7.0;
            lambda(2, 0, U, 4)    = -1.0;
            lambda(2, 0, V, 0)    =  1.0;

            lambda(2, 1, U, 0)    = -1.0;
            lambda(2, 1, V, 2)    =  1.0;

            lambda(2, 2, U, 3)    =  1.0;
            lambda(2, 2, V, 2)    =  1.0;

            GLdouble u_min        =  7.0 * PI / 2.0,
                     u_max        = 57.0 * PI / 8.0;
            GLint    row_count    =  5;

            GLdouble tau          =  2.0;

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

                SnailUECSpace ET(u, u + u_step,
                                 check_for_ill_conditioned_matrices,
                                 expected_correct_significant_digits);

                for (GLint l = 0; l < column_count; l++)
                {
                    GLdouble v = v_min + l * v_step;

                    SnailVECSpace T(v, v + v_step,
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
                            &sp_isoparametric_lines = (*matrix_of_isoparametric_lines[direction])(k, l);

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

            _apply_reflection_lines                 = false;
            _show_patches                           = true;
            _show_control_nets                      = true;
            _show_u_isoparametric_lines             = false;
            _show_v_isoparametric_lines             = false;
            _show_tangents_of_u_isoparametric_lines = false;
            _show_tangents_of_v_isoparametric_lines = false;
            _transparency                           = false;

            glClearColor(1.0, 1.0, 1.0, 1.0);          // set the background color
            glEnable(GL_DEPTH_TEST);                   // enable depth testing
            glEnable(GL_LINE_SMOOTH);                  // enable the anti-aliasing of line primitives
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
            glEnable(GL_POINT_SMOOTH);                 // enable the anti-aliasing of point primitives
            glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        }
        catch (Exception &e)
        {
            cout << e << endl;
            QMessageBox::critical(this, "An exception occured...", QString::fromStdString(e.reason()));
            exit(0);
        }
    }

    GLboolean GLWidget::_updateTransformationMatrices()
    {
        if (_P && _V)
        {
            _M              = _Rx * _Ry * _Rz * _T * _S;
            _VM             = (*_V) * _M;
            _PVM            = (*_P) * _VM;

            bool invertible = false;
            _tN             = _VM.inverse(&invertible);

            return invertible;
        }

        return GL_FALSE;
    }

    //-----------------------
    // the rendering function
    //-----------------------
    void GLWidget::paintGL()
    {
        // clear the color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // revert to the original transformation matrices in case of the color shader program
        _color_shader.enable();
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
        _color_shader.disable();

        // At first, we render non-transparent geometries like control nets, spheres that represent control points,
        // isoparametric lines and their tangent vectors:

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
            // try to render the isoparametric lines in the selected direction\ldots
            if (show_isoparametric_lines[direction] &&
                matrix_of_isoparametric_lines[direction])
            {
                const Matrix<SP<RowMatrix< SP<GenericCurve3>::Default> >::Default>
                        &matrix = *matrix_of_isoparametric_lines[direction];

                _color_shader.enable();

                    glLineWidth(2.0);

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

            // try to render the tangent vectors of the isoparametric lines in the selected direction\ldots
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

        // try to render the control nets of the B-surface patches\ldots
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

                        _renderControlNet(surface, colors::gray, false);
                    }
                }
            }
        }

        // Second, we render possible transparent geometries like the triangle mesh images of all B-surfaces:
        if (_show_patches)
        {
            if (_show_u_isoparametric_lines || _show_v_isoparametric_lines)
            {
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(30, 1.0);
            }

            // select the constant reference of either the two-sided lighting or the reflection line generating shader program\ldots
            const ShaderProgram &surface_shader = _apply_reflection_lines ?
                                                  _reflection_lines : _two_sided_lighting;

            // revert to the original transformation matrices\ldots
            surface_shader.enable();

                surface_shader.setUniformMatrix4fv("VM",  1, GL_FALSE, _VM.address());
                surface_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
                surface_shader.setUniformMatrix4fv("N",   1, GL_TRUE,  _tN.address());

            surface_shader.disable();

            // in order to obtain the best coloring effects, color schemes different than the DEFAULT_NULL_FRAGMENT
            // should be used together with black uniform (color) materials\ldots
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
                                default_color_scheme ? colors::baby_blue : colors::black,
                                default_color_scheme ? colors::light_blue : colors::black,
                                _transparency);
                        }
                        else
                        {
                            _renderTransparentMesh(
                                surface_shader, mesh,
                                default_color_scheme ? colors::purple : colors::black,
                                default_color_scheme ? colors::light_purple : colors::black,
                                _transparency);
                        }
                    }
                }
            }

            if (_show_u_isoparametric_lines || _show_v_isoparametric_lines)
            {
                glDisable(GL_POLYGON_OFFSET_FILL);
            }
        }
    }


    void GLWidget::_renderSpheresAtControlPoints(const BSurface3 &surface, const Color4 &front_color_material, const Color4 &back_color_material) const
    {
        _two_sided_lighting.enable();

            _two_sided_lighting.setUniformColorMaterial("front_material", front_color_material);
            _two_sided_lighting.setUniformColorMaterial("back_material", back_color_material);

            for (int r = 0; r < surface.rowCount(); r++)
            {
                for (int c = 0; c < surface.columnCount(); c++)
                {
                    const Cartesian3 &reference = surface(r, c);

                    Translate sphere_T(reference[0], reference[1], reference[2]);

                    GLTransformation sphere_VM  = _VM * sphere_T * _sphere_S;
                    GLTransformation sphere_PVM = (*_P) * sphere_VM;

                    _two_sided_lighting.setUniformMatrix4fv("VM",  1, GL_FALSE, sphere_VM.address());
                    _two_sided_lighting.setUniformMatrix4fv("PVM", 1, GL_FALSE, sphere_PVM.address());
                    _two_sided_lighting.setUniformMatrix4fv("N",   1, GL_TRUE,  sphere_VM.inverse().address());

                    _sphere.render(_two_sided_lighting);
                }
            }

        _two_sided_lighting.disable();
    }

    void GLWidget::_renderControlNet(const BSurface3 &surface, const Color4 &color, bool use_dashed_line) const
    {
//        glLineWidth(2.0);
        _color_shader.enable();

            _color_shader.setUniformColor("color", color);

            if (use_dashed_line)
            {
                glEnable(GL_LINE_STIPPLE);
                    glLineStipple(1, 0xf0f0);
                    surface.renderData(_color_shader);
                glDisable(GL_LINE_STIPPLE);
            }
            else
            {
                surface.renderData(_color_shader);
            }

        _color_shader.disable();
//        glLineWidth(1.0);
    }

    void GLWidget::_renderTransparentMesh(
            const ShaderProgram &shader, const TriangleMesh3 &mesh,
            const Color4 &front_color_material, const Color4 &back_color_material,
            GLfloat transparency) const
    {
        shader.enable();

            if (!shader.setUniformColorMaterial("front_material", front_color_material))
            {
                throw Exception("Currently active shader program: could not initialize the "
                                "uniform variable \"front_material\"!");
            }

            if (!shader.setUniformColorMaterial("back_material",  back_color_material))
            {
                throw Exception("Currently active shader program: could not initialize the "
                                "uniform variable \"back_material\"!");
            }

            if (transparency)
            {
                // in order to ensure transparency on white background, we render the geometry in two steps:
                // at first we cull the front faces and render the back ones,
                // then we cull the back faces and render the front ones
                glEnable(GL_BLEND);

                    glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
                    glDepthMask(GL_FALSE);
                    glEnable(GL_CULL_FACE);

                        glCullFace(GL_FRONT);

                        if (!mesh.render(shader))
                        {
                            throw Exception("YourGLWidget:_renderTransparentMesh : could not "
                                            "render the given triangle mesh!");
                        }

                        glCullFace(GL_BACK);
                        mesh.render(shader);

                    glDisable(GL_CULL_FACE);
                    glDepthMask(GL_TRUE);

                glDisable(GL_BLEND);
            }
            else
            {
                if (!shader.setUniformValue1f("transparency", 0.0f))
                {
                    throw Exception("Currently active shader program: could not initialize the "
                                    "uniform variable \"transparency\"!");
                }

                if (!mesh.render(shader))
                {
                    throw Exception("YourGLWidget:_renderTransparentMesh : could not "
                                    "render the given triangle mesh!");
                }
            }

        shader.disable();
    }
    
    bool GLWidget::_renderSpheresAndConesAtEndpointsOfTangentVectors(const GenericCurve3 &curve, const Color4 &front_color_material, const Color4 &back_color_material) const
    {
        if (curve.maximumOrderOfDerivatives() < 1)
        {
            return false;
        }
        
        _two_sided_lighting.enable();

            _two_sided_lighting.setUniformColorMaterial("front_material", front_color_material);
            _two_sided_lighting.setUniformColorMaterial("back_material", back_color_material);

            for (int i = 0; i < curve.pointCount(); i++)
            {
                    const Cartesian3 &point = curve(0, i);

                    Translate sphere_T(point[0], point[1], point[2]);

                    GLTransformation sphere_VM  = _VM * sphere_T * _cone_S;
                    GLTransformation sphere_PVM = (*_P) * sphere_VM;

                    _two_sided_lighting.setUniformMatrix4fv("VM",  1, GL_FALSE, sphere_VM.address());
                    _two_sided_lighting.setUniformMatrix4fv("PVM", 1, GL_FALSE, sphere_PVM.address());
                    _two_sided_lighting.setUniformMatrix4fv("N",   1, GL_TRUE,  sphere_VM.inverse().address());

                    _sphere.render(_two_sided_lighting);
                    
                    const Cartesian3 &tangent = curve(1, i);
                    
                    Cartesian3 sum = point;
                    sum += tangent;

                    Cartesian3 K = tangent, I((GLdouble)rand() / (GLdouble)RAND_MAX,
                                              (GLdouble)rand() / (GLdouble)RAND_MAX,
                                              (GLdouble)rand() / (GLdouble)RAND_MAX);
                    K.normalize();
                    I.normalize();

                    Cartesian3 J = K ^ I;
                    J.normalize();

                    I = J ^ K;
                    I.normalize();

                    GLTransformation C;

                    C[ 0] =   I[0], C[ 1] =   I[1], C[ 2] =   I[2], C[ 3] = 0.0f;
                    C[ 4] =   J[0], C[ 5] =   J[1], C[ 6] =   J[2], C[ 7] = 0.0f;
                    C[ 8] =   K[0], C[ 9] =   K[1], C[10] =   K[2], C[11] = 0.0f;
                    C[12] = sum[0], C[13] = sum[1], C[14] = sum[2], C[15] = 1.0f;

                    GLTransformation cone_VM  = _VM * C * _cone_S;
                    GLTransformation cone_PVM = (*_P) * cone_VM;

                    _two_sided_lighting.setUniformMatrix4fv("VM",  1, GL_FALSE, cone_VM.address());
                    _two_sided_lighting.setUniformMatrix4fv("PVM", 1, GL_FALSE, cone_PVM.address());
                    _two_sided_lighting.setUniformMatrix4fv("N",   1, GL_TRUE,  cone_VM.inverse().address());

                    _cone.render(_two_sided_lighting);
            }

        _two_sided_lighting.disable();
        
        return true;
    }

    //----------------------------------------------------------------------------
    // when the main window is resized one needs to redefine the projection matrix
    //----------------------------------------------------------------------------
    void GLWidget::resizeGL(int w, int h)
    {
        // setting the new size of the rendering context
        glViewport(0, 0, w, h);

        if (_P)
        {
            _P->setAspectRatio((float)w / (float)h);

            _PVM = (*_P) * _VM;
            _tN = _VM.inverse();
        }

        updateGL();
    }

    //-----------------------------------
    // implementation of the public slots
    //-----------------------------------

    void GLWidget::setAngleX(int value)
    {
        _Rx.setAngle(value * DEG_TO_RADIAN);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setAngleY(int value)
    {
        _Ry.setAngle(value * DEG_TO_RADIAN);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setAngleZ(int value)
    {
        _Rz.setAngle(value * DEG_TO_RADIAN);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setZoomFactor(double value)
    {
        _S.setScalingFactors(value, value, value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setTransX(double value)
    {
        _T.setXDirectionalUnits(value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setTransY(double value)
    {
        _T.setYDirectionalUnits(value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setTransZ(double value)
    {
        _T.setZDirectionalUnits(value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::toggleReflectionLines(bool enabled)
    {
        if (_apply_reflection_lines != enabled)
        {
            _apply_reflection_lines = enabled;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfPatches(bool value)
    {
        if (_show_patches != value)
        {
            _show_patches = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfControlNets(bool value)
    {
        if (_show_control_nets != value)
        {
            _show_control_nets = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfUIsoparametricLines(bool value)
    {
        if (_show_u_isoparametric_lines != value)
        {
            _show_u_isoparametric_lines = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfVIsoparametricLines(bool value)
    {
        if (_show_v_isoparametric_lines != value)
        {
            _show_v_isoparametric_lines = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfTangentsOfUIsoparametricLines(bool value)
    {
        if (_show_tangents_of_u_isoparametric_lines != value)
        {
            _show_tangents_of_u_isoparametric_lines = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfTangentsOfVIsoparametricLines(bool value)
    {
        if (_show_tangents_of_v_isoparametric_lines != value)
        {
            _show_tangents_of_v_isoparametric_lines = value;
            updateGL();
        }
    }

    void GLWidget::setTransparency(bool value)
    {
        if (_transparency != value)
        {
            _transparency = value;

            if (_transparency)
            {
                _two_sided_lighting.enable();
                    _two_sided_lighting.setUniformValue1f("transparency", 0.5f);
                _two_sided_lighting.disable();

                _reflection_lines.enable();
                    _reflection_lines.setUniformValue1f("transparency", 0.5f);
                _reflection_lines.disable();
            }
            else
            {
                _two_sided_lighting.enable();
                    _two_sided_lighting.setUniformValue1f("transparency", 0.0f);
                _two_sided_lighting.disable();

                _reflection_lines.enable();
                    _reflection_lines.setUniformValue1f("transparency", 0.0f);
                _reflection_lines.disable();
            }

            updateGL();
        }
    }
}
