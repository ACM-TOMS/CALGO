#if defined(_MSC_VER)
#pragma warning(disable:4503)
#endif

#include "GLWidget.h"

#include <QTime>
#include <QMessageBox>

#include <iostream>
#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Math/Constants.h>
#include <Core/Utilities.h>
#include <Core/Geometry/Surfaces/Materials.h>
#include <Core/Math/RealMatrixDecompositions.h>
#include "../Spaces/SpecializedECSpaces.h"
#include <Core/Math/RealMatrices.h>
#include <fstream>

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
            // try to initialize the OpenGL Extension Wrangler library}\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:start
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
            _Rx.setAngle((GLfloat)(-27.0 * DEG_TO_RADIAN));
            _Rz.setAngle((GLfloat)(50.0 * DEG_TO_RADIAN));
            _S.setScalingFactors(1.9f, 1.9f, 1.9f);
            _T.setZDirectionalUnits(6.9f);

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
            } //(*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_color_shader:end

            if (!_two_sided_lighting.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX, "../../00_Dependencies/Shaders/two_sided_lighting.vert",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the vertex shader of two "
                                "sided lighting shader program!");
            }

            if (!_two_sided_lighting.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT, "../../00_Dependencies/Shaders/two_sided_lighting.frag",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the fragment shader of "
                                "two sided lighting shader program!");
            }

            if (!_two_sided_lighting.linkAttachedShaders(logging_is_enabled, shader_log))
            {
                throw Exception("Could not link the attached shaders of the "
                                "two sided lighting shader program!");
            } //(*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_two_sided_lighting:end

            if (!_reflection_lines.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX, "../../00_Dependencies/Shaders/reflection_lines.vert",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the vertex shader of the "
                                "reflection lines shader program!");
            }

            if (!_reflection_lines.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT, "../../00_Dependencies/Shaders/reflection_lines.frag",
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

            if (!_reflection_lines.setUniformValue1f("transparency", 0.5f))
            {
                throw Exception("Two-sided per pixed lighting: could not initialize the "
                                "uniform variable \"transparency\"!");
            }

            if (!_reflection_lines.setUniformValue1f("scale_factor", 9.7f))
            {
                throw Exception("Two-sided per pixed lighting: could not initialize the "
                                "uniform variable \"transparency\"!");
            }

            if (!_reflection_lines.setUniformValue1f("smoothing", 1.0f))
            {
                throw Exception("Reflection lines: could not initialize the "
                                "uniform variable \"smoothing\"!");
            }

            if (!_reflection_lines.setUniformValue1f("shading", 0.1f))
            {
                throw Exception("Reflection lines: could not initialize the "
                                "uniform variable \"shading\"!");
            }

            _reflection_lines.disable();

            if (!_sphere.loadFromOFF("../../00_Dependencies/Models/sphere.off"))
            {
                throw Exception("Could not load model file!");
            }

            _control_point_radius = 0.05;
            _sphere_S.setScalingFactors(_control_point_radius, _control_point_radius, _control_point_radius);

            if (!_sphere.updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBO of the model!");
            }

            GLboolean check_for_ill_conditioned_matrices  = GL_FALSE;
            GLint     expected_correct_significant_digits = 3;

            _alpha.resize(2);
            _beta.resize(2);
            _n.resize(2);
            _ratio.resize(2);
            _subdivision.resize(2);
            _img_subdivision.resize(2);
            _div_point_count.resize(2);

            _alpha[U]           =  0.0;
            _beta[U]            =  1.0;
            _n[U]               =  4;
            _ratio[U]           =  0.5;
            _div_point_count[U] = 50;

            _alpha[V]           =  -PI / 3.0;
            _beta[V]            =   PI / 3.0;
            _n[V]               =   2;
            _ratio[V]           =   0.5;
            _div_point_count[V] = 100;

            _space.resize(2);

            _space[U] = SP<ECSpace>::Default(
                        new (nothrow) PolynomialECSpace(
                            _alpha[U], _beta[U], _n[U],
                            check_for_ill_conditioned_matrices, expected_correct_significant_digits));

            if (!_space[U])
            {
                throw Exception("Could not create the EC space associated with the direction u!");
            }

            _space[V] = SP<ECSpace>::Default(
                        new (nothrow) TrigonometricECSpace(
                            _alpha[V], _beta[V], _n[V],
                            check_for_ill_conditioned_matrices, expected_correct_significant_digits));

            if (!_space[V])
            {
                throw Exception("Could not create the EC space associated with the direction v!");
            }

            _bsurface = SP<BSurface3>::Default(new (nothrow) BSurface3(*_space[U], *_space[V]));

            if (!_bsurface)
            {
                throw Exception("Could not create the EC B-surface!");
            }

            // define the control points
            GLdouble u_min = -3.0, u_max = 3.0, u_step = (u_max - u_min) / (_bsurface->rowCount() - 1);
            GLdouble v_min = -2.0, v_max = 2.0, v_step = (v_max - v_min) / (_bsurface->columnCount() - 1);
            GLdouble z_min = -3.5, z_max = 3.5, length = z_max - z_min;

            for (int r = 0; r < _bsurface->rowCount(); r++)
            {
                for (int c = 0; c < _bsurface->columnCount(); c++)
                {
                    Cartesian3 &reference = (*_bsurface)(r, c);

                    reference[0] = min(u_min + r * u_step, u_max);
                    reference[1] = min(v_min + c * v_step, v_max);
                    reference[2] = +(z_min + length * (GLdouble)rand() / (GLdouble)RAND_MAX);
                }
            }

            for (int r = 0; r < _bsurface->rowCount(); r++)
            {
                for (int c = 0; c < _bsurface->columnCount(); c++)
                {
                    Cartesian3 &reference = (*_bsurface)(r, c);

                    reference[2] += (GLdouble)rand() / (GLdouble)RAND_MAX < 0.5 ? -0.8 * (GLdouble)rand() / (GLdouble)RAND_MAX : +0.8 * (GLdouble)rand() / (GLdouble)RAND_MAX;
                }
            }

            // update the VBO of the control net
            if (!_bsurface->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the control net!");
            }

            // generate the image of the B-surface

            _color_scheme = TensorProductSurface3::DEFAULT_NULL_FRAGMENT;
            //_color_scheme = TensorProductSurface3::LOG_UMBILIC_DEVIATION_ENERGY_FRAGMENT;

            _img_bsurface = SP<TriangleMesh3>::Default(
                        _bsurface->generateImage(
                            _div_point_count[U], _div_point_count[V],
                            _color_scheme));

            if (!_img_bsurface)
            {
                throw Exception("Could not generate the image of the B-surface!");
            }

            if (!_img_bsurface->updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBOs of the B-surface's image!");
            }

            // perform order elevation in direction $u$
            CharacteristicPolynomial::Zero u_zero(0.0, 0.0, 9);

            _oe_bsurface = SP<BSurface3>::Default(
                        _bsurface->performOrderElevation(
                            U,  u_zero,
                            check_for_ill_conditioned_matrices, expected_correct_significant_digits));

            if (!_oe_bsurface)
            {
                throw Exception("Could not perform order elevation in direction u!");
            }

            // perform order elevation in direction $v$
            CharacteristicPolynomial::Zero v_zero(0.0, 1.0, 3);

            _oe_bsurface = SP<BSurface3>::Default(
                        _oe_bsurface->performOrderElevation(
                            V,  v_zero,
                            check_for_ill_conditioned_matrices, expected_correct_significant_digits));

            if (!_oe_bsurface)
            {
                throw Exception("Could not perform order elevation in direction v!");
            }

            // update the VBO of the control net of the order elevated B-surface
            if (!_oe_bsurface->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the control net of the order elevated B-surface!");
            }

            _oe_img_bsurface = SP<TriangleMesh3>::Default(
                        _oe_bsurface->generateImage(
                            _div_point_count[U], _div_point_count[V],
                            _color_scheme));

            if (!_oe_img_bsurface)
            {
                throw Exception("Could not generate the image of the order elevated B-surface!");
            }

            if (!_oe_img_bsurface->updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBOs of the order elevated B-surface's image!");
            }

            // perform subdivision
            for (GLint direction = U; direction <= V; direction++)
            {
                GLdouble gamma = (1.0 - _ratio[direction]) * _alpha[direction] +
                                        _ratio[direction]  * _beta[direction];

                _subdivision[direction] = SP< RowMatrix<SP<BSurface3>::Default> >::Default(
                            _bsurface->performSubdivision(
                                (variable::Type)direction, gamma,
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

                if (!_subdivision[direction])
                {
                    throw Exception("Could not perform subdivision in direction " +
                                    string(direction ? "v" : "u") + ".!");
                }

                _img_subdivision[direction].resizeColumns(2);

                for (GLint part = 0; part < 2; part++)
                {
                    if (!(*_subdivision[direction])[part]->updateVertexBufferObjectsOfData())
                    {
                        throw Exception("Could not update the VBOs of the control net of the " +
                                        string(part ? "right" : "left") + " part of the " +
                                        string(direction ? "v" : "u") + "-directional subdivisions!");
                    }

                    _img_subdivision[direction][part] = SP<TriangleMesh3>::Default(
                                (*_subdivision[direction])[part]->generateImage(
                                    _div_point_count[U], _div_point_count[V],
                                    _color_scheme));

                    if (!_img_subdivision[direction][part])
                    {
                        throw Exception("Could not generate the image of the " +
                                        string(part ? "right" : "left") + " part of the " +
                                        string(direction ? "v" : "u") + "-directional subdivisions!");
                    }

                    if (!_img_subdivision[direction][part]->updateVertexBufferObjects())
                    {
                        throw Exception("Could not update the VBOs of the image of the " +
                                        string(part ? "right" : "left") + " part of the " +
                                        string(direction ? "v" : "u") + "-directional subdivisions!");
                    }
                }
            }

            _apply_reflection_lines                    = false;
            _show_randomly_generated_initial_B_surface = true;
            _show_order_elevated_B_surface             = false;
            _show_u_subdivided_B_surfaces              = false;
            _show_v_subdivided_B_surfaces              = false;
            _compare_control_nets                      = false;

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
        // clears the color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // revert to the original transformation matrices in case of the color shader program
        _color_shader.enable();
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
        _color_shader.disable();

        // At first, we render non-transparent geometries like control nets and spheres
        // representing control points.

        if (_show_order_elevated_B_surface && _oe_bsurface)
        {
            // render small spheres at the positions of the control points of the order elevated B-surface
            _renderSpheresAtControlPoints(*_oe_bsurface, colors::light_purple, colors::silver);

            // render the control net of the order elevated B-surface
            _renderControlNet(*_oe_bsurface, colors::gray);
        }

        // render small spheres at the positions of the control points of the $u$- and $v$-directional subdivisions
        bool show_subdivision[2] = {_show_u_subdivided_B_surfaces, _show_v_subdivided_B_surfaces};
        const Color4 * const control_point_color[2] = {&colors::ice_blue, &colors::purple};

        if (_subdivision.size() == 2)
        {
            for (GLint direction = U; direction <= V; direction++)
            {
                if (show_subdivision[direction] && _subdivision[direction])
                {
                    for (GLint part = 0; part < 2; part++)
                    {
                        _renderSpheresAtControlPoints(*(*_subdivision[direction])[part],
                                                      *control_point_color[part], colors::silver);
                        _renderControlNet(*(*_subdivision[direction])[part], colors::gray);
                    }
                }
            }
        }

        // render small spheres at the positions of the initial control points
        if ((_show_randomly_generated_initial_B_surface || _compare_control_nets) &&
                _bsurface)
        {
            _renderSpheresAtControlPoints(*_bsurface, colors::red, colors::silver);
            _renderControlNet(*_bsurface, colors::gray, _compare_control_nets);
        }

        const ShaderProgram &surface_shader = _apply_reflection_lines ?
                                              _reflection_lines : _two_sided_lighting;

        // revert to the original transformation matrices
        surface_shader.enable();

            surface_shader.setUniformMatrix4fv("VM",  1, GL_FALSE, _VM.address());
            surface_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
            surface_shader.setUniformMatrix4fv("N",   1, GL_TRUE,  _tN.address());

        surface_shader.disable();

        bool default_color_scheme = (_color_scheme == TensorProductSurface3::DEFAULT_NULL_FRAGMENT);

        if (_show_randomly_generated_initial_B_surface && _img_bsurface)
        {
            _renderTransparentMesh(
                    surface_shader, *_img_bsurface,
                    default_color_scheme ? colors::light_blue : colors::black,
                    default_color_scheme ? colors::silver : colors::black);
        }

        if (_show_order_elevated_B_surface && _oe_bsurface && _oe_img_bsurface)
        {
            _renderTransparentMesh(
                    surface_shader, *_oe_img_bsurface,
                    default_color_scheme ? colors::light_blue : colors::black,
                    default_color_scheme ? colors::silver : colors::black);
        }

        const Color4 * const front_surface_color[2] = {&colors::light_blue, &colors::light_purple};
        const Color4 * const back_surface_color[2]  = {&colors::baby_blue,  &colors::dark_purple};

        if (_img_subdivision.size() == 2)
        {
            for (GLint direction = U; direction <= V; direction++)
            {
                if (show_subdivision[direction])
                {
                    for (GLint part = 0; part < 2; part++)
                    {
                        // render the images of the left and right patches
                        if (_img_subdivision[direction][part])
                        {
                            _renderTransparentMesh(
                                surface_shader, *_img_subdivision[direction][part],
                                default_color_scheme ? *front_surface_color[part] : colors::black,
                                default_color_scheme ? *back_surface_color[part]  : colors::black);
                        }
                    }
                }
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

            if (!shader.setUniformValue1f("transparency", 0.5f))
            {
                throw Exception("Currently active shader program: could not initialize the "
                                "uniform variable \"transparency\"!");
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
                mesh.render(shader);
            }

        shader.disable();
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

    void GLWidget::toggleReflectionLines(bool value)
    {
        if (_apply_reflection_lines != value)
        {
            _apply_reflection_lines = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfInitialBSurface(bool value)
    {
        static bool compare_control_nets = _compare_control_nets;

        if (_show_randomly_generated_initial_B_surface != value)
        {
            _show_randomly_generated_initial_B_surface = value;

            if (_show_randomly_generated_initial_B_surface)
            {
                compare_control_nets = _compare_control_nets;
                _compare_control_nets = false;
            }
            else
            {
                _compare_control_nets = compare_control_nets;
            }

            updateGL();
        }
    }

    void GLWidget::setVisibilityOfOrderElevatedBSurface(bool value)
    {
        if (_show_order_elevated_B_surface != value)
        {
            _show_order_elevated_B_surface = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfSubdivisionsInDirectionU(bool value)
    {
        if (_show_u_subdivided_B_surfaces != value)
        {
            _show_u_subdivided_B_surfaces = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfSubdivisionsInDirectionV(bool value)
    {
        if (_show_v_subdivided_B_surfaces != value)
        {
            _show_v_subdivided_B_surfaces = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfInitialControlNet(bool value)
    {
        if (_compare_control_nets != value)
        {
            _compare_control_nets = value;
            updateGL();
        }
    }
}
