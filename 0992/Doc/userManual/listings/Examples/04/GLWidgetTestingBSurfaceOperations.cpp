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
        _color_shader.enable();

        if (!_color_shader.setUniformColor("color", colors::gray))
        {
            throw Exception("Color shader program: could not initialize the "
                            "uniform variable \"color\"!");
        }

        _color_shader.disable();

        _two_sided_lighting.enable();

        if (!_two_sided_lighting.setUniformDirectionalLight(
                "light_source[0]", *_light))
        {
            throw Exception("Two-sided per pixed lighting: could not initialize the "
                            "uniform variable \"light_source[0]\"!");
        }

        if (!_two_sided_lighting.setUniformValue1i(
                "light_source[0].enabled", GL_TRUE))
        {
            throw Exception("Two-sided per pixed lighting: could not initialize the "
                            "uniform variable \"light_source[0].enabled\"!");
        }

        _two_sided_lighting.disable();

        _reflection_lines.enable();

        if (!_reflection_lines.setUniformDirectionalLight(
                "light_source[0]", *_light))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"light_source[0]\"!");
        }

        if (!_reflection_lines.setUniformValue1i("light_source[0].enabled", GL_TRUE))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"light_source[0].enabled\"!");
        }

        if (!_reflection_lines.setUniformValue1f("scale_factor", 9.7f))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"scale_factor\"!");
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

        (*@\Green{// Loading a triangulated unit sphere centered at the origin. Using translation and scaling transformations,}@*)
        (*@\Green{// it will be rendered multiple times at the positions of the control points.}@*)
        if (!_sphere.loadFromOFF("Models/sphere.off"))
        {
            throw Exception("Could not load model file!");
        }

        if (!_sphere.updateVertexBufferObjects())
        {
            throw Exception("Could not update the VBOs of the model!");
        }

        (*@\Green{// Determines the common scaling factor of the unit sphere that has to be rendered at the positions}@*)
        (*@\Green{// of the control points.}@*)
        _control_point_radius = 0.05;

        (*@\Green{// Determines the scaling transformation of the unit sphere that has to be rendered at the positions}@*)
        (*@\Green{// of the control points.}@*)
        _sphere_S.setScalingFactors(
                _control_point_radius, _control_point_radius, _control_point_radius);

        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:start}--\mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:end} of Listing \mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// Memory allocations and parameter settings\ldots}@*)
        _alpha.resize(2);
        _beta.resize(2);
        _n.resize(2);
        _ratio.resize(2);
        _subdivision.resize(2);
        _img_subdivision.resize(2);
        _div_point_count.resize(2);

        _alpha[U]           = 0.0;
        _beta[U]            = 1.0;
        _n[U]               = 4;
        _ratio[U]           = 0.5;
        _div_point_count[U] = 50;

        _alpha[V]           = -PI / 3.0;
        _beta[V]            =  PI / 3.0;
        _n[V]               =  2;
        _ratio[V]           =  0.5;
        _div_point_count[V] =  70;

        _space.resize(2);

        (*@\Green{// Try to create different EC spaces:}@*)
        _space[U] = SP<ECSpace>::Default(
                    new (nothrow) PolynomialECSpace(
                        _alpha[U], _beta[U], _n[U],
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits));

        if (!_space[U])
        {
            throw Exception("Could not create the EC space associated with the "
                            "direction u!");
        }

        _space[V] = SP<ECSpace>::Default(
                    new (nothrow) TrigonometricECSpace(
                        _alpha[V], _beta[V], _n[V],
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits));

        if (!_space[V])
        {
            throw Exception("Could not create the EC space associated with the "
                            "direction v!");
        }

        (*@\Green{// Try to create a B-surface:}@*)
        _bsurface = SP<BSurface3>::Default(
                new (nothrow) BSurface3(*_space[U], *_space[V]));

        if (!_bsurface)
        {
            throw Exception("Could not create the B-surface!");
        }

        (*@\Green{// Generate randomly the control points of the B-surface:}@*)
        GLdouble u_min = -3.0, u_max = 3.0,
                 u_step = (u_max - u_min) / (_bsurface->rowCount() - 1);
        GLdouble v_min = -2.0, v_max = 2.0,
                 v_step = (v_max - v_min) / (_bsurface->columnCount() - 1);
        GLdouble z_min = -3.5, z_max = 3.5, length = z_max - z_min;

        for (int r = 0; r < _bsurface->rowCount(); r++)
        {
            for (int c = 0; c < _bsurface->columnCount(); c++)
            {
                Cartesian3 &reference = (*_bsurface)(r, c);

                reference[0] = min(u_min + r * u_step, u_max);
                reference[1] = min(v_min + c * v_step, v_max);
                reference[2] = z_min + length * (GLdouble)rand()/(GLdouble)RAND_MAX;
            }
        }

        (*@\Green{// Try to update the VBO of the initial control net:}@*)
        if (!_bsurface->updateVertexBufferObjectsOfData())
        {
            throw Exception("Could not update the VBO of the initial control net!");
        }

        (*@\Green{// Select a color scheme:}@*)
        _color_scheme = TensorProductSurface3::DEFAULT_NULL_FRAGMENT; (*@\label{src:GLWidgetTestingBSurfaceOperations:color_scheme:initialization}@*)

        (*@\Green{// Try to generate the image of the initial B-surface:}@*)
        _img_bsurface = SP<TriangleMesh3>::Default(
                    _bsurface->generateImage(
                        _div_point_count[U], _div_point_count[V],
                        _color_scheme));

        if (!_img_bsurface)
        {
            throw Exception("Could not generate the image of the initial B-surface!");
        }

        (*@\Green{// Try to update the VBOs of initial B-surface's image:}@*)
        if (!_img_bsurface->updateVertexBufferObjects())
        {
            throw Exception("Could not update the VBOs of the initial B-surface's "
                            "image!");
        }

        (*@\Green{// Try to perform order elevation in direction $u$:}@*)
        CharacteristicPolynomial::Zero u_zero(0.0, 0.0, 9);

        _oe_bsurface = SP<BSurface3>::Default(
                    _bsurface->performOrderElevation(
                        U,  u_zero,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits));

        if (!_oe_bsurface)
        {
            throw Exception("Could not perform order elevation in direction u!");
        }

        (*@\Green{// Try to perform order elevation in direction $v$:}@*)
        CharacteristicPolynomial::Zero v_zero(0.0, 1.0, 3);

        _oe_bsurface = SP<BSurface3>::Default(
                    _oe_bsurface->performOrderElevation(
                        V,  v_zero,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits));

        if (!_oe_bsurface)
        {
            throw Exception("Could not perform order elevation in direction v!");
        }

        (*@\Red{\textbf{// Remark}}@*)
        (*@\Red{// After successful $u$- and $v$-directional order elevations, the obtained surface is described by means}@*)
        (*@\Red{// of the unique normalized B-bases of the EC spaces $\mathbb{P}_{8}^{0,1}=\operatorname{span}\left\{u^i:u\in\left[0,1\right]\right\}_{i=0}^{8}$ and}@*)
        (*@\Red{// $\mathbb{AT}_{8}^{-\frac{\pi}{3},\frac{\pi}{3}}=\operatorname{span}\left\{1,\cos\left(v\right),\sin\left(v\right),v \cos\left(v\right), v \sin\left(v\right), v^{2} \cos\left(v\right), v^{2} \sin\left(v\right), \cos\left(2 v\right), \sin\left(2 v\right) : v\in\left[-\frac{\pi}{3},\frac{\pi}{3}\right]\right\}$.}@*)

        (*@\Green{// Try to update the VBO of the control net of the order elevated B-surface:}@*)
        if (!_oe_bsurface->updateVertexBufferObjectsOfData())
        {
            throw Exception("Could not update the VBOs of the control net of the "
                            "order elevated B-surface!");
        }

        (*@\Green{// Try to generate the image of the order elevated B-surface:}@*)
        _oe_img_bsurface = SP<TriangleMesh3>::Default(
                    _oe_bsurface->generateImage(
                        _div_point_count[U], _div_point_count[V],
                        _color_scheme));

        if (!_oe_img_bsurface)
        {
            throw Exception("Could not generate the image of the order elevated "
                            "B-surface!");
        }

        (*@\Green{// Try to update the VBOs of the order elevated B-surface's image:}@*)
        if (!_oe_img_bsurface->updateVertexBufferObjects())
        {
            throw Exception("Could not update the VBOs of the order elevated "
                            "B-surface's image!");
        }

        (*@\Green{// Try to perform subdivision on the initial B-surface in both directions:}@*)
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
                    throw Exception("Could not update the VBOs of the control net of "
                                    "the " + string(part ? "right" : "left") +
                                    " part of the " + string(direction ? "v" : "u") +
                                    "-directional subdivisions!");
                }

                _img_subdivision[direction][part] = SP<TriangleMesh3>::Default(
                        (*_subdivision[direction])[part]->generateImage(
                                _div_point_count[U], _div_point_count[V],
                                _color_scheme));

                if (!_img_subdivision[direction][part])
                {
                    throw Exception("Could not generate the image of the " +
                                    string(part ? "right" : "left") + " part of the " +
                                    string(direction ? "v" : "u") +
                                    "-directional subdivisions!");
                }

                if (!_img_subdivision[direction][part]->updateVertexBufferObjects())
                {
                    throw Exception("Could not update the VBOs of the image of the " +
                                    string(part ? "right" : "left") + " part of the " +
                                    string(direction ? "v" : "u") +
                                    "-directional subdivisions!");
                }
            }
        }

        (*@\Green{// decides whether the two-sided lighting or the reflection lines shader program should be used}@*)
        (*@\Green{// during rendering}@*)
        _apply_reflection_lines                    = false; (*@\label{src:GLWidgetTestingBSurfaceOperations:boolean_values:initialization:start}@*)

        (*@\Green{// visibility flags that should be synchronized with radio buttons and check boxes of your graphical}@*)
        (*@\Green{// user interface}@*)
        _show_randomly_generated_initial_B_surface = true;  (*@\Green{// synchronize it with a radio button}@*)
        _show_order_elevated_B_surface             = false; (*@\Green{// synchronize it with a radio button}@*)
        _show_u_subdivided_B_surfaces              = false; (*@\Green{// synchronize it with a radio button}@*)
        _show_v_subdivided_B_surfaces              = false; (*@\Green{// synchronize it with a radio button}@*)
        _compare_control_nets                      = false; (*@\Green{// synchronize it with a check box}@*) (*@\label{src:GLWidgetTestingBSurfaceOperations:boolean_values:initialization:end}@*)

        glClearColor(1.0, 1.0, 1.0, 1.0);        (*@\Green{// set the background color}@*)
        glEnable(GL_DEPTH_TEST);                 (*@\Green{// enable depth testing}@*)
        glEnable(GL_LINE_SMOOTH);                (*@\Green{// enable the anti-aliasing of line primitives}@*)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_POINT_SMOOTH);               (*@\Green{// enable the anti-aliasing of point primitives}@*)
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

    (*@\Green{// At first, we render non-transparent geometries like control nets and spheres that represent control points:}@*)

    if (_show_order_elevated_B_surface && _oe_bsurface)
    {
        (*@\Green{// render small spheres at the positions of the control points of the order elevated B-surface\ldots}@*)
        _renderSpheresAtControlPoints(*_oe_bsurface,
                                      colors::light_purple, colors::silver);

        (*@\Green{// render the control net of the order elevated B-surface\ldots}@*)
        _renderControlNet(*_oe_bsurface, colors::gray);
    }

    (*@\Green{// auxiliary arrays that store boolean values or pointers to existing color objects\ldots}@*)
    (*@\Green{// these values will be used in the next for-loop\ldots}@*)
    bool show_subdivision[2] = {_show_u_subdivided_B_surfaces,
                                _show_v_subdivided_B_surfaces};

    const Color4 * const control_point_color[2] = {&colors::ice_blue, &colors::purple};

    for (GLint direction = U; direction <= V; direction++)
    {
        if (show_subdivision[direction] && _subdivision[direction])
        {
            for (GLint part = 0; part < 2; part++)
            {
                (*@\Green{// render small spheres at the positions of the control points of the $u$- and $v$-directional}@*)
                (*@\Green{// subdivisions\ldots}@*)
                _renderSpheresAtControlPoints(
                        *(*_subdivision[direction])[part],
                        *control_point_color[part], colors::silver);

                (*@\Green{// render the control nets of the $u$- and $v$-directional subdivisions\ldots}@*)
                _renderControlNet(*(*_subdivision[direction])[part], colors::gray);
            }
        }
    }

    if ((_show_randomly_generated_initial_B_surface || _compare_control_nets) &&
        _bsurface)
    {
        (*@\Green{// render small spheres at the positions of the control points of the initial B-surface\ldots}@*)
        _renderSpheresAtControlPoints(*_bsurface, colors::red, colors::silver);

        (*@\Green{// render the control net of the initial B-surface... based on the true or false values of the variable}@*)
        (*@\Green{// \_compare\_control\_nets we will use either dashed or continuous line styles, respectively\ldots}@*)
        _renderControlNet(*_bsurface, colors::gray, _compare_control_nets);
    }

    (*@\Green{// Second, we render transparent geometries like the triangle mesh images of all B-surfaces:}@*)

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

    if (_show_randomly_generated_initial_B_surface && _img_bsurface)
    {
        (*@\Green{// render the triangle mesh of the initial B-surface\ldots}@*)
        _renderTransparentMesh(
                surface_shader, *_img_bsurface,
                default_color_scheme ? colors::light_blue : colors::black,
                default_color_scheme ? colors::silver : colors::black);
    }

    if (_show_order_elevated_B_surface && _oe_bsurface && _oe_img_bsurface)
    {
        (*@\Green{// render the triangle mesh of the order elevated B-surface\ldots}@*)
        _renderTransparentMesh(
                surface_shader, *_oe_img_bsurface,
                default_color_scheme ? colors::light_blue : colors::black,
                default_color_scheme ? colors::silver : colors::black);
    }

    (*@\Green{// an auxiliary array of pointers to existing color objects that will be used in the next foor-loop\ldots}@*)
    const Color4 * const surface_color[2] = {&colors::light_blue, &colors::light_purple};

    for (GLint direction = U; direction <= V; direction++)
    {
        if (show_subdivision[direction])
        {
            for (GLint part = 0; part < 2; part++)
            {
                (*@\Green{// render the triangle meshes of the left and right B-patches obtained during $u$- and $v$-directional}@*)
                (*@\Green{// subdivisions\ldots}@*)
                if (_img_subdivision[direction][part])
                {
                    _renderTransparentMesh(
                            surface_shader, *_img_subdivision[direction][part],
                            default_color_scheme ? *surface_color[part] : colors::black,
                            default_color_scheme ? colors::silver : colors::black);
                }
            }
        }
    }

    (*@\Green{// render other geometries\ldots}@*)
}

(*@\Green{// auxiliary private method used for control point rendering}@*)
void YourGLWidget::_renderSpheresAtControlPoints(
        const BSurface3 &surface,
        const Color4 &front_color_material, const Color4 &back_color_material) const
{
    _two_sided_lighting.enable(); (*@\label{src:GLWidgetTestingBSurfaceOperations:_renderSpheresAtControlPoints:implementation:start}@*)

        _two_sided_lighting.setUniformColorMaterial(
                "front_material", front_color_material);

        _two_sided_lighting.setUniformColorMaterial(
                "back_material", back_color_material);

        for (int r = 0; r < surface.rowCount(); r++)
        {
            for (int c = 0; c < surface.columnCount(); c++)
            {
                const Cartesian3 &reference = surface(r, c);

                Translate sphere_T(reference[0], reference[1], reference[2]);

                GLTransformation sphere_VM  = _VM * sphere_T * _sphere_S;
                GLTransformation sphere_PVM = (*_P) * sphere_VM;

                _two_sided_lighting.setUniformMatrix4fv(
                        "VM",  1, GL_FALSE, sphere_VM.address());
                _two_sided_lighting.setUniformMatrix4fv(
                        "PVM", 1, GL_FALSE, sphere_PVM.address());
                _two_sided_lighting.setUniformMatrix4fv(
                        "N",   1, GL_TRUE,  sphere_VM.inverse().address());

                _sphere.render(_two_sided_lighting);
            }
        }

    _two_sided_lighting.disable(); (*@\label{src:GLWidgetTestingBSurfaceOperations:_renderSpheresAtControlPoints:implementation:end}@*)
}

(*@\Green{// auxiliary private method used for control net rendering}@*)
void YourGLWidget::_renderControlNet(
        const BSurface3 &surface, const Color4 &color, bool use_dashed_line) const
{
    glLineWidth(2.0); (*@\label{src:GLWidgetTestingBSurfaceOperations:_renderControlNet:implementation:start}@*)
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
    glLineWidth(1.0); (*@\label{src:GLWidgetTestingBSurfaceOperations:_renderControlNet:implementation:end}@*)
}

(*@\Green{// auxiliary private method used for transparent triangle mesh rendering}@*)
void YourGLWidget::_renderTransparentMesh(
        const ShaderProgram &shader, const TriangleMesh3 &mesh,
        const Color4 &front_color_material, const Color4 &back_color_material,
        GLfloat transparency) const
{
    shader.enable(); (*@\label{src:GLWidgetTestingBSurfaceOperations:_renderTransparentMesh:implementation:start}@*)

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
            (*@\Green{// in order to ensure transparency on white background, we render the geometry in two steps:}@*)
            glEnable(GL_BLEND);

                glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
                glDepthMask(GL_FALSE);
                glEnable(GL_CULL_FACE);

                    (*@\Green{// at first we cull the front faces and render the back ones,}@*)
                    glCullFace(GL_FRONT);

                    if (!mesh.render(shader))
                    {
                        throw Exception("YourGLWidget:_renderTransparentMesh : "
                                        "could not render the given triangle mesh!");
                    }

                    (*@\Green{// then we cull the back faces and render the front ones}@*)
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

    shader.disable(); (*@\label{src:GLWidgetTestingBSurfaceOperations:_renderTransparentMesh:implementation:end}@*)
}

(*@\Red{// The following comments specify the parameter settings that have to be applied in the constructor of the class}@*)
(*@\Red{// YourGLWidget in order to obtain the images presented in Figs.\ \mref{fig:B_surface_order_elevation}(\textit{a})--(\textit{f}), \mref{fig:color_schemes}(\textit{a})--($\ell$), \mref{fig:B_surface_subdivision_u}(\textit{a})--(\textit{d}) and}@*)
(*@\Red{// \mref{fig:B_surface_subdivision_v}(\textit{a})--(\textit{d}) of the current subsection. One has to modify the lines \mref{src:GLWidgetTestingBSurfaceOperations:color_scheme:initialization} and \mref{src:GLWidgetTestingBSurfaceOperations:boolean_values:initialization:start}--\mref{src:GLWidgetTestingBSurfaceOperations:boolean_values:initialization:end} of the constructor.}@*)

(*@\Green{// Figs.\ \mref{fig:B_surface_order_elevation}(\textit{a}), \mref{fig:B_surface_subdivision_u}(\textit{a}) and \mref{fig:B_surface_subdivision_v}(\textit{a}) are the same and illustrate the randomly generated initial B-surface.}@*)
(*@\Green{// The image was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Figs.\ \mref{fig:B_surface_order_elevation}(\textit{b})--(\textit{f}) and \mref{fig:color_schemes}(\textit{a})--($\ell$) illustrate the same order elevated B-surface.}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_order_elevation}(\textit{b}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = \Red{true};}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_order_elevation}(\textit{c}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_order_elevation}(\textit{d}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_order_elevation}(\textit{e}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::LOG\_UMBILIC\_DEVIATION\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_order_elevation}(\textit{f}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::LOG\_UMBILIC\_DEVIATION\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Figs.\ \mref{fig:color_schemes}(\textit{a})--($\ell$) were obtained by using the common parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)
(*@\Green{// and by applying the color schemes:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::X\_VARIATION\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::Y\_VARIATION\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::Z\_VARIATION\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::NORMAL\_LENGTH\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::GAUSSIAN\_CURVATURE\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::MEAN\_CURVATURE\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::WILLMORE\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::LOG\_WILLMORE\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::UMBILIC\_DEVIATION\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::TOTAL\_CURVATURE\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::LOG\_TOTAL\_CURVATURE\_ENERGY\_FRAGMENT};}@*)
(*@\Green{// respectively.}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_subdivision_u}(\textit{b}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = \Red{true};}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_subdivision_u}(\textit{c}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_subdivision_u}(\textit{d}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_subdivision_v}(\textit{b}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = \Red{true};}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_subdivision_v}(\textit{c}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)

(*@\Green{// Fig.\ \mref{fig:B_surface_subdivision_v}(\textit{d}) was obtained by using the parameter settings:}@*)
(*@\Green{// \hspace{1cm}\_color\_scheme = \Red{TensorProductSurface3::DEFAULT\_NULL\_FRAGMENT};}@*)
(*@\Green{// \hspace{1cm}\_apply\_reflection\_lines = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_show\_randomly\_generated\_initial\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_order\_elevated\_B\_surface = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_u\_subdivided\_B\_surfaces = false;}@*)
(*@\Green{// \hspace{1cm}\_show\_v\_subdivided\_B\_surfaces = \Red{true};}@*)
(*@\Green{// \hspace{1cm}\_compare\_control\_nets = false;}@*)
