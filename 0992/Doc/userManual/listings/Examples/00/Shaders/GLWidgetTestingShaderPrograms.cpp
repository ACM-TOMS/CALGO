#include "YourGLWidget.h"

#include <iostream>
#include <fstream>

#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Geometry/Coordinates/Colors4.h>
#include <Core/Geometry/Surfaces/Materials.h>

using namespace cagd;
using namespace std;

(*@\Green{// your default/special constructor}@*)
YourGLWidget::YourGLWidget(/*...*/)
{
    try
    {
        (*@\Green{// try to initialize the OpenGL Extension Wrangler library}\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:start}@*)
        if (glewInit() != GLEW_OK)
        {
            throw Exception("Could not initialize the "
                            "OpenGL Extension Wrangler Library!");
        }

        (*@\Green{// test whether your platform is supported from the perspective of the proposed function library}@*)
        if (!platformIsSupported())
        {
            throw Exception("The platform is not supported!");
        }

        (*@\Green{// creating an orthogonal perspective projection matrix}@*)
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

        (*@\Green{// creating a view (or world) transformation matrix}@*)
        _eye[0]    = _eye[1]    = 0.0, _eye[2]    = 6.0;
        _center[0] = _center[1] =      _center[2] = 0.0;
        _up[0]     = _up[2]     = 0.0, _up[1]     = 1.0;

        _V = SP<LookAt>::Default(new (nothrow) LookAt(_eye, _center, _up));

        if (!_V)
        {
            throw Exception("Could not create the view/world transformation matrix!");
        }

        (*@\Green{// specifying the axes of rotation matrices}@*)
        _Rx.setDirection(Cartesian3(1.0, 0.0, 0.0));
        _Ry.setDirection(Cartesian3(0.0, 1.0, 0.0));
        _Rz.setDirection(Cartesian3(0.0, 0.0, 1.0));

        (*@\Green{// By default, all rotation angles, translation units and scaling factors are set to 0, 0 and 1, respectively.}@*)

        (*@\Green{// In what follows, we assume that the vertex and fragment shader files color.vert/frag,}@*)
        (*@\Green{// two\_sided\_lighting.vert/frag and reflection\_lines.vert/frag are downloaded and placed}@*)
        (*@\Green{// in a folder named Shaders that is created along side of your executable.}@*)

        (*@\Green{// If the variable logging\_is\_enabled is set to true, one obtains logging information at run-time about}@*)
        (*@\Green{// creating, loading, compiling, attaching and linking different types of shaders.}@*)
        (*@\Green{// If one is convinced that the applied shaders do not contain bugs, this logging mechanism can be deactivated.}@*)
        GLboolean logging_is_enabled = GL_TRUE; (*@\label{src:GLWidgetTestingShaderPrograms:loading_shaders:start}@*)

        (*@\Green{// By default, logging information appear in the standard console output, but they can be redirected}@*)
        (*@\Green{// to any output streams as it is shown in the following lines.}@*)
        fstream   shader_log("shader.log", ios_base::out); (*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:end}@*)

        if (!_color_shader.attachNewShaderFromSourceFile( (*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_color_shader:start}@*)
                    ShaderProgram::Shader::VERTEX, "Shaders/color.vert",
                    logging_is_enabled, shader_log))
        {
            throw Exception("Could not attach the vertex shader of the "
                            "color shader program!");
        }

        if (!_color_shader.attachNewShaderFromSourceFile(
                    ShaderProgram::Shader::FRAGMENT, "Shaders/color.frag",
                    logging_is_enabled, shader_log))
        {
            throw Exception("Could not attach the fragment shader of the "
                            "color shader program!");
        }

        if (!_color_shader.linkAttachedShaders(logging_is_enabled, shader_log))
        {
            throw Exception("Could not link the attached shaders of the "
                            "color shader program!");
        } (*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_color_shader:end}@*)

        if (!_two_sided_lighting.attachNewShaderFromSourceFile( (*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_two_sided_lighting:start}@*)
                    ShaderProgram::Shader::VERTEX, "Shaders/two_sided_lighting.vert",
                    logging_is_enabled, shader_log))
        {
            throw Exception("Could not attach the vertex shader of two "
                            "sided lighting shader program!");
        }

        if (!_two_sided_lighting.attachNewShaderFromSourceFile(
                    ShaderProgram::Shader::FRAGMENT, "Shaders/two_sided_lighting.frag",
                    logging_is_enabled, shader_log))
        {
            throw Exception("Could not attach the fragment shader of "
                            "two sided lighting shader program!");
        }

        if (!_two_sided_lighting.linkAttachedShaders(logging_is_enabled, shader_log))
        {
            throw Exception("Could not link the attached shaders of the "
                            "two sided lighting shader program!");
        } (*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_two_sided_lighting:end}@*)

        if (!_reflection_lines.attachNewShaderFromSourceFile( (*@\label{src:GLWidgetTestingShaderPrograms:loading_reflection_lines:start}@*)
                    ShaderProgram::Shader::VERTEX, "Shaders/reflection_lines.vert",
                    logging_is_enabled, shader_log))
        {
            throw Exception("Could not attach the vertex shader of the "
                            "reflection lines shader program!");
        }

        if (!_reflection_lines.attachNewShaderFromSourceFile(
                    ShaderProgram::Shader::FRAGMENT, "Shaders/reflection_lines.frag",
                    logging_is_enabled, shader_log))
        {
            throw Exception("Could not attach the fragment shader of the "
                            "reflection lines shader program!");
        }

        if (!_reflection_lines.linkAttachedShaders(logging_is_enabled, shader_log))
        {
            throw Exception("Could not link the attached shaders of the "
                            "reflection lines shader program!");
        } (*@\label{src:GLWidgetTestingShaderPrograms:loading_reflection_lines:end}@*)

        shader_log.close(); (*@\label{src:GLWidgetTestingShaderPrograms:loading_shaders:end}@*)

        (*@\Green{// Try to update all required transformation matrices (the method \_updateTransformationMatrices() has to}\label{src:GLWidgetTestingShaderPrograms:constructor:first_call_of_updateTransformationMatrices:start}@*)
        (*@\Green{// be called after every successful transformation related event handling).}@*)
        if (!_updateTransformationMatrices())
        {
            throw Exception("Could not update all transformation matrices!");
        } (*@\label{src:GLWidgetTestingShaderPrograms:constructor:first_call_of_updateTransformationMatrices:end}@*)

        (*@\Green{// creating a directional light object}\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:light:start}@*)

        Cartesian3 direction(0.0, 0.0, 1.0);
        Cartesian3 eye = _eye;
        eye.normalize();

        Cartesian3 half_vector = direction;
        half_vector += eye;

        half_vector.normalize();

        _light = SP<DirectionalLight>::Default(
                    new DirectionalLight(
                        Homogeneous3(0.0f, 0.0f, 1.0f, 0.0f), (*@\Green{// direction vector of the light}@*)
                        Homogeneous3(half_vector),            (*@\Green{// homogeneous half vector}@*)
                        Color4(0.4f, 0.4f, 0.4f, 1.0f),       (*@\Green{// ambient light intensity}@*)
                        Color4(0.8f, 0.8f, 0.8f, 1.0f),       (*@\Green{// diffuse light intensity}@*)
                        Color4(1.0f, 1.0f, 1.0f, 1.0f)));     (*@\Green{// specular light intensity}@*)

        if (!_light)
        {
            throw Exception("Could not create the directional light object!");
        } (*@\label{src:GLWidgetTestingShaderPrograms.cpp:constructor:light:end}@*)

        (*@\Green{// examples for communicating with our shader programs via uniform variables}@*)

        _two_sided_lighting.enable(); (*@\label{src:GLWidgetTestingShaderPrograms:initializing_shaders:start}@*)

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

        if (!_two_sided_lighting.setUniformMaterial("front_material", materials::brass))
        {
            throw Exception("Two-sided per pixed lighting: could not initialize the "
                            "uniform variable \"front_material\"!");
        }

        if (!_two_sided_lighting.setUniformMaterial("back_material", materials::chrome))
        {
            throw Exception("Two-sided per pixed lighting: could not initialize the "
                            "uniform variable \"back_material\"!");
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

        if (!_reflection_lines.setUniformColorMaterial("front_material", colors::green))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"front_material\"!");
        }

        if (!_reflection_lines.setUniformColorMaterial("back_material", colors::orange))
        {
            throw Exception("Reflection lines: could not initialize the "
                            "uniform variable \"back_material\"!");
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

        _reflection_lines.disable(); (*@\label{src:GLWidgetTestingShaderPrograms:initializing_shaders:end}@*)

        (*@\Red{// \textbf{Remark}}@*)
        (*@\Red{// Lines \mref{src:ShaderProgram:setUniformVariables:declaration:start}--\mref{src:ShaderProgram:setUniformVariables:declaration:end} of Listing \mref{src:ShaderPrograms.h} declare many other methods that can be used to initialize}@*)
        (*@\Red{// possible uniform variables. Here we used only a few of them.}@*)

        (*@\Red{// \textbf{Remark}}@*)
        (*@\Red{// If one does not want to use our class ShaderProgram, then lines \mref{src:GLWidgetTestingShaderPrograms:loading_shaders:start}--\mref{src:GLWidgetTestingShaderPrograms:loading_shaders:end}  and \mref{src:GLWidgetTestingShaderPrograms:initializing_shaders:start}--\mref{src:GLWidgetTestingShaderPrograms:initializing_shaders:end} of}@*)
        (*@\Red{// the current listing should be replaced by using other data types and programming techniques.}@*)

        glEnable(GL_DEPTH_TEST);              (*@\Green{// enable depth testing}@*)
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f); (*@\Green{// set the background color}@*)

        (*@\Green{// try to create your renderable geometry, by instantiating, e.g., our classes ECSpace, BCurve3, BSurface3,}@*)
        (*@\Green{// GenericCurve3 and TriangleMesh3\ldots}@*)
    }
    catch (Exception &e)
    {
        cout << e << endl;
    }
}

(*@\Green{// a private method that calculates the transformationmatrices \_M, \_VM, \_PVM and \_tN}@*)
GLboolean YourGLWidget::_updateTransformationMatrices()
{
    if (_P && _V) (*@\label{src:GLWidgetTestingShaderPrograms:_updateTransformationMatrices:implementation:start}@*)
    {
        _M              = _Rx * _Ry * _Rz * _T * _S;
        _VM             = (*_V) * _M;
        _PVM            = (*_P) * _VM;

        bool invertible = false;
        _tN             = _VM.inverse(&invertible);

        return invertible;
    }

    return GL_FALSE; (*@\label{src:GLWidgetTestingShaderPrograms:_updateTransformationMatrices:implementation:end}@*)
}

(*@\Green{// your rendering method}@*)
void YourGLWidget::render()
{
    (*@\Green{// clears the color and depth buffers}@*)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    (*@\Green{// For the sake of simplicity, in this method we have omitted the try--catch mechanism that would detect }@*)
    (*@\Green{// possible errors that originate from initializations of non-existent uniform variables.}@*)
    (*@\Green{// In what follows we outline three possibilities to use the provided shader programs.}@*)

    (*@\Green{// $1$st possibility}@*)
    if ((*@\Green{/* your curve points, derivatives, control polygons/nets that have to be rendered exist and are renderable*/}@*))
    {
        _color_shader.enable();

        _color_shader.setUniformColor("color", colors::black);
        _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

        (*@\Green{// Render your line geometry with the color shader program, e.g., you may pass the constant reference}@*)
        (*@\Green{// of the variable \_color\_shader as the first input variable to the rendering methods:}@*)
        (*@\Green{//}@*)
        (*@\Green{// GLboolean GenericCurve3::renderDerivatives(}@*)
        (*@\Green{// \hspace{0.5cm}const ShaderProgram \&program, GLint order, GLenum render\_mode) const;}@*)
        (*@\Green{//}@*)
        (*@\Green{// GLboolean BCurve3::renderData(}@*)
        (*@\Green{// \hspace{0.5cm}const ShaderProgram \&program, GLenum render\_mode = GL\_LINE\_STRIP) const;}@*)
        (*@\Green{//}@*)
        (*@\Green{// GLboolean BSurface3::renderData(}@*)
        (*@\Green{// \hspace{0.5cm}const ShaderProgram \&program,}@*)
        (*@\Green{// \hspace{0.5cm}GLenum u\_render\_mode = GL\_LINE\_STRIP, GLenum v\_render\_mode = GL\_LINE\_STRIP) const.}@*)

        (*@\Red{// \textbf{Remark}}@*)
        (*@\Red{// If you do not use our class ShaderProgram, then instead of the previously listed three rendering methods}@*)
        (*@\Red{// you should use the member functions:}@*)
        (*@\Red{//}@*)
        (*@\Red{// GLboolean GenericCurve3::renderDerivatives(}@*)
        (*@\Red{// \hspace{0.5cm}GLint order, GLenum render\_mode, GLint vec3\_position\_location = 0) const;}@*)
        (*@\Red{//}@*)
        (*@\Red{// GLboolean BCurve3::renderData(}@*)
        (*@\Red{// \hspace{0.5cm}GLenum render\_mode = GL\_LINE\_STRIP, GLint vec3\_position\_location = 0) const;}@*)
        (*@\Red{//}@*)
        (*@\Red{// GLboolean BSurface3::renderData(}@*)
        (*@\Red{// \hspace{0.5cm}GLenum u\_render\_mode = GL\_LINE\_STRIP, GLenum v\_render\_mode = GL\_LINE\_STRIP,}@*)
        (*@\Red{// \hspace{0.5cm}GLint vec3\_position\_location = 0) const,}@*)
        (*@\Red{//}@*)
        (*@\Red{// but in these cases you should provide a valid and active attribute location that is associated with positions}@*)
        (*@\Red{// of type vec3 in your custom shader programs. (At the same time it is your job to handle possible}@*)
        (*@\Red{// uniform variables.)}@*)

        _color_shader.disable();
    }

    (*@\Green{// $2$nd possibility}@*)
    if ((*@\Green{/* your triangles meshes that have to be rendered exist and are renderable*/}@*))
    {
        _two_sided_lighting.enable();

        _two_sided_lighting.setUniformMatrix4fv("VM",  1, GL_FALSE, _VM.address());
        _two_sided_lighting.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
        _two_sided_lighting.setUniformMatrix4fv("N",   1, GL_TRUE,  _tN.address());

        (*@\Green{// Render your triangles meshes with the two sided lighting shader program, e.g., you may pass the constant}@*)
        (*@\Green{// reference of the variable \_two\_sided\_lighting as the first input variable to the rendering method:}@*)
        (*@\Green{//}@*)
        (*@\Green{// GLboolean TriangleMesh3::render(}@*)
        (*@\Green{// \hspace{0.5cm}const ShaderProgram \&program, GLenum render\_mode = GL\_TRIANGLES) const.}@*)

        (*@\Red{// \textbf{Remark}}@*)
        (*@\Red{// If you do not use our class ShaderProgram, then instead of the previous rendering method you should}@*)
        (*@\Red{// use the member function:}@*)
        (*@\Red{//}@*)
        (*@\Red{// GLboolean TriangleMesh3::render(}@*)
        (*@\Red{// \hspace{0.5cm}GLenum render\_mode = GL\_TRIANGLES,}@*)
        (*@\Red{// \hspace{0.5cm}GLint vec3\_position\_location = 0, GLint vec3\_normal\_location = 1,}@*)
        (*@\Red{// \hspace{0.5cm}GLint vec4\_color\_location = 2, GLint vec4\_texture\_location = 3) const,}@*)
        (*@\Red{//}@*)
        (*@\Red{// but in this case you should provide valid and active attribute locations that are associated in your}@*)
        (*@\Red{// custom shader programs with positions, normals, colors and texture coordinates of type vec3, vec3,}@*)
        (*@\Red{// vec4 and vec4, respectively. (At the same time it is your job to handle possible uniform variables.)}@*)

        _two_sided_lighting.disable();
    }

    (*@\Green{// $3$rd possibility}@*)
    if ((*@\Green{/* your triangles meshes that have to be rendered exist and are renderable*/}@*))
    {
        _reflection_lines.enable();

        _reflection_lines.setUniformMatrix4fv("VM",  1, GL_FALSE, _VM.address());
        _reflection_lines.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());
        _reflection_lines.setUniformMatrix4fv("N",   1, GL_TRUE,  _tN.address());

        (*@\Green{// Render your triangles meshes with the reflection lines shader program, e.g., you may pass the constant}@*)
        (*@\Green{// reference of the variable \_reflection\_lines as the first input variable to the rendering method:}@*)
        (*@\Green{//}@*)
        (*@\Green{// GLboolean TriangleMesh3::render(}@*)
        (*@\Green{// \hspace{0.5cm}const ShaderProgram \&program, GLenum render\_mode = GL\_TRIANGLES) const.}@*)

        (*@\Red{// \textbf{Remark}}@*)
        (*@\Red{// If you do not use our class ShaderProgram, then instead of the previous rendering method you should}@*)
        (*@\Red{// use the member function:}@*)
        (*@\Red{//}@*)
        (*@\Red{// GLboolean TriangleMesh3::render(}@*)
        (*@\Red{// \hspace{0.5cm}GLenum render\_mode = GL\_TRIANGLES,}@*)
        (*@\Red{// \hspace{0.5cm}GLint vec3\_position\_location = 0, GLint vec3\_normal\_location = 1,}@*)
        (*@\Red{// \hspace{0.5cm}GLint vec4\_color\_location = 2, GLint vec4\_texture\_location = 3) const,}@*)
        (*@\Red{//}@*)
        (*@\Red{// but in this case you should provide valid and active attribute locations that are associated in your}@*)
        (*@\Red{// custom shader programs with positions, normals, colors and texture coordinates of type vec3, vec3,}@*)
        (*@\Red{// vec4 and vec4, respectively. (At the same time it is your job to handle possible uniform variables.)}@*)

        _reflection_lines.disable();
    }

    (*@\Red{// \textbf{Remark}}@*)
    (*@\Red{// Note that, the constant memory addresses of all transformation matrices have to be passed in this rendering}@*)
    (*@\Red{// method of your application to the corresponding uniform variables. The reason of this is the fact that the}@*)
    (*@\Red{// underlying matrices may change in your transformation related event handling methods. (Do not forget to}@*)
    (*@\Red{// call the method \_updateTransformationMatrices() before leaving your event handling slots.)}@*)
}
