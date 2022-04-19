#ifndef YOURGLWIDGET_H
#define YOURGLWIDGET_H

#include <GL/glew.h> (*@\Green{// in order to handle vertex buffer and shader program objects, we rely on GLEW}@*)
                     (*@\Green{// we assume that this external dependency is added to your project}@*)

#include <Core/Geometry/Coordinates/Cartesians3.h>
#include <Core/Geometry/Surfaces/Lights.h>
#include <Core/Math/SpecialGLTransformations.h>
#include <Core/Shaders/ShaderPrograms.h>
#include <Core/SmartPointers/SpecializedSmartPointers.h>

namespace cagd (*@\Green{// all data structures in the proposed function library are defined under the namespace cagd}@*)
{
    class YourGLWidget: public ...
    {
    private:
        (*@\Green{// variables defining the orthogonal projection matrix}\label{src:GLWidgetTestingShaderPrograms.h:part_1:start}@*)
        GLfloat _aspect;         (*@\Green{// aspect ratio of the rendering window}@*)
        GLfloat _left,   _right; (*@\Green{// minimum and maximum values of $x$-coordinates}@*)
        GLfloat _bottom, _top;   (*@\Green{// minimum and maximum values of $y$-coordinates}@*)
        GLfloat _near,   _far;   (*@\Green{// minimum and maximum values of $z$-coordinates}@*)
        SP<OrthogonalProjection>::Default _P; (*@\Green{// smart pointer to an orthogonal projection matrix}@*)

        (*@\Green{// variables defining the view matrix}@*)
        Cartesian3          _eye, _center, _up;
        SP<LookAt>::Default _V; (*@\Green{// smart pointer to the view or world matrix}@*)

        (*@\Green{// transformation matrices}@*)
        Rotate           _Rx, _Ry, _Rz; (*@\Green{// rotation matrices around axis $x$, $y$ and $z$, respectively}@*)
        Translate        _T;            (*@\Green{// translation}@*)
        Scale            _S;            (*@\Green{// scaling}@*)
        GLTransformation _M,            (*@\Green{// model matrix, i.e., the product \_Rx * \_Ry * \_Rz * \_T * \_S}@*)
                         _VM,           (*@\Green{// product of view and model matrices}@*)
                         _PVM,          (*@\Green{// product of projection, view and model matrices}@*)
                         _tN;           (*@\Green{// transposed normal matrix, i.e., inverse of \_VM}@*)

        (*@\Green{// a private method that calculates the transformationmatrices \_M, \_VM, \_PVM and \_tN}@*)
        GLboolean        _updateTransformationMatrices(); (*@\label{src:GLWidgetTestingShaderPrograms.h:part_1:end}@*)

        (*@\Green{// shader program objects}@*)
        ShaderProgram    _color_shader;       (*@\label{src:GLWidgetTestingShaderPrograms.h:color_shader:declaration}@*)
        ShaderProgram    _two_sided_lighting; (*@\label{src:GLWidgetTestingShaderPrograms.h:two_sided_lighting:declaration}@*)
        ShaderProgram    _reflection_lines;   (*@\label{src:GLWidgetTestingShaderPrograms.h:reflection_lines:declaration}@*)

        (*@\Green{// a smart pointer to a directional light object}@*)
        SP<DirectionalLight>::Default _light; (*@\label{src:GLWidgetTestingShaderPrograms.h:directional_light:declaration}@*)

        (*@\Green{// your other declarations, e.g.\ one may declare smart pointers to objects of type ECSpace, BCurve3,}@*)
        (*@\Green{// BSurface3, GenericCurve3 and TriangleMesh3...}@*)

    public:
        (*@\Green{// your default/special constructor}@*)
        YourGLWidget(/*...*/);

        (*@\Green{// your rendering method}@*)
        void render();

        (*@\Green{// your other methods...}@*)
    };
}

#endif (*@\Green{// YOURGLWIDGET\_H}@*)
