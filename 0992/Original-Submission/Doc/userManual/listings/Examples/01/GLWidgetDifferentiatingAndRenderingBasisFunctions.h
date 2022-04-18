#ifndef YOURGLWIDGET_H
#define YOURGLWIDGET_H

#include <GL/glew.h> (*@\Green{// in order to handle vertex buffer and shader program objects, we rely on GLEW}@*)
                     (*@\Green{// we assume that this external dependency is added to your project}@*)

#include <Core/Geometry/Curves/GenericCurves3.h>
#include <Core/Math/SpecialGLTransformations.h>
#include <Core/Shaders/ShaderPrograms.h>
#include <Core/SmartPointers/SpecializedSmartPointers.h>
#include <EC/ECSpaces.h>

namespace cagd (*@\Green{// all data structures in the proposed function library are defined under the namespace cagd}@*)
{
    class YourGLWidget: public ...
    {
    private:
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.h:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.h:color_shader:declaration} of Listing \mref{src:GLWidgetTestingShaderPrograms.h}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// Parameters of different EC spaces:}@*)
        GLint                             _space_count;  (*@\Green{// number of considered EC spaces;}@*)
        std::vector<GLdouble>             _alpha, _beta; (*@\Green{// endpoints of definition domains;}@*)
        std::vector<GLint>                _n;            (*@\Green{// orders of considered EC spaces;}@*)
        std::vector<SP<ECSpace>::Default> _space;        (*@\Green{// smart pointers to different EC spaces.}@*)

        (*@\Green{// Parameters used for the image generation of normalized B-basis functions:}@*)
        (*@\Green{// the common number of subdivision points in each definition domain;}@*)
        GLint _div_point_count;
        (*@\Green{// the maximum order of derivatives that have to be evaluated along the basis functions;}@*)
        GLint _maximum_order_of_derivatives;
        (*@\Green{// an array of smart pointers pointing to row matrices that store the image smart pointers of all}@*)
        (*@\Green{// normalized B-basis functions of each EC space.}@*)
        std::vector< SP< RowMatrix<SP<GenericCurve3>::Default> >::Default > _img_B_basis;

        (*@\Green{// your other declarations\ldots}@*)

    public:
        (*@\Green{// your default/special constructor}@*)
        YourGLWidget(/*...*/);

        (*@\Green{// your rendering method}@*)
        void render();

        (*@\Green{// your other methods...}@*)
    };
}

#endif (*@\Green{// YOURGLWIDGET\_H}@*)
