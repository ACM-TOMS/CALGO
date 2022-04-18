#ifndef YOURGLWIDGET_H
#define YOURGLWIDGET_H

#include <GL/glew.h> (*@\Green{// in order to handle vertex buffer and shader program objects, we rely on GLEW}@*)
                     (*@\Green{// we assume that this external dependency is added to your project}@*)

#include <Core/Geometry/Curves/GenericCurves3.h>
#include <Core/Math/SpecialGLTransformations.h>
#include <Core/Shaders/ShaderPrograms.h>
#include <Core/SmartPointers/SpecializedSmartPointers.h>
#include <EC/ECSpaces.h>
#include <EC/BCurves3.h>

namespace cagd (*@\Green{// all data structures in the proposed function library are under the namespace cagd}@*)
{
    class YourGLWidget: public ...
    {
    private:
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.h:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.h:color_shader:declaration} of Listing \mref{src:GLWidgetTestingShaderPrograms.h}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// EC space parameters:}@*)
        GLdouble                                _alpha, _beta; (*@\Green{// endpoints;}@*)
        GLdouble                                _omega;        (*@\Green{// shape parameter;}@*)
        SP<ECSpace>::Default                    _space;        (*@\Green{// smart pointer to the EC space.}@*)

        (*@\Green{// number of uniform subdivision points in the definition domains, it is used for image generation}@*)
        GLint                                   _maximum_order_of_derivatives;
        (*@\Green{// determines the maximum order of derivatives that have to be evaluated along the B-curves}@*)
        GLint                                   _div_point_count;

        (*@\Green{// number of arcs along the logarithmic spireal (\mref{eq:logarithmic_spiral})}@*)
        GLint                                   _arc_count;
        (*@\Green{// a vector of smart pointers to dynamically allocated B-curves}@*)
        std::vector<SP<BCurve3>::Default>       _bcurve;
        (*@\Green{// a vector of smart pointers to dynamically allocated B-curve images (i.e., GenericCurve3 objects)}@*)
        std::vector<SP<GenericCurve3>::Default> _img_bcurve;

    public:
        (*@\Green{// your default/special constructor}@*)
        YourGLWidget(/*...*/);

        (*@\Green{// your rendering method}@*)
        void render();

        (*@\Green{// your other methods...}@*)
    };
}

#endif (*@\Green{// YOURGLWIDGET\_H}@*)
