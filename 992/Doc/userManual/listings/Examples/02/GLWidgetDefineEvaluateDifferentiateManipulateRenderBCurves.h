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
        std::vector<SP<ECSpace>::Default> _space;        (*@\Green{// smart pointers to different EC spaces;}@*)
        std::vector<Color4>               _color;        (*@\Green{// colors associated with different EC spaces.}@*)

        (*@\Green{// number of uniform subdivision points in the definition domains, it is used for image generation}@*)
        GLint                                   _div_point_count;
        (*@\Green{// determines the maximum order of derivatives that have to be evaluated along the B-curves}@*)
        GLint                                   _maximum_order_of_derivatives;

        (*@\Green{// a vector of smart pointers to dynamically allocated B-curves}@*)
        std::vector<SP<BCurve3>::Default>       _bcurve;
        (*@\Green{// a vector of smart pointers to dynamically allocated B-curve images (i.e., GenericCurve3 objects)}@*)
        std::vector<SP<GenericCurve3>::Default> _img_bcurve;

        (*@\Green{// a vector of smart pointers to dynamically allocated order elevated B-curves}@*)
        std::vector<SP<BCurve3>::Default>       _oe_bcurve;
        (*@\Green{// a vector of smart pointers to dynamically allocated order elevated B-curve images}@*)
        std::vector<SP<GenericCurve3>::Default> _img_oe_bcurve;

        (*@\Green{// determines the subdivision point of the definition domains, where the initial B-curves have to be subdivided}@*)
        GLdouble                                                      _ratio;
        (*@\Green{// a vector of smart pointers to row matrices that store two smart pointers to the subdivided arcs}@*)
        std::vector< SP< RowMatrix<SP<BCurve3>::Default> >::Default > _subdivision;
        (*@\Green{// a vector of row matrices that store two smart pointers to the images of the subdivided arcs}@*)
        std::vector< RowMatrix<SP<GenericCurve3>::Default> >          _img_subdivision;

    public:
        (*@\Green{// your default/special constructor}@*)
        YourGLWidget(/*...*/);

        (*@\Green{// your rendering method}@*)
        void render();

        (*@\Green{// your other methods...}@*)
    };
}

#endif (*@\Green{// YOURGLWIDGET\_H}@*)
