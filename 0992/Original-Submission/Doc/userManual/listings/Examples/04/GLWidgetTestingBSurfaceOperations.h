#ifndef YOURGLWIDGET_H
#define YOURGLWIDGET_H

#include <GL/glew.h>

#include <Core/SmartPointers/SpecializedSmartPointers.h>
#include <Core/Geometry/Surfaces/Lights.h>
#include <Core/Geometry/Surfaces/TriangleMeshes3.h>
#include <Core/Math/SpecialGLTransformations.h>
#include <Core/Shaders/ShaderPrograms.h>
#include <EC/BSurfaces3.h>
#include <EC/ECSpaces.h>

namespace cagd (*@\Green{// all data structures in the proposed function library are defined under the namespace cagd}@*)
{
    class YourGLWidget: public ...
    {
    private:
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.h:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.h:directional_light:declaration} of Listing \mref{src:GLWidgetTestingShaderPrograms.h}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// A triangle mesh that will store a triangulated unit sphere centered at the origin. Using translation and}@*)
        (*@\Green{// scaling transformations, it will be rendered multiple times at the positions of the control points.}@*)
        TriangleMesh3                     _sphere;

        (*@\Green{// Determines the common scaling factor of the unit sphere that has to be rendered at the positions}@*)
        (*@\Green{// of the control points.}@*)
        GLdouble                          _control_point_radius;

        (*@\Green{// Determines the scaling transformation of the unit sphere that has to be rendered at the positions}@*)
        (*@\Green{// of the control points.}@*)
        Scale                             _sphere_S;

        (*@\Green{// Parameters of the EC spaces associated with directions $u$ and $v$:}@*)
        std::vector<GLdouble>             _alpha, _beta; (*@\Green{// endpoints of definition domains;}@*)
        std::vector<GLint>                _n;            (*@\Green{// orders of the corresponding EC spaces;}@*)
        std::vector<SP<ECSpace>::Default> _space;        (*@\Green{// an array of smart pointers to EC spaces.}@*)

        (*@\Green{// an array of numbers of uniform subdivision points in the corresponding definition domains,}@*)
        (*@\Green{// they will be used for image generation}@*)
        std::vector<GLint>                _div_point_count;

        (*@\Green{// Determines the color scheme of the images (i.e., triangle meshes) of the B-surfaces that have to be}@*)
        (*@\Green{// rendered. The user can choose from $13$ possible color schemes such as DEFAULT\_NULL\_FRAGMENT,}@*)
        (*@\Green{// X\_VARIATION\_FRAGMENT, Y\_VARIATION\_FRAGMENT, Z\_VARIATION\_FRAGMENT,}@*)
        (*@\Green{// NORMAL\_LENGTH\_FRAGMENT, GAUSSIAN\_CURVATURE\_FRAGMENT,}@*)
        (*@\Green{// MEAN\_CURVATURE\_FRAGMENT, WILLMORE\_ENERGY\_FRAGMENT,}@*)
        (*@\Green{// LOG\_WILLMORE\_ENERGY\_FRAGMENT, UMBILIC\_DEVIATION\_ENERGY\_FRAGMENT,}@*)
        (*@\Green{// LOG\_UMBILIC\_DEVIATION\_ENERGY\_FRAGMENT, TOTAL\_CURVATURE\_ENERGY\_FRAGMENT,}@*)
        (*@\Green{// LOG\_TOTAL\_CURVATURE\_ENERGY\_FRAGMENT. For more details see the lines \mref{src:TensorProductSurface3.h:ImageColorScheme:start}--\mref{src:TensorProductSurface3.h:ImageColorScheme:end} of}@*)
        (*@\Green{// Listing \mref{src:TensorProductSurfaces3.h} and Fig.\ \mref{fig:color_schemes}. \Red{(With the exception of DEFAULT\_NULL\_FRAGMENT all remaining}}@*)
        (*@\Green{// \Red{color schemes should be used together with uniform black (color) materials.)}}@*)
        TensorProductSurface3::ImageColorScheme                       _color_scheme;

        (*@\Green{// smart pointer to a dynamically allocated randomly generated B-surface}@*)
        SP<BSurface3>::Default                                        _bsurface;

        (*@\Green{// smart pointer to a dynamically allocated B-surface image (i.e., TriangleMesh3 object)}@*)
        SP<TriangleMesh3>::Default                                    _img_bsurface;

        (*@\Green{// smart pointer to dynamically allocated order elevated B-surface}@*)
        SP<BSurface3>::Default                                        _oe_bsurface;

        (*@\Green{// smart pointer to the image of the dynamically allocated order elevated B-surface}@*)
        SP<TriangleMesh3>::Default                                    _oe_img_bsurface;

        (*@\Green{// an array of percentages that determine the subdivision points of the $u$- and $v$-directional definition}@*)
        (*@\Green{// domains, where the initial B-surfaces have to be subdivided}@*)
        std::vector<GLdouble>                                         _ratio;

        (*@\Green{// an array of smart pointers to row matrices that store two smart pointers to the $u$- and $v$-directional}@*)
        (*@\Green{// subdivided B-surfaces}@*)
        std::vector<SP< RowMatrix<SP<BSurface3>::Default> >::Default> _subdivision;

        (*@\Green{// an array of row matrices that store two smart pointers to the images of the $u$- and $v$-directional}@*)
        (*@\Green{// subdivided B-surfaces}@*)
        std::vector< RowMatrix<SP<TriangleMesh3>::Default> >          _img_subdivision;

        (*@\Green{// decides whether the two-sided lighting or the reflection lines shader program should be used}@*)
        (*@\Green{// during rendering}@*)
        bool _apply_reflection_lines;                    (*@\Green{// synchronize it with a check box}@*)

        (*@\Green{// visibility flags that should be synchronized with radio buttons and check boxes of your graphical}@*)
        (*@\Green{// user interface}@*)
        bool _show_randomly_generated_initial_B_surface; (*@\Green{// synchronize it with a radio button}@*)
        bool _show_order_elevated_B_surface;             (*@\Green{// synchronize it with a radio button}@*)
        bool _show_u_subdivided_B_surfaces;              (*@\Green{// synchronize it with a radio button}@*)
        bool _show_v_subdivided_B_surfaces;              (*@\Green{// synchronize it with a radio button}@*)
        bool _compare_control_nets;                      (*@\Green{// synchronize it with a check box}@*)

        (*@\Green{// auxiliary private methods used for control point, control net and transparent triangle mesh rendering}@*)
        void _renderSpheresAtControlPoints(
            const BSurface3 &surface, const Color4 &front_color_material) const;

        void _renderControlNet(
            const BSurface3 &surface, const Color4 &color,
            bool use_dashed_line = false) const;

        void _renderTransparentMesh(
            const ShaderProgram &shader, const TriangleMesh3 &mesh,
            const Color4 &front_color_material, const Color4 &back_color_material,
            GLfloat transparency = 0.5f) const;

    public:
        (*@\Green{// your default/special constructor}@*)
        YourGLWidget(/*...*/);

        (*@\Green{// your rendering method}@*)
        void render();

        (*@\Green{// your other methods...}@*)
    };
}

#endif (*@\Green{// YOURGLWIDGET\_H}@*)
