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

        (*@\Green{// A triangle mesh that will store a triangulated right circular cone with unit base radius and appex $(0, 0, 2)$.}@*)
        (*@\Green{// Using translation and scaling transformation, this cone will be rendered multiple times at the endpoints}@*)
        (*@\Green{// of the tangent vectors of the isoparametric lines.}@*)
        TriangleMesh3                     _cone;

        (*@\Green{// Determines the scaling transformation of the cone that has to be rendered at the endpoints of the}@*)
        (*@\Green{// tangent vectors of the isoparametric lines.}@*)
        Scale                             _cone_S;

        (*@\Green{// Parameters of the EC spaces associated with directions $u$ and $v$:}@*)
        std::vector<GLdouble>             _alpha, _beta; (*@\Green{// endpoints of definition domains;}@*)
        std::vector<GLint>                _n;            (*@\Green{// orders of the corresponding EC spaces;}@*)
        std::vector<SP<ECSpace>::Default> _space;        (*@\Green{// an array of smart pointers to EC spaces.}@*)

        (*@\Green{// a 2-element array that stores the dimensions of the EC spaces applied in direction $u$ and $v$}@*)
        std::vector<GLint>                _dimension;

        (*@\Green{// a 2-element array of numbers of uniform subdivision points in the corresponding definition domains,}@*)
        (*@\Green{// they will be used for surface image (i.e., TriangleMesh3) generation}@*)
        std::vector<GLint>                _surface_div_point_count;

        (*@\Green{// a 2-element array of isoparametric line counts in the corresponding directions}@*)
        std::vector<GLint>                _isoparametric_line_count;

        (*@\Green{// a 2-element array of maximum differentiation orders in the corresponding directions}@*)
        std::vector<GLint>                _maximum_order_of_derivatives;

        (*@\Green{// a 2-element array of numbers of uniform subdivision points in the corresponding definition domains,}@*)
        (*@\Green{// they will be used for curve image (i.e., GenericCurve3) generation}@*)
        std::vector<GLint>                _curve_div_point_count;

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
        TensorProductSurface3::ImageColorScheme _color_scheme;

        (*@\Green{// a matrix of smart pointers to dynamically allocated B-surface patches}@*)
        Matrix<SP<BSurface3>::Default>          _patches;

        (*@\Green{// a matrix of smart pointers to dynamically allocated B-surface images (i.e., TriangleMesh3 object)}@*)
        Matrix<SP<TriangleMesh3>::Default>      _img_patches;

        (*@\Green{// a matrix of smart pointers to dynamically allocated row matrices that store smart pointers to the }@*)
        (*@\Green{// dynamically allocated images (i.e., GenericCurve3 objects) of the $u$-directional isoparametric lines}@*)
        (*@\Green{// of all B-surface patches}@*)
        Matrix<SP< RowMatrix<SP<GenericCurve3>::Default> >::Default>
                _u_isoparametric_lines;

        (*@\Green{// a matrix of smart pointers to dynamically allocated row matrices that store smart pointers to the }@*)
        (*@\Green{// dynamically allocated images (i.e., GenericCurve3 objects) of the $v$-directional isoparametric lines}@*)
        (*@\Green{// of all B-surface patches}@*)
        Matrix<SP< RowMatrix<SP<GenericCurve3>::Default> >::Default>
                _v_isoparametric_lines;

        (*@\Green{// decides whether the two-sided lighting or the reflection lines shader program should be used}@*)
        (*@\Green{// during rendering}@*)
        bool    _apply_reflection_lines;                 (*@\Green{// synchronize it with a check box}@*)

        (*@\Green{// visibility flags that should be synchronized with check boxes of your graphical user interface}@*)
        bool    _show_patches;                           (*@\Green{// synchronize it with a check box}@*)
        bool    _show_control_nets;                      (*@\Green{// synchronize it with a check box}@*)
        bool    _show_u_isoparametric_lines;             (*@\Green{// synchronize it with a check box}@*)
        bool    _show_v_isoparametric_lines;             (*@\Green{// synchronize it with a check box}@*)
        bool    _show_tangents_of_u_isoparametric_lines; (*@\Green{// synchronize it with a check box}@*)
        bool    _show_tangents_of_v_isoparametric_lines; (*@\Green{// synchronize it with a check box}@*)

        (*@\Green{// determines the transparency of the rendered triangle meshes, its values should be in the interval $[0, 1]$}@*)
        GLfloat     _transparency;

        (*@\Green{// auxiliary private methods used for control point, control net, triangle mesh and tangent vector rendering}@*)
        void _renderSpheresAtControlPoints(
            const BSurface3 &surface, const Color4 &front_color_material) const;

        void _renderControlNet(
            const BSurface3 &surface, const Color4 &color,
            bool use_dashed_line = false) const;

        void _renderTransparentMesh(
            const ShaderProgram &shader, const TriangleMesh3 &mesh,
            const Color4 &front_color_material, const Color4 &back_color_material,
            GLfloat transparency = 0.5f) const;

        bool _renderSpheresAndConesAtEndpointsOfTangentVectors(
            const GenericCurve3 &curve,
            const Color4 &front_color_material,
            const Color4 &back_color_material) const;

    public:
        (*@\Green{// your default/special constructor}@*)
        YourGLWidget(/*...*/);

        (*@\Green{// your rendering method}@*)
        void render();

        (*@\Green{// your other methods...}@*)
    };
}

#endif (*@\Green{// YOURGLWIDGET\_H}@*)
