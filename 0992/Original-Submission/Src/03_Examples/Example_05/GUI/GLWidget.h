#ifndef GLWIDGET_H
#define GLWIDGET_H

#if defined(_MSC_VER)
#pragma warning(disable:4503)
#endif

#include <GL/glew.h>

#include <QtOpenGL/QGLWidget>
#include <QtOpenGL/QGLFormat>

#include <Core/SmartPointers/SpecializedSmartPointers.h>
#include <Core/Geometry/Surfaces/Lights.h>
#include <Core/Geometry/Surfaces/TriangleMeshes3.h>
#include <Core/Math/SpecialGLTransformations.h>
#include <Core/Shaders/ShaderPrograms.h>
#include <EC/BSurfaces3.h>
#include <EC/ECSpaces.h>

namespace cagd
{
    class GLWidget: public QGLWidget
    {
        Q_OBJECT

    private:
        // variables defining the orthogonal projection matrix
        GLfloat                             _aspect;         // aspect ratio of the rendering window
        GLfloat                             _left,   _right; // minimum and maximum values of $x$-coordinates
        GLfloat                             _bottom, _top;   // minimum and maximum values of $y$-coordinates
        GLfloat                             _near,   _far;   // minimum and maximum values of $z$-coordinates
        SP<OrthogonalProjection>::Default   _P; // smart pointer to an orthogonal projection matrix

        // variables defining the view matrix
        Cartesian3                          _eye, _center, _up;
        SP<LookAt>::Default                 _V; // smart pointer to the view or world matrix

        // transformation matrices
        Rotate                              _Rx, _Ry, _Rz; // rotation matrices around axis $x$, $y$ and $z$, respectively
        Translate                           _T;            // translation
        Scale                               _S;            // scaling
        GLTransformation                    _M,            // model matrix, i.e., the product _Rx * _Ry * _Rz * _T * _S
                                            _VM,           // product of view and model matrices
                                            _PVM,          // product of projection, view and model matrices
                                            _tN;           // transposed normal matrix, i.e., inverse of _VM

        // a private method that calculates the transformationmatrices _M, _VM, _PVM and _tN
        GLboolean                           _updateTransformationMatrices();

        // shader program objects
        ShaderProgram                       _color_shader;
        ShaderProgram                       _two_sided_lighting;
        ShaderProgram                       _reflection_lines;

        // a smart pointer to a directional light object
        SP<DirectionalLight>::Default       _light;

        // A triangle mesh that will store a triangulated unit sphere centered at the origin.
        // Using translation and scaling transformation, this sphere will be rendered multiple
        // times at the positions of the control points.
        TriangleMesh3                      _sphere;

        // Determines the common scaling factor of the unit sphere that has to be rendered at
        // the positions of the control points.
        GLdouble                           _control_point_radius;

        // Determines the scaling transformation of the unit sphere that has to be rendered at
        // the positions of the control points.
        Scale                              _sphere_S;

        // A triangle mesh that will store a triangulated right circular cone with unit base
        // radius and appex $(0, 0, 2)$.
        // Using translation and scaling transformation, this cone will be rendered multiple
        // times at the endpoints of the tangent vectors of the isoparametric lines.
        TriangleMesh3                      _cone;

        // Determines the scaling transformation of the cone that has to be rendered at the
        // endpoints of the tangent vectors of the isoparametric lines.
        Scale                              _cone_S;

        // a 2-element array that stores the dimensions of the EC spaces applied in direction $u$ and $v$
        std::vector<GLint>                 _dimension;

        // a 2-element array of numbers of uniform subdivision points in the corresponding definition domains,
        // they will be used for surface image (i.e., TriangleMesh3) generation
        std::vector<GLint>                 _surface_div_point_count;

        // a 2-element array of isoparametric line counts in the corresponding directions
        std::vector<GLint>                 _isoparametric_line_count;

        // a 2-element array of maximum differentiation orders in the corresponding directions
        std::vector<GLint>                 _maximum_order_of_derivatives;

        // a 2-element array of numbers of uniform subdivision points in the corresponding definition domains,
        // they will be used for curve image (i.e., GenericCurve3) generation
        std::vector<GLint>                 _curve_div_point_count;

        // Determines the color scheme of the images (i.e., triangle meshes) of the B-surfaces that have to be
        // rendered. The user can choose from $13$ possible color schemes such as
        // DEFAULT_NULL_FRAGMENT,
        // X_VARIATION_FRAGMENT, Y_VARIATION_FRAGMENT, Z_VARIATION_FRAGMENT,
        // NORMAL_LENGTH_FRAGMENT,
        // GAUSSIAN_CURVATURE_FRAGMENT, MEAN_CURVATURE_FRAGMENT,
        // WILLMORE_ENERGY_FRAGMENT, LOG_WILLMORE_ENERGY_FRAGMENT,
        // UMBILIC_DEVIATION_ENERGY_FRAGMENT, LOG_UMBILIC_DEVIATION_ENERGY_FRAGMENT,
        // TOTAL_CURVATURE_ENERGY_FRAGMENT, LOG_TOTAL_CURVATURE_ENERGY_FRAGMENT.
        TensorProductSurface3::ImageColorScheme _color_scheme;

        // a matrix of smart pointers to dynamically allocated B-surfaces
        Matrix<SP<BSurface3>::Default>     _patches;

        // a matrix of smart pointers to dynamically allocated B-surface images (i.e., TriangleMesh3 object)
        Matrix<SP<TriangleMesh3>::Default> _img_patches;

        // a matrix of smart pointers to dynamically allocated row matrices that store smart pointers
        // to $u$-directional isoparametric lines
        Matrix<SP< RowMatrix<SP<GenericCurve3>::Default> >::Default> _u_isoparametric_lines;

        // a matrix of smart pointers to dynamically allocated row matrices that store smart pointers
        // to $v$-directional isoparametric lines
        Matrix<SP< RowMatrix<SP<GenericCurve3>::Default> >::Default> _v_isoparametric_lines;

        // decides whether the two-sided lighting or the reflection lines shader program should be used during rendering
        bool                               _apply_reflection_lines;

        // visibility flags
        bool                               _show_patches;
        bool                               _show_control_nets;
        bool                               _show_u_isoparametric_lines;
        bool                               _show_v_isoparametric_lines;
        bool                               _show_tangents_of_u_isoparametric_lines;
        bool                               _show_tangents_of_v_isoparametric_lines;
        bool                               _transparency;

        // auxiliar private rendering methods
        void _renderSpheresAtControlPoints(
                const BSurface3 &surface,
                const Color4 &front_color_material,
                const Color4 &back_color_material) const;

        void _renderControlNet(
                const BSurface3 &surface,
                const Color4 &color,
                bool use_dashed_line = false) const;

        void _renderTransparentMesh(
                const ShaderProgram &shader,
                const TriangleMesh3 &mesh,
                const Color4 &front_color_material,
                const Color4 &back_color_material,
                GLfloat transparency = 0.5f) const;

        bool _renderSpheresAndConesAtEndpointsOfTangentVectors(
                const GenericCurve3 &curve,
                const Color4 &front_color_material,
                const Color4 &back_color_material) const;

    public:
        // special and default constructor
        // the format specifies the properties of the rendering window
        GLWidget(QWidget* parent = 0, const QGLFormat& format = QGL::Rgba | QGL::DepthBuffer | QGL::DoubleBuffer);

    protected:
        // redeclared virtual functions
        void initializeGL();
        void paintGL();
        void resizeGL(int w, int h);

    public slots:
        // public event handling methods/slots
        void setAngleX(int value);
        void setAngleY(int value);
        void setAngleZ(int value);

        void setZoomFactor(double value);

        void setTransX(double value);
        void setTransY(double value);
        void setTransZ(double value);

        void toggleReflectionLines(bool enabled);
        void setVisibilityOfPatches(bool value);
        void setVisibilityOfControlNets(bool value);
        void setVisibilityOfUIsoparametricLines(bool value);
        void setVisibilityOfVIsoparametricLines(bool value);
        void setVisibilityOfTangentsOfUIsoparametricLines(bool value);
        void setVisibilityOfTangentsOfVIsoparametricLines(bool value);
        void setTransparency(bool value);
    };
}

#endif // GLWIDGET_H
