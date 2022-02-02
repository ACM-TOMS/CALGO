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
        GLfloat _aspect;         // aspect ratio of the rendering window
        GLfloat _left,   _right; // minimum and maximum values of $x$-coordinates
        GLfloat _bottom, _top;   // minimum and maximum values of $y$-coordinates
        GLfloat _near,   _far;   // minimum and maximum values of $z$-coordinates
        SP<OrthogonalProjection>::Default _P; // smart pointer to an orthogonal projection matrix

        // variables defining the view matrix
        Cartesian3          _eye, _center, _up;
        SP<LookAt>::Default _V; // smart pointer to the view or world matrix

        // transformation matrices
        Rotate           _Rx, _Ry, _Rz; // rotation matrices around axis $x$, $y$ and $z$, respectively
        Translate        _T;            // translation
        Scale            _S;            // scaling
        GLTransformation _M,            // model matrix, i.e., the product _Rx * _Ry * _Rz * _T * _S
                         _VM,           // product of view and model matrices
                         _PVM,          // product of projection, view and model matrices
                         _tN;           // transposed normal matrix, i.e., inverse of _VM

        // a private method that calculates the transformationmatrices _M, _VM, _PVM and _tN
        GLboolean        _updateTransformationMatrices();

        // shader program objects
        ShaderProgram    _color_shader;
        ShaderProgram    _two_sided_lighting;
        ShaderProgram    _reflection_lines;

        // a smart pointer to a directional light object
        SP<DirectionalLight>::Default _light;

        // A triangle mesh that will store a triangulated unit sphere centered at the origin.
        // Using translation and scaling transformation, this sphere will be rendered multiple times at the positions of the control points.
        TriangleMesh3                                   _sphere;

        // Determines the common scaling factor of the unit sphere that has to be rendered at the positions of the control points.
        GLdouble                                        _control_point_radius;

        // Determines the scaling transformation of the unit sphere that has to be rendered at the positions of the control points.
        Scale                                            _sphere_S;

        // Parameters of the EC spaces associated with directions $u$ and $v$:
        std::vector<GLdouble>                            _alpha, _beta;   // endpoints of definition domains
        std::vector<GLint>                               _n;              // orders of the corresponding EC spaces
        std::vector<SP<ECSpace>::Default>                _space;          // an array of smart pointers to EC spaces

        // an array of numbers of uniform subdivision points in the corresponding definition domains, they will be used for image generation
        std::vector<GLint>                               _div_point_count;


        // smart pointer to a dynamically allocated randomly generated B-surface
        SP<BSurface3>::Default                           _bsurface;

        TensorProductSurface3::ImageColorScheme          _color_scheme;

        // smart pointer to a dynamically allocated B-surface image (i.e., TriangleMesh3 object)
        SP<TriangleMesh3>::Default                       _img_bsurface;

        // smart pointer to dynamically allocated order elevated B-surface
        SP<BSurface3>::Default                           _oe_bsurface;

        // smart pointer to dynamically allocated order elevated B-surface images
        SP<TriangleMesh3>::Default                       _oe_img_bsurface;

        // an array of percentages that determine the subdivision points of the $u$- and $v$-directional definition domains, where the initial B-surfaces have to be subdivided
        std::vector<GLdouble>                            _ratio;

        // an array of smart pointers to row matrices that store two smart pointers to the $u$- and $v$-directional subdivided B-surfaces
        std::vector<SP< RowMatrix<SP<BSurface3>::Default> >::Default> _subdivision;

        // an array of row matrices that store two smart pointers to the images of the $u$- and $v$-directional subdivided B-surfaces
        std::vector< RowMatrix<SP<TriangleMesh3>::Default> >           _img_subdivision;

        // decides whether the two-sided lighting or the reflection lines shader program should be used during rendering
        bool                                             _apply_reflection_lines;

        // visibility flags
        bool                                             _show_randomly_generated_initial_B_surface;
        bool                                             _show_order_elevated_B_surface;
        bool                                             _show_u_subdivided_B_surfaces;
        bool                                             _show_v_subdivided_B_surfaces;
        bool                                             _compare_control_nets;

        // auxiliar private rendering methods
        void _renderSpheresAtControlPoints(const BSurface3 &surface, const Color4 &front_color_material, const Color4 &back_color_material) const;
        void _renderControlNet(const BSurface3 &surface, const Color4 &color, bool use_dashed_line = false) const;
        void _renderTransparentMesh(const ShaderProgram &shader, const TriangleMesh3 &mesh, const Color4 &front_color_material, const Color4 &back_color_material, GLfloat transparency = 0.5f) const;

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

        void toggleReflectionLines(bool value);
        void setVisibilityOfInitialBSurface(bool value);
        void setVisibilityOfOrderElevatedBSurface(bool value);
        void setVisibilityOfSubdivisionsInDirectionU(bool value);
        void setVisibilityOfSubdivisionsInDirectionV(bool value);
        void setVisibilityOfInitialControlNet(bool value);
    };
}

#endif // GLWIDGET_H
