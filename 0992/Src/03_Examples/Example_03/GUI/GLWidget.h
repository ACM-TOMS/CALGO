#ifndef GLWIDGET_H
#define GLWIDGET_H

#if defined(_MSC_VER)
#pragma warning(disable:4503)
#endif

#include <GL/glew.h>

#include <QtOpenGL/QGLWidget>
#include <QtOpenGL/QGLFormat>

#include <Core/Geometry/Curves/GenericCurves3.h>
#include <Core/Math/SpecialGLTransformations.h>
#include <Core/Shaders/ShaderPrograms.h>
#include <Core/SmartPointers/SpecializedSmartPointers.h>
#include <EC/ECSpaces.h>
#include <EC/BCurves3.h>

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
        GLboolean                               _updateTransformationMatrices();

        // color shader program object
        ShaderProgram                           _color_shader;

        // EC space parameters
        GLdouble                                _alpha, _beta;
        GLdouble                                _omega;
        SP<ECSpace>::Default                    _space;

        // number of uniform subdivision points in the definition domains, it is used for image generation
        GLint                                   _maximum_order_of_derivatives;
        // determines the maximum order of derivatives that have to be evaluated along the B-curves
        GLint                                   _div_point_count;

        // number of arcs along the logarithmic spireal (\ref{eq:logarithmic_spiral})
        GLint                                   _arc_count;
        // a vector of smart pointers to dynamically allocated B-curves
        std::vector<SP<BCurve3>::Default>       _bcurve;
        // a vector of smart pointers to dynamically allocated B-curve images (i.e., GenericCurve3 objects)
        std::vector<SP<GenericCurve3>::Default> _img_bcurve;

        // visibility flags
        bool                                    _show_arcs;
        bool                                    _show_control_polygons;
        bool                                    _show_tangent_vectors;

    public:
        // special and default constructor
        // the format specifies the properties of the rendering window
        GLWidget(QWidget* parent = 0, const QGLFormat& format = QGL::Rgba | QGL::DepthBuffer | QGL::DoubleBuffer);

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

        void setVisibilityOfArcs(bool value);
        void setVisibilityOfControlPolygons(bool value);
        void setVisibilityOfTangentVectors(bool value);
    };
}

#endif // GLWIDGET_H
