#ifndef GLWIDGET_H
#define GLWIDGET_H

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
        GLboolean                           _updateTransformationMatrices();

        // color shader program object
        ShaderProgram                       _color_shader;

        // parameters of different EC spaces
        GLint                               _space_count;  // number of considered EC spaces
        std::vector<GLdouble>               _alpha, _beta; // definition domain endpoints for each EC space
        std::vector<GLint>                  _n;            // order of each EC space
        std::vector<SP<ECSpace>::Default>   _space;        // smart pointers to different EC spaces

        GLint _div_point_count;              // number of subdivision points for basis function image generations
        GLint _maximum_order_of_derivatives; // specifies the maximum order of derivatives that have to be evaluated along the basis functions

        // an array of smart pointers that point to row matrices that store the image pointers of all non-negative normalized B-basis functions of each EC spaces
        std::vector< SP< RowMatrix<SP<GenericCurve3>::Default> >::Default > _img_B_basis;

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
    };
}

#endif // GLWIDGET_H
