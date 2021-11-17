#if defined(_MSC_VER)
#pragma warning(disable:4503)
#endif

#include "GLWidget.h"

#include <QMessageBox>

#include <fstream>
#include <iostream>

#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Math/Constants.h>
#include <Core/Utilities.h>

#include "../Spaces/SpecializedECSpaces.h"

using namespace std;

namespace cagd
{
    //--------------------------------
    // special and default constructor
    //--------------------------------
    GLWidget::GLWidget(QWidget *parent, const QGLFormat &format): QGLWidget(format, parent)
    {
    }

    //--------------------------------------------------------------------------------------
    // this virtual function is called once before the first call to paintGL() or resizeGL()
    //--------------------------------------------------------------------------------------
    void GLWidget::initializeGL()
    {
        try
        {
            // initializing the OpenGL Extension Wrangler library
            if (glewInit() != GLEW_OK)
            {
                throw Exception("Could not initialize the OpenGL Extension Wrangler Library!");
            }

            if (!platformIsSupported())
            {
                throw Exception("The platform is not supported!\n\n"
                                "We assume that the user has: \n"
                                "- a multi-core CPU; and \n"
                                "- a graphics adapter managed by a driver that supports at least OpenGL 3.0.");
            }


            // creating an orthogonal perspective projection matrix
            _aspect = (float)width() / (float)height();
            _left   = _bottom = -10.0f;
            _right  = _top    = +10.0f;
            _near   = -20.0f;
            _far    = +20.0f;

            _P = SP<OrthogonalProjection>::Default(
                        new (nothrow) OrthogonalProjection(_aspect,
                                                 _left, _right,
                                                 _bottom, _top,
                                                 _near, _far));

            if (!_P)
            {
                throw Exception("Could not create orthogonal projection matrix!");
            }

            // creating a view (or world) transformation matrix
            _eye[0]    = _eye[1]    = 0.0, _eye[2]    = 6.0;
            _center[0] = _center[1] =      _center[2] = 0.0;
            _up[0]     = _up[2]     = 0.0, _up[1]     = 1.0;

            _V = SP<LookAt>::Default(new (nothrow) LookAt(_eye, _center, _up));

            if (!_V)
            {
                throw Exception("Could not create the view/world transformation matrix!");
            }

            // specifying the axes of rotation matrices
            _Rx.setDirection(Cartesian3(1.0, 0.0, 0.0));
            _Ry.setDirection(Cartesian3(0.0, 1.0, 0.0));
            _Rz.setDirection(Cartesian3(0.0, 0.0, 1.0));

            // If the variable logging_is_enabled is set to true, one obtains logging information at run-time about
            // creating, loading, compiling, attaching and linking different types of shaders.
            // If one is convinced that the applied shaders do not contain bugs, this logging mechanism can be deactivated.
            GLboolean logging_is_enabled = GL_FALSE;

            // By default, logging information appear in the standard console output, but they can be redirected
            // to any output streams as it is shown in the following lines.
            fstream   shader_log;

            if (logging_is_enabled)
            {
                shader_log.open("shader.log", ios_base::out);
            }

            if (!_color_shader.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX, "../../00_Dependencies/Shaders/color.vert",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the vertex shader of the "
                                "color shader program!");
            }

            if (!_color_shader.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT, "../../00_Dependencies/Shaders/color.frag",
                        logging_is_enabled, shader_log))
            {
                throw Exception("Could not attach the fragment shader of the "
                                "color shader program!");
            }

            if (!_color_shader.linkAttachedShaders(logging_is_enabled, shader_log))
            {
                throw Exception("Could not link the attached shaders of the "
                                "color shader program!");
            }

            if (logging_is_enabled)
            {
                shader_log.close();
            }

            // Try to update all required transformation matrices (the method _updateTransformationMatrices() has to
            // be called after every successful transformation related event handling).
            if (!_updateTransformationMatrices())
            {
                throw Exception("Could not update all transformation matrices!");
            }

            // Remark
            //
            // Several methods of classes ECSpace, BCurve3 and BSurface3 expect a boolean flag named
            // check_for_ill_conditioned_matrices and a non-negative integer named expected_correct_significant_digits.
            // These variables appear e.g.\ in the input argument lists of the constructor of EC spaces and of the methods
            // that perform either order elevation or subdivision on B-curves and B-surfaces.
            //
            // If the flag check_for_ill_conditioned_matrices is set to true, these methods will calculate the condition
            // number of each matrix that appears in the construction of the unique normalized B-basis of the underlying
            // EC space. Using singular value decomposition, each condition number is determined as the as the ratio
            // of the largest and smallest singular values of the corresponding matrices.
            //
            // If at least one of the obtained condition numbers is too large, i.e., when the number of estimated
            // correct significant digits is less than number of expected ones, these methods will throw an exception
            // that states that one of the systems of linear equations is ill-conditioned and therefore its solution
            // may be not accurate.
            //
            // If the user catches such an exception, one can try:
            //      1) to lower the number of expected correct significant digits;
            //      2) to decrease the dimension of the underlying EC space;
            //      3) to change the endpoints of the definition domain $\left[\alpha, \beta\right]$;
            //      4) to run the code without testing for ill-conditioned matrices and hope for the best.
            //
            // Note that the standard condition number may lead to an overly pessimistic estimate for the overall error
            // and, by activating this boolean flag, the run-time of these methods will increase. Several numerical tests
            // show that ill-conditioned matrices appear only when one defines EC spaces with relatively big dimensions.
            // Considering that, in practice, curves and surfaces are mostly composed of smoothly joined lower order arcs
            // and patches, by default we opted for speed, i.e., initially the flag check_for_ill_conditioned_matrices is set
            // to false. Naturally, if one obtains mathematically or geometrically unexpected results, then one should
            // (also) study the condition numbers mentioned above.
            GLboolean check_for_ill_conditioned_matrices  = GL_FALSE;
            GLint     expected_correct_significant_digits = 3;

            // Try to create different EC spaces:
            _alpha = 0.0;
            _beta  = 5.0 * PI / 6.0;
            _space = SP<ECSpace>::Default(
                        new ECSpace(
                            _alpha, _beta,
                            check_for_ill_conditioned_matrices, expected_correct_significant_digits));

            _omega = 1.0 / 3.0 / PI;

            // the order of the following zero insertions is important, since it determines the order of the ordinary
            // basis functions and consequently the order of their coefficients in case of B-representations
            _space->insertZero(-_omega, 1.0, 1, false);
            _space->insertZero(+_omega, 1.0, 1, false);

            if (!_space->updateBothBases(check_for_ill_conditioned_matrices,
                                         expected_correct_significant_digits))
            {
                throw Exception("Could not update the normalized B-basis of the EC space!");
            }

            // memory allocations
            _arc_count = 7;
            _bcurve.resize(_arc_count);
            _img_bcurve.resize(_arc_count);

            // we will evaluate the zeroth and first order derivatives at 30 uniform subdivision points
            _maximum_order_of_derivatives =  2;
            _div_point_count              = 30;

            // multi-threaded operations performed on the CPU-side
            bool    B_representations_aborted = false;
            string  reason;

            #pragma omp parallel for
            for (GLint k = 0; k < _arc_count; k++)
            {
                #pragma omp flush (B_representations_aborted)
                if (!B_representations_aborted)
                {
                    string th;

                    switch (k % 5)
                    {
                    case 1:  th = "st"; break;
                    case 2:  th = "nd"; break;
                    case 3:  th = "rd"; break;
                    default: th = "th"; break;
                    }

                    // a smart pointer to the B-representation of the $k$th arc of the logarithmic spiral (\ref{eq:logarithmic_spiral})
                    _bcurve[k] = SP<BCurve3>::Default(new BCurve3(*_space));

                    if (!_bcurve[k])
                    {
                        B_representations_aborted = true;
                        #pragma omp flush (B_representations_aborted)
                        reason = "Could not create the " + toString(k) + th + " B-arc!";
                    }
                    else
                    {
                        // calculating the coefficients of the ordinary basis functions cf.\ formula (\ref{eq:logarithmic_spiral})
                        GLdouble scale = exp(k * _omega * _beta);

                        RowMatrix<Cartesian3> lambda(_space->dimension());
                        lambda[3][0]  =  scale * cos(k * _beta);
                        lambda[3][1]  =  scale * sin(k * _beta);
                        lambda[4][0]  = -lambda[3][1];
                        lambda[4][1]  =  lambda[3][0];

                        // updating the control points for B-representation
                        if (!_bcurve[k]->updateControlPointsForExactDescription(lambda))
                        {
                            B_representations_aborted = true;
                            #pragma omp flush (B_representations_aborted)
                            reason = "Could not perform the " + toString(k) + th +
                                     " B-representation!";
                        }
                        else
                        {
                            // generating the image of the B-curve and updating its VBOs
                            _img_bcurve[k] = SP<GenericCurve3>::Default(
                                    _bcurve[k]->generateImage(
                                            _maximum_order_of_derivatives, _div_point_count));

                            if (!_img_bcurve[k])
                            {
                                B_representations_aborted = true;
                                #pragma omp flush (B_representations_aborted)
                                reason = "Could not generate the image of the " +
                                         toString(k) + th + " B-representation!";
                            }
                        }
                    }
                }
            }

            if (B_representations_aborted)
            {
                throw Exception(reason);
            }

            // sequential operations performed on the GPU-side
            for (GLint k = 0; k < _arc_count; k++)
            {
                string th;

                switch (k % 5)
                {
                case 1:  th = "st"; break;
                case 2:  th = "nd"; break;
                case 3:  th = "rd"; break;
                default: th = "th"; break;
                }

                // updating the VBO of the $k$th control polygon
                if (!_bcurve[k]->updateVertexBufferObjectsOfData())
                {
                    throw Exception("Could not update the VBO of the control polygon of the "
                                    + toString(k) + th + " B-curve!");
                }

                // updating the VBOs of the $k$th B-curve's image
                if (!_img_bcurve[k]->updateVertexBufferObjects())
                {
                    throw Exception("Could not update the VBOs of the image of the " +
                                    toString(k) + th + " B-curve!");
                }
            }

            glClearColor(1.0, 1.0, 1.0, 1.0);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_LINE_SMOOTH);
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
            glEnable(GL_POINT_SMOOTH);
            glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

            // default values of visibility flags
            _show_arcs             = true;
            _show_control_polygons = true;
            _show_tangent_vectors  = false;
        }
        catch (Exception &e)
        {
            cout << e << endl;
            QMessageBox::critical(this, "An exception occured...", QString::fromStdString(e.reason()));
            exit(0);
        }
    }

    GLboolean GLWidget::_updateTransformationMatrices()
    {
        if (_P && _V)
        {
            _M              = _Rx * _Ry * _Rz * _T * _S;
            _VM             = (*_V) * _M;
            _PVM            = (*_P) * _VM;

            bool invertible = false;
            _tN             = _VM.inverse(&invertible);

            return invertible;
        }

        return GL_FALSE;
    }

    //-----------------------
    // the rendering function
    //-----------------------
    void GLWidget::paintGL()
    {
        // clear the color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        _color_shader.enable();

        // revert to the original projection-view-model matrix
        _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

        // rendering the control polygons and images of the B-representations
//        glLineWidth(2.0);
        glPointSize(8.0);
        for (GLint k = 0; k < _arc_count; k++)
        {
            if (_bcurve[k] && _img_bcurve[k])
            {
                _color_shader.setUniformColor(
                        "color", coldToHotColormap(k, 0, _arc_count - 1));

                if (_show_control_polygons)
                {
                    _bcurve[k]->renderData(_color_shader, GL_LINE_STRIP);
                    _bcurve[k]->renderData(_color_shader, GL_POINTS);
                }

                if (_show_arcs)
                {
                    _img_bcurve[k]->renderDerivatives(_color_shader, 0, GL_LINE_STRIP);
                }

                if (_show_tangent_vectors)
                {
                    _img_bcurve[k]->renderDerivatives(_color_shader, 1, GL_LINES);
                }
            }
        }
        glPointSize(1.0);
//        glLineWidth(1.0);

        _color_shader.disable();
    }

    //----------------------------------------------------------------------------
    // when the main window is resized one needs to redefine the projection matrix
    //----------------------------------------------------------------------------
    void GLWidget::resizeGL(int w, int h)
    {
        // setting the new size of the rendering context
        glViewport(0, 0, w, h);

        if (_P)
        {
            _P->setAspectRatio((float)w / (float)h);

            _PVM = (*_P) * _VM;
            _tN = _VM.inverse();
        }

        updateGL();
    }

    //-----------------------------------
    // implementation of the public slots
    //-----------------------------------

    void GLWidget::setAngleX(int value)
    {
        _Rx.setAngle(value * DEG_TO_RADIAN);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setAngleY(int value)
    {
        _Ry.setAngle(value * DEG_TO_RADIAN);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setAngleZ(int value)
    {
        _Rz.setAngle(value * DEG_TO_RADIAN);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setZoomFactor(double value)
    {
        _S.setScalingFactors(value, value, value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setTransX(double value)
    {
        _T.setXDirectionalUnits(value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setTransY(double value)
    {
        _T.setYDirectionalUnits(value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setTransZ(double value)
    {
        _T.setZDirectionalUnits(value);
        _updateTransformationMatrices();
        updateGL();
    }

    void GLWidget::setVisibilityOfArcs(bool value)
    {
        if (_show_arcs != value)
        {
            _show_arcs = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfControlPolygons(bool value)
    {
        if (_show_control_polygons != value)
        {
            _show_control_polygons = value;
            updateGL();
        }
    }

    void GLWidget::setVisibilityOfTangentVectors(bool value)
    {
        if (_show_tangent_vectors != value)
        {
            _show_tangent_vectors = value;
            updateGL();
        }
    }
}
