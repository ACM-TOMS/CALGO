#if defined(_MSC_VER)
#pragma warning(disable:4503)
#endif

#include "GLWidget.h"
#include <QTime>
#include <QMessageBox>
#include <iostream>

#include <Core/Exceptions.h>
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

            _S.setScalingFactors(3.2f, 3.2f, 3.2f);

            if (!_color_shader.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::VERTEX, "../../00_Dependencies/Shaders/color.vert"))
            {
                throw Exception("Could not attach the vertex shader of the color shader program!");
            }

            if (!_color_shader.attachNewShaderFromSourceFile(
                        ShaderProgram::Shader::FRAGMENT, "../../00_Dependencies/Shaders/color.frag"))
            {
                throw Exception("Could not attach the fragment shader of the color shader program!");
            }

            if (!_color_shader.linkAttachedShaders())
            {
                throw Exception("Could not link the attached shaders of the color shader program!");
            }

            // Try to update all required transformation matrices (the method _updateTransformationMatrices() has to
            // be called after every successful transformation related event handling).
            if (!_updateTransformationMatrices())
            {
                throw Exception("Could not update all transformation matrices!");
            }

            // creating different types of EC spaces and evaluating their normalized B-basis functions
            _space_count = 6;

            _space.resize(_space_count);
            _alpha.resize(_space_count);
            _beta.resize(_space_count);
            _n.resize(_space_count);
            _img_B_basis.resize(_space_count);

            // Test these parameters without checking for ill conditioned matrices.
            // Note that, in certain cases the standard condition number may lead to an overly
            // pessimistic estimate for the overall error.
            _n[0] = 10;
            _n[1] =  5;
            _n[2] =  7;
            _n[3] =  4;
            _n[4] =  3;
            _n[5] = 10;

            _alpha[0] = 0.0,        _beta[0] = 1.0;
            _alpha[1] = 0.0,        _beta[1] = HALF_PI;
            _alpha[2] = 0.0,        _beta[2] = PI;
            _alpha[3] = 0.0,        _beta[3] = TWO_PI;
            _alpha[4] = -1.25 * PI, _beta[4] = 1.25 * PI;
            _alpha[5] = 0.0,        _beta[5] = 16.694941067922716 / 2.0;

            // additional shape parameters of the mixed algebraic-exponential-trigonometric EC space $\mathbb{M}_{6}^{0,16.694941067922716}$
            double a  = 1.0, b = 0.2;

            bool check_for_ill_conditioned_matrices = false;
            int expected_correct_significant_digits = 1;

            _space[0] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{P}_{10}^{0,1}$;
                        new (nothrow) PolynomialECSpace(
                                _alpha[0], _beta[0], _n[0],
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

            _space[1] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{T}_{10}^{0,\frac{\pi}{2}}$;
                        new (nothrow) TrigonometricECSpace(
                                _alpha[1], _beta[1], _n[1],
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

            _space[2] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{H}_{14}^{0,\pi}$;
                        new (nothrow) HyperbolicECSpace(
                                _alpha[2], _beta[2], _n[2],
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

            _space[3] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{AT}_{24}^{0,2\pi}$;
                        new (nothrow) ATECSpace(
                                _alpha[3], _beta[3], _n[3],
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

            _space[4] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{AET}_{27}^{-5\pi/4,5\pi/4}$;
                        new (nothrow) AETECSpace(
                                _alpha[4], _beta[4], _n[4],
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

            _space[5] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{M}_{14,1,0.2}^{0,8.347470533961358}$
                        new (nothrow) BrilleaudMazureSpace(
                                _alpha[5], _beta[5], _n[5], a, b,
                                check_for_ill_conditioned_matrices,
                                expected_correct_significant_digits));

            _maximum_order_of_derivatives =   1;
            _div_point_count              = 200;

            for (GLuint i = 0; i < _space.size(); i++)
            {
                string th;

                switch (i % 5)
                {
                case 1:  th = "st"; break;
                case 2:  th = "nd"; break;
                case 3:  th = "rd"; break;
                default: th = "th"; break;
                }

                if (!_space[i])
                {
                    throw Exception("The " + toString(i) + th + " EC space could not be created!");
                }

                cout << "Dimension of the " << i << th << " EC space: " << _space[i]->dimension() << endl;
                cout << "Ordinary basis functions of the " << i << th << " EC space: " << endl;

                string expression;
                for (GLint j = 0; j < _space[i]->dimension(); j++)
                {
                    _space[i]->LaTeXExpression(j, expression);
                    cout << "\t" << expression << endl;
                }
                cout << endl;

                _img_B_basis[i] = SP< RowMatrix<SP<GenericCurve3>::Default> >::Default(
                            _space[i]->generateImagesOfAllBasisFunctions(
                            ECSpace::B_BASIS,
                            _maximum_order_of_derivatives, _div_point_count));

                if (!_img_B_basis[i])
                {
                    throw Exception("Could not generate the image of all normalized "
                                    "B-basis functions of the " + toString(i) +
                                    th + " EC space!");
                }

                RowMatrix<SP<GenericCurve3>::Default> &reference = *_img_B_basis[i];
                for (GLint j = 0; j < reference.columnCount(); j++)
                {
                    if (reference[j])
                    {
                        if (!reference[j]->updateVertexBufferObjects())
                        {
                            throw Exception("Could not update the VBOs of a normalized B-basis function!");
                        }
                    }
                }
            }

            // calculating the second order derivatives of the nth trigonometric ordinary and B-basis functions
            // at the midpoint of their definition domain
            GLdouble ordinary = (*_space[1])(ECSpace::ORDINARY_BASIS, _n[1], 2, (_alpha[1] + _beta[1]) / 2.0);
            GLdouble B_basis  = (*_space[1])(ECSpace::B_BASIS,        _n[1], 2, (_alpha[1] + _beta[1]) / 2.0);

            cout << "Second order derivative of the " << _n[1] << "th ordinary "
                    "trigonometric basis function: " << ordinary << endl;
            cout << "Second order derivative of the " << _n[1] << "th normalized "
                    "trigonometric B-basis function: " << B_basis << endl;

            glClearColor(1.0, 1.0, 1.0, 1.0);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_LINE_SMOOTH);
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
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
        // clears the color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        _color_shader.enable();
        _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

        for (GLint i = 0; i < (GLint)_space.size(); i++)
        {
            if (_img_B_basis[i])
            {
                const RowMatrix<SP<GenericCurve3>::Default> &reference = *_img_B_basis[i];

                Translate local_T(-(_alpha[i] + _beta[i]) / 2.0, 1.7 - 0.75 * i, 0.0);

                GLTransformation local_PVM = _PVM * local_T;
                _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

                for (GLint j = 0; j < reference.columnCount(); j++)
                {
                    _color_shader.setUniformColor("color", coldToHotColormap(j, 0, reference.columnCount() - 1));

                    reference[j]->renderDerivatives(0, GL_LINE_STRIP);
                    //reference[j]->renderDerivatives(1, GL_LINES);
                }
            }
        }

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
}
