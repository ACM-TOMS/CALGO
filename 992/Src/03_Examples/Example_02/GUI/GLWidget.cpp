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

using namespace cagd;
using namespace std;

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

        _S.setScalingFactors(2.5f, 2.5f, 2.5f);

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

        // Try to update all required transformation matrices (the method _updateTransformationMatrices() has to}\label{src:GLWidgetTestingShaderPrograms:constructor:first_call_of_updateTransformationMatrices:start
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
        _space_count = 5;                   // number of considered EC spaces;

        _space.resize(_space_count);        // memory allocations;
        _alpha.resize(_space_count);
        _beta.resize(_space_count);
        _n.resize(_space_count);

        _color.resize(_space_count);

        _bcurve.resize(_space_count);
        _img_bcurve.resize(_space_count);

        _oe_bcurve.resize(_space_count);
        _img_oe_bcurve.resize(_space_count);

        _subdivision.resize(_space_count);
        _img_subdivision.resize(_space_count);

        // orders of the EC spaces
        _n[0] = 6; _n[1] = 3; _n[2] = 3; _n[3] = 2; _n[4] = 2;

        // endpoints of the corresponding definition domains;
        _alpha[0] = 0.0, _beta[0] = 1.0;
        _alpha[1] = 0.0, _beta[1] = HALF_PI;
        _alpha[2] = 0.0, _beta[2] = PI;
        _alpha[3] = 0.0, _beta[3] = 4.0 * PI / 3.0;
        _alpha[4] = 0.0, _beta[4] = 16.694941067922716;

        // additional shape parameters of the mixed algebraic-exponential-trigonometric EC space $\mathbb{M}_{6}^{0,16.694941067922716}$
        double a  = 1.0, b = 0.2;

        _space[0] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{P}_{6}^{0,1}$;
            new (nothrow) PolynomialECSpace(
              _alpha[0], _beta[0], _n[0],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[1] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{T}_{6}^{0,\frac{\pi}{2}}$;
            new (nothrow) TrigonometricECSpace(
              _alpha[1], _beta[1], _n[1],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[2] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{H}_{6}^{0,\pi}$;
            new (nothrow) HyperbolicECSpace(
              _alpha[2], _beta[2], _n[2],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[3] = SP<ECSpace>::Default(   // dynamical allocation of the EC space $\mathbb{AT}_{8}^{0,\frac{4\pi}{3}}$;
            new (nothrow) ATECSpace(
              _alpha[3], _beta[3], _n[3],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[4] = SP<ECSpace>::Default(  // dynamical allocation of the EC space $\mathbb{M}_{6,1,0.2}^{0,16.694941067922716}$
            new (nothrow) BrilleaudMazureSpace(
                    _alpha[4], _beta[4], _n[4], a, b,
                    check_for_ill_conditioned_matrices, expected_correct_significant_digits));

        // colors associated with different EC spaces
        for (int i = 0; i < _space_count; i++)
        {
            _color[i] = coldToHotColormap(i * 1.1, 0, _space_count - 1);
        }

        // we will evaluate the zeroth, first and second order derivatives at 50 uniform subdivision points
        _maximum_order_of_derivatives = 2;
        _div_point_count              = 50;

        // determines the subdivision point of the definition domains, where the initial B-curves have to be subdivided
         _ratio                       = 0.5;

        // labels associated with the left and right arcs of a subdivided B-curve
        string direction[2] = {"left", "right"};

        // generating B-curves, performing order elevations and subdivisions on them
        for (GLint i = 0; i < _space_count; i++)
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
                throw Exception("The " + toString(i) + th +
                                " EC space could not be created!");
            }

            _bcurve[i] = SP<BCurve3>::Default(new (nothrow) BCurve3(*_space[i]));

            if (!_bcurve[i])
            {
                throw Exception("Could not genereate the " + toString(i) + th +
                                " EC B-curve!");
            }

            // the control points of each B-curve will be obtained by uniform sampling the unit circle
            GLdouble step = TWO_PI / _bcurve[i]->dataCount();

            for (GLint j = 0; j < _bcurve[i]->dataCount(); j++)
            {
                GLdouble u = j * step;

                Cartesian3 &reference = (*_bcurve[i])[j];

                reference[0] = -cos(u);
                reference[1] =  sin(u);
            }

            // updating the VBOs of the $i$th control polygon
            if (!_bcurve[i]->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the " +
                                toString(i) + th + " control polygon!");
            }

            // generating the image of the $i$th B-curve and updating its VBOs
            _img_bcurve[i] = SP<GenericCurve3>::Default(
                    _bcurve[i]->generateImage(
                            _maximum_order_of_derivatives, _div_point_count));

            if (!_img_bcurve[i])
            {
                throw Exception("Could not generate the image of the " +
                                toString(i) + th + " B-curve!");
            }

            if (!_img_bcurve[i]->updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBOs of the " +
                                toString(i) + th + " B-curve's image!");
            }

            // generating the $i$th order elevated B-curve: we insert the complex zero $z = 0$ of mulitplicity $n_i + 2$
            _oe_bcurve[i] = SP<BCurve3>::Default(
                    _bcurve[i]->performOrderElevation(
                            0.0, 0.0, _n[i] + 2,
                            check_for_ill_conditioned_matrices,
                            expected_correct_significant_digits));

            if (!_oe_bcurve[i])
            {
                throw Exception("Could not perform order elevation on the " +
                                toString(i) + th + " B-curve!");
            }

            // updating the VBOs of the $i$th order elevated control polygon
            if (!_oe_bcurve[i]->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the " + toString(i) + th +
                                " order elevated B-curve's control polygon!");
            }

            // generating the image of the $i$th order elevated B-curve and updating its VBOs
            _img_oe_bcurve[i] = SP<GenericCurve3>::Default(
                    _oe_bcurve[i]->generateImage(
                            _maximum_order_of_derivatives, _div_point_count));

            if (!_img_oe_bcurve[i])
            {
                throw Exception("Could not generate the image of the " +
                                toString(i) + th + " order elevated B-curve!");
            }

            if (!_img_oe_bcurve[i]->updateVertexBufferObjects())
            {
                throw Exception("Could not update the VBOs of the "  +
                                toString(i) + th + " order elevated B-curve's image!");
            }

            // subdivision of the initial B-curves
            _subdivision[i] = SP< RowMatrix<SP<BCurve3>::Default> >::Default(
                    _bcurve[i]->performSubdivision(
                        _alpha[i] + _ratio * (_beta[i] - _alpha[i]),
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits));

            if (!_subdivision[i])
            {
                throw Exception("Could not perform subdivision on the " +
                                toString(i) + th + " B-curve!");
            }

            _img_subdivision[i].resizeColumns(2);

            for (GLint j = 0; j < 2; j++)
            {
                if (!(*_subdivision[i])[j]->updateVertexBufferObjectsOfData())
                {
                    throw Exception("Could not update the VBO of the control polygon "
                                    "of the " + direction[j] + " arc of the " +
                                    toString(i) + th + " subdivided B-curve!");
                }

                _img_subdivision[i][j] = SP<GenericCurve3>::Default(
                        (*_subdivision[i])[j]->generateImage(
                            _maximum_order_of_derivatives, _div_point_count));

                if (!_img_subdivision[i][j])
                {
                    throw Exception("Could not generate the image of the " +
                                    direction[j] + " arc of the " + toString(i) + th +
                                    " subdivided B-curve!");
                }

                if (!_img_subdivision[i][j]->updateVertexBufferObjects())
                {
                    throw Exception("Could not update the VBOs of the image of the " +
                                    direction[j] + " arc of the " + toString(i) + th +
                                    " subdivided B-curve!");
                }
            }
        }

        glClearColor(1.0, 1.0, 1.0, 1.0);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_POINT_SMOOTH);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
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

    GLfloat x_min = -4.5f, x_step = 11.0f / _space_count;

//    glLineWidth(2.0);
    glPointSize(8.0);
    for (GLint i = 0; i < _space_count; i++)
    {
        if (_bcurve[i] && _img_bcurve[i] &&
            _oe_bcurve[i] && _img_oe_bcurve[i] &&
            _subdivision[i])
        {
            // rendering the control polygons and images of the initial B-curves
            _color_shader.setUniformColor("color", _color[i]);

            Translate        local_T(x_min + i * x_step, +2.5f, 0.0f);
            GLTransformation local_PVM = _PVM * local_T;
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            _bcurve[i]->renderData(GL_LINE_STRIP);
            _bcurve[i]->renderData(GL_POINTS);

            _img_bcurve[i]->renderDerivatives(0, GL_LINE_STRIP);

            // uncomment the line below in order to see the tangent vectors of the $i$th B-curve
            //_img_bcurve[i]->renderDerivatives(1, GL_LINES);

            // uncomment the line below in order to see the acceleration vectors of the $i$th B-curve
            //_img_bcurve[i]->renderDerivatives(2, GL_LINES);

            // rendering the control polygons and images of the order elevated B-curves
            local_T.setYDirectionalUnits(0.0f);
            local_PVM = _PVM * local_T;

            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            _oe_bcurve[i]->renderData(GL_LINE_STRIP);
            _oe_bcurve[i]->renderData(GL_POINTS);

            _img_oe_bcurve[i]->renderDerivatives(0, GL_LINE_STRIP);

            // for the sake of comparison, we re-render the control polygons of the initial B-curves
            _color_shader.setUniformColor("color", colors::gray);
            glEnable(GL_LINE_STIPPLE);
                glLineStipple(1, 0xf0f0);
                _bcurve[i]->renderData(GL_LINE_STRIP);
                _bcurve[i]->renderData(GL_POINTS);
            glDisable(GL_LINE_STIPPLE);

            // rendering the control polygons and images of the subdivided arcs of the initial B-curves
            local_T.setYDirectionalUnits(-2.5f);
            local_PVM = _PVM * local_T;
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            for (GLint j = 0; j < 2; j++)
            {
                _color_shader.setUniformColor("color", _color[(i + j) % _space_count]);

                (*_subdivision[i])[j]->renderData(GL_LINE_STRIP);
                (*_subdivision[i])[j]->renderData(GL_POINTS);

                _img_subdivision[i][j]->renderDerivatives(0, GL_LINE_STRIP);
            }

            // for the sake of comparison, we re-render the control polygons of the initial B-curves
            _color_shader.setUniformColor("color", colors::gray);
            glEnable(GL_LINE_STIPPLE);
                glLineStipple(1, 0xf0f0);
                _bcurve[i]->renderData(GL_LINE_STRIP);
                _bcurve[i]->renderData(GL_POINTS);
            glDisable(GL_LINE_STIPPLE);
        }
    }
    glPointSize(1.0);
//    glLineWidth(1.0);

    // revert to the original projection-view-model matrix
    _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

    // \ldots
    // render other geometries
    // \ldots

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
