#pragma warning(disable:4503)
#include "YourGLWidget.h"

#include <iostream>
#include <fstream>

#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Math/Constants.h>

#include "../Spaces/SpecializedECSpaces.h"

using namespace cagd;
using namespace std;

(*@\Green{// your default/special constructor}@*)
void YourGLWidget::YourGLWidget(/*...*/)
{
    try
    {
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_color_shader:end} and \mref{src:GLWidgetTestingShaderPrograms:loading_shaders:end}--\mref{src:GLWidgetTestingShaderPrograms:constructor:first_call_of_updateTransformationMatrices:end} of Listing \mref{src:GLWidgetTestingShaderPrograms.cpp}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Red{// \textbf{Remark}} \label{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:start}@*)
        (*@\Red{//}@*)
        (*@\Red{// Several methods of classes ECSpace, BCurve3 and BSurface3 expect a boolean flag named}@*)
        (*@\Red{// check\_for\_ill\_conditioned\_matrices and a non-negative integer named expected\_correct\_significant\_digits.}@*)
        (*@\Red{// These variables appear, e.g.,\ in the input argument lists of the constructor of EC spaces and of the}@*)
        (*@\Red{// methods that perform either order elevation or subdivision on B-curves and B-surfaces.}@*)
        (*@\Red{//}@*)
        (*@\Red{// If the flag check\_for\_ill\_conditioned\_matrices is set to true, these methods will calculate the condition}@*)
        (*@\Red{// number of each matrix that appears in the construction of the unique normalized B-basis of the underlying}@*)
        (*@\Red{// EC space. Using singular value decomposition, each condition number is determined as the ratio of the}@*)
        (*@\Red{// largest and smallest singular values of the corresponding matrices.}@*)
        (*@\Red{//}@*)
        (*@\Red{// If at least one of the obtained condition numbers is too large, i.e., when the number of estimated}@*)
        (*@\Red{// correct significant digits is less than number of expected ones, these methods will throw an exception}@*)
        (*@\Red{// that states that one of the systems of linear equations is ill-conditioned and therefore its solution}@*)
        (*@\Red{// may be not accurate.}@*)
        (*@\Red{//}@*)
        (*@\Red{// If the user catches such an exception, one can try:}@*)
        (*@\Red{// \ \ 1) to lower the number of expected correct significant digits;}@*)
        (*@\Red{// \ \ 2) to decrease the dimension of the underlying EC space;}@*)
        (*@\Red{// \ \ 3) to change the endpoints of the definition domain $\left[\alpha, \beta\right]$;}@*)
        (*@\Red{// \ \ 4) to run the code without testing for ill-conditioned matrices and hope for the best.}@*)
        (*@\Red{//}@*)
        (*@\Red{// Note that the standard condition number may lead to an overly pessimistic estimate for the overall error}@*)
        (*@\Red{// and, by activating this boolean flag, the run-time of these methods will increase. Several numerical tests}@*)
        (*@\Red{// show that ill-conditioned matrices appear only when one defines EC spaces with relatively big dimensions.}@*)
        (*@\Red{// Considering that, in practice, curves and surfaces are mostly composed of smoothly joined lower order arcs}@*)
        (*@\Red{// and patches, by default we opted for speed, i.e., initially the flag check\_for\_ill\_conditioned\_matrices is set}@*)
        (*@\Red{// to false. Naturally, if one obtains mathematically or geometrically unexpected results, then one should}@*)
        (*@\Red{// (also) study the condition numbers mentioned above.}@*)
        GLboolean check_for_ill_conditioned_matrices  = GL_FALSE;
        GLint     expected_correct_significant_digits = 3; (*@\label{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:end}@*)

        (*@\Green{// Try to create different EC spaces:}@*)
        _space_count = 5;                   (*@\Green{// number of considered EC spaces;}@*)

        _space.resize(_space_count);        (*@\Green{// memory allocations;}@*)
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

        (*@\Green{// orders of the EC spaces (\mref{eq:pure_polynomial}), (\mref{eq:pure_trigonometric}), (\mref{eq:pure_hyperbolic}), (\mref{eq:algebraic_trigonometric}) and (\mref{eq:Mazure_vector_space}), respectively;}@*)
        _n[0] = 6; _n[1] = 3; _n[2] = 3; _n[3] = 2; _n[4] = 2;

        (*@\Green{// endpoints of the corresponding definition domains;}@*)
        _alpha[0] = 0.0, _beta[0] = 1.0;
        _alpha[1] = 0.0, _beta[1] = HALF_PI;
        _alpha[2] = 0.0, _beta[2] = PI;
        _alpha[3] = 0.0, _beta[3] = 4.0 * PI / 3.0;
        _alpha[4] = 0.0, _beta[4] = 16.694941067922716;

        (*@\Green{// additional shape parameters of the mixed algebraic-exponential-trigonometric EC space (\mref{eq:Mazure_vector_space});}@*)
        double a  = 1.0, b = 0.2;

        _space[0] = SP<ECSpace>::Default(   (*@\Green{// dynamical allocation of the EC space $\mathbb{P}_{6}^{0,1}$;}@*)
            new (nothrow) PolynomialECSpace(
              _alpha[0], _beta[0], _n[0],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[1] = SP<ECSpace>::Default(   (*@\Green{// dynamical allocation of the EC space $\mathbb{T}_{6}^{0,\frac{\pi}{2}}$;}@*)
            new (nothrow) TrigonometricECSpace(
              _alpha[1], _beta[1], _n[1],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[2] = SP<ECSpace>::Default(   (*@\Green{// dynamical allocation of the EC space $\mathbb{H}_{6}^{0,\pi}$;}@*)
            new (nothrow) HyperbolicECSpace(
              _alpha[2], _beta[2], _n[2],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[3] = SP<ECSpace>::Default(   (*@\Green{// dynamical allocation of the EC space $\mathbb{AT}_{8}^{0,\frac{4\pi}{3}}$;}@*)
            new (nothrow) ATECSpace(
              _alpha[3], _beta[3], _n[3],
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));
        _space[4] = SP<ECSpace>::Default(  (*@\Green{// dynamical allocation of the EC space $\mathbb{M}_{6,1,0.2}^{0,16.694941067922716}$.}@*)
          new (nothrow) BrilleaudMazureSpace(
              _alpha[4], _beta[4], _n[4], a, b,
              check_for_ill_conditioned_matrices, expected_correct_significant_digits));

        (*@\Green{// colors associated with different EC spaces}@*)
        for (GLint i = 0; i < _space_count; i++)
        {
            _color[i] = coldToHotColormap(i, 0, _space_count - 1);
        }

        (*@\Green{// we will evaluate the zeroth, first and second order derivatives at 50 uniform subdivision points}@*)
        _maximum_order_of_derivatives =  2;
        _div_point_count              = 50;

        (*@\Green{// determines the subdivision point of the definition domains, where the initial B-curves have to be subdivided}@*)
         _ratio                       = 0.5;

        (*@\Green{// labels associated with the left and right arcs of a subdivided B-curve}@*)
        string direction[2] = {"left", "right"};

        (*@\Green{// generating B-curves, performing order elevations and subdivisions on them}@*)
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

            (*@\Green{// the control points of each B-curve will be obtained by uniform sampling the unit circle}@*)
            GLdouble step = TWO_PI / _bcurve[i]->dataCount();

            for (GLint j = 0; j < _bcurve[i]->dataCount(); j++)
            {
                GLdouble u = j * step;

                Cartesian3 &reference = (*_bcurve[i])[j];

                reference[0] = -cos(u);
                reference[1] =  sin(u);
            }

            (*@\Green{// updating the VBOs of the $i$th control polygon}@*)
            if (!_bcurve[i]->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the " +
                                toString(i) + th + " control polygon!");
            }

            (*@\Green{// generating the image of the $i$th B-curve and updating its VBOs}@*)
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

            (*@\Green{// generating the $i$th order elevated B-curve: we insert the complex zero $z = 0$ of mulitplicity $n_i + 2$}@*)
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

            (*@\Green{// updating the VBOs of the $i$th order elevated control polygon}@*)
            if (!_oe_bcurve[i]->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the " + toString(i) + th +
                                " order elevated B-curve's control polygon!");
            }

            (*@\Green{// generating the image of the $i$th order elevated B-curve and updating its VBOs}@*)
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

            (*@\Green{// subdivision of the initial B-curves}@*)
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

        glClearColor(1.0, 1.0, 1.0, 1.0);       (*@\Green{// set the background color}@*)
        glEnable(GL_DEPTH_TEST);                (*@\Green{// enable depth testing}@*)
        glEnable(GL_LINE_SMOOTH);               (*@\Green{// enable the anti-aliasing of line primitives}@*)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_POINT_SMOOTH);              (*@\Green{// enable the anti-aliasing of point primitives}@*)
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    }
    catch (Exception &e)
    {
        cout << e << endl;
    }
}

GLboolean YourGLWidget::_updateTransformationMatrices()
{
    (*@\Green{// \ldots}@*)
    (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms:_updateTransformationMatrices:implementation:start}--\mref{src:GLWidgetTestingShaderPrograms:_updateTransformationMatrices:implementation:end} of Listing \mref{src:GLWidgetTestingShaderPrograms.cpp}}@*)
    (*@\Green{// \ldots}@*)
}

(*@\Green{// your rendering method}@*)
void YourGLWidget::render()
{
    (*@\Green{// clear the color and depth buffers}@*)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    _color_shader.enable();

    GLfloat x_min = -4.5f, x_step = 11.0f / _space_count;

    glLineWidth(2.0);
    glPointSize(8.0);
    for (GLint i = 0; i < _space_count; i++)
    {
        if (_bcurve[i] && _img_bcurve[i] &&
            _oe_bcurve[i] && _img_oe_bcurve[i] &&
            _subdivision[i])
        {
            (*@\Green{// rendering the control polygons and images of the initial B-curves}@*)
            _color_shader.setUniformColor("color", _color[i]);

            Translate        local_T(x_min + i * x_step, +2.5f, 0.0f);
            GLTransformation local_PVM = _PVM * local_T;
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            _bcurve[i]->renderData(_color_shader, GL_LINE_STRIP);
            _bcurve[i]->renderData(_color_shader, GL_POINTS);

            _img_bcurve[i]->renderDerivatives(_color_shader, 0, GL_LINE_STRIP);

            (*@\Green{// uncomment the line below in order to see the tangent vectors of the $i$th B-curve}@*)
            (*@\Green{// \_img\_bcurve[i]-$>$renderDerivatives(\_color\_shader, 1, GL\_LINES);}@*)

            (*@\Green{// uncomment the line below in order to see the acceleration vectors of the $i$th B-curve}@*)
            (*@\Green{// \_img\_bcurve[i]-$>$renderDerivatives(\_color\_shader, 2, GL\_LINES);}@*)

            (*@\Green{// rendering the control polygons and images of the order elevated B-curves}@*)
            local_T.setYDirectionalUnits(0.0f);
            local_PVM = _PVM * local_T;

            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            _oe_bcurve[i]->renderData(_color_shader, GL_LINE_STRIP);
            _oe_bcurve[i]->renderData(_color_shader, GL_POINTS);

            _img_oe_bcurve[i]->renderDerivatives(_color_shader, 0, GL_LINE_STRIP);

            (*@\Green{// for the sake of comparison, we re-render the control polygons of the initial B-curves}@*)
            _color_shader.setUniformColor("color", colors::gray);
            glEnable(GL_LINE_STIPPLE);
                glLineStipple(1, 0xf0f0);
                _bcurve[i]->renderData(_color_shader, GL_LINE_STRIP);
                _bcurve[i]->renderData(_color_shader, GL_POINTS);
            glDisable(GL_LINE_STIPPLE);

            (*@\Green{// rendering the control polygons and images of the subdivided arcs of the initial B-curves}@*)
            local_T.setYDirectionalUnits(-2.5f);
            local_PVM = _PVM * local_T;
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            for (GLint j = 0; j < 2; j++)
            {
                _color_shader.setUniformColor("color", _color[(i + j) % _space_count]);

                (*_subdivision[i])[j]->renderData(_color_shader, GL_LINE_STRIP);
                (*_subdivision[i])[j]->renderData(_color_shader, GL_POINTS);

                _img_subdivision[i][j]->renderDerivatives(
                        _color_shader, 0, GL_LINE_STRIP);
            }

            (*@\Green{// for the sake of comparison, we re-render the control polygons of the initial B-curves}@*)
            _color_shader.setUniformColor("color", colors::gray);
            glEnable(GL_LINE_STIPPLE);
                glLineStipple(1, 0xf0f0);
                _bcurve[i]->renderData(_color_shader, GL_LINE_STRIP);
                _bcurve[i]->renderData(_color_shader, GL_POINTS);
            glDisable(GL_LINE_STIPPLE);
        }
    }
    glPointSize(1.0);
    glLineWidth(1.0);

    (*@\Green{// revert to the original projection-view-model matrix}@*)
    _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

    (*@\Green{// render other geometries\ldots}@*)

    _color_shader.disable();

    (*@\Green{// Fig.\ \mref{fig:B-curve_operations} illustrates the image obtained at run-time.}@*)
}
