#pragma warning(disable:4503)
#include "YourGLWidget.h"

#include <fstream>
#include <iostream>


#include <Core/Exceptions.h>
#include <Core/Utilities.h>
#include <Core/Math/Constants.h>
#include <Core/Utilities.h>

using namespace cagd;
using namespace std;

(*@\Green{// your default/special constructor}@*)
YourGLWidget::YourGLWidget(/*...*/)
{
    try
    {
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_color_shader:end} and \mref{src:GLWidgetTestingShaderPrograms:loading_shaders:end}--\mref{src:GLWidgetTestingShaderPrograms:constructor:first_call_of_updateTransformationMatrices:end} of Listing \mref{src:GLWidgetTestingShaderPrograms.cpp}}@*)
        (*@\Green{// \ldots}@*)

        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:start}--\mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp:constructor:remark:end} of Listing \mref{src:GLWidgetDefineEvaluateDifferentiateManipulateRenderBCurves.cpp}}@*)
        (*@\Green{// \ldots}@*)

        _alpha = 0.0;
        _beta  = 5.0 * PI / 6.0;
        _space = SP<ECSpace>::Default(new ECSpace(_alpha, _beta,
                                                  check_for_ill_conditioned_matrices,
                                                  expected_correct_significant_digits));

        _omega = 1.0 / 3.0 / PI;

        (*@\Green{// the order of the following zero insertions is important, since it determines the order of the ordinary}@*)
        (*@\Green{// basis functions and consequently the order of their coefficients in case of B-representations}@*)
        _space->insertZero(-_omega, 1.0, 1, false);
        _space->insertZero(+_omega, 1.0, 1, false);

        if (!_space->updateBothBases(check_for_ill_conditioned_matrices,
                                     expected_correct_significant_digits))
        {
            throw Exception("Could not update the normalized B-basis of the EC space!");
        }

        (*@\Green{// If no exception occured so far, one should have the ordinary basis functions $\varphi_{4,0}\left(u\right) \equiv 1,$}@*)
        (*@\Green{// $\varphi_{4,1}\left(u\right)=e^{-\omega u} \cos\left(u\right), ~\varphi_{4,2}\left(u\right)=e^{-\omega u} \sin\left(u\right),~\varphi_{4,3}\left(u\right)=e^{\omega u} \cos\left(u\right),~\varphi_{4,4}\left(u\right)=e^{\omega u} \sin\left(u\right),$}@*)
        (*@\Green{// where $u \in \left[0,\beta\right]$.}@*)

        (*@\Green{// memory allocations}@*)
        _arc_count = 7;
        _bcurve.resize(_arc_count);
        _img_bcurve.resize(_arc_count);

        (*@\Green{// we will evaluate the zeroth and first order derivatives at 30 uniform subdivision points}@*)
        _maximum_order_of_derivatives = 1;
        _div_point_count = 30;

        (*@\Green{// multi-threaded operations performed on the CPU-side (if one is convinced that there are no errors}@*)
        (*@\Green{// in the lines \mref{src:GLWigetLogarithmicSpiral:constructor:B_representation:start}--\mref{src:GLWigetLogarithmicSpiral:constructor:B_representation:end} below, variables B\_representations\_aborted and reason can be eliminated from}@*)
        (*@\Green{// the code in order to further increase the performance)}@*)
        bool    B_representations_aborted = false; (*@\label{src:GLWigetLogarithmicSpiral:constructor:B_representation:start}@*)
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

                (*@\Green{// a smart pointer to the B-representation of the $k$th arc of the logarithmic spiral (\mref{eq:logarithmic_spiral})}@*)
                _bcurve[k] = SP<BCurve3>::Default(new BCurve3(*_space));

                if (!_bcurve[k])
                {
                    B_representations_aborted = true;
                    #pragma omp flush (B_representations_aborted)
                    reason = "Could not create the " + toString(k) + th + " B-arc!";
                }
                else
                {
                    (*@\Green{// calculating the coefficients of the ordinary basis functions cf.\ formula (\mref{eq:logarithmic_spiral}), i.e.,}@*)
                    (*@\Green{// $\boldsymbol{\lambda}_i = \left[\begin{array}{c}0 \\ 0 \end{array}\right], ~i=0,1,2,~\boldsymbol{\lambda}_3 = e^{k\omega\beta}\left[\begin{array}{c}\cos\left(k\beta\right) \\ \sin\left(k\beta\right) \end{array}\right], ~\boldsymbol{\lambda}_4 = e^{k\omega\beta}\left[\begin{array}{r}-\sin\left(k\beta\right) \\ \cos\left(k\beta\right) \end{array}\right]$}@*)
                    GLdouble scale = exp(k * _omega * _beta);

                    RowMatrix<Cartesian3> lambda(_space->dimension());
                    lambda[3][0]  =  scale * cos(k * _beta);
                    lambda[3][1]  =  scale * sin(k * _beta);
                    lambda[4][0]  = -lambda[3][1];
                    lambda[4][1]  =  lambda[3][0];

                    (*@\Green{// updating the control points for B-representation}@*)
                    if (!_bcurve[k]->updateControlPointsForExactDescription(lambda))
                    {
                        B_representations_aborted = true;
                        #pragma omp flush (B_representations_aborted)
                        reason = "Could not perform the " + toString(k) + th +
                                 " B-representation!";
                    }
                    else
                    {
                        (*@\Green{// generating the image of the B-curve and updating its VBOs}@*)
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
        } (*@(\label{src:GLWigetLogarithmicSpiral:constructor:B_representation:end}@*)

        (*@\Green{// sequential operations performed on the GPU-side}@*)
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

            (*@\Green{// updating the VBO of the $k$th control polygon}@*)
            if (!_bcurve[k]->updateVertexBufferObjectsOfData())
            {
                throw Exception("Could not update the VBO of the control polygon of the "
                                + toString(k) + th + " B-curve!");
            }

            (*@\Green{// updating the VBOs of the $k$th B-curve's image}@*)
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

    (*@\Green{// set the projection-view-model matrix}@*)
    _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

    (*@\Green{// render the control polygons and images of the B-representations}@*)
    glLineWidth(2.0);
    glPointSize(8.0);
    for (GLint k = 0; k < _arc_count; k++)
    {
        if (_bcurve[k] && _img_bcurve[k])
        {
            _color_shader.setUniformColor(
                    "color", coldToHotColormap(k, 0, _arc_count - 1));

            _bcurve[k]->renderData(_color_shader, GL_LINE_STRIP);
            _bcurve[k]->renderData(_color_shader, GL_POINTS);

            _img_bcurve[k]->renderDerivatives(_color_shader, 0, GL_LINE_STRIP);
            _img_bcurve[k]->renderDerivatives(_color_shader, 1, GL_LINES);
        }
    }
    glPointSize(1.0);
    glLineWidth(1.0);

    (*@\Green{// render other geometries\ldots}@*)

    _color_shader.disable();

    (*@\Green{// Fig.\ \mref{fig:logarithmic_spiral} illustrates the image obtained at run-time.}@*)
}
