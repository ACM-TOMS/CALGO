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
YourGLWidget::YourGLWidget()
{
    try
    {
        (*@\Green{// \ldots}@*)
        (*@\Green{// these lines coincide with the lines \mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:part_1:start}--\mref{src:GLWidgetTestingShaderPrograms.cpp:constructor:loading_color_shader:end} and \mref{src:GLWidgetTestingShaderPrograms:loading_shaders:end}--\mref{src:GLWidgetTestingShaderPrograms:constructor:first_call_of_updateTransformationMatrices:end} of Listing \mref{src:GLWidgetTestingShaderPrograms.cpp}}@*)
        (*@\Green{// \ldots}@*)


        (*@\Green{// Try to create different EC spaces:}@*)
        _space_count = 6;                   (*@\Green{// number of considered EC spaces;}@*)

        _space.resize(_space_count);        (*@\Green{// memory allocations;}@*)
        _alpha.resize(_space_count);
        _beta.resize(_space_count);
        _n.resize(_space_count);
        _img_B_basis.resize(_space_count);

        _n[0] = 10;             (*@\Green{// degree of the pure polynomial EC space (\mref{eq:pure_polynomial});}@*)
        _n[1] =  5;             (*@\Green{// order of the pure trigonometric EC space (\mref{eq:pure_trigonometric});}@*)
        _n[2] =  7;             (*@\Green{// order of the pure hyperbolic EC space (\mref{eq:pure_hyperbolic});}@*)
        _n[3] =  4;             (*@\Green{// order of the mixed algebraic-trigonometric EC space (\mref{eq:algebraic_trigonometric});}@*)
        _n[4] =  3;             (*@\Green{// order of the mixed algebraic-exponential-trigonometric EC space (\mref{eq:algebraic-exponential-trigonometric});}@*)
        _n[5] = 10;             (*@\Green{// order of the mixed algebraic-exponential-trigonometric EC space (\mref{eq:Mazure_vector_space});}@*)

        _alpha[0] = 0.0,        _beta[0] = 1.0;
        _alpha[1] = 0.0,        _beta[1] = HALF_PI;
        _alpha[2] = 0.0,        _beta[2] = PI;       (*@\Green{// endpoints of corresponding definition domains;}@*)
        _alpha[3] = 0.0,        _beta[3] = TWO_PI;
        _alpha[4] = -1.25 * PI, _beta[4] = 1.25 * PI;
        _alpha[4] = 0.0,        _beta[4] = 16.694941067922716 / 2.0;

        (*@\Green{// additional shape parameters of the mixed algebraic-exponential-trigonometric EC space (\mref{eq:Mazure_vector_space});}@*)
        double a  = 1.0, b = 0.2;

        _space[0] = SP<ECSpace>::Default( (*@\Green{// dynamical allocation of the EC space $\mathbb{P}_{10}^{0,1}$;}@*)
                new (nothrow) PolynomialECSpace(_alpha[0], _beta[0], _n[0]));
        _space[1] = SP<ECSpace>::Default( (*@\Green{// dynamical allocation of the EC space $\mathbb{T}_{10}^{0,\frac{\pi}{2}}$;}@*)
                new (nothrow) TrigonometricECSpace(_alpha[1], _beta[1], _n[1]));
        _space[2] = SP<ECSpace>::Default( (*@\Green{// dynamical allocation of the EC space $\mathbb{H}_{14}^{0,\pi}$;}@*)
                new (nothrow) HyperbolicECSpace(_alpha[2], _beta[2], _n[2]));
        _space[3] = SP<ECSpace>::Default( (*@\Green{// dynamical allocation of the EC space $\mathbb{AT}_{24}^{0,2\pi}$;}@*)
                new (nothrow) ATECSpace(_alpha[3], _beta[3], _n[3]));
        _space[4] = SP<ECSpace>::Default( (*@\Green{// dynamical allocation of the EC space $\mathbb{AET}_{27}^{-5\pi/4,\,5\pi/4}$;}@*)
                new (nothrow) AETECSpace(_alpha[4], _beta[4], _n[4]));
        _space[5] = SP<ECSpace>::Default( (*@\Green{// dynamical allocation of the EC space $\mathbb{M}_{14,1,0.2}^{0,8.347470533}$.}@*)
                new (nothrow) BrilleaudMazureSpace(_alpha[5], _beta[5], _n[5], a, b));

        (*@\Green{// Parameters used for the image generation of the normalized B-basis functions:}@*)

        (*@\Green{// the common number of subdivision points in each definition domain;}@*)
        _div_point_count = 200;

        (*@\Green{// the maximum order of derivatives that have to be evaluated along the basis functions;}@*)
        _maximum_order_of_derivatives = 1;

        (*@\Green{// generate dynamicaly the row matrices that store the image smart pointers of all normalized}@*)
        (*@\Green{// B-basis functions of each EC space.}@*)
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
                throw Exception("The " + toString(i) + th +
                                " EC space could not be created!");
            }

            _img_B_basis[i] = SP< RowMatrix<SP<GenericCurve3>::Default> >::Default(
                _space[i]->generateImagesOfAllBasisFunctions(
                    ECSpace::B_BASIS, _maximum_order_of_derivatives, _div_point_count)); (*@\label{src:GLWidgetDifferentiatingAndRenderingBasisFunctions.cpp:generateImagesOfAllBasisFunctions:B_BASIS}@*)

            (*@\Red{// \textbf{Remark}}@*)
            (*@\Red{// If instead of the normalized B-basis functions one intends to generate the images of all}@*)
            (*@\Red{// ordinary basis functions, then in line \mref{src:GLWidgetDifferentiatingAndRenderingBasisFunctions.cpp:generateImagesOfAllBasisFunctions:B_BASIS} of the current listing one should use the constant}@*)
            (*@\Red{// ECSpace::ORDINARY\_BASIS instead of ECSpace::B\_BASIS.}@*)

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
                        throw Exception("Could not update the VBOs of a normalized "
                                        "B-basis function!");
                    }
                }
            }
        }

        (*@\Green{// If required, one can differentiate individual ordinary or B-basis functions, e.g.,}@*)
        GLint    space_index           = 1;
        GLint    function_index        = _n[space_index];
        GLint    differentiation_order = 2;
        GLdouble parameter             = (_alpha[space_index]+_beta[space_index])/2.0;

        GLdouble d_ordinary = (*_space[space_index])(
            ECSpace::ORDINARY_BASIS, function_index, differentiation_order, parameter);
        GLdouble d_B_basis  = (*_space[space_index])(
            ECSpace::B_BASIS, function_index, differentiation_order, parameter);
        (*@\Green{// Do something nice and meaningful with the obtained derivatives d\_ordinary and d\_B\_basis\ldots}@*)
        (*@\Green{// For verification purposes, variables d\_order and d\_B\_basis should contain the numerical values of}@*)
        (*@\Green{// the constants $\left.\frac{\mathrm{d}^2}{\mathrm{d}u^2}\cos\left(3u\right)\right|_{u=\frac{\pi}{4}}=\frac{9\sqrt{2}}{2}\approx 6.36396$ and $\left. \frac{\mathrm{d}^2}{\mathrm{d}u^2}\left[c_{10,5}\sin^5\left(\frac{\beta - u}{2}\right) \sin^5\left(\frac{u}{2}\right)\right]\right|_{\beta =\frac{\pi}{2}, u=\frac{\pi}{4}} = $}@*)
        (*@\Green{// $2220-\frac{3145 \sqrt{2}}{2} \approx -3.85083$, respectively, where $c_{10,5}=2368\sqrt{2}$.}@*)

        (*@\Green{// For debugging purposes one can also list the \LaTeX{} expressions of the ordinary basis functions, e.g.,}@*)
        string expression;
        for (GLint i = 0; i < _space[space_index]->dimension(); i++)
        {
            if (_space[space_index]->LaTeXExpression(i, expression))
            {
                cout << expression << endl;
            }
        }

        glClearColor(1.0, 1.0, 1.0, 1.0);       (*@\Green{// set the background color}@*)
        glEnable(GL_DEPTH_TEST);                (*@\Green{// enable depth testing}@*)
        glEnable(GL_LINE_SMOOTH);               (*@\Green{// enable the anti-aliasing of line primitives}@*)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
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

    (*@\Green{// revert to the original projection-view-model matrix}@*)
    _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, _PVM.address());

    for (GLuint i = 0; i < _space.size(); i++)
    {
        if (_img_B_basis[i])
        {
            (*@\Green{// constant reference to the $i$th normalized B-basis function of the current EC space}@*)
            const RowMatrix<SP<GenericCurve3>::Default> &reference = *_img_B_basis[i];

            (*@\Green{// in order to render each system of B-basis functions at different positions, we will use local}@*)
            (*@\Green{// translation matrices, based on which we also temporarily update the projection-view-model}@*)
            (*@\Green{// matrix of the color shader program}@*)
            Translate        local_T(-(_alpha[i] + _beta[i]) / 2.0, 1.7 - 0.75 * i, 0.0);
            GLTransformation local_PVM = _PVM * local_T;
            _color_shader.setUniformMatrix4fv("PVM", 1, GL_FALSE, local_PVM.address());

            for (GLint j = 0; j < reference.columnCount(); j++)
            {
                (*@\Green{// each B-basis function will be rendered by using different colors}@*)
                _color_shader.setUniformColor(
                        "color", coldToHotColormap(j, 0, reference.columnCount() - 1));
                reference[j]->renderDerivatives(_color_shader, 0, GL_LINE_STRIP);

                (*@\Green{// uncomment the line below in order to see the tangent vectors of the basis functions}@*)
                (*@\Green{// reference[j]-$>$renderDerivatives(\_color\_shader, 1, GL\_LINES);}@*)
            }
        }
    }

    _color_shader.disable();

    (*@\Green{// Fig.\ \mref{fig:differentiation_and_rendering_of_B_bases} illustrates the image obtained at run-time.}@*)
}
