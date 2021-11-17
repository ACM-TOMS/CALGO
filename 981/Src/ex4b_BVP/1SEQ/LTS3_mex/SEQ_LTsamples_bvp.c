/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

            U"   = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
            U(0) = s/(s^2+9)^2

                   s*cos(3*L)-3*sin(3*L)   (s*L+1)*sin(3*L)+3*L*cos(3*L)
            U(L) = --------------------- + -----------------------------
                         (s^2+9)^2                    6*(s^2+9)

      The analytical solution is:
                     s*cos(3*x)-3*sin(3*x)   (s*x+1)*sin(3*x)+3*x*cos(3*x)
            U(x,s) = --------------------- + -----------------------------
                           (s^2+9)^2                   6*(s^2+9)

        THE ODE PROBLEM IS SOLVED BY MATLAB bvp5c.m
 **/
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "MEX_complexArray.h"


double complex *SEQ_LTsamples_bvp (unsigned int NXval, double X[], unsigned int NOPTS, double complex S[], double tol)
{
    unsigned int h, k; /* for loop indices */

    /* ALLOCATE INPUT/OUTPUT PARAMETERS TO LT_samples MATLAB FUNCTION */
    mxArray *plhs[2], *prhs[3];     /* pointers to left and right side MATLAB parameters */

    /* BUILD INPUT PARAMETERS TO LT_samples.m function */
    double  *vr, *vi;
    /* Xval = [X(1), X(2), ..., X(NXval)]: user defined */
    prhs[0] = mxCreateDoubleMatrix(1,NXval,mxREAL);      /* 1st param.: X */
    vr = mxGetPr(prhs[0]);
    for (k=0; k<NXval; k++)
        vr[k] = X[k];

    /* Compute sample points on Talbot's contour */
    prhs[1] = mxCreateDoubleMatrix(1, NOPTS, mxCOMPLEX); /* 2nd param.: S (row-wise array) */
    vr = mxGetPr(prhs[1]); /* pointer to the double array REAL(prhs[0]) */
    vi = mxGetPi(prhs[1]); /* pointer to the double array IMAG(prhs[0]) */
    for (k=0; k<NOPTS; k++)
    {   vr[k] = creal(S[k]);
        vi[k] = cimag(S[k]);
    }

    prhs[2] = mxCreateDoubleScalar(tol);                 /* 3rd param.: tol */

    /* CALL THE MATLAB FUNCTION LT_samples:
                    [x, U] = LT_samples(Xval, S, tol)
       BY MEANS OF mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName) */
    mexCallMATLAB(2, plhs, 3, prhs, "LT_samples");

    /* THE 1st OUTPUT PARAMETER x OF LT_samples() IS NOT EXTRACTED
       Here the function is called with a user defined array Xspan and not
       ODE solver defined. */


    /* EXTRACT THE 2nd OUTPUT PARAMETER OF LT_samples():
        it is the solution of ODE problems and the output array */
    double complex *FS;
    if ( mxIsComplex(plhs[1]) )
        FS = cplxArray_fromMATLABtoC (NXval*NOPTS, 1, (const mxArray **)plhs);
    else
        FS = realArray_fromMATLABtoC (NXval*NOPTS, 1, (const mxArray **)plhs);

    return FS;
}

