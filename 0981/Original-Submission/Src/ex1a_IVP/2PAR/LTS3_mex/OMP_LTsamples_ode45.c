/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

               U'    = s*U - x
               U(X0) = X0/s + 1/s^2

        THE ODE PROBLEM IS SOLVED BY ode45.m.

        OPENMP-BASED PARALLEL VERSION, WITHOUT MATLAB PARALLEL COMPUTING TOOLBOX.
 **/

#include <stdlib.h>
#include <complex.h>

#include "OMP_MEX_complexArray.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


double complex *OMP_LTsamples_ode45 (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix FS of LT samples computed at (Xval[h], S[k])
    by solving ODE problems.
    FS has to be used row by row in the summation step
**************************************************************************/
{
    unsigned int h, k; /* for loop indices */

    /* ********************************************************************** */
    /*  ALLOCATE INPUT/OUTPUT PARAMETERS TO LT_samples.m FUNCTION  */
    mxArray *plhs[2], *prhs[3]; /* pointers to left and right MATLAB params */

    /*  BUILD INPUT PARAMETERS TO LT_samples.m  */
    double  *vr, *vi;
    /* Xval = [X(1), X(2), ..., X(NXval)]: user defined */
    prhs[0] = mxCreateDoubleMatrix(1,NXval,mxREAL);      /* 1st param.: X */
    vr = mxGetPr(prhs[0]); /* pointer to the double array prhs[0] */

    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
        for (k=0; k<NXval; k++)
            vr[k] = Xval[k];

    /* SAMPLE POINTS ON TALBOT'S CONTOUR */
    prhs[1] = mxCreateDoubleMatrix(1, NOPTS, mxCOMPLEX); /* 2nd param.: S
                                                            (row-wise array) */
    vr = mxGetPr(prhs[1]); /* pointer to the double array REAL(prhs[0]) */
    vi = mxGetPi(prhs[1]); /* pointer to the double array IMAG(prhs[0]) */

    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
        for (k=0; k<NOPTS; k++)
        {   vr[k] = creal(S[k]);
            vi[k] = cimag(S[k]);
        }

    prhs[2] = mxCreateDoubleScalar(tol);                 /* 3rd param.: tol */

    /* CALL THE MATLAB FUNCTION LT_samples:
                    [x, U] = LT_samples(Xval, S, tol)
        BY MEANS OF:
        mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName) */
        mexCallMATLAB(  2,  plhs,   3,  prhs, "LT_samples");

    /* THE 1st OUTPUT PARAMETER x OF LT_samples() IS NOT EXTRACTED
        Here the function is called with a user defined array Xval and not
        ODE solver defined. */


    /* EXTRACT THE 2nd OUTPUT PARAMETER OF LT_samples():
        it is the solution of ODE problems and the output array */
    double complex *FS;
    if ( mxIsComplex(plhs[1]) ) /* 2nd output parameter */
        FS = OMP_cplxArray_fromMATLABtoC (NXval*NOPTS, 1, (const mxArray **)plhs, THREADS);
    else
        FS = OMP_realArray_fromMATLABtoC (NXval*NOPTS, 1, (const mxArray **)plhs, THREADS);

    return FS;
}

