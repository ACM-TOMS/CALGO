/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

               U' = exp(-x)/(s^2+1),        x>0
               U(0) = 0

        THE ODE PROBLEM IS SOLVED BY MATLAB ode45.m + PCT.

        OPENMP-BASED PARALLEL VERSION.
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
    mxArray *plhs[2], *prhs[4]; /* pointers to left and right MATLAB params */


    /*  BUILD INPUT PARAMETERS TO LT_samples.m  */

    /*  1st param.: X
                    Xval = [X(1), X(2), ..., X(NXval)]: user defined */
    double  *vr, *vi;
    prhs[0] = mxCreateDoubleMatrix(1,NXval,mxREAL);
    vr = mxGetPr(prhs[0]); /* pointer to the double array prhs[0] */
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
        for (k=0; k<NXval; k++)
            vr[k] = Xval[k];

    /*  2nd param.: S (row-wise array)
                    ARRAY OF POINTS ON TALBOT'S CONTOUR */
    prhs[1] = mxCreateDoubleMatrix(1, NOPTS, mxCOMPLEX);
    vr = mxGetPr(prhs[1]); /* pointer to the double array REAL(prhs[0]) */
    vi = mxGetPi(prhs[1]); /* pointer to the double array IMAG(prhs[0]) */
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
        for (k=0; k<NOPTS; k++)
        {   vr[k] = creal(S[k]);
            vi[k] = cimag(S[k]);
        }

    /* 3rd param.: tol */
    prhs[2] = mxCreateDoubleScalar(tol);

    /* 4th param.: THREADS */
    prhs[3] = mxCreateDoubleScalar((double)THREADS);


    /*  CALL THE MATLAB FUNCTION LT_samples:   [x, U] = LT_samples(Xval, S, tol)
        BY MEANS OF:
        mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName) */
        mexCallMATLAB(  2,  plhs,   4,  prhs, "LT_samples");


    /*  EXTRACT THE OUTPUT PARAMETERS FROM LT_samples()

        THE 1st OUTPUT PARAMETER x IS NOT EXTRACTED
        Here the function is called with a user defined array x=Xval and not
        with a ODE solver defined array. */

    /*  EXTRACT THE 2nd OUTPUT PARAMETER FS OF LT_samples():
        it is the solution of ODE problems and the output array.
        FS is a row-wise matrix of size (NXval,NOPTS) */
    double complex *FS;
    if ( mxIsComplex(plhs[1]) )
        FS = OMP_cplxArray_fromMATLABtoC (NXval*NOPTS, 1, (const mxArray **)plhs,THREADS);
    else
        FS = OMP_realArray_fromMATLABtoC (NXval*NOPTS, 1, (const mxArray **)plhs,THREADS);

    return FS;
}

