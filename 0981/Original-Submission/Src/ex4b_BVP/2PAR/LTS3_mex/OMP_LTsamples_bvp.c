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

        THE ODE PROBLEM IS SOLVED BY MATLAB bvp5c.m + PCT
 **/

#include <stdlib.h>
#include <complex.h>

#include "OMP_MEX_complexArray.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


double complex *OMP_LTsamples_bvp (unsigned int NXval, double X[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix FS of LT samples computed at (Xval[h], S[k])
    by solving ODE problems.
    FS has to be used row by row in the summation step
**************************************************************************/
{
    unsigned int h, k; /* for loop indices */

    /* ALLOCATE INPUT/OUTPUT PARAMETERS TO LT_samples MATLAB FUNCTION */
    mxArray *plhs[2], *prhs[4];     /* pointers to left and right side MATLAB parameters */


    /* BUILD INPUT PARAMETERS TO LT_samples.m function */
    double  *vr, *vi;
    /* Xval = [X(1), X(2), ..., X(NXval)]: user defined */
    prhs[0] = mxCreateDoubleMatrix(1,NXval,mxREAL);      /* 1st param.: X */
    vr = mxGetPr(prhs[0]); /* pointer to the double array prhs[0] */
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
    for (k=0; k<NXval; k++)
        vr[k] = X[k];

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


    /*  CALL THE MATLAB FUNCTION LT_samples:   [x, U] = LT_samples(Xval, S, tol, thrds)
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

