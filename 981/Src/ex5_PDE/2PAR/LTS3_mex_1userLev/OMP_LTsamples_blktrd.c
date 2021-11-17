/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

                U_xx + U_yy - s*U = - [x*(x-1) + y*(y-1)],      0 < x,y < 1
                U(0,y) = U(1,y) = 4/s^2 + y*(y-1)/s
                U(x,0) = U(x,1) = 4/s^2 + x*(x-1)/s

        THE ODE PROBLEM IS SOLVED BY MATLAB backslash operator.

        OPEN MP-BASED IMPLEMENTATION.
 **/

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "OMP_MEX_complexArray.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


double complex *OMP_LTsamples_blktrd (unsigned int Nrows, double XY[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix FS of LT samples computed at (Xval[h], S[k])
    by solving ODE problems.
    FS has to be used row by row in the summation step
  *************************************************************************
    XY is a (col-wise) matrix of Nrows rows and two columns, where Nrows = (NXYval-2)^2,
    containing the cartesian coordinates of the internal mesh points.
    They have to be passed to the MATLAB function LT_samples().
**************************************************************************/
{
    unsigned int h, k; /* for loop indices */

    /* Nrows = (NXYval-2)*(NXYval-2) */
    unsigned int NXYval = (unsigned int)sqrt((double)Nrows) + 2;


    /* ALLOCATE INPUT/OUTPUT PARAMETERS TO LT_samples MATLAB FUNCTION */
    mxArray *plhs[1], *prhs[6];     /* pointers to left and right side MATLAB parameters */
    double  *vr, *vi;

    /* BUILD INPUT PARAMETERS TO LT_samples.m function */
    prhs[0] = mxCreateDoubleScalar((double)NXYval);      /* 1st param.: NXYval */
    prhs[1] = mxCreateDoubleMatrix(Nrows, 2, mxREAL);    /* 2nd param.: XY (col-wise matrix) */
    vr = mxGetPr(prhs[1]);
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
    for (k=0; k<2*Nrows; k++)
        vr[k] = XY[k];
    prhs[2] = mxCreateDoubleScalar((double)NOPTS);       /* 3rd param.: NOPTS */
    prhs[3] = mxCreateDoubleMatrix(1, NOPTS, mxCOMPLEX); /* 4th param.: S (row-wise array) */
    vr = mxGetPr(prhs[3]);
    vi = mxGetPi(prhs[3]);
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
    for (k=0; k<NOPTS; k++)
    {   vr[k] = creal(S[k]);
        vi[k] = cimag(S[k]);
    }
    prhs[4] = mxCreateDoubleScalar(tol);                 /* 6th param.: tol */
    prhs[5] = mxCreateDoubleScalar((double)THREADS);     /* 7th param.: thrds */

    /*  CALL THE MATLAB FUNCTION LT_samples:  U = LTsamples (NXYval,XY,NOPTS,S,tol,thrds)
        BY MEANS OF:
        mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName) */
        mexCallMATLAB(  1,  plhs,   6,  prhs, "LT_samples");

    /*  EXTRACT THE 1st OUTPUT PARAMETER FS OF LT_samples():
        it is the solution of ODE problems and the output array.
        FS is a row-wise matrix of size (Nrows, NOPTS) */
    double complex *FS;
    if ( mxIsComplex(plhs[0]) )
        FS = OMP_cplxArray_fromMATLABtoC (Nrows*NOPTS, 0, (const mxArray **)plhs,THREADS);
    else
        FS = OMP_realArray_fromMATLABtoC (Nrows*NOPTS, 0, (const mxArray **)plhs,THREADS);

    return FS;
}

