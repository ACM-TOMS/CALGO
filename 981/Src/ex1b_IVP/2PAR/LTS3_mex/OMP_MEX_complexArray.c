/**   OMP_MEX_complexArray.c    -    AUXILIARY FUNCTIONS

    Convert MATLAB arrays into C double complex arrays

    >>>>>>>>>>>>>        VERSION 4.0     Jul 5th, 2015       <<<<<<<<<<
                            by Mariarosaria Rizzardi
 **/

#include <stdlib.h> /* for NULL */
#include <stdlib.h>
#include <complex.h>

#include "mex.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


double complex* OMP_cplxArray_fromMATLABtoC (int Nv, mwIndex Jprhs, const mxArray *prhs[], int THREADS)
/** Convert a complex array from MATLAB to C **/
{
    double  *vr, *vi;
    vr = mxGetPr(prhs[Jprhs]);
    vi = mxGetPi(prhs[Jprhs]);

    double complex  *v;
    v = (double complex *)malloc( Nv*sizeof(double complex) );
    if ( v == NULL )
        mexErrMsgTxt("\n***   ERROR IN OMP_cplxArray_fromMATLABtoC: DYNAMIC ALLOCATION OF v IS FAILED. ***\n");

    int  k;
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
    for (k=0; k<Nv; k++)
        *(v+k) = *(vr+k) + I* *(vi+k);

    return v;
}

double complex* OMP_realArray_fromMATLABtoC (int Nv, mwIndex Jprhs, const mxArray *prhs[], int THREADS)
/** Convert a real array from MATLAB to a C complex one **/
{
    double  *vr;
    vr = mxGetPr(prhs[Jprhs]);

    double complex  *v;
    v = (double complex *)malloc( Nv*sizeof(double complex) );
    if ( v == NULL )
        mexErrMsgTxt("\n***   ERROR IN OMP_realArray_fromMATLABtoC: DYNAMIC ALLOCATION OF v IS FAILED. ***\n");

    int  k;
    #pragma omp parallel for    default   (shared)    \
                                private   (k)         \
                                num_threads (THREADS)
    for (k=0; k<Nv; k++)
        *(v+k) = *(vr+k) + I*0.0;

    return v;
}

