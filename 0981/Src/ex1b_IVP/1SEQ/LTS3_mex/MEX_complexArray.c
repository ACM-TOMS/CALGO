/**   MEX_complexArray.c    -    AUXILIARY FUNCTIONS

    Convert MATLAB arrays into C double complex arrays

    >>>>>>>>>>>>>        VERSION 2.0     Jul 13th, 2013       <<<<<<<<<<
                            by Mariarosaria Rizzardi
 **/

#include "mex.h"
#include <stdlib.h>
#include <complex.h>


double complex* cplxArray_fromMATLABtoC (int Nv, mwIndex Jprhs, const mxArray *prhs[])
/** Convert a complex array from MATLAB to C **/
{
    double  *vr, *vi;
    vr = mxGetPr(prhs[Jprhs]);
    vi = mxGetPi(prhs[Jprhs]);

    double complex  *v;
    v = (double complex *)malloc(Nv*sizeof(double complex));
    /* v = (double complex *)mxMalloc(Nv*sizeof(double complex)); */
    if ( v == NULL )
        mexErrMsgTxt("\n***   ERROR IN cplxArray_fromMATLABtoC: DYNAMIC ALLOCATION OF v IS FAILED. ***\n");

    int  k;
    for (k=0; k<Nv; k++)
        *(v+k) = *(vr+k) + I* *(vi+k);

    return v;
}

double complex* realArray_fromMATLABtoC (int Nv, mwIndex Jprhs, const mxArray *prhs[])
/** Convert a real array from MATLAB to a C complex one **/
{
    double  *vr;
    vr = mxGetPr(prhs[Jprhs]);

    double complex  *v;
    v = (double complex *)malloc(Nv*sizeof(double complex));
    /* v = (double complex *)mxMalloc(Nv*sizeof(double complex)); */
    if ( v == NULL )
        mexErrMsgTxt("\n***   ERROR IN cplxArray_fromMATLABtoC: DYNAMIC ALLOCATION OF v IS FAILED. ***\n");

    int  k;
    for (k=0; k<Nv; k++)
        *(v+k) = *(vr+k) + I*0.0;

    return v;
}

