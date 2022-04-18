/**   MEX_complexArray.h    -    AUXILIARY HEADER

    header file for MEX_complexArray.c

    >>>>>>>>>>>>>        VERSION 2.0     Jul 13th, 2013       <<<<<<<<<<
                            by Mariarosaria Rizzardi
 **/


#include "mex.h"
#include <complex.h>


/** Convert a complex array from MATLAB to C **/
double complex *cplxArray_fromMATLABtoC (int Nv, mwIndex Jprhs, const mxArray *prhs[]);

/** Convert a real array from MATLAB to a C complex one **/
double complex *realArray_fromMATLABtoC (int Nv, mwIndex Jprhs, const mxArray *prhs[]);
