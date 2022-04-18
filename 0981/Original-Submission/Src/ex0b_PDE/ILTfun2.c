/** Example 2.2.6 (pag.100)
    Dean G. Duffy - Transform Methods for solving partial differential equations.
    Chapman & Hall/CRC, 2004
 */

#include <math.h>
#include "mex.h"

double ILTfun2(double eta, double t)
/* Inverse Laplace Transform function u(r,t) */
{
    double ft;

    /* u(eta,0+) = 0 */
    if ( fabs(t) < 1.0e-8 )
        return 0.0;

    /* u(0,0+) = f(t) = t*sin(3*t)/6 */
    if ( fabs(eta) < 1.0e-8 )
        return t*sin(3*t)/6;

    /* ALLOCATE INPUT/OUTPUT PARAMETERS TO ILT_fun MATLAB FUNCTION */
    mxArray *plhs[1], *prhs[2]; /* pointers to left and right side MATLAB parameters */
    prhs[0] = mxCreateDoubleScalar(eta); /* 1st par */
    prhs[1] = mxCreateDoubleScalar(t);   /* 2nd par */

    /* CALL THE MATLAB FUNCTION ILT_fun:
                    ft = ILT_fun(eta,t)
       BY MEANS OF mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName) */
    mexCallMATLAB(1, plhs, 2, prhs, "ILT_fun");
    ft = mxGetScalar( (const mxArray *) plhs[0] );

    return ft;
}

