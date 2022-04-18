#include "mex.h"

int gcd(int x, int y)
{
    if (y == 0)
        return (x);
    else
        return (gcd(y, x % y));
}

void mexFunction(
         int nlhs,       mxArray *plhs[],
         int nrhs, const mxArray *prhs[]
         )
{
    /* Check for proper number of input and output arguments */
    if (nrhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:colgcd:invalidNumInputs",
                "One input argument required.");
    }
    if (nlhs > 2){
        mexErrMsgIdAndTxt( "MATLAB:colgcd:maxlhs",
                "Too many output arguments.");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0]))){
        mexErrMsgIdAndTxt( "MATLAB:colgcd:inputNotDouble",
                "Input argument must be of type double.");
    }
    if (!(mxIsSparse(prhs[0]))){
        mexErrMsgIdAndTxt( "MATLAB:colgcd:inputNotSparse",
                "Input argument must be of format sparse.");
    }
    if (mxGetNumberOfDimensions(prhs[0]) != 2){
        mexErrMsgIdAndTxt( "MATLAB:colgcd:inputNot2D",
                "Input argument must be two dimensional\n");
    }
    
    /* Get the size and pointers to input data */
    mwSize m  = mxGetM(prhs[0]);
    mwSize n  = mxGetN(prhs[0]);
    double *pr = mxGetPr(prhs[0]);
    //mwIndex *irs = mxGetIr(prhs[0]);
    mwIndex *jcs = mxGetJc(prhs[0]);

    double *pi = mxGetPi(prhs[0]);
    if (pi != NULL){
        mexErrMsgIdAndTxt( "MATLAB:colgcd:inputNotReal",
                "Input argument must be real.");
    }
    
    /* Allocate space for output vector */
    plhs[0] = mxCreateNumericMatrix(1, n, mxDOUBLE_CLASS, mxREAL);
    double *pl = mxGetPr(plhs[0]);
    
    /*  */
    for (mwSize j=0; j<n; j++) {
        mwSize top = jcs[j];
        mwSize end = jcs[j+1];
        int val = 0;
        for(mwSize i=top; i<end; i++){
            val = gcd((int)pr[i],val);
        }
        if (val==0){
            pl[j] = 1;
        }else{
            pl[j] = (double)val;
        }
    }
}