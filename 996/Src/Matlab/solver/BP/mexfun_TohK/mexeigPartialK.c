/***********************************************************************
* mexeigPartial.c : C mex file for the routine DSYEVX in LAPACK
*
* mex -O -largeArrayDims -lmwlapack -lmwblas mexeigPartial.c
*
*  [V,D,nv] = mexeigPartial(A, vl, vu, options);
*
* options.jobz = 1 (default), compute both eigenvectors and eigenvalues;
*              = 2, compute only eigenvalues.
* options.range = 1 (default), compute all;
*               = 2, compute  eigenvalues between [vl, vu].
*               = 3, compute  vl-th to vu-th eigenvalues.
* options.abstol = tolerance (default = 1e-8)
***********************************************************************/

#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <blas.h>
#include <lapack.h>
#include <string.h> /* needed for memcpy() */

#if !defined(_WIN32)
#define dsyevd dsyevd_
#define dsyevx dsyevx_
#define dgesdd dgesdd_
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{    double   *A, *AA, *nV, *V, *D, VL, VU, abstol, *work, *work2, *tmp;  

     ptrdiff_t   m, n, IL, IU, lwork, *iwork, *ifail, info, options; 
     char     *jobz="V";
     char     *uplo="U"; 
     char     *range="A";

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 4){
      mexErrMsgTxt("mexeigPartial: requires at most 4 input arguments."); }
   if (nlhs > 3){ 
      mexErrMsgTxt("mexeigPartial: requires at most 3 output argument."); }   

/* CHECK THE DIMENSIONS */

    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]); 
    if (m != n) { 
       mexErrMsgTxt("mexeigPartial: matrix must be square."); }
    if (mxIsSparse(prhs[0])) {
       mexErrMsgTxt("mexeigPartial: sparse matrix not allowed."); }   
    A = mxGetPr(prhs[0]);     

    options = 1;    
    abstol = 1e-8;
    /*Get fields from options */
    tmp = mxGetPr(mxGetField(prhs[3], 0, "abstol"));  abstol  = *tmp;
    tmp = mxGetPr(mxGetField(prhs[3], 0, "jobz"));    options = (ptrdiff_t)*tmp;
    tmp = mxGetPr(mxGetField(prhs[3], 0, "range"));   options = (ptrdiff_t)*tmp;
    if (options==1) { 
       range="A";  } 
    else if (options==2) { 
       range="V"; VL = *mxGetPr(prhs[1]); VU = *mxGetPr(prhs[2]);  }
    else if (options == 3) { 
       range="I"; 
       IL = (ptrdiff_t)*mxGetPr(prhs[1]); 
       IU = (ptrdiff_t)*mxGetPr(prhs[2]);
    }   
    /***** create return argument *****/
    plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL); 
    V = mxGetPr(plhs[0]);  
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    D = mxGetPr(plhs[1]);  
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); 
    nV = mxGetPr(plhs[2]);

    /***** Do the computations in a subroutine *****/
    lwork  = 8*n;  
    work   =  mxCalloc(lwork,sizeof(double)); 
    iwork  =  mxCalloc(5*n, sizeof(ptrdiff_t));
    ifail  =  mxCalloc(n, sizeof(ptrdiff_t));  
    
    AA =  mxCalloc(n*n,sizeof(double)); 
    memcpy(AA,mxGetPr(prhs[0]),(n*n)*sizeof(double));
        
    dsyevx(jobz, range, uplo, &n, AA, &n, &VL, &VU, &IL, &IU,
           &abstol, &m, D, V, &n, work, &lwork, iwork, ifail, &info);
    
    *nV = m;    
    return;
 }
/**********************************************************/
