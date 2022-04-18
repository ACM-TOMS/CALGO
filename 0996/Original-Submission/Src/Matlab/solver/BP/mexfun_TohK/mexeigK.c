/***********************************************************************
* mexeig.c : C mex file 
*
* mex -O -largeArrayDims -lmwlapack -lmwblas mexeig.c
*
* [V,D] = mexeig(A,options); 
* options = 1 (default), compute both eigenvectors and eigenvalues;
*         = 0, compute only eigenvalues.
***********************************************************************/

#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <string.h> /* needed for memcpy() */

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

#if !defined(_WIN32)
#define dsyevd dsyevd_
#define dgesdd dgesdd_
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{    double   *A, *V, *D, *work, *work2;  

     mwIndex  subs[2];
     mwSize   nsubs=2;
     mwIndex  *irD, *jcD;  
     mwSize   m, n, lwork, lwork2, info, j, k, jn, options; 
     char     *jobz="V";
     char     *uplo="U"; 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 2){
      mexErrMsgTxt("mexeig: requires at most 2 input arguments."); }
   if (nlhs > 2){ 
      mexErrMsgTxt("mexeig: requires at most 2 output argument."); }   

/* CHECK THE DIMENSIONS */

    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]); 
    if (m != n) { 
       mexErrMsgTxt("mexeig: matrix must be square."); }
    if (mxIsSparse(prhs[0])) {
       mexErrMsgTxt("mexeig: sparse matrix not allowed."); }   
    A = mxGetPr(prhs[0]);     
    options = 1; 
    if (nrhs==2) { options = (int)*mxGetPr(prhs[1]); } 
    if (options==1) { jobz="V"; } else { jobz="N"; } 

    /***** create return argument *****/
    plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL); 
    V = mxGetPr(plhs[0]);  
    plhs[1] = mxCreateSparse(n,n,n,mxREAL); 
    D   = mxGetPr(plhs[1]); 
    irD = mxGetIr(plhs[1]); 
    jcD = mxGetJc(plhs[1]);

    /***** Do the computations in a subroutine *****/
    lwork  = 1+6*n+2*n*n;  
    work   = mxCalloc(lwork,sizeof(double)); 
    lwork2 = 3 + 5*n; 
    work2  = mxCalloc(lwork2,sizeof(double)); 

    memcpy(mxGetPr(plhs[0]),mxGetPr(prhs[0]),(m*n)*sizeof(double));
    dsyevd(jobz,uplo,&n, V,&n, D, work,&lwork, work2,&lwork2, &info); 

    for (k=0; k<n; k++) { irD[k] = k; }
    jcD[0] = 0;
    for (k=1; k<=n; k++) { jcD[k] = k; }  
    return;
 }
/**********************************************************/
