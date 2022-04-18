/***********************************************************************
* mexeig.c : C mex file 
*
* mex -O -largeArrayDims -lmwlapack -lmwblas mexeig.c
*
* [Vvec, d] = blkEigMex(Xvec, PSDcone, options); 
* options = 1 (default), compute both eigenvectors and eigenvalues;
*         = 0, compute only eigenvalues.
*
* Xvec: m times 0  ->  Vvec: m times 0,  d: empty
* Xvec: 0 times n  ->  Vvec: 0 times n,  d: empty
* Xvec: 0 times 0  ->  Vvec: empty,  d: empty
***********************************************************************/

#include <math.h>
#include "mex.h"
#include "lapack.h"
#include <matrix.h>
#include <string.h> /* needed for memcpy() */

#if !defined(_WIN32) && !defined(__APPLE__)
#define dsyevd dsyevd_
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )
{    
    double   *X, *V, *D, *work, *d_PSDcone;

    size_t   m, n, numblk, numel, numeig, maxnk, options, k, Vpt, Dpt;
    ptrdiff_t  nk, lwork, lwork2, info;
    ptrdiff_t  *work2;
    char     jobz = 'V';
    char     uplo = 'U';

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

    if (nrhs > 3) { mexErrMsgTxt("mexeig: requires at most 3 input arguments."); }
    if (nlhs > 2) { mexErrMsgTxt("mexeig: requires at most 2 output argument."); }

/* Check the proper input type */

    if (mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]))  {mexErrMsgTxt("mexeig: sparse matrix not allowed.");}
    if (nrhs == 3 && !mxIsScalar(prhs[2])) {mexErrMsgTxt("mexeig: options must be 0 or 1");}

/* Read inputs */

    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]); 
    X = mxGetPr(prhs[0]);

    d_PSDcone = mxGetPr(prhs[1]);
    numblk = mxGetNumberOfElements(prhs[1]);
    numeig = 0;
    maxnk = 0;
    numel = 0;
    for (k = 0; k < numblk; k++)
    {
        nk = (ptrdiff_t)d_PSDcone[k];
        numel += nk * nk;
        numeig += nk;
        if (maxnk < nk) { maxnk = nk; }
    }
    options = 1; 
    if (nrhs == 3) { options = (size_t) *mxGetPr(prhs[2]); } 
    if (options == 1) { jobz = 'V'; } else { jobz = 'N'; } 

/* Check the input */

    if (numel != m * n) { mexErrMsgTxt("mexeig: dimensions do not match"); }
    if (m * n == 0)
    {
        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); 
        plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
        return;
    }

/* Create return argument */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); 
    V = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(numeig, 1, mxREAL); 
    D = mxGetPr(plhs[1]); 

/* Do the computations in a subroutine */
    lwork  = 1 + 6 * maxnk + 2 * maxnk * maxnk;
    work   = (double *)mxCalloc(lwork, sizeof(double)); 
    lwork2 = 3 + 5 * maxnk; 
    work2  = (ptrdiff_t *)mxCalloc(lwork2, sizeof(ptrdiff_t)); 

    memcpy(V, X, numel * sizeof(double));

    Vpt = 0;
    Dpt = 0;
    for (k = 0; k < numblk; k++)
    {
        nk = (ptrdiff_t)d_PSDcone[k];
        if (nk == 0) { continue; }
        dsyevd(&jobz, &uplo, &nk, V + Vpt, &nk, D + Dpt, work, &lwork, work2, &lwork2, &info); 
        Vpt += nk * nk;
        Dpt += nk; 
    }
    mxFree(work);
    mxFree(work2);
    return;
 }
/**********************************************************/
