/**********************************************************************
* [y, B] = mexProxMonotonicNew(b, w, I)
* Use the pool-adjacent-violators algorithm
* to compute
* y = argmin   sum(w.*(y-b).^2)
*     sub. to  y(0+1) >= y(0+2) >= .... >= y(I[1]),
*              y(I[1]+1) >= ....   >= y(I[2]),
*              ...
*              y(I[m-1]+1) >= .... >= y(I[m]).
* w = positive wight
* B(k,1) = k-th distinct value in y
* B(k,2) = total weight for the k-th distinct value
* 
* mex -O  -largeArrayDims  mexProxMonotonicNew.c
*********************************************************************/
#include <mex.h>
#include <math.h>
#include <matrix.h> 
#include <stdlib.h>

#ifndef MWSIZE_MAX
    #define  mwIndex        size_t
    #define  mwSignedIndex  int
    #define  mwSize         size_t
#endif

/*************************************************************
*   PROCEDURE mexFunction - Entry for Matlab
**************************************************************/
 void mexFunction(const int nlhs, mxArray *plhs[],
                  const int nrhs, const mxArray *prhs[])
{
   mwIndex *irb, *jcb, *I; 
   mwSize   n, m, j, k, kend, cnt, precnt, *index, tmpn;
   double  *b, *y, *w, *btmp, *summ, *B, *Itmp;
   double  tmpsumm; 

   if (nrhs > 3) {
      mexErrMsgTxt("mexProxMonotonic requires at most 3 input argument."); }
   if (nlhs > 2) {
      mexErrMsgTxt("mexProxMonotonic generates at most 2 output argument."); }
 
   n = mxGetM(prhs[0]); 
   if (mxGetN(prhs[0])!= 1) {
      mexErrMsgTxt("mexProxMonotonic: b should be a column vector."); }   
   if (mxIsSparse(prhs[0])) {
      btmp = mxGetPr(prhs[0]);
      irb = mxGetIr(prhs[0]); jcb = mxGetJc(prhs[0]); 
      b = (double *)mxCalloc(n, sizeof(double));
      kend = jcb[1]; 
      for (k=0; k<kend; k++) { b[irb[k]] = btmp[k]; } 
   } else {
      b = mxGetPr(prhs[0]); 
   }
   if (nrhs >= 2) {
      if (mxGetM(prhs[1]) != n) {
          mexErrMsgTxt("mexProxMonotonic: size of b and w mismatch."); }
      if (mxIsSparse(prhs[1])) {
         mexErrMsgTxt("mexProxMonotonic: w cannot be sparse."); }
      w = mxGetPr(prhs[1]);
   } else {
      w = (double *)mxCalloc(n, sizeof(double));
      for (j=0; j<n; j++) { w[j]=1; }
   }
   if (nrhs == 3) {
      if (mxIsSparse(prhs[2])) {
         mexErrMsgTxt("mexProxMonotonic: I cannot be sparse."); }
      m = mxGetM(prhs[2]) - 1;
      Itmp = mxGetPr(prhs[2]);
      if (Itmp[m] != n) { 
          mexErrMsgTxt("mexProxMonotonic: Index of I mismatch."); }
      I = (mwIndex *)mxCalloc(m+1, sizeof(mwIndex));
      for (j=0; j<=m; j++) { I[j] = (mwIndex)Itmp[j]; }
   } else  {
      m = 1;
      I = (mwIndex *)mxCalloc(2, sizeof(mwIndex));
      I[0] = 0; I[1] = n;
   }
   /************************************************/
   plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);  
   y = mxGetPr(plhs[0]);  
   summ = (double *)mxCalloc(n, sizeof(double));
   index = (mwSize *)mxCalloc(n, sizeof(mwSize));   
   /************************************************/
   
   cnt = 0; 
   precnt = 0;
   for (k = 0; k < m; k++){
       for (j=I[k]; j<I[k+1]; j++) {
           index[cnt] = j;
           y[cnt] = b[j];
           summ[cnt] = w[j];
           while ((cnt > precnt) && (y[cnt-1] <= y[cnt])) {
               tmpsumm = summ[cnt-1]+summ[cnt];
               y[cnt-1] = (summ[cnt-1]*y[cnt-1]+summ[cnt]*y[cnt])/tmpsumm;
               summ[cnt-1] = tmpsumm;
               cnt--;
           }
           cnt++;
       }
       precnt = cnt;
   }
   cnt--;
   
   if (nlhs==2) {
       plhs[1] = mxCreateDoubleMatrix(cnt+1, 2, mxREAL);
       B = mxGetPr(plhs[1]);
       for (j=0; j<=cnt; j++) {
           B[j] = y[j]; 
           B[j+cnt+1] = summ[j];
       }
   }
   
   tmpn = n;
   for (; ; cnt--) { 
       for (j = index[cnt]; j < tmpn; j++) {
           y[j] = y[cnt];
       }
       tmpn = index[cnt];
       if (cnt==0) { break; }
   }

   if (mxIsSparse(prhs[0])) { mxFree(b); }
   if (nrhs < 2) { mxFree(w); }
   if (nrhs < 3) { mxFree(I); }
   mxFree(summ);
   mxFree(index);
   return;   
}
/*************************************************************/
