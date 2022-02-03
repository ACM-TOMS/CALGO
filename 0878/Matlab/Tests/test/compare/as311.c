//AS311  Exact likelihood of VARMA model
//  This is a C-gateway function to the Fortran 77 subroutine ELF1 of [1], which
//  calculates the exact log-likelihood value of a vector ARMA time series 
//  model.
//
//       [LOGELF,F1,F2,A,IFAULT,QQ] = AS311(W,PHI,THETA,QQ,MU,ATF,SIGMA2,XITOL)
//
//  The Fortran sorce code and a PDF-copy of the journal article [1] are
//  available at the authors web site [2]. The article gives a detailed
//  description of the parameters. Here is an excerpt:
//
//     W      an m×n matrix containing observed data (m = dimension of each
//            observation, n = number of observations)
//     PHI    m×(p·m) matrix with autoregressive coefficients; PHI = 
//            [PHI1 PHI2...PHIp]
//     THETA  m×(q·m) matrix with autoregressive coefficients; THETA = 
//            [THETA1 THETA2...THETAq]
//     QQ     m×m covariance matrix of the driving shock series. On entry only
//            lower triangle need be set, but on exit the upper triangle is a
//            copy of the lower one.
//     MU     m-vector of means of W(t) (subtracted from each column of W
//            at start of calculation). Use MU=[] to obtain mean = 0
//     ATF    Set to 1 to return the shock series in the m×n matrix A, or to 0 
//            to return the intermediate quantities eta (see [1]) in A
//     SIGMA2 Scalar multiplier that multiplies QQ.
//     XITOL  Scalar convergence tolerance; see [1]
//     LOGELF The calculated log-likeliood (scalar)
//     F1,F2  Scalars with intermediate results (see [1])
//     A      See ATF above
//     IFAULT A fault indicator; 0 on success, otherwise nonzero (see [1] for
//            details)
//
//  The time-series model is
//
//     w(t) = PHI1·w(t-1)+...+PHIp·w(t-p) + a(t)+THETA1·a(t-1)+...+THETAq·a(t-q)
//
//  where w(t) is the m-vector observed at time t and a(t) is the random shock
//  at time t (t-th columns of W and A). The shocks are multivariate normal with
//  mean zero and covariance QQ.
//
//  If "mexing" has been set up in matlab, the mex-file as311.mexw32 (or
//  as311.dll in Matlab-versions before 7.1) may be created with the commands:
//
//           mex -f fortran-mexopts-file -c jam197.f
//           mex -f C-mexopts-file as311.c jam197.obj
//
//  [1] Mauricio, J. A. (1997). Algorithm AS 311: The exact likelihood function
//      of a vector autoregressive moving average model. J. R. Statist. Soc. C:
//      Applied Statistics, 46 no. 1, 157–171.
//
//  [2] Mauricio, J. A. (1997). JAM197.FOR (source code for [1]). Web page
//      http://www.ucm.es/info/ecocuan/jam/jam197.for

#include "mex.h"
#include <string.h>
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
  if (nin != 8) mexErrMsgTxt("Wrong number of input arguments");
  if (nout != 6) mexErrMsgTxt("Wrong number of output arguments");
  int m, n, im, p, q, atf, ismu, g, b1, b2, b3, nws, niws, *iws, int_ifault;
  double *w, *phi, *theta, *mu, *qq, *a, *sigma2, *xitol, *logelf, *f1, *f2, 
         *ws, *ifault;
  extern void elf1_();
  
  // w:
  w = mxGetPr(in[0]);
  m = mxGetM(in[0]);
  n = mxGetN(in[0]);
  im = m;
  //mexPrintf("m=%d n=%d\n",m,n);
  
  // phi:
  p = mxGetN(in[1])/mxGetM(in[3]);
  if (mxIsEmpty(in[1])) {
    p = 0;
    phi = mxCalloc(m, sizeof(double));
  }
  else {
    p = mxGetN(in[1])/mxGetM(in[3]);
    phi = mxCalloc(p*m*(m+1),sizeof(double));
    memcpy(phi, mxGetPr(in[1]), p*m*m*sizeof(double));
  }
  
  // theta:
  if (mxIsEmpty(in[2])) {
    q = 0;
    theta = mxCalloc(m, sizeof(double));
  }
  else {
    q = mxGetN(in[2])/mxGetM(in[3]);
    theta = mxCalloc(q*m*(m+1),sizeof(double));
    memcpy(theta, mxGetPr(in[2]), q*m*m*sizeof(double));
  }
  //mexPrintf("p=%d q=%d\n",p,q);
  
  // qq:
  out[5] = mxDuplicateArray(in[3]);
  qq = mxGetPr(out[5]);
  
  // ismu & mu:
  if (mxIsEmpty(in[4])) {
    ismu = 0; mu = 0;
  }
  else {
    ismu = 1; mu = mxGetPr(in[4]);
  }
  
  //atf & a:
  atf = *mxGetPr(in[5]); // convert to integer
  out[3] = mxCreateDoubleMatrix(m,n,mxREAL);
  a = mxGetPr(out[3]);
  
  //simga2:
  sigma2 = mxGetPr(in[6]);
  
  //xitol:
  xitol = mxGetPr(in[7]);
  
  //logelf:
  out[0] = mxCreateDoubleMatrix(1,1,mxREAL); 
  logelf = mxGetPr(out[0]);
  
  //f1 & f2:
  out[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  f1 = mxGetPr(out[1]);
  out[2] = mxCreateDoubleMatrix(1,1,mxREAL);
  f2 = mxGetPr(out[2]);

  //ws:
  g = p>q ? p : q;
  b1 = p>0 ? m*(m+1)/2+m*m*(p-1) : 1;
  b2 = b1 > g*m ? b1 : g*m;
  b3 = n > q ? n : q;
  nws = m*m*(3+3*g*g+(p+q)*g+b3)+b1*b1+b2+m;
  //mexPrintf("nws=%d\n",nws);
  ws = mxCalloc(nws, sizeof(double));
  //mexPrintf("ws=%d\n",ws);
  
  //iws:
  niws = b1;
  //mexPrintf("niws=%d\n",niws);
  iws = mxCalloc(niws, sizeof(int));
  //mexPrintf("iws=%d\n",iws);
  
  //call fortran subroutine:
  
  elf1_(&m,&im,&p,&q,&n,w,phi,theta,qq,&ismu,mu,&atf,a,sigma2,xitol,logelf,f1,
        f2,ws,&nws,iws,&niws,&int_ifault);
  
  //ifault:
  out[4] = mxCreateDoubleMatrix(1,1,mxREAL);
  ifault = mxGetPr(out[4]);
  *ifault = int_ifault; //change ifault to double
}
