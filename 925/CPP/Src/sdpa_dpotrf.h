/*-----------------------------------------------
  sdpa_dpotrf.cpp
  modification of ATL_dpotrfL
  int rATL_dpotrfL(int N, double *A,int lda)
  $Id: rsdpa_dpotrf.h,v 1.2 2004/09/01 06:34:12 makoto Exp $
-----------------------------------------------*/

#ifndef __sdpa_dpotrf_h__
#define __sdpa_dpotrf_h__

#ifdef __cplusplus
namespace sdpa {
  extern "C" int rATL_dpotrfL(int N, double *A,int lda);
  extern "C" void rdpotrfl_(int* N, double *A,int* lda,int* info);
} // end of namespace 'sdpa'

#else
  extern int rATL_dpotrfL(int N, double *A,int lda);
  void rdpotrfl_(int* N, double *A,int* lda,int* info);
#endif

#endif // __sdpa_dpotrf_h__


