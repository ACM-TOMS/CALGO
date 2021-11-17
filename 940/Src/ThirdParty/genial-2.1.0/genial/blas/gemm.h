//GENIAL - GENeric Image & Array Library
//Copyright (C) 2006  Patrick LAURENT
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#ifndef GEMM_H
#define GEMM_H

#ifdef __cplusplus

#include "blas/copy.h"

#if defined(BLAS_THREADING)
#include "threads.h"
using namespace gmt;
#endif



//Group = Linear Algebra

#ifndef BLAS_M
#define BLAS_M  84
#endif
#ifndef BLAS_N
#define BLAS_N  84
#endif
#ifndef BLAS_K
#define BLAS_K  44
#endif
#ifndef BLAS_K2
#define BLAS_K2 24
#endif


template<template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void real_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
  F<value_type,value_type,value_type> f;
  Func<value_type,value_type> func;

  int m=A.nrows(), n=B.nrows(), k=A.ncols(); assert(C.nrows()==m && C.ncols()==n && B.ncols()==k);


  int i,j=0;
  for (i=0; i<m; ++i)
  {
    j=0;
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);
    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);
    for (int jmax=n-2; j<=jmax; j+=2)
    {
      value_type a,s0,s1;
      typename MatrixRow<const Matrix<G2> >::self B0=row(B,j), B1=row(B,j+1);
      {  a=Ai[0]; s0 =f(a,B0[0]); s1 =f(a,B1[0]); }
      for (int h=1; h<k; ++h) { a=Ai[h]; s0+=f(a,B0[h]); s1+=f(a,B1[h]); }
      adapt(j  ,func(s0));
      adapt(j+1,func(s1));
    }
  }
  if (j<n)
  for (int i=0; i<m; ++i)
  {
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);
    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);
    {
      value_type s0, a;
      typename MatrixRow<const Matrix<G2> >::self B0=row(B,j);
      { a=Ai[0]; s0 =f(a,B0[0]); }
      for (int h=1; h<k; ++h) { a=Ai[h]; s0+=f(a,B0[h]); }
      adapt(j  ,func(s0));
    }
  }
}

template<template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void complex_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
  F<value_type,value_type,value_type> f;
  Func<value_type,value_type> func;

  int m=A.nrows(), n=B.nrows(), k=A.ncols(), k2=k/2; assert(C.nrows()==m && C.ncols()==n && B.ncols()==k);

  for (int i=0; i<m; ++i)
  {
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);
    typename SubArray<typename MatrixRow<const Matrix<G1> >::self>::self A0=sub(Ai, 0,k2);
    typename SubArray<typename MatrixRow<const Matrix<G1> >::self>::self A1=sub(Ai,k2,k2);

    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);
    for (int j=0; j<n; ++j)
    {
      value_type s0,s1, a;
      typename MatrixRow<const Matrix<G2> >::self Bj=row(B,j);
      typename SubArray<typename MatrixRow<const Matrix<G2> >::self>::self B0=sub(Bj, 0,k2);
      typename SubArray<typename MatrixRow<const Matrix<G2> >::self>::self B1=sub(Bj,k2,k2);

      { a=A0[0]; s0 =f(a,B0[0]); s1 =f(a,B1[0]); }
      for (int h=1; h<k2; ++h) { a=A0[h]; s0+=f(a,B0[h]); s1+=f(a,B1[h]); }
      for (int h=0; h<k2; ++h) { a=A1[h]; s1+=f(a,B0[h]); s0-=f(a,B1[h]); }
      adapt(j,typename Matrix<G3>::value_type(func(s0),func(s1)));
    }
  }
}



#define MATMUL_0           { a=Ai[0]; s0 =f(a,Bj[0]); s1 =f(a,Bj[K  ]); s2 =f(a,Bj[2*K  ]); s3 =f(a,Bj[3*K  ]); s4 =f(a,Bj[4*K  ]); s5 =f(a,Bj[5*K  ]); }
#define MATMUL(k) if (k<K) { a=Ai[k]; s0+=f(a,Bj[k]); s1+=f(a,Bj[K+k]); s2+=f(a,Bj[2*K+k]); s3+=f(a,Bj[3*K+k]); s4+=f(a,Bj[4*K+k]); s5+=f(a,Bj[5*K+k]); }

template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void real_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
  F<value_type,value_type,value_type> f;
  Func<value_type,value_type> func;

  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);

  int i=0,j=0;
  for (; i<m; ++i)
  {
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);
    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);
    j=0;
    for (int jmax=n-6; j<=jmax; j+=6)
    {
      typename DataVector<const typename Matrix<G2>::value_type>::self Bj(6*K,&B(j,0));
      value_type s0,s1,s2,s3,s4,s5, a;
      MATMUL_0   MATMUL( 1) MATMUL( 2) MATMUL( 3) MATMUL( 4) MATMUL( 5) MATMUL( 6) MATMUL( 7) MATMUL( 8) MATMUL( 9)
      MATMUL(10) MATMUL(11) MATMUL(12) MATMUL(13) MATMUL(14) MATMUL(15) MATMUL(16) MATMUL(17) MATMUL(18) MATMUL(19)
      MATMUL(20) MATMUL(21) MATMUL(22) MATMUL(23) MATMUL(24)
      adapt(j+0,func(s0)); adapt(j+1,func(s1)); adapt(j+2,func(s2)); adapt(j+3,func(s3)); adapt(j+4,func(s4)); adapt(j+5,func(s5));
    }
  }
  if (j<n)
  {
    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
    real_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
  }
}

template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void sse_real_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
  F<value_type,value_type,value_type> f;
  Func<value_type,value_type> func;

  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);

  int i=0,j=0;
  for (; i<m; ++i)
  {
    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);
    j=0;
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);
    const typename Matrix<G2>::value_type *pb=&B(0,0);
    for (int jmax=n-6; j<=jmax; j+=6, pb+=6*K)
    {
      typename DataVector<const typename Matrix<G2>::value_type>::self Bj(6*K,pb);
      value_type s0,s1,s2,s3,s4,s5, a;
      MATMUL_0   MATMUL( 1) MATMUL( 2) MATMUL( 3) MATMUL( 4) MATMUL( 5) MATMUL( 6) MATMUL( 7) MATMUL( 8) MATMUL( 9)
      MATMUL(10) MATMUL(11) MATMUL(12) MATMUL(13) MATMUL(14) MATMUL(15) MATMUL(16) MATMUL(17) MATMUL(18) MATMUL(19)
      MATMUL(20) MATMUL(21) MATMUL(22) MATMUL(23) MATMUL(24)
      adapt(j  ,func(s0,s1,s2,s3));
      adapt(j+4,func(s4,s5));
    }
  }
  if (j<n)
  {
    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
    real_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
  }
}



//#define MATMUL2x4_0           { a0=Ai[0]; a1=Ai[K  ]; s0 =hadd(f(a0,Bj[0]),f(a1,Bj[0])); s1 =hadd(f(a0,Bj[K  ]),f(a1,Bj[K  ])); s2 =hadd(f(a0,Bj[2*K  ]),f(a1,Bj[2*K  ])); s3 =hadd(f(a0,Bj[3*K  ]),f(a1,Bj[3*K  ])); }
//#define MATMUL2x4(k) if (k<K) { a0=Ai[k]; a1=Ai[K+k]; s0+=hadd(f(a0,Bj[k]),f(a1,Bj[k])); s1+=hadd(f(a0,Bj[K+k]),f(a1,Bj[K+k])); s2+=hadd(f(a0,Bj[2*K+k]),f(a1,Bj[2*K+k])); s3+=hadd(f(a0,Bj[3*K+k]),f(a1,Bj[3*K+k])); }
//
//template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
//void sse3_float_real_2x4_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
//{
//  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
//  F<value_type,value_type,value_type> f;
//  Func<value_type,value_type> func;
//
//  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);
//
//  int i=0,j=0;
//  for (int imax=m-2; i<=imax; i+=2)
//  {
//    j=0;
//    typename DataVector<const typename Matrix<G1>::value_type>::self Ai(2*K,&A(i,0));
//    typename MatrixRow<Matrix<G3> >::self C0(C,i), C1(C,i+1);
//    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt0=adaptor(C0), adapt1=adaptor(C1);
//    for (int jmax=n-4; j<=jmax; j+=4)
//    {
//      typename DataVector<const typename Matrix<G2>::value_type>::self Bj(4*K,&B(j,0));
//      value_type s0,s1,s2,s3, a0,a1;
//      MATMUL2x4_0   MATMUL2x4( 1) MATMUL2x4( 2) MATMUL2x4( 3) MATMUL2x4( 4) MATMUL2x4( 5) MATMUL2x4( 6) MATMUL2x4( 7) MATMUL2x4( 8) MATMUL2x4( 9)
//      MATMUL2x4(10) MATMUL2x4(11) MATMUL2x4(12) MATMUL2x4(13) MATMUL2x4(14) MATMUL2x4(15) MATMUL2x4(16) MATMUL2x4(17) MATMUL2x4(18) MATMUL2x4(19)
//      MATMUL2x4(20) MATMUL2x4(21) MATMUL2x4(22) MATMUL2x4(23) MATMUL2x4(24)
//      s0=hadd(s0,s1);
//      s1=hadd(s2,s3);
//      adapt0(j,shuffle<_MM_SHUFFLE(2,0,2,0)>(s0,s1));
//      adapt1(j,shuffle<_MM_SHUFFLE(3,1,3,1)>(s0,s1));
//    }
//  }
//  if (j<n)
//  {
//    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
//    real_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
//  }
//}
//
//template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
//void sse3_double_real_2x4_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
//{
//  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
//  F<value_type,value_type,value_type> f;
//  Func<value_type,value_type> func;
//
//  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);
//
//  int i=0,j=0;
//  for (int imax=m-2; i<=imax; i+=2)
//  {
//    j=0;
//    typename DataVector<const typename Matrix<G1>::value_type>::self Ai(2*K,&A(i,0));
//    typename MatrixRow<Matrix<G3> >::self C0(C,i), C1(C,i+1);
//    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt0=adaptor(C0), adapt1=adaptor(C1);
//    for (int jmax=n-4; j<=jmax; j+=4)
//    {
//      typename DataVector<const typename Matrix<G2>::value_type>::self Bj(4*K,&B(j,0));
//      value_type s0,s1,s2,s3, a0,a1;
//      MATMUL2x4_0   MATMUL2x4( 1) MATMUL2x4( 2) MATMUL2x4( 3) MATMUL2x4( 4) MATMUL2x4( 5) MATMUL2x4( 6) MATMUL2x4( 7) MATMUL2x4( 8) MATMUL2x4( 9)
//      MATMUL2x4(10) MATMUL2x4(11) MATMUL2x4(12) MATMUL2x4(13) MATMUL2x4(14) MATMUL2x4(15) MATMUL2x4(16) MATMUL2x4(17) MATMUL2x4(18) MATMUL2x4(19)
//      MATMUL2x4(20) MATMUL2x4(21) MATMUL2x4(22) MATMUL2x4(23) MATMUL2x4(24)
//      adapt0(j,unpacklo(s0,s1)); adapt0(j+2,unpacklo(s2,s3));
//      adapt1(j,unpackhi(s0,s1)); adapt1(j+2,unpackhi(s2,s3));
//    }
//  }
//  if (j<n)
//  {
//    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
//    real_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
//  }
//}


//#define MATMUL2x2_0           { a0=Ai[0]; a1=Ai[  K]; s00 =f(a0,Bj[0]); s10 =f(a1,Bj[0]); s01 =f(a0,Bj[K  ]); s11 =f(a1,Bj[K  ]); }
//#define MATMUL2x2(k) if (k<K) { a0=Ai[k]; a1=Ai[K+k]; s00+=f(a0,Bj[k]); s10+=f(a1,Bj[k]); s01+=f(a0,Bj[K+k]); s11+=f(a1,Bj[K+k]); }
//
//template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
//void sse_real_2x2_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
//{
//  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
//  F<value_type,value_type,value_type> f;
//  Func<value_type,value_type> func;
//
//  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);
//
//  int i=0,j=0;
//  for (int imax=m-2; i<=imax; i+=2)
//  {
//    j=0;
//    typename DataVector<const typename Matrix<G1>::value_type>::self Ai(2*K,&A(i,0));
//    typename MatrixRow<Matrix<G3> >::self C0(C,i), C1(C,i+1);
//    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt0=adaptor(C0), adapt1=adaptor(C1);
//    for (int jmax=n-2; j<=jmax; j+=2)
//    {
//      typename DataVector<const typename Matrix<G2>::value_type>::self Bj(2*K,&B(j,0));
//      value_type s00,s01,s10,s11, a0,a1;
//      MATMUL2x2_0   MATMUL2x2( 1) MATMUL2x2( 2) MATMUL2x2( 3) MATMUL2x2( 4) MATMUL2x2( 5) MATMUL2x2( 6) MATMUL2x2( 7) MATMUL2x2( 8) MATMUL2x2( 9)
//      MATMUL2x2(10) MATMUL2x2(11) MATMUL2x2(12) MATMUL2x2(13) MATMUL2x2(14) MATMUL2x2(15) MATMUL2x2(16) MATMUL2x2(17) MATMUL2x2(18) MATMUL2x2(19)
//      MATMUL2x2(20) MATMUL2x2(21) MATMUL2x2(22) MATMUL2x2(23) MATMUL2x2(24)
//      adapt0(j,func(s00,s01));
//      adapt1(j,func(s10,s11));
//    }
//  }
//  if (j<n)
//  {
//    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
//    real_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
//  }
//}

#define MATMULR_0             { value_type a=Ai[    0]; s0 =f(a,Bj[0]); s1 =f(a,Bj[K/2  ]); s2 =f(a,Bj[2*K/2  ]); s3 =f(a,Bj[3*K/2  ]); s4 =f(a,Bj[4*K/2  ]); s5 =f(a,Bj[5*K/2  ]); }
#define MATMULR(k) if (k<K/2) { value_type a=Ai[    k]; s0+=f(a,Bj[k]); s1+=f(a,Bj[K/2+k]); s2+=f(a,Bj[2*K/2+k]); s3+=f(a,Bj[3*K/2+k]); s4+=f(a,Bj[4*K/2+k]); s5+=f(a,Bj[5*K/2+k]); }
#define MATMULI(k) if (k<K/2) { value_type a=Ai[K/2+k]; s1+=f(a,Bj[k]); s0-=f(a,Bj[K/2+k]); s3+=f(a,Bj[2*K/2+k]); s2-=f(a,Bj[3*K/2+k]); s5+=f(a,Bj[4*K/2+k]); s4-=f(a,Bj[5*K/2+k]); }

template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void complex_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
  F<value_type,value_type,value_type> f;
  Func<value_type,value_type> func;

  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);

  int i,j;
  for (i=0; i<m; ++i)
  {
    j=0;
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);
    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);

    for (int jmax=n-3; j<=jmax; j+=3)
    {
      value_type s0,s1,s2,s3,s4,s5;
      {
        typename DataVector<const typename Matrix<G2>::value_type>::self Bj(3*K,&B(j,0));
        MATMULR_0   MATMULR( 1) MATMULR( 2) MATMULR( 3) MATMULR( 4) MATMULR( 5) MATMULR( 6) MATMULR( 7) MATMULR( 8) MATMULR( 9)
        MATMULR(10) MATMULR(11) MATMULR(12) MATMULR(13) MATMULR(14) MATMULR(15) MATMULR(16) MATMULR(17) MATMULR(18) MATMULR(19)
        MATMULR(20) MATMULR(21) MATMULR(22) MATMULR(23) MATMULR(24) MATMULR(25) MATMULR(26) MATMULR(27) MATMULR(28) MATMULR(29)
        MATMULR(30) MATMULR(31) MATMULR(32) MATMULR(33) MATMULR(34) MATMULR(35) MATMULR(36) MATMULR(37) MATMULR(38) MATMULR(39)
        MATMULR(40) MATMULR(41) MATMULR(42) MATMULR(43) MATMULR(44) MATMULR(45) MATMULR(46) MATMULR(47) MATMULR(48) MATMULR(49)
      }
      {
        typename DataVector<const typename Matrix<G2>::value_type>::self Bj(3*K,&B(j,0));
        MATMULI( 0) MATMULI( 1) MATMULI( 2) MATMULI( 3) MATMULI( 4) MATMULI( 5) MATMULI( 6) MATMULI( 7) MATMULI( 8) MATMULI( 9)
        MATMULI(10) MATMULI(11) MATMULI(12) MATMULI(13) MATMULI(14) MATMULI(15) MATMULI(16) MATMULI(17) MATMULI(18) MATMULI(19)
        MATMULI(20) MATMULI(21) MATMULI(22) MATMULI(23) MATMULI(24) MATMULI(25) MATMULI(26) MATMULI(27) MATMULI(28) MATMULI(29)
        MATMULI(30) MATMULI(31) MATMULI(32) MATMULI(33) MATMULI(34) MATMULI(35) MATMULI(36) MATMULI(37) MATMULI(38) MATMULI(39)
        MATMULI(40) MATMULI(41) MATMULI(42) MATMULI(43) MATMULI(44) MATMULI(45) MATMULI(46) MATMULI(47) MATMULI(48) MATMULI(49)
      }
      adapt(j  ,typename Matrix<G3>::value_type(func(s0),func(s1)));
      adapt(j+1,typename Matrix<G3>::value_type(func(s2),func(s3)));
      adapt(j+2,typename Matrix<G3>::value_type(func(s4),func(s5)));
    }
  }
  if (j<n)
  {
    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
    complex_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
  }
}


template<int K,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void sse_complex_matmul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type;
  F<value_type,value_type,value_type> f;
  Func<value_type,value_type> func;

  int m=A.nrows(), n=B.nrows(); assert(C.nrows()==m && C.ncols()==n && A.ncols()==K && B.ncols()==K);

  int i=0,j=0;
  for (; i<m; ++i)
  {
    typename MatrixRow<      Matrix<G3> >::self Ci=row(C,i);
    typename Adapt::template result_rebind<typename MatrixRow<Matrix<G3> >::self>::other adapt=adaptor(Ci);
    j=0;
    typename MatrixRow<const Matrix<G1> >::self Ai=row(A,i);

    for (int jmax=n-3; j<=jmax; j+=3)
    {
      value_type s0,s1,s2,s3,s4,s5;
      {
        typename DataVector<const typename Matrix<G2>::value_type>::self Bj(3*K,&B(j,0));
        MATMULR_0   MATMULR( 1) MATMULR( 2) MATMULR( 3) MATMULR( 4) MATMULR( 5) MATMULR( 6) MATMULR( 7) MATMULR( 8) MATMULR( 9)
        MATMULR(10) MATMULR(11) MATMULR(12) MATMULR(13) MATMULR(14) MATMULR(15) MATMULR(16) MATMULR(17) MATMULR(18) MATMULR(19)
        MATMULR(20) MATMULR(21) MATMULR(22) MATMULR(23) MATMULR(24) MATMULR(25)
      }
      {
        typename DataVector<const typename Matrix<G2>::value_type>::self Bj(3*K,&B(j,0));
        MATMULI( 0) MATMULI( 1) MATMULI( 2) MATMULI( 3) MATMULI( 4) MATMULI( 5) MATMULI( 6) MATMULI( 7) MATMULI( 8) MATMULI( 9)
        MATMULI(10) MATMULI(11) MATMULI(12) MATMULI(13) MATMULI(14) MATMULI(15) MATMULI(16) MATMULI(17) MATMULI(18) MATMULI(19)
        MATMULI(20) MATMULI(21) MATMULI(22) MATMULI(23) MATMULI(24) MATMULI(25)
      }
      typename SimdVector<4,typename Matrix<G3>::value_type::value_type>::self a=func(s0,s1,s2,s3);
      typename SimdVector<2,typename Matrix<G3>::value_type::value_type>::self b=func(s4,s5);
      //Vector<tiny_vector_generator<4,typename Matrix<G3>::value_type::value_type> > a=func(s0,s1,s2,s3);
      //Vector<tiny_vector_generator<2,typename Matrix<G3>::value_type::value_type> > b=func(s4,s5);
      adapt(j  ,reinterpret_cast<typename SimdVector<2,typename Matrix<G3>::value_type>::self &>(a));
      adapt(j+2,reinterpret_cast<typename SimdVector<1,typename Matrix<G3>::value_type>::self &>(b));
    }
  }
  if (j<n)
  {
    typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,m,n-j);
    complex_matmul<F,Func>(sub(A,0,0,m,K),sub(B,j,0,n-j,K),C2,adaptor);
  }
}

template<int M,int N,int K,int L,class T,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function    { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return real_matmul   <K,F,Func>(X,Y,C,adapt); } };
template<int M,int N,int K,int L,class T,template<class,class,class> class F,template<class,class> class Func> struct complex_matmul_function { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return complex_matmul<K,F,Func>(X,Y,C,adapt); } };

template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function   <M,N,K,4,float ,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse_real_matmul    <K,F,Func>(X,Y,C,adapt); } };
template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function   <M,N,K,2,double,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse_real_matmul    <K,F,Func>(X,Y,C,adapt); } };

//template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function   <M,N,K,4,float ,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse3_float_real_2x4_matmul <K,F,Func>(X,Y,C,adapt); } };
//template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function   <M,N,K,2,double,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse3_double_real_2x4_matmul<K,F,Func>(X,Y,C,adapt); } };

//template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function   <M,N,K,4,float ,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse_real_2x2_matmul<K,F,Func>(X,Y,C,adapt); } };
//template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct real_matmul_function   <M,N,K,2,double,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse_real_2x2_matmul<K,F,Func>(X,Y,C,adapt); } };


template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct complex_matmul_function<M,N,K,4,float ,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse_complex_matmul <K,F,Func>(X,Y,C,adapt); } };
template<int M,int N,int K,template<class,class,class> class F,template<class,class> class Func> struct complex_matmul_function<M,N,K,2,double,F,Func>{ template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return sse_complex_matmul <K,F,Func>(X,Y,C,adapt); } };


template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_function { };
template<int M,int N,            int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_function<M,N,0,0,L,T,F,Func> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return real_matmul   <  F,Func>(X,Y,C,adapt); } };
template<int M,int N,            int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_function<M,N,0,1,L,T,F,Func> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return complex_matmul<  F,Func>(X,Y,C,adapt); } };
template<int M,int N,int K,      int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_function<M,N,K,0,L,T,F,Func> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return real_matmul_function   <M,N,K,L,T,F,Func>()(X,Y,C,adapt); } };
template<int M,int N,int K,      int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_function<M,N,K,1,L,T,F,Func> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return complex_matmul_function<M,N,K,L,T,F,Func>()(X,Y,C,adapt); } };



template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3>
void real_matmul_switch(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const id_adaptor_function &adapt)
{
  assert(C.nrows()==A.nrows() && C.ncols()==B.nrows() && B.ncols()==A.ncols());

  switch (A.ncols())
  {
    case  0: { return; }
    case  1: if ( 1<=K) return matmul_function<M,N, 1,R,L,T,F,Func>()(A,B,C,adapt);
    case  2: if ( 2<=K) return matmul_function<M,N, 2,R,L,T,F,Func>()(A,B,C,adapt);
    case  3: if ( 3<=K) return matmul_function<M,N, 3,R,L,T,F,Func>()(A,B,C,adapt);
    case  4: if ( 4<=K) return matmul_function<M,N, 4,R,L,T,F,Func>()(A,B,C,adapt);
    case  5: if ( 5<=K) return matmul_function<M,N, 5,R,L,T,F,Func>()(A,B,C,adapt);
    case  6: if ( 6<=K) return matmul_function<M,N, 6,R,L,T,F,Func>()(A,B,C,adapt);
    case  7: if ( 7<=K) return matmul_function<M,N, 7,R,L,T,F,Func>()(A,B,C,adapt);
    case  8: if ( 8<=K) return matmul_function<M,N, 8,R,L,T,F,Func>()(A,B,C,adapt);
    case  9: if ( 9<=K) return matmul_function<M,N, 9,R,L,T,F,Func>()(A,B,C,adapt);
    case 10: if (10<=K) return matmul_function<M,N,10,R,L,T,F,Func>()(A,B,C,adapt);
    case 11: if (11<=K) return matmul_function<M,N,11,R,L,T,F,Func>()(A,B,C,adapt);
    case 12: if (12<=K) return matmul_function<M,N,12,R,L,T,F,Func>()(A,B,C,adapt);
    case 13: if (13<=K) return matmul_function<M,N,13,R,L,T,F,Func>()(A,B,C,adapt);
    case 14: if (14<=K) return matmul_function<M,N,14,R,L,T,F,Func>()(A,B,C,adapt);
    case 15: if (15<=K) return matmul_function<M,N,15,R,L,T,F,Func>()(A,B,C,adapt);
    case 16: if (16<=K) return matmul_function<M,N,16,R,L,T,F,Func>()(A,B,C,adapt);
    case 17: if (17<=K) return matmul_function<M,N,17,R,L,T,F,Func>()(A,B,C,adapt);
    case 18: if (18<=K) return matmul_function<M,N,18,R,L,T,F,Func>()(A,B,C,adapt);
    case 19: if (19<=K) return matmul_function<M,N,19,R,L,T,F,Func>()(A,B,C,adapt);
    case 20: if (20<=K) return matmul_function<M,N,20,R,L,T,F,Func>()(A,B,C,adapt);
    case 21: if (21<=K) return matmul_function<M,N,21,R,L,T,F,Func>()(A,B,C,adapt);
    case 22: if (22<=K) return matmul_function<M,N,22,R,L,T,F,Func>()(A,B,C,adapt);
    case 23: if (23<=K) return matmul_function<M,N,23,R,L,T,F,Func>()(A,B,C,adapt);
    case 24: if (24<=K) return matmul_function<M,N,24,R,L,T,F,Func>()(A,B,C,adapt);
    case 25: if (25<=K) return matmul_function<M,N,25,R,L,T,F,Func>()(A,B,C,adapt);
    default: { return matmul_function<M,N,0,R,L,T,F,Func>()(A,B,C,adapt); }
  }
}

template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3>
void complex_matmul_switch(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const id_adaptor_function &adapt)
{
  assert(C.nrows()==A.nrows() && C.ncols()==B.nrows() && B.ncols()==A.ncols());

  switch (A.ncols())
  {
    case  0: { return; }
    case  2: if ( 2<=K) return matmul_function<M,N, 2,R,L,T,F,Func>()(A,B,C,adapt);
    case  4: if ( 4<=K) return matmul_function<M,N, 4,R,L,T,F,Func>()(A,B,C,adapt);
    case  6: if ( 6<=K) return matmul_function<M,N, 6,R,L,T,F,Func>()(A,B,C,adapt);
    case  8: if ( 8<=K) return matmul_function<M,N, 8,R,L,T,F,Func>()(A,B,C,adapt);
    case 10: if (10<=K) return matmul_function<M,N,10,R,L,T,F,Func>()(A,B,C,adapt);
    case 12: if (12<=K) return matmul_function<M,N,12,R,L,T,F,Func>()(A,B,C,adapt);
    case 14: if (14<=K) return matmul_function<M,N,14,R,L,T,F,Func>()(A,B,C,adapt);
    case 16: if (16<=K) return matmul_function<M,N,16,R,L,T,F,Func>()(A,B,C,adapt);
    case 18: if (18<=K) return matmul_function<M,N,18,R,L,T,F,Func>()(A,B,C,adapt);
    case 20: if (20<=K) return matmul_function<M,N,20,R,L,T,F,Func>()(A,B,C,adapt);
    case 22: if (22<=K) return matmul_function<M,N,22,R,L,T,F,Func>()(A,B,C,adapt);
    case 24: if (24<=K) return matmul_function<M,N,24,R,L,T,F,Func>()(A,B,C,adapt);
    case 26: if (26<=K) return matmul_function<M,N,26,R,L,T,F,Func>()(A,B,C,adapt);
    default: { return matmul_function<M,N,0,R,L,T,F,Func>()(A,B,C,adapt); }
  }
}

template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_switch_function { };
template<int M,int N,int K,      int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_switch_function<M,N,K,0,L,T,F,Func> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return real_matmul_switch   <M,N,K,0,L,T,F,Func>(X,Y,C,adapt); } };
template<int M,int N,int K,      int L,class T,template<class,class,class> class F,template<class,class> class Func> struct matmul_switch_function<M,N,K,1,L,T,F,Func> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const Adapt &adapt) { return complex_matmul_switch<M,N,K,1,L,T,F,Func>(X,Y,C,adapt); } };

template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3>
void inline matmul_switch(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const id_adaptor_function &adapt)
{
  matmul_switch_function<M,N,K,R,L,T,F,Func>()(A,B,C,adapt);
}

template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt        > void matmul_switch(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const add_adaptor_function  <Adapt  > &adapt)
{
  typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
  array_type Z(C.size());
  typename SubArray<array_type>::self Z2=sub(Z,0,0,Z.nrows(),Z.ncols());
  matmul_switch_function<M,N,K,R,L,T,F,Func>()(data(X),data(Y),Z2,id_adapt());
  copy<id_function>(Z,C,adapt);
}

template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3                    > void matmul_switch(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const neg_adaptor_function            &adapt)
{
  typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
  array_type Z(C.size());
  typename SubArray<array_type>::self Z2=sub(Z,0,0,Z.nrows(),Z.ncols());
  matmul_switch_function<M,N,K,R,L,T,F,Func>()(data(X),data(Y),Z2,id_adapt());
  copy<id_function>(Z,C,adapt);
}
template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,            class V> void matmul_switch(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const ax_adaptor_function   <V      > &adapt)
{
  typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
  array_type Z(C.size());
  typename SubArray<array_type>::self Z2=sub(Z,0,0,Z.nrows(),Z.ncols());
  matmul_switch_function<M,N,K,R,L,T,F,Func>()(data(X),data(Y),Z2,id_adapt());
  copy<id_function>(Z,C,adapt);
}
template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt        > void matmul_switch(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const sub_adaptor_function  <Adapt  > &adapt)
{
  typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
  array_type Z(C.size());
  typename SubArray<array_type>::self Z2=sub(Z,0,0,Z.nrows(),Z.ncols());
  matmul_switch_function<M,N,K,R,L,T,F,Func>()(data(X),data(Y),Z2,id_adapt());
  copy<id_function>(Z,C,adapt);
}
template<int M,int N,int K,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt,class V> void matmul_switch(const Matrix<G1> &X,const Matrix<G2> &Y,Matrix<G3> &C,const scale_adaptor_function<Adapt,V> &adapt)
{
  typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
  array_type Z(C.size());
  typename SubArray<array_type>::self Z2=sub(Z,0,0,Z.nrows(),Z.ncols());
  matmul_switch_function<M,N,K,R,L,T,F,Func>()(data(X),data(Y),Z2,id_adapt());
  copy<id_function>(Z,C,adapt);
}



template<int M,int N,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const Adapt &adapt)
{
  typedef typename Matrix<G3>::value_type value_type;
  int m=C.nrows(), n=C.ncols(), k=X.size()/m; assert(Y.size()/n==k);
  int dk=k%K2;

  int h=0;
  if (dk)            { matmul_switch  <M,N,K2,R,L,T,F,Func>  (data_matrix(m,dk,&X[0  ]),data_matrix(n,dk,&Y[0  ]),C,id_adapt ()); h+=dk; }
  else               { matmul_function<M,N,K2,R,L,T,F,Func>()(data_matrix(m,K2,&X[0  ]),data_matrix(n,K2,&Y[0  ]),C,id_adapt ()); h+=K2; }
  for (; h<k; h+=K2) { matmul_function<M,N,K2,R,L,T,F,Func>()(data_matrix(m,K2,&X[h*m]),data_matrix(n,K2,&Y[h*n]),C,add_adapt()); }
}

template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3>
void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C)
{
  typedef typename Matrix<G3>::value_type value_type;
  int m=C.nrows(), n=C.ncols(), k=X.size()/m; assert(Y.size()/n==k);
  int dk=k%K;

  int h=0;
  if (dk)           { mul            <M,N,K2,R,L,T,F,Func>  (data_vector(m*dk,&X[0  ]),data_vector(n*dk,&Y[0  ]),C,id_adapt ()); h+=dk; }
  else              { matmul_function<M,N,K ,R,L,T,F,Func>()(data_matrix(m, K,&X[0  ]),data_matrix(n, K,&Y[0  ]),C,id_adapt ()); h+= K; }
  for (; h<k; h+=K) { matmul_function<M,N,K ,R,L,T,F,Func>()(data_matrix(m, K,&X[h*m]),data_matrix(n, K,&Y[h*n]),C,add_adapt()); }
}


template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
struct adapt_mul_body
{
  typedef Matrix<G3> array_type;
  const Vector<G1> &X; const Vector<G2> &Y; Matrix<G3> &C; const Adapt &adapt;
  inline adapt_mul_body(const Vector<G1> &x, const Vector<G2> &y, Matrix<G3> &c, const Adapt &a) : X(x),Y(y),C(c),adapt(a) {}
  inline void operator()(int i,int imax) const
  {
    int m=C.nrows(), n=C.ncols(), k=X.size()/m, dm=(m-1)%M+1, dn=n%N; assert(Y.size()/m==k);
    i=(i>0)*(M*(i-1)+dm);
    imax=M*(imax-1)+dm;
    dm=(imax-i)%M;
    typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
    array_type S(M,N);
    if (dm)
    {
      typename SubArray<const Vector<G1> >::self X2 =sub(X,0,dm*k);
      if (dn)
      {
        typename SubArray<array_type >::self S2=sub(S,0,0,dm,dn);
        mul<0,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),S2);
        typename SubArray<Matrix<G3> >::self C2=sub(C,0,0,dm,dn);
        copy<id_function>(S2,C2,adapt);
      }
      {
        typename SubArray<array_type >::self S2=sub(S,0,0,dm, N);
        for (int j=dn; j<n; j+=N)
        {
          mul<0,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),S2);
          typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,dm, N);
          copy<id_function>(S2,C2,adapt);
        }
      }
    }

    for (int i=dm; i<m; i+=M)
    {
      typename SubArray<const Vector<G1> >::self X2 =sub(X,i*k,M*k);
      if (dn)
      {
        typename SubArray<array_type >::self S2=sub(S,0,0, M,dn);
        mul<M,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),S2);
        typename SubArray<Matrix<G3> >::self C2=sub(C,i,0, M,dn);
        copy<id_function>(S2,C2,adapt);
      }
      {
        typename SubArray<array_type >::self S2=sub(S,0,0, M, N);
        for (int j=dn; j<n; j+=N)
        {
          mul<M,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),S2);
          typename SubArray<Matrix<G3> >::self C2=sub(C,i,j, M, N);
          copy<id_function>(S2,C2,adapt);
        }
      }
    }
  }
};

template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void adapt_mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const Adapt &adapt)
{
#if defined(BLAS_THREADING)
  adapt_mul_body<M,N,K,K2,R,L,T,F,Func,G1,G2,G3,Adapt> func(X,Y,C,adapt);
  int n=(C.nrows()-1)/M+1;
  parallel_for(0,n,func,(n<=2)+1);
#else
  int m=C.nrows(), n=C.ncols(), k=X.size()/m, dm=m%M, dn=n%N; assert(Y.size()/m==k);
  typedef typename Matrix<G3>::value_type value_type;
  typedef typename DenseMatrix<typename Matrix<G3>::value_type>::self array_type;
  array_type S(M,N);

  if (dm)
  {
    typename SubArray<const Vector<G1> >::self X2 =sub(X,0,dm*k);
    if (dn)
    {
      typename SubArray<array_type >::self S2=sub(S,0,0,dm,dn);
      mul<0,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),S2);
      typename SubArray<Matrix<G3> >::self C2=sub(C,0,0,dm,dn);
      copy<id_function>(S2,C2,adapt);
    }
    {
      typename SubArray<array_type >::self S2=sub(S,0,0,dm, N);
      for (int j=dn; j<n; j+=N)
      {
        mul<0,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),S2);
        typename SubArray<Matrix<G3> >::self C2=sub(C,0,j,dm, N);
        copy<id_function>(S2,C2,adapt);
      }
    }
  }

  for (int i=dm; i<m; i+=M)
  {
    typename SubArray<const Vector<G1> >::self X2 =sub(X,i*k,M*k);
    if (dn)
    {
      typename SubArray<array_type >::self S2=sub(S,0,0, M,dn);
      mul<M,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),S2);
      typename SubArray<Matrix<G3> >::self C2=sub(C,i,0, M,dn);
      copy<id_function>(S2,C2,adapt);
    }
    {
      typename SubArray<array_type >::self S2=sub(S,0,0, M, N);
      for (int j=dn; j<n; j+=N)
      {
        mul<M,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),S2);
        typename SubArray<Matrix<G3> >::self C2=sub(C,i,j, M, N);
        copy<id_function>(S2,C2,adapt);
      }
    }
  }
#endif
}

template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
struct mul_body
{
  typedef Matrix<G3> array_type;
  const Vector<G1> &X; const Vector<G2> &Y; Matrix<G3> &C; const Adapt &adapt;
  inline mul_body(const Vector<G1> &x, const Vector<G2> &y, Matrix<G3> &c, const Adapt &a) : X(x),Y(y),C(c),adapt(a) {}
  inline void operator()(int i,int imax) const
  {
    int m=C.nrows(), n=C.ncols(), k=X.size()/m, dm=(m-1)%M+1, dn=n%N; assert(Y.size()/m==k);
    i=(i>0)*(M*(i-1)+dm);
    imax=M*(imax-1)+dm;
    dm=(imax-i)%M;
    if (dm)
    {
      typename SubArray<const Vector<G1> >::self X2 =sub(X,i*k,dm*k);
      if (dn)                   { typename SubArray<array_type>::self C2=sub(C,i,0,dm,dn); mul<0,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),C2); }
      for (int j=dn; j<n; j+=N) { typename SubArray<array_type>::self C2=sub(C,i,j,dm, N); mul<0,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),C2); }
    }
    for (i+=dm; i<imax; i+=M)
    {
      typename SubArray<const Vector<G1> >::self X2 =sub(X,i*k,M*k);
      if (dn)                   { typename SubArray<array_type>::self C2=sub(C,i,0, M,dn); mul<M,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),C2); }
      for (int j=dn; j<n; j+=N) { typename SubArray<array_type>::self C2=sub(C,i,j, M, N); mul<M,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),C2); }
    }
  }
};

template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3>
void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const id_adaptor_function &adapt)
{
#if defined(BLAS_THREADING)
  mul_body<M,N,K,K2,R,L,T,F,Func,G1,G2,G3,id_adaptor_function> func(X,Y,C,adapt);
  int n=(C.nrows()-1)/M+1;
  parallel_for(0,n,func,(n<=2)+1);
#else
  typedef Matrix<G3> array_type;
  int m=C.nrows(), n=C.ncols(), k=X.size()/m, dm=m%M, dn=n%N; assert(Y.size()/m==k);
  if (dm)
  {
    typename SubArray<const Vector<G1> >::self X2 =sub(X,0,dm*k);
    if (dn)                   { typename SubArray<array_type>::self C2=sub(C,0,0,dm,dn); mul<0,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),C2); }
    for (int j=dn; j<n; j+=N) { typename SubArray<array_type>::self C2=sub(C,0,j,dm, N); mul<0,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),C2); }
  }
  for (int i=dm; i<m; i+=M)
  {
    typename SubArray<const Vector<G1> >::self X2 =sub(X,i*k,M*k);
    if (dn)                   { typename SubArray<array_type>::self C2=sub(C,i,0, M,dn); mul<M,0,K,K2,R,L,T,F,Func>(X2,sub(Y,0  ,dn*k),C2); }
    for (int j=dn; j<n; j+=N) { typename SubArray<array_type>::self C2=sub(C,i,j, M, N); mul<M,N,K,K2,R,L,T,F,Func>(X2,sub(Y,j*k,N *k),C2); }
  }
#endif
}

template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3                    > void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const neg_adaptor_function            &adapt) { return adapt_mul<M,N,K,K2,R,L,T,F,Func>(X,Y,C,adapt); }
template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,            class V> void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const ax_adaptor_function   <V      > &adapt) { return adapt_mul<M,N,K,K2,R,L,T,F,Func>(X,Y,C,adapt); }
template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt        > void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const add_adaptor_function  <Adapt  > &adapt) { return adapt_mul<M,N,K,K2,R,L,T,F,Func>(X,Y,C,adapt); }
template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt        > void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const sub_adaptor_function  <Adapt  > &adapt) { return adapt_mul<M,N,K,K2,R,L,T,F,Func>(X,Y,C,adapt); }
template<int M,int N,int K,int K2,int R,int L,class T,template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt,class V> void mul(const Vector<G1> &X,const Vector<G2> &Y,Matrix<G3> &C,const scale_adaptor_function<Adapt,V> &adapt) { return adapt_mul<M,N,K,K2,R,L,T,F,Func>(X,Y,C,adapt); }


template<class T,int Copy=0>
struct complex_mul_copy_adaptor : public two_arrays_generator<T,T,Copy>
{
  typedef complex_mul_copy_adaptor self;
  typedef two_arrays_generator<T,T,Copy> base;
  TWO_ARRAYS_BASE_TYPES

  typedef id_adaptor<T> adaptor_type;

  template<class T2,int C2=Copy> struct rebind { typedef complex_mul_copy_adaptor<T2,C2> other; };

  using base::first_array;
  using base::second_array;
  inline complex_mul_copy_adaptor(first_array_type &x, second_array_type &y) : base(x,y) {}
  inline complex_mul_copy_adaptor(const self &x) : base(x) {}

  template<class V> inline void operator()(const index_type &i, const V &x) { first_array()[i]=real(x); second_array()[i]=imag(x); }

  typedef complex_mul_copy_adaptor<typename SubArray<T>::self,1> sub_type;
  inline sub_type sub(index_type i, size_type n) { typename SubArray<T>::self A=::sub(first_array (),i,n), B=::sub(second_array(),i,n); return sub_type(A,B); }
};

template<class T,int C=0> struct ComplexMulCopyAdaptor { typedef complex_mul_copy_adaptor<T,C> self; };


struct complex_mul_copy_adaptor_function
{
  template<class T,int C=0> struct result_rebind { typedef typename ComplexMulCopyAdaptor<typename SubArray<T,C>::self,1>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { typename SubArray<Vector<G> >::self X1=sub(X,0,X.size()/2), X2=sub(X,X.size()/2,X.size()/2); return typename result_rebind<Vector<G> >::other(X1,X2); }
};
inline complex_mul_copy_adaptor_function complex_mul_copy_adapt() { return complex_mul_copy_adaptor_function(); }




template<class G1,class G2,class Adapt>
void mul_copy(const Matrix<G1> &A, Vector<G2> &X, const Adapt &adapt, int M,int d1,int d2,int a=1)
{
  typedef typename Vector<G2>::value_type value_type;

  typedef typename DataMatrix<typename Vector<G2>::value_type>::self array_type;

  int m=A.nrows(), dm=m%M;
  int n1=A.ncols(), n2=n1%d1;
  int n3=n2%d2;
  int d3=n3;

  int d=0;
  if (dm)
  {
    if (n3)                     { array_type B(dm,a*d3,&X[d]); d+=B.nelms(); copy<id_function>(sub(A,0,0,dm,n3),B,adapt); }
    for (int j=n3; j<n2; j+=d2) { array_type B(dm,a*d2,&X[d]); d+=B.nelms(); copy<id_function>(sub(A,0,j,dm,d2),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B(dm,a*d1,&X[d]); d+=B.nelms(); copy<id_function>(sub(A,0,j,dm,d1),B,adapt); }
  }
  for (int i=dm; i<m; i+=M)
  {
    if (n3)                     { array_type B( M,a*d3,&X[d]); d+=B.nelms(); copy<id_function>(sub(A,i,0, M,n3),B,adapt); }
    for (int j=n3; j<n2; j+=d2) { array_type B( M,a*d2,&X[d]); d+=B.nelms(); copy<id_function>(sub(A,i,j, M,d2),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B( M,a*d1,&X[d]); d+=B.nelms(); copy<id_function>(sub(A,i,j, M,d1),B,adapt); }
  }
}

template<int L,template<class,class> class F,class G1,class G2,class Adapt>
void mul_fill_copy(const Matrix<G1> &A, Vector<G2> &X, const Adapt &adapt, int M,int d1,int d2,int a=1)
{
  typedef typename Vector<G1>::value_type value_type;

  typedef typename DataMatrix<typename Vector<G2>::value_type>::self array_type;

  int m=A.nrows(), dm=m%M;
  int n1=A.ncols(), n2=n1%d1;
  if (d1-n2<L) d2=d1;
  int n3=n2%d2;
  int d3=((n3+(L-1))/L)*L;

  int d=0;
  if (dm)
  {
    if (n3)                     { array_type B(dm,a*d3,&X[d]); d+=B.nelms(); fill_copy<F>(sub(A,0,0,dm,n3),B,adapt,value_type(0)); }
    for (int j=n3; j<n2; j+=d2) { array_type B(dm,a*d2,&X[d]); d+=B.nelms(); copy     <F>(sub(A,0,j,dm,d2),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B(dm,a*d1,&X[d]); d+=B.nelms(); copy     <F>(sub(A,0,j,dm,d1),B,adapt); }
  }
  for (int i=dm; i<m; i+=M)
  {
    if (n3)                     { array_type B( M,a*d3,&X[d]); d+=B.nelms(); fill_copy<F>(sub(A,i,0, M,n3),B,adapt,value_type(0)); }
    for (int j=n3; j<n2; j+=d2) { array_type B( M,a*d2,&X[d]); d+=B.nelms(); copy     <F>(sub(A,i,j, M,d2),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B( M,a*d1,&X[d]); d+=B.nelms(); copy     <F>(sub(A,i,j, M,d1),B,adapt); }
  }
}

template<int L,class G1,class G2,class Adapt>
void mul_trn(const Matrix<G1> &A, Vector<G2> &X, const Adapt &adapt, int M,int d1,int d2,int a=1)
{
  typedef typename Vector<G2>::value_type value_type;

  typedef typename DataMatrix<typename Vector<G2>::value_type>::self array_type;

  int m=L*A.ncols(), dm=m%M;
  int n1=A.nrows(), n2=n1%d1;
  int n3=n2%d2;
  int d3=n3;

  int d=0;
  if (dm)
  {
    if (n3)                     { array_type B(dm,a*d3/L,&X[d]); d+=B.nelms(); trn_function<L,id_function>()(sub(A,0,0  ,n3,dm/L),B,adapt); }
    for (int j=n3; j<n2; j+=d2) { array_type B(dm,a*d2/L,&X[d]); d+=B.nelms(); trn_function<L,id_function>()(sub(A,j,0  ,d2,dm/L),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B(dm,a*d1/L,&X[d]); d+=B.nelms(); trn_function<L,id_function>()(sub(A,j,0  ,d1,dm/L),B,adapt); }
  }
  for (int i=dm; i<m; i+=M)
  {
    if (n3)                     { array_type B( M,a*d3/L,&X[d]); d+=B.nelms(); trn_function<L,id_function>()(sub(A,0,i/L,n3, M/L),B,adapt); }
    for (int j=n3; j<n2; j+=d2) { array_type B( M,a*d2/L,&X[d]); d+=B.nelms(); trn_function<L,id_function>()(sub(A,j,i/L,d2, M/L),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B( M,a*d1/L,&X[d]); d+=B.nelms(); trn_function<L,id_function>()(sub(A,j,i/L,d1, M/L),B,adapt); }
  }
}

template<int L,template<class,class> class F,class G1,class G2,class Adapt>
void mul_fill_trn(const Matrix<G1> &A, Vector<G2> &X, const Adapt &adapt, int M,int d1,int d2,int a=1)
{
  typedef typename Vector<G1>::value_type value_type;

  typedef typename DataMatrix<typename Vector<G2>::value_type>::self array_type;

  int m=A.ncols(), dm=m%M;
  int n1=A.nrows(), n2=n1%d1;
  if (d1-n2<L) d2=d1;
  int n3=n2%d2;
  int d3=((n3+(L-1))/L)*L;

  int d=0;
  if (dm)
  {
    if (n3)                     { array_type B(dm,a*d3,&X[d]); d+=B.nelms(); fill_trn<F>(sub(A,0,0,n3,dm),B,adapt,value_type(0)); }
    for (int j=n3; j<n2; j+=d2) { array_type B(dm,a*d2,&X[d]); d+=B.nelms(); trn     <F>(sub(A,j,0,d2,dm),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B(dm,a*d1,&X[d]); d+=B.nelms(); trn     <F>(sub(A,j,0,d1,dm),B,adapt); }
  }
  for (int i=dm; i<m; i+=M)
  {
    if (n3)                     { array_type B( M,a*d3,&X[d]); d+=B.nelms(); fill_trn<F>(sub(A,0,i,n3, M),B,adapt,value_type(0)); }
    for (int j=n3; j<n2; j+=d2) { array_type B( M,a*d2,&X[d]); d+=B.nelms(); trn     <F>(sub(A,j,i,d2, M),B,adapt); }
    for (int j=n2; j<n1; j+=d1) { array_type B( M,a*d1,&X[d]); d+=B.nelms(); trn     <F>(sub(A,j,i,d1, M),B,adapt); }
  }
}


template<int TransA,int TransB,int M,int N,int K,int K2,int L,class T,template<class,class,class> class F,class G1,class G2,class G3,class Adapt>
void simd_real_mul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt)
{
  int m=(!TransA)?A.nrows():A.ncols(), k=(!TransA)?A.ncols():A.nrows(), n=(!TransB)?A.ncols():A.nrows();
  assert( ((!TransB)?B.nrows():B.ncols())==k && C.nrows()==m && C.ncols()==n);
  int l=L*((k-1)/L+1);

  if (k%L==0 && k/L<=K2 && k*(m+n)<=L*K*(M+N) )
  {
    if (!TransA && !TransB)
    {
      typedef typename DenseMatrix<T>::self array_type;
      array_type B2(n,k); typename SimdBlock<array_type>::self B3(B2); trn_function<L,id_function>()(simd_block(B),B3,id_adapt());
      return matmul_switch<0,0,K2,0,L,T,F,sum_function>(simd_block(dense(A)),B3,C,adapt);
    }
    else if (!TransA &&  TransB) return matmul_switch<0,0,K2,0,L,T,F,sum_function>(simd_block(data(dense(A))),simd_block(data(dense(B))),C,adapt);
    else if ( TransA && !TransB)
    {
      typedef typename DenseMatrix<T>::self array_type;
      array_type A2(m,k); typename SimdBlock<array_type>::self A3(A2); trn_function<L,id_function>()(simd_block(A),A3,id_adapt());
      array_type B2(n,k); typename SimdBlock<array_type>::self B3(B2); trn_function<L,id_function>()(simd_block(B),B3,id_adapt());
      return matmul_switch<0,0,K2,0,L,T,F,sum_function>(A3,B3,C,adapt);
    }
    else if ( TransA &&  TransB)
    {
      typedef typename DenseMatrix<T>::self array_type;
      array_type A2(m,k); typename SimdBlock<array_type>::self A3(A2); trn_function<L,id_function>()(simd_block(A),A3,id_adapt());
      return matmul_switch<0,0,K2,0,L,T,F,sum_function>(A3,simd_block(dense(B)),C,adapt);
    }
  }

  typedef typename DenseVector<T>::self array_type;
  array_type X(m*l), Y(n*l);
  if (m%L==0 && n%L==0 && k%L==0)
  {
    typename SimdBlock<array_type>::self X2(X), Y2(Y);
    if (!TransA) mul_copy(simd_block(A),X2,id_adapt(),M,K,K2); else mul_trn<L>(simd_block(A),X2,id_adapt(),M,L*K,L*K2);
    if ( TransB) mul_copy(simd_block(B),Y2,id_adapt(),N,K,K2); else mul_trn<L>(simd_block(B),Y2,id_adapt(),N,L*K,L*K2);
  }
  else
  {
    if (!TransA) mul_fill_copy<L,id_function>(A,X,id_adapt(),M,L*K,L*K2); else mul_fill_trn<L,id_function>(A,X,id_adapt(),M,L*K,L*K2);
    if ( TransB) mul_fill_copy<L,id_function>(B,Y,id_adapt(),N,L*K,L*K2); else mul_fill_trn<L,id_function>(B,Y,id_adapt(),N,L*K,L*K2);
  }
  return mul<M,N,K,K2,0,L,T,F,sum_function>(simd_block(X),simd_block(Y),C,adapt);
}

template<int TransA,int TransB,int M,int N,int K,int K2,int L,class T,template<class,class,class> class F,class G1,class G2,class G3,class Adapt>
void simd_complex_mul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt)
{
  int m=(!TransA)?A.nrows():A.ncols(), k=(!TransA)?A.ncols():A.nrows(), n=(!TransB)?A.ncols():A.nrows();
  assert( ((!TransB)?B.nrows():B.ncols())==k && C.nrows()==m && C.ncols()==n);
  int l=4*L*((k-1)/(2*L)+1);

  typedef typename DenseVector<T>::self array_type;
  array_type X(m*l), Y(n*l);
       if (TransA==0) mul_fill_copy<2*L,  id_function>(A,X,complex_mul_copy_adapt(),M,2*L*(K/2),2*L*((K2+1)/2),2);
  else if (TransA==1) mul_fill_trn <2*L,  id_function>(A,X,complex_mul_copy_adapt(),M,2*L*(K/2),2*L*((K2+1)/2),2);
  else if (TransA==2) mul_fill_trn <2*L,conj_function>(A,X,complex_mul_copy_adapt(),M,2*L*(K/2),2*L*((K2+1)/2),2);
       if (TransB==0) mul_fill_trn <2*L,  id_function>(B,Y,complex_mul_copy_adapt(),N,2*L*(K/2),2*L*((K2+1)/2),2);
  else if (TransB==1) mul_fill_copy<2*L,  id_function>(B,Y,complex_mul_copy_adapt(),N,2*L*(K/2),2*L*((K2+1)/2),2);
  else if (TransB==2) mul_fill_copy<2*L,conj_function>(B,Y,complex_mul_copy_adapt(),N,2*L*(K/2),2*L*((K2+1)/2),2);
  return mul<M,N,2*(K/2),2*((K2+1)/2),1,2*L,T,F,sum_function>(simd_block(X),simd_block(Y),C,adapt);
}

template<int TransA,int TransB,int M,int N,int K,int K2,int L,class T,template<class,class,class> class F> struct simd_mul_function2                                        { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) const { return simd_real_mul   <TransA!=0,TransB!=0,M,N,K,K2,L,T,F>(A,B,C,adapt); } };
template<int TransA,int TransB,int M,int N,int K,int K2,int L,class T,template<class,class,class> class F> struct simd_mul_function2<TransA,TransB,M,N,K,K2,L,complex<T>,F> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) const { return simd_complex_mul<TransA,TransB,M,N,K,K2,L,T,F>(A,B,C,adapt); } };
template<int TransA,int TransB,int M,int N,int K,int K2,int L,class T,template<class,class,class> class F> struct simd_mul_function                                         { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) const { return simd_mul_function2<TransA,TransB,M,N,K,K2,L,T,F>()(A,B,C,adapt); } };
template<int TransA,int TransB,int M,int N,int K,int K2,      class T,template<class,class,class> class F>
struct simd_mul_function<TransA,TransB,M,N,K,K2,0,T,F>
{
  template<class G1,class G2,class G3,class Adapt>
  inline void operator()(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) const
  {
    int m=(!TransA)?A.nrows():A.ncols(), k=(!TransA)?A.ncols():A.nrows(), n=(!TransB)?A.ncols():A.nrows();
    assert( ((!TransB)?B.nrows():B.ncols())==k && C.nrows()==m && C.ncols()==n);

    typedef typename DenseVector<T>::self array_type;
    array_type X(m*k), Y(n*k);
         if (TransA==0) mul_fill_copy<1,  id_function>(A,X,id_adapt(),M,K,K2);
    else if (TransA==1) mul_fill_trn <1,  id_function>(A,X,id_adapt(),M,K,K2);
    else if (TransA==2) mul_fill_trn <1,conj_function>(A,X,id_adapt(),M,K,K2);
         if (TransB==0) mul_fill_trn <1,  id_function>(B,Y,id_adapt(),N,K,K2);
    else if (TransB==1) mul_fill_copy<1,  id_function>(B,Y,id_adapt(),N,K,K2);
    else if (TransB==2) mul_fill_copy<1,conj_function>(B,Y,id_adapt(),N,K,K2);
    mul<M,N,K,K2,0,0,T,F,id_function>(X,Y,C,adapt);
  }
};


template<int TransA,int TransB,class G1,class G2,class G3,class Adapt,class V> inline void aux_mul(const Matrix<G1> &A,const Matrix<G2> &B,Matrix<G3> &C,const Adapt &adapt,const V &,const V &) { return simd_mul_function<TransA,TransB,BLAS_M,BLAS_N,BLAS_K,BLAS_K2,SimdLength<V>::RET,V,dot_function>()(A,B,C,adapt); }

template<class G1,class G2,class G3,class Adapt> void aux_mul00(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<0,0>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul01(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<0,1>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul02(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<0,2>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul10(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<1,0>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul11(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<1,1>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul12(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<1,2>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul20(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<2,0>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul21(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<2,1>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul22(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt) { return aux_mul<2,2>(A,B,C,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<int TransA,int TransB,class G1,class G2,class G3,class Adapt> void aux_mul(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt)
{
       if (TransA==0) { if (TransB==0) return aux_mul00(A,B,C,adapt); else if (TransB==1) return aux_mul01(A,B,C,adapt); else if (TransB==2) return aux_mul02(A,B,C,adapt); }
  else if (TransA==1) { if (TransB==0) return aux_mul10(A,B,C,adapt); else if (TransB==1) return aux_mul11(A,B,C,adapt); else if (TransB==2) return aux_mul12(A,B,C,adapt); }
  else if (TransA==2) { if (TransB==0) return aux_mul20(A,B,C,adapt); else if (TransB==1) return aux_mul21(A,B,C,adapt); else if (TransB==2) return aux_mul22(A,B,C,adapt); }
}

template<int TransA,int TransB,class G1,class G2,class G3,class Adapt>
inline void mul(const Matrix<G1> &A, const Matrix<G2> &B,Matrix<G3> &C, const Adapt &adapt)
{
  typename ArrayData<Matrix<G3> >::self C2 = data(C);
  return aux_mul<TransA,TransB>(data(A),data(B),C2,adapt);
}

template<int TransA,int TransB,        class G1,class G2,class G3> inline void mul    (           const Matrix<G1> &A,const Matrix<G2> &B,           Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,id_adapt       (   )); }
template<int TransA,int TransB,        class G1,class G2,class G3> inline void neg_mul(           const Matrix<G1> &A,const Matrix<G2> &B,           Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,neg_adapt      (   )); }
template<int TransA,int TransB,class V,class G1,class G2,class G3> inline void mul    (const V &a,const Matrix<G1> &A,const Matrix<G2> &B,           Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,ax_adapt       (  a)); }

template<int TransA,int TransB,        class G1,class G2,class G3> inline void add_mul(           const Matrix<G1> &A,const Matrix<G2> &B,           Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,add_adapt      (   )); }
template<int TransA,int TransB,        class G1,class G2,class G3> inline void sub_mul(           const Matrix<G1> &A,const Matrix<G2> &B,           Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,sub_adapt      (   )); }
template<int TransA,int TransB,class V,class G1,class G2,class G3> inline void add_mul(const V &a,const Matrix<G1> &A,const Matrix<G2> &B,           Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,add_adapt      (  a)); }

template<int TransA,int TransB,class V,class G1,class G2,class G3> inline void add_mul(           const Matrix<G1> &A,const Matrix<G2> &B,const V &s,Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,scale_add_adapt(s  )); }
template<int TransA,int TransB,class V,class G1,class G2,class G3> inline void sub_mul(           const Matrix<G1> &A,const Matrix<G2> &B,const V &s,Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,scale_sub_adapt(s  )); }
template<int TransA,int TransB,class V,class G1,class G2,class G3> inline void add_mul(const V &a,const Matrix<G1> &A,const Matrix<G2> &B,const V &s,Matrix<G3> &C) { return mul<TransA,TransB>(A,B,C,scale_add_adapt(s,a)); }


//{unsecret}
//Description: C=a*mul(A,B)+b*C
//Remarks:
//  TransA and TransB are template arguments to transpose A and B.
//  (0=No Transpose, 1=Transpose, 2=Conjugate Transpose)
//Example:
//  gemm     (  A,B,  C);  // C =  mul(A,     B ), equivalent to gemm<0,0>(1,A,B,0,C)
//  gemm     (2,A,B,  C);  // C =2*mul(A,     B ), equivalent to gemm<0,0>(2,A,B,0,C)
//  gemm     (2,A,B,1,C);  // C+=2*mul(A,     B ), equivalent to gemm<0,0>(2,A,B,1,C)
//  gemm<0,1>(  A,B,  C);  // C =  mul(A,trn (B)), equivalent to gemm<0,1>(1,A,B,0,C)
//  gemm<0,1>(2,A,B,  C);  // C =2*mul(A,trn (B)), equivalent to gemm<0,1>(2,A,B,0,C)
//  gemm<0,1>(2,A,B,1,C);  // C+=2*mul(A,trn (B))
//  gemm<0,2>(  A,B,  C);  // C =  mul(A,htrn(B)), equivalent to gemm<0,2>(1,A,B,0,C)
//  gemm<0,2>(2,A,B,  C);  // C =2*mul(A,htrn(B)), equivalent to gemm<0,2>(2,A,B,0,C)
//  gemm<0,2>(2,A,B,1,C);  // C+=2*mul(A,htrn(B))
template<int TransA,int TransB,class G1,class G2,class G3> void gemm(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Matrix<G2> &B, const typename Matrix<G3>::value_type &b, Matrix<G3> &C)
{
  typedef typename Matrix<G3>::value_type V;
       if (a==V( 1)) { if (b==V(0)) mul    <TransA,TransB>(  A,B,C); else if (b==V(1)) add_mul<TransA,TransB>(  A,B,C); else add_mul<TransA,TransB>(  A,B,b,C); }
  else if (a==V(-1)) { if (b==V(0)) neg_mul<TransA,TransB>(  A,B,C); else if (b==V(1)) sub_mul<TransA,TransB>(  A,B,C); else sub_mul<TransA,TransB>(  A,B,b,C); }
  else if (a!=V( 0)) { if (b==V(0)) mul    <TransA,TransB>(a,A,B,C); else if (b==V(1)) add_mul<TransA,TransB>(a,A,B,C); else add_mul<TransA,TransB>(a,A,B,b,C); }
  else scal(b,C);
}
template<int TransA,int TransB,class G1,class G2,class G3> void gemm(                                          const Matrix<G1> &A, const Matrix<G2> &B, const typename Matrix<G3>::value_type &b, Matrix<G3> &C) { typedef typename Matrix<G3>::value_type V; if (b==V(0)) mul<TransA,TransB>(A,B,C); else if (b==V( 1)) add_mul<TransA,TransB>(A,B,C); else add_mul<TransA,TransB>(A,B,b,C); }
template<int TransA,int TransB,class G1,class G2,class G3> void gemm(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Matrix<G2> &B,                                           Matrix<G3> &C) { typedef typename Matrix<G3>::value_type V; if (a==V(1)) mul<TransA,TransB>(A,B,C); else if (a==V(-1)) neg_mul<TransA,TransB>(A,B,C); else if (a!=V(0)) mul<TransA,TransB>(a,A,B,C); else fill(C,typename Matrix<G3>::value_type(0)); }
template<int TransA,int TransB,class G1,class G2,class G3> void gemm(                                          const Matrix<G1> &A, const Matrix<G2> &B,                                           Matrix<G3> &C) { typedef typename Matrix<G3>::value_type V; mul<TransA,TransB>(A,B,C); }

template<class G1,class G2,class G3> void gemm(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Matrix<G2> &B, const typename Matrix<G3>::value_type &b, Matrix<G3> &C) { return gemm<0,0>(a,A,B,b,C); }
template<class G1,class G2,class G3> void gemm(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Matrix<G2> &B,                                           Matrix<G3> &C) { return gemm<0,0>(a,A,B,  C); }
template<class G1,class G2,class G3> void gemm(                                          const Matrix<G1> &A, const Matrix<G2> &B, const typename Matrix<G3>::value_type &b, Matrix<G3> &C) { return gemm<0,0>(  A,B,b,C); }
template<class G1,class G2,class G3> void gemm(                                          const Matrix<G1> &A, const Matrix<G2> &B,                                           Matrix<G3> &C) { return gemm<0,0>(  A,B,  C); }


template<class V> void gemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, V a, const V *pA, int lda, const V *pB, int ldb, V b, V *pC, int ldc)
{
  V *pX, *pY, *pZ;
  typename DenseMatrix<V>::self X,Y,Z;
  if (lda==k) pX=(V *)pA; else { X.resize(m,k); pX=&X(0,0); copy(sub(data_matrix(m,lda,(V*)pA),0,0,m,k),X); }
  if (ldb==n) pY=(V *)pB; else { Y.resize(k,n); pY=&Y(0,0); copy(sub(data_matrix(k,ldb,(V*)pB),0,0,k,n),Y); }
  if (ldc==n) pZ=(V *)pC; else { Z.resize(m,n); pZ=&Z(0,0); }
  typename DataMatrix<V>::self A(m,k,pX), B(k,n,pY), C(m,n,pZ);

  if (Order==CblasRowMajor)
  {
         if (TransA==CblasNoTrans  ) { if (TransB==CblasNoTrans) gemm<0,0>(a,A,B,b,C); else if (TransB==CblasTrans) gemm<0,1>(a,A,B,b,C); else if (TransB==CblasConjTrans) gemm<0,2>(a,A,B,b,C); }
    else if (TransA==CblasTrans    ) { if (TransB==CblasNoTrans) gemm<1,0>(a,A,B,b,C); else if (TransB==CblasTrans) gemm<1,1>(a,A,B,b,C); else if (TransB==CblasConjTrans) gemm<1,2>(a,A,B,b,C); }
    else if (TransA==CblasConjTrans) { if (TransB==CblasNoTrans) gemm<2,0>(a,A,B,b,C); else if (TransB==CblasTrans) gemm<2,1>(a,A,B,b,C); else if (TransB==CblasConjTrans) gemm<2,2>(a,A,B,b,C); }
  }
  else
  {
         if (TransB==CblasNoTrans  ) { if (TransA==CblasNoTrans) gemm<0,0>(a,B,A,b,C); else if (TransA==CblasTrans) gemm<0,1>(a,B,A,b,C); else if (TransA==CblasConjTrans) gemm<0,2>(a,B,A,b,C); }
    else if (TransB==CblasTrans    ) { if (TransA==CblasNoTrans) gemm<1,0>(a,B,A,b,C); else if (TransA==CblasTrans) gemm<1,1>(a,B,A,b,C); else if (TransA==CblasConjTrans) gemm<1,2>(a,B,A,b,C); }
    else if (TransB==CblasConjTrans) { if (TransA==CblasNoTrans) gemm<2,0>(a,B,A,b,C); else if (TransA==CblasTrans) gemm<2,1>(a,B,A,b,C); else if (TransA==CblasConjTrans) gemm<2,2>(a,B,A,b,C); }
  }

  if (ldc==n) return;
  else
  {
    typename DataMatrix<V>::self C2(m,ldc,pC);
    typename SubArray<typename DataMatrix<V>::self>::self C3=sub(C2,0,0,m,n);
    copy(C,C3);
  }
}

#ifdef AUTO_BLAS
template<class G,class T1,class T2> inline void       affect    (Matrix<G> &Y, const Matrix<multiply_matrix_generator<T1,T2> > &X) { mul    (X.generator().first_array(),X.generator().second_array(),Y); }
template<class G,class T1,class T2> inline Matrix<G> &operator+=(Matrix<G> &Y, const Matrix<multiply_matrix_generator<T1,T2> > &X) { add_mul(X.generator().first_array(),X.generator().second_array(),Y); return Y; }
template<class G,class T1,class T2> inline Matrix<G> &operator-=(Matrix<G> &Y, const Matrix<multiply_matrix_generator<T1,T2> > &X) { sub_mul(X.generator().first_array(),X.generator().second_array(),Y); return Y; }

template<class G,class T1,class T2,class Arg,          class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,neg_function    <Arg,    Res> > > &X) { neg_mul(                                 X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,val_mul_function<Arg,Val,Res> > > &X) { mul    (X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,mul_val_function<Arg,Val,Res> > > &X) { mul    (X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator+=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,val_mul_function<Arg,Val,Res> > > &X) { add_mul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator-=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,val_mul_function<Arg,Val,Res> > > &X) { sub_mul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator+=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,mul_val_function<Arg,Val,Res> > > &X) { add_mul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator-=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,T2> >,mul_val_function<Arg,Val,Res> > > &X) { sub_mul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array(),Y); return Y; }


template<class G,class T1,class T2> inline void       affect    (Matrix<G> &Y, const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > > &X) { tmul(X.generator().first_array(),X.generator().second_array().generator().array(),Y); }
template<class G,class T1,class T2> inline Matrix<G> &operator+=(Matrix<G> &Y, const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > > &X) { tmul(X.generator().first_array(),X.generator().second_array().generator().array(),Y); return Y; }
template<class G,class T1,class T2> inline Matrix<G> &operator-=(Matrix<G> &Y, const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > > &X) { tmul(X.generator().first_array(),X.generator().second_array().generator().array(),Y); return Y; }

template<class G,class T1,class T2,class Arg,          class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,neg_function    <Arg,    Res> > > &X) { neg_tmul(                                 X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,val_mul_function<Arg,Val,Res> > > &X) { tmul    (X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,mul_val_function<Arg,Val,Res> > > &X) { tmul    (X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator+=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,val_mul_function<Arg,Val,Res> > > &X) { add_tmul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator-=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,val_mul_function<Arg,Val,Res> > > &X) { sub_tmul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator+=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,mul_val_function<Arg,Val,Res> > > &X) { add_tmul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Matrix<G> &operator-=(Matrix<G> &Y, const Matrix<function_array_generator<const Matrix<multiply_matrix_generator<T1,const Matrix<trans_matrix_generator<T2> > > >,mul_val_function<Arg,Val,Res> > > &X) { sub_tmul(X.generator().function().value(),X.generator().array().generator().first_array(),X.generator().array().generator().second_array().generator().array(),Y); return Y; }
#endif


#ifdef BLAS_PRECOMPILE

#define DECLARE_BLAS_GEMM(V) \
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);


#define DEFINE_BLAS_GEMM(V) \
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul00(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<0,0>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul01(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<0,1>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul02(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<0,2>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul10(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<1,0>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul11(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<1,1>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul12(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<1,2>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul20(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<2,0>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul21(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<2,1>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const id_adaptor_function                                &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const neg_adaptor_function                               &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const ax_adaptor_function   <                        V > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }\
void aux_mul22(const Matrix<data_matrix_generator<V > > &A,const Matrix<data_matrix_generator<V > > &B,Matrix<data_matrix_generator<V > > &C,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul<2,2>(A,B,C,adapt,V(),V()); }




DECLARE_BLAS_GEMM(float)
DECLARE_BLAS_GEMM(double)
DECLARE_BLAS_GEMM(complex<float>)
DECLARE_BLAS_GEMM(complex<double>)

#endif //BLAS_PRECOMPILE

#endif //__cplusplus



#ifdef BLAS_PRECOMPILE
#ifdef __cplusplus
extern "C" {
#endif

//Group = CBLAS Level 3

//{unsecret}
GENIAL_API void BLAS_NAME(sgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const float  a, const float  *pA, int lda, const float  *pB, int ldb, const float  b, float  *pC, int ldc);
//{unsecret}
GENIAL_API void BLAS_NAME(dgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const double a, const double *pA, int lda, const double *pB, int ldb, const double b, double *pC, int ldc);
//{unsecret}
GENIAL_API void BLAS_NAME(cgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const void  *a, const void   *pA, int lda, const void   *pB, int ldb, const void  *b, void   *pC, int ldc);
//{unsecret}
GENIAL_API void BLAS_NAME(zgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const void  *a, const void   *pA, int lda, const void   *pB, int ldb, const void  *b, void   *pC, int ldc);

#ifdef __cplusplus
}
#endif
#endif //BLAS_PRECOMPILE

#endif //GEMM_H
