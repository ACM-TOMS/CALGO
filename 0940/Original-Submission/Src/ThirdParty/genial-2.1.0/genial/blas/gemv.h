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


#ifndef GEMV_H
#define GEMV_H

#ifdef __cplusplus

#include "copy.h"
#include "dot.h"
#include "array/matrix.h"


//Group = Linear Algebra

template<template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void mul(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Vector<G2>::value_type) value_type;

  typedef typename Adapt::template result_rebind<Vector<G3> >::other adaptor_type;
  adaptor_type adapt = adaptor(Y);

  int m=A.nrows();
  assert(A.ncols()==X.size() && A.nrows()==Y.size());

  int i=0;
  for (int imax=m-2; i<=imax; i+=2)
  {
    typename Vector<G3>::value_type y0,y1;
    dot<F,Func>(row(A,i),row(A,i+1),X,y0,y1);
    adapt(i,y0); adapt(i+1,y1);
  }
  for (; i<m; ++i)
  {
    typename Vector<G3>::value_type y;
    dot<F,Func>(row(A,i),X, y);
    adapt(i,y);
  }
}

template<template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class G3,class Adapt>
void tmul(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adaptor)
{
  typedef PROMOTE2(typename Matrix<G1>::value_type,typename Vector<G2>::value_type) value_type;

  F<value_type,value_type,value_type> f;
  Func<value_type,typename Vector<G3>::value_type> func;
  typedef typename Adapt::template result_rebind<Vector<G3> >::other adaptor_type;
  adaptor_type adapt = adaptor(Y);

  int m=A.nrows(), n=A.ncols();
  typename DenseVector<value_type>::self S(n);

  {
    value_type x(X[0]);
    const typename MatrixRow<const Matrix<G1> >::self A0(A,0);

    int j=0;
    for (int jmax=n-16; j<=jmax; j+=16)
    {
      typename SubArray<typename DenseVector<value_type>::self >::self s(S,j,16);
      typename SubArray<const typename MatrixRow<const Matrix<G1> >::self>::self a(A0,j,16);
      s[ 0] =f(a[ 0],x); s[ 1] =f(a[ 1],x); s[ 2] =f(a[ 2],x); s[ 3] =f(a[ 3],x); s[ 4] =f(a[ 4],x); s[ 5] =f(a[ 5],x); s[ 6] =f(a[ 6],x); s[ 7] =f(a[ 7],x); s[ 8] =f(a[ 8],x); s[ 9] =f(a[ 9],x); s[10] =f(a[10],x); s[11] =f(a[11],x); s[12] =f(a[12],x); s[13] =f(a[13],x); s[14] =f(a[14],x); s[15] =f(a[15],x);
    }
    for (; j<n; ++j) S[j]=f(A0[j],x);
  }

  for (int i=1, imax=m-1; i<imax; ++i)
  {
    value_type x(X[i]);
    const typename MatrixRow<const Matrix<G1> >::self A0(A,i);

    int j=0;
    for (int jmax=n-16; j<=jmax; j+=16)
    {
      typename SubArray<typename DenseVector<value_type>::self >::self s(S,j,16);
      typename SubArray<const typename MatrixRow<const Matrix<G1> >::self>::self a(A0,j,16);
      s[ 0]+=f(a[ 0],x); s[ 1]+=f(a[ 1],x); s[ 2]+=f(a[ 2],x); s[ 3]+=f(a[ 3],x); s[ 4]+=f(a[ 4],x); s[ 5]+=f(a[ 5],x); s[ 6]+=f(a[ 6],x); s[ 7]+=f(a[ 7],x); s[ 8]+=f(a[ 8],x); s[ 9]+=f(a[ 9],x); s[10]+=f(a[10],x); s[11]+=f(a[11],x); s[12]+=f(a[12],x); s[13]+=f(a[13],x); s[14]+=f(a[14],x); s[15]+=f(a[15],x);
    }
    for (; j<n; ++j) S[j]+=f(A0[j],x);
  }

  {
    value_type x(X[m-1]);
    const typename MatrixRow<const Matrix<G1> >::self A0(A,m-1);

    int j=0;
    for (int jmax=n-16; j<=jmax; j+=16)
    {
      typename adaptor_type::sub_type ad=adapt.sub(j,16);
      typename SubArray<typename DenseVector<value_type>::self >::self s(S,j,16);
      typename SubArray<const typename MatrixRow<const Matrix<G1> >::self>::self a(A0,j,16);
      ad( 0,func(s[ 0]+f(a[ 0],x))); ad( 1,func(s[ 1]+f(a[ 1],x))); ad( 2,func(s[ 2]+f(a[ 2],x))); ad( 3,func(s[ 3]+f(a[ 3],x))); ad( 4,func(s[ 4]+f(a[ 4],x))); ad( 5,func(s[ 5]+f(a[ 5],x))); ad( 6,func(s[ 6]+f(a[ 6],x))); ad( 7,func(s[ 7]+f(a[ 7],x))); ad( 8,func(s[ 8]+f(a[ 8],x))); ad( 9,func(s[ 9]+f(a[ 9],x))); ad(10,func(s[10]+f(a[10],x))); ad(11,func(s[11]+f(a[11],x))); ad(12,func(s[12]+f(a[12],x))); ad(13,func(s[13]+f(a[13],x))); ad(14,func(s[14]+f(a[14],x))); ad(15,func(s[15]+f(a[15],x)));
    }
    for (; j<n; ++j) adapt(j,func(S[j]+f(A0[j],x)));
  }
}

template<template<class,class> class Func,class G1,class G2,class G3,class Adapt> inline void mul (const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt) { mul <dot_function,Func>(A,X,Y,adapt); }
template<template<class,class> class Func,class G1,class G2,class G3,class Adapt> inline void tmul(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt) { tmul<dot_function,Func>(A,X,Y,adapt); }


template<int N,template<class,class,class> class F>
struct simd_mulv_function
{
  template<class G1,class G2,class G3,class Adapt>
  inline void operator()(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt) const
  {
    int n=A.ncols();
    if (n%N==0) return mul <F,sum_function>(simd_block(A),simd_block(X),Y,adapt);
    return mul <F,id_function>(A,X,Y,adapt);
  }
};

template<int N,template<class,class,class> class F>
struct simd_tmulv_function
{
  template<class G1,class G2,class G3,class Adapt>
  inline void operator()(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt) const
  {
    int n=A.ncols();
    if (n%N==0) { typename SimdBlock<Vector<G3> >::self Y2=simd_block(Y); return tmul<F,id_function>(simd_block(A),X,Y2,adapt); }
    return tmul<F,id_function>(A,X,Y,adapt);
  }
};

template<template<class,class,class> class F> struct simd_mulv_function <0,F> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt) const { return mul <F,id_function>(A,X,Y,adapt); } };
template<template<class,class,class> class F> struct simd_tmulv_function<0,F> { template<class G1,class G2,class G3,class Adapt> inline void operator()(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt) const { return tmul<F,id_function>(A,X,Y,adapt); } };

template<template<class,class,class> class F,class G1,class G2,class G3,class Adapt,class V> inline void aux_mul (const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt,const V &, const V &) { return simd_mulv_function <SimdLength<V>::RET,F>()(A,X,Y,adapt); }
template<template<class,class,class> class F,class G1,class G2,class G3,class Adapt,class V> inline void aux_tmul(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt,const V &, const V &) { return simd_tmulv_function<SimdLength<V>::RET,F>()(A,X,Y,adapt); }

template<class G1,class G2,class G3,class Adapt> void aux_mul0(const Matrix<G1> &A, const Vector<G2> &X,Vector<G3> &Y, const Adapt &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul1(const Matrix<G1> &A, const Vector<G2> &X,Vector<G3> &Y, const Adapt &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul2(const Matrix<G1> &A, const Vector<G2> &X,Vector<G3> &Y, const Adapt &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }
template<class G1,class G2,class G3,class Adapt> void aux_mul3(const Matrix<G1> &A, const Vector<G2> &X,Vector<G3> &Y, const Adapt &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<int TransA,class G1,class G2,class G3,class Adapt>
void aux_mul(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt)
{
       if (TransA==0) return aux_mul0(A,X,Y,adapt);
  else if (TransA==1) return aux_mul1(A,X,Y,adapt);
  else if (TransA==2) return aux_mul2(A,X,Y,adapt);
  else if (TransA==3) return aux_mul3(A,X,Y,adapt);
}

template<int TransA,class G1,class G2,class G3,class Adapt>
inline void mul(const Matrix<G1> &A,const Vector<G2> &X,Vector<G3> &Y,const Adapt &adapt)
{
  typename ArrayData<Vector<G3> >::self Y2=data(Y);
  return aux_mul<TransA>(data(A),data(X),Y2,adapt);
}

template<int TransA,        class G1,class G2,class G3> inline void mul    (           const Matrix<G1> &A,const Vector<G2> &X,           Vector<G3> &Y) { return mul<TransA>(A,X,Y,id_adapt       (   )); }
template<int TransA,        class G1,class G2,class G3> inline void neg_mul(           const Matrix<G1> &A,const Vector<G2> &X,           Vector<G3> &Y) { return mul<TransA>(A,X,Y,neg_adapt      (   )); }
template<int TransA,class V,class G1,class G2,class G3> inline void mul    (const V &a,const Matrix<G1> &A,const Vector<G2> &X,           Vector<G3> &Y) { return mul<TransA>(A,X,Y,ax_adapt       (  a)); }

template<int TransA,        class G1,class G2,class G3> inline void add_mul(           const Matrix<G1> &A,const Vector<G2> &X,           Vector<G3> &Y) { return mul<TransA>(A,X,Y,add_adapt      (   )); }
template<int TransA,        class G1,class G2,class G3> inline void sub_mul(           const Matrix<G1> &A,const Vector<G2> &X,           Vector<G3> &Y) { return mul<TransA>(A,X,Y,sub_adapt      (   )); }
template<int TransA,class V,class G1,class G2,class G3> inline void add_mul(const V &a,const Matrix<G1> &A,const Vector<G2> &X,           Vector<G3> &Y) { return mul<TransA>(A,X,Y,add_adapt      (  a)); }

template<int TransA,class V,class G1,class G2,class G3> inline void add_mul(           const Matrix<G1> &A,const Vector<G2> &X,const V &s,Vector<G3> &Y) { return mul<TransA>(A,X,Y,scale_add_adapt(s  )); }
template<int TransA,class V,class G1,class G2,class G3> inline void sub_mul(           const Matrix<G1> &A,const Vector<G2> &X,const V &s,Vector<G3> &Y) { return mul<TransA>(A,X,Y,scale_sub_adapt(s  )); }
template<int TransA,class V,class G1,class G2,class G3> inline void add_mul(const V &a,const Matrix<G1> &A,const Vector<G2> &X,const V &s,Vector<G3> &Y) { return mul<TransA>(A,X,Y,scale_add_adapt(s,a)); }


//{unsecret}
//Description: Y=a*mul(A,X)+b*Y
//Remarks:
//  TransA is a template argument to transpose A.
//  (0=No Transpose, 1=Transpose, 2=Conjugate Transpose)
//Example:
//  gemv   (  A,X,  Y);  // Y =  mul(     A ,X), equivalent to gemv<0>(1,A,X,0,Y)
//  gemv   (2,A,X,  Y);  // Y =2*mul(     A ,X), equivalent to gemv<0>(2,A,X,0,Y)
//  gemv   (2,A,X,1,Y);  // Y+=2*mul(     A ,X), equivalent to gemv<0>(2,A,X,1,Y)
//  gemv<1>(  A,X,  Y);  // Y =  mul( trn(A),X), equivalent to gemv<1>(1,A,X,0,Y)
//  gemv<1>(2,A,X,  Y);  // Y =2*mul( trn(A),X), equivalent to gemv<1>(2,A,X,0,Y)
//  gemv<1>(2,A,X,1,Y);  // Y+=2*mul( trn(A),X)
//  gemv<2>(  A,X,  Y);  // Y =  mul(htrn(A),X), equivalent to gemv<2>(1,A,X,0,Y)
//  gemv<2>(2,A,X,  Y);  // Y =2*mul(htrn(A),X), equivalent to gemv<2>(2,A,X,0,Y)
//  gemv<2>(2,A,X,1,Y);  // Y+=2*mul(htrn(A),X)
template<int TransA,class G1,class G2,class G3> void gemv(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Vector<G2> &X, const typename Vector<G3>::value_type &b, Vector<G3> &Y)
{
  typedef typename Matrix<G1>::value_type V;
  if      (a==V( 1)) { if (b==V(0)) mul    <TransA>(  A,X,Y); else if (b==V(1)) add_mul<TransA>(  A,X,Y); else add_mul<TransA>(  A,X,b,Y); }
  else if (a==V(-1)) { if (b==V(0)) neg_mul<TransA>(  A,X,Y); else if (b==V(1)) sub_mul<TransA>(  A,X,Y); else sub_mul<TransA>(  A,X,b,Y); }
  else if (a!=V( 0)) { if (b==V(0)) mul    <TransA>(a,A,X,Y); else if (b==V(1)) add_mul<TransA>(a,A,X,Y); else add_mul<TransA>(a,A,X,b,Y); }
  else scal(b,Y);
}

template<int TransA,class G1,class G2,class G3> void gemv(                                          const Matrix<G1> &A, const Vector<G2> &X, const typename Vector<G3>::value_type &b, Vector<G3> &Y) { typedef typename Matrix<G1>::value_type V; if (b==V(0)) mul<TransA>(A,X,Y); else if (b==V( 1)) add_mul<TransA>(A,X,Y); else add_mul<TransA>(A,X,b,Y); }
template<int TransA,class G1,class G2,class G3> void gemv(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Vector<G2> &X,                                           Vector<G3> &Y) { typedef typename Matrix<G1>::value_type V; if (a==V(1)) mul<TransA>(A,X,Y); else if (a==V(-1)) neg_mul<TransA>(A,X,Y); else if (a!=V(0)) mul<TransA>(a,A,X,Y); else fill(Y,typename Matrix<G3>::value_type(0)); }
template<int TransA,class G1,class G2,class G3> void gemv(                                          const Matrix<G1> &A, const Vector<G2> &X,                                           Vector<G3> &Y) { mul    <TransA>(  A,X,Y); }

template<class G1,class G2,class G3> void gemv(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Vector<G2> &X, const typename Vector<G3>::value_type &b, Vector<G3> &Y) { return gemv<0>(a,A,X,b,Y); }
template<class G1,class G2,class G3> void gemv(const typename Matrix<G1>::value_type &a, const Matrix<G1> &A, const Vector<G2> &X,                                           Vector<G3> &Y) { return gemv<0>(a,A,X,  Y); }
template<class G1,class G2,class G3> void gemv(                                          const Matrix<G1> &A, const Vector<G2> &X, const typename Vector<G3>::value_type &b, Vector<G3> &Y) { return gemv<0>(  A,X,b,Y); }
template<class G1,class G2,class G3> void gemv(                                          const Matrix<G1> &A, const Vector<G2> &X,                                           Vector<G3> &Y) { return gemv<0>(  A,X,  Y); }


template<class V> void gemv(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, V a, const V *pA, int lda, const V *pX, int dx, V b, const V *pY, int dy)
{
  int m2,n2; bool trans;
  if (Order==CblasRowMajor) { m2=m; n2=n; trans=TransA!=CblasNoTrans; } else { m2=n; n2=m; trans=TransA==CblasNoTrans; }
  typename DataMatrix<V>::self A(m2,lda, (V *)pA);
  
  typename DenseVector<V>::self Y2;
  if (!trans)
  {
    typename DenseVector<V>::self X2; if (dx==1); else { X2.resize(n); copy(stride(data_vector(n*dx,pX),dx),X2); pX=&X2[0]; }
    typename DataVector<V>::self X(n,(V *)pX), Y(m,(V *)pY);
    
    if (dy!=1) { Y2.resize(m); Y.generator().data=&Y2[0]; }
    if (TransA!=CblasConjTrans) gemv<0>(a,sub(A,0,0,m2,n2),X,b,Y ); else gemv<3>(a,sub(A,0,0,m2,n2),X,b,Y );
    if (dy!=1) { typename DataVector<V>::self Z(m*dy,(V *)pY); typename StrideArray<typename DataVector<V>::self>::self Z2=stride(Z,dy); Z2=Y; }
  }
  else
  {
    typename DenseVector<V>::self X2; if (dx==1); else { X2.resize(n); copy(stride(data_vector(n*dx,pX),dx),X2); pX=&X2[0]; }
    typename DataVector<V>::self X(n,(V *)pX), Y(m*dy,(V *)pY);
    
    if (dy!=1) { Y2.resize(m); Y.generator().data=&Y2[0]; }
    if (TransA!=CblasConjTrans) gemv<1>(a,sub(A,0,0,m2,n2),X,b,Y ); else  gemv<2>(a,sub(A,0,0,m2,n2),X,b,Y );
    if (dy!=1) { typename DataVector<V>::self Z(m*dy,(V *)pY); typename StrideArray<typename DataVector<V>::self>::self Z2=stride(Z,dy); Z2=Y; }
  }
}

#ifdef AUTO_BLAS
template<class G,class T1,class T2> inline void       affect    (Vector<G> &Y, const Vector<multiply_matrix_vector_generator<T1,T2> > &X) { mul    (X.generator().second_array(),X.generator().first_array(),Y); }
template<class G,class T1,class T2> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<multiply_matrix_vector_generator<T1,T2> > &X) { add_mul(X.generator().second_array(),X.generator().first_array(),Y); return Y; }
template<class G,class T1,class T2> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<multiply_matrix_vector_generator<T1,T2> > &X) { sub_mul(X.generator().second_array(),X.generator().first_array(),Y); return Y; }

template<class G,class T1,class T2,class Arg,          class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,neg_function    <Arg,    Res> > > &X) { neg_mul(                                 X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,val_mul_function<Arg,Val,Res> > > &X) { mul    (X.generator().function().value(),X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,mul_val_function<Arg,Val,Res> > > &X) { mul    (X.generator().function().value(),X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,val_mul_function<Arg,Val,Res> > > &X) { add_mul(X.generator().function().value(),X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,val_mul_function<Arg,Val,Res> > > &X) { sub_mul(X.generator().function().value(),X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,mul_val_function<Arg,Val,Res> > > &X) { add_mul(X.generator().function().value(),X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<T1,T2> >,mul_val_function<Arg,Val,Res> > > &X) { sub_mul(X.generator().function().value(),X.generator().array().generator().second_array(),X.generator().array().generator().first_array(),Y); return Y; }


template<class G,class T1,class T2> inline void       affect    (Vector<G> &Y, const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> > &X) { tmul(X.generator().second_array().generator().array(),X.generator().first_array(),Y); }
template<class G,class T1,class T2> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> > &X) { tmul(X.generator().second_array().generator().array(),X.generator().first_array(),Y); return Y; }
template<class G,class T1,class T2> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> > &X) { tmul(X.generator().second_array().generator().array(),X.generator().first_array(),Y); return Y; }

template<class G,class T1,class T2,class Arg,          class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,neg_function    <Arg,    Res> > > &X) { neg_tmul(                                 X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,val_mul_function<Arg,Val,Res> > > &X) { tmul    (X.generator().function().value(),X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,mul_val_function<Arg,Val,Res> > > &X) { tmul    (X.generator().function().value(),X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,val_mul_function<Arg,Val,Res> > > &X) { add_tmul(X.generator().function().value(),X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,val_mul_function<Arg,Val,Res> > > &X) { sub_tmul(X.generator().function().value(),X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,mul_val_function<Arg,Val,Res> > > &X) { add_tmul(X.generator().function().value(),X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); return Y; }
template<class G,class T1,class T2,class Arg,class Val,class Res> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<function_array_generator<const Vector<multiply_matrix_vector_generator<const Matrix<trans_matrix_generator<T1> >,T2> >,mul_val_function<Arg,Val,Res> > > &X) { sub_tmul(X.generator().function().value(),X.generator().array().generator().second_array().generator().array(),X.generator().array().generator().first_array(),Y); return Y; }
#endif


#ifdef BLAS_PRECOMPILE

#define DECLARE_BLAS_GEMV(V) \
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt);\
GENIAL_API void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt);

#define DEFINE_BLAS_GEMV(V) \
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul0(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul < dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul1(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_tmul< dot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul2(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_tmul<cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const neg_adaptor_function                               &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<id_adaptor_function    ,V > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function<V >   > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<neg_adaptor_function   ,V > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }\
void aux_mul3(const Matrix<data_matrix_generator<V > > &A,const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const scale_adaptor_function<ax_adaptor_function<V >,V > &adapt) { return aux_mul <cdot_function>(A,X,Y,adapt,V(),V()); }

DECLARE_BLAS_GEMV(float)
DECLARE_BLAS_GEMV(double)
DECLARE_BLAS_GEMV(complex<float>)
DECLARE_BLAS_GEMV(complex<double>)

#endif //BLAS_PRECOMPILE

#endif //__cplusplus



#ifdef BLAS_PRECOMPILE
#ifdef __cplusplus
extern "C" {
#endif

//Group = CBLAS Level 2

//{unsecret}
GENIAL_API void BLAS_NAME(sgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const float  a, const float  *pA, int lda, const float  *pX, int dx, const float  b, float  *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(dgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const double a, const double *pA, int lda, const double *pX, int dx, const double b, double *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(cgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const void  *a, const void   *pA, int lda, const void   *pX, int dx, const void  *b, void   *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(zgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const void  *a, const void   *pA, int lda, const void   *pX, int dx, const void  *b, void   *pY, int dy);

#ifdef __cplusplus
}
#endif
#endif //BLAS_PRECOMPILE

#endif //GEMV_H
