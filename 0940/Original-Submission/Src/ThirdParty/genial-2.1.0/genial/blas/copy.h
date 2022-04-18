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


#ifndef COPY_H
#define COPY_H

#ifdef __cplusplus

#include "blas.h"


//Group = Linear Algebra

template<template<class,class> class F,class G1,class G2,class Adapt>
void copy(const Vector<G1> &X, Vector<G2> &Y, const Adapt &adaptor)
{
  F<typename Vector<G1>::value_type,typename Vector<G2>::value_type> f;
  typedef typename Adapt::template result_rebind<Vector<G2> >::other adaptor_type;
  adaptor_type adapt = adaptor(Y);
  int n=adapt.size();
  
  int i=0;
  for (int imax=n-16; i<=imax; i+=16)
  {
    typename SubArray<const Vector<G1> >::self x=sub(X,i,16);
    typename adaptor_type::sub_type ad=adapt.sub(i,16);
    ad( 0,f(x[ 0])); ad( 1,f(x[ 1])); ad( 2,f(x[ 2])); ad( 3,f(x[ 3])); ad( 4,f(x[ 4])); ad( 5,f(x[ 5])); ad( 6,f(x[ 6])); ad( 7,f(x[ 7])); ad( 8,f(x[ 8])); ad( 9,f(x[ 9])); ad(10,f(x[10])); ad(11,f(x[11])); ad(12,f(x[12])); ad(13,f(x[13])); ad(14,f(x[14])); ad(15,f(x[15]));
  }
  for (; i<n; ++i)
    adapt(i,f(X[i]));
}
template<template<class,class> class F,class G1,class G2,class Adapt>
void copy(const Matrix<G1> &X, Matrix<G2> &Y, const Adapt &adapt)
{
  for (int i=0, m=X.nrows(); i<m; ++i)
  {
    typename MatrixRow<Matrix<G2> >::self Bi = row(Y,i);
    copy<F>(row(X,i),Bi,adapt);
  }
}

template<template<class,class> class F,class G1,class G2> inline void copy(const Vector<G1> &X, Vector<G2> &Y) { return copy<F>(X,Y,id_adapt()); }
template<template<class,class> class F,class G1,class G2> inline void copy(const Matrix<G1> &X, Matrix<G2> &Y) { return copy<F>(X,Y,id_adapt()); }

template<int N> 
struct simd_copy_function    
{
  template<class G1,class G2,class Adapt> inline void operator()(const Vector<G1> &X,Vector<G2> &Y, const Adapt &adapt) const { int n=X.size (); if (n%N) return copy<id_function>(X,Y,adapt); typename SimdBlock<Vector<G2> >::self B2(Y); return copy<id_function>(simd_block(X),B2,adapt); } 
  template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &X,Matrix<G2> &Y, const Adapt &adapt) const { int n=X.ncols(); if (n%N) return copy<id_function>(X,Y,adapt); typename SimdBlock<Matrix<G2> >::self B2(Y); return copy<id_function>(simd_block(X),B2,adapt); } 
};
template<> struct simd_copy_function<0> 
{ 
  template<class G1,class G2,class Adapt> inline void operator()(const Vector<G1> &X,Vector<G2> &Y, const Adapt &adapt) const { return copy<id_function>(X,Y,adapt); }
  template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &X,Matrix<G2> &Y, const Adapt &adapt) const { return copy<id_function>(X,Y,adapt); }
};

template<class V          ,class Adapt> inline void aux_copy(const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<V> > &Y, const Adapt &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }
template<class V          ,class Adapt> inline void aux_copy(const Matrix<data_matrix_generator<V> > &X, Matrix<data_matrix_generator<V> > &Y, const Adapt &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }
template<class G1,class G2,class Adapt> inline void aux_copy(const Vector<G1>                        &X, Vector<G2>                        &Y, const Adapt &adapt) { return simd_copy_function<0                       >()(X,Y,adapt); }
template<class G1,class G2,class Adapt> inline void aux_copy(const Matrix<G1>                        &X, Matrix<G2>                        &Y, const Adapt &adapt) { return simd_copy_function<0                       >()(X,Y,adapt); }

template<class G1,class G2,class Adapt> inline void _copy(const Vector<G1> &X, Vector<G2> &Y, const Adapt &adapt) { typename ArrayData<Vector<G2> >::self B2=data(Y); return aux_copy(data(X),B2,adapt); }
template<class G1,class G2,class Adapt> inline void _copy(const Matrix<G1> &X, Matrix<G2> &Y, const Adapt &adapt) { typename ArrayData<Matrix<G2> >::self B2=data(Y); return aux_copy(data(X),B2,adapt); }

//{unsecret}
//Summary: Y=X
template<class G1,class G2> inline void copy(const Vector<G1> &X, Vector<G2> &Y) { return _copy(X,Y,id_adapt()); }
//{unsecret}
template<class G1,class G2> inline void copy(const Matrix<G1> &X, Matrix<G2> &Y) { return _copy(X,Y,id_adapt()); }

template<class V,class G1,class G2> inline void copy    (const V &a,const Vector<G1> &X, Vector<G2> &Y) { return _copy(X,Y,ax_adapt <typename Vector<G2>::value_type>(a)); }
template<class V,class G1,class G2> inline void copy    (const V &a,const Matrix<G1> &X, Matrix<G2> &Y) { return _copy(X,Y,ax_adapt <typename Matrix<G2>::value_type>(a)); }
template<        class G1,class G2> inline void add_copy(           const Vector<G1> &X, Vector<G2> &Y) { return _copy(X,Y,add_adapt()); }
template<        class G1,class G2> inline void add_copy(           const Matrix<G1> &X, Matrix<G2> &Y) { return _copy(X,Y,add_adapt()); }
template<class V,class G1,class G2> inline void add_copy(const V &a,const Vector<G1> &X, Vector<G2> &Y) { return _copy(X,Y,add_adapt<typename Vector<G2>::value_type>(a)); }
template<class V,class G1,class G2> inline void add_copy(const V &a,const Matrix<G1> &X, Matrix<G2> &Y) { return _copy(X,Y,add_adapt<typename Vector<G2>::value_type>(a)); }
template<        class G1,class G2> inline void sub_copy(           const Vector<G1> &X, Vector<G2> &Y) { return _copy(X,Y,sub_adapt()); }
template<        class G1,class G2> inline void sub_copy(           const Matrix<G1> &X, Matrix<G2> &Y) { return _copy(X,Y,sub_adapt()); }
template<class V,class G1,class G2> inline void sub_copy(const V &a,const Vector<G1> &X, Vector<G2> &Y) { return _copy(X,Y,sub_adapt<typename Vector<G2>::value_type>(a)); }
template<class V,class G1,class G2> inline void sub_copy(const V &a,const Matrix<G1> &X, Matrix<G2> &Y) { return _copy(X,Y,sub_adapt<typename Vector<G2>::value_type>(a)); }




template<template<class,class> class F,class G,class Adapt,class V>
void fill(Vector<G> &X, const Adapt &adaptor, const V &val)
{
  F<V,typename Vector<G>::value_type> f;
  typedef typename Adapt::template result_rebind<Vector<G> >::other adaptor_type;
  adaptor_type adapt = adaptor(X);
  int n=adapt.size();
  
  int i=0;
  for (int imax=n-16; i<=imax; i+=16)
  {  
    typename adaptor_type::sub_type ad=adapt.sub(i,16);
    ad( 0,f(val)); ad( 1,f(val)); ad( 2,f(val)); ad( 3,f(val)); ad( 4,f(val)); ad( 5,f(val)); ad( 6,f(val)); ad( 7,f(val)); ad( 8,f(val)); ad( 9,f(val)); ad(10,f(val)); ad(11,f(val)); ad(12,f(val)); ad(13,f(val)); ad(14,f(val)); ad(15,f(val));
  }  
  for (; i<n; ++i) adapt(i,f(val));
}

template<template<class,class> class F,class G,class Adapt,class V>
void fill(Matrix<G> &X, const Adapt &adaptor, const V &val)
{
  for (int i=0, m=X.nrows(); i<m; ++i)
  {
    typename MatrixRow<Matrix<G> >::self Ai = row(X,i);
    fill<F>(Ai,adaptor,val);
  }
}

template<template<class,class> class F,class G> inline void fill(Vector<G> &X, const typename Vector<G>::value_type &val) { return fill<F>(X,id_adapt(),val); }
template<template<class,class> class F,class G> inline void fill(Matrix<G> &X, const typename Matrix<G>::value_type &val) { return fill<F>(X,id_adapt(),val); }

template<int N> 
struct simd_fill_function    
{
  template<class G,class Adapt> inline void operator()(Vector<G> &X, const Adapt &adapt, const typename Vector<G>::value_type &val) const { int n=X.size (); if (n%N) return fill<id_function>(X,adapt,val); typename SimdBlock<Vector<G> >::self A2(X); return fill<id_function>(A2,adapt,typename SimdBlock<Matrix<G> >::self::value_type(val)); }
  template<class G,class Adapt> inline void operator()(Matrix<G> &X, const Adapt &adapt, const typename Matrix<G>::value_type &val) const { int n=X.ncols(); if (n%N) return fill<id_function>(X,adapt,val); typename SimdBlock<Matrix<G> >::self A2(X); return fill<id_function>(A2,adapt,typename SimdBlock<Matrix<G> >::self::value_type(val)); } 
};
template<> struct simd_fill_function<0> 
{
  template<class G,class Adapt> inline void operator()(Vector<G> &X, const Adapt &adapt, const typename Vector<G>::value_type &val) const { return fill<id_function>(X,adapt,val); }
  template<class G,class Adapt> inline void operator()(Matrix<G> &X, const Adapt &adapt, const typename Matrix<G>::value_type &val) const { return fill<id_function>(X,adapt,val); }
};

template<class V,class Adapt        > inline void aux_fill(Vector<data_vector_generator<V> > &X, const Adapt &adapt, const V &val) { return simd_fill_function<SimdLength<V >::RET>()(X,adapt,val); }
template<class V,class Adapt        > inline void aux_fill(Matrix<data_matrix_generator<V> > &X, const Adapt &adapt, const V &val) { return simd_fill_function<SimdLength<V >::RET>()(X,adapt,val); }
template<class G,class Adapt,class V> inline void aux_fill(Vector<G>                         &X, const Adapt &adapt, const V &val) { return simd_fill_function<0                       >()(X,adapt,val); }
template<class G,class Adapt,class V> inline void aux_fill(Matrix<G>                         &X, const Adapt &adapt, const V &val) { return simd_fill_function<0                       >()(X,adapt,val); }

template<class G,class Adapt> inline void fill(Vector<G> &X, const Adapt &adapt, const typename Vector<G>::value_type &val) { typename ArrayData<Vector<G> >::self A2=data(X); return aux_fill(A2,adapt,val); }
template<class G,class Adapt> inline void fill(Matrix<G> &X, const Adapt &adapt, const typename Matrix<G>::value_type &val) { typename ArrayData<Matrix<G> >::self A2=data(X); return aux_fill(A2,adapt,val); }

template<class G> void fill(Vector<G> &X, const typename Vector<G>::value_type &val) { fill(X,id_adapt(),val); }
template<class G> void fill(Matrix<G> &X, const typename Matrix<G>::value_type &val) { fill(X,id_adapt(),val); }




template<template<class,class> class F,class G1,class G2,class Adapt>
void fill_copy(const Vector<G1> &A, Vector<G2> &B, const Adapt &adaptor, const typename Vector<G1>::value_type &val)
{
  F<typename Vector<G1>::value_type,typename Vector<G2>::value_type> f;
  typedef typename Adapt::template result_rebind<Vector<G2> >::other adaptor_type;
  adaptor_type adapt = adaptor(B);
  int n1=A.size(), n2=adapt.size();  assert(n1<=n2);
  
  int i=0;
  {
    for (int imax=n1-16; i<=imax; i+=16)
    {
      typename SubArray<const Vector<G1> >::self x=sub(A,i,16);
      typename adaptor_type::sub_type ad=adapt.sub(i,16);
      ad( 0,f(x[ 0])); ad( 1,f(x[ 1])); ad( 2,f(x[ 2])); ad( 3,f(x[ 3])); ad( 4,f(x[ 4])); ad( 5,f(x[ 5])); ad( 6,f(x[ 6])); ad( 7,f(x[ 7])); ad( 8,f(x[ 8])); ad( 9,f(x[ 9])); ad(10,f(x[10])); ad(11,f(x[11])); ad(12,f(x[12])); ad(13,f(x[13])); ad(14,f(x[14])); ad(15,f(x[15]));
    }
    for (; i<n1; ++i)
      adapt(i,f(A[i]));
  }
  
  {  
    for (int imax=n2-16; i<=imax; i+=16)
    {
      typename adaptor_type::sub_type ad=adapt.sub(i,16);
      ad( 0,f(val)); ad( 1,f(val)); ad( 2,f(val)); ad( 3,f(val)); ad( 4,f(val)); ad( 5,f(val)); ad( 6,f(val)); ad( 7,f(val)); ad( 8,f(val)); ad( 9,f(val)); ad(10,f(val)); ad(11,f(val)); ad(12,f(val)); ad(13,f(val)); ad(14,f(val)); ad(15,f(val));
    }
    for (; i<n2; ++i)
      adapt(i,f(val));
  }
}

template<template<class,class> class F,class G1,class G2,class Adapt>
void fill_copy(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor, const typename Matrix<G1>::value_type &val)
{
  int n1=A.nrows(), n2=B.nrows(); assert(n1<=n2);
  for (int i=0 ; i<n1; ++i) { typename MatrixRow<Matrix<G2> >::self Bi=row(B,i); fill_copy<F>(row(A,i),Bi,adaptor,val); }
  for (int i=n1; i<n2; ++i) { typename MatrixRow<Matrix<G2> >::self Bi=row(B,i); fill     <F>(         Bi,adaptor,val); }
}

template<template<class,class>class F,class G1,class G2> inline void fill_copy(const Vector<G1> &A, Vector<G2> &B, const typename Vector<G1>::value_type &val) { return fill_copy<F>(A,B,id_adapt(),val); }
template<template<class,class>class F,class G1,class G2> inline void fill_copy(const Matrix<G1> &A, Matrix<G2> &B, const typename Matrix<G1>::value_type &val) { return fill_copy<F>(A,B,id_adapt(),val); }

template<int N> 
struct simd_fill_copy_function    
{
  template<class G1,class G2,class Adapt> inline void operator()(const Vector<G1> &A,Vector<G2> &B, const Adapt &adapt1, const typename Vector<G1>::value_type &val) const { int n=A.size (); if (n%N) return fill_copy<id_function>(A,B,adapt1,val); typename SimdBlock<Vector<G2> >::self B2(B); return fill_copy<id_function>(simd_block(A),B2,adapt1,val); }
  template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt1, const typename Matrix<G1>::value_type &val) const { int n=A.ncols(); if (n%N) return fill_copy<id_function>(A,B,adapt1,val); typename SimdBlock<Matrix<G2> >::self B2(B); return fill_copy<id_function>(simd_block(A),B2,adapt1,val); } 
};
template<> struct simd_fill_copy_function<0> 
{ 
  template<class G1,class G2,class Adapt> inline void operator()(const Vector<G1> &A,Vector<G2> &B, const Adapt &adapt1, const typename Vector<G1>::value_type &val) const { return fill_copy<id_function>(A,B,adapt1,val); } 
  template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt1, const typename Matrix<G1>::value_type &val) const { return fill_copy<id_function>(A,B,adapt1,val); } 
};

template<class G1,class G2,class Adapt,class V1,class V2> inline void aux_fill_copy(const Vector<G1> &A,Vector<G2> &B, const Adapt &adapt1, const typename Vector<G1>::value_type &val, const V1 &,const V2 &) { return simd_fill_copy_function<0                       >()(A,B,adapt1,val); } 
template<class G1,class G2,class Adapt,class V1,class V2> inline void aux_fill_copy(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt1, const typename Matrix<G1>::value_type &val, const V1 &,const V2 &) { return simd_fill_copy_function<0                       >()(A,B,adapt1,val); } 
template<class G1,class G2,class Adapt,class V1         > inline void aux_fill_copy(const Vector<G1> &A,Vector<G2> &B, const Adapt &adapt1, const typename Vector<G1>::value_type &val, const V1 &,const V1 &) { return simd_fill_copy_function<SimdLength<V1>::RET>()(A,B,adapt1,val); }
template<class G1,class G2,class Adapt,class V1         > inline void aux_fill_copy(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt1, const typename Matrix<G1>::value_type &val, const V1 &,const V1 &) { return simd_fill_copy_function<SimdLength<V1>::RET>()(A,B,adapt1,val); } 

template<class G1,class G2,class Adapt> inline void fill_copy(const Vector<G1> &A, Vector<G2> &B, const Adapt &adapt1, const typename Vector<G1>::value_type &val) { typename ArrayData<Vector<G2> >::self B2=data(B); return aux_fill_copy(data(A),B2,adapt1,val,typename Vector<G1>::value_type(),typename Vector<G2>::value_type()); }
template<class G1,class G2,class Adapt> inline void fill_copy(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt1, const typename Matrix<G1>::value_type &val) { typename ArrayData<Matrix<G2> >::self B2=data(B); return aux_fill_copy(data(A),B2,adapt1,val,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<class G1,class G2> inline void fill_copy(const Vector<G1> &A, Vector<G2> &B, const typename Vector<G1>::value_type &val) { return fill_copy(A,B,id_adapt(),val); }
template<class G1,class G2> inline void fill_copy(const Matrix<G1> &A, Matrix<G2> &B, const typename Matrix<G1>::value_type &val) { return fill_copy(A,B,id_adapt(),val); }






template<template<class,class> class F,class G1,class G2,class Adapt>
void trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  for (int i=0, m=A.ncols(); i<m; ++i)
  {
    typename MatrixRow<Matrix<G2> >::self Bi = row(B,i);
    copy<F>(col(A,i),Bi,adaptor);
  }
}

template<template<class,class> class F,class G1,class G2> inline void trn(const Matrix<G1> &A, Matrix<G2> &B) { return trn<F>(A,B,id_adapt()); }

template<template<class,class> class F,class G1,class G2,class Adapt>
void trn2(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  
  int m=2*A.ncols(), n=A.nrows();
    
  for (int i=0; i<m; i+=2)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    
    col_type Ai=col(A,i/2);
    row_type B0=row(B,i),B1=row(B,i+1);
    adaptor_type adapt0=adaptor(B0),adapt1=adaptor(B1);
    
    typename Vector<G1>::value_type a,b;
    int j=0;
    for (int jmax=n/2; j<jmax; ++j)
    {
      typename SubArray<col_type>::self x=sub(Ai,2*j,2);
      a=x[ 0]; b=x[ 1]; trn(a,b); adapt0(j,f(a)); adapt1(j,f(b)); 
    }    
  }
}

template<template<class,class> class F,class G1,class G2,class Adapt>
void trn4(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  
  int m=4*A.ncols(), n=A.nrows();
    
  for (int i=0; i<m; i+=4)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    
    col_type Ai=col(A,i/4);
    row_type B0=row(B,i),B1=row(B,i+1),B2=row(B,i+2),B3=row(B,i+3);
    adaptor_type adapt0=adaptor(B0),adapt1=adaptor(B1),adapt2=adaptor(B2),adapt3=adaptor(B3);
    
    typename Vector<G1>::value_type a,b,c,d;
    int j=0;
    for (int jmax=n/4; j<jmax; ++j)
    {
      typename SubArray<col_type>::self x=sub(Ai,4*j,4);
      a=x[ 0]; b=x[ 1]; c=x[ 2]; d=x[ 3]; trn(a,b,c,d); adapt0(j,f(a)); adapt1(j,f(b)); adapt2(j,f(c)); adapt3(j,f(d)); 
    }    
  }
}

template<int N,template<class,class> class F> struct trn_function      { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return trn <F>(A,B,adapt); } };
template<      template<class,class> class F> struct trn_function<2,F> { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return trn2<F>(A,B,adapt); } };
template<      template<class,class> class F> struct trn_function<4,F> { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return trn4<F>(A,B,adapt); } };


template<int N,template<class,class> class F> 
struct simd_trn_function    
{
  template<class G1,class G2,class Adapt> 
  inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt) const 
  { 
    int m=A.nrows(), n=A.ncols(); 
    if (n%N==0 && m%N==0)
    {
      typename SimdBlock<Matrix<G2> >::self B2(B); 
      return trn_function<N,F>()(simd_block(A),B2,adapt); 
    }
    return trn<F>(A,B,adapt); 
  }
};
template<template<class,class> class F> 
struct simd_trn_function<0,F> 
{
  template<class G1,class G2,class Adapt> 
  inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt) const 
  { 
    return trn<F>(A,B,adapt); 
  }
};

template<class G1,class G2,class Adapt,class V1,class V2> inline void aux_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const V1 &,const V2 &) { return simd_trn_function<0                  ,id_function>()(A,B,adapt); } 
template<class G1,class G2,class Adapt,class V1         > inline void aux_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const V1 &,const V1 &) { return simd_trn_function<SimdLength<V1>::RET,id_function>()(A,B,adapt); }

template<class G1,class G2,class Adapt> inline void trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { typename ArrayData<Matrix<G2> >::self B2=data(B); return aux_trn(data(A),B2,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<class G1,class G2> void trn(const Matrix<G1> &A, Matrix<G2> &B) { trn(A,B,id_adapt()); }




template<template<class,class> class F,class G1,class G2,class Adapt>
void up_trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  int m=A.ncols(), n=A.nrows();
  for (int i=0; i<m; ++i)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    col_type Ai=col(A,i);
    row_type Bi=row(B,i);
    adaptor_type adapt=adaptor(Bi);
    for (int j=0; j<=i; ++j)
      adapt(j,f(Ai[j]));
  }
}

template<template<class,class> class F,class G1,class G2> inline void up_trn(const Matrix<G1> &A, Matrix<G2> &B) { return up_trn<F>(A,B,id_adapt()); }

template<template<class,class> class F,class G1,class G2,class Adapt>
void up_trn2(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  int m=A.ncols(), n=A.nrows();
  for (int i=0; i<m; ++i)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    
    col_type Ai=col(A,i);
    row_type B0=row(B,2*i),B1=row(B,2*i+1);
    adaptor_type adapt0=adaptor(B0),adapt1=adaptor(B1);
    
    for (int j=0; j<i; ++j)
    {
      typename SubArray<col_type>::self x=sub(Ai,2*j,2);
      typename Vector<G1>::value_type a=x[0], b=x[1]; trn(a,b); adapt0(j,f(a)); adapt1(j,f(b)); 
    }
    {
      typename SubArray<col_type>::self x=sub(Ai,2*i,2);
      typename Vector<G1>::value_type a=x[0], b=x[1]; up_trn(a,b); adapt0(i,f(a)); adapt1(i,f(b)); 
    }
  }
}

template<template<class,class> class F,class G1,class G2,class Adapt>
void up_trn4(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  int m=A.ncols(), n=A.nrows();
    
  for (int i=0; i<m; ++i)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    
    col_type Ai=col(A,i);
    row_type B0=row(B,4*i),B1=row(B,4*i+1),B2=row(B,4*i+2),B3=row(B,4*i+3);
    adaptor_type adapt0=adaptor(B0),adapt1=adaptor(B1),adapt2=adaptor(B2),adapt3=adaptor(B3);
    
    for (int j=0; j<i; ++j)
    {
      typename SubArray<col_type>::self x=sub(Ai,4*j,4);
      typename Vector<G1>::value_type a=x[0], b=x[1], c=x[2], d=x[3]; trn(a,b,c,d); adapt0(j,f(a)); adapt1(j,f(b)); adapt2(j,f(c)); adapt3(j,f(d)); 
    }    
    {
      typename SubArray<col_type>::self x=sub(Ai,4*i,4);
      typename Vector<G1>::value_type a=x[0], b=x[1], c=x[2], d=x[3]; up_trn(a,b,c,d); adapt0(i,f(a)); adapt1(i,f(b)); adapt2(i,f(c)); adapt3(i,f(d)); 
    }    
  }
}

template<int N,template<class,class> class F> struct up_trn_function      { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return up_trn <F>(A,B,adapt); } };
template<      template<class,class> class F> struct up_trn_function<2,F> { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return up_trn2<F>(A,B,adapt); } };
template<      template<class,class> class F> struct up_trn_function<4,F> { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return up_trn4<F>(A,B,adapt); } };


template<int N,template<class,class> class F> 
struct simd_up_trn_function    
{
  template<class G1,class G2,class Adapt> 
  inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt) const 
  { 
    int m=A.nrows(), n=A.ncols(); 
    if (n%N==0 && m%N==0)
    {
      typename SimdBlock<Matrix<G2> >::self B2(B); 
      return up_trn_function<N,F>()(simd_block(A),B2,adapt); 
    }
    return up_trn<F>(A,B,adapt); 
  }
};
template<template<class,class> class F> 
struct simd_up_trn_function<0,F> 
{
  template<class G1,class G2,class Adapt> 
  inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt) const 
  { 
    return up_trn<F>(A,B,adapt); 
  }
};

template<class G1,class G2,class Adapt,class V1,class V2> inline void aux_up_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const V1 &,const V2 &) { return simd_up_trn_function<0                  ,id_function>()(A,B,adapt); } 
template<class G1,class G2,class Adapt,class V1         > inline void aux_up_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const V1 &,const V1 &) { return simd_up_trn_function<SimdLength<V1>::RET,id_function>()(A,B,adapt); }

template<class G1,class G2,class Adapt> inline void up_trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { typename ArrayData<Matrix<G2> >::self B2=data(B); return aux_up_trn(data(A),B2,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<class G1,class G2> void up_trn(const Matrix<G1> &A, Matrix<G2> &B) { up_trn(A,B,id_adapt()); }






template<template<class,class> class F,class G1,class G2,class Adapt>
void low_trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  int m=A.ncols(), n=A.nrows();
  for (int i=0; i<m; ++i)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    col_type Ai=col(A,i);
    row_type Bi=row(B,i);
    adaptor_type adapt=adaptor(Bi);
    for (int j=i; j<m; ++j)
      adapt(j,f(Ai[j]));
  }
}

template<template<class,class> class F,class G1,class G2> inline void low_trn(const Matrix<G1> &A, Matrix<G2> &B) { return low_trn<F>(A,B,id_adapt()); }

template<template<class,class> class F,class G1,class G2,class Adapt>
void low_trn2(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  int m=A.ncols(), n=A.nrows();
  for (int i=0; i<m; ++i)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    
    col_type Ai=col(A,i);
    row_type B0=row(B,2*i),B1=row(B,2*i+1);
    adaptor_type adapt0=adaptor(B0),adapt1=adaptor(B1);
    
    {
      typename SubArray<col_type>::self x=sub(Ai,2*i,2);
      typename Vector<G1>::value_type a=x[0], b=x[1]; low_trn(a,b); adapt0(i,f(a)); adapt1(i,f(b)); 
    }
    for (int j=i+1; j<m; ++j)
    {
      typename SubArray<col_type>::self x=sub(Ai,2*j,2);
      typename Vector<G1>::value_type a=x[0], b=x[1]; trn(a,b); adapt0(j,f(a)); adapt1(j,f(b)); 
    }
  }
}

template<template<class,class> class F,class G1,class G2,class Adapt>
void low_trn4(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor)
{
  F<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type> f;
  int m=A.ncols(), n=A.nrows();
    
  for (int i=0; i<m; ++i)
  {
    typedef typename MatrixCol<const Matrix<G1> >::self col_type;
    typedef typename MatrixRow<      Matrix<G2> >::self row_type;
    typedef typename Adapt::template result_rebind<row_type>::other adaptor_type;
    
    col_type Ai=col(A,i);
    row_type B0=row(B,4*i),B1=row(B,4*i+1),B2=row(B,4*i+2),B3=row(B,4*i+3);
    adaptor_type adapt0=adaptor(B0),adapt1=adaptor(B1),adapt2=adaptor(B2),adapt3=adaptor(B3);
    
    {
      typename SubArray<col_type>::self x=sub(Ai,4*i,4);
      typename Vector<G1>::value_type a=x[0], b=x[1], c=x[2], d=x[3]; low_trn(a,b,c,d); adapt0(i,f(a)); adapt1(i,f(b)); adapt2(i,f(c)); adapt3(i,f(d)); 
    }    
    for (int j=i+1; j<m; ++j)
    {
      typename SubArray<col_type>::self x=sub(Ai,4*j,4);
      typename Vector<G1>::value_type a=x[0], b=x[1], c=x[2], d=x[3]; trn(a,b,c,d); adapt0(j,f(a)); adapt1(j,f(b)); adapt2(j,f(c)); adapt3(j,f(d)); 
    }    
  }
}

template<int N,template<class,class> class F> struct low_trn_function      { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return low_trn <F>(A,B,adapt); } };
template<      template<class,class> class F> struct low_trn_function<2,F> { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return low_trn2<F>(A,B,adapt); } };
template<      template<class,class> class F> struct low_trn_function<4,F> { template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { return low_trn4<F>(A,B,adapt); } };

template<int N,template<class,class> class F> 
struct simd_low_trn_function    
{
  template<class G1,class G2,class Adapt> 
  inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt) const 
  { 
    int m=A.nrows(), n=A.ncols(); 
    if (n%N==0 && m%N==0)
    {
      typename SimdBlock<Matrix<G2> >::self B2(B); 
      return low_trn_function<N,F>()(simd_block(A),B2,adapt); 
    }
    return low_trn<F>(A,B,adapt); 
  }
};
template<template<class,class> class F> 
struct simd_low_trn_function<0,F> 
{
  template<class G1,class G2,class Adapt> 
  inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt) const 
  { 
    return low_trn<F>(A,B,adapt); 
  }
};

template<class G1,class G2,class Adapt,class V1,class V2> inline void aux_low_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const V1 &,const V2 &) { return simd_low_trn_function<0                  ,id_function>()(A,B,adapt); } 
template<class G1,class G2,class Adapt,class V1         > inline void aux_low_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const V1 &,const V1 &) { return simd_low_trn_function<SimdLength<V1>::RET,id_function>()(A,B,adapt); }

template<class G1,class G2,class Adapt> inline void low_trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt) { typename ArrayData<Matrix<G2> >::self B2=data(B); return aux_low_trn(data(A),B2,adapt,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<class G1,class G2> void low_trn(const Matrix<G1> &A, Matrix<G2> &B) { low_trn(A,B,id_adapt()); }





template<template<class,class> class F,class G1,class G2,class Adapt>
void fill_trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adaptor,const typename Matrix<G1>::value_type &val)
{
  int m1=A.ncols(), m2=B.nrows(); assert(m1<=m2);
  for (int i= 0; i<m1; ++i) { typename MatrixRow<Matrix<G2> >::self Bi=row(B,i); fill_copy<F>(col(A,i),Bi,adaptor,val); }
  for (int i=m1; i<m2; ++i) { typename MatrixRow<Matrix<G2> >::self Bi=row(B,i); fill     <F>(         Bi,adaptor,val); }
}

template<template<class,class>class F,class G1,class G2> inline void fill_trn(const Matrix<G1> &A, Matrix<G2> &B, const typename Matrix<G2>::value_type &val) { return fill_trn<F>(A,B,id_adapt(),val); }

template<int N> 
struct simd_fill_trn_function    
{
  template<class G1,class G2,class Adapt> inline void operator()(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt,const typename Matrix<G1>::value_type &val) const 
  { 
    return fill_trn<id_function>(A,B,adapt,val); 
  } 
};

template<class G1,class G2,class Adapt,class V1,class V2> inline void aux_fill_trn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const typename Matrix<G1>::value_type &val, const V1 &,const V2 &) { return simd_fill_trn_function<0                       >()(A,B,adapt,val); } 
template<class G1,class G2,class Adapt,class V1         > inline void aux_fill_tn(const Matrix<G1> &A,Matrix<G2> &B, const Adapt &adapt, const typename Matrix<G1>::value_type &val, const V1 &,const V1 &) { return simd_fill_trn_function<SimdLength<V1>::RET>()(A,B,adapt,val); } 

template<class G1,class G2,class Adapt> inline void fill_trn(const Matrix<G1> &A, Matrix<G2> &B, const Adapt &adapt, const typename Matrix<G1>::value_type &val) { typename ArrayData<Matrix<G2> >::self B2=data(B); return aux_fill_trn(data(A),B2,adapt,val,typename Matrix<G1>::value_type(),typename Matrix<G2>::value_type()); }

template<class G1,class G2> inline void fill_trn(const Matrix<G1> &A, Matrix<G2> &B, const typename Matrix<G1>::value_type &val) { return fill_trn(A,B,id_adapt(),val); }




template<class V> void copy(int n, const V *pX, int dx, V *pY, int dy) 
{
  typedef typename DataVector<V>::self array_type;
  array_type X(n*dx, (V *)pX);
  array_type Y(n*dy, (V *)pY);
  
  if (dy==1)
    if (dx==1) return copy(X,Y); 
          else return copy(stride(X,dx),Y);
  else
  {
    typename StrideArray<array_type>::self Y2=stride(Y,dy);
    if (dx==1) return copy(X,Y2); 
          else return copy(stride(X,dx),Y2);
  }        
}

//{unsecret}
//Summary: Y=a*X+Y
template<class V,class G1,class G2> void axpy(const V &a, const Vector<G1> &X, Vector<G2> &Y)
{
  if (a!=V(0) && a!=V(1)) return add_copy(a,X,Y);
  if (a==V(1)) return add_copy(X,Y);
}
//{unsecret}
template<class V,class G1,class G2> void axpy(const V &a, const Matrix<G1> &X, Matrix<G2> &Y)
{
  if (a!=V(0) && a!=V(1)) return add_copy(a,X,Y);
  if (a==V(1)) return add_copy(X,Y);
}

template<class V> void axpy(int n, const V &a, const V *pX, int dx, V *pY, int dy) 
{
  typedef typename DataVector<V>::self array_type;
  array_type X(n*dx, (V *)pX);
  array_type Y(n*dy, (V *)pY);
  
  if (dy==1) 
    if (dx==1) return axpy(a,X,Y);
          else return axpy(a,stride(X,dx),Y);
  else
  {
    typename StrideArray<array_type>::self Y2=stride(Y,dy);
    if (dx==1) return axpy(a,X,Y2);
          else return axpy(a,stride(X,dx),Y2);
  }       
}



//{unsecret}
// X=a*X
template<class V,class G1> void scal(V a, Vector<G1> &X)
{
  if (a==V(1)) return;
  if (a==V(0)) return fill(X,typename Vector<G1>::value_type(0));
  return copy(a,X,X);
}
//{unsecret}
template<class V,class G1> void scal(V a, Matrix<G1> &X)
{
  if (a!=V(0) && a!=V(1)) return copy(a,X,X);
  if (a==V(0)) return fill(X,typename Matrix<G1>::value_type(0));
}


template<class V> void scal(int n, V a, V *pX, int dx) 
{
  typedef typename DataVector<V>::self array_type;
  array_type X(n*dx, pX);
  if (dx==1) scal(a,X); 
  else 
  {
    typename StrideArray<array_type>::self X2=stride(X,dx);
    scal(a,X2);
  }
}



#ifdef AUTO_BLAS
template<class G,class V,class A> inline void affect(Vector<G> &Y, const Vector<dense_vector_generator<V,A> > &X) { copy(X,Y); }
template<class G,class V,class A> inline void affect(Matrix<G> &Y, const Matrix<dense_vector_generator<V,A> > &X) { copy(X,Y); }

template<class G> inline Vector<G> &operator*=(Vector<G> &Y, const typename Vector<G>::value_type &val) { scal(val,X); return Y; }
template<class G> inline Vector<G> &operator*=(Matrix<G> &Y, const typename Matrix<G>::value_type &val) { scal(val,X); return Y; }

template<class G,class T1,class Arg,class Val,class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<T1,mul_val_function<Arg,Val,Res> > > &X) { copy(X.generator().function().value(),X.generator().array(),Y); return Y; }
template<class G,class T1,class Arg,class Val,class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<T1,mul_val_function<Arg,Val,Res> > > &X) { copy(X.generator().function().value(),X.generator().array(),Y); return Y; }
template<class G,class T1,class Arg,class Val,class Res> inline void       affect    (Vector<G> &Y, const Vector<function_array_generator<T1,val_mul_function<Arg,Val,Res> > > &X) { copy(X.generator().function().value(),X.generator().array(),Y); return Y; }
template<class G,class T1,class Arg,class Val,class Res> inline void       affect    (Matrix<G> &Y, const Matrix<function_array_generator<T1,val_mul_function<Arg,Val,Res> > > &X) { copy(X.generator().function().value(),X.generator().array(),Y); return Y; }

template<class G,class T1,class Arg,class Val,class Res> inline Vector<G> &operator+=(Vector<G> &Y, const Vector<function_array_generator<T1,mul_val_function<Arg,Val,Res> > > &X) { add_copy(X.generator().function().value(),X.generator().array(),Y); return Y; }
template<class G,class T1,class Arg,class Val,class Res> inline Vector<G> &operator+=(Matrix<G> &Y, const Matrix<function_array_generator<T1,mul_val_function<Arg,Val,Res> > > &X) { add_copy(X.generator().function().value(),X.generator().array(),Y); return Y; }
template<class G,class T1,class Arg,class Val,class Res> inline Vector<G> &operator-=(Vector<G> &Y, const Vector<function_array_generator<T1,mul_val_function<Arg,Val,Res> > > &X) { sub_copy(X.generator().function().value(),X.generator().array(),Y); return Y; }
template<class G,class T1,class Arg,class Val,class Res> inline Vector<G> &operator-=(Matrix<G> &Y, const Matrix<function_array_generator<T1,mul_val_function<Arg,Val,Res> > > &X) { sub_copy(X.generator().function().value(),X.generator().array(),Y); return Y; }

template<class G,class T1> inline void affect(Matrix<G> &Y, const Matrix<trans_matrix_generator<T1> > &X) { trn(X,Y); }
#endif




#ifdef BLAS_PRECOMPILE

#define DECLARE_BLAS_COPY(V) \
GENIAL_API void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const id_adaptor_function                                &adapt);\
GENIAL_API void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt);\
GENIAL_API void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt);\
GENIAL_API void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function <V >  > &adapt);\
GENIAL_API void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function <V >  > &adapt);

#define DEFINE_BLAS_COPY(V) \
void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const id_adaptor_function                                &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const id_adaptor_function                                &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const ax_adaptor_function   <                        V > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const add_adaptor_function  <id_adaptor_function       > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const sub_adaptor_function  <id_adaptor_function       > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Vector<data_vector_generator<V > > &X,Vector<data_vector_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function <V >  > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }\
void aux_copy(const Matrix<data_matrix_generator<V > > &X,Matrix<data_matrix_generator<V > > &Y,const add_adaptor_function  <ax_adaptor_function <V >  > &adapt) { return simd_copy_function<SimdLength<V >::RET>()(X,Y,adapt); }

DECLARE_BLAS_COPY(float)
DECLARE_BLAS_COPY(double)
DECLARE_BLAS_COPY(complex<float>)
DECLARE_BLAS_COPY(complex<double>)


#define DECLARE_BLAS_FILL(V) \
GENIAL_API void aux_fill(Vector<data_vector_generator<V > > &X,const id_adaptor_function &adapt, const V &val);\
GENIAL_API void aux_fill(Matrix<data_matrix_generator<V > > &X,const id_adaptor_function &adapt, const V &val);

#define DEFINE_BLAS_FILL(V) \
void aux_fill(Vector<data_vector_generator<V > > &X,const id_adaptor_function &adapt, const V &val) { return simd_fill_function<SimdLength<V >::RET>()(X,adapt,val); }\
void aux_fill(Matrix<data_matrix_generator<V > > &X,const id_adaptor_function &adapt, const V &val) { return simd_fill_function<SimdLength<V >::RET>()(X,adapt,val); }

DECLARE_BLAS_FILL(float)
DECLARE_BLAS_FILL(double)
DECLARE_BLAS_FILL(complex<float>)
DECLARE_BLAS_FILL(complex<double>)

#endif //BLAS_PRECOMPILE

#endif //__cplusplus



#ifdef BLAS_PRECOMPILE
#ifdef __cplusplus
extern "C" {
#endif

//Group = CBLAS Level 1

//{unsecret}
GENIAL_API void BLAS_NAME(sscal)(int n, const float  a, float  *pX, int dx);
//{unsecret}
GENIAL_API void BLAS_NAME(dscal)(int n, const double a, double *pX, int dx);
//{unsecret}
GENIAL_API void BLAS_NAME(cscal)(int n, const void  *a, void   *pX, int dx);
//{unsecret}
GENIAL_API void BLAS_NAME(zscal)(int n, const void  *a, void   *pX, int dx);

//{unsecret}
GENIAL_API void BLAS_NAME(scopy)(int n, const float  *pX, int dx, float  *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(dcopy)(int n, const double *pX, int dx, double *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(ccopy)(int n, const void   *pX, int dx, void   *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(zcopy)(int n, const void   *pX, int dx, void   *pY, int dy);

//{unsecret}
GENIAL_API void BLAS_NAME(saxpy)(int n, const float  a, const float  *pX, int dx, float  *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(daxpy)(int n, const double a, const double *pX, int dx, double *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(caxpy)(int n, const void  *a, const void   *pX, int dx, void   *pY, int dy);
//{unsecret}
GENIAL_API void BLAS_NAME(zaxpy)(int n, const void  *a, const void   *pX, int dx, void   *pY, int dy);


#ifdef __cplusplus
}
#endif
#endif //BLAS_PRECOMPILE


#endif // COPY_H
