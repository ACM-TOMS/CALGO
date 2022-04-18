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


#ifndef DOT_H
#define DOT_H

#ifdef __cplusplus

#include "blas.h"
#include "array/vector.h"


//Group = Linear Algebra

template<template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class V>
inline void dot(const Vector<G1> &X1,const Vector<G2> &Y,V &z)
{
  typedef typename Vector<G1>::value_type value_type;
    
  F<value_type,typename Vector<G2>::value_type,value_type> f;
  Func<value_type,V> func;
   
  int i=0;
  int n = Y.size();
  value_type s1(0);
  for (int imax=n-16; i<=imax; i+=16)
  {
    typename SubArray<const Vector<G2> >::self y(Y,i,16);
    typename SubArray<const Vector<G1> >::self x1(X1,i,16);  
    s1+=f(x1[ 0],y[ 0]); s1+=f(x1[ 1],y[ 1]); s1+=f(x1[ 2],y[ 2]); s1+=f(x1[ 3],y[ 3]); s1+=f(x1[ 4],y[ 4]); s1+=f(x1[ 5],y[ 5]); s1+=f(x1[ 6],y[ 6]); s1+=f(x1[ 7],y[ 7]); s1+=f(x1[ 8],y[ 8]); s1+=f(x1[ 9],y[ 9]); s1+=f(x1[10],y[10]); s1+=f(x1[11],y[11]); s1+=f(x1[12],y[12]); s1+=f(x1[13],y[13]); s1+=f(x1[14],y[14]); s1+=f(x1[15],y[15]); //s1+=f(x1[16],y[16]); s1+=f(x1[17],y[17]); s1+=f(x1[18],y[18]); s1+=f(x1[19],y[19]); s1+=f(x1[20],y[20]); s1+=f(x1[21],y[21]); s1+=f(x1[22],y[22]); s1+=f(x1[23],y[23]); s1+=f(x1[24],y[24]); s1+=f(x1[25],y[25]); s1+=f(x1[26],y[26]); s1+=f(x1[27],y[27]); s1+=f(x1[28],y[28]); s1+=f(x1[29],y[29]); s1+=f(x1[30],y[30]); s1+=f(x1[31],y[31]);
  }  
  for (; i<n; ++i) s1+=f(X1[i],Y[i]);
  z=func(s1);
}

template<template<class,class,class> class F,template<class,class> class Func,class G1,class G2,class V>
inline void dot(const Vector<G1> &X1,const Vector<G1> &X2,const Vector<G2> &Y,V &z1,V &z2)
{
  typedef typename Vector<G1>::value_type value_type;
    
  F<value_type,typename Vector<G2>::value_type,value_type> f;
  Func<value_type,V> func;

  int i=0;
  int n = Y.size();
  value_type s1(0),s2(0);
  for (int imax=n-16; i<=imax; i+=16)
  {
    typename SubArray<const Vector<G2> >::self y(Y,i,16);
    typename SubArray<const Vector<G1> >::self x1(X1,i,16),x2(X2,i,16);
    { value_type b=y[ 0]; s1+=f(x1[ 0],b); s2+=f(x2[ 0],b); }
    { value_type b=y[ 1]; s1+=f(x1[ 1],b); s2+=f(x2[ 1],b); }
    { value_type b=y[ 2]; s1+=f(x1[ 2],b); s2+=f(x2[ 2],b); }
    { value_type b=y[ 3]; s1+=f(x1[ 3],b); s2+=f(x2[ 3],b); }
    { value_type b=y[ 4]; s1+=f(x1[ 4],b); s2+=f(x2[ 4],b); }
    { value_type b=y[ 5]; s1+=f(x1[ 5],b); s2+=f(x2[ 5],b); }
    { value_type b=y[ 6]; s1+=f(x1[ 6],b); s2+=f(x2[ 6],b); }
    { value_type b=y[ 7]; s1+=f(x1[ 7],b); s2+=f(x2[ 7],b); }
    { value_type b=y[ 8]; s1+=f(x1[ 8],b); s2+=f(x2[ 8],b); }
    { value_type b=y[ 9]; s1+=f(x1[ 9],b); s2+=f(x2[ 9],b); }
    { value_type b=y[10]; s1+=f(x1[10],b); s2+=f(x2[10],b); }
    { value_type b=y[11]; s1+=f(x1[11],b); s2+=f(x2[11],b); }
    { value_type b=y[12]; s1+=f(x1[12],b); s2+=f(x2[12],b); }
    { value_type b=y[13]; s1+=f(x1[13],b); s2+=f(x2[13],b); }
    { value_type b=y[14]; s1+=f(x1[14],b); s2+=f(x2[14],b); }
    { value_type b=y[15]; s1+=f(x1[15],b); s2+=f(x2[15],b); }
  }     
  for (; i<n; ++i) { value_type y=Y[i]; s1+=f(X1[i],y); s2+=f(X2[i],y); }
  z1=func(s1); z2=func(s2);
}

template<template<class,class> class Func,class G1,class G2,class V>
inline void dot(const Vector<G1> &X,const Vector<G2> &Y, V &z)
{
  return dot<dot_function,Func>(X,Y,z);
}


template<int N,template<class,class,class> class F> 
struct simd_dot_function    
{
  template<class G1,class G2,class V> 
  inline void operator()(const Vector<G1> &X,const Vector<G2> &Y,V &v) const
  {
    int n=X.size(); 
    if (n%N==0) return dot<F,sum_function>(simd_block(X),simd_block(Y),v); 
    return dot<F,id_function>(X,Y,v); 
  } 
};
template<template<class,class,class> class F> struct simd_dot_function<0,F> { template<class G1,class G2,class V> inline void operator()(const Vector<G1> &X,const Vector<G2> &Y,V &v) const { return dot<F,id_function>(X,Y,v); } };

template<                  class V> inline void aux_dot (const Vector<data_vector_generator<V> > &X,const Vector<data_vector_generator<V> > &Y,V &v) { return simd_dot_function<SimdLength<V>::RET, dot_function>()(X,Y,v); }
template<                  class V> inline void aux_cdot(const Vector<data_vector_generator<V> > &X,const Vector<data_vector_generator<V> > &Y,V &v) { return simd_dot_function<SimdLength<V>::RET,cdot_function>()(X,Y,v); }
template<class G1,class G2,class V> inline void aux_dot (const Vector<G1                       > &X,const Vector<G2                       > &Y,V &v) { return simd_dot_function<0                 , dot_function>()(X,Y,v); } 
template<class G1,class G2,class V> inline void aux_cdot(const Vector<G1                       > &X,const Vector<G2                       > &Y,V &v) { return simd_dot_function<0                 ,cdot_function>()(X,Y,v); } 

template<class G1,class G2,class V> inline void dot (const Vector<G1> &X, const Vector<G2> &Y, V &v) { return aux_dot (data(X),data(Y),v); }
template<class G1,class G2,class V> inline void cdot(const Vector<G1> &X, const Vector<G2> &Y, V &v) { return aux_cdot(data(X),data(Y),v); }

//{unsecret}
//Summary: X.Y
template<class G1,class G2> inline PROMOTE2(typename Vector<G1>::value_type,typename Vector<G2>::value_type) dot (const Vector<G1> &X, const Vector<G2> &Y) { PROMOTE2(typename Vector<G1>::value_type,typename Vector<G2>::value_type) val; dot (X,Y,val); return val; }
//{unsecret}
//Summary: conj(X).Y
template<class G1,class G2> inline PROMOTE2(typename Vector<G1>::value_type,typename Vector<G2>::value_type) cdot(const Vector<G1> &X, const Vector<G2> &Y) { PROMOTE2(typename Vector<G1>::value_type,typename Vector<G2>::value_type) val; cdot(X,Y,val); return val; }


template<class V> V dot(int n, const V *pX, int dx, const V *pY, int dy) 
{  
  typename DataVector<V>::self X(n*dx, (V *)pX);
  typename DataVector<V>::self Y(n*dy, (V *)pY);
  if (dx==1) 
    if (dy==1) return dot(X,Y);
          else return dot(X,stride(Y,dy));
  else
    if (dy==1) return dot(Y,stride(X,dx));
          else return dot(stride(X,dx),stride(Y,dy));
}

template<class V> V cdot(int n, const V *pX, int dx, const V *pY, int dy) 
{  
  typename DataVector<V>::self X(n*dx, (V *)pX);
  typename DataVector<V>::self Y(n*dy, (V *)pY);
  if (dx==1) 
    if (dy==1) return cdot(X,Y);
          else return cdot(X,stride(Y,dy));
  else
    if (dy==1) return cdot(Y,stride(X,dx));
          else return cdot(stride(X,dx),stride(Y,dy));
}


#ifdef BLAS_PRECOMPILE

#define DECLARE_BLAS_DOT(V)\
GENIAL_API void aux_dot (const Vector<data_vector_generator<V > > &X, const Vector<data_vector_generator<V > > &Y, V &v);\
GENIAL_API void aux_cdot(const Vector<data_vector_generator<V > > &X, const Vector<data_vector_generator<V > > &Y, V &v);

#define DEFINE_BLAS_DOT(V)\
void aux_dot (const Vector<data_vector_generator<V > > &X, const Vector<data_vector_generator<V > > &Y, V &v) { return simd_dot_function<SimdLength<V >::RET,dot_function >()(X,Y,v); }\
void aux_cdot(const Vector<data_vector_generator<V > > &X, const Vector<data_vector_generator<V > > &Y, V &v) { return simd_dot_function<SimdLength<V >::RET,dot_function >()(X,Y,v); }

DECLARE_BLAS_DOT(float)
DECLARE_BLAS_DOT(double)
DECLARE_BLAS_DOT(complex<float>)
DECLARE_BLAS_DOT(complex<double>)

#endif //BLAS_PRECOMPILE

#endif //__cplusplus



#ifdef BLAS_PRECOMPILE
#ifdef __cplusplus
extern "C" {
#endif

//Group = CBLAS Level 1

//{unsecret}
GENIAL_API float  BLAS_NAME(sdot)     (int n, const float  *pX, int dx, const float  *pY, int dy            );
//{unsecret}
GENIAL_API double BLAS_NAME(ddot)     (int n, const double *pX, int dx, const double *pY, int dy            );
//{unsecret}
GENIAL_API void   BLAS_NAME(cdotu_sub)(int n, const void   *pX, int dx, const void   *pY, int dy, void *dotu);
//{unsecret}
GENIAL_API void   BLAS_NAME(zdotu_sub)(int n, const void   *pX, int dx, const void   *pY, int dy, void *dotu);
//{unsecret}
GENIAL_API void   BLAS_NAME(cdotc_sub)(int n, const void   *pX, int dx, const void   *pY, int dy, void *dotc);
//{unsecret}
GENIAL_API void   BLAS_NAME(zdotc_sub)(int n, const void   *pX, int dx, const void   *pY, int dy, void *dotc);

#ifdef __cplusplus
}
#endif
#endif //BLAS_PRECOMPILE


#endif // DOT_H







