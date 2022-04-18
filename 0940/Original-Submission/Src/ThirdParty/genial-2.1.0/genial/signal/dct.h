//GENIAL - GENeric Image Array Library
//Copyright (C) 2007  Patrick LAURENT
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

#ifndef DCT_H
#define DCT_H

#ifdef __cplusplus

//namespace genial
//{

#include "signal/fft.h"

//Group = DCT

#ifndef FFT_TWIDDLES
#define FFT_TWIDDLES 1
#endif

template<class V1,class V2,class V3> class ColorSpace;
template<class V1,class V2,class V3> class RGBColor;
template<class V1,class V2,class V3> class YCbCrColor;

template<class G1,class G2,class V>
void aux_real_dct(const Vector<G1> &X, Vector<G2> &Y, const V &)
{
  typedef Vector<G1> array_type;
  typedef PROMOTE2(float         ,typename array_type::const_value_type) value_type;
  typedef PROMOTE2(float         ,typename array_type::const_value_type) real_value_type;
  typedef PROMOTE2(complex<float>,typename array_type::const_value_type) complex_value_type;
  typedef typename value_traits<real_value_type   >::value_type real_type;
  typedef complex<real_type> complex_type;
  
  int n=X.size();

  static __thread struct tls_twiddle { int n; void *z; } twiddles[FFT_TWIDDLES]={0,NULL};
  tls_twiddle w=twiddles[0]; if (w.n==n) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.n==n) { twiddles[0]=w; goto found; } }
  {
  aligned_free(w.z); w.n=n; 
  w.z=aligned_malloc((n/2+1)*sizeof(complex_type)); 
  twiddles[0]=w; 
  typename DataVector<complex_type>::self W(n/2+1, (complex_type *)w.z); 
  W=2*real_type(inv(sqrt(real_type(2*n))))*generate_vector(n/2+1,scale(polar_function<real_type>(),-PI/(2*n))); 
  }
  found: typename DataVector<complex_type>::self W(n/2+1, (complex_type *)w.z);

  typename DenseVector<   real_value_type>::self T(n);
  typename DenseVector<complex_value_type>::self U(n);

	T[0] = X[0];
	for (int i=1, imax=(n+1)/2; i<imax; ++i) { T[i]=X[2*i]; T[n-i]=X[2*i-1]; }
	if (is_even(n)) T[n/2]=X[n-1];
  
  half_real_fft(T,U);

	Y[0] = real(W[0]) * real(U[0]) / SQRT2;
	for (int i=1, imax=(n+1)/2; i<imax; ++i) 
	{
	  complex_value_type z = W[i] * U[i];
	  Y[i  ] =  real(z);
	  Y[n-i] = -imag(z);
	}
	if (is_even(n)) Y[n/2]=real(W[n/2])*real(U[n/2]);
}

template<class G1,class G2>
inline void aux_real_dct(const Vector<G1> &X, Vector<G2> &Y)
{
  aux_real_dct(X,Y, typename Vector<G1>::value_type());
}

template<class G1,class G2>
inline void real_dct(const Vector<G1> &X, Vector<G2> &Y)
{
  typename ArrayData<Vector<G2> >::self Y2 = data(Y);
  aux_real_dct(data(X),Y2);
}

template<class G> typename DenseVector<PROMOTE2(float,typename Vector<G>::value_type)>::self real_dct(const Vector<G> &X) 
{
  typename DenseVector<PROMOTE2(float,typename Vector<G>::value_type)>::self Y(X.size());
  real_dct(X,Y);
  return Y;
}



template<class G1,class G2,class V> void aux_dct(const Vector<G1> &X, Vector<G2> &Y, const V &) { return real_dct(X,Y); }
template<class G1,class G2,class V1,class V2,class V3> void aux_dct(const Vector<G1> &X, Vector<G2> &Y, const ColorSpace<V1,V2,V3> &) 
{
  typename firstComponentArray <Vector<G2> >::self Y1=first (Y); dct(first (X),Y1);
  typename secondComponentArray<Vector<G2> >::self Y2=second(Y); dct(second(X),Y2);
  typename thirdComponentArray <Vector<G2> >::self Y3=third (Y); dct(third (X),Y3);
}
template<class G1,class G2,class V1,class V2,class V3> inline void aux_dct(const Vector<G1> &X, Vector<G2> &Y, const RGBColor  <V1,V2,V3> &x) { return aux_dct(X,Y,static_cast<const ColorSpace<V1,V2,V3> &>(x)); }
template<class G1,class G2,class V1,class V2,class V3> inline void aux_dct(const Vector<G1> &X, Vector<G2> &Y, const YCbCrColor<V1,V2,V3> &x) { return aux_dct(X,Y,static_cast<const ColorSpace<V1,V2,V3> &>(x)); }

//{unsecret}
//Summary: Discrete Cosine Transform
//Arguments: 
//  X - A signal
//Output:
//  Y - The DCT
//Example:
//  DenseVector<float>::self X(8, 1);
//  DenseVector<float>::self Y = dct(X);
//
//  DenseVector<double>::self X(8, 1);
//  dct(X,X); // in place
template<class G1,class G2> void dct(const Vector<G1> &X, Vector<G2> &Y)
{
  return aux_dct(X,Y,typename Vector<G1>::const_value_type());
}

//{unsecret}
//Return: Dense array
template<class G> PROMOTE2(float,Vector<G>) dct(const Vector<G> &X) 
{
  PROMOTE2(float,Vector<G>) Y(X.size());
  dct(X,Y);
  return Y;
}


template<class G1,class G2,class V> inline void aux_dct(const Matrix<G1> &X, Matrix<G2> &Y, const V &)
{
  Y.resize(X.size()); Y.set_lower_bound(X.lower_bound());
  for (int i=X.row_lower_bound(), imax=X.row_upper_bound()+1; i<imax; ++i) 
  {
    typename MatrixRow<Matrix<G2> >::self T = rows(Y)[i];
    dct(rows(X)[i],T);
  }
  for (int j=Y.col_lower_bound(), jmax=Y.col_upper_bound()+1; j<jmax; ++j) 
  {
    typename MatrixCol<Matrix<G2> >::self T = cols(Y)[j];
    dct(T,T);
  }
}

template<class G1,class G2>
inline void aux_dct(const Matrix<G1> &X, Matrix<G2> &Y)
{
  aux_dct(X,Y, typename Matrix<G1>::value_type());
}

//{unsecret}
template<class G1,class G2> inline void dct(const Matrix<G1> &X, Matrix<G2> &Y)
{
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_dct(data(X),Y2);
}

//{unsecret}
template<class G> PROMOTE2(float,Matrix<G>) dct(const Matrix<G> &X) 
{
  PROMOTE2(float,Matrix<G>) Y(X.size());
  dct(X,Y);
  return Y;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////


template<class G1,class G2,class V>
void aux_real_idct(const Vector<G1> &X, Vector<G2> &Y, const V &)
{
  typedef Vector<G1> array_type;
  typedef PROMOTE2(float         ,typename array_type::const_value_type) value_type;
  typedef PROMOTE2(float         ,typename array_type::const_value_type) real_value_type;
  typedef PROMOTE2(complex<float>,typename array_type::const_value_type) complex_value_type;
  typedef typename value_traits<real_value_type   >::value_type real_type;
  typedef complex<real_type> complex_type;

  int n=X.size();

  static __thread struct tls_twiddle { int n; void *z; } twiddles[FFT_TWIDDLES]={0,NULL};
  tls_twiddle w=twiddles[0]; if (w.n==n) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.n==n) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.n=n; w.z=aligned_malloc((n/2+1)*sizeof(complex_type)); twiddles[0]=w; typename DataVector<complex_type>::self W(n/2+1, (complex_type *)w.z); W=real_type(inv(sqrt(real_type(2*n))))*generate_vector(n/2+1,scale(polar_function<real_type>(),-PI/(2*n))); }
  found: typename DataVector<complex_type>::self W(n/2+1, (complex_type *)w.z);

  typename DenseVector<        value_type>::self T(n);
  typename DenseVector<complex_value_type>::self U(n);

	T[0] = SQRT2*real(W[0])*X[0];
	for (int i=1; i<(n+1)/2; ++i) 
	{
	  real_value_type a=X[i], b=X[n-i];
	  complex_value_type z(a-b,a+b);
	  z *= W[i];
	  T[i  ] = real(z);
	  T[n-i] = imag(z);
	}
	if (is_even(n)) T[n/2]=2*real(W[n/2])*X[n/2];

  half_real_fft(T,U);
 
	Y[0] = real(U[0]);
	for (int i=1, imax=(n+1)/2; i<imax; ++i)
	{
	  complex_value_type z=U[i];
	  real_value_type a=real(z), b=imag(z);
	  Y[2*i-1] = a-b;
	  Y[2*i]   = a+b;
	}
	if (is_even(n)) Y[n-1]=real(U[n/2]);
}

template<class G1,class G2>
inline void aux_real_idct(const Vector<G1> &X, Vector<G2> &Y)
{
  aux_real_idct(X,Y, typename Vector<G1>::value_type());
}

template<class G1,class G2>
inline void real_idct(const Vector<G1> &X, Vector<G2> &Y)
{
  typename ArrayData<Vector<G2> >::self Y2 = data(Y);
  aux_real_idct(data(X),Y2);
}

template<class G> typename DenseVector<PROMOTE2(float,typename Vector<G>::value_type)>::self real_idct(const Vector<G> &X) 
{
  typename DenseVector<PROMOTE2(float,typename Vector<G>::value_type)>::self Y(X.size());
  real_idct(X,Y);
  return Y;
}

template<class G1,class G2,class V> void aux_idct(const Vector<G1> &X, Vector<G2> &Y, const V &) { return real_idct(X,Y); }
template<class G1,class G2,class V1,class V2,class V3> void aux_idct(const Vector<G1> &X, Vector<G2> &Y, const ColorSpace<V1,V2,V3> &) 
{
  typename firstComponentArray <Vector<G2> >::self Y1=first (Y); idct(first (X),Y1);
  typename secondComponentArray<Vector<G2> >::self Y2=second(Y); idct(second(X),Y2);
  typename thirdComponentArray <Vector<G2> >::self Y3=third (Y); idct(third (X),Y3);
}
template<class G1,class G2,class V1,class V2,class V3> inline void aux_idct(const Vector<G1> &X, Vector<G2> &Y, const RGBColor  <V1,V2,V3> &x) { return aux_idct(X,Y,static_cast<const ColorSpace<V1,V2,V3> &>(x)); }
template<class G1,class G2,class V1,class V2,class V3> inline void aux_idct(const Vector<G1> &X, Vector<G2> &Y, const YCbCrColor<V1,V2,V3> &x) { return aux_idct(X,Y,static_cast<const ColorSpace<V1,V2,V3> &>(x)); }

//{unsecret}
//Summary: Inverse Discrete Cosine Transform
//Arguments: 
//  X - A signal
//Output:
//  Y - The IDCT
template<class G1,class G2> void idct(const Vector<G1> &X, Vector<G2> &Y)
{
  return aux_idct(X,Y,typename Vector<G1>::const_value_type());
}

//{unsecret}
//Return: Dense array
template<class G> PROMOTE2(float,Vector<G>) idct(const Vector<G> &X) 
{
  PROMOTE2(float,Vector<G>) Y(X.size());
  idct(X,Y);
  return Y;
}

template<class G1,class G2,class V> inline void aux_idct(const Matrix<G1> &X, Matrix<G2> &Y, const V &)
{
  Y.resize(X.size()); Y.set_lower_bound(X.lower_bound());
  for (int i=X.row_lower_bound(), imax=X.row_upper_bound()+1; i<imax; ++i) 
  {
    typename MatrixRow<Matrix<G2> >::self T = rows(Y)[i];
    idct(rows(X)[i],T);
  }
  for (int j=Y.col_lower_bound(), jmax=Y.col_upper_bound()+1; j<jmax; ++j) 
  {
    typename MatrixCol<Matrix<G2> >::self T = cols(Y)[j];
    idct(T,T);
  }
}

template<class G1,class G2>
inline void aux_idct(const Matrix<G1> &X, Matrix<G2> &Y)
{
  aux_idct(X,Y, typename Matrix<G1>::value_type());
}

//{unsecret}
template<class G1,class G2> inline void idct(const Matrix<G1> &X, Matrix<G2> &Y)
{
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_idct(data(X),Y2);
}

//{unsecret}
template<class G> PROMOTE2(float,Matrix<G>) idct(const Matrix<G> &X) 
{
  PROMOTE2(float,Matrix<G>) Y;
  idct(X,Y);
  return Y;
}




#ifdef FFT_PRECOMPILE

#define DECLARE_DCT(V)\
GENIAL_API void aux_real_dct (const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<V> > &Y);\
GENIAL_API void aux_real_idct(const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<V> > &Y);\
GENIAL_API void aux_dct2     (const Matrix<data_matrix_generator<V> > &X, Matrix<data_matrix_generator<V> > &Y);\
GENIAL_API void aux_idct2    (const Matrix<data_matrix_generator<V> > &X, Matrix<data_matrix_generator<V> > &Y);

#define DEFINE_DCT(V)\
void aux_real_dct (const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<V> > &Y) { aux_real_dct (X,Y,V()); }\
void aux_real_idct(const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<V> > &Y) { aux_real_idct(X,Y,V()); }\
void aux_dct2     (const Matrix<data_matrix_generator<V> > &X, Matrix<data_matrix_generator<V> > &Y) { aux_dct      (X,Y,V()); }\
void aux_idct2    (const Matrix<data_matrix_generator<V> > &X, Matrix<data_matrix_generator<V> > &Y) { aux_idct     (X,Y,V()); }

DECLARE_DCT(float)
DECLARE_DCT(double)

#endif

//} // namespace

#endif //__cplusplus




#ifdef __cplusplus
extern "C" {
#endif

#ifdef FFT_PRECOMPILE

//Group = DCT C Interface

//{unsecret}
GENIAL_API void sdct  (int n, const float  *px, float  *py);
//{unsecret}
GENIAL_API void sdct2 (int m,int n, const float  *px, float  *py);
//{unsecret}
GENIAL_API void sidct (int n, const float  *px, float  *py);
//{unsecret}
GENIAL_API void sidct2(int m,int n, const float  *px, float  *py);

//{unsecret}
GENIAL_API void ddct  (int n, const double *px, double *py);
//{unsecret}
GENIAL_API void ddct2 (int m,int n, const double *px, double *py);
//{unsecret}
GENIAL_API void didct (int n, const double *px, double *py);
//{unsecret}
GENIAL_API void didct2(int m,int n, const double *px, double *py);


#endif // FFT_PRECOMPILE

#ifdef __cplusplus
}
#endif



#endif //DCT_H

