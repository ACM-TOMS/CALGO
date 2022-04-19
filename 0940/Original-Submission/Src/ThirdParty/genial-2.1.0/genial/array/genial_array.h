//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
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

#ifndef ARRAY_H
#define ARRAY_H

#include "array/arraygenerator.h"
#include "array/arrayfunction.h"
#include "numeric.h"

//namespace genial
//{

using namespace std;

template<class G> class Vector;
template<class G> class Matrix;


template<int D, class G> class Array { };

template<class V,class A> struct DenseVector;
template<class V,class A> struct DenseMatrix;
template<int N,class V> struct TinyVector;
template<int M,int N,class V> struct TinyMatrix;

template<class V,class A> class dense_vector_generator;
template<class V,class A> class dense_matrix_generator;
template<int N,class V> class tiny_vector_generator;
template<int N,class V> class simd_vector_generator;
template<int M,int N,class V> class tiny_matrix_generator;
template<int N,class T,int C> class tiny_row_dense_matrix_generator;
template<      class T,int C> class row_dense_matrix_generator;

template<int D,class V> struct DenseArray { };
template<class V> struct DenseArray<1,V> { typedef Array<1,dense_vector_generator<V> > self; };
template<class V> struct DenseArray<2,V> { typedef Array<2,dense_matrix_generator<V> > self; };

template<int D,class V> struct shiftDenseArray { typedef typename shiftArray<typename DenseArray<D,V>::self,1>::self self; };


template<class G>       struct promotion_traits<Vector<G>  > { typedef typename Vector<G> ::template rebind<typename Vector<G> ::const_value_type>::other value_type; };
template<class G>       struct promotion_traits<Matrix<G>  > { typedef typename Matrix<G> ::template rebind<typename Matrix<G> ::const_value_type>::other value_type; };
template<int D,class G> struct promotion_traits<Array<D,G> > { typedef typename Array<D,G>::template rebind<typename Array<D,G>::const_value_type>::other value_type; };




template<class G ,class V > struct promotion2_traits<Vector<G> , V          > { typedef typename Vector<G>::template rebind<PROMOTE2(typename Vector<G>::const_value_type,V)>::other value_type; };
template<class V ,class G > struct promotion2_traits<V         , Vector<G>  > { typedef typename Vector<G>::template rebind<PROMOTE2(V,typename Vector<G>::const_value_type)>::other value_type; };
template<class G1,class G2> struct promotion2_traits<Vector<G1>, Vector<G2> > 
{
  typedef typename Vector<G1>::template rebind<typename Vector<G1>::const_value_type>::other other1;
  typedef typename Vector<G2>::template rebind<typename Vector<G2>::const_value_type>::other other2;
  typedef PROMOTE2(other1,other2) value_type;
};

template<       class V1,       class V2> struct promotion2_traits<Vector<dense_vector_generator<   V1> >,Vector<dense_vector_generator<   V2> > > { typedef typename DenseVector<   typename promotion2_traits<typename Vector<dense_vector_generator<   V1> >::const_value_type,typename Vector<dense_vector_generator<   V2> >::const_value_type>::value_type>::self value_type; };
template<       class V1,int N2,class V2> struct promotion2_traits<Vector<dense_vector_generator<   V1> >,Vector<tiny_vector_generator <N2,V2> > > { typedef typename TinyVector <N2,typename promotion2_traits<typename Vector<dense_vector_generator<   V1> >::const_value_type,typename Vector<tiny_vector_generator <N2,V2> >::const_value_type>::value_type>::self value_type; };
template<int N1,class V1,       class V2> struct promotion2_traits<Vector<tiny_vector_generator <N1,V1> >,Vector<dense_vector_generator<   V2> > > { typedef typename TinyVector <N1,typename promotion2_traits<typename Vector<tiny_vector_generator <N1,V1> >::const_value_type,typename Vector<dense_vector_generator<   V2> >::const_value_type>::value_type>::self value_type; };
template<int N1,class V1,int N2,class V2> struct promotion2_traits<Vector<tiny_vector_generator <N1,V1> >,Vector<tiny_vector_generator <N2,V2> > > { typedef typename TinyVector <MAX(N1,N2) ,typename promotion2_traits<typename Vector<tiny_vector_generator <N1,V1> >::const_value_type,typename Vector<tiny_vector_generator <N2,V2> >::const_value_type>::value_type>::self value_type; };

template<int N,class V> struct promotion2_traits<Vector<simd_vector_generator <N,V> >,Vector<simd_vector_generator <N,V> > > { typedef Vector<simd_vector_generator <N,V> > value_type; };


template<class G, class V > struct promotion2_traits<Matrix<G>,  V          > { typedef typename Matrix<G>::template rebind<PROMOTE2(typename Matrix<G>::const_value_type,V)>::other value_type; };
template<class V, class G > struct promotion2_traits<V,          Matrix<G>  > { typedef typename Matrix<G>::template rebind<PROMOTE2(V,typename Matrix<G>::const_value_type)>::other value_type; };
template<class G1,class G2> struct promotion2_traits<Matrix<G1>, Matrix<G2> > { typedef typename shiftDenseMatrix<PROMOTE2(typename Matrix<G1>::const_value_type, typename Matrix<G2>::const_value_type)>::self value_type; };

template<int D,class G,class V>   struct promotion2_traits<Array<D,G >,V           > { typedef typename Array<D,G>::template rebind<typename promotion2_traits<typename Array<D,G>::const_value_type,V>::value_type>::other value_type; };
template<class V,int D,class G>   struct promotion2_traits<V          ,Array<D,G > > { typedef typename Array<D,G>::template rebind<typename promotion2_traits<V,typename Array<D,G>::const_value_type>::value_type>::other value_type; };
template<int D,class G1,class G2> struct promotion2_traits<Array<D,G1>,Array<D,G2> > { typedef typename DenseArray<D,typename promotion2_traits<typename Array<D,G1>::const_value_type,typename Array<D,G2>::const_value_type>::value_type>::self value_type; };

template<class T1,class T2> struct promotion2_traits;
template<class T,class V> struct promotion2_traits<complex<T>,V > { typedef complex<typename promotion2_traits<T,V>::value_type> value_type; };
template<class V,class T> struct promotion2_traits<V,complex<T> > { typedef complex<typename promotion2_traits<V,T>::value_type> value_type; };
template<class T,class V> struct promotion2_traits<complex<T>,const V > { typedef complex<typename promotion2_traits<T,const V>::value_type> value_type; };
template<class V,class T> struct promotion2_traits<const V,complex<T> > { typedef complex<typename promotion2_traits<const V,T>::value_type> value_type; };
template<class T1,class T2> struct promotion2_traits<complex<T1>,complex<T2> > { typedef complex<typename promotion2_traits<T1,T2>::value_type> value_type; };


template<class G> struct value_type_traits<      Vector<G> > { typedef typename Vector<G>::value_type value_type; };
template<class G> struct value_type_traits<const Vector<G> > { typedef typename Vector<G>::value_type value_type; };
template<class G> struct value_type_traits<      Matrix<G> > { typedef typename Matrix<G>::value_type value_type; };
template<class G> struct value_type_traits<const Matrix<G> > { typedef typename Matrix<G>::value_type value_type; };


//group=Arrays functions

template<class V> struct ArrayDense          { typedef       V &self; };
template<class V> struct ArrayDense<const V> { typedef const V &self; };
template<class G> struct ArrayDense<      Vector<G> > { typedef       typename DenseVector<typename Vector<G>::const_value_type>::self self; };
template<class G> struct ArrayDense<const Vector<G> > { typedef const typename DenseVector<typename Vector<G>::const_value_type>::self self; };
template<class G> struct ArrayDense<      Matrix<G> > { typedef       typename DenseMatrix<typename Matrix<G>::const_value_type>::self self; };
template<class G> struct ArrayDense<const Matrix<G> > { typedef const typename DenseMatrix<typename Matrix<G>::const_value_type>::self self; };
template<class V            > struct ArrayDense<      Vector<data_vector_generator          <V    > > > { typedef       Vector<data_vector_generator          <V    > > &self; };
template<class V            > struct ArrayDense<const Vector<data_vector_generator          <V    > > > { typedef const Vector<data_vector_generator          <V    > > &self; };
template<class V            > struct ArrayDense<      Matrix<data_matrix_generator          <V    > > > { typedef       Matrix<data_matrix_generator          <V    > > &self; };
template<class V            > struct ArrayDense<const Matrix<data_matrix_generator          <V    > > > { typedef const Matrix<data_matrix_generator          <V    > > &self; };
template<class V,class A    > struct ArrayDense<      Vector<dense_vector_generator         <V,A  > > > { typedef       Vector<dense_vector_generator         <V,A  > > &self; };
template<class V,class A    > struct ArrayDense<const Vector<dense_vector_generator         <V,A  > > > { typedef const Vector<dense_vector_generator         <V,A  > > &self; };
template<class V,class A    > struct ArrayDense<      Matrix<dense_matrix_generator         <V,A  > > > { typedef       Matrix<dense_matrix_generator         <V,A  > > &self; };
template<class V,class A    > struct ArrayDense<const Matrix<dense_matrix_generator         <V,A  > > > { typedef const Matrix<dense_matrix_generator         <V,A  > > &self; };
template<int N,      class V> struct ArrayDense<      Vector<tiny_vector_generator          <N  ,V> > > { typedef       Vector<tiny_vector_generator          <N  ,V> > &self; };
template<int N,      class V> struct ArrayDense<const Vector<tiny_vector_generator          <N  ,V> > > { typedef const Vector<tiny_vector_generator          <N  ,V> > &self; };
template<int M,int N,class V> struct ArrayDense<      Matrix<tiny_matrix_generator          <M,N,V> > > { typedef       Matrix<tiny_matrix_generator          <M,N,V> > &self; };
template<int M,int N,class V> struct ArrayDense<const Matrix<tiny_matrix_generator          <M,N,V> > > { typedef const Matrix<tiny_matrix_generator          <M,N,V> > &self; };
template<      class T,int C> struct ArrayDense<      Vector<sub_dense_vector_generator     <T  ,C> > > { typedef       Vector<sub_dense_vector_generator     <T  ,C> > &self; }; 
template<      class T,int C> struct ArrayDense<const Vector<sub_dense_vector_generator     <T  ,C> > > { typedef const Vector<sub_dense_vector_generator     <T  ,C> > &self; }; 
template<      class T,int C> struct ArrayDense<      Vector<row_dense_matrix_generator     <T  ,C> > > { typedef       Vector<row_dense_matrix_generator     <T  ,C> > &self; }; 
template<      class T,int C> struct ArrayDense<const Vector<row_dense_matrix_generator     <T  ,C> > > { typedef const Vector<row_dense_matrix_generator     <T  ,C> > &self; }; 
template<int N,class T,int C> struct ArrayDense<      Vector<tiny_sub_dense_vector_generator<N,T,C> > > { typedef       Vector<tiny_sub_dense_vector_generator<N,T,C> > &self; }; 
template<int N,class T,int C> struct ArrayDense<const Vector<tiny_sub_dense_vector_generator<N,T,C> > > { typedef const Vector<tiny_sub_dense_vector_generator<N,T,C> > &self; }; 
template<int N,class T,int C> struct ArrayDense<      Vector<tiny_row_dense_matrix_generator<N,T,C> > > { typedef       Vector<tiny_row_dense_matrix_generator<N,T,C> > &self; }; 
template<int N,class T,int C> struct ArrayDense<const Vector<tiny_row_dense_matrix_generator<N,T,C> > > { typedef const Vector<tiny_row_dense_matrix_generator<N,T,C> > &self; }; 


//{unsecret}
//{noAutoLink}
//Summary: Converts (if necessary) an array in a temporary dense array
template<class G> inline typename ArrayDense<      Vector<G> >::self dense(      Vector<G> &X) { return typename ArrayDense<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename ArrayDense<const Vector<G> >::self dense(const Vector<G> &X) { return typename ArrayDense<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename ArrayDense<      Matrix<G> >::self dense(      Matrix<G> &X) { return typename ArrayDense<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename ArrayDense<const Matrix<G> >::self dense(const Matrix<G> &X) { return typename ArrayDense<const Matrix<G> >::self(X); }


template<class T> struct dense_vector_function : public unary_value_function<T,typename T::template dense_rebind<typename T::const_value_type>::other> { typedef unary_value_function<T,typename T::template dense_rebind<typename T::const_value_type>::other> base; UNARY_FUNCTION_BASE_TYPES; const_reference operator()(argument_type &x) const { return dense(x); } };
template<class T> struct dense_matrix_function : public unary_value_function<T,typename T::template dense_rebind<typename T::const_value_type>::other> { typedef unary_value_function<T,typename T::template dense_rebind<typename T::const_value_type>::other> base; UNARY_FUNCTION_BASE_TYPES; const_reference operator()(argument_type &x) const { return dense(x); } };

template<class G> typename UnaryFunctionArray<const Vector<G>,const dense_vector_function<const typename Vector<G>::const_reference> >::self ele_dense(const Vector<G> &X) { return apply(X,dense_vector_function<const typename Vector<G>::const_reference>()); }
template<class G> typename UnaryFunctionArray<const Matrix<G>,const dense_matrix_function<const typename Matrix<G>::const_reference> >::self ele_dense(const Matrix<G> &X) { return apply(X,dense_matrix_function<const typename Matrix<G>::const_reference>()); }




template<class V> struct ArrayData          { typedef       V &self; };
template<class V> struct ArrayData<const V> { typedef const V &self; };
template<class G> struct ArrayData<      Vector<G> > { typedef       Vector<G> &self; };
template<class G> struct ArrayData<const Vector<G> > { typedef const Vector<G> &self; };
template<class G> struct ArrayData<      Matrix<G> > { typedef       Matrix<G> &self; };
template<class G> struct ArrayData<const Matrix<G> > { typedef const Matrix<G> &self; };
template<class V            > struct ArrayData<      Vector<data_vector_generator          <V    > > > { typedef       Vector<data_vector_generator<V> > &self; };
template<class V            > struct ArrayData<const Vector<data_vector_generator          <V    > > > { typedef const Vector<data_vector_generator<V> > &self; };
template<class V            > struct ArrayData<      Matrix<data_matrix_generator          <V    > > > { typedef       Matrix<data_matrix_generator<V> > &self; };
template<class V            > struct ArrayData<const Matrix<data_matrix_generator          <V    > > > { typedef const Matrix<data_matrix_generator<V> > &self; };
template<class V,class A    > struct ArrayData<      Vector<dense_vector_generator         <V,A  > > > { typedef       typename DataVector<V>::self self; };
template<class V,class A    > struct ArrayData<const Vector<dense_vector_generator         <V,A  > > > { typedef const typename DataVector<V>::self self; };
template<class V,class A    > struct ArrayData<      Matrix<dense_matrix_generator         <V,A  > > > { typedef       typename DataMatrix<V>::self self; };
template<class V,class A    > struct ArrayData<const Matrix<dense_matrix_generator         <V,A  > > > { typedef const typename DataMatrix<V>::self self; };
template<int N,      class V> struct ArrayData<      Vector<tiny_vector_generator          <N  ,V> > > { typedef       typename DataVector<V>::self self; };
template<int N,      class V> struct ArrayData<const Vector<tiny_vector_generator          <N  ,V> > > { typedef const typename DataVector<V>::self self; };
template<int M,int N,class V> struct ArrayData<      Matrix<tiny_matrix_generator          <M,N,V> > > { typedef       typename DataMatrix<V>::self self; };
template<int M,int N,class V> struct ArrayData<const Matrix<tiny_matrix_generator          <M,N,V> > > { typedef const typename DataMatrix<V>::self self; };
template<      class T,int C> struct ArrayData<      Vector<sub_dense_vector_generator     <T  ,C> > > { typedef       typename DataVector<typename T::const_value_type>::self self; }; 
template<      class T,int C> struct ArrayData<const Vector<sub_dense_vector_generator     <T  ,C> > > { typedef const typename DataVector<typename T::const_value_type>::self self; }; 
template<      class T,int C> struct ArrayData<      Vector<row_dense_matrix_generator     <T  ,C> > > { typedef       typename DataVector<typename T::const_value_type>::self self; }; 
template<      class T,int C> struct ArrayData<const Vector<row_dense_matrix_generator     <T  ,C> > > { typedef const typename DataVector<typename T::const_value_type>::self self; }; 
template<int N,class T,int C> struct ArrayData<      Vector<tiny_sub_dense_vector_generator<N,T,C> > > { typedef       typename DataVector<typename T::const_value_type>::self self; }; 
template<int N,class T,int C> struct ArrayData<const Vector<tiny_sub_dense_vector_generator<N,T,C> > > { typedef const typename DataVector<typename T::const_value_type>::self self; }; 
template<int N,class T,int C> struct ArrayData<      Vector<tiny_row_dense_matrix_generator<N,T,C> > > { typedef       typename DataVector<typename T::const_value_type>::self self; }; 
template<int N,class T,int C> struct ArrayData<const Vector<tiny_row_dense_matrix_generator<N,T,C> > > { typedef const typename DataVector<typename T::const_value_type>::self self; }; 


//{unsecret}
//{noAutoLink}
//Summary: Converts (if possible) an array in a temporary data-array
template<class G> inline typename ArrayData<      Vector<G> >::self data(      Vector<G> &X) { return typename ArrayData<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename ArrayData<const Vector<G> >::self data(const Vector<G> &X) { return typename ArrayData<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename ArrayData<      Matrix<G> >::self data(      Matrix<G> &X) { return typename ArrayData<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename ArrayData<const Matrix<G> >::self data(const Matrix<G> &X) { return typename ArrayData<const Matrix<G> >::self(X); }



//{unsecret}
//Summary: Writes all the elements in a string
template<class G> string to_string(const Vector<G> &X) { return to_string(X.begin(), X.end()); }
//{unsecret}
template<class G> string to_string(const Matrix<G> &X) { return to_string(X.begin(), X.end()); }
//{unsecret}
template<int D,class G> string to_string(const Array<D,G> &X) { return to_string(X.begin(), X.end()); }

//{unsecret}
template<class G> void operator>>(const string &s, Vector<G> &X) { istringstream is(s.c_str()); for (typename Vector<G>::iterator i=X.begin(); i!=X.end(); ++i) is>>*i; return; }
//{unsecret}
template<class G> void operator>>(const string &s, Matrix<G> &X) { istringstream is(s.c_str()); for (typename Matrix<G>::iterator i=X.begin(); i!=X.end(); ++i) is>>*i; return; }

//{unsecret}
//Summary: Swaps the 2 arrays
template<class G> inline void swap(Vector<G> &X, Vector<G> &Y) { X.swap(Y); }
//{unsecret}
template<class G> inline void swap(Matrix<G> &X, Matrix<G> &Y) { X.swap(Y); }

//{unsecret}
//Summary: Swaps the contents of 2 arrays
template<class G> inline void swap_ranges(Vector<G> &X, Vector<G> &Y) { swap_ranges_n(X.nelms(),X.begin(),Y.begin()); }
//{unsecret}
template<class G> inline void swap_ranges(Matrix<G> &X, Matrix<G> &Y) { swap_ranges_n(X.nelms(),X.begin(),Y.begin()); }


//{unsecret}
template<class V ,class G > inline Vector<G > &operator+=(Vector<G > &X, const V          &v) { X=X+v; return X; }
//{unsecret}
template<class V ,class G > inline Matrix<G > &operator+=(Matrix<G > &X, const V          &v) { X=X+v; return X; }
//{unsecret}
template<class G1,class G2> inline Vector<G1> &operator+=(Vector<G1> &X, const Vector<G2> &Y) { X=X+Y; return X; }
//{unsecret}
template<class G1,class G2> inline Matrix<G1> &operator+=(Matrix<G1> &X, const Matrix<G2> &Y) { X=X+Y; return X; }

//{unsecret}
template<class V ,class G > inline Vector<G > &operator-=(Vector<G > &X, const V          &v) { X=X-v; return X; }
//{unsecret}
template<class V ,class G > inline Matrix<G > &operator-=(Matrix<G > &X, const V          &v) { X=X-v; return X; }
//{unsecret}
template<class G1,class G2> inline Vector<G1> &operator-=(Vector<G1> &X, const Vector<G2> &Y) { X=X-Y; return X; }
//{unsecret}
template<class G1,class G2> inline Matrix<G1> &operator-=(Matrix<G1> &X, const Matrix<G2> &Y) { X=X-Y; return X; }

//{unsecret}
template<class V, class G> inline Vector<G> &operator*=(Vector<G> &X, const V &v) { X=X*v; return X; }
//{unsecret}
template<class V, class G> inline Matrix<G> &operator*=(Matrix<G> &X, const V &v) { X=X*v; return X; }

//{unsecret}
template<class V, class G> inline Vector<G> &operator/=(Vector<G> &X, const V &v) { X=X/v; return X; }
//{unsecret}
template<class V, class G> inline Matrix<G> &operator/=(Matrix<G> &X, const V &v) { X=X/v; return X; }

//{unsecret}
//Summary: Lowest position of the lowest value
//Return: An iterator describing the position
template<class G> inline typename Vector<G>::iterator       min_element(      Vector<G> &X) { return min_element(X.begin(), X.end()); }
//{unsecret}
template<class G> inline typename Matrix<G>::iterator       min_element(      Matrix<G> &X) { return min_element(X.begin(), X.end()); }
//{unsecret}
template<class G> inline typename Vector<G>::const_iterator min_element(const Vector<G> &X) { return min_element(X.begin(), X.end()); }
//{unsecret}
template<class G> inline typename Matrix<G>::const_iterator min_element(const Matrix<G> &X) { return min_element(X.begin(), X.end()); }

//{unsecret}
//Summary: Lowest position of the largest value
//Return: An iterator describing the position
template<class G> inline typename Vector<G>::iterator       max_element(      Vector<G> &X) { return max_element(X.begin(), X.end()); }
//{unsecret}
template<class G> inline typename Matrix<G>::iterator       max_element(      Matrix<G> &X) { return max_element(X.begin(), X.end()); }
//{unsecret}
template<class G> inline typename Vector<G>::const_iterator max_element(const Vector<G> &X) { return max_element(X.begin(), X.end()); }
//{unsecret}
template<class G> inline typename Matrix<G>::const_iterator max_element(const Matrix<G> &X) { return max_element(X.begin(), X.end()); }


//{unsecret}
//Summary: Lowest value
template<class G> inline typename Vector<G>::const_value_type min(const Vector<G> &X) { return min_n(X.nelms(), X.begin()); }
//{unsecret}
template<class G> inline typename Matrix<G>::const_value_type min(const Matrix<G> &X) { return min_n(X.nelms(), X.begin()); }


//{unsecret}
//Summary: Largest value
template<class G> inline typename Vector<G>::const_value_type max(const Vector<G> &X) { return max_n(X.nelms(), X.begin()); }
//{unsecret}
template<class G> inline typename Matrix<G>::const_value_type max(const Matrix<G> &X) { return max_n(X.nelms(), X.begin()); }


template<      int N,class G> inline typename Vector<G>::const_value_type tiny_sum(const Vector<G> &X) { assert(X.nelms()!=0); return accumulate_n<N-1>(++X.begin(), typename Vector<G>::const_value_type(*X.begin())); }
template<int M,int N,class G> inline typename Matrix<G>::const_value_type tiny_sum(const Matrix<G> &X) { return sum(ele_sum(rows(X))); }


template<int N,class V      > inline typename Vector<tiny_vector_generator          <N,V  > >::const_value_type sum(const Vector<tiny_vector_generator          <N,V  > > &X) { return tiny_sum<N>(X); }
template<int N,class T,int C> inline typename Vector<tiny_row_dense_matrix_generator<N,T,C> >::const_value_type sum(const Vector<tiny_row_dense_matrix_generator<N,T,C> > &X) { return tiny_sum<N>(X); }
template<int N,class T,int C> inline typename Vector<tiny_sub_vector_generator      <N,T,C> >::const_value_type sum(const Vector<tiny_sub_vector_generator      <N,T,C> > &X) { return tiny_sum<N>(X); }
template<int N,class T,int C> inline typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::const_value_type sum(const Vector<tiny_sub_dense_vector_generator<N,T,C> > &X) { return tiny_sum<N>(X); }

template<int M,int N,class V      > inline typename Matrix<tiny_matrix_generator          <M,N,V  > >::const_value_type sum(const Matrix<tiny_matrix_generator          <M,N,V  > > &X) { return tiny_sum<M,N>(X); }
template<int M,int N,class T,int C> inline typename Matrix<tiny_sub_matrix_generator      <M,N,T,C> >::const_value_type sum(const Matrix<tiny_sub_matrix_generator      <M,N,T,C> > &X) { return tiny_sum<M,N>(X); }

//{unsecret}
//{noAutoLink}
//summary: Sum of all values
template<class G> inline typename Vector<G>::const_value_type sum(const Vector<G> &X) { return accumulate_n(X.nelms(), X.begin(), typename Vector<G>::const_value_type(0)); }
//{unsecret}
template<class G> inline typename Matrix<G>::const_value_type sum(const Matrix<G> &X) { return accumulate_n(X.nelms(), X.begin(), typename Matrix<G>::const_value_type(0)); }


template<class G> inline pair<typename Vector<G>::const_value_type,typename Vector<G>::const_value_type> square_sum(const Vector<G> &X) { assert(X.nelms()!=0); return square_accumulate(++X.begin(), X.end(), typename Vector<G>::const_value_type(*X.begin()), sqr(typename Vector<G>::const_value_type(*X.begin()))); }
template<class G> inline pair<typename Matrix<G>::const_value_type,typename Matrix<G>::const_value_type> square_sum(const Matrix<G> &X) { assert(X.nelms()!=0); return square_accumulate(++X.begin(), X.end(), typename Matrix<G>::const_value_type(*X.begin()), sqr(typename Matrix<G>::const_value_type(*X.begin()))); }


//{unsecret}
//{noAutoLink}
//Summary: Average
template<class G> inline PROMOTE2(float,typename Vector<G>::const_value_type) mean(const Vector<G> &X) { assert(X.nelms()!=0); return PROMOTE2(float,typename Vector<G>::const_value_type)(sum(X))/X.nelms(); }
//{unsecret}
template<class G> inline PROMOTE2(float,typename Matrix<G>::const_value_type) mean(const Matrix<G> &X) { assert(X.nelms()!=0); return PROMOTE2(float,typename Matrix<G>::const_value_type)(sum(X))/X.nelms(); }

//{unsecret}
//{noAutoLink}
//Summary: Centroid/Barycenter of the values
template<class G> PROMOTE2(float,typename Matrix<G>::index_type)
centroid(const Matrix<G> &X)
{
  typedef PROMOTE2(float,typename Matrix<G>::index_type) centroid_type;
  typedef typename centroid_type::value_type value_type;
  centroid_type g(0);
  value_type s=0;
  for (int i=0, imax=X.nrows(); i<imax; ++i)
    for (int j=0, jmax=X.ncols(); j<jmax; ++j)
    {
      value_type v(X(i,j));
      s+=v;
      g+=v*centroid_type(i,j);
    }
  return g/s;
}

//{unsecret}
//Summary: Sum of the absolute differences
template<class G1,class G2> inline typename abs_function<PROMOTE2(typename Vector<G1>::value_type, typename Vector<G2>::value_type)>::result_type dist(const Vector<G1> &A, const Vector<G2> &B) { return sum(abs(A-B)); }
//{unsecret}
template<class G1,class G2> inline typename abs_function<PROMOTE2(typename Matrix<G1>::value_type, typename Matrix<G2>::value_type)>::result_type dist(const Matrix<G1> &A, const Matrix<G2> &B) { return sum(abs(A-B)); }

//{unsecret}
//Summary: Standard deviation
//Arguments:
//  X - The array
//  m - Sum of the squared values of X, if already known
template<class G> inline typename promotion2_traits<float,typename Vector<G>::const_value_type>::value_type std_dev(const Vector<G> &X, const typename promotion2_traits<float,typename Vector<G>::const_value_type>::value_type &m) { return sqrt(promotion2_traits<float,typename Vector<G>::const_value_type>::value_type(sum(ele_sqr(X)))/X.nelms()-sqr(m)); }
//{unsecret}
template<class G> inline typename promotion2_traits<float,typename Matrix<G>::const_value_type>::value_type std_dev(const Matrix<G> &X, const typename promotion2_traits<float,typename Matrix<G>::const_value_type>::value_type &m) { return sqrt(promotion2_traits<float,typename Matrix<G>::const_value_type>::value_type(sum(ele_sqr(X)))/X.nelms()-sqr(m)); }
//{unsecret}
template<class G> inline typename promotion2_traits<float,typename Vector<G>::const_value_type>::value_type std_dev(const Vector<G> &X) { typedef typename Vector<G>::const_value_type const_value_type; typedef typename promotion2_traits<float,const_value_type>::value_type value_type; pair<const_value_type,const_value_type> sum = square_sum(X); return sqrt(value_type(sum.second)/X.nelms()-sqr(value_type(sum.first)/X.nelms())); }
//{unsecret}
template<class G> inline typename promotion2_traits<float,typename Matrix<G>::const_value_type>::value_type std_dev(const Matrix<G> &X) { typedef typename Matrix<G>::const_value_type const_value_type; typedef typename promotion2_traits<float,const_value_type>::value_type value_type; pair<const_value_type,const_value_type> sum = square_sum(X); return sqrt(value_type(sum.second)/X.nelms()-sqr(value_type(sum.first)/X.nelms())); }

//{unsecret}
//Summary: Average and standard deviation at once
//Return: A pair whose first value contains the mean, and whose second value contains the standard deviation
template<class G> inline pair<PROMOTE2(float,typename Vector<G>::const_value_type),PROMOTE2(float,typename Vector<G>::const_value_type)> mean_std_dev(const Vector<G> &X) { typedef PROMOTE2(float,typename Vector<G>::const_value_type) value_type; pair<value_type,value_type> sum = square_sum(X); return make_pair(sum.first/X.nelms(),sqrt(sum.second/X.nelms()-sqr(sum.first/X.nelms()))); }
//{unsecret}
template<class G> inline pair<PROMOTE2(float,typename Matrix<G>::const_value_type),PROMOTE2(float,typename Matrix<G>::const_value_type)> mean_std_dev(const Matrix<G> &X) { typedef PROMOTE2(float,typename Matrix<G>::const_value_type) value_type; pair<value_type,value_type> sum = square_sum(X); return make_pair(sum.first/X.nelms(),sqrt(sum.second/X.nelms()-sqr(sum.first/X.nelms()))); }

//{unsecret}
//Summary: Sum of the squared differences
template<class G1,class G2> inline typename promotion2_traits<typename Vector<G1>::const_value_type,typename Vector<G2>::const_value_type>::value_type square_error(const Vector<G1> &X, const Vector<G2> &Y) { return sum(ele_sqr(X-Y)); }
//{unsecret}
template<class G1,class G2> inline typename promotion2_traits<typename Matrix<G1>::const_value_type,typename Matrix<G2>::const_value_type>::value_type square_error(const Matrix<G1> &X, const Matrix<G2> &Y) { return sum(ele_sqr(X-Y)); }

//{unsecret}
//Summary: Average of the squared differences
template<class G1,class G2> inline typename promotion3_traits<float,typename Vector<G1>::const_value_type,typename Vector<G2>::const_value_type>::value_type mse(const Vector<G1> &X, const Vector<G2> &Y) { return mean(ele_sqr(X-Y)); }
//{unsecret}
template<class G1,class G2> inline typename promotion3_traits<float,typename Matrix<G1>::const_value_type,typename Matrix<G2>::const_value_type>::value_type mse(const Matrix<G1> &X, const Matrix<G2> &Y) { return mean(ele_sqr(X-Y)); }

//{unsecret}
//{noAutoLink}
//Summary: Initalizes an array with a function
//Remarks:
//  The function /f/ is called for each position of /X/
//
//    - for vectors
//
//          X[i]   = f(i);   
//
//    - for matrices
//
//          X(i,j) = f(i,j); 
//See: ^dense_vector_generate^, ^dense_matrix_generate^
template<class G, class F> Vector<G> &generate(Vector<G> &X, const F &f) { typedef typename Vector<G>::int_type int_type; for (int_type i=X.lower_bound(),imax=X.upper_bound(); i<=imax; ++i) X[i]=f(i); return X; }
//{unsecret}
template<class G, class F> Matrix<G> &generate(Matrix<G> &X, const F &f) { typedef typename Matrix<G>::int_type int_type; for (int_type i=X.row_lower_bound(), imax=X.row_upper_bound(); i<=imax; ++i) for (int_type j=X.col_lower_bound(), jmax=X.col_upper_bound(); j<=jmax; ++j) X(i,j)=f(j,i); return X; }

//{unsecret}
//{Group:Vectors functions}
//Summary: Creates an dense vector and initalizes it with a function
//Arguments:
//  n - size of the dense vector
//  f - the function to call
//Remarks:
//  The function /f/ is called for each position of /X/
//
//      X[i] = f(i);
//See: ^generate^, ^dense_matrix_generate^, ^DenseVector^
template<class F>             inline typename DenseVector<typename F::result_type>::self dense_vector_generate(int n, const F &f)    { typename DenseVector<typename F::result_type>::self X(n); return generate(X,f); }
//{unsecret}
template<class Arg,class Res> inline typename DenseVector<Res                    >::self dense_vector_generate(int n, Res (*f)(Arg)) { return dense_vector_generate(n,ptr_fun(f)); }

//{unsecret}
//{Group:Matrices functions}
//Summary: Creates an dense matrix and initalizes it with a function
//Arguments:
//  m - Height of the dense vector
//  n - Width of the dense vector
//  f - The function to call
//Remarks:
//  The function /f/ is called for each position of /X/
//
//      X(i,j) = f(i,j);
//See: ^generate^, ^dense_vector_generate^, ^DenseVector^
template<class F                        > inline typename DenseMatrix<typename F::result_type>::self dense_matrix_generate(int m, int n, const F &f)          { typename DenseMatrix<typename F::result_type>::self X(m,n); return generate(X,f); }
//{unsecret}
template<class Arg1,class Arg2,class Res> inline typename DenseMatrix<Res                    >::self dense_matrix_generate(int m, int n, Res (*f)(Arg1,Arg2)) { return dense_matrix_generate(m,n,ptr_fun(f)); }

    
template<class G> inline void affect(Vector<G> &X, const typename Vector<G>::value_type &v) { fill(X.begin(), X.end(), v); }
template<class G> inline void affect(Matrix<G> &X, const typename Matrix<G>::value_type &v) { fill(X.begin(), X.end(), v); }

template<class G1,class G2> inline void affect(Vector<G1> &X, const Vector<G2>  &Y) { typedef typename Vector<G1>::template const_iterator_rebind<Vector<G2> >::other const_iterator; copyn(X.nelms(), const_iterator(Y), X.begin()); }
template<class G1,class G2> inline void affect(Vector<G1> &X, const Array<1,G2> &Y) { typedef typename Vector<G1>::template const_iterator_rebind<Vector<G2> >::other const_iterator; copyn(X.nelms(), const_iterator(Y), X.begin()); }
//template<class G1,class G2> inline void affect(Vector<G1> &X, const Vector<G2>  &Y) { typedef typename Vector<G1>::template const_iterator_rebind<Vector<G2> >::other const_iterator; copy(const_iterator(Y), const_iterator(Y)+X.nelms(), X.begin()); }
//template<class G1,class G2> inline void affect(Vector<G1> &X, const Array<1,G2> &Y) { typedef typename Vector<G1>::template const_iterator_rebind<Vector<G2> >::other const_iterator; copy(const_iterator(Y), const_iterator(Y)+X.nelms(), X.begin()); }

template<class G1,class G2> inline void affect(Matrix<G1> &X, const Matrix<G2>  &Y) { typedef typename Matrix<G1>::template const_iterator_rebind<Matrix<G2> >::other const_iterator; copyn(X.nelms(), const_iterator(Y), X.begin()); }
template<class G1,class G2> inline void affect(Matrix<G1> &X, const Array<2,G2> &Y) { typedef typename Matrix<G1>::template const_iterator_rebind<Matrix<G2> >::other const_iterator; copyn(X.nelms(), const_iterator(Y), X.begin()); }
//template<class G1,class G2> inline void affect(Matrix<G1> &X, const Matrix<G2>  &Y) { typedef typename Matrix<G1>::template const_iterator_rebind<Matrix<G2> >::other const_iterator; copy(const_iterator(Y), const_iterator(Y)+X.nelms(), X.begin()); }
//template<class G1,class G2> inline void affect(Matrix<G1> &X, const Array<2,G2> &Y) { typedef typename Matrix<G1>::template const_iterator_rebind<Matrix<G2> >::other const_iterator; copy(const_iterator(Y), const_iterator(Y)+X.nelms(), X.begin()); }

//}

#endif

