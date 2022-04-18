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


#ifndef BLAS_H
#define BLAS_H

#include "array/matrix.h"

#ifndef BLAS_NAME
#define BLAS_NAME(name) cblas_ ## name
#endif
//RegExpHook BLAS_NAME\(([^)]*)\) cblas_\1


#ifndef CBLAS_ENUM_DEFINED_H
#define CBLAS_ENUM_DEFINED_H
enum CBLAS_ORDER     { CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE { CblasNoTrans =111, CblasTrans   =112, CblasConjTrans=113 };
enum CBLAS_UPLO      { CblasUpper   =121, CblasLower   =122 };
enum CBLAS_DIAG      { CblasNonUnit =131, CblasUnit    =132 };
enum CBLAS_SIDE      { CblasLeft    =141, CblasRight   =142 };
#endif


template<class T, int Copy=0> 
struct id_adaptor : public array_generator<T,Copy>
{
  typedef id_adaptor self;
  typedef array_generator<T,Copy> base;
  ARRAY_BASE_TYPES
  
  template<class T2,int C2=Copy> struct rebind { typedef id_adaptor<T2,C2> other; };
            
  using base::array;
  inline id_adaptor(const self &x) : base(static_cast<const base &>(x)) {}  
  inline id_adaptor(array_type &x) : base(x) {}
  
  template<class V2> inline V2 eval(const index_type &,const V2 &x) const { return x; } 
  template<class V2> inline V2 eval(int_type,int_type, const V2 &x) const { return x; } 
  
  inline void operator()(const index_type &i) { } 
  inline void operator()(const index_type &i,const value_type &x) { array()[i]=eval(i,x); } 
  template<int N> inline void operator()(const index_type &i,const Vector<simd_vector_generator<N,value_type> > &x) 
  {
    value_type *p=&array()[i]; 
    //_mm_prefetch((char*)p,_MM_HINT_T0);
    storeu(p,eval(i,x)); 
  } 

  typedef id_adaptor<typename SubArray<T>::self,1> sub_type;
  inline sub_type sub(index_type i, size_type n) { typename SubArray<T>::self A=::sub(array (),i,n); return sub_type(A); }
};

template<class Adapt> 
struct add_adaptor : public array_traits<typename Adapt::array_type>
{
  typedef add_adaptor self;
  typedef array_traits<typename Adapt::array_type> base;
  typedef Adapt adaptor_type;
  ARRAY_BASE_TYPES  
    
  template<class Adapt2> struct rebind { typedef add_adaptor<Adapt2> other; };
    
  adaptor_type adapt; 
    
  inline add_adaptor(const self &x) : adapt(x.adapt) {}
  inline add_adaptor(const adaptor_type &a) : adapt(a) {}
  template<class A> inline add_adaptor(A &a) : adapt(a) {}
  template<class A,class B> inline add_adaptor(A &a,const B &b) : adapt(a,b) {}
  
  inline size_type         size () const { return adapt.size (); }
  inline array_type       &array()       { return adapt.array(); }
  inline const array_type &array() const { return adapt.array(); }
    
  inline void operator()(const index_type &i,const value_type &x) 
  { 
    reference r=array()[i];
    r+=adapt.eval(i,x); 
  } 
  template<int N> inline void operator()(const index_type &i,const Vector<simd_vector_generator<N,value_type> > &x) 
  {
    typename SimdVector<N,value_type>::self a;
    value_type *p=&array()[i]; 
    //_mm_prefetch((char*)p,_MM_HINT_T0);
    loadu(a,p); 
    a+=adapt.eval(i,x); 
    storeu(p,a);
  }

  typedef add_adaptor<typename Adapt::sub_type> sub_type;
  inline sub_type sub(index_type i, size_type n) { return sub_type(adapt.sub(i,n)); }  
};

template<class T, int Copy=0> 
struct neg_adaptor : public array_generator<T,Copy>
{
  typedef neg_adaptor self;
  typedef array_generator<T,Copy> base;
  ARRAY_BASE_TYPES
  
  template<class T2,int C2=Copy> struct rebind { typedef neg_adaptor<T2,C2> other; };
            
  using base::array;
  inline neg_adaptor(const self &x) : base(static_cast<const base &>(x)) {}  
  inline neg_adaptor(array_type &x) : base(x) {}

  template<class V2> inline V2 eval(const index_type &,const V2 &x) const { return -x; } 
  template<class V2> inline V2 eval(int_type,int_type, const V2 &x) const { return -x; } 

  inline void operator()(const index_type &i) { reference r=array()[i]; r=-r; } 
  inline void operator()(const index_type &i,const value_type &x) { array()[i]=eval(i,x); } 
  template<int N> inline void operator()(const index_type &i,const Vector<simd_vector_generator<N,value_type> > &x) { storeu(&array()[i],eval(i,x)); } 

  typedef neg_adaptor<typename SubArray<T>::self,1> sub_type;
  inline sub_type sub(index_type i, size_type n) { typename SubArray<T>::self A=::sub(array (),i,n); return sub_type(A); }
};

template<class T, int Copy=0> 
struct ax_adaptor : public array_generator<T,Copy>
{
  typedef ax_adaptor self;
  typedef array_generator<T,Copy> base;
  ARRAY_BASE_TYPES
  
  template<class T2,int C2=Copy> struct rebind { typedef neg_adaptor<T2,C2> other; };
  
  value_type val;
            
  using base::array;
  inline ax_adaptor(const self &x) : base(static_cast<const base &>(x)), val(x.val) {}
  inline ax_adaptor(array_type &x,const value_type &alpha) : base(x), val(alpha) {}

  template<class V2> inline V2 eval(const index_type &,const V2 &x) const { return val*x; } 
  template<class V2> inline V2 eval(int_type,int_type, const V2 &x) const { return val*x; } 
    
  inline void operator()(const index_type &i) { reference r=array()[i]; r=val*r; } 
  inline void operator()(const index_type &i,const value_type &x) { array()[i]=eval(i,x); } 
  template<int N> inline void operator()(const index_type &i,const Vector<simd_vector_generator<N,value_type> > &x) { storeu(&array()[i],eval(i,x)); } 
    
  typedef ax_adaptor<typename SubArray<T>::self,1> sub_type;
  inline sub_type sub(index_type i, size_type n) { typename SubArray<T>::self A=::sub(array (),i,n); return sub_type(A,val); }  
};


template<class Adapt> 
struct sub_adaptor : public array_traits<typename Adapt::array_type>
{
  typedef sub_adaptor self;
  typedef array_traits<typename Adapt::array_type> base;
  typedef Adapt adaptor_type;
  ARRAY_BASE_TYPES  
    
  template<class Adapt2> struct rebind { typedef sub_adaptor<Adapt2> other; };
    
  adaptor_type adapt; 
    
  inline sub_adaptor(const self &x) : adapt(x.adapt) {}
  inline sub_adaptor(const adaptor_type &a) : adapt(a) {}
  template<class A> inline sub_adaptor(A &a) : adapt(a) {}
  template<class A,class B> inline sub_adaptor(A &a,const B &b) : adapt(a,b) {}

  inline size_type         size () const { return adapt.size (); }
  inline array_type       &array()       { return adapt.array(); }
  inline const array_type &array() const { return adapt.array(); }
  
  inline void operator()(const index_type &i,const value_type &x) { reference r=array()[i]; r-=adapt.eval(i,x); } 
  template<int N> inline void operator()(const index_type &i,const Vector<simd_vector_generator<N,value_type> > &x) { value_type *p=&array()[i]; typename SimdVector<N,value_type>::self a; loadu(a,p); a-=adapt.eval(i,x); storeu(p,a); }

  typedef sub_adaptor<typename Adapt::sub_type> sub_type;
  inline sub_type sub(index_type i, size_type n) { return sub_type(adapt.sub(i,n)); }  
};

template<class Adapt> 
struct scale_adaptor : public array_traits<typename Adapt::array_type>
{
  typedef scale_adaptor self;
  typedef array_traits<typename Adapt::array_type> base;
  typedef Adapt adaptor_type;
  ARRAY_BASE_TYPES  
    
  template<class Adapt2> struct rebind { typedef scale_adaptor<Adapt2> other; };
    
  adaptor_type adapt; 
  value_type val;
    
  inline scale_adaptor(const self &x) : adapt(x.adapt), val(x.val) {}
  inline scale_adaptor(const adaptor_type &a, const value_type &v) : adapt(a), val(v) {}
  template<class A> inline scale_adaptor(A &a, const value_type &v) : adapt(a), val(v) {}
  template<class A,class B> inline scale_adaptor(A &a,const B &b, const value_type &v) : adapt(a,b), val(v) {}

  inline size_type         size () const { return adapt.size (); }
  inline array_type       &array()       { return adapt.array(); }
  inline const array_type &array() const { return adapt.array(); }
    
  inline void operator()(const index_type &i) { reference r=array()[i]; r*=val; } 
  inline void operator()(const index_type &i,const value_type &x) { reference r=array()[i]; r*=val; r+=adapt.eval(i,x); } 
  template<int N> inline void operator()(const index_type &i,const Vector<simd_vector_generator<N,value_type> > &x) { value_type *p=&array()[i]; typename SimdVector<N,value_type>::self a; loadu(a,p); a*=val; a+=adapt.eval(i,x); storeu(p,a); }
  
  typedef scale_adaptor<typename Adapt::sub_type> sub_type;
  inline sub_type sub(index_type i, size_type n) { return sub_type(adapt.sub(i,n),val); }    
};


template<class T,int C=0> struct IdAdaptor  { typedef id_adaptor <T,C> self; };
template<class T,int C=0> struct NegAdaptor { typedef neg_adaptor<T,C> self; };
template<class T,int C=0> struct AxAdaptor  { typedef ax_adaptor <T,C> self; };

template<class Adapt> struct AddAdaptor   { typedef add_adaptor  <Adapt> self; };
template<class Adapt> struct SubAdaptor   { typedef sub_adaptor  <Adapt> self; };
template<class Adapt> struct ScaleAdaptor { typedef scale_adaptor<Adapt> self; };


struct id_adaptor_function
{
  template<class T> struct result_rebind { typedef typename IdAdaptor <T,0>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { return typename result_rebind<Vector<G> >::other(X); }
  template<class G> inline typename result_rebind<Matrix<G> >::other operator()(Matrix<G> &X) const { return typename result_rebind<Matrix<G> >::other(X); }
};
struct neg_adaptor_function
{
  template<class T> struct result_rebind { typedef typename NegAdaptor<T,0>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { return typename result_rebind<Vector<G> >::other(X); }
  template<class G> inline typename result_rebind<Matrix<G> >::other operator()(Matrix<G> &X) const { return typename result_rebind<Matrix<G> >::other(X); }
};
template<class V>
struct ax_adaptor_function
{
  typedef V value_type;
  value_type val;
  inline ax_adaptor_function(const value_type &v) : val(v) {}
  template<class T> struct result_rebind { typedef typename AxAdaptor <T,0>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { return typename result_rebind<Vector<G> >::other(X,typename Vector<G>::value_type(val)); }
  template<class G> inline typename result_rebind<Matrix<G> >::other operator()(Matrix<G> &X) const { return typename result_rebind<Matrix<G> >::other(X,typename Matrix<G>::value_type(val)); }
};

template<class Adapt>
struct add_adaptor_function
{
  typedef add_adaptor_function self;
  typedef Adapt adaptor_type;
  adaptor_type adapt;
  
  inline add_adaptor_function() : adapt() {}
  inline add_adaptor_function(const self &x) : adapt(x.adapt) {}  
  template<class A> inline add_adaptor_function(const A &a) : adapt(a) {}
  
  template<class T> struct result_rebind { typedef typename AddAdaptor<typename adaptor_type::template result_rebind<T>::other>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { return typename result_rebind<Vector<G> >::other(adapt(X)); }
  template<class G> inline typename result_rebind<Matrix<G> >::other operator()(Matrix<G> &X) const { return typename result_rebind<Matrix<G> >::other(adapt(X)); }
};

template<class Adapt>
struct sub_adaptor_function
{
  typedef sub_adaptor_function self;
  typedef Adapt adaptor_type;
  adaptor_type adapt;
  
  inline sub_adaptor_function() : adapt() {}
  inline sub_adaptor_function(const self &x) : adapt(x.adapt) {}
  template<class A> inline sub_adaptor_function(const A &a) : adapt(a) {}
  
  template<class T> struct result_rebind { typedef typename SubAdaptor<typename adaptor_type::template result_rebind<T>::other>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { return typename result_rebind<Vector<G> >::other(adapt(X)); }
  template<class G> inline typename result_rebind<Matrix<G> >::other operator()(Matrix<G> &X) const { return typename result_rebind<Matrix<G> >::other(adapt(X)); }
};

template<class Adapt,class V>
struct scale_adaptor_function
{
  typedef scale_adaptor_function self;
  typedef Adapt adaptor_type;
  typedef V value_type;
  adaptor_type adapt;
  value_type val;
  
  inline scale_adaptor_function(const value_type &v) : adapt(), val(v) {}
  inline scale_adaptor_function(const self &x) : adapt(x.adapt), val(x.val) {}
  template<class A> inline scale_adaptor_function(const A &a,const value_type &v) : adapt(a), val(v) {}
  
  template<class T> struct result_rebind { typedef typename ScaleAdaptor<typename adaptor_type::template result_rebind<T>::other>::self other; };
  template<class G> inline typename result_rebind<Vector<G> >::other operator()(Vector<G> &X) const { return typename result_rebind<Vector<G> >::other(adapt(X),typename Vector<G>::value_type(val)); }
  template<class G> inline typename result_rebind<Matrix<G> >::other operator()(Matrix<G> &X) const { return typename result_rebind<Matrix<G> >::other(adapt(X),typename Matrix<G>::value_type(val)); }
};

                  inline id_adaptor_function     id_adapt (          ) { return id_adaptor_function    ( ); }
                  inline neg_adaptor_function    neg_adapt(          ) { return neg_adaptor_function   ( ); }
template<class V> inline ax_adaptor_function <V> ax_adapt (const V &a) { return ax_adaptor_function <V>(a); }

                  inline add_adaptor_function  <id_adaptor_function     > add_adapt      (                     ) { return add_adaptor_function  <id_adaptor_function     >(    ); }
                  inline sub_adaptor_function  <id_adaptor_function     > sub_adapt      (                     ) { return sub_adaptor_function  <id_adaptor_function     >(    ); }
template<class V> inline scale_adaptor_function<id_adaptor_function   ,V> scale_add_adapt(const V &s           ) { return scale_adaptor_function<id_adaptor_function   ,V>(   s); }
template<class V> inline scale_adaptor_function<neg_adaptor_function  ,V> scale_sub_adapt(const V &s           ) { return scale_adaptor_function<neg_adaptor_function  ,V>(   s); }
template<class V> inline add_adaptor_function  <ax_adaptor_function<V>  > add_adapt      (           const V &a) { return add_adaptor_function  <ax_adaptor_function<V>  >( a  ); }
template<class V> inline add_adaptor_function  <ax_adaptor_function<V>  > sub_adapt      (           const V &a) { return add_adaptor_function  <ax_adaptor_function<V>  >(-a  ); }
template<class V> inline scale_adaptor_function<ax_adaptor_function<V>,V> scale_add_adapt(const V &s,const V &a) { return scale_adaptor_function<ax_adaptor_function<V>,V>( a,s); }
template<class V> inline scale_adaptor_function<ax_adaptor_function<V>,V> scale_sub_adapt(const V &s,const V &a) { return scale_adaptor_function<ax_adaptor_function<V>,V>(-a,s); }




template<class Arg1,class Arg2,class Res=PROMOTE2(Arg1,Arg2)>
struct dot_function : public binary_function<Arg1,Arg2,Res>
{
  typedef binary_function<Arg1,Arg2,Res> base;
  typedef typename base::first_argument_type  first_argument_type;
  typedef typename base::second_argument_type second_argument_type;
  typedef typename base::result_type          result_type;
  template<class T1,class T2> inline typename dot_function<T1,T2>::result_type operator()(const T1 &x, const T2 &y) { return x*y; }
  //inline result_type operator()(const first_argument_type &x, const second_argument_type &y) const { return x*y; }
  //inline void operator()(result_type &r, const first_argument_type &x, const second_argument_type &y) const { r=x*y; }
};

template<class Arg1,class Arg2,class Res=PROMOTE2(Arg1,Arg2)>
struct cdot_function : public binary_function<Arg1,Arg2,Res>
{
  typedef binary_function<Arg1,Arg2,Res> base;
  typedef typename base::first_argument_type  first_argument_type;
  typedef typename base::second_argument_type second_argument_type;
  typedef typename base::result_type          result_type;
  template<class T1,class T2> inline typename cdot_function<T1,T2>::result_type operator()(const T1 &x, const T2 &y) { return cmul(x,y); }
  //inline result_type operator()(const first_argument_type &x, const second_argument_type &y) const { return cmul(x,y); } // conj(x)*y;
  //inline void operator()(result_type &r, const first_argument_type &x, const second_argument_type &y) const { r=cmul(x,y); } // r=conj(x)*y;
};

#endif // BLAS_H
