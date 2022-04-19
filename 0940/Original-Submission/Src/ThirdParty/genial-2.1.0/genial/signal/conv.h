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

#ifndef CONV_H
#define CONV_H

#include "array/matrix.h"

//namespace genial
//{

//Group = Convolution

template<class A1,class A2=A1,class V=typename promotion2_traits<typename A1::value_type,typename A2::value_type>::value_type>
class conv_vector_generator : public two_arrays_value_generator<const A1,const A2,V>
{
  public:
    typedef conv_vector_generator self;
    typedef two_arrays_value_generator<const A1,const A2,V> base;
    TWO_ARRAYS_BASE_TYPES
    typedef typename base::first_const_iterator  first_const_iterator;
    typedef typename base::second_const_iterator second_const_iterator;
    typedef typename base::first_const_reverse_iterator  first_const_reverse_iterator;
    typedef typename base::second_const_reverse_iterator second_const_reverse_iterator;
    
  protected:
    using base::X;
    using base::Y;  
    
  public:
    inline conv_vector_generator(first_array_type &x, second_array_type &y) : base(x,y) {}

    inline const_reference operator[](const index_type &i) const 
    {
      index_type i0=i-lower_bound();
      const size_type n1=X.size(), n2=Y.size();
      first_const_iterator  it1=X.begin(), end1=it1+i0+1;
      second_const_iterator it2=Y.begin();
      assert(n2<=n1);

      value_type r;
      if (i0>=n2 && i0<n1) 
      { 
        it1+=i0+1-n2; it2+=n2-1; 
        value_type t = *it1 * *it2;
        r = inner_product_n(n2-1,++it1, second_const_reverse_iterator(it2), t);
      }
      else      
      if (i0<n2) 
      {
        it2+=i0;  
        value_type t= *it1 * *it2;
        r = inner_product_n(i0, ++it1, second_const_reverse_iterator(it2), t); 
      }
      else 
      { 
        it1+=i0+1-n2; it2+=n2-1; 
        value_type t= *it1 * *it2;
        int n = n1-i0-2+n2;
        r = inner_product_n(n, ++it1, second_const_reverse_iterator(it2), t); 
      }
      return r;
    }
    
    inline size_type lower_bound() const { return X.lower_bound()+Y.lower_bound(); }
    inline size_type size() const { return X.size()+Y.size()-1; }
    inline void resize(const size_type &d) { assert(d==size()); }
};

template<class T1,class T2=T1,class V=typename promotion2_traits<typename T1::value_type,typename T2::value_type>::value_type>
struct convVector  
{
  typedef T1 first_array_type;
  typedef T2 second_array_type;
  typedef V value_type; 
  typedef conv_vector_generator<first_array_type, second_array_type, value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
}; 

//{unsecret}
//Summary: Convolution
//Arguments:
//  X,Y - The arrays
//Return: An array representing the convolution
template<class G1,class G2> inline typename convVector<const Vector<G1>,const Vector<G2> >::self conv(const Vector<G1> &X, const Vector<G2> &Y)
{
  return typename convVector<const Vector<G1>,const Vector<G2> >::self(X,Y);
}



template<class A,class V=typename A::value_type>
class sum_cyclic_vector_generator : public array_value_generator<const A,V>
{
  public:
    typedef sum_cyclic_vector_generator self;
    typedef array_value_generator<const A,V> base;
    ARRAY_BASE_TYPES

  protected:
    size_type sz;
    using base::X;

  public:
    sum_cyclic_vector_generator(array_type &x, const size_type &s) : base(x), sz(s) { }

    const_reference operator[](const index_type &i) const 
    {    
      value_type sum(0);
      index_type j;
      for (j=i       ; j>=X.lower_bound(); j-=size()) sum+=X[j];
      for (j=i+size(); j<=X.upper_bound(); j+=size()) sum+=X[j];
      return sum;
    }

    size_type size() const { return sz; }
    index_type lower_bound() const { return 0; }
    void resize(const size_type &s) { sz=s; }
};

template<class A,class V=typename A::value_type>
struct sumCyclicVector  
{
  typedef A array_type;
  typedef V value_type; 
  typedef sum_cyclic_vector_generator<array_type, value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
}; 


template<class G>
inline typename sumCyclicVector<const Vector<G> >::self sum_cyclic(const Vector<G> &X, const typename Vector<G>::size_type &n)
{
  return typename sumCyclicVector<const Vector<G> >::self(X,n);
}


template<class T1,class T2=T1,class V=typename promotion2_traits<typename T1::value_type,typename T2::value_type>::value_type,int Copy=0>
class cyclic_conv_vector_generator : public two_arrays_value_generator<const T1,const T2,V,Copy>
{
  public:
    typedef cyclic_conv_vector_generator self;
    typedef two_arrays_value_generator<const T1,const T2,V,Copy> base;
    TWO_ARRAYS_BASE_TYPES
    
  protected:
    using base::X;
    using base::Y;

  public:
    cyclic_conv_vector_generator(first_array_type &x, second_array_type &y) : base(x,y) {}
    
    template<class T3,class T4,class V2,int C2> cyclic_conv_vector_generator(const        cyclic_conv_vector_generator<T3,T4,V2,C2>   &x) : base(x.            first_array(),x.            second_array()) {}
    template<class T3,class T4,class V2,int C2> cyclic_conv_vector_generator(const Vector<cyclic_conv_vector_generator<T3,T4,V2,C2> > &x) : base(x.generator().first_array(),x.generator().second_array()) {}
    
    const_reference operator[](const index_type &i) const { return sum_cyclic(conv(X,Y),size())[i]; }
      
    size_type  size       () const { return X.size(); }
    index_type lower_bound() const { return 0; } 
};

template<class T1,class T2=T1,class V=typename promotion2_traits<typename T1::value_type,typename T2::value_type>::value_type,int Copy=0>
struct cyclicConvVector  
{
  typedef T1 first_array_type;
  typedef T2 second_array_type;
  typedef V value_type; 
  typedef cyclic_conv_vector_generator<first_array_type, second_array_type, value_type,Copy> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
}; 

template<class G1, class G2>
inline typename cyclicConvVector<const Vector<G1>,const Vector<G2> >::self cyclic_conv(const Vector<G1> &X, const Vector<G2> &Y)
{
  return typename cyclicConvVector<const Vector<G1>,const Vector<G2> >::self(X,Y);
}

template<class Arg1, class Arg2>
struct cyclic_conv_function : public binary_value_function<const Arg1,const Arg2,typename cyclicConvVector<const Arg1,const Arg2,PROMOTE2(typename Arg1::value_type,typename Arg2::value_type),1>::self>
{
  typedef binary_value_function<const Arg1,const Arg2,typename cyclicConvVector<const Arg1,const Arg2,PROMOTE2(typename Arg1::value_type,typename Arg2::value_type),1>::self> base;
  BINARY_FUNCTION_BASE_TYPES
  const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return cyclic_conv(x,y); }
};




#if defined(MMX) || defined (SSE2)

template<int S,class G1,class G2,class G3,int N>
void aux_simd_conv7(const Vector<G1> &X, const Vector<G2> &Y, Vector<G3> &Z, Vector<simd_vector_generator<N,unsigned char> >, short, Vector<simd_vector_generator<N,unsigned char> >) 
{
  int n=X.size();
  assert(n%8==0 && n>=8);
  assert(Y.size()==7);

  typename TinyVector<N/2,short>::self y0(Y[0]),y1(Y[1]),y2(Y[2]),y3(Y[3]),y4(Y[4]),y5(Y[5]),y6(Y[6]);
  typename TinyVector<N,unsigned char>::self x;
  typename TinyVector<N/2,short>::self x0l,x1l,x2l,x3l,x4l,x5l(0),x6l(0),x7l(0);
  typename TinyVector<N/2,short>::self x0h,x1h,x2h,x3h,x4h,x5h(0),x6h(0),x7h(0);
  typename TinyVector<N/2,short>::self zl, zh;
  
  x=X[0]; x0l=unpacklo(x); x0h=unpackhi(x);
  x=X[1]; x1l=unpacklo(x); x1h=unpackhi(x);
  x=X[2]; x2l=unpacklo(x); x2h=unpackhi(x);
  int i=0;
  for (int imax=n-8; i<imax; i+=8)
  {    
    x=X[i+3]; x3l=unpacklo(x); x3h=unpackhi(x);
    zl =mullo(y6,x5l); zl+=mullo(y5,x6l); zl+=mullo(y4,x7l); zl+=mullo(y3,x0l); zl+=mullo(y2,x1l); zl+=mullo(y1,x2l); zl+=mullo(y0,x3l);
    zh =mullo(y6,x5h); zh+=mullo(y5,x6h); zh+=mullo(y4,x7h); zh+=mullo(y3,x0h); zh+=mullo(y2,x1h); zh+=mullo(y1,x2h); zh+=mullo(y0,x3h);
    Z[i  ]=packus(shiftra<S>(zl),shiftra<S>(zh));
 
    x=X[i+4]; x4l=unpacklo(x); x4h=unpackhi(x);
    zl =mullo(y6,x6l); zl+=mullo(y5,x7l); zl+=mullo(y4,x0l); zl+=mullo(y3,x1l); zl+=mullo(y2,x2l); zl+=mullo(y1,x3l); zl+=mullo(y0,x4l);
    zh =mullo(y6,x6h); zh+=mullo(y5,x7h); zh+=mullo(y4,x0h); zh+=mullo(y3,x1h); zh+=mullo(y2,x2h); zh+=mullo(y1,x3h); zh+=mullo(y0,x4h);
    Z[i+1]=packus(shiftra<S>(zl),shiftra<S>(zh));
   
    x=X[i+5]; x5l=unpacklo(x); x5h=unpackhi(x);
    zl =mullo(y6,x7l); zl+=mullo(y5,x0l); zl+=mullo(y4,x1l); zl+=mullo(y3,x2l); zl+=mullo(y2,x3l); zl+=mullo(y1,x4l); zl+=mullo(y0,x5l);
    zh =mullo(y6,x7h); zh+=mullo(y5,x0h); zh+=mullo(y4,x1h); zh+=mullo(y3,x2h); zh+=mullo(y2,x3h); zh+=mullo(y1,x4h); zh+=mullo(y0,x5h);
    Z[i+2]=packus(shiftra<S>(zl),shiftra<S>(zh));
 
    x=X[i+6]; x6l=unpacklo(x); x6h=unpackhi(x);
    zl =mullo(y6,x0l); zl+=mullo(y5,x1l); zl+=mullo(y4,x2l); zl+=mullo(y3,x3l); zl+=mullo(y2,x4l); zl+=mullo(y1,x5l); zl+=mullo(y0,x6l);
    zh =mullo(y6,x0h); zh+=mullo(y5,x1h); zh+=mullo(y4,x2h); zh+=mullo(y3,x3h); zh+=mullo(y2,x4h); zh+=mullo(y1,x5h); zh+=mullo(y0,x6h);
    Z[i+3]=packus(shiftra<S>(zl),shiftra<S>(zh));

    x=X[i+7]; x7l=unpacklo(x); x7h=unpackhi(x);
    zl =mullo(y6,x1l); zl+=mullo(y5,x2l); zl+=mullo(y4,x3l); zl+=mullo(y3,x4l); zl+=mullo(y2,x5l); zl+=mullo(y1,x6l); zl+=mullo(y0,x7l);
    zh =mullo(y6,x1h); zh+=mullo(y5,x2h); zh+=mullo(y4,x3h); zh+=mullo(y3,x4h); zh+=mullo(y2,x5h); zh+=mullo(y1,x6h); zh+=mullo(y0,x7h);
    Z[i+4]=packus(shiftra<S>(zl),shiftra<S>(zh));

    x=X[i+8]; x0l=unpacklo(x); x0h=unpackhi(x);
    zl =mullo(y6,x2l); zl+=mullo(y5,x3l); zl+=mullo(y4,x4l); zl+=mullo(y3,x5l); zl+=mullo(y2,x6l); zl+=mullo(y1,x7l); zl+=mullo(y0,x0l);
    zh =mullo(y6,x2h); zh+=mullo(y5,x3h); zh+=mullo(y4,x4h); zh+=mullo(y3,x5h); zh+=mullo(y2,x6h); zh+=mullo(y1,x7h); zh+=mullo(y0,x0h);
    Z[i+5]=packus(shiftra<S>(zl),shiftra<S>(zh));

    x=X[i+9]; x1l=unpacklo(x); x1h=unpackhi(x);
    zl =mullo(y6,x3l); zl+=mullo(y5,x4l); zl+=mullo(y4,x5l); zl+=mullo(y3,x6l); zl+=mullo(y2,x7l); zl+=mullo(y1,x0l); zl+=mullo(y0,x1l);
    zh =mullo(y6,x3h); zh+=mullo(y5,x4h); zh+=mullo(y4,x5h); zh+=mullo(y3,x6h); zh+=mullo(y2,x7h); zh+=mullo(y1,x0h); zh+=mullo(y0,x1h);
    Z[i+6]=packus(shiftra<S>(zl),shiftra<S>(zh));

    x=X[i+10]; x2l=unpacklo(x); x2h=unpackhi(x);
    zl =mullo(y6,x4l); zl+=mullo(y5,x5l); zl+=mullo(y4,x6l); zl+=mullo(y3,x7l); zl+=mullo(y2,x0l); zl+=mullo(y1,x1l); zl+=mullo(y0,x2l);
    zh =mullo(y6,x4h); zh+=mullo(y5,x5h); zh+=mullo(y4,x6h); zh+=mullo(y3,x7h); zh+=mullo(y2,x0h); zh+=mullo(y1,x1h); zh+=mullo(y0,x2h);
    Z[i+7]=packus(shiftra<S>(zl),shiftra<S>(zh));
  }
  
  x=X[i+3]; x3l=unpacklo(x); x3h=unpackhi(x);
  zl =mullo(y6,x5l); zl+=mullo(y5,x6l); zl+=mullo(y4,x7l); zl+=mullo(y3,x0l); zl+=mullo(y2,x1l); zl+=mullo(y1,x2l); zl+=mullo(y0,x3l);
  zh =mullo(y6,x5h); zh+=mullo(y5,x6h); zh+=mullo(y4,x7h); zh+=mullo(y3,x0h); zh+=mullo(y2,x1h); zh+=mullo(y1,x2h); zh+=mullo(y0,x3h);
  Z[i  ]=packus(shiftra<S>(zl),shiftra<S>(zh));

  x=X[i+4]; x4l=unpacklo(x); x4h=unpackhi(x);
  zl =mullo(y6,x6l); zl+=mullo(y5,x7l); zl+=mullo(y4,x0l); zl+=mullo(y3,x1l); zl+=mullo(y2,x2l); zl+=mullo(y1,x3l); zl+=mullo(y0,x4l);
  zh =mullo(y6,x6h); zh+=mullo(y5,x7h); zh+=mullo(y4,x0h); zh+=mullo(y3,x1h); zh+=mullo(y2,x2h); zh+=mullo(y1,x3h); zh+=mullo(y0,x4h);
  Z[i+1]=packus(shiftra<S>(zl),shiftra<S>(zh));
  
  x=X[i+5]; x5l=unpacklo(x); x5h=unpackhi(x);
  zl =mullo(y6,x7l); zl+=mullo(y5,x0l); zl+=mullo(y4,x1l); zl+=mullo(y3,x2l); zl+=mullo(y2,x3l); zl+=mullo(y1,x4l); zl+=mullo(y0,x5l);
  zh =mullo(y6,x7h); zh+=mullo(y5,x0h); zh+=mullo(y4,x1h); zh+=mullo(y3,x2h); zh+=mullo(y2,x3h); zh+=mullo(y1,x4h); zh+=mullo(y0,x5h);
  Z[i+2]=packus(shiftra<S>(zl),shiftra<S>(zh));

  x=X[i+6]; x6l=unpacklo(x); x6h=unpackhi(x);
  zl =mullo(y6,x0l); zl+=mullo(y5,x1l); zl+=mullo(y4,x2l); zl+=mullo(y3,x3l); zl+=mullo(y2,x4l); zl+=mullo(y1,x5l); zl+=mullo(y0,x6l);
  zh =mullo(y6,x0h); zh+=mullo(y5,x1h); zh+=mullo(y4,x2h); zh+=mullo(y3,x3h); zh+=mullo(y2,x4h); zh+=mullo(y1,x5h); zh+=mullo(y0,x6h);
  Z[i+3]=packus(shiftra<S>(zl),shiftra<S>(zh));

  x=X[i+7]; x7l=unpacklo(x); x7h=unpackhi(x);
  zl =mullo(y6,x1l); zl+=mullo(y5,x2l); zl+=mullo(y4,x3l); zl+=mullo(y3,x4l); zl+=mullo(y2,x5l); zl+=mullo(y1,x6l); zl+=mullo(y0,x7l);
  zh =mullo(y6,x1h); zh+=mullo(y5,x2h); zh+=mullo(y4,x3h); zh+=mullo(y3,x4h); zh+=mullo(y2,x5h); zh+=mullo(y1,x6h); zh+=mullo(y0,x7h);
  Z[i+4]=packus(shiftra<S>(zl),shiftra<S>(zh));

  zl =mullo(y6,x2l); zl+=mullo(y5,x3l); zl+=mullo(y4,x4l); zl+=mullo(y3,x5l); zl+=mullo(y2,x6l); zl+=mullo(y1,x7l);
  zh =mullo(y6,x2h); zh+=mullo(y5,x3h); zh+=mullo(y4,x4h); zh+=mullo(y3,x5h); zh+=mullo(y2,x6h); zh+=mullo(y1,x7h);
  Z[i+5]=packus(shiftra<S>(zl),shiftra<S>(zh));

  zl =mullo(y6,x3l); zl+=mullo(y5,x4l); zl+=mullo(y4,x5l); zl+=mullo(y3,x6l); zl+=mullo(y2,x7l);
  zh =mullo(y6,x3h); zh+=mullo(y5,x4h); zh+=mullo(y4,x5h); zh+=mullo(y3,x6h); zh+=mullo(y2,x7h);
  Z[i+6]=packus(shiftra<S>(zl),shiftra<S>(zh));

  zl =mullo(y6,x4l); zl+=mullo(y5,x5l); zl+=mullo(y4,x6l); zl+=mullo(y3,x7l);
  zh =mullo(y6,x4h); zh+=mullo(y5,x5h); zh+=mullo(y4,x6h); zh+=mullo(y3,x7h);
  Z[i+7]=packus(shiftra<S>(zl),shiftra<S>(zh));
}

template<int S,class G1,class G2,class G3,int N>
void aux_simd_sample2_conv7(const Vector<G1> &X, const Vector<G2> &Y, Vector<G3> &Z, Vector<simd_vector_generator<N,unsigned char> >, short, Vector<simd_vector_generator<N,unsigned char> >) 
{
  int n=X.size();
  assert(n%8==0 && n>=8);
  assert(Y.size()==7);

  typename TinyVector<N/2,short>::self y0(Y[0]),y1(Y[1]),y2(Y[2]),y3(Y[3]),y4(Y[4]),y5(Y[5]),y6(Y[6]);
  typename TinyVector<N,unsigned char>::self x;
  typename TinyVector<N/2,short>::self x0l, x1l, x2l, x3l, x4l, x5l, x6l, x7l;
  typename TinyVector<N/2,short>::self x0h, x1h, x2h, x3h, x4h, x5h, x6h, x7h;
  typename TinyVector<N/2,short>::self zl, zh;
  
  x=X[4]; x4l=0; x4h=0;
  x=X[5]; x5l=0; x5h=0;
  x=X[6]; x6l=0; x6h=0;
  x=X[7]; x7l=0; x7h=0;
  x=X[0]; x0l=unpacklo(x); x0h=unpackhi(x);
  x=X[1]; x1l=unpacklo(x); x1h=unpackhi(x);
  int i=0, j=0;
  for (int imax=n/2-4; i<imax; i+=4, j+=8)
  {    
    x=X[j+2]; x2l=unpacklo(x); x2h=unpackhi(x);
    x=X[j+3]; x3l=unpacklo(x); x3h=unpackhi(x);
    zl =mullo(y6,x5l); zl+=mullo(y5,x6l); zl+=mullo(y4,x7l); zl+=mullo(y3,x0l); zl+=mullo(y2,x1l); zl+=mullo(y1,x2l); zl+=mullo(y0,x3l);
    zh =mullo(y6,x5h); zh+=mullo(y5,x6h); zh+=mullo(y4,x7h); zh+=mullo(y3,x0h); zh+=mullo(y2,x1h); zh+=mullo(y1,x2h); zh+=mullo(y0,x3h);
    Z[i]=packus(shiftra<S>(zl),shiftra<S>(zh));
    
    x=X[j+4]; x4l=unpacklo(x); x4h=unpackhi(x);
    x=X[j+5]; x5l=unpacklo(x); x5h=unpackhi(x);
    zl =mullo(y6,x7l); zl+=mullo(y5,x0l); zl+=mullo(y4,x1l); zl+=mullo(y3,x2l); zl+=mullo(y2,x3l); zl+=mullo(y1,x4l); zl+=mullo(y0,x5l);
    zh =mullo(y6,x7h); zh+=mullo(y5,x0h); zh+=mullo(y4,x1h); zh+=mullo(y3,x2h); zh+=mullo(y2,x3h); zh+=mullo(y1,x4h); zh+=mullo(y0,x5h);
    Z[i+1]=packus(shiftra<S>(zl),shiftra<S>(zh));
    
    x=X[j+6]; x6l=unpacklo(x); x6h=unpackhi(x);
    x=X[j+7]; x7l=unpacklo(x); x7h=unpackhi(x);
    zl =mullo(y6,x1l); zl+=mullo(y5,x2l); zl+=mullo(y4,x3l); zl+=mullo(y3,x4l); zl+=mullo(y2,x5l); zl+=mullo(y1,x6l); zl+=mullo(y0,x7l);
    zh =mullo(y6,x1h); zh+=mullo(y5,x2h); zh+=mullo(y4,x3h); zh+=mullo(y3,x4h); zh+=mullo(y2,x5h); zh+=mullo(y1,x6h); zh+=mullo(y0,x7h);
    Z[i+2]=packus(shiftra<S>(zl),shiftra<S>(zh));
    
    x=X[j+8]; x0l=unpacklo(x); x0h=unpackhi(x);
    x=X[j+9]; x1l=unpacklo(x); x1h=unpackhi(x);
    zl =mullo(y6,x3l); zl+=mullo(y5,x4l); zl+=mullo(y4,x5l); zl+=mullo(y3,x6l); zl+=mullo(y2,x7l); zl+=mullo(y1,x0l); zl+=mullo(y0,x1l);
    zh =mullo(y6,x3h); zh+=mullo(y5,x4h); zh+=mullo(y4,x5h); zh+=mullo(y3,x6h); zh+=mullo(y2,x7h); zh+=mullo(y1,x0h); zh+=mullo(y0,x1h);
    Z[i+3]=packus(shiftra<S>(zl),shiftra<S>(zh));
  }
  
  x=X[j+2]; x2l=unpacklo(x); x2h=unpackhi(x);
  x=X[j+3]; x3l=unpacklo(x); x3h=unpackhi(x);
  zl =mullo(y6,x5l); zl+=mullo(y5,x6l); zl+=mullo(y4,x7l); zl+=mullo(y3,x0l); zl+=mullo(y2,x1l); zl+=mullo(y1,x2l); zl+=mullo(y0,x3l);
  zh =mullo(y6,x5h); zh+=mullo(y5,x6h); zh+=mullo(y4,x7h); zh+=mullo(y3,x0h); zh+=mullo(y2,x1h); zh+=mullo(y1,x2h); zh+=mullo(y0,x3h);
  Z[i]=packus(shiftra<S>(zl),shiftra<S>(zh));
  
  x=X[j+4]; x4l=unpacklo(x); x4h=unpackhi(x);
  x=X[j+5]; x5l=unpacklo(x); x5h=unpackhi(x);
  zl =mullo(y6,x7l); zl+=mullo(y5,x0l); zl+=mullo(y4,x1l); zl+=mullo(y3,x2l); zl+=mullo(y2,x3l); zl+=mullo(y1,x4l); zl+=mullo(y0,x5l);
  zh =mullo(y6,x7h); zh+=mullo(y5,x0h); zh+=mullo(y4,x1h); zh+=mullo(y3,x2h); zh+=mullo(y2,x3h); zh+=mullo(y1,x4h); zh+=mullo(y0,x5h);
  Z[i+1]=packus(shiftra<S>(zl),shiftra<S>(zh));
  
  x=X[j+6]; x6l=unpacklo(x); x6h=unpackhi(x);
  x=X[j+7]; x7l=unpacklo(x); x7h=unpackhi(x);
  zl =mullo(y6,x1l); zl+=mullo(y5,x2l); zl+=mullo(y4,x3l); zl+=mullo(y3,x4l); zl+=mullo(y2,x5l); zl+=mullo(y1,x6l); zl+=mullo(y0,x7l);
  zh =mullo(y6,x1h); zh+=mullo(y5,x2h); zh+=mullo(y4,x3h); zh+=mullo(y3,x4h); zh+=mullo(y2,x5h); zh+=mullo(y1,x6h); zh+=mullo(y0,x7h);
  Z[i+2]=packus(shiftra<S>(zl),shiftra<S>(zh));
  
  zl =mullo(y6,x3l); zl+=mullo(y5,x4l); zl+=mullo(y4,x5l); zl+=mullo(y3,x6l); zl+=mullo(y2,x7l);
  zh =mullo(y6,x3h); zh+=mullo(y5,x4h); zh+=mullo(y4,x5h); zh+=mullo(y3,x6h); zh+=mullo(y2,x7h);
  Z[i+3]=packus(shiftra<S>(zl),shiftra<S>(zh));
}
#endif //defined(MMX) || defined (SSE2)


#if defined(SSE2)

template<int S>
void sse2_conv7(unsigned char *p, const unsigned char *p1, const short *p2, int n) 
{
  assert(n>=16 && n%16==0);
  SimdVector<16,unsigned char>::self a;
  SimdVector< 8,short        >::self x2(0),x0,x1, y, zl,zh;
  SimdVector< 4,int          >::self z0, z1;
  loadr(y,p2); y=shiftl<1>(y); y=shiftr<1>(y); load(a,p1); x0=unpacklo(a); x1=unpackhi(a);
  for(n-=16; n>0; n-=16, p+=16)
	{
    z0=sum( madd(shiftl<4>(x2)+shiftr<4>(x0),y), madd(shiftl<5>(x2)+shiftr<3>(x0),y), madd(shiftl<6>(x2)+shiftr<2>(x0),y), madd(shiftl<7>(x2)+shiftr<1>(x0),y) ); 
    z1=sum( madd(         (x0)              ,y), madd(shiftl<1>(x0)+shiftr<7>(x1),y), madd(shiftl<2>(x0)+shiftr<6>(x1),y), madd(shiftl<3>(x0)+shiftr<5>(x1),y) ); 
    zl=packs(shiftra<S>(z0),shiftra<S>(z1));
    z0=sum( madd(shiftl<4>(x0)+shiftr<4>(x1),y), madd(shiftl<5>(x0)+shiftr<3>(x1),y), madd(shiftl<6>(x0)+shiftr<2>(x1),y), madd(shiftl<7>(x0)+shiftr<1>(x1),y) ); x2=x1; load(a,p1+=16); x0=unpacklo(a); x1=unpackhi(a);
    z1=sum( madd(         (x2)              ,y), madd(shiftl<1>(x2)+shiftr<7>(x0),y), madd(shiftl<2>(x2)+shiftr<6>(x0),y), madd(shiftl<3>(x2)+shiftr<5>(x0),y) ); zh=packs(shiftra<S>(z0),shiftra<S>(z1));
	  store(p,packus(zl,zh));
	}
  z0=sum( madd(shiftl<4>(x2)+shiftr<4>(x0),y), madd(shiftl<5>(x2)+shiftr<3>(x0),y), madd(shiftl<6>(x2)+shiftr<2>(x0),y), madd(shiftl<7>(x2)+shiftr<1>(x0),y) );
  z1=sum( madd(         (x0)              ,y), madd(shiftl<1>(x0)+shiftr<7>(x1),y), madd(shiftl<2>(x0)+shiftr<6>(x1),y), madd(shiftl<3>(x0)+shiftr<5>(x1),y) ); zl=packs(shiftra<S>(z0),shiftra<S>(z1));
  z0=sum( madd(shiftl<4>(x0)+shiftr<4>(x1),y), madd(shiftl<5>(x0)+shiftr<3>(x1),y), madd(shiftl<6>(x0)+shiftr<2>(x1),y), madd(shiftl<7>(x0)+shiftr<1>(x1),y) ); x2=x1;
  z1=sum( madd(         (x2)              ,y), madd(shiftl<1>(x2)              ,y), madd(shiftl<2>(x2)              ,y), madd(shiftl<3>(x2)              ,y) ); zh=packs(shiftra<S>(z0),shiftra<S>(z1));
	store(p,packus(zl,zh));
}

template<int S>
void sse2_sample2_conv7(unsigned char *p, const unsigned char *p1, const short *p2, int n) 
{
  assert(n>=32 && n%32==0);
  SimdVector<16,unsigned char>::self a;
  SimdVector< 8,short        >::self x0,x1,x2,x3(0), y, zl,zh;
  SimdVector< 4,int          >::self z0, z1;
  loadr(y,p2); y=shiftl<1>(y); y=shiftr<1>(y); load(a,p1); x0=unpacklo(a); x1=unpackhi(a);
  for(n-=32; n>0; n-=32, p+=16)
	{
    z0=sum( madd(shiftl<4>(x3)+shiftr<4>(x0),y), madd(shiftl<6>(x3)+shiftr<2>(x0),y), madd(x0,y), madd(shiftl<2>(x0)+shiftr<6>(x1),y) ); load(a,p1+=16); x2=unpacklo(a); x3=unpackhi(a);
    z1=sum( madd(shiftl<4>(x0)+shiftr<4>(x1),y), madd(shiftl<6>(x0)+shiftr<2>(x1),y), madd(x1,y), madd(shiftl<2>(x1)+shiftr<6>(x2),y) ); zl=packs(shiftra<S>(z0),shiftra<S>(z1));
    z0=sum( madd(shiftl<4>(x1)+shiftr<4>(x2),y), madd(shiftl<6>(x1)+shiftr<2>(x2),y), madd(x2,y), madd(shiftl<2>(x2)+shiftr<6>(x3),y) ); load(a,p1+=16); x0=unpacklo(a); x1=unpackhi(a);
    z1=sum( madd(shiftl<4>(x2)+shiftr<4>(x3),y), madd(shiftl<6>(x2)+shiftr<2>(x3),y), madd(x3,y), madd(shiftl<2>(x3)+shiftr<6>(x0),y) ); zh=packs(shiftra<S>(z0),shiftra<S>(z1));
	  store(p,packus(zl,zh));
	}
  z0=sum( madd(shiftl<4>(x3)+shiftr<4>(x0),y), madd(shiftl<6>(x3)+shiftr<2>(x0),y), madd(x0,y), madd(shiftl<2>(x0)+shiftr<6>(x1),y) ); load(a,p1+=16); x2=unpacklo(a); x3=unpackhi(a);
  z1=sum( madd(shiftl<4>(x0)+shiftr<4>(x1),y), madd(shiftl<6>(x0)+shiftr<2>(x1),y), madd(x1,y), madd(shiftl<2>(x1)+shiftr<6>(x2),y) ); zl=packs(shiftra<S>(z0),shiftra<S>(z1));
  z0=sum( madd(shiftl<4>(x1)+shiftr<4>(x2),y), madd(shiftl<6>(x1)+shiftr<2>(x2),y), madd(x2,y), madd(shiftl<2>(x2)+shiftr<6>(x3),y) ); 
  z1=sum( madd(shiftl<4>(x2)+shiftr<4>(x3),y), madd(shiftl<6>(x2)+shiftr<2>(x3),y), madd(x3,y), madd(shiftl<2>(x3)              ,y) ); zh=packs(shiftra<S>(z0),shiftra<S>(z1));
	store(p,packus(zl,zh));
}

template<int S,class G1,class V,class G2> void aux_conv        (const Vector<G1> &X, const Vector<tiny_vector_generator<7,V> > &Y, Vector<G2> &Z,unsigned char,short,unsigned char) { return sse2_conv7<S>(&*Z.begin(), &*X.begin(), &*Y.begin(), X.size()); }
template<int S,class G1,class V,class G2> void aux_sample2_conv(const Vector<G1> &X, const Vector<tiny_vector_generator<7,V> > &Y, Vector<G2> &Z,unsigned char,short,unsigned char) { return sse2_sample2_conv7<S>(&*Z.begin(), &*X.begin(), &*Y.begin(), X.size()); }

template<int S,class G1,class V,class G2> void aux_conv        (const Vector<G1> &X, const Vector<tiny_vector_generator<7,V> > &Y, Vector<G2> &Z, const TinyVector<16,unsigned char>::self &, short, const TinyVector<16,unsigned char>::self &) { return aux_simd_conv7        <S>(X,Y,Z,typename Vector<G1>::value_type(),V(),typename Vector<G2>::value_type()); }
template<int S,class G1,class V,class G2> void aux_sample2_conv(const Vector<G1> &X, const Vector<tiny_vector_generator<7,V> > &Y, Vector<G2> &Z, const TinyVector<16,unsigned char>::self &, short, const TinyVector<16,unsigned char>::self &) { return aux_simd_sample2_conv7<S>(X,Y,Z,typename Vector<G1>::value_type(),V(),typename Vector<G2>::value_type()); }

#endif //SSE2


template<int S,class G1,class G2,class G3>
void conv(const Vector<G1> &X, const Vector<G2> &Y, Vector<G3> &Z) 
{
  typename ArrayData<Vector<G3> >::self Z2 = data(Z);
  return aux_conv<S>(data(X),Y,Z2,typename Vector<G1>::value_type(),typename Vector<G2>::value_type(),typename Vector<G3>::value_type());
}

template<int S,class G1,class G2,class G3>
void sample2_conv(const Vector<G1> &X, const Vector<G2> &Y, Vector<G3> &Z) 
{
  typename ArrayData<Vector<G3> >::self Z2 = data(Z);
  return aux_sample2_conv<S>(data(X),Y,Z2,typename Vector<G1>::value_type(),typename Vector<G2>::value_type(),typename Vector<G3>::value_type());
}


template<int S,class G1,class G2,class G3>
void conv(const Matrix<G1> &X, const Vector<G2> &Y, Matrix<G3> &Z) 
{
  typedef const Matrix<G1> first_array_type;
  typedef Matrix<G3> third_array_type;
  typedef typename third_array_type::value_type value_type;
  typedef typename DenseMatrix<value_type>::self array_type;
  
  array_type T(X.nrows(),X.ncols());

  
  typename RowSplit<third_array_type>::self T1 = rows(T);
  for (int i=0, imax=T.nrows(); i<imax; ++i)
  { 
    typename RowSplit<third_array_type>::self::reference T2 = T1[i];  
    conv<S>(rows(X)[i],Y,T2);
  }
  
  typename SimdBlock<array_type>::self Z1 = simd_block(Z);
  typename ColSplit<typename SimdBlock<array_type>::self>::self Z2 = cols(Z1);
  for (int j=0, jmax=Z1.ncols(); j<jmax; ++j)
  {
    typename ColSplit<typename SimdBlock<array_type>::self>::self::reference Z3 = Z2[j];
    conv<S>(cols(simd_block(T))[j],Y,Z3);
  }

    
  //typename SimdBlock<array_type>::self T1 = simd_block(T);
  //typename ColSplit<typename SimdBlock<array_type>::self>::self T2 = cols(T1);
  //for (int j=0, jmax=T1.ncols(); j<jmax; ++j)
  //{
  //  typename ColSplit<typename SimdBlock<array_type>::self>::self::reference T3 = T2[j];
  //  conv<S>(cols(simd_block(X))[j],Y,T3);
  //}
  //
  //typename RowSplit<third_array_type>::self Z1 = rows(Z);
  //for (int i=0, imax=T.nrows(); i<imax; ++i)
  //{ 
  //  typename RowSplit<third_array_type>::self::reference Z2 = Z1[i];  
  //  conv<S>(rows(T)[i],Y,Z2);
  //}
}

template<int S,class G1,class G2,class G3>
void sample2_conv(const Matrix<G1> &X, const Vector<G2> &Y, Matrix<G3> &Z) 
{
  typedef const Matrix<G1> first_array_type;
  typedef Matrix<G3> third_array_type;
  typedef typename third_array_type::value_type value_type;
  typedef typename DenseMatrix<value_type>::self array_type;
  
  array_type T(X.nrows()/2,X.ncols());
  
  typename SimdBlock<array_type>::self T1 = simd_block(T);
  typename ColSplit<typename SimdBlock<array_type>::self>::self T2 = cols(T1);
  for (int j=0, jmax=T1.ncols(); j<jmax; ++j)
  {
    typename ColSplit<typename SimdBlock<array_type>::self>::self::reference T3 = T2[j];
    sample2_conv<S>(cols(simd_block(X))[j],Y,T3);
  }
    
  typename RowSplit<third_array_type>::self Z1 = rows(Z);
  for (int i=0, imax=T.nrows(); i<imax; ++i)
  { 
    typename RowSplit<third_array_type>::self::reference Z2 = Z1[i];  
    sample2_conv<S>(rows(T)[i],Y,Z2);
  }
}


//}

#endif



//template<int N,class A1,class A2=A1,class V=PROMOTE2(typename A1::value_type,typename A2::value_type)>
//class sse_conv_vector_generator : public two_arrays_value_generator<const A1,const A2,V>
//{
//  public:
//	  typedef sse_conv_vector_generator self;
//    typedef two_arrays_value_generator<const A1,const A2,V> base;
//  public:
//    sse_conv_vector_generator(first_array_type &x, second_array_type &y) : base(x,y) {}
//    const_reference operator[](const index_type &i) const 
//    {
//      index_type i0=i-lower_bound();
//      const size_type n1=X.size(), n2=Y.size();
//      const value_type *it1 = &*X.begin();
//      const value_type *it2 = &*Y.begin(), *end2=&*Y.begin();
// 
//	    value_type result=0;
//
//      if (n1<=n2) { if (i0<n1) { it1+=i0; end2+=(i0+1); } else if (i0<n2) { it1+=i0;                 end2+=n1;     } else { it1+=n2-1; it2+=i0-n2+1; end2+=n1; } }
//             else { if (i0<n2) { it1+=i0; end2+=(i0+1); } else if (i0<n1) { it1+=n2-1; it2+=i0-n2+1; end2+=(i0+1); } else { it1+=n2-1; it2+=i0-n2+1; end2+=n1; } }   
//
//
//      //Problem if n2<4
//	    int rest=( ((int)(&*it2)-(16/N))%16 )/(16/N)+1;
//      
//      const value_type *end = it2+((16/N)-rest);
//
//      for (; it2!=end;++it2,--it1) 
//        result+=*it1 * *it2; 
//
//	    rest=( ((int)(&*end2))%16 )/(16/N);
//      end = end2-rest;
//      
//      if (it2!=end)
//      {
//	      TinyVector<N,value_type>::self a,b,t;
//
//        it1-= (N-1);
//        loadr(b,&*it2);
//		    loadu(a,&*it1);
//		    t=a*b;
//        for(it2+=4,it1-=N; it2!=end;it2+=N,it1-=N)
//	      {
//          loadr(b,&*it2);
//		      loadu(a,&*it1);
//		      t+=a*b;
//	      }
//	      it1+=(N-1);
//        result += sum(t);
//      }
//
//      for (;it2!=end2; --it1,++it2) 
//        result+= *it1 * *it2;
//
//	    return result;
//	  }
//
//    size_type lower_bound() const { return X.lower_bound()+Y.lower_bound(); }
//    size_type size() const { return X.size()+Y.size()-1; }
//    void resize(const size_type &d) { assert(d==size()); }
//};
//
//#ifdef SSE
//template<class A1,class A2>
//struct convVector<const Vector<dense_vector_generator<float,A1> >,const Vector<dense_vector_generator<float,A2> >,float>
//{
//  typedef const Vector<dense_vector_generator<float,A1> > first_array_type;
//  typedef const Vector<dense_vector_generator<float,A2> > second_array_type;
//  typedef float value_type; 
//  typedef sse_conv_vector_generator<4,first_array_type,second_array_type,value_type> generator_type; 
//  typedef typename GeneratorArray<generator_type>::self self;
//};
//template<int N1,int N2>
//struct convVector<const Vector<simd_vector_generator<N1,float> >,const Vector<simd_vector_generator<N2,float> >,float>
//{
//  typedef const Vector<simd_vector_generator<N1,float> > first_array_type;
//  typedef const Vector<simd_vector_generator<N2,float> > second_array_type;
//  typedef float value_type; 
//  typedef sse_conv_vector_generator<4,first_array_type,second_array_type,value_type> generator_type; 
//  typedef typename GeneratorArray<generator_type>::self self;
//};
//#endif //MMX
//
//#ifdef SSE2
//template<class A1,class A2>
//struct convVector<const Vector<dense_vector_generator<double,A1> >,const Vector<dense_vector_generator<double,A2> >,double>
//{
//  typedef const Vector<dense_vector_generator<double,A1> > first_array_type;
//  typedef const Vector<dense_vector_generator<double,A2> > second_array_type;
//  typedef double value_type; 
//  typedef sse_conv_vector_generator<2,first_array_type,second_array_type,value_type> generator_type; 
//  typedef typename GeneratorArray<generator_type>::self self;
//}; 
//template<int N1,int N2>
//struct convVector<const Vector<simd_vector_generator<N1,double> >,const Vector<simd_vector_generator<N2,double> >,double>
//{
//  typedef const Vector<simd_vector_generator<N1,double> > first_array_type;
//  typedef const Vector<simd_vector_generator<N2,double> > second_array_type;
//  typedef double value_type; 
//  typedef sse_conv_vector_generator<2,first_array_type,second_array_type,value_type> generator_type; 
//  typedef typename GeneratorArray<generator_type>::self self;
//};
//#endif //SSE

