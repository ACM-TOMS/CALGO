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

#ifndef VECTORITERATOR_H
#define VECTORITERATOR_H

#include "array/arrayiterator.h"

//namespace genial
//{

//template<class V> class dense_vector_generator;
template<class V,class A> class dense_vector_generator;
template<int N,class V> class tiny_vector_generator;

using namespace std;

template<class A>
class index_vector_iterator : public index_array_iterator<A>
{
  public:
    typedef index_vector_iterator self;
    typedef index_array_iterator<A> base;
    ARRAY_BASE_TYPES

  protected:
    inline explicit index_vector_iterator(array_type &x) : base(x,0) {}
    inline index_vector_iterator(array_type &x, const index_type &i) : base(x,i) {}
    
    inline index_vector_iterator(const self &x) : base(x) {}
    template<class T2> inline index_vector_iterator(const index_vector_iterator<T2> &x) : base(x) {}
};


template<class A>
class raster_vector_iterator : public index_vector_iterator<A>
{
  public:
    typedef raster_vector_iterator self;
    typedef index_vector_iterator<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;
    
    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  protected:
    using base::X;

  public:
    inline explicit raster_vector_iterator(array_type &x) : base(x) {}
    inline raster_vector_iterator(array_type &x, const index_type &n) : base(x,n) {}
    inline raster_vector_iterator(const self &x) : base(x) {}
    template<class T2> inline raster_vector_iterator(const raster_vector_iterator<T2> &x) : base(x) {}

    inline self &operator=(const self &r) { base::operator=(r); return *this; }
    
    using base::array;
    using base::pos;

    inline self &operator++() { ++pos(); return *this; }
    inline self &operator--() { --pos(); return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    inline self  operator+ (const difference_type &n) const { return self(X,pos()+n); }
    inline self  operator- (const difference_type &n) const { return self(X,pos()-n); }
    inline self &operator+=(const difference_type &n)       { pos()+=n; return *this; }
    inline self &operator-=(const difference_type &n)       { pos()-=n; return *this; }

    inline difference_type operator-(const self &x) const { return pos()-x.pos(); }
};

template<class G> struct raster_iterator_rebind<      Vector<G> > { typedef raster_vector_iterator<      Vector<G> > other; };
template<class G> struct raster_iterator_rebind<const Vector<G> > { typedef raster_vector_iterator<const Vector<G> > other; };


template<class A>
class raster_dense_vector_iterator : public array_traits<A>
{
  public:
    typedef raster_dense_vector_iterator self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;

    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  public:
    pointer data;

  public:
    inline explicit raster_dense_vector_iterator(array_type &x) : data(&x[0]) {}
    inline raster_dense_vector_iterator(array_type &x, const index_type &n) : data(&x[n]) {}
    inline raster_dense_vector_iterator(const pointer &p) : data(p) {}

    inline raster_dense_vector_iterator(const self &x) : data(x.data) {}
    
    inline self &operator=(const self &x) { data=x.data; return *this; }

    inline reference       operator*()       { return *data; }
    inline const_reference operator*() const { return *data; }

    inline self &operator++() { ++data; return *this; }
    inline self &operator--() { --data; return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    inline self  operator+ (const difference_type &n) const { return self(data+n); }
    inline self  operator- (const difference_type &n) const { return self(data-n); }
    inline self &operator+=(const difference_type &n)       { data+=n; return *this; }
    inline self &operator-=(const difference_type &n)       { data-=n; return *this; }

    inline difference_type operator-(const self &x) const { return &**this-&*x; }

    inline bool operator==(const self &r) const { return data==r.data; }
    inline bool operator!=(const self &r) const { return data!=r.data; }
    
    inline bool operator<(const self &x) { return data<x.data; }    
};

template<class A>
class raster_stride_dense_vector_iterator : public array_traits<A>
{
  public:
    typedef raster_stride_dense_vector_iterator self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;

    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  private:
    pointer data;
    size_type is;

  public:
    inline explicit raster_stride_dense_vector_iterator(array_type &x) : data(&x[0]), is(x.generator().stride()) {}
    inline raster_stride_dense_vector_iterator(array_type &x, const index_type &n) : data(&x[n]), is(x.generator().stride()) {}
    inline raster_stride_dense_vector_iterator(pointer p, const size_type &s) : data(p), is(s) {}

    inline raster_stride_dense_vector_iterator(const self &x) : data(x.data), is(x.is) {}
    
    inline self &operator=(const self &x) { data=x.data; is=x.is; return *this; }

    inline reference       operator*()       { return *data; }
    inline const_reference operator*() const { return *data; }

    inline self &operator++() { data+=stride(); return *this; }
    inline self &operator--() { data-=stride(); return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    inline self  operator+ (const difference_type &n) const { return self(data+n*stride(),stride()); }
    inline self  operator- (const difference_type &n) const { return self(data-n*stride(),stride()); }
    inline self &operator+=(const difference_type &n)       { data+=n*stride(); return *this; }
    inline self &operator-=(const difference_type &n)       { data-=n*stride(); return *this; }

    inline difference_type operator-(const self &x) const { return (&**this-&*x)/stride(); }

    inline bool operator==(const self &r) const { return data==r.data; }
    inline bool operator!=(const self &r) const { return data!=r.data; }
    
    inline size_type stride() const { return is; }
    
    inline bool operator<(const self &x) { return data<x.data; }        
};

template<class T,int N>
class raster_simd_block_dense_vector_iterator : public array_traits<T>
{
  public:
    typedef raster_simd_block_dense_vector_iterator self;
    typedef array_traits<T> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;

    template<class T2> struct rebind { typedef typename raster_iterator_rebind<T2>::other other; };
    
    typedef typename array_type::generator_type::template_reference template_reference;

  private:
    typename array_traits<reference>::pointer ptr;

  public:
    inline explicit raster_simd_block_dense_vector_iterator(array_type &x) : ptr(&x.generator().array()[0]) {}
    inline raster_simd_block_dense_vector_iterator(array_type &x, const index_type &n) : ptr(&x.generator().array()[n*N]) {}
    inline raster_simd_block_dense_vector_iterator(typename array_traits<reference>::pointer p) : ptr(p) {}

    inline raster_simd_block_dense_vector_iterator(const self &x) : ptr(x.ptr) {}
    
    inline self &operator=(const self &x) { ptr=x.ptr; return *this; }

    inline reference       operator*()       { return template_reference((typename template_reference::pointer)ptr); }
    inline const_reference operator*() const { return template_reference((typename template_reference::pointer)ptr); }

    inline self &operator++() { ptr+=N; return *this; }
    inline self &operator--() { ptr-=N; return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    inline self  operator+ (const difference_type &n) const { return self(ptr+n*N); }
    inline self  operator- (const difference_type &n) const { return self(ptr-n*N); }
    inline self &operator+=(const difference_type &n)       { ptr+=n*N; return *this; }
    inline self &operator-=(const difference_type &n)       { ptr-=n*N; return *this; }

    inline difference_type operator-(const self &x) const { return &**this-&*x; }

    inline bool operator==(const self &r) const { return ptr==r.ptr; }
    inline bool operator!=(const self &r) const { return ptr!=r.ptr; }
    
    inline bool operator<(const self &x) { return ptr<x.ptr; }            
};

#ifdef NDEBUG
template<class V        > struct raster_iterator_rebind<      Vector<data_vector_generator     <V  > > > { typedef raster_dense_vector_iterator<      Vector<data_vector_generator     <V  > > > other; };
template<class V        > struct raster_iterator_rebind<const Vector<data_vector_generator     <V  > > > { typedef raster_dense_vector_iterator<const Vector<data_vector_generator     <V  > > > other; };
template<class V,class A> struct raster_iterator_rebind<      Vector<dense_vector_generator    <V,A> > > { typedef raster_dense_vector_iterator<      Vector<dense_vector_generator    <V,A> > > other; };
template<class V,class A> struct raster_iterator_rebind<const Vector<dense_vector_generator    <V,A> > > { typedef raster_dense_vector_iterator<const Vector<dense_vector_generator    <V,A> > > other; };
template<class T,int C  > struct raster_iterator_rebind<      Vector<sub_dense_vector_generator<T,C> > > { typedef raster_dense_vector_iterator<      Vector<sub_dense_vector_generator<T,C> > > other; };
template<class T,int C  > struct raster_iterator_rebind<const Vector<sub_dense_vector_generator<T,C> > > { typedef raster_dense_vector_iterator<const Vector<sub_dense_vector_generator<T,C> > > other; };
template<class T,int C  > struct raster_iterator_rebind<      Vector<row_dense_matrix_generator<T,C> > > { typedef raster_dense_vector_iterator<      Vector<row_dense_matrix_generator<T,C> > > other; };
template<class T,int C  > struct raster_iterator_rebind<const Vector<row_dense_matrix_generator<T,C> > > { typedef raster_dense_vector_iterator<const Vector<row_dense_matrix_generator<T,C> > > other; };
template<int N  ,class V> struct raster_iterator_rebind<      Vector<tiny_vector_generator     <N,V> > > { typedef raster_dense_vector_iterator<      Vector<tiny_vector_generator     <N,V> > > other; };
template<int N  ,class V> struct raster_iterator_rebind<const Vector<tiny_vector_generator     <N,V> > > { typedef raster_dense_vector_iterator<const Vector<tiny_vector_generator     <N,V> > > other; };

template<class T,int C  > struct raster_iterator_rebind<      Vector<stride_dense_vector_generator<T,C> > > { typedef raster_stride_dense_vector_iterator<      Vector<stride_dense_vector_generator<T,C> > > other; };
template<class T,int C  > struct raster_iterator_rebind<const Vector<stride_dense_vector_generator<T,C> > > { typedef raster_stride_dense_vector_iterator<const Vector<stride_dense_vector_generator<T,C> > > other; };
template<class T,int C  > struct raster_iterator_rebind<      Vector<col_dense_matrix_generator   <T,C> > > { typedef raster_stride_dense_vector_iterator<      Vector<col_dense_matrix_generator   <T,C> > > other; };
template<class T,int C  > struct raster_iterator_rebind<const Vector<col_dense_matrix_generator   <T,C> > > { typedef raster_stride_dense_vector_iterator<const Vector<col_dense_matrix_generator   <T,C> > > other; };

//Pas certain que ca marche
template<class V> struct raster_iterator_rebind<      Vector<shift_array_generator<      Vector<dense_vector_generator<V> > > > > { typedef raster_dense_vector_iterator<      Vector<shift_array_generator<      Vector<dense_vector_generator<V> > > > > other; };
template<class V> struct raster_iterator_rebind<      Vector<shift_array_generator<const Vector<dense_vector_generator<V> > > > > { typedef raster_dense_vector_iterator<      Vector<shift_array_generator<const Vector<dense_vector_generator<V> > > > > other; };
template<class V> struct raster_iterator_rebind<const Vector<shift_array_generator<      Vector<dense_vector_generator<V> > > > > { typedef raster_dense_vector_iterator<const Vector<shift_array_generator<      Vector<dense_vector_generator<V> > > > > other; };
template<class V> struct raster_iterator_rebind<const Vector<shift_array_generator<const Vector<dense_vector_generator<V> > > > > { typedef raster_dense_vector_iterator<const Vector<shift_array_generator<const Vector<dense_vector_generator<V> > > > > other; };

//Pas certain que ca marche
template<int N,class V> struct raster_iterator_rebind<      Vector<shift_array_generator<      Vector<tiny_vector_generator<N,V> > > > > { typedef raster_dense_vector_iterator<      Vector<shift_array_generator<      Vector<tiny_vector_generator<N,V> > > > > other; };
template<int N,class V> struct raster_iterator_rebind<      Vector<shift_array_generator<const Vector<tiny_vector_generator<N,V> > > > > { typedef raster_dense_vector_iterator<      Vector<shift_array_generator<const Vector<tiny_vector_generator<N,V> > > > > other; };
template<int N,class V> struct raster_iterator_rebind<const Vector<shift_array_generator<      Vector<tiny_vector_generator<N,V> > > > > { typedef raster_dense_vector_iterator<const Vector<shift_array_generator<      Vector<tiny_vector_generator<N,V> > > > > other; };
template<int N,class V> struct raster_iterator_rebind<const Vector<shift_array_generator<const Vector<tiny_vector_generator<N,V> > > > > { typedef raster_dense_vector_iterator<const Vector<shift_array_generator<const Vector<tiny_vector_generator<N,V> > > > > other; };

template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef,int C> struct raster_iterator_rebind<      Vector<simd_block_dense_vector_generator<T,N,Cast,TempRef,C> > > { typedef raster_simd_block_dense_vector_iterator<      Vector<simd_block_dense_vector_generator<T,N,Cast,TempRef,C> >,N> other; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef,int C> struct raster_iterator_rebind<const Vector<simd_block_dense_vector_generator<T,N,Cast,TempRef,C> > > { typedef raster_simd_block_dense_vector_iterator<const Vector<simd_block_dense_vector_generator<T,N,Cast,TempRef,C> >,N> other; };
#endif



template<class A>
class shift_raster_vector_iterator : public index_vector_iterator<A>
{
  public:
    typedef shift_raster_vector_iterator self;
    typedef index_vector_iterator<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;
    
    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  protected:
    using base::X;

  public:
    inline explicit shift_raster_vector_iterator(array_type &x) : base(x,x.lower_bound()) {}
    inline shift_raster_vector_iterator(array_type &x, const index_type &n) : base(x,n) {}
    inline shift_raster_vector_iterator(const self &x) : base(x) {}
    template<class T2> inline shift_raster_vector_iterator(const shift_raster_vector_iterator<T2> &x) : base(x) {}

    inline self &operator=(const self &r) { base::operator=(r); return *this; }
    
    using base::array;
    using base::pos;

    inline self &operator++() { ++pos(); return *this; }
    inline self &operator--() { --pos(); return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    inline self  operator+ (const difference_type &n) const { return self(X,pos()+n); }
    inline self  operator- (const difference_type &n) const { return self(X,pos()-n); }
    inline self &operator+=(const difference_type &n)       { pos()+=n; return *this; }
    inline self &operator-=(const difference_type &n)       { pos()-=n; return *this; }

    inline difference_type operator-(const self &x) const { return pos()-x.pos(); }
};

template<class G> struct shift_raster_iterator_rebind<      Vector<G> > { typedef shift_raster_vector_iterator<      Vector<G> > other; };
template<class G> struct shift_raster_iterator_rebind<const Vector<G> > { typedef shift_raster_vector_iterator<const Vector<G> > other; };


#endif
