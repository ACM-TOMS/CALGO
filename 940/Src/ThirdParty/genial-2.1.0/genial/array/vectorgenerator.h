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

#ifndef VECTORGENERATOR_H
#define VECTORGENERATOR_H

#include <algorithm>
#include <sstream>

//namespace genial
//{

//Group=Vectors functions

using namespace std;


template<class G> class Vector;

template<class V>
struct value_vector_generator_traits
{
  typedef vector_array_tag array_category;

  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef value_type reference;
  typedef value_type const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;

  //typedef size_t int_type;
  //typedef ptrdiff_t difference_type;
  typedef int int_type;
  typedef int difference_type;
  typedef int_type index_type;
  typedef int_type size_type;

  template<class V2> struct array_rebind { typedef Vector<dense_vector_generator<V2> > other; };

  template<class A> struct iterator_rebind       { typedef typename raster_iterator_rebind<      A>::other other; };
  template<class A> struct const_iterator_rebind { typedef typename raster_iterator_rebind<const A>::other other; };
};

template<class V, class Ref=V &, class ConstRef=const V &>
struct reference_vector_generator_traits
{
  typedef vector_array_tag array_category;

  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;

  //typedef size_t int_type;
  //typedef ptrdiff_t difference_type;
  typedef int int_type;
  typedef int difference_type;
  typedef int_type index_type;
  typedef int_type size_type;

  template<class V2> struct array_rebind { typedef Vector<dense_vector_generator<V2> > other; };

  template<class A> struct iterator_rebind       { typedef typename raster_iterator_rebind<      A>::other other; };
  template<class A> struct const_iterator_rebind { typedef typename raster_iterator_rebind<const A>::other other; };
};


template<class V>
class value_vector_generator : public value_vector_generator_traits<V>
{
  public:
    typedef value_vector_generator self;
    typedef value_vector_generator_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

  public:
    inline value_vector_generator() {}

    inline index_type lower_bound() const { return 0; }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline size_type size() const { return 0; }
    inline void resize(const size_type &s) { assert(s==size()); }
};


template<class V,class Ref,class ConstRef>
class reference_vector_generator : public reference_vector_generator_traits<V,Ref,ConstRef>
{
  public:
    typedef reference_vector_generator self;
    typedef reference_vector_generator_traits<V,Ref,ConstRef> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

  public:
    inline reference_vector_generator() {}

    inline index_type lower_bound() const { return 0; }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline size_type size() const { return 0; }
    inline void resize(const size_type &s) { assert(s==size()); }
};

template<class V>
class sized_reference_vector_generator : public reference_vector_generator_traits<V>
{
  public:
    typedef sized_reference_vector_generator self;
    typedef reference_vector_generator_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

  protected:
    size_type si;

  public:
    inline sized_reference_vector_generator(const size_type &s) : si(s) {}

    inline index_type lower_bound() const { return 0; }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline size_type size() const { return si; }
    inline int_type  nelms() const { return si; }
    inline void resize(const size_type &s) { si=s; }

    inline void swap(self &x) { std::swap(si,x.si); }
};

template<class V>
class sized_value_vector_generator : public value_vector_generator_traits<V>
{
  public:
    typedef sized_value_vector_generator self;
    typedef value_vector_generator_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

  protected:
    size_type si;

  public:
    inline sized_value_vector_generator(const size_type &s) : si(s) {}

    inline index_type lower_bound() const { return 0; }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline size_type size () const { return si; }
    inline int_type  nelms() const { return si; }
    inline void resize(const size_type &s) { si=s; }

    inline void swap(self &x) { std::swap(si,x.si); }
};


template<class V>
class data_vector_generator : public sized_reference_vector_generator<V>
{
  public:
    typedef data_vector_generator self;
    typedef sized_reference_vector_generator<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

  //protected:
  public:
    value_type *data;

  public:
    inline data_vector_generator(int_type n, value_type *p) : base(n), data(p) { }
    template<class G> inline data_vector_generator(const Vector<G> &X) : base((X.size()*sizeof(typename Vector<G>::value_type))/sizeof(value_type)), data((pointer)&*X.begin()) {}

    using base::size;
    using base::lower_bound;
    using base::nelms;

    inline reference       operator[](const index_type &i)       { assert((i>=lower_bound()) && (i-lower_bound()<size())); return data[i]; }
    inline const_reference operator[](const index_type &i) const { assert((i>=lower_bound()) && (i-lower_bound()<size())); return data[i]; }

    inline void resize(const size_type &s) { assert(s==size()); }

    inline pointer       begin()       { return data; }
    inline const_pointer begin() const { return data; }
    inline pointer       end  ()       { return data+size(); }
    inline const_pointer end  () const { return data+size(); }

    inline void swap(self &x) { base::swap(x); std::swap(data,x.data); }
};

template<class V> inline void swap(data_vector_generator<V> &x, data_vector_generator<V> &y) { x.swap(y); }

template<class V>
struct DataVector : public reference_vector_generator_traits<V>
{
  typedef V value_type;
  typedef data_vector_generator<value_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class V> typename DataVector<      V>::self data_vector(int n,       V *p) { return typename DataVector<      V>::self(n,p); }
//template<class V> typename DataVector<const V>::self data_vector(int n, const V *p) { return typename DataVector<const V>::self(n,p); }



template<class V,class A>
class dense_vector_generator : public data_vector_generator<V>
{
  public:
    typedef dense_vector_generator self;
    typedef data_vector_generator<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

    typedef A allocator_type;

  protected:
    allocator_type alloc;

  public:
    using base::size;
    using base::lower_bound;
    using base::nelms;
    using base::begin;
    using base::end;

    inline dense_vector_generator() : base(0,NULL) {}
    inline explicit dense_vector_generator(int_type n            ) : base(n,alloc.allocate(n)) { uninitialized_fill(begin(),end()); }
    inline dense_vector_generator(int_type n, const value_type &v) : base(n,alloc.allocate(n)) { uninitialized_fill(begin(),end(),v           ); }
    dense_vector_generator(int_type n, const char       *s) : base(n,alloc.allocate(n)) { istringstream iss(s        ); for (pointer i=begin();i!=end();++i) iss>>*i; }
    dense_vector_generator(int_type n, const string     &s) : base(n,alloc.allocate(n)) { istringstream iss(s.c_str()); for (pointer i=begin();i!=end();++i) iss>>*i; }

    inline dense_vector_generator                           (const self       &X) : base(X.size(),alloc.allocate(X.nelms())) { copyn(X.size(),X.begin(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit dense_vector_generator(const Vector<G>  &X) : base(X.size(),alloc.allocate(X.nelms())) { copyn(X.size(),X.begin(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit dense_vector_generator(const Array<1,G> &X) : base(X.size(),alloc.allocate(X.nelms())) { copyn(X.size(),X.begin(),raw_storage_iterator<pointer,value_type>(begin())); }

    template<class V2,class A2> inline explicit dense_vector_generator(const list<V2,A2> &X) : base(X.size(),alloc.allocate(X.nelms())) { copy(X.begin(),X.end(),raw_storage_iterator<pointer,value_type>(begin())); }

    inline ~dense_vector_generator() { alloc.deallocate(this->data,nelms()); }

    inline void swap(self &x) { base::swap(x); std::swap(alloc,x.alloc); }

    inline reference       operator[](const index_type &i)       { return base::operator[](i); }
    inline const_reference operator[](const index_type &i) const { return base::operator[](i); }

    inline void inv(const value_type &y, const int_type &n) { assert(n<size()); (*this)[n]=y; }

    inline void resize(const int_type &n) { if (size()==n) return; alloc.deallocate(this->data,nelms()); this->si=n; this->data=alloc.allocate(nelms()); uninitialized_fill(begin(),end()); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }

    inline bool operator==(const self &x) const { return size()==x.size() &&  equal(begin(), end(), x.begin()); }
    inline bool operator!=(const self &x) const { return size()!=x.size() || !equal(begin(), end(), x.begin()); }
};

template<class V,class A> inline void swap(dense_vector_generator<V,A> &x, dense_vector_generator<V,A> &y) { x.swap(y); }

//{unsecret}
//{group:Vectors Interfaces}
//summary: Dense vector type
//Arguments:
//  V - Type of the elements
//  A - Optional allocator
//example:
//  DenseVector<int>::self X(4,"1 2 3 4");
//  typedef DenseVector<complex<float> >::self MyVector;
template<class V,class A>
struct DenseVector : public reference_vector_generator_traits<V>
{
  typedef V value_type;
  typedef A allocator_type;
  typedef dense_vector_generator<value_type,allocator_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};


template<int N, class V>
class tiny_vector_generator : public reference_vector_generator_traits<V>
{
  public:
    typedef tiny_vector_generator self;
    typedef reference_vector_generator_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;

    template<class V2> struct array_rebind { typedef Vector<tiny_vector_generator<N,V2> > other; };

  protected:
    __aligned value_type data[N];

  public:
    inline tiny_vector_generator() { uninitialized_fill(begin(),end()); }
    inline explicit tiny_vector_generator(const value_type &v) { uninitialized_fill(begin(),end(),v); }
    inline tiny_vector_generator(const value_type &v0,const value_type &v1) { (*this)[0]=v0; (*this)[1]=v1; }
    inline tiny_vector_generator(const value_type &v0,const value_type &v1,const value_type &v2) { (*this)[0]=v0; (*this)[1]=v1; (*this)[2]=v2; }
    inline tiny_vector_generator(const value_type &v0,const value_type &v1,const value_type &v2,const value_type &v3) { (*this)[0]=v0; (*this)[1]=v1; (*this)[2]=v2; (*this)[3]=v3; }
    tiny_vector_generator(const char   *s) { istringstream iss(s        ); for (pointer i=begin();i!=end();++i) iss>>*i; }
    tiny_vector_generator(const string &s) { istringstream iss(s.c_str()); for (pointer i=begin();i!=end();++i) iss>>*i; }

    template<class G> inline tiny_vector_generator(const Vector<G> &x) { assert(size()==x.size()); copy(x.begin(),x.end(),begin()); }

    inline tiny_vector_generator(const self &x) { assert(size()==x.size()); copy(x.begin(),x.end(),begin()); }

    inline void swap(self &x) { swap_ranges_n(size(),begin(),x.begin()); }

    inline reference       operator[](const index_type &n)       { assert(n<size()); return data[n]; }
    inline const_reference operator[](const index_type &n) const { assert(n<size()); return data[n]; }

    inline size_type  size       () const { return N; }
    inline index_type lower_bound() const { return 0; }
    inline int_type   nelms      () const { return N; }

    inline pointer       begin()       { return data; }
    inline const_pointer begin() const { return data; }
    inline pointer       end  ()       { return begin()+nelms(); }
    inline const_pointer end  () const { return begin()+nelms(); }

    inline void resize(const size_type &s) { assert(size()==s); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }

    inline bool operator==(const self &r) const { return size()==r.size() &&  equal(begin(), end(), r.begin()); }
    inline bool operator!=(const self &r) const { return !(*this==r); }
};

template<int N,class V> inline void swap(tiny_vector_generator<N,V> &x, tiny_vector_generator<N,V> &y) { x.swap(y); }



template<int N,class V,class G2> inline void affect(Vector<tiny_vector_generator<N,V> > &X, const Vector<G2>  &Y) { for (int i=0; i<N; ++i) X[i]=Y[i]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 0,V> > &X, const Vector<G2>  &Y) {  }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 1,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 2,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 3,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 4,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 5,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 6,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 7,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 8,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator< 9,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<10,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<11,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<12,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; X[11]=Y[11]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<13,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; X[11]=Y[11]; X[12]=Y[12]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<14,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; X[11]=Y[11]; X[12]=Y[12]; X[13]=Y[13]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<15,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; X[11]=Y[11]; X[12]=Y[12]; X[13]=Y[13]; X[14]=Y[14]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<16,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; X[11]=Y[11]; X[12]=Y[12]; X[13]=Y[13]; X[14]=Y[14]; X[15]=Y[15]; }
template<class V,class G2> inline void affect(Vector<tiny_vector_generator<32,V> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; X[4]=Y[4]; X[5]=Y[5]; X[6]=Y[6]; X[7]=Y[7]; X[8]=Y[8]; X[9]=Y[9]; X[10]=Y[10]; X[11]=Y[11]; X[12]=Y[12]; X[13]=Y[13]; X[14]=Y[14]; X[15]=Y[15]; X[16]=Y[16]; X[17]=Y[17]; X[18]=Y[18]; X[19]=Y[19]; X[20]=Y[20]; X[21]=Y[21]; X[22]=Y[22]; X[23]=Y[23]; X[24]=Y[24]; X[25]=Y[25]; X[26]=Y[26]; X[27]=Y[27]; X[28]=Y[28]; X[29]=Y[29]; X[30]=Y[30]; X[31]=Y[31]; }


template<int V,class T>
struct simd_data
{
};

#ifdef MMX
template<> struct simd_data< 8,char            > { typedef m64c  type; };
template<> struct simd_data< 8,unsigned char   > { typedef m64b  type; };
template<> struct simd_data< 4,short           > { typedef m64s  type; };
template<> struct simd_data< 4,unsigned short  > { typedef m64w  type; };
template<> struct simd_data< 2,int             > { typedef m64i  type; };
#endif // MMX

#ifdef SSE
template<> struct simd_data< 2,float           > { typedef m64f   type; };
template<> struct simd_data< 4,float           > { typedef m128f  type; };
template<> struct simd_data< 1,complex<float > > { typedef m64cf  type; };
template<> struct simd_data< 2,complex<float > > { typedef m128cf type; };
template<> struct simd_data< 4,complex<float > > { typedef m256cf type; };
#endif //SSE

#ifdef SSE2
template<> struct simd_data< 1,double          > { typedef m64d   type; };
template<> struct simd_data< 2,double          > { typedef m128d  type; };
template<> struct simd_data< 4,double          > { typedef m256d  type; };
template<> struct simd_data< 1,complex<double> > { typedef m128cd type; };
template<> struct simd_data< 2,complex<double> > { typedef m256cd type; };
template<> struct simd_data<16,char            > { typedef m128c  type; };
template<> struct simd_data<16,unsigned char   > { typedef m128b  type; };
template<> struct simd_data< 8,short           > { typedef m128s  type; };
template<> struct simd_data< 8,unsigned short  > { typedef m128w  type; };
template<> struct simd_data< 4,int             > { typedef m128i  type; };
#endif //SSE2

template<int N,class V>
class simd_vector_generator : public reference_vector_generator_traits<V>
{
  public:
    typedef simd_vector_generator self;
    typedef reference_vector_generator_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;
    typedef typename base::index_type      index_type;
    typedef typename base::size_type       size_type;
    typedef typename base::difference_type difference_type;
    typedef typename base::int_type        int_type;
    
    typedef typename simd_data<N,V>::type data_type;

    template<class V2> struct array_rebind { typedef typename SimdVector<N,V2>::self other; };    
    //template<class V2> struct array_rebind { typedef Vector<self> other; };    

  public:
    data_type data;

  public:
    inline simd_vector_generator() { }
    inline explicit simd_vector_generator(const data_type &d) : data(d) {}
    inline explicit simd_vector_generator(const value_type &v) : data(v) { }
    inline simd_vector_generator(const value_type &v0,const value_type &v1) : data(v0,v1) { }
    inline simd_vector_generator(const value_type &v0,const value_type &v1,const value_type &v2,const value_type &v3) : data(v0,v1,v2,v3) {}
    inline simd_vector_generator(const value_type &v0,const value_type &v1,const value_type &v2,const value_type &v3, const value_type &v4,const value_type &v5,const value_type &v6,const value_type &v7) : data(v0,v1,v2,v3,v4,v5,v6,v7) {}
    explicit simd_vector_generator(const char   *s) { istringstream iss(s        ); for (pointer i=begin();i!=end();++i) iss>>*i; }
    explicit simd_vector_generator(const string &s) { istringstream iss(s.c_str()); for (pointer i=begin();i!=end();++i) iss>>*i; }

    template<class G> inline simd_vector_generator(const Vector<G> &x) { assert(size()==x.size()); copy(x.begin(),x.end(),begin()); }

    inline simd_vector_generator(const self &x) { assert(size()==x.size()); data=x.data; }

    inline void swap(self &x) { data_type t=data; data=x.data; x.data=t; }

    inline size_type  size       () const { return N; }
    inline index_type lower_bound() const { return 0; }
    inline int_type   nelms      () const { return N; }

   	inline reference       operator[](const index_type &n)       { assert(n<size()); return data[n]; };
   	inline const_reference operator[](const index_type &n) const { assert(n<size()); return data[n]; };

    inline pointer       begin()       { return (value_type *)&data; }
    inline const_pointer begin() const { return (value_type *)&data; }
    inline pointer       end  ()       { return begin()+nelms(); }
    inline const_pointer end  () const { return begin()+nelms(); }

    inline void resize(const size_type &s) { assert(size()==s); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }

    inline bool operator==(const self &r) const { return size()==r.size() &&  equal(begin(), end(), r.begin()); }
    inline bool operator!=(const self &r) const { return !(*this==r); }
};

//{unsecret}
//{group:Vectors Interfaces}
//Summary: Tiny dense vector type
//Arguments:
//  N - Size
//  V - Type of the elements
//Example:
//  TinyVector<4,int>::self X("1 2 3 4");
//  typedef TinyDenseVector<4,complex<float> >::self MyVector;
template<int N, class V>
struct TinyVector : public reference_vector_generator_traits<V>
{
  typedef V value_type;
  typedef tiny_vector_generator<N,value_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int N,class V> struct SimdVector
{ 
  typedef V value_type;
  typedef simd_vector_generator<N,value_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

#ifdef MMX
template<> struct TinyVector    < 8,char            > { enum {N= 8}; typedef char            value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 8,unsigned char   > { enum {N= 8}; typedef unsigned char   value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 4,short           > { enum {N= 4}; typedef short           value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 4,unsigned short  > { enum {N= 4}; typedef unsigned short  value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 2,int             > { enum {N= 2}; typedef int             value_type; typedef SimdVector<N,value_type>::self self; };
#endif //MMX

#ifdef SSE
template<> struct TinyVector    < 2,float           > { enum {N= 2}; typedef float           value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 4,float           > { enum {N= 4}; typedef float           value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 1,complex<float > > { enum {N= 1}; typedef complex<float > value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 2,complex<float > > { enum {N= 2}; typedef complex<float > value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 4,complex<float > > { enum {N= 4}; typedef complex<float > value_type; typedef SimdVector<N,value_type>::self self; };
#endif //SSE

#ifdef SSE2
template<> struct TinyVector    < 2,double          > { enum {N= 2}; typedef double          value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 4,double          > { enum {N= 4}; typedef double          value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 1,complex<double> > { enum {N= 1}; typedef complex<double> value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 2,complex<double> > { enum {N= 2}; typedef complex<double> value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    <16,char            > { enum {N=16}; typedef char            value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    <16,unsigned char   > { enum {N=16}; typedef unsigned char   value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 8,short           > { enum {N= 8}; typedef short           value_type; typedef SimdVector<N,value_type>::self self; };
template<> struct TinyVector    < 4,int             > { enum {N= 4}; typedef int             value_type; typedef SimdVector<N,value_type>::self self; };
#endif //SSE2



template<class E,class Tr,int N,class V> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const Vector<simd_vector_generator<N,V> > &X) { return os << X.generator().data; }

template<int N,class V> inline void affect(Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { X.generator().data = Y.generator().data; }
template<int N,class V,class T> inline void affect(Vector<simd_vector_generator<N,V> > &X, const Vector<dense_vector_generator<T>  > &Y) { load(X.generator().data,&Y[0]); }
template<int N,class V,class T> inline void affect(Vector<simd_vector_generator<N,V> > &X, const Vector<tiny_vector_generator<N,T> > &Y) { load(X.generator().data,&Y[0]); }
template<int N,class V,class T,int C> inline void affect(Vector<simd_vector_generator<N,V> > &X, const Vector<tiny_sub_dense_vector_generator<N,T,C> > &Y) { load(X.generator().data,&Y[0]); }
template<class T,int C,int N,class V> inline void affect(Vector<tiny_sub_dense_vector_generator<2,T,C> > &X, const Vector<simd_vector_generator<N,V> > &Y) { store(&X[0],Y.generator().data); }

#if !defined(_MSC_VER)
template<int N,class V> inline void affect(Vector<simd_vector_generator<N,V> > &X, const V &v) { X.generator().data=v; }
#else
template<int N,class V> inline void affect(Vector<simd_vector_generator<N,V> > &X, const float           &v) { X.generator().data=v; }
template<int N,class V> inline void affect(Vector<simd_vector_generator<N,V> > &X, const complex<float > &v) { X.generator().data=v; }
template<int N,class V> inline void affect(Vector<simd_vector_generator<N,V> > &X, const double          &v) { X.generator().data=v; }
template<int N,class V> inline void affect(Vector<simd_vector_generator<N,V> > &X, const complex<double> &v) { X.generator().data=v; }
#endif

template<int N,class V> inline void load   (Vector<simd_vector_generator<N,V> > &X, const V *p) { load   (X.generator().data, p); }
template<int N,class V> inline void loadu  (Vector<simd_vector_generator<N,V> > &X, const V *p) { loadu  (X.generator().data, p); }
template<int N,class V> inline void loadr  (Vector<simd_vector_generator<N,V> > &X, const V *p) { loadr  (X.generator().data, p); }
template<int N,class V> inline void loadur (Vector<simd_vector_generator<N,V> > &X, const V *p) { loadur (X.generator().data, p); }
template<int N,class V> inline void loadl  (Vector<simd_vector_generator<N,V> > &X, const V *p) { loadl  (X.generator().data, p); }
template<int N,class V> inline void loadh  (Vector<simd_vector_generator<N,V> > &X, const V *p) { loadh  (X.generator().data, p); }
template<int N,class V> inline void loadlh (Vector<simd_vector_generator<N,V> > &X, const V *p) { loadlh (X.generator().data, p); }
template<int N,class V> inline void loadlh (Vector<simd_vector_generator<N,V> > &X, const V *p1, const V *p2) { loadlh (X.generator().data, p1,p2); }
template<int N,class V> inline void store  (V *p, const Vector<simd_vector_generator<N,V> > &X) { store  (p, X.generator().data); }
template<int N,class V> inline void storeu (V *p, const Vector<simd_vector_generator<N,V> > &X) { storeu (p, X.generator().data); }
template<int N,class V> inline void storer (V *p, const Vector<simd_vector_generator<N,V> > &X) { storer (p, X.generator().data); }
template<int N,class V> inline void storeur(V *p, const Vector<simd_vector_generator<N,V> > &X) { storeur(p, X.generator().data); }
template<int N,class V> inline void storel (V *p, const Vector<simd_vector_generator<N,V> > &X) { storel (p, X.generator().data); }
template<int N,class V> inline void storeh (V *p, const Vector<simd_vector_generator<N,V> > &X) { storeh (p, X.generator().data); }
template<int N,class V> inline void stream (V *p, const Vector<simd_vector_generator<N,V> > &X) { stream (p, X.generator().data); }
template<int N,class V> inline void storelh(V *p1, V *p2, const Vector<simd_vector_generator<N,V> > &X) { storelh(p1,p2,X.generator().data); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator- (                                              const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(                  -Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator+ (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data+Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator- (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data-Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data*Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data/Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > mullo     (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(mullo(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > mulhi     (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(mulhi(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > inv       (                                              const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(inv(Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator& (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data&Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator| (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data|Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator^ (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data^Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > and_not   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(and_not(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator+=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X+Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator-=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X-Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator*=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X*Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator/=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X/Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator*=(      Vector<simd_vector_generator<N,V> > &X, const V                                   &Y) { return X=X*Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator/=(      Vector<simd_vector_generator<N,V> > &X, const V                                   &Y) { return X=X/Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator&=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X&Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator|=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X|Y; }
template<int N,class V> inline Vector<simd_vector_generator<N,V> >&operator^=(      Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return X=X^Y; }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const float                               &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 *Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const complex<float>                      &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 *Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const double                              &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 *Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const complex<double>                     &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 *Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const float                               &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 /Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const complex<float>                      &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 /Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const double                              &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 /Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const complex<double>                     &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 /Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const float                               &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data*v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const complex<float>                      &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data*v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const double                              &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data*v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const complex<double>                     &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data*v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const float                               &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data/v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const complex<float>                      &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data/v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const double                              &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data/v                 ); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const complex<double>                     &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data/v                 ); }

//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const V                                   &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 *Y.generator().data); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const V                                   &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 /Y.generator().data); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const V                                   &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data*v                 ); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const V                                   &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data/v                 ); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const typename V::value_type              &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 *Y.generator().data); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const typename V::value_type              &v, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(v                 /Y.generator().data); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator* (const Vector<simd_vector_generator<N,V> > &X, const typename V::value_type              &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data*v                 ); }
//template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator/ (const Vector<simd_vector_generator<N,V> > &X, const typename V::value_type              &v) { return Vector<simd_vector_generator<N,V> >(X.generator().data/v                 ); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator==(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data==Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator!=(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data!=Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator> (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data> Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator>=(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data>=Y.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator< (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(Y.generator().data> X.generator().data); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > operator<=(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(X.generator().data<=Y.generator().data); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > min       (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(min(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > max       (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(max(X.generator().data,Y.generator().data)); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > sqrt (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(sqrt (X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > rcp  (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(rcp  (X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > rsqrt(const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(rsqrt(X.generator().data)); }

template<int n,int N,class V> inline Vector<simd_vector_generator<N,V> > shiftl (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(shiftl <n>(X.generator().data)); }
template<int n,int N,class V> inline Vector<simd_vector_generator<N,V> > shiftr (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(shiftr <n>(X.generator().data)); }
template<int n,int N,class V> inline Vector<simd_vector_generator<N,V> > shiftll(const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(shiftll<n>(X.generator().data)); }
template<int n,int N,class V> inline Vector<simd_vector_generator<N,V> > shiftrl(const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(shiftrl<n>(X.generator().data)); }
template<int n,int N,class V> inline Vector<simd_vector_generator<N,V> > shiftla(const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(shiftla<n>(X.generator().data)); }
template<int n,int N,class V> inline Vector<simd_vector_generator<N,V> > shiftra(const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(shiftra<n>(X.generator().data)); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > conj     (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(conj  (X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > imul     (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(imul  (X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > mimul    (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(mimul (X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > i1mul    (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(i1mul (X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > mi1mul   (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(mi1mul(X.generator().data)); }

template<int N,class V> inline Vector<simd_vector_generator<N,complex<V> > > conj     (const V &v, const Vector<simd_vector_generator<N,complex<V> > > &X) { return Vector<simd_vector_generator<N,complex<V> > >(conj(v,X.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,complex<V> > > imul     (const V &v, const Vector<simd_vector_generator<N,complex<V> > > &X) { return Vector<simd_vector_generator<N,complex<V> > >(imul(v,X.generator().data)); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > cmul     (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(cmul    (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > icmul    (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(icmul   (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > mimul    (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(mimul   (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > imul     (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(imul    (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > add_imul (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(add_imul(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > sub_imul (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(sub_imul(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > addsub   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(addsub  (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > subadd   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(subadd  (X.generator().data,Y.generator().data)); }

template<int i                      ,int N,class V> inline Vector<simd_vector_generator<N,V> > shuffle(const Vector<simd_vector_generator<N,V> > &X                                                    ) { return Vector<simd_vector_generator<N,V> >(shuffle<i          >(X.generator().data)); }
template<int i0,int i1              ,int N,class V> inline Vector<simd_vector_generator<N,V> > shuffle(const Vector<simd_vector_generator<N,V> > &X                                                    ) { return Vector<simd_vector_generator<N,V> >(shuffle<i0,i1      >(X.generator().data)); }
template<int i0,int i1,int i2,int i3,int N,class V> inline Vector<simd_vector_generator<N,V> > shuffle(const Vector<simd_vector_generator<N,V> > &X                                                    ) { return Vector<simd_vector_generator<N,V> >(shuffle<i0,i1,i2,i3>(X.generator().data)); }
template<int i                      ,int N,class V> inline Vector<simd_vector_generator<N,V> > shuffle(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(shuffle<i          >(X.generator().data,Y.generator().data)); }
template<int i0,int i1              ,int N,class V> inline Vector<simd_vector_generator<N,V> > shuffle(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(shuffle<i0,i1      >(X.generator().data,Y.generator().data)); }
template<int i0,int i1,int i2,int i3,int N,class V> inline Vector<simd_vector_generator<N,V> > shuffle(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(shuffle<i0,i1,i2,i3>(X.generator().data,Y.generator().data)); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > flip   (const Vector<simd_vector_generator<N,V> > &X) { return Vector<simd_vector_generator<N,V> >(flip (X.generator().data)); }

template<int N,class V> inline typename SimdVector<  N,typename complex_traits<V>::value_type>::self norm(const Vector<simd_vector_generator<N,V> > &X                                             ) { return typename SimdVector<  N,typename complex_traits<V>::value_type>::self(norm(X.generator().data)); }
template<int N,class V> inline typename SimdVector<  N,typename complex_traits<V>::value_type>::self abs (const Vector<simd_vector_generator<N,V> > &X                                             ) { return typename SimdVector<  N,typename complex_traits<V>::value_type>::self(abs (X.generator().data)); }
template<int N,class V> inline typename SimdVector<  N,typename complex_traits<V>::value_type>::self real(const Vector<simd_vector_generator<N,V> > &X                                             ) { return typename SimdVector<  N,typename complex_traits<V>::value_type>::self(real(X.generator().data)); }
template<int N,class V> inline typename SimdVector<  N,typename complex_traits<V>::value_type>::self imag(const Vector<simd_vector_generator<N,V> > &X                                             ) { return typename SimdVector<  N,typename complex_traits<V>::value_type>::self(imag(X.generator().data)); }
template<int N,class V> inline typename SimdVector<2*N,typename complex_traits<V>::value_type>::self norm(const Vector<simd_vector_generator<N,V> > &X,const Vector<simd_vector_generator<N,V> > &Y) { return typename SimdVector<2*N,typename complex_traits<V>::value_type>::self(norm(X.generator().data,Y.generator().data)); }
template<int N,class V> inline typename SimdVector<2*N,typename complex_traits<V>::value_type>::self abs (const Vector<simd_vector_generator<N,V> > &X,const Vector<simd_vector_generator<N,V> > &Y) { return typename SimdVector<2*N,typename complex_traits<V>::value_type>::self(abs (X.generator().data,Y.generator().data)); }
template<int N,class V> inline typename SimdVector<2*N,typename complex_traits<V>::value_type>::self real(const Vector<simd_vector_generator<N,V> > &X,const Vector<simd_vector_generator<N,V> > &Y) { return typename SimdVector<2*N,typename complex_traits<V>::value_type>::self(real(X.generator().data,Y.generator().data)); }
template<int N,class V> inline typename SimdVector<2*N,typename complex_traits<V>::value_type>::self imag(const Vector<simd_vector_generator<N,V> > &X,const Vector<simd_vector_generator<N,V> > &Y) { return typename SimdVector<2*N,typename complex_traits<V>::value_type>::self(imag(X.generator().data,Y.generator().data)); }

template<int N,class V> inline typename TinyVector<N/2,complex<V> >::self complexlo(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return typename TinyVector<N/2,complex<V> >::self(complexlo(X.generator().data,Y.generator().data)); }
template<int N,class V> inline typename TinyVector<N/2,complex<V> >::self complexhi(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return typename TinyVector<N/2,complex<V> >::self(complexhi(X.generator().data,Y.generator().data)); }

template<int N,class V> inline Vector<simd_vector_generator<N,V> > unpacklo (const Vector<simd_vector_generator<N,V> > &X                                                    ) { return Vector<simd_vector_generator<N,V> >(unpacklo(X.generator().data                   )); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > unpackhi (const Vector<simd_vector_generator<N,V> > &X                                                    ) { return Vector<simd_vector_generator<N,V> >(unpackhi(X.generator().data                   )); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > unpacklo (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(unpacklo(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > unpackhi (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(unpackhi(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > movehl   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(movehl  (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > movelh   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(movelh  (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > packlo   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(packlo  (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > packhi   (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(packhi  (X.generator().data,Y.generator().data)); }



template<int N,class V> inline Vector<simd_vector_generator<N,V> > hadd (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(hadd (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > hsub (const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(hsub (X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > h2add(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(h2add(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > h2sub(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(h2sub(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > h3add(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(h3add(X.generator().data,Y.generator().data)); }
template<int N,class V> inline Vector<simd_vector_generator<N,V> > h3sub(const Vector<simd_vector_generator<N,V> > &X, const Vector<simd_vector_generator<N,V> > &Y) { return Vector<simd_vector_generator<N,V> >(h3sub(X.generator().data,Y.generator().data)); }

template<int N,class V> inline void trn    (Vector<simd_vector_generator<N,V> > &X0, Vector<simd_vector_generator<N,V> > &X1) { trn    (X0.generator().data,X1.generator().data); }
template<int N,class V> inline void up_trn (Vector<simd_vector_generator<N,V> > &X0, Vector<simd_vector_generator<N,V> > &X1) { up_trn (X0.generator().data,X1.generator().data); }
template<int N,class V> inline void low_trn(Vector<simd_vector_generator<N,V> > &X0, Vector<simd_vector_generator<N,V> > &X1) { low_trn(X0.generator().data,X1.generator().data); }
template<int N,class V> inline void trn    (Vector<simd_vector_generator<N,V> > &X0, Vector<simd_vector_generator<N,V> > &X1, Vector<simd_vector_generator<N,V> > &X2,Vector<simd_vector_generator<N,V> > &X3) { trn    (X0.generator().data,X1.generator().data,X2.generator().data,X3.generator().data); }
template<int N,class V> inline void up_trn (Vector<simd_vector_generator<N,V> > &X0, Vector<simd_vector_generator<N,V> > &X1, Vector<simd_vector_generator<N,V> > &X2,Vector<simd_vector_generator<N,V> > &X3) { up_trn (X0.generator().data,X1.generator().data,X2.generator().data,X3.generator().data); }
template<int N,class V> inline void low_trn(Vector<simd_vector_generator<N,V> > &X0, Vector<simd_vector_generator<N,V> > &X1, Vector<simd_vector_generator<N,V> > &X2,Vector<simd_vector_generator<N,V> > &X3) { low_trn(X0.generator().data,X1.generator().data,X2.generator().data,X3.generator().data); }

template<int N,class V> inline V min  (const Vector<simd_vector_generator<N,V> > &X) { return get0(min  (X.generator().data)); }
template<int N,class V> inline V max  (const Vector<simd_vector_generator<N,V> > &X) { return get0(max  (X.generator().data)); }
template<int N,class V> inline V sum  (const Vector<simd_vector_generator<N,V> > &X) { return get0(sum  (X.generator().data)); }
template<int N,class V> inline V sumlo(const Vector<simd_vector_generator<N,V> > &X) { return get0(sumlo(X.generator().data)); }
template<int N,class V> inline V sumhi(const Vector<simd_vector_generator<N,V> > &X) { return get0(sumhi(X.generator().data)); }
template<int N,class V> inline typename TinyVector<2,V>::self sum  (const Vector<simd_vector_generator<N,V> > &X0, const Vector<simd_vector_generator<N,V> > &X1) { return typename TinyVector<2,V>::self((typename TinyVector<2,V>::self::generator_type::data_type)sum(X0.generator().data,X1.generator().data)); }
template<int N,class V> inline typename TinyVector<4,V>::self sum  (const Vector<simd_vector_generator<N,V> > &X0, const Vector<simd_vector_generator<N,V> > &X1, const Vector<simd_vector_generator<N,V> > &X2, const Vector<simd_vector_generator<N,V> > &X3) { return typename TinyVector<4,V>::self(sum(X0.generator().data,X1.generator().data,X2.generator().data,X3.generator().data)); }

template<      int N,class V> inline V get_sad(const Vector<simd_vector_generator<N,V> > &X) { return get_sad   (X.generator().data); }
template<int n,int N,class V> inline V get_sad(const Vector<simd_vector_generator<N,V> > &X) { return get_sad<n>(X.generator().data); }


//template<> struct norm_function<Vector<simd_vector_generator<2,complex<float> > > > : public unary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,float> > > { template<class V> typename norm_function<V>::const_reference operator()(const V &x) const { return norm(x); } };
//template<> struct abs_function <Vector<simd_vector_generator<2,complex<float> > > > : public unary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,float> > > { template<class V> typename abs_function <V>::const_reference operator()(const V &x) const { return abs (x); } };
//template<> struct arg_function <Vector<simd_vector_generator<2,complex<float> > > > : public unary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,float> > > { template<class V> typename arg_function <V>::const_reference operator()(const V &x) const { return arg (x); } };
//template<> struct real_function<Vector<simd_vector_generator<2,complex<float> > > > : public unary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,float> > > { template<class V> typename real_function<V>::const_reference operator()(const V &x) const { return real(x); } };
//template<> struct imag_function<Vector<simd_vector_generator<2,complex<float> > > > : public unary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,float> > > { template<class V> typename imag_function<V>::const_reference operator()(const V &x) const { return imag(x); } };
//template<> struct norm_binary_function<Vector<simd_vector_generator<2,complex<float> > > > : public binary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<4,float> > > { template<class V1,class V2> typename norm_binary_function<V1,V2>::const_reference operator()(const V1 &x,const V2 &y) const { return norm(x,y); } };
//template<> struct abs_binary_function <Vector<simd_vector_generator<2,complex<float> > > > : public binary_value_function<Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<2,complex<float> > >,Vector<simd_vector_generator<4,float> > > { template<class V1,class V2> typename abs_binary_function <V1,V2>::const_reference operator()(const V1 &x,const V2 &y) const { return abs (x,y); } };


template<int N,class V> struct norm_function<Vector<simd_vector_generator<N,V> > > : public unary_value_function<Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<N,typename complex_traits<V>::value_type> > > { template<class T> typename norm_function<V>::const_reference operator()(const T &x) const { return norm(x); } };
template<int N,class V> struct abs_function <Vector<simd_vector_generator<N,V> > > : public unary_value_function<Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<N,typename complex_traits<V>::value_type> > > { template<class T> typename abs_function <V>::const_reference operator()(const T &x) const { return abs (x); } };
template<int N,class V> struct real_function<Vector<simd_vector_generator<N,V> > > : public unary_value_function<Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<N,typename complex_traits<V>::value_type> > > { template<class T> typename real_function<V>::const_reference operator()(const T &x) const { return real(x); } };
template<int N,class V> struct imag_function<Vector<simd_vector_generator<N,V> > > : public unary_value_function<Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<N,typename complex_traits<V>::value_type> > > { template<class T> typename imag_function<V>::const_reference operator()(const T &x) const { return imag(x); } };
template<int N,class V> struct norm_binary_function<Vector<simd_vector_generator<N,V> > > : public binary_value_function<Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<2*N,typename complex_traits<V>::value_type> > > { template<class T1,class T2> typename norm_binary_function<V>::const_reference operator()(const T1 &x,const T2 &y) const { return norm(x,y); } };
template<int N,class V> struct abs_binary_function <Vector<simd_vector_generator<N,V> > > : public binary_value_function<Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<N,V> >,Vector<simd_vector_generator<2*N,typename complex_traits<V>::value_type> > > { template<class T1,class T2> typename abs_binary_function <V>::const_reference operator()(const T1 &x,const T2 &y) const { return abs (x,y); } };





template<class T> struct matrixrow_rebind {};
template<class V,class A>          struct matrixrow_rebind<Matrix<dense_matrix_generator<V,A  > > > { typedef typename DenseVector<V,A>::self other; };
template<int M,int N,class V>      struct matrixrow_rebind<Matrix< tiny_matrix_generator<M,N,V> > > { typedef typename TinyVector<N,V>::self other; };
template<class V,class A,int Copy> struct matrixrow_rebind<Matrix<shift_array_generator      <Matrix<dense_matrix_generator<V,A> >,Copy> > > { typedef typename shiftDenseVector<V,A>::self other; };
template<class V,class A,int Copy> struct matrixrow_rebind<Matrix<shift_dense_array_generator<Matrix<dense_matrix_generator<V,A> >,Copy> > > { typedef typename shiftDenseVector<V,A>::self other; };
template<int M,int N,class V,int Copy> struct matrixrow_rebind<Matrix<shift_array_generator      <Matrix<tiny_matrix_generator<M,N,V> >,Copy> > > { typedef typename shiftTinyVector<N,V>::self other; };
template<int M,int N,class V,int Copy> struct matrixrow_rebind<Matrix<shift_dense_array_generator<Matrix<tiny_matrix_generator<M,N,V> >,Copy> > > { typedef typename shiftTinyVector<N,V>::self other; };

template<class A,int C=0>
class matrixrow_generator : public array_generator<A,C>
{
  public:
    typedef matrixrow_generator self;
    typedef array_generator<A> base;
    typedef typename base::array_type      array_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::int_type int_type;
    typedef typename base::int_type size_type;
    typedef typename base::int_type index_type;

    typedef vector_array_tag array_category;

    template<class V2> struct array_rebind { typedef typename matrixrow_rebind<typename array_type::template rebind<V2>::other>::other other; };

  private:
    int_type pos;

  public:
    inline matrixrow_generator(array_type &x, const int_type &i) : base(x), pos(i) {}
    inline matrixrow_generator(const self &r) : base(r), pos(r.pos) {}

    inline reference       operator[](const index_type &n)       { return array()(pos,n); }
    inline const_reference operator[](const index_type &n) const { return array()(pos,n); }

    inline void inv(const value_type &y, const int_type &n) { (*this)[n]=y; }

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }
    inline index_type lower_bound() const { return array().col_lower_bound(); }
    inline size_type size () const { return array().ncols(); }
    inline int_type  nelms() const { return array().ncols(); }
    inline void resize(int_type n) { assert(size()==n); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
};

template<class T,int C>
class row_dense_matrix_generator : public matrixrow_generator<T,C>
{
  public:
    typedef row_dense_matrix_generator self;
    typedef matrixrow_generator<T> base;
    ARRAY_BASE_TYPES

    typedef vector_array_tag array_category;

  protected:
    pointer data;

  public:
    inline row_dense_matrix_generator(array_type &x, const int_type &i) : base(x,i), data(&x(i,0)) {}
    inline row_dense_matrix_generator(const self &x) : base((base &)x), data(x.data) {}

    inline reference       operator[](const index_type &n)       { return data[n]; }
    inline const_reference operator[](const index_type &n) const { return data[n]; }
};

template<int N,class T,int C>
class tiny_row_dense_matrix_generator : public row_dense_matrix_generator<T,C>
{
  public:
    typedef tiny_row_dense_matrix_generator self;
    typedef row_dense_matrix_generator<T,C> base;
    ARRAY_BASE_TYPES

    typedef vector_array_tag array_category;

  public:
    inline tiny_row_dense_matrix_generator(array_type &x, const int_type &i) : base(x,i) {}
    inline tiny_row_dense_matrix_generator(const self &x) : base((base &)x) {}

    inline size_type size () const { return N; }
};

template<class T,int N,int C>
class row_simd_dense_matrix_generator : public matrixrow_generator<T,C>
{
  public:
    typedef row_simd_dense_matrix_generator self;
    typedef matrixrow_generator<T> base;
    ARRAY_BASE_TYPES

    typedef vector_array_tag array_category;
    typedef typename array_type::generator_type::template_reference template_reference;

  public:
    typename array_traits<reference>::pointer data;

  public:
    inline row_simd_dense_matrix_generator(array_type &x, const int_type &i) : base(x,i), data(x.generator().data+i*x.generator().row_stride()) {}
    inline row_simd_dense_matrix_generator(const self &x) : base((base &)x), data(x.data) {}

    inline reference       operator[](const index_type &i)       { return template_reference((typename template_reference::pointer)(data+N*i)); }
    inline const_reference operator[](const index_type &i) const { return template_reference((typename template_reference::pointer)(data+N*i)); }
};

template<class T,int C> struct MatrixRow { typedef T array_type; typedef matrixrow_generator<array_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

#ifdef NDEBUG
template<class V             ,int C> struct MatrixRow<      Matrix<data_matrix_generator          <V     > >,C> { typedef typename GeneratorArray<row_dense_matrix_generator     <      Matrix<data_matrix_generator          <V     > >  ,C> >::self self; };
template<class V             ,int C> struct MatrixRow<const Matrix<data_matrix_generator          <V     > >,C> { typedef typename GeneratorArray<row_dense_matrix_generator     <const Matrix<data_matrix_generator          <V     > >  ,C> >::self self; };
template<class V,class A     ,int C> struct MatrixRow<      Matrix<dense_matrix_generator         <V,A   > >,C> { typedef typename GeneratorArray<row_dense_matrix_generator     <      Matrix<dense_matrix_generator         <V,A   > >  ,C> >::self self; };
template<class V,class A     ,int C> struct MatrixRow<const Matrix<dense_matrix_generator         <V,A   > >,C> { typedef typename GeneratorArray<row_dense_matrix_generator     <const Matrix<dense_matrix_generator         <V,A   > >  ,C> >::self self; };
template<class T      ,int C2,int C> struct MatrixRow<      Matrix<sub_dense_matrix_generator     <T  ,C2> >,C> { typedef typename GeneratorArray<row_dense_matrix_generator     <      Matrix<sub_dense_matrix_generator     <T  ,C2> >  ,C> >::self self; };
template<class T      ,int C2,int C> struct MatrixRow<const Matrix<sub_dense_matrix_generator     <T  ,C2> >,C> { typedef typename GeneratorArray<row_dense_matrix_generator     <const Matrix<sub_dense_matrix_generator     <T  ,C2> >  ,C> >::self self; };
template<class T,int N,int C2,int C> struct MatrixRow<      Matrix<sub_simd_dense_matrix_generator<T,N,C2> >,C> { typedef typename GeneratorArray<row_simd_dense_matrix_generator<      Matrix<sub_simd_dense_matrix_generator<T,N,C2> >,N,C> >::self self; };
template<class T,int N,int C2,int C> struct MatrixRow<const Matrix<sub_simd_dense_matrix_generator<T,N,C2> >,C> { typedef typename GeneratorArray<row_simd_dense_matrix_generator<const Matrix<sub_simd_dense_matrix_generator<T,N,C2> >,N,C> >::self self; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct MatrixRow<      Matrix<simd_sub_dense_matrix_generator  <T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<row_simd_dense_matrix_generator<      Matrix<simd_sub_dense_matrix_generator  <T,N,Cast,R,C2> >,N,C> >::self self; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct MatrixRow<const Matrix<simd_sub_dense_matrix_generator  <T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<row_simd_dense_matrix_generator<const Matrix<simd_sub_dense_matrix_generator  <T,N,Cast,R,C2> >,N,C> >::self self; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct MatrixRow<      Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<row_simd_dense_matrix_generator<      Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,N,C> >::self self; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct MatrixRow<const Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<row_simd_dense_matrix_generator<const Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,N,C> >::self self; };
#endif

#ifdef NDEBUG
template<int M,int N,class V       ,int C> struct MatrixRow<      Matrix<tiny_matrix_generator          <M,N,V   > >,C> { typedef typename GeneratorArray<tiny_row_dense_matrix_generator<N,      Matrix<tiny_matrix_generator          <M,N,V   > >,C> >::self self; };
template<int M,int N,class V       ,int C> struct MatrixRow<const Matrix<tiny_matrix_generator          <M,N,V   > >,C> { typedef typename GeneratorArray<tiny_row_dense_matrix_generator<N,const Matrix<tiny_matrix_generator          <M,N,V   > >,C> >::self self; };
template<int M,int N,class T,int C2,int C> struct MatrixRow<      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> { typedef typename GeneratorArray<tiny_row_dense_matrix_generator<N,      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> >::self self; };
template<int M,int N,class T,int C2,int C> struct MatrixRow<const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> { typedef typename GeneratorArray<tiny_row_dense_matrix_generator<N,const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> >::self self; };
#endif

//{unsecret}
//{noAutoLink}
//{group: Matrices functions}
//Summary : Isolates a row
//Parameters:
//  X : The matrix
//  n : Position of the row
//Return: A vector representing the n-th row of X
//Example:
//  DenseMatrix<int>::self X(2,2,"1 2 3 4");
//  cout << row(X,0) << endl; // [1 2]
//See: ^rows^, ^col^, ^cols^
template<class G> inline typename MatrixRow<      Matrix<G> >::self row(      Matrix<G> &X, typename Matrix<G>::int_type n) { return typename MatrixRow<      Matrix<G> >::self(X,n); }
//{unsecret}
template<class G> inline typename MatrixRow<const Matrix<G> >::self row(const Matrix<G> &X, typename Matrix<G>::int_type n) { return typename MatrixRow<const Matrix<G> >::self(X,n); }



template<class T> struct matrixcol_rebind {};
template<class V,class A>          struct matrixcol_rebind<Matrix<dense_matrix_generator<V,A  > > > { typedef typename DenseVector<V,A>::self other; };
template<int M,int N,class V>      struct matrixcol_rebind<Matrix< tiny_matrix_generator<M,N,V> > > { typedef typename TinyVector<M,V>::self other; };
template<class V,class A,int Copy> struct matrixcol_rebind<Matrix<shift_array_generator      <Matrix<dense_matrix_generator<V,A> >,Copy> > > { typedef typename shiftDenseVector<V,A>::self other; };
template<class V,class A,int Copy> struct matrixcol_rebind<Matrix<shift_dense_array_generator<Matrix<dense_matrix_generator<V,A> >,Copy> > > { typedef typename shiftDenseVector<V,A>::self other; };
template<int M,int N,class V,int Copy> struct matrixcol_rebind<Matrix<shift_array_generator      <Matrix<tiny_matrix_generator<M,N,V> >,Copy> > > { typedef typename shiftTinyVector<M,V>::self other; };
template<int M,int N,class V,int Copy> struct matrixcol_rebind<Matrix<shift_dense_array_generator<Matrix<tiny_matrix_generator<M,N,V> >,Copy> > > { typedef typename shiftTinyVector<M,V>::self other; };

template<class T,int C=0>
class matrixcol_generator : public array_generator<T,C>
{
  public:
    typedef matrixcol_generator self;
    typedef array_generator<T,C> base;
    typedef typename base::array_type      array_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::int_type int_type;
    typedef typename base::int_type size_type;
    typedef typename base::int_type index_type;

    typedef vector_array_tag array_category;

    template<class V2> struct array_rebind { typedef typename matrixcol_rebind<typename array_type::template rebind<V2>::other>::other other; };

  protected:
    int_type ind;

  public:
    inline matrixcol_generator(array_type &x, const int_type &i) : base(x), ind(i) {}
    inline matrixcol_generator(array_type &x, const self &r) : base(x,r), ind(r.ind) {}
    inline matrixcol_generator(const self &r) : base(r), ind(r.ind) {}


    inline reference       operator[](const index_type &n)       { return array()(n,ind); }
    inline const_reference operator[](const index_type &n) const { return array()(n,ind); }

    inline void inv(const value_type &y, const int_type &n) { (*this)[n]=y; }

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }

    inline index_type lower_bound() const { return array().row_lower_bound(); }
    inline size_type size () const { return array().nrows(); }
    inline int_type  nelms() const { return array().nrows(); }
    inline void resize(const size_type &s) { assert(size()==s); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    index_type pos() const { return ind; }
};

template<class T,int C=0>
struct MatrixCol : public array_traits<T>
{
  typedef T array_type;
  typedef matrixcol_generator<array_type,C> generator_type;
  typedef Vector<generator_type> self;
};

//{unsecret}
//{group: Matrices functions}
//Summary : Isolates a column
//Parameters:
//  X : The matrix
//  n : Position of the column
//Return: A vector representing the n-th column of X
//Example:
//  DenseMatrix<int>::self X(2,2,"1 2 3 4");
//  cout << col(X,0) << endl; // [1 3]
//See: ^cols^, ^row^, ^rows^
template<class G> inline typename MatrixCol<      Matrix<G> >::self col(      Matrix<G> &X, typename Matrix<G>::int_type n) { return typename MatrixCol<      Matrix<G> >::self(X,n); }
//{unsecret}
template<class G> inline typename MatrixCol<const Matrix<G> >::self col(const Matrix<G> &X, typename Matrix<G>::int_type n) { return typename MatrixCol<const Matrix<G> >::self(X,n); }


template<class T,int C>
class col_dense_matrix_generator : public matrixcol_generator<T,C>
{
  public:
    typedef col_dense_matrix_generator self;
    typedef matrixcol_generator<T,C> base;
    ARRAY_BASE_TYPES

    typedef vector_array_tag array_category;

  protected:
    pointer data;
    int_type s;

  public:
    inline col_dense_matrix_generator(array_type &x, const int_type &i) : base(x,i), data(&x(0,i)), s(array().generator().row_stride()) {}
    inline col_dense_matrix_generator(array_type &x,const self &r) : base(x,(base &)r), data(&x(0,r.ind)), s(array().generator().row_stride()) {}
    inline col_dense_matrix_generator(const self &r) : base((base &)r), data(r.data), s(r.s) {}

    using base::array;

    inline reference       operator[](const index_type &n)       { return data[stride()*n]; }
    inline const_reference operator[](const index_type &n) const { return data[stride()*n]; }

    inline int_type stride() const { return array().generator().row_stride(); }
    //inline int_type stride() const { return s; }
};

#ifdef NDEBUG
template<            class V        ,int C> struct MatrixCol<      Matrix<data_matrix_generator          <    V   > >,C> { typedef typename GeneratorArray<col_dense_matrix_generator<      Matrix<data_matrix_generator          <    V   > >,C> >::self self; };
template<            class V        ,int C> struct MatrixCol<const Matrix<data_matrix_generator          <    V   > >,C> { typedef typename GeneratorArray<col_dense_matrix_generator<const Matrix<data_matrix_generator          <    V   > >,C> >::self self; };
template<            class V,class A,int C> struct MatrixCol<      Matrix<dense_matrix_generator         <    V,A > >,C> { typedef typename GeneratorArray<col_dense_matrix_generator<      Matrix<dense_matrix_generator         <    V,A > >,C> >::self self; };
template<            class V,class A,int C> struct MatrixCol<const Matrix<dense_matrix_generator         <    V,A > >,C> { typedef typename GeneratorArray<col_dense_matrix_generator<const Matrix<dense_matrix_generator         <    V,A > >,C> >::self self; };
template<            class T,int C2 ,int C> struct MatrixCol<      Matrix<sub_dense_matrix_generator     <    T,C2> >,C> { typedef typename GeneratorArray<col_dense_matrix_generator<      Matrix<sub_dense_matrix_generator     <    T,C2> >,C> >::self self; };
template<            class T,int C2 ,int C> struct MatrixCol<const Matrix<sub_dense_matrix_generator     <    T,C2> >,C> { typedef typename GeneratorArray<col_dense_matrix_generator<const Matrix<sub_dense_matrix_generator     <    T,C2> >,C> >::self self; };
#endif



template<int N,class T,int C>
class tiny_col_dense_matrix_generator : public col_dense_matrix_generator<T,C>
{
  public:
    typedef tiny_col_dense_matrix_generator self;
    typedef col_dense_matrix_generator<T,C> base;
    ARRAY_BASE_TYPES

    typedef vector_array_tag array_category;

  public:
    inline tiny_col_dense_matrix_generator(array_type &x, const int_type &i) : base(x,i) {}
    inline tiny_col_dense_matrix_generator(const self &x) : base((base &)x) {}

    inline size_type size () const { return N; }
};

#ifdef NDEBUG
template<int M,int N,class V       ,int C> struct MatrixCol<      Matrix<tiny_matrix_generator          <M,N,V   > >,C> { typedef typename GeneratorArray<tiny_col_dense_matrix_generator<M,      Matrix<tiny_matrix_generator          <M,N,V   > >,C> >::self self; };
template<int M,int N,class V       ,int C> struct MatrixCol<const Matrix<tiny_matrix_generator          <M,N,V   > >,C> { typedef typename GeneratorArray<tiny_col_dense_matrix_generator<M,const Matrix<tiny_matrix_generator          <M,N,V   > >,C> >::self self; };
template<int M,int N,class T,int C2,int C> struct MatrixCol<      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> { typedef typename GeneratorArray<tiny_col_dense_matrix_generator<M,      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> >::self self; };
template<int M,int N,class T,int C2,int C> struct MatrixCol<const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> { typedef typename GeneratorArray<tiny_col_dense_matrix_generator<M,const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,C> >::self self; };
#endif



template<int N,class T> struct tiny_block_matrix_col_generator_traits            { typedef typename TinySubVector<N,      typename MatrixRow<      T>::self,1>::self value_type; };
template<int N,class T> struct tiny_block_matrix_col_generator_traits<N,const T> { typedef typename TinySubVector<N,/*const*/ typename MatrixRow<const T>::self,1>::self value_type; };

template<int N,class T>
class tiny_block_matrix_col_generator : public array_reference_generator<T,typename tiny_block_matrix_col_generator_traits<N,T>::value_type,typename tiny_block_matrix_col_generator_traits<N,T>::value_type,typename tiny_block_matrix_col_generator_traits<N,const T>::value_type>
{
  public:
    typedef tiny_block_matrix_col_generator self;
    typedef array_reference_generator<T,typename tiny_block_matrix_col_generator_traits<N,T>::value_type,typename tiny_block_matrix_col_generator_traits<N,T>::value_type,typename tiny_block_matrix_col_generator_traits<N,const T>::value_type> base;
    typedef typename base::array_type      array_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::int_type int_type;
    typedef typename base::int_type size_type;
    typedef typename base::int_type index_type;

    typedef vector_array_tag array_category;

    template<class V2> struct array_rebind { typedef typename matrixrow_rebind<typename array_type::template rebind<V2>::other>::other other; };

  private:
    int_type ind;

  public:
    inline tiny_block_matrix_col_generator(array_type &x, const int_type &i) : base(x), ind(i) {}
    inline tiny_block_matrix_col_generator(const self &x) : base(x), ind(x.ind) {}

    inline reference       operator[](const index_type &n)       { return sub<N>(row(array(),n),ind*N); }
    inline const_reference operator[](const index_type &n) const { return sub<N>(row(array(),n),ind*N); }

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }

    inline index_type lower_bound() const { return array().row_lower_bound(); }
    inline size_type size () const { return array().nrows(); }
    inline int_type  nelms() const { return array().nrows(); }
    inline void resize(const size_type &s) { assert(size()==s); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
};

template<int N,class T>
struct tinyBlockMatrixCol
{
  typedef tiny_block_matrix_col_generator<N,T> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int N, class G> typename tinyBlockMatrixCol<N,      Matrix<G> >::self block_col(      Matrix<G> &X,typename Matrix<G>::int_type j) { return typename tinyBlockMatrixCol<N,      Matrix<G> >::self(X,j); }
template<int N, class G> typename tinyBlockMatrixCol<N,const Matrix<G> >::self block_col(const Matrix<G> &X,typename Matrix<G>::int_type j) { return typename tinyBlockMatrixCol<N,const Matrix<G> >::self(X,j); }


template<class A,class V=typename A::value_type>
class rotate_vector_generator : public array_value_generator<const A,V>
{
  public:
    typedef rotate_vector_generator self;
    typedef array_value_generator<const A,V> base;
    ARRAY_BASE_TYPES

  private:
    int n;

  public:
    inline rotate_vector_generator(array_type &x, int r) : base(x), n(r) {}

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }

    inline const_reference operator[](const index_type &i) const { return array()[positive_mod(i-n-array().lower_bound(),array().size())+array().lower_bound()]; }
};

template<class A,class V=typename A::value_type>
struct rotateVector
{
  typedef A array_type;
  typedef V value_type;
  typedef rotate_vector_generator<array_type, value_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename rotateVector<Vector<G>  >::self rotate(const Vector<G>  &X, int n) { return typename rotateVector<Vector<G>  >::self(X,n); }
template<class G> inline typename rotateVector<Array<1,G> >::self rotate(const Array<1,G> &X, int n) { return typename rotateVector<Array<1,G> >::self(X,n); }

template<class T, class V=typename T::value_type>
class alternate_generator : public array_value_generator<T,V>
{
  public:
    typedef alternate_generator self;
    typedef array_value_generator<T,V> base;
    ARRAY_BASE_TYPES

  public:
    inline alternate_generator(array_type &x) : base(x) {}

    using base::array;

    inline const_reference operator[](const index_type &i) const { return (i%2)?-array()[i]:array()[i]; }
};

template<class T, class V=typename T::value_type>
struct alternateVector : public array_value_traits<T,V>
{
  typedef T array_type;
  typedef V value_type;
  typedef alternate_generator<array_type,value_type> generator_type;
  typedef Vector<generator_type> self;
};

template<class G>
inline typename alternateVector<const Vector<G> >::self alternate(const Vector<G> &X)
{
  return typename alternateVector<const Vector<G> >::self(X);
}


//}

#endif
