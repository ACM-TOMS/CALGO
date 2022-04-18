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

#ifndef VECTOR_H
#define VECTOR_H

#include "genial_config.h"

#include <assert.h>
#include <numeric>
#include <iterator>
#include <iostream>
#include <algorithm>

#ifdef HAS_SSE3
#include "simd/sse3.h"
#endif
#ifdef HAS_SSE2
#include "simd/sse2.h"
#endif
#ifdef HAS_SSE
#include "simd/sse.h"
#endif
#ifdef HAS_MMX
#include "simd/mmx.h"
#endif


#include "util.h"
#include "array/genial_array.h"

#include "array/arraygenerator.h"
#include "array/arrayfunction.h"

#include "array/vectoriterator.h"
#include "array/vectorgenerator.h"



//namespace genial
//{

using namespace std;



//{unsecret}
//{group:Vectors}
//Summary:Class that represents one-dimentional arrays
//Include: vector.h
//Example:
//  template<class G>
//  typename Vector<G>::size_type size(const Vector<G> &X)
//  {
//    return X.size();
//  }
template<class G>
class Vector : public generator_traits<G>
{
  public:
    typedef Vector self;
    typedef generator_traits<G> base;
    typedef typename base::generator_type generator_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::index_type index_type;
    typedef typename base::size_type  size_type;
    typedef typename base::int_type   int_type;

    template<class V> struct rebind { typedef typename generator_type::template array_rebind<V>::other other; };

    template<class V> struct dense_rebind { typedef typename DenseVector<V>::self other; };
    template<class V> struct shift_dense_rebind { typedef typename shiftDenseVector<V>::self other; };

    typedef typename generator_type::template       iterator_rebind<self>::other iterator;
    typedef typename generator_type::template const_iterator_rebind<self>::other const_iterator;

    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<      iterator> reverse_iterator;

  protected:
    generator_type gen;

  public:
    inline Vector() : gen() {}
    inline explicit Vector(const generator_type &g) : gen(g) {}

    inline Vector(      self &x) : gen(const_cast<const generator_type &>(x.gen)) {}
    inline Vector(const self &x) : gen(x.gen) {}
    template<class G2> inline Vector(      Vector<G2> &x) : gen(x) {}
    template<class G2> inline Vector(const Vector<G2> &x) : gen(x) {}

    template<class A> inline Vector(      A &a,      self &x) : gen(a,const_cast<const generator_type &>(x.gen)) {}
    template<class A> inline Vector(      A &a,const self &x) : gen(a,x.gen) {}
    template<class A> inline Vector(const A &a,      self &x) : gen(a,const_cast<const generator_type &>(x.gen)) {}
    template<class A> inline Vector(const A &a,const self &x) : gen(a,x.gen) {}

    inline Vector(      Array<1,generator_type> &x) : gen(const_cast<const generator_type &>(x.generator())) {}
    inline Vector(const Array<1,generator_type> &x) : gen(x.generator()) {}
    template<class G2> inline Vector(      Array<1,G2> &x) : gen(x) {}
    template<class G2> inline Vector(const Array<1,G2> &x) : gen(x) {}

    template<class A> inline explicit Vector(A &a) : gen(a) {}
    template<class A,class B> inline Vector(A &a, const B &b) : gen(a,b) { }
    template<class A,class B,class C> inline Vector(A &a, const B &b, const C &c) : gen(a,b,c) {}
    template<class A,class B,class C,class D> inline Vector(A &a, const B &b, const C &c, const D &d) : gen(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> inline Vector(A &a, const B &b, const C &c, const D &d, const E &e) : gen(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> inline Vector(A &a, const B &b, const C &c, const D &d, const E &e, const F &f) : gen(a,b,c,d,e,f) {}

    template<class A> inline explicit Vector(const A &a) : gen(a) {}
    template<class A,class B> inline Vector(const A &a, const B &b) : gen(a,b) { }
    template<class A,class B,class C> inline Vector(const A &a, const B &b, const C &c) : gen(a,b,c) {}
    template<class A,class B,class C,class D> inline Vector(const A &a, const B &b, const C &c, const D &d) : gen(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> inline Vector(const A &a, const B &b, const C &c, const D &d, const E &e) : gen(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> inline Vector(const A &a, const B &b, const C &c, const D &d, const E &e, const F &f) : gen(a,b,c,d,e,f) {}

    inline self                    &operator=(const self        &Y) { resize(Y.size()); set_lower_bound(Y.lower_bound()); affect(*this,Y); return *this; }
    template<class G2> self inline &operator=(const Vector<G2>  &Y) { resize(Y.size()); set_lower_bound(Y.lower_bound()); affect(*this,Y); return *this; }
    template<class G2> self inline &operator=(const Array<1,G2> &Y) { resize(Y.size()); set_lower_bound(Y.lower_bound()); affect(*this,Y); return *this; }
    
    inline self &operator=(const value_type &v) { affect(*this,v); return *this; }

    template<class V,class A> inline explicit Vector(      list<V,A> &x) { resize(x.size()); copy(x.begin(), x.end(), begin()); }
    template<class V,class A> inline explicit Vector(const list<V,A> &x) { resize(x.size()); copy(x.begin(), x.end(), begin()); }
    template<class V,class A> inline self &operator=(      list<V,A> &x) { resize(x.size()); copy(x.begin(), x.end(), begin()); return *this; }
    template<class V,class A> inline self &operator=(const list<V,A> &x) { resize(x.size()); copy(x.begin(), x.end(), begin()); return *this; }

    //Summary: Swap the contents
    template<class G2> inline void swap(Vector<G2> &x) { ::swap(generator(),x.generator()); }

    //{secret}
    inline generator_type       &generator()       { return gen; }
    //{secret}
    inline const generator_type &generator() const { return gen; }

    //Summary: Iterator that points at the first element
    inline iterator       begin()       { return iterator(*this); }
    inline const_iterator begin() const { return const_iterator(*this); }
    //Summary: Iterator that points just beyond the last element
    inline iterator       end  ()       { return iterator(*this)+nelms(); }
    inline const_iterator end  () const { return const_iterator(*this)+nelms(); }
    //Summary: Reverse iterator that points at the last element
    inline reverse_iterator       rbegin()       { return reverse_iterator(end  ()); }
    inline const_reverse_iterator rbegin() const { return const_reverse_iterator(end  ()); }
    //Summary: Reverse iterator that points just before the first element
    inline reverse_iterator       rend  ()       { return reverse_iterator(begin()); }
    inline const_reverse_iterator rend  () const { return const_reverse_iterator(begin()); }

    //Summary: First element access
    inline reference       front()       { return *begin(); }
    inline const_reference front() const { return *begin(); }
    //Summary: Last element access
    inline reference       back ()       { return *--end(); }
    inline const_reference back () const { return *--end(); }

    //Summary: Element access
    inline reference       operator[](const index_type &n)       { return gen[n]; }
    inline const_reference operator[](const index_type &n) const { return gen[n]; }

    //{noAutoLink}
    //Summary: Number of elements
    inline size_type size () const { return gen.size(); }
    //Summary: Number of elements, same result than size
    inline int_type  nelms() const { return size(); }
    //Summary: Lowest valid position
    inline index_type lower_bound() const { return gen.lower_bound(); }
    //Summary: Greatest valid position
    inline index_type upper_bound() const { return size()+(lower_bound()-1); }

    //Summary: Tests the validity of the position
    inline bool withinbounds  (const index_type &i) const { return (i>=lower_bound()) && (i<=upper_bound()); }

    //inline int_type nrows() const { return size(); }
    //inline int_type ncols() const { return 1; }

    //Summary: Changes the size
    inline void resize(const size_type &s) { gen.resize(s); }
    //Summary: Changes the lowest valid position
    //Remarks: Only for shifted array
    inline void set_lower_bound(const index_type &i) { gen.set_lower_bound(i); }

    //inline bool operator==(const self &r) const { return gen==r.gen; }
    //inline bool operator!=(const self &r) const { return gen!=r.gen; }
};

template<class G>
class Array<1,G> : public Vector<G>
{
  public:
    typedef Array self;
    typedef Vector<G> base;
    typedef typename base::generator_type generator_type;

  public:
    Array() : base() {}
    explicit Array(const generator_type &g) : base(g) {}

    Array(const self &r) : base(r) {}

    template<class A> explicit Array(A &a) : base(a) {}
    template<class A,class B> Array(A &a, const B &b) : base(a,b) {}
    template<class A,class B,class C> Array(A &a, const B &b, const C &c) : base(a,b,c) {}
    template<class A,class B,class C,class D> Array(A &a, const B &b, const C &c, const D &d) : base(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> Array(A &a, const B &b, const C &c, const D &d, const E &e) : base(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> Array(A &a, const B &b, const C &c, const D &d, const E &e,const F &f) : base(a,b,c,d,e,f) {}

    template<class A> explicit Array(const A &a) : base(a) {}
    template<class A,class B> Array(const A &a, const B &b) : base(a,b) {}
    template<class A,class B,class C> Array(const A &a, const B &b, const C &c) : base(a,b,c) {}
    template<class A,class B,class C,class D> Array(const A &a, const B &b, const C &c, const D &d) : base(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> Array(const A &a, const B &b, const C &c, const D &d,const E &e) : base(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> Array(const A &a, const B &b, const C &c, const D &d,const E &e,const F &f) : base(a,b,c,d,e,f) {}

    self &operator=(const self &r) { base::operator=(r); return *this; }
    template<class A> self &operator=(const A &a) { base::operator=(a); return *this; }
};

//Group = Array functions

//{unsecret}
//Summary: Conversions
template<class E,class Tr,class G> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const Vector<G> &X)
{
  typedef typename Vector<G>::const_iterator const_iterator;

  os << "% " << X.size();
  if (X.lower_bound()!=0) os << "(" << X.lower_bound() << ")";
  os << endl;
  for (const_iterator it=X.begin(); it!=X.end(); ++it)
    os << *it << " ";
  os << endl;
  return os;
}

//{unsecret}
//Summary: Conversions
template<class E,class Tr,class G> basic_istream<E,Tr> &operator>>(basic_istream<E,Tr> &is, Vector<G> &X)
{
  typedef typename Vector<G>::int_type int_type;
  typedef typename Vector<G>::size_type size_type;

  char c; is >> c && c=='%';

  size_type d; is >> d; X.resize(d);
  for (int_type i=0; i<X.size(); ++i) is >> X[i];
  return is;
}

template<class E,class T,int N,class V> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const Vector<tiny_vector_generator<N,V> > &X)
{
  os << '[';
  for (int i=0; i<N; ++i)
    os << X[i] << " ";
  os << ']';
  return os;
}

template<class E,class T,int N,class V> basic_istream<E,T> &operator>>(basic_istream<E,T> &is, Vector<tiny_vector_generator<N,V> > &X)
{
  char c;
  is >> c && c=='[';
  for (int i=0; i<N; ++i) is >> X[i];
  is >> c && c==']';
  return is;
}

template<class FwdIt,class V,class A> inline void uninitialized_fill_aux(FwdIt begin, FwdIt end, const Vector<dense_vector_generator<V,A> > &v) { uninitialized_fill(begin,end,v); }


//Group = Vectors functions

//{unsecret}
//Summary: Inner product
template<class G1,class G2> inline typename promotion2_traits<typename Vector<G1>::value_type, typename Vector<G2>::value_type>::value_type inner_product(const Vector<G1> &X, const Vector<G2> &Y)
{
  assert(X.size()==Y.size());
  return inner_product_n(X.size(), X.begin(), Y.begin(), PROMOTE2(typename Vector<G1>::value_type, typename Vector<G2>::value_type)(0));
}


//Group = Arrays functions

#ifdef MMX
inline SimdVector<4,short>::self unpacklo(const SimdVector<8,unsigned char>::self &X) { return SimdVector<4,short>::self(unpacklo(X.generator().data)); }
inline SimdVector<4,short>::self unpackhi(const SimdVector<8,unsigned char>::self &X) { return SimdVector<4,short>::self(unpackhi(X.generator().data)); }
inline SimdVector<2,int  >::self madd    (const SimdVector<4,short        >::self &X,const SimdVector<4,short        >::self &Y) { return SimdVector<2,int  >::self(madd (X.generator().data,Y.generator().data)); }
inline SimdVector<4,short>::self sad     (const SimdVector<8,unsigned char>::self &X,const SimdVector<8,unsigned char>::self &Y) { return SimdVector<4,short>::self(sad  (X.generator().data,Y.generator().data)); }

inline SimdVector<8,char         >::self packs (const SimdVector<4,short>::self &X,const SimdVector<4,short>::self &Y) { return SimdVector<8,char         >::self(packs (X.generator().data,Y.generator().data)); }
inline SimdVector<8,unsigned char>::self packus(const SimdVector<4,short>::self &X,const SimdVector<4,short>::self &Y) { return SimdVector<8,unsigned char>::self(packus(X.generator().data,Y.generator().data)); }
inline SimdVector<4,short        >::self packs (const SimdVector<2,int  >::self &X,const SimdVector<2,int  >::self &Y) { return SimdVector<4,short        >::self(packs (X.generator().data,Y.generator().data)); }

#endif //MMX

#ifdef SSE2
inline SimdVector<8,short>::self unpacklo(const SimdVector<16,unsigned char>::self &X) { return SimdVector<8,short>::self(unpacklo(X.generator().data)); }
inline SimdVector<8,short>::self unpackhi(const SimdVector<16,unsigned char>::self &X) { return SimdVector<8,short>::self(unpackhi(X.generator().data)); }
inline SimdVector<4,int  >::self madd    (const SimdVector< 8,short        >::self &X,const SimdVector< 8,short        >::self &Y) { return SimdVector<4,int  >::self(madd (X.generator().data,Y.generator().data)); }
inline SimdVector<8,short>::self sad     (const SimdVector<16,unsigned char>::self &X,const SimdVector<16,unsigned char>::self &Y) { return SimdVector<8,short>::self(sad  (X.generator().data,Y.generator().data)); }

inline SimdVector<16,char         >::self packs (const SimdVector<8,short>::self &X,const SimdVector<8,short>::self &Y) { return SimdVector<16,char         >::self(packs (X.generator().data,Y.generator().data)); }
inline SimdVector<16,unsigned char>::self packus(const SimdVector<8,short>::self &X,const SimdVector<8,short>::self &Y) { return SimdVector<16,unsigned char>::self(packus(X.generator().data,Y.generator().data)); }
inline SimdVector< 8,short        >::self packs (const SimdVector<4,int  >::self &X,const SimdVector<4,int  >::self &Y) { return SimdVector< 8,short        >::self(packs (X.generator().data,Y.generator().data)); }
#endif // SSE2


//} //namespace genial

#endif

