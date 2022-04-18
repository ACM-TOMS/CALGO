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

#ifndef MATRIX_H
#define MATRIX_H

#include "array/vector.h"
#include "array/matrixiterator.h"
#include "array/matrixindex.h"
#include "array/matrixsize.h"
#include "array/matrixgenerator.h"

//namespace genial
//{

using namespace std;

//{unsecret}
//{group:Matrices}
//Summary:Class that represents two-dimentional arrays
//Includes: matrix
//Example:
//  template<class G>
//  typename Matrix<G>::size_type size(const Matrix<G> &X)
//  {
//    return X.size();
//  }
template<class G>
class Matrix : public generator_traits<G>
{
  public:
    typedef Matrix self;
    typedef generator_traits<G> base;
    typedef typename base::generator_type  generator_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::index_type index_type;
    typedef typename base::size_type  size_type;
    typedef typename base::int_type   int_type;

    template<class V> struct rebind { typedef typename generator_type::template array_rebind<V>::other other; };

    template<class V> struct dense_rebind { typedef typename DenseMatrix<V>::self other; };
    template<class V> struct shift_dense_rebind { typedef typename shiftDenseMatrix<V>::self other; };

    typedef typename generator_type::template       iterator_rebind<self>::other iterator;
    typedef typename generator_type::template const_iterator_rebind<self>::other const_iterator;

    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<      iterator> reverse_iterator;

  protected:
    generator_type gen;

  public:
    inline Matrix() : gen() {}
    inline  explicit Matrix(const generator_type &g) : gen(g) {}

    inline Matrix(      self &x) : gen(const_cast<const generator_type &>(x.gen)) {}
    inline Matrix(const self &x) : gen(x.gen) {}
    template<class G2> inline Matrix(      Matrix<G2>  &x) : gen(x) {}
    template<class G2> inline Matrix(const Matrix<G2>  &x) : gen(x) {}
    
    inline Matrix(      Array<2,generator_type> &x) : gen(const_cast<const generator_type &>(x.generator())) {}
    inline Matrix(const Array<2,generator_type> &x) : gen(x.generator()) {}
    template<class G2> inline Matrix(      Array<2,G2> &x) : gen(x) {}
    template<class G2> inline Matrix(const Array<2,G2> &x) : gen(x) {}

    template<class A> inline explicit Matrix(A &a) : gen(a) {}
    template<class A,class B> inline Matrix(A &a, const B &b) : gen(a,b) {}
    template<class A,class B,class C> inline Matrix(A &a, const B &b, const C &c) : gen(a,b,c) {}
    template<class A,class B,class C,class D> inline Matrix(A &a, const B &b, const C &c, const D &d) : gen(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> inline Matrix(A &a, const B &b, const C &c, const D &d, const E &e) : gen(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> inline Matrix(A &a, const B &b, const C &c, const D &d, const E &e, const F &f) : gen(a,b,c,d,e,f) {}

    template<class A> inline explicit Matrix(const A &a) : gen(a) {}
    template<class A,class B> inline Matrix(const A &a, const B &b) : gen(a,b) {}
    template<class A,class B,class C> inline Matrix(const A &a, const B &b, const C &c) : gen(a,b,c) {}
    template<class A,class B,class C,class D> inline Matrix(const A &a, const B &b, const C &c, const D &d) : gen(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> inline Matrix(const A &a, const B &b, const C &c, const D &d, const E &e) : gen(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> inline Matrix(const A &a, const B &b, const C &c, const D &d, const E &e, const F &f) : gen(a,b,c,d,e,f) {}

    inline self                    &operator=(const self        &Y) { resize(Y.size()); set_lower_bound(Y.lower_bound()); affect(*this,Y); return *this; }
    template<class G2> inline self &operator=(const Matrix<G2>  &Y) { resize(Y.size()); set_lower_bound(Y.lower_bound()); affect(*this,Y); return *this; }
    template<class G2> inline self &operator=(const Array<2,G2> &Y) { resize(Y.size()); set_lower_bound(Y.lower_bound()); affect(*this,Y); return *this; }

    inline self &operator=(const value_type &v) { affect(*this,v); return *this; }

    //Summary: Swap the contents
    inline void swap(self &x) { ::swap(generator(),x.generator()); }

    //{secret}
    inline generator_type       &generator()       { return gen; }
    //{secret}
    inline const generator_type &generator() const { return gen; }

    //Summary: Iterator that points at the first element
    inline iterator       begin()       { return iterator      (*this); }
    inline const_iterator begin() const { return const_iterator(*this); }
    //Summary: Iterator that points just beyond the last element
    inline iterator       end  ()       { return iterator      (*this)+nelms(); }
    inline const_iterator end  () const { return const_iterator(*this)+nelms(); }
    //Summary: Reverse iterator that points at the last element
    inline reverse_iterator       rbegin()       { return reverse_iterator(end  ()); }
    inline const_reverse_iterator rbegin() const { return const_reverse_iterator(end  ()); }
    //Summary: Reverse iterator that points just before the the first element
    inline reverse_iterator       rend  ()       { return reverse_iterator(begin()); }
    inline const_reverse_iterator rend  () const { return const_reverse_iterator(begin()); }

    //Summary: First element access
    inline reference       front()       { return *begin(); }
    inline const_reference front() const { return *begin(); }
    //Summary: Last element access
    inline reference       back ()       { return *--end(); }
    inline const_reference back () const { return *--end(); }

    //{noAutoLink}
    //Summary: Two-dimensional size, height & width
    inline size_type size () const { return gen.size(); }
    //Summary: Number of rows, height
    inline int_type nrows () const { return size().nrows(); }
    inline int_type height() const { return size().nrows(); }
    //Summary: Number of columns, width
    inline int_type ncols () const { return size().ncols(); }
    inline int_type width () const { return size().ncols(); }
    //Summary: Number of elements
    inline int_type nelms () const { return size().nelms(); }

    //Summary: Lowest valid position
    inline index_type   lower_bound() const { return gen.lower_bound(); }
    //Summary: Lowest valid row position
    inline int_type row_lower_bound() const { return lower_bound().i; }
    //Summary: Lowest valid comlumn position
    inline int_type col_lower_bound() const { return lower_bound().j; }
    //Summary: Greatest valid position
    inline index_type   upper_bound() const { return (lower_bound()-1)+size(); }
    //Summary: Greatest valid row position
    inline int_type row_upper_bound() const { return upper_bound().i; }
    //Summary: Greatest valid comlumn position
    inline int_type col_upper_bound() const { return upper_bound().j; }

    //Summary: Tests the validity of the position
    inline bool withinbounds  (const index_type &i)    const { return size().withinbounds  (i-lower_bound()); }
    inline bool withinbounds  (int_type i, int_type j) const { return withinbounds  (index_type(i,j)); }
    //Summary: Tests if the position is on the left
    inline bool leftboundary  (const index_type &i)    const { return size().leftboundary  (i-lower_bound()); }
    inline bool leftboundary  (int_type i, int_type j) const { return leftboundary  (index_type(i,j)); }
    //Summary: Tests if the position is on the right
    inline bool rightboundary (const index_type &i)    const { return size().rightboundary (i-lower_bound()); }
    inline bool rightboundary (int_type i, int_type j) const { return rightboundary (index_type(i,j)); }
    //Summary: Tests if the position is on the top
    inline bool topboundary   (const index_type &i)    const { return size().topboundary   (i-lower_bound()); }
    inline bool topboundary   (int_type i, int_type j) const { return topboundary   (index_type(i,j)); }
    //Summary: Tests if the position is on the bottom
    inline bool bottomboundary(const index_type &i)    const { return size().bottomboundary(i-lower_bound()); }
    inline bool bottomboundary(int_type i, int_type j) const { return bottomboundary(index_type(i,j)); }
    //Summary: Tests if the position is on any boundary
    inline bool boundary      (const index_type &i)    const { return size().boundary      (i-lower_bound()); }
    inline bool boundary      (int_type i, int_type j) const { return boundary      (index_type(i,j)); }

    //Summary: Element access
    inline reference       operator()(int_type i, int_type j)       { return (*this)[index_type(i,j)]; }
    inline const_reference operator()(int_type i, int_type j) const { return (*this)[index_type(i,j)]; }
    //Summary: Element access
    inline reference       operator[](const index_type &p)       { return gen[p]; }
    inline const_reference operator[](const index_type &p) const { return gen[p]; }

    //Summary: Changes the size
    inline void resize(int_type m, int_type n) { resize(size_type(m,n)); }
    inline void resize(const size_type &d)     { gen.resize(d); }
    //Summary: Changes the lowest valid position
    //Remarks: Only for shifted array
    inline void set_lower_bound(int_type m, int_type n) { gen.set_lower_bound(index_type(m,n)); }
    inline void set_lower_bound(const index_type &i)    { gen.set_lower_bound(i); }

    //inline bool operator==(const self &r) const { return gen==r.gen; }
    //inline bool operator!=(const self &r) const { return gen!=r.gen; }
};

template<class T>
struct matrix_traits
{
  typedef T array_type;
  typedef typename array_type::value_type      value_type;
  typedef typename array_type::reference       reference;
  typedef typename array_type::const_reference const_reference;
  typedef typename array_type::pointer         pointer;
  typedef typename array_type::const_pointer   const_pointer;
};

template<class G>
struct matrix_traits<const Matrix<G> >
{
  typedef const Matrix<G> array_type;
  typedef typename array_type::value_type      value_type;
  typedef typename array_type::const_reference reference;
  typedef typename array_type::const_reference const_reference;
  typedef typename array_type::const_pointer   pointer;
  typedef typename array_type::const_pointer   const_pointer;
};

template<class G>
class Array<2,G> : public Matrix<G>
{
  public:
    typedef Array self;
    typedef Matrix<G> base;
    typedef typename base::generator_type generator_type;

  public:
    Array() : base() {}
    Array(const generator_type &g) : base(g) {}

    Array(const self &r) : base(r) {}

    template<class A> explicit Array(A &a) : base(a) {}
    template<class A,class B> Array(A &a, const B &b) : base(a,b) {}
    template<class A,class B,class C> Array(A &a, const B &b, const C &c) : base(a,b,c) {}
    template<class A,class B,class C,class D> Array(A &a, const B &b, const C &c, const D &d) : base(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> Array(A &a, const B &b, const C &c, const D &d, const E &e) : base(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> Array(A &a, const B &b, const C &c, const D &d, const E &e, const F &f) : base(a,b,c,d,e,f) {}

    template<class A> explicit Array(const A &a) : base(a) {}
    template<class A,class B> Array(const A &a, const B &b) : base(a,b) {}
    template<class A,class B,class C> Array(const A &a, const B &b, const C &c) : base(a,b,c) {}
    template<class A,class B,class C,class D> Array(const A &a, const B &b, const C &c, const D &d) : base(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> Array(const A &a, const B &b, const C &c, const D &d, const E &e) : base(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F> Array(const A &a, const B &b, const C &c, const D &d, const E &e, const F &f) : base(a,b,c,d,e,f) {}

    self &operator=(const self &r) { base::operator=(r); return *this; }
    template<class A> self &operator=(const A &a) { base::operator=(a); return *this; }
};

//Group=Arrays functions

//{unsecret}
template<class E,class Tr, class G> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const Matrix<G> &X)
{
  typedef typename Matrix<G>::index_type index_type;
  os << "% " << X.nrows() << " " << X.ncols();
  if (X.lower_bound()!=index_type(0,0)) os << X.lower_bound();
  os << endl;
  for (row_matrix_iterator<const Matrix<G> > it=row_it(X); it!=row_it(X)+X.nrows(); ++it)
    os << to_string(*it) << endl;
  return os;
}

//{unsecret}
template<class E,class Tr,class G> basic_istream<E,Tr> &operator>>(basic_istream<E,Tr> &is, Matrix<G> &X)
{
  typedef typename Matrix<G>::int_type int_type;
  typedef typename Matrix<G>::size_type size_type;

  char c; is >> c && c=='%';
	int_type m,n; is >> m >> n; X.resize(m,n);
  for (int_type i=0; i<X.nrows(); ++i)
    for (int_type j=0; j<X.ncols(); ++j)
      is >> X(i,j);
  return is;
}

template<class FwdIt,class V,class A> inline void uninitialized_fill_aux(FwdIt begin, FwdIt end, const Matrix<dense_matrix_generator<V,A> > &v) { uninitialized_fill(begin,end,v); }

//}

#endif




