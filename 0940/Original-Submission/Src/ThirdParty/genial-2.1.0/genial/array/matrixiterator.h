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

#ifndef MATRIXITERATOR_H
#define MATRIXITERATOR_H


#include "array/arrayiterator.h"

//namespace genial
//{

using namespace std;

//Group = Matrices functions


template<class A>
class matrix_iterator : public array_iterator<A>
{
  public:
    typedef matrix_iterator self;
    typedef array_iterator<A> base;
    ARRAY_BASE_TYPES
        
  protected:
    explicit inline matrix_iterator(array_type &x) : base(x) {}
     
    using base::array;
    using base::size;    
 
    inline int_type nrows() const { return array().nrows(); }
    inline int_type ncols() const { return array().ncols(); }
    inline int_type nelms() const { return array().nelms(); }
};

template<class A>
class index_matrix_iterator : public index_array_iterator<A>
{
  public:
    typedef index_matrix_iterator self;
    typedef index_array_iterator<A> base;
    ARRAY_BASE_TYPES

  protected:
    explicit inline index_matrix_iterator(array_type &x) : base(x, index_type(0,0)) {}
    inline index_matrix_iterator(array_type &x, const index_type &p) : base(x,p) {}
    inline index_matrix_iterator(array_type &x, int_type i, int_type j) : base(x,index_type(i,j)) {}
    
    inline index_matrix_iterator(const self &x) : base(x) {}
    template<class T2> inline index_matrix_iterator(const index_matrix_iterator<T2> &x) : base(x) {}
    
  public:
    using base::array;
    using base::size;
    using base::pos;
  
    inline int_type nrows() const { return array().nrows(); }
    inline int_type ncols() const { return array().ncols(); }
    inline int_type nelms() const { return array().nelms(); }
  
    inline bool leftboundary () const { return array().leftboundary (pos()); }
    inline bool rightboundary() const { return array().rightboundary(pos()); }
  
    inline int_type row_lower_bound() { return array().row_lower_bound(); }
    inline int_type col_lower_bound() { return array().col_lower_bound(); }
    inline int_type row_upper_bound() { return array().row_upper_bound(); }
    inline int_type col_upper_bound() { return array().col_upper_bound(); }
    
    inline difference_type absolute_pos() const { return array().ncols()*pos().i+pos().j; }
    
    inline bool operator==(const self &x) const { return base::operator==(x); }
    inline bool operator!=(const self &x) const { return base::operator!=(x); }    
};


template<class A>
class raster_matrix_iterator : public index_matrix_iterator<A>
{
  public:
    typedef raster_matrix_iterator self;
    typedef index_matrix_iterator<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;    

    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  protected:
    using base::X;

  public:
    explicit inline raster_matrix_iterator(array_type &x) : base(x) {} // penser à changer raster_dense_matrix_iterator
    inline raster_matrix_iterator(array_type &x, const index_type &p) : base(x,p) {}
    inline raster_matrix_iterator(array_type &x, int_type i, int_type j) : base(x,i,j) {}
    inline raster_matrix_iterator(const self &x) : base(x) {}
    template<class T2> inline raster_matrix_iterator(const raster_matrix_iterator<T2> &x) : base(x) {}

    using base::array;
    using base::size;
    using base::nrows;
    using base::ncols;
    using base::nelms;
    using base::pos;
    using base::absolute_pos;

    inline self &operator++() { if (!size().rightboundary(pos())) ++pos().j; else { pos().j=0        ; ++pos().i; } return *this; }   
    inline self &operator--() { if (!size().rightboundary(pos())) --pos().j; else { pos().j=ncols()-1; --pos().i; } return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    
    inline self  operator+ (const difference_type &n) const { difference_type i=absolute_pos()+n; return self(X,i/ncols(),i%ncols()); }
    inline self  operator- (const difference_type &n) const { difference_type i=absolute_pos()-n; return self(X,i/ncols(),i%ncols()); }
    inline self &operator+=(const difference_type &n)       { difference_type i=absolute_pos()+n; pos()=index_type(i/ncols(),i%ncols()); return *this; }
    inline self &operator-=(const difference_type &n)       { difference_type i=absolute_pos()-n; pos()=index_type(i/ncols(),i%ncols()); return *this; }

    inline difference_type operator-(const self &x) const { return absolute_pos() - x.absolute_pos(); }
    
    inline bool operator==(const self &x) const { return base::operator==(x); }
    inline bool operator!=(const self &x) const { return base::operator!=(x); }    
};

template<class G> struct raster_iterator_rebind<      Matrix<G> > { typedef raster_matrix_iterator<      Matrix<G> > other; };
template<class G> struct raster_iterator_rebind<const Matrix<G> > { typedef raster_matrix_iterator<const Matrix<G> > other; };


template<class A>
class raster_dense_matrix_iterator : public array_traits<A>
{
  public:
    typedef raster_dense_matrix_iterator self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES    

    typedef random_access_iterator_tag iterator_category;    

    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  private:
    pointer ptr;

  public:
    explicit inline raster_dense_matrix_iterator(array_type &x) : ptr(&x(0,0)) {}
    inline raster_dense_matrix_iterator(array_type &x, const index_type &i) : ptr(&x[i]) {}
    inline raster_dense_matrix_iterator(const pointer &p) : ptr(p) {}

    inline raster_dense_matrix_iterator(const self &r) : ptr(r.ptr) {}
    template<class T2> inline raster_dense_matrix_iterator(const raster_dense_matrix_iterator<T2> &x) : ptr(&*x) {}

    inline self &operator=(const self &r) { ptr=r.ptr; return *this; }

    inline reference       operator*()       { return *ptr; }
    inline const_reference operator*() const { return *ptr; }

    inline self &operator++() { ++ptr; return *this; }
    inline self &operator--() { --ptr; return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    inline self  operator+ (const difference_type &n) const { return self(ptr+n); }
    inline self  operator- (const difference_type &n) const { return self(ptr-n); }
    inline self &operator+=(const difference_type &n) { ptr+=n; return *this; }
    inline self &operator-=(const difference_type &n) { ptr-=n; return *this; }    

    inline bool operator==(const self &r) const { return ptr==r.ptr; }
    inline bool operator!=(const self &r) const { return ptr!=r.ptr; }
    
    inline difference_type operator-(const self &x) const { return &**this-&*x; }
    
    inline bool operator<(const self &x) { return ptr<x.ptr; }        
};

#ifdef NDEBUG
template<class V> struct raster_iterator_rebind<      Matrix<dense_matrix_generator<V> > > { typedef raster_dense_matrix_iterator<      Matrix<dense_matrix_generator<V> > > other; };
template<class V> struct raster_iterator_rebind<const Matrix<dense_matrix_generator<V> > > { typedef raster_dense_matrix_iterator<const Matrix<dense_matrix_generator<V> > > other; };

//pas certain que ça marche

template<class V> struct raster_iterator_rebind<      Matrix<shift_array_generator<      Matrix<dense_matrix_generator<V> > > > > { typedef raster_dense_matrix_iterator<      Matrix<shift_array_generator<      Matrix<dense_matrix_generator<V> > > > > other; };
template<class V> struct raster_iterator_rebind<      Matrix<shift_array_generator<const Matrix<dense_matrix_generator<V> > > > > { typedef raster_dense_matrix_iterator<      Matrix<shift_array_generator<const Matrix<dense_matrix_generator<V> > > > > other; };
template<class V> struct raster_iterator_rebind<const Matrix<shift_array_generator<      Matrix<dense_matrix_generator<V> > > > > { typedef raster_dense_matrix_iterator<const Matrix<shift_array_generator<      Matrix<dense_matrix_generator<V> > > > > other; };
template<class V> struct raster_iterator_rebind<const Matrix<shift_array_generator<const Matrix<dense_matrix_generator<V> > > > > { typedef raster_dense_matrix_iterator<const Matrix<shift_array_generator<const Matrix<dense_matrix_generator<V> > > > > other; };
#endif


template<class A>
class shift_raster_matrix_iterator : public index_matrix_iterator<A>
{
  public:
    typedef shift_raster_matrix_iterator self;
    typedef index_matrix_iterator<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;    

    template<class A2> struct rebind { typedef typename raster_iterator_rebind<A2>::other other; };

  protected:
    using base::X;

  public:
    explicit inline shift_raster_matrix_iterator(array_type &x) : base(x,x.lower_bound()) {}
    inline shift_raster_matrix_iterator(array_type &x, const index_type &p) : base(x,p) {}
    inline shift_raster_matrix_iterator(array_type &x, int_type i, int_type j) : base(x,i,j) {}
    inline shift_raster_matrix_iterator(const self &x) : base(x) {}
    template<class T2> inline shift_raster_matrix_iterator(const shift_raster_matrix_iterator<T2> &x) : base(x) {}

    using base::array;
    using base::size;
    using base::nrows;
    using base::ncols;
    using base::nelms;
    using base::pos;
    using base::absolute_pos;
    using base::leftboundary;
    using base::rightboundary;
    using base::row_lower_bound;
    using base::col_lower_bound;
    using base::row_upper_bound;
    using base::col_upper_bound;

    inline self &operator++() { if (!rightboundary()) ++pos().j; else { pos().j=col_lower_bound(); ++pos().i; } return *this; }   
    inline self &operator--() { if (!leftboundary ()) --pos().j; else { pos().j=col_upper_bound(); --pos().i; } return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    
    inline self  operator+ (const difference_type &n) const { difference_type i=absolute_pos()+n; return self(X,i/ncols(),i%ncols()); }
    inline self  operator- (const difference_type &n) const { difference_type i=absolute_pos()-n; return self(X,i/ncols(),i%ncols()); }
    inline self &operator+=(const difference_type &n)       { difference_type i=absolute_pos()+n; pos()=index_type(i/ncols(),i%ncols()); return *this; }
    inline self &operator-=(const difference_type &n)       { difference_type i=absolute_pos()-n; pos()=index_type(i/ncols(),i%ncols()); return *this; }

    inline difference_type operator-(const self &x) const { return absolute_pos() - x.absolute_pos(); }
    
    inline bool operator==(const self &x) const { return base::operator==(x); }
    inline bool operator!=(const self &x) const { return base::operator!=(x); }    
};

template<class G> struct shift_raster_iterator_rebind<      Matrix<G> > { typedef shift_raster_matrix_iterator<      Matrix<G> > other; };
template<class G> struct shift_raster_iterator_rebind<const Matrix<G> > { typedef shift_raster_matrix_iterator<const Matrix<G> > other; };


template<class A>
class lower_triangle_matrix_iterator : public index_matrix_iterator<A>
{
  public:
    typedef lower_triangle_matrix_iterator self;
    typedef index_matrix_iterator<A> base;
    ARRAY_BASE_TYPES

    typedef random_access_iterator_tag iterator_category;    

  protected:
    using base::X;

  public:
    inline explicit lower_triangle_matrix_iterator(array_type &x) : base(x) {}
    inline lower_triangle_matrix_iterator(array_type &x, const index_type &p) : base(x,p) {}
    inline lower_triangle_matrix_iterator(array_type &x, int_type i, int_type j) : base(x,i,j) {}
    inline lower_triangle_matrix_iterator(const self &x) : base(x) {}
    template<class A2> inline lower_triangle_matrix_iterator(const lower_triangle_matrix_iterator<A2> &x) : base(x) {}
  
    using base::array;
    using base::size;
    using base::nrows;
    using base::ncols;
    using base::nelms;
    using base::pos;
    using base::leftboundary;
    using base::rightboundary;
    using base::row_lower_bound;
    using base::col_lower_bound;
    using base::row_upper_bound;
    using base::col_upper_bound;

    inline difference_type absolute_pos() const { return (pos().i*(pos().i+1))/2+pos().j; } 
    inline void set_absolute_pos(difference_type n) { pos().i=(sqrt(1+8*n)-1)/2; pos().j=n-(pos().i*(pos().i+1)/2); }  

    inline self &operator++() { if (pos().j<pos().i) ++pos().j; else { ++pos().i; pos().j=0;       } return *this; }   
    inline self &operator--() { if (pos().j>0      ) --pos().j; else { --pos().i; pos().j=pos().i; } return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    
    inline self  operator+ (const difference_type &n) const { return self(X,absolute_pos()+n); }
    inline self  operator- (const difference_type &n) const { return self(X,absolute_pos()-n); }
    inline self &operator+=(const difference_type &n)       { set_absolut_pos(absolute_pos()+n); return *this; }
    inline self &operator-=(const difference_type &n)       { set_absolut_pos(absolute_pos()-n); return *this; }

    inline difference_type operator-(const self &x) const { return absolute_pos() - x.absolute_pos(); }
 
    inline bool operator==(const self &x) const { return base::operator==(x); }
    inline bool operator!=(const self &x) const { return base::operator!=(x); }    
   
  private:
    inline lower_triangle_matrix_iterator(array_type &x, difference_type n) : base(x) { set_absolut_pos(n); }
};

template<class G> struct lower_triangle_iterator_rebind<      Matrix<G> > { typedef lower_triangle_matrix_iterator<      Matrix<G> > other; };
template<class G> struct lower_triangle_iterator_rebind<const Matrix<G> > { typedef lower_triangle_matrix_iterator<const Matrix<G> > other; };


template<class A>
class matrix_zigzag_iterator : public index_matrix_iterator<A>
{
  public:
    typedef matrix_zigzag_iterator self;
    typedef index_matrix_iterator<A> base;
    ARRAY_BASE_TYPES

  protected:
    using base::X;

  public:
    inline explicit matrix_zigzag_iterator(array_type &x) : base(x, index_type(0,0)) {}
    inline matrix_zigzag_iterator(array_type &x, const index_type &p) : base(x,p) {}
    inline matrix_zigzag_iterator(const self &r) : base(r) {}

    using base::array;
    using base::size;
    using base::nrows;
    using base::ncols;
    using base::nelms;
    using base::pos;
    using base::leftboundary;
    using base::rightboundary;
    using base::row_lower_bound;
    using base::col_lower_bound;
    using base::row_upper_bound;
    using base::col_upper_bound;

    inline self &operator++() { if ((pos().j+pos().i)&1) if (size().leftboundary (pos())) if (size().bottomboundary(pos())) ++pos().j; else ++pos().i; else if (size().bottomboundary(pos())) ++pos().j; else { ++pos().i; --pos().j; } else if (size().topboundary   (pos())) if (size().rightboundary(pos())) ++pos().i; else ++pos().j; else if (size().rightboundary(pos())) ++pos().i; else { ++pos().j; --pos().i; } return *this; } 
    inline self &operator--() { if ((pos().j+pos().i)&1) if (size().rightboundary(pos())) if (size().topboundary   (pos())) --pos().j; else --pos().i; else if (size().topboundary   (pos())) --pos().j; else { --pos().i; ++pos().j; } else if (size().bottomboundary(pos())) if (size().leftboundary (pos())) --pos().i; else --pos().j; else if (size().leftboundary (pos())) --pos().i; else { --pos().j; ++pos().i; } return *this; } 

    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    
    inline bool operator==(const self &x) const { return base::operator==(x); }
    inline bool operator!=(const self &x) const { return base::operator!=(x); }    
};

//{unsecret}
//Summary : Zigzag iterator on a matrix
//Arguments:
//  X - The matrix to iterate on
//  p - Start position
template<class G> matrix_zigzag_iterator<      Matrix<G> > zigzag_iterator(      Matrix<G> &X                                         ) { return matrix_zigzag_iterator<      Matrix<G> >(X  ); }
//{unsecret}
template<class G> matrix_zigzag_iterator<      Matrix<G> > zigzag_iterator(      Matrix<G> &X, const typename Matrix<G>::index_type &p) { return matrix_zigzag_iterator<      Matrix<G> >(X,p); }
//{unsecret}
template<class G> matrix_zigzag_iterator<const Matrix<G> > zigzag_iterator(const Matrix<G> &X                                         ) { return matrix_zigzag_iterator<const Matrix<G> >(X  ); }
//{unsecret}
template<class G> matrix_zigzag_iterator<const Matrix<G> > zigzag_iterator(const Matrix<G> &X, const typename Matrix<G>::index_type &p) { return matrix_zigzag_iterator<const Matrix<G> >(X,p); }


template<class A>
class diamond_matrix_iterator : public index_matrix_iterator<A>
{
  public:
    typedef diamond_matrix_iterator self;
    typedef index_matrix_iterator<A> base;
    ARRAY_BASE_TYPES

  public:
    inline explicit diamond_matrix_iterator(array_type &x) : base(x, index_type(0,0)) {}
    inline diamond_matrix_iterator(array_type &x, const index_type &i) : base(x,i) {}
    inline diamond_matrix_iterator(const self &r) : base(r) {}

    using base::pos;

    inline self &operator++() { if (pos().x> 0 && pos().y< 0) { ++pos().x; ++pos().y; } else if (pos().x> 0 && pos().y>=0) { --pos().x; ++pos().y; } else if (pos().x<=0 && pos().y> 0) { --pos().x; --pos().y; } else if (pos().x< 0 && pos().y<=0) { ++pos().x; --pos().y; } else pos().x+=1; return *this; }
    inline self &operator--() { if (pos().x> 1 && pos().y<=0) { --pos().x; --pos().y; } else if (pos().x>=0 && pos().y> 0) { ++pos().x; --pos().y; } else if (pos().x< 0 && pos().y>=0) { ++pos().x; ++pos().y; } else if (pos().x<=0 && pos().y< 0) { --pos().x; ++pos().y; } else pos().x-=1; return *this; }

    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    
    inline bool operator==(const self &x) const { return base::operator==(x); }
    inline bool operator!=(const self &x) const { return base::operator!=(x); }    
};

template<class G> diamond_matrix_iterator<      Matrix<G> > diamond_iterator(      Matrix<G> &X                                         ) { return diamond_matrix_iterator<      Matrix<G> >(X  ); }
template<class G> diamond_matrix_iterator<      Matrix<G> > diamond_iterator(      Matrix<G> &X, const typename Matrix<G>::index_type &i) { return diamond_matrix_iterator<      Matrix<G> >(X,i); }
template<class G> diamond_matrix_iterator<const Matrix<G> > diamond_iterator(const Matrix<G> &X                                         ) { return diamond_matrix_iterator<const Matrix<G> >(X  ); }
template<class G> diamond_matrix_iterator<const Matrix<G> > diamond_iterator(const Matrix<G> &X, const typename Matrix<G>::index_type &i) { return diamond_matrix_iterator<const Matrix<G> >(X,i); }


template<class A>
class row_matrix_iterator : public matrix_iterator<A>
{
  public:
    typedef row_matrix_iterator self;
    typedef matrix_iterator<A> base;
    typedef typename base::array_type array_type;
    typedef typename base::int_type index_type;
    typedef typename base::int_type int_type;
    typedef typename base::difference_type difference_type;  

    typedef random_access_iterator_tag iterator_category;  
    
  private:
    using base::X;
    index_type ind;

  public:
    inline explicit row_matrix_iterator(array_type &x) : base(x), ind(x.row_lower_bound()) {}
    inline row_matrix_iterator(array_type &x, index_type n) : base(x), ind(n) {}
    inline row_matrix_iterator(const self &x) : base(x), ind(x.ind) {}
    template<class T2> inline row_matrix_iterator(const row_matrix_iterator<T2> &x) : base(x), ind(x.pos()) {}

    inline self &operator=(const self &x) { ind=x.ind; return *this; }

    using base::array;
    using base::nrows;
    using base::ncols;

    inline index_type        &pos  ()       { return ind; }
    inline const index_type  &pos  () const { return ind; }
    
    inline int_type     nelms      () const { return nrows(); }
    
    inline index_type   lower_bound() const { return array().row_lower_bound(); }
    inline index_type   upper_bound() const { return array().row_upper_bound(); }

    inline typename MatrixRow<      array_type>::self operator*()       { return row(array(),ind); }
    inline typename MatrixRow<const array_type>::self operator*() const { return row(array(),ind); }

    inline self &operator++() { ++ind; return *this; }   
    inline self &operator--() { --ind; return *this; }
    inline self  operator++(int) { self t=*this; ++*this; return t; }
    inline self  operator--(int) { self t=*this; --*this; return t; }
    
    inline self  operator+ (const difference_type &n) const { return self(X,pos()+n); }
    inline self  operator- (const difference_type &n) const { return self(X,pos()-n); }
    inline self &operator+=(const difference_type &n)       { ind+=n; return *this; }
    inline self &operator-=(const difference_type &n)       { ind-=n; return *this; }

    inline difference_type operator-(const self &x) const { return ind-x.ind; }

    inline bool operator==(const self &x) const { return ind==x.ind; }
    inline bool operator!=(const self &x) const { return ind!=x.ind; }
};

//{unsecret}
//Summary: Iterator on each row
//Arguments:
//  X - The matrix to iterate on
//  n - Start row position
template<class G> row_matrix_iterator<      Matrix<G> > row_it(      Matrix<G> &X                                       ) { return row_matrix_iterator<      Matrix<G> >(X  ); }
//{unsecret}
template<class G> row_matrix_iterator<      Matrix<G> > row_it(      Matrix<G> &X, const typename Matrix<G>::int_type &n) { return row_matrix_iterator<      Matrix<G> >(X,n); }
//{unsecret}
template<class G> row_matrix_iterator<const Matrix<G> > row_it(const Matrix<G> &X                                       ) { return row_matrix_iterator<const Matrix<G> >(X  ); }
//{unsecret}
template<class G> row_matrix_iterator<const Matrix<G> > row_it(const Matrix<G> &X, const typename Matrix<G>::int_type &n) { return row_matrix_iterator<const Matrix<G> >(X,n); }


//}


#endif

