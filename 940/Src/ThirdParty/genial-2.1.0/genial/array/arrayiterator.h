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

#ifndef ARRAYITERATOR_H
#define ARRAYITERATOR_H

//namespace genial
//{

using namespace std;

template<class A> struct array_traits;


template<class A> struct raster_iterator_rebind { };
template<class A> struct shift_raster_iterator_rebind { };

template<class A> struct lower_triangle_iterator_rebind { typedef typename raster_iterator_rebind<A>::other other; };

template<class It1,class It2> struct iterator_promotion_traits { typedef It1 iterator; };
template<class It1,class It2> struct iterator_promotion_traits<const It1,      It2> { typedef const typename iterator_promotion_traits<It1,It2>::iterator iterator; };
template<class It1,class It2> struct iterator_promotion_traits<      It1,const It2> { typedef const typename iterator_promotion_traits<It1,It2>::iterator iterator; };
template<class It1,class It2> struct iterator_promotion_traits<const It1,const It2> { typedef const typename iterator_promotion_traits<It1,It2>::iterator iterator; };



template<class A>
class array_iterator : public array_traits<A>
{
  public:
    typedef array_iterator self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES

    typedef bidirectional_iterator_tag iterator_category;

  protected:
    array_type &X;

  protected:
    inline explicit array_iterator(array_type &x) : X(x) {}

    inline array_iterator(const self &x) : X(x.X) {}
    template<class T2> inline array_iterator(const array_iterator<T2> &x) : X(x.array()) {}

    inline self &operator=(const self &r) { return *this; }
  
  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }

    inline size_type size() const { return X.size(); }

    inline index_type lower_bound() const { return X.lower_bound(); }
    inline index_type upper_bound() const { return X.upper_bound(); }
};


template<class A>
class index_array_iterator : public array_iterator<A>
{
  public:
    typedef index_array_iterator self;
    typedef array_iterator<A> base;
    ARRAY_BASE_TYPES

  protected:
    index_type ind;

  protected:
    inline explicit index_array_iterator(array_type &x) : base(x), ind() {}
    inline index_array_iterator(array_type &x, const index_type &i) : base(x), ind(i) {}

    inline index_array_iterator(const self &x) : base(x), ind(x.pos()) {}
    template<class T2> inline index_array_iterator(const index_array_iterator<T2> &x) : base(x), ind(x.pos()) {}

   inline  self &operator=(const self &x) { ind=x.ind; return *this; }
  
  public:
    using base::array;
    
    inline index_type       &pos()       { return ind; }
    inline const index_type &pos() const { return ind; }

    inline reference       operator*()       { return array()[ind]; }
    inline const_reference operator*() const { return array()[ind]; }

    inline bool operator==(const self &x) const { return ind==x.ind; }
    inline bool operator!=(const self &x) const { return ind!=x.ind; }
    
    inline bool operator<(const self &x) { return pos()<x.pos(); }     
};


template<class A> class shift_array_iterator;
template<class A> class shift_dense_array_iterator;


template<class It> struct shift_array_iterator_traits { /*typedef array_traits<It>::array_type array_type;*/ typedef It iterator; };
template<class A> struct shift_array_iterator_traits<      shift_array_iterator      <A> > { typedef       A array_type; typedef typename array_traits<array_type>::generator_type generator_type; typedef typename shift_array_iterator_traits<typename array_generator_traits<generator_type>::      iterator>::iterator iterator; };
template<class A> struct shift_array_iterator_traits<const shift_array_iterator      <A> > { typedef const A array_type; typedef typename array_traits<array_type>::generator_type generator_type; typedef typename shift_array_iterator_traits<typename array_generator_traits<generator_type>::const_iterator>::iterator iterator; };
template<class A> struct shift_array_iterator_traits<      shift_dense_array_iterator<A> > { typedef       A array_type; typedef typename array_traits<array_type>::generator_type generator_type; typedef typename shift_array_iterator_traits<typename array_generator_traits<generator_type>::      iterator>::iterator iterator; };
template<class A> struct shift_array_iterator_traits<const shift_dense_array_iterator<A> > { typedef const A array_type; typedef typename array_traits<array_type>::generator_type generator_type; typedef typename shift_array_iterator_traits<typename array_generator_traits<generator_type>::const_iterator>::iterator iterator; };


template<class A>
class shift_array_iterator : public array_iterator<A>
{
  public:
    typedef shift_array_iterator self;
    typedef array_iterator<A> base;
    ARRAY_BASE_TYPES
    
    typedef shift_array_iterator_traits<typename array_type::iterator> iter_traits;
    typedef typename iter_traits::iterator::template rebind<array_type>::other iter_type;

    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;

    template<class A2> struct rebind { typedef typename shift_iterator_rebind<A2>::other other; };

  protected:
    iter_type it;
 
  public:
    explicit shift_array_iterator(array_type &x) : base(x), it(x) {}
    shift_array_iterator(array_type &x, const iter_type  &i) : base(x), it(i) {}
    shift_array_iterator(array_type &x, const index_type &i) : base(x), it(x,i) {}
    shift_array_iterator(const self &x) : base(x), it(x.it) {}

    iter_type       &iter()       { return it; }
    const iter_type &iter() const { return it; }
          
    index_type pos() const { return it.pos()+this->lower_bound(); }
    
    reference       operator* ()       { return  (this->array())[pos()]; }
    const_reference operator* () const { return  (this->array())[pos()]; }
    pointer         operator->()       { return &(this->array())[pos()]; } 
    const_pointer   operator->() const { return &(this->array())[pos()]; }
    
    self &operator++() { ++it; return *this; }
    self &operator--() { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(this->array(),it+n); }
    self  operator- (const difference_type &n) const { return self(this->array(),it-n); }
    self &operator+=(const difference_type &n) { it+=n; return *this; }
    self &operator-=(const difference_type &n) { it-=n; return *this; }

    difference_type operator-(const self &x) const { return iter()-x.iter(); }

    bool operator==(const self &x) const { return it==x.it; }
    bool operator!=(const self &x) const { return it!=x.it; } 
};

template<class A ,class It> struct iterator_promotion_traits<shift_array_iterator<A>,It > { typedef shift_array_iterator<A> iterator; };
template<class It,class A > struct iterator_promotion_traits<It,shift_array_iterator<A> > { typedef shift_array_iterator<A> iterator; };
template<class A1,class A2> struct iterator_promotion_traits<shift_array_iterator<A1>,shift_array_iterator<A2> > { typedef shift_array_iterator<A1> iterator; };

template<class A> struct shift_iterator_rebind          { typedef shift_array_iterator<A> other; };



template<class A>
class shift_dense_array_iterator : public array_traits<A>
{
  public:
    typedef shift_dense_array_iterator self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES
    typedef typename base::generator_type generator_type;
    
    typedef typename array_traits<typename array_generator_traits<generator_type>::array_type>::iterator iter_type;
    
    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;
       
    template<class A2> struct rebind { typedef typename shift_iterator_rebind<A2>::other other; };
                     
  protected:
    iter_type it;
    
  public:
    explicit shift_dense_array_iterator(array_type &x) : it(x.generator().array()) {}
    explicit shift_dense_array_iterator(const iter_type &i) : it(i) {}
    shift_dense_array_iterator(const self &x) : it(x.it) {}
 
    self &operator=(const self &x) { it=x.it; return *this; }
    template<class B> self &operator=(const shift_dense_array_iterator<B> &x) { it=x.iter(); return *this; }
     
    iter_type       &iter()       { return it; }
    const iter_type &iter() const { return it; }
    
    reference       operator* ()       { return  *it; }
    const_reference operator* () const { return  *it; }
    pointer         operator->()       { return &*it; } 
    const_pointer   operator->() const { return &*it; }
    
    self &operator++() { ++it; return *this; }
    self &operator--() { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n) { it+=n; return *this; }
    self &operator-=(const difference_type &n) { it-=n; return *this; }

    difference_type operator-(const self &x) const { return iter()-x.iter(); }

    bool operator==(const self &x) const { return it==x.it; }
    bool operator!=(const self &x) const { return it!=x.it; } 
};

template<class A ,class It> struct iterator_promotion_traits<shift_dense_array_iterator<A>,It > { typedef shift_array_iterator<A> iterator; };
template<class It,class A > struct iterator_promotion_traits<It,shift_dense_array_iterator<A> > { typedef shift_array_iterator<A> iterator; };
template<class A1,class A2> struct iterator_promotion_traits<shift_dense_array_iterator<A1>,shift_dense_array_iterator<A2> > { typedef shift_array_iterator<A1> iterator; };

template<class A1,class A2> struct iterator_promotion_traits<shift_array_iterator      <A1>,shift_dense_array_iterator<A2> > { typedef shift_array_iterator<A1> iterator; };
template<class A1,class A2> struct iterator_promotion_traits<shift_dense_array_iterator<A1>,shift_array_iterator      <A2> > { typedef shift_array_iterator<A1> iterator; };

#ifdef NDEBUG
template<      class A,int C> struct shift_iterator_rebind<      Vector< shift_dense_array_generator<A,C> > > { typedef shift_dense_array_iterator<      Vector< shift_dense_array_generator<A,C> > > other; };
template<      class A,int C> struct shift_iterator_rebind<const Vector< shift_dense_array_generator<A,C> > > { typedef shift_dense_array_iterator<const Vector< shift_dense_array_generator<A,C> > > other; };
template<      class A,int C> struct shift_iterator_rebind<      Matrix< shift_dense_array_generator<A,C> > > { typedef shift_dense_array_iterator<      Matrix< shift_dense_array_generator<A,C> > > other; };
template<      class A,int C> struct shift_iterator_rebind<const Matrix< shift_dense_array_generator<A,C> > > { typedef shift_dense_array_iterator<const Matrix< shift_dense_array_generator<A,C> > > other; };
#endif


//}

#endif
