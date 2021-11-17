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

#ifndef ITERATOR_ADAPTOR_H
#define ITERATOR_ADAPTOR_H

#include <iterator>
#include "functional.h"

//namespace genial
//{

using namespace std;


template <class It, class Dist> inline It &inc(It &it, Dist n) { advance(it, n); return it; }
template <class It, class Dist> inline It &inc(It &it        ) { return inc(it, 1); }
template <class It, class Dist> inline It &dec(It &it        ) { return inc(it,-1); }
template <class It, class Dist> inline It  adv(It &it        ) { return adv(it, 1); }

template <class It, class Dist> inline It  adv(const It &it, Dist n) { It it2=it; return inc(it2,n); }


template<class IS,class V=typename IS::char_type>
class is_iterator
{
  public:
    typedef is_iterator self;
    
    typedef IS istream_type;
    
    typedef input_iterator_tag iterator_category;
    
    typedef V                 value_type;
    typedef const value_type &reference;
    typedef const value_type &const_reference;
    typedef const value_type *pointer;
    typedef const value_type *const_pointer;
    typedef ptrdiff_t difference_type;
    //typedef int difference_type;
  
  private:
    istream_type *p;
    value_type val;
  
  public:
    is_iterator() : p(NULL) {}
    explicit is_iterator(istream_type &is) : p(&is) { ++*this; }
                
    self &operator++() 
    { 
      *p>>val; 
      if (p->fail()) 
        p=NULL; 
      return *this; 
    }
    
    const_reference operator*() const { return val; }
    
    bool operator==(const self &x) const { return p==x.p; }
    bool operator!=(const self &x) const { return !(*this==x); }
};


template<class It, class Pred>
class filter_iterator
{
  public:
    typedef filter_iterator self;

    typedef It iter_type;
    typedef Pred predicate_type;

    typedef forward_iterator_tag iterator_category;

    typedef typename iterator_traits<iter_type>::value_type      value_type;
    typedef typename iterator_traits<iter_type>::reference       reference;
    typedef typename iterator_traits<iter_type>::const_reference const_reference;
    typedef typename iterator_traits<iter_type>::pointer         pointer;
    typedef typename iterator_traits<iter_type>::const_pointer   const_pointer;
    typedef typename iterator_traits<iter_type>::difference_type difference_type;

  protected:
    iter_type it;
    iter_type end;
    predicate_type pred;

  public:
    filter_iterator() {};
    filter_iterator(iter_type b                                       ) : it(b), end( ), pred(  ) { if (it!=end && !pred(*it)) ++(*this); }
    filter_iterator(iter_type b, iter_type e                          ) : it(b), end(e), pred(  ) { if (it!=end && !pred(*it)) ++(*this); }
    filter_iterator(iter_type b,              const predicate_type &pr) : it(b), end( ), pred(pr) { if (it!=end && !pred(*it)) ++(*this); }
    filter_iterator(iter_type b, iter_type e, const predicate_type &pr) : it(b), end(e), pred(pr) { if (it!=end && !pred(*it)) ++(*this); }

    iter_type       &base()       { return it; }
    const iter_type &base() const { return it; }

    reference       operator* ()       { return  *it; }
    const_reference operator* () const { return  *it; }
    pointer         operator->()       { return &*it; } 
    const_pointer   operator->() const { return &*it; }

    self &operator++()    { do ++it; while (it!=end && !pred(*it)); return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }

    bool operator==(const self &x) const { return it==x.it; }
    bool operator!=(const self &x) const { return !(*this==x); } 
};

template <class It,class Pred> inline filter_iterator<It,Pred> filter_it(It b,       const Pred& pr) { return filter_iterator<It,Pred>(b,  pr); }
template <class It,class Pred> inline filter_iterator<It,Pred> filter_it(It b, It e, const Pred& pr) { return filter_iterator<It,Pred>(b,e,pr); }


template<class It>
class cyclic_iterator
{ 
  public:
    typedef cyclic_iterator self;

    typedef It iter_type;

    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;

    typedef typename iterator_traits<iter_type>::value_type      value_type;
    typedef typename iterator_traits<iter_type>::reference       reference;
    typedef const value_type                                    &const_reference;
    typedef typename iterator_traits<iter_type>::pointer         pointer;
    typedef const value_type                                    *const_pointer;
    typedef typename iterator_traits<iter_type>::difference_type difference_type;

  protected:
    iter_type it;
    iter_type begin;
    iter_type end;

  public:
    cyclic_iterator(iter_type b, iter_type e) : it(b), begin(b), end(e) {};
    cyclic_iterator(iter_type b, iter_type e, iter_type &i) : it(i), begin(b), end(e) {};

    iter_type       &base()       { return it; }
    const iter_type &base() const { return it; }

    reference       operator* ()       { return  *it; }
    const_reference operator* () const { return  *it; }
    pointer         operator->()       { return &*it; } 
    const_pointer   operator->() const { return &*it; }

    self &operator++()    { ++it; if (it==end) it=begin; return *this; }
    self &operator--()    { if (it==begin) it=end; --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }

    //manque random access:modulo?

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; } 
};

template<class T > inline cyclic_iterator<typename T::      iterator> cyclic_it(      T &x                        ) { return cyclic_iterator<typename T::      iterator>(x.begin(),x.end()); }
template<class T > inline cyclic_iterator<typename T::      iterator> cyclic_it(      T &x, typename T::iterator i) { return cyclic_iterator<typename T::      iterator>(x.begin(),x.end(),i); }
template<class It> inline cyclic_iterator<It                        > cyclic_it(It b, It e                        ) { return cyclic_iterator<It                        >(b,e); }



template<class It>
class adjacency_iterator
{ 
  public:
    typedef adjacency_iterator self;

    typedef It iter_type;

    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;

    typedef typename iterator_traits<iter_type>::value_type      iter_value_type;

    typedef pair<iter_value_type,iter_value_type> value_type;
    typedef value_type        reference;
    typedef const value_type  const_reference;
    typedef value_type       *pointer;
    typedef const value_type *const_pointer;

    typedef typename iterator_traits<iter_type>::difference_type difference_type;

  protected:
    iter_type it;

  public:
    explicit adjacency_iterator(iter_type a) : it(a) {};

    iter_type       &base()       { return it; }
    const iter_type &base() const { return it; }

    const_reference operator*() const { return value_type(*it,*++iter_type(it)); }

    self &operator++()    { ++it; return *this; }
    self &operator--()    { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n)       { it+=n; return *this; }
    self &operator-=(const difference_type &n)       { it-=n; return *this; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; } 
};

template<class It> inline adjacency_iterator<It> adjacency_it(It it) { return adjacency_iterator<It>(it); }



template<class It>
class iterator_iterator
{
  public:
    typedef iterator_iterator self;

    typedef It iter_type;

    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;

    typedef typename iterator_traits<iter_type>::value_type iter_value_type;
    typedef typename iterator_traits<iter_value_type>::value_type value_type;
    typedef typename iterator_traits<iter_value_type>::reference reference;
    typedef const value_type &const_reference;
    typedef typename iterator_traits<iter_value_type>::pointer pointer;
    typedef const value_type *const_pointer;

    typedef typename iterator_traits<iter_type>::difference_type difference_type;

  protected:
    iter_type it;

  public:
    iterator_iterator() : it() {}
    explicit iterator_iterator(const iter_type &i) : it(i) {}
    iterator_iterator(const self &x) : it(x.it) {}

    self &operator=(const self &x) { it=x.it; return *this; }
    template<class A> self &operator=(const iterator_iterator<A> &x) { it=x.base(); return *this; }

    iter_type       &base()       { return it; }
    const iter_type &base() const { return it; }

    reference       operator* ()       { return  **it; }
    const_reference operator* () const { return  **it; }
    pointer         operator->()       { return &**it; } 
    const_pointer   operator->() const { return &**it; }

    self &operator++() { ++it; return *this; }
    self &operator--() { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n) { it+=n; return *this; }
    self &operator-=(const difference_type &n) { it-=n; return *this; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; } 
    
    bool operator!() const { return !it; } 
};



template<class It>
class buffer_iterator
{
  public:
    typedef buffer_iterator self;

    typedef It iter_type;

    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;

    typedef typename iterator_traits<iter_type>::value_type value_type;
    typedef value_type reference;
    typedef const value_type const_reference;
    typedef typename iterator_traits<iter_type>::pointer pointer;
    typedef const value_type *const_pointer;

    typedef typename iterator_traits<iter_type>::difference_type difference_type;

  protected:
    iter_type it;
    mutable bool b;
    mutable value_type val;

  public:
    buffer_iterator() : it(), b(true), val() {}
    explicit buffer_iterator(const iter_type &i) : it(i), b(true), val() {}
    buffer_iterator(const self &r) : it(r.it), b(r.b), val(r.val) {}

    self &operator=(const self &r) { it=r.it; b=r.b; val=r.val; return *this; }
    template<class A> self &operator=(const buffer_iterator<A> &r) { it=r.base(); b=r.empty(); val=r.value(); return *this; }

    iter_type       &base()       { return it; }
    const iter_type &base() const { return it; }

    bool            empty() const { return b; }

    reference       value()       { return val; }
    const_reference value() const { return val; }

    typename iter_type::index_type       &pos()       { return base().pos(); }
    const typename iter_type::index_type &pos() const { return base().pos(); }

//    reference       operator* ()       { return  value(); }

    const_reference operator* () const { if (empty()) update(); return value(); }
//    pointer         operator->()       { return &value(); } 
//    const_pointer   operator->() const { return &value(); }

    self &operator++() { ++it; b=true; return *this; }
    self &operator--() { --it; b=true; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n) { it+=n; b=true; return *this; }
    self &operator-=(const difference_type &n) { it-=n; b=true; return *this; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; } 
    
    bool operator!() const { return !it; } 

  private:
    void update() const { b=false; val=*it; }
};

template<class It> inline buffer_iterator<It> buffer_it(const It &it) { return buffer_iterator<It>(it); }



template<class It, class F>
class function_iterator
{
  public:
    typedef function_iterator self;

    typedef It iter_type;
    typedef F function_type;

    typedef typename iterator_traits<iter_type>::iterator_category iterator_category;
    typedef typename iterator_traits<iter_type>::difference_type difference_type;

    typedef typename function_type::argument_type argument_type;
    typedef typename function_type::value_type value_type;
    typedef typename function_type::reference reference;
    typedef typename function_type::const_reference const_reference;
    typedef typename function_type::pointer pointer;
    typedef typename function_type::const_pointer const_pointer;

  protected:
    iter_type it;
    mutable function_type func;

  public:
    function_iterator() : it(), func() {}
    explicit function_iterator(const iter_type &i) : it(i), func() {}
    template<class A> explicit function_iterator(const A &a) : it(a), func() {}
    function_iterator(const self &r) : it(r.it), func(r.func) {}

    self &operator=(const self &r) { it=r.it; func=r.func; return *this; }

    iter_type       &base()       { return it; }
    const iter_type &base() const { return it; }

    reference       operator *()       { typename iterator_traits<iter_type>::reference       ref=*it; return  func(ref); }
    const_reference operator *() const { typename iterator_traits<iter_type>::const_reference ref=*it; return  func(ref); }
    pointer         operator->()       { typename iterator_traits<iter_type>::reference       ref=*it; return &func(ref); }
    const_pointer   operator->() const { typename iterator_traits<iter_type>::const_reference ref=*it; return &func(ref); }

    self &operator++() { ++it; return *this; }
    self &operator--() { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n) { it+=n; return *this; }
    self &operator-=(const difference_type &n) { it-=n; return *this; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; }  
};


template<class Arg> struct name_function;
template<class Arg> struct to_string_function;
template<class It> inline function_iterator<It, name_function     <const typename It::value_type> > name_it     (const It &i) { return function_iterator<It, name_function     <const typename It::value_type> >(i); }
template<class It> inline function_iterator<It, to_string_function<const typename It::value_type> > to_string_it(const It &i) { return function_iterator<It, to_string_function<const typename It::value_type> >(i); }


//}

#endif
