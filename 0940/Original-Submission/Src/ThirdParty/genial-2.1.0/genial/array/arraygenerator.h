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

#ifndef ARRAYGENERATOR_H
#define ARRAYGENERATOR_H

//namespace genial
//{

using namespace std;

template<class G> class Vector;
template<class G> class Matrix;
template<int D, class G> class Array;

template<class V,class A=alignment_allocator<V> > struct DenseVector;
template<class V,class A=alignment_allocator<V> > struct DenseMatrix;

template<int N,      class V> struct TinyVector;
template<int M,int N,class V> struct TinyMatrix;

template<int N,class V> struct SimdVector;
template<int N,class V> struct SaturatedSimdTinyVector;

template<class I> class MatrixSize;
template<class I> class MatrixIndex;

struct vector_array_tag {};
struct matrix_array_tag {};

template<class V,class A=alignment_allocator<V> > class dense_vector_generator;
template<class V,class A=alignment_allocator<V> > class dense_matrix_generator;
template<      int N,class V> class tiny_vector_generator;
template<int M,int N,class V> class tiny_matrix_generator;
template<class V> class data_vector_generator;
template<class V> class data_matrix_generator;
template<      class T,int C> class row_dense_matrix_generator;
template<int N,class T,int C> class tiny_row_dense_matrix_generator;

template<int N, class V> struct TinyVector;
template<class V> struct DataVector;
template<class V> struct DataMatrix;


template<class G, class C=typename G::array_category> struct GeneratorArray { };
//#if defined(__ICL) 
//template<class G> struct GeneratorArray<G, vector_array_tag> { typedef G generator_type; typedef Array<1,generator_type> self; };
//template<class G> struct GeneratorArray<G, matrix_array_tag> { typedef G generator_type; typedef Array<2,generator_type> self; };
//#else
template<class G> struct GeneratorArray<G, vector_array_tag> { typedef G generator_type; typedef Vector<generator_type> self; };
template<class G> struct GeneratorArray<G, matrix_array_tag> { typedef G generator_type; typedef Matrix<generator_type> self; };
//#endif

#ifdef HIDE_FROM_DOCJET
#else
#define ARRAY_BASE_TYPES \
  typedef typename base::array_type      array_type;      \
  typedef typename base::value_type      value_type;      \
  typedef typename base::reference       reference;       \
  typedef typename base::const_reference const_reference; \
  typedef typename base::pointer         pointer;         \
  typedef typename base::const_pointer   const_pointer;   \
  typedef typename base::index_type      index_type;      \
  typedef typename base::size_type       size_type;       \
  typedef typename base::difference_type difference_type; \
  typedef typename base::int_type        int_type;

#define TWO_ARRAYS_BASE_TYPES \
  typedef typename base::first_array_type  first_array_type;  \
  typedef typename base::second_array_type second_array_type; \
  typedef typename base::value_type      value_type;          \
  typedef typename base::reference       reference;           \
  typedef typename base::const_reference const_reference;     \
  typedef typename base::pointer         pointer;             \
  typedef typename base::const_pointer   const_pointer;       \
  typedef typename base::index_type      index_type;          \
  typedef typename base::size_type       size_type;           \
  typedef typename base::difference_type difference_type;     \
  typedef typename base::int_type        int_type;

#define THREE_ARRAYS_BASE_TYPES \
  typedef typename base::first_array_type  first_array_type;  \
  typedef typename base::second_array_type second_array_type; \
  typedef typename base::third_array_type  third_array_type;  \
  typedef typename base::value_type      value_type;          \
  typedef typename base::reference       reference;           \
  typedef typename base::const_reference const_reference;     \
  typedef typename base::pointer         pointer;             \
  typedef typename base::const_pointer   const_pointer;       \
  typedef typename base::index_type      index_type;          \
  typedef typename base::size_type       size_type;           \
  typedef typename base::difference_type difference_type;     \
  typedef typename base::int_type        int_type;

#endif


//{unsecret}
//{group:Arrays}
template<class G>
struct generator_traits
{
  typedef G generator_type;
  typedef typename generator_type::array_category array_category;

  typedef typename generator_type::value_type value_type; 
  typedef typename generator_type::const_value_type const_value_type;
  typedef typename generator_type::reference reference;
  typedef typename generator_type::const_reference const_reference;
  typedef typename generator_type::pointer pointer;
  typedef typename generator_type::const_pointer const_pointer;   
  
  typedef typename generator_type::int_type int_type;
  typedef typename generator_type::difference_type difference_type;
  typedef typename generator_type::index_type index_type;
  typedef typename generator_type::size_type size_type;  

  template<class A2> struct array_rebind          { typedef typename generator_type::template array_rebind      <A2>::other other; };
    
  template<class A2> struct iterator_rebind       { typedef typename generator_type::template iterator_rebind      <A2>::other other; };
  template<class A2> struct const_iterator_rebind { typedef typename generator_type::template const_iterator_rebind<A2>::other other; };
};

template<class G>
struct generator_traits<const G>
{
  typedef const G generator_type;
  typedef typename generator_type::array_category array_category;

  typedef typename generator_type::value_type value_type;
  typedef typename generator_type::const_value_type const_value_type;
  typedef typename generator_type::const_reference reference;
  typedef typename generator_type::const_reference const_reference;
  typedef typename generator_type::const_pointer pointer;
  typedef typename generator_type::const_pointer const_pointer;  

  typedef typename generator_type::int_type int_type;
  typedef typename generator_type::difference_type difference_type;
  typedef typename generator_type::index_type index_type;
  typedef typename generator_type::size_type size_type;     
  
  template<class A2> struct array_rebind          { typedef typename generator_type::template array_rebind      <A2>::other other; };
  
  template<class A2> struct iterator_rebind       { typedef typename generator_type::template const_iterator_rebind<A2>::other other; };
  template<class A2> struct const_iterator_rebind { typedef typename generator_type::template const_iterator_rebind<A2>::other other; };
};

template<class G>
struct array_generator_traits : public generator_traits<G>
{
  typedef generator_traits<G> base;
  typedef typename base::generator_type generator_type;

  typedef typename generator_type::array_type array_type;
  
  typedef typename generator_type::iterator       iterator;
  typedef typename generator_type::const_iterator const_iterator;
};

template<class G>
struct array_generator_traits<const G> : public generator_traits<const G>
{
  typedef generator_traits<G> base;
  typedef typename base::generator_type generator_type;

  typedef const typename generator_type::array_type array_type;
  
  typedef typename generator_type::const_iterator iterator;
  typedef typename generator_type::const_iterator const_iterator;
};


template<class A>
struct array_traits
{
  typedef A array_type;
  typedef typename array_type::array_category array_category;
  
  typedef typename array_type::generator_type generator_type;

  typedef typename array_type::value_type value_type; 
  typedef typename array_type::const_value_type const_value_type;
  typedef typename array_type::reference reference;
  typedef typename array_type::const_reference const_reference;
  typedef typename array_type::pointer pointer;
  typedef typename array_type::const_pointer const_pointer;   
  
  typedef typename array_type::int_type int_type;
  typedef typename array_type::difference_type difference_type;
  typedef typename array_type::index_type index_type;
  typedef typename array_type::size_type size_type;  
  
  typedef typename array_type::iterator       iterator;
  typedef typename array_type::const_iterator const_iterator;  
  typedef typename array_type::reverse_iterator       reverse_iterator;
  typedef typename array_type::const_reverse_iterator const_reverse_iterator;  
  
  template<class V2> struct rebind          { typedef typename array_type::template rebind<V2>::other other; };
  template<class V2> struct array_rebind    { typedef typename array_type::template rebind<V2>::other other; };
  
  template<class A2> struct iterator_rebind       { typedef typename array_type::template iterator_rebind      <A2>::other other; };
  template<class A2> struct const_iterator_rebind { typedef typename array_type::template const_iterator_rebind<A2>::other other; };
};

template<class A>
struct array_traits<const A>
{
  typedef const A array_type;
  typedef typename array_type::array_category array_category;
  
  typedef const typename array_type::generator_type generator_type;  

  typedef typename array_type::value_type value_type;
  typedef typename array_type::const_value_type const_value_type;
  typedef typename array_type::const_reference reference;
  typedef typename array_type::const_reference const_reference;
  typedef typename array_type::const_pointer pointer;
  typedef typename array_type::const_pointer const_pointer;  

  typedef typename array_type::int_type int_type;
  typedef typename array_type::difference_type difference_type;
  typedef typename array_type::index_type index_type;
  typedef typename array_type::size_type size_type;     

  typedef typename array_type::const_iterator iterator;
  typedef typename array_type::const_iterator const_iterator;  
  typedef typename array_type::const_reverse_iterator reverse_iterator;
  typedef typename array_type::const_reverse_iterator const_reverse_iterator;  
  
  template<class V2> struct rebind          { typedef typename array_type::template rebind<V2>::other other; };
  template<class V2> struct array_rebind    { typedef typename array_type::template rebind<V2>::other other; };
  
  template<class A2> struct iterator_rebind       { typedef typename array_type::template const_iterator_rebind<A2>::other other; };
  template<class A2> struct const_iterator_rebind { typedef typename array_type::template const_iterator_rebind<A2>::other other; };
};

template<class It1,class It2> struct iterator_promotion_traits;

template<class A1,class A2=A1>
struct two_arrays_traits: public array_traits<A1>
{
  typedef two_arrays_traits self;

  typedef A1 first_array_type;
  typedef A2 second_array_type;
 
  typedef typename first_array_type::array_category array_category;

  typedef PROMOTE2(typename array_traits<first_array_type>::value_type,typename array_traits<second_array_type>::value_type) value_type;

  typedef typename array_traits<first_array_type >::iterator       first_iterator;
  typedef typename array_traits<first_array_type >::const_iterator first_const_iterator;
  typedef typename array_traits<second_array_type>::iterator       second_iterator;
  typedef typename array_traits<second_array_type>::const_iterator second_const_iterator;
  
  typedef typename array_traits<first_array_type >::reverse_iterator       first_reverse_iterator;
  typedef typename array_traits<first_array_type >::const_reverse_iterator first_const_reverse_iterator;
  typedef typename array_traits<second_array_type>::reverse_iterator       second_reverse_iterator;
  typedef typename array_traits<second_array_type>::const_reverse_iterator second_const_reverse_iterator;
  
  template<class A> struct first_iterator_rebind        { typedef typename array_traits<first_array_type >::template       iterator_rebind<A>::other other; };
  template<class A> struct first_const_iterator_rebind  { typedef typename array_traits<first_array_type >::template const_iterator_rebind<A>::other other; };
  template<class A> struct second_iterator_rebind       { typedef typename array_traits<second_array_type>::template       iterator_rebind<A>::other other; };
  template<class A> struct second_const_iterator_rebind { typedef typename array_traits<second_array_type>::template const_iterator_rebind<A>::other other; };
  
  typedef typename iterator_promotion_traits<first_iterator      , second_iterator      >::iterator iterator;
  typedef typename iterator_promotion_traits<first_const_iterator, second_const_iterator>::iterator const_iterator;
  
  template<class B> struct iterator_rebind       { typedef typename iterator      ::template rebind<      B>::other other; };
  template<class B> struct const_iterator_rebind { typedef typename const_iterator::template rebind<const B>::other other; };

  template<class V> struct first_array_rebind  { typedef typename array_traits<first_array_type >::template rebind<V>::other other; };
  template<class V> struct second_array_rebind { typedef typename array_traits<second_array_type>::template rebind<V>::other other; };
    
  template<class V> struct array_rebind { typedef PROMOTE2(typename first_array_rebind<V>::other,typename second_array_rebind<V>::other) other; };
};

template<class A1,class A2=A1,class A3=A2>
struct three_arrays_traits: public array_traits<A1>
{
  typedef three_arrays_traits self;

  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef A3 third_array_type;
 
  typedef typename first_array_type::array_category array_category;

  typedef PROMOTE3(typename array_traits<first_array_type>::value_type,typename array_traits<second_array_type>::value_type,typename array_traits<third_array_type>::value_type) value_type;

  typedef typename array_traits<first_array_type >::iterator       first_iterator;
  typedef typename array_traits<first_array_type >::const_iterator first_const_iterator;
  typedef typename array_traits<second_array_type>::iterator       second_iterator;
  typedef typename array_traits<second_array_type>::const_iterator second_const_iterator;
  typedef typename array_traits<third_array_type >::iterator       third_iterator;
  typedef typename array_traits<third_array_type >::const_iterator third_const_iterator;
  
  typedef typename array_traits<first_array_type >::reverse_iterator       first_reverse_iterator;
  typedef typename array_traits<first_array_type >::const_reverse_iterator first_const_reverse_iterator;
  typedef typename array_traits<second_array_type>::reverse_iterator       second_reverse_iterator;
  typedef typename array_traits<second_array_type>::const_reverse_iterator second_const_reverse_iterator;
  typedef typename array_traits<third_array_type >::reverse_iterator       third_reverse_iterator;
  typedef typename array_traits<third_array_type >::const_reverse_iterator third_const_reverse_iterator;

  template<class A> struct first_iterator_rebind        { typedef typename array_traits<first_array_type >::template       iterator_rebind<A>::other other; };
  template<class A> struct first_const_iterator_rebind  { typedef typename array_traits<first_array_type >::template const_iterator_rebind<A>::other other; };
  template<class A> struct second_iterator_rebind       { typedef typename array_traits<second_array_type>::template       iterator_rebind<A>::other other; };
  template<class A> struct second_const_iterator_rebind { typedef typename array_traits<second_array_type>::template const_iterator_rebind<A>::other other; };
  template<class A> struct third_iterator_rebind        { typedef typename array_traits<third_array_type >::template       iterator_rebind<A>::other other; };
  template<class A> struct third_const_iterator_rebind  { typedef typename array_traits<third_array_type >::template const_iterator_rebind<A>::other other; };
  
  typedef typename iterator_promotion_traits<first_iterator      , typename iterator_promotion_traits<second_iterator      , third_iterator      >::iterator>::iterator iterator;
  typedef typename iterator_promotion_traits<first_const_iterator, typename iterator_promotion_traits<second_const_iterator, third_const_iterator>::iterator>::iterator const_iterator;
  
  template<class B> struct iterator_rebind       { typedef typename iterator      ::template rebind<      B>::other other; };
  template<class B> struct const_iterator_rebind { typedef typename const_iterator::template rebind<const B>::other other; };
  
  template<class V> struct first_array_rebind  { typedef typename array_traits<first_array_type >::template rebind<V>::other other; };
  template<class V> struct second_array_rebind { typedef typename array_traits<second_array_type>::template rebind<V>::other other; };
  template<class V> struct third_array_rebind  { typedef typename array_traits<third_array_type >::template rebind<V>::other other; };
  
  template<class V> struct array_rebind { typedef PROMOTE3(typename first_array_rebind<V>::other,typename second_array_rebind<V>::other,typename third_array_rebind<V>::other) other; };
};

template<class T, class V=typename T::value_type>
struct array_value_traits : public array_traits<T>
{
  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef value_type reference;
  typedef const value_type const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class A1,class A2=A1,class V=typename two_arrays_traits<A1,A2>::value_type>
struct two_arrays_value_traits : public two_arrays_traits<A1,A2>
{
  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef value_type reference;
  typedef const value_type const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class T, class V=typename T::value_type, class Ref=V &, class ConstRef=const V &>
struct array_reference_traits : public array_traits<T>
{
  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class T, class V, class Ref, class ConstRef>
struct array_reference_traits<const T,V,Ref,ConstRef> : public array_traits<const T>
{
  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef ConstRef reference;
  typedef ConstRef const_reference;
  typedef const value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class A1,class A2=A1,class V=typename two_arrays_traits<A1,A2>::value_type,class Ref=V &,class ConstRef=const V &>
struct two_arrays_reference_traits : public two_arrays_traits<A1,A2>
{
  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class A1,class A2=A1,class A3=A2,class V=typename three_arrays_traits<A1,A2,A3>::value_type,class Ref=V &,class ConstRef=const V &>
struct three_arrays_reference_traits : public three_arrays_traits<A1,A2,A3>
{
  typedef V value_type;
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class V, class Ref=V&,class ConstRef=const V &> class reference_vector_generator;
template<class V, class Ref=V&,class ConstRef=const V &> class reference_matrix_generator;
template<int D,class V,class Ref=V &,class ConstRef=const V &> class reference_array_generator;
template<class V,class Ref,class ConstRef> class reference_array_generator<1,V,Ref,ConstRef> : public reference_vector_generator<V,Ref,ConstRef> { typedef reference_vector_generator<V,Ref,ConstRef> base; };
template<class V,class Ref,class ConstRef> class reference_array_generator<2,V,Ref,ConstRef> : public reference_matrix_generator<V,Ref,ConstRef> { typedef reference_matrix_generator<V,Ref,ConstRef> base; };

template<class T,int Copy=0>
class array_generator : public array_traits<T>
{
  public:
    typedef array_traits<T> base;
    typedef array_generator self;
    ARRAY_BASE_TYPES

  protected:
    array_type &X;

  protected:
    inline array_generator(array_type &x) : X(x) {}
    inline array_generator(const self &r) : X(r.X) {}

  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }
  
    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); } 
    inline int_type  nrows() const { return X.nrows(); } 
    inline int_type  ncols() const { return X.ncols(); } 
    inline void resize(const size_type &s) { assert(s==size()); }       
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); } 
    inline void reset() const { X=value_type(0); }
};

template<class T>
class array_generator<T,1> : public array_traits<T>
{
  public:
    typedef array_traits<T> base;
    typedef array_generator self;
    ARRAY_BASE_TYPES

  protected:
    array_type X;

  protected:
    inline array_generator(const array_type &x) : X(x) {}
    inline array_generator(      self &r) : X(r.X) {}
    inline array_generator(const self &r) : X(r.X) {}

    inline array_generator() : X() {}
    template<class A> inline array_generator(const A &a) : X(a) {}
    template<class A,class B> inline array_generator(const A &a, const B &b) : X(a,b) {}
    template<class A,class B,class C> inline array_generator(const A &a, const B &b, const C &c) : X(a,b,c) {}

  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }
  
    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline int_type  nrows() const { return X.nrows(); } 
    inline int_type  ncols() const { return X.ncols(); } 
    inline void resize(const size_type &s) { assert(s==size()); } 
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline void reset() const { X=value_type(0); }    
};


template<class A1,class A2=A1,int Copy=0>
class two_arrays_generator : public two_arrays_traits<A1,A2>
{
  public:
    typedef two_arrays_traits<A1,A2> base;
    typedef two_arrays_generator self;
    TWO_ARRAYS_BASE_TYPES

  protected:
    first_array_type &X;
    second_array_type &Y;

  protected:
    inline two_arrays_generator(first_array_type &x, second_array_type &y) : X(x), Y(y) {}
    inline two_arrays_generator(const self &r) : X(r.X), Y(r.Y) {}

  public:    
    inline first_array_type        &first_array ()       { return X; }
    inline const first_array_type  &first_array () const { return X; }
    inline second_array_type       &second_array()       { return Y; }
    inline const second_array_type &second_array() const { return Y; }
    
    inline size_type size () const { return X.size (); }      
};

template<class A1,class A2>
class two_arrays_generator<A1,A2,1> : public two_arrays_traits<A1,A2>
{
  public:
    typedef two_arrays_traits<A1,A2> base;
    typedef two_arrays_generator self;
    TWO_ARRAYS_BASE_TYPES    

  protected:
    first_array_type X;
    second_array_type Y;

  protected:
    inline two_arrays_generator(first_array_type &x, second_array_type &y) : X(x), Y(y) {}
    inline two_arrays_generator(const self &r) : X(r.X), Y(r.Y) {}

    inline two_arrays_generator() : X(), Y() {}

  public:    
    inline first_array_type        &first_array ()       { return X; }
    inline const first_array_type  &first_array () const { return X; }
    inline second_array_type       &second_array()       { return Y; }
    inline const second_array_type &second_array() const { return Y; }
    
    inline size_type size () const { return X.size (); }          
};


template<class T,class V=typename T::value_type,int Copy=0>
class array_value_generator : public array_value_traits<T,V>
{
  public:
    typedef array_value_traits<T,V> base;
    typedef array_value_generator self;
    ARRAY_BASE_TYPES    

  protected:
    array_type &X;

  protected:
    inline array_value_generator(array_type &x) : X(x) {}
    inline array_value_generator(const self &r) : X(r.X) {}

  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }

    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &s) { X.resize(s); } // { assert(s==size()); }  
    inline void set_lower_bound(const index_type &i) { X.set_lower_bound(i); }
};

template<class T,class V>
class array_value_generator<T,V,1> : public array_value_traits<T,V>
{
  public:
    typedef array_value_traits<T,V> base;
    typedef array_value_generator self;
    ARRAY_BASE_TYPES    

  protected:
    array_type X;

  protected:
    inline array_value_generator() : X() {}
    inline array_value_generator(array_type &x) : X(x) {}
    inline array_value_generator(const self &r) : X(r.X) {}
    
    //template<class A> inline array_value_generator(const A &a) : X(a) {}
    //template<class A,class B> inline array_value_generator(const A &a, const B &b) : X(a,b) {}

  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }

    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &s) { X.resize(s); } // { assert(s==size()); }    
    inline void set_lower_bound(const index_type &i) { X.set_lower_bound(i); }
};


template<class T,class V=typename T::value_type, class Ref=V &, class ConstRef=const V &, int Copy=0>
class array_reference_generator : public array_reference_traits<T,V,Ref,ConstRef>
{
  public:
    typedef array_reference_traits<T,V,Ref,ConstRef> base;
    typedef array_reference_generator self;
    ARRAY_BASE_TYPES    

  protected:
    array_type &X;

  protected:
    inline array_reference_generator(array_type &x) : X(x) {}
    inline array_reference_generator(const self &r) : X(r.X) {}

  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }

    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &s) { X.resize(s); } // { assert(s==size()); }    
    inline void set_lower_bound(const index_type &i) { X.set_lower_bound(i); }
};

template<class T,class V, class Ref, class ConstRef>
class array_reference_generator<T,V,Ref,ConstRef,1> : public array_reference_traits<T,V,Ref,ConstRef>
{
  public:
    typedef array_reference_traits<T,V,Ref,ConstRef> base;
    typedef array_reference_generator self;
    ARRAY_BASE_TYPES    

  protected:
    array_type X;

  protected:
    inline array_reference_generator(array_type &x) : X(x) {}
    inline array_reference_generator(const self &r) : X(r.X) {}

    inline array_reference_generator() : X() {}
    template<class A> inline array_reference_generator(const A &a) : X(a) {}
    template<class A,class B> inline array_reference_generator(const A &a, const B &b) : X(a,b) {}
    template<class A,class B,class C> inline array_reference_generator(const A &a, const B &b, const C &c) : X(a,b,c) {}

  public:
    inline array_type       &array()       { return X; }
    inline const array_type &array() const { return X; }

    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &s) { X.resize(s); } // { assert(s==size()); }   
    inline void set_lower_bound(const index_type &i) { X.set_lower_bound(i); }
};


template<class A1,class A2=A1,class V=typename two_arrays_traits<A1,A2>::value_type,int Copy=0>
class two_arrays_value_generator : public two_arrays_value_traits<A1,A2,V>
{
  public:
    typedef two_arrays_value_traits<A1,A2,V> base;
    typedef two_arrays_value_generator self;
    TWO_ARRAYS_BASE_TYPES    

  protected:
    first_array_type  &X;
    second_array_type &Y;

  protected:
    inline two_arrays_value_generator(first_array_type &x, second_array_type &y) : X(x), Y(y) {}
    inline two_arrays_value_generator(const self &r) : X(r.X), Y(r.Y) {}

  public:
    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &n) { X.resize(n); Y.resize(n); }  
    
    inline first_array_type        &first_array ()       { return X; }
    inline const first_array_type  &first_array () const { return X; }
    inline second_array_type       &second_array()       { return Y; }
    inline const second_array_type &second_array() const { return Y; }
};

template<class A1,class A2,class V>
class two_arrays_value_generator<A1,A2,V,1> : public two_arrays_value_traits<A1,A2,V>
{
  public:
    typedef two_arrays_value_traits<A1,A2,V> base;
    typedef two_arrays_value_generator self;
    TWO_ARRAYS_BASE_TYPES

  protected:
    first_array_type  X;
    second_array_type Y;

  protected:
    inline two_arrays_value_generator(first_array_type &x, second_array_type &y) : X(x), Y(y) {}
    inline two_arrays_value_generator(const self &r) : X(r.X), Y(r.Y) {}

    inline two_arrays_value_generator() : X(), Y() {}
    template<class T1,class T2> inline two_arrays_value_generator(const T1 &x, const T2 &y) : X(x), Y(y) {}

  public:
    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &n) { X.resize(n); Y.resize(n); }   

    inline first_array_type        &first_array ()       { return X; }
    inline const first_array_type  &first_array () const { return X; }
    inline second_array_type       &second_array()       { return Y; }
    inline const second_array_type &second_array() const { return Y; }
};

template<class A1,class A2=A1,class V=typename two_arrays_traits<A1,A2>::value_type,class Ref=V &,class ConstRef=const V &>
class two_arrays_reference_generator : public two_arrays_reference_traits<A1,A2,V,Ref,ConstRef>
{
  public:
    typedef two_arrays_reference_traits<A1,A2,V,Ref,ConstRef> base;
    typedef two_arrays_reference_generator self;
    TWO_ARRAYS_BASE_TYPES
    
  protected:
    first_array_type  &X;
    second_array_type &Y;

  protected:
    inline two_arrays_reference_generator(first_array_type &x, second_array_type &y) : X(x), Y(y) {}
    inline two_arrays_reference_generator(const self &r) : X(r.X), Y(r.Y) {}

  public:
    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &n) { X.resize(n); Y.resize(n); }   

    inline first_array_type        &first_array ()       { return X; }
    inline const first_array_type  &first_array () const { return X; }
    inline second_array_type       &second_array()       { return Y; }
    inline const second_array_type &second_array() const { return Y; }
};


template<class A1,class A2=A1,class A3=A2,class V=typename three_arrays_traits<A1,A2,A3>::value_type,class Ref=V &,class ConstRef=const V &>
class three_arrays_reference_generator : public three_arrays_reference_traits<A1,A2,A3,V,Ref,ConstRef>
{
  public:
    typedef three_arrays_reference_traits<A1,A2,A3,V,Ref,ConstRef> base;
    typedef three_arrays_reference_generator self;
    THREE_ARRAYS_BASE_TYPES    

  protected:
    first_array_type  &X;
    second_array_type &Y;
    third_array_type  &Z;

  protected:
    inline three_arrays_reference_generator(first_array_type &x, second_array_type &y, third_array_type &z) : X(x), Y(y), Z(z) {}
    inline three_arrays_reference_generator(const self &r) : X(r.X), Y(r.Y), Z(r.Z) {}

  public:
    inline index_type lower_bound() const { return X.lower_bound(); }
    inline size_type size () const { return X.size (); }  
    inline void resize(const size_type &n) { X.resize(n); Y.resize(n); Z.resize(n); }   

    inline first_array_type        &first_array ()       { return X; }
    inline const first_array_type  &first_array () const { return X; }
    inline second_array_type       &second_array()       { return Y; }
    inline const second_array_type &second_array() const { return Y; }    
    inline third_array_type        &third_array ()       { return Z; }
    inline const third_array_type  &third_array () const { return Z; }
};


//Group = Arrays functions


template<class T, int C=0>
class sub_array_generator : public array_generator<T,C>
{
  public:
    typedef array_generator<T,C> base;
    typedef sub_array_generator self;
    ARRAY_BASE_TYPES

  private:
    index_type ind;
    size_type si;

  public:
    inline sub_array_generator(array_type &x, const index_type &p) : base(x), ind(p), si() {}
    inline sub_array_generator(array_type &x, const index_type &p, const size_type &s) : base(x), ind(p), si(s) {}
    inline sub_array_generator(const self &r) : base(static_cast<const base &>(r)), ind(r.ind), si(r.si) {}
    template<class T2,int C2> inline sub_array_generator(const Vector<sub_array_generator<T2,C2> > &x) : base(x.generator().array()), ind(x.generator().pos()), si(x.generator().size()) {}

    using base::array;

    inline index_type lower_bound() const { return index_type(0); }
    inline size_type size() const { return si; }
    inline index_type pos() const { return ind; }
     
    inline reference       operator[](const index_type &i)       { return array()[i+pos()]; }
    inline const_reference operator[](const index_type &i) const { return array()[i+pos()]; }
    
    inline void resize(const size_type &s) { si=s; }
};

template<class T,int Copy> class shift_array_generator;
template<class A,int Copy> class shift_dense_array_generator;

template<      class T,int C=0> struct SubArray { typedef typename GeneratorArray<sub_array_generator<T,C> >::self self; };
template<      class T,int C> struct SubArray<      Vector< sub_array_generator        <T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<      class T,int C> struct SubArray<const Vector< sub_array_generator        <T,C> > > { typedef typename SubArray<const T,C>::self self; };
template<      class T,int C> struct SubArray<      Matrix< sub_array_generator        <T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<      class T,int C> struct SubArray<const Matrix< sub_array_generator        <T,C> > > { typedef typename SubArray<const T,C>::self self; };
template<int D,class T,int C> struct SubArray<      Array<D,sub_array_generator        <T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<int D,class T,int C> struct SubArray<const Array<D,sub_array_generator        <T,C> > > { typedef typename SubArray<const T,C>::self self; };

template<      class T,int C> struct SubArray<      Vector< shift_array_generator      <T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<      class T,int C> struct SubArray<const Vector< shift_array_generator      <T,C> > > { typedef typename SubArray<const T,C>::self self; };
template<      class T,int C> struct SubArray<      Matrix< shift_array_generator      <T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<      class T,int C> struct SubArray<const Matrix< shift_array_generator      <T,C> > > { typedef typename SubArray<const T,C>::self self; };
template<int D,class T,int C> struct SubArray<      Array<D,shift_array_generator      <T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<int D,class T,int C> struct SubArray<const Array<D,shift_array_generator      <T,C> > > { typedef typename SubArray<const T,C>::self self; };

template<      class T,int C> struct SubArray<      Vector< shift_dense_array_generator<T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<      class T,int C> struct SubArray<const Vector< shift_dense_array_generator<T,C> > > { typedef typename SubArray<const T,C>::self self; };
template<      class T,int C> struct SubArray<      Matrix< shift_dense_array_generator<T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<      class T,int C> struct SubArray<const Matrix< shift_dense_array_generator<T,C> > > { typedef typename SubArray<const T,C>::self self; };
template<int D,class T,int C> struct SubArray<      Array<D,shift_dense_array_generator<T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<int D,class T,int C> struct SubArray<const Array<D,shift_dense_array_generator<T,C> > > { typedef typename SubArray<const T,C>::self self; };

//{unsecret}
//summary: gives a part of an array
//Arguments:
//  X - The array to cut out
//  p - Position index
//  i - Vertical position
//  j - Hoirizontal position
//  s - Size
//  m - Height
//  n - Width
//Return: An array representing the part of X
//example:
//  DenseVector<int>::self X(8,"0 1 2 3 4 5 6 7");
//  cout << sub(X,4,2) << endl; // [4 5]
//  cout << sub<2>(X,4) << endl; // [4 5] size fixed to 2
//  DenseMatrix<int>::self Y(4,4,"0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15");
//  typedef DenseMatrix<int>::index_type index_type;
//  typedef DenseMatrix<int>::size_type  size_type;
//  cout << sub(Y,index_type(1,0),size_type(2,2)) << endl; // [4 5; 8 9]
//  cout << sub<2,2>(Y,1,0) << endl; // [4 5; 8 9] size fixed to (2,2)
//See: ^extract^
template<class G> inline typename SubArray<      Vector<G> >::self sub(      Vector<G> &X, const typename Vector<G>::index_type &p)                                         { return typename SubArray<      Vector<G> >::self(X,p  ); }
//{unsecret}
template<class G> inline typename SubArray<const Vector<G> >::self sub(const Vector<G> &X, const typename Vector<G>::index_type &p)                                         { return typename SubArray<const Vector<G> >::self(X,p  ); }
//{unsecret}
template<class G> inline typename SubArray<      Matrix<G> >::self sub(      Matrix<G> &X, const typename Matrix<G>::index_type &p)                                         { return typename SubArray<      Matrix<G> >::self(X,p  ); }
//{unsecret}
template<class G> inline typename SubArray<const Matrix<G> >::self sub(const Matrix<G> &X, const typename Matrix<G>::index_type &p)                                         { return typename SubArray<const Matrix<G> >::self(X,p  ); }
//{unsecret}
template<class G> inline typename SubArray<      Vector<G> >::self sub(      Vector<G> &X, const typename Vector<G>::index_type &p, const typename Vector<G>::size_type &s) { return typename SubArray<      Vector<G> >::self(X,p,s); }
//{unsecret}
template<class G> inline typename SubArray<const Vector<G> >::self sub(const Vector<G> &X, const typename Vector<G>::index_type &p, const typename Vector<G>::size_type &s) { return typename SubArray<const Vector<G> >::self(X,p,s); }
//{unsecret}
template<class G> inline typename SubArray<      Matrix<G> >::self sub(      Matrix<G> &X, const typename Matrix<G>::index_type &p, const typename Matrix<G>::size_type &s) { return typename SubArray<      Matrix<G> >::self(X,p,s); }
//{unsecret}
template<class G> inline typename SubArray<const Matrix<G> >::self sub(const Matrix<G> &X, const typename Matrix<G>::index_type &p, const typename Matrix<G>::size_type &s) { return typename SubArray<const Matrix<G> >::self(X,p,s); }
//{unsecret}
template<class G> inline typename SubArray<      Matrix<G> >::self sub(      Matrix<G> &X, typename Matrix<G>::index_type::int_type i, typename Matrix<G>::index_type::int_type j, typename Matrix<G>::size_type::int_type m,typename Matrix<G>::size_type::int_type n) { return sub(X,typename Matrix<G>::index_type(i,j),typename Matrix<G>::size_type(m,n)); }
//{unsecret}
template<class G> inline typename SubArray<const Matrix<G> >::self sub(const Matrix<G> &X, typename Matrix<G>::index_type::int_type i, typename Matrix<G>::index_type::int_type j, typename Matrix<G>::size_type::int_type m,typename Matrix<G>::size_type::int_type n) { return sub(X,typename Matrix<G>::index_type(i,j),typename Matrix<G>::size_type(m,n)); }

template<class T,int C> inline typename SubArray<      Vector<sub_array_generator  <T,C> > >::self sub(      Vector<sub_array_generator  <T,C> > &X, const typename Vector<sub_array_generator  <T,C> >::index_type &p) { return typename SubArray<      Vector<sub_array_generator  <T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<const Vector<sub_array_generator  <T,C> > >::self sub(const Vector<sub_array_generator  <T,C> > &X, const typename Vector<sub_array_generator  <T,C> >::index_type &p) { return typename SubArray<const Vector<sub_array_generator  <T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<      Matrix<sub_array_generator  <T,C> > >::self sub(      Matrix<sub_array_generator  <T,C> > &X, const typename Matrix<sub_array_generator  <T,C> >::index_type &p) { return typename SubArray<      Matrix<sub_array_generator  <T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<const Matrix<sub_array_generator  <T,C> > >::self sub(const Matrix<sub_array_generator  <T,C> > &X, const typename Matrix<sub_array_generator  <T,C> >::index_type &p) { return typename SubArray<const Matrix<sub_array_generator  <T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<      Vector<sub_array_generator  <T,C> > >::self sub(      Vector<sub_array_generator  <T,C> > &X, const typename Vector<sub_array_generator  <T,C> >::index_type &p, const typename Vector<sub_array_generator  <T,C> >::size_type &s) { return typename SubArray<      Vector<sub_array_generator  <T,C> > >::self(X.generator().array(),X.generator().pos()+p,s); }
template<class T,int C> inline typename SubArray<const Vector<sub_array_generator  <T,C> > >::self sub(const Vector<sub_array_generator  <T,C> > &X, const typename Vector<sub_array_generator  <T,C> >::index_type &p, const typename Vector<sub_array_generator  <T,C> >::size_type &s) { return typename SubArray<const Vector<sub_array_generator  <T,C> > >::self(X.generator().array(),X.generator().pos()+p,s); }
template<class T,int C> inline typename SubArray<      Matrix<sub_array_generator  <T,C> > >::self sub(      Matrix<sub_array_generator  <T,C> > &X, const typename Matrix<sub_array_generator  <T,C> >::index_type &p, const typename Vector<sub_array_generator  <T,C> >::size_type &s) { return typename SubArray<      Matrix<sub_array_generator  <T,C> > >::self(X.generator().array(),X.generator().pos()+p,s); }
template<class T,int C> inline typename SubArray<const Matrix<sub_array_generator  <T,C> > >::self sub(const Matrix<sub_array_generator  <T,C> > &X, const typename Matrix<sub_array_generator  <T,C> >::index_type &p, const typename Vector<sub_array_generator  <T,C> >::size_type &s) { return typename SubArray<const Matrix<sub_array_generator  <T,C> > >::self(X.generator().array(),X.generator().pos()+p,s); }

template<class T,int C> inline typename SubArray<      Vector<shift_array_generator<T,C> > >::self sub(      Vector<shift_array_generator<T,C> > &X, const typename Vector<shift_array_generator<T,C> >::index_type &p) { return typename SubArray<      Vector<shift_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<const Vector<shift_array_generator<T,C> > >::self sub(const Vector<shift_array_generator<T,C> > &X, const typename Vector<shift_array_generator<T,C> >::index_type &p) { return typename SubArray<const Vector<shift_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<      Matrix<shift_array_generator<T,C> > >::self sub(      Matrix<shift_array_generator<T,C> > &X, const typename Matrix<shift_array_generator<T,C> >::index_type &p) { return typename SubArray<      Matrix<shift_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<const Matrix<shift_array_generator<T,C> > >::self sub(const Matrix<shift_array_generator<T,C> > &X, const typename Matrix<shift_array_generator<T,C> >::index_type &p) { return typename SubArray<const Matrix<shift_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<      Vector<shift_array_generator<T,C> > >::self sub(      Vector<shift_array_generator<T,C> > &X, const typename Vector<shift_array_generator<T,C> >::index_type &p, const typename Vector<shift_array_generator<T,C> >::size_type &s) { return typename SubArray<      Vector<shift_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }
template<class T,int C> inline typename SubArray<const Vector<shift_array_generator<T,C> > >::self sub(const Vector<shift_array_generator<T,C> > &X, const typename Vector<shift_array_generator<T,C> >::index_type &p, const typename Vector<shift_array_generator<T,C> >::size_type &s) { return typename SubArray<const Vector<shift_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }
template<class T,int C> inline typename SubArray<      Matrix<shift_array_generator<T,C> > >::self sub(      Matrix<shift_array_generator<T,C> > &X, const typename Matrix<shift_array_generator<T,C> >::index_type &p, const typename Vector<shift_array_generator<T,C> >::size_type &s) { return typename SubArray<      Matrix<shift_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }
template<class T,int C> inline typename SubArray<const Matrix<shift_array_generator<T,C> > >::self sub(const Matrix<shift_array_generator<T,C> > &X, const typename Matrix<shift_array_generator<T,C> >::index_type &p, const typename Vector<shift_array_generator<T,C> >::size_type &s) { return typename SubArray<const Matrix<shift_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }

template<class T,int C> inline typename SubArray<      Vector<shift_dense_array_generator<T,C> > >::self sub(      Vector<shift_dense_array_generator<T,C> > &X, const typename Vector<shift_dense_array_generator<T,C> >::index_type &p) { return typename SubArray<      Vector<shift_dense_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<const Vector<shift_dense_array_generator<T,C> > >::self sub(const Vector<shift_dense_array_generator<T,C> > &X, const typename Vector<shift_dense_array_generator<T,C> >::index_type &p) { return typename SubArray<const Vector<shift_dense_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<      Matrix<shift_dense_array_generator<T,C> > >::self sub(      Matrix<shift_dense_array_generator<T,C> > &X, const typename Matrix<shift_dense_array_generator<T,C> >::index_type &p) { return typename SubArray<      Matrix<shift_dense_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<const Matrix<shift_dense_array_generator<T,C> > >::self sub(const Matrix<shift_dense_array_generator<T,C> > &X, const typename Matrix<shift_dense_array_generator<T,C> >::index_type &p) { return typename SubArray<const Matrix<shift_dense_array_generator<T,C> > >::self(X.generator().array(), p-X.generator().index()); }
template<class T,int C> inline typename SubArray<      Vector<shift_dense_array_generator<T,C> > >::self sub(      Vector<shift_dense_array_generator<T,C> > &X, const typename Vector<shift_dense_array_generator<T,C> >::index_type &p, const typename Vector<shift_dense_array_generator<T,C> >::size_type &s) { return typename SubArray<      Vector<shift_dense_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }
template<class T,int C> inline typename SubArray<const Vector<shift_dense_array_generator<T,C> > >::self sub(const Vector<shift_dense_array_generator<T,C> > &X, const typename Vector<shift_dense_array_generator<T,C> >::index_type &p, const typename Vector<shift_dense_array_generator<T,C> >::size_type &s) { return typename SubArray<const Vector<shift_dense_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }
template<class T,int C> inline typename SubArray<      Matrix<shift_dense_array_generator<T,C> > >::self sub(      Matrix<shift_dense_array_generator<T,C> > &X, const typename Matrix<shift_dense_array_generator<T,C> >::index_type &p, const typename Vector<shift_dense_array_generator<T,C> >::size_type &s) { return typename SubArray<      Matrix<shift_dense_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }
template<class T,int C> inline typename SubArray<const Matrix<shift_dense_array_generator<T,C> > >::self sub(const Matrix<shift_dense_array_generator<T,C> > &X, const typename Matrix<shift_dense_array_generator<T,C> >::index_type &p, const typename Vector<shift_dense_array_generator<T,C> >::size_type &s) { return typename SubArray<const Matrix<shift_dense_array_generator<T,C> > >::self(X.generator().array(),p-X.generator().index(),s); }

template<class T,int C=0>
class sub_dense_vector_generator : public sub_array_generator<T,C>
{
  public:
    typedef sub_array_generator<T,C> base;
    typedef sub_dense_vector_generator self;
    ARRAY_BASE_TYPES
        
  protected:
    pointer data; 
  
  public:
    inline sub_dense_vector_generator(array_type &x, const index_type &p) : base(x,p), data(&x[p]) {}
    inline sub_dense_vector_generator(array_type &x, const index_type &p, const size_type &s) : base(x,p,s), data(&x[p]) {}
         
    template<class T2> inline sub_dense_vector_generator(const Vector<sub_dense_vector_generator<T2> > &x) : base((base &)(x.generator())), data(&x[0]) {}
         
    inline reference       operator[](const index_type &i)       { return data[i]; }
    inline const_reference operator[](const index_type &i) const { return data[i]; }
};

template<class T,int N,int Copy=0>
class sub_simd_dense_vector_generator : public sub_array_generator<T,Copy>
{
  public:
    typedef sub_array_generator<T,Copy> base;
    typedef sub_simd_dense_vector_generator self;
    ARRAY_BASE_TYPES
    
    typedef typename array_type::generator_type::template_reference template_reference;    
        
  public:
    typename array_traits<reference>::pointer data; 
  
  public:
    inline sub_simd_dense_vector_generator(array_type &x, const index_type &i                    ) : base(x,i  ), data(x.generator().data+N*i) {} // data(&x.generator().array()[N*i])
    inline sub_simd_dense_vector_generator(array_type &x, const index_type &i, const size_type &s) : base(x,i,s), data(x.generator().data+N*i) {} // data(&x.generator().array()[N*i]) 

    inline reference       operator[](const index_type &i)       { return template_reference((typename template_reference::pointer)(data+N*i)); }
    inline const_reference operator[](const index_type &i) const { return template_reference((typename template_reference::pointer)(data+N*i)); }
};


#ifdef NDEBUG
template<        class V,int C> struct SubArray<      Vector<data_vector_generator <V  > >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<      Vector<data_vector_generator <V  > >,C> >::self self; };
template<        class V,int C> struct SubArray<const Vector<data_vector_generator <V  > >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<const Vector<data_vector_generator <V  > >,C> >::self self; };
template<class V,class A,int C> struct SubArray<      Vector<dense_vector_generator<V,A> >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<      Vector<dense_vector_generator<V,A> >,C> >::self self; };
template<class V,class A,int C> struct SubArray<const Vector<dense_vector_generator<V,A> >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<const Vector<dense_vector_generator<V,A> >,C> >::self self; };
template<int   N,class V,int C> struct SubArray<      Vector< tiny_vector_generator<N,V> >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<      Vector< tiny_vector_generator<N,V> >,C> >::self self; };
template<int   N,class V,int C> struct SubArray<const Vector< tiny_vector_generator<N,V> >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<const Vector< tiny_vector_generator<N,V> >,C> >::self self; };

template<class T,int C> class row_dense_matrix_generator;
template<class T,int C2,int C> struct SubArray<      Vector<row_dense_matrix_generator<T,C2> >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<      Vector<row_dense_matrix_generator<T,C2> >,C> >::self self; };
template<class T,int C2,int C> struct SubArray<const Vector<row_dense_matrix_generator<T,C2> >,C> { typedef typename GeneratorArray<sub_dense_vector_generator<const Vector<row_dense_matrix_generator<T,C2> >,C> >::self self; };

template<int N,class T,template<int,class> class Cast> class simd_block_dense_reference;
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef=simd_block_dense_reference,int C=0> class simd_block_dense_vector_generator;
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct SubArray<      Vector<simd_block_dense_vector_generator<T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<sub_simd_dense_vector_generator<      Vector<simd_block_dense_vector_generator<T,N,Cast,R,C2> >,N,C> >::self self; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct SubArray<const Vector<simd_block_dense_vector_generator<T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<sub_simd_dense_vector_generator<const Vector<simd_block_dense_vector_generator<T,N,Cast,R,C2> >,N,C> >::self self; };

template<class T,int N,int C> class row_simd_dense_matrix_generator;
template<class T,int N,int C2,int C> struct SubArray<      Vector<row_simd_dense_matrix_generator<T,N,C2> >,C> { typedef typename GeneratorArray<sub_simd_dense_vector_generator<      Vector<row_simd_dense_matrix_generator<T,N,C2> >,N,C> >::self self; };
template<class T,int N,int C2,int C> struct SubArray<const Vector<row_simd_dense_matrix_generator<T,N,C2> >,C> { typedef typename GeneratorArray<sub_simd_dense_vector_generator<const Vector<row_simd_dense_matrix_generator<T,N,C2> >,N,C> >::self self; };
#endif




template<int N,class T, int Copy=0>
class tiny_sub_vector_generator : public array_generator<T,Copy>
{
  public:
    typedef array_generator<T,Copy> base;  
    typedef tiny_sub_vector_generator self; 
    ARRAY_BASE_TYPES

    template<class V2> struct array_rebind { typedef typename TinyVector<N,V2>::self other; };

  private:
    index_type ind;

  public:
    inline tiny_sub_vector_generator(array_type &x, const index_type &p) : base(x), ind(p) {}
    inline tiny_sub_vector_generator(const self &x) : base(static_cast<const base &>(x)), ind(x.ind) {}
    
    template<class T2,int C2> inline tiny_sub_vector_generator(const Vector<tiny_sub_vector_generator<N,T2,C2> > &x) : base(x.generator().array()), ind(x.generator().pos()) {}

    using base::array;
    inline index_type lower_bound() const { return index_type(0); }
    inline size_type size() const { return N; }
    inline index_type pos() const { return ind; }
     
    inline reference       operator[](const index_type &i)       { return array()[i+pos()]; }
    inline const_reference operator[](const index_type &i) const { return array()[i+pos()]; }
    
    inline void resize(const size_type &s) { assert(s==N); }
};

template<class T,int C,class G2> inline void affect(Vector<tiny_sub_vector_generator<1,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0];  }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_vector_generator<2,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_vector_generator<3,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_vector_generator<4,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; }

template<int N,       class T,     int C=0> struct TinySubVector { typedef typename GeneratorArray<tiny_sub_vector_generator<N,T,C> >::self self; inline self operator()(T &X, const typename T::index_type &p) const {return self(X,p); }  };
template<int N,       class T,int C2,int C> struct TinySubVector<N,      Vector<sub_array_generator        <   T,C2> >,C> { typedef typename TinySubVector<N,      T,C>::self self; inline self operator()(      Vector<sub_array_generator        <   T,C2> > &X, const typename Vector<sub_array_generator        <   T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,       class T,int C2,int C> struct TinySubVector<N,const Vector<sub_array_generator        <   T,C2> >,C> { typedef typename TinySubVector<N,const T,C>::self self; inline self operator()(const Vector<sub_array_generator        <   T,C2> > &X, const typename Vector<sub_array_generator        <   T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,int N2,class T,int C2,int C> struct TinySubVector<N,      Vector<tiny_sub_vector_generator  <N2,T,C2> >,C> { typedef typename TinySubVector<N,      T,C>::self self; inline self operator()(      Vector<tiny_sub_vector_generator  <N2,T,C2> > &X, const typename Vector<tiny_sub_vector_generator  <N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,int N2,class T,int C2,int C> struct TinySubVector<N,const Vector<tiny_sub_vector_generator  <N2,T,C2> >,C> { typedef typename TinySubVector<N,const T,C>::self self; inline self operator()(const Vector<tiny_sub_vector_generator  <N2,T,C2> > &X, const typename Vector<tiny_sub_vector_generator  <N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,       class T,int C2,int C> struct TinySubVector<N,      Vector<shift_array_generator      <   T,C2> >,C> { typedef typename TinySubVector<N,      T,C>::self self; inline self operator()(      Vector<shift_array_generator      <   T,C2> > &X, const typename Vector<shift_array_generator      <   T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };
template<int N,       class T,int C2,int C> struct TinySubVector<N,const Vector<shift_array_generator      <   T,C2> >,C> { typedef typename TinySubVector<N,const T,C>::self self; inline self operator()(const Vector<shift_array_generator      <   T,C2> > &X, const typename Vector<shift_array_generator      <   T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };
template<int N,       class T,int C2,int C> struct TinySubVector<N,      Vector<shift_dense_array_generator<   T,C2> >,C> { typedef typename TinySubVector<N,      T,C>::self self; inline self operator()(      Vector<shift_dense_array_generator<   T,C2> > &X, const typename Vector<shift_dense_array_generator<   T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };
template<int N,       class T,int C2,int C> struct TinySubVector<N,const Vector<shift_dense_array_generator<   T,C2> >,C> { typedef typename TinySubVector<N,const T,C>::self self; inline self operator()(const Vector<shift_dense_array_generator<   T,C2> > &X, const typename Vector<shift_dense_array_generator<   T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };

template<int N,class G> inline typename TinySubVector<N,      Vector<G> >::self sub(      Vector<G> &X, const typename Vector<G>::index_type &p) { return TinySubVector<N,      Vector<G> >()(X,p); }
template<int N,class G> inline typename TinySubVector<N,const Vector<G> >::self sub(const Vector<G> &X, const typename Vector<G>::index_type &p) { return TinySubVector<N,const Vector<G> >()(X,p); }


template<int N,class T,int Copy=0>
class tiny_sub_dense_vector_generator : public tiny_sub_vector_generator<N,T,Copy>
{
  public:
    typedef tiny_sub_vector_generator<N,T,Copy> base;
    typedef tiny_sub_dense_vector_generator self;
    ARRAY_BASE_TYPES
        
    template<class V2> struct array_rebind { typedef typename TinyVector<N,V2>::self other; };
        
  protected:
    pointer data; 
  
  public:
    inline tiny_sub_dense_vector_generator(array_type &x, const index_type &p) : base(x,p), data(&x[p]) {}
    template<class T2,int C2> inline tiny_sub_dense_vector_generator(const Vector<tiny_sub_dense_vector_generator<N,T2,C2> > &x) : base((base &)(x.generator())), data(const_cast<pointer>(&x[0])) {}

    inline index_type lower_bound() const { return index_type(0); }
    inline size_type size() const { return N; }
    inline void resize(const size_type &s) { assert(s==N); }

    inline reference       operator[](const index_type &i)       
    {
      return data[i]; 
    }
    inline const_reference operator[](const index_type &i) const { return data[i]; }
};

template<class T,int C,class G2> inline void affect(Vector<tiny_sub_dense_vector_generator<1,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_dense_vector_generator<2,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_dense_vector_generator<3,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_dense_vector_generator<4,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; }

template<int N,       class T,int C2,int C> struct TinySubVector<N,      Vector<sub_dense_vector_generator     <   T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      T,C> >::self self; inline self operator()(      Vector<sub_dense_vector_generator     <   T,C2> > &X, const typename Vector<sub_dense_vector_generator     <   T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,       class T,int C2,int C> struct TinySubVector<N,const Vector<sub_dense_vector_generator     <   T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const T,C> >::self self; inline self operator()(const Vector<sub_dense_vector_generator     <   T,C2> > &X, const typename Vector<sub_dense_vector_generator     <   T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,int N2,class T,int C2,int C> struct TinySubVector<N,      Vector<tiny_sub_dense_vector_generator<N2,T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      T,C> >::self self; inline self operator()(      Vector<tiny_sub_dense_vector_generator<N2,T,C2> > &X, const typename Vector<tiny_sub_dense_vector_generator<N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int N,int N2,class T,int C2,int C> struct TinySubVector<N,const Vector<tiny_sub_dense_vector_generator<N2,T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const T,C> >::self self; inline self operator()(const Vector<tiny_sub_dense_vector_generator<N2,T,C2> > &X, const typename Vector<tiny_sub_dense_vector_generator<N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };

#ifdef NDEBUG
template<int N,class V,        int C> struct TinySubVector<N,      Vector<data_vector_generator   <V  > >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      Vector<data_vector_generator   <V  > >,C> >::self self; inline self operator()(      Vector<data_vector_generator  <V  > > &X, const typename Vector<dense_vector_generator  <V  > >::index_type &p) const { return self(X,p); } };
template<int N,class V,        int C> struct TinySubVector<N,const Vector<data_vector_generator   <V  > >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const Vector<data_vector_generator   <V  > >,C> >::self self; inline self operator()(const Vector<data_vector_generator  <V  > > &X, const typename Vector<dense_vector_generator  <V  > >::index_type &p) const { return self(X,p); } };
template<int N,class V,class A,int C> struct TinySubVector<N,      Vector<dense_vector_generator  <V,A> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      Vector<dense_vector_generator  <V,A> >,C> >::self self; inline self operator()(      Vector<dense_vector_generator <V,A> > &X, const typename Vector<dense_vector_generator  <V,A> >::index_type &p) const { return self(X,p); } };
template<int N,class V,class A,int C> struct TinySubVector<N,const Vector<dense_vector_generator  <V,A> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const Vector<dense_vector_generator  <V,A> >,C> >::self self; inline self operator()(const Vector<dense_vector_generator <V,A> > &X, const typename Vector<dense_vector_generator  <V,A> >::index_type &p) const { return self(X,p); } };
template<int N,int N2 ,class V,int C> struct TinySubVector<N,      Vector<tiny_vector_generator  <N2,V> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      Vector<tiny_vector_generator  <N2,V> >,C> >::self self; inline self operator()(      Vector<tiny_vector_generator <N2,V> > &X, const typename Vector<tiny_vector_generator  <N2,V> >::index_type &p) const { return self(X,p); } };
template<int N,int N2 ,class V,int C> struct TinySubVector<N,const Vector<tiny_vector_generator  <N2,V> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const Vector<tiny_vector_generator  <N2,V> >,C> >::self self; inline self operator()(const Vector<tiny_vector_generator <N2,V> > &X, const typename Vector<tiny_vector_generator  <N2,V> >::index_type &p) const { return self(X,p); } };
template<      class T,int C=0> class row_dense_matrix_generator;
template<int N,class T,int C=0> class tiny_row_dense_matrix_generator;
template<int N,       class T, int C> struct TinySubVector<N,      Vector<row_dense_matrix_generator     <   T> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      Vector<row_dense_matrix_generator     <   T> >,C> >::self self; inline self operator()(      Vector<row_dense_matrix_generator     <T   > >&X, const typename Vector<row_dense_matrix_generator     <   T> >::index_type &p) const { return self(X,p); } };
template<int N,       class T, int C> struct TinySubVector<N,const Vector<row_dense_matrix_generator     <   T> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const Vector<row_dense_matrix_generator     <   T> >,C> >::self self; inline self operator()(const Vector<row_dense_matrix_generator     <T   > >&X, const typename Vector<row_dense_matrix_generator     <   T> >::index_type &p) const { return self(X,p); } };
template<int N,int N2,class T, int C> struct TinySubVector<N,      Vector<tiny_row_dense_matrix_generator<N2,T> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,      Vector<tiny_row_dense_matrix_generator<N2,T> >,C> >::self self; inline self operator()(      Vector<tiny_row_dense_matrix_generator<N2,T> >&X, const typename Vector<tiny_row_dense_matrix_generator<N2,T> >::index_type &p) const { return self(X,p); } };
template<int N,int N2,class T, int C> struct TinySubVector<N,const Vector<tiny_row_dense_matrix_generator<N2,T> >,C> { typedef typename GeneratorArray<tiny_sub_dense_vector_generator<N,const Vector<tiny_row_dense_matrix_generator<N2,T> >,C> >::self self; inline self operator()(const Vector<tiny_row_dense_matrix_generator<N2,T> >&X, const typename Vector<tiny_row_dense_matrix_generator<N2,T> >::index_type &p) const { return self(X,p); } };
#endif



template<class T,int C> struct SubArray<      Vector<sub_dense_vector_generator<T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<class T,int C> struct SubArray<const Vector<sub_dense_vector_generator<T,C> > > { typedef typename SubArray<const T,C>::self self; };

template<int N,class T,int C> struct SubArray<      Vector<tiny_sub_dense_vector_generator<N,T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<int N,class T,int C> struct SubArray<const Vector<tiny_sub_dense_vector_generator<N,T,C> > > { typedef typename SubArray<const T,C>::self self; };

template<class T,int N,int C> struct SubArray<      Vector<sub_simd_dense_vector_generator<T,N,C> > > { typedef typename SubArray<      T,C>::self self; };
template<class T,int N,int C> struct SubArray<const Vector<sub_simd_dense_vector_generator<T,N,C> > > { typedef typename SubArray<const T,C>::self self; };


template<class T,int C> inline typename SubArray<      Vector<sub_dense_vector_generator<T,C> > >::self sub(      Vector<sub_dense_vector_generator<T,C> > &X, const typename Vector<sub_dense_vector_generator<T,C> >::index_type &p)                                                                        { return typename SubArray<      Vector<sub_dense_vector_generator<T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<const Vector<sub_dense_vector_generator<T,C> > >::self sub(const Vector<sub_dense_vector_generator<T,C> > &X, const typename Vector<sub_dense_vector_generator<T,C> >::index_type &p)                                                                        { return typename SubArray<const Vector<sub_dense_vector_generator<T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<      Vector<sub_dense_vector_generator<T,C> > >::self sub(      Vector<sub_dense_vector_generator<T,C> > &X, const typename Vector<sub_dense_vector_generator<T,C> >::index_type &p, const typename Vector<sub_dense_vector_generator<T,C> >::size_type &s) { return typename SubArray<      Vector<sub_dense_vector_generator<T,C> > >::self(X.generator().array(), X.generator().pos()+p,s); }
template<class T,int C> inline typename SubArray<const Vector<sub_dense_vector_generator<T,C> > >::self sub(const Vector<sub_dense_vector_generator<T,C> > &X, const typename Vector<sub_dense_vector_generator<T,C> >::index_type &p, const typename Vector<sub_dense_vector_generator<T,C> >::size_type &s) { return typename SubArray<const Vector<sub_dense_vector_generator<T,C> > >::self(X.generator().array(), X.generator().pos()+p,s); }

template<int N,class T,int C> inline typename SubArray<      Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self sub(      Vector<tiny_sub_dense_vector_generator<N,T,C> > &X, const typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::index_type &p)                                                                               { return typename SubArray<      Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<int N,class T,int C> inline typename SubArray<const Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self sub(const Vector<tiny_sub_dense_vector_generator<N,T,C> > &X, const typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::index_type &p)                                                                               { return typename SubArray<const Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<int N,class T,int C> inline typename SubArray<      Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self sub(      Vector<tiny_sub_dense_vector_generator<N,T,C> > &X, const typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::index_type &p, const typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::size_type &s) { return typename SubArray<      Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self(X.generator().array(), X.generator().pos()+p,s); }
template<int N,class T,int C> inline typename SubArray<const Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self sub(const Vector<tiny_sub_dense_vector_generator<N,T,C> > &X, const typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::index_type &p, const typename Vector<tiny_sub_dense_vector_generator<N,T,C> >::size_type &s) { return typename SubArray<const Vector<tiny_sub_dense_vector_generator<N,T,C> > >::self(X.generator().array(), X.generator().pos()+p,s); }

template<class T,int N,int C> inline typename SubArray<      Vector<sub_simd_dense_vector_generator<T,N,C> > >::self sub(      Vector<sub_simd_dense_vector_generator<T,N,C> > &X, const typename Vector<sub_simd_dense_vector_generator<T,N,C> >::index_type &p)                                                                               { return typename SubArray<      Vector<sub_simd_dense_vector_generator<T,N,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int N,int C> inline typename SubArray<const Vector<sub_simd_dense_vector_generator<T,N,C> > >::self sub(const Vector<sub_simd_dense_vector_generator<T,N,C> > &X, const typename Vector<sub_simd_dense_vector_generator<T,N,C> >::index_type &p)                                                                               { return typename SubArray<const Vector<sub_simd_dense_vector_generator<T,N,C> > >::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int N,int C> inline typename SubArray<      Vector<sub_simd_dense_vector_generator<T,N,C> > >::self sub(      Vector<sub_simd_dense_vector_generator<T,N,C> > &X, const typename Vector<sub_simd_dense_vector_generator<T,N,C> >::index_type &p, const typename Vector<sub_simd_dense_vector_generator<T,N,C> >::size_type &s) { return typename SubArray<      Vector<sub_simd_dense_vector_generator<T,N,C> > >::self(X.generator().array(), X.generator().pos()+p,s); }
template<class T,int N,int C> inline typename SubArray<const Vector<sub_simd_dense_vector_generator<T,N,C> > >::self sub(const Vector<sub_simd_dense_vector_generator<T,N,C> > &X, const typename Vector<sub_simd_dense_vector_generator<T,N,C> >::index_type &p, const typename Vector<sub_simd_dense_vector_generator<T,N,C> >::size_type &s) { return typename SubArray<const Vector<sub_simd_dense_vector_generator<T,N,C> > >::self(X.generator().array(), X.generator().pos()+p,s); }


template<class T,int Copy=0>
class sub_dense_matrix_generator : public sub_array_generator<T,Copy>
{
  public:
    typedef sub_array_generator<T,Copy> base;
    typedef sub_dense_matrix_generator self;
    ARRAY_BASE_TYPES
        
  private:
    pointer data; 
  
  public:
    inline sub_dense_matrix_generator(array_type &x, const index_type &p) : base(x,p), data(&x[p]) {}
    inline sub_dense_matrix_generator(array_type &x, const index_type &p, const size_type &s) : base(x,p,s), data(&x[p]) {}
         
    template<class T2> inline sub_dense_matrix_generator(const Matrix<sub_dense_matrix_generator<T2> > &x) : base((base &)(x.generator())), data(&x[0]) {}
    
    using base::array;

    int_type row_stride() const { return array().generator().row_stride(); }
         
    inline reference       operator[](const index_type &i)       { return data[array().ncols()*i.i+i.j]; }
    inline const_reference operator[](const index_type &i) const { return data[array().ncols()*i.i+i.j]; }
};

template<class T,int C> struct SubArray<      Matrix< sub_dense_matrix_generator<T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<class T,int C> struct SubArray<const Matrix< sub_dense_matrix_generator<T,C> > > { typedef typename SubArray<const T,C>::self self; };

template<int M,int N,class T,int C> class tiny_sub_dense_matrix_generator;
template<int M,int N,class T,int C> struct SubArray<      Matrix< tiny_sub_dense_matrix_generator<M,N,T,C> > > { typedef typename SubArray<      T,C>::self self; };
template<int M,int N,class T,int C> struct SubArray<const Matrix< tiny_sub_dense_matrix_generator<M,N,T,C> > > { typedef typename SubArray<const T,C>::self self; };

#ifdef NDEBUG
template<      class V,        int C> struct SubArray<      Matrix<data_matrix_generator <V  > >,C> { typedef typename GeneratorArray<sub_dense_matrix_generator<      Matrix<data_matrix_generator <  V  > >,C> >::self self; };
template<      class V,        int C> struct SubArray<const Matrix<data_matrix_generator <V  > >,C> { typedef typename GeneratorArray<sub_dense_matrix_generator<const Matrix<data_matrix_generator <  V  > >,C> >::self self; };
template<      class V,class A,int C> struct SubArray<      Matrix<dense_matrix_generator<V,A> >,C> { typedef typename GeneratorArray<sub_dense_matrix_generator<      Matrix<dense_matrix_generator<  V,A> >,C> >::self self; };
template<      class V,class A,int C> struct SubArray<const Matrix<dense_matrix_generator<V,A> >,C> { typedef typename GeneratorArray<sub_dense_matrix_generator<const Matrix<dense_matrix_generator<  V,A> >,C> >::self self; };
template<int M,int   N,class V,int C> struct SubArray<      Matrix<tiny_matrix_generator <M,N,V> >,C> { typedef typename GeneratorArray<sub_dense_matrix_generator<      Matrix< tiny_matrix_generator<M,N,V> >,C> >::self self; };
template<int M,int   N,class V,int C> struct SubArray<const Matrix<tiny_matrix_generator <M,N,V> >,C> { typedef typename GeneratorArray<sub_dense_matrix_generator<const Matrix< tiny_matrix_generator<M,N,V> >,C> >::self self; };
#endif



template<class T,int N,int Copy=0>
class sub_simd_dense_matrix_generator : public sub_array_generator<T,Copy>
{
  public:
    typedef sub_array_generator<T,Copy> base;
    typedef sub_simd_dense_matrix_generator self;
    ARRAY_BASE_TYPES

    typedef typename array_type::generator_type::template_reference template_reference;    
        
  public:
    typename array_traits<reference>::pointer data; 
  
  public:
    inline sub_simd_dense_matrix_generator(array_type &x, const index_type &i                    ) : base(x,i  ), data(&x.generator().array()(i.i,N*i.j)) {}
    inline sub_simd_dense_matrix_generator(array_type &x, const index_type &i, const size_type &s) : base(x,i,s), data(&x.generator().array()(i.i,N*i.j)) {}
         
    using base::array;
    
    inline reference       operator[](const index_type &i)       { return template_reference((typename template_reference::pointer)(data+array().ncols()*i.i+N*i.j)); }
    inline const_reference operator[](const index_type &i) const { return template_reference((typename template_reference::pointer)(data+array().ncols()*i.i+N*i.j)); }
    
    int_type stride() const { return array().generator().stride(); }
};

template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef,int C> struct simd_block_dense_matrix_generator;
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct SubArray<      Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<sub_simd_dense_matrix_generator<      Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,N,C> >::self self; };
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class R,int C2,int C> struct SubArray<const Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,C> { typedef typename GeneratorArray<sub_simd_dense_matrix_generator<const Matrix<simd_block_dense_matrix_generator<T,N,Cast,R,C2> >,N,C> >::self self; };




template<class T,int C> inline typename SubArray<      T>::self sub(      Matrix<sub_dense_matrix_generator<T,C> > &X, const typename Matrix<sub_dense_matrix_generator<T,C> >::index_type &p) { return typename SubArray<      T>::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<const T>::self sub(const Matrix<sub_dense_matrix_generator<T,C> > &X, const typename Matrix<sub_dense_matrix_generator<T,C> >::index_type &p) { return typename SubArray<const T>::self(X.generator().array(), X.generator().pos()+p); }
template<class T,int C> inline typename SubArray<      T>::self sub(      Matrix<sub_dense_matrix_generator<T,C> > &X, const typename Matrix<sub_dense_matrix_generator<T,C> >::index_type &p, const typename Matrix<sub_dense_matrix_generator<T,C> >::size_type &s) { return typename SubArray<      T>::self(X.generator().array(), X.generator().pos()+p,s); }
template<class T,int C> inline typename SubArray<const T>::self sub(const Matrix<sub_dense_matrix_generator<T,C> > &X, const typename Matrix<sub_dense_matrix_generator<T,C> >::index_type &p, const typename Matrix<sub_dense_matrix_generator<T,C> >::size_type &s) { return typename SubArray<const T>::self(X.generator().array(), X.generator().pos()+p,s); }

template<int M,int N,class T,int C> inline typename SubArray<      T>::self sub(      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> > &X, const typename Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> >::index_type &p) { return typename SubArray<      T>::self(X.generator().array(), X.generator().pos()+p); }
template<int M,int N,class T,int C> inline typename SubArray<const T>::self sub(const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> > &X, const typename Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> >::index_type &p) { return typename SubArray<const T>::self(X.generator().array(), X.generator().pos()+p); }
template<int M,int N,class T,int C> inline typename SubArray<      T>::self sub(      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> > &X, const typename Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> >::index_type &p, const typename Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> >::size_type &s) { return typename SubArray<      T>::self(X.generator().array(), X.generator().pos()+p,s); }
template<int M,int N,class T,int C> inline typename SubArray<const T>::self sub(const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> > &X, const typename Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> >::index_type &p, const typename Matrix<tiny_sub_dense_matrix_generator<M,N,T,C> >::size_type &s) { return typename SubArray<const T>::self(X.generator().array(), X.generator().pos()+p,s); }


template<int M,int N,class T, int Copy=0>
class tiny_sub_matrix_generator : public array_generator<T,Copy>
{
  public:
    typedef array_generator<T,Copy> base;
    typedef tiny_sub_matrix_generator self; 
    ARRAY_BASE_TYPES

    template<class V2> struct array_rebind { typedef Matrix<tiny_matrix_generator<M,N,V2> > other; };

  private:
    index_type ind;

  public:
    using base::array;
    
    inline tiny_sub_matrix_generator(array_type &x, const index_type &p) : base(x), ind(p) {}
    inline tiny_sub_matrix_generator(const self &x) : base(x), ind(x.ind) {}
    
    template<class T2,int C2> inline tiny_sub_matrix_generator(const Matrix<tiny_sub_matrix_generator<M,N,T2,C2> > &x) : base(x.generator().array()), ind(x.generator().pos()) {}

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return size_type(M,N); }
    inline index_type pos() const { return ind; }
     
    inline reference       operator[](const index_type &i)       { return array()[i+pos()]; }
    inline const_reference operator[](const index_type &i) const { return array()[i+pos()]; }
    
    inline void resize(const size_type &s) { assert(s==size_type(M,N)); }
};

template<int M,int N,              class T,     int C=0> struct TinySubMatrix { typedef typename GeneratorArray<tiny_sub_matrix_generator<M,N,T,C> >::self self; inline self operator()(T &X, const typename T::index_type &p) const {return self(X,p); }  };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,      Matrix<sub_array_generator        <      T,C2> >,C> { typedef typename TinySubMatrix<M,N,      T,C>::self self; inline self operator()(      Matrix<sub_array_generator        <      T,C2> > &X, const typename Matrix<sub_array_generator        <      T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,const Matrix<sub_array_generator        <      T,C2> >,C> { typedef typename TinySubMatrix<M,N,const T,C>::self self; inline self operator()(const Matrix<sub_array_generator        <      T,C2> > &X, const typename Matrix<sub_array_generator        <      T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,int M2,int N2,class T,int C2,int C> struct TinySubMatrix<M,N,      Matrix<tiny_sub_matrix_generator  <M2,N2,T,C2> >,C> { typedef typename TinySubMatrix<M,N,      T,C>::self self; inline self operator()(      Matrix<tiny_sub_matrix_generator  <M2,N2,T,C2> > &X, const typename Matrix<tiny_sub_matrix_generator  <M2,N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,int M2,int N2,class T,int C2,int C> struct TinySubMatrix<M,N,const Matrix<tiny_sub_matrix_generator  <M2,N2,T,C2> >,C> { typedef typename TinySubMatrix<M,N,const T,C>::self self; inline self operator()(const Matrix<tiny_sub_matrix_generator  <M2,N2,T,C2> > &X, const typename Matrix<tiny_sub_matrix_generator  <M2,N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,      Matrix<shift_array_generator      <      T,C2> >,C> { typedef typename TinySubMatrix<M,N,      T,C>::self self; inline self operator()(      Matrix<shift_array_generator      <      T,C2> > &X, const typename Matrix<shift_array_generator      <      T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,const Matrix<shift_array_generator      <      T,C2> >,C> { typedef typename TinySubMatrix<M,N,const T,C>::self self; inline self operator()(const Matrix<shift_array_generator      <      T,C2> > &X, const typename Matrix<shift_array_generator      <      T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,      Matrix<shift_dense_array_generator<      T,C2> >,C> { typedef typename TinySubMatrix<M,N,      T,C>::self self; inline self operator()(      Matrix<shift_dense_array_generator<      T,C2> > &X, const typename Matrix<shift_dense_array_generator<      T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,const Matrix<shift_dense_array_generator<      T,C2> >,C> { typedef typename TinySubMatrix<M,N,const T,C>::self self; inline self operator()(const Matrix<shift_dense_array_generator<      T,C2> > &X, const typename Matrix<shift_dense_array_generator<      T,C2> >::index_type &p) const { return self(X.generator().array(),p-X.generator().index()); } };

template<int M,int N,class G> inline typename TinySubMatrix<M,N,      Matrix<G> >::self sub(      Matrix<G> &X, const typename Matrix<G>::index_type &p) { return TinySubMatrix<M,N,      Matrix<G> >()(X,p); }
template<int M,int N,class G> inline typename TinySubMatrix<M,N,const Matrix<G> >::self sub(const Matrix<G> &X, const typename Matrix<G>::index_type &p) { return TinySubMatrix<M,N,const Matrix<G> >()(X,p); }
template<int M,int N,class G> inline typename TinySubMatrix<M,N,      Matrix<G> >::self sub(      Matrix<G> &X, const typename Matrix<G>::int_type &i, const typename Matrix<G>::int_type &j) { return sub<M,N>(X,typename Matrix<G>::index_type(i,j)); }
template<int M,int N,class G> inline typename TinySubMatrix<M,N,const Matrix<G> >::self sub(const Matrix<G> &X, const typename Matrix<G>::int_type &i, const typename Matrix<G>::int_type &j) { return sub<M,N>(X,typename Matrix<G>::index_type(i,j)); }

template<int M,int N,class T,int C=0>
class tiny_sub_dense_matrix_generator : public tiny_sub_matrix_generator<M,N,T,C>
{
  public:
    typedef tiny_sub_dense_matrix_generator self;
    typedef tiny_sub_matrix_generator<M,N,T,C> base;
    ARRAY_BASE_TYPES    
        
    template<class V2> struct array_rebind { typedef Matrix<tiny_matrix_generator<M,N,V2> > other; };
        
  private:
    pointer data; 
  
  public:
    using base::array;
    inline tiny_sub_dense_matrix_generator(array_type &x, const index_type &p) : base(x,p), data(&x[p]) {}
         
    template<class T2> inline tiny_sub_dense_matrix_generator(const Matrix<tiny_sub_dense_matrix_generator<M,N,T2> > &x) : base((base &)(x.generator())), data(&x[0]) {}

    int_type row_stride() const { return array().generator().row_stride(); }
         
		inline reference       operator[](const index_type &i)       { return data[array().ncols()*i.i+i.j]; }
    inline const_reference operator[](const index_type &i) const { return data[array().ncols()*i.i+i.j]; }
};

template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,      Matrix<sub_dense_matrix_generator     <      T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,      T,C> >::self self; inline self operator()(      Matrix<sub_dense_matrix_generator     <      T,C2> > &X, const typename Matrix<sub_dense_matrix_generator     <      T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,              class T,int C2,int C> struct TinySubMatrix<M,N,const Matrix<sub_dense_matrix_generator     <      T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,const T,C> >::self self; inline self operator()(const Matrix<sub_dense_matrix_generator     <      T,C2> > &X, const typename Matrix<sub_dense_matrix_generator     <      T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,int M2,int N2,class T,int C2,int C> struct TinySubMatrix<M,N,      Matrix<tiny_sub_dense_matrix_generator<M2,N2,T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,      T,C> >::self self; inline self operator()(      Matrix<tiny_sub_dense_matrix_generator<M2,N2,T,C2> > &X, const typename Matrix<tiny_sub_dense_matrix_generator<M2,N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
template<int M,int N,int M2,int N2,class T,int C2,int C> struct TinySubMatrix<M,N,const Matrix<tiny_sub_dense_matrix_generator<M2,N2,T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,const T,C> >::self self; inline self operator()(const Matrix<tiny_sub_dense_matrix_generator<M2,N2,T,C2> > &X, const typename Matrix<tiny_sub_dense_matrix_generator<M2,N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };

#ifdef NDEBUG
template<int M,int N,class V,        int C>        struct TinySubMatrix<M,N,      Matrix<data_matrix_generator    <V  > >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,      Matrix<data_matrix_generator    <V  > >,C> >::self self; inline self operator()(      Matrix<data_matrix_generator    <V  > > &X, const typename Matrix<dense_matrix_generator   <V  > >::index_type &p) const { return self(X,p); } };
template<int M,int N,class V,        int C>        struct TinySubMatrix<M,N,const Matrix<data_matrix_generator    <V  > >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,const Matrix<data_matrix_generator    <V  > >,C> >::self self; inline self operator()(const Matrix<data_matrix_generator    <V  > > &X, const typename Matrix<dense_matrix_generator   <V  > >::index_type &p) const { return self(X,p); } };
template<int M,int N,class V,class A,int C>        struct TinySubMatrix<M,N,      Matrix<dense_matrix_generator   <V,A> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,      Matrix<dense_matrix_generator   <V,A> >,C> >::self self; inline self operator()(      Matrix<dense_matrix_generator   <V,A> > &X, const typename Matrix<dense_matrix_generator   <V,A> >::index_type &p) const { return self(X,p); } };
template<int M,int N,class V,class A,int C>        struct TinySubMatrix<M,N,const Matrix<dense_matrix_generator   <V,A> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,const Matrix<dense_matrix_generator   <V,A> >,C> >::self self; inline self operator()(const Matrix<dense_matrix_generator   <V,A> > &X, const typename Matrix<dense_matrix_generator   <V,A> >::index_type &p) const { return self(X,p); } };
template<int M,int N,int M2,int N2 ,class V,int C> struct TinySubMatrix<M,N,      Matrix<tiny_matrix_generator<M2,N2,V> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,      Matrix<tiny_matrix_generator<M2,N2,V> >,C> >::self self; inline self operator()(      Matrix<tiny_matrix_generator<M2,N2,V> > &X, const typename Matrix<tiny_matrix_generator<M2,N2,V> >::index_type &p) const { return self(X,p); } };
template<int M,int N,int M2,int N2 ,class V,int C> struct TinySubMatrix<M,N,const Matrix<tiny_matrix_generator<M2,N2,V> >,C> { typedef typename GeneratorArray<tiny_sub_dense_matrix_generator<M,N,const Matrix<tiny_matrix_generator<M2,N2,V> >,C> >::self self; inline self operator()(const Matrix<tiny_matrix_generator<M2,N2,V> > &X, const typename Matrix<tiny_matrix_generator<M2,N2,V> >::index_type &p) const { return self(X,p); } };
#endif



template<class T, int Copy=0>
class sub_stride_dense_vector_generator : public sub_dense_vector_generator<T,Copy>
{
  public:
    typedef sub_dense_vector_generator<T,Copy> base;
    typedef sub_stride_dense_vector_generator self;
    ARRAY_BASE_TYPES
        
  protected:
    using base::data;
    size_type dim; 
  
  public:
    inline sub_stride_dense_vector_generator(array_type &x, const index_type &p) : base(x,p), dim(x.generator().stride()) {}
    inline sub_stride_dense_vector_generator(array_type &x, const index_type &p, const size_type &s) : base(x,p,s), dim(x.generator().stride()) {}
         
    template<class T2> inline sub_stride_dense_vector_generator(const Vector<sub_stride_dense_vector_generator<T2> > &x) : base((base &)(x.generator())), dim(x.generator().stride()) {}
    
    inline size_type stride() const { return dim; }
         
    inline reference       operator[](const index_type &i)       { return data[i*stride()]; }
    inline const_reference operator[](const index_type &i) const { return data[i*stride()]; }
};

#ifdef NDEBUG
template<class T,int C> class col_dense_matrix_generator;
template<class T,int C1,int C> struct SubArray<      Vector<col_dense_matrix_generator<T,C1> >,C> { typedef typename GeneratorArray<sub_stride_dense_vector_generator<      Vector<col_dense_matrix_generator<T,C1> >,C> >::self self; };
template<class T,int C1,int C> struct SubArray<const Vector<col_dense_matrix_generator<T,C1> >,C> { typedef typename GeneratorArray<sub_stride_dense_vector_generator<const Vector<col_dense_matrix_generator<T,C1> >,C> >::self self; };
#endif



template<int N,class T,int Copy=0>
class tiny_sub_stride_dense_vector_generator : public tiny_sub_dense_vector_generator<N,T,Copy>
{
  public:
    typedef tiny_sub_dense_vector_generator<N,T,Copy> base;
    typedef tiny_sub_stride_dense_vector_generator self;
    ARRAY_BASE_TYPES
        
    template<class V2> struct array_rebind { typedef typename TinyVector<N,V2>::self other; };
        
  protected:
    using base::data;
    size_type dim; 
  
  public:
    inline tiny_sub_stride_dense_vector_generator(array_type &x, const index_type &p) : base(x,p), dim(x.generator().stride()) {}
    template<class T2,int C2> inline tiny_sub_stride_dense_vector_generator(const Vector<tiny_sub_stride_dense_vector_generator<N,T2,C2> > &x) : base((base &)(x.generator())), dim(x.generator().stride()) {}

    inline size_type stride() const { return dim; }
         
    inline reference       operator[](const index_type &i)       { return data[i*stride()]; }
    inline const_reference operator[](const index_type &i) const { return data[i*stride()]; }
};

template<class T,int C,class G2> inline void affect(Vector<tiny_sub_stride_dense_vector_generator<1,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_stride_dense_vector_generator<2,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_stride_dense_vector_generator<3,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; }
template<class T,int C,class G2> inline void affect(Vector<tiny_sub_stride_dense_vector_generator<4,T,C> > &X, const Vector<G2>  &Y) { X[0]=Y[0]; X[1]=Y[1]; X[2]=Y[2]; X[3]=Y[3]; }

//template<int N,       class T,int C2,int C> struct TinySubVector<N,      Vector<sub_stride_dense_vector_generator     <   T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_stride_dense_vector_generator<N,      T,C> >::self self; inline self operator()(      Vector<sub_stride_dense_vector_generator     <   T,C2> > &X, const typename Vector<sub_stride_dense_vector_generator     <   T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
//template<int N,       class T,int C2,int C> struct TinySubVector<N,const Vector<sub_stride_dense_vector_generator     <   T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_stride_dense_vector_generator<N,const T,C> >::self self; inline self operator()(const Vector<sub_stride_dense_vector_generator     <   T,C2> > &X, const typename Vector<sub_stride_dense_vector_generator     <   T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
//template<int N,int N2,class T,int C2,int C> struct TinySubVector<N,      Vector<tiny_sub_stride_dense_vector_generator<N2,T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_stride_dense_vector_generator<N,      T,C> >::self self; inline self operator()(      Vector<tiny_sub_stride_dense_vector_generator<N2,T,C2> > &X, const typename Vector<tiny_sub_stride_dense_vector_generator<N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };
//template<int N,int N2,class T,int C2,int C> struct TinySubVector<N,const Vector<tiny_sub_stride_dense_vector_generator<N2,T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_stride_dense_vector_generator<N,const T,C> >::self self; inline self operator()(const Vector<tiny_sub_stride_dense_vector_generator<N2,T,C2> > &X, const typename Vector<tiny_sub_stride_dense_vector_generator<N2,T,C2> >::index_type &p) const { return self(X.generator().array(),X.generator().pos()+p); } };

#ifdef NDEBUG
template<class T,int C> class stride_dense_vector_generator;
template<int N,class T,int C2,int C> struct TinySubVector<N,      Vector<stride_dense_vector_generator  <T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_stride_dense_vector_generator<N,      Vector<stride_dense_vector_generator  <T,C2> >,C> >::self self; inline self operator()(      Vector<stride_dense_vector_generator  <T,C2> > &X, const typename Vector<stride_dense_vector_generator  <T,C2> >::index_type &i) const { return self(X,i); } };
template<int N,class T,int C2,int C> struct TinySubVector<N,const Vector<stride_dense_vector_generator  <T,C2> >,C> { typedef typename GeneratorArray<tiny_sub_stride_dense_vector_generator<N,const Vector<stride_dense_vector_generator  <T,C2> >,C> >::self self; inline self operator()(const Vector<stride_dense_vector_generator  <T,C2> > &X, const typename Vector<stride_dense_vector_generator  <T,C2> >::index_type &i) const { return self(X,i); } };
#endif




template<class T, int C=0>
class stride_array_generator : public array_generator<T,C>
{
  public:
    typedef array_generator<T,C> base;
    typedef stride_array_generator self;
    ARRAY_BASE_TYPES

  private:
    size_type m;
    index_type ind;

  public:
    inline stride_array_generator(array_type &x, const size_type &n) : base(x), m(n), ind(x.lower_bound()) {}
    inline stride_array_generator(array_type &x, const size_type &n, const index_type &i) : base(x), m(n), ind(i) {}
    inline stride_array_generator(const self &r) : base(static_cast<const base &>(r)), m(r.m), ind(r.ind) {}
    template<class T2,int C2> inline stride_array_generator(const Vector<sub_array_generator<T2,C2> > &x) : base(x.generator().array()), m(x.generator().stride()), ind(x.generator().pos()) {}

    using base::array;

    inline index_type lower_bound() const { return array().lower_bound()/stride(); }
    inline size_type  size       () const { return array().size()/stride(); }
    inline index_type stride     () const { return m; }
    inline index_type pos        () const { return ind; }
     
    inline reference       operator[](const index_type &i)       { return array()[i*stride()+pos()]; }
    inline const_reference operator[](const index_type &i) const { return array()[i*stride()+pos()]; }
    
    inline void resize(const size_type &s) { assert(size()==s); }
};

template<class T,int C=0> struct StrideArray { typedef typename GeneratorArray<stride_array_generator<T,C> >::self self; };

template<class G> inline typename StrideArray<      Vector<G> >::self stride(      Vector<G> &X, const typename Vector<G>::size_type &n                                         ) { return typename StrideArray<      Vector<G> >::self(X,n  ); }
template<class G> inline typename StrideArray<const Vector<G> >::self stride(const Vector<G> &X, const typename Vector<G>::size_type &n                                         ) { return typename StrideArray<const Vector<G> >::self(X,n  ); }
template<class G> inline typename StrideArray<      Vector<G> >::self stride(      Vector<G> &X, const typename Vector<G>::size_type &n, const typename Vector<G>::index_type &i) { return typename StrideArray<      Vector<G> >::self(X,n,i); }
template<class G> inline typename StrideArray<const Vector<G> >::self stride(const Vector<G> &X, const typename Vector<G>::size_type &n, const typename Vector<G>::index_type &i) { return typename StrideArray<const Vector<G> >::self(X,n,i); }
template<class G> inline typename StrideArray<      Matrix<G> >::self stride(      Matrix<G> &X, const typename Matrix<G>::size_type &n                                         ) { return typename StrideArray<      Matrix<G> >::self(X,n  ); }
template<class G> inline typename StrideArray<const Matrix<G> >::self stride(const Matrix<G> &X, const typename Matrix<G>::size_type &n                                         ) { return typename StrideArray<const Matrix<G> >::self(X,n  ); }
template<class G> inline typename StrideArray<      Matrix<G> >::self stride(      Matrix<G> &X, const typename Matrix<G>::size_type &n, const typename Matrix<G>::index_type &i) { return typename StrideArray<      Matrix<G> >::self(X,n,i); }
template<class G> inline typename StrideArray<const Matrix<G> >::self stride(const Matrix<G> &X, const typename Matrix<G>::size_type &n, const typename Matrix<G>::index_type &i) { return typename StrideArray<const Matrix<G> >::self(X,n,i); }


template<      class T,int C> struct StrideArray<      Vector< stride_array_generator        <T,C> > > { typedef typename StrideArray<      T,C>::self self; };
template<      class T,int C> struct StrideArray<const Vector< stride_array_generator        <T,C> > > { typedef typename StrideArray<const T,C>::self self; };

template<class T,int C> inline typename StrideArray<      Vector<stride_array_generator  <T,C> > >::self stride(      Vector<stride_array_generator  <T,C> > &X, const typename Vector<stride_array_generator  <T,C> >::size_type &n, const typename Vector<stride_array_generator  <T,C> >::index_type &i) { return typename StrideArray<      Vector<stride_array_generator  <T,C> > >::self(X.generator().array(),n*X.generator().stride(),X.generator().pos()+i*X.generator().stride()); }
template<class T,int C> inline typename StrideArray<const Vector<stride_array_generator  <T,C> > >::self stride(const Vector<stride_array_generator  <T,C> > &X, const typename Vector<stride_array_generator  <T,C> >::size_type &n, const typename Vector<stride_array_generator  <T,C> >::index_type &i) { return typename StrideArray<const Vector<stride_array_generator  <T,C> > >::self(X.generator().array(),n*X.generator().stride(),X.generator().pos()+i*X.generator().stride()); }
template<class T,int C> inline typename StrideArray<      Matrix<stride_array_generator  <T,C> > >::self stride(      Matrix<stride_array_generator  <T,C> > &X, const typename Matrix<stride_array_generator  <T,C> >::size_type &n, const typename Vector<stride_array_generator  <T,C> >::index_type &i) { return typename StrideArray<      Matrix<stride_array_generator  <T,C> > >::self(X.generator().array(),n*X.generator().stride(),X.generator().pos()+i*X.generator().stride()); }
template<class T,int C> inline typename StrideArray<const Matrix<stride_array_generator  <T,C> > >::self stride(const Matrix<stride_array_generator  <T,C> > &X, const typename Matrix<stride_array_generator  <T,C> >::size_type &n, const typename Vector<stride_array_generator  <T,C> >::index_type &i) { return typename StrideArray<const Matrix<stride_array_generator  <T,C> > >::self(X.generator().array(),n*X.generator().stride(),X.generator().pos()+i*X.generator().stride()); }



template<class T,int Copy=0>
class stride_dense_vector_generator : public stride_array_generator<T,Copy>
{
  public:
    typedef stride_array_generator<T,Copy> base;
    typedef stride_dense_vector_generator self;
    ARRAY_BASE_TYPES
        
  protected:
    pointer data; 
  
  public:
    inline stride_dense_vector_generator(array_type &x, const size_type &n                     ) : base(x,n,x.lower_bound()), data(&*x.begin()) {}
    inline stride_dense_vector_generator(array_type &x, const size_type &n, const index_type &i) : base(x,n,i), data(&x[i]) {}
         
    template<class T2> inline stride_dense_vector_generator(const Vector<stride_dense_vector_generator<T2> > &x) : base((base &)(x.generator())), data(&x[0]) {}
       
    using base::stride;
      
    inline reference       operator[](const index_type &i)       { return data[i*stride()]; }
    inline const_reference operator[](const index_type &i) const { return data[i*stride()]; }
};


#ifdef NDEBUG
template<        class V,int C> struct StrideArray<      Vector<data_vector_generator     <V  > >,C > { typedef typename GeneratorArray<stride_dense_vector_generator<      Vector<data_vector_generator     <V  > >,C > >::self self; };
template<        class V,int C> struct StrideArray<const Vector<data_vector_generator     <V  > >,C > { typedef typename GeneratorArray<stride_dense_vector_generator<const Vector<data_vector_generator     <V  > >,C > >::self self; };
template<class V,class A,int C> struct StrideArray<      Vector<dense_vector_generator    <V,A> >,C > { typedef typename GeneratorArray<stride_dense_vector_generator<      Vector<dense_vector_generator    <V,A> >,C > >::self self; };
template<class V,class A,int C> struct StrideArray<const Vector<dense_vector_generator    <V,A> >,C > { typedef typename GeneratorArray<stride_dense_vector_generator<const Vector<dense_vector_generator    <V,A> >,C > >::self self; };
template<int   N,class V,int C> struct StrideArray<      Vector<tiny_vector_generator     <N,V> >,C > { typedef typename GeneratorArray<stride_dense_vector_generator<      Vector<tiny_vector_generator     <N,V> >,C > >::self self; };
template<int   N,class V,int C> struct StrideArray<const Vector<tiny_vector_generator     <N,V> >,C > { typedef typename GeneratorArray<stride_dense_vector_generator<const Vector<tiny_vector_generator     <N,V> >,C > >::self self; };
template<class T,int C, int C1> struct StrideArray<      Vector<row_dense_matrix_generator<T,C> >,C1> { typedef typename GeneratorArray<stride_dense_vector_generator<      Vector<row_dense_matrix_generator<T,C> >,C1> >::self self; };
template<class T,int C, int C1> struct StrideArray<const Vector<row_dense_matrix_generator<T,C> >,C1> { typedef typename GeneratorArray<stride_dense_vector_generator<const Vector<row_dense_matrix_generator<T,C> >,C1> >::self self; };
#endif

template<class T,int C> struct StrideArray<      Vector<stride_dense_vector_generator<T,C> > > { typedef typename StrideArray<      T,C>::self self; };
template<class T,int C> struct StrideArray<const Vector<stride_dense_vector_generator<T,C> > > { typedef typename StrideArray<const T,C>::self self; };
template <class T,int C> class col_dense_matrix_generator;
template<class T,int C, int C1> struct StrideArray<      Vector<col_dense_matrix_generator<T,C> >,C1> { typedef typename GeneratorArray<stride_dense_vector_generator<      Vector<col_dense_matrix_generator<T,C> >,C1> >::self self; };
template<class T,int C, int C1> struct StrideArray<const Vector<col_dense_matrix_generator<T,C> >,C1> { typedef typename GeneratorArray<stride_dense_vector_generator<const Vector<col_dense_matrix_generator<T,C> >,C1> >::self self; };

template<class T,int C> inline typename StrideArray<      Vector<stride_dense_vector_generator<T,C> > >::self stride(      Vector<stride_dense_vector_generator<T,C> > &X, const typename Vector<stride_dense_vector_generator<T,C> >::size_type &n, const typename Vector<stride_dense_vector_generator<T,C> >::index_type &i) { return typename StrideArray<      Vector<stride_dense_vector_generator<T,C> > >::self(X.generator().array(),n*X.generator().stride(),X.generator().pos()+i*X.generator().stride()); }
template<class T,int C> inline typename StrideArray<const Vector<stride_dense_vector_generator<T,C> > >::self stride(const Vector<stride_dense_vector_generator<T,C> > &X, const typename Vector<stride_dense_vector_generator<T,C> >::size_type &n, const typename Vector<stride_dense_vector_generator<T,C> >::index_type &i) { return typename StrideArray<const Vector<stride_dense_vector_generator<T,C> > >::self(X.generator().array(),n*X.generator().stride(),X.generator().pos()+i*X.generator().stride()); }
template<class T,int C> inline typename StrideArray<      Vector<col_dense_matrix_generator   <T,C> > >::self stride(      Vector<col_dense_matrix_generator   <T,C> > &X, const typename Vector<col_dense_matrix_generator   <T,C> >::size_type &n, const typename Vector<col_dense_matrix_generator   <T,C> >::index_type &i) { return typename StrideArray<      Vector<col_dense_matrix_generator   <T,C> > >::self(X,n*X.generator().stride(),i); }
template<class T,int C> inline typename StrideArray<const Vector<col_dense_matrix_generator   <T,C> > >::self stride(const Vector<col_dense_matrix_generator   <T,C> > &X, const typename Vector<col_dense_matrix_generator   <T,C> >::size_type &n, const typename Vector<col_dense_matrix_generator   <T,C> >::index_type &i) { return typename StrideArray<const Vector<col_dense_matrix_generator   <T,C> > >::self(X,n*X.generator().stride(),i); }



template<class A> struct shift_iterator_rebind;
template<class A,int Copy> class shift_dense_array_generator;

template<class A,int C=0> struct shiftArray;

template<class T,int Copy=0>
class shift_array_generator : public array_generator<T,Copy>
{
  public:
    typedef array_generator<T,Copy> base;  
    typedef shift_array_generator self;
    ARRAY_BASE_TYPES
  
    template<class A2> struct iterator_rebind       { typedef typename shift_iterator_rebind<      A2>::other other; };
    template<class A2> struct const_iterator_rebind { typedef typename shift_iterator_rebind<const A2>::other other; };
    
    template<class V2> struct array_rebind { typedef typename shiftArray<typename array_type::template rebind<V2>::other,1>::self other; };
    
  protected:
    index_type ind;

  public:
    inline shift_array_generator(array_type &x, const index_type &i) : base(x), ind(i) {}
    template<class A2,int Copy2> inline shift_array_generator(const Vector<shift_array_generator       <A2,Copy2> > &X) : base(X.generator().array()), ind(X.generator().index()) { }
    template<class A2,int Copy2> inline shift_array_generator(const Vector<shift_dense_array_generator<A2,Copy2> > &X) : base(X.generator().array()), ind(X.generator().index()) { }

    template<class G> inline explicit shift_array_generator(      Vector<G>  &x) : base(x), ind(x.lower_bound()) {}
    template<class G> inline explicit shift_array_generator(const Vector<G>  &x) : base(x), ind(x.lower_bound()) {}
    template<class G> inline explicit shift_array_generator(      Matrix<G>  &x) : base(x), ind(x.lower_bound()) {}
    template<class G> inline explicit shift_array_generator(const Matrix<G>  &x) : base(x), ind(x.lower_bound()) {}
    template<int D,class G> inline explicit shift_array_generator(      Array<D,G> &x) : base(x), ind(x.lower_bound()) {}
    template<int D,class G> inline explicit shift_array_generator(const Array<D,G> &x) : base(x), ind(x.lower_bound()) {}

    inline shift_array_generator() : base(), ind(0) {}
    template<class A> inline shift_array_generator(const A &a) : base(a), ind(0) {}
    template<class A,class B> inline shift_array_generator(const A &a, const B &b) : base(a,b), ind(0) {}
    template<class A,class B,class C> inline shift_array_generator(const A &a, const B &b, const C &c) : base(a,b,c), ind(0) {}

    inline shift_array_generator(const index_type &i) : base(), ind(i) {}
    template<class A> inline shift_array_generator(const A &a, const index_type &i) : base(a), ind(i) {}
    template<class A,class B> inline shift_array_generator(const A &a, const B &b, const index_type &i) : base(a,b), ind(i) {}
    template<class A,class B,class C> inline shift_array_generator(const A &a, const B &b, const C &c, const index_type &i) : base(a,b,c), ind(i) {}

    inline self                         &operator=(const self       &x) { index()=x.lower_bound(); this->array()=no_shift(x); return *this; }
    template<class G>       inline self &operator=(const Vector<G>  &x) { index()=x.lower_bound(); this->array()=no_shift(x); return *this; }
    template<class G>       inline self &operator=(const Matrix<G>  &x) { index()=x.lower_bound(); this->array()=no_shift(x); return *this; }
    template<int D,class G> inline self &operator=(const Array<D,G> &x) { index()=x.lower_bound(); this->array()=no_shift(x); return *this; }
    
    inline void swap(self &g) { ::swap(this->array(),(&g)->array()); std::swap(ind,g.ind); }
    template<class G> inline void swap(G &g) { ::swap(this->array().generator(),g); }

    inline reference       operator[](const index_type &i)       { return this->array()[i-ind]; }
    inline const_reference operator[](const index_type &i) const { return this->array()[i-ind]; }
    
    inline index_type       &index()       { return ind; }
    inline const index_type &index() const { return ind; }

    inline void set_lower_bound(const index_type &i) { index()=i; }

    inline index_type lower_bound() const { return ind+this->array().lower_bound(); }

    inline void resize(const size_type &s) { this->array().resize(s); }  
};

template<class G,int C> inline void swap(shift_array_generator<Vector<G>,C> &x, G &y) { x.swap(y); }
template<class G,int C> inline void swap(G &y, shift_array_generator<Vector<G>,C> &x) { x.swap(y); }


template<      class A,int C> struct shiftArray { typedef typename GeneratorArray<shift_array_generator<A,C> >::self self; };
template<      class A,int C> struct shiftArray<      Vector< shift_array_generator      <A,C> > > { typedef typename shiftArray<      A,C>::self self; };
template<      class A,int C> struct shiftArray<const Vector< shift_array_generator      <A,C> > > { typedef typename shiftArray<const A,C>::self self; };
template<      class A,int C> struct shiftArray<      Matrix< shift_array_generator      <A,C> > > { typedef typename shiftArray<      A,C>::self self; };
template<      class A,int C> struct shiftArray<const Matrix< shift_array_generator      <A,C> > > { typedef typename shiftArray<const A,C>::self self; };
template<int D,class A,int C> struct shiftArray<      Array<D,shift_array_generator      <A,C> > > { typedef typename shiftArray<      A,C>::self self; };
template<int D,class A,int C> struct shiftArray<const Array<D,shift_array_generator      <A,C> > > { typedef typename shiftArray<const A,C>::self self; };

//{unsecret}
//Summary: Shifts the indexes
//Arguments:
//  X - The array to shift
//  p - The offset to use
//  i - The vertical offset to use
//  j - The horizontal offset to use
//Return: An array representing the shifted array
//Example:
//  DenseVector<int>::self X(4, "0 1 2 3");
//  cout << shift(X,-1)[-1] << endl; // 0
//  cout << shift(X,-1)[ 0] << endl; // 1
//  DenseMatrix<int>::self X(4,4, "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15");
//  cout << shift(X,-1,-1)(-1,-1) << endl; // 0
//  cout << shift(X,-1,-1)( 0, 0) << endl; // 5
//See: ^middle_shift^, ^no_shift^
template<class G> inline typename shiftArray<      Vector<G> >::self shift(      Vector<G> &X, const typename Vector<G>::index_type &p) { return typename shiftArray<      Vector<G> >::self(X,p); }
//{unsecret}
template<class G> inline typename shiftArray<const Vector<G> >::self shift(const Vector<G> &X, const typename Vector<G>::index_type &p) { return typename shiftArray<const Vector<G> >::self(X,p); }
//{unsecret}
template<class G> inline typename shiftArray<      Matrix<G> >::self shift(      Matrix<G> &X, const typename Matrix<G>::index_type &p) { return typename shiftArray<      Matrix<G> >::self(X,p); }
//{unsecret}
template<class G> inline typename shiftArray<const Matrix<G> >::self shift(const Matrix<G> &X, const typename Matrix<G>::index_type &p) { return typename shiftArray<const Matrix<G> >::self(X,p); }
//{unsecret}
template<class G> inline typename shiftArray<      Matrix<G> >::self shift(      Matrix<G> &X, const typename Matrix<G>::int_type   &i) { return typename shiftArray<      Matrix<G> >::self(X,typename Matrix<G>::index_type(i,i)); }
//{unsecret}
template<class G> inline typename shiftArray<const Matrix<G> >::self shift(const Matrix<G> &X, const typename Matrix<G>::int_type   &i) { return typename shiftArray<const Matrix<G> >::self(X,typename Matrix<G>::index_type(i,i)); }
//{unsecret}
template<class G> inline typename shiftArray<      Matrix<G> >::self shift(      Matrix<G> &X, const typename Matrix<G>::int_type   &i, const typename Matrix<G>::int_type   &j) { return typename shiftArray<      Matrix<G> >::self(X,typename Matrix<G>::index_type(i,j)); }
//{unsecret}
template<class G> inline typename shiftArray<const Matrix<G> >::self shift(const Matrix<G> &X, const typename Matrix<G>::int_type   &i, const typename Matrix<G>::int_type   &j) { return typename shiftArray<const Matrix<G> >::self(X,typename Matrix<G>::index_type(i,j)); }

template<class A,int C> inline typename shiftArray<      A>::self shift(      Vector<shift_array_generator      <A,C> > &X, const typename Vector<shift_array_generator      <A,C> >::index_type &i) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Vector<shift_array_generator      <A,C> > &X, const typename Vector<shift_array_generator      <A,C> >::index_type &i) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<      A>::self shift(      Matrix<shift_array_generator      <A,C> > &X, const typename Matrix<shift_array_generator      <A,C> >::index_type &i) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Matrix<shift_array_generator      <A,C> > &X, const typename Matrix<shift_array_generator      <A,C> >::index_type &i) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<      A>::self shift(      Matrix<shift_array_generator      <A,C> > &X, const typename Matrix<shift_array_generator      <A,C> >::int_type   &i) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Matrix<shift_array_generator      <A,C> > &X, const typename Matrix<shift_array_generator      <A,C> >::int_type   &i) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<      A>::self shift(      Matrix<shift_array_generator      <A,C> > &X, const typename Matrix<shift_array_generator      <A,C> >::int_type   &i, const typename Matrix<shift_array_generator      <A,C> >::int_type   &j) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+typename Matrix<shift_array_generator<A,C> >::index_type(i,j)); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Matrix<shift_array_generator      <A,C> > &X, const typename Matrix<shift_array_generator      <A,C> >::int_type   &i, const typename Matrix<shift_array_generator      <A,C> >::int_type   &j) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+typename Matrix<shift_array_generator<A,C> >::index_type(i,j)); }

//{unsecret}
//Summary: Shifts the indexes of the half size
//Arguments: X - The array to shift
//Return: An array representing the shifted array
//Example:
//  DenseVector<int>::self X(4, "0 1 2 3");
//  cout << middle_shift(X)[-1] << endl; // 1
//  cout << middle_shift(X)[ 0] << endl; // 2
//  DenseMatrix<int>::self X(4,4, "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15");
//  cout << middle_shift(X)(-1,-1) << endl; // 5
//  cout << middle_shift(X)( 0, 0) << endl; // 10
//See: ^shift^, ^no_shift^
template<class G> inline typename shiftArray<      Vector<G> >::self middle_shift(      Vector<G> &X) { return shift(X,-X.lower_bound()-X.size()/2); }
//{unsecret}
template<class G> inline typename shiftArray<const Vector<G> >::self middle_shift(const Vector<G> &X) { return shift(X,-X.lower_bound()-X.size()/2); }
//{unsecret}
template<class G> inline typename shiftArray<      Matrix<G> >::self middle_shift(      Matrix<G> &X) { return shift(X,-X.lower_bound()-X.size()/2); }
//{unsecret}
template<class G> inline typename shiftArray<const Matrix<G> >::self middle_shift(const Matrix<G> &X) { return shift(X,-X.lower_bound()-X.size()/2); }

//{unsecret}
//Summary: Suppresses the offset on a shifted array
//Arguments: X - The array to shift
//Return: An array representing the shifted array
//Example:
//  shiftDenseVector<int>::self X(4, "0 1 2 3", -2);
//  cout << X[0] << endl; // 2
//  cout << no_shift(X)[0] << endl; // 0
//See: ^shift^, ^middle_shift^
template<class G> inline typename shiftArray<      Vector<G> >::self no_shift(      Vector<G> &X) { return shift(X,-X.lower_bound()); }
//{unsecret}
template<class G> inline typename shiftArray<const Vector<G> >::self no_shift(const Vector<G> &X) { return shift(X,-X.lower_bound()); }
//{unsecret}
template<class G> inline typename shiftArray<      Matrix<G> >::self no_shift(      Matrix<G> &X) { return shift(X,-X.lower_bound()); }
//{unsecret}
template<class G> inline typename shiftArray<const Matrix<G> >::self no_shift(const Matrix<G> &X) { return shift(X,-X.lower_bound()); }

template<class V,class A,int C> inline       Vector<dense_vector_generator<V,A> > &no_shift(      Vector<shift_array_generator      <      Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }
template<class V,class A,int C> inline const Vector<dense_vector_generator<V,A> > &no_shift(      Vector<shift_array_generator      <const Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }
template<class V,class A,int C> inline const Vector<dense_vector_generator<V,A> > &no_shift(const Vector<shift_array_generator      <      Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }
template<class V,class A,int C> inline const Vector<dense_vector_generator<V,A> > &no_shift(const Vector<shift_array_generator      <const Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }

template<class V,class A,int C> inline       Vector<dense_vector_generator<V,A> > &no_shift(      Vector<shift_dense_array_generator<      Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }
template<class V,class A,int C> inline const Vector<dense_vector_generator<V,A> > &no_shift(      Vector<shift_dense_array_generator<const Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }
template<class V,class A,int C> inline const Vector<dense_vector_generator<V,A> > &no_shift(const Vector<shift_dense_array_generator<      Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }
template<class V,class A,int C> inline const Vector<dense_vector_generator<V,A> > &no_shift(const Vector<shift_dense_array_generator<const Vector<dense_vector_generator<V,A> >,C> > &X) { return X.generator().array(); }


template<class T,int Copy=0>
class shift_dense_array_generator : public shift_array_generator<T,Copy>
{
  public:
    typedef shift_dense_array_generator self;
    typedef shift_array_generator<T,Copy> base;
    ARRAY_BASE_TYPES    
        
  private:
    pointer data; 
  
  public:
    shift_dense_array_generator(array_type &x, const index_type &i) : base(x,i), data(&x[-i]) {}
    template<      class A2,int Copy2> shift_dense_array_generator(const Vector< shift_array_generator      <A2,Copy2> > &x) : base(x.generator().array(),x.generator().index()) { data = &(this->array())[-this->index()]; }
    template<      class A2,int Copy2> shift_dense_array_generator(const Vector< shift_dense_array_generator<A2,Copy2> > &x) : base(x.generator().array(),x.generator().index()) { data = &(this->array())[-this->index()]; }
    template<      class A2,int Copy2> shift_dense_array_generator(const Matrix< shift_array_generator      <A2,Copy2> > &x) : base(x.generator().array(),x.generator().index()) { data = &(this->array())[-this->index()]; }
    template<      class A2,int Copy2> shift_dense_array_generator(const Matrix< shift_dense_array_generator<A2,Copy2> > &x) : base(x.generator().array(),x.generator().index()) { data = &(this->array())[-this->index()]; }
    template<int D,class A2,int Copy2> shift_dense_array_generator(const Array<D,shift_array_generator      <A2,Copy2> > &x) : base(x.generator().array(),x.generator().index()) { data = &(this->array())[-this->index()]; }
    template<int D,class A2,int Copy2> shift_dense_array_generator(const Array<D,shift_dense_array_generator<A2,Copy2> > &x) : base(x.generator().array(),x.generator().index()) { data = &(this->array())[-this->index()]; }

    shift_dense_array_generator() : base() { data=&(this->array())[-this->index()]; }
    template<class A> shift_dense_array_generator(A &a) : base(a) { data=&(this->array())[-this->index()]; }
    template<class A,class B> shift_dense_array_generator(A &a, B &b) : base(a,b) { data=&(this->array())[-this->index()]; }
    template<class A,class B,class C> shift_dense_array_generator(A &a, B &b, C &c) : base(a,b,c) { data=&(this->array())[-this->index()]; }

    shift_dense_array_generator(const index_type &i) : base(i) { data=&(this->array())[-this->index()]; }
    template<class A> shift_dense_array_generator(A &a, const index_type &i) : base(a,i) { data=&(this->array())[-this->index()]; }
    template<class A,class B> shift_dense_array_generator(A &a, B &b, const index_type &i) : base(a,b,i) { data=&(this->array())[-this->index()]; }
    template<class A,class B,class C> shift_dense_array_generator(A &a, B &b, C &c, const index_type &i) : base(a,b,c,i) { data=&(this->array())[-this->index()]; }

    self                         &operator=(const self       &x) { static_cast<base &>(*this)=x; data=&(this->array())[-this->index()]; return *this; }
    template<      class G> self &operator=(const Vector<G>  &x) { static_cast<base &>(*this)=x; data=&(this->array())[-this->index()]; return *this; }
    template<      class G> self &operator=(const Matrix<G>  &x) { static_cast<base &>(*this)=x; data=&(this->array())[-this->index()]; return *this; }
    template<int D,class G> self &operator=(const Array<D,G> &x) { static_cast<base &>(*this)=x; data=&(this->array())[-this->index()]; return *this; }

    void swap(self &g) { std::swap(data,g.data); base::swap(g); }
    template<class G> void swap(G &g) { base::swap(g); }

    void set_lower_bound(const index_type &i) { this->index()=i; data=&(this->array())[-this->index()]; }

    reference       operator[](const int &i)       { return data[i]; }
    const_reference operator[](const int &i) const { return data[i]; }
    template<class I> reference       operator[](const MatrixIndex<I> &p)       { return data[p.i*this->size().ncols()+p.j]; }
    template<class I> const_reference operator[](const MatrixIndex<I> &p) const { return data[p.i*this->size().ncols()+p.j]; }

    void resize(const size_type &s) { (this->array()).resize(s); data=&(this->array())[-this->index()]; }  
};

template<class G,int C> void swap(shift_dense_array_generator<Vector<G>,C> &x, G &y) { x.swap(y); }
template<class G,int C> void swap(G &y, shift_dense_array_generator<Vector<G>,C> &x) { x.swap(y); }


template<      class A,int C> struct shiftArray<      Vector< shift_dense_array_generator      <A,C> > > { typedef typename shiftArray<      A,C>::self self; };
template<      class A,int C> struct shiftArray<const Vector< shift_dense_array_generator      <A,C> > > { typedef typename shiftArray<const A,C>::self self; };
template<      class A,int C> struct shiftArray<      Matrix< shift_dense_array_generator      <A,C> > > { typedef typename shiftArray<      A,C>::self self; };
template<      class A,int C> struct shiftArray<const Matrix< shift_dense_array_generator      <A,C> > > { typedef typename shiftArray<const A,C>::self self; };
template<int D,class A,int C> struct shiftArray<      Array<D,shift_dense_array_generator      <A,C> > > { typedef typename shiftArray<      A,C>::self self; };
template<int D,class A,int C> struct shiftArray<const Array<D,shift_dense_array_generator      <A,C> > > { typedef typename shiftArray<const A,C>::self self; };

#ifdef NDEBUG
template<    class V        ,int C> struct shiftArray<      Vector<data_vector_generator  <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Vector<data_vector_generator  <V    > >,C> >::self self; };
template<    class V        ,int C> struct shiftArray<const Vector<data_vector_generator  <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Vector<data_vector_generator  <V    > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<      Vector<dense_vector_generator <V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Vector<dense_vector_generator <V,A  > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<const Vector<dense_vector_generator <V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Vector<dense_vector_generator <V,A  > >,C> >::self self; };
template<      int N,class V,int C> struct shiftArray<      Vector<tiny_vector_generator  <N,V  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Vector< tiny_vector_generator <N,V  > >,C> >::self self; };
template<      int N,class V,int C> struct shiftArray<const Vector<tiny_vector_generator  <N,V  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Vector< tiny_vector_generator <N,V  > >,C> >::self self; };

template<    class V        ,int C> struct shiftArray<      Matrix<data_matrix_generator  <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Matrix<data_matrix_generator  <V    > >,C> >::self self; };
template<    class V        ,int C> struct shiftArray<const Matrix<data_matrix_generator  <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Matrix<data_matrix_generator  <V    > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<      Matrix<dense_matrix_generator <V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Matrix<dense_matrix_generator <V,A  > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<const Matrix<dense_matrix_generator <V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Matrix<dense_matrix_generator <V,A  > >,C> >::self self; };
template<int M,int N,class V,int C> struct shiftArray<      Matrix<tiny_matrix_generator  <M,N,V> >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Matrix< tiny_matrix_generator <M,N,V> >,C> >::self self; };
template<int M,int N,class V,int C> struct shiftArray<const Matrix<tiny_matrix_generator  <M,N,V> >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Matrix< tiny_matrix_generator <M,N,V> >,C> >::self self; };

template<    class V        ,int C> struct shiftArray<      Array<1,data_vector_generator <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Array<1,data_vector_generator <V    > >,C> >::self self; };
template<    class V        ,int C> struct shiftArray<const Array<1,data_vector_generator <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Array<1,data_vector_generator <V    > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<      Array<1,dense_vector_generator<V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Array<1,dense_vector_generator<V,A  > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<const Array<1,dense_vector_generator<V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Array<1,dense_vector_generator<V,A  > >,C> >::self self; };
template<      int N,class V,int C> struct shiftArray<      Array<1,tiny_vector_generator <N,V  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Array<1,tiny_vector_generator <N,V  > >,C> >::self self; };
template<      int N,class V,int C> struct shiftArray<const Array<1,tiny_vector_generator <N,V  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Array<1,tiny_vector_generator <N,V  > >,C> >::self self; };

template<    class V        ,int C> struct shiftArray<      Array<2,data_matrix_generator <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Array<2,data_matrix_generator <V    > >,C> >::self self; };
template<    class V        ,int C> struct shiftArray<const Array<2,data_matrix_generator <V    > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Array<2,data_matrix_generator <V    > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<      Array<2,dense_matrix_generator<V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Array<2,dense_matrix_generator<V,A  > >,C> >::self self; };
template<    class V,class A,int C> struct shiftArray<const Array<2,dense_matrix_generator<V,A  > >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Array<2,dense_matrix_generator<V,A  > >,C> >::self self; };
template<int M,int N,class V,int C> struct shiftArray<      Array<2,tiny_matrix_generator <M,N,V> >,C> { typedef typename GeneratorArray<shift_dense_array_generator<      Array<2,tiny_matrix_generator <M,N,V> >,C> >::self self; };
template<int M,int N,class V,int C> struct shiftArray<const Array<2,tiny_matrix_generator <M,N,V> >,C> { typedef typename GeneratorArray<shift_dense_array_generator<const Array<2,tiny_matrix_generator <M,N,V> >,C> >::self self; };
#endif

template<class A,int C> inline typename shiftArray<      A>::self shift(      Vector<shift_dense_array_generator<A,C> > &X, const typename Vector<shift_dense_array_generator<A,C> >::index_type &i) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Vector<shift_dense_array_generator<A,C> > &X, const typename Vector<shift_dense_array_generator<A,C> >::index_type &i) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<      A>::self shift(      Matrix<shift_dense_array_generator<A,C> > &X, const typename Matrix<shift_dense_array_generator<A,C> >::index_type &i) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Matrix<shift_dense_array_generator<A,C> > &X, const typename Matrix<shift_dense_array_generator<A,C> >::index_type &i) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<      A>::self shift(      Matrix<shift_dense_array_generator<A,C> > &X, const typename Matrix<shift_dense_array_generator<A,C> >::int_type   &i) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Matrix<shift_dense_array_generator<A,C> > &X, const typename Matrix<shift_dense_array_generator<A,C> >::int_type   &i) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+i); }
template<class A,int C> inline typename shiftArray<      A>::self shift(      Matrix<shift_dense_array_generator<A,C> > &X, const typename Matrix<shift_dense_array_generator<A,C> >::int_type   &i, const typename Matrix<shift_dense_array_generator<A,C> >::int_type   &j) { return typename shiftArray<      A>::self(X.generator().array(), X.generator().lower_bound()+typename Matrix<shift_dense_array_generator<A,C> >::index_type(i,j) ); }
template<class A,int C> inline typename shiftArray<const A>::self shift(const Matrix<shift_dense_array_generator<A,C> > &X, const typename Matrix<shift_dense_array_generator<A,C> >::int_type   &i, const typename Matrix<shift_dense_array_generator<A,C> >::int_type   &j) { return typename shiftArray<const A>::self(X.generator().array(), X.generator().lower_bound()+typename Matrix<shift_dense_array_generator<A,C> >::index_type(i,j) ); }

//{unsecret}
//{group:Vectors Interfaces}
//summary: Shifted dense vector type
//Arguments:
//  V - Type of the elements
//  A - Optional allocator
//example:
//  shiftDenseVector<int>::self X(4, "1 2 3 4", -2); // validity domain of indexes: -2..1
//  DenseVector<int>::self Y(4,"1 2 3 4");
//  shiftDenseVector<int>::self Z = shift(X,-2) // shift and copy, validity domain of Y indexes: -2..1
//See: ^DenseVector^, ^shiftDenseMatrix^, ^shiftTinyVector^, ^shift^
template<class V,class A=alignment_allocator<V> > struct shiftDenseVector { typedef typename shiftArray<typename DenseVector<V,A>::self,1>::self self; };

//{unsecret}
//{group:Matrices Interfaces}
//summary: Shifted dense matrix type
//Arguments:
//  V - Type of the elements
//  A - Optional allocator
//example:
//  shiftDenseMatrix<int>::self X(3,3,"1 2 3 4 5 6 7 8 9", MatrixIndex<int>(-1,-1)); // validity domain of indexes: -1..1, -1..1
//  DenseMatrix<int>::self Y(3,3,"1 2 3 4 5 6 7 8 9");
//  shiftDenseMatrix<int>::self Z = shift(X,-1,-1) // shift and copy, validity domain of Y indexes: -1..1, -1..1
//See: ^DenseMatrix^, ^shiftDenseVector^, ^shiftTinyMatrix^, ^shift^
template<class V,class A=alignment_allocator<V> > struct shiftDenseMatrix { typedef typename shiftArray<typename DenseMatrix<V,A>::self,1>::self self; };

//{unsecret}
//{group:Vectors Interfaces}
//summary:Shifted dense vector
//Arguments:
//  N - Size of the vector
//  V - Type of the elements
//See: ^TinyVector^, ^shiftDenseVector^, ^shiftTinyMatrix^, ^shift^
template<int N,class V> struct shiftTinyVector { typedef typename shiftArray<typename TinyVector<N,V>::self,1>::self self; };

//{unsecret}
//{group:Matrices Interfaces}
//summary: Shifted dense matrix type
//Arguments:
//  M - Number of rows
//  N - Number of columns
//  V - Type of the elements
//See: ^TinyMatrix^, ^shiftDenseMatrix^, ^shiftTinyVector^, ^shift^
template<int M,int N,class V> struct shiftTinyMatrix { typedef typename shiftArray<typename TinyMatrix<M,N,V>::self,1>::self self; };

//{unsecret}
//Summary: Copies an array in a temporary shifted dense array
template<class G> inline typename shiftDenseVector<typename Vector<G>::const_value_type>::self shift_dense(const Vector<G> &X) { typename shiftDenseVector<typename Vector<G>::const_value_type>::self Y; Y=X; return Y; }
//{unsecret}
template<class G> inline typename shiftDenseMatrix<typename Matrix<G>::const_value_type>::self shift_dense(const Matrix<G> &X) { typename shiftDenseMatrix<typename Matrix<G>::const_value_type>::self Y; Y=X; return Y; }

template<class T>
struct shift_dense_vector_function : public unary_value_function<T,typename T::template shift_dense_rebind<typename T::const_value_type>::other>
{
  typedef unary_value_function<T,typename T::template shift_dense_rebind<typename T::const_value_type>::other> base;
  UNARY_FUNCTION_BASE_TYPES
  const_reference operator()(argument_type &x) const { return shift_dense(x); }
};

template<class A,class F,int Copy=0> struct UnaryFunctionArray;

template<class G>
typename UnaryFunctionArray<const Vector<G>,const shift_dense_vector_function<const typename Vector<G>::value_type> >::self
ele_shift_dense(const Vector<G> &X)
{
  return apply(X,shift_dense_vector_function<const typename Vector<G>::value_type>());
}

template<class T>
struct shift_dense_matrix_function : public unary_value_function<T,typename T::template shift_dense_rebind<typename T::const_value_type>::other>
{
  typedef unary_value_function<T,typename T::template shift_dense_rebind<typename T::const_value_type>::other> base;
  UNARY_FUNCTION_BASE_TYPES
  const_reference operator()(argument_type &x) const { return shift_dense(x); }
};

template<class G>
typename UnaryFunctionArray<const Matrix<G>,const shift_dense_matrix_function<const typename Matrix<G>::value_type> >::self
ele_shift_dense(const Matrix<G> &X)
{
  return apply(X,shift_dense_matrix_function<const typename Matrix<G>::value_type>());
}



//template<class A> struct block_rebind;
//template<class V>          struct block_rebind<Vector<dense_vector_generator     <V     > > > { typedef typename DenseVector<typename DenseVector<V>::self>::self other; };
//template<class V>          struct block_rebind<Matrix<dense_matrix_generator     <V     > > > { typedef typename DenseMatrix<typename DenseMatrix<V>::self>::self other; };
//template<class T,int Copy> struct block_rebind<Vector<shift_array_generator      <T,Copy> > > { typedef typename shiftDenseVector<typename shiftDenseVector<typename T::const_value_type>::self>::self other; };
//template<class T,int Copy> struct block_rebind<Matrix<shift_array_generator      <T,Copy> > > { typedef typename shiftDenseMatrix<typename shiftDenseMatrix<typename T::const_value_type>::self>::self other; };
//template<class T,int Copy> struct block_rebind<Vector<shift_dense_array_generator<T,Copy> > > { typedef typename shiftDenseVector<typename shiftDenseVector<typename T::const_value_type>::self>::self other; };
//template<class T,int Copy> struct block_rebind<Matrix<shift_dense_array_generator<T,Copy> > > { typedef typename shiftDenseMatrix<typename shiftDenseMatrix<typename T::const_value_type>::self>::self other; };


template<class T,int C=0>
class block_array_generator : public array_reference_generator<T,typename SubArray<T>::self,typename SubArray<T>::self,typename SubArray<const T>::self,C>
{
  public:
    typedef block_array_generator self;
    typedef array_reference_generator<T,typename SubArray<T>::self,typename SubArray<T>::self,typename SubArray<const T>::self,C> base;
    ARRAY_BASE_TYPES

    //template<class V2> struct array_rebind { typedef typename block_rebind<typename array_type::rebind<V2>::other>::other other; };

  private:
    size_type dim;

  public:
    inline block_array_generator(array_type &x, const size_type &d) : base(x), dim(d) {}
    inline block_array_generator(const self &r) : base(r), dim(r.dim) {}
    template<class T2,int C2> inline block_array_generator(const Vector<block_array_generator<T2,C2> > &x) : base(x.generator().array()), dim(x.generator().stride()) {}    
    
    inline index_type lower_bound() const { return this->array().lower_bound()/dim; }
    inline size_type size() const { return this->array().size()/dim; }
    inline size_type stride() const { return dim; }

    inline reference       operator[](const index_type &p)       { return sub(this->array(),p*dim, dim); }
    inline const_reference operator[](const index_type &p) const { return sub(this->array(),p*dim, dim); }

    inline void resize(const size_type &d) { this->array().resize(d*dim); }
};

template<class T,int C=0>
struct blockArray
{
  typedef T array_type;
  typedef block_array_generator<array_type,C> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
//Summary: Decomposition in blocks
//Arguments:
//  X - The array to decompose
//  s - Size of each block
//  m - Height of each block
//  n - Width of each block
//Result: An array representing every blocks
//Example:
//  DenseVector<int>::self X(9, "0 1 2 3 4 5 6 7 8");
//  cout << block(X,3)[1] << endl; // [3 4 5]
//  DenseMatrix<int>::self Y(16, "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15");
//  cout << block(Y,2,2)(1,1) << endl; // [10 11; 14 15]
template<class G> inline typename blockArray<      Vector<G> >::self block(      Vector<G> &X, const typename Vector<G>::size_type &s) { return typename blockArray<      Vector<G> >::self(X,s); }
//{unsecret}
template<class G> inline typename blockArray<const Vector<G> >::self block(const Vector<G> &X, const typename Vector<G>::size_type &s) { return typename blockArray<const Vector<G> >::self(X,s); }
//{unsecret}
template<class G> inline typename blockArray<      Matrix<G> >::self block(      Matrix<G> &X, const typename Matrix<G>::size_type &s) { return typename blockArray<      Matrix<G> >::self(X,s); }
//{unsecret}
template<class G> inline typename blockArray<const Matrix<G> >::self block(const Matrix<G> &X, const typename Matrix<G>::size_type &s) { return typename blockArray<const Matrix<G> >::self(X,s); }
//{unsecret}
template<class G> inline typename blockArray<      Matrix<G> >::self block(      Matrix<G> &X, const typename Matrix<G>::size_type::size_type &m, const typename Matrix<G>::size_type::size_type &n) { return block(X,typename Matrix<G>::size_type(m,n)); }
//{unsecret}
template<class G> inline typename blockArray<const Matrix<G> >::self block(const Matrix<G> &X, const typename Matrix<G>::size_type::size_type &m, const typename Matrix<G>::size_type::size_type &n) { return block(X,typename Matrix<G>::size_type(m,n)); }
//{unsecret}
template<class G> inline typename blockArray<      Matrix<G> >::self block(      Matrix<G> &X, const typename Matrix<G>::size_type::size_type &n) { return block(X,n,n); }
//{unsecret}
template<class G> inline typename blockArray<const Matrix<G> >::self block(const Matrix<G> &X, const typename Matrix<G>::size_type::size_type &n) { return block(X,n,n); }


template<class T,int Copy=0>
class block_dense_vector_generator : public block_array_generator<T,Copy>
{
  public:
    typedef block_dense_vector_generator self;
    typedef block_array_generator<T,Copy> base;
    ARRAY_BASE_TYPES
          
  public:
    inline block_dense_vector_generator(array_type &x, const size_type &d) : base(x,d) {}
    inline block_dense_vector_generator(const self &r) : base(r) {}
         
    template<class T2> inline block_dense_vector_generator(const Vector<block_dense_vector_generator<T2> > &x) : base((base &)(x.generator())) {}
};

#ifdef NDEBUG
template<      class V,        int C> struct blockArray<      Vector<data_vector_generator          <V     > >,C> { typedef typename GeneratorArray<block_dense_vector_generator<      Vector<data_vector_generator          <V     > >,C> >::self self; };
template<      class V,        int C> struct blockArray<const Vector<data_vector_generator          <V     > >,C> { typedef typename GeneratorArray<block_dense_vector_generator<const Vector<data_vector_generator          <V     > >,C> >::self self; };
template<      class V,class A,int C> struct blockArray<      Vector<dense_vector_generator         <V,A   > >,C> { typedef typename GeneratorArray<block_dense_vector_generator<      Vector<dense_vector_generator         <V,A   > >,C> >::self self; };
template<      class V,class A,int C> struct blockArray<const Vector<dense_vector_generator         <V,A   > >,C> { typedef typename GeneratorArray<block_dense_vector_generator<const Vector<dense_vector_generator         <V,A   > >,C> >::self self; };
template<int N,class V,        int C> struct blockArray<      Vector<tiny_vector_generator          <N,V   > >,C> { typedef typename GeneratorArray<block_dense_vector_generator<      Vector<tiny_vector_generator          <N,V   > >,C> >::self self; };
template<int N,class V,        int C> struct blockArray<const Vector<tiny_vector_generator          <N,V   > >,C> { typedef typename GeneratorArray<block_dense_vector_generator<const Vector<tiny_vector_generator          <N,V   > >,C> >::self self; };
template<      class T,int C1, int C> struct blockArray<      Vector<sub_dense_vector_generator     <T  ,C1> >,C> { typedef typename GeneratorArray<block_dense_vector_generator<      Vector<sub_dense_vector_generator     <T  ,C1> >,C> >::self self; };
template<      class T,int C1, int C> struct blockArray<const Vector<sub_dense_vector_generator     <T  ,C1> >,C> { typedef typename GeneratorArray<block_dense_vector_generator<const Vector<sub_dense_vector_generator     <T  ,C1> >,C> >::self self; };
template<int N,class T,int C1, int C> struct blockArray<      Vector<tiny_sub_dense_vector_generator<N,T,C1> >,C> { typedef typename GeneratorArray<block_dense_vector_generator<      Vector<tiny_sub_dense_vector_generator<N,T,C1> >,C> >::self self; };
template<int N,class T,int C1, int C> struct blockArray<const Vector<tiny_sub_dense_vector_generator<N,T,C1> >,C> { typedef typename GeneratorArray<block_dense_vector_generator<const Vector<tiny_sub_dense_vector_generator<N,T,C1> >,C> >::self self; };
#endif


template<class T,int Copy=0>
class block_stride_dense_vector_generator : public block_array_generator<T,Copy>
{
  public:
    typedef block_stride_dense_vector_generator self;
    typedef block_array_generator<T,Copy> base;
    ARRAY_BASE_TYPES
          
  public:
    inline block_stride_dense_vector_generator(array_type &x, const size_type &d) : base(x,d) {}
    inline block_stride_dense_vector_generator(const self &r) : base(r) {}
         
    template<class T2> inline block_stride_dense_vector_generator(const Vector<block_stride_dense_vector_generator<T2> > &x) : base((base &)(x.generator())) {}
};

#ifdef NDEBUG
template<class T,int C> class col_dense_matrix_generator;
template<class T,int C1, int C> struct blockArray<      Vector<col_dense_matrix_generator       <T,C1> >,C> { typedef typename GeneratorArray<block_stride_dense_vector_generator<      Vector<col_dense_matrix_generator       <T,C1> >,C> >::self self; };
template<class T,int C1, int C> struct blockArray<const Vector<col_dense_matrix_generator       <T,C1> >,C> { typedef typename GeneratorArray<block_stride_dense_vector_generator<const Vector<col_dense_matrix_generator       <T,C1> >,C> >::self self; };
template<class T,int C1, int C> struct blockArray<      Vector<sub_stride_dense_vector_generator<T,C1> >,C> { typedef typename GeneratorArray<block_stride_dense_vector_generator<      Vector<sub_stride_dense_vector_generator<T,C1> >,C> >::self self; };
template<class T,int C1, int C> struct blockArray<const Vector<sub_stride_dense_vector_generator<T,C1> >,C> { typedef typename GeneratorArray<block_stride_dense_vector_generator<const Vector<sub_stride_dense_vector_generator<T,C1> >,C> >::self self; };
#endif



template<class T,int N,int Copy=0>
class tiny_block_vector_generator : public array_reference_generator<T,typename TinySubVector<N,T>::self,typename TinySubVector<N,T>::self,typename TinySubVector<N,const T>::self,Copy>
{
  public:
    typedef tiny_block_vector_generator self;
    typedef array_reference_generator<T,typename TinySubVector<N,T>::self,typename TinySubVector<N,T>::self,typename TinySubVector<N,const T>::self,Copy> base;
    ARRAY_BASE_TYPES

  public:
    inline tiny_block_vector_generator(array_type &x) : base(x) {}
    inline tiny_block_vector_generator(const self &r) : base(r) {}
    
    inline index_type lower_bound() const { return this->array().lower_bound()/N; }
    inline size_type size() const { return this->array().size()/N; }
    inline size_type stride() const { return N; }    

    inline reference       operator[](const index_type &p)       { return sub<N>(this->array(),p*N); }
    inline const_reference operator[](const index_type &p) const { return sub<N>(this->array(),p*N); }

    inline void resize(const size_type &d) { this->array().resize(d*N); }
};

template<class T,int N,int Copy=0>
struct tinyBlockVector
{
  typedef T array_type;
  typedef tiny_block_vector_generator<array_type,N,Copy> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
template<int N,class G> inline typename tinyBlockVector<      Vector<G>,N>::self block(      Vector<G> &X) { return typename tinyBlockVector<      Vector<G>,N>::self(X); }
//{unsecret}
template<int N,class G> inline typename tinyBlockVector<const Vector<G>,N>::self block(const Vector<G> &X) { return typename tinyBlockVector<const Vector<G>,N>::self(X); }


template<class T,int N,int Copy=0>
class tiny_block_dense_vector_generator : public tiny_block_vector_generator<T,N,Copy>
{
  public:
    typedef tiny_block_dense_vector_generator self;
    typedef tiny_block_vector_generator<T,N,Copy> base;
    ARRAY_BASE_TYPES
          
  public:
    inline tiny_block_dense_vector_generator(array_type &x) : base(x) {}
    inline tiny_block_dense_vector_generator(const self &r) : base(r) {}
         
    template<class T2> tiny_block_dense_vector_generator(const Vector<tiny_block_dense_vector_generator<T2,N> > &x) : base((base &)(x.generator())) {}
};

#ifdef NDEBUG
template<      class V        ,int N,int C> struct tinyBlockVector<      Vector<data_vector_generator          <  V   > >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<      Vector<data_vector_generator          <  V   > >,N,C> >::self self; };
template<      class V        ,int N,int C> struct tinyBlockVector<const Vector<data_vector_generator          <  V   > >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<const Vector<data_vector_generator          <  V   > >,N,C> >::self self; };
template<      class V,class A,int N,int C> struct tinyBlockVector<      Vector<dense_vector_generator         <  V,A > >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<      Vector<dense_vector_generator         <  V,A > >,N,C> >::self self; };
template<      class V,class A,int N,int C> struct tinyBlockVector<const Vector<dense_vector_generator         <  V,A > >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<const Vector<dense_vector_generator         <  V,A > >,N,C> >::self self; };
template<int M,class V        ,int N,int C> struct tinyBlockVector<      Vector<tiny_vector_generator          <M,V   > >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<      Vector< tiny_vector_generator         <M,V   > >,N,C> >::self self; };
template<int M,class V        ,int N,int C> struct tinyBlockVector<const Vector<tiny_vector_generator          <M,V   > >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<const Vector< tiny_vector_generator         <M,V   > >,N,C> >::self self; };
template<      class T,int C1 ,int N,int C> struct tinyBlockVector<      Vector<sub_dense_vector_generator     <  T,C1> >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<      Vector<sub_dense_vector_generator     <  T,C1> >,N,C> >::self self; };
template<      class T,int C1 ,int N,int C> struct tinyBlockVector<const Vector<sub_dense_vector_generator     <  T,C1> >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<const Vector<sub_dense_vector_generator     <  T,C1> >,N,C> >::self self; };
template<int M,class T,int C1 ,int N,int C> struct tinyBlockVector<      Vector<tiny_sub_dense_vector_generator<M,T,C1> >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<      Vector<tiny_sub_dense_vector_generator<M,T,C1> >,N,C> >::self self; };
template<int M,class T,int C1 ,int N,int C> struct tinyBlockVector<const Vector<tiny_sub_dense_vector_generator<M,T,C1> >,N,C> { typedef typename GeneratorArray<tiny_block_dense_vector_generator<const Vector<tiny_sub_dense_vector_generator<M,T,C1> >,N,C> >::self self; };
#endif




template<int N,class T,template<int,class> class Cast>
class simd_block_reference : public array_generator<T,1>
{
  public:
    typedef simd_block_reference self;
    typedef array_generator<T,1> base;
    ARRAY_BASE_TYPES

    typedef typename Cast<N,value_type>::self cast_type;

  public:       
    inline simd_block_reference(array_type &x) : base(x) {}
    
    using base::array;
    inline int_type size() const { return N; }
    inline operator cast_type() const { return cast_type(array()); }
    template<class G> inline self &operator=(const Vector<G> &Y) { array()=Y; return *this; }
};

template<int N,class T,template<int,class> class Cast,class A> inline simd_block_reference<N,T,Cast> &operator+=(simd_block_reference<N,T,Cast> &x, const A &y) { x.array()+=y; return x; }
template<int N,class T,template<int,class> class Cast,class A> inline simd_block_reference<N,T,Cast> &operator-=(simd_block_reference<N,T,Cast> &x, const A &y) { x.array()-=y; return x; }
template<int N,class T,template<int,class> class Cast,class A> inline simd_block_reference<N,T,Cast> &operator*=(simd_block_reference<N,T,Cast> &x, const A &y) { x.array()*=y; return x; }
template<int N,class T,template<int,class> class Cast,class A> inline simd_block_reference<N,T,Cast> &operator/=(simd_block_reference<N,T,Cast> &x, const A &y) { x.array()/=y; return x; }

template<int N,class T,template<int,class> class Cast>
class simd_block_dense_reference : public array_traits<T>
{
  public:
    typedef simd_block_dense_reference self;
    typedef array_traits<T> base;
    ARRAY_BASE_TYPES

    typedef typename Cast<N,value_type>::self cast_type;

  public:
    pointer data;
    
  public:       
    inline simd_block_dense_reference(array_type &x) : data(&x.front()) {}
    inline simd_block_dense_reference(const pointer &p) : data(p) {}
    
    inline int_type size() const { return N; }
    inline operator cast_type() const { cast_type r; load(r,data); return r; }
    template<class G> inline self &operator=(const Vector<G> &Y) { store(data,Y); return *this; }
};

template<int N,class T,template<int,class> class Cast>
class simd_unaligned_block_dense_reference : public array_traits<T>
{
  public:
    typedef simd_unaligned_block_dense_reference self;
    typedef array_traits<T> base;
    ARRAY_BASE_TYPES

    typedef typename Cast<N,value_type>::self cast_type;
    
  public:
    pointer data;

  public:       
    inline simd_unaligned_block_dense_reference(array_type &x) : data(&x.front()) {}
    inline simd_unaligned_block_dense_reference(const pointer &p) : data(p) {}
    
    inline int_type size() const { return N; }
    inline operator cast_type() const { cast_type r; loadu(r,data); return r; }
    template<class G> inline self &operator=(const Vector<G> &Y) { storeu(data,Y); return *this; }
};


template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef=simd_block_reference,int C=0>
struct simd_block_vector_generator : public array_reference_generator<T,typename Cast<N,typename T::const_value_type>::self,TempRef<N,typename TinySubVector<N,T>::self,Cast>,const typename Cast<N,typename T::const_value_type>::self,C>
{
  public:
    typedef simd_block_vector_generator self;
    typedef TempRef<N,typename TinySubVector<N,T>::self,Cast> template_reference;
    typedef array_reference_generator<T,typename Cast<N,typename T::const_value_type>::self,template_reference,const typename Cast<N,typename T::const_value_type>::self,C> base;
    ARRAY_BASE_TYPES

  public:
    inline simd_block_vector_generator(array_type &x) : base(x) {}
    inline simd_block_vector_generator(const self &r) : base(r) {}
    
    using base::array;
    
    inline size_type stride() const { return N; }    

    inline index_type lower_bound() const { return array().lower_bound()/N; }
    inline void set_lower_bound(const index_type &i) const { assert(i==lower_bound()); }

    inline size_type size() const { return array().size()/N; }
    inline void resize(const size_type &s) { assert(s==size()); }

    inline reference       operator[](const index_type &i)       { typename TinySubVector<N,      T>::self r=sub<N>(array(),i*N); return reference      (r); }
    inline const_reference operator[](const index_type &i) const { typename TinySubVector<N,const T>::self r=sub<N>(array(),i*N); return const_reference(r); }
};

//template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef=simd_block_dense_reference,int C=0>
template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef,int C>
class simd_block_dense_vector_generator : public simd_block_vector_generator<T,N,Cast,TempRef,C>
{
  public:
    typedef simd_block_dense_vector_generator self;
    typedef simd_block_vector_generator<T,N,Cast,TempRef,C> base;
    ARRAY_BASE_TYPES
    
    typedef typename base::template_reference template_reference;

  protected:
    typename array_traits<reference>::pointer data;

  public:
    inline simd_block_dense_vector_generator(array_type &x) : base(x), data(&x[0]) {}
    inline simd_block_dense_vector_generator(const self &r) : base(r), data(r.data) {}

    inline reference       operator[](const index_type &i)       { return template_reference((typename template_reference::pointer)(data+N*i)); }
    inline const_reference operator[](const index_type &i) const { return template_reference((typename template_reference::pointer)(data+N*i)); }
};


template<class T,int C=0> struct MatrixRow;

template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef=simd_block_reference,int C=0>
struct simd_block_matrix_generator : public array_reference_generator<T,typename Cast<N,typename T::const_value_type>::self,TempRef<N,typename TinySubVector<N,typename MatrixRow<T>::self,1>::self,Cast>,const typename Cast<N,typename T::const_value_type>::self,C>
{
  public:
    typedef simd_block_matrix_generator self;
    typedef TempRef<N,typename TinySubVector<N,typename MatrixRow<T>::self,1>::self,Cast> template_reference;    
    typedef array_reference_generator<T,typename Cast<N,typename T::const_value_type>::self,template_reference,const typename Cast<N,typename T::const_value_type>::self,C> base;
    ARRAY_BASE_TYPES

  public:
    inline simd_block_matrix_generator(array_type &x) : base(x) {}
    inline simd_block_matrix_generator(const self &r) : base(r) {}
    
    using base::array;
    
    inline int_type stride() const { return N; }    

    inline index_type lower_bound() const { return index_type(array().row_lower_bound(),array().col_lower_bound()/N); }
    inline void set_lower_bound(const index_type &i) const { assert(i==lower_bound()); }

    inline size_type size() const { return size_type(array().nrows(),array().ncols()/N); }
    inline void resize(const size_type &s) { assert(s==size()); }

    inline reference       operator[](const index_type &p)       { typename TinySubVector<N,      typename MatrixRow<      T>::self,1>::self r=sub<N>(row(array(),p.i),p.j*N); return reference      (r); }
    inline const_reference operator[](const index_type &p) const { typename TinySubVector<N,const typename MatrixRow<const T>::self,1>::self r=sub<N>(row(array(),p.i),p.j*N); return const_reference(r); }
};

template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef=simd_block_dense_reference,int C=0>
struct simd_block_dense_matrix_generator : public simd_block_matrix_generator<T,N,Cast,TempRef,C>
{
  public:
    typedef simd_block_dense_matrix_generator self;
    typedef simd_block_matrix_generator<T,N,Cast,TempRef,C> base;
    ARRAY_BASE_TYPES
  
    typedef typename base::template_reference template_reference;
    
  public:
    typename array_traits<reference>::pointer data;

  public:
    inline simd_block_dense_matrix_generator(array_type &x) : base(x), data(&x(0,0)) {}
    inline simd_block_dense_matrix_generator(const self &r) : base(r), data(r.data) {}
    
    using base::array;
        
    inline reference       operator[](const index_type &i)       { return template_reference((typename template_reference::pointer)(data+i.i*array().ncols()+N*i.j)); }
    inline const_reference operator[](const index_type &i) const { return template_reference((typename template_reference::pointer)(data+i.i*array().ncols()+N*i.j)); }
    
    int_type stride() const { return array().ncols(); }
};

template<class T,int N,template<int,class> class Cast,template<int,class,template<int,class> class> class TempRef=simd_block_dense_reference,int C=0>
struct simd_sub_dense_matrix_generator : public simd_block_matrix_generator<T,N,Cast,TempRef,C>
{
  public:
    typedef simd_sub_dense_matrix_generator self;
    typedef simd_block_matrix_generator<T,N,Cast,TempRef,C> base;
    ARRAY_BASE_TYPES
  
    typedef typename base::template_reference template_reference;
    
  public:
    typename array_traits<reference>::pointer data;

  public:
    inline simd_sub_dense_matrix_generator(array_type &x) : base(x), data(&x(0,0)) {}
    inline simd_sub_dense_matrix_generator(const self &r) : base(r), data(r.data) {}
    
    using base::array;
        
    inline reference       operator[](const index_type &i)       { return template_reference((typename template_reference::pointer)(data+i.i*array().generator().row_stride()+N*i.j)); }
    inline const_reference operator[](const index_type &i) const { return template_reference((typename template_reference::pointer)(data+i.i*array().generator().row_stride()+N*i.j)); }
    
    inline int_type row_stride() const { return array().generator().row_stride(); }
};

template<class V> struct SimdLength { enum { RET=0}; };
template<class V> struct SimdLength<const V> { enum { RET=SimdLength<V>::RET}; };
#if defined(MMX) && !defined(SSE2)
template<> struct SimdLength<unsigned char   > { enum { RET=8  }; };
template<> struct SimdLength<         char   > { enum { RET=8  }; };
template<> struct SimdLength<unsigned short  > { enum { RET=4  }; };
template<> struct SimdLength<         short  > { enum { RET=4  }; };
template<> struct SimdLength<         int    > { enum { RET=2  }; };
#endif
#if defined(SSE)
template<> struct SimdLength<        float   > { enum { RET=4  }; };
template<> struct SimdLength<complex<float>  > { enum { RET=2  }; };
#endif
#if defined(SSE2)
template<> struct SimdLength<unsigned char   > { enum { RET=16 }; };
template<> struct SimdLength<         char   > { enum { RET=16 }; };
template<> struct SimdLength<unsigned short  > { enum { RET=8  }; };
template<> struct SimdLength<         short  > { enum { RET=8  }; };
template<> struct SimdLength<         int    > { enum { RET=4  }; };
template<> struct SimdLength<        double  > { enum { RET=2  }; };
template<> struct SimdLength<complex<double> > { enum { RET=1  }; };
#endif

template<class V> struct SimdBlockLength
{
  enum { RET=(SimdLength<V>::RET==0)?1:SimdLength<V>::RET };
};

template<class T,int L,template<int,class> class Cast,int C=0> struct AuxSimdBlock { };
template<class G,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<G>,L,Cast,C> { typedef simd_block_vector_generator<      Vector<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class G,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<G>,L,Cast,C> { typedef simd_block_vector_generator<const Vector<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class G,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Matrix<G>,L,Cast,C> { typedef simd_block_matrix_generator<      Matrix<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class G,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Matrix<G>,L,Cast,C> { typedef simd_block_matrix_generator<const Matrix<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<            class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<data_vector_generator          <V  >      >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<data_vector_generator          <V  >    >::value_type>::self>::self self; };
template<            class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<data_vector_generator          <V  >      >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<data_vector_generator          <V  >    >::value_type>::self>::self self; };
template<class A    ,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<dense_vector_generator         <V,A>      >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<dense_vector_generator         <V,A>    >::value_type>::self>::self self; };
template<class A    ,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<dense_vector_generator         <V,A>      >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<dense_vector_generator         <V,A>    >::value_type>::self>::self self; };
template<      class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<sub_dense_vector_generator     <  T,C1>   >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<sub_dense_vector_generator     <  T,C1> >::value_type>::self>::self self; };
template<      class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<sub_dense_vector_generator     <  T,C1>   >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<sub_dense_vector_generator     <  T,C1> >::value_type>::self>::self self; };
template<      class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<row_dense_matrix_generator     <  T,C1>   >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<row_dense_matrix_generator     <  T,C1> >::value_type>::self>::self self; };
template<      class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<row_dense_matrix_generator     <  T,C1>   >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<row_dense_matrix_generator     <  T,C1> >::value_type>::self>::self self; };
template<int N      ,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<tiny_vector_generator          <N,  V>    >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<tiny_vector_generator          <N,  V>  >::value_type>::self>::self self; };
template<int N      ,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<tiny_vector_generator          <N,  V>    >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<tiny_vector_generator          <N,  V>  >::value_type>::self>::self self; };
template<int N,class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<tiny_sub_dense_vector_generator<N,T,C1>   >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<tiny_sub_dense_vector_generator<N,T,C1> >::value_type>::self>::self self; };
template<int N,class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<tiny_sub_dense_vector_generator<N,T,C1>   >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<tiny_sub_dense_vector_generator<N,T,C1> >::value_type>::self>::self self; };
template<int N,class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Vector<tiny_row_dense_matrix_generator<N,T,C1>   >,L,Cast,C> { typedef       typename DataVector<typename Cast<L,typename Vector<tiny_row_dense_matrix_generator<N,T,C1> >::value_type>::self>::self self; };
template<int N,class T,int C1      ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Vector<tiny_row_dense_matrix_generator<N,T,C1>   >,L,Cast,C> { typedef const typename DataVector<typename Cast<L,typename Vector<tiny_row_dense_matrix_generator<N,T,C1> >::value_type>::self>::self self; };
template<            class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Matrix<data_matrix_generator          <V  >      >,L,Cast,C> { typedef       typename DataMatrix<typename Cast<L,typename Matrix<data_matrix_generator          <V  >    >::value_type>::self>::self self; };
template<            class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Matrix<data_matrix_generator          <V  >      >,L,Cast,C> { typedef const typename DataMatrix<typename Cast<L,typename Matrix<data_matrix_generator          <V  >    >::value_type>::self>::self self; };
template<int M,int N,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Matrix<tiny_matrix_generator          <M,N,V>    >,L,Cast,C> { typedef       typename DataMatrix<typename Cast<L,typename Matrix<tiny_matrix_generator          <M,N,V>  >::value_type>::self>::self self; };
template<int M,int N,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Matrix<tiny_matrix_generator          <M,N,V>    >,L,Cast,C> { typedef const typename DataMatrix<typename Cast<L,typename Matrix<tiny_matrix_generator          <M,N,V>  >::value_type>::self>::self self; };
template<class A    ,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Matrix<dense_matrix_generator         <V,A>      >,L,Cast,C> { typedef       typename DataMatrix<typename Cast<L,typename Matrix<dense_matrix_generator         <V,A>    >::value_type>::self>::self self; };
template<class A    ,class V       ,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Matrix<dense_matrix_generator         <V,A>      >,L,Cast,C> { typedef const typename DataMatrix<typename Cast<L,typename Matrix<dense_matrix_generator         <V,A>    >::value_type>::self>::self self; };
template<            class T,int C2,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Matrix<     sub_dense_matrix_generator<    T,C2> >,L,Cast,C> { typedef simd_sub_dense_matrix_generator<      Matrix<     sub_dense_matrix_generator<    T,C2> >,L,Cast,simd_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<            class T,int C2,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Matrix<     sub_dense_matrix_generator<    T,C2> >,L,Cast,C> { typedef simd_sub_dense_matrix_generator<const Matrix<     sub_dense_matrix_generator<    T,C2> >,L,Cast,simd_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T,int C2,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,L,Cast,C> { typedef simd_sub_dense_matrix_generator<      Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,L,Cast,simd_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T,int C2,int L,template<int,class> class Cast,int C> struct AuxSimdBlock<const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,L,Cast,C> { typedef simd_sub_dense_matrix_generator<const Matrix<tiny_sub_dense_matrix_generator<M,N,T,C2> >,L,Cast,simd_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };




template<class T             ,int L,template<int,class> class Cast,int C=0> struct AuxSimdUnalignedBlock { };
template<class G             ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<G>,L,Cast,C> { typedef simd_block_vector_generator<      Vector<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class G             ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<G>,L,Cast,C> { typedef simd_block_vector_generator<const Vector<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class G             ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Matrix<G>,L,Cast,C> { typedef simd_block_matrix_generator<      Matrix<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class G             ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Matrix<G>,L,Cast,C> { typedef simd_block_matrix_generator<const Matrix<G>,L,Cast,simd_block_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<            class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<data_vector_generator          <V  >    >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<data_vector_generator          <V  >    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<            class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<data_vector_generator          <V  >    >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<data_vector_generator          <V  >    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class A    ,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<dense_vector_generator         <V,A>    >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<dense_vector_generator         <V,A>    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class A    ,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<dense_vector_generator         <V,A>    >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<dense_vector_generator         <V,A>    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N      ,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<tiny_vector_generator          <N,  V>  >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<tiny_vector_generator          <N,  V>  >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N      ,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<tiny_vector_generator          <N,  V>  >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<tiny_vector_generator          <N,  V>  >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<sub_dense_vector_generator     <  T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<sub_dense_vector_generator     <  T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<sub_dense_vector_generator     <  T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<sub_dense_vector_generator     <  T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<row_dense_matrix_generator     <  T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<row_dense_matrix_generator     <  T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<row_dense_matrix_generator     <  T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<row_dense_matrix_generator     <  T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<tiny_sub_dense_vector_generator<N,T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<tiny_sub_dense_vector_generator<N,T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<tiny_sub_dense_vector_generator<N,T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<tiny_sub_dense_vector_generator<N,T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Vector<tiny_row_dense_matrix_generator<N,T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<      Vector<tiny_row_dense_matrix_generator<N,T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,int C1,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Vector<tiny_row_dense_matrix_generator<N,T,C1> >,L,Cast,C> { typedef simd_block_dense_vector_generator<const Vector<tiny_row_dense_matrix_generator<N,T,C1> >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<            class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Matrix<data_matrix_generator          <V  >    >,L,Cast,C> { typedef simd_block_dense_matrix_generator<      Matrix<data_matrix_generator          <V  >    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<            class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Matrix<data_matrix_generator          <V  >    >,L,Cast,C> { typedef simd_block_dense_matrix_generator<const Matrix<data_matrix_generator          <V  >    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class A    ,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Matrix<dense_matrix_generator         <V,A>    >,L,Cast,C> { typedef simd_block_dense_matrix_generator<      Matrix<dense_matrix_generator         <V,A>    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class A    ,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Matrix<dense_matrix_generator         <V,A>    >,L,Cast,C> { typedef simd_block_dense_matrix_generator<const Matrix<dense_matrix_generator         <V,A>    >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<      Matrix<tiny_matrix_generator          <M,N,V>  >,L,Cast,C> { typedef simd_block_dense_matrix_generator<      Matrix<tiny_matrix_generator          <M,N,V>  >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V ,int L,template<int,class> class Cast,int C> struct AuxSimdUnalignedBlock<const Matrix<tiny_matrix_generator          <M,N,V>  >,L,Cast,C> { typedef simd_block_dense_matrix_generator<const Matrix<tiny_matrix_generator          <M,N,V>  >,L,Cast,simd_unaligned_block_dense_reference,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<class T,int L,template<int,class> class Cast,int C> struct Aux2SimdBlock                      { typedef typename AuxSimdBlock         <T,L,Cast,C>::self self; typedef typename self::value_type value_type; };
template<class T,int L,template<int,class> class Cast,int C> struct Aux2SimdUnalignedBlock             { typedef typename AuxSimdUnalignedBlock<T,L,Cast,C>::self self; typedef typename self::value_type value_type; };
template<class T,      template<int,class> class Cast,int C> struct Aux2SimdBlock         <T,0,Cast,C> { typedef T &self; typedef typename T::value_type value_type; };
template<class T      ,template<int,class> class Cast,int C> struct Aux2SimdUnalignedBlock<T,0,Cast,C> { typedef T &self; typedef typename T::value_type value_type; };

template<class T,int L=SimdLength<typename T::value_type>::RET,int C=0> struct SimdBlock                   { typedef typename Aux2SimdBlock         <T,L,SimdVector             ,C>::self self; typedef typename Aux2SimdBlock         <T,L,SimdVector             ,C>::value_type value_type; };
template<class T,int L=SimdLength<typename T::value_type>::RET,int C=0> struct SaturatedSimdBlock          { typedef typename Aux2SimdBlock         <T,L,SaturatedSimdTinyVector,C>::self self; typedef typename Aux2SimdBlock         <T,L,SaturatedSimdTinyVector,C>::value_type value_type; };
template<class T,int L=SimdLength<typename T::value_type>::RET,int C=0> struct SimdUnalignedBlock          { typedef typename Aux2SimdUnalignedBlock<T,L,SimdVector             ,C>::self self; typedef typename Aux2SimdUnalignedBlock<T,L,SimdVector             ,C>::value_type value_type; };
template<class T,int L=SimdLength<typename T::value_type>::RET,int C=0> struct SaturatedSimdUnalignedBlock { typedef typename Aux2SimdUnalignedBlock<T,L,SaturatedSimdTinyVector,C>::self self; typedef typename Aux2SimdUnalignedBlock<T,L,SaturatedSimdTinyVector,C>::value_type value_type; };

template<int L,class G> inline typename SimdBlock<      Vector<G>,L>::self simd_block(      Vector<G> &X) { return typename SimdBlock<      Vector<G>,L>::self(X); }
template<int L,class G> inline typename SimdBlock<const Vector<G>,L>::self simd_block(const Vector<G> &X) { return typename SimdBlock<const Vector<G>,L>::self(X); }
template<int L,class G> inline typename SimdBlock<      Matrix<G>,L>::self simd_block(      Matrix<G> &X) { return typename SimdBlock<      Matrix<G>,L>::self(X); }
template<int L,class G> inline typename SimdBlock<const Matrix<G>,L>::self simd_block(const Matrix<G> &X) { return typename SimdBlock<const Matrix<G>,L>::self(X); }


//{unsecret}
//Summary: Decomposition in blocks for parallel calculation with SIMD instructions
//Arguments:
//  X - The array to decompose
//Remarks:
//  Does not decompose the array if no SIMD header file was previously included.
//  The array is assumed to be aligned in memory.
//Result: An array representing every blocks that fit in a SIMD register
//Example:
//  DenseVector<short int>::self X(16, "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15");
//  cout << simd_block(X)[0] << endl; [0 1 2 3 4 5 6 7] with SSE2
//  simd_block(Y)=simd_block(X) + simd_block(Y); // quicker than Y=X+Y;
//See: ^saturated_simd_block^, ^simd_unaligned_block^, ^saturated_simd_unaligned_block^
template<class G> inline typename SimdBlock<      Vector<G> >::self simd_block(      Vector<G> &X) { return typename SimdBlock<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SimdBlock<const Vector<G> >::self simd_block(const Vector<G> &X) { return typename SimdBlock<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SimdBlock<      Matrix<G> >::self simd_block(      Matrix<G> &X) { return typename SimdBlock<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename SimdBlock<const Matrix<G> >::self simd_block(const Matrix<G> &X) { return typename SimdBlock<const Matrix<G> >::self(X); }


//{unsecret}
//Summary: Decomposition in blocks for parallel calculation with saturated SIMD instructions
//Arguments:
//  X - The array to decompose
//Remarks:
//  Only a few operations (addition, substraction) differ from unsatarated operations,
//  and only for a few data types (char, unsigned char, short, unsigned short).
//Result: An array representing every blocks that fit in a SIMD register
//Example:
//  DenseVector<unsigned char>::self X(16, 250);
//  DenseVector<unsigned char>::self Y(16, 10);
//  cout << simd_block(X)+10 << endl; // 4 4 4...
//  cout << saturated_simd_block(X)+10 << endl; // 255 255 255...
//See: ^simd_block^, ^simd_unaligned_block^, ^saturated_simd_unaligned_block^
template<class G> inline typename SaturatedSimdBlock<      Vector<G> >::self saturated_simd_block(      Vector<G> &X) { return typename SaturatedSimdBlock<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SaturatedSimdBlock<const Vector<G> >::self saturated_simd_block(const Vector<G> &X) { return typename SaturatedSimdBlock<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SaturatedSimdBlock<      Matrix<G> >::self saturated_simd_block(      Matrix<G> &X) { return typename SaturatedSimdBlock<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename SaturatedSimdBlock<const Matrix<G> >::self saturated_simd_block(const Matrix<G> &X) { return typename SaturatedSimdBlock<const Matrix<G> >::self(X); }


//{unsecret}
//Summary: Decomposition in blocks for parallel calculation with SIMD instructions
//Arguments:
//  X - The array to decompose
//Remarks:
//  No assumption is made about the alignment in memory, resulting in slower calculations than ^simd_block^.
//Result: An array representing every blocks that fit in a SIMD register
//See: ^simd_block^, ^saturated_simd_block^, ^saturated_simd_unaligned_block^
template<class G> inline typename SimdUnalignedBlock<      Vector<G> >::self simd_unaligned_block(      Vector<G> &X) { return typename SimdUnalignedBlock<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SimdUnalignedBlock<const Vector<G> >::self simd_unaligned_block(const Vector<G> &X) { return typename SimdUnalignedBlock<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SimdUnalignedBlock<      Matrix<G> >::self simd_unaligned_block(      Matrix<G> &X) { return typename SimdUnalignedBlock<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename SimdUnalignedBlock<const Matrix<G> >::self simd_unaligned_block(const Matrix<G> &X) { return typename SimdUnalignedBlock<const Matrix<G> >::self(X); }


//{unsecret}
//Summary: Decomposition in blocks for parallel calculation with saturated SIMD instructions
//Arguments:
//  X - The array to decompose
//Result: An array representing every blocks that fit in a SIMD register
//See: ^simd_block^, ^saturated_simd_block^, ^simd_unaligned_block^
template<class G> inline typename SaturatedSimdUnalignedBlock<      Vector<G> >::self saturated_simd_unaligned_block(      Vector<G> &X) { return typename SaturatedSimdUnalignedBlock<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SaturatedSimdUnalignedBlock<const Vector<G> >::self saturated_simd_unaligned_block(const Vector<G> &X) { return typename SaturatedSimdUnalignedBlock<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename SaturatedSimdUnalignedBlock<      Matrix<G> >::self saturated_simd_unaligned_block(      Matrix<G> &X) { return typename SaturatedSimdUnalignedBlock<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename SaturatedSimdUnalignedBlock<const Matrix<G> >::self saturated_simd_unaligned_block(const Matrix<G> &X) { return typename SaturatedSimdUnalignedBlock<const Matrix<G> >::self(X); }


////{unsecret}
//template<class G> inline typename SubArray<      Vector<G> >::self simd_block_rest(      Vector<G> &X) { int n=X.size(); const int N=SimdLength<typename Vector<G>::value_type>::RET; int i=N*(n/N); return typename SubArray<      Vector<G> >::self(X,i,X.size()-i); }
////{unsecret}
//template<class G> inline typename SubArray<const Vector<G> >::self simd_block_rest(const Vector<G> &X) { int n=X.size(); const int N=SimdLength<typename Vector<G>::value_type>::RET; int i=N*(n/N); return typename SubArray<const Vector<G> >::self(X,i,X.size()-i); }
////{unsecret}
//template<class G> inline typename SubArray<      Matrix<G> >::self simd_block_rest(      Matrix<G> &X) { int n=X.ncols(); const int N=SimdLength<typename Matrix<G>::value_type>::RET; int i=N*(n/N); return sub(X,0,i,X.nrows(),X.ncols()-i); }
////{unsecret}
//template<class G> inline typename SubArray<const Matrix<G> >::self simd_block_rest(const Matrix<G> &X) { int n=X.ncols(); const int N=SimdLength<typename Matrix<G>::value_type>::RET; int i=N*(n/N); return sub(X,0,i,X.nrows(),X.ncols()-i); }



template<class T>
class unblock_generator : public array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference>
{
  public:
    typedef unblock_generator self;
    typedef array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference> base;
    ARRAY_BASE_TYPES

  private:
    size_type dim;
    
  public:
    inline unblock_generator(array_type &x) : base(x), dim(x[0].size()) {}

    inline reference       operator[](const index_type &n)       { return (this->array())[n/dim][mod(n,dim)]; }
    inline const_reference operator[](const index_type &n) const { return (this->array())[n/dim][mod(n,dim)]; }

    inline size_type size () const { return dim*(this->array()).size(); }

    inline void resize(size_type n) { assert(size()==n); }
};

template<class A>
struct UnblockArray
{
  typedef A array_type;
  typedef unblock_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> typename UnblockArray<      Vector<G> >::self unblock(      Vector<G> &X) { return typename UnblockArray<      Vector<G> >::self(X); }
template<class G> typename UnblockArray<const Vector<G> >::self unblock(const Vector<G> &X) { return typename UnblockArray<const Vector<G> >::self(X); }
template<class G> typename UnblockArray<      Matrix<G> >::self unblock(      Matrix<G> &X) { return typename UnblockArray<      Matrix<G> >::self(X); }
template<class G> typename UnblockArray<const Matrix<G> >::self unblock(const Matrix<G> &X) { return typename UnblockArray<const Matrix<G> >::self(X); }



template<class T>
class block_function : public unary_reference_function<T,typename blockArray<T,1>::self,typename blockArray<T,1>::self,typename blockArray<const T,1>::self>
{
  public:
    typedef block_function self;
    typedef unary_reference_function<T,typename blockArray<T,1>::self,typename blockArray<T,1>::self,typename blockArray<const T,1>::self> base;
    UNARY_FUNCTION_BASE_TYPES
    
    typedef T array_type;
    typedef typename array_type::index_type index_type;
    typedef typename array_type::size_type  size_type;
        
  private:
    size_type dim;
    
  public:
    inline block_function(const size_type &s) : dim(s) {}
    
    inline reference       operator()(      argument_type &x)       { return reference      (x,dim); }
    inline const_reference operator()(const argument_type &x) const { return const_reference(x,dim); }
};

template<class A>
struct eleblockArray
{
  typedef A array_type;
  typedef typename array_traits<array_type>::value_type value_type;
  typedef block_function<value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<class A>
struct eleblockArray<const A>
{
  typedef const A array_type;
  typedef const typename array_type::value_type value_type;
  typedef const block_function<value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<class G> inline typename eleblockArray<      Vector<G> >::self ele_block(      Vector<G> &X, const typename Vector<G>::value_type::size_type &s) { return apply(X,typename eleblockArray<      Vector<G> >::function_type(s)); }
template<class G> inline typename eleblockArray<const Vector<G> >::self ele_block(const Vector<G> &X, const typename Vector<G>::value_type::size_type &s) { return apply(X,typename eleblockArray<const Vector<G> >::function_type(s)); }
template<class G> inline typename eleblockArray<      Matrix<G> >::self ele_block(      Matrix<G> &X, const typename Matrix<G>::value_type::size_type &s) { return apply(X,typename eleblockArray<      Matrix<G> >::function_type(s)); }
template<class G> inline typename eleblockArray<const Matrix<G> >::self ele_block(const Matrix<G> &X, const typename Matrix<G>::value_type::size_type &s) { return apply(X,typename eleblockArray<const Matrix<G> >::function_type(s)); }
template<class G> inline typename eleblockArray<      Matrix<G> >::self ele_block(      Matrix<G> &X, const typename Matrix<G>::value_type::size_type::value_type &m, const typename Matrix<G>::value_type::size_type::value_type &n) { return ele_block(X,typename Matrix<G>::value_type::size_type(m,n)); }
template<class G> inline typename eleblockArray<const Matrix<G> >::self ele_block(const Matrix<G> &X, const typename Matrix<G>::value_type::size_type::value_type &m, const typename Matrix<G>::value_type::size_type::value_type &n) { return ele_block(X,typename Matrix<G>::value_type::size_type(m,n)); }


template<int N,class T>
class tiny_block_vector_function : public unary_reference_function<T,typename tinyBlockVector<T,N,1>::self,typename tinyBlockVector<T,N,1>::self,typename tinyBlockVector<const T,N,1>::self>
{
  public:
    typedef tiny_block_vector_function self;
    typedef unary_reference_function<T,typename tinyBlockVector<T,N,1>::self,typename tinyBlockVector<T,N,1>::self,typename tinyBlockVector<const T,N,1>::self> base;
    UNARY_FUNCTION_BASE_TYPES
        
    typedef T array_type;
    typedef typename array_type::index_type index_type;
    typedef typename array_type::size_type  size_type;
            
  public:
    inline tiny_block_vector_function() {}
    
    inline reference       operator()(argument_type       &x)       { return reference      (x); }
    inline const_reference operator()(const argument_type &x) const { return const_reference(x); }
};

template<int N,class T>
struct eleTinyBlockVector
{
  typedef T array_type;
  typedef typename array_traits<array_type>::value_type value_type;
  typedef tiny_block_vector_function<N,value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<int N,class T>
struct eleTinyBlockVector<N,const T>
{
  typedef const T array_type;
  typedef const typename array_type::value_type value_type;
  typedef const tiny_block_vector_function<N,value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<int N,class G> inline typename eleTinyBlockVector<N,      Vector<G> >::self ele_block(      Vector<G> &X) { return apply(X,typename eleTinyBlockVector<N,      Vector<G> >::function_type()); }
template<int N,class G> inline typename eleTinyBlockVector<N,const Vector<G> >::self ele_block(const Vector<G> &X) { return apply(X,typename eleTinyBlockVector<N,const Vector<G> >::function_type()); }



template<class I>
inline I sample_lower_bound(const I &p, const I &d)
{
  return (p>=0) ? p/d : (p-1)/d+1; // pas sûr du tout
}

template<class I,class S>
inline MatrixIndex<I> sample_lower_bound(const MatrixIndex<I> &p, const S &d)
{
  return MatrixIndex<I>(sample_lower_bound(p.i,d),sample_lower_bound(p.j,d));
}

template<class T>
class sample_generator : public array_generator<T>
{
  public:
    typedef sample_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES 

  private:
    difference_type dp;

  public:
    inline sample_generator(array_type &x, const difference_type &p) : base(x), dp(p) {}
    inline sample_generator(const self &r) : base(r), dp(r.dp) {}
    
    inline size_type size() const { return (this->array()).size()/dp; } // pas bon pour n'importe quel lower_bound
    inline index_type lower_bound() const { return sample_lower_bound((this->array()).lower_bound(),dp); }
     
    inline reference       operator[](const index_type &p)       { return (this->array())[dp*p]; }
    inline const_reference operator[](const index_type &p) const { return (this->array())[dp*p]; }
    
    inline void inv(const value_type &y, const index_type &p) { (*this)[p]=y; }

    inline void resize(const size_type &d) { assert(d==size()); }     
};

template<class T> 
struct sampleArray : public array_traits<T>
{
  typedef array_traits<T> base;
  ARRAY_BASE_TYPES
  typedef sample_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
};

template<class G> inline typename sampleArray<Vector<G> >::self sample(Vector<G> &X, const typename G::int_type &n) { return typename sampleArray<Vector<G> >::self(X, n); }
template<class G> inline typename sampleArray<Matrix<G> >::self sample(Matrix<G> &X, const typename G::int_type &n) { return typename sampleArray<Matrix<G> >::self(X, n); }
template<class G> inline typename sampleArray<const Vector<G> >::self sample(const Vector<G> &X, const typename G::int_type &n) { return typename sampleArray<const Vector<G> >::self(X, n); }
template<class G> inline typename sampleArray<const Matrix<G> >::self sample(const Matrix<G> &X, const typename G::int_type &n) { return typename sampleArray<const Matrix<G> >::self(X, n); }


template<class T>
class upsample_generator : public array_value_generator<const T>
{
  public:
    typedef upsample_generator self;
    typedef array_value_generator<const T> base;
    ARRAY_BASE_TYPES

  private:
    difference_type dp;

  public:
    inline upsample_generator(const array_type &x, const difference_type &p) : base(x), dp(p) {}
    inline upsample_generator(const self &r) : base(r), dp(r.dp) {}
    
    inline size_type size() const { return dp*((this->array()).size()-1)+1; }
     
    inline const_reference operator[](const index_type &p) const { return (p%dp)?0:(this->array())[p/dp]; }

    inline void resize(const size_type &d) { assert(d==size()); }     
};

template<class T> 
struct upsampleArray : public array_traits<T>
{
  typedef array_traits<T> base;
  ARRAY_BASE_TYPES
  typedef upsample_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
};

template<class G> inline typename upsampleArray<const Vector<G> >::self upsample(const Vector<G> &A, const typename G::int_type &n) { return typename upsampleArray<const Vector<G> >::self(A, n); }
template<class G> inline typename upsampleArray<const Matrix<G> >::self upsample(const Matrix<G> &A, const typename G::int_type &n) { return typename upsampleArray<const Matrix<G> >::self(A, n); }


template<class T>
class flip_array_generator : public array_generator<T>
{
  public:
    typedef flip_array_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

  public:
    inline flip_array_generator(array_type &x) : base(x) {}
     
    inline reference       operator[](const index_type &i)       { return (this->array())[(this->array()).upper_bound()+(this->array()).lower_bound()-i]; }
    inline const_reference operator[](const index_type &i) const { return (this->array())[(this->array()).upper_bound()+(this->array()).lower_bound()-i]; }
};

template<class T>
struct flipArray 
{
  typedef T array_type;
  typedef flip_array_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
}; 

template<class G> inline typename flipArray<const Vector<G> >::self flip(const Vector<G> &X) { return typename flipArray<const Vector<G> >::self(X); }
template<class G> inline typename flipArray<const Matrix<G> >::self flip(const Matrix<G> &X) { return typename flipArray<const Matrix<G> >::self(X); }


template<class T>
class mirror_array_generator : public array_generator<T>
{
  public:
    typedef mirror_array_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

    template<class A2> struct iterator_rebind       { typedef typename shift_iterator_rebind<      A2>::other other; };
    template<class A2> struct const_iterator_rebind { typedef typename shift_iterator_rebind<const A2>::other other; };

  public:
    inline mirror_array_generator(array_type &x) : base(x) {}
     
    inline reference       operator[](const index_type &p)       { return (this->array())[-p]; }
    inline const_reference operator[](const index_type &p) const { return (this->array())[-p]; }

    inline index_type lower_bound() const { return -(this->array()).upper_bound(); }
};

template<class T>
struct mirrorArray 
{
  typedef T array_type;
  typedef mirror_array_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
}; 

template<class G> inline typename mirrorArray<      Vector<G> >::self mirror(      Vector<G> &X) { return typename mirrorArray<      Vector<G> >::self(X); }
template<class G> inline typename mirrorArray<const Vector<G> >::self mirror(const Vector<G> &X) { return typename mirrorArray<const Vector<G> >::self(X); }
template<class G> inline typename mirrorArray<      Matrix<G> >::self mirror(      Matrix<G> &X) { return typename mirrorArray<      Matrix<G> >::self(X); }
template<class G> inline typename mirrorArray<const Matrix<G> >::self mirror(const Matrix<G> &X) { return typename mirrorArray<const Matrix<G> >::self(X); }


template<class A>
class buffer_vector_generator : public array_value_generator<A>
{
  public:
    typedef buffer_vector_generator self;
    typedef array_value_generator<A> base;
    ARRAY_BASE_TYPES

    typedef pair<bool, value_type> buffer_value_type;
    typedef typename DenseVector<buffer_value_type>::self buffer_type;

  private:
    mutable buffer_type buff;

  public:
    inline buffer_vector_generator(array_type &x) : base(x), buff(x.size(), buffer_value_type(false, value_type())) {}
    inline buffer_vector_generator(const self &r) : base(r), buff(r.buff) {}

    inline const_reference operator[](const index_type &i) const { if (!buff[i].first) buff[i]=buffer_value_type(true,(this->array())[i]); return buff[i].second; }
};

template<class A>
class buffer_matrix_generator : public array_value_generator<A>
{
  public:
    typedef buffer_matrix_generator self;
    typedef array_value_generator<A> base;
    ARRAY_BASE_TYPES

    typedef pair<bool, value_type> buffer_value_type;
    typedef typename DenseMatrix<buffer_value_type>::self buffer_type;

  private:
    mutable buffer_type buff;

  public:
    inline buffer_matrix_generator(array_type &x) : base(x), buff(x.size(), buffer_value_type(false, value_type())) {}
    inline buffer_matrix_generator(const self &r) : base(r), buff(r.buff) {}

    inline const_reference operator[](const index_type &i) const { if (!buff[i].first) buff[i]=buffer_value_type(true,(this->array())[i]); return buff[i].second; }
};

template<class T>
struct bufferVector : public array_value_traits<T>
{
  typedef array_value_traits<T> base;
  ARRAY_BASE_TYPES
  typedef buffer_vector_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class T>
struct bufferMatrix : public array_value_traits<T>
{
  typedef array_value_traits<T> base;
  ARRAY_BASE_TYPES
  typedef buffer_matrix_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename bufferVector<const Vector<G> >::self buffer(const Vector<G> &X) { return typename bufferVector<const Vector<G> >::self(X); }
template<class G> inline typename bufferMatrix<const Matrix<G> >::self buffer(const Matrix<G> &X) { return typename bufferMatrix<const Matrix<G> >::self(X); }


template<class T>
class zone_array_generator : public array_value_generator<T, typename SubArray<T>::self>
{
  public:
    typedef zone_array_generator self;
    typedef array_value_generator<T, typename SubArray<T>::self> base;
    ARRAY_BASE_TYPES

  private:
    size_type dim;

  public:
    inline zone_array_generator(array_type &x, const size_type &d) : base(x), dim(d) {}
    inline zone_array_generator(const self &r) : base(r), dim(r.dim) {}

    using base::array;

    inline size_type size() const { return array().size()-dim+1; }

    inline const_reference operator[](const index_type &i) const { return sub(array(),i,dim); }

    inline void resize(const size_type &s) { assert(s==size()); }
};

template<class T>
struct zoneArray
{
  typedef T array_type;
  typedef zone_array_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename zoneArray<Vector<G> >::self zone(Vector<G> &X, const typename Vector<G>::size_type &d) { return zoneArray<Vector<G> >::self(X,d); }
template<class G> inline typename zoneArray<Matrix<G> >::self zone(Matrix<G> &X, const typename Matrix<G>::size_type &d) { return zoneArray<Matrix<G> >::self(X,d); }
template<class G> inline typename zoneArray<const Vector<G> >::self zone(const Vector<G> &X, const typename Vector<G>::size_type &d) { return zoneArray<const Vector<G> >::self(X,d); }
template<class G> inline typename zoneArray<const Matrix<G> >::self zone(const Matrix<G> &X, const typename Matrix<G>::size_type &d) { return zoneArray<const Matrix<G> >::self(X,d); }

template<int M,int N,class T>
class tiny_zone_matrix_generator : public array_value_generator<T, typename TinySubMatrix<M,N,T>::self>
{
  public:
    typedef tiny_zone_matrix_generator self;
    typedef array_value_generator<T, typename TinySubMatrix<M,N,T>::self> base;
    ARRAY_BASE_TYPES

  public:
    using base::array;
  
    inline tiny_zone_matrix_generator(array_type &x) : base(x) {}
    inline tiny_zone_matrix_generator(const self &r) : base(r) {}

    inline size_type size() const { return (array()).size()-size_type(M-1,N-1); }

    inline const_reference operator[](const index_type &i) const { return sub<M,N>(array(),i); }

    inline void resize(const size_type &s) { assert(s==size()); }
};

template<int M,int N,class T>
struct tinyZoneMatrix
{
  typedef T array_type;
  typedef tiny_zone_matrix_generator<M,N,array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int M,int N,class G> inline typename tinyZoneMatrix<M,N,      Matrix<G> >::self zone(      Matrix<G> &X) { return typename tinyZoneMatrix<M,N,      Matrix<G> >::self(X); }
template<int M,int N,class G> inline typename tinyZoneMatrix<M,N,const Matrix<G> >::self zone(const Matrix<G> &X) { return typename tinyZoneMatrix<M,N,const Matrix<G> >::self(X); }



template<class A1,class A2>
class mask_array_generator : public array_value_generator<const A1>
{
  public:
    typedef mask_array_generator self;
    typedef array_value_generator<const A1> base;
    ARRAY_BASE_TYPES

    typedef A2 mask_type;

  private:
    mask_type mask;

  public:
    inline mask_array_generator(array_type &x, mask_type &m) : base(x), mask(m) {}
    inline mask_array_generator(const self &x) : base(x), mask(x.mask) {}

    inline const_reference operator[](const index_type &i) const { if (mask[i]) return (this->array())[i]; else return value_type(); }
};


template<class A1,class A2>
struct maskArray : public array_traits<A1>
{
  typedef array_traits<A1> base;
  ARRAY_BASE_TYPES
  typedef A2 mask_type;
  typedef mask_array_generator<array_type,mask_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G1,class G2> inline typename maskArray<const Vector<G1>,const Vector<G2> >::self mask(const Vector<G1> &X, const Vector<G2> &M) { return maskArray<const Vector<G1>,const Vector<G2> >::self(X,M); }
template<class G1,class G2> inline typename maskArray<const Matrix<G1>,const Matrix<G2> >::self mask(const Matrix<G1> &X, const Matrix<G2> &M) { return maskArray<const Matrix<G1>,const Matrix<G2> >::self(X,M); }


template<class A>
class prolong_array_generator : public array_value_generator<const A>
{
  public:
    typedef prolong_array_generator self; 
    typedef array_value_generator<const A> base;
    ARRAY_BASE_TYPES

  private:
    size_type si;

  public:
    inline prolong_array_generator(array_type &x) : base(x), si((this->array()).size()) {}
    inline prolong_array_generator(array_type &x, const size_type &s) : base(x), si(s) {}

    inline const_reference operator[](const index_type &i) const { if ((this->array()).withinbounds(i)) return (this->array())[i]; else return value_type(0); }

    inline size_type size() const { return si; }
    inline void resize(const size_type &s) const { si=s; }
};

template<class A>
struct prolongArray
{
  typedef A array_type;
  typedef prolong_array_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename prolongArray<const Vector<G> >::self prolong(const Vector<G> &X) { return prolongArray<const Vector<G> >::self(X); }
template<class G> inline typename prolongArray<const Matrix<G> >::self prolong(const Matrix<G> &X) { return prolongArray<const Matrix<G> >::self(X); }
template<class G> inline typename prolongArray<const Vector<G> >::self prolong(const Vector<G> &X, const typename Vector<G>::size_type &s) { return prolongArray<const Vector<G> >::self(X,s); }
template<class G> inline typename prolongArray<const Matrix<G> >::self prolong(const Matrix<G> &X, const typename Vector<G>::size_type &s) { return prolongArray<const Matrix<G> >::self(X,s); }



template<class A>
class symmetric_prolong_array_generator : public array_generator<A>
{
  public:
    typedef symmetric_prolong_array_generator self; 
    typedef array_generator<A> base;
    ARRAY_BASE_TYPES

    template<class A2> struct iterator_rebind       { typedef typename shift_iterator_rebind<      A2>::other other; };
    template<class A2> struct const_iterator_rebind { typedef typename shift_iterator_rebind<const A2>::other other; };

  public:
    inline symmetric_prolong_array_generator(array_type &x) : base(x) {}

    // ne marche que pour les vecteurs
    inline reference       operator[](const index_type &i)       { if ((this->array()).withinbounds(i)) return (this->array())[i]; else return (this->array())[-i]; }
    inline const_reference operator[](const index_type &i) const { if ((this->array()).withinbounds(i)) return (this->array())[i]; else return (this->array())[-i]; }
    
    inline void resize(const size_type &s) { (this->array()).resize(s); }
    inline size_type size() const { return 2*(this->array()).size()-1; }
    inline index_type lower_bound() const { return -(this->array()).upper_bound(); }
};

template<class A>
struct symmetricProlongArray
{
  typedef A array_type;
  typedef symmetric_prolong_array_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename symmetricProlongArray<      Vector<G> >::self symmetric_prolong(      Vector<G> &X) { return symmetricProlongArray<      Vector<G> >::self(X); }
template<class G> inline typename symmetricProlongArray<const Vector<G> >::self symmetric_prolong(const Vector<G> &X) { return symmetricProlongArray<const Vector<G> >::self(X); }


template<class A>
class periodize_array_generator : public array_generator<A>
{
  public:
    typedef periodize_array_generator self; 
    typedef array_generator<A> base;
    ARRAY_BASE_TYPES

  protected:
    using base::X;

  public:
    inline periodize_array_generator(array_type &x) : base(x) {}
    
    using base::array;
    using base::size;
        
    inline reference       operator[](const index_type &i)       { return array()[mod(i-X.lower_bound(),size()) + X.lower_bound()]; }
    inline const_reference operator[](const index_type &i) const { return array()[mod(i-X.lower_bound(),size()) + X.lower_bound()]; }
    
    //void resize(const size_type &s) { X.resize(s); }
    //size_type size() const { return X.size(); }
    //index_type lower_bound() const { return index_type(0); }
};

template<class A>
struct periodizeArray
{
  typedef A array_type;
  typedef periodize_array_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename periodizeArray<      Vector<G> >::self periodize(      Vector<G> &X) { return typename periodizeArray<      Vector<G> >::self(X); }
template<class G> inline typename periodizeArray<const Vector<G> >::self periodize(const Vector<G> &X) { return typename periodizeArray<const Vector<G> >::self(X); }
template<class G> inline typename periodizeArray<      Matrix<G> >::self periodize(      Matrix<G> &X) { return typename periodizeArray<      Matrix<G> >::self(X); }
template<class G> inline typename periodizeArray<const Matrix<G> >::self periodize(const Matrix<G> &X) { return typename periodizeArray<const Matrix<G> >::self(X); }


template<class T,int Copy=0>
class extract_generator : public array_generator<T,Copy>
{
  public:
    typedef extract_generator self; 
    typedef array_generator<T,Copy> base;
    ARRAY_BASE_TYPES

    template<class A2> struct iterator_rebind       { typedef typename shift_iterator_rebind<      A2>::other other; };
    template<class A2> struct const_iterator_rebind { typedef typename shift_iterator_rebind<const A2>::other other; };

  private:
    index_type ind;
    size_type si;

  public:
//    inline extract_generator(array_type &x, const index_type &p) : base(x), ind(p), si(x.size()-(p-x.lower_bound())) {}
    inline extract_generator(array_type &x, const index_type &p) : base(x), ind(p), si(x.size()) {}
    inline extract_generator(array_type &x, const index_type &p, const size_type &s) : base(x), ind(p), si(s) {}

    inline extract_generator(const index_type &p) : base(), ind(p), si(this->array().size()) {}
    template<class A> inline extract_generator(const A &a, const index_type &p) : base(a), ind(p), si(this->array().size()) {}
    template<class A,class B> inline extract_generator(const A &a, const B &b, const index_type &p) : base(a,b), ind(p), si(this->array().size()) {}
    template<class A,class B,class C> inline extract_generator(const A &a, const B &b, const C &c, const index_type &p) : base(a,b,c), ind(p), si(this->array().size()) {}

    template<class T2,int C2> inline extract_generator(const        extract_generator<T2,C2>   &x) : base(x.            array()), ind(x.            lower_bound()), si(x.            size()) {}
    template<class T2,int C2> inline extract_generator(const Vector<extract_generator<T2,C2> > &x) : base(x.generator().array()), ind(x.generator().lower_bound()), si(x.generator().size()) {}
    template<class T2,int C2> inline extract_generator(const Matrix<extract_generator<T2,C2> > &x) : base(x.generator().array()), ind(x.generator().lower_bound()), si(x.generator().size()) {}

    inline index_type lower_bound() const { return ind; }
    inline size_type size() const { return si; }
     
    inline reference       operator[](const index_type &i)       { return (this->array())[i]; }
    inline const_reference operator[](const index_type &i) const { return (this->array())[i]; }
    
    inline void resize(const size_type &s) { si=s; }
};

template<class T,int Copy=0>
struct extractArray : public array_traits<T>
{
  typedef array_traits<T> base;
  ARRAY_BASE_TYPES
  typedef extract_generator<array_type,Copy> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
//Summary: Gives a part of an array
//Remarks:
//  This function differs from the function /sub/,
//  so that each element of the extracted array has the same index position of its origin array.
//  Because of this, the extracted array has to be handled carefully.
//  In some cases the use of /extract/ can be quicker than /sub/,
//  because the library does not recalculate the indexes.
//Arguments:
//  X - The array to cut out
//  p - Start position of the part
//  pi - Vertical start position of the part
//  pj - Horizontal start position of the part
//  s - Size of the part
//  si - Height of the part
//  sj - Width of the part
//Return: An array representing the part.
//See: ^sub^
//Example:
//  DenseVector<int>::self X(5,"1 2 3 4 5");
//  extract(X,2,3)[2] = 0; // X=[1 2 0 4 5];
template<class G> inline typename extractArray<      Vector<G> >::self extract(      Vector<G> &X)                                                                                  { return typename extractArray<      Vector<G> >::self(X,X.lower_bound(),X.size()); }
//{unsecret}
template<class G> inline typename extractArray<      Matrix<G> >::self extract(      Matrix<G> &X)                                                                                  { return typename extractArray<      Matrix<G> >::self(X,X.lower_bound(),X.size()); }
//{unsecret}
template<class G> inline typename extractArray<const Vector<G> >::self extract(const Vector<G> &X)                                                                                  { return typename extractArray<const Vector<G> >::self(X,X.lower_bound(),X.size()); }
//{unsecret}
template<class G> inline typename extractArray<const Matrix<G> >::self extract(const Matrix<G> &X)                                                                                  { return typename extractArray<const Matrix<G> >::self(X,X.lower_bound(),X.size()); }
//{unsecret}
template<class G> inline typename extractArray<      Vector<G> >::self extract(      Vector<G> &X, const typename Vector<G>::index_type &p)                                         { return typename extractArray<      Vector<G> >::self(X,p              ,X.size()); }
//{unsecret}
template<class G> inline typename extractArray<      Matrix<G> >::self extract(      Matrix<G> &X, const typename Matrix<G>::index_type &p)                                         { return typename extractArray<      Matrix<G> >::self(X,p              ,X.size()); }
//{unsecret}
template<class G> inline typename extractArray<const Vector<G> >::self extract(const Vector<G> &X, const typename Vector<G>::index_type &p)                                         { return typename extractArray<const Vector<G> >::self(X,p              ,X.size()); }
//{unsecret}
template<class G> inline typename extractArray<const Matrix<G> >::self extract(const Matrix<G> &X, const typename Matrix<G>::index_type &p)                                         { return typename extractArray<const Matrix<G> >::self(X,p              ,X.size()); }

//{unsecret}
template<class G> inline typename extractArray<      Vector<G> >::self extract(      Vector<G> &X, const typename Vector<G>::index_type &p, const typename Vector<G>::size_type &s) { return typename extractArray<      Vector<G> >::self(X,p              ,s       ); }
//{unsecret}
template<class G> inline typename extractArray<      Matrix<G> >::self extract(      Matrix<G> &X, const typename Matrix<G>::index_type &p, const typename Matrix<G>::size_type &s) { return typename extractArray<      Matrix<G> >::self(X,p              ,s       ); }
//{unsecret}
template<class G> inline typename extractArray<const Vector<G> >::self extract(const Vector<G> &X, const typename Vector<G>::index_type &p, const typename Vector<G>::size_type &s) { return typename extractArray<const Vector<G> >::self(X,p              ,s       ); }
//{unsecret}
template<class G> inline typename extractArray<const Matrix<G> >::self extract(const Matrix<G> &X, const typename Matrix<G>::index_type &p, const typename Matrix<G>::size_type &s) { return typename extractArray<const Matrix<G> >::self(X,p              ,s       ); }

//{unsecret}
template<class G> inline typename extractArray<      Matrix<G> >::self extract(      Matrix<G> &X, typename Matrix<G>::int_type pi, typename Matrix<G>::int_type pj)                { return typename extractArray<      Matrix<G> >::self(X,typename Vector<G>::index_type(pi,pj),X.size()); }
//{unsecret}
template<class G> inline typename extractArray<const Matrix<G> >::self extract(const Matrix<G> &X, typename Matrix<G>::int_type pi, typename Matrix<G>::int_type pj)                { return typename extractArray<const Matrix<G> >::self(X,typename Vector<G>::index_type(pi,pj),X.size()); }
//{unsecret}
template<class G> inline typename extractArray<      Matrix<G> >::self extract(      Matrix<G> &X, typename Matrix<G>::int_type pi, typename Matrix<G>::int_type pj, typename Matrix<G>::size_type s) { return typename extractArray<      Matrix<G> >::self(X,typename Vector<G>::index_type(pi,pj),s); }
//{unsecret}
template<class G> inline typename extractArray<const Matrix<G> >::self extract(const Matrix<G> &X, typename Matrix<G>::int_type pi, typename Matrix<G>::int_type pj, typename Matrix<G>::size_type s) { return typename extractArray<const Matrix<G> >::self(X,typename Vector<G>::index_type(pi,pj),s); }
//{unsecret}
template<class G> inline typename extractArray<      Matrix<G> >::self extract(      Matrix<G> &X, typename Matrix<G>::int_type pi, typename Matrix<G>::int_type pj, typename Matrix<G>::int_type si, typename Matrix<G>::int_type sj) { return typename extractArray<      Matrix<G> >::self(X,typename Vector<G>::index_type(pi,pj),typename Vector<G>::size_type(si,sj) ); }
//{unsecret}
template<class G> inline typename extractArray<const Matrix<G> >::self extract(const Matrix<G> &X, typename Matrix<G>::int_type pi, typename Matrix<G>::int_type pj, typename Matrix<G>::int_type si, typename Matrix<G>::int_type sj) { return typename extractArray<const Matrix<G> >::self(X,typename Vector<G>::index_type(pi,pj),typename Vector<G>::size_type(si,sj) ); }

template<class T>
class extract_function : public unary_value_function<T, typename extractArray<T,1>::self>
{
  public:
    typedef unary_value_function<T, typename extractArray<T,1>::self> base;
    UNARY_FUNCTION_BASE_TYPES
    
    typedef T array_type;
    typedef typename array_type::index_type index_type;
    typedef typename array_type::size_type  size_type;
    
  private:
    index_type p1;
    size_type  s1;
    
  public:
    inline extract_function(const index_type &p, const size_type &s) : p1(p), s1(s) {}
    
    inline const_reference operator()(argument_type &x) const { return const_reference(x,p1,s1); }
};

template<class G>
typename UnaryFunctionArray<const Vector<G>,const extract_function<const typename Vector<G>::value_type> >::self
inline ele_extract(const Vector<G> &X, const typename Vector<G>::value_type::index_type &p, const typename Vector<G>::value_type::size_type &s)
{
  return apply(X,extract_function<const typename Vector<G>::value_type>(p,s));
}

template<class G>
inline typename UnaryFunctionArray<const Matrix<G>,const extract_function<const typename Matrix<G>::value_type> >::self
ele_extract(const Matrix<G> &X, const typename Matrix<G>::value_type::index_type &p, const typename Matrix<G>::value_type::size_type &s)
{
  return apply(X,extract_function<const typename Matrix<G>::value_type>(p,s));
}

template<class G>
inline typename UnaryFunctionArray<const Matrix<G>,const extract_function<const typename Matrix<G>::value_type> >::self
ele_extract(const Matrix<G> &X, const typename Matrix<G>::value_type::index_type::value_type &p1, const typename Matrix<G>::value_type::index_type::value_type &p2, const typename Matrix<G>::value_type::size_type::value_type &m, const typename Matrix<G>::value_type::size_type::value_type &n)
{
  return ele_extract(X,typename Matrix<G>::value_type::index_type(p1,p2),typename Matrix<G>::value_type::size_type(m,n));
}


template<class T,int Copy=0>
class resize_array_generator : public array_generator<T,Copy>
{
  public:
    typedef resize_array_generator self; 
    typedef array_generator<T,Copy> base;
    ARRAY_BASE_TYPES

  private:
    size_type si;

  public:
    resize_array_generator(array_type &x, const size_type &s) : base(x), si(s) {}
    template<class A> resize_array_generator(const A &a, const size_type &s) : base(a), si(s) {}
    
    template<class T2,int C2> resize_array_generator(const resize_array_generator<T2,C2> &x) : base(x.array()), si(x.size()) {}
    template<class T2,int C2> resize_array_generator(const Vector<resize_array_generator<T2,C2> > &x) : base(x.generator().array()), si(x.generator().size()) {}
    template<class T2,int C2> resize_array_generator(const Matrix<resize_array_generator<T2,C2> > &x) : base(x.generator().array()), si(x.generator().size()) {}

    inline reference       operator[](const index_type &i)       { return (this->array())[i]; }
    inline const_reference operator[](const index_type &i) const { return (this->array())[i]; }

    inline size_type size() const { return si; }
         
    inline void resize(const size_type &s) { si=s; }
};

template<class T,int Copy=0>
struct resizeArray : public array_traits<T>
{
  typedef array_traits<T> base;
  ARRAY_BASE_TYPES
  typedef resize_array_generator<array_type,Copy> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename resizeArray<      Vector<G> >::self resize(      Vector<G> &X, const typename Vector<G>::size_type &s) { return typename resizeArray<      Vector<G> >::self(X,s); }
template<class G> inline typename resizeArray<      Matrix<G> >::self resize(      Matrix<G> &X, const typename Matrix<G>::size_type &s) { return typename resizeArray<      Matrix<G> >::self(X,s); }
template<class G> inline typename resizeArray<const Vector<G> >::self resize(const Vector<G> &X, const typename Vector<G>::size_type &s) { return typename resizeArray<const Vector<G> >::self(X,s); }
template<class G> inline typename resizeArray<const Matrix<G> >::self resize(const Matrix<G> &X, const typename Matrix<G>::size_type &s) { return typename resizeArray<const Matrix<G> >::self(X,s); }
template<class G> inline typename resizeArray<      Matrix<G> >::self resize(      Matrix<G> &X, typename Matrix<G>::int_type m, typename Matrix<G>::int_type n) { return typename resizeArray<      Matrix<G> >::self(X,Matrix<G>::size_type(m,n)); }
template<class G> inline typename resizeArray<const Matrix<G> >::self resize(const Matrix<G> &X, typename Matrix<G>::int_type m, typename Matrix<G>::int_type n) { return typename resizeArray<const Matrix<G> >::self(X,Matrix<G>::size_type(m,n)); }


template<int D, class F>
class function_generated_array_generator : public reference_array_generator<D,typename F::value_type, typename F::reference, typename F::const_reference>
{
  public:
    typedef function_generated_array_generator self;
    typedef reference_array_generator<D,typename F::value_type, typename F::reference, typename F::const_reference> base;
    
    typedef F function_type;
    
    typedef typename base::value_type value_type;
    typedef typename base::reference reference;
    typedef typename base::const_reference const_reference;
   
  protected:
    function_type func;  

  public: 
    inline function_generated_array_generator() : base(), func() { }
    inline function_generated_array_generator(const function_type &f) : base(), func(f) { }

    inline function_generated_array_generator(const self &r) : base(r), func(r.func) {}
    
    template<class F2> inline function_generated_array_generator(const Vector<function_generated_array_generator<D,F2> > &x) : base(x.generator()), func(x.generator().function()) {}
    template<class F2> inline function_generated_array_generator(const Matrix<function_generated_array_generator<D,F2> > &x) : base(x.generator()), func(x.generator().function()) {}
   
    inline function_type       &function()       { return func; }
    inline const function_type &function() const { return func; }
   
    inline reference       operator[](const int i)       { return func(i); }
    inline const_reference operator[](const int i) const { return func(i); }
    
    template<class I> inline reference       operator[](const MatrixIndex<I> &i)       { return func(i.x,i.y); }
    template<class I> inline const_reference operator[](const MatrixIndex<I> &i) const { return func(i.x,i.y); }
};

template<int D,class F>
struct FunctionGeneratedArray
{
  typedef F function_type;
  typedef function_generated_array_generator<D,function_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int D,class F>
struct resizedFunctionGeneratedArray
{
  typedef typename FunctionGeneratedArray<D,F>::self array_type;
  typedef typename resizeArray<array_type,1>::self self;
};

template<int D,class F>
struct extractedFunctionGeneratedArray
{
  typedef typename FunctionGeneratedArray<D,F>::self array_type;
  typedef typename extractArray<array_type,1>::self self;
}; 


template<int D, class F> inline typename FunctionGeneratedArray<D,const F>::self generate_array(const F &f) { return typename FunctionGeneratedArray<D,const F>::self(f); }

template<class F> inline typename FunctionGeneratedArray<1,const F>::self generate_vector(const F &f) { return generate_array<1>(f); }
template<class F> inline typename FunctionGeneratedArray<2,const F>::self generate_matrix(const F &f) { return generate_array<2>(f); }

template<class F> inline typename resizedFunctionGeneratedArray<1,const F>::self generate_vector(int n, const F &f) { return resize(generate_vector(f),n); }
template<class F> inline typename resizedFunctionGeneratedArray<2,const F>::self generate_matrix(const MatrixSize<int> &s, const F &f) { return resize(generate_matrix(f),s); }
template<class F> inline typename resizedFunctionGeneratedArray<2,const F>::self generate_matrix(int m, int n, const F &f) { return generate_matrix(MatrixSize<int>(m,n),f); }

template<class F> inline typename extractedFunctionGeneratedArray<1,const F>::self generate_shift_vector(int p, int n, const F &f) { return extract(FunctionGeneratedArray<1,const F>::self(f),p,n); }
template<class F> inline typename extractedFunctionGeneratedArray<2,const F>::self generate_shift_matrix(const MatrixIndex<int> &p, const MatrixSize<int> &s, const F &f) { return extract(FunctionGeneratedArray<2,const F>::self(f),p,s); }
template<class F> inline typename extractedFunctionGeneratedArray<2,const F>::self generate_shift_matrix(int pi, int pj, int m, int n, const F &f) { return generate_shift_matrix(MatrixIndex<int>(pi,pj),MatrixSize<int>(m,n),f); }


template<class F>
class function_vector_function : public unary_value_function<double, typename FunctionGeneratedArray<1,F>::self>
{
  public:
    typedef function_vector_function self;
    typedef unary_value_function<double, typename FunctionGeneratedArray<1,F>::self> base;
    UNARY_FUNCTION_BASE_TYPES
    
    typedef F function_type;
    
  public:
    function_type func;

  public:
    inline function_vector_function(const function_type &f) : func(f) {}
    template<class A> inline function_vector_function(const A &a) : func(a) {}
        
    inline const_reference operator()(const argument_type &x) const 
    {
      return generate_vector(function_type(func,x));
    }
};

template<class F>
function_vector_function<const F> fun_vec_fun(const F &f)
{
  return function_vector_function<const F>(f);
}


template<class F>
class function_matrix_function : public binary_value_function<double,double,typename FunctionGeneratedArray<2,F>::self>
{
  public:
    typedef function_matrix_function self;
    typedef binary_value_function<double,double,typename FunctionGeneratedArray<2,F>::self> base;
    BINARY_FUNCTION_BASE_TYPES
    
    typedef F function_type;
    
  public:
    function_type func;

  public:
    inline function_matrix_function(const function_type &f) : func(f) {}
    template<class A> inline function_matrix_function(const A &a) : func(a) {}
    
    
    inline const_reference operator()(const first_argument_type &i) const 
    {
      return generate_matrix(function_type(func,i));
    }
        
    inline const_reference operator()(const first_argument_type &j,const second_argument_type &i) const 
    {
      return generate_matrix(function_type(func,i,j));
    }
};

template<class F>
function_matrix_function<const F> fun_mat_fun(const F &f)
{
  return function_matrix_function<const F>(f);
}



template<class T,int C=0>
class isolate_array_generator : public array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference,C>
{
  public:
    typedef isolate_array_generator self;
    typedef array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference,C> base;
    ARRAY_BASE_TYPES
    
  protected:
    typename T::value_type::index_type p;

  public:
    inline isolate_array_generator(array_type &x, const typename T::value_type::index_type &i) : base(x), p(i) {}
    //inline isolate_array_generator(array_type &x, const self &r) : base(x), p(r.p) {}
    inline isolate_array_generator(const self &r) : base(r), p(r.p) {}
    template<class T2,int C2> inline isolate_array_generator(const Vector<isolate_array_generator<T2,C2> > &x) : base(x.generator().array()), p(x.generator().pos()) {}

    inline reference       operator[](const index_type &i)       { return (this->array())[i][pos()]; }
    inline const_reference operator[](const index_type &i) const { return (this->array())[i][pos()]; }
    
    inline typename T::value_type::index_type pos() const { return p; }
};

template<class T,int C=0>
struct IsolateArray : public array_reference_traits<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference>
{
  typedef array_reference_traits<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference> base;
  ARRAY_BASE_TYPES
  typedef isolate_array_generator<array_type,C> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename IsolateArray<      Vector<G> >::self isolate(      Vector<G> & X, const typename Vector<G>::value_type::index_type &i) { return typename IsolateArray<      Vector<G> >::self(X,i); }
template<class G> inline typename IsolateArray<const Vector<G> >::self isolate(const Vector<G> & X, const typename Vector<G>::value_type::index_type &i) { return typename IsolateArray<const Vector<G> >::self(X,i); }
template<class G> inline typename IsolateArray<      Matrix<G> >::self isolate(      Matrix<G> & X, const typename Matrix<G>::value_type::index_type &i) { return typename IsolateArray<      Matrix<G> >::self(X,i); }
template<class G> inline typename IsolateArray<const Matrix<G> >::self isolate(const Matrix<G> & X, const typename Matrix<G>::value_type::index_type &i) { return typename IsolateArray<const Matrix<G> >::self(X,i); }
template<class G> inline typename IsolateArray<      Vector<G> >::self isolate(      Vector<G> & X, typename Vector<G>::value_type::int_type m, typename Vector<G>::value_type::int_type n) { return isolate(X,typename Vector<G>::value_type::index_type(m,n)); }
template<class G> inline typename IsolateArray<const Vector<G> >::self isolate(const Vector<G> & X, typename Vector<G>::value_type::int_type m, typename Vector<G>::value_type::int_type n) { return isolate(X,typename Vector<G>::value_type::index_type(m,n)); }
template<class G> inline typename IsolateArray<      Matrix<G> >::self isolate(      Matrix<G> & X, typename Vector<G>::value_type::int_type m, typename Vector<G>::value_type::int_type n) { return isolate(X,typename Vector<G>::value_type::index_type(m,n)); }
template<class G> inline typename IsolateArray<const Matrix<G> >::self isolate(const Matrix<G> & X, typename Vector<G>::value_type::int_type m, typename Vector<G>::value_type::int_type n) { return isolate(X,typename Vector<G>::value_type::index_type(m,n)); }


template<class T>
class isolate_block_dense_vector_generator : public array_reference_generator<T>
{
  public:
    typedef isolate_block_dense_vector_generator self;
    typedef array_reference_generator<T> base;
    ARRAY_BASE_TYPES
    
  private:
    pointer data; 
    size_type dim;
  
  public:
    inline isolate_block_dense_vector_generator(array_type &x, const index_type &i, const size_type &d) : base(x), data((pointer)&x[i]), dim(d) {}
    inline isolate_block_dense_vector_generator(const self &r) : base(r), data(r.data), dim(r.dim) {}

    inline reference       operator[](const index_type &i)       { return data[stride()*i]; }
    inline const_reference operator[](const index_type &i) const { return data[stride()*i]; }
    
    inline size_type stride() const { return dim; }
    inline size_type size() const { return this->array().size()/stride(); }   
    inline void resize(const size_type &s) { assert(size()==s); }
    
    inline int pos() const { return data - (pointer)&(this->array())[0]; }
};

template<int N,class T> class isolate_tiny_block_dense_vector_generator;
template<      class T,int C> struct IsolateArray<      Vector<block_dense_vector_generator       <      T                                                      ,C> > > { typedef       T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C> struct IsolateArray<const Vector<block_dense_vector_generator       <const T                                                      ,C> > > { typedef const T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C> struct IsolateArray<      Vector<block_array_generator              <      Vector<isolate_block_dense_vector_generator     <  T> >,C> > > { typedef       T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int C> struct IsolateArray<const Vector<block_array_generator              <const Vector<isolate_block_dense_vector_generator     <  T> >,C> > > { typedef const T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,int C> struct IsolateArray<      Vector<block_array_generator              <      Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > > { typedef       T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,int C> struct IsolateArray<const Vector<block_array_generator              <const Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > > { typedef const T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<      class T,int C> inline typename IsolateArray<      Vector<block_dense_vector_generator       <T                                                            ,C> > >::self isolate(      Vector<block_dense_vector_generator       <T                                                            ,C> > &X, const typename Vector<block_dense_vector_generator       <T                                                            ,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<block_dense_vector_generator       <T                                                             ,C> > >::self(X.generator().array(),i,X.generator().stride()); }
template<      class T,int C> inline typename IsolateArray<const Vector<block_dense_vector_generator       <T                                                            ,C> > >::self isolate(const Vector<block_dense_vector_generator       <T                                                            ,C> > &X, const typename Vector<block_dense_vector_generator       <T                                                            ,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<block_dense_vector_generator       <T                                                             ,C> > >::self(X.generator().array(),i,X.generator().stride()); }
template<      class T,int C> inline typename IsolateArray<      Vector<block_array_generator              <      Vector<isolate_block_dense_vector_generator     <  T> >,C> > >::self isolate(      Vector<block_array_generator              <      Vector<isolate_block_dense_vector_generator     <  T> >,C> > &X, const typename Vector<block_array_generator              <      Vector<isolate_block_dense_vector_generator     <  T> >,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<block_array_generator               <      Vector<isolate_block_dense_vector_generator     <  T> >,C> > >::self(X.generator().array().generator().array(), X.generator().array().generator().stride()*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]), X.generator().stride()*X.generator().array().generator().stride() ); }
template<      class T,int C> inline typename IsolateArray<const Vector<block_array_generator              <const Vector<isolate_block_dense_vector_generator     <  T> >,C> > >::self isolate(const Vector<block_array_generator              <const Vector<isolate_block_dense_vector_generator     <  T> >,C> > &X, const typename Vector<block_array_generator              <const Vector<isolate_block_dense_vector_generator     <  T> >,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<block_array_generator               <const Vector<isolate_block_dense_vector_generator     <  T> >,C> > >::self(X.generator().array().generator().array(), X.generator().array().generator().stride()*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]), X.generator().stride()*X.generator().array().generator().stride() ); }
template<int N,class T,int C> inline typename IsolateArray<      Vector<block_array_generator              <      Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > >::self isolate(      Vector<block_array_generator              <      Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > &X, const typename Vector<block_array_generator              <      Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<block_array_generator               <      Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > >::self(X.generator().array().generator().array(), N*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]), X.generator().stride()*N ); }
template<int N,class T,int C> inline typename IsolateArray<const Vector<block_array_generator              <const Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > >::self isolate(const Vector<block_array_generator              <const Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > &X, const typename Vector<block_array_generator              <const Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<block_array_generator               <const Vector<isolate_tiny_block_dense_vector_generator<N,T> >,C> > >::self(X.generator().array().generator().array(), N*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]), X.generator().stride()*N ); }

//template<int N,class T,int C> struct IsolateArray<      Vector<tiny_block_vector_generator<N,      Vector<isolate_block_dense_vector_generator<T> >,C> > > { typedef       T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
//template<int N,class T,int C> struct IsolateArray<const Vector<tiny_block_vector_generator<N,const Vector<isolate_block_dense_vector_generator<T> >,C> > > { typedef const T array_type; typedef isolate_block_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
//template<int N,class T,int C> inline typename IsolateArray<      Vector<tiny_block_vector_generator<N,      Vector<isolate_block_dense_vector_generator<T> >,C> > >::self isolate(      Vector<tiny_block_vector_generator<N,      Vector<isolate_block_dense_vector_generator<T> >,C> > &X, const typename Vector<tiny_block_vector_generator<N,      Vector<isolate_block_dense_vector_generator<T> >,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<tiny_block_vector_generator<N,      Vector<isolate_block_dense_vector_generator<T> >,C> > >::self (X.generator().array().generator().array(), X.generator().array().generator().stride()*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]), N*X.generator().array().generator().stride() ); }
//template<int N,class T,int C> inline typename IsolateArray<const Vector<tiny_block_vector_generator<N,const Vector<isolate_block_dense_vector_generator<T> >,C> > >::self isolate(const Vector<tiny_block_vector_generator<N,const Vector<isolate_block_dense_vector_generator<T> >,C> > &X, const typename Vector<tiny_block_vector_generator<N,const Vector<isolate_block_dense_vector_generator<T> >,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<tiny_block_vector_generator<N,const Vector<isolate_block_dense_vector_generator<T> >,C> > >::self (X.generator().array().generator().array(), X.generator().array().generator().stride()*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]), N*X.generator().array().generator().stride() ); }


template<class T>
class isolate_block_stride_dense_vector_generator : public array_reference_generator<T>
{
  public:
    typedef isolate_block_stride_dense_vector_generator self;
    typedef array_reference_generator<T> base;
    ARRAY_BASE_TYPES
    
  private:
    pointer data; 
    size_type dim;
  
  public:
    inline isolate_block_stride_dense_vector_generator(array_type &x, const index_type &i, const size_type &d) : base(x), data((pointer)&x[i]), dim(x.generator().stride()*d) {}
    inline isolate_block_stride_dense_vector_generator(const self &r) : base(r), data(r.data), dim(r.dim) {}
  
    using base::array;
  
    inline reference       operator[](const index_type &i)       { return data[stride()*i]; }
    inline const_reference operator[](const index_type &i) const { return data[stride()*i]; }
    
    inline size_type stride() const { return dim; }
    inline size_type size() const { return array().size()/(stride()/array().generator().stride()); }   
    inline void resize(const size_type &s) { assert(size()==s); }
};

template<class T,int C> struct IsolateArray<      Vector<block_stride_dense_vector_generator<      T,C> > > { typedef       T array_type; typedef isolate_block_stride_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<class T,int C> struct IsolateArray<const Vector<block_stride_dense_vector_generator<const T,C> > > { typedef const T array_type; typedef isolate_block_stride_dense_vector_generator<array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<class T,int C> inline typename IsolateArray<      Vector<block_stride_dense_vector_generator<T,C> > >::self isolate(      Vector<block_stride_dense_vector_generator<T,C> > &X, const typename Vector<block_stride_dense_vector_generator<T,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<block_stride_dense_vector_generator<T,C> > >::self(X.generator().array(),i,X.generator().stride()); }
template<class T,int C> inline typename IsolateArray<const Vector<block_stride_dense_vector_generator<T,C> > >::self isolate(const Vector<block_stride_dense_vector_generator<T,C> > &X, const typename Vector<block_stride_dense_vector_generator<T,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<block_stride_dense_vector_generator<T,C> > >::self(X.generator().array(),i,X.generator().stride()); }



template<int N,class T>
class isolate_tiny_block_dense_vector_generator : public array_reference_generator<T>
{
  public:
    typedef isolate_tiny_block_dense_vector_generator self;
    typedef array_reference_generator<T> base;
    ARRAY_BASE_TYPES
    
  private:
    pointer data; 
  
  public:
    inline isolate_tiny_block_dense_vector_generator(array_type &x, const index_type &i) : base(x), data((pointer)&x[i]) {}
    inline isolate_tiny_block_dense_vector_generator(const self &r) : base(r), data(r.data) {}

    inline reference       operator[](const index_type &i)       { return data[stride()*i]; }
    inline const_reference operator[](const index_type &i) const { return data[stride()*i]; }
    
    inline size_type stride() const { return N; }
    inline size_type size() const { return this->array().size()/stride(); }   
    inline index_type pos() const { return data - (pointer)&(this->array())[0]; }        
};

#ifdef NDEBUG
template<      class T,int N,int C> struct IsolateArray<      Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> > > { typedef       T array_type; typedef isolate_tiny_block_dense_vector_generator<N  ,array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<      class T,int N,int C> struct IsolateArray<const Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> > > { typedef const T array_type; typedef isolate_tiny_block_dense_vector_generator<N  ,array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,class T,int N,int C> struct IsolateArray<      Vector<tiny_block_vector_generator      <      Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > > { typedef       T array_type; typedef isolate_tiny_block_dense_vector_generator<N*M,array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,class T,int N,int C> struct IsolateArray<const Vector<tiny_block_vector_generator      <const Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > > { typedef const T array_type; typedef isolate_tiny_block_dense_vector_generator<N*M,array_type> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<      class T,int N,int C> inline typename IsolateArray<      Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> > >::self isolate(      Vector<tiny_block_dense_vector_generator<T                                                             ,N,C> > &X, const typename Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> > >::self(X.generator().array(),i); }
template<      class T,int N,int C> inline typename IsolateArray<const Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> > >::self isolate(const Vector<tiny_block_dense_vector_generator<T                                                             ,N,C> > &X, const typename Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<tiny_block_dense_vector_generator<T                                                            ,N,C> > >::self(X.generator().array(),i); }
template<int M,class T,int N,int C> inline typename IsolateArray<      Vector<tiny_block_vector_generator      <      Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > >::self isolate(      Vector<tiny_block_vector_generator       <      Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > &X, const typename Vector<tiny_block_vector_generator      <      Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> >::value_type::index_type &i) { return typename IsolateArray<      Vector<tiny_block_vector_generator      <      Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > >::self(X.generator().array().generator().array(), M*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]) ); }
template<int M,class T,int N,int C> inline typename IsolateArray<const Vector<tiny_block_vector_generator      <const Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > >::self isolate(const Vector<tiny_block_vector_generator       <const Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > &X, const typename Vector<tiny_block_vector_generator      <const Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> >::value_type::index_type &i) { return typename IsolateArray<const Vector<tiny_block_vector_generator      <const Vector<isolate_tiny_block_dense_vector_generator<M,T> >,N,C> > >::self(X.generator().array().generator().array(), M*i+(&X.generator().array()[0]-&X.generator().array().generator().array()[0]) ); }
#endif


template<class A,int Copy=0>
class split_array_generator : public array_reference_generator<A,typename IsolateArray<A>::self,typename IsolateArray<A>::self,typename IsolateArray<const A>::self,Copy>
{
  public:
    typedef split_array_generator self;
    typedef array_reference_generator<A,typename IsolateArray<A>::self,typename IsolateArray<A>::self,typename IsolateArray<const A>::self,Copy> base;

    typedef A array_type;
    typedef typename array_type::value_type::array_category array_category;

    typedef typename array_type::value_type::int_type   int_type;
    typedef typename array_type::value_type::size_type  size_type;
    typedef typename array_type::value_type::index_type index_type;
    typedef typename base::reference reference;
    typedef typename base::const_reference const_reference;
  
  protected:
    using base::X;
    
  public:
    inline split_array_generator(array_type &x) : base(x) {}
    
    using base::array;
    
    inline split_array_generator(const self &x) : base(x) {}
    template<class T2> inline split_array_generator(const split_array_generator<T2> &x) : base(x.array()) {}
    
    inline reference       operator[](const index_type &i)       { return isolate(array(),i); }
    inline const_reference operator[](const index_type &i) const { return isolate(array(),i); }
    
    inline index_type lower_bound() const { return array().front().lower_bound(); }
    inline size_type  size ()       const { return array().front().size       (); }
    inline int_type   nelms()       const { return array().front().nelms      (); }  
    inline void resize(size_type s) { assert(size()==s); }    
};

template<class A,int Copy=0>
struct splitArray : public array_value_traits<A,typename IsolateArray<A>::self>
{
  typedef array_value_traits<A,typename IsolateArray<A>::self> base;
  ARRAY_BASE_TYPES
  typedef split_array_generator<array_type,Copy> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> typename splitArray<      Vector<G> >::self split(      Vector<G> &X) { return typename splitArray<      Vector<G> >::self(X); }
template<class G> typename splitArray<const Vector<G> >::self split(const Vector<G> &X) { return typename splitArray<const Vector<G> >::self(X); }
template<class G> typename splitArray<      Matrix<G> >::self split(      Matrix<G> &X) { return typename splitArray<      Matrix<G> >::self(X); }
template<class G> typename splitArray<const Matrix<G> >::self split(const Matrix<G> &X) { return typename splitArray<const Matrix<G> >::self(X); }

template<class T>
class split_function : public unary_value_function<T, typename splitArray<T,1>::self>
{
  public:
    typedef unary_value_function<T, typename splitArray<T,1>::self> base;
    UNARY_FUNCTION_BASE_TYPES
    
    typedef T array_type;
    typedef typename array_type::index_type index_type;
    typedef typename array_type::size_type  size_type;
        
  public:
    split_function() {}
    
    const_reference operator()(argument_type &x) const { return const_reference(x); }
};

template<class G>
typename UnaryFunctionArray<const Vector<G>,const split_function<const typename Vector<G>::value_type> >::self
inline ele_split(const Vector<G> &X)
{
  return apply(X,split_function<const typename Vector<G>::value_type>());
}

template<class G>
typename UnaryFunctionArray<const Matrix<G>,const split_function<const typename Matrix<G>::value_type> >::self
inline ele_split(const Matrix<G> &X)
{
  return apply(X,split_function<const typename Matrix<G>::value_type>());
}

//}

#endif








