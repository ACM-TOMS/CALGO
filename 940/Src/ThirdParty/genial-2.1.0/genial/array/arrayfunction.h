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

#ifndef ARRAYFUNCTION_H
#define ARRAYFUNCTION_H

#include <functional>
#include <complex>

//namespace genial
//{

//Group=Arrays functions

template<class G> class Vector;
template<class G> class Matrix;
template<class I> class MatrixIndex;
template<class I> class MatrixSize;

#ifdef HIDE_FROM_DOCJET
#else
#define FUNCTION_ARRAY_BASE_TYPES \
  typedef typename base::function_type   function_type;   \
  typedef typename base::array_type      array_type;      \
  typedef typename base::value_type      value_type;      \
  typedef typename base::reference       reference;       \
  typedef typename base::const_reference const_reference; \
  typedef typename base::pointer         pointer;         \
  typedef typename base::const_pointer   const_pointer;   \
  typedef typename base::index_type index_type;           \
  typedef typename base::size_type  size_type;            \
  typedef typename base::int_type   int_type;

#define BINARY_FUNCTION_ARRAY_BASE_TYPES \
  typedef typename base::function_type   function_type;       \
  typedef typename base::first_array_type  first_array_type;  \
  typedef typename base::second_array_type second_array_type; \
  typedef typename base::value_type      value_type;          \
  typedef typename base::reference       reference;           \
  typedef typename base::const_reference const_reference;     \
  typedef typename base::pointer         pointer;             \
  typedef typename base::const_pointer   const_pointer;       \
  typedef typename base::index_type index_type;               \
  typedef typename base::size_type  size_type;                \
  typedef typename base::int_type   int_type;
#endif

template<class T, class F, int Copy=0>
class function_array_generator : public array_reference_generator<T, typename F::value_type, typename F::reference, typename F::const_reference,Copy>
{
  public:
    typedef function_array_generator self;
    typedef array_reference_generator<T, typename F::value_type, typename F::reference, typename F::const_reference,Copy> base;
    ARRAY_BASE_TYPES

    typedef F function_type;

  protected:
    function_type func;

  public:
    inline explicit function_array_generator(array_type &x) : base(x), func() {}
    inline function_array_generator(array_type &x, const function_type &f) : base(x), func(f) { }

    template<class A> inline function_array_generator(array_type &x, const A &a) : base(x), func(a) {}
    template<class A, class B> inline function_array_generator(array_type &x, const A &a, const B &b) : base(x), func(a,b) {}

    template<class A> inline function_array_generator(const A &a) : base(a), func() {}
    template<class A, class B> inline function_array_generator(const A &a, const B &b) : base(a,b), func() {}
    template<class A, class B, class C> inline function_array_generator(const A &a, const B &b, const C &c) : base(a,b,c), func() {}

    inline function_array_generator(const self &r) : base(r.X), func(r.func) {}

    inline function_type       &function()       { return func; }
    inline const function_type &function() const { return func; }

 // attention: si func est un adaptateur d'array, l'adaptation doit copier l'array
 // à moins éventuellement que X[i] renvoye une vrai référence, mais c'est risqué pour la compatibilité future
    inline reference       operator[](const index_type &i)
    {
      //typename array_type::reference r=X[i];
      typename array_traits<array_type>::reference r=(this->array())[i];
      return func(r);
     }
    inline const_reference operator[](const index_type &i) const
    {
      typename array_traits<array_type>::const_reference r=(this->array())[i];
      return func(r);
    }
};

template<class A,class F,int Copy>
struct UnaryFunctionArray : public array_reference_traits<A,typename F::value_type,typename F::reference,typename F::const_reference>
{
  typedef array_reference_traits<A,typename F::value_type,typename F::reference,typename F::const_reference> base;
  ARRAY_BASE_TYPES
  typedef F function_type;
  typedef function_array_generator<array_type,function_type,Copy> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int N,class T,class F,int Copy=0>
class tiny_function_vector_generator : public function_array_generator<T,F,Copy>
{
  public:
    typedef tiny_function_vector_generator self;
    typedef function_array_generator<T,F,Copy> base;
    FUNCTION_ARRAY_BASE_TYPES;

  public:
    inline explicit tiny_function_vector_generator(array_type &x) : base(x) {}
    inline tiny_function_vector_generator(array_type &x, const function_type &f) : base(x,f) { }

    template<class A> inline tiny_function_vector_generator(array_type &x, const A &a) : base(x,a) {}
    template<class A, class B> inline tiny_function_vector_generator(array_type &x, const A &a, const B &b) : base(x,a,b) {}

    template<class A> inline tiny_function_vector_generator(const A &a) : base(a) {}
    template<class A, class B> inline tiny_function_vector_generator(const A &a, const B &b) : base(a,b) {}
    template<class A, class B, class C> inline tiny_function_vector_generator(const A &a, const B &b, const C &c) : base(a,b,c) {}

    inline tiny_function_vector_generator(const self &r) : base(r) {}
};

#ifdef NDEBUG
template<int N,class V,class F,int C> struct UnaryFunctionArray<      Vector<tiny_vector_generator<N,V> >,F,C> : public array_reference_traits<      Vector<tiny_vector_generator<N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<      Vector<tiny_vector_generator<N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_vector_generator<N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V,class F,int C> struct UnaryFunctionArray<const Vector<tiny_vector_generator<N,V> >,F,C> : public array_reference_traits<const Vector<tiny_vector_generator<N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<const Vector<tiny_vector_generator<N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_vector_generator<N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,class F2,int C2,class F,int C> struct UnaryFunctionArray<      Vector<tiny_function_vector_generator<N,T,F2,C2> >,F,C> : public array_reference_traits<      Vector<tiny_function_vector_generator<N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<      Vector<tiny_function_vector_generator<N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_vector_generator<N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T,class F2,int C2,class F,int C> struct UnaryFunctionArray<const Vector<tiny_function_vector_generator<N,T,F2,C2> >,F,C> : public array_reference_traits<const Vector<tiny_function_vector_generator<N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<const Vector<tiny_function_vector_generator<N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_vector_generator<N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
#endif

template<int M, int N,class T,class F,int Copy=0>
class tiny_function_matrix_generator : public function_array_generator<T,F,Copy>
{
  public:
    typedef tiny_function_matrix_generator self;
    typedef function_array_generator<T,F,Copy> base;
    FUNCTION_ARRAY_BASE_TYPES

  public:
    inline explicit tiny_function_matrix_generator(array_type &x) : base(x) {}
    inline tiny_function_matrix_generator(array_type &x, const function_type &f) : base(x,f) { }

    template<class A> inline tiny_function_matrix_generator(array_type &x, const A &a) : base(x,a) {}
    template<class A, class B> inline tiny_function_matrix_generator(array_type &x, const A &a, const B &b) : base(x,a,b) {}

    template<class A> inline tiny_function_matrix_generator(const A &a) : base(a) {}
    template<class A, class B> inline tiny_function_matrix_generator(const A &a, const B &b) : base(a,b) {}
    template<class A, class B, class C> inline tiny_function_matrix_generator(const A &a, const B &b, const C &c) : base(a,b,c) {}

    inline tiny_function_matrix_generator(const self &r) : base(r) {}
};

#ifdef NDEBUG
template<int M,int N,class V,class F,int C> struct UnaryFunctionArray<      Matrix<tiny_matrix_generator<M,N,V> >,F,C> : public array_reference_traits<      Matrix<tiny_matrix_generator<M,N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<      Matrix<tiny_matrix_generator<M,N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_matrix_generator<M,N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V,class F,int C> struct UnaryFunctionArray<const Matrix<tiny_matrix_generator<M,N,V> >,F,C> : public array_reference_traits<const Matrix<tiny_matrix_generator<M,N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<const Matrix<tiny_matrix_generator<M,N,V> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_matrix_generator<M,N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T,class F2,int C2,class F,int C> struct UnaryFunctionArray<      Matrix<tiny_function_matrix_generator<M,N,T,F2,C2> >,F,C> : public array_reference_traits<      Matrix<tiny_function_matrix_generator<M,N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<      Matrix<tiny_function_matrix_generator<M,N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_matrix_generator<M,N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T,class F2,int C2,class F,int C> struct UnaryFunctionArray<const Matrix<tiny_function_matrix_generator<M,N,T,F2,C2> >,F,C> : public array_reference_traits<const Matrix<tiny_function_matrix_generator<M,N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> { typedef array_reference_traits<const Matrix<tiny_function_matrix_generator<M,N,T,F2,C2> >,typename F::value_type,typename F::reference,typename F::const_reference> base; ARRAY_BASE_TYPES typedef F function_type; typedef tiny_function_matrix_generator<M,N,array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
#endif

template<class G,class F> inline typename UnaryFunctionArray<      Vector<G>,      F>::self apply(      Vector<G> &X, const F &f) { return typename UnaryFunctionArray<      Vector<G>,      F>::self(X,f); }
template<class G,class F> inline typename UnaryFunctionArray<const Vector<G>,const F>::self apply(const Vector<G> &X, const F &f) { return typename UnaryFunctionArray<const Vector<G>,const F>::self(X,f); }
template<class G,class F> inline typename UnaryFunctionArray<      Matrix<G>,      F>::self apply(      Matrix<G> &X, const F &f) { return typename UnaryFunctionArray<      Matrix<G>,      F>::self(X,f); }
template<class G,class F> inline typename UnaryFunctionArray<const Matrix<G>,const F>::self apply(const Matrix<G> &X, const F &f) { return typename UnaryFunctionArray<const Matrix<G>,const F>::self(X,f); }

template<class G,class Arg,class Res> inline typename UnaryFunctionArray<const Vector<G>,const pointer_to_unary_function<Arg,Res> >::self _apply(const Vector<G> &X, Res (*f)(Arg)) { return apply(X,ptr_fun(f)); }

template<template<class> class F,class G> inline typename UnaryFunctionArray<const Vector<G>,const F<const typename Vector<G>::const_reference> >::self apply(const Vector<G> &X) { return apply(X,F<const typename Vector<G>::const_reference>()); }
template<template<class> class F,class G> inline typename UnaryFunctionArray<const Matrix<G>,const F<const typename Matrix<G>::const_reference> >::self apply(const Matrix<G> &X) { return apply(X,F<const typename Matrix<G>::const_reference>()); }
//template<template<class> class F,class G> inline typename UnaryFunctionArray<const Vector<G>,const F<typename Vector<G>::value_type> >::self apply(const Vector<G> &X) { return apply(X,F<typename Vector<G>::value_type>()); }
//template<template<class> class F,class G> inline typename UnaryFunctionArray<const Matrix<G>,const F<typename Matrix<G>::value_type> >::self apply(const Matrix<G> &X) { return apply(X,F<typename Matrix<G>::value_type>()); }

//template<template<class,class> class F,class G,class A> inline typename UnaryFunctionArray<const Vector<G>,const F<typename Vector<G>::const_reference,typename Vector<G>::const_value_type> >::self apply(const Vector<G> &X, const A &a) { return apply(X,F<typename Vector<G>::const_reference,typename Vector<G>::const_value_type>(a)); }
//template<template<class,class> class F,class G,class A> inline typename UnaryFunctionArray<const Matrix<G>,const F<typename Matrix<G>::const_reference,typename Vector<G>::const_value_type> >::self apply(const Matrix<G> &X, const A &a) { return apply(X,F<typename Matrix<G>::const_reference,typename Vector<G>::const_value_type>(a)); }

template<class G,class F,class V> inline typename UnaryFunctionArray<const Vector<G>,const binder1st<F> >::self distrib1st(const Vector<G> &X, const F &f, const V &v) { return apply(X, bind1st(f,v) ); }
template<class G,class F,class V> inline typename UnaryFunctionArray<const Vector<G>,const binder2nd<F> >::self distrib2nd(const Vector<G> &X, const F &f, const V &v) { return apply(X, bind2nd(f,v) ); }

template<template<class,class> class F, class G, class V> inline typename UnaryFunctionArray<const Vector<G>,const binder1st<F<const V,const typename Vector<G>::value_type> > >::self distrib1st(const Vector<G> &X,const V &v) { return apply(X,bind1st(F<const V,const typename Vector<G>::value_type>(),v)); }
template<template<class,class> class F, class G, class V> inline typename UnaryFunctionArray<const Vector<G>,const binder2nd<F<const typename Vector<G>::value_type,const V> > >::self distrib2nd(const Vector<G> &X,const V &v) { return apply(X,bind2nd(F<const typename Vector<G>::value_type,const V>(),v)); }
template<template<class,class,class> class F, class Res, class G, class V> inline typename UnaryFunctionArray<const Vector<G>,const binder1st<F<const V,const typename Vector<G>::value_type,Res> > >::self distrib1st(const Vector<G> &X,const V &v) { return apply(X,bind1st(F<const V,const typename Vector<G>::value_type,Res>(),v)); }
template<template<class,class,class> class F, class Res, class G, class V> inline typename UnaryFunctionArray<const Vector<G>,const binder2nd<F<const V,const typename Vector<G>::value_type,Res> > >::self distrib2nd(const Vector<G> &X,const V &v) { return apply(X,bind2nd(F<const V,const typename Vector<G>::value_type,Res>(),v)); }

//template parameter Copy not used yet.
template<class T1, class T2, class F,int Copy=0>
class binary_function_array_generator : public two_arrays_reference_generator<T1,T2,typename F::value_type,typename F::reference,typename F::const_reference>
{
  public:
    typedef binary_function_array_generator self;
    typedef two_arrays_reference_generator<T1,T2,typename F::value_type,typename F::reference,typename F::const_reference> base;
    TWO_ARRAYS_BASE_TYPES

    typedef F function_type;

  protected:
    function_type func;
    using base::X;
    using base::Y;

  public:
    inline binary_function_array_generator(first_array_type &x, second_array_type &y) : base(x,y), func() {}
    template<class A> inline binary_function_array_generator(first_array_type &x, second_array_type &y, const A &a) : base(x,y), func(a) {}
    template<class A,class B> inline binary_function_array_generator(first_array_type &x, second_array_type &y, const A &a, const B &b) : base(x,y), func(a,b) {}

    inline binary_function_array_generator(const self &r) : base(r), func(r.func) {}

    inline reference       operator[](const index_type &p)       { return func((this->first_array())[p],(this->second_array())[p]); }
    inline const_reference operator[](const index_type &p) const { return func((this->first_array())[p],(this->second_array())[p]); }

    inline size_type size() const { return (this->first_array()).size(); }
    inline void resize(const size_type &d) { assert(size()==d); }
};

template<class A1,class A2,class F,int Copy=0>
struct BinaryFunctionArray : public two_arrays_reference_traits<A1,A2,typename F::value_type,typename F::reference,typename F::const_reference>
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef F function_type;
  typedef binary_function_array_generator<first_array_type,second_array_type,function_type,Copy> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int N,class T1, class T2, class F,int Copy=0>
class tiny_binary_function_vector_generator : public binary_function_array_generator<T1,T2,F>
{
  public:
    typedef tiny_binary_function_vector_generator self;
    typedef binary_function_array_generator<T1,T2,F> base;
    BINARY_FUNCTION_ARRAY_BASE_TYPES

  public:
    inline tiny_binary_function_vector_generator(first_array_type &x, second_array_type &y) : base(x,y) { }
    template<class A> inline tiny_binary_function_vector_generator(first_array_type &x, second_array_type &y, const A &a) : base(x,y,a) { }
    template<class A,class B> inline tiny_binary_function_vector_generator(first_array_type &x, second_array_type &y, const A &a, const B &b) : base(x,y,a,b) { }
};

#ifdef NDEBUG
template<int N,class V,class T,class F,int C> struct BinaryFunctionArray<      Vector<tiny_vector_generator<N,V> >,T,F,C> { typedef       Vector<tiny_vector_generator<N,V> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V,class T,class F,int C> struct BinaryFunctionArray<const Vector<tiny_vector_generator<N,V> >,T,F,C> { typedef const Vector<tiny_vector_generator<N,V> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V,class T,class F,int C> struct BinaryFunctionArray<T,      Vector<tiny_vector_generator<N,V> >,F,C> { typedef       Vector<tiny_vector_generator<N,V> > second_array_type; typedef T first_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V,class T,class F,int C> struct BinaryFunctionArray<T,const Vector<tiny_vector_generator<N,V> >,F,C> { typedef const Vector<tiny_vector_generator<N,V> > second_array_type; typedef T first_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<      Vector<tiny_vector_generator<N,V1> >,      Vector<tiny_vector_generator<N,V2> >,F,C> { typedef       Vector<tiny_vector_generator<N,V1> > first_array_type; typedef       Vector<tiny_vector_generator<N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<      Vector<tiny_vector_generator<N,V1> >,const Vector<tiny_vector_generator<N,V2> >,F,C> { typedef       Vector<tiny_vector_generator<N,V1> > first_array_type; typedef const Vector<tiny_vector_generator<N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<const Vector<tiny_vector_generator<N,V1> >,      Vector<tiny_vector_generator<N,V2> >,F,C> { typedef const Vector<tiny_vector_generator<N,V1> > first_array_type; typedef       Vector<tiny_vector_generator<N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<const Vector<tiny_vector_generator<N,V1> >,const Vector<tiny_vector_generator<N,V2> >,F,C> { typedef const Vector<tiny_vector_generator<N,V1> > first_array_type; typedef const Vector<tiny_vector_generator<N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };


template<int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<      Vector<tiny_function_vector_generator<N,T1,F1,C1> >,T,F,C> { typedef       Vector<tiny_function_vector_generator<N,T1,F1,C1> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<const Vector<tiny_function_vector_generator<N,T1,F1,C1> >,T,F,C> { typedef const Vector<tiny_function_vector_generator<N,T1,F1,C1> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<T,      Vector<tiny_function_vector_generator<N,T1,F1,C1> >,F,C> { typedef       Vector<tiny_function_vector_generator<N,T1,F1,C1> > second_array_type; typedef T first_array_type;  typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<T,const Vector<tiny_function_vector_generator<N,T1,F1,C1> >,F,C> { typedef const Vector<tiny_function_vector_generator<N,T1,F1,C1> > second_array_type; typedef T first_array_type;  typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<      Vector<tiny_function_vector_generator<N,T1,F1,C1> >,      Vector<tiny_function_vector_generator<N,T2,F2,C2> >,F,C> { typedef       Vector<tiny_function_vector_generator<N,T1,F1,C1> > first_array_type; typedef       Vector<tiny_function_vector_generator<N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<      Vector<tiny_function_vector_generator<N,T1,F1,C1> >,const Vector<tiny_function_vector_generator<N,T2,F2,C2> >,F,C> { typedef       Vector<tiny_function_vector_generator<N,T1,F1,C1> > first_array_type; typedef const Vector<tiny_function_vector_generator<N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<const Vector<tiny_function_vector_generator<N,T1,F1,C1> >,      Vector<tiny_function_vector_generator<N,T2,F2,C2> >,F,C> { typedef const Vector<tiny_function_vector_generator<N,T1,F1,C1> > first_array_type; typedef       Vector<tiny_function_vector_generator<N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<const Vector<tiny_function_vector_generator<N,T1,F1,C1> >,const Vector<tiny_function_vector_generator<N,T2,F2,C2> >,F,C> { typedef const Vector<tiny_function_vector_generator<N,T1,F1,C1> > first_array_type; typedef const Vector<tiny_function_vector_generator<N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_vector_generator<N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
#endif

template<int M,int N,class T1, class T2, class F,int Copy=0>
class tiny_binary_function_matrix_generator : public binary_function_array_generator<T1,T2,F>
{
  public:
    typedef tiny_binary_function_matrix_generator self;
    typedef binary_function_array_generator<T1,T2,F> base;
    BINARY_FUNCTION_ARRAY_BASE_TYPES

  public:
    inline tiny_binary_function_matrix_generator(first_array_type &x, second_array_type &y) : base(x,y) { }
    template<class A> inline tiny_binary_function_matrix_generator(first_array_type &x, second_array_type &y, const A &a) : base(x,y,a) { }
    template<class A,class B> inline tiny_binary_function_matrix_generator(first_array_type &x, second_array_type &y, const A &a, const B &b) : base(x,y,a,b) { }
};

#ifdef NDEBUG
template<int M,int N,class V,class T,class F,int C> struct BinaryFunctionArray<      Matrix<tiny_matrix_generator<M,N,V> >,T,F,C> { typedef       Matrix<tiny_matrix_generator<M,N,V> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V,class T,class F,int C> struct BinaryFunctionArray<const Matrix<tiny_matrix_generator<M,N,V> >,T,F,C> { typedef const Matrix<tiny_matrix_generator<M,N,V> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V,class T,class F,int C> struct BinaryFunctionArray<T,      Matrix<tiny_matrix_generator<M,N,V> >,F,C> { typedef       Matrix<tiny_matrix_generator<M,N,V> > second_array_type; typedef T first_array_type;  typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V,class T,class F,int C> struct BinaryFunctionArray<T,const Matrix<tiny_matrix_generator<M,N,V> >,F,C> { typedef const Matrix<tiny_matrix_generator<M,N,V> > second_array_type; typedef T first_array_type;  typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<int M,int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<      Matrix<tiny_matrix_generator<M,N,V1> >,      Matrix<tiny_matrix_generator<M,N,V2> >,F,C> { typedef       Matrix<tiny_matrix_generator<M,N,V1> > first_array_type; typedef       Matrix<tiny_matrix_generator<M,N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<      Matrix<tiny_matrix_generator<M,N,V1> >,const Matrix<tiny_matrix_generator<M,N,V2> >,F,C> { typedef       Matrix<tiny_matrix_generator<M,N,V1> > first_array_type; typedef const Matrix<tiny_matrix_generator<M,N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<const Matrix<tiny_matrix_generator<M,N,V1> >,      Matrix<tiny_matrix_generator<M,N,V2> >,F,C> { typedef const Matrix<tiny_matrix_generator<M,N,V1> > first_array_type; typedef       Matrix<tiny_matrix_generator<M,N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class V1,class V2,class F,int C> struct BinaryFunctionArray<const Matrix<tiny_matrix_generator<M,N,V1> >,const Matrix<tiny_matrix_generator<M,N,V2> >,F,C> { typedef const Matrix<tiny_matrix_generator<M,N,V1> > first_array_type; typedef const Matrix<tiny_matrix_generator<M,N,V2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };


template<int M,int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<      Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,T,F,C> { typedef       Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,T,F,C> { typedef const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > first_array_type;  typedef T second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<T,      Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,F,C> { typedef       Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > second_array_type; typedef T first_array_type;  typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T1,class F1,int C1,class T,class F,int C> struct BinaryFunctionArray<T,const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,F,C> { typedef const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > second_array_type; typedef T first_array_type;  typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };

template<int M,int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<      Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,      Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> >,F,C> { typedef       Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > first_array_type; typedef       Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<      Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,const Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> >,F,C> { typedef       Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > first_array_type; typedef const Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,      Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> >,F,C> { typedef const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > first_array_type; typedef       Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
template<int M,int N,class T1,class F1,int C1,class T2,class F2,int C2,class F,int C> struct BinaryFunctionArray<const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> >,const Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> >,F,C> { typedef const Matrix<tiny_function_matrix_generator<M,N,T1,F1,C1> > first_array_type; typedef const Matrix<tiny_function_matrix_generator<M,N,T2,F2,C2> > second_array_type; typedef F function_type; typedef tiny_binary_function_matrix_generator<M,N,first_array_type,second_array_type,function_type,C> generator_type; typedef typename GeneratorArray<generator_type>::self self; };
#endif

template<class G1,class G2,class F> inline typename BinaryFunctionArray<const Vector<G1>,const Vector<G2>,const F>::self apply(const Vector<G1> &X, const Vector<G2> &Y, const F &f) { return typename BinaryFunctionArray<const Vector<G1>,const Vector<G2>,const F>::self(X,Y,f); }
template<class G1,class G2,class F> inline typename BinaryFunctionArray<const Matrix<G1>,const Matrix<G2>,const F>::self apply(const Matrix<G1> &X, const Matrix<G2> &Y, const F &f) { return typename BinaryFunctionArray<const Matrix<G1>,const Matrix<G2>,const F>::self(X,Y,f); }


template<class T1, class T2, class T3, class F>
class tertiary_function_generator : public three_arrays_reference_generator<T1,T2,T3,typename F::value_type,typename F::reference,typename F::const_reference>
{
  public:
    typedef tertiary_function_generator self;
    typedef three_arrays_reference_generator<T1,T2,T3,typename F::value_type,typename F::reference,typename F::const_reference> base;
    THREE_ARRAYS_BASE_TYPES

    typedef F function_type;

  protected:
    function_type func;
    using base::X;
    using base::Y;
    using base::Z;

  public:
    inline tertiary_function_generator(first_array_type &x, second_array_type &y, third_array_type &z) : base(x,y,z), func() { assert(X.size()==Y.size() && X.size()==Z.size()); }
    template<class A> inline tertiary_function_generator(first_array_type &x, second_array_type &y, third_array_type &z, const A &a) : base(x,y,z), func(a) { assert(X.size()==Y.size() && X.size()==Z.size()); }
    template<class A,class B> inline tertiary_function_generator(first_array_type &x, second_array_type &y, third_array_type &z, const A &a, const B &b) : base(x,y,z), func(a,b) { assert(X.size()==Y.size() && X.size()==Z.size()); }

    inline tertiary_function_generator(const self &r) : base(r), func(r.func) {}

    inline reference       operator[](const index_type &p)       { return func((this->first_array())[p],(this->second_array())[p],(this->third_array())[p]); }
    inline const_reference operator[](const index_type &p) const { return func((this->first_array())[p],(this->second_array())[p],(this->third_array())[p]); }

    inline size_type size() const { return (this->first_array()).size(); }
    inline void resize(const size_type &d) { assert(size()==d); }
};

template<class A1,class A2,class A3,class F>
struct TertiaryFunctionArray
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef A3 third_array_type;
  typedef F function_type;
  typedef tertiary_function_generator<first_array_type,second_array_type,third_array_type,function_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G1,class G2,class G3,class F> typename TertiaryFunctionArray<const Vector<G1>,const Vector<G2>,const Vector<G3>,const F>::self apply(const Vector<G1> &X, const Vector<G2> &Y, const Vector<G3> &Z, const F &f) { return typename TertiaryFunctionArray<const Vector<G1>,const Vector<G2>,const Vector<G3>,const F>::self(X,Y,Z,f); }
template<class G1,class G2,class G3,class F> typename TertiaryFunctionArray<const Matrix<G1>,const Matrix<G2>,const Vector<G3>,const F>::self apply(const Matrix<G1> &X, const Matrix<G2> &Y, const Vector<G3> &Z, const F &f) { return typename TertiaryFunctionArray<const Matrix<G1>,const Matrix<G2>,const Vector<G3>,const F>::self(X,Y,Z,f); }


template<class T, class F>
class predicate_generator : public array_value_generator<T, typename F::result_type>
{
  public:
    typedef predicate_generator self;
    typedef array_value_generator<T, typename F::result_type> base;
    ARRAY_BASE_TYPES

    typedef F function_type;

  protected:
    function_type func;

  public:
    explicit predicate_generator(array_type &x) : base(x), func() {}
    template<class A> predicate_generator(array_type &x, const A &a) : base(x), func(a) {}
    template<class A, class B> predicate_generator(array_type &x, const A &a, const B &b) : base(x), func(a,b) {}

    predicate_generator(const self &r) : base(r.X), func(r.func) {}

    const_reference operator[](const index_type &i) const { return func((this->array())[i]); }
};

template<class T1, class T2, class F>
class binary_predicate_generator : public two_arrays_value_generator<T1,T2,typename F::result_type>
{
  public:
    typedef binary_predicate_generator self;
    typedef two_arrays_value_generator<T1,T2,typename F::result_type> base;
    TWO_ARRAYS_BASE_TYPES

    typedef F function_type;

  protected:
    function_type func;
    using base::X;
    using base::Y;

  public:
    binary_predicate_generator(first_array_type &x, second_array_type &y) : base(x,y), func() { assert(X.size()==Y.size()); }
    template<class A> binary_predicate_generator(first_array_type &x, second_array_type &y, const A &a) : base(x,y), func(a) { assert(X.size()==Y.size()); }
    template<class A,class B> binary_predicate_generator(first_array_type &x, second_array_type &y, const A &a, const B &b) : base(x,y), func(a,b) { assert(X.size()==Y.size()); }

    binary_predicate_generator(const self &r) : base(r), func(r.func) {}

    const_reference operator[](const index_type &p) const { return func((this->first_array())[p],(this->second_array())[p]); }

    size_type size() const { return (this->first_array()).size(); }
    void resize(const size_type &d) { assert(size()==d); }
};


template<class Res>
struct value_cast_function : public unary_value_function<Res,Res>
{
  typedef unary_value_function<Res,Res> base;
  UNARY_FUNCTION_BASE_TYPES
  template<class V> result_type operator()(const V &x) const { return result_type(x); }
};

template<class T, class V,int C=0>
struct valuecastArray
{
  typedef T array_type;
  typedef V value_type;
  typedef value_cast_function<value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type,C>::self self;
};

//{unsecret}
//Summary: Casts each element
//Example:
//Return: An array representing the array with casted elements
//example:
//  DenseVector<float>::self X(3, "1.1 2.2 3.3");
//  cout << value_cast<int>(X)[0] << endl; // 1
//example:
//  DenseVector<int>::self X(4, "1 2 3 4");
//  DenseVector<float>::self Y=inv(value_cast<float>(X));
//
//  RGBImage Y(512,512,RGBA(12,54,67));
//  cout << value_cast<YUV>(Y)(0,0) << endl;
template<class V,class G> inline typename valuecastArray<const Vector<G>,V>::self value_cast(const Vector<G> &X) { return typename valuecastArray<const Vector<G>,V>::self(X); }
//{unsecret}
template<class V,class G> inline typename valuecastArray<const Matrix<G>,V>::self value_cast(const Matrix<G> &X) { return typename valuecastArray<const Matrix<G>,V>::self(X); }



template<int M,int N,class V> class tiny_matrix_generator;

template<class V,class A> inline const Vector<dense_vector_generator<V,A> > &value_cast(const Vector<dense_vector_generator<V,A> > &X) { return X; }
template<class V,class A> inline const Matrix<dense_matrix_generator<V,A> > &value_cast(const Matrix<dense_matrix_generator<V,A> > &X) { return X; }
template<int N,      class V> inline const Vector<tiny_vector_generator <N,  V> > &value_cast(const Vector<tiny_vector_generator <N,  V> > &X) { return X; }
template<int M,int N,class V> inline const Matrix<tiny_matrix_generator <M,N,V> > &value_cast(const Matrix<tiny_matrix_generator <M,N,V> > &X) { return X; }


#ifdef HIDE_FROM_DOCJET
#else

#define DECLARE_FUNCTION(FUNC,F)                                                                           \
template<class Arg, class Res=Arg>                                                                         \
struct F : public unary_value_function<Arg,Res>                                                            \
{                                                                                                          \
  typedef unary_value_function<Arg,Res> base; UNARY_FUNCTION_BASE_TYPES                                    \
  template<class V> struct rebind { typedef F<V> other; };                                                 \
  template<class V> inline typename F<V>::reference       operator()(const V &x)       { return FUNC(x); } \
  template<class V> inline typename F<V>::const_reference operator()(const V &x) const { return FUNC(x); } \
};

#define DECLARE_REFERENCE_FUNCTION(FUNC,F)                                                      \
template<class Arg, class Res=Arg, class Ref=Res &, class ConstRef=const Res&>                  \
struct F : public unary_reference_function<Arg,Res,Ref,ConstRef>                                \
{                                                                                               \
  typedef unary_reference_function<Arg,Res,Ref,ConstRef> base; UNARY_FUNCTION_BASE_TYPES        \
  template<class V> inline reference       operator()(      V &x)       { return FUNC(x); }     \
  template<class V> inline const_reference operator()(const V &x) const { return FUNC(x); }     \
};

#define DECLARE_FUNCTION2(FUNC,F)                                                                                                   \
template<class Arg1,class Arg2=Arg1,class Res=PROMOTE2(Arg1,Arg2)>                                                                  \
struct F : public binary_value_function<Arg1,Arg2,Res>                                                                              \
{                                                                                                                                   \
  typedef binary_value_function<Arg1,Arg2,Res> base; BINARY_FUNCTION_BASE_TYPES                                                     \
  template<class V1,class V2> inline typename F<V1,V2>::reference       operator()(const V1 &x1,const V2 &x2)       { typename F<V1,V2>::reference       r; r = FUNC(x1,x2); return r; } \
  template<class V1,class V2> inline typename F<V1,V2>::const_reference operator()(const V1 &x1,const V2 &x2) const { typename F<V1,V2>::const_reference r; r = FUNC(x1,x2); return r; } \
};

#define DECLARE_OPERATOR(OP,F)                                                                                 \
template<class Arg, class Res=Arg>                                                                             \
struct F : public unary_value_function<Arg,Res>                                                                \
{                                                                                                              \
  typedef unary_value_function<Arg,Res> base; UNARY_FUNCTION_BASE_TYPES                                        \
  template<class V> inline const_reference operator()(const V &x) const { result_type r; r = OP x; return r; } \
};

#define DECLARE_OPERATOR2(OP,F)                                                                                                        \
template<class Arg1,class Arg2=Arg1,class Res=PROMOTE2(Arg1,Arg2)>                                                                     \
struct F : public binary_value_function<Arg1,Arg2,Res>                                                                                 \
{                                                                                                                                      \
  typedef binary_value_function<Arg1,Arg2,Res>  base; BINARY_FUNCTION_BASE_TYPES                                                       \
  template<class V1,class V2> inline result_type operator()(const V1 &x1,const V2 &x2) const { result_type r; r = x1 OP x2; return r; }\
};


#define DECLARE_FUNCTION_ARRAY(CLASS_NAME,F)                                        \
template<class T, class V=typename F<typename T::const_value_type>::result_type>    \
struct CLASS_NAME : public array_value_traits<T,V>                                  \
{                                                                                   \
  typedef array_value_traits<T,V> base; ARRAY_BASE_TYPES                            \
  typedef F<typename array_type::const_value_type> function_type;                   \
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;         \
};

#define DECLARE_FUNCTION2_ARRAY(CLASS_NAME,F)                                                                               \
template<class T1, class T2, class V=typename F<typename T1::const_value_type,typename T2::const_value_type>::result_type>  \
struct CLASS_NAME : public two_arrays_value_traits<T1,T2,V>                                                                 \
{                                                                                                                           \
  typedef two_arrays_value_traits<T1,T2,V>  base; TWO_ARRAYS_BASE_TYPES                            \
  typedef F<typename T1::const_value_type,typename T2::const_value_type> function_type;                                     \
  typedef typename BinaryFunctionArray<first_array_type,second_array_type,function_type>::self self;                        \
};

#define DECLARE_BIND1ST_ARRAY(CLASS_NAME,F) \
template<class T, class Arg2=typename T::const_value_type, class V=typename F<typename T::const_value_type,Arg2>::result_type> \
struct CLASS_NAME : public array_value_traits<T,V>                          \
{                                                                           \
  typedef array_value_traits<T,V> base;                                     \
  typedef typename base::array_type array_type;                             \
  typedef typename base::value_type value_type;                             \
  typedef typename T::value_type argument_type;                             \
  typedef Arg2 second_param_type;                                           \
  typedef F<argument_type, second_param_type, value_type> function_type;    \
  typedef reference_binder1st<function_type> binder_type;                   \
  typedef typename UnaryFunctionArray<array_type,binder_type>::self self;   \
};

#define DECLARE_BIND2ND_ARRAY(CLASS_NAME,F) \
template<class T, class Arg2=typename T::const_value_type, class V=typename F<typename T::const_value_type,Arg2>::result_type> \
struct CLASS_NAME : public array_value_traits<T,V>                          \
{                                                                           \
  typedef array_value_traits<T,V> base;                                     \
  typedef typename base::array_type array_type;                             \
  typedef typename base::value_type value_type;                             \
  typedef typename T::value_type argument_type;                             \
  typedef Arg2 second_param_type;                                           \
  typedef F<argument_type, second_param_type, value_type> function_type;    \
  typedef reference_binder2nd<function_type> binder_type;                   \
  typedef typename UnaryFunctionArray<array_type,binder_type>::self self;   \
};

#define DECLARE_VECTOR_FUNCTION(FUNC,C) template<class G> inline const typename C<const Vector<G> >::self FUNC(const Vector<G> &A) { return typename C<const Vector<G> >::self(A); }
#define DECLARE_MATRIX_FUNCTION(FUNC,C) template<class G> inline const typename C<const Matrix<G> >::self FUNC(const Matrix<G> &A) { return typename C<const Matrix<G> >::self(A); }
#define DECLARE_VECTOR_OPERATOR(OP,C)   template<class G> inline const typename C<const Vector<G> >::self operator OP(const Vector<G> &A) { return typename C<const Vector<G> >::self(A); }
#define DECLARE_MATRIX_OPERATOR(OP,C)   template<class G> inline const typename C<const Matrix<G> >::self operator OP(const Matrix<G> &A) { return typename C<const Matrix<G> >::self(A); }

#define DECLARE_VECTOR_FUNCTION2(FUNC,C) template<class G1,class G2> inline const typename C<const Vector<G1>, const Vector<G2> >::self FUNC(const Vector<G1> &A, const Vector<G2> &B) { return typename C<const Vector<G1>, const Vector<G2> >::self(A,B); }
#define DECLARE_MATRIX_FUNCTION2(FUNC,C) template<class G1,class G2> inline const typename C<const Matrix<G1>, const Matrix<G2> >::self FUNC(const Matrix<G1> &A, const Matrix<G2> &B) { return typename C<const Matrix<G1>, const Matrix<G2> >::self(A,B); }
#define DECLARE_VECTOR_OPERATOR2(OP,C)   template<class G1,class G2> inline const typename C<const Vector<G1>, const Vector<G2> >::self operator OP(const Vector<G1> &A, const Vector<G2> &B) { return typename C<const Vector<G1>, const Vector<G2> >::self(A,B); }
#define DECLARE_MATRIX_OPERATOR2(OP,C)   template<class G1,class G2> inline const typename C<const Matrix<G1>, const Matrix<G2> >::self operator OP(const Matrix<G1> &A, const Matrix<G2> &B) { return typename C<const Matrix<G1>, const Matrix<G2> >::self(A,B); }

#define DECLARE_ARRAY_FUNCTION(FUNC,C)  DECLARE_VECTOR_FUNCTION (FUNC,C) DECLARE_MATRIX_FUNCTION (FUNC,C) template<int D,class G> inline const typename C<const Array<D,G> >::self FUNC       (const Array<D,G> &A) { return typename C<const Array<D,G> >::self(A); }
#define DECLARE_ARRAY_OPERATOR(OP,C)    DECLARE_VECTOR_OPERATOR (OP,  C) DECLARE_MATRIX_OPERATOR (OP,  C) template<int D,class G> inline const typename C<const Array<D,G> >::self operator OP(const Array<D,G> &A) { return typename C<const Array<D,G> >::self(A); }
#define DECLARE_ARRAY_FUNCTION2(FUNC,C) DECLARE_VECTOR_FUNCTION2(FUNC,C) DECLARE_MATRIX_FUNCTION2(FUNC,C) template<int D,class G1,class G2> inline const typename C<const Array<D,G1>,const Array<D,G2> >::self FUNC       (const Array<D,G1> &A, const Array<D,G2> &B) { return typename C<const Array<D,G1>,const Array<D,G2> >::self(A,B); }
#define DECLARE_ARRAY_OPERATOR2(OP,C)   DECLARE_VECTOR_OPERATOR2(OP,  C) DECLARE_MATRIX_OPERATOR2(OP,  C) template<int D,class G1,class G2> inline const typename C<const Array<D,G1>,const Array<D,G2> >::self operator OP(const Array<D,G1> &A, const Array<D,G2> &B) { return typename C<const Array<D,G1>,const Array<D,G2> >::self(A,B); }


#define DECLARE_NON_CONST_VECTOR_FUNCTION(FUNC,C) template<class G> inline typename C<Vector<G> >::self FUNC       (Vector<G> &A) { return typename C<Vector<G> >::self(A); }
#define DECLARE_NON_CONST_MATRIX_FUNCTION(FUNC,C) template<class G> inline typename C<Matrix<G> >::self FUNC       (Matrix<G> &A) { return typename C<Matrix<G> >::self(A); }
#define DECLARE_NON_CONST_VECTOR_OPERATOR(OP,C)   template<class G> inline typename C<Vector<G> >::self operator OP(Vector<G> &A) { return typename C<Vector<G> >::self(A); }
#define DECLARE_NON_CONST_MATRIX_OPERATOR(OP,C)   template<class G> inline typename C<Matrix<G> >::self operator OP(Matrix<G> &A) { return typename C<Matrix<G> >::self(A); }

#define DECLARE_NON_CONST_ARRAY_FUNCTION(FUNC,C) DECLARE_NON_CONST_VECTOR_FUNCTION (FUNC,C) DECLARE_NON_CONST_MATRIX_FUNCTION (FUNC,C) template<int D,class G> inline typename C<Array<D,G> >::self FUNC       (Array<D,G> &A) { return typename C<Array<D,G> >::self(A); }
#define DECLARE_NON_CONST_ARRAY_OPERATOR(OP,C)   DECLARE_NON_CONST_VECTOR_OPERATOR (OP  ,C) DECLARE_NON_CONST_MATRIX_OPERATOR (OP  ,C) template<int D,class G> inline typename C<Array<D,G> >::self operator OP(Array<D,G> &A) { return typename C<Array<D,G> >::self(A); }



//DocJet RegExpHook
//
//  DECLARE_ARRAY_FUNCTION([ ]*)\(([^,]*),([^)]*)\)
//  template<class G> inline const typename \3<const Vector<G> >::self \2(const Vector<G> &A); \n  template<class G> inline const typename \3<const Matrix<G> >::self \2(const Matrix<G> &A); \n  template<int D,class G> inline const typename \3<const Array<D,G> >::self \2(const Array<D,G> &A);
//
//  DECLARE_ARRAY_OPERATOR([ ]*)\(([^,]*),([^)]*)\)
//  template<class G> inline const typename \3<const Vector<G> >::self operator\2(const Vector<G> &A); \n  template<class G> inline const typename \3<const Matrix<G> >::self operator\2(const Matrix<G> &A); \n  template<int D,class G> inline const typename \3<const Array<D,G> >::self operator\2(const Array<D,G> &A);
//
//  DECLARE_ARRAY_FUNCTION2([ ]*)\(([^,]*),([^)]*)\)
//  template<class G1,class G2> inline const typename \3<const Vector<G1>, const Vector<G2> >::self \2(const Vector<G1> &A, const Vector<G2> &B); \n  template<class G1,class G2> inline const typename \3<const Matrix<G1>, const Matrix<G2> >::self \2(const Matrix<G1> &A, const Matrix<G2> &B); \n  template<int D,class G1,class G2> inline const typename \3<const Array<D,G1>,const Array<D,G2> >::self \2(const Array<D,G1> &A, const Array<D,G2> &B);
//
//  DECLARE_ARRAY_OPERATOR2([ ]*)\(([^,]*),([^)]*)\)
//  template<class G1,class G2> inline const typename \3<const Vector<G1>, const Vector<G2> >::self operator\2(const Vector<G1> &A, const Vector<G2> &B); \n  template<class G1,class G2> inline const typename \3<const Matrix<G1>, const Matrix<G2> >::self operator\2(const Matrix<G1> &A, const Matrix<G2> &B); \n  template<int D,class G1,class G2> inline const typename \3<const Array<D,G1>,const Array<D,G2> >::self operator\2(const Array<D,G1> &A, const Array<D,G2> &B);
//
//  DECLARE_NON_CONST_ARRAY_FUNCTION([ ]*)\(([^,]*),([^)]*)\)
//  template<class G> inline typename \3<Vector<G> >::self \2(Vector<G> &A); \n  template<class G> inline typename \3<Matrix<G> >::self \2(Matrix<G> &A); \n  template<int D,class G> inline typename \3<Array<D,G> >::self \2(Array<D,G> &A);
//
//  DECLARE_NON_CONST_ARRAY_OPERATOR([ ]*)\(([^,]*),([^)]*)\)
//  template<class G> inline typename \3<Vector<G> >::self operator\2(Vector<G> &A); \n  template<class G> inline typename \3<Matrix<G> >::self operator\2(Matrix<G> &A); \n  template<int D,class G> inline typename \3<Array<D,G> >::self operator\2(Array<D,G> &A);
//
#endif

DECLARE_OPERATOR       (-,neg_function)
DECLARE_FUNCTION_ARRAY (negArray,neg_function)
//{unsecret}
//Summary: Negation/Substraction
DECLARE_ARRAY_OPERATOR (-,negArray)


DECLARE_OPERATOR2      (+,add_function)

//template<class Arg1,class Arg2=Arg1,class Res=PROMOTE2(Arg1,Arg2)>                                                                     
//struct add_function : public binary_value_function<Arg1,Arg2,Res>                                                                                 
//{                                                                                                                                      
//  typedef binary_value_function<Arg1,Arg2,Res>  base; BINARY_FUNCTION_BASE_TYPES                                                       
//  template<class V1,class V2> inline result_type operator()(const V1 &x1,const V2 &x2) const 
//  {
//    result_type r; 
//    r = x1 + x2; 
//    return r; 
//  }
//};

DECLARE_FUNCTION2_ARRAY(addArray,add_function)
DECLARE_BIND1ST_ARRAY  (valaddArray,add_function)
DECLARE_BIND2ND_ARRAY  (addvalArray,add_function)

//template<class T, class Arg2=typename T::const_value_type, class V=typename add_function<typename T::const_value_type,Arg2>::result_type> 
//struct addvalArray : public array_value_traits<T,V>                          
//{                                                                           
//  typedef array_value_traits<T,V> base;                                     
//  typedef typename base::array_type array_type;                             
//  typedef typename base::value_type value_type;                             
//  typedef typename T::value_type argument_type;                             
//  typedef Arg2 second_param_type;                                           
//  typedef add_function<argument_type, second_param_type, value_type> function_type;    
//  typedef reference_binder2nd<function_type> binder_type;                   
//  typedef typename UnaryFunctionArray<array_type,binder_type>::self self;   
//};

//{unsecret}
//Summary: Addition
//See: ^operator-^
DECLARE_ARRAY_OPERATOR2(+,addArray)
template<class G> inline typename valaddArray<const Vector<G>,typename Vector<G>::value_type>::self operator+(const typename Vector<G>::value_type &v, const Vector<G> &X) { return typename valaddArray<const Vector<G>,typename Vector<G>::value_type>::self(X,v); }
template<class G> inline typename valaddArray<const Matrix<G>,typename Matrix<G>::value_type>::self operator+(const typename Matrix<G>::value_type &v, const Matrix<G> &X) { return typename valaddArray<const Matrix<G>,typename Matrix<G>::value_type>::self(X,v); }
template<class G> inline typename addvalArray<const Vector<G>,typename Vector<G>::value_type>::self operator+(const Vector<G> &X, const typename Vector<G>::value_type &v) { return typename addvalArray<const Vector<G>,typename Vector<G>::value_type>::self(X,v); }
template<class G> inline typename addvalArray<const Matrix<G>,typename Matrix<G>::value_type>::self operator+(const Matrix<G> &X, const typename Matrix<G>::value_type &v) { return typename addvalArray<const Matrix<G>,typename Matrix<G>::value_type>::self(X,v); }


DECLARE_OPERATOR2      (-,minus_function)
DECLARE_FUNCTION2_ARRAY(minusArray,minus_function)
DECLARE_BIND1ST_ARRAY  (valminusArray,minus_function)
DECLARE_BIND2ND_ARRAY  (minusvalArray,minus_function)
//{unsecret}
DECLARE_ARRAY_OPERATOR2(-,minusArray)
template<class G> inline typename valminusArray<const Vector<G>,typename Vector<G>::value_type>::self operator-(const typename Vector<G>::value_type &v, const Vector<G> &X) { return typename valminusArray<const Vector<G>,typename Vector<G>::value_type>::self(X,v); }
template<class G> inline typename valminusArray<const Matrix<G>,typename Matrix<G>::value_type>::self operator-(const typename Matrix<G>::value_type &v, const Matrix<G> &X) { return typename valminusArray<const Matrix<G>,typename Matrix<G>::value_type>::self(X,v); }
template<class G> inline typename minusvalArray<const Vector<G>,typename Vector<G>::value_type>::self operator-(const Vector<G> &X, const typename Vector<G>::value_type &v) { return typename minusvalArray<const Vector<G>,typename Vector<G>::value_type>::self(X,v); }
template<class G> inline typename minusvalArray<const Matrix<G>,typename Matrix<G>::value_type>::self operator-(const Matrix<G> &X, const typename Matrix<G>::value_type &v) { return typename minusvalArray<const Matrix<G>,typename Matrix<G>::value_type>::self(X,v); }


DECLARE_OPERATOR2      (*, mul_function)
DECLARE_FUNCTION2_ARRAY(elemulArray,mul_function)
DECLARE_BIND1ST_ARRAY  (valmulArray,mul_function)
DECLARE_BIND2ND_ARRAY  (mulvalArray,mul_function)
////{unsecret}
////Summary: Multiplication (element-wise)
////See: ^sqr^, ^mul^, ^operator/^
DECLARE_ARRAY_OPERATOR2(*,elemulArray)
template<class G> inline typename mulvalArray<const Vector<G>,typename Vector<G>::value_type            >::self operator*(const typename Vector<G>::value_type             &x, const Vector<G> &X) { return typename mulvalArray<const Vector<G>,typename Vector<G>::value_type            >::self(X,x); }
template<class G> inline typename mulvalArray<const Matrix<G>,typename Vector<G>::value_type            >::self operator*(const typename Matrix<G>::value_type             &x, const Matrix<G> &X) { return typename mulvalArray<const Matrix<G>,typename Matrix<G>::value_type            >::self(X,x); }
template<class G> inline typename mulvalArray<const Vector<G>,typename Vector<G>::value_type::value_type>::self operator*(const typename Vector<G>::value_type::value_type &x, const Vector<G> &X) { return typename mulvalArray<const Vector<G>,typename Vector<G>::value_type::value_type>::self(X,x); }
template<class G> inline typename mulvalArray<const Matrix<G>,typename Vector<G>::value_type::value_type>::self operator*(const typename Matrix<G>::value_type::value_type &x, const Matrix<G> &X) { return typename mulvalArray<const Matrix<G>,typename Matrix<G>::value_type::value_type>::self(X,x); }
template<class G> inline typename mulvalArray<const Vector<G>,typename Vector<G>::value_type            >::self operator*(const Vector<G> &X, const typename Vector<G>::value_type             &x) { return typename mulvalArray<const Vector<G>,typename Vector<G>::value_type            >::self(X,x); }
template<class G> inline typename mulvalArray<const Matrix<G>,typename Vector<G>::value_type            >::self operator*(const Matrix<G> &X, const typename Matrix<G>::value_type             &x) { return typename mulvalArray<const Matrix<G>,typename Matrix<G>::value_type            >::self(X,x); }
template<class G> inline typename mulvalArray<const Vector<G>,typename Vector<G>::value_type::value_type>::self operator*(const Vector<G> &X, const typename Vector<G>::value_type::value_type &x) { return typename mulvalArray<const Vector<G>,typename Vector<G>::value_type::value_type>::self(X,x); }
template<class G> inline typename mulvalArray<const Matrix<G>,typename Vector<G>::value_type::value_type>::self operator*(const Matrix<G> &X, const typename Matrix<G>::value_type::value_type &x) { return typename mulvalArray<const Matrix<G>,typename Matrix<G>::value_type::value_type>::self(X,x); }


DECLARE_OPERATOR2      (/,div_function)
DECLARE_FUNCTION2_ARRAY(eledivArray,div_function)
DECLARE_BIND2ND_ARRAY  (divvalArray,div_function)
//{unsecret}
//Summary: Division (element-wise)
//See: ^inv^, ^operator*^
DECLARE_ARRAY_OPERATOR2(/,eledivArray)
template<class G> inline typename divvalArray<const Vector<G>,typename Vector<G>::value_type            >::self operator/(const Vector<G> &X, const typename Vector<G>::value_type             &x) { return typename divvalArray<const Vector<G>,typename Vector<G>::value_type            >::self(X,x); }
template<class G> inline typename divvalArray<const Matrix<G>,typename Matrix<G>::value_type            >::self operator/(const Matrix<G> &X, const typename Matrix<G>::value_type             &x) { return typename divvalArray<const Matrix<G>,typename Matrix<G>::value_type            >::self(X,x); }
template<class G> inline typename divvalArray<const Vector<G>,typename Vector<G>::value_type::value_type>::self operator/(const Vector<G> &X, const typename Vector<G>::value_type::value_type &x) { return typename divvalArray<const Vector<G>,typename Vector<G>::value_type::value_type>::self(X,x); }
template<class G> inline typename divvalArray<const Matrix<G>,typename Matrix<G>::value_type::value_type>::self operator/(const Matrix<G> &X, const typename Matrix<G>::value_type::value_type &x) { return typename divvalArray<const Matrix<G>,typename Matrix<G>::value_type::value_type>::self(X,x); }


DECLARE_FUNCTION       (sqr,sqr_function)
DECLARE_FUNCTION_ARRAY (sqrArray,sqr_function)
//{unsecret}
//Summary: Square (element-wise)
//See: ^operator*^
DECLARE_ARRAY_FUNCTION(sqr,sqrArray)

DECLARE_FUNCTION2      (dist,dist_function)
template<class G1,class G2> struct dist_function <Vector<G1>,Vector<G2> > : public binary_value_function<Vector<G1>,Vector<G2>,PROMOTE2(typename Vector<G1>::value_type,typename Vector<G2>::value_type)> { typedef PROMOTE2(typename Vector<G1>::value_type,typename Vector<G2>::value_type) value_type; typedef value_type const_reference; template<class V1,class V2> const_reference operator()(const V1 &x,const V2 &y) const { return dist(x,y); } };
template<class G1,class G2> struct dist_function <Matrix<G1>,Matrix<G2> > : public binary_value_function<Matrix<G1>,Matrix<G2>,PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type)> { typedef PROMOTE2(typename Matrix<G1>::value_type,typename Matrix<G2>::value_type) value_type; typedef value_type const_reference; template<class V1,class V2> const_reference operator()(const V1 &x,const V2 &y) const { return dist(x,y); } };
DECLARE_FUNCTION2_ARRAY(eledistArray,dist_function)
//{unsecret}
//Summary: Distance (element-wise)
//Remarks: The elements should be arrays
//See: ^dist^
DECLARE_ARRAY_FUNCTION2(ele_dist,eledistArray)

template<class Arg1,class Arg2>
struct sad_function : public binary_value_function<Arg1,Arg2,PROMOTE2(Arg1,Arg2)>
{
  typedef Arg1 first_argument;
  typedef Arg2 second_argument;
  typedef PROMOTE2(first_argument,second_argument) result_type;
  template<class V1,class V2> result_type operator()(const V1 &x1,const V2 &x2) const { result_type r; r = sad(x1,x2); return r; }
};
template<class G1,class G2> struct sad_function<Vector<G1>,Vector<G2> > : public binary_value_function<Vector<G1>,Vector<G2>,PROMOTE2(Vector<G1>,Vector<G2>)> {	typedef PROMOTE2(Vector<G1>,Vector<G2>) const_reference; template<class V1,class V2> const_reference operator()(const V1 &x,const V2 &y) const { return sad(x,y); } };
DECLARE_FUNCTION2_ARRAY(elesadArray, sad_function)
DECLARE_ARRAY_FUNCTION2(ele_sad, elesadArray)

DECLARE_FUNCTION       (inv,inv_function)
DECLARE_FUNCTION_ARRAY (invArray,inv_function)
//{unsecret}
//Summary: Inverse (element-wise)
//See: ^operator/^
DECLARE_ARRAY_FUNCTION (inv,invArray)

DECLARE_FUNCTION       (log,log_function)
DECLARE_FUNCTION_ARRAY (logArray,log_function)
//{unsecret}
//Summary: Logarithm (element-wise), base /e/
//See: ^log10^, ^exp^, ^exp10^
DECLARE_ARRAY_FUNCTION (log,logArray)

DECLARE_FUNCTION       (exp,exp_function)
DECLARE_FUNCTION_ARRAY (expArray,exp_function)
//{unsecret}
//Summary: Exponential (element-wise), base /e/
//See: ^log^, ^exp10^, ^log10^
DECLARE_ARRAY_FUNCTION (exp,expArray)

DECLARE_FUNCTION       (log10,log10_function)
DECLARE_FUNCTION_ARRAY (log10Array,log10_function)
//{unsecret}
//Summary: Logarithm (element-wise), base 10
//See: ^exp10^, ^log^, ^exp^
DECLARE_ARRAY_FUNCTION (log10,log10Array)

DECLARE_FUNCTION       (exp10,exp10_function)
DECLARE_FUNCTION_ARRAY (exp10Array,exp10_function)
//{unsecret}
//Summary: Exponential (element-wise), base 10
//See: ^log10^, ^exp^, ^log^
DECLARE_ARRAY_FUNCTION (exp10,exp10Array)

DECLARE_FUNCTION       (norm,norm_function)
DECLARE_FUNCTION2      (norm,norm_binary_function)
DECLARE_FUNCTION_ARRAY (normArray,norm_function)
//{unsecret}
//Summary: Squared magnitude (element-wise)
//See: ^abs^, ^arg^
DECLARE_ARRAY_FUNCTION (norm,normArray)

DECLARE_FUNCTION       (abs,abs_function)
DECLARE_FUNCTION2      (abs,abs_binary_function)
DECLARE_FUNCTION_ARRAY (absArray,abs_function)
//{unsecret}
//Summary: Magnitude (element-wise)
//See: ^norm^, ^arg^
DECLARE_ARRAY_FUNCTION (abs,absArray)

DECLARE_FUNCTION       (conj,conj_function)
DECLARE_FUNCTION_ARRAY (conjArray,conj_function)
//{unsecret}
//Summary: Conjugate (element-wise)
//See: ^imul^
DECLARE_ARRAY_FUNCTION (conj,conjArray)

DECLARE_FUNCTION       (arg,arg_function)
DECLARE_FUNCTION_ARRAY (argArray,arg_function)
//{unsecret}
//Summary: Phase angle (element-wise)
//See: ^abs^, ^norm^
DECLARE_ARRAY_FUNCTION (arg,argArray)

DECLARE_FUNCTION       (imul,imul_function)
DECLARE_FUNCTION_ARRAY (imulArray,imul_function)
//{unsecret}
//Summary: Multiplication by /i/ (element-wise)
//See: ^conj^
DECLARE_ARRAY_FUNCTION (imul,imulArray)

DECLARE_REFERENCE_FUNCTION       (real,real_function)
DECLARE_FUNCTION_ARRAY           (realArray,real_function)
DECLARE_ARRAY_FUNCTION           (real,realArray)
//{unsecret}
//{noAutoLink}
//Summary: Real part (element-wise)
//See: ^imag^
DECLARE_NON_CONST_ARRAY_FUNCTION (real,realArray)

DECLARE_REFERENCE_FUNCTION       (imag,imag_function)
DECLARE_FUNCTION_ARRAY           (imagArray,imag_function)
DECLARE_ARRAY_FUNCTION           (imag,imagArray)
//{unsecret}
//Summary: Imaginary part (element-wise)
//See: ^real^
DECLARE_NON_CONST_ARRAY_FUNCTION (imag,imagArray)


DECLARE_FUNCTION(round,round_function)
DECLARE_FUNCTION_ARRAY (roundArray, round_function)
//{unsecret}
//Summary: Rounded value (element-wise)
DECLARE_ARRAY_FUNCTION (round, roundArray)


template<class Res> struct abs_function <complex<Res> > : public unary_value_function<complex<Res>,Res> { typedef Res const_reference; template<class V> const_reference operator()(const V &x) const { return abs(x); } };
template<class Res> struct arg_function <complex<Res> > : public unary_value_function<complex<Res>,Res> { typedef Res const_reference; template<class V> const_reference operator()(const V &x) const { return arg(x); } };
template<class Res> struct norm_function<complex<Res> > : public unary_value_function<complex<Res>,Res> { typedef Res const_reference; template<class V> const_reference operator()(const V &x) const { return norm(x); } };
template<class Res> struct real_function<complex<Res> > : public unary_reference_function<complex<Res>,Res,Res&,Res> { typedef Res &reference; typedef Res const_reference; template<class V> reference operator()(V &x) { return real_ref(x); } template<class V> const_reference operator()(const V &x) const { return real(x); } };
template<class Res> struct imag_function<complex<Res> > : public unary_reference_function<complex<Res>,Res,Res&,Res> { typedef Res &reference; typedef Res const_reference; template<class V> reference operator()(V &x) { return imag_ref(x); } template<class V> const_reference operator()(const V &x) const { return imag(x); } };

template<class I> struct abs_function <MatrixIndex<I> > : public unary_value_function<MatrixIndex<I>,double> { typedef double const_reference; template<class V> const_reference operator()(const V &x) const { return abs (x); } };
template<class I> struct norm_function<MatrixIndex<I> > : public unary_value_function<MatrixIndex<I>,double> { typedef double const_reference; template<class V> const_reference operator()(const V &x) const { return norm(x); } };

template<class I> struct abs_function <MatrixSize <I> > : public unary_value_function<MatrixSize <I>,double> { typedef double const_reference; template<class V> const_reference operator()(const V &x) const { return abs (x); } };
template<class I> struct norm_function<MatrixSize <I> > : public unary_value_function<MatrixSize <I>,double> { typedef double const_reference; template<class V> const_reference operator()(const V &x) const { return norm(x); } };


                                                                                   
DECLARE_FUNCTION2      (max,max_function)
DECLARE_FUNCTION2      (min,min_function)
DECLARE_FUNCTION2_ARRAY(maxArray,max_function)
DECLARE_FUNCTION2_ARRAY(minArray,min_function)
DECLARE_BIND1ST_ARRAY  (bindmaxArray,max_function)
DECLARE_BIND1ST_ARRAY  (bindminArray,min_function)
//{unsecret}
//Summary: Maximum value (element-wise)
//See: ^ele_min^
DECLARE_ARRAY_FUNCTION2(ele_max,maxArray)
//{unsecret}
//Summary: Minimum value (element-wise)
//See: ^ele_max^
DECLARE_ARRAY_FUNCTION2(ele_min,minArray)
//{unsecret}
template<class V,class G> inline typename bindmaxArray<const Vector<G> >::self ele_max(const V &a, const Vector<G> &A) { return typename bindmaxArray<const Vector<G> >::self(A,a); }
//{unsecret}
template<class V,class G> inline typename bindmaxArray<const Matrix<G> >::self ele_max(const V &a, const Matrix<G> &A) { return typename bindmaxArray<const Matrix<G> >::self(A,a); }
//{unsecret}
template<class V,class G> inline typename bindminArray<const Vector<G> >::self ele_min(const V &a, const Vector<G> &A) { return typename bindminArray<const Vector<G> >::self(A,a); }
//{unsecret}
template<class V,class G> inline typename bindminArray<const Matrix<G> >::self ele_min(const V &a, const Matrix<G> &A) { return typename bindminArray<const Matrix<G> >::self(A,a); }
//{unsecret}
template<class G,class V> inline typename bindmaxArray<const Vector<G> >::self ele_max(const Vector<G> &A, const V &a) { return typename bindmaxArray<const Vector<G> >::self(A,a); }
//{unsecret}
template<class G,class V> inline typename bindmaxArray<const Matrix<G> >::self ele_max(const Matrix<G> &A, const V &a) { return typename bindmaxArray<const Matrix<G> >::self(A,a); }
//{unsecret}
template<class G,class V> inline typename bindminArray<const Vector<G> >::self ele_min(const Vector<G> &A, const V &a) { return typename bindminArray<const Vector<G> >::self(A,a); }
//{unsecret}
template<class G,class V> inline typename bindminArray<const Matrix<G> >::self ele_min(const Matrix<G> &A, const V &a) { return typename bindminArray<const Matrix<G> >::self(A,a); }




DECLARE_FUNCTION2    (clip_low     ,clip_low_function )
DECLARE_FUNCTION2    (clip_high    ,clip_high_function)
DECLARE_BIND2ND_ARRAY(cliplowArray ,clip_low_function)
DECLARE_BIND2ND_ARRAY(cliphighArray,clip_high_function)
template<class G,class V> inline typename cliplowArray <const Vector<G> >::self clip_low (const Vector<G> &X, const V &v) { return typename cliplowArray <const Vector<G> >::self(X,v); }
template<class G,class V> inline typename cliplowArray <const Matrix<G> >::self clip_low (const Matrix<G> &X, const V &v) { return typename cliplowArray <const Matrix<G> >::self(X,v); }
template<class G,class V> inline typename cliphighArray<const Vector<G> >::self clip_high(const Vector<G> &X, const V &v) { return typename cliphighArray<const Vector<G> >::self(X,v); }
template<class G,class V> inline typename cliphighArray<const Matrix<G> >::self clip_high(const Matrix<G> &X, const V &v) { return typename cliphighArray<const Matrix<G> >::self(X,v); }


//DECLARE_FUNCTION(sum,sum_function)
template<class Arg,class Res=Arg> struct sum_function : public unary_value_function<Arg,Res>
{
  typedef unary_value_function<Arg,Res> base; UNARY_FUNCTION_BASE_TYPES;
  template<class V> struct rebind { typedef sum_function<V> other; };
  template<class V> inline typename sum_function<V>::reference       operator()(const V &x)       { return sum(x); }
  template<class V> inline typename sum_function<V>::const_reference operator()(const V &x) const { return sum(x); }
};
template<class T> struct sum_function<complex<T> > : public unary_value_function<complex<T>,complex<typename sum_function<T>::value_type> >
{
  typedef unary_value_function<complex<T>,complex<typename sum_function<T>::value_type> > base; UNARY_FUNCTION_BASE_TYPES;
  template<class V> struct rebind { typedef sum_function<V> other; };
  template<class V> inline typename sum_function<V>::reference       operator()(const V &x)       { return sum(x); }
  template<class V> inline typename sum_function<V>::const_reference operator()(const V &x) const { return sum(x); }
};
template<class T> struct array_sum_function : public unary_value_function<T,typename T::const_value_type>
{
  typedef unary_value_function<T,typename T::const_value_type> base; UNARY_FUNCTION_BASE_TYPES;
  template<class V> typename sum_function<V>::reference       operator()(const V &x)       { return sum(x); }
  template<class V> typename sum_function<V>::const_reference operator()(const V &x) const { return sum(x); }
  template<class V> typename TinyVector<2,typename sum_function<V>::reference      >::self operator()(const V &x0,const V &x1)       { return sum(x0,x1); }
  template<class V> typename TinyVector<2,typename sum_function<V>::const_reference>::self operator()(const V &x0,const V &x1) const { return sum(x0,x1); }
  template<class V> typename TinyVector<4,typename sum_function<V>::reference      >::self operator()(const V &x0,const V &x1,const V &x2,const V &x3)       { return sum(x0,x1,x2,x3); }
  template<class V> typename TinyVector<4,typename sum_function<V>::const_reference>::self operator()(const V &x0,const V &x1,const V &x2,const V &x3) const { return sum(x0,x1,x2,x3); }
};
template<       class G> struct sum_function<Vector<  G> > : public array_sum_function<Vector<  G> > {};
template<       class G> struct sum_function<Matrix<  G> > : public array_sum_function<Matrix<  G> > {};
template<int D, class G> struct sum_function<Array <D,G> > : public array_sum_function<Array <D,G> > {};
DECLARE_FUNCTION_ARRAY (elesumArray, sum_function)
//{unsecret}
//Summary: Sum (element-wise)
//Remarks: The elements should be arrays
//See: ^sum^
DECLARE_ARRAY_FUNCTION (ele_sum, elesumArray)

DECLARE_FUNCTION(mean,mean_function)
template<class T> struct array_mean_function : public unary_value_function<T,PROMOTE2(float,typename T::const_value_type)> { typedef unary_value_function<T,PROMOTE2(float,typename T::const_value_type)> base; UNARY_FUNCTION_BASE_TYPES; template<class V>inline  reference       operator()(const V &x)       { return mean(x); } template<class V> inline const_reference operator()(const V &x) const { return mean(x); } };
template<       class G> struct mean_function<Vector<  G> > : public array_mean_function<Vector<  G> > {};
template<       class G> struct mean_function<Matrix<  G> > : public array_mean_function<Matrix<  G> > {};
template<int D, class G> struct mean_function<Array <D,G> > : public array_mean_function<Array <D,G> > {};
DECLARE_FUNCTION_ARRAY (elemeanArray, mean_function)
//{unsecret}
//Summary: Mean (element-wise)
//Remarks: The elements should be arrays
//See: ^mean^
DECLARE_ARRAY_FUNCTION (ele_mean, elemeanArray)


DECLARE_OPERATOR(*,dereference_function);
DECLARE_FUNCTION_ARRAY(dereferenceArray,dereference_function);
//{unsecret}
//Summary: Dereference (element-wise)
DECLARE_ARRAY_OPERATOR          (*,dereferenceArray)
DECLARE_NON_CONST_ARRAY_OPERATOR(*,dereferenceArray)


//{unsecret}
//Summary: Logical not (element-wise)
template<class G> inline Vector<predicate_generator<const Vector<G>,logical_not<typename Vector<G>::value_type> > > ele_logical_not(const Vector<G> &X) { return Vector<predicate_generator<const Vector<G>,logical_not<typename Vector<G>::value_type> > >(X); }
//{unsecret}
template<class G> inline Vector<predicate_generator<const Vector<G>,logical_not<typename Matrix<G>::value_type> > > ele_logical_not(const Matrix<G> &X) { return Vector<predicate_generator<const Matrix<G>,logical_not<typename Vector<G>::value_type> > >(X); }

template<class G> inline Vector<predicate_generator<const Vector<G>,binder1st<equal_to    <typename Vector<G>::value_type> > > > ele_equal_to    (const typename Vector<G>::value_type &a, const Vector<G> &X) { return Vector<predicate_generator<const Vector<G>,binder1st<equal_to    <typename Vector<G>::value_type> > > >(X,a); }
template<class G> inline Vector<predicate_generator<const Vector<G>,binder1st<not_equal_to<typename Vector<G>::value_type> > > > ele_not_equal_to(const typename Vector<G>::value_type &a, const Vector<G> &X) { return Vector<predicate_generator<const Vector<G>,binder1st<not_equal_to<typename Vector<G>::value_type> > > >(X,a); }
template<class G> inline Matrix<predicate_generator<const Matrix<G>,binder1st<equal_to    <typename Matrix<G>::value_type> > > > ele_equal_to    (const typename Matrix<G>::value_type &a, const Matrix<G> &X) { return Matrix<predicate_generator<const Matrix<G>,binder1st<equal_to    <typename Matrix<G>::value_type> > > >(X,a); }
template<class G> inline Matrix<predicate_generator<const Matrix<G>,binder1st<not_equal_to<typename Matrix<G>::value_type> > > > ele_not_equal_to(const typename Matrix<G>::value_type &a, const Matrix<G> &X) { return Matrix<predicate_generator<const Matrix<G>,binder1st<not_equal_to<typename Matrix<G>::value_type> > > >(X,a); }

template<class G> inline Vector<predicate_generator<const Vector<G>,binder2nd<equal_to    <typename Vector<G>::value_type> > > > ele_equal_to    (const Vector<G> &X, const typename Vector<G>::value_type &a) { return Vector<predicate_generator<const Vector<G>,binder2nd<equal_to    <typename Vector<G>::value_type> > > >(X,a); }
template<class G> inline Vector<predicate_generator<const Vector<G>,binder2nd<not_equal_to<typename Vector<G>::value_type> > > > ele_not_equal_to(const Vector<G> &X, const typename Vector<G>::value_type &a) { return Vector<predicate_generator<const Vector<G>,binder2nd<not_equal_to<typename Vector<G>::value_type> > > >(X,a); }
template<class G> inline Matrix<predicate_generator<const Matrix<G>,binder2nd<equal_to    <typename Matrix<G>::value_type> > > > ele_equal_to    (const Matrix<G> &X, const typename Matrix<G>::value_type &a) { return Matrix<predicate_generator<const Matrix<G>,binder2nd<equal_to    <typename Matrix<G>::value_type> > > >(X,a); }
template<class G> inline Matrix<predicate_generator<const Matrix<G>,binder2nd<not_equal_to<typename Matrix<G>::value_type> > > > ele_not_equal_to(const Matrix<G> &X, const typename Matrix<G>::value_type &a) { return Matrix<predicate_generator<const Matrix<G>,binder2nd<not_equal_to<typename Matrix<G>::value_type> > > >(X,a); }

template<class G1,class G2> inline Vector<binary_predicate_generator<const Vector<G1>,const Vector<G2>,equal_to    <typename promotion2_traits<typename Vector<G1>::value_type,typename Vector<G2>::value_type>::value_type> > > ele_equal_to    (const Vector<G1> &X, const Vector<G2> &Y) { return Vector<binary_predicate_generator<const Vector<G1>,const Vector<G2>,equal_to    <typename promotion2_traits<typename Vector<G1>::value_type,typename Vector<G2>::value_type>::value_type> > >(X,Y); }
template<class G1,class G2> inline Vector<binary_predicate_generator<const Vector<G1>,const Vector<G2>,not_equal_to<typename promotion2_traits<typename Vector<G1>::value_type,typename Vector<G2>::value_type>::value_type> > > ele_not_equal_to(const Vector<G1> &X, const Vector<G2> &Y) { return Vector<binary_predicate_generator<const Vector<G1>,const Vector<G2>,not_equal_to<typename promotion2_traits<typename Vector<G1>::value_type,typename Vector<G2>::value_type>::value_type> > >(X,Y); }
template<class G1,class G2> inline Matrix<binary_predicate_generator<const Matrix<G1>,const Matrix<G2>,equal_to    <typename promotion2_traits<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type>::value_type> > > ele_equal_to    (const Matrix<G1> &X, const Matrix<G2> &Y) { return Matrix<binary_predicate_generator<const Matrix<G1>,const Matrix<G2>,equal_to    <typename promotion2_traits<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type>::value_type> > >(X,Y); }
template<class G1,class G2> inline Matrix<binary_predicate_generator<const Matrix<G1>,const Matrix<G2>,not_equal_to<typename promotion2_traits<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type>::value_type> > > ele_not_equal_to(const Matrix<G1> &X, const Matrix<G2> &Y) { return Matrix<binary_predicate_generator<const Matrix<G1>,const Matrix<G2>,not_equal_to<typename promotion2_traits<typename Matrix<G1>::value_type,typename Matrix<G2>::value_type>::value_type> > >(X,Y); }


template<class Arg1,class Arg2,class Res=Arg2>
class if_then_function : public binary_value_function<Arg1,Arg2,Res>
{
  public:
    typedef binary_value_function<Arg1,Arg2,Res> base;
    BINARY_FUNCTION_BASE_TYPES

  private:
    value_type val;

  public:
    inline if_then_function() : val() {}
    inline explicit if_then_function(const value_type &v) :val(v) {}

    template<class V1,class V2> inline const_reference operator()(const V1 &x1, const V2 &x2) const { if (x1) return x2; else return val; }
};

template<class A1, class A2, class V=typename if_then_function<typename A1::const_value_type,typename A2::const_value_type>::result_type>
struct if_then_array : public array_value_traits<A2,V>
{
  typedef array_value_traits<A2,V> base;
  ARRAY_BASE_TYPES
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef typename first_array_type ::const_value_type first_argument_type;
  typedef typename second_array_type::const_value_type second_argument_type;

  typedef if_then_function<first_argument_type,second_argument_type,value_type> function_type;
  typedef typename BinaryFunctionArray<first_array_type,second_array_type,function_type>::self self;
};

template<class G1,class G2> inline typename if_then_array<const Vector<G1>, const Vector<G2> >::self if_then(const Vector<G1> &X, const Vector<G2> &Y, typename if_then_function<const typename Vector<G1>::value_type, const typename Vector<G2>::value_type >::value_type &v) { return typename if_then_array<const Vector<G1>, const Vector<G2> >::self(X,Y,v); }
template<class G1,class G2> inline typename if_then_array<const Matrix<G1>, const Matrix<G2> >::self if_then(const Matrix<G1> &X, const Matrix<G2> &Y, typename if_then_function<const typename Matrix<G1>::value_type, const typename Matrix<G2>::value_type >::value_type &v) { return typename if_then_array<const Matrix<G1>, const Matrix<G2> >::self(X,Y,v); }


template<class Arg1,class Arg2,class Arg3,class Res=Arg3>
class if_then_else_function : public tertiary_value_function<Arg1,Arg2,Arg3,Res>
{
  public:
    typedef tertiary_value_function<Arg1,Arg2,Arg3,Res> base;
    TERTIARY_FUNCTION_BASE_TYPES

  public:
    inline if_then_else_function() {}

    template<class V1,class V2,class V3> inline const_reference operator()(const V1 &x, const V2 &y, const V3 &z) const { if (x) return y; else return z; }
};

template<class A1,class A2,class A3=A2,
   class V       =typename if_then_else_function<typename A1::const_value_type,typename A2::const_value_type,typename A3::const_value_type>::result_type,
   class Ref     =typename if_then_else_function<typename A1::const_value_type,typename A2::const_value_type,typename A3::const_value_type>::reference,
   class ConstRef=typename if_then_else_function<typename A1::const_value_type,typename A2::const_value_type,typename A3::const_value_type>::const_reference>
struct if_then_else_array : public three_arrays_reference_traits<A1,A2,A3,V,Ref,ConstRef>
{
  typedef three_arrays_reference_traits<A1,A2,A3,V,Ref,ConstRef> base;
  THREE_ARRAYS_BASE_TYPES

  typedef typename first_array_type ::const_value_type first_argument_type;
  typedef typename second_array_type::const_value_type second_argument_type;
  typedef typename third_array_type ::const_value_type third_argument_type;

  typedef if_then_else_function<first_argument_type,second_argument_type,third_argument_type,value_type> function_type;
  typedef tertiary_function_generator<first_array_type, second_array_type, third_array_type, function_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G1,class G2,class G3> inline typename if_then_else_array<const Vector<G1>, const Vector<G2>, const Vector<G3> >::self if_then_else(const Vector<G1> &X, const Vector<G2> &Y, const Vector<G3> &Z) { return typename if_then_else_array<const Vector<G1>,const Vector<G2>,const Vector<G3> >::self(X,Y,Z); }
template<class G1,class G2,class G3> inline typename if_then_else_array<const Matrix<G1>, const Matrix<G2>, const Vector<G3> >::self if_then_else(const Matrix<G1> &X, const Matrix<G2> &Y, const Matrix<G3> &Z) { return typename if_then_else_array<const Matrix<G1>,const Matrix<G2>,const Matrix<G3> >::self(X,Y,Z); }


template<class Arg>
struct first_component_function : public unary_reference_function<Arg,typename Arg::first_type>
{
  typedef unary_reference_function<Arg,typename Arg::first_type> base;
  UNARY_FUNCTION_BASE_TYPES
  template<class Arg1> inline reference       operator()(      Arg1 &x)       { return first(x); }
  template<class Arg1> inline const_reference operator()(const Arg1 &x) const { return first(x); }
};

template<class Arg>
struct second_component_function : public unary_reference_function<Arg,typename Arg::second_type>
{
  typedef unary_reference_function<Arg,typename Arg::second_type> base;
  UNARY_FUNCTION_BASE_TYPES
  template<class Arg1> inline reference       operator()(      Arg1 &x)       { return second(x); }
  template<class Arg1> inline const_reference operator()(const Arg1 &x) const { return second(x); }
};

template<class Arg>
struct third_component_function : public unary_reference_function<Arg,typename Arg::third_type>
{
  typedef unary_reference_function<Arg,typename Arg::third_type> base;
  UNARY_FUNCTION_BASE_TYPES
  template<class Arg1> inline reference       operator()(      Arg1 &x)       { return third(x); }
  template<class Arg1> inline const_reference operator()(const Arg1 &x) const { return third(x); }
};

template<class T>
struct firstComponentArray : public array_traits<T>
{
  typedef array_traits<T> base; ARRAY_BASE_TYPES
  typedef first_component_function<value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<class T>
struct secondComponentArray : public array_traits<T>
{
  typedef array_traits<T> base; ARRAY_BASE_TYPES
  typedef second_component_function<value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<class T>
struct thirdComponentArray : public array_traits<T>
{
  typedef array_traits<T> base; ARRAY_BASE_TYPES
  typedef third_component_function<value_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

//{unsecret}
//{noAutoLink}
//Summary: First part (element-wise)
//Remarks:
// This function supposes the existence of a function /first/ applicable on each element.
// A function /first/ is defined for color spaces of the library and for the STL /pair/ structur.
// The function /first/ has to be overloaded for any other classes,
// and should give a reference back, if a use as left-value is expected.
//Return: An array representing the first part of X
//See: ^second^, ^third^
//Example:
//  DenseVector<pair<int,int> >::self X(3);
//  X[0]=pair<int,float>(1,2); X[0]=pair<int,float>(3,4); X[2]=pair<int,float>(5,6);
//  cout << first(X) << endl; // [1 3 5]
//  first(X) = 0;
//  cout << first(X) << endl; // [0 0 0]
//
//  RGBImage X;
//  PPMFile fin("x.ppm");
//  fin >> X; // load a RGB image from a PPM file
//  ucharImage Y = first(value_cast<YUV>(X)); // Greyscale image of X
template<class G> inline typename firstComponentArray <      Vector<G> >::self first (      Vector<G> &X) { return typename firstComponentArray <      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename firstComponentArray <const Vector<G> >::self first (const Vector<G> &X) { return typename firstComponentArray <const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename firstComponentArray <      Matrix<G> >::self first (      Matrix<G> &X) { return typename firstComponentArray <      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename firstComponentArray <const Matrix<G> >::self first (const Matrix<G> &X) { return typename firstComponentArray <const Matrix<G> >::self(X); }

//{unsecret}
//{noAutoLink}
//Summary: Second part (element-wise)
//Remarks:
// This function supposes the existence of a function /second/ applicable on each element.
// A function /second/ is defined for color spaces of the library and for the STL /pair/ structur.
// The function /second/ has to be overloaded for any other classes,
// and should give a reference back, if a use as left-value is expected.
//Return: An array representing the second part of X
//See: ^first^, ^third^
template<class G> inline typename secondComponentArray<      Vector<G> >::self second(      Vector<G> &X) { return typename secondComponentArray<      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename secondComponentArray<const Vector<G> >::self second(const Vector<G> &X) { return typename secondComponentArray<const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename secondComponentArray<      Matrix<G> >::self second(      Matrix<G> &X) { return typename secondComponentArray<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename secondComponentArray<const Matrix<G> >::self second(const Matrix<G> &X) { return typename secondComponentArray<const Matrix<G> >::self(X); }

//{unsecret}
//{noAutoLink}
//Summary: Third part (element-wise)
//Remarks:
// This function supposes the existence of a function /third/ applicable on each element.
// A function /third/ is defined for color spaces of the library.
// The function /third/ has to be overloaded for any other classes,
// and should give a reference back, if a use as left-value is expected.
//Return: An array representing the third part of X
//See: ^first^, ^second^
template<class G> inline typename thirdComponentArray <      Vector<G> >::self third (      Vector<G> &X) { return typename thirdComponentArray <      Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename thirdComponentArray <const Vector<G> >::self third (const Vector<G> &X) { return typename thirdComponentArray <const Vector<G> >::self(X); }
//{unsecret}
template<class G> inline typename thirdComponentArray <      Matrix<G> >::self third (      Matrix<G> &X) { return typename thirdComponentArray <      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename thirdComponentArray <const Matrix<G> >::self third (const Matrix<G> &X) { return typename thirdComponentArray <const Matrix<G> >::self(X); }

//}

#endif

