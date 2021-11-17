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

#ifndef MEMORY_H
#define MEMORY_H

#include <memory>

//namespace genial;
//{

using namespace std;

// Thread-Local Storage
#if defined(__ICL) || defined(_MSC_VER)
  //#define __thread
  #define __thread __declspec( thread )
#elif defined(__CYGWIN__) || defined(__MINGW32__) || defined(__APPLE__) // do not support TLS
  #define __thread
#elif defined(__GNUC__)
  //#define __thread
#else
  #define __thread
#endif


#if defined(__ICL)
  #define __aligned __declspec(align(16))
#elif defined(_MSC_VER)
  #define __aligned __declspec(align(16))
#elif defined (__GNUC__)
  #define __aligned __attribute__((__aligned__(16)))
#else
  #define __aligned
#endif

#if defined(__GNUC__) && defined(__i386__) && !(defined(__MACOSX__) || defined(__APPLE__))
  #define __aligned_stack __builtin_alloca(16); __asm__ __volatile__ ("andl $-16, %esp");
#elif defined(__ICC)
  #define __aligned_stack  alloca(16);
#else
  #define __aligned_stack 
#endif

#if (defined(__ICL) || defined(_MSC_VER) || defined(__ICC))
  #include <fvec.h>
  inline void *aligned_malloc (size_t size, size_t align=16) { return _mm_malloc(size+align,align); }
  inline void  aligned_free   (void *p)                      { return _mm_free(p); }
#elif defined (__CYGWIN__)
  #include <xmmintrin.h>
  inline void *aligned_malloc (size_t size, size_t align=16) { return _mm_malloc(size+align,align);  }
  inline void  aligned_free   (void *p)                      { return _mm_free(p); }
#elif defined(__MINGW32__)
  #include <malloc.h>
  inline void *aligned_malloc (size_t size, size_t align=16) { return __mingw_aligned_malloc(size+align,align);  }
  inline void  aligned_free   (void *p)                      { return __mingw_aligned_free(p);             }
#elif defined(__FreeBSD__) // Many thanks to Myles Rossiya
  #include <stdlib.h>
  inline void* aligned_malloc (size_t size, size_t align=16) { return malloc(size+align); }
  inline void  aligned_free   (void *p)                      { return free(p); }
#elif (defined(__MACOSX__) || defined(__APPLE__)) // Many thanks to Epsilon68
  #include <stdlib.h>
  inline void* aligned_malloc (size_t size, size_t align=16) { return malloc(size+align); }
  inline void  aligned_free   (void *p)                      { return free(p); }
#else
  #include <malloc.h>
  inline void* aligned_malloc (size_t size, size_t align=16) { return memalign(align,size+align); }
  inline void  aligned_free   (void *p)                      { return free(p); }
#endif


template<class T, int N=16> class alignment_allocator
{
  public:
    typedef T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    typedef T* pointer;
    typedef const T* const_pointer;

    typedef T& reference;
    typedef const T& const_reference;

  public:
    inline alignment_allocator() throw() {}
    template <class T2> inline alignment_allocator(const alignment_allocator<T2,N>&) throw() {}

    inline ~alignment_allocator() throw() {}

    inline pointer       address(reference       r)       { return &r; }
    inline const_pointer address(const_reference r) const { return &r; }

    inline pointer allocate(size_type n) { return (pointer)aligned_malloc(n*sizeof(value_type),N); }
    inline void deallocate(pointer p, size_type) { aligned_free(p); }

    inline void construct (pointer p,const value_type& val)  { new (p) value_type(val); }
    inline void destroy   (pointer p                      )  { p->~value_type();         }

    inline size_type max_size() const throw() { return size_type(-1)/sizeof(value_type); }

    template<class T2> struct rebind { typedef alignment_allocator<T2,N> other; };
};



template<class FwdIt,class V> inline void uninitialized_fill_aux(FwdIt begin, FwdIt end, const V &v) { }
template<class FwdIt> inline void uninitialized_fill(FwdIt begin, FwdIt end)
{
  uninitialized_fill_aux(begin,end,typename iterator_traits<FwdIt>::value_type());
}


template<class V> inline V *copy(const V *p) { return new V(*p); }
//template<class V> inline V *copy(const V *p) { return p?new V(*p):NULL; }
template<class V> inline void destroy(V *p) { delete p; }
template<class V> inline V *isolate(V *p) { return p; }


template<class V>
struct value_traits
{
  typedef V                 value_type;
  typedef value_type       &reference;
  typedef const value_type &const_reference;
  typedef value_type       *pointer;
  typedef const value_type *const_pointer;
};

template<class V>
struct value_traits<const V>
{
  typedef const V           value_type;
  typedef const value_type &reference;
  typedef const value_type &const_reference;
  typedef const value_type *pointer;
  typedef const value_type *const_pointer;
};


template<class V>
struct counter_value : public value_traits<V>
{
  public:
    typedef value_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;

    typedef counter_value self;

    typedef int counter_type;

  private:
    value_type val;
    mutable counter_type count;

  public:
    inline counter_value() : val(), count(1) {}
    inline counter_value(const value_type &v) : val(v), count(1) {}

    inline counter_value(const self &r) : val(r.val), count(1) { }

//    counter_value &operator=(const value_type &v) { val=v; return *this; }
//    counter_value &operator=(const self &r) { val=r; return *this; }

    ~counter_value() {}

    counter_type counter() const { return count; }
    counter_type inc() const { return ++count; }
    counter_type dec() const { return --count; }

    reference       value()       { return val; }
    const_reference value() const { return val; }

//    operator reference      ()       { return value(); }
//    operator const_reference() const { return value(); }
};

template<class E, class T, class V> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const counter_value<V> &x) { os << x.value(); return os; }

template<class V> inline counter_value<V> *copy(const counter_value<V> *p) { p->inc(); return const_cast<counter_value<V> *>(p); }
template<class V> inline void destroy(counter_value<V> *p) { if (p && !p->dec()) delete p; }
template<class V> inline counter_value<V> *isolate(counter_value<V> *p) { if (p->counter()==1) return p; p->dec(); return new counter_value<V>(*p); }



template<class V>
struct auto_pointer : public value_traits<V>
{
  public:
    typedef auto_pointer self;
    typedef value_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;

    typedef random_access_iterator_tag iterator_category;
    typedef ptrdiff_t difference_type;

  private:
    pointer ptr;

  public:
    inline auto_pointer() : ptr(NULL) {}
    inline auto_pointer(pointer p) : ptr(p) {}
    template<class V2> inline explicit auto_pointer(V2 *p) : ptr(p) {}

    inline auto_pointer(const self &p) : ptr(p.get()?copy(p.get()):NULL) {}
    template<class T2> inline explicit auto_pointer(const auto_pointer<T2> &p) : ptr(p.get()?copy(p.get()):NULL) {}

    inline self &operator=(pointer p) { destroy(); ptr=p; return *this; }
    inline self &operator=(const self &p) { destroy(); ptr=p.get()?copy(p.get()):NULL; return *this; }
    template<class T2> inline self &operator=(T2 *p) { destroy(); ptr=p; return *this; }
    template<class T2> inline self &operator=(const auto_pointer<T2> &p) { destroy(); ptr=p.get()?copy(p.get()):NULL; return *this; }

    inline ~auto_pointer() { destroy(); }

    inline pointer       get()       { return  ptr; }
    inline const_pointer get() const { return  ptr; }
    inline operator pointer      ()       { return ptr; }
    inline operator const_pointer() const { return ptr; }

    inline reference       operator* ()       { return *ptr; }
    inline const_reference operator* () const { return *ptr; }
    inline pointer         operator->()       { return  ptr; }
    inline const_pointer   operator->() const { return  ptr; }

  private:
    inline void destroy() { ::destroy(ptr); }
};

template<class V> auto_pointer<      V> autopointer(      V *p) { return auto_pointer<      V>(p); }
template<class V> auto_pointer<const V> autopointer(const V *p) { return auto_pointer<const V>(p); }


template<class V>
struct shared_pointer : public value_traits<V>
{
  public:
    typedef shared_pointer self;
    typedef value_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;

    typedef random_access_iterator_tag iterator_category;
    typedef ptrdiff_t difference_type;

  private:
    counter_value<value_type> *ptr;

  public:
    inline shared_pointer() : ptr(NULL) {}
    inline shared_pointer(counter_value<value_type> *p) : ptr(p) {}
    template<class V2> inline explicit shared_pointer(counter_value<V2> *p) : ptr(p) {}

    inline shared_pointer(const self &p) : ptr(p.get()?copy(p.get()):NULL) {}
    template<class T2> inline explicit shared_pointer(const shared_pointer<T2> &p) : ptr(p.get()?copy(p.get()):NULL) {}

    inline self &operator=(pointer p) { destroy(); ptr=p; return *this; }
    inline self &operator=(const self &p) { destroy(); ptr=p.get()?copy(p.get()):NULL; return *this; }
    template<class T2> inline self &operator=(T2 *p) { destroy(); ptr=p; return *this; }
    template<class T2> inline self &operator=(const shared_pointer<T2> &p) { destroy(); ptr=p.get()?copy(p.get()):NULL; return *this; }

    inline ~shared_pointer() { destroy(); }

    inline counter_value<value_type>       *get()       { return  ptr; }
    inline const counter_value<value_type> *get() const { return  ptr; }

    inline operator pointer      ()       { return &ptr->value(); }
    inline operator const_pointer() const { return &ptr->value(); }

    inline reference       operator* ()       { return ptr->value(); }
    inline const_reference operator* () const { return ptr->value(); }

    inline pointer         operator->()       { return &ptr->value(); }
    inline const_pointer   operator->() const { return &ptr->value(); }

  private:
    inline void destroy() { ::destroy(ptr); }
};

template<class V> shared_pointer<      V> sharedpointer(      V *p) { return shared_pointer<      V>(p); }
template<class V> shared_pointer<const V> sharedpointer(const V *p) { return shared_pointer<const V>(p); }



template<class V>
class auto_value : public value_traits<V>
{
  public:
    typedef auto_value self;
    typedef value_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;

  private:
    pointer ptr;

  public:
    inline auto_value() : ptr(new value_type()) {}
    inline explicit auto_value(const value_type &v) : ptr(copy(&v)) {}
    template<class V2> inline explicit auto_value(const V2 &v) : ptr(copy(&v)) {}

    inline auto_value(const self &r) : ptr(copy(r.get())) { }
    template<class V2> inline explicit auto_value(const auto_value<V2> &r) : ptr(copy(r.get())) { }

    inline auto_value &operator=(const value_type &v) { destroy(); ptr=copy(&v); return *this; }
    inline auto_value &operator=(const self &r) { destroy(); ptr=copy(r.get()); return *this; }
    template<class V2> inline auto_value &operator=(V2 &v) { destroy(); ptr=copy(&v); return *this; }
    template<class V2> inline auto_value &operator=(const auto_value<V2> &r) { destroy(); ptr=copy(r.get()); return *this; }

    inline ~auto_value() { destroy(); }

    inline reference       value()       { return *ptr; }
    inline const_reference value() const { return *ptr; }

    inline pointer       get()       { return  ptr; }
    inline const_pointer get() const { return  ptr; }
    inline operator reference      ()       { return value(); }
    inline operator const_reference() const { return value(); }

  private:
    inline void destroy() { ::destroy(ptr); }
};

template<class E, class T, class V> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const auto_value<V> &x) { os << x.value(); return os; }

template<class Arg> struct name_function;
template<class Arg> struct name_function<      auto_value<Arg> > : public name_function<      Arg> {};
template<class Arg> struct name_function<const auto_value<Arg> > : public name_function<const Arg> {};


template<class V>
struct value_traits<auto_value<V> >
{
  typedef typename value_traits<V>::value_type      value_type;
  typedef typename value_traits<V>::reference       reference;
  typedef typename value_traits<V>::const_reference const_reference;
  typedef typename value_traits<V>::pointer         pointer;
  typedef typename value_traits<V>::const_pointer   const_pointer;
};

template<class V>
struct value_traits<const auto_value<V> >
{
  typedef const typename value_traits<V>::value_type      value_type;
  typedef const typename value_traits<V>::reference       reference;
  typedef const typename value_traits<V>::const_reference const_reference;
  typedef const typename value_traits<V>::pointer         pointer;
  typedef const typename value_traits<V>::const_pointer   const_pointer;
};

template<class V>
class shared_value : public value_traits<V>
{
  public:
    typedef shared_value self;
    typedef value_traits<V> base;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::pointer         pointer;
    typedef typename base::const_pointer   const_pointer;

  private:
    counter_value<value_type> *ptr;

  public:
    inline shared_value() : ptr(new counter_value<value_type>()) {}
    inline explicit shared_value(const value_type &v) : ptr(new counter_value<value_type>(v)) {}

    inline shared_value(const self &r) : ptr(copy(r.ptr)) {}
    template<class V2> inline explicit shared_value(const shared_value<V2> &r) : ptr(copy(r.ptr)) {}

    inline ~shared_value() { destroy(); }

    inline self &operator=(const value_type &v) { destroy(); ptr=new counter_value<value_type>(v); return *this; }
    inline self &operator=(const self &r) { destroy(); ptr=copy(r.get()); return *this; }
    template<class V2> inline self &operator=(const V2 &v) { destroy(); ptr=copy(&v); return *this; }
    template<class V2> inline self &operator=(const shared_value<V2> &r) { destroy(); ptr=copy(r.get()); return *this; }

    inline counter_value<value_type>       *get()       { return ptr; }
    inline const counter_value<value_type> *get() const { return ptr; }

    inline reference       value()       { isolate(); return ptr->value(); }
    inline const_reference value() const { return ptr->value(); }

    inline operator reference      ()       { return value(); }
    inline operator const_reference() const { return value(); }

  private:
    inline void destroy() { ::destroy(ptr); }
    inline void isolate() { ptr=::isolate(ptr); }
};

template<class E, class T, class V> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const shared_value<V> &x) { os << x.value(); return os; }

template<class Arg> struct name_function<      shared_value<Arg> > : public name_function<      Arg> {};
template<class Arg> struct name_function<const shared_value<Arg> > : public name_function<const Arg> {};

template<class V>
struct value_traits<shared_value<V> >
{
  typedef typename value_traits<V>::value_type      value_type;
  typedef typename value_traits<V>::reference       reference;
  typedef typename value_traits<V>::const_reference const_reference;
  typedef typename value_traits<V>::pointer         pointer;
  typedef typename value_traits<V>::const_pointer   const_pointer;
};

template<class V>
struct value_traits<const shared_value<V> >
{
  typedef const typename value_traits<V>::value_type      value_type;
  typedef const typename value_traits<V>::reference       reference;
  typedef const typename value_traits<V>::const_reference const_reference;
  typedef const typename value_traits<V>::pointer         pointer;
  typedef const typename value_traits<V>::const_pointer   const_pointer;
};


//}

#endif
