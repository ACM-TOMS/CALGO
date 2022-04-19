//-*- C++ -*-
//
// Tiny 1 Dimemsional Expression Templates of the Constant Size
// (Glommable and Disambiguated)
//
// ref. Computers in Physics, 11, 263 (1997)
//      Computers in Physics, 10, 552 (1996)

#ifndef XPR1_H_
#define XPR1_H_

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>

#include "cmplxtype.h"

#include "xpr1.h"

// Template Metaprogram for Loop Unrolling
template <int N, class I, class A>
class Metaloop1 {
  static inline void assign(I* me, const A& rhs) {
    Metaloop1<N-1,I,A>::assign(me,rhs);
    me->operator()(N-1) = rhs(N-1);
  }
};

#define XXX(scalar) \
template <int N, class I> \
class Metaloop1<N,I,scalar> {\
  static inline void assign(I* me, scalar rhs) {\
    Metaloop1<N-1,I,scalar>::assign(me,rhs);\
    me->operator()(N-1) = rhs;\
  }\
};
XXX(int)
XXX(double)
XXX(Complex )
#undef XXX

template <class I, class A>
class Metaloop1<1,I,A> {
  static inline void assign(I* me, const A& rhs) {
    me->operator()(0) = rhs(0);
  }
};

#define XXX(scalar) \
template <class I> \
class Metaloop1<1,I,scalar> {\
  static inline void assign(I* me, scalar rhs) {\
    me->operator()(0) = rhs;\
  }\
};
XXX(int)
XXX(double)
XXX(Complex )
#undef XXX

// Tiny 1 Dimensional Array Base Class
// for  Glommable Expression Templates and Template Metaprograms
template <int N, class P, class I>
class TDim1 {
public:
  explicit TDim1() {}
  int size() const { return N; }
  P operator() (int n) const {
    return static_cast<const I*>(this)->operator()(n);
  }
  template <class E> I& assignFrom(const Xpr1<P,E>& x) {
    I *me = static_cast<I*>(this);
    Metaloop1<N,I,Xpr1<P,E> >::assign(me,x);
    return *me;
  }
  template <class V> I& assignFrom(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
    Metaloop1<N,I,Dim1<P,V> >::assignFrom(me,x);
    return *me;
  }
  template <class V> I& assignFrom(const TDim1<N,P,V>& x) {
    I *me = static_cast<I*>(this);
    Metaloop1<N,I,Dim1<P,V> >::assignFrom(me,x);
    return *me;
  }
  I& assignFrom(P x) {
    I *me = static_cast<I*>(this);
    Metaloop1<N,I,P>::assignFrom(me,x);
    return *me;
  }
};

template <class P,class I> template <int N,class V>
I& Dim1<P,I>::assignFrom(const TDim1<N,P,V>& x) {
  I *me = static_cast<I*>(this);
  for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
  return *me;
}

#endif // XPR1_H_
