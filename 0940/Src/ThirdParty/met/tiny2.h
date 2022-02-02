//-*- C++ -*-
//
// Disambiguated Glommarble Expression Templates
// for Matrix Calculations
//
// ref. Computers in Physics, 11, 263 (1997)

#ifndef XPR2_H_
#define XPR2_H_

#include "xpr1.h"
#include "tiny1.h"
#include "xpr2.h"

// Template Metaprogram for Inner (Column) Loop Unrolling
template <int Nr, int Nc, class I, class A>
class Metaloop1of2 {
  static inline void assign(I* me, const A& rhs) {
    Metaloop1of2<Nr,Nc-1,I,A>::assign(me,rhs);
    me->operator()(Nr-1,Nc-1) = rhs(Nr-1,Nc-1);
  }
};

#define XXX(scalar) \
template <int Nr, int Nc, class I> \
class Metaloop1of2<Nr,Nc,I,scalar> {\
  static inline void assign(I* me, scalar rhs) {\
    Metaloop1of2<Nr,Nc-1,I,scalar>::assign(me,rhs);\
    me->operator()(Nr-1,Nc-1) = scalar(0);\
  }\
};
XXX(int)
XXX(double)
XXX(Complex )
#undef XXX

#define XXX(scalar) \
template <int N, class I> \
class Metaloop1of2<N,N,I,scalar> {\
  static inline void assign(I* me, scalar rhs) {\
    Metaloop1of2<N,N-1,I,scalar>::assign(me,rhs);\
    me->operator()(N-1,N-1) = rhs;\
  }\
};
XXX(int)
XXX(double)
XXX(Complex )
#undef XXX

template <int Nr, class I, class A>
class Metaloop1of2<Nr,1,I,A> {
  static inline void assign(I* me, const A& rhs) {
    me->operator()(Nr-1,0) = rhs(Nr-1,0);
  }
};

#define XXX(scalar) \
template <int Nr,class I> \
class Metaloop1of2<Nr,1,I,scalar> {\
  static inline void assign(I* me, scalar rhs) {\
    me->operator()(Nr-1,0) = scalar(0);\
  }\
};
XXX(int)
XXX(double)
XXX(Complex )
#undef XXX

#define XXX(scalar) \
template <class I> \
class Metaloop1of2<1,1,I,scalar> {\
  static inline void assign(I* me, scalar rhs) {\
    me->operator()(0,0) = rhs;\
  }\
};
XXX(int)
XXX(double)
XXX(Complex )
#undef XXX

// Template Metaprogram for Outer (Row) Loop Unrolling
template <int Nr, int Nc, class I, class A>
class Metaloop2 {
  static inline void assign(I* me, const A& rhs) {
    Metaloop2<Nr-1,Nc,I,A>::assign(me,rhs);
    Metaloop1of2<Nr,Nc,I,A>::assign(me,rhs);
  }
};

template <int Nc, class I, class A>
class Metaloop2<1,Nc,I,A> {
  static inline void assign(I* me, const A& rhs) {
    Metaloop1of2<1,Nc,I,A>::assign(me,rhs);
  }
};

template <int Nr, int Nc, class P, class I>
class TDim2 {
public:
  explicit TDim2() {}
  int size() const { return Nr*Nc; }
  int rows() const { return Nr; }
  int cols() const { return Nc; }
  //P operator() (int n) const {
  //  return static_cast<const I*>(this)->operator()(n);
  //}
  P operator() (int i, int j) const {
    return static_cast<const I*>(this)->operator()(i,j);
  }
  template <class X> I& assignFrom(const Xpr2<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    Metaloop2<Nr,Nc,I,Xpr2<P,X> >::assign(me,rhs);
    return *me;
  }
  template <class M> I& assignFrom(const Dim2<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    Metaloop2<Nr,Nc,I,Dim2<P,M> >::assign(me,rhs);
    return *me;
  }
  I& assignFrom(P rhs) {
    I *me = static_cast<I*>(this);
    Metaloop2<Nr,Nc,I,P>::assign(me,rhs);
    return *me;
  }
};

template <class P,class I> template <int Nr,int Nc,class M>
I& Dim2<P,I>::assignFrom(const TDim2<Nr,Nc,P,M>& rhs) {
  I *me = static_cast<I*>(this);
  for (int i=0; i<me->rows(); i++) {
    for (int j=0; j<me->cols(); j++) me->operator()(i,j) = rhs(i,j);
  }
  return *me;
}


// Template Metaprogram for (Coloumn-Row) Reduction Loop Unrolling
template <int Nc,class T,class D2, class D1>
class MetaReduct3to1 {
  static inline T sum(const D2& a, const D1& b, int i) {
    return MetaReduct3to1<Nc-1,T,D2,D1>::sum(a,b,i) + a(i,Nc-1)*b(Nc-1);
  }
};

template <class T,class D2, class D1>
class MetaReduct3to1<1,T,D2,D1> {
  static inline T sum(const D2& a, const D1& b, int i) {
    return a(i,0)*b(0);
  }
};

// Loop Unrolled Matrix Vector Multiplication
template <int N, class P, class A, class B>
class TXpr1Reduct {
  A a;
  B b;
public:
  TXpr1Reduct(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i) const {
    return MetaReduct3to1<N,P,A,B>::sum(a,b,i);
  }
};



#endif // XPR2_H_
