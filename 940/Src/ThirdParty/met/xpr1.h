//-*- C++ -*-
//
// 1 Dimensional Expression Templates
// (Glommable and Disambiguated)
//
// ref. Computers in Physics, 11, 263 (1997)
//      Computers in Physics, 10, 552 (1996)
//
//     Copyright (C) 2000 Masakatsu Ito
//     Institute for Molecular Science
//     Myodaiji  Okazaki  444-0867 Japan
//     E-mail mito@ims.ac.jp

//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 

#ifndef XPR1_H_
#define XPR1_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "error.h"
#include "cmplxtype.h"

class Range {
  const int sz, sta, ste;
public:
  Range(int sta_, int end, int ste_=1) : // [ sta_, end )
  sz( (end - sta_)/ste_), sta(sta_), ste(ste_) {}
  Range(const Range& sl) : sz(sl.sz), sta(sl.sta), ste(sl.ste) {}
  //Range& operator=(const Range& x) { sz = x.sz; sta = x.sta; ste=x.ste; }
  Range operator-() const { return Range(sta, sta-sz*ste, -ste); }
  Range operator+(int x) const { return Range(sta+x, sta+sz*ste+x, ste); }
  Range operator-(int x) const { return Range(sta-x, sta+sz*ste-x, ste); }
  Range operator*(int x) const { 
    return Range(sta*x, (sta+sz*ste)*x, ste*x); 
  }
//   Range operator/(int x) const { 
//     if (sta) {
//       return Range(sta/x, (sta+sz*ste)/x, ste/x); 
//     } else {
//       return Range(sta, (sta+sz*ste)/x, ste/x); 
//     }
//   }
  int operator()(int i) const { return sta + ste*i ; }
  int size() const { return sz; }
  int start() const { return sta; }
  int step() const { return ste; }
};

template <class P>
class Identity {
public:
  static inline P apply(const P& a) { return a; }
};

template <class P>
class UnaryMinus {
public:
  static inline P apply(const P& a) { return -a; }
};

template <class P>
class Exp {
public:
  static inline P apply(const P& a) { return P( exp(a) ); }
};

template <class P, class A, class F>
class Xpr1Func {
  A a;
public:
  Xpr1Func(const A& a_) : a(a_) {}
  P operator()(int n) const { return F::apply( a(n) ); }
// #ifdef SIZECHECK
//   int size() const { return a.size(); }
// #endif
};

template <class P, class F>
class Xpr1Func<P,Range,F> {
  //const Range sl;
  const Range& sl;
public:
  Xpr1Func(const Range& sl_) : sl(sl_) {}
  P operator()(int n) const { return F::apply( P(sl(n)) ); }
// #ifdef SIZECHECK
//   int size() const { return sl.size(); }
// #endif
};


#define XXX(ap,op) \
template <class P> \
class ap {\
public:\
  static inline P apply(const P& a, const P& b) { return a op b; }\
};
XXX(OpAdd,+)
XXX(OpSub,-)
XXX(OpMul,*)
XXX(OpDiv,/)
#undef XXX

template <class P, class A, class B, class Op>
class Xpr1BinOp {
  A a;
  B b;
public:
  Xpr1BinOp(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator() (int n) const { return Op::apply( a(n), b(n) ); }
// #ifdef SIZECHECK
//   int size() const { 
//     if ( a.size() != b.size() ) error("The size of Xpr1BinOp operands should be the same!");
//     return a.size(); 
//   }
// #endif
};

// Specialziation for the all operations with Scalar
template <class P, class A, class Op>
class Xpr1OpScalar {
  A a;
  P b;
public:
  Xpr1OpScalar(const A& a_, const P& b_) : a(a_), b(b_) {}
  P operator() (int n) const { return Op::apply( a(n), b ); }
// #ifdef SIZECHECK
//   int size() const { return a.size(); }
// #endif
};

// Specialziation for the +,- and * operations with Scalar
template <class P, class B, class Op>
class Xpr1ScalarOp {
  P a;
  B b;
public:
  Xpr1ScalarOp(const P& a_, const B& b_) : a(a_), b(b_) {}
  P operator() (int n) const { return Op::apply( a, b(n) ); }
// #ifdef SIZECHECK
//   int size() const { return b.size(); }
// #endif
};

template <class P, class V>
class ConstRef1 {
  const V& v;
public:
  ConstRef1(const V& v_) : v(v_) {}
  P operator()(int n) const { return v(n); }
// #ifdef SIZECHECK
//   int size() const { return v.size(); }
// #endif
};


template <class P, class E>
class Xpr1 {
private:
  E e;
public:
  Xpr1(const E& e_) : e(e_) {}
  P operator() (int n) const { return e(n); }
// #ifdef SIZECHECK
//   int size() const { return e.size(); }
// #endif
};

//template <int N, class P, class I> class TDim1;

// 1 Dimensional Array Base Class
// for Glommable Expression Templates
template <class P, class I>
class Dim1 {
private:
  void error(const char *msg) const {
      std::cerr << "Dim1 error: " << msg << std::endl;
    exit(1);
  }

public:
  explicit Dim1() {}
  int size() const {
    return static_cast<const I*>(this)->size();
  }
  P operator() (int n) const {
    return static_cast<const I*>(this)->operator()(n);
  }
  template <class E> I& assignFrom(const Xpr1<P,E>& x) {
    I *me = static_cast<I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for = operator.");
// #endif
    for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  template <class V> I& assignFrom(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for = operator.");
// #endif
    for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  //  template <int N,class V> I& assignFrom(const TDim1<N,P,V>& x);
  I& assignFrom(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) = x;
    return *me;
  }
  template <class E> Dim1<P,I>& operator+=(const Xpr1<P,E>& x) {
    I* me = static_cast<I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for += operator.");
// #endif
    for (int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    return *me;
  }
  template <class V> Dim1<P,I>& operator+=(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for += operator.");
// #endif
    for (int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    return *me;
  }
  Dim1<P,I>& operator+=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) += x;
    return *me;
  }
  template <class E> Dim1<P,I>& operator-=(const Xpr1<P,E>& x) {
    I* me = static_cast<I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for -= operator.");
// #endif
    for (int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    return *me;
  }
  template <class V> Dim1<P,I>& operator-=(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for -= operator.");
// #endif
    for (int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    return *me;
  }
  Dim1<P,I>& operator-=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) -= x;
    return *me;
  }
  Dim1<P,I>& operator*=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) *= x;
    return *me;
  }
  Dim1<P,I>& operator/=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) /= x;
    return *me;
  }
  
  template <class E> double in(const Xpr1<P,E>& x) const {
    double sum = 0.0;
    const I* me = static_cast<const I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for innner product.");
// #endif
    for (int i=0; i < me->size(); i++) sum += me->operator()(i) * x(i);
    return sum;
  }
  template <class V> double in(const Dim1<P,V>& x) const {
    double sum = 0.0;
    const I* me = static_cast<const I*>(this);
// #ifdef SIZECHECK
//     if ( me->size() != x.size() ) error("The size should be the same for innner product.");
// #endif
    for (int i=0; i < me->size(); i++) sum += me->operator()(i) * x(i);
    return sum;
  }
};

template <class T,class A>
std::ostream& operator<<(std::ostream& s, const Dim1<T,A>& a) {
  for (int i=0; i< a.size(); i++)
      s << std::setw(6) << std::setprecision(3) << a(i) << std::endl;
  return s;
}



// Functions of Dim1
#define XXX(f,ap) \
template <class P,class A> \
Xpr1<P, Xpr1Func<P, ConstRef1<P,Dim1<P,A> >, ap<P> > > \
f(const Dim1<P,A>& a) \
{\
   typedef Xpr1Func<P, ConstRef1<P, Dim1<P,A> >, ap<P> > ExprT;\
   return Xpr1<P,ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a))); \
}
XXX(ident, Identity)
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX

// Functions of Range
#define XXX(P,f,ap) \
Xpr1<P, Xpr1Func<P, Range, ap<P> > > \
f(const Range& a) \
{\
   typedef Xpr1Func<P, Range, ap<P> > ExprT;\
   return Xpr1<P,ExprT>(ExprT(a));\
}
XXX(double,dble, Identity)
XXX(double, dexp, Exp)
XXX(Complex ,cmplx, Identity)
XXX(Complex ,cexp, Exp)
#undef XXX


// Functions of Xpr1
#define XXX(f,ap) \
template <class P, class A> \
Xpr1<P, Xpr1Func<P, Xpr1<P,A>, ap<P> > > \
f(const Xpr1<P,A>& a) \
{\
   typedef Xpr1Func<P, Xpr1<P,A> , ap<P> > ExprT;\
   return Xpr1<P,ExprT>(ExprT(a)); \
}
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX

//Binary operations between Two Dim1s
#define XXX(op,ap) \
template <class P,class A,class B> \
Xpr1<P, Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, ConstRef1<P,Dim1<P,B> >, ap<P> > >\
op (const Dim1<P,A>& a, const Dim1<P,B>& b) {\
  typedef \
    Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, ConstRef1<P, Dim1<P,B> >, ap<P> > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a),\
			      ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Dim1 and Scalar
#define XXX(op,ap) \
template <class P,class A>\
Xpr1<P, Xpr1OpScalar<P, ConstRef1<P, Dim1<P,A> >, ap<P> > > \
op (const Dim1<P,A>& a, const P& b) {\
  typedef Xpr1OpScalar<P, ConstRef1<P, Dim1<P,A> >, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a),b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Scalar and Dim1
#define XXX(op,ap) \
template <class P,class B> \
Xpr1<P, Xpr1ScalarOp<P, ConstRef1<P, Dim1<P,B> >, ap<P> > > \
op (const P& a, const Dim1<P,B>& b) {\
  typedef Xpr1ScalarOp<P, ConstRef1<P, Dim1<P,B> >, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Xpr1 and Scalar
#define XXX(op,ap) \
template <class P,class A>\
Xpr1<P, Xpr1OpScalar<P, Xpr1<P,A>, ap<P> > > \
op (const Xpr1<P,A>& a, const P& b) {\
  typedef Xpr1OpScalar<P, Xpr1<P,A>, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Scalar and Xpr1
#define XXX(op,ap) \
template <class P,class B> \
Xpr1<P, Xpr1ScalarOp<P, Xpr1<P,B>, ap<P> > > \
op (const P& a, const Xpr1<P,B>& b) {\
  typedef Xpr1ScalarOp<P, Xpr1<P,B>, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Dim1 and Xpr1
#define XXX(op,ap) \
template <class P,class A,class B>\
Xpr1<P, Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, Xpr1<P,B>, ap<P> > > \
op (const Dim1<P,A>& a, const Xpr1<P,B>& b) {\
  typedef Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, Xpr1<P,B>, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a), b));\
}

XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Xpr1 and Dim1
#define XXX(op,ap) \
template <class P,class A,class B>\
Xpr1<P, Xpr1BinOp<P, Xpr1<P,A>, ConstRef1<P, Dim1<P,B> >, ap<P> > > \
op (const Xpr1<P,A>& a, const Dim1<P,B>& b) {\
  typedef Xpr1BinOp<P, Xpr1<P,A>, ConstRef1<P, Dim1<P,B> >, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b)));\
}

XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Two Xpr1's
#define XXX(op,ap) \
template <class P, class A, class B>\
Xpr1<P, Xpr1BinOp<P, Xpr1<P,A>, Xpr1<P,B>, ap<P> > >\
op (const Xpr1<P,A>& a, const Xpr1<P,B>& b) {\
  typedef Xpr1BinOp<P, Xpr1<P,A>, Xpr1<P,B>, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a, b));\
}

XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

#endif // XPR1_H_
