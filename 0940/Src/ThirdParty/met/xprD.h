//-*- C++ -*-
//
// Disambiguated Glommarble Expression Templates
// for Digaonal Matrix Calculations
//
// ref. Computers in Physics, 11, 263 (1997)
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

#ifndef XPRD_H_
#define XPRD_H_

#include "xpr1.h"

class DiagRange {
  const int sz, sta, ste;
public:
  DiagRange(int sta_, int end, int ste_=1) : // [ sta_, end )
  sz( (end - sta_)/ste_), sta(sta_), ste(ste_) {}
  DiagRange(const DiagRange& sl) : sz(sl.sz), sta(sl.sta), ste(sl.ste) {}
  //DiagRange& operator=(const DiagRange& x) { sz = x.sz; sta = x.sta; ste=x.ste; }
  DiagRange operator-() const { return DiagRange(sta, sta-sz*ste, -ste); }
  DiagRange operator+(int x) const { return DiagRange(sta+x, sta+sz*ste+x, ste); }
  DiagRange operator-(int x) const { return DiagRange(sta-x, sta+sz*ste-x, ste); }
  DiagRange operator*(int x) const { 
    return DiagRange(sta*x, (sta+sz*ste)*x, ste*x); 
  }
//   DiagRange operator/(int x) const { 
//     if (sta) {
//       return DiagRange(sta/x, (sta+sz*ste)/x, ste/x); 
//     } else {
//       return DiagRange(sta, (sta+sz*ste)/x, ste/x); 
//     }
//   }
  int operator()(int i) const { return sta + ste*i ; }
  int size() const { return sz; }
  int rows() const { return sz; }
  int cols() const { return sz; }
  int start() const { return sta; }
  int step() const { return ste; }
};

template <class P, class A, class F>
class XprDFunc {
  A a;
public:
  XprDFunc(const A& a_) : a(a_) {}
  P operator()(int n) const { return F::apply( a(n)); }
};

template <class P, class F>
class XprDFunc< P, DiagRange, F> {
  const DiagRange& sl;
public:
  XprDFunc(const DiagRange& sl_) : sl(sl_) {}
  P operator()(int n) const { return F::apply( P(sl(n)) ); }
};

template <class P, class A, class B, class Op>
class XprDBinOp {
  A a;
  B b;
public:
  XprDBinOp(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int n) const { return Op::apply( a(n), b(n) ); }
};

// Specialziation for Scalar Addition and Substractin
// (Not Scalar Multiplication or Division)
template <class P, class A, class Op>
class XprDOpScalar {
  A a;
  P b;
public:
  XprDOpScalar(const A& a_, const P& b_) : a(a_), b(b_) {}
  P operator()(int n) const { return Op::apply( a(n), b ); }
};

template <class P, class B, class Op>
class XprDScalarOp {
  P a;
  B b;
public:
  XprDScalarOp(const P& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int n) const { return Op::apply( a, b(n) ); }
};

template <class P, class M>
class ConstRefD {
  const M& m;
public:
  ConstRefD(const M& m_) : m(m_) {}
  P operator()(int n) const { return m(n); }
};

template <class P, class E>
class XprD {
private:
  E e;
public:
  XprD(const E& e_) : e(e_) {}
  P operator()(int n) const { return e(n); }
};

//template <int Nr,int Nc,class P,class I> class TDimD;

template <class P, class I>
class DimD {
public:
  explicit DimD() {}
  int size() const { return static_cast<const I*>(this)->size(); }
  int rows() const { return static_cast<const I*>(this)->size(); }
  int cols() const { return static_cast<const I*>(this)->size(); }
  //P operator() (int n) const {
  //  return static_cast<const I*>(this)->operator()(n);
  //}
  P operator()(int n) const {
    return static_cast<const I*>(this)->operator()(n);
  }
  template <class X> I& assignFrom(const XprD<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->size(); i++) me->operator()(i) = rhs(i);
    return *me;
  }
  template <class M> I& assignFrom(const DimD<P,M>& x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  //template <int Nr,int Nc,class M> I& assignFrom(const TDimD<Nr,Nc,P,M>& x);
  I& assignFrom(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) = x;
    return *me;
  }
  template <class X> DimD<P,I>& operator+=(const XprD<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) += rhs(i);
    return *me;
  }
  template <class M> DimD<P,I>& operator+=(const DimD<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) += rhs(i);
    return *me;
  }
  DimD<P,I>& operator+=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) += x;
    return *me;
  }
  template <class X> DimD<P,I>& operator-=(const XprD<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) -= rhs(i);
    return *me;
  }
  template <class M> DimD<P,I>& operator-=(const DimD<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) -= rhs(i);
    return *me;
  }
  DimD<P,I>& operator-=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) -= x;
    return *me;
  }
  template <class X> DimD<P,I>& operator*=(const XprD<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) *= rhs(i);
    return *me;
  }
  template <class M> DimD<P,I>& operator*=(const DimD<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) *= rhs(i);
    return *me;
  }
  DimD<P,I>& operator*=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) *= x;
    return *me;
  }
  template <class X> DimD<P,I>& operator/=(const XprD<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) /= rhs(i);
    return *me;
  }
  template <class M> DimD<P,I>& operator/=(const DimD<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) /= rhs(i);
    return *me;
  }
  DimD<P,I>& operator/=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) /= x;
    return *me;
  }
};

template <class T,class A>
std::ostream& operator<<(std::ostream& s, const DimD<T,A>& a) {
  for (int i=0; i< a.size(); i++)
    s << std::setw(6) << std::setprecision(3) << a(i) << std::endl;
  return s;
}


// DiagonalMatrix Vector Multiplication
// DimD * Dim1
// DimD * Xpr1
// XprD * Dim1
template <class P, class A, class B>
class Xpr1ReductD {
  A a;
  B b;
public:
  Xpr1ReductD(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i) const { return a(i)*b(i) ; }
};


// Functions of DimD
#define XXX(f,ap) \
template <class P, class A> \
XprD<P, XprDFunc<P, ConstRefD<P,DimD<P,A> >, ap<P> > > \
f(const DimD<P,A>& a) \
{\
   typedef XprDFunc<P, ConstRefD<P,DimD<P,A> >, ap<P> > ExprT;\
   return XprD<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a))); \
}
//XXX(ident, Identity)
XXX(operator- , UnaryMinus)
//XXX(exp, Exp)
#undef XXX

// Functions of DiagRange 
#define XXX(P,r,f,ap) \
XprD< P, XprDFunc<P, r, ap<P> > > \
f(const r& a) \
{\
   typedef XprDFunc<P, r, ap<P> > ExprT;\
   return XprD<P,ExprT>(ExprT(a)); \
}
XXX(double, DiagRange, dble, Identity)
XXX(double, DiagRange, dexp, Exp)
XXX(Complex , DiagRange, cmplx, Identity)
XXX(Complex , DiagRange, cexp, Exp)
#undef XXX

// Functions of XprD
#define XXX(f,ap) \
template <class P, class E> \
XprD<P, XprDFunc<P, XprD<P,E>, ap<P> > > \
f(const XprD<P,E>& a) \
{\
   typedef XprDFunc<P, XprD<P,E>, ap<P> > ExprT;\
   return XprD<P,ExprT>(ExprT(a)); \
}
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX

//Binary operations between Two DimDs
#define XXX(op,ap) \
template <class P,class A,class B> \
XprD<P, XprDBinOp<P, ConstRefD<P, DimD<P,A> >, ConstRefD<P,DimD<P,B> >, ap<P> > >\
op (const DimD<P,A>& a, const DimD<P,B>& b) {\
  typedef \
    XprDBinOp<P, ConstRefD<P,DimD<P,A> >, ConstRefD<P,DimD<P,B> >, ap<P> > \
      ExprT;\
  return XprD<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),\
			      ConstRefD<P,DimD<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX


//Binary operations between DimD and Scalar
#define XXX(op,ap) \
template <class P,class A> \
XprD<P, XprDOpScalar<P, ConstRefD<P, DimD<P,A> >, ap<P> > >\
op (const DimD<P,A>& a, P& b) {\
  typedef \
    XprDOpScalar<P, ConstRefD<P,DimD<P,A> >, ap<P> > \
      ExprT;\
  return XprD<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Scalar and DimD
#define XXX(op,ap) \
template <class P,class B> \
XprD<P, XprDScalarOp<P, ConstRefD<P,DimD<P,B> >, ap<P> > >\
op (const P& a, const DimD<P,B>& b) {\
  typedef \
    XprDScalarOp<P, ConstRefD<P,DimD<P,B> >, ap<P> > \
      ExprT;\
  return XprD<P,ExprT>(ExprT(a,ConstRefD<P,DimD<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

// Multiplication with DimD and Dim1
#define XXX(op) \
template <class P, class A, class B> \
Xpr1<P, Xpr1ReductD<P, ConstRefD<P, DimD<P,A> >, ConstRef1<P,Dim1<P,B> > > >\
op (const DimD<P,A>& a, const Dim1<P,B>& b) {\
  typedef \
    Xpr1ReductD<P, ConstRefD<P,DimD<P,A> >, ConstRef1<P,Dim1<P,B> > > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),\
			      ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

// Muitiplicatino with  DimD and Xpr1
#define XXX(op) \
template <class P,class A,class B> \
Xpr1<P, Xpr1ReductD<P, ConstRefD<P, DimD<P,A> >, Xpr1<P,B> > >\
op (const DimD<P,A>& a, const Xpr1<P,B>& b) {\
  typedef \
    Xpr1ReductD<P, ConstRefD<P, DimD<P,A> >, Xpr1<P,B> > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),b));\
}
XXX(operator*)
#undef XXX

//Binary operations between XprD and Scalar
#define XXX(op,ap) \
template <class P, class A> \
XprD<P, XprDOpScalar<P, XprD<P,A>, ap<P> > >\
op (const XprD<P,A>& a, const P& b) {\
  typedef XprDOpScalar<P, XprD<P,A>, ap<P> > ExprT;\
  return XprD<P,ExprT>(ExprT(a,b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Scalar and XprD
#define XXX(op,ap) \
template <class P, class B> \
XprD<P, XprDScalarOp<P, XprD<P,B>, ap<P> > >\
op (const P& a, const XprD<P,B>& b) {\
  typedef XprDScalarOp<P, XprD<P,B>, ap<P> > ExprT;\
  return XprD<P,ExprT>(ExprT(a,b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

// Multiplication with XprD and Dim1
#define XXX(op) \
template <class P,class A,class B> \
Xpr1<P, Xpr1ReductD<P, XprD<P,A>, ConstRef1<P,Dim1<P,B> > > >\
op (const XprD<P,A>& a, const Dim1<P,B>& b) {\
  typedef Xpr1ReductD<P, XprD<P,A>, ConstRef1<P,Dim1<P,B> > > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

//Binary operations between XprD and DimD
#define XXX(op,ap) \
template <class P,class A,class B> \
XprD<P, XprDBinOp<P, XprD<P,A>, ConstRefD<P,DimD<P,B> >, ap<P> > >\
op (const XprD<P,A>& a, const DimD<P,B>& b) {\
  typedef XprDBinOp<P, XprD<P,A>, ConstRefD<P,DimD<P,B> >, ap<P> > ExprT;\
  return XprD<P,ExprT>(ExprT(a,ConstRefD<P,DimD<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between DimD and XprD
#define XXX(op,ap) \
template <class P,class A,class B> \
XprD<P, XprDBinOp<P, ConstRefD<P, DimD<P,A> >, XprD<P,B>, ap<P> > >\
op (const DimD<P,A>& a, const XprD<P,B>& b) {\
  typedef XprDBinOp<P, ConstRefD<P,DimD<P,A> >, XprD<P,B>, ap<P> > ExprT;\
  return XprD<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between XprD and Xpr1
#define XXX(op) \
template <class P,class A,class B> \
Xpr1<P, Xpr1ReductD<P, XprD<P,A>,Xpr1<P,B> > >\
op (const XprD<P,A>& a, const Xpr1<P,B>& b) {\
  typedef Xpr1ReductD<P, XprD<P,A>,Xpr1<P,B> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,b));\
}
XXX(operator*)
#undef XXX

//Binary operations between Two XprDs
#define XXX(op,ap) \
template <class P, class A, class B> \
XprD<P, XprDBinOp<P, XprD<P,A>, XprD<P,B>, ap<P> > >\
op (const XprD<P,A>& a, const XprD<P,B>& b) {\
  typedef XprDBinOp<P, XprD<P,A>, XprD<P,B>, ap<P> > ExprT;\
  return XprD<P,ExprT>(ExprT(XprD<P,A>(a),XprD<P,B>(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

#endif // XPRD_H_
