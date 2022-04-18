//-*- C++ -*-
//
// Disambiguated Glommarble Expression Templates
// for Matrix Calculations
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

#ifndef XPR2_H_
#define XPR2_H_

#include "xpr1.h"
#include "xprD.h"

class RowRange {
  const int sz, sta, ste;
public:
  RowRange(int sta_, int end, int ste_=1) : // [ sta_, end )
  sz( (end - sta_)/ste_), sta(sta_), ste(ste_) {}
  RowRange(const RowRange& sl) : sz(sl.sz), sta(sl.sta), ste(sl.ste) {}
  //RowRange& operator=(const Range& x) { sz = x.sz; sta = x.sta; ste=x.ste; }
  RowRange operator-() const { return RowRange(sta, sta-sz*ste, -ste); }
  RowRange operator+(int x) const { return RowRange(sta+x, sta+sz*ste+x, ste); }
  RowRange operator-(int x) const { return RowRange(sta-x, sta+sz*ste-x, ste); }
  RowRange operator*(int x) const { 
    return RowRange(sta*x, (sta+sz*ste)*x, ste*x); 
  }
//   RowRange operator/(int x) const { 
//     if (sta) {
//       return RowRange(sta/x, (sta+sz*ste)/x, ste/x); 
//     } else {
//       return RowRange(sta, (sta+sz*ste)/x, ste/x); 
//     }
//   }
  int operator()(int i) const { return sta + ste*i ; }
  int size() const { return sz; }
  int rows() const { return sz; }
  int start() const { return sta; }
  int step() const { return ste; }
};

class ColRange {
  const int sz, sta, ste;
public:
  ColRange(int sta_, int end, int ste_=1) : // [ sta_, end )
  sz( (end - sta_)/ste_), sta(sta_), ste(ste_) {}
  ColRange(const ColRange& sl) : sz(sl.sz), sta(sl.sta), ste(sl.ste) {}
  //ColRange& operator=(const Range& x) { sz = x.sz; sta = x.sta; ste=x.ste; }
  ColRange operator-() const { return ColRange(sta, sta-sz*ste, -ste); }
  ColRange operator+(int x) const { return ColRange(sta+x, sta+sz*ste+x, ste); }
  ColRange operator-(int x) const { return ColRange(sta-x, sta+sz*ste-x, ste); }
  ColRange operator*(int x) const { 
    return ColRange(sta*x, (sta+sz*ste)*x, ste*x); 
  }
//   ColRange operator/(int x) const { 
//     if (sta) {
//       return ColRange(sta/x, (sta+sz*ste)/x, ste/x); 
//     } else {
//       return ColRange(sta, (sta+sz*ste)/x, ste/x); 
//     }
//   }
  int operator()(int i) const { return sta + ste*i ; }
  int size() const { return sz; }
  int cols() const { return sz; }
  int start() const { return sta; }
  int step() const { return ste; }
};

template <class P, class A, class F>
class Xpr2Func {
  A a;
public:
  Xpr2Func(const A& a_) : a(a_) {}
  P operator()(int i, int j) const { return F::apply( a(i,j) ); }
// #ifdef SIZECHECK
//   int rows() const { return a.rows(); }
//   int cols() const { return a.cols(); }
// #endif
};

template <class P, class F>
class Xpr2Func< P, RowRange, F> {
  const RowRange& r;
public:
  Xpr2Func(const RowRange& r_) : r(r_) {}
  P operator()(int i, int j) const { return F::apply( P(r(i)) ); }
};

template <class P, class F>
class Xpr2Func<P, ColRange, F> {
  const ColRange& c;
public:
  Xpr2Func(const ColRange& c_) : c(c_) {}
  P operator()(int i, int j) const { return F::apply( P(c(j)) ); }
};

template <class P, class A, class B, class Op>
class Xpr2BinOp {
  A a;
  B b;
public:
  Xpr2BinOp(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const { return Op::apply( a(i,j), b(i,j) ); }
};

// Specialziation for Scalar Addition and Substractin
// (Not Scalar Multiplication or Division)
template <class P, class A, class Op>
class Xpr2OpScalar {
  A a;
  P b;
public:
  Xpr2OpScalar(const A& a_, const P& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const { return Op::apply( a(i,j), b ); }
};

template <class P, class B, class Op>
class Xpr2ScalarOp {
  P a;
  B b;
public:
  Xpr2ScalarOp(const P& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const { return Op::apply( a, b(i,j) ); }
};


template <class P, class M>
class ConstRef2 {
  const M& m;
public:
  ConstRef2(const M& m_) : m(m_) {}
  P operator()(int i, int j) const { return m(i,j); }
  int cols() const { return m.cols(); }
};

template <class P, class E>
class Xpr2 {
private:
  E e;
public:
  Xpr2(const E& e_) : e(e_) {}
  P operator()(int i, int j) const { return e(i,j); }
};

//template <int Nr,int Nc,class P,class I> class TDim2;

template <class P, class I>
class Dim2 {
private:
  void error(const char *msg) const {
      std::cerr << "Dim2 error: " << msg << std::endl;
    exit(1);
  }

public:
  explicit Dim2() {}
  int size() const { return static_cast<const I*>(this)->size(); }
  int rows() const { return static_cast<const I*>(this)->rows(); }
  int cols() const { return static_cast<const I*>(this)->cols(); }
  //P operator() (int n) const {
  //  return static_cast<const I*>(this)->operator()(n);
  //}
  P operator() (int i, int j) const {
    return static_cast<const I*>(this)->operator()(i,j);
  }
  template <class X> I& assignFrom(const Xpr2<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) = rhs(i,j);
    }
    return *me;
  }
  template <class M> I& assignFrom(const Dim2<P,M>& x) {
    I *me = static_cast<I*>(this);
    //for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) = x(i,j);
    }
    return *me;
  }
  //template <int Nr,int Nc,class M> I& assignFrom(const TDim2<Nr,Nc,P,M>& x);
  I& assignFrom(P x) {  // for Square Matrix Only
    I *me = static_cast<I*>(this);
    int i,j, rsz = me->rows(), csz = me->cols(),
      //min = ( rsz < csz ) ? rsz : csz ;
      min = rsz;
    for (i=0; i<min; i++) {
      for (j=0; j<i; j++) me->operator()(i,j) = P(0);
      me->operator()(i,i) = x;
      for (j=i+1; j<csz; j++) me->operator()(i,j) = P(0);
    }
    //for (i=min; i<rsz; i++) {
    //  for (j=0; j<csz; j++) me->operator()(i,j) = P(0);
    //}
    return *me;
  }
  template <class X> I& assignFrom(const XprD<P,X>& rhs) {  // for Square Matrix Only
    I *me = static_cast<I*>(this);
    int i,j, rsz = me->rows(), csz = me->cols(),
      //min = ( rsz < csz ) ? rsz : csz ;
      min = rsz;
    for (i=0; i<min; i++) {
      for (j=0; j<i; j++) me->operator()(i,j) = P(0);
      me->operator()(i,i) = rhs(i);
      for (j=i+1; j<csz; j++) me->operator()(i,j) = P(0);
    }
    //for (i=min; i<rsz; i++) {
    //  for (j=0; j<csz; j++) me->operator()(i,j) = T(0);
    //}
    return *me;
  }
  template <class M> I& assignFrom(const DimD<P,M>& rhs) {  // for Square Matrix Only
    I *me = static_cast<I*>(this);
    int i,j, rsz = me->rows(), csz = me->cols(),
      //min = ( rsz < csz ) ? rsz : csz ;
      min = rsz;
    for (i=0; i<min; i++) {
      for (j=0; j<i; j++) me->operator()(i,j) = P(0);
      me->operator()(i,i) = rhs(i);
      for (j=i+1; j<csz; j++) me->operator()(i,j) = P(0);
    }
    //for (i=min; i<rsz; i++) {
    //  for (j=0; j<csz; j++) me->operator()(i,j) = T(0);
    //}
    return *me;
  }
  template <class X> Dim2<P,I>& operator+=(const Xpr2<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) += rhs(i,j);
    }
    return *me;
  }
  template <class M> Dim2<P,I>& operator+=(const Dim2<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    //for (int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) += rhs(i,j);
    }
    return *me;
  }
  Dim2<P,I>& operator+=(P x) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    //int rsz = me->rows(), csz = me->cols(),
    //min = ( rsz < csz ) ? rsz : csz ;
    for (int i=0; i<me->rows(); i++) me->operator()(i,i) += x;
    return *me;
  }
  template <class X> Dim2<P,I>& operator+=(const XprD<P,X>& rhs) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    //int rsz = me->rows(), csz = me->cols(),
    //min = ( rsz < csz ) ? rsz : csz ;
    for (int i=0; i<me->rows(); i++) me->operator()(i,i) += rhs(i);
    return *me;
  }
  template <class M> Dim2<P,I>& operator+=(const DimD<P,M>& rhs) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    //int rsz = me->rows(), csz = me->cols(),
    //min = ( rsz < csz ) ? rsz : csz ;
    for (int i=0; i<me->rows(); i++) me->operator()(i,i) += rhs(i);
    return *me;
  }
  template <class X> Dim2<P,I>& operator-=(const Xpr2<P,X>& rhs) {
      I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) -= rhs(i,j);
    }
    return *me;
  }
  template <class M> Dim2<P,I>& operator-=(const Dim2<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    //for (int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) -= rhs(i,j);
    }
    return *me;
  }
  Dim2<P,I>& operator-=(P x) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    //int rsz = me->rows(), csz = me->cols(),
    //min = ( rsz < csz ) ? rsz : csz ;
    for (int i=0; i<me->rows(); i++) me->operator()(i,i) -= x;
    return *me;
  }
  template <class X> Dim2<P,I>& operator-=(const XprD<P,X>& rhs) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    //int rsz = me->rows(), csz = me->cols(),
    //min = ( rsz < csz ) ? rsz : csz ;
    for (int i=0; i<me->rows(); i++) me->operator()(i,i) -= rhs(i);
    return *me;
  }
  template <class M> Dim2<P,I>& operator-=(const DimD<P,M>& rhs) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    //int rsz = me->rows(), csz = me->cols(),
    //min = ( rsz < csz ) ? rsz : csz ;
    for (int i=0; i<me->rows(); i++) me->operator()(i,i) -= rhs(i);
    return *me;
  }
  Dim2<P,I>& operator*=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) *= x;
    }
    return *me;
  }
  template <class X> Dim2<P,I>& operator*=(const XprD<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) *= rhs(j);
    }
    return *me;
  }
  template <class M> Dim2<P,I>& operator*=(const DimD<P,M>& rhs) { // for Square Matrix Only
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) *= rhs(j);
    }
    return *me;
  }
  Dim2<P,I>& operator/=(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i<me->rows(); i++) {
      for (int j=0; j<me->cols(); j++) me->operator()(i,j) /= x;
    }
    return *me;
  }
};

template <class T,class A>
std::ostream& operator<<(std::ostream& s, const Dim2<T,A>& a) {
  for (int i=0; i< a.rows(); i++) {
    for (int j=0; j< a.cols(); j++)
        s << std::setw(6) << std::setprecision(3) << a(i,j);
    s << std::endl;
  }
  return s;
}


// Matrix Vector Multiplication
// Dim2 * Dim1
// Dim2 * Xpr1
// Xpr2 * Dim1
template <class P, class A, class B>
class Xpr1Reduct {
  A a;
  B b;
public:
  Xpr1Reduct(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i) const {
    P sum = a(i,0)*b(0);
    for (int j=1; j < a.cols(); j++) sum += a(i,j)*b(j);
    return sum;
  }
};

// Matrix Matrix Multiplication
// Dim2 * Dim2
// Dim2 * Xpr2
// Xpr2 * Dim2
// Xpr2 * Xpr2
template <class P, class A, class B>
class Xpr2Reduct {
  A a;
  B b;
public:
  Xpr2Reduct(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const {
    P sum = a(i,0)*b(0,j);
    for (int k=1; k < a.cols(); k++) sum += a(i,k)*b(k,j);
    return sum;
  }
};

// DiagonalMatrix Matrix Multiplication
// DimD * Dim2
// DimD * Xpr2
// XprD * Dim2
// XprD * Xpr2
template <class P, class A, class B>
class Xpr2ReductD2 {
  A a;
  B b;
public:
  Xpr2ReductD2(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const { return a(i)*b(i,j); }
};

// Matrix DiagonalMatrix Multiplication
// Dim2 * DimD
// Dim2 * XprD
// Xpr2 * DimD
// Xpr2 * XprD
template <class P, class A, class B>
class Xpr2Reduct2D {
  A a;
  B b;
public:
  Xpr2Reduct2D(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const { return a(i,j)*b(j); }
};


// Functions of Dim2
#define XXX(f,ap) \
template <class P, class A> \
Xpr2<P, Xpr2Func<P, ConstRef2<P,Dim2<P,A> >, ap<P> > > \
f(const Dim2<P,A>& a) \
{\
   typedef Xpr2Func<P, ConstRef2<P,Dim2<P,A> >, ap<P> > ExprT;\
   return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a))); \
}
//XXX(ident, Identity)
XXX(operator- , UnaryMinus)
//XXX(exp, Exp)
#undef XXX

// Functions of RowRange and ColRange
#define XXX(P,r,f,ap) \
Xpr2< P, Xpr2Func<P, r, ap<P> > > \
f(const r& a) \
{\
   typedef Xpr2Func<P, r, ap<P> > ExprT;\
   return Xpr2<P,ExprT>(ExprT(a)); \
}
XXX(double, RowRange, dble, Identity)
XXX(double, RowRange, dexp, Exp)
XXX(Complex , RowRange, cmplx, Identity)
XXX(Complex , RowRange, cexp, Exp)
XXX(double, ColRange, dble, Identity)
XXX(double, ColRange, dexp, Exp)
XXX(Complex , ColRange, cmplx, Identity)
XXX(Complex , ColRange, cexp, Exp)
#undef XXX

// Functions of Xpr2
#define XXX(f,ap) \
template <class P, class E> \
Xpr2<P, Xpr2Func<P, Xpr2<P,E>, ap<P> > > \
f(const Xpr2<P,E>& a) \
{\
   typedef Xpr2Func<P, Xpr2<P,E>, ap<P> > ExprT;\
   return Xpr2<P,ExprT>(ExprT(a)); \
}
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX

//Binary operations between Two Dim2s
#define XXX(op,ap) \
template <class P,class A,class B> \
Xpr2<P, Xpr2BinOp<P, ConstRef2<P, Dim2<P,A> >, ConstRef2<P,Dim2<P,B> >, ap<P> > >\
op (const Dim2<P,A>& a, const Dim2<P,B>& b) {\
  typedef \
    Xpr2BinOp<P, ConstRef2<P,Dim2<P,A> >, ConstRef2<P,Dim2<P,B> >, ap<P> > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),\
			      ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

// Multiplication with Two Dim2s
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct<P, ConstRef2<P, Dim2<P,A> >, ConstRef2<P,Dim2<P,B> > > >\
op (const Dim2<P,A>& a, const Dim2<P,B>& b) {\
  typedef \
    Xpr2Reduct<P, ConstRef2<P,Dim2<P,A> >, ConstRef2<P,Dim2<P,B> > > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),\
			      ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

//Binary operations between Dim2 and Xpr2
#define XXX(op,ap) \
template <class P,class A,class B> \
Xpr2<P, Xpr2BinOp<P, ConstRef2<P, Dim2<P,A> >, Xpr2<P,B>, ap<P> > >\
op (const Dim2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2BinOp<P, ConstRef2<P,Dim2<P,A> >, Xpr2<P,B>, ap<P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

/* Not Efficient !! // Muitiplicatino with  Dim2 and Xpr2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct<P, ConstRef2<P, Dim2<P,A> >, Xpr2<P,B> > >\
op (const Dim2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef \
    Xpr2Reduct<P, ConstRef2<P, Dim2<P,A> >, Xpr2<P,B> > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
XXX(operator*)
#undef XXX */

//Binary operations between Dim2 and Scalar
#define XXX(op,ap) \
template <class P,class A> \
Xpr2<P, Xpr2OpScalar<P, ConstRef2<P, Dim2<P,A> >, ap<P> > >\
op (const Dim2<P,A>& a, P& b) {\
  typedef \
    Xpr2OpScalar<P, ConstRef2<P,Dim2<P,A> >, ap<P> > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
  //XXX(operator+, OpAdd) //impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

// impossible ! // Binary Operations between Dim2 and Dim1

// Multiplication with Dim2 and Dim1
#define XXX(op) \
template <class P, class A, class B> \
Xpr1<P, Xpr1Reduct<P, ConstRef2<P, Dim2<P,A> >, ConstRef1<P,Dim1<P,B> > > >\
op (const Dim2<P,A>& a, const Dim1<P,B>& b) {\
  typedef \
    Xpr1Reduct<P, ConstRef2<P,Dim2<P,A> >, ConstRef1<P,Dim1<P,B> > > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),\
			      ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

// impossible ! // Binary Operations between Dim2 and Xpr1

/* Not Efficient !! // Muitiplicatino with  Dim2 and Xpr1
#define XXX(op) \
template <class P,class A,class B> \
Xpr1<P, Xpr1Reduct<P, ConstRef2<P, Dim2<P,A> >, Xpr1<P,B> > >\
op (const Dim2<P,A>& a, const Xpr1<P,B>& b) {\
  typedef \
    Xpr1Reduct<P, ConstRef2<P, Dim2<P,A> >, Xpr1<P,B> > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
XXX(operator*)
#undef XXX */

// impossible ! //Binary Operations between Dim2 and DimD

//Multiplication with Dim2 and DimD
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct2D<P, ConstRef2<P, Dim2<P,A> >, ConstRefD<P,DimD<P,B> > > >\
op (const Dim2<P,A>& a, const DimD<P,B>& b) {\
  typedef \
    Xpr2Reduct2D<P, ConstRef2<P,Dim2<P,A> >, ConstRefD<P,DimD<P,B> > > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),\
			      ConstRefD<P,DimD<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

// impossible ! //Binary Operations between Dim2 and XprD

//Multiplication with Dim2 and XprD
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct2D<P, ConstRef2<P, Dim2<P,A> >, XprD<P,B> > > \
op (const Dim2<P,A>& a, const XprD<P,B>& b) {\
  typedef \
    Xpr2Reduct2D<P, ConstRef2<P,Dim2<P,A> >, XprD<P,B> > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
XXX(operator*)
#undef XXX

/* Not Efficient! //Multiplication with Xpr2 and Dim2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> > > >\
op (const Xpr2<P,A>& a, const Dim2<P,B>& b) {\
  typedef \
    Xpr2Reduct<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> > > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(a, ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator*)
#undef XXX */

//Binary operations between Two Xpr2s
#define XXX(op,ap) \
template <class P, class A, class B> \
Xpr2<P, Xpr2BinOp<P, Xpr2<P,A>, Xpr2<P,B>, ap<P> > >\
op (const Xpr2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2BinOp<P, Xpr2<P,A>, Xpr2<P,B>, ap<P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(Xpr2<P,A>(a),Xpr2<P,B>(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

/* Not Efficient! //Multiplication with Two Xpr2s
#define XXX(op,ap) \
template <class P, class A, class B> \
Xpr2<P, Xpr2Reduct<P, Xpr2<P,A>, Xpr2<P,B> > > >\
op (const Xpr2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2Reduct<P, Xpr2<P,A>, Xpr2<P,B> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(Xpr2<P,A>(a),Xpr2<P,B>(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX */

//Binary operations between Xpr2 and Scalar
#define XXX(op,ap) \
template <class P, class A> \
Xpr2<P, Xpr2OpScalar<P, Xpr2<P,A>, ap<P> > >\
op (const Xpr2<P,A>& a, const P& b) {\
  typedef Xpr2OpScalar<P, Xpr2<P,A>, ap<P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a,b));\
}
  //XXX(operator+, OpAdd) // impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

// impossible! //Binary Operations between Xpr2 and Dim1
// impossible! //Multiplications with Xpr2 and Dim1

// impossible! //Binary Operations between Xpr2 and Xpr1
// impossible! //Multiplications with Xpr2 and Xpr1


/* Not Efficient !! // Multiplication with Xpr2 and Dim1
#define XXX(op) \
template <class P,class A,class B> \
Xpr1<P, Xpr1Reduct<P, Xpr2<P,A>, ConstRef1<P,Dim1<P,B> > > >\
op (const Xpr2<P,A>& a, const Dim1<P,B>& b) {\
  typedef Xpr1Reduct<P, Xpr2<P,A>, ConstRef1<P,Dim1<P,B> > > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator*)
#undef XXX */

//Binary operations between Xpr2 and Dim2
#define XXX(op,ap) \
template <class P,class A,class B> \
Xpr2<P, Xpr2BinOp<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> >, ap<P> > >\
op (const Xpr2<P,A>& a, const Dim2<P,B>& b) {\
  typedef Xpr2BinOp<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> >, ap<P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a, ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

/* Not Efficient! //Mutliplication with Xpr2 and Dim2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> > > >\
op (const Xpr2<P,A>& a, const Dim2<P,B>& b) {\
  typedef Xpr2Reduct<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> > > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a,ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator*)
#undef XXX */

// impossible ! // Binary Operations between Xpr2 and Dim1

/* Not Efficient // Multiplication with Xpr2 and Dim1
#define XXX(op) \
template <class P, class A, class B> \
Xpr1<P, Xpr1Reduct<P, Xpr2<P,A>, ConstRef1<P,Dim1<P,B> > > >\
op (const Xpr2<P,A>& a, const Dim1<P,B>& b) {\
  typedef \
    Xpr1Reduct<P, Xpr2<P,A>, ConstRef1<P,Dim1<P,B> > > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(a, ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator*)
#undef XXX */

// impossible ! // Binary Operations between Xpr2 and Xpr1

/* Not Efficient !! // Muitiplicatino with  Xpr2 and Xpr1
#define XXX(op) \
template <class P,class A,class B> \
Xpr1<P, Xpr1Reduct<P, Xpr2<P,A>, Xpr1<P,B> > >\
op (const Xpr2<P,A>& a, const Xpr1<P,B>& b) {\
  typedef \
    Xpr1Reduct<P, Xpr2<P,A>, Xpr1<P,B> > \
      ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,b));\
}
XXX(operator*)
#undef XXX */

// impossible ! //Binary Operations between Xpr2 and DimD

//Multiplication with Xpr2 and DimD
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct2D<P, Xpr2<P,A>, ConstRefD<P,DimD<P,B> > > >\
op (const Xpr2<P,A>& a, const DimD<P,B>& b) {\
  typedef \
    Xpr2Reduct2D<P, Xpr2<P,A>, ConstRefD<P,DimD<P,B> > > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(a, ConstRefD<P,DimD<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

// impossible ! //Binary Operations between Xpr2 and XprD

//Multiplication with Xpr2 and XprD
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct2D<P, Xpr2<P,A>, XprD<P,B> > > \
op (const Xpr2<P,A>& a, const XprD<P,B>& b) {\
  typedef Xpr2Reduct2D<P, Xpr2<P,A>, XprD<P,B> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a,b));\
}
XXX(operator*)
#undef XXX

//Binary operations between Scalar and Dim2
#define XXX(op,ap) \
template <class P,class B> \
Xpr2<P, Xpr2ScalarOp<P, ConstRef2<P,Dim2<P,B> >, ap<P> > >\
op (const P& a, const Dim2<P,B>& b) {\
  typedef \
    Xpr2ScalarOp<P, ConstRef2<P,Dim2<P,B> >, ap<P> > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(a,ConstRef2<P,Dim2<P,B> >(b)));\
}
  //XXX(operator+, OpAdd) //impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Scalar and Xpr2
#define XXX(op,ap) \
template <class P, class B> \
Xpr2<P, Xpr2ScalarOp<P, Xpr2<P,B>, ap<P> > >\
op (const P& a, const Xpr2<P,B>& b) {\
  typedef Xpr2ScalarOp<P, Xpr2<P,B>, ap<P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a,b));\
}
  //XXX(operator+, OpAdd) // impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

// ALREADY DEFINED! //Binary operations between Scalars
// Already Defined! //Binary operations between Scalar and Dim1
// Already Defined! //Binary operations between Scalar and Xpr1
// Already Defined! //Binary operations between Scalar and DimD
// Already Defined! //Binary operations between Scalar and XprD

// impossible ! //Binary operations between Dim1 and Dim2
// impossible ! //Multiplication with Dim1 and Dim2
// impossible ! //Binary operations between Dim1 and Xpr2
// impossible ! //Multiplication with Dim1 and Xpr2
// Already Defined! //Binary operations between Dim1 and Scalar
// Already Defined! //Binary operations between Dim1 and Dim1
// Already Defined! //Binary operations between Dim1 and Xpr1
// Already Defined! //Binary operations between Dim1 and DimD
// Already Defined! //Binary operations between Dim1 and XprD

// impossible ! //Binary operations between Xpr1 and Dim2
// impossible ! //Multiplication with Xpr1 and Dim2
// impossible ! //Binary operations between Xpr1 and Xpr2
// impossible ! //Multiplication with Xpr1 and Xpr2
// Already Defined! //Binary operations between Xpr1 and Scalar
// Already Defined! //Binary operations between Xpr1 and Dim1
// Already Defined! //Binary operations between Xpr1 and Xpr1
// Already Defined! //Binary operations between Xpr1 and DimD
// Already Defined! //Binary operations between Xpr1 and XprD


// impossible ! //Binary Operations between DimD and Dim2

//Multiplication with DimD and Dim2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2ReductD2<P, ConstRefD<P, DimD<P,A> >, ConstRef2<P,Dim2<P,B> > > >\
op (const DimD<P,A>& a, const Dim2<P,B>& b) {\
  typedef \
    Xpr2ReductD2<P, ConstRefD<P,DimD<P,A> >, ConstRef2<P,Dim2<P,B> > > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),\
			      ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

// impossible ! //Binary Operations between DimD and Xpr2

//Multiplication with DimD and Xpr2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2ReductD2<P, ConstRefD<P, DimD<P,A> >, Xpr2<P,B> > > \
op (const DimD<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2ReductD2<P, ConstRefD<P,DimD<P,A> >, Xpr2<P,B> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRefD<P,DimD<P,A> >(a),b));\
}
XXX(operator*)
#undef XXX

// Already Defined! //Binary operations between DimD and Scalar

// impossible ! //Binary Operations between DimD and Dim1
// impossible ! //Binary Operations between DimD and Xpr1
// Already Defined! //Binary operations between DimD and DimD
// Already Defined! //Binary operations between DimD and XprD


// impossible ! //Binary Operations between XprD and Dim2

//Multiplication with XprD and Dim2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2ReductD2<P, XprD<P,A>, ConstRef2<P,Dim2<P,B> > > >\
op (const XprD<P,A>& a, const Dim2<P,B>& b) {\
  typedef Xpr2ReductD2<P, XprD<P,A>, ConstRef2<P,Dim2<P,B> > > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a, ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator*)
#undef XXX

// impossible ! //Binary Operations between XprD and Xpr2

//Multiplication with XprD and Xpr2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2ReductD2<P, XprD<P,A>, Xpr2<P,B> > >\
op (const XprD<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2ReductD2<P, XprD<P,A>, Xpr2<P,B> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a,b));\
}
XXX(operator*)
#undef XXX

// Already Defined! //Binary operations between XprD and Scalar

// impossible ! //Binary Operations between XprD and Dim1
// impossible ! //Binary Operations between XprD and Xpr1
// Already Defined! //Binary operations between XprD and DimD
// Already Defined! //Binary operations between XprD and XprD

#endif // XPR2_H_

