//-*- C++ -*-
//
// Expression Templates
// (Not Glommable nor Disambiguated)
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
  Range operator*(int x) const { return Range(sta, sta+sz*ste, ste*x); }
  Range operator/(int x) const { return Range(sta, sta+sz*ste, ste/x); }
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
};

template <class P, class F>
class Xpr1Func<P,Range,F> {
  //const Range sl;
  const Range& sl;
public:
  Xpr1Func(const Range& sl_) : sl(sl_) {}
  P operator()(int n) const { return F::apply( P(sl(n)) ); }
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
};

// Specialziation for the all operations with Scalar
template <class P, class A, class Op>
class Xpr1BinOp<P,A,P,Op> {
  A a;
  P b;
public:
  Xpr1BinOp(const A& a_, const P& b_) : a(a_), b(b_) {}
  P operator() (int n) const { return Op::apply( a(n), b ); }
};

// Specialziation for the +,- and * operations with Scalar 
template <class P, class B, class Op>
class Xpr1BinOp<P,P,B,Op> {
  P a;
  B b;
public:
  Xpr1BinOp(const P& a_, const B& b_) : a(a_), b(b_) {}
  P operator() (int n) const { return Op::apply( a, b(n) ); }
};

template <class P, class V>
class ConstRef1 {
  const V& v;
public:
  ConstRef1(const V& v_) : v(v_) {}
  P operator()(int n) const { return v(n); }
};


template <class P, class E>
class Xpr1 {
private:
  E e;
public:
  Xpr1(const E& e_) : e(e_) {}
  P operator() (int n) const { return e(n); }
};

template <int N, class P, class I> class TDim1;

// 1 Dimensional Array Base Class
// for Glommable Expression Templates
template <class P, class I>
class Dim1 { 
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
    for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  template <class V> I& assignFrom(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  I& assignFrom(P x) {
    I *me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) = x;
    return *me;
  }
  template <int N,class V> I& assignFrom(const TDim1<N,P,V>& x);
  /* template <class E> Dim1<P,I>& operator+=(const Xpr1<P,E>& x) {
    I* me = static_cast<I*>(this);
    for (int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    return *me;
  }
  template <class V> Dim1<P,I>& operator+=(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
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
    for (int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    return *me;
  }
  template <class V> Dim1<P,I>& operator-=(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
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
  } */
};

template <class T,class A>
ostream& operator<<(ostream& s, const Dim1<T,A>& a) {
  for (int i=0; i< a.size(); i++) 
    s << setw(6) << setprecision(3) << a(i) << endl;
  return s;
}

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

/* // Tiny 1 Dimensional Array Base Class
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
    Metaloop1<N,I,Xpr<P,E> >::assign(me,x);
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
} */

/* // Functions of Dim1
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
Xpr1<P, Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, P, ap<P> > > \
op (const Dim1<P,A>& a, const P& b) {\
  typedef Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, P, ap<P> > ExprT;\
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
Xpr1<P, Xpr1BinOp<P, P, ConstRef1<P, Dim1<P,B> >, ap<P> > > \
op (const P& a, const Dim1<P,B>& b) {\
  typedef Xpr1BinOp<P, P, ConstRef1<P, Dim1<P,B> >, ap<P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Xpr1 and Scalar
#define XXX(op,ap) \
template <class P,class A>\
Xpr1<P, Xpr1BinOp<P, Xpr1<P,A>, P, ap<P> > > \
op (const Xpr1<P,A>& a, const P& b) {\
  typedef Xpr1BinOp<P, Xpr1<P,A>, P, ap<P> > ExprT;\
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
Xpr1<P, Xpr1BinOp<P, P, Xpr1<P,B>, ap<P> > > \
op (const P& a, const Xpr1<P,B>& b) {\
  typedef Xpr1BinOp<P, P, Xpr1<P,B>, ap<P> > ExprT;\
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
#undef XXX */

#endif // XPR1_H_

#ifndef VECMAT_H_
#define VECMAT_H_

#include <assert.h>

//#include "xpr1.h"
//#include "xpr2.h"

template <class T> class RangeVec;

template <class T>
class Vec : public Dim1<T,Vec<T> > {
private:
  const int sz;
  T* v;
public:
  Vec(int sz_) : Dim1<T,Vec<T> >(), sz(sz_), v(new T[sz]) {}
  ~Vec() { delete [] v; }
  template <class X> Vec<T>& operator=(const Xpr1<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> Vec<T>& operator=(const Dim1<T,V>& rhs) {
    return assignFrom(rhs);
  }
  Vec<T>& operator=(const Vec<T>& rhs) { 
    return assignFrom(rhs); 
  }
  Vec<T>& operator=(T rhs) { 
    return assignFrom(rhs); 
  }
  template <int N,class V> Vec<T>& operator=(const TDim1<N,T,V>& rhs) {
    return assignFrom(rhs);
  }  
  template <class Closure> Vec<T>& operator=(const Closure& rhs) { 
    rhs.assignTo(v);
    return *this;
  }
  Vec(const Vec<T>& rhs) : Dim1<T,Vec<T> >(), sz(x.sz), v(new T[sz]) { 
    (*this) = x; 
  }
  int size() const { return sz; }
  T& operator()(int i) { return v[i]; }
  T operator()(int i) const { return v[i]; }
  // RangeVec<T> operator()(Range sl) const;
};

/* template <class T>
class RangeVec : public Dim1<T,RangeVec<T> > {
  const Range sl;
  T* const dat;
public:
  RangeVec(Range sl_,T* dat_) : Dim1<T,RangeVec<T> >(), sl(sl_), dat(dat_) {}
  template <class X> RangeVec<T>& operator=(const Xpr1<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> RangeVec<T>& operator=(const Dim1<T,V>& rhs) {
    return assignFrom(rhs);
  }
  template <int N,class V> RangeVec<T>& assignFrom(const TDim1<N,T,V>& rhs) {
    return assignFrom(rhs);
  }
  RangeVec<T>& operator=(T rhs) { return assignFrom(rhs); }
  int size() const { return sl.size(); }
  T& operator()(int i) { return dat[sl(i)]; }
  T operator()(int i) const { return dat[sl(i)]; }
};

template <class T>
RangeVec<T> Vec<T>::operator()(Range sl) const { 
  return RangeVec<T>(sl,v); 
} */

/* template <class T> class RangeMat;

template <class T>
class Mat : public Dim2<T,Mat<T> > {
private:
  const int rsz, csz, tsz;
  T* const dat;
  T** m;
public:
  Mat(int rsz_, int csz_) : Dim2<T,Mat<T> >(),
  rsz(rsz_), csz(csz_), tsz(rsz*csz), 
  dat(new T[tsz]), m(new T*[rsz]) {
    for (int i=0; i<rsz; i++) m[i] = dat + i*csz;
  }
  ~Mat() { 
    delete [] m;
    delete [] dat; 
  }
  template <class X> Mat<T>& operator=(const Xpr2<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class M> Mat<T>& operator=(const Dim2<T,M>& rhs) {
    return assignFrom(rhs);
  }
  template <int Nr,int Nc,class M> 
    Mat<T>& assignFrom(const TDim2<Nr,Nc,T,M>& rhs) { return assignFrom(rhs); }
  Mat<T>& operator=(T rhs) { return assignFrom(rhs); }
  template <class Closure> Mat<T>& operator=(const Closure& rhs) { 
    rhs.assignTo(m);
    return *this;
  }
  int size() const { return tsz; }
  int rows() const { return rsz; }
  int cols() const { return csz; }
  //T& operator()(int n) { return dat[n]; }
  //T operator()(int n) const { return dat[n]; }
  T& operator()(int i, int j) { return m[i][j]; }
  T operator()(int i, int j) const { return m[i][j]; }
  RangeVec<T> operator()(RowRange rsl, int j) const {
    return RangeVec<T>(Range(rsl.start()*csz+j,rsl.size()*rsl.step()*csz+j,
			     rsl.step()*csz),
		       v);
  }
  RangeVec<T> operator()(int i, ColRange csl) const {
    return RangeVec<T>(Range(i*csz+csl.start(),i*csz+csl.size()*csl.step(),
			     csl.step()),
		       v);
  }
  RangeMat<T> operator()(RowRange rsl, ColRange csl) const;
};

template <class T>
class RangeMat : public Dim2<T,RangeMat<T> > {
private:
  const Range rsl ,csl;
  T** const  m;
public:
  RangeMat(RowRange rsl_, ColRange csl_, T** m_) :
  Dim2<T,RangeMat<T> >(), rsl(rsl_), csl(csl_), m(m_) {}
  template <class X> RangeMat<T>& operator=(const Xpr2<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class M> RangeMat<T>& operator=(const Dim2<T,M>& rhs) {
    return assignFrom(rhs);
  }
  template <int Nr,int Nc,class M> 
    RangeMat<T>& assignFrom(const TDim2<Nr,Nc,T,M>& rhs) { 
      return assignFrom(rhs); 
    }
  RangeVec<T>& operator=(T rhs) { return assignFrom(rhs); }
  int size() const { return sl.size(); }
  T& operator()(int i, int j) { return m[rsl(i)][csl(j)]; }
  T operator()(int i, int j) const { return m[rsl(i)][csl(j)]; }
};

template <class T>
RangeMat<T> Mat<T>::operator()(RowRange rsl, ColRange csl) const {
  return RangeMat<T>(rsl,csl,m);
}; */

#endif // VECMAT_H_

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

//#include "vecmat.h"

int main(int argc, char * argv[]) {
  Vec<double> a(7), b(7), c(7), x(7);
  Range i(0,7,1);
  double err, de;
  int j;

  a = 3.0;
  /* a = b + c;

  a(i) = dble(-i + 7);
  cout << "a = " << endl;
  cout << a;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = double(-j+7) - a(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  b(i) = -dble(i*4)*2.0 + dexp(i);
  //b = -a*3.0;
  cout << "b = " << endl;
  cout << b;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = -j*4.0*2.0 + exp(j) - b(j);
    //de = b(j) + a(j)*3.0;
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  c = dble(i*3);
  cout << "c = " << endl;
  cout << c;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = j*3.0 - c(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  x = -a + b - c; // -2*i
  cout << " x = -a + b - c =" << endl;
  cout << x;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = -double(-j+7) - j*4.0*2.0 + exp(j) - j*3.0 - x(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  x = exp(a);
  cout << "x = exp(a);" << endl;
  cout << x;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = exp(a(j)) - x(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  Range sl02(0,2);
  x(sl02) = dexp( sl02 );
  cout << "exp( Range(1.0,2.0) ) = " << endl;
  cout << x(sl02);
  err = 0.0;
  for (j=0; j<2; j++) {
    de = exp(double(sl02(j))) - x(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10); */

  return 0;
}
