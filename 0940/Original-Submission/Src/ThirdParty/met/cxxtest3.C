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
  template <int N,class V> Vec<T>& operator=(const TDim1<N,T,V>& rhs) {
    return assignFrom(rhs);
  }  
  Vec<T>& operator=(T rhs) { return assignFrom(rhs); }
  int size() const { return sz; }
  T& operator()(int i) { return v[i]; }
  T operator()(int i) const { return v[i]; }
  // RangeVec<T> operator()(Range sl) const;
};

int main(int argc, char * argv[]) {
  Vec<double> a(7);

  a = 3.0;

  return 0;
}
