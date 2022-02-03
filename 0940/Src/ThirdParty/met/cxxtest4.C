#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

class A {};

template <class P, class A, class B>
class Bin {
  A a;
  B b;
public:
  Bin(const A& a_, const B& b_) : a(a_), b(b_) {}
  P operator()(int i, int j) const { return P(i+j); }
};

/* DIGITAL C++ cannot compile these partical specialization!!!

cxx: Error: cxxtest4.C, line 9: mangled name collision for function types "Bin<double, A, double> &(Bin<double, A, double> *, const Bin<double, A,
          double> &)" and "Bin<double, double, A> &(Bin<double, double, A> *, const
          Bin<double, double, A> &)"
class Bin {
------^

template <class B>
class Bin<double,double,B> {
  double a;
  B b;
public:
  Bin(const double& a_, const B& b_) : a(a_), b(b_) {}
  double operator()(int i, int j) const { return double(i+j); }
};

template <class A>
class Bin<double,A,double> {
  A a;
  double b;
public:
  Bin(const A& a_, const double& b_) : a(a_), b(b_) {}
  double operator()(int i, int j) const { return double(i+j); }
};

*/

template <class P>
class Foo {
  P dat;
public:
  Foo() {}
  template <class A, class B> Foo<P>& operator=(const Bin<P,A,B>& rhs) {
    dat = rhs(0,0);
    return *this;
  }
};

int main(int argc, char * argv[]) {
  Foo<double> f;
  double d;
  A a;
  /* Bin<double,A,double> ad(a,d);
  Bin<double,double,A> da(d,a);

  f = ad;
  f = da; */

  return 0;
}
