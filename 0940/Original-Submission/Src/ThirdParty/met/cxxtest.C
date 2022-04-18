#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

template <class T, class A>
class Hoge {};

template <int N, class T, class I> class Foo1 {
public:
  explicit Foo1() {}
  template <class X> I& assignFrom(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    int i = N;
    return *me;
  }
};
template <class T, class I> class Foo2 {
public:
  explicit Foo2() {}
  template <class X> I& assignFrom(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
};
template <class T, class I> class Foo3 {
public:
  explicit Foo3() {}
  template <class X> I& assignFrom(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
};
template <class T, class I> class Foo4 {
public:
  explicit Foo4() {}
  template <class X> I& assignFrom(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
};

template <class T, class I>
class Foo {
public:
  explicit Foo() {}
  template <class X> I& assignFrom(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <int N,class V> I& assignFrom(const Foo1<N,T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> I& assignFrom(const Foo2<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> I& assignFrom(const Foo3<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> I& assignFrom(const Foo4<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> I& assignFrom(const Foo<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class X> Foo<T,I>& operator+=(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <int N,class V> Foo<T,I>& operator+=(const Foo1<N,T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator+=(const Foo2<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator+=(const Foo3<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator+=(const Foo4<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator+=(const Foo<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class X> Foo<T,I>& operator-=(const Hoge<T,X>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <int N,class V> Foo<T,I>& operator-=(const Foo1<N,T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator-=(const Foo2<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator-=(const Foo3<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator-=(const Foo4<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
  template <class V> Foo<T,I>& operator-=(const Foo<T,V>& rhs) {
    I *me = static_cast<I*>(this);
    return *me;
  }
};

template <class T>
class Bar : public Foo<T,Bar<T> > {
public:
  Bar() : Foo<T,Bar<T> >() {}
  template <class X> Bar<T>& operator=(const Hoge<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> Bar<T>& operator=(const Foo<T,V>& rsh) {
    return assignFrom(rhs);
  }
  template <int N,class V> Bar<T>& operator=(const Foo1<N,T,V>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> Bar<T>& operator=(const Foo2<T,V>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> Bar<T>& operator=(const Foo3<T,V>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> Bar<T>& operator=(const Foo4<T,V>& rhs) {
    return assignFrom(rhs);
  }
};

int main(int argc, char * argv[]) {
  Hoge<double, Foo<double, Bar<double> > > h;
  Bar<double> b;

  b = h;

  return 0;
}


