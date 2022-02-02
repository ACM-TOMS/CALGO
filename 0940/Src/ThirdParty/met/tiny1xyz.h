//-*- C++ -*-

#ifdef TINYXYZ_H_
#define TINYXYZ_H_

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>

#include "error.h"
#include "cmplxtype.h"

template <class P, class I> class BaseMat3;

// XYZ Component Vector Base Class (Glommable)
template <class P, class I>
class BaseVec3 {
public:
  P x,y,z;

protected:
  BaseVec3(P x_, P y_, P z_) : x(x_), y(y_), z(z_) {}
  template <class V> void assignFrom(const BaseVec3<P,V>& rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
  }
  void assignFrom(P rhs) {
    x = rhs; y = rhs; z = rhs;
  }

public:
  template <class V>
  I& operator+=(const BaseVec3<P,V>& rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *(static_cast<I*>(this));
  }
  I& operator+=(P rhs) {
    x += rhs; y += rhs; z += rhs;
    return *(static_cast<I*>(this));
  }
  I& operator-=(const BaseVec3<P,V>& rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;
    return *(static_cast<I*>(this));
  }
  I& operator-=(P rhs) {
    x -= rhs; y -= rhs; z -= rhs;
    return *(static_cast<I*>(this));
  }
  I& operator*=(P rhs) {
    x *= rhs; y *= rhs; z *= rhs;
    return *(static_cast<I*>(this));
  }
  I& operator/=(P rhs) {
    x /= rhs; y /= rhs; z /= rhs;
    return *(static_cast<I*>(this));
  }

  template <class V>
  BaseVec3<P,I> operator+(const BaseVec3<P,V>& v) const { 
    return BaseVec3<P,I>(x+v.x, y+v.y, z+v.z); 
  }
  BaseVec3<P,I> operator+(double c) const { 
    return BaseVec3<P,I>(x+c, y+c, z+c); 
  }
  template <class V>
  BaseVec3<P,I> operator-(const BaseVec<P,V>3& v) const { 
    return BaseVec3<P,I>(x-v.x, y-v.y, z-v.z); 
  }
  BaseVec3<P,I> operator-(double c) const { 
    return BaseVec3<P,I>(x-c, y-c, z-c); 
  }
  BaseVec3<P,I> operator-() const { 
    return BaseVec3<P,I>(-x,-y,-z); 
  }
  BaseVec3<P,I> operator*(double c) const { 
    return BaseVec3<P,I>(x*c, y*c, z*c); 
  }

  // transpose( transpose(Mat3) * Vec3 )

  template <class M>
  BaseVec3<P,I> operator*(const BaseMat3<P,M>& m) const;
  BaseVec3<P,I> operator/(double c) const { 
    return BaseVec3<P,I>(x/c, y/c, z/c); 
  }

  template <class V>
  double in(const BaseVec3<P,V>& v) { return ( x*v.x + y*v.y + z*v.z ); }

  template <class V>
  BaseVec3<P,I> out(const BaseVec3<P,V>& v) const { 
    return BaseVec3<P,I>(y*v.z-z*v.y,
			 z*v.x-x*v.z,
			 x*v.y-y*v.x); }
  Mat3 tns(const Vec3& v) const ;
  double norm() const { return (x*x+y*y+z*z); }
  double abs() const { return sqrt(x*x+y*y+z*z); }
  Vec3& normalize() { double l=abs();
                      x/=l; y/=l; z/=l;
                      return *this; }

 void print(const char *msg = "") const {
    if (*msg) cout << msg << ":" << endl;
    cout << setw(6) << setprecision(3) << x << 
      setw(6) << setprecision(3) << y << 
      setw(6) << setprecision(3) << z << endl;
  }
};

template <class T, class A>
ostream& operator<<(ostream& s, const BaseVec3<T,A>& a) {
  s << setw(6) << setprecision(3) << a.x;
  s << setw(6) << setprecision(3) << a.y;
  s << setw(6) << setprecision(3) << a.z << endl;
  return s;
}

//    Input Format
//    FLOAT
//    ( FLOAT )
//    ( FLOAT, FLOAT )
//    ( FLOAT, FLOAT, FLOAT)
template <class T, class A>
ostream& operator>>(ostream& s, const BaseVec3<T,A>& a) {
  double x=0.0, y=0.0, z=0.0;
  char c=0;
  s >> c;
  if ( c=='(' ) {
    s >> x >> c;
    if ( c==',' ) {
      s >> y >> c;
      if ( c==',' ) s >> z >> c;
      if ( c!=')' ) s.clear(ios::badbit);
    }
    if ( c!=')' ) s.clear(ios::badbit);
  } else {
    s.putback(c);
    s >> x;
  }
  if (s) {
    a.x = x; a.y = y; a.z = z;
  }
  return s;
}

#endif // TINYXYZ_H_
