//-*- C++ -*-
/*
   2D Vector and Matrix
*/
#ifndef MAT2_H_
#define MAT2_H_

#include <math.h>
#include <iostream.h>

#include "cmplxtype.h"

class Mat2;

class Vec2 {
public:
  double x,y;

  Vec2(double cx=0,double cy=0) : x(cx),y(cy) {}
  Vec2(const Vec2& v) { x = v.x; y = v.y; }
  const Vec2& operator=(const Vec2& v) { x = v.x; y = v.y; return *this; }
  Vec2 operator+(const Vec2& v) const { return Vec2(x+v.x, y+v.y); }
  Vec2 operator+(double c) const { return Vec2(x+c, y+c); }
  Vec2 operator-(const Vec2& v) const { return Vec2(x-v.x, y-v.y); }
  Vec2 operator-(double c) const { return Vec2(x-c, y-c); }
  Vec2 operator-() const { return Vec2(-x,-y); }
  Vec2 operator*(double c) const { return Vec2(x*c, y*c); }
  // transpose( transpose(Mat2) * Vec2 )
  Vec2 operator*(const Mat2& m) const;
  Vec2 operator/(double c) const { return Vec2(x/c, y/c); }
  Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
  Vec2& operator+=(double c) { x += c; y += c; return *this; }
  Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
  Vec2& operator*=(double c) { x *= c; y *= c; return *this; }
  Vec2& operator/=(double c) { x /= c; y /= c; return *this; }
  double in(const Vec2& v) const { return ( x*v.x+y*v.y ); }
  Mat2 tns(const Vec2& v) const ;
  double norm() const { return (x*x+y*y); }
  double abs() const { return sqrt(x*x+y*y); }
  Vec2& normalize() { double l=abs();
                      x/=l; y/=l;
                      return *this; }

  void print(const char *msg="") const;
};

ostream& operator<<(ostream& s, const Vec2& v);
istream& operator>>(istream& s, Vec2& v);

class Mat2 {
public:
  double xx,xy;
  double yx,yy;
  Mat2(double cxx=0, double cxy=0,
       double cyx=0, double cyy=0)
    : xx(cxx),xy(cxy),
      yx(cyx),yy(cyy) {}
  Mat2(const Mat2& m) {
    xx=m.xx; xy=m.xy;
    yx=m.yx; yy=m.yy;
  }
  const Mat2& operator=(const Mat2& m) {
    xx=m.xx; xy=m.xy;
    yx=m.yx; yy=m.yy;
    return *this;
  }
  Mat2 operator+(const Mat2& m) const { 
    return Mat2(xx+m.xx, xy+m.xy,
		yx+m.yx, yy+m.yy); 
  }
  Mat2 operator+(double c) const { return Mat2(xx+c, xy  ,
					       yx  , yy+c); }
  Mat2 operator-(const Mat2& m) const { 
    return Mat2(xx-m.xx, xy-m.xy,
		yx-m.yx, yy-m.yy); 
  }
  Mat2 operator-(double c) const { return Mat2(xx-c, xy  ,
					       yx  , yy-c); }
  Mat2 operator-() const { return Mat2(-xx,-xy,
				       -yx,-yy); }
  Mat2 operator*(const Mat2& m) const { 
    return Mat2(xx*m.xx + xy*m.yx,
		xx*m.xy + xy*m.yy,
		yx*m.xx + yy*m.yx,
		yx*m.xy + yy*m.yy); 
  }
  Vec2 operator*(const Vec2& v) const { 
    return Vec2(xx*v.x + xy*v.y,
		yx*v.x + yy*v.y); 
  }
  Mat2 operator*(double c) const { return Mat2(xx*c, xy*c,
					       yx*c, yy*c); }
  Mat2 operator/(double c) const { return Mat2(xx/c, xy/c,
					       yx/c, yy/c); }
  double abs() const { return sqrt(xx*xx+xy*xy+
				 yx*yx+yy*yy); }
  void print(const char *msg="") const;
};

inline Vec2 Vec2::operator*(const Mat2& m) const { 
  return Vec2(x*m.xx+y*m.yx,
	      x*m.xy+y*m.yy);
}

inline Mat2 Vec2::tns(const Vec2& v) const { 
  return Mat2(x*v.x, x*v.y,
	      y*v.x, y*v.y); 
}

//
// Complex Vector and Matrix
//

class CMat2;

class CVec2 {
public:
  Complex x,y;

  CVec2(Complex cx=0,Complex cy=0) : x(cx),y(cy) {}
  CVec2(const CVec2& v) { x = v.x; y = v.y; }
  CVec2& operator=(const CVec2& v) { x = v.x; y = v.y; return *this; }
  CVec2 operator+(const CVec2& v) const { return CVec2(x+v.x, y+v.y); }
  CVec2 operator+(Complex c) const { return CVec2(x+c, y+c); }
  CVec2 operator-(const CVec2& v) const { return CVec2(x-v.x, y-v.y); }
  CVec2 operator-(Complex c) const { return CVec2(x-c, y-c); }
  CVec2 operator-() const { return CVec2(-x,-y); }
  CVec2 operator*(Complex c) const { return CVec2(x*c, y*c); }
  // transpose( transpose(CMat2) * CVec2 )
  CVec2 operator*(const CMat2& m) const;
  CVec2 operator/(Complex c) const { return CVec2(x/c, y/c); }
  CVec2& operator+=(const CVec2& v) { x += v.x; y += v.y; return *this; }
  CVec2& operator+=(Complex c) { x += c; y += c; return *this; }
  CVec2& operator-=(const CVec2& v) { x -= v.x; y -= v.y; return *this; }
  CVec2& operator*=(Complex c) { x *= c; y *= c; return *this; }
  CVec2& operator/=(Complex c) { x /= c; y /= c; return *this; }
  Complex in(const CVec2& v) const { return ( x*v.x+y*v.y ); }
  CMat2 tns(const CVec2& v) const ;
  Complex abs() const { return sqrt(x*x+y*y); }
  CVec2& normalize() { Complex l=abs();
                      x/=l; y/=l;
                      return *this; }

  void print(const char *msg="") const;
};

ostream& operator<<(ostream& s, const CVec2& v);
istream& operator>>(istream& s, CVec2& v);

class CMat2 {
public:
  Complex xx,xy;
  Complex yx,yy;
  CMat2(Complex cxx=0, Complex cxy=0,
       Complex cyx=0, Complex cyy=0)
    : xx(cxx),xy(cxy),
      yx(cyx),yy(cyy) {}
  CMat2(const CMat2& m) {
    xx=m.xx; xy=m.xy;
    yx=m.yx; yy=m.yy;
  }
  CMat2& operator=(const CMat2& m) {
    xx=m.xx; xy=m.xy;
    yx=m.yx; yy=m.yy;
    return *this;
  }
  CMat2 operator+(const CMat2& m) const { 
    return CMat2(xx+m.xx, xy+m.xy,
		yx+m.yx, yy+m.yy); 
  }
  CMat2 operator+(Complex c) const { return CMat2(xx+c, xy  ,
					       yx  , yy+c); }
  CMat2 operator-(const CMat2& m) const { 
    return CMat2(xx-m.xx, xy-m.xy,
		yx-m.yx, yy-m.yy); 
  }
  CMat2 operator-(Complex c) const { return CMat2(xx-c, xy  ,
					       yx  , yy-c); }
  CMat2 operator-() const { return CMat2(-xx,-xy,
				       -yx,-yy); }
  CMat2 operator*(const CMat2& m) const { 
    return CMat2(xx*m.xx + xy*m.yx,
		xx*m.xy + xy*m.yy,
		yx*m.xx + yy*m.yx,
		yx*m.xy + yy*m.yy); 
  }
  CVec2 operator*(const CVec2& v) const { 
    return CVec2(xx*v.x + xy*v.y,
		yx*v.x + yy*v.y); 
  }
  CMat2 operator*(Complex c) const { return CMat2(xx*c, xy*c,
					       yx*c, yy*c); }
  CMat2 operator/(Complex c) const { return CMat2(xx/c, xy/c,
					       yx/c, yy/c); }
  Complex abs() const { return sqrt(xx*xx+xy*xy+
				 yx*yx+yy*yy); }
  void print(const char *msg="") const;
};

inline CVec2 CVec2::operator*(const CMat2& m) const { 
  return CVec2(x*m.xx+y*m.yx,
	      x*m.xy+y*m.yy);
}

inline CMat2 CVec2::tns(const CVec2& v) const { 
  return CMat2(x*v.x, x*v.y,
	      y*v.x, y*v.y); 
}

#endif // MAT2_H_
