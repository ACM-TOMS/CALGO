//-*- C++ -*-
/*
   3D Vector and Matrix
*/
#ifndef MAT3_H_
#define MAT3_H_

#include <math.h>
#include <iostream.h>

class Mat3;
class DiagMat3;

class Vec3 {
public:
  double x,y,z;

  Vec3(double cx=0,double cy=0,double cz=0) : x(cx),y(cy),z(cz) {}
  Vec3(const Vec3& v) { x = v.x; y = v.y; z = v.z; }
  void assignFrom(const Vec3& v) {  x = v.x; y = v.y; z = v.z; }
  Vec3& operator=(const Vec3& v) { 
    x = v.x; y = v.y; z = v.z; return *this; 
  }
  Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
  Vec3 operator+(double c) const { return Vec3(x+c, y+c, z+c); }
  Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
  Vec3 operator-(double c) const { return Vec3(x-c, y-c, z-c); }
  Vec3 operator-() const { return Vec3(-x,-y,-z); }
  Vec3 operator*(double c) const { return Vec3(x*c, y*c, z*c); }
  // transpose( transpose(Mat3) * Vec3 )
  // as a row vector
  Vec3 operator*(const Mat3& m) const;
  Vec3 operator*(const DiagMat3& m) const;
  Vec3 operator/(double c) const { return Vec3(x/c, y/c, z/c); }
  Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
  Vec3& operator+=(double c) { x += c; y += c; z += c; return *this; }
  Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
  Vec3& operator*=(double c) { x *= c; y *= c; z *= c; return *this; }
  Vec3& operator/=(double c) { x /= c; y /= c; z /= c; return *this; }
  double in(const Vec3& v) const { return ( x*v.x+y*v.y+z*v.z ); }
  Vec3 out(const Vec3& v) const { return Vec3(y*v.z-z*v.y,
					      z*v.x-x*v.z,
					      x*v.y-y*v.x); }
  Mat3 tns(const Vec3& v) const ;
  double norm() const { return (x*x+y*y+z*z); }
  double abs() const { return sqrt(x*x+y*y+z*z); }
  Vec3& normalize() { double l=abs();
                      x/=l; y/=l; z/=l;
                      return *this; }

  void print(const char *msg="") const;
};

ostream& operator<<(ostream& s, const Vec3& v);
istream& operator>>(istream& s, Vec3& v);

class Mat3 {
public:
  double xx,xy,xz;
  double yx,yy,yz;
  double zx,zy,zz;
  Mat3(double cxx=0, double cxy=0, double cxz=0,
       double cyx=0, double cyy=0, double cyz=0,
       double czx=0, double czy=0, double czz=0)
    : xx(cxx),xy(cxy),xz(cxz),
      yx(cyx),yy(cyy),yz(cyz),
      zx(czx),zy(czy),zz(czz) {}
  Mat3(const Mat3& m) {
    xx=m.xx; xy=m.xy; xz=m.xz;
    yx=m.yx; yy=m.yy; yz=m.yz;
    zx=m.zx; zy=m.zy; zz=m.zz;
  }
  void assignFrom(const Mat3& m) {
    xx=m.xx; xy=m.xy; xz=m.xz;
    yx=m.yx; yy=m.yy; yz=m.yz;
    zx=m.zx; zy=m.zy; zz=m.zz;
  }
  Mat3& operator=(const Mat3& m) {
    xx=m.xx; xy=m.xy; xz=m.xz;
    yx=m.yx; yy=m.yy; yz=m.yz;
    zx=m.zx; zy=m.zy; zz=m.zz;
    return *this;
  }
  Mat3 operator+(const Mat3& m) const { 
    return Mat3(xx+m.xx, xy+m.xy, xz+m.xz,
		yx+m.yx, yy+m.yy, yz+m.yz,
		zx+m.zx, zy+m.zy, zz+m.zz); 
  }
  Mat3 operator+(double c) const { return Mat3(xx+c, xy  , xz  ,
					       yx  , yy+c, yz  ,
					       zx  , zy  , zz+c); }
  Mat3 operator-(const Mat3& m) const { 
    return Mat3(xx-m.xx, xy-m.xy, xz-m.xz,
		yx-m.yx, yy-m.yy, yz-m.yz,
		zx-m.zx, zy-m.zy, zz-m.zz); 
  }
  Mat3 operator-(double c) const { return Mat3(xx-c, xy  , xz  ,
					       yx  , yy-c, yz  ,
					       zx  , zy  , zz-c); }
  Mat3 operator-() const { return Mat3(-xx,-xy,-xz,
				       -yx,-yy,-yz,
				       -zx,-zy,-zz); }
  Mat3 operator*(const Mat3& m) const { 
    return Mat3(xx*m.xx + xy*m.yx + xz*m.zx,
		xx*m.xy + xy*m.yy + xz*m.zy,
		xx*m.xz + xy*m.yz + xz*m.zz,
		yx*m.xx + yy*m.yx + yz*m.zx,
		yx*m.xy + yy*m.yy + yz*m.zy,
		yx*m.xz + yy*m.yz + yz*m.zz,
		zx*m.xx + zy*m.yx + zz*m.zx,
		zx*m.xy + zy*m.yy + zz*m.zy,
		zx*m.xz + zy*m.yz + zz*m.zz); 
  }
  Vec3 operator*(const Vec3& v) const { 
    return Vec3(xx*v.x + xy*v.y + xz*v.z,
		yx*v.x + yy*v.y + yz*v.z,
		zx*v.x + zy*v.y + zz*v.z); 
  }
  Mat3 operator*(double c) const { return Mat3(xx*c, xy*c, xz*c,
					       yx*c, yy*c, yz*c,
					       zx*c, zy*c, zz*c); }
  Mat3 operator/(double c) const { return Mat3(xx/c, xy/c, xz/c,
					       yx/c, yy/c, yz/c,
					       zx/c, zy/c, zz/c); }
  double abs() const { return sqrt(xx*xx+xy*xy+xz*xz+
				 yx*yx+yy*yy+yz*yz+
				 zx*zx+zy*zy+zz*zz); }
  void print(const char *msg="") const;
};

// as a row vector
inline Vec3 Vec3::operator*(const Mat3& m) const { 
  return Vec3(x*m.xx+y*m.yx+z*m.zx,
	      x*m.xy+y*m.yy+z*m.zy,
	      x*m.xz+y*m.yz+z*m.zz);
}

inline Mat3 Vec3::tns(const Vec3& v) const { 
  return Mat3(x*v.x, x*v.y, x*v.z,
	      y*v.x, y*v.y, y*v.z,
	      z*v.x, z*v.y, z*v.z); 
}

class DiagMat3 {
public:
  double xx,yy,zz;
  DiagMat3(double cxx=0, double cyy=0, double czz=0)
    : xx(cxx),yy(cyy),zz(czz) {}
  DiagMat3(const DiagMat3& m) {
    xx=m.xx; yy=m.yy; zz=m.zz;
  }
  void assignFrom(const DiagMat3& m) {
    xx=m.xx; yy=m.yy; zz=m.zz;
  }
  DiagMat3& operator=(const DiagMat3& m) {
    xx=m.xx; yy=m.yy; zz=m.zz;
    return *this;
  }
  DiagMat3 operator+(const DiagMat3& m) const { 
    return DiagMat3(xx+m.xx, yy+m.yy, zz+m.zz); 
  }
  DiagMat3 operator+(double c) const { return DiagMat3(xx+c, yy+c, zz+c); }
  DiagMat3 operator-(const DiagMat3& m) const { 
    return DiagMat3(xx-m.xx, yy-m.yy, zz-m.zz); 
  }
  DiagMat3 operator-(double c) const { return DiagMat3(xx-c, yy-c, zz-c); }
  DiagMat3 operator-() const { return DiagMat3(-xx,-yy,-zz); }
  DiagMat3 operator*(const DiagMat3& m) const { 
    return DiagMat3(xx*m.xx, yy*m.yy, zz*m.zz); 
  }
  Vec3 operator*(const Vec3& v) const { return Vec3(xx*v.x, yy*v.y, zz*v.z); }
  DiagMat3 operator*(double c) const { return DiagMat3(xx*c, yy*c, zz*c); }
  DiagMat3 operator/(double c) const { return DiagMat3(xx/c, yy/c, zz/c); }
  double abs() const { return sqrt(xx*xx+yy*yy+zz*zz); }
  void print(const char *msg="") const;
};

inline Vec3 Vec3::operator*(const DiagMat3& m) const { 
  return Vec3(x*m.xx, y*m.yy, z*m.zz);
}

#endif // MAT3_H_
