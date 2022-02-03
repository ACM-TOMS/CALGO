/////////////////////////////////////////////////////////
//                                                     //
//    Cubpack++                                        //
//                                                     //
//        A Package For Automatic Cubature             //
//                                                     //
//        Authors : Ronald Cools                       //
//                  Dirk Laurie                        //
//                  Luc Pluym                          //
//                                                     //
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// File : point2D.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(file name changed)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Point_2D
// -------------------------
//
// BASECLASSES:
//   None
//
// PURPOSE:
//   Implements 2-dimensional points.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Point_2D()
//     ----------
//     default constructor.
//     the point must be initialized afterwards with an
//     assignment.
//
//     2) Point_2D(const Point_2D&)
//     ----------------------
//     copy constructor
//
//     3) Point_2D(const real, const real)
//     ------------------------------------
//     constructor. Two Cartesian coordinates are given.
//
//     4) Point_2D (int l)
//     ---------------------
//     constructs  a point of dimension l, l must be 2
//
//   SELECTORS:
//     1) real X()const
//     --------------------
//     returns the first component of the Point_2D.
//
//     2) real Y() const
//     --------------------
//     returns the second component of the Point_2D.
//
//     3) real R()const
//     --------------------
//     returns the first component of the Point_2D.
//     not to be confused with Length()
//
//     4) real Theta() const
//     ------------------------
//     returns the second component of the Point_2D.
//     not to be confused with Angle()
//
//     5) real operator[] (int i) const
//     --------------------------------
//     returns the i-th component. 0<=i<=1.
//
//     6)  real Length() const
//     ----------------------
//     returns the length of the vector from (0,0) to point.
//
//     7) real Angle() const
//     -----------------------
//     returns the angle (in [0, 2pi) ) between the
//     vector and positive X-axis
//
//     8) Point_2D Proj(const Point_2D& d) const
//     -----------------------------------------
//     returns the projection of the point on a
//     line with direction d
//
//     9) int Size() const
//     --------------------
//     returns the dimension of the point. (=2)
//
//   MODIFIERS:
//     1) real& X()
//     ------------
//     returns the first component of the Point_2D.
//
//     2) real& Y()
//     ------------
//     returns the second component of the Point_2D.
//
//     3) real& R()
//     -----------
//     returns the first component of the Point_2D.
//     not to be confused with Length()
//
//     4) real& Theta()
//     ----------------
//     returns the second component of the Point_2D.
//     not to be confused with Angle()
//
//     5) real& operator[] (int i)
//     ----------------------------
//     returns the i-th component. 0<=i<=1.
//
//   OPERATORS:
//     1)Point_2D& operator=(const Point_2D&);
//     --------------------------------
//     assignment operator
//
//     2)Point_2D operator+(const Point_2D&) const;
//     --------------------------------
//     point addition
//
//     3)Point_2D operator-(const Point_2D&) const;
//     --------------------------------
//     point subtraction
//
//     4)Point_2D operator-() const;
//     --------------------------------
//     multiplication by -1
//
//     5)Point_2D operator*( const real) const;
//     --------------------------------
//     multiplication with scalar
//
//     6)Point_2D operator/(const real) const;
//     --------------------------------
//     division by scalar
//
//     7)real operator*(const Point_2D& ) const;
//     --------------------------------
//     dot-product of two point
//
//     8)Point_2D& operator+=(const Point_2D&);
//     --------------------------------
//     addition and assignment
//
//     9)Point_2D& operator-=(const Point_2D&);
//     --------------------------------
//     subtraction and assignment
//
//     10)Point_2D& operator*=(real);
//     --------------------------------
//     scalar multiply and assignment
//
//     11)Point_2D& operator/=(real);
//     --------------------------------
//     division and assignment
//
//     12)Boolean operator==(const Point_2D&);
//     ---------------------------------------
//
//     13)Boolean operator!=(const Point_2D&);
//     ---------------------------------------
//
//   SPECIAL:
//     1) friend Point_2D operator*(real d,const Point_2D& p)
//     ----------------------------------------------------
//     product of the point with the scalar d
//
//     2) friend ostream& operator<<(ostream&,const Point_2D&p)
//     --------------------------------------------------------
//     output function
///////////////////////////////////////////////////////////
#ifndef _2DPOINT_H
#define _2DPOINT_H

/////////////////////////////////////////////////////////
// IMPORT FOR INTERFACE
/////////////////////////////////////////////////////////

#include <boolean.h>
#include <tools.h>

// introduction
class Point_2D;
  static Point_2D operator*(real,const Point_2D&);
  static ostream& operator<<(ostream&,const Point_2D& );

/////////////////////////////////////////////////////////
// EXPORT
/////////////////////////////////////////////////////////
class Point_2D
  {

  friend

  Point_2D operator*(real,const Point_2D&);

  friend

  ostream& operator<<(ostream&,const Point_2D& );

  public:

  Point_2D(const Point_2D&);
  Point_2D(const real, const real);
  Point_2D(int  n);
  Point_2D();
  ~Point_2D();
  Point_2D& operator=(const Point_2D&);
  Point_2D operator+(const Point_2D&) const;
  Point_2D operator-(const Point_2D&) const;
  Point_2D operator-() const;
  Point_2D operator*( const real) const;
  Point_2D operator/(const real) const;
  real operator*(const Point_2D& ) const;
  Point_2D& operator+=(const Point_2D&);
  Point_2D& operator-=(const Point_2D&);
  Point_2D& operator*=(real);
  Point_2D& operator/=(real);
  real X() const;
  real Y() const;
  real R() const;
  real Theta() const;
  real operator[](int i) const;
  real& X() ;
  real& Y() ;
  real& R() ;
  real& Theta() ;
  real& operator[](int i) ;
  real/* in[0,2pi] */ Angle() const;
  real Length() const;
  Boolean operator == (const Point_2D&) const;
  Boolean operator != (const Point_2D&) const;
  int Size() const;
  Point_2D Proj(const Point_2D& p)const;


  private:

  real x,y;
  };


/////////////////////////////////////////////
//INLINE DEFINITIONS
/////////////////////////////////////////////

#include <math.h>
#include <iostream.h>
#include <error.h>

//////////////////////////////////////////////////////

inline
Point_2D&
Point_2D::operator=(const Point_2D& v)
  {
  x = v.x;
  y = v.y;
  return *this;
  }
///////////////////////////////////////////////

inline
Point_2D
operator*(real d ,const Point_2D& v)
  {
  Point_2D r(v);
  r.x *= d;
  r.y *= d;
  return r;
  }
/////////////////////////////////////////////////
inline
real
Point_2D::X()
const
  {
  return x;
  }

////////////////////////////////////////////////////
inline
real
Point_2D::Y()
const
  {
  return y;
  }
////////////////////////////////////////////////////
inline
real
Point_2D::R()
const
  {
  return x;
  }
////////////////////////////////////////////////////
inline
real
Point_2D::Theta()
const
  {
  return y;
  }
/////////////////////////////////////////////////
inline
real&
Point_2D::X()
  {
  return x;
  }

////////////////////////////////////////////////////
inline
real&
Point_2D::Y()
  {
  return y;
  }
////////////////////////////////////////////////////
inline
real&
Point_2D::R()
  {
  return x;
  }
////////////////////////////////////////////////////
inline
real&
Point_2D::Theta()
  {
  return y;
  }
////////////////////////////////////////////////////
inline
real
Point_2D::Angle()
const
  {
  real Return;
  if (x==0)
    {
    Return = M_PI-M_PI/2*sign(y);
    };
  if (x<0)
    {
    Return = M_PI + atan(y/x);
    }
  else
    {
    real a=2*M_PI +atan(y/x);
    Return = (a>=2*M_PI) ? a-2*M_PI : a;
    };
  return Return ;
  }
/////////////////////////////////////////////////////
inline
ostream&
operator<< (ostream& os,const Point_2D& p)
  {
  os<< "Point_2D(" <<p.x<<","<<p.y<<")";
  return os;
  }
//////////////////////////////////////////////////////////
inline
Point_2D
Point_2D::operator+(const Point_2D& v)
const
  {
  Point_2D r(*this);
  r.x += v.x;
  r.y += v.y;
  return r;
  }

/////////////////////////////////////////////////////////
inline
Point_2D
Point_2D::operator-(const Point_2D& v)
const
  {
  Point_2D r(*this);
  r.x -= v.x;
  r.y -= v.y;
  return r;
  }

////////////////////////////////////////////////////
inline
Point_2D
Point_2D::operator-()
const
  {
  Point_2D r(-x,-y);
  return r;
  }

////////////////////////////////////////////////////
inline
Point_2D
Point_2D::operator*(const real d)
const
  {
  Point_2D r(*this);
  r.x *= d;
  r.y *= d;
  return r;
  }

////////////////////////////////////////////////////
inline
Point_2D
Point_2D::operator/(const real d)
const
  {
  Error(d==0,"Point_2D:division by zero error");
  Point_2D r(*this);
  r.x/= d;
  r.y/= d;
  return r;
  }

////////////////////////////////////////////////////
inline
real
Point_2D::Length()
const
  {
  return sqrt(x*x+y*y);
  //if (y>x)
    //{
    //return y*sqrt((x/y)*(x/y)+1);
    //}
  //else
    //{
    //return x*sqrt((y/x)*(y/x)+1);
    //};
  }
////////////////////////////////////////////////////
inline
real
Point_2D::operator *(const Point_2D& v)
const
  {
  real inprod =x*v.x+y*v.y;
  return inprod;
  }
////////////////////////////////////////////////////
inline
Point_2D&
Point_2D::operator +=(const Point_2D& v)
  {
  x += v.x;
  y += v.y;
  return (*this);
  }
///////////////////////////////////////////////////
inline
Point_2D&
Point_2D::operator -=(const Point_2D& v)
  {
  x -= v.x;
  y -= v.y;
  return (*this);
  }
///////////////////////////////////////////////////
inline
Point_2D&
Point_2D::operator *=(real d)
  {
  x *=d;
  y *=d;
  return (*this);
  }
///////////////////////////////////////////////////
inline
Point_2D&
Point_2D::operator /=(real d)
  {
  Error(d==0,"Point_2D:division by zero error");
  x /= d;
  y /= d;
  return (*this);
  }
///////////////////////////////////////////////////
inline
real
Point_2D::operator[](int i)
const
  {
  return i ? y : x;
  }
///////////////////////////////////////////////
inline
real&
Point_2D::operator[](int i)
  {
  return i ? y : x;
  }
///////////////////////////////////////////////
inline
Boolean
Point_2D::operator==(const Point_2D& v)
const
  {
    int b =  (x == v.x)&&(y == v.y);
    return (Boolean) b;
  }
////////////////////////////////////////////////////
inline
Boolean
Point_2D::operator!=(const Point_2D& v)
const
  {
    int b =  (x != v.x)||(y != v.y);
    return (Boolean) b;
  }
////////////////////////////////////////////////////
inline
int
Point_2D::Size() const
  {
  return 2;
  }
//////////////////////////////////////////////////
inline
Point_2D
Point_2D::Proj(const Point_2D& v)
const
  {
  Point_2D r(v);
  r *= (*this) *v;
  r /= (*this) *  (*this);
  return r;
  }
//////////////////////////////////////////////////
inline
Point_2D::Point_2D(int n)
  {
  Error(n  != 2, "a Point_2D should be two-dimensional");
  }
/////////////////////////////////////////////////////
#endif
