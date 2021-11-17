//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifndef COLOR_H
#define COLOR_H

#include "util.h"

//namespace genial
//{

using namespace std;

template<class T1,class T2,class T3> class ColorSpace;
template<class T1,class T2,class T3> class RGBColor;
template<class T1,class T2,class T3> class BGRColor;
template<class T1,class T2,class T3> class YUVColor;
template<class T1,class T2,class T3> class YCbCrColor;
template<class T1,class T2,class T3,class T4> class RGBAColor;
template<class T1,class T2,class T3,class T4> class YUVAColor;

template<class T, class Q> struct quantizerIndexArray;
template<class T, class Q> struct quantizerValueArray;
template<class T, class Q> struct quantizerQuantizeArray;

template<class T1,class T2,class T3,class V> struct promotion2_traits<ColorSpace<T1,T2,T3>, V > { typedef ColorSpace<PROMOTE2(T1,V),PROMOTE2(T2,V),PROMOTE2(T2,V) > value_type; };
template<class V,class T1,class T2,class T3> struct promotion2_traits<V, ColorSpace<T1,T2,T3> > { typedef ColorSpace<PROMOTE2(V,T1),PROMOTE2(V,T2),PROMOTE2(V,T3) > value_type; };
template<class T1,class T2,class T3,class U1,class U2,class U3> struct promotion2_traits<ColorSpace<T1,T2,T3>, ColorSpace<U1,U2,U3> > { typedef ColorSpace<PROMOTE2(T1,U1),PROMOTE2(T2,U2),PROMOTE2(T2,U3)> value_type; };

template<class T1,class T2,class T3,class V> struct promotion2_traits<RGBColor<T1,T2,T3>, V > { typedef RGBColor<PROMOTE2(T1,V),PROMOTE2(T2,V),PROMOTE2(T2,V) > value_type; };
template<class V,class T1,class T2,class T3> struct promotion2_traits<V, RGBColor<T1,T2,T3> > { typedef RGBColor<PROMOTE2(V,T1),PROMOTE2(V,T2),PROMOTE2(V,T3) > value_type; };
template<class T1,class T2,class T3,class U1,class U2,class U3> struct promotion2_traits<RGBColor<T1,T2,T3>, RGBColor<U1,U2,U3> > { typedef RGBColor<PROMOTE2(T1,U1),PROMOTE2(T2,U2),PROMOTE2(T2,U3) > value_type; };

template<class T1,class T2,class T3,class V> struct promotion2_traits<YUVColor<T1,T2,T3>, V > { typedef YUVColor<PROMOTE2(T1,V),PROMOTE2(T2,V),PROMOTE2(T2,V) > value_type; };
template<class V,class T1,class T2,class T3> struct promotion2_traits<V, YUVColor<T1,T2,T3> > { typedef YUVColor<PROMOTE2(V,T1),PROMOTE2(V,T2),PROMOTE2(V,T3) > value_type; };
template<class T1,class T2,class T3,class U1,class U2,class U3> struct promotion2_traits<YUVColor<T1,T2,T3>, YUVColor<U1,U2,U3> > { typedef YUVColor<PROMOTE2(T1,U1),PROMOTE2(T2,U2),PROMOTE2(T2,U3) > value_type; };

template<class T1,class T2,class T3,class V> struct promotion2_traits<YCbCrColor<T1,T2,T3>, V > { typedef YCbCrColor<PROMOTE2(T1,V),PROMOTE2(T2,V),PROMOTE2(T2,V) > value_type; };
template<class V,class T1,class T2,class T3> struct promotion2_traits<V, YCbCrColor<T1,T2,T3> > { typedef YCbCrColor<PROMOTE2(V,T1),PROMOTE2(V,T2),PROMOTE2(V,T3) > value_type; };
template<class T1,class T2,class T3,class U1,class U2,class U3> struct promotion2_traits<YCbCrColor<T1,T2,T3>, YCbCrColor<U1,U2,U3> > { typedef YCbCrColor<PROMOTE2(T1,U1),PROMOTE2(T2,U2),PROMOTE2(T2,U3) > value_type; };

//{unsecret}
//{group:Images Color Spaces Classes}
//Summary: Base class for color spaces with 3 components
//Arguments:
//  T1 - First component type
//  T2 - Second component type
//  T3 - Third component type
template<class T1=unsigned char, class T2=T1, class T3=T2>
class ColorSpace
{
  public:
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

  public:
    first_type  x;
    second_type y;
    third_type  z;
    
  public:
    inline ColorSpace() : x(), y(), z() {}
    inline ColorSpace(first_type a, second_type b, third_type c) : x(a), y(b), z(c) {}
    template<class U1,class U2,class U3> inline ColorSpace(U1 a, U2 b, U3 c) : x(bound<T1>(a)), y(bound<T2>(b)), z(bound<T3>(c)) {}

    inline ColorSpace(const ColorSpace &r) : x(r.x), y(r.y), z(r.z) {}
    template<class U1,class U2,class U3> inline ColorSpace(const ColorSpace<U1,U2,U3> &r) : x(bound<T1>(r.x)), y(bound<T2>(r.y)), z(bound<T3>(r.z)) {}

    inline ColorSpace &operator=(const ColorSpace &r) { x=r.x; y=r.y; z=r.z; return *this; }

    inline bool operator==(const ColorSpace &r) const { return x==r.x && y==r.y && z==r.z; }
    inline bool operator!=(const ColorSpace &r) const { return x!=r.x || y!=r.y || z!=r.z; }
};

//{unsecret}
//{group:Images Color Spaces External Functions}
//Summary: First color component
template<class T1,class T2,class T3> inline T1       &first (      ColorSpace<T1,T2,T3> &x) { return x.x; }
//{unsecret}
template<class T1,class T2,class T3> inline const T1 &first (const ColorSpace<T1,T2,T3> &x) { return x.x; }
//{unsecret}
//{group:Images Color Spaces External Functions}
//Summary: Second color component
template<class T1,class T2,class T3> inline T2       &second(      ColorSpace<T1,T2,T3> &x) { return x.y; }
//{unsecret}
template<class T1,class T2,class T3> inline const T2 &second(const ColorSpace<T1,T2,T3> &x) { return x.y; }
//{unsecret}
//{group:Images Color Spaces External Functions}
//Summary: Third color component
template<class T1,class T2,class T3> inline T3       &third (      ColorSpace<T1,T2,T3> &x) { return x.z; }
//{unsecret}
template<class T1,class T2,class T3> inline const T3 &third (const ColorSpace<T1,T2,T3> &x) { return x.z; }


template<class E,class Tr,class T1,class T2,class T3> basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const ColorSpace<           T1,           T2,           T3> &X) { return os << "(" <<      X.x << "," <<      X.y << "," <<      X.z << ")"; }
template<class E,class Tr>                            basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const ColorSpace<unsigned char,unsigned char,unsigned char> &X) { return os << "(" << (int)X.x << "," << (int)X.y << "," << (int)X.z << ")"; }
template<class E,class Tr>                            basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const ColorSpace<unsigned char,         char,         char> &X) { return os << "(" << (int)X.x << "," << (int)X.y << "," << (int)X.z << ")"; }
template<class E,class Tr,class T1,class T2,class T3> basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       ColorSpace<           T1,           T2,           T3> &X) { E c; is>>c && c=='(' && is>>first(X)>>c && c==',' && is>>second(X)>>c && c==',' && is>>third(X)>>c && c==')'; return is; }
template<class E,class Tr>                            basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       ColorSpace<unsigned char,unsigned char,unsigned char> &X) { E c; int i0=0,i1=0,i2=0; is>>c && c=='(' && is>>i0>>c && c==',' && is>>i1>>c && c==',' && is>>i2>>c && c==')'; first(X)=i0; second(X)=i1; third(X)=i2; return is; }
template<class E,class Tr>                            basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       ColorSpace<unsigned char,         char,         char> &X) { E c; int i0=0,i1=0,i2=0; is>>c && c=='(' && is>>i0>>c && c==',' && is>>i1>>c && c==',' && is>>i2>>c && c==')'; first(X)=i0; second(X)=i1; third(X)=i2; return is; }

template<class T1,class T2,class T3> inline ColorSpace<T1,T2,T3> operator+(const ColorSpace<T1,T2,T3> &X, const ColorSpace<T1,T2,T3> &Y) { return ColorSpace<T1,T2,T3>(X.x+Y.x, X.y+Y.y, X.z+Y.z); }
template<class T1,class T2,class T3> inline ColorSpace<T1,T2,T3> operator*(const ColorSpace<T1,T2,T3> &X, const ColorSpace<T1,T2,T3> &Y) { return ColorSpace<T1,T2,T3>(X.x*Y.x, X.y*Y.y, X.z*Y.z); }

template<class T1,class T2,class T3> inline ColorSpace<PROMOTE2(T1,T1),PROMOTE2(T1,T2),PROMOTE2(T1,T3)> operator*(const T1 &v, const ColorSpace<T1,T2,T3> &X) { return ColorSpace<PROMOTE2(T1,T1),PROMOTE2(T1,T2),PROMOTE2(T1,T3)>(v*X.x, v*X.y, v*X.z); }
template<class T1,class T2,class T3> inline ColorSpace<PROMOTE2(T1,T1),PROMOTE2(T2,T1),PROMOTE2(T3,T1)> operator*(const ColorSpace<T1,T2,T3> &X, const T1 &v) { return ColorSpace<PROMOTE2(T1,T1),PROMOTE2(T2,T1),PROMOTE2(T3,T1)>(X.x*v, X.y*v, X.z*v); }

template<class T1,class T2,class T3,class V> inline ColorSpace<PROMOTE2(T1,V),PROMOTE2(T2,V),PROMOTE2(T3,V)> operator/(const ColorSpace<T1,T2,T3> &X, const V &v) { return ColorSpace<PROMOTE2(T1,V),PROMOTE2(T2,V),PROMOTE2(T3,V)>(X.x/v, X.y/v, X.z/v); }

template<class T1, class T2, class T3> inline ColorSpace<T1,T2,T3> round(const ColorSpace<T1,T2,T3> &X) { return ColorSpace<T1,T2,T3>(round(X.x), round(X.y), round(X.z)); }

//{unsecret}
//{group:Images Color Spaces Classes}
//Summary: Base class for color spaces with 3 components and an alpha value
//Arguments:
//  T1 - First component type
//  T2 - Second component type
//  T3 - Third component type
//  T4 - Alpha type
template<class T1=unsigned char, class T2=T1, class T3=T2,class T4=unsigned char>
class AlphaColorSpace : public ColorSpace<T1,T2,T3>
{
  public:
    typedef AlphaColorSpace self;
    typedef ColorSpace<T1,T2,T3> base;
    
    typedef T4 alpha_type;

  public:
    alpha_type alpha;
    
  public:
    inline AlphaColorSpace() : base(), alpha()  {}
    template<class U1,class U2,class U3> inline AlphaColorSpace(U1 a, U2 b, U3 c) : base(a,b,c), alpha(0) {}
    template<class U1,class U2,class U3> inline AlphaColorSpace(U1 a, U2 b, U3 c, alpha_type d) : base(a,b,c), alpha(d) {}

    inline AlphaColorSpace(const AlphaColorSpace &r) : base(r), alpha(r.alpha) {}
    template<class U1,class U2,class U3,class U4> inline AlphaColorSpace(const AlphaColorSpace<U1,U2,U3,U4> &r) : base(r), alpha(r.alpha) {}
};

template<class E,class Tr,class T1,class T2,class T3,class T4> basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const AlphaColorSpace<           T1,           T2,           T3,           T4> &X) { return os << "(" <<      X.x << "," <<      X.y << "," <<      X.z << "," <<      X.alpha << ")"; }
template<class E,class Tr>                                     basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const AlphaColorSpace<unsigned char,unsigned char,unsigned char,unsigned char> &X) { return os << "(" << (int)X.x << "," << (int)X.y << "," << (int)X.z << "," << (int)X.alpha << ")"; }
template<class E,class Tr>                                     basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const AlphaColorSpace<unsigned char,         char,         char,unsigned char> &X) { return os << "(" << (int)X.x << "," << (int)X.y << "," << (int)X.z << "," << (int)X.alpha << ")"; }
template<class E,class Tr,class T1,class T2,class T3,class T4> basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       AlphaColorSpace<           T1,           T2,           T3,           T4> &X) { E c; is>>c && c=='(' && is>>first(X)>>c && c==',' && is>>second(X)>>c && c==',' && is>>third(X)>>c && c==',' && is>>X.alpha>>c && c==')'; return is; }
template<class E,class Tr>                                     basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       AlphaColorSpace<unsigned char,unsigned char,unsigned char,unsigned char> &X) { E c; int i0=0,i1=0,i2=0,i3=0; is>>c && c=='(' && is>>i0>>c && c==',' && is>>i1>>c && c==',' && is>>i2>>c && c==',' && is>>i3>>c && c==')'; first(X)=i0; second(X)=i1; third(X)=i2; X.alpha=i3; return is; }
template<class E,class Tr>                                     basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       AlphaColorSpace<unsigned char,         char,         char,unsigned char> &X) { E c; int i0=0,i1=0,i2=0,i3=0; is>>c && c=='(' && is>>i0>>c && c==',' && is>>i1>>c && c==',' && is>>i2>>c && c==',' && is>>i3>>c && c==')'; first(X)=i0; second(X)=i1; third(X)=i2; X.alpha=i3; return is; }


//{unsecret}
//{group:Images Color Spaces Classes}
//Summary: RGB color space
//Arguments:
//  T1 - First component type
//  T2 - Second component type
//  T3 - Third component type
template<class T1=unsigned char, class T2=T1, class T3=T2>
class RGBColor : public ColorSpace<T1,T2,T3>
{
  public:
    typedef ColorSpace<T1,T2,T3> base;

    using base::x;
    using base::y;
    using base::z;
    
  public:
    inline RGBColor() : base() {}
    inline RGBColor(T1 R, T2 G, T3 B) : base(R,G,B) {}

    template<class U1,class U2,class U3> inline RGBColor(const ColorSpace<U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline RGBColor(const RGBColor  <U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline RGBColor(const BGRColor  <U1,U2,U3> &r) : base(r.z,r.y,r.x) {}
    template<class U1,class U2,class U3> inline RGBColor(const YUVColor  <U1,U2,U3> &r) : base(r.x+1.140 *r.z, r.x-0.396 *r.y-0.581 *r.z, r.x+2.029 *r.y) { }
    template<class U1,class U2,class U3> inline RGBColor(const YCbCrColor<U1,U2,U3> &r) : base(r.x+1.4022*r.z, r.x-0.3456*r.y-0.7145*r.z, r.x+1.7710*r.y) { }

    template<class U1,class U2,class U3> inline RGBColor &operator=(const YUVColor<U1,U2,U3> &r) { *this=RGBColor(r); return *this; }
    
    inline operator unsigned char() { return 0.299*x+0.587*y+0.114*z; }
};


template<class T1=unsigned char, class T2=T1, class T3=T2>
class BGRColor : public ColorSpace<T1,T2,T3>
{
  public:
    typedef ColorSpace<T1,T2,T3> base;

    using base::x;
    using base::y;
    using base::z;

  public:
    inline BGRColor() : base() {}
    inline BGRColor(T1 B, T2 G, T3 R) : base(B,G,R) {}

    template<class U1,class U2,class U3> inline BGRColor(const ColorSpace<U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline BGRColor(const BGRColor  <U1,U2,U3> &r) : base(r) {}

    template<class U1,class U2,class U3> inline BGRColor(const RGBColor  <U1,U2,U3> &r) : base(r.z,r.y,r.x) { }

    inline operator unsigned char() { return 0.114*x+0.587*y+0.299*z; }
};

//{unsecret}
//{group:Images Color Spaces Classes}
//Summary: YUV color space
//Arguments:
//  T1 - First component type
//  T2 - Second component type
//  T3 - Third component type
template<class T1=unsigned char, class T2=T1, class T3=T2>
class YUVColor : public ColorSpace<T1,T2,T3>
{
  public:
    typedef ColorSpace<T1,T2,T3> base;

  public:
    inline YUVColor() : base() {}
    inline YUVColor(T1 Y, T2 U, T3 V) : base(Y,U,V) {}

    template<class U1,class U2,class U3> inline YUVColor(const ColorSpace<U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline YUVColor(const YUVColor  <U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline YUVColor(const RGBColor  <U1,U2,U3> &r) : base(0.299*r.x+0.587*r.y+0.114*r.z, -0.147*r.x-0.289*r.y+0.436*r.z, 0.615*r.x-0.515*r.y-0.100*r.z) { }

    template<class U1,class U2,class U3> inline YUVColor &operator=(const RGBColor<U1,U2,U3> &r) { *this=YUVColor(r); return *this; }
};

//{unsecret}
//{group:Images Color Spaces Classes}
//Summary: YCbCr color space
//Arguments:
//  T1 - First component type
//  T2 - Second component type
//  T3 - Third component type
template<class T1=unsigned char, class T2=T1, class T3=T2>
class YCbCrColor : public ColorSpace<T1,T2,T3>
{
  public:
    typedef ColorSpace<T1,T2,T3> base;

  public:
    inline YCbCrColor() : base() {}
    inline YCbCrColor(T1 Y, T2 Cb, T3 Cr) : base(Y,Cb,Cr) {}

    template<class U1,class U2,class U3> inline YCbCrColor(const ColorSpace<U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline YCbCrColor(const YCbCrColor<U1,U2,U3> &r) : base(r) {}
    template<class U1,class U2,class U3> inline YCbCrColor(const RGBColor  <U1,U2,U3> &r) : base(0.2989*r.x+0.5866*r.y+0.1145*r.z, -0.1688*r.x-0.3312*r.y+0.5*r.z, 0.5*r.x-0.4184*r.y-0.0816*r.z) { }

    template<class U1,class U2,class U3> inline YCbCrColor &operator=(const RGBColor  <U1,U2,U3> &r) { *this=YCbCrColor(r); return *this; }
};


//group = Images Color Spaces Type Definitions

//{unsecret}
//Remarks:
//  The /RGB/ type is not defined because it would collapse with another definition in some Windows SDK header files.
typedef RGBAColor<unsigned char, unsigned char, unsigned char,unsigned char> RGBA;
//{unsecret}
typedef YUVColor <unsigned char, char, char> YUV;
//{unsecret}
//{PartOf:YUV}
typedef YUVAColor<unsigned char, char, char,unsigned char> YUVA;
//{unsecret}
typedef YCbCrColor<unsigned char, char, char> YCbCr;

typedef BGRColor<unsigned char> BGR;



template<class T1=unsigned char, class T2=T1, class T3=T2, class T4=unsigned char>
class RGBAColor : public AlphaColorSpace<T1,T2,T3,T4>
{
  public:
	  typedef RGBAColor self;
	  typedef AlphaColorSpace<T1,T2,T3,T4> base;
	  
	  typedef typename base::first_type  first_type;
	  typedef typename base::second_type second_type;
	  typedef typename base::third_type  third_type;
	  typedef typename base::alpha_type  alpha_type;

    using base::x;
    using base::y;
    using base::z;
	   
  public:
    inline RGBAColor() : base() {}
    inline RGBAColor(first_type R, second_type G, third_type B) : base(R,G,B) {}
    inline RGBAColor(first_type R, second_type G, third_type B, alpha_type A) : base(R,G,B,A) {}

    template<class U1,class U2,class U3> inline RGBAColor(const RGBColor<U1,U2,U3> &r) : base(r.x,r.y,r.z,0) {}
    template<class U1,class U2,class U3> inline RGBAColor(const BGRColor<U1,U2,U3> &r) : base(r.z,r.y,r.x,0) {}
    
    inline RGBAColor(const YUVA &r);
    inline RGBAColor &operator=(const YUVA &r) { *this=self(r); return *this; }
    
    inline operator unsigned char() { return 0.299*x+0.587*y+0.114*z; }
};

template<class T1=unsigned char, class T2=char, class T3=T2, class T4=unsigned char>
class YUVAColor : public AlphaColorSpace<T1,T2,T3,T4>
{
  public:
	  typedef YUVAColor self;
	  typedef AlphaColorSpace<T1,T2,T3,T4> base;    

	  typedef typename base::first_type  first_type;
	  typedef typename base::second_type second_type;
	  typedef typename base::third_type  third_type;
	  typedef typename base::alpha_type  alpha_type;

    using base::x;
    using base::y;
    using base::z;
  
  public:
    inline YUVAColor() : base() {}
    inline YUVAColor(first_type Y, second_type U, third_type V) : base(Y,U,V) {}
    inline YUVAColor(first_type Y, second_type U, third_type V, alpha_type A) : base(Y,U,V,A) {}

    template<class U1,class U2,class U3> inline YUVAColor(const RGBColor<U1,U2,U3> &r) : base(0.299*r.x+0.587*r.y+0.114*r.z, -0.147*r.x-0.289*r.y+0.436*r.z, 0.615*r.x-0.515*r.y-0.100*r.z) { }

    inline YUVAColor(const RGBA &r);
    inline YUVAColor &operator=(const RGBA &r) { *this=self(r); return *this; }
};

#if defined(MMX)
template<> inline YUVA::YUVAColor(const RGBA &x)
{	
  m64b m(m64i((int&)x));
  m64s m0 = unpacklo(m);

//0.299*r.x+0.587*r.y+0.114*r.z
//-0.147*r.x-0.289*r.y+0.437*r.z
//0.615*r.x-0.515*r.y-0.100*r.z

  m64i Y0 = madd(m0,m64s( 9798,19235,3736,0));
  m64i Y1 = madd(m0,m64s(-4817,-9470,14320,0));
  m64i Y01 = shiftra<15>(unpacklo(Y0,Y1)+unpackhi(Y0,Y1));
  
  m64i Y2 = madd(m0,m64s(10076,-8438,-1638,    0));
  m64i Y3 = madd(m0,m64s(    0,    0,    0,16384));
  m64i Y23 = shiftra<14>(unpacklo(Y2,Y3)+unpackhi(Y2,Y3));
  
  m64s Y = packs(Y01,Y23)+m64s(-128,0,0,-128);
  m64c n = packs(Y,Y)+m64c(-128,0,0,-128);
  
  (int &)(*this) = _mm_cvtsi64_si32(n);
}

template<> inline RGBA::RGBAColor(const YUVA &X)
{	
  m64b m(m64i((int&)X));
  m = m+m64b(0,128,128,0);

  m64s m0 = unpacklo(m);
  m0 = m0-m64s(0,128,128,0);
  
//r.x           +1.140 *r.z
//r.x-0.396 *r.y-0.581 *r.z
//r.x+2.029 *r.y

  m64i Y0 = madd(m0,m64s(16384,    0,18678,0));
  m64i Y1 = madd(m0,m64s(16384,-6488,-9519,0));
  m64i Y01 = shiftra<14>(unpacklo(Y0,Y1)+unpackhi(Y0,Y1));
  
  m64i Y23 = madd(m0,m64s( 8192,16622,    0, 8192));
  Y23 = shiftra<13>(Y23);

  m64s Y = packs(Y01,Y23);
  m64b n = packus(Y,Y);
  
  (int &)(*this) = _mm_cvtsi64_si32(n);
}
#endif // MMX





template<class V, class I, class Q1, class Q2=Q1, class Q3=Q2>
class ColorQuantizer
{
  public:
    typedef ColorQuantizer self;

    typedef V argument_type;
    typedef I value_type;

  private:
    Q1 q1;
    Q2 q2;
    Q3 q3;

  public:
    inline ColorQuantizer() : q1(), q2(), q3() {}
    inline ColorQuantizer(const Q1 &r) : q1(r), q2(r), q3(r) {}
    inline ColorQuantizer(const Q1 &r1, const Q2 &r2, const Q3 &r3) : q1(r1), q2(r2), q3(r3) {}

    inline operator value_type() const { return q1*q2*q3; }

    inline value_type operator()(const argument_type &v) const { return q3*(q2*q1(v.x)+q2(v.y))+q3(v.z); }
    inline argument_type operator[](value_type i) const { return argument_type(q1[i/(q2*q3)], q2[(i/q3)%q1], q3[i%(q1*q2)]); }
    inline argument_type quantize(const argument_type &v) const { return argument_type(q1.quantize(v.x),q2.quantize(v.y),q3.quantize(v.z)); }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return typename quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return typename quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return typename quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return typename quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return typename quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return typename quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class V, class I, class Q1,class Q2=Q1, class Q3=Q2>
class ColorColorQuantizer
{
  public:
    typedef ColorColorQuantizer self;

    typedef V argument_type;
    typedef I value_type;

  private:
    Q1 q1;
    Q2 q2;
    Q3 q3;

  public:
    inline ColorColorQuantizer() : q1(), q2(), q3() {}
    inline ColorColorQuantizer(const Q1 &r1                            ) : q1(r1), q2(r1), q3(r1) {}
    inline ColorColorQuantizer(const Q1 &r1, const Q2 &r2              ) : q1(r1), q2(r2), q3(r2) {}
    inline ColorColorQuantizer(const Q1 &r1, const Q2 &r2, const Q3 &r3) : q1(r1), q2(r2), q3(r3) {}

    inline operator value_type() const { return value_type(q1,q2,q3); }

    inline value_type    operator()(const argument_type &v) const { return value_type   (q1(v.x),q2(v.y),q3(v.z)); }
    inline argument_type operator[](const value_type    &n) const { return argument_type(q1[n.x],q2[n.y],q3[n.z]); }
    inline argument_type quantize  (const argument_type &v) const { return argument_type(q1.quantize(v.x),q2.quantize(v.y),q3.quantize(v.z)); }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return typename quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return typename quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return typename quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return typename quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return typename quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return typename quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};


//}

#endif

