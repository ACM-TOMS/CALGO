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

#ifndef UTIL_H
#define UTIL_H



#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <iostream>
#include <numeric>
#include <limits>
#include <bitset>

#include "functional.h"
#include "iterator.h"
#include "memory.h"
#include "exception.h"
#include "utility.h"
#include "cmath.h"
#include "ctime.h"
#include "complex.h"


//namespace genial
//{

using namespace std;

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#if !defined(__ICL) && !defined(_MSC_VER)
template<class T> inline T __min(const T &x, const T &y) { return min(x,y); }
template<class T> inline T __max(const T &x, const T &y) { return max(x,y); }
#else
inline char           max(char           x,char           y) { return __max(x,y); }
inline unsigned char  max(unsigned char  x,unsigned char  y) { return __max(x,y); }
inline short          max(short          x,short          y) { return __max(x,y); }
inline unsigned short max(unsigned short x,unsigned short y) { return __max(x,y); }
inline int            max(int            x,int            y) { return __max(x,y); }
inline unsigned int   max(unsigned int   x,unsigned int   y) { return __max(x,y); }
inline float          max(float          x,float          y) { return __max(x,y); }
inline double         max(double         x,double         y) { return __max(x,y); }
inline long double    max(long double    x,long double    y) { return __max(x,y); }

inline char           min(char           x,char           y) { return __min(x,y); }
inline unsigned char  min(unsigned char  x,unsigned char  y) { return __min(x,y); }
inline short          min(short          x,short          y) { return __min(x,y); }
inline unsigned short min(unsigned short x,unsigned short y) { return __min(x,y); }
inline int            min(int            x,int            y) { return __min(x,y); }
inline unsigned int   min(unsigned int   x,unsigned int   y) { return __min(x,y); }
inline float          min(float          x,float          y) { return __min(x,y); }
inline double         min(double         x,double         y) { return __min(x,y); }
inline long double    min(long double    x,long double    y) { return __min(x,y); }
#endif


#define JOIN(Sa,Sb) JOIN_(Sa,Sb)
#define JOIN2(Sa,Sb) JOIN_(Sa,Sb)
#define JOIN_(Sa,Sb) Sa ## Sb
#define JOIN3(Sa,Sb,Sc) JOIN2(JOIN2(Sa,Sb),Sc)
#define JOIN4(Sa,Sb,Sc,Sd) JOIN2(JOIN3(Sa,Sb,Sc),Sd)
#define JOIN5(Sa,Sb,Sc,Sd,Se) JOIN3(JOIN3(Sa,Sb,Sc),Sd,Se)


template<int m,int n> struct min_traits { enum { RET=(m<n)?m:n }; };
template<int m,int n> struct max_traits { enum { RET=(m>n)?m:n }; };

#define MIN(m,n) min_traits<m,n>::RET
#define MAX(m,n) max_traits<m,n>::RET


template<int cond,class A, class B> struct IF {};
template<class A,class B> struct IF<0, A, B> { typedef B RET; };
template<class A,class B> struct IF<1, A, B> { typedef A RET; };

template<int m,int n> struct EQUAL         { enum { RET=(m==n)?1:0 }; };
template<int m,int n> struct NOT_EQUAL     { enum { RET=(m!=n)?1:0 }; };
template<int m,int n> struct LOWER         { enum { RET=(m< n)?1:0 }; };
template<int m,int n> struct GREATER       { enum { RET=(m> n)?1:0 }; };
template<int m,int n> struct LOWER_EQUAL   { enum { RET=(m<=n)?1:0 }; };
template<int m,int n> struct GREATER_EQUAL { enum { RET=(m>=n)?1:0 }; };




template<class T> struct promotion {  };

#define DECLARE_PROMOTION(TYPE,RANK) template<> struct promotion< TYPE > { enum { rank = RANK }; }; template<> struct promotion<const TYPE > { enum { rank = RANK }; };

DECLARE_PROMOTION(unsigned char,         1)
DECLARE_PROMOTION(char,                  2)
DECLARE_PROMOTION(unsigned short,        3)
DECLARE_PROMOTION(short,                 4)
DECLARE_PROMOTION(unsigned int,          5)
DECLARE_PROMOTION(int,                   6)
DECLARE_PROMOTION(unsigned long,         7)
DECLARE_PROMOTION(long,                  8)
DECLARE_PROMOTION(float,                 9)
DECLARE_PROMOTION(double,               10)
DECLARE_PROMOTION(long double,          11)

template<class T>
struct promotion_traits
{
  typedef T value_type;
};

template<class T>
struct promotion_traits<const T>
{
  typedef typename promotion_traits<T>::value_type value_type;
};

template<class T1,class T2>
struct promotion2_traits
{
  typedef typename IF<LOWER<promotion<T1>::rank,promotion<T2>::rank>::RET, T2, T1>::RET value_type;
//  typedef T1 value_type;
};

template<class T1,class T2> struct promotion2_traits<const T1,      T2> { typedef typename promotion2_traits<T1,T2>::value_type value_type; };
template<class T1,class T2> struct promotion2_traits<      T1,const T2> { typedef typename promotion2_traits<T1,T2>::value_type value_type; };
template<class T1,class T2> struct promotion2_traits<const T1,const T2> { typedef typename promotion2_traits<T1,T2>::value_type value_type; };

template<class G> class Vector;
template<class G,class V> struct promotion2_traits<Vector <G>,complex<V> > { typedef typename Vector<G>::template rebind<typename promotion2_traits<typename Vector<G>::const_value_type,complex<V> >::value_type>::other value_type; };
template<class V,class G> struct promotion2_traits<complex<V>,Vector <G> > { typedef typename Vector<G>::template rebind<typename promotion2_traits<complex<V>,typename Vector<G>::const_value_type >::value_type>::other value_type; };

template<class G> class Matrix;
template<class G,class V> struct promotion2_traits<Matrix <G>,complex<V> > { typedef typename Matrix<G>::template rebind<typename promotion2_traits<typename Matrix<G>::const_value_type,complex<V> >::value_type>::other value_type; };
template<class V,class G> struct promotion2_traits<complex<V>,Matrix <G> > { typedef typename Matrix<G>::template rebind<typename promotion2_traits<complex<V>,typename Matrix<G>::const_value_type >::value_type>::other value_type; };

template<int D,class G> class Array;
template<int D,class G,class V>   struct promotion2_traits<Array<D,G >,complex<V>  > { typedef typename Array<D,G>::template rebind<typename promotion2_traits<typename Array<D,G>::const_value_type,complex<V> >::value_type>::other value_type; };
template<class V,int D,class G>   struct promotion2_traits<complex<V> ,Array<D,G > > { typedef typename Array<D,G>::template rebind<typename promotion2_traits<complex<V> ,typename Array<D,G>::const_value_type>::value_type>::other value_type; };


template<class T1,class T2,class T3>
struct promotion3_traits
{
  typedef typename promotion2_traits<typename promotion2_traits<T1,T2>::value_type, T3>::value_type value_type;
//  typedef T1 value_type;
};

#define PROMOTE(T) typename promotion_traits< T >::value_type
//#define PROMOTE2(T1,T2) T1
#define PROMOTE2(T1,T2) typename promotion2_traits< T1,T2 >::value_type
//#define PROMOTE3(T1,T2,T3) T1
#define PROMOTE3(T1,T2,T3) typename promotion3_traits< T1,T2,T3 >::value_type



template<class V> struct value_type_traits { typedef V value_type; };

template<class T> struct value_type_traits;
template<class V> struct value_type_traits<      complex<V> > { typedef typename complex<V>::value_type       value_type; };
template<class V> struct value_type_traits<const complex<V> > { typedef const typename complex<V>::value_type value_type; };

#define VALUE_TYPE(T) typename value_type_traits< T >::value_type


template<class T,class U> inline T bound(U x)
{
  if (x<numeric_limits<T>::min()) return numeric_limits<T>::min();
  if (x>numeric_limits<T>::max()) return numeric_limits<T>::max();
  return (T)x;
}

template<class T> inline T clip_low (const T &x, const T &y            ) { return (y< x)*x; }
template<class T> inline T clip_low (const T &x, const T &y,const T &x0) { return (y< x)?x:x0; }
template<class T> inline T clip_high(const T &x, const T &y            ) { return (y> x)?x:numeric_limits<T>::max(); }
template<class T> inline T clip_high(const T &x, const T &y,const T &x0) { return (y> x)?x:x0; }



inline int to_int(const string &str) { int i; istringstream(str.c_str()) >> i; return i; }

template<class E, class T, class A, class X> inline void operator>>(const basic_string<E,T,A> &s, X &x) { istringstream(s)>>x; return; }

template<class T> string to_string(const T &x) { ostringstream oss; oss <<  x; return oss.str(); }
inline string to_string(bool          x) { ostringstream oss; oss <<               x; return oss.str(); }
inline string to_string(         int  x) { ostringstream oss; oss <<               x; return oss.str(); }
inline string to_string(unsigned int  x) { ostringstream oss; oss <<               x; return oss.str(); }
inline string to_string(         char x) { ostringstream oss; oss << (         int)x; return oss.str(); }
inline string to_string(unsigned char x) { ostringstream oss; oss << (unsigned int)x; return oss.str(); }


template<class InIt> string to_string(InIt begin, InIt end)
{
  ostringstream oss;
  if (begin!=end)
  {
    oss<<*begin;
    for (++begin; begin!=end; ++begin) 
      oss<<" "<<*begin;
  }
  return oss.str();
}

template<size_t N> string to_string(const bitset<N> &x) { return x.to_string(); }

inline string to_string(         int  x, int w, char c=' ') { ostringstream oss; oss.width(w); oss.fill(c); oss <<               x; return oss.str(); }
inline string to_string(unsigned int  x, int w, char c=' ') { ostringstream oss; oss.width(w); oss.fill(c); oss <<               x; return oss.str(); }
inline string to_string(         char x, int w, char c=' ') { ostringstream oss; oss.width(w); oss.fill(c); oss << (         int)x; return oss.str(); }
inline string to_string(unsigned char x, int w, char c=' ') { ostringstream oss; oss.width(w); oss.fill(c); oss << (unsigned int)x; return oss.str(); }
inline string to_string(float         x, int w, char c=' ') { ostringstream oss; oss.width(w); oss.fill(c); oss <<               x; return oss.str(); }
inline string to_string(double        x, int w, char c=' ') { ostringstream oss; oss.width(w); oss.fill(c); oss <<               x; return oss.str(); }

template<class E,class Tr,class T,class U> inline basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const pair<T,U> &x) { os << "(" << x.first << ", " << x.second << ")"; return os; }
template<class E,class Tr,class T,class A> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &s, const vector<T,A> &v) { /*s << v.size();*/ s << "[ "; for (typename vector<T,A>::const_iterator i=v.begin(); i!=v.end(); ++i) s << *i << " "; s << "]"; return s; }
template<class E,class Tr,class T,class A> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &s, const list  <T,A> &l) { /*s << l.size();*/ s << "{ "; for (typename list  <T,A>::const_iterator i=l.begin(); i!=l.end(); ++i) s << *i << " "; s << "}"; return s; }


class File
{
  private:
    string fname;

  public:
    File() : fname() {}
    File(const char *s) : fname(s) {}
    File(const string &s) : fname(s) {}
  
    string       &name()       { return fname; }
    const string &name() const { return fname; }
};



//}

#endif
