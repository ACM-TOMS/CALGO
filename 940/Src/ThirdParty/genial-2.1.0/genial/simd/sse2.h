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

#ifndef SSE2_H
#define SSE2_H

#ifndef NO_SIMD

#ifndef SSE2
#define SSE2
#endif

#include <emmintrin.h>
#include "sse.h"

using namespace std;

class m128d;
class m256d;
class m128cd;

class m128c;
class m128b;
class m128s;
class m128w;
class m128i;

//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 2 doubles
class m64d
{
  public:
    typedef m64d self;
	  typedef double value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=1 };

  public:
    __aligned __m128d vec;

  public:
    inline m64d(        ) {}
    explicit inline m64d(__m128d m) : vec(m) {}
    inline explicit m64d(double d) : vec(_mm_set1_pd(d))	{}
    inline m64d(const self &x) : vec(x.vec) {}

    inline m64d &operator=(double  d) { vec = _mm_set1_pd(d); return *this; }
    inline m64d &operator=(__m128d x) { vec = x;              return *this; }
    inline m64d &operator=(const self &x) { vec=x.vec; return *this; }

    inline operator __m128d() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};

template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64d &x)	 { double *p = (double*)&x; os << '(' << *p << ')'; return os; }

inline void store  (double *p, const m64d &x) { _mm_storel_pd(p,x.vec); }
inline void storeu (double *p, const m64d &x) { store(p,x); }
inline void storer (double *p, const m64d &x) { store(p,x); }
inline void storeur(double *p, const m64d &x) { store(p,x); }
inline void store0 (double *p, const m64d &x) { store(p,x); }
inline void storelh(double *p1, double *p2, const m64d &x) { store(p1,x); store(p2,x); }

inline m64d operator- (               const m64d &y) { return m64d(_mm_sub_pd(_mm_setzero_pd(),y.vec)); }
inline m64d operator+ (const m64d &x, const m64d &y) { return m64d(_mm_add_pd(x.vec,y.vec)); }
inline m64d operator- (const m64d &x, const m64d &y) { return m64d(_mm_sub_pd(x.vec,y.vec)); }
inline m64d operator* (const m64d &x, const m64d &y) { return m64d(_mm_mul_pd(x.vec,y.vec)); }
inline m64d operator/ (const m64d &x, const m64d &y) { return m64d(_mm_div_pd(x.vec,y.vec)); }
inline m64d operator* (const m64d &x,      double y) { return m64d(_mm_mul_pd(x.vec,_mm_set1_pd(y))); }
inline m64d operator/ (const m64d &x,      double y) { return m64d(_mm_div_pd(x.vec,_mm_set1_pd(y))); }
inline m64d operator* (     double x, const m64d &y) { return m64d(_mm_mul_pd(_mm_set1_pd(x),y.vec)); }
inline m64d operator/ (     double x, const m64d &y) { return m64d(_mm_div_pd(_mm_set1_pd(x),y.vec)); }
inline m64d inv       (const m64d &x               ) { return double(1)/x; }
inline m64d operator& (const m64d &x, const m64d &y) { return m64d(_mm_and_pd   (x.vec,y.vec)); }
inline m64d operator| (const m64d &x, const m64d &y) { return m64d(_mm_or_pd    (x.vec,y.vec)); }
inline m64d operator^ (const m64d &x, const m64d &y) { return m64d(_mm_xor_pd   (x.vec,y.vec)); }
inline m64d and_not   (const m64d &x, const m64d &y) { return m64d(_mm_andnot_pd(x.vec,y.vec)); }
inline m64d&operator+=(      m64d &x, const m64d &y) { return x=x+y; }
inline m64d&operator-=(      m64d &x, const m64d &y) { return x=x-y; }
inline m64d&operator*=(      m64d &x, const m64d &y) { return x=x*y; }
inline m64d&operator/=(      m64d &x, const m64d &y) { return x=x/y; }
inline m64d&operator*=(      m64d &x,      double y) { return x=x*y; }
inline m64d&operator/=(      m64d &x,      double y) { return x=x/y; }
inline m64d&operator&=(      m64d &x, const m64d &y) { return x=x&y; }
inline m64d&operator|=(      m64d &x, const m64d &y) { return x=x|y; }
inline m64d&operator^=(      m64d &x, const m64d &y) { return x=x^y;}

inline m64d operator==(const m64d &x, const m64d &y) { return m64d(_mm_cmpeq_pd (x.vec,y.vec)); }
inline m64d operator!=(const m64d &x, const m64d &y) { return m64d(_mm_cmpneq_pd(x.vec,y.vec)); }
inline m64d operator> (const m64d &x, const m64d &y) { return m64d(_mm_cmpgt_pd (x.vec,y.vec)); }
inline m64d operator>=(const m64d &x, const m64d &y) { return m64d(_mm_cmpge_pd (x.vec,y.vec)); }
inline m64d operator< (const m64d &x, const m64d &y) { return m64d(_mm_cmplt_pd (x.vec,y.vec)); }
inline m64d operator<=(const m64d &x, const m64d &y) { return m64d(_mm_cmple_pd (x.vec,y.vec)); }

inline m64d sqrt (const m64d &x) { return m64d(_mm_sqrt_pd (x.vec)); }

inline m64d conj (const m64d &x) { return x; }

inline m64d min(const m64d &x, const m64d &y) { return m64d(_mm_min_pd(x.vec,y.vec)); }
inline m64d max(const m64d &x, const m64d &y) { return m64d(_mm_max_pd(x.vec,y.vec)); }

//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 2 doubles
class m128d
{
  public:
    typedef m128d self;
	  typedef double value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=2 };

  public:
    __aligned __m128d vec;

  public:
    inline m128d(        ) {}
    explicit inline m128d(__m128d m) : vec(m) {}
    inline m128d(double d0, double d1) : vec(_mm_set_pd(d1,d0)) {}
    inline explicit m128d(double d) : vec(_mm_set1_pd(d))	{}
    inline m128d(const self &x) : vec(x.vec) {}

    inline m128d &operator=(double  d) { vec = _mm_set1_pd(d); return *this; }
    inline m128d &operator=(__m128d x) { vec = x;              return *this; }
    inline m128d &operator=(const self &x) { vec=x.vec; return *this; }

    inline operator __m128d() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};

template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m128d &x)	 { double *dp = (double*)&x; os << '(' << *dp << ',' << *(dp+1) << ')'; return os; }

inline m128d operator- (                const m128d &y) { return m128d(_mm_sub_pd(_mm_setzero_pd(),y.vec)); }
inline m128d operator+ (const m128d &x, const m128d &y) { return m128d(_mm_add_pd(x.vec,y.vec)); }
inline m128d operator- (const m128d &x, const m128d &y) { return m128d(_mm_sub_pd(x.vec,y.vec)); }
inline m128d operator* (const m128d &x, const m128d &y) { return m128d(_mm_mul_pd(x.vec,y.vec)); }
inline m128d operator/ (const m128d &x, const m128d &y) { return m128d(_mm_div_pd(x.vec,y.vec)); }
inline m128d operator* (const m128d &x,       double y) { return m128d(_mm_mul_pd(x.vec,_mm_set1_pd(y))); }
inline m128d operator/ (const m128d &x,       double y) { return m128d(_mm_div_pd(x.vec,_mm_set1_pd(y))); }
inline m128d operator* (      double x, const m128d &y) { return m128d(_mm_mul_pd(_mm_set1_pd(x),y.vec)); }
inline m128d operator/ (      double x, const m128d &y) { return m128d(_mm_div_pd(_mm_set1_pd(x),y.vec)); }
inline m128d inv       (const m128d &x                ) { return double(1)/x; }
inline m128d operator& (const m128d &x, const m128d &y) { return m128d(_mm_and_pd   (x.vec,y.vec)); }
inline m128d operator| (const m128d &x, const m128d &y) { return m128d(_mm_or_pd    (x.vec,y.vec)); }
inline m128d operator^ (const m128d &x, const m128d &y) { return m128d(_mm_xor_pd   (x.vec,y.vec)); }
inline m128d and_not   (const m128d &x, const m128d &y) { return m128d(_mm_andnot_pd(x.vec,y.vec)); }
inline m128d&operator+=(      m128d &x, const m128d &y) { return x=x+y; }
inline m128d&operator-=(      m128d &x, const m128d &y) { return x=x-y; }
inline m128d&operator*=(      m128d &x, const m128d &y) { return x=x*y; }
inline m128d&operator/=(      m128d &x, const m128d &y) { return x=x/y; }
inline m128d&operator*=(      m128d &x,       double y) { return x=x*y; }
inline m128d&operator/=(      m128d &x,       double y) { return x=x/y; }
inline m128d&operator&=(      m128d &x, const m128d &y) { return x=x&y; }
inline m128d&operator|=(      m128d &x, const m128d &y) { return x=x|y; }
inline m128d&operator^=(      m128d &x, const m128d &y) { return x=x^y;}

inline m128d operator==(const m128d &x, const m128d &y) { return m128d(_mm_cmpeq_pd (x.vec,y.vec)); }
inline m128d operator!=(const m128d &x, const m128d &y) { return m128d(_mm_cmpneq_pd(x.vec,y.vec)); }
inline m128d operator> (const m128d &x, const m128d &y) { return m128d(_mm_cmpgt_pd (x.vec,y.vec)); }
inline m128d operator>=(const m128d &x, const m128d &y) { return m128d(_mm_cmpge_pd (x.vec,y.vec)); }
inline m128d operator< (const m128d &x, const m128d &y) { return m128d(_mm_cmplt_pd (x.vec,y.vec)); }
inline m128d operator<=(const m128d &x, const m128d &y) { return m128d(_mm_cmple_pd (x.vec,y.vec)); }

inline m128d flip   (const m128d &x) { return m128d(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,1))); }
inline m128d conj   (const m128d &x) { return x; }

inline m128d  unpacklo (const m128d &X, const m128d &Y) { return m128d (_mm_unpacklo_pd(X.vec,Y.vec)); }
inline m128d  unpackhi (const m128d &X, const m128d &Y) { return m128d (_mm_unpackhi_pd(X.vec,Y.vec)); }

inline void trn    (m128d &x0, m128d &x1) { m128d a=unpacklo(x0,x1); x1=unpackhi(x0,x1); x0=a; }
inline void up_trn (m128d &x0, m128d &x1) { x1=unpackhi(x0,x1); }
inline void low_trn(m128d &x0, m128d &x1) { x0=unpacklo(x0,x1); }

inline m128d sqrt (const m128d &x) { return m128d(_mm_sqrt_pd (x.vec)); }

template<int i        > inline m128d shuffle(const m128d &x                ) { return m128d(_mm_shuffle_pd(x.vec,x.vec,i)); }
template<int i0,int i1> inline m128d shuffle(const m128d &x                ) { return shuffle<_MM_SHUFFLE2(i0,i1)>(x); }
template<int i        > inline m128d shuffle(const m128d &x, const m128d &y) { return m128d(_mm_shuffle_pd(x.vec,y.vec,i)); }
template<int i0,int i1> inline m128d shuffle(const m128d &x, const m128d &y) { return shuffle<_MM_SHUFFLE2(i0,i1)>(x,y); }
template<int i> inline m128d repeat (const m128d &x) { return shuffle<i,i>(x,x); }

inline void loadl (m128d &x, const double *p) { x.vec=_mm_loadl_pd (x.vec,p); }
inline void loadh (m128d &x, const double *p) { x.vec=_mm_loadh_pd (x.vec,p); }
inline void load0 (m128d &x, const double *p) { loadl(x,p); }
inline void load1 (m128d &x, const double *p) { loadh(x,p); }
inline void load  (m128d &x, const double *p) { x.vec = _mm_load_pd (p); }
inline void loadu (m128d &x, const double *p) { x.vec = _mm_loadu_pd(p); }
inline void loadr (m128d &x, const double *p) { x.vec = _mm_loadr_pd(p); }
inline void loadlh(m128d &x, const double *p1, const double *p2) { loadl(x,p1); loadh(x,p2); }

inline void storel (double *p, const m128d &x) { _mm_storel_pd(p,x.vec); }
inline void storeh (double *p, const m128d &x) { _mm_storeh_pd(p,x.vec); }
inline void store0 (double *p, const m128d &x) { storel(p,x); }
inline void store1 (double *p, const m128d &x) { storeh(p,x); }
inline void store  (double *p, const m128d &x) { _mm_store_pd (p, x.vec); }
inline void stream (double *p, const m128d &x) { _mm_stream_pd(p, x.vec); }
inline void storeu (double *p, const m128d &x) { _mm_storeu_pd(p, x.vec); }
inline void storer (double *p, const m128d &x) { _mm_storer_pd(p, x.vec); }
inline void storeur(double *p, const m128d &x) { _mm_storeu_pd(p, flip(x).vec); }
inline void storelh(double *p1, double *p2, const m128d &x) { storel(p1,x); storeh(p2,x); }

inline double get0(const m128d &x) { return *(double*)&(x.vec); }
inline double get1(const m128d &x) { m128d a=shuffle<0,1>(x,x); return get0(a); }

inline m128d min(const m128d &x, const m128d &y) { return m128d(_mm_min_pd(x.vec,y.vec)); }
inline m128d max(const m128d &x, const m128d &y) { return m128d(_mm_max_pd(x.vec,y.vec)); }

inline m128d min(const m128d &x) { return m128d(min(x,shuffle<_MM_SHUFFLE2(0,1)>(x,x))); }
inline m128d max(const m128d &x) { return m128d(max(x,shuffle<_MM_SHUFFLE2(0,1)>(x,x))); }

inline m128d hadd(const m128d &x, const m128d &y);
inline m128d hsub(const m128d &x, const m128d &y);
inline m128d sum  (const m128d &x);
inline m128d sumlo(const m128d &x) { return x; }
inline m128d sumhi(const m128d &x) { return unpackhi(x,x); }
inline m128d sum  (const m128d &x, const m128d &y);

#ifndef SSE3
inline m128d sum (const m128d &x) { return x+shuffle<_MM_SHUFFLE2(0,1)>(x,x); }
inline m128d sum (const m128d &x, const m128d &y) { return m128d(unpacklo(x,y)+unpackhi(x,y)); }
inline m128d hadd(const m128d &x, const m128d &y) { return m128d(unpacklo(x,y)+unpackhi(x,y)); }
inline m128d hsub(const m128d &x, const m128d &y) { return m128d(unpacklo(x,y)-unpackhi(x,y)); }
#endif



//{secret}
//{group:SSE2}
//Summary: Two SSE registers repesenting 4 doubles
class m256d
{
  public:
    typedef m256d self;
    typedef double value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=4 };    

  public:
    m128d vec0;
    m128d vec1;

  public:
    inline m256d() {};
	  inline m256d(const __m128d &m0,const __m128d &m1) : vec0(m0),vec1(m1)  {};
	  inline m256d(const   m128d &m0,const   m128d &m1) : vec0(m0),vec1(m1)  {};
    inline explicit m256d(const value_type &x) { vec0=m128d(x); vec1=vec0; }
    inline m256d(const value_type &x0, const value_type &x1, const value_type &x2, const value_type &x3) : vec0(x0,x1), vec1(x2,x3) {}
    inline m256d(const self &x) : vec0(x.vec0), vec1(x.vec1) {}

    inline self& operator=(const value_type &x) { vec0=m128d(x); vec1=vec0; return *this; }
    inline self& operator=(const self &x) { vec0=x.vec0; vec1=x.vec1; return *this; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec0)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec0)+i); };
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m256d &x)	{ double *p = (double*) &x; return os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ')'; }

inline void load (m256d &x, const double *p) { load (x.vec0,p); load (x.vec1,p+2); }
inline void loadu(m256d &x, const double *p) { loadu(x.vec0,p); loadu(x.vec1,p+2); }
inline void loadr(m256d &x, const double *p) { loadr(x.vec1,p); loadr(x.vec0,p+2); }

inline void store  (double *p, const m256d &x) { store (p,x.vec0); store (p+2,x.vec1); }
inline void stream (double *p, const m256d &x) { stream(p,x.vec0); stream(p+2,x.vec1); }
inline void storeu (double *p, const m256d &x) { storeu(p,x.vec0); storeu(p+2,x.vec1); }

inline m256d operator- (                const m256d &y) { return m256d(      -y.vec0,      -y.vec1); }
inline m256d operator+ (const m256d &x, const m256d &y) { return m256d(x.vec0+y.vec0,x.vec1+y.vec1); }
inline m256d operator- (const m256d &x, const m256d &y) { return m256d(x.vec0-y.vec0,x.vec1-y.vec1); }
inline m256d operator* (const m256d &x, const m256d &y) { return m256d(x.vec0*y.vec0,x.vec1*y.vec1); }
inline m256d operator/ (const m256d &x, const m256d &y) { return m256d(x.vec0/y.vec0,x.vec1/y.vec1); }
inline m256d operator* (const m256d &x,       double y) { return m256d(x.vec0*y     ,x.vec1*y     ); }
inline m256d operator/ (const m256d &x,       double y) { return m256d(x.vec0/y     ,x.vec1/y     ); }
inline m256d operator* (      double x, const m256d &y) { return y*x; }
inline m256d operator& (const m256d &x, const m256d &y) { return m256d(x.vec0&y.vec0,x.vec1&y.vec1); }
inline m256d operator| (const m256d &x, const m256d &y) { return m256d(x.vec0|y.vec0,x.vec1|y.vec1); }
inline m256d operator^ (const m256d &x, const m256d &y) { return m256d(x.vec0^y.vec0,x.vec1^y.vec1); }
inline m256d&operator+=(      m256d &x, const m256d &y) { return x=x+y; }
inline m256d&operator-=(      m256d &x, const m256d &y) { return x=x-y; }
inline m256d&operator*=(      m256d &x, const m256d &y) { return x=x*y; }
inline m256d&operator/=(      m256d &x, const m256d &y) { return x=x/y; }
inline m256d&operator*=(      m256d &x,       double y) { return x=x*y; }
inline m256d&operator/=(      m256d &x,       double y) { return x=x/y; }
inline m256d&operator&=(      m256d &x, const m256d &y) { return x=x&y; }
inline m256d&operator|=(      m256d &x, const m256d &y) { return x=x|y; }
inline m256d&operator^=(      m256d &x, const m256d &y) { return x=x^y; }

inline m256d  sum  (const m128d &x0, const m128d &x1, const m128d &x2, const m128d &x3) { return m256d(sum(x0,x1),sum(x2,x3)); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 1 double precision complex
class m128cd
{
  public:
	  typedef m128cd self;
	  typedef complex<double> value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=1 };

  public:
	  __aligned __m128d vec;

  public:
	  inline m128cd() {};
	  explicit inline m128cd(__m128d m) : vec(m)	{};
	  inline m128cd(double d0, double d1) : vec(_mm_set_pd(d1,d0))	{}
	  inline m128cd(const value_type &x) : vec(_mm_set_pd(x.imag(),x.real()))	{}
	  inline m128cd(const self &x) : vec(x.vec)	{}

	  inline self& operator=(const value_type &x) { vec=_mm_loadu_pd((double*)&x);	return *this;}
	  inline self& operator=(const self       &x) { vec=x.vec;                     return *this;}

	  inline operator __m128d() const { return vec; }
	  inline operator value_type() const { value_type z; _mm_storeu_pd((double*)&z, vec); return z; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m128cd &x) { double *dp = (double*)&x; os << '(' << *dp << ',' << *(dp+1) << ')'; return os; }

inline m128cd imul(const m128cd &x);
inline m128cd flip_ri(const m128cd &x)	{return m128cd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,1))); }

inline m128cd operator-(                 const m128cd &y) { return m128cd(_mm_sub_pd(_mm_setzero_pd(),y.vec)); }
inline m128cd operator+(const m128cd &x, const m128cd &y)	{ return m128cd(_mm_add_pd(x.vec,y.vec)); }
inline m128cd operator-(const m128cd &x, const m128cd &y)	{ return m128cd(_mm_sub_pd(x.vec,y.vec)); }
inline m128cd operator*(const m128cd &x, const m128cd &y);
inline m128cd operator*(const m128cd &x, const complex<double> &y)	{ __m128d a=_mm_mul_pd(x.vec,_mm_set1_pd(y.real())); __m128d b=_mm_mul_pd((flip_ri(x)).vec,_mm_set_pd(y.imag(),-y.imag())); return m128cd(_mm_add_pd(a,b)); }
inline m128cd operator*(const m128cd &x,       double           y)	{ return m128cd(_mm_mul_pd(x.vec,_mm_set1_pd(y))); }
inline m128cd operator/(const m128cd &x,       double           y)	{ return m128cd(_mm_div_pd(x.vec,_mm_set1_pd(y))); }
inline m128cd operator*(const complex<double> &x, const m128cd &y)	{ return y*x;}
inline m128cd operator*(      double           x, const m128cd &y)	{ return y*x; }
inline m128cd&operator+=(      m128cd &x, const m128cd &y)	{ return x=x+y; }
inline m128cd&operator-=(      m128cd &x, const m128cd &y)	{ return x=x-y; }
inline m128cd&operator*=(      m128cd &x, const m128cd &y)  { return x=x*y; }
inline m128cd&operator*=(      m128cd &x, const complex<double> &y)	{ return x=x*y; }
inline m128cd&operator*=(      m128cd &x,       double           y)	{ return x=x*y; }
inline m128cd&operator/=(      m128cd &x,       double           y)	{ return x=x/y; }

inline m128cd flip   (const m128cd &x)	{return x; }
inline m128cd conj(const m128cd &x)	{return m128cd(_mm_mul_pd(_mm_set_pd(-1,1),x.vec));}

inline complex<double> get0(const m128cd &x) { return *(complex<double>*)&(x.vec); }

inline void store (complex<double> *p, const m128cd &x) { _mm_store_pd ((double *)p, x.vec); }
inline void storeu(complex<double> *p, const m128cd &x) { _mm_storeu_pd((double *)p, x.vec); }
inline void storer(complex<double> *p, const m128cd &x) { _mm_storer_pd((double *)p, x.vec); }
inline void store0(complex<double> *p, const m128cd &x) { store(p,x); }
inline void storel(complex<double> *p, const m128cd &x) { _mm_storel_pd((double *)p,x.vec); }
inline void storeh(complex<double> *p, const m128cd &x) { _mm_storeh_pd((double *)p,x.vec); }
inline void stream(complex<double> *p, const m128cd &x) { _mm_stream_pd((double *)p, x.vec); }
inline void storelh(complex<double> *p1, complex<double> *p2, const m128cd &x) { storel(p1,x); storeh(p2,x); }

inline void loadl (m128cd &x, const complex<double> *p) { x.vec=_mm_loadl_pd (x.vec,(double *)p); }
inline void loadh (m128cd &x, const complex<double> *p) { x.vec=_mm_loadh_pd (x.vec,(double *)p); }
inline void loadlh(m128cd &x, const complex<double> *p1, const complex<double> *p2) { loadl(x,p1); loadh(x,p2); }


inline void load (m128cd &x, const complex<double> *p) { x.vec = _mm_load_pd ((double *)p); }
inline void loadu(m128cd &x, const complex<double> *p) { x.vec = _mm_loadu_pd((double *)p); }
inline void loadr(m128cd &x, const complex<double> *p) { x.vec = _mm_loadr_pd((double *)p); }

inline m128cd operator*(const m128cd &x, const m128cd &y);
inline m128cd cmul(const m128cd &x, const m128cd &y);
inline m128cd icmul(const m128cd &x, const m128cd &y);
inline m128cd imul    (const m128cd &x);
inline m128cd mimul   (const m128cd &x);
inline m128cd mimul(const m128cd &x, const m128cd &y);
inline m128cd imul (const m128cd &x, const m128cd &y);
inline m128cd add_imul (const m128cd &a, const m128cd &x);
inline m128cd sub_imul (const m128cd &a, const m128cd &x);
inline m128cd addsub   (const m128cd &x, const m128cd &y);
inline m128cd subadd   (const m128cd &x, const m128cd &y);

inline m64d  norm(const m128cd &x);
inline m128d norm(const m128cd &x,const m128cd &y);
inline m64d  abs (const m128cd &x) { return sqrt(norm(x));}
inline m128d norm(const m128cd &x,const m128cd &y);
inline m128d abs (const m128cd &x,const m128cd &y) { return sqrt(norm(x,y));}
inline m128d real(const m128cd &x) { return m128d(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,0))); }
inline m128d imag(const m128cd &x) { return m128d(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1))); }

inline m128cd conj   (double a, const m128cd &x) { return m128cd(_mm_mul_pd(_mm_set_pd(-a,a),x.vec)); }
inline m128cd imul   (double a, const m128cd &x) { return flip_ri(conj(a,x)); }

inline m128cd complexlo(const m128d &X, const m128d &Y) { return m128cd(_mm_unpacklo_pd(X.vec,Y.vec)); }
inline m128cd complexhi(const m128d &X, const m128d &Y) { return m128cd(_mm_unpackhi_pd(X.vec,Y.vec)); }


#if !defined(SSE3)
inline m64d norm(const m128cd &x)
{
  __m128d a=_mm_mul_pd(x.vec,x.vec);
  return m64d(_mm_add_pd(_mm_unpacklo_pd(a,a),_mm_unpackhi_pd(a,a)));
}
inline m128d norm(const m128cd &x,const m128cd &y)
{
  __m128d a0=_mm_mul_pd(x.vec,x.vec);
  __m128d a1=_mm_mul_pd(y.vec,y.vec);
  return m128d(_mm_add_pd(_mm_unpacklo_pd(a0,a1),_mm_unpackhi_pd(a0,a1)));
}

inline m128cd operator*(const m128cd &x, const m128cd &y)
{
  __m128d a=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,0)),y.vec);
  __m128d b=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1)),imul(y).vec);
  return m128cd(_mm_add_pd(a,b));
}
inline m128cd cmul(const m128cd &x, const m128cd &y)
{
  __m128d a=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,0)),y.vec);
  __m128d b=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1)),imul(y).vec);
  return m128cd(_mm_sub_pd(a,b));
}
inline m128cd icmul(const m128cd &x, const m128cd &y)
{
  __m128d a=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1)),y.vec);
  __m128d b=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,0)),imul(y).vec);
  return m128cd(_mm_add_pd(a,b));
}

inline m128cd mimul(const m128cd &x, const m128cd &y)
{
  __m128d a=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1)),y.vec);
  __m128d b=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,0)),imul(y).vec);
  return m128cd(_mm_sub_pd(a,b));
}

inline m128cd imul (const m128cd &x, const m128cd &y)
{
  __m128d a=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1)),y.vec);
  __m128d b=_mm_mul_pd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(0,0)),imul(y).vec);
  return m128cd(_mm_sub_pd(b,a));
}


inline m128cd imul    (const m128cd &x)	{return flip_ri(conj(x));}
inline m128cd mimul   (const m128cd &x)	{return flip_ri(x);}
inline m128cd add_imul (const m128cd &a, const m128cd &x) { return a+imul (x); }
inline m128cd sub_imul (const m128cd &a, const m128cd &x) { return a-imul (x); }
inline m128cd addsub   (const m128cd &x, const m128cd &y) { return x+conj(y); }
inline m128cd subadd   (const m128cd &x, const m128cd &y) { return x-conj(y); }
#endif


//{secret}
//{group:SSE2}
//Summary: Two SSE registers repesenting 2 double precision complexes
class m256cd
{
  public:
    typedef m256cd self;
    typedef complex<double> value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=2 };

  public:
    m128cd vec0;
    m128cd vec1;

  public:
    inline m256cd() {};
	  inline m256cd(const __m128d  &m0,const __m128d  &m1) : vec0(m0),vec1(m1)  {};
	  inline m256cd(const   m128cd &m0,const   m128cd &m1) : vec0(m0),vec1(m1)  {};
    inline explicit m256cd(const value_type &x) { vec0=m128cd(x); vec1=vec0; }
    inline m256cd(const value_type &x0, const value_type &x1) : vec0(x0), vec1(x1) {}
    inline m256cd(const self &x) : vec0(x.vec0), vec1(x.vec1) {}

    inline self& operator=(const value_type &x) { vec0=m128cd(x); vec1=vec0; return *this; }
    inline self& operator=(const self &x) { vec0=x.vec0; vec1=x.vec1; return *this; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec0)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec0)+i); };
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m256cd &x)	{ double *p = (double*) &x; return os << '(' << *p << ',' << *(p+1) << ')'; }

inline void load (m256cd &x, const complex<double> *p) { load (x.vec0,p); load (x.vec1,p+1); }
inline void loadu(m256cd &x, const complex<double> *p) { loadu(x.vec0,p); loadu(x.vec1,p+1); }
inline void loadr(m256cd &x, const complex<double> *p) { loadr(x.vec1,p); loadr(x.vec0,p+1); }

inline void store (complex<double> *p, const m256cd &x) { store (p,x.vec0); store (p+1,x.vec1); }
inline void stream(complex<double> *p, const m256cd &x) { stream(p,x.vec0); stream(p+1,x.vec1); }
inline void storeu(complex<double> *p, const m256cd &x) { storeu(p,x.vec0); storeu(p+1,x.vec1); }

inline m256cd operator- (                 const m256cd &y) { return m256cd(      -y.vec0,      -y.vec1); }
inline m256cd operator+ (const m256cd &x, const m256cd &y) { return m256cd(x.vec0+y.vec0,x.vec1+y.vec1); }
inline m256cd operator- (const m256cd &x, const m256cd &y) { return m256cd(x.vec0-y.vec0,x.vec1-y.vec1); }
inline m256cd operator* (const m256cd &x, const m256cd &y) { return m256cd(x.vec0*y.vec0,x.vec1*y.vec1); }
inline m256cd operator* (const m256cd &x,        double y) { return m256cd(x.vec0*y     ,x.vec1*y     ); }
inline m256cd operator/ (const m256cd &x,        double y) { return m256cd(x.vec0/y     ,x.vec1/y     ); }
inline m256cd operator* (      double  x, const m256cd &y) { return y*x; }
inline m256cd&operator+=(      m256cd &x, const m256cd &y) { return x=x+y; }
inline m256cd&operator-=(      m256cd &x, const m256cd &y) { return x=x-y; }
inline m256cd&operator*=(      m256cd &x, const m256cd &y) { return x=x*y; }
inline m256cd&operator*=(      m256cd &x,        double y) { return x=x*y; }
inline m256cd&operator/=(      m256cd &x,        double y) { return x=x/y; }

inline m256cd addsub    (const m256cd &x, const m256cd &y) { return m256cd(addsub(x.vec0,y.vec0),addsub(x.vec1,y.vec1)); }
inline m256cd subadd    (const m256cd &x, const m256cd &y) { return m256cd(subadd(x.vec0,y.vec0),subadd(x.vec1,y.vec1)); }




//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 16 signed chars
class m128c
{
  public:
	  typedef m128c self;
	  typedef char value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=16 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128c()  {}
    inline m128c(__m128i m) { vec = m; }
    inline m128c(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7, value_type v8, value_type v9, value_type v10, value_type v11, value_type v12, value_type v13, value_type v14, value_type v15) { vec= _mm_set_epi8(v15,v14,v13,v12,v11,v10,v9,v8,v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128c(value_type v)	{ vec = _mm_set1_epi8(v); }
    inline m128c(const self &x) : vec(x.vec) { }

    inline m128c& operator=(value_type  v) { vec=_mm_set1_epi8(v); return *this; }
    inline m128c& operator=(const self &x) { vec=x.vec;            return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m128c &a){ char *p = (char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) << ',' << (int)*(p+8) << ',' << (int)*(p+9) << ',' << (int)*(p+10) << ',' << (int)*(p+11) << ',' << (int)*(p+12) << ',' << (int)*(p+13) << ',' << (int)*(p+14) << ',' << (int)*(p+15) <<')'; return os; }

inline m128c operator+ (const m128c &x, const m128c &y) { return m128c(_mm_add_epi8(x.vec,y.vec)); }
inline m128c operator- (const m128c &x, const m128c &y) { return m128c(_mm_sub_epi8(x.vec,y.vec)); }
inline m128c operator- (                const m128c &y) { return m128c(_mm_setzero_si128())-y; }
inline m128c&operator+=(      m128c &x, const m128c &y) { return x=x+y; }
inline m128c&operator-=(      m128c &x, const m128c &y) { return x=x-y; }

inline void load	  (m128c &x, const char *p) { x.vec = *(__m128i *)p; }
inline void store	  (char *p, const m128c &x) { *(__m128i *)p = x.vec; }

inline m128c unpacklo (const m128c &a, const m128c &b) { return m128c(_mm_unpacklo_epi8(a.vec,b.vec)); }
inline m128c unpackhi (const m128c &a, const m128c &b) { return m128c(_mm_unpackhi_epi8(a.vec,b.vec)); }
inline m128c flip     (const m128c &a) { char *p = (char *)&a.vec; return m128c(*(p+15),*(p+14),*(p+13),*(p+12),*(p+11),*(p+10),*(p+9),*(p+8),*(p+7),*(p+6),*(p+5),*(p+4),*(p+3),*(p+2),*(p+1),*p); }

template<int N> inline m128c shiftl   (const m128c &a) { return m128c(_mm_srli_si128(a.vec,N)); }
template<int N> inline m128c shiftr   (const m128c &a) { return m128c(_mm_slli_si128(a.vec,N)); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 16 signed chars with saturation
class m128cs
{
  public:
	  typedef m128cs self;
	  typedef char value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=16 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128cs()  {}
    inline m128cs(__m128i m) { vec = m; }
    inline m128cs(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7, value_type v8, value_type v9, value_type v10, value_type v11, value_type v12, value_type v13, value_type v14, value_type v15) { vec= _mm_set_epi8(v15,v14,v13,v12,v11,v10,v9,v8,v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128cs(value_type v)	{ vec = _mm_set1_epi8(v); }
    inline m128cs(const self &x) : vec(x.vec) { }

    inline m128cs& operator=(value_type  v) { vec=_mm_set1_epi8(v); return *this; }
    inline m128cs& operator=(const self &x) { vec=x.vec;            return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m128cs &a){ char *p = (char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) << ',' << (int)*(p+8) << ',' << (int)*(p+9) << ',' << (int)*(p+10) << ',' << (int)*(p+11) << ',' << (int)*(p+12) << ',' << (int)*(p+13) << ',' << (int)*(p+14) << ',' << (int)*(p+15) <<')'; return os; }

inline m128cs operator+ (const m128cs &x, const m128cs &y) { return m128cs(_mm_adds_epi8(x.vec,y.vec)); }
inline m128cs operator- (const m128cs &x, const m128cs &y) { return m128cs(_mm_subs_epi8(x.vec,y.vec)); }
inline m128cs operator- (                 const m128cs &y) { return m128cs(_mm_setzero_si128())-y; }
inline m128cs&operator+=(      m128cs &x, const m128cs &y) { return x=x+y; }
inline m128cs&operator-=(      m128cs &x, const m128cs &y) { return x=x-y; }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 16 unsigned chars
class m128b
{
  public:
	  typedef m128b self;
	  typedef unsigned char value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=16 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128b() {}
    inline m128b(__m128i m) { vec = m; }
    inline m128b(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7, value_type v8, value_type v9, value_type v10, value_type v11, value_type v12, value_type v13, value_type v14, value_type v15) { vec= _mm_set_epi8(v15,v14,v13,v12,v11,v10,v9,v8,v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128b(value_type v) { vec = _mm_set1_epi8(v); }
    inline m128b(const self &x) : vec(x.vec) {}

    inline operator __m128i() const { return vec; }

    inline m128b& operator=(const value_type &v) { vec=_mm_set1_epi8(v); return *this; }
    
   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };    
};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m128b &a) { unsigned char *p = (unsigned char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) << ',' << (int)*(p+8) << ',' << (int)*(p+9) << ',' << (int)*(p+10) << ',' << (int)*(p+11) << ',' << (int)*(p+12) << ',' << (int)*(p+13) << ',' << (int)*(p+14) << ',' << (int)*(p+15) <<')'; return os; }

inline m128b operator+ (const m128b &a, const m128b &b) { return m128b(_mm_add_epi8(a.vec,b.vec)); }
inline m128b operator- (const m128b &a, const m128b &b) { return m128b(_mm_sub_epi8(a.vec,b.vec)); }
inline m128b operator- (const m128b &a                ) { return m128b(_mm_sub_epi8(_mm_setzero_si128(),a.vec)); }
inline m128b&operator+=(      m128b &x, const m128b &y) { return x=x+y; }
inline m128b&operator-=(      m128b &x, const m128b &y) { return x=x-y; }

inline m128b min       (const m128b &a,const m128b &b) { return m128b(_mm_min_epu8   (a.vec,b.vec)); }
inline m128b max       (const m128b &a,const m128b &b) { return m128b(_mm_max_epu8   (a.vec,b.vec)); }

inline void load (m128b &x, const unsigned char *p) { x.vec = _mm_load_si128 ((__m128i *)p); }
#if !defined(SSE3)
inline void loadu(m128b &x, const unsigned char *p) { x.vec = _mm_loadu_si128((__m128i *)p); }
#endif

inline void store (unsigned char *p, const m128b &x) { _mm_store_si128 ((__m128i *)p, x.vec); }
inline void storeu(unsigned char *p, const m128b &x) { _mm_storeu_si128((__m128i *)p, x.vec); }

inline void loadl (m128b &x, const unsigned char *p) { __m128 m=(__m128&)x.vec; m=_mm_loadl_pi (m,(__m64 *)p); x.vec=(__m128i&)m; }
inline void loadh (m128b &x, const unsigned char *p) { __m128 m=(__m128&)x.vec; m=_mm_loadh_pi (m,(__m64 *)p); x.vec=(__m128i&)m; }
inline void loadlh(m128b &x, const unsigned char *p) { x.vec=_mm_loadl_epi64((const __m128i *)p); x.vec= _mm_unpacklo_epi64(x.vec,x.vec); }
inline void loadlh(m128b &x, const unsigned char *p1, const unsigned char *p2) { x.vec=_mm_loadl_epi64((const __m128i *)p1); loadh(x,p2); }

inline void storel(unsigned char *p, const m128b &x) { _mm_storel_pi((__m64 *)p,(__m128&)x.vec); }
inline void storeh(unsigned char *p, const m128b &x) { _mm_storeh_pi((__m64 *)p,(__m128&)x.vec); }

template<int N> inline m128b shiftl   (const m128b &a) { return m128b(_mm_srli_si128(a.vec,N)); }
template<int N> inline m128b shiftr   (const m128b &a) { return m128b(_mm_slli_si128(a.vec,N)); }

inline m128b unpacklo (const m128b &a, const m128b &b) { return m128b(_mm_unpacklo_epi8(a.vec,b.vec)); }
inline m128b unpackhi (const m128b &a, const m128b &b) { return m128b(_mm_unpackhi_epi8(a.vec,b.vec)); }
inline m128s unpacklo(const m128b &a);
inline m128s unpackhi(const m128b &a);



//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 16 unsigned chars with saturation
class m128bs
{
  public:
	  typedef m128bs self;
	  typedef unsigned char value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=16 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128bs() {}
    inline m128bs(__m128i m) { vec = m; }
    inline m128bs(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7, value_type v8, value_type v9, value_type v10, value_type v11, value_type v12, value_type v13, value_type v14, value_type v15) { vec= _mm_set_epi8(v15,v14,v13,v12,v11,v10,v9,v8,v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128bs(value_type v) { vec=_mm_set1_epi8(v); }
    inline m128bs(const self   &x) : vec(x.vec) { }

    inline m128bs& operator=(value_type  v) { vec=_mm_set1_epi8(v); return *this; }
    inline m128bs& operator=(const self &x) { vec=x.vec;            return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m128bs &a)	{ unsigned char *p = (unsigned char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) << ',' << (int)*(p+8) << ',' << (int)*(p+9) << ',' << (int)*(p+10) << ',' << (int)*(p+11) << ',' << (int)*(p+12) << ',' << (int)*(p+13) << ',' << (int)*(p+14) << ',' << (int)*(p+15) <<')'; return os; }

inline m128bs operator+ (const m128bs &a, const m128bs &b) { return m128bs(_mm_adds_epu8(a.vec,b.vec)); }
inline m128bs operator- (const m128bs &a, const m128bs &b) { return m128bs(_mm_subs_epu8(a.vec,b.vec)); }
inline m128bs operator- (const m128bs &a                 ) { return m128bs(_mm_subs_epu8(_mm_setzero_si128(),a.vec)); }
inline m128bs&operator+=(      m128bs &x, const m128bs &y) { return x=x+y; }
inline m128bs&operator-=(      m128bs &x, const m128bs &y) { return x=x-y; }

inline m128bs min       (const m128bs &a,const m128bs &b) { return m128bs(_mm_min_epu8   (a.vec,b.vec)); }
inline m128bs max       (const m128bs &a,const m128bs &b) { return m128bs(_mm_max_epu8   (a.vec,b.vec)); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 8 signed short integers
class m128s
{
  public:
	  typedef m128s self;
	  typedef short value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=8 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128s() {}
    inline m128s(__m128i m) { vec = m; }
    inline m128s(value_type v0,value_type v1, value_type v2, value_type v3) { vec= _mm_set_epi16(v3,v2,v1,v0,v3,v2,v1,v0); }
    inline m128s(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7) { vec= _mm_set_epi16(v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128s(value_type v) { vec=_mm_set1_epi16(v); }
    inline m128s(const self &x) : vec(x.vec) {}

    inline m128s &operator=(value_type  v) { vec=_mm_set1_epi16(v); return *this; }
    inline m128s &operator=(const self &x) { vec=x.vec;             return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> & operator<<(basic_ostream<E,Tr> &os,const m128s &a)	{ short *p = (short*)&a; os << '(' << *p     << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ',' << *(p+4) << ',' << *(p+5) << ',' << *(p+6) << ',' << *(p+7) <<')'; return os; }

inline m128s operator- (const m128s &a                ) { return m128s(_mm_sub_epi16(_mm_setzero_si128(),a.vec)); }
inline m128s operator+ (const m128s &a, const m128s &b) { return m128s(_mm_add_epi16(a.vec,b.vec)); }
inline m128s operator- (const m128s &a, const m128s &b) { return m128s(_mm_sub_epi16(a.vec,b.vec)); }
inline m128s operator* (const m128s &a, const m128s &b) { return m128s(_mm_mullo_epi16   (a.vec,b.vec)); }
inline m128s&operator+=(      m128s &x, const m128s &y) { return x=x+y; }
inline m128s&operator-=(      m128s &x, const m128s &y) { return x=x-y; }
inline m128s mullo    (const m128s &a, const m128s &b) { return m128s(_mm_mullo_epi16   (a.vec,b.vec)); }
inline m128s mulhi    (const m128s &a, const m128s &b) { return m128s(_mm_mulhi_epi16   (a.vec,b.vec)); }

inline m128s sad(const m128b &a,const m128b &b) { return m128s(_mm_sad_epu8(a.vec,b.vec)); }
inline m128s shufflelh (const m128s &a,const m128s &b) { return m128s(_mm_unpacklo_epi64(a.vec,_mm_shuffle_epi32(b.vec,14))); }

inline m128s operator==(const m128s &a,const m128s &b) { return m128s(_mm_cmpeq_epi16 (a.vec,b.vec)); }
inline m128s operator> (const m128s &a,const m128s &b) { return m128s(_mm_cmpgt_epi16 (a.vec,b.vec)); }
inline m128s operator< (const m128s &a,const m128s &b) { return m128s(_mm_cmplt_epi16 (a.vec,b.vec)); }
inline m128s operator& (const m128s &a,const m128s &b) { return m128s(_mm_and_si128   (a.vec,b.vec)); }
inline m128s operator| (const m128s &a,const m128s &b) { return m128s(_mm_or_si128    (a.vec,b.vec)); }
inline m128s operator^ (const m128s &a,const m128s &b) { return m128s(_mm_xor_si128   (a.vec,b.vec)); }

inline m128s and_not   (const m128s &a,const m128s &b) { return m128s(_mm_andnot_si128(a.vec,b.vec)); }
inline m128s min       (const m128s &a,const m128s &b) { return m128s(_mm_min_epi16   (a.vec,b.vec)); }
inline m128s max       (const m128s &a,const m128s &b) { return m128s(_mm_max_epi16   (a.vec,b.vec)); }
inline m128s abs       (const m128s &a) { return max(a,-a); }

template<int i> inline m128s shufflelo(const m128s &x) { return m128s(_mm_shufflelo_epi16(x.vec,i)); }
template<int i> inline m128s shufflehi(const m128s &x) { return m128s(_mm_shufflehi_epi16(x.vec,i)); }

inline m128s flip    (const m128s &a) { return m128s(_mm_shuffle_epi32( shufflelo<_MM_SHUFFLE(0,1,2,3)>(shufflehi<_MM_SHUFFLE(0,1,2,3)>(a)).vec ,_MM_SHUFFLE(1,0,3,2) )); }
inline m128s unpacklo(const m128b &a) { return m128s(_mm_unpacklo_epi8 (a.vec,_mm_setzero_si128())); }
inline m128s unpackhi(const m128b &a) { return m128s(_mm_unpackhi_epi8 (a.vec,_mm_setzero_si128())); }
inline m128s unpacklo(const m128s &a, const m128s &b) { return m128s(_mm_unpacklo_epi16(a.vec,b.vec)); }
inline m128s unpackhi(const m128s &a, const m128s &b) { return m128s(_mm_unpackhi_epi16(a.vec,b.vec)); }
inline m128c packs   (const m128s &a, const m128s &b) { return m128c (_mm_packs_epi16   (a.vec,b.vec)); }
inline m128b packus  (const m128s &a, const m128s &b) { return m128b(_mm_packus_epi16  (a.vec,b.vec)); }

inline void load  (m128s &x, const short *p) { x.vec = _mm_load_si128 ((__m128i *)p); }
inline void loadu	(m128s &x, const short *p) { x.vec = _mm_loadu_si128((__m128i *)p); }
inline void loadr	(m128s &x, const short *p) { load(x,p); x=flip(x); }
inline void loadl (m128s &x, const short *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadh (m128s &x, const short *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadlh(m128s &x, const short *p1, const short *p2) { loadl(x,p1); loadh(x,p2); }

inline void store (short *p, const m128s &x) { _mm_store_si128 ((__m128i *)p,x.vec); }
inline void storeu(short *p, const m128s &x) { _mm_storeu_si128((__m128i *)p,x.vec); }
inline void storer(short *p, const m128s &x) { store(p,flip(x)); }
inline void storel(short *p, const m128s &x) { _mm_storel_pd((double *)p, (__m128d &)x.vec); }
inline void storeh(short *p, const m128s &x) { _mm_storeh_pd((double *)p, (__m128d &)x.vec); }
inline void storelh(short *p1, short *p2, const m128s &x) { storel(p1,x); storeh(p2,x); }

template<int N> inline m128s shiftl (const m128s &x) { return m128s(_mm_srli_si128(x.vec,2*N)); }
template<int N> inline m128s shiftr (const m128s &x) { return m128s(_mm_slli_si128(x.vec,2*N)); }
template<int N> inline m128s shiftll(const m128s &x) { return m128s(_mm_slli_epi16(x.vec,N)); }
template<int N> inline m128s shiftrl(const m128s &x) { return m128s(_mm_srli_epi16(x.vec,N)); }
template<int N> inline m128s shiftla(const m128s &x) { return m128s(_mm_slli_epi16(x.vec,N)); }
template<int N> inline m128s shiftra(const m128s &x) { return m128s(_mm_srai_epi16(x.vec,N)); }

inline short get_sad(const m128s &a) { m128s b(_mm_shuffle_epi32(a.vec,78)); m128s c=a+b; return *(short *)&c; }
template<int N> inline short get_sad(const m128s &x) { return sad_result(shiftl<N>(x)); }

inline m128s sum(const m128s &x) { __m128i a=_mm_add_epi16(x.vec,_mm_shuffle_epi32(x.vec,14)); __m128i b=_mm_add_epi16(a,_mm_shufflelo_epi16(a,14)); return m128s(_mm_add_epi16(b,_mm_shufflelo_epi16(b,1))); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 8 signed short integers with saturation
class m128ss
{
  public:
	  typedef m128ss self;
	  typedef short value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=8 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128ss() {}
    inline m128ss(__m128i m) { vec = m; }
    inline m128ss(value_type v0,value_type v1, value_type v2, value_type v3) { vec= _mm_set_epi16(v3,v2,v1,v0,v3,v2,v1,v0); }
    inline m128ss(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7) { vec= _mm_set_epi16(v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128ss(value_type v) { vec=_mm_set1_epi16(v); }
    inline m128ss(const self &x) : vec(x.vec) {}

    inline m128ss &operator=(value_type  v) { vec=_mm_set1_epi16(v); return *this; }
    inline m128ss &operator=(const self &x) { vec=x.vec;             return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> & operator<<(basic_ostream<E,Tr> &os,const m128ss &a)	{ short *p = (short*)&a; os << '(' << *p     << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ',' << *(p+4) << ',' << *(p+5) << ',' << *(p+6) << ',' << *(p+7) <<')'; return os; }

template<int i> inline m128ss shufflelo(const m128ss &x) { return m128ss(_mm_shufflelo_epi16(x.vec,i)); }
template<int i> inline m128ss shufflehi(const m128ss &x) { return m128ss(_mm_shufflehi_epi16(x.vec,i)); }

inline m128ss flip    (const m128ss &x) { return m128ss(_mm_shuffle_epi32( shufflelo<_MM_SHUFFLE(0,1,2,3)>(shufflehi<_MM_SHUFFLE(0,1,2,3)>(x)).vec ,_MM_SHUFFLE(1,0,3,2) )); }

inline void load  (m128ss &x, const short *p) { x.vec = _mm_load_si128 ((__m128i *)p); }
inline void loadu	(m128ss &x, const short *p) { x.vec = _mm_loadu_si128((__m128i *)p); }
inline void loadr	(m128ss &x, const short *p) { load(x,p); x=flip(x); }
inline void loadl (m128ss &x, const short *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadh (m128ss &x, const short *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadlh(m128ss &x, const short *p1, const short *p2) { loadl(x,p1); loadh(x,p2); }

inline void store (short *p, const m128ss &x) { _mm_store_si128 ((__m128i *)p,x.vec); }
inline void storeu(short *p, const m128ss &x) { _mm_storeu_si128((__m128i *)p,x.vec); }
inline void storer(short *p, const m128ss &x) { store(p,flip(x)); }
inline void storel(short *p, const m128ss &x) { _mm_storel_pd((double *)p, (__m128d &)x.vec); }
inline void storeh(short *p, const m128ss &x) { _mm_storeh_pd((double *)p, (__m128d &)x.vec); }
inline void storelh(short *p1, short *p2, const m128ss &x) { storel(p1,x); storeh(p2,x); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 8 unsigned short integers
class m128w
{
  public:
	  typedef m128w self;
	  typedef unsigned short value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=8 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128w() {}
    inline m128w(__m128i m) { vec = m; }
    inline m128w(value_type v0,value_type v1, value_type v2, value_type v3) { vec= _mm_set_epi16(v3,v2,v1,v0,v3,v2,v1,v0); }
    inline m128w(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7) { vec= _mm_set_epi16(v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128w(value_type v) { vec=_mm_set1_epi16(v); }
    inline m128w(const self &x) : vec(x.vec) {}

    inline m128w &operator=(value_type  v) { vec=_mm_set1_epi16(v); return *this; }
    inline m128w &operator=(const self &x) { vec=x.vec;             return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> & operator<<(basic_ostream<E,Tr> &os,const m128w &a)	{ unsigned short *p = (unsigned short*)&a; os << '(' << *p     << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ',' << *(p+4) << ',' << *(p+5) << ',' << *(p+6) << ',' << *(p+7) <<')'; return os; }

inline m128w operator+ (const m128w &x, const m128w &y) { return m128w(_mm_add_epi16(x.vec,y.vec)); }
inline m128w operator- (const m128w &x, const m128w &y) { return m128w(_mm_sub_epi16(x.vec,y.vec)); }
inline m128w operator- (                const m128w &y) { return m128w(_mm_setzero_si128())-y; }
inline m128w operator* (const m128w &x, const m128w &y) { return m128w(_mm_mullo_epi16(x.vec,y.vec)); }
inline m128w&operator+=(      m128w &x, const m128w &y) { return x=x+y; }
inline m128w&operator-=(      m128w &x, const m128w &y) { return x=x-y; }

template<int i> inline m128w shufflelo(const m128w &x) { return m128w(_mm_shufflelo_epi16(x.vec,i)); }
template<int i> inline m128w shufflehi(const m128w &x) { return m128w(_mm_shufflehi_epi16(x.vec,i)); }

inline m128w flip    (const m128w &x) { return m128w(_mm_shuffle_epi32( shufflelo<_MM_SHUFFLE(0,1,2,3)>(shufflehi<_MM_SHUFFLE(0,1,2,3)>(x)).vec ,_MM_SHUFFLE(1,0,3,2) )); }

inline void load  (m128w &x, const unsigned short *p) { x.vec = _mm_load_si128 ((__m128i *)p); }
inline void loadu	(m128w &x, const unsigned short *p) { x.vec = _mm_loadu_si128((__m128i *)p); }
inline void loadr	(m128w &x, const unsigned short *p) { load(x,p); x=flip(x); }
inline void loadl (m128w &x, const unsigned short *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadh (m128w &x, const unsigned short *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadlh(m128w &x, const unsigned short *p1, const unsigned short *p2) { loadl(x,p1); loadh(x,p2); }

inline void store (unsigned short *p, const m128w &x) { _mm_store_si128 ((__m128i *)p,x.vec); }
inline void storeu(unsigned short *p, const m128w &x) { _mm_storeu_si128((__m128i *)p,x.vec); }
inline void storer(unsigned short *p, const m128w &x) { store(p,flip(x)); }
inline void storel(unsigned short *p, const m128w &x) { _mm_storel_pd((double *)p, (__m128d &)x.vec); }
inline void storeh(unsigned short *p, const m128w &x) { _mm_storeh_pd((double *)p, (__m128d &)x.vec); }
inline void storelh(unsigned short *p1, unsigned short *p2, const m128w &x) { storel(p1,x); storeh(p2,x); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 8 unsigned short integers with saturation
class m128ws
{
  public:
	  typedef m128ws self;
	  typedef unsigned short value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=8 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128ws() {}
    inline m128ws(__m128i m) { vec = m; }
    inline m128ws(value_type v0,value_type v1, value_type v2, value_type v3) { vec= _mm_set_epi16(v3,v2,v1,v0,v3,v2,v1,v0); }
    inline m128ws(value_type v0,value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7) { vec= _mm_set_epi16(v7,v6,v5,v4,v3,v2,v1,v0); }
    inline explicit m128ws(value_type v) { vec=_mm_set1_epi16(v); }
    inline m128ws(const self &x) : vec(x.vec) {}

    inline m128ws &operator=(value_type  v) { vec=_mm_set1_epi16(v); return *this; }
    inline m128ws &operator=(const self &x) { vec=x.vec;             return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class Tr> basic_ostream<E,Tr> & operator<<(basic_ostream<E,Tr> &os,const m128ws &a)	{ unsigned short *p = (unsigned short*)&a; os << '(' << *p     << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ',' << *(p+4) << ',' << *(p+5) << ',' << *(p+6) << ',' << *(p+7) <<')'; return os; }

inline m128ws operator+ (const m128ws &x, const m128ws &y) { return m128ws(_mm_adds_epu16(x.vec,y.vec)); }
inline m128ws operator- (const m128ws &x, const m128ws &y) { return m128ws(_mm_subs_epu16(x.vec,y.vec)); }
inline m128ws operator- (                 const m128ws &y) { return m128ws(_mm_setzero_si128())-y; }
inline m128ws operator* (const m128ws &x, const m128ws &y) { return m128ws(_mm_mullo_epi16(x.vec,y.vec)); }
inline m128ws&operator+=(      m128ws &x, const m128ws &y) { return x=x+y; }
inline m128ws&operator-=(      m128ws &x, const m128ws &y) { return x=x-y; }

template<int i> inline m128ws shufflelo(const m128ws &x) { return m128ws(_mm_shufflelo_epi16(x.vec,i)); }
template<int i> inline m128ws shufflehi(const m128ws &x) { return m128ws(_mm_shufflehi_epi16(x.vec,i)); }

inline m128ws flip    (const m128ws &x) { return m128ws(_mm_shuffle_epi32( shufflelo<_MM_SHUFFLE(0,1,2,3)>(shufflehi<_MM_SHUFFLE(0,1,2,3)>(x)).vec ,_MM_SHUFFLE(1,0,3,2) )); }

inline void load  (m128ws &x, const unsigned short *p) { x.vec = _mm_load_si128 ((__m128i *)p); }
inline void loadu	(m128ws &x, const unsigned short *p) { x.vec = _mm_loadu_si128((__m128i *)p); }
inline void loadr	(m128ws &x, const unsigned short *p) { load(x,p); x=flip(x); }
inline void loadl (m128ws &x, const unsigned short *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadh (m128ws &x, const unsigned short *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadlh(m128ws &x, const unsigned short *p1, const unsigned short *p2) { loadl(x,p1); loadh(x,p2); }

inline void store (unsigned short *p, const m128ws &x) { _mm_store_si128 ((__m128i *)p,x.vec); }
inline void storeu(unsigned short *p, const m128ws &x) { _mm_storeu_si128((__m128i *)p,x.vec); }
inline void storer(unsigned short *p, const m128ws &x) { store(p,flip(x)); }
inline void storel(unsigned short *p, const m128ws &x) { _mm_storel_pd((double *)p, (__m128d &)x.vec); }
inline void storeh(unsigned short *p, const m128ws &x) { _mm_storeh_pd((double *)p, (__m128d &)x.vec); }
inline void storelh(unsigned short *p1, unsigned short *p2, const m128ws &x) { storel(p1,x); storeh(p2,x); }


//{unsecret}
//{group:SSE2}
//Summary: SSE2 register repesenting 4 signed integers
class m128i
{
  public:
	  typedef m128i self;
	  typedef int value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=4 };

  public:
	  __aligned __m128i vec;

  public:
    inline m128i() {}
    inline m128i(__m128i m) { vec = m; }
    inline m128i(value_type v0,value_type v1) { vec= _mm_set_epi32(v1,v0,v1,v0); }
    inline m128i(value_type v0,value_type v1,value_type v2,value_type v3) { vec= _mm_set_epi32(v3,v2,v1,v0); }
    inline explicit m128i(value_type v) { vec = _mm_set1_epi32(v); }
    inline m128i(const self &x) : vec(x.vec) { }

    inline self &operator=(value_type  v) { vec=_mm_set1_epi32(v); return *this; }
    inline self &operator=(const self &x) { vec=x.vec;             return *this; }

    inline operator __m128i() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};
template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> & os,const m128i &a) { int *p=(int*)&a; os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3)  <<')'; return os; }

inline m128i operator- (const m128i &a                ) { return m128i(_mm_sub_epi32(_mm_setzero_si128(),a.vec)); }
inline m128i operator+ (const m128i &a, const m128i &b) { return m128i(_mm_add_epi32(a.vec,b.vec)); }
inline m128i operator- (const m128i &a, const m128i &b) { return m128i(_mm_sub_epi32(a.vec,b.vec)); }
inline m128i&operator+=(      m128i &x, const m128i &y) { return x=x+y; }
inline m128i&operator-=(      m128i &x, const m128i &y) { return x=x-y; }
inline m128i madd      (const m128s &x, const m128s &y) { return m128i(_mm_madd_epi16(x.vec,y.vec)); }

template<int i> inline m128i shuffle(const m128i &x) { return m128i(_mm_shuffle_epi32(x.vec,i)); }

inline m128i flip     (const m128i &a) { return shuffle<_MM_SHUFFLE(0,1,2,3)>(a); }
inline m128i unpacklo (const m128s &a) { return m128i(_mm_unpacklo_epi16 (a.vec,_mm_setzero_si128())); }
inline m128i unpackhi (const m128s &a) { return m128i(_mm_unpackhi_epi16 (a.vec,_mm_setzero_si128())); }
inline m128i unpacklo (const m128i &a, const m128i &b) { return m128i(_mm_unpacklo_epi32(a.vec,b.vec)); }
inline m128i unpackhi (const m128i &a, const m128i &b) { return m128i(_mm_unpackhi_epi32(a.vec,b.vec)); }
inline m128s packs    (const m128i &a, const m128i &b) { return m128s(_mm_packs_epi32  (a.vec,b.vec)); }
inline m128i movehl   (const m128i &a, const m128i &b) { return m128i( _mm_unpackhi_epi64( a.vec, b.vec) ); }
inline m128i movelh   (const m128i &a, const m128i &b) { return m128i( _mm_unpacklo_epi64( a.vec, b.vec) ); }

inline void load  (m128i &x, const int *p) { x.vec = _mm_load_si128 ((__m128i *)p); }
inline void loadu	(m128i &x, const int *p) { x.vec = _mm_loadu_si128((__m128i *)p); }
inline void loadr	(m128i &x, const int *p) { load(x,p); x=flip(x); }
inline void loadl (m128i &x, const int *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadh (m128i &x, const int *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128i &)a; }
inline void loadlh(m128i &x, const int *p1, const int *p2) { loadl(x,p1); loadh(x,p2); }

inline void store (int *p, const m128i &x) { _mm_store_si128 ((__m128i *)p,x.vec); }
inline void storeu(int *p, const m128i &x) { _mm_storeu_si128((__m128i *)p,x.vec); }
inline void storer(int *p, const m128i &x) { store(p,flip(x)); }
inline void storel(int *p, const m128i &x) { _mm_storel_pd((double *)p, (__m128d &)x.vec); }
inline void storeh(int *p, const m128i &x) { _mm_storeh_pd((double *)p, (__m128d &)x.vec); }
inline void storelh(int *p0, int *p1, const m128i &x) { storel(p0,x); storeh(p1,x); }

template<int N> inline m128i shiftl (const m128i &a) { return m128i(_mm_srli_si128(a.vec,4*N)); }
template<int N> inline m128i shiftr (const m128i &a) { return m128i(_mm_slli_si128(a.vec,4*N)); }
template<int N> inline m128i shiftll(const m128i &a) { return m128i(_mm_slli_epi32(a.vec,N)); } 
template<int N> inline m128i shiftrl(const m128i &a) { return m128i(_mm_srli_epi32(a.vec,N)); }
template<int N> inline m128i shiftla(const m128i &a) { return m128i(_mm_slli_epi32(a.vec,N)); } 
template<int N> inline m128i shiftra(const m128i &a) { return m128i(_mm_srai_epi32(a.vec,N)); }

inline int get0(const m128i &x) { return *(int*)&(x.vec); }
inline int get1(const m128i &x) { m128i a=shuffle<_MM_SHUFFLE(3,1,3,1)>(x); return *(int*)&(a.vec); }
inline int get2(const m128i &x) { m128i a=movehl(x,x); return *(int*)&(a.vec); }
inline int get3(const m128i &x) { m128i a=shuffle<_MM_SHUFFLE(2,3,2,3)>(x); return *(int*)&(a.vec); }

inline m128i h2add(const m128i &x, const m128i &y) { return movelh(x,y)+movehl(x,y); }
inline m128i h2sub(const m128i &x, const m128i &y) { return movelh(x,y)-movehl(y,x); }
inline m128i h3add(const m128i &x, const m128i &y) { return unpacklo(x,y)+unpackhi(x,y); }
inline m128i h3sub(const m128i &x, const m128i &y) { return unpacklo(x,y)-unpackhi(x,y); }
inline m128i hadd (const m128i &x, const m128i &y) { return h2add(shuffle<_MM_SHUFFLE(3,1,2,0)>(x),shuffle<_MM_SHUFFLE(3,1,2,0)>(y)); }
inline m128i hsub (const m128i &x, const m128i &y) { return h2sub(shuffle<_MM_SHUFFLE(3,1,2,0)>(x),shuffle<_MM_SHUFFLE(3,1,2,0)>(y)); }

inline m128i sum  (const m128i &x0) { m128i a=h3add(x0,x0); return h2add(a,a); }
inline m128i sum  (const m128i &x0, const m128i &x1) { m128i a=h3add(x0,x1); return h2add(a,a); }
inline m128i sum  (const m128i &x0, const m128i &x1, const m128i &x2, const m128i &x3) { return h2add(h3add(x0,x1),h3add(x2,x3)); }
inline m128i sumlo(const m128i &x0) { return x0           +shuffle<_MM_SHUFFLE(3,1,3,1)>(x0); }
inline m128i sumhi(const m128i &x0) { return movehl(x0,x0)+shuffle<_MM_SHUFFLE(2,3,2,3)>(x0); }


#endif // NO_SIMD

#endif

