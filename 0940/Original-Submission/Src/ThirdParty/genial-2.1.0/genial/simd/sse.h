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

#ifndef SSE_H
#define SSE_H

#ifndef NO_SIMD

#ifndef SSE
#define SSE
#endif

#include <xmmintrin.h>


#include "memory.h"
//#include "simd/mmx.h"
#include <complex>

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

using namespace std;

class m64f;
class m128f;
class m128cf;


//{unsecret}
//{group:SSE}
//Summary: SSE register repesenting 4 floats
class m64f
{
  public:
    typedef m64f self;
    typedef float value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;
    
    enum { N=2 };

  public:
    __aligned __m128 vec;

  public:
    inline m64f(        ) {}
    explicit inline m64f(__m128 m) : vec(m) {}
    inline explicit m64f(value_type x) : vec(_mm_set_ps1(x))	{}
    inline m64f(value_type x0, value_type x1) : vec(_mm_set_ps(x1,x0,x1,x0)) {}
    inline m64f(const self &x) : vec(x.vec) {}
    
    inline m64f &operator=(value_type  x) { vec=_mm_set_ps1(x); return *this; }
    inline m64f &operator=(const self &x) { vec=x.vec;          return *this; }

    inline operator __m128() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m64f &x)	{ float *p = (float*)&x; return os << '(' << *p << ',' << *(p+1) << ')'; }

inline m64f operator- (               const m64f &y) { return m64f(_mm_sub_ps(_mm_setzero_ps(),y.vec)); }
inline m64f operator+ (const m64f &x, const m64f &y) { return m64f(_mm_add_ps(x.vec,y.vec)); }
inline m64f operator- (const m64f &x, const m64f &y) { return m64f(_mm_sub_ps(x.vec,y.vec)); }
inline m64f operator* (const m64f &x, const m64f &y) { return m64f(_mm_mul_ps(x.vec,y.vec)); }
inline m64f operator/ (const m64f &x, const m64f &y) { return m64f(_mm_div_ps(x.vec,y.vec)); }
inline m64f operator* (const m64f &x, const float y) { return m64f(_mm_mul_ps(x.vec,_mm_set_ps1(y))); }
inline m64f operator/ (const m64f &x, const float y) { return m64f(_mm_div_ps(x.vec,_mm_set_ps1(y))); }
inline m64f operator* (const float x, const m64f &y) { return m64f(_mm_mul_ps(_mm_set_ps1(x),y.vec)); }
inline m64f operator/ (const float x, const m64f &y) { return m64f(_mm_div_ps(_mm_set_ps1(x),y.vec)); }
inline m64f inv       (const m64f &x               ) { return float(1)/x; }
inline m64f operator& (const m64f &x, const m64f &y) { return m64f(_mm_and_ps   (x.vec,y.vec)); }
inline m64f operator| (const m64f &x, const m64f &y) { return m64f(_mm_or_ps    (x.vec,y.vec)); }
inline m64f operator^ (const m64f &x, const m64f &y) { return m64f(_mm_xor_ps   (x.vec,y.vec)); }
inline m64f and_not   (const m64f &x, const m64f &y) { return m64f(_mm_andnot_ps(x.vec,y.vec)); }
inline m64f&operator+=(      m64f &x, const m64f &y) { return x=x+y; }
inline m64f&operator-=(      m64f &x, const m64f &y) { return x=x-y; }
inline m64f&operator*=(      m64f &x, const m64f &y) { return x=x*y; }
inline m64f&operator/=(      m64f &x, const m64f &y) { return x=x/y; }
inline m64f&operator*=(      m64f &x, const float y) { return x=x*y; }
inline m64f&operator/=(      m64f &x, const float y) { return x=x/y; }
inline m64f&operator&=(      m64f &x, const m64f &y) { return x=x&y; }
inline m64f&operator|=(      m64f &x, const m64f &y) { return x=x|y; }
inline m64f&operator^=(      m64f &x, const m64f &y) { return x=x^y; }

inline m64f operator==(const m64f &x, const m64f &y) { return m64f(_mm_cmpeq_ps (x.vec,y.vec)); }
inline m64f operator!=(const m64f &x, const m64f &y) { return m64f(_mm_cmpneq_ps(x.vec,y.vec)); }
inline m64f operator> (const m64f &x, const m64f &y) { return m64f(_mm_cmpgt_ps (x.vec,y.vec)); }
inline m64f operator>=(const m64f &x, const m64f &y) { return m64f(_mm_cmpge_ps (x.vec,y.vec)); }
inline m64f operator< (const m64f &x, const m64f &y) { return m64f(_mm_cmplt_ps (x.vec,y.vec)); }
inline m64f operator<=(const m64f &x, const m64f &y) { return m64f(_mm_cmple_ps (x.vec,y.vec)); }

inline m64f min       (const m64f &x, const m64f &y) { return m64f(_mm_min_ps(x.vec,y.vec)); }
inline m64f max       (const m64f &x, const m64f &y) { return m64f(_mm_max_ps(x.vec,y.vec)); }

inline m64f sqrt (const m64f &x) { return m64f(_mm_sqrt_ps (x.vec)); }
inline m64f rcp  (const m64f &x) { return m64f(_mm_rcp_ps  (x.vec)); }
inline m64f rsqrt(const m64f &x) { return m64f(_mm_rsqrt_ps(x.vec)); }

template<int i                      > inline m64f shuffle(const m64f &x               ) { return m64f(_mm_shuffle_ps(x.vec,x.vec,i)); }
template<int i0,int i1,int i2,int i3> inline m64f shuffle(const m64f &x               ) { return shuffle<_MM_SHUFFLE(i0,i1,i2,i3)>(x); }
template<int i                      > inline m64f shuffle(const m64f &x, const m64f &y) { return m64f(_mm_shuffle_ps(x.vec,y.vec,i)); }
template<int i0,int i1,int i2,int i3> inline m64f shuffle(const m64f &x, const m64f &y) { return shuffle<_MM_SHUFFLE(i0,i1,i2,i3)>(x,y); }

inline m64f alternate(const m64f &x) { return m64f(1,-1)*x; }
inline m64f flip   (const m64f &x) { return shuffle<2,3,0,1>(x); }
inline m64f flip_ri(const m64f &x) { return shuffle<2,3,0,1>(x); }

inline m64f conj   (const m64f &x) { return x; }
inline m64f real   (const m64f &x) { return x; }
inline m64f imag   (const m64f &x) { return m64f(_mm_setzero_ps()); }

inline float get0(const m64f &x) { return *(float*)&(x.vec); }
inline float get1(const m64f &x) { m64f a=shuffle<_MM_SHUFFLE(3,1,3,1)>(x); return *(float*)&(a.vec); }

inline void load0 (m64f &x, const float *p) { x[0]=*p; }
inline void load1 (m64f &x, const float *p) { x[1]=*p; }
inline void loadl (m64f &x, const float *p) { load0(x,p); }
inline void loadh (m64f &x, const float *p) { load1(x,p); }

inline void store0(float *p, const m64f &x) { *p=x[0]; }
inline void store1(float *p, const m64f &x) { *p=x[1]; }
inline void storel(float *p, const m64f &x) { store0(p,x); }
inline void storeh(float *p, const m64f &x) { store1(p,x); }
#ifdef SSE2
inline void load   (m64f &x, const float *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128 &)a; }
inline void loadu  (m64f &x, const float *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128 &)a; }
inline void store  (float *p, const m64f &x) { _mm_storel_pd((double *)p,(__m128d &)x.vec); }
inline void storeu (float *p, const m64f &x) { _mm_storel_pd((double *)p,(__m128d &)x.vec); }
#else
inline void load   (m64f &x, const float *p) { x.vec=_mm_loadl_pi (x.vec,(__m64 *)p); }
inline void loadu  (m64f &x, const float *p) { x.vec=_mm_loadl_pi (x.vec,(__m64 *)p); }
inline void store  (float *p, m64f &x) { _mm_storel_pi (x.vec,(__m64 *)p); }
inline void storeu (float *p, m64f &x) { _mm_storel_pi (x.vec,(__m64 *)p); }
#endif
inline void loadr  (m64f &x, const float *p) { load(x,p); x=flip(x); }
inline void loadur (m64f &x, const float *p) { load(x,p); x=flip(x); }
inline void loadlh (m64f &x, const float *p1, const float *p2) { loadl (x,p1); loadh (x,p2); }
inline void storer (float *p,const m64f &x) { store (p,flip(x)); }
inline void storeur(float *p,const m64f &x) { storeu(p,flip(x)); }
inline void storelh(float *p1, float *p2, const m64f &x) { storel(p1,x); storeh(p2,x); }


//{unsecret}
//{group:SSE}
//Summary: SSE register repesenting 4 floats
class m128f
{
  public:
    typedef m128f self;
    typedef float value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;
    
    enum { N=4 };

  public:
    __aligned __m128 vec;

  public:
    inline m128f(        ) {}
    explicit inline m128f(__m128 m) : vec(m) {}
    inline explicit m128f(value_type x) : vec(_mm_set_ps1(x))	{}
    inline m128f(value_type x0, value_type x1, value_type x2, value_type x3) : vec(_mm_set_ps(x3,x2,x1,x0)) {}
    inline m128f(const self &x) : vec(x.vec) {}
    
    inline m128f &operator=(value_type  x) { vec=_mm_set_ps1(x); return *this; }
    inline m128f &operator=(const self &x) { vec=x.vec;          return *this; }

    inline operator __m128() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m128f &x)	{ float *p = (float*)&x; return os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ')'; }

inline m128f operator- (                const m128f &y) { return m128f(_mm_sub_ps(_mm_setzero_ps(),y.vec)); }
inline m128f operator+ (const m128f &x, const m128f &y) { return m128f(_mm_add_ps(x.vec,y.vec)); }
inline m128f operator- (const m128f &x, const m128f &y) { return m128f(_mm_sub_ps(x.vec,y.vec)); }
inline m128f operator* (const m128f &x, const m128f &y) { return m128f(_mm_mul_ps(x.vec,y.vec)); }
inline m128f operator/ (const m128f &x, const m128f &y) { return m128f(_mm_div_ps(x.vec,y.vec)); }
inline m128f operator* (const m128f &x, const float  y) { return m128f(_mm_mul_ps(x.vec,_mm_set_ps1(y))); }
inline m128f operator/ (const m128f &x, const float  y) { return m128f(_mm_div_ps(x.vec,_mm_set_ps1(y))); }
inline m128f operator* (const float  x, const m128f &y) { return m128f(_mm_mul_ps(_mm_set_ps1(x),y.vec)); }
inline m128f operator/ (const float  x, const m128f &y) { return m128f(_mm_div_ps(_mm_set_ps1(x),y.vec)); }
inline m128f inv       (                const m128f &y) { return float(1)/y; }
inline m128f operator& (const m128f &x, const m128f &y) { return m128f(_mm_and_ps   (x.vec,y.vec)); }
inline m128f operator| (const m128f &x, const m128f &y) { return m128f(_mm_or_ps    (x.vec,y.vec)); }
inline m128f operator^ (const m128f &x, const m128f &y) { return m128f(_mm_xor_ps   (x.vec,y.vec)); }
inline m128f and_not   (const m128f &x, const m128f &y) { return m128f(_mm_andnot_ps(x.vec,y.vec)); }
inline m128f&operator+=(      m128f &x, const m128f &y) { return x=x+y; }
inline m128f&operator-=(      m128f &x, const m128f &y) { return x=x-y; }
inline m128f&operator*=(      m128f &x, const m128f &y) { return x=x*y; }
inline m128f&operator/=(      m128f &x, const m128f &y) { return x=x/y; }
inline m128f&operator*=(      m128f &x, const float  y) { return x=x*y; }
inline m128f&operator/=(      m128f &x, const float  y) { return x=x/y; }
inline m128f&operator&=(      m128f &x, const m128f &y) { return x=x&y; }
inline m128f&operator|=(      m128f &x, const m128f &y) { return x=x|y; }
inline m128f&operator^=(      m128f &x, const m128f &y) { return x=x^y; }

inline m128f operator==(const m128f &x, const m128f &y) { return m128f(_mm_cmpeq_ps (x.vec,y.vec)); }
inline m128f operator!=(const m128f &x, const m128f &y) { return m128f(_mm_cmpneq_ps(x.vec,y.vec)); }
inline m128f operator> (const m128f &x, const m128f &y) { return m128f(_mm_cmpgt_ps (x.vec,y.vec)); }
inline m128f operator>=(const m128f &x, const m128f &y) { return m128f(_mm_cmpge_ps (x.vec,y.vec)); }
inline m128f operator< (const m128f &x, const m128f &y) { return m128f(_mm_cmplt_ps (x.vec,y.vec)); }
inline m128f operator<=(const m128f &x, const m128f &y) { return m128f(_mm_cmple_ps (x.vec,y.vec)); }

inline m128f min       (const m128f &x, const m128f &y) { return m128f(_mm_min_ps(x.vec,y.vec)); }
inline m128f max       (const m128f &x, const m128f &y) { return m128f(_mm_max_ps(x.vec,y.vec)); }

inline m128f sqrt (const m128f &x) { return m128f(_mm_sqrt_ps (x.vec)); }
inline m128f rcp  (const m128f &x) { return m128f(_mm_rcp_ps  (x.vec)); }
inline m128f rsqrt(const m128f &x) { return m128f(_mm_rsqrt_ps(x.vec)); }

template<int i                      > inline m128f shuffle(const m128f &x                ) { return m128f(_mm_shuffle_ps(x.vec,x.vec,i)); }
template<int i0,int i1,int i2,int i3> inline m128f shuffle(const m128f &x                ) { return shuffle<_MM_SHUFFLE(i0,i1,i2,i3)>(x); }
template<int i                      > inline m128f shuffle(const m128f &x, const m128f &y) { return m128f(_mm_shuffle_ps(x.vec,y.vec,i)); }
template<int i0,int i1,int i2,int i3> inline m128f shuffle(const m128f &x, const m128f &y) { return shuffle<_MM_SHUFFLE(i0,i1,i2,i3)>(x,y); }

inline m128f alternate(const m128f &x) { return m128f(1,-1,1,-1)*x; }
inline m128f flip   (const m128f &x) { return shuffle<0,1,2,3>(x); }
inline m128f flip_ri(const m128f &x) { return shuffle<2,3,0,1>(x); }

inline m128f conj   (const m128f &x) { return x; }
inline m128f real   (const m128f &x) { return x; }
inline m128f imag   (const m128f &x) { return m128f(_mm_setzero_ps()); }

inline m128f movehl   (const m128f &x, const m128f &y) { return m128f(_mm_movehl_ps(x.vec,y.vec)); }
inline m128f movelh   (const m128f &x, const m128f &y) { return m128f(_mm_movelh_ps(x.vec,y.vec)); }
inline m128f unpacklo (const m128f &x, const m128f &y) { return m128f(_mm_unpacklo_ps(x.vec,y.vec)); }
inline m128f unpackhi (const m128f &x, const m128f &y) { return m128f(_mm_unpackhi_ps(x.vec,y.vec)); }
inline m128f packlo   (const m128f &x, const m128f &y) { return shuffle<2,0,2,0>(x,y); }
inline m128f packhi   (const m128f &x, const m128f &y) { return shuffle<3,1,3,1>(x,y); }

inline void trn    (m128f &x0, m128f &x1, m128f &x2,m128f &x3) { m128f a=movelh(x0,x1); m128f b=movehl(x1,x0); m128f c=movelh(x2,x3); m128f d=movehl(x3,x2); x0=shuffle<0x88>(a,c); x1=shuffle<0xDD>(a, c); x2=shuffle<0x88>(b, d); x3=shuffle<0xDD>(b,d); }
inline void up_trn (m128f &x0, m128f &x1, m128f &x2,m128f &x3) { m128f a=movelh(x0,x1); m128f b=movehl(x1,x0);                        m128f d=movehl(x3,x2);                        x1=shuffle<0xED>(a,x1); x2=shuffle<0xE8>(b,x2); x3=shuffle<0xDD>(b,d); }
inline void low_trn(m128f &x0, m128f &x1, m128f &x2,m128f &x3) { m128f a=movelh(x0,x1); m128f c=movelh(x2,x3);                        m128f d=movehl(x3,x2); x0=shuffle<0x88>(a,c); x1=shuffle<0xD4>(x1,c); x2=shuffle<0x84>(x2,d);                        }


inline float get0(const m128f &x) { return *(float*)&(x.vec); }
inline float get1(const m128f &x) { m128f a=shuffle<_MM_SHUFFLE(3,1,3,1)>(x,x); return get0(a); }
inline float get2(const m128f &x) { m128f a=movehl(x,x); return get0(a); }
inline float get3(const m128f &x) { m128f a=shuffle<_MM_SHUFFLE(2,3,2,3)>(x,x); return get0(a); }

inline void load0 (m128f &x, const float *p) { x[0]=*p; }
inline void load1 (m128f &x, const float *p) { x[1]=*p; }
inline void load2 (m128f &x, const float *p) { x[2]=*p; }
inline void load3 (m128f &x, const float *p) { x[3]=*p; }
inline void store0(float *p, const m128f &x) { _mm_store_ss (p, x.vec); }
inline void store1(float *p, const m128f &x) { *p = x[1]; }
inline void store2(float *p, const m128f &x) { *p = x[2]; }
inline void store3(float *p, const m128f &x) { *p = x[3]; }

#ifdef SSE2
inline void loadl (m128f &x, const float *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128 &)a; }
inline void loadh (m128f &x, const float *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128 &)a; }
inline void storel(float *p, const m128f &x) { _mm_storel_pd((double *)p, (__m128d &)x.vec); }
inline void storeh(float *p, const m128f &x) { _mm_storeh_pd((double *)p, (__m128d &)x.vec); }
#else
inline void loadl (m128f &x, const float *p) { x.vec=_mm_loadl_pi (x.vec,(__m64 *)p); }
inline void loadh (m128f &x, const float *p) { x.vec=_mm_loadh_pi (x.vec,(__m64 *)p); }
inline void storel(float *p, const m128f &x) { _mm_storel_pi((__m64 *)p,x.vec); }
inline void storeh(float *p, const m128f &x) { _mm_storeh_pi((__m64 *)p,x.vec); }
#endif

inline void load   (m128f &x, const float *p) { x.vec = _mm_load_ps ((float *)p); }
inline void loadu  (m128f &x, const float *p) { x.vec = _mm_loadu_ps((float *)p); }
inline void loadr  (m128f &x, const float *p) { x.vec = _mm_loadr_ps((float *)p); }
inline void loadur (m128f &x, const float *p) { loadu(x,p); x=flip(x); }
inline void loadlh (m128f &x, const float *p1, const float *p2) { loadl(x,p1); loadh(x,p2); }
inline void store  (float *p, const m128f &x) { _mm_store_ps (p, x.vec); }
inline void storeu (float *p, const m128f &x) { _mm_storeu_ps(p, x.vec); }
inline void storer (float *p, const m128f &x) { _mm_storer_ps(p, x.vec); }
inline void storeur(float *p, const m128f &x) { _mm_storeu_ps(p, flip(x).vec); }
inline void stream (float *p, const m128f &x) { _mm_stream_ps(p, x.vec); }
inline void storelh(float *p1, float *p2, const m128f &x) { storel(p1,x); storeh(p2,x); }

inline void store01(float *p0, float *p1, const m128f &x) { store0(p0,x); store1(p1,x); }
inline void store23(float *p0, float *p1, const m128f &x) { store2(p0,x); store3(p1,x); }


inline m128f min  (const m128f &x) { m128f a=m128f(min(unpacklo(x,x),unpackhi(x,x))); return m128f(min(movelh(a,a),movehl(a,a))); }
inline m128f max  (const m128f &x) { m128f a=m128f(max(unpacklo(x,x),unpackhi(x,x))); return m128f(max(movelh(a,a),movehl(a,a))); }

inline m128f hadd (const m128f &x, const m128f &y);
inline m128f hsub (const m128f &x, const m128f &y);
inline m128f h2add(const m128f &x, const m128f &y) { return m128f(movelh(x,y)+movehl(y,x)); }
inline m128f h2sub(const m128f &x, const m128f &y) { return m128f(movelh(x,y)-movehl(y,x)); }
inline m128f h3add(const m128f &x, const m128f &y) { return m128f(unpacklo(x,y)+unpackhi(x,y)); }
inline m128f h3sub(const m128f &x, const m128f &y) { return m128f(unpacklo(x,y)-unpackhi(x,y)); }

inline m128f sum  (const m128f &x);
inline m128f sumlo(const m128f &x);
inline m128f sumhi(const m128f &x);
inline m128f sum  (const m128f &x0, const m128f &x1);
inline m128f sum  (const m128f &x0, const m128f &x1, const m128f &x2, const m128f &x3);

#ifndef SSE3
inline m128f sum  (const m128f &x0) { m128f a=h3add(x0,x0); return h2add(a,a); }
inline m128f sumlo(const m128f &x0) { return x0           +shuffle<_MM_SHUFFLE(3,1,3,1)>(x0,x0); }
inline m128f sumhi(const m128f &x0) { return movehl(x0,x0)+shuffle<_MM_SHUFFLE(2,3,2,3)>(x0,x0); }
inline m128f sum  (const m128f &x0, const m128f &x1) { m128f a=h3add(x0,x1); return h2add(a,a); }
inline m128f sum  (const m128f &x0, const m128f &x1, const m128f &x2, const m128f &x3) { return h2add(h3add(x0,x1),h3add(x2,x3)); }
inline m128f hadd (const m128f &x, const m128f &y) { return shuffle<_MM_SHUFFLE(2,0,2,0)>(x,y)+shuffle<_MM_SHUFFLE(3,1,3,1)>(x,y); }
inline m128f hsub (const m128f &x, const m128f &y) { return shuffle<_MM_SHUFFLE(2,0,2,0)>(x,y)-shuffle<_MM_SHUFFLE(3,1,3,1)>(x,y); }
#endif



//{unsecret}
//{group:SSE}
//Summary: SSE register repesenting 2 single precision complexes
class m128cf
{
  public:
    typedef m128cf self;
    typedef complex<float> value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;
    
    enum { N=2 };    

  public:
    __aligned __m128 vec;

  public:
    inline m128cf() {};
	  explicit inline m128cf(__m128 m) : vec(m)  {};
    inline explicit m128cf(const value_type &x) : vec(_mm_set_ps(x.imag(),x.real(),x.imag(),x.real())) {}
    inline m128cf(float x0, float x1, float x2, float x3) : vec(_mm_set_ps(x3,x2,x1,x0)) {}
    inline m128cf(const value_type &x0, const value_type &x1) : vec(_mm_set_ps(x1.imag(),x1.real(),x0.imag(),x0.real())) {}
    inline m128cf(const self &x) : vec(x) {}

    inline self& operator=(const value_type &x) { vec=_mm_set_ps(x.imag(),x.real(),x.imag(),x.real()); return *this; }
    inline self& operator=(const self       &x) { vec=x.vec;                                           return *this; }

    inline operator __m128() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m128cf &x)	{ complex<float> *p = (complex<float>*) &x; return os << '(' << *p << ',' << *(p+1) << ')'; }

inline m128cf imul   (const m128cf &x);

inline m128cf operator&(const m128cf &x, const m128cf &y) { return m128cf(_mm_and_ps(x.vec,y.vec)); }
inline m128cf operator|(const m128cf &x, const m128cf &y) { return m128cf(_mm_or_ps (x.vec,y.vec)); }
inline m128cf operator^(const m128cf &x, const m128cf &y) { return m128cf(_mm_xor_ps(x.vec,y.vec)); }

inline m128cf operator-(const m128cf &y) { return m128cf(_mm_sub_ps(_mm_setzero_ps(),y.vec)); }
inline m128cf operator+(const m128cf &x, const m128cf &y) { return m128cf(_mm_add_ps(x.vec,y.vec)); }
inline m128cf operator-(const m128cf &x, const m128cf &y) { return m128cf(_mm_sub_ps(x.vec,y.vec)); }
inline m128cf operator*(const m128cf &x, const m128cf &y);
inline m128cf operator*(const m128cf &x, const complex<float> &y) { __m128 a=_mm_mul_ps(x.vec,_mm_set_ps1(y.real())); __m128 b=_mm_mul_ps((imul(x)).vec,_mm_set_ps1(y.imag())); return m128cf(_mm_add_ps(a,b)); }
inline m128cf operator*(const m128cf &x,       float           y) { return m128cf(_mm_mul_ps(x.vec,_mm_set_ps1(y))); }
inline m128cf operator/(const m128cf &x,       float           y) { return m128cf(_mm_div_ps(x.vec,_mm_set_ps1(y))); }
inline m128cf operator*(const complex<float> &x, const m128cf &y) { return y*x; }
inline m128cf operator*(      float           x, const m128cf &y) { return y*x; }

inline m128cf &operator+=(m128cf &x, const m128cf &y) { return x=x+y; }
inline m128cf &operator-=(m128cf &x, const m128cf &y) { return x=x-y; }
inline m128cf &operator*=(m128cf &x, const m128cf &y) { return x=x*y; }
inline m128cf &operator*=(m128cf &x, const complex<float> &y) { return x=x*y; }
inline m128cf &operator*=(m128cf &x,       float           y) { return x=x*y; }
inline m128cf &operator/=(m128cf &x,       float           y) { return x=x/y; }

inline m128cf cmul (const m128cf &x, const m128cf &y);
inline m128cf icmul(const m128cf &x, const m128cf &y);
inline m128cf mimul(const m128cf &x, const m128cf &y);
inline m128cf imul (const m128cf &x, const m128cf &y);

inline m128cf imul     (const m128cf &x);
inline m128cf mimul    (const m128cf &x);
inline m128cf add_imul (const m128cf &a, const m128cf &x);
inline m128cf sub_imul (const m128cf &a, const m128cf &x);
inline m128cf addsub   (const m128cf &x, const m128cf &y);
inline m128cf subadd   (const m128cf &x, const m128cf &y);

inline m128cf flip   (const m128cf &x) { return m128cf(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(1,0,3,2))); }
inline m128cf flip_ri(const m128cf &x) { return m128cf(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,3,0,1))); }
inline m128cf conj   (const m128cf &x) { return m128cf(_mm_mul_ps(_mm_set_ps(-1,1,-1,1),x.vec)); }

inline m128cf complexlo(const m128f &x, const m128f &y) { return m128cf(_mm_unpacklo_ps(x.vec,y.vec)); }
inline m128cf complexhi(const m128f &x, const m128f &y) { return m128cf(_mm_unpackhi_ps(x.vec,y.vec)); }

inline complex<float> get0(const m128cf &x) { return *(complex<float>*)&(x.vec); }
inline complex<float> get1(const m128cf &x) { m128cf a=flip(x); return get0(a); }

#ifdef SSE2
inline void loadl  (m128cf &x, const complex<float> *p) { __m128d a=_mm_loadl_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128 &)a; }
inline void loadh  (m128cf &x, const complex<float> *p) { __m128d a=_mm_loadh_pd ((__m128d &)x.vec,(double *)p); x.vec=(__m128 &)a; }
#else
inline void loadl  (m128cf &x, const complex<float> *p) { x.vec=_mm_loadl_pi (x.vec,(__m64 *)p); }
inline void loadh  (m128cf &x, const complex<float> *p) { x.vec=_mm_loadh_pi (x.vec,(__m64 *)p); }
#endif
inline void loadlh (m128cf &x, const complex<float> *p1, const complex<float> *p2) { loadl(x,p1); loadh(x,p2); }

inline void storel (complex<float> *p, const m128cf &x) { _mm_storel_pi((__m64 *)p,x.vec); }
inline void storeh (complex<float> *p, const m128cf &x) { _mm_storeh_pi((__m64 *)p,x.vec); }
inline void storelh(complex<float> *p1, complex<float> *p2, const m128cf &x) { storel(p1,x); storeh(p2,x); }

inline void load   (m128cf &x, const complex<float> *p) { x.vec = _mm_load_ps ((float *)p); }
inline void loadu  (m128cf &x, const complex<float> *p) { x.vec = _mm_loadu_ps((float *)p); }
inline void loadr  (m128cf &x, const complex<float> *p) { x.vec = _mm_loadr_ps((float *)p); }
inline void store  (complex<float> *p, const m128cf &x) { _mm_store_ps ((float *)p, x.vec); }
inline void stream (complex<float> *p, const m128cf &x) { _mm_stream_ps((float *)p, x.vec); }
inline void storeu (complex<float> *p, const m128cf &x) { _mm_storeu_ps((float *)p, x.vec); }
inline void storer (complex<float> *p, const m128cf &x) { _mm_storer_ps((float *)p, x.vec); }
inline void storeur(complex<float> *p, const m128cf &x) { storeu(p, flip(x)); }
inline void store0 (complex<float> *p, const m128cf &x) { storel(p,x); }
inline void store1 (complex<float> *p, const m128cf &x) { storeh(p,x); }

inline m128cf movehl(const m128cf &x, const m128cf &y) { return m128cf(_mm_movehl_ps(x.vec,y.vec)); }
inline m128cf movelh(const m128cf &x, const m128cf &y) { return m128cf(_mm_movelh_ps(x.vec,y.vec)); }

inline m128cf i1mul  (const m128cf &x) { return m128cf(_mm_mul_ps( _mm_set_ps(1,-1,1,1),_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,3,1,0)) )); }
inline m128cf mi1mul (const m128cf &x) { return m128cf(_mm_mul_ps( _mm_set_ps(-1,1,1,1),_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,3,1,0)) )); }

inline m128cf conj   (float a, const m128cf &x) { return m128cf(_mm_mul_ps(_mm_set_ps(-a,a,-a,a),x.vec)); }
inline m128cf imul   (float a, const m128cf &x) { return flip_ri(conj(a,x)); }

inline m64f  norm(const m128cf &x);
inline m64f  abs (const m128cf &x)	{ return sqrt(norm(x));}
inline m128f norm(const m128cf &x,const m128cf &y);
inline m128f abs (const m128cf &x,const m128cf &y)	{ return sqrt(norm(x,y));}

inline m128f real(const m128cf &x) { return m128f(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,0,2,0))); }
inline m128f imag(const m128cf &x) { return m128f(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(3,1,3,1))); }
inline m128f real(const m128cf &x,const m128cf &y) { return m128f(_mm_shuffle_ps(x.vec,y.vec,_MM_SHUFFLE(2,0,2,0))); }
inline m128f imag(const m128cf &x,const m128cf &y) { return m128f(_mm_shuffle_ps(x.vec,y.vec,_MM_SHUFFLE(3,1,3,1))); }

inline void trn    (m128cf &x0, m128cf &x1) { m128cf a=movelh(x0,x1); x1=movehl(x1,x0); x0=a; }
inline void up_trn (m128cf &x0, m128cf &x1) { x1=movehl(x1,x0); }
inline void low_trn(m128cf &x0, m128cf &x1) { x0=movelh(x0,x1); }

inline m128cf sum(const m128cf &x) { return x+movehl(x,x); }

#ifndef SSE3
inline m64f  norm(const m128cf &x) { __m128 a=_mm_mul_ps(x,x); return m64f(_mm_add_ps(_mm_shuffle_ps(a,a,_MM_SHUFFLE(2,0,2,0)),_mm_shuffle_ps(a,a,_MM_SHUFFLE(3,1,3,1)))); }
inline m128f norm(const m128cf &x,const m128cf &y)
{
  __m128 a0=_mm_mul_ps(x.vec,x.vec);
  __m128 a1=_mm_mul_ps(y.vec,y.vec);
  return m128f(_mm_add_ps(_mm_shuffle_ps(a0,a1,_MM_SHUFFLE(2,0,2,0)),_mm_shuffle_ps(a0,a1,_MM_SHUFFLE(3,1,3,1))));
}

inline m128cf operator*(const m128cf &x, const m128cf &y)
{
  __m128 a=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,2,0,0)),y.vec);
  __m128 b=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(3,3,1,1)),(imul(y)).vec);
  return m128cf(_mm_add_ps(a,b));
}

inline m128cf cmul(const m128cf &x, const m128cf &y)
{
  __m128 a=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,2,0,0)),y.vec);
  __m128 b=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(3,3,1,1)),(imul(y)).vec);
  return m128cf(_mm_sub_ps(a,b));
}
inline m128cf icmul(const m128cf &x, const m128cf &y)
{
  __m128 a=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(3,3,1,1)),y.vec);
  __m128 b=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,2,0,0)),(imul(y)).vec);
  return m128cf(_mm_add_ps(a,b));
}

inline m128cf mimul(const m128cf &x, const m128cf &y)
{
  __m128 a=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(3,3,1,1)),y.vec);
  __m128 b=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,2,0,0)),(imul(y)).vec);
  return m128cf(_mm_sub_ps(a,b));
}

inline m128cf imul (const m128cf &x, const m128cf &y)
{
  __m128 a=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(3,3,1,1)),y.vec);
  __m128 b=_mm_mul_ps(_mm_shuffle_ps(x.vec,x.vec,_MM_SHUFFLE(2,2,0,0)),(imul(y)).vec);
  return m128cf(_mm_sub_ps(b,a));
}

inline m128cf imul     (const m128cf &x) { return flip_ri(conj(x)); }
inline m128cf mimul    (const m128cf &x) { return flip_ri(x); }
inline m128cf add_imul (const m128cf &a, const m128cf &x) { return a+imul (x); }
inline m128cf sub_imul (const m128cf &a, const m128cf &x) { return a-imul (x); }
inline m128cf addsub   (const m128cf &x, const m128cf &y) { return m128cf(x+conj(y)); }
inline m128cf subadd   (const m128cf &x, const m128cf &y) { return m128cf(x-conj(y)); }
#endif


//{unsecret}
//{group:SSE}
//Summary: SSE register repesenting 2 single precision complexes
class m64cf
{
  public:
    typedef m64cf self;
    typedef complex<float> value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;
    
    enum { N=1 };    

  public:
    __aligned __m128 vec;

  public:
    inline m64cf() {};
	  explicit inline m64cf(__m128 m) : vec(m)  {};
    inline explicit m64cf(const value_type &x) : vec(_mm_set_ps(x.imag(),x.real(),x.imag(),x.real())) {}
    inline m64cf(float x0, float x1) : vec(_mm_set_ps(x1,x0,x1,x0)) {}
    inline m64cf(const self &x) : vec(x) {}

    inline self& operator=(const value_type &x) { vec=_mm_set_ps(x.imag(),x.real(),x.imag(),x.real()); return *this; }
    inline self& operator=(const self       &x) { vec=x.vec;                                           return *this; }

    inline operator __m128() const { return vec; }

   	inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m64cf &x)	{ complex<float> *p = (complex<float>*) &x; return os << '(' << *p << ')'; }

inline m64cf operator&(const m64cf &x, const m64cf &y) { return (m64cf)((m128cf&)x&(m128cf&)y); }
inline m64cf operator|(const m64cf &x, const m64cf &y) { return (m64cf)((m128cf&)x|(m128cf&)y); }
inline m64cf operator^(const m64cf &x, const m64cf &y) { return (m64cf)((m128cf&)x^(m128cf&)y); }

inline m64cf imul     (const m64cf &x) { return (m64cf)flip_ri(conj((m128cf&)x)); }

inline m64cf operator-(                const m64cf          &y) { return (m64cf)(          -(m128cf&)y); }
inline m64cf operator+(const m64cf &x, const m64cf          &y) { return (m64cf)((m128cf&)x+(m128cf&)y); }
inline m64cf operator-(const m64cf &x, const m64cf          &y) { return (m64cf)((m128cf&)x-(m128cf&)y); }
inline m64cf operator*(const m64cf &x, const m64cf          &y) { return (m64cf)((m128cf&)x*(m128cf&)y); }
inline m64cf operator*(const m64cf &x, const complex<float> &y) { return (m64cf)((m128cf&)x*         y); }
inline m64cf operator*(const m64cf &x,       float           y) { return (m64cf)((m128cf&)x*         y); }
inline m64cf operator/(const m64cf &x,       float           y) { return (m64cf)((m128cf&)x/         y); }
inline m64cf operator*(const complex<float> &x, const m64cf &y) { return y*x; }
inline m64cf operator*(      float           x, const m64cf &y) { return y*x; }

inline m64cf &operator+=(m64cf &x, const m64cf          &y) { return x=x+y; }
inline m64cf &operator-=(m64cf &x, const m64cf          &y) { return x=x-y; }
inline m64cf &operator*=(m64cf &x, const m64cf          &y) { return x=x*y; }
inline m64cf &operator*=(m64cf &x, const complex<float> &y) { return x=x*y; }
inline m64cf &operator*=(m64cf &x,       float           y) { return x=x*y; }
inline m64cf &operator/=(m64cf &x,       float           y) { return x=x/y; }

inline m64cf flip   (const m64cf &x) { return x; }
inline m64cf flip_ri(const m64cf &x) { return (m64cf)flip_ri((m128cf&)x); }
inline m64cf conj   (const m64cf &x) { return (m64cf)conj   ((m128cf&)x); }

inline void load   (m64cf &x, const complex<float> *p) { loadl((m128cf&)x,p); }
inline void loadu  (m64cf &x, const complex<float> *p) { load(x,p); }
inline void loadr  (m64cf &x, const complex<float> *p) { load(x,p); }

inline void store  (complex<float> *p, const m64cf &x) { storel(p,(m128cf&)x); }
inline void storeu (complex<float> *p, const m64cf &x) { store(p,x); }
inline void storer (complex<float> *p, const m64cf &x) { store(p,x); }
inline void storeur(complex<float> *p, const m64cf &x) { store(p,x); }
inline void store0 (complex<float> *p, const m64cf &x) { store(p,x); }



//{secret}
//{group:SSE}
//Summary: Two SSE registers repesenting 4 single precision complexes
class m256cf
{
  public:
    typedef m256cf self;
    typedef complex<float> value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;

    enum { N=4 };

  public:
    m128cf vec0;
    m128cf vec1;

  public:
    inline m256cf() {};
	  inline m256cf(const __m128 &m0,const __m128 &m1) : vec0(m0),vec1(m1)  {};
	  inline m256cf(const m128cf &m0,const m128cf &m1) : vec0(m0),vec1(m1)  {};
    inline explicit m256cf(const value_type &x) { vec0=m128cf(x); vec1=vec0; }
    inline m256cf(float x0, float x1, float x2, float x3,float x4, float x5, float x6, float x7) : vec0(x0,x1,x2,x3), vec1(x4,x5,x6,x7) {}
    inline m256cf(const value_type &x0, const value_type &x1, const value_type &x2, const value_type &x3) : vec0(x0,x1), vec1(x2,x3) {}
    inline m256cf(const self &x) : vec0(x.vec0), vec1(x.vec1) {}

   	inline reference       operator[](int i)       { return *(((value_type *)&vec0)+i); };
   	inline const_reference operator[](int i) const { return *(((value_type *)&vec0)+i); };

    inline self& operator=(const value_type &x) { vec0=m128cf(x); vec1=vec0;   return *this; }
    inline self& operator=(const self       &x) { vec0=x.vec0;    vec1=x.vec1; return *this; }
};


template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m256cf &x)	{ complex<float> *p = (complex<float>*) &x; return os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) << ')'; }

inline m256cf operator&(const m256cf &x, const m256cf &y) { return m256cf(x.vec0&y.vec0,x.vec1&y.vec1); }
inline m256cf operator|(const m256cf &x, const m256cf &y) { return m256cf(x.vec0|y.vec0,x.vec1|y.vec1); }
inline m256cf operator^(const m256cf &x, const m256cf &y) { return m256cf(x.vec0^y.vec0,x.vec1^y.vec1); }

inline m256cf operator-(                 const m256cf &y) { return m256cf(      -y.vec0,      -y.vec1); }
inline m256cf operator+(const m256cf &x, const m256cf &y) { return m256cf(x.vec0+y.vec0,x.vec1+y.vec1); }
inline m256cf operator-(const m256cf &x, const m256cf &y) { return m256cf(x.vec0-y.vec0,x.vec1-y.vec1); }
inline m256cf operator*(const m256cf &x, const m256cf &y) { return m256cf(x.vec0*y.vec0,x.vec1*y.vec1); }
inline m256cf operator*(const m256cf &x, const complex<float> &y) { return m256cf(x.vec0*y,x.vec1*y); }
inline m256cf operator*(const m256cf &x,       float           y) { return m256cf(x.vec0*y,x.vec1*y); }
inline m256cf operator/(const m256cf &x,       float           y) { return m256cf(x.vec0/y,x.vec1/y); }
inline m256cf operator*(const complex<float> &x, const m256cf &y) { return y*x; }
inline m256cf operator*(      float           x, const m256cf &y) { return y*x; }

inline m256cf &operator+=(m256cf &x, const m256cf &y) { x.vec0+=y.vec0; x.vec1+=y.vec1; return x; }
inline m256cf &operator-=(m256cf &x, const m256cf &y) { x.vec0-=y.vec0; x.vec1-=y.vec1; return x; }
inline m256cf &operator*=(m256cf &x, const m256cf &y) { x.vec0*=y.vec0; x.vec1*=y.vec1; return x; }
inline m256cf &operator*=(m256cf &x, const complex<float> &y) { x.vec0*=y; x.vec1*=y; return x; }
inline m256cf &operator*=(m256cf &x,       float           y) { x.vec0*=y; x.vec1*=y; return x; }
inline m256cf &operator/=(m256cf &x,       float           y) { x.vec0/=y; x.vec1/=y; return x; }

inline m256cf cmul (const m256cf &x, const m256cf &y) { return m256cf(cmul (x.vec0,y.vec0),cmul (x.vec1,y.vec1)); }
inline m256cf icmul(const m256cf &x, const m256cf &y) { return m256cf(icmul(x.vec0,y.vec0),icmul(x.vec1,y.vec1)); }
inline m256cf mimul(const m256cf &x, const m256cf &y) { return m256cf(mimul(x.vec0,y.vec0),mimul(x.vec1,y.vec1)); }
inline m256cf imul (const m256cf &x, const m256cf &y) { return m256cf(imul (x.vec0,y.vec0),imul (x.vec1,y.vec1)); }

inline m256cf imul     (const m256cf &x) { return m256cf(imul (x.vec0),imul (x.vec1)); };
inline m256cf mimul    (const m256cf &x) { return m256cf(mimul(x.vec0),mimul(x.vec1)); }
inline m256cf add_imul (const m256cf &x, const m256cf &y) { return m256cf(add_imul(x.vec0,y.vec0),add_imul(x.vec1,y.vec1)); }
inline m256cf sub_imul (const m256cf &x, const m256cf &y) { return m256cf(sub_imul(x.vec0,y.vec0),sub_imul(x.vec1,y.vec1)); }
inline m256cf addsub   (const m256cf &x, const m256cf &y) { return m256cf(addsub  (x.vec0,y.vec0),addsub  (x.vec1,y.vec1)); }
inline m256cf subadd   (const m256cf &x, const m256cf &y) { return m256cf(subadd  (x.vec0,y.vec0),subadd  (x.vec1,y.vec1)); }

inline m256cf flip   (const m256cf &x) { return m256cf(flip   (x.vec0),flip   (x.vec1)); }
inline m256cf flip_ri(const m256cf &x) { return m256cf(flip_ri(x.vec0),flip_ri(x.vec1)); }
inline m256cf conj   (const m256cf &x) { return m256cf(conj   (x.vec0),conj   (x.vec1)); }

inline void load  (m256cf &x, const complex<float> *p) { load (x.vec0,p); load (x.vec1,p+2); }
inline void loadu (m256cf &x, const complex<float> *p) { loadu(x.vec0,p); loadu(x.vec1,p+2); }
inline void loadr (m256cf &x, const complex<float> *p) { loadr(x.vec1,p); loadr(x.vec0,p+2); }
inline void loadl (m256cf &x, const complex<float> *p) { load(x.vec0,p); }
inline void loadh (m256cf &x, const complex<float> *p) { load(x.vec1,p); }
inline void loadlh(m256cf &x, const complex<float> *p1, const complex<float> *p2) { loadl(x,p1); loadh(x,p2); }

inline void store0 (complex<float> *p, const m256cf &x) { storel(p,x.vec0); }
inline void store1 (complex<float> *p, const m256cf &x) { storeh(p,x.vec0); }
inline void store  (complex<float> *p, const m256cf &x) { store (p,x.vec0); store (p+2,x.vec1); }
inline void storeu (complex<float> *p, const m256cf &x) { storeu(p,x.vec0); storeu(p+2,x.vec1); }
inline void storer (complex<float> *p, const m256cf &x) { storer(p,x.vec1); storer(p+2,x.vec0); }
inline void storeur(complex<float> *p, const m256cf &x) { storeu(p, flip(x)); }
inline void storel (complex<float> *p, const m256cf &x) { store (p,x.vec0); }
inline void storeh (complex<float> *p, const m256cf &x) { store (p,x.vec1); }
inline void storelh(complex<float> *p1, complex<float> *p2, const m256cf &x) { storel(p1,x); storeh(p2,x); }
inline void stream (complex<float> *p, const m256cf &x) { stream(p,x.vec0); stream(p+2,x.vec1); }

inline m256cf conj   (float a, const m256cf &x) { return m256cf(conj(a,x.vec0),conj(a,x.vec1)); }
inline m256cf imul   (float a, const m256cf &x) { return m256cf(imul(a,x.vec0),imul(a,x.vec1)); }

inline m128f abs (const m256cf &x) { return abs (x.vec0,x.vec1);}
inline m128f norm(const m256cf &x) { return norm(x.vec0,x.vec1); }
inline m128f real(const m256cf &x) { return real(x.vec0,x.vec1); }
inline m128f imag(const m256cf &x) { return imag(x.vec0,x.vec1); }

inline m128cf sum(const m256cf &x) { return sum(x.vec0+x.vec1); }



#endif // NO_SIMD

#endif // SSE_H
