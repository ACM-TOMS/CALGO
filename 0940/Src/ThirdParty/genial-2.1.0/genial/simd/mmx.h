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

#ifndef MMX_H
#define MMX_H

#ifndef NO_SIMD

#ifndef MMX
#define MMX
#endif

#include <mmintrin.h>
#include <iostream>

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

using namespace std;

class m64c;
class m64b;
class m64s;
class m64w;
class m64i;


//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 8 signed chars
class m64c
{
  public:
	  typedef m64c self;
  	
	  typedef char value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=8 };

  public:
    __m64 vec;

  public:
    inline m64c(       ) {}
    inline m64c(__m64 m) : vec(m) {}
    inline explicit m64c(value_type v) : vec(_mm_set1_pi8(v))	{}		
    inline m64c(value_type v0, value_type v1, value_type v2, value_type v3) : vec(_mm_set_pi8(v3,v2,v1,v0,v3,v2,v1,v0))  {}	
    inline m64c(value_type v0, value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7) : vec(_mm_set_pi8(v7,v6,v5,v4,v3,v2,v1,v0)) {}	

    inline m64c& operator=(value_type v) { vec=_mm_set1_pi8(v); return *this; }

    inline operator __m64() const { return vec; }
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m64c &a)	 { char *p=(char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) <<')'; return os; }

inline m64c operator- (               const m64c &y) { return m64c(_mm_sub_pi8(_mm_setzero_si64(),y.vec)); } 
inline m64c operator+ (const m64c &x, const m64c &y) { return m64c(_mm_add_pi8(x.vec,y.vec)); }
inline m64c operator- (const m64c &x, const m64c &y) { return m64c(_mm_sub_pi8(x.vec,y.vec)); } 

inline void load (m64c &x, const char *p) { x.vec = *(__m64 *)p; }
inline void loadu(m64c &x, const char *p) { x.vec = *(__m64 *)p; }
inline void store(char *p, const m64c &x) { *(__m64 *)p = x.vec; }

inline m64c unpacklo(const m64c &a, const m64c &b) { return m64c(_mm_unpacklo_pi8 (a.vec,b.vec)); }

template<int N> inline m64c shiftl   (const m64c &a) { return m64c(_mm_srli_si64(a.vec,8*N)); } 
template<int N> inline m64c shiftr   (const m64c &a) { return m64c(_mm_slli_si64(a.vec,8*N)); } 


//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 8 signed chars with saturation
class m64cs
{
  public:
	  typedef m64cs self;
  	
	  typedef char value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=8 };

  public:
    __m64 vec;

  public:
    inline m64cs(       ) {}
    inline m64cs(__m64 m) : vec(m) {}
    inline explicit m64cs(value_type v) : vec(_mm_set1_pi8(v))	{}		
    inline m64cs(value_type v0, value_type v1, value_type v2, value_type v3) : vec(_mm_set_pi8(v3,v2,v1,v0,v3,v2,v1,v0))  {}	
    inline m64cs(value_type v0, value_type v1, value_type v2, value_type v3, value_type v4, value_type v5, value_type v6, value_type v7) : vec(_mm_set_pi8(v7,v6,v5,v4,v3,v2,v1,v0)) {}	

    inline m64cs& operator=(value_type v) { vec=_mm_set1_pi8(v); return *this; }

    inline operator __m64() const { return vec; }
};

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os,const m64cs &x)	 { char *p=(char*)&x; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) <<')'; return os; }

inline m64cs operator- (                const m64cs &y) { return m64cs(_mm_subs_pi8(_mm_setzero_si64(),y.vec)); } 
inline m64cs operator+ (const m64cs &x, const m64cs &y) { return m64cs(_mm_adds_pi8(x.vec,y.vec)); }
inline m64cs operator- (const m64cs &x, const m64cs &y) { return m64cs(_mm_subs_pi8(x.vec,y.vec)); } 


//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 8 unsigned chars
class m64b
{
  public:
	  typedef m64b self;

    typedef unsigned char value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=8 };
    
  public:
    __m64 vec;

  public:
    inline m64b() {}
    inline m64b(__m64 m) : vec(m) {}
    inline explicit m64b(value_type v) : vec(_mm_set1_pi8(v))	{}		
    inline m64b(value_type v0,value_type v1,value_type v2,value_type v3) : vec(_mm_set_pi8(v3,v2,v1,v0,v3,v2,v1,v0)) {}	
    inline m64b(value_type v0,value_type v1,value_type v2,value_type v3,value_type v4,value_type v5,value_type v6,value_type v7) : vec(_mm_set_pi8(v7,v6,v5,v4,v3,v2,v1,v0)) {}	

    inline m64b& operator=(value_type v) { vec=_mm_set1_pi8(v); return *this; }

    inline operator __m64() const { return vec; }
};

template<class E,class Tr>  basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64b &a)	{ unsigned char *p=(unsigned char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) <<')'; return os; }

inline  m64b operator+(const m64b &a , const m64b  &b) { return m64b(_mm_add_pi8(a.vec,b.vec)); }
inline  m64b operator+(const m64c &a, const m64b  &b) { return m64b(_mm_add_pi8(a.vec,b.vec)); }
inline  m64b operator+(const m64b &a , const m64c &b) { return b+a; }
inline  m64b operator-(const m64b &a,  const m64b  &b) { return m64b(_mm_sub_pi8(a.vec,b.vec)); } 
inline  m64b operator-(const m64b &a                 ) { return m64b(_mm_sub_pi8(_mm_setzero_si64(),a.vec)); } 

inline  void load (m64b &x, const unsigned char *p) { x.vec = *(__m64 *)p; }
inline  void loadu(m64b &x, const unsigned char *p) { x.vec = *(__m64 *)p; }
inline  void store(unsigned char *p, const m64b &x) { *(__m64 *)p = x.vec; }

inline  m64b unpacklo(const m64b &a, const m64b &b) { return m64b(_mm_unpacklo_pi8 (a.vec,b.vec)); }
inline  m64s unpacklo(const m64b &a);
inline  m64s unpackhi(const m64b &a);
inline  m64b shiftr1 (const m64b &a               ) { return m64b(_mm_and_si64(_mm_srli_pi16(a.vec,1),_mm_setr_pi8(127,-128,127,-128,127,-128,127,-128))); }
inline  m64b abs_diff(const m64b &a,  const m64b  &b) { return m64b(_mm_or_si64(_mm_subs_pu8(a.vec,b.vec),_mm_subs_pu8(b.vec,a.vec))); }

template<int N> inline  m64b shiftl   (const m64b &a) { return m64b(_mm_srli_si64(a.vec,8*N)); } 
template<int N> inline  m64b shiftr   (const m64b &a) { return m64b(_mm_slli_si64(a.vec,8*N)); } 



//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 8 unsigned chars with saturation 
class m64bs
{
  public:
	  typedef m64bs self;
  	
	  typedef unsigned char value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=8 };

  public:
    __m64 vec;

  public:
    inline m64bs(       ) {}
    inline m64bs(__m64 m) : vec(m) {}
    inline explicit m64bs(value_type v) : vec(_mm_set1_pi8(v)) {}		
    inline m64bs(value_type v0,value_type v1,value_type v2,value_type v3)  : vec(_mm_set_pi8(v3,v2,v1,v0,v3,v2,v1,v0)) {}	
    inline m64bs(value_type v0,value_type v1,value_type v2,value_type v3,value_type v4,value_type v5,value_type v6,value_type v7)  : vec(_mm_set_pi8(v7,v6,v5,v4,v3,v2,v1,v0)) {}	

    inline m64bs& operator=(value_type v) { vec=_mm_set1_pi8(v);  return *this; }

    inline operator __m64() const { return vec; }
};
template<class E,class Tr>  basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64bs &a)	{ unsigned char *p=(unsigned char*)&a; os << '(' << (int)*p     << ',' << (int)*(p+1) << ',' << (int)*(p+2) << ',' << (int)*(p+3) << ',' << (int)*(p+4) << ',' << (int)*(p+5) << ',' << (int)*(p+6) << ',' << (int)*(p+7) <<')'; return os; }

inline  m64bs operator+(const m64bs &a, const m64bs &b) { return m64bs(_mm_adds_pu8(a.vec,b.vec)); }
inline  m64bs operator+(const m64c  &a, const m64bs &b) { return m64bs(_mm_adds_pu8(a.vec,b.vec)); }
inline  m64bs operator+(const m64bs &a, const m64c  &b) { return b+a; }
inline  m64bs operator-(const m64bs &a, const m64bs &b) { return m64bs(_mm_subs_pu8(a.vec,b.vec)); } 

inline  void load (m64bs &x, const unsigned char *p) { x.vec = *(__m64 *)p; }
inline  void loadu(m64bs &x, const unsigned char *p) { x.vec = *(__m64 *)p; }
inline  void store(unsigned char *p, const m64bs &x) { *(__m64 *)p = x.vec; }

inline  m64bs unpacklo(const m64bs &a, const m64bs &b) { return m64bs(_mm_unpacklo_pi8 (a.vec,b.vec)); }
inline  m64bs abs_diff(const m64bs &a, const m64bs &b) { return m64bs(_mm_or_si64(_mm_subs_pu8(a.vec,b.vec),_mm_subs_pu8(b.vec,a.vec))); }

template<int N> inline  m64bs shiftl   (const m64bs &a) { return m64bs(_mm_srli_si64(a.vec,8*N)); } 
template<int N> inline  m64bs shiftr   (const m64bs &a) { return m64bs(_mm_slli_si64(a.vec,8*N)); } 


//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 4 signed short integers
class m64s
{
  public:
	  typedef m64s self;
	  
	  typedef short value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;
	  
	  enum { N=4 };

  public:
    __m64 vec;

  public:
    inline m64s(       ) {}
    inline m64s(__m64 m) : vec(m) {}
    inline explicit m64s(value_type v) : vec(_mm_set1_pi16(v)) {}		
    inline m64s(value_type v0,value_type v1,value_type v2,value_type v3) : vec(_mm_set_pi16 (v3,v2,v1,v0)) {}

    inline m64s& operator=(value_type v) { vec = _mm_set1_pi16(v);  return *this; }

    inline operator __m64() const { return vec; }
    
    inline reference       operator[](int i)       { return *(((value_type *)&vec)+i); };    
    inline const_reference operator[](int i) const { return *(((value_type *)&vec)+i); };    
};
template<class E,class Tr>  basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64s &a)	{ short *p = (short *)&a; os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) <<')'; return os; }

inline m64s operator-(               const m64s &b) { return m64s(_mm_sub_pi16     (_mm_setzero_si64(),b.vec)); } 
inline m64s operator+(const m64s &a, const m64s &b) { return m64s(_mm_add_pi16     (a.vec,b.vec)); }
inline m64s operator-(const m64s &a, const m64s &b) { return m64s(_mm_sub_pi16     (a.vec,b.vec)); } 
inline m64s operator*(const m64s &a, const m64s &b) { return m64s(_mm_mullo_pi16   (a.vec,b.vec)); }
inline m64s mullo    (const m64s &a, const m64s &b) { return m64s(_mm_mullo_pi16   (a.vec,b.vec)); }
inline m64s mulhi    (const m64s &a, const m64s &b) { return m64s(_mm_mulhi_pi16   (a.vec,b.vec)); }

inline m64s operator==(const m64s &a, const m64s &b) { return m64s(_mm_cmpeq_pi16   (a.vec,b.vec)); } 
inline m64s operator> (const m64s &a, const m64s &b) { return m64s(_mm_cmpgt_pi16   (a.vec,b.vec)); } 
inline m64s operator< (const m64s &a, const m64s &b) { return m64s(_mm_cmpgt_pi16   (b.vec,a.vec)); } 
inline m64s operator& (const m64s &a, const m64s &b) { return m64s(_mm_and_si64     (a.vec,b.vec)); } 
inline m64s operator| (const m64s &a, const m64s &b) { return m64s(_mm_or_si64      (a.vec,b.vec)); } 
inline m64s and_not   (const m64s &a, const m64s &b) { return m64s(_mm_andnot_si64  (a.vec,b.vec)); } 

inline void load (m64s &x, const short *p) { x.vec = *(__m64 *)p; }
inline void loadu(m64s &x, const short *p) { x.vec = *(__m64 *)p; }
inline void store(short *p, const m64s &x) { *(__m64 *)p = x.vec; }

inline m64s unpacklo (const m64s &a, const m64s &b) { return m64s(_mm_unpacklo_pi16(a.vec,b.vec)); }
inline m64s unpackhi (const m64s &a, const m64s &b) { return m64s(_mm_unpackhi_pi16(a.vec,b.vec)); }
inline m64c packs    (const m64s &a, const m64s &b) { return m64c(_mm_packs_pi16   (a.vec,b.vec)); }
inline m64b packus   (const m64s &a, const m64s &b) { return m64b(_mm_packs_pu16   (a.vec,b.vec)); }


inline m64i madd     (const m64s &a, const m64s &b);

template<int N> inline m64s shiftra(const m64s &a) { return m64s(_mm_srai_pi16(a.vec,N)); }
template<int N> inline m64s shiftla(const m64s &a) { return m64s(_mm_slli_pi16(a.vec,N)); } 
template<int N> inline m64s shiftrl(const m64s &a) { return m64s(_mm_srli_pi16(a.vec,N)); }
template<int N> inline m64s shiftll(const m64s &a) { return m64s(_mm_slli_pi16(a.vec,N)); } 

template<int N> inline m64s shiftl   (const m64s &a) { return m64s(_mm_srli_si64(a.vec,16*N)); } 
template<int N> inline m64s shiftr   (const m64s &a) { return m64s(_mm_slli_si64(a.vec,16*N)); } 

inline m64s sum(const m64s &a) { __m64 c=_mm_adds_pu16(a.vec,_mm_srli_si64(a.vec,32)); __m64 b=_mm_adds_pu16(c,_mm_srli_si64(c,16)); return m64s(b); }

inline short get_sad(const m64s &a) { return *(short *)&a; }


#ifndef SSE
inline m64s min(const m64s &a,const m64s &b) { m64s c=a>b; return and_not(c,a)|(c&b); }
inline m64s max(const m64s &a,const m64s &b) { m64s c=b>a; return and_not(c,a)|(c&b); }
inline m64s sad(const m64b &a,const m64b &b)
{ 
  __m64 c=_mm_or_si64(_mm_subs_pu8(a.vec,b.vec),_mm_subs_pu8(b.vec,a.vec));
  __m64 d=_mm_adds_pu16(_mm_unpacklo_pi8(c,_mm_setzero_si64()),_mm_unpackhi_pi8(c,_mm_setzero_si64()));
        c=_mm_adds_pu16(d,_mm_srli_si64(d,32));
        d=_mm_adds_pu16(c,_mm_srli_si64(c,16));

  return m64s(_mm_and_si64(d,_mm_set_pi16(0,0,0,-1)));
}
inline m64s flip(const m64s &a) { short *p = (short *)&a; return m64s(*(p+3),*(p+2),*(p+1),*p); }	      
#else
inline short get0(const m64s &x) { return _mm_extract_pi16(x.vec,0); }
inline m64s min (const m64s &a,const m64s &b) { return m64s(_mm_min_pi16(a.vec,b.vec)); }
inline m64s max (const m64s &a,const m64s &b) { return m64s(_mm_max_pi16(a.vec,b.vec)); }
inline m64s sad (const m64b &a,const m64b &b) { return m64s(_mm_sad_pu8(a.vec,b.vec)); }
inline m64s flip(const m64s &a) { return m64s(_mm_shuffle_pi16(a.vec,_MM_SHUFFLE(0,1,2,3))); }	  
#endif
inline m64s abs (const m64s &a) { return max(a,-a); }



//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 4 signed short integers with saturation
class m64ss
{
  public:
	  typedef m64ss self;
	  
	  typedef short value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=4 };

  public:
    __m64 vec;

  public:
    inline m64ss(       ) {}
    inline m64ss(__m64 m) : vec(m) {}
    inline explicit m64ss(value_type v) : vec(_mm_set1_pi16(v)) {}		
    inline m64ss(value_type v0,value_type v1,value_type v2,value_type v3) : vec(_mm_set_pi16 (v3,v2,v1,v0)) {}

    inline m64ss& operator=(value_type v) { vec = _mm_set1_pi16(v);  return *this; }

    inline operator __m64() const { return vec; }
};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64ss &a)	{ short *p = (short *)&a; os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) <<')'; return os; }

inline m64ss operator+(const m64ss &a, const m64ss &b) { return m64ss(_mm_adds_pi16     (a.vec,b.vec)); }
inline m64ss operator-(const m64ss &a, const m64ss &b) { return m64ss(_mm_subs_pi16     (a.vec,b.vec)); } 
inline m64ss operator+(const m64ss &a, const m64s  &b) { return m64ss(_mm_adds_pi16     (a.vec,b.vec)); }
inline m64ss operator-(const m64ss &a, const m64s  &b) { return m64ss(_mm_subs_pi16     (a.vec,b.vec)); } 
inline m64ss operator+(const m64s  &a, const m64ss &b) { return m64ss(_mm_adds_pi16     (a.vec,b.vec)); }
inline m64ss operator-(const m64s  &a, const m64ss &b) { return m64ss(_mm_subs_pi16     (a.vec,b.vec)); } 

inline void load (m64ss &x, const short *p) { x.vec = *(__m64 *)p; }
inline void loadu(m64ss &x, const short *p) { x.vec = *(__m64 *)p; }
inline void store(short *p, const m64ss &x) { *(__m64 *)p = x.vec; }

inline m64s unpacklo(const m64b &a) { return m64s(_mm_unpacklo_pi8 (a.vec,_mm_setzero_si64())); }
inline m64s unpackhi(const m64b &a) { return m64s(_mm_unpackhi_pi8 (a.vec,_mm_setzero_si64())); }

#ifndef SSE
inline m64ss flip(const m64ss &a) { short *p = (short *)&a; return m64ss(*(p+3),*(p+2),*(p+1),*p); }	      
#else
inline m64ss flip(const m64ss &a) { return m64ss(_mm_shuffle_pi16(a.vec,_MM_SHUFFLE(0,1,2,3))); }	  
#endif



//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 4 unsigned short integers
class m64w
{
  public:
	  typedef m64w self;
  	
	  typedef unsigned short value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=4 };

  public:
    __m64 vec;

  public:
    inline m64w(       ) {}
    inline m64w(__m64 m) : vec(m) {}
    inline explicit m64w(value_type v) : vec(_mm_set1_pi16(v)) {}		
    inline m64w(value_type v0,value_type v1,value_type v2,value_type v3) : vec(_mm_set_pi16 (v3,v2,v1,v0)) {}

    inline self &operator=(value_type v) { vec = _mm_set1_pi16(v);  return *this; }

    inline operator __m64() const { return vec; }

};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64w &a)	{ unsigned short *p=(unsigned short *)&a; os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) <<')'; return os; }

inline m64w operator+(const m64w &a, const m64w &b) { return m64w(_mm_add_pi16    (a.vec,b.vec)); }
inline m64w operator-(const m64w &a, const m64w &b) { return m64w(_mm_sub_pi16    (a.vec,b.vec)); } 

inline void load (m64w &x, const unsigned short *p) { x.vec = *(__m64 *)p; }
inline void loadu(m64w &x, const unsigned short *p) { x.vec = *(__m64 *)p; }
inline void store(unsigned short *p, const m64w &x) { *(__m64 *)p = x.vec; }

inline m64b packs(const m64w &a, const m64w &b) { return m64b(_mm_packs_pu16   (a.vec,b.vec)); }

#ifndef SSE
inline m64w flip(const m64w &a) { unsigned short *p = (unsigned short *)&a; return m64w(*(p+3),*(p+2),*(p+1),*p); }	      
#else
inline m64w flip(const m64w &a) { return m64w(_mm_shuffle_pi16(a.vec,_MM_SHUFFLE(0,1,2,3))); }	  
#endif



//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 4 unsigned short integers with saturation
class m64ws
{
  public:
	  typedef m64ws self;
  	
	  typedef unsigned short value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=4 };

  public:
    __m64 vec;

  public:
    inline m64ws(       ) {}
    inline m64ws(__m64 m) : vec(m) {}
    inline explicit m64ws(value_type v) : vec(_mm_set1_pi16(v)) {}		
    inline m64ws(value_type v0,value_type v1,value_type v2,value_type v3) : vec(_mm_set_pi16 (v3,v2,v1,v0)) {}

    inline self &operator=(value_type v) { vec = _mm_set1_pi16(v);  return *this; }

    inline operator __m64() const { return vec; }
};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64ws &a)	{ unsigned short *p=(unsigned short *)&a; os << '(' << *p << ',' << *(p+1) << ',' << *(p+2) << ',' << *(p+3) <<')'; return os; }

inline m64ws operator+(const m64ws &a, const m64ws &b) { return m64ws(_mm_adds_pu16    (a.vec,b.vec)); }
inline m64ws operator-(const m64ws &a, const m64ws &b) { return m64ws(_mm_subs_pu16    (a.vec,b.vec)); } 

inline void load (m64ws &x, const unsigned short *p) { x.vec = *(__m64 *)p; }
inline void loadu(m64ws &x, const unsigned short *p) { x.vec = *(__m64 *)p; }
inline void store(unsigned short *p, const m64ws &x) { *(__m64 *)p = x.vec; }

#ifndef SSE
inline m64ws flip(const m64ws &a) { unsigned short *p = (unsigned short *)&a; return m64ws(*(p+3),*(p+2),*(p+1),*p); }	      
#else
inline m64ws flip(const m64ws &a) { return m64ws(_mm_shuffle_pi16(a.vec,_MM_SHUFFLE(0,1,2,3))); }	  
#endif



//{unsecret}
//{group:MMX}
//Summary: MMX register repesenting 2 signed integers
class m64i
{
  public:
	  typedef m64i self;
  	
	  typedef int value_type;
	  typedef value_type       &reference;
	  typedef const value_type &const_reference;

    enum { N=2 };

  public:
    __m64 vec;

  public:
    inline m64i(       ) {}
    inline m64i(__m64 m) : vec(m) {}
    inline explicit m64i(value_type v) : vec(_mm_set1_pi32(v)) {}		
    inline m64i(value_type v0,value_type v1) : vec(_mm_set_pi32 (v1,v0)) {}

    inline self &operator=(value_type v) { vec = _mm_set1_pi32(v);  return *this; }

    inline operator __m64() const { return vec; }

};
template<class E,class Tr> basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os,const m64i &a)	{ int *p=(int *)&a; os << '(' << *p << ',' << *(p+1) <<')'; return os; }

inline m64i operator+ (const m64i &a, const m64i &b) { return m64i(_mm_add_pi32     (a.vec,b.vec)); }
inline m64i operator- (const m64i &a, const m64i &b) { return m64i(_mm_sub_pi32     (a.vec,b.vec)); } 

inline void load (m64i &x, const int *p) { x.vec = *(__m64 *)p; }
inline void loadu(m64i &x, const int *p) { x.vec = *(__m64 *)p; }
inline void store(int *p, const m64i &x) { *(__m64 *)p = x.vec; }

inline m64i unpacklo  (const m64i &a, const m64i &b) { return m64i(_mm_unpacklo_pi32(a.vec,b.vec)); }
inline m64i unpackhi  (const m64i &a, const m64i &b) { return m64i(_mm_unpackhi_pi32(a.vec,b.vec)); }
inline m64s packs     (const m64i &a, const m64i &b) { return m64s(_mm_packs_pi32   (a.vec,b.vec)); }

template<int N> inline m64i shiftl (const m64i &a) { return m64i(_mm_srli_si64(a.vec,32*N)); } 
template<int N> inline m64i shiftr (const m64i &a) { return m64i(_mm_slli_si64(a.vec,32*N)); }
template<int N> inline m64i shiftll(const m64i &a) { return m64i(_mm_slli_pi32(a.vec,N)); } 
template<int N> inline m64i shiftrl(const m64i &a) { return m64i(_mm_srli_pi32(a.vec,N)); }
template<int N> inline m64i shiftla(const m64i &a) { return m64i(_mm_slli_pi32(a.vec,N)); } 
template<int N> inline m64i shiftra(const m64i &a) { return m64i(_mm_srai_pi32(a.vec,N)); }

inline m64i madd     (const m64s &a, const m64s &b) { return m64i(_mm_madd_pi16    (a.vec,b.vec)); }




#ifdef SSE
template<int i> inline m64c  shuffle(const m64c  &a) { return m64c (_mm_shuffle_pi16(a.vec,i)); }
template<int i> inline m64b  shuffle(const m64b  &a) { return m64b (_mm_shuffle_pi16(a.vec,i)); }
template<int i> inline m64bs shuffle(const m64bs &a) { return m64bs(_mm_shuffle_pi16(a.vec,i)); }
template<int i> inline m64s  shuffle(const m64s  &a) { return m64s (_mm_shuffle_pi16(a.vec,i)); }
#endif



#endif // NO_SIMD

#endif


