//GENIAL - GENeric Image & Array Library
//Copyright (C) 2007  Patrick LAURENT
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

#ifndef FFT_H
#define FFT_H

#ifdef __cplusplus

#include "array/matrix.h"
#if defined(FFT_THREADING)
#include "threads.h"
using namespace gmt;
#endif

//group=FFT

#ifndef FFT_LEVEL
#define FFT_LEVEL 8
#endif

#ifndef FFT_TWIDDLES
#define FFT_TWIDDLES 1
#endif

#define SQRT2      real_type(1.414213562373095048801688724209698078569671875)
#define SQRT_5_16  real_type(0.559016994374947424102293417182819058860154590)


#define SIN00 real_type(0)
#define SIN01 real_type(0.049067674327418014254954976942682658314745363)
#define SIN02 real_type(0.098017140329560601994195563888641845861136673)
#define SIN03 real_type(0.146730474455361751658850129646717819706215317)
#define SIN04 real_type(0.195090322016128267848284868477022240927691618)
#define SIN05 real_type(0.242980179903263889948274162077471118320990783)
#define SIN06 real_type(0.290284677254462367636192375817395274691476278)
#define SIN07 real_type(0.336889853392220050689253212619147570477766780)
#define SIN08 real_type(0.382683432365089771728459984030398866761344562)
#define SIN09 real_type(0.427555093430282094320966856888798534304578629)
#define SIN10 real_type(0.471396736825997648556387625905254377657460319)
#define SIN11 real_type(0.514102744193221726593693838968815772608049120)
#define SIN12 real_type(0.555570233019602224742830813948532874374937191)
#define SIN13 real_type(0.595699304492433343467036528829969889511926338)
#define SIN15 real_type(0.671558954847018400625376850427421803228750632)
#define SIN14 real_type(0.634393284163645498215171613225493370675687095)
#define SIN16 real_type(0.707106781186547524400844362104849039284835938)
#define SIN17 real_type(0.740951125354959091175616897495162729728955309)
#define SIN18 real_type(0.773010453362736960810906609758469800971041293)
#define SIN19 real_type(0.803207531480644909806676512963141923879569427)
#define SIN20 real_type(0.831469612302545237078788377617905756738560812)
#define SIN21 real_type(0.857728610000272069902269984284770137042490799)
#define SIN22 real_type(0.881921264348355029712756863660388349508442621)
#define SIN23 real_type(0.903989293123443331586200297230537048710132025)
#define SIN24 real_type(0.923879532511286756128183189396788286822416626)
#define SIN25 real_type(0.941544065183020778412509402599502357185589796)
#define SIN26 real_type(0.956940335732208864935797886980269969482849206)
#define SIN27 real_type(0.970031253194543992603984207286100251456865962)
#define SIN28 real_type(0.980785280403230449126182236134239036973933731)
#define SIN29 real_type(0.989176509964780973451673738016243063983689533)
#define SIN30 real_type(0.995184726672196886244836953109479921575474869)
#define SIN31 real_type(0.998795456205172392714771604759100694443203615)
#define SIN32 real_type(1)

#define SIN01_03 real_type(0.866025403784438646763723170752936183471402627)
#define SIN01_05 real_type(0.587785252292473129168705954639072768597652438)
#define SIN02_05 real_type(0.951056516295153572116439333379382143405698634)
#define SIN01_10 real_type(0.309016994374947424102293417182819058860154590)
#define SIN03_10 real_type(0.809016994374947424102293417182819058860154590)


inline bool is_pow2(int n)
{
  return (n&(n-1))==0;
}

template<class G>
PROMOTE(Vector<G>)
fft_to_ifft(const Vector<G> X)
{
  PROMOTE(Vector<G>) T;
  if (X.lower_bound())
    T=extract(periodize(mirror(X)),X.lower_bound());
  else
    T=rotate(flip(X),1);
  return T;
}


template<class A>
class fft_base_adaptor : public array_traits<A>
{
  public:
    typedef fft_base_adaptor self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES

  public:
    mutable array_type &X;

  public:
    inline fft_base_adaptor(array_type &x) : X(x) {}

    inline array_type &array()       { return X; }
    inline array_type &array() const { return X; }

    template<class V> static inline void sfirst0   (V *p, const V &v) { *p=v; }
    template<class V> static inline void sfirst    (V *p, const V &v) { *p=v; }
    template<class V> static inline void sfirst_ru (V *p, const V &v) { sfirst (p,v); }
    template<class V> static inline void ssecond   (V *p, const V &v) { sfirst (p,v); }
    template<class V> static inline void ssecond_ru(V *p, const V &v) { ssecond(p,v); }
    template<class V> static inline void sfirst_ru (Vector<simd_vector_generator<1,V> > *p, const Vector<simd_vector_generator<1,V> > &a) { store((V *)p,a); }
    template<class V> static inline void ssecond_ru(Vector<simd_vector_generator<1,V> > *p, const Vector<simd_vector_generator<1,V> > &a) { store((V *)p,a); }
    template<int M,class V> static inline void sfirst_ru (Vector<simd_vector_generator<M,V> > *p, const Vector<simd_vector_generator<M,V> > &a) { storeur((V *)p-(M-1),a); }
    template<int M,class V> static inline void ssecond_ru(Vector<simd_vector_generator<M,V> > *p, const Vector<simd_vector_generator<M,V> > &a) { sfirst_ru(p,a); }

    inline void first0   (              const value_type &v) const { X[0]=v; }
    inline void first    (index_type i, const value_type &v) const { X[i]=v; }
    inline void second   (index_type i, const value_type &v) const { first(i,v); }
    inline void first_u  (index_type i, const value_type &v) const { first(i,v); }
    inline void second_ru(index_type i, const value_type &v) const { first(i,v); }

#if !defined(_MSC_VER)
    template<      class V> inline void first0   (              const Vector<simd_vector_generator<1,V> > &a) const { store(&X[0],a); }
    template<int M,class V> inline void first0   (              const Vector<simd_vector_generator<M,V> > &a) const { sub<M>(X,0)=a; }
    template<      class V> inline void first    (index_type i, const Vector<simd_vector_generator<1,V> > &a) const { store(&X[i],a); }
    template<int M,class V> inline void first    (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { sub<M>(X,i)=a; }
    template<      class V> inline void second   (index_type i, const Vector<simd_vector_generator<1,V> > &a) const { first(i,a); }
    template<int M,class V> inline void second   (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { first(i,a); }

    template<class V> inline void first_u  (index_type i, const Vector<simd_vector_generator<1,V> > &a) const { store(&X[i],a); }
    template<class V> inline void second_ru(index_type i, const Vector<simd_vector_generator<1,V> > &a) const { store(&X[i],a); }
    template<class V> inline void first_u  (index_type i, const Vector<simd_vector_generator<2,V> > &a) const { storelh(&X[i],&X[i+1],a); }
    template<class V> inline void second_ru(index_type i, const Vector<simd_vector_generator<2,V> > &a) const { storelh(&X[i],&X[i-1],a); }
 
    template<class V> inline void second_ru(index_type i, index_type j,const Vector<simd_vector_generator<2,V> > &a) const { storelh(&X[i],&X[j-1],a); }
#else
    template<int M,class V> inline void first0   (              const Vector<simd_vector_generator<M,V> > &a) const { if (M==1) store(&X[0],a); else sub<M>(X,0)=a; }
    template<int M,class V> inline void first    (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { if (M==1) store(&X[i],a); else sub<M>(X,i)=a; }
    template<int M,class V> inline void second   (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { first(i,a); }

    template<int M,class V> inline void first_u  (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { if (M==1) store(&X[i],a); else storelh(&X[i],&X[i+1],a); }
    template<int M,class V> inline void second_ru(index_type i, const Vector<simd_vector_generator<M,V> > &a) const { if (M==1) store(&X[i],a); else storelh(&X[i],&X[i-1],a); }
 
    template<int M, class V> inline void second_ru(index_type i, index_type j,const Vector<simd_vector_generator<M,V> > &a) const { assert(M==2); storelh(&X[i],&X[j-1],a); }
#endif
};


template<class A>
class fft_adaptor : public fft_base_adaptor<A>
{
  public:
    typedef fft_adaptor self;
    typedef fft_base_adaptor<A> base;
    ARRAY_BASE_TYPES

  public:
    inline fft_adaptor(array_type &x) : base(x) {}
};


template<class A>
class simd_fft_adaptor : public array_traits<A>
{
  public:
    typedef simd_fft_adaptor self;
    typedef array_traits<A> base;
    ARRAY_BASE_TYPES

  public:
    mutable array_type &X;

  public:
    inline simd_fft_adaptor(array_type &x) : X(x) {}

    inline array_type &array()       { return X; }
    inline array_type &array() const { return X; }

    inline void first0(              const value_type &x) const { X[0]=x; }
    inline void first (index_type i, const value_type &x) const { X[i]=x; }
    inline void second(index_type i, const value_type &x) const { first(i,x); }

    template<int M, class V> inline void first0   (              const Vector<simd_vector_generator<M,V> > &a) const { store(&X[0],a); }
    template<int M, class V> inline void first    (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { store(&X[i],a); }
    template<int M, class V> inline void second   (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { first(i,a); }
    template<int M, class V> inline void first_u  (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { storeu (&X[i],a); }
    template<int M, class V> inline void second_r (index_type i, const Vector<simd_vector_generator<M,V> > &a) const { storer (&X[i]-(M-1),a); }
    template<int M, class V> inline void second_ru(index_type i, const Vector<simd_vector_generator<M,V> > &a) const { storeur(&X[i]-(M-1),a); }

    inline void second_ru(index_type i, const value_type &x) const { X[i]=x; }
    template<class V> inline void second_ru(index_type i, index_type j, const Vector<simd_vector_generator<2,V> > &x) const { storelh(&X[i],&X[j]-1,x); }
};

template<int N,class S>
class pair_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef pair_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adaptx;
    adaptor_type adapty;

  public:
    inline pair_fft_adaptor(array_type &x, array_type &y) : adaptx(x), adapty(y) {}

    template<class T> inline void first0 (              const T &v) const { adaptx.first0   (             v ); }
    template<class T> inline void first  (index_type i, const T &v) const { adaptx.first    (i      ,     v ); }
    template<class T> inline void second (index_type i, const T &v) const { adapty.second_ru((N-1)-i,conj(v)); }
};

template<int N,class S>
class full_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef full_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adaptx;
    adaptor_type adapty;
    const void *p0;

  public:
    inline full_fft_adaptor(array_type &x, array_type &y, const void *p) : adaptx(x), adapty(y), p0(p) {}

    inline void first0(              const value_type &x) const
    {
      adaptx.first0(  x);
      value_type *p=&adaptx.array()[0];
      if (p!=p0) adapty.second_ru( N-1   ,conj(x));
    }
    
#if !defined(_MSC_VER)
    template<class V> inline void first0(const Vector<simd_vector_generator<1,V> > &x) const
    {
      adaptx.first0(  x);
      value_type *p=&adaptx.array()[0];
      if (p!=p0) adapty.second_ru( N-1   ,conj(x));
    }
    template<class V> inline void first0(const Vector<simd_vector_generator<2,V> > &x) const
    {
      adaptx.first0(  x);
      value_type *p=&adaptx.array()[0];
      if (p!=p0) adapty.second_ru( N-1   ,conj(x));
      else adapty.second_ru(-1,N-1,conj(x));
    }
#else   
    template<int M,class V> inline void first0(const Vector<simd_vector_generator<M,V> > &x) const
    {
      if (M==1)
      {
        adaptx.first0(  x);
        value_type *p=&adaptx.array()[0];
        if (p!=p0) adapty.second_ru( N-1   ,conj(x));
      }
      else
      {
        adaptx.first0(  x);
        value_type *p=&adaptx.array()[0];
        if (p!=p0) adapty.second_ru( N-1   ,conj(x));
        else adapty.second_ru(-1,N-1,conj(x));
      }
    }    
#endif
    
    template<class T> inline void first (index_type i, const T &v) const { adaptx.first (i,v); adapty.second_ru((N-1)-i,conj(v)); }
    template<class T> inline void second(index_type i, const T &v) const { first(i,v); }
};


template<int N,class Adapt>
class real_fft_adaptor : public array_traits<typename Adapt::array_type>
{
  public:
    typedef Adapt adaptor_type;
    typedef real_fft_adaptor self;
    typedef array_traits<typename Adapt::array_type> base;
    ARRAY_BASE_TYPES
    
  public:
    adaptor_type adapt;

  public:
    inline real_fft_adaptor(array_type &x) : adapt(adaptor_type(x)) {}
    inline real_fft_adaptor(adaptor_type &s) : adapt(s) {}

    inline void first0(              const value_type &x) const { adapt.first0(    x); }
    inline void firstN(              const value_type &x) const { adapt.first (N/2,x); }
    inline void first (index_type i, const value_type &x) const { adapt.first (i  ,x); }
    inline void second(index_type i, const value_type &x) const { adapt.second(i  ,x); }

#if !defined(_MSC_VER)
    template<class V> inline void first0 (              const Vector<simd_vector_generator<1,V> > &x) const { adapt.first(0    ,     x ); }
    template<class V> inline void firstN (              const Vector<simd_vector_generator<1,V> > &x) const { adapt.first(N/2  ,     x ); }
    template<class V> inline void secondN(              const Vector<simd_vector_generator<1,V> > &x) const { adapt.first_u(N/2,conj(x)); }
    template<class V> inline void first  (index_type i, const Vector<simd_vector_generator<1,V> > &x) const { adapt.first(i    ,     x ); }
    template<class V> inline void second (index_type i, const Vector<simd_vector_generator<1,V> > &x) const { adapt.first(N-i  ,conj(x)); }

    template<int M,class V> inline void first0 (              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (0          ,     x ); }
    template<int M,class V> inline void firstN (              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (N/2        ,   x[0]); }
    template<int M,class V> inline void secondN(              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first_u((N+1-M)-N/2,conj(x)); }
    template<int M,class V> inline void first  (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (i          ,     x ); }
    template<int M,class V> inline void second (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { adapt.first_u((N+1-M)-i  ,conj(x)); }
#else
    template<int M,class V> inline void first0 (              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (0          ,     x ); }
    template<int M,class V> inline void firstN (              const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) adapt.first(N/2  ,     x ); else adapt.first  (N/2        ,   x[0]); }
    template<int M,class V> inline void secondN(              const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) adapt.first_u(N/2,conj(x)); else adapt.first_u((N+1-M)-N/2,conj(x)); }
    template<int M,class V> inline void first  (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (i          ,     x ); }
    template<int M,class V> inline void second (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { adapt.first_u((N+1-M)-i  ,conj(x)); }
#endif

    template<class V> inline void first0(              const V &x, const V &y) const { adapt.first0(  value_type(x,y)); }
    template<class V> inline void first (index_type i, const V &x, const V &y) const { adapt.first (i,value_type(x,y)); }
    template<class V> inline void second(index_type i, const V &x, const V &y) const { adapt.second(i,value_type(x,y)); }

    template<int M,class V> inline void first0(              const Vector<simd_vector_generator<M,V> > &x, const Vector<simd_vector_generator<M,V> > &y) const { adapt.first0(  complexlo(x,y)); adapt.first (  M/2*(N/M+1),complexhi(x,y)); }
    template<int M,class V> inline void first (index_type i, const Vector<simd_vector_generator<M,V> > &x, const Vector<simd_vector_generator<M,V> > &y) const { adapt.first (i,complexlo(x,y)); adapt.first (i+M/2*(N/M+1),complexhi(x,y)); }
    template<int M,class V> inline void second(index_type i, const Vector<simd_vector_generator<M,V> > &x, const Vector<simd_vector_generator<M,V> > &y) const { adapt.second(i,complexlo(x,y)); adapt.second(i+M/2*(N/M+1),complexhi(x,y)); }
};


template<int N,class S>
class full_real_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef full_real_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;

  public:
    inline full_real_fft_adaptor(array_type &x) : adapt(adaptor_type(x)) {}
    inline full_real_fft_adaptor(adaptor_type &s) : adapt(s) {}

    inline void firstN(              const value_type &z) const { adapt.second(N/2,conj(z)); adapt.first (N/2,z); }
    inline void first0(              const value_type &z) const {                            adapt.first0(  z); }
    inline void first (index_type i, const value_type &z) const { adapt.second(N-i,conj(z)); adapt.first (i,z); }
    inline void second(index_type i, const value_type &z) const { adapt.first (N-i,conj(z)); adapt.second(i,z); }

#if !defined(_MSC_VER)
    template<class V> inline void first0 (              const Vector<simd_vector_generator<1,V> > &x) const { adapt.first0 (         x ); }
    template<class V> inline void firstN (              const Vector<simd_vector_generator<1,V> > &x) const { adapt.first  (N/2,     x ); }
    template<class V> inline void secondN(              const Vector<simd_vector_generator<1,V> > &x) const { adapt.first_u(N/2,conj(x)); adapt.second(N/2,     x ); }
    template<class V> inline void first  (index_type i, const Vector<simd_vector_generator<1,V> > &x) const { adapt.first  (i  ,     x ); adapt.second(N-i,conj(x)); }
    template<class V> inline void second (index_type i, const Vector<simd_vector_generator<1,V> > &x) const { adapt.first  (N-i,conj(x)); adapt.second(i  ,     x ); }

    template<int M,class V> inline void first0 (              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first0 (               x    ); adapt.second_ru(N-1,conj(x[1])); }
    template<int M,class V> inline void firstN (              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (N/2      ,     x [0]); }
    template<int M,class V> inline void secondN(              const Vector<simd_vector_generator<M,V> > &x) const { adapt.first_u((N/2+1-M),conj(x)[0]); adapt.second   (N/2,flip(x   )); }
    template<int M,class V> inline void first  (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { adapt.first  (i,             x    ); adapt.second_ru(N-i,conj(x   )); }
    template<int M,class V> inline void second (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { adapt.first_u((N+1-M)-i,conj(x   )); adapt.second   (i  ,flip(x   )); }
#else
    template<int M,class V> inline void first0 (              const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) { adapt.first0 (         x );                            } else { adapt.first0 (               x    ); adapt.second_ru(N-1,conj(x[1])); } }
    template<int M,class V> inline void firstN (              const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) { adapt.first  (N/2,     x );                            } else { adapt.first  (N/2      ,     x [0]); } }
    template<int M,class V> inline void secondN(              const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) { adapt.first_u(N/2,conj(x)); adapt.second(N/2,     x ); } else { adapt.first_u((N/2+1-M),conj(x)[0]); adapt.second   (N/2,flip(x   )); } }
    template<int M,class V> inline void first  (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) { adapt.first  (i  ,     x ); adapt.second(N-i,conj(x)); } else { adapt.first  (i,             x    ); adapt.second_ru(N-i,conj(x   )); } }
    template<int M,class V> inline void second (index_type i, const Vector<simd_vector_generator<M,V> > &x) const { if (M==1) { adapt.first  (N-i,conj(x)); adapt.second(i  ,     x ); } else { adapt.first_u((N+1-M)-i,conj(x   )); adapt.second   (i  ,flip(x   )); } }
#endif

    template<class V> inline void first0(              const V &x, const V &y) const { value_type z(x,y); adapt.first0(  z); }
    template<class V> inline void first (index_type i, const V &x, const V &y) const { value_type z(x,y); adapt.first (i,z); adapt.first (N-i,conj(z)); }
    template<class V> inline void second(index_type i, const V &x, const V &y) const { value_type z(x,y); adapt.second(i,z); adapt.second(N-i,conj(z)); }
  
    template<int M,class V> inline void first (index_type i, const Vector<simd_vector_generator<M,V> > &x, const Vector<simd_vector_generator<M,V> > &y) const { Vector<simd_vector_generator<M/2,value_type> > z; z=complexlo(x,y); adapt.first (i,z); adapt.first (N-i,conj(z)); z=complexhi(x,y); adapt.first (N+i,z); adapt.first (2*N-i,conj(z)); }
    template<int M,class V> inline void second(index_type i, const Vector<simd_vector_generator<M,V> > &x, const Vector<simd_vector_generator<M,V> > &y) const { Vector<simd_vector_generator<M/2,value_type> > z; z=complexlo(x,y); adapt.second(i,z); adapt.second(N-i,conj(z)); z=complexhi(x,y); adapt.second(N+i,z); adapt.second(2*N-i,conj(z)); }
};


template <int N> struct aux_fft_adaptor;
struct make_fft_adaptor
{
  template<int N> struct rebind { typedef aux_fft_adaptor<N> other; };
  template<class G> inline make_fft_adaptor operator()(Vector<G> & ) const { return make_fft_adaptor(); }
  template<class G> inline make_fft_adaptor operator()(Matrix<G> & ) const { return make_fft_adaptor(); }
  template<class G> inline int normalization(const Vector<G> &) const { return 1; }
};
template<int N> struct aux_fft_adaptor
{
  inline aux_fft_adaptor(const make_fft_adaptor &) {}

  template<class G> struct rebind           { typedef                         fft_adaptor<Vector<G> >   other; };
  template<class G> struct pair_rebind      { typedef pair_fft_adaptor<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct real_rebind      { typedef real_fft_adaptor<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct simd_rebind      { typedef                    simd_fft_adaptor<Vector<G> >   other; };
  template<class G> struct pair_simd_rebind { typedef pair_fft_adaptor<N,simd_fft_adaptor<Vector<G> > > other; };

  template<class G> inline typename rebind          <G>::other make0         (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x,Vector<G> & ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename real_rebind     <G>::other real_make     (Vector<G> &x             ) const { return typename real_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_rebind     <G>::other pair_make     (Vector<G> &x,Vector<G> &y) const { return typename pair_rebind     <G>::other(x,y); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x             ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x,Vector<G> & ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_simd_rebind<G>::other simd_pair_make(Vector<G> &x,Vector<G> &y) const { return typename pair_simd_rebind<G>::other(x,y); }
};


template<int N,class S>
class full_matrix_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef full_matrix_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;
    value_type *p1;
    value_type *p2;

  public:
    inline full_matrix_fft_adaptor(array_type &s, value_type *px, value_type *py) : adapt(s),p1(px),p2(py) {}

    inline array_type &array()       { return adapt.array(); }
    inline array_type &array() const { return adapt.array(); }

    //template<class V> inline void first (index_type i, const V &v) const
    //{
    //  value_type *p=&array()[i];
    //  adaptor_type::sfirst (p,v);
    //  adaptor_type::ssecond_ru(p2-(p-p1),conj(v));
    //}
    //template<class V> inline void second(index_type i, const V &v) const
    //{
    //  value_type *p=&array()[i];
    //  adaptor_type::ssecond(p,v);
    //  value_type *q=p2-(p-p1);
    //  adaptor_type::sfirst_ru(q,conj(v));
    //}
    
    inline void first (index_type i, const value_type &v) const
    {
      value_type *p=&array()[i];
      adaptor_type::sfirst (p,v);
      adaptor_type::ssecond_ru(p2-(p-p1),conj(v));
    }
    inline void second(index_type i, const value_type &v) const
    {
      value_type *p=&array()[i];
      adaptor_type::ssecond(p,v);
      value_type *q=p2-(p-p1);
      adaptor_type::sfirst_ru(q,conj(v));
    }

 #ifdef NDEBUG
    inline void first0(              const value_type &v) const
    {
      value_type *p=&array()[0]; adaptor_type::sfirst (p,v);
      p=&array()[(p==p1)*N];     adaptor_type::ssecond_ru(p2-(p-p1),conj(v));
    }
    template<class V> inline void first0 (              const Vector<simd_vector_generator<2,V> > &a) const
    {
      value_type *p;
      Vector<simd_vector_generator<2,V> > x=conj(a);
      p=&array()[0]; adaptor_type::sfirst (p,a[0]); p=&array()[(p==p1)*N]; adaptor_type::ssecond_ru(p2-(p-p1),x[0]);
      p=&array()[1]; adaptor_type::sfirst (p,a[1]);                        adaptor_type::ssecond_ru(p2-(p-p1),x[1]);
    }
 #else
    inline void first0(              const value_type &v) const
    {
      value_type *p=&array()[0];      adaptor_type::sfirst (p,v);
      p+=((p==p1)*N)*(&array()[1]-p); adaptor_type::ssecond_ru(p2-(p-p1),conj(v));
    }
    template<class V> inline void first0 (              const Vector<simd_vector_generator<2,V> > &a) const
    {
      value_type *p;
      Vector<simd_vector_generator<2,V> > x=conj(a);
      p=&array()[0]; adaptor_type::sfirst (p,a[0]); p+=((p==p1)*N)*(&array()[1]-p); adaptor_type::ssecond_ru(p2-(p-p1),x[0]);
      p=&array()[1]; adaptor_type::sfirst (p,a[1]);                                 adaptor_type::ssecond_ru(p2-(p-p1),x[1]);
    }
 #endif

    template<class V> inline void first  (index_type i, const Vector<simd_vector_generator<2,V> > &a) const
    {
      value_type *p;
      Vector<simd_vector_generator<2,V> > x=conj(a);
      p=&array()[i  ]; adaptor_type::sfirst (p,a[0]); adaptor_type::ssecond_ru(p2-(p-p1),x[0]);
      p=&array()[i+1]; adaptor_type::sfirst (p,a[1]); adaptor_type::ssecond_ru(p2-(p-p1),x[1]);
    }
    template<class V> inline void second (index_type i, const Vector<simd_vector_generator<2,V> > &a) const
    {
      value_type *p;
      Vector<simd_vector_generator<2,V> > x=conj(a);
      p=&array()[i  ]; adaptor_type::ssecond(p,a[0]); adaptor_type::sfirst_ru(p2-(p-p1),x[0]);
      p=&array()[i+1]; adaptor_type::ssecond(p,a[1]); adaptor_type::sfirst_ru(p2-(p-p1),x[1]);
    }
};

template<int N,class F,class T> struct aux_full_matrix_fft_adaptor;

template<class F,class T>
struct make_full_matrix_fft_adaptor : public array_traits<T>
{
  typedef array_traits<T> base;
  ARRAY_BASE_TYPES;

  value_type *p1;
  value_type *p2;

#ifdef NDEBUG
  inline make_full_matrix_fft_adaptor(array_type &x, array_type &y) : p1(&x[0]), p2(&*y.end()) {}
#else
  inline make_full_matrix_fft_adaptor(array_type &x, array_type &y) : p1(&x[0]) { p2=&y.back()+(&y[1]-&y[0]); }
#endif
  template<int N> struct rebind { typedef aux_full_matrix_fft_adaptor<N,F,T> other; };
  template<class G> inline int normalization(const Vector<G> &) const { return 1; }
};

template<int N,class F,class T>
struct aux_full_matrix_fft_adaptor : public array_traits<T>
{
  typedef array_traits<T> base;
  ARRAY_BASE_TYPES;

  value_type *p1;
  value_type *p2;

  inline aux_full_matrix_fft_adaptor(const make_full_matrix_fft_adaptor<F,T> &x) : p1(x.p1), p2(x.p2) {}

// a priori il faudrait appliquer un rebind sur l'adaptateur plutôt que d'utiliser les types fft_adaptor et simd_fft_adaptor
  template<class T2> struct rebind       { typedef full_matrix_fft_adaptor<N,     fft_adaptor<T2> > other; };
  template<class T2> struct simd_rebind  { typedef full_matrix_fft_adaptor<N,simd_fft_adaptor<T2> > other; };

  template<class G> inline typename rebind          <Vector<G> >::other make0    (Vector<G> &x             ) const { return typename rebind          <Vector<G> >::other(x,p1,p2); }\
  template<class G> inline typename rebind          <Vector<G> >::other make     (Vector<G> &x             ) const { return typename rebind          <Vector<G> >::other(x,p1,p2); }\
  template<class G> inline typename rebind          <Vector<G> >::other make     (Vector<G> &x,Vector<G> & ) const { return typename rebind          <Vector<G> >::other(x,p1,p2); }\
  template<class G> inline typename simd_rebind     <Vector<G> >::other simd_make(Vector<G> &x             ) const { return typename simd_rebind     <Vector<G> >::other(x,p1,p2); }\
  template<class G> inline typename simd_rebind     <Vector<G> >::other simd_make(Vector<G> &x,Vector<G> & ) const { return typename simd_rebind     <Vector<G> >::other(x,p1,p2); }
};



template<int N> struct aux_full_fft_adaptor;
struct make_full_fft_adaptor
{
  const void *p0;
  template<int N> struct rebind { typedef aux_full_fft_adaptor<N> other; };
  template<class G> inline make_full_fft_adaptor(Vector<G> &X) : p0(&X[0]) {}
  template<class G> inline make_full_fft_adaptor(Matrix<G> & ) : p0(NULL) {}
  template<class G> inline make_fft_adaptor operator()(Vector<G> &) const { return make_fft_adaptor(); }
  template<class G> inline make_fft_adaptor operator()(Matrix<G> &) const { return make_fft_adaptor(); }
  template<class G> inline make_full_matrix_fft_adaptor<make_fft_adaptor,Vector<G> > operator()(Vector<G> &X,Vector<G> &Y) const { return make_full_matrix_fft_adaptor<make_fft_adaptor,Vector<G> >(X,Y); }
};
template<int N> struct aux_full_fft_adaptor
{
  const void *p0;

  inline aux_full_fft_adaptor(const make_full_fft_adaptor &x) : p0(x.p0) {}

  template<class G> struct      rebind      { typedef                              fft_adaptor<Vector<G> >   other; };
  template<class G> struct pair_rebind      { typedef full_fft_adaptor     <N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct real_rebind      { typedef full_real_fft_adaptor<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct simd_rebind      { typedef                         simd_fft_adaptor<Vector<G> >   other; };
  template<class G> struct pair_simd_rebind { typedef full_fft_adaptor     <N,simd_fft_adaptor<Vector<G> > > other; };

  template<class G> inline typename rebind          <G>::other make0         (Vector<G> &x             ) const { return typename rebind          <G>::other(x     ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x             ) const { return typename rebind          <G>::other(x     ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x,Vector<G> & ) const { return typename rebind          <G>::other(x     ); }\
  template<class G> inline typename real_rebind     <G>::other real_make     (Vector<G> &x             ) const { return typename real_rebind     <G>::other(x     ); }\
  template<class G> inline typename pair_rebind     <G>::other pair_make     (Vector<G> &x,Vector<G> &y) const { return typename pair_rebind     <G>::other(x,y,p0); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x             ) const { return typename simd_rebind     <G>::other(x     ); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x,Vector<G> & ) const { return typename simd_rebind     <G>::other(x     ); }\
  template<class G> inline typename pair_simd_rebind<G>::other simd_pair_make(Vector<G> &x,Vector<G> &y) const { return typename pair_simd_rebind<G>::other(x,y,p0); }
};



template<int N,class S>
class in_place_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef in_place_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;

  public:
    inline in_place_fft_adaptor(array_type &x) : adapt(x) {}
    inline in_place_fft_adaptor(adaptor_type &s) : adapt(s) {}

    inline void first0(              const complex<value_type> &x) const { adapt.first0(     real(x)); }
    inline void first (index_type i, const complex<value_type> &x) const { adapt.second(N-i, imag(x)); adapt.first(i,  real(x)); }
    inline void second(index_type i, const complex<value_type> &x) const { adapt.second(i,  -imag(x)); adapt.first(N-i,real(x)); } //ordre important

    inline void first0(              const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.first0(real(x)); }
    inline void first (index_type i, const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.second(N-i,imag(x)); adapt.first(i,real(x)); }
    inline void second(index_type i, const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.first(i,-imag(x));  adapt.second(N-i,real(x)); }

    inline void first0(              const Vector<simd_vector_generator<2,complex<value_type> > > &x) const { adapt.second(N-1,imag(x[1])); adapt.first0(real(x)); }
    inline void first (index_type i, const Vector<simd_vector_generator<2,complex<value_type> > > &x) const { adapt.second_ru(N-i,imag(x)); adapt.first(i,real(x)); }
    inline void second(index_type i, const Vector<simd_vector_generator<2,complex<value_type> > > &x) const { adapt.first(i,-imag(x));  adapt.second_ru(N-i,real(x)); }
};

template<int N,class S>
class pair_in_place_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef pair_in_place_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type setx;
    adaptor_type sety;

  public:
    inline pair_in_place_fft_adaptor(array_type &x, array_type &y) : setx(x), sety(y) {}

    inline void first0(              const complex<value_type> &x) const { sety.first(N-1  ,imag(x)); setx.first0(   real(x)); }
    inline void first (index_type i, const complex<value_type> &x) const { sety.first(N-1-i,imag(x)); setx.first (i, real(x)); }
    inline void second(index_type i, const complex<value_type> &x) const { sety.first(N-1-i,real(x)); setx.first (i,-imag(x));  }

    inline void first0(              const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { sety.second(N-1  ,imag(x)); setx.first0(   real(x)); }
    inline void first (index_type i, const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { sety.second(N-i-1,imag(x)); setx.first (i, real(x)); }
    inline void second(index_type i, const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { sety.second(N-i-1,real(x)); setx.first (i,-imag(x)); }

    inline void first0(              const Vector<simd_vector_generator<2,complex<value_type> > > &x) const { sety.second_ru(N-1  ,imag(x)); setx.first0(   real(x)); }
    inline void first (index_type i, const Vector<simd_vector_generator<2,complex<value_type> > > &x) const { sety.second_ru(N-i-1,imag(x)); setx.first (i, real(x)); }
    inline void second(index_type i, const Vector<simd_vector_generator<2,complex<value_type> > > &x) const { sety.second_ru(N-i-1,real(x)); setx.first (i,-imag(x)); }
};

template<int N,class S>
class real_in_place_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef real_in_place_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;

  public:
    inline real_in_place_fft_adaptor(array_type &x) : adapt(x) {}

    inline void first0(              const value_type &x                     ) const { adapt.first0 (  x); }
    inline void firstN(              const value_type &x                     ) const { adapt.first  (N/2,x); }
    inline void first (index_type i, const value_type &x                     ) const { adapt.first  (i,x); }
    inline void second(index_type i, const value_type &x                     ) const { adapt.second (i,x); }
    inline void first (index_type i, const value_type &x, const value_type &y) const { adapt.first  (i,x); adapt.first(N-i,y); }
    inline void second(index_type i, const value_type &x, const value_type &y) const { }

    inline void firstN(              const complex<value_type>             &z) const { firstN(  real(z)        ); }
    inline void first (index_type i, const complex<value_type>             &z) const { first (i,real(z),imag(z)); }
    inline void second(index_type i, const complex<value_type>             &z) const { second(i,real(z),imag(z)); }

    inline void first0(              const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.first0(    real(x)); }
    inline void firstN(              const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.first (N/2,real(x)); }
    inline void first (index_type i, const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.first (i  ,real(x)); adapt.second(N-i,imag(x)); }
    inline void second(index_type i, const Vector<simd_vector_generator<1,complex<value_type> > > &x) const { adapt.first (N-i,real(x)); adapt.second(i  ,-imag(x)); }

    template<int M> inline void first0 (              const Vector<simd_vector_generator<M,complex<value_type> > > &z) const { adapt.first0 (        real(z)); adapt.second(N-M+1,imag(z[1])); }
    template<int M> inline void firstN (              const Vector<simd_vector_generator<M,complex<value_type> > > &z) const { adapt.first_u(N/2-M+1,real(z)); adapt.second(N/2+M-1,-imag(z[0])); }
    template<int M> inline void secondN(              const Vector<simd_vector_generator<M,complex<value_type> > > &z) const { adapt.first_u(N/2-M+1  ,real(z)); adapt.second   (N/2,-flip(imag(z))); }
    template<int M> inline void first  (index_type i, const Vector<simd_vector_generator<M,complex<value_type> > > &z) const { adapt.first  (i        ,real(z)); adapt.second_ru(N-i,imag(z)); }
    template<int M> inline void second (index_type i, const Vector<simd_vector_generator<M,complex<value_type> > > &z) const { adapt.first_u((N-M+1)-i,real(z)); adapt.second   (i  ,-flip(imag(z))); }
};

template<int N> struct aux_in_place_fft_adaptor;
struct make_in_place_fft_adaptor { template<int N> struct rebind { typedef aux_in_place_fft_adaptor<N> other; }; };
template<int N> struct aux_in_place_fft_adaptor
{
  inline aux_in_place_fft_adaptor(const make_in_place_fft_adaptor &) {}

  template<class G> struct      rebind      { typedef      in_place_fft_adaptor<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct pair_rebind      { typedef pair_in_place_fft_adaptor<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct real_rebind      { typedef real_in_place_fft_adaptor<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct simd_rebind      { typedef      in_place_fft_adaptor<N,simd_fft_adaptor<Vector<G> > > other; };
  template<class G> struct pair_simd_rebind { typedef pair_in_place_fft_adaptor<N,simd_fft_adaptor<Vector<G> > > other; };

  template<class G> inline typename rebind          <G>::other make0         (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x,Vector<G> & ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename real_rebind     <G>::other real_make     (Vector<G> &x             ) const { return typename real_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_rebind     <G>::other pair_make     (Vector<G> &x,Vector<G> &y) const { return typename pair_rebind     <G>::other(x,y); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x             ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x,Vector<G> & ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_simd_rebind<G>::other simd_pair_make(Vector<G> &x,Vector<G> &y) const { return typename pair_simd_rebind<G>::other(x,y); }
};


struct fft_norm_function
{
  template<class V> inline typename norm_function       <V>::const_reference operator()(const V &x           ) const { return norm(x); }
  template<class V> inline typename norm_binary_function<V>::const_reference operator()(const V &x,const V &y) const { return norm(x,y); }
};

struct fft_abs_function
{
  template<class V> inline typename abs_function        <V>::const_reference operator()(const V &x           ) const { return abs (x); }
  template<class V> inline typename abs_binary_function <V>::const_reference operator()(const V &x,const V &y) const { return abs (x,y); }
};

template<int N,class F,class S>
class function_fft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef F function_type;
    typedef function_fft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;
    function_type func;

  public:
    inline explicit function_fft_adaptor(array_type &x) : adapt(x),func() {}
    inline function_fft_adaptor(array_type &x, array_type &y) : adapt(x,y),func() {}
    template<class A> inline function_fft_adaptor(array_type &x, array_type &y, const A &a) : adapt(x,y,a),func() {}
    inline explicit function_fft_adaptor(adaptor_type &s) : adapt(s),func() {}

    template<class T> inline void first0 (              const T &x) const { adapt.first0 (  func(x)); }
    template<class T> inline void firstN (              const T &x) const { adapt.firstN (  func(x)); }
    template<class T> inline void secondN(              const T &x) const { adapt.secondN(  func(x)); }
    template<class T> inline void first  (index_type i, const T &x) const { adapt.first  (i,func(x)); }
    template<class T> inline void second (index_type i, const T &x) const { adapt.second (i,func(x)); }

    template<class T> inline void first0 (              const T &x, const T &y) const { adapt.first0 (  func(x,y)); }
    template<class T> inline void firstN (              const T &x, const T &y) const { adapt.firstN (  func(x,y)); }
    template<class T> inline void secondN(              const T &x, const T &y) const { adapt.secondN(  func(x,y)); }
    template<class T> inline void first  (index_type i, const T &x, const T &y) const { adapt.first  (i,func(x,y)); }
    template<class T> inline void second (index_type i, const T &x, const T &y) const { adapt.second (i,func(x,y)); }
};


template<int N,class F> struct aux_function_fft_adaptor;
template<      class F> struct make_function_fft_adaptor
{
  template<class G> inline int normalization(const Vector<G> &) const { return 1; }
  template<int N> struct rebind { typedef aux_function_fft_adaptor<N,F> other; };
};
template<int N,class F> struct aux_function_fft_adaptor
{
  inline aux_function_fft_adaptor(const make_function_fft_adaptor<F> &) {}

  template<class G> struct      rebind      { typedef function_fft_adaptor<N,F,                        fft_adaptor<Vector<G> >   > other; };
  template<class G> struct pair_rebind      { typedef function_fft_adaptor<N,F,pair_fft_adaptor<N,     fft_adaptor<Vector<G> > > > other; };
  template<class G> struct real_rebind      { typedef function_fft_adaptor<N,F,real_fft_adaptor<N,     fft_adaptor<Vector<G> > > > other; };
  template<class G> struct simd_rebind      { typedef function_fft_adaptor<N,F,                   simd_fft_adaptor<Vector<G> >   > other; };
  template<class G> struct pair_simd_rebind { typedef function_fft_adaptor<N,F,pair_fft_adaptor<N,simd_fft_adaptor<Vector<G> > > > other; };

  template<class G> inline typename rebind          <G>::other make0         (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x,Vector<G> & ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename real_rebind     <G>::other real_make     (Vector<G> &x             ) const { return typename real_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_rebind     <G>::other pair_make     (Vector<G> &x,Vector<G> &y) const { return typename pair_rebind     <G>::other(x,y); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x             ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x,Vector<G> & ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_simd_rebind<G>::other simd_pair_make(Vector<G> &x,Vector<G> &y) const { return typename pair_simd_rebind<G>::other(x,y); }
};


template<int N,class F> struct aux_function_full_fft_adaptor;
template<      class F> struct make_function_full_fft_adaptor
{
  const void *p0;
  template<class G> inline make_function_full_fft_adaptor(Vector<G> &X) : p0(&X[0]  ) {}
  template<class G> inline make_function_full_fft_adaptor(Matrix<G> &X) : p0(&X(0,0)) {}
  template<int N> struct rebind { typedef aux_function_full_fft_adaptor<N,F> other; };
};
template<int N,class F> struct aux_function_full_fft_adaptor
{
  const void *p0;
  inline aux_function_full_fft_adaptor(const make_function_full_fft_adaptor<F> &x) : p0(x.p0) {}

  template<class G> struct      rebind      { typedef function_fft_adaptor<N,F,                             fft_adaptor<Vector<G> >   > other; };
  template<class G> struct pair_rebind      { typedef function_fft_adaptor<N,F,     full_fft_adaptor<N,     fft_adaptor<Vector<G> > > > other; };
  template<class G> struct real_rebind      { typedef function_fft_adaptor<N,F,full_real_fft_adaptor<N,     fft_adaptor<Vector<G> > > > other; };
  template<class G> struct simd_rebind      { typedef function_fft_adaptor<N,F,                        simd_fft_adaptor<Vector<G> >   > other; };
  template<class G> struct pair_simd_rebind { typedef function_fft_adaptor<N,F,     full_fft_adaptor<N,simd_fft_adaptor<Vector<G> > > > other; };

  template<class G> inline typename rebind          <G>::other make0         (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x             ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename rebind          <G>::other make          (Vector<G> &x,Vector<G> & ) const { return typename rebind          <G>::other(x  ); }\
  template<class G> inline typename real_rebind     <G>::other real_make     (Vector<G> &x             ) const { return typename real_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_rebind     <G>::other pair_make     (Vector<G> &x,Vector<G> &y) const { return typename pair_rebind     <G>::other(x,y,p0); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x             ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename simd_rebind     <G>::other simd_make     (Vector<G> &x,Vector<G> & ) const { return typename simd_rebind     <G>::other(x  ); }\
  template<class G> inline typename pair_simd_rebind<G>::other simd_pair_make(Vector<G> &x,Vector<G> &y) const { return typename pair_simd_rebind<G>::other(x,y,p0); }
};

typedef make_function_fft_adaptor     <fft_abs_function > make_abs_fft_adaptor;
typedef make_function_full_fft_adaptor<fft_abs_function > make_abs_full_fft_adaptor;
typedef make_function_fft_adaptor     <fft_norm_function> make_norm_fft_adaptor;
typedef make_function_full_fft_adaptor<fft_norm_function> make_norm_full_fft_adaptor;


template<int N,class S>
class ifft_adaptor0 : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef ifft_adaptor0 self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES
    typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
    typedef typename complex_traits<complex_type>::value_type real_type;

  public:
    adaptor_type adapt;

  public:
    inline explicit ifft_adaptor0(array_type &x) : adapt(x) {}
    inline ifft_adaptor0(array_type &x, array_type &y) : adapt(x,y) {}
    inline explicit ifft_adaptor0(adaptor_type &s) : adapt(s) {}

    inline void first0(const value_type &x) const { adapt.first0   (    x*(1/real_type(N))); }
    template<class T> inline void first0(const Vector<simd_vector_generator<1,T> > &x) const { adapt.first0(x*(1/real_type(N))); }
    template<class T> inline void first0(const Vector<simd_vector_generator<2,T> > &x) const { Vector<simd_vector_generator<2,T> > y=x*(1/real_type(N)); adapt.first0(y[0]); adapt.second(N-1,y[1]); }

    inline void first (index_type i, const value_type &x) const { adapt.second_ru(N-i,x*(1/real_type(N))); }
    inline void second(index_type i, const value_type &x) const { adapt.second_ru(N-i,x*(1/real_type(N))); }
    template<int M,class T> inline void first (index_type i, const Vector<simd_vector_generator<M,T> > &x) const { adapt.second_ru(N-i,x*(1/real_type(N))); }
    template<int M,class T> inline void second(index_type i, const Vector<simd_vector_generator<M,T> > &x) const { adapt.second_ru(N-i,x*(1/real_type(N))); }
};

template<int N,class S>
class ifft_adaptor1 : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef ifft_adaptor1 self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;

  public:
    inline explicit ifft_adaptor1(array_type &x) : adapt(x) {}
    inline ifft_adaptor1(array_type &x, array_type &y) : adapt(x,y) {}
    inline explicit ifft_adaptor1(adaptor_type &s) : adapt(s) {}

    inline void first0(const value_type &x) const { adapt.first0   (    x); }
    template<class T> inline void first0(const Vector<simd_vector_generator<1,T> > &x) const { adapt.first0(x); }
    template<class T> inline void first0(const Vector<simd_vector_generator<2,T> > &x) const { adapt.first0(x[0]); adapt.second(N-1,x[1]); }

    inline void first (index_type i, const value_type &x) const { adapt.second_ru(N-i,x); }
    inline void second(index_type i, const value_type &x) const { adapt.second_ru(N-i,x); }
    template<int M,class T> inline void first (index_type i, const Vector<simd_vector_generator<M,T> > &x) const { adapt.second_ru(N-i,x); }
    template<int M,class T> inline void second(index_type i, const Vector<simd_vector_generator<M,T> > &x) const { adapt.second_ru(N-i,x); }
};

template<int N,class S>
class ifft_adaptor : public array_traits<typename S::array_type>
{
  public:
    typedef S adaptor_type;
    typedef ifft_adaptor self;
    typedef array_traits<typename S::array_type> base;
    ARRAY_BASE_TYPES

  public:
    adaptor_type adapt;
    const void *p0;

  public:
    inline explicit ifft_adaptor(array_type &x, const void *p) : adapt(x), p0(p) {}
    inline ifft_adaptor(array_type &x, array_type &y, const void *p) : adapt(x,y), p0(p) {}
    inline explicit ifft_adaptor(adaptor_type &s, const void *p) : adapt(s), p0(p) {}

    inline array_type &array()       { return adapt.array(); }
    inline array_type &array() const { return adapt.array(); }

#ifdef NDEBUG
                      inline void first0(const value_type                          &x) const { value_type *p=&array()[-1]; adapt.second_ru((p0!=p)*N-1,x); }
    template<class T> inline void first0(const Vector<simd_vector_generator<1,T> > &x) const { value_type *p=&array()[-1]; adapt.second_ru((p0!=p)*N-1,x); }
    template<class T> inline void first0(const Vector<simd_vector_generator<2,T> > &x) const { value_type *p=&array()[-1]; if (p0!=p) adapt.second_ru(N-1,x); else adapt.second_ru(-1,N-1,x); }
#else
    template<class T> inline void first0(const T &x) const { adapt.second_ru(N-1,x); }
#endif
    inline void first  (index_type i, const value_type &x) const { adapt.second_ru(N-1-i,x); }
    inline void second (index_type i, const value_type &x) const { adapt.second_ru(N-1-i,x); }
    template<int M,class T> inline void first (index_type i, const Vector<simd_vector_generator<M,T> > &x) const { adapt.second_ru(N-1-i,x); }
    template<int M,class T> inline void second(index_type i, const Vector<simd_vector_generator<M,T> > &x) const { adapt.second_ru(N-1-i,x); }
};

template<int N> struct aux_ifft_adaptor;

struct make_ifft_adaptor
{
  const void *p0;

  template<int N> struct rebind { typedef aux_ifft_adaptor<N> other; };
  template<class G> inline make_ifft_adaptor(Vector<G> &X) : p0(&X[0]  ) {}
  template<class G> inline make_ifft_adaptor(Matrix<G> &X) : p0(&X(0,0)) {}
  template<class G> inline int normalization(const Vector<G> &X) const { return X.size(); }
  template<class G> inline make_ifft_adaptor operator()(Vector<G> &x) const { return make_ifft_adaptor(x); }
  template<class G> inline make_ifft_adaptor operator()(Matrix<G> &x) const { return make_ifft_adaptor(x); }
};
template<int N> struct aux_ifft_adaptor
{
  const void *p0;

  inline aux_ifft_adaptor(const make_ifft_adaptor &x) : p0(x.p0) {}
  template<class G> struct rebind0     { typedef ifft_adaptor0<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct rebind1     { typedef ifft_adaptor1<N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct rebind      { typedef ifft_adaptor <N,     fft_adaptor<Vector<G> > > other; };
  template<class G> struct simd_rebind { typedef ifft_adaptor <N,simd_fft_adaptor<Vector<G> > > other; };

  template<class G> inline typename rebind0     <G>::other make0    (Vector<G> &x             ) const { return typename rebind0    <G>::other(x); }
  template<class G> inline typename rebind1     <G>::other make     (Vector<G> &x             ) const { return typename rebind1    <G>::other(x); }
  template<class G> inline typename rebind      <G>::other make     (Vector<G> & ,Vector<G> &y) const { return typename rebind     <G>::other(y,p0); }
  template<class G> inline typename simd_rebind <G>::other simd_make(Vector<G> & ,Vector<G> &y) const { return typename simd_rebind<G>::other(y,p0); }
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class G1,class S>
inline void complex_fft_1(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;

  value_type t0;
  t0=X[0];
	Y.first0(  t0);
}

template<class G1,class S>
inline void complex_fft_2(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;

  value_type t0,t1;
  t0=X[0]; t1=X[1];
	Y.first0(  t0+t1);
	Y.second(1,t0-t1);
}

template<class G1,class S>
inline void complex_fft_3(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3;
  {
    value_type v0, v1;
    t0=X[0];
    v0=X[1]; v1=X[2]; t1=v0+v1;
    t2=t0-real_type(0.5)*t1;
    t3=imul(SIN01_03,v1-v0);
  }
	Y.first0(  t0+t1);
	Y.first (1,t2+t3);
	Y.second(2,t2-t3);
}

template<class G1,class S>
inline void complex_fft_4(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3;
  { value_type v0,v1; v0=X[0]; v1=X[2]; t0=v0+v1; t1=v0-v1; }
  { value_type v0,v1; v0=X[1]; v1=X[3]; t2=v0+v1; t3=imul(v0-v1); }

  Y.first0(  t0+t2); Y.second(2,t0-t2);
  Y.first (1,t1-t3); Y.second(3,t1+t3);
}

template<class G1,class S>
inline void complex_fft_5(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3,t4;

  t0=X[0];
  {
    value_type v0,v1;
    {
      value_type v2,v3;
      v2=X[1]; v3=X[4]; v0=v2+v3; t2=imul(v2-v3);
      v2=X[2]; v3=X[3]; v1=v2+v3; t4=imul(v2-v3);
    }
    t1=v0+v1; t3=real_type(SQRT_5_16)*(v0-v1);
  }
	Y.first0(t0+t1);
	t0-=real_type(0.25)*t1;
	{ value_type v0,v1; v0=SIN02_05*t2+SIN01_05*t4; v1=t0+t3; Y.first (1,v1-v0); Y.second(4,v1+v0); }
	{ value_type v0,v1; v0=SIN01_05*t2-SIN02_05*t4; v1=t0-t3; Y.first (2,v1-v0); Y.second(3,v1+v0); }
}

template<class G1,class S>
inline void complex_fft_8(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3,t4,t5,t6,t7;

  { value_type v0,v1; v0=X[0],v1=X[4]; t0=v0+v1; t1=     v0-v1 ; }
  { value_type v0,v1; v0=X[2],v1=X[6]; t2=v0+v1; t3=imul(v0-v1); }
  { value_type v0,v1; v0=X[7],v1=X[3]; t4=v0+v1; v0=v0-v1; t5=add_imul(v0,v0); }
  { value_type v0,v1; v0=X[1],v1=X[5]; t6=v0+v1; v0=v0-v1; t7=sub_imul(v0,v0); }

  { value_type v0,v1; v0=     t0+t2 ; v1=t6+t4; Y.first0(  v0+v1); Y.second(4,v0-v1); }
  { value_type v0,v1; v0=imul(t4-t6); v1=t0-t2; Y.first (2,v1+v0); Y.second(6,v1-v0); }
  { value_type v0,v1; v0=     t1-t3 ; v1=     SIN16*(t7+t5); Y.first (1,v0+v1); Y.second(5,v0-v1); }
  { value_type v0,v1; v0=     t1+t3 ; v1=imul(SIN16, t5-t7); Y.first (3,v0+v1); Y.second(7,v0-v1); }
}

template<class G,class S>
inline void complex_fft_16(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15;
  { value_type v0,v1; v0=X[ 0]; v1=X[ 8]; t01=v0+v1; t03=v0-v1; v0=X[ 4]; v1=X[12]; t02=v0+v1; v0=imul(v0-v1); t00=t01+t02; t01=     t01-t02 ; t02=t03+v0; t03-=v0; }
  { value_type v0,v1; v0=X[15]; v1=X[ 7]; t05=v0+v1; t07=v0-v1; v0=X[ 3]; v1=X[11]; t06=v0+v1; v0=imul(v0-v1); t04=t05+t06; t05=     t05-t06 ; t06=t07+v0; t07-=v0; }
  { value_type v0,v1; v0=X[ 6]; v1=X[14]; t09=v0+v1; t11=v1-v0; v0=X[10]; v1=X[ 2]; t10=v0+v1; v0=     v0-v1 ; t08=t09+t10; t09=imul(t09-t10); t10=t11-v0; t11=imul(t11+v0); }
  { value_type v0,v1; v0=X[ 1]; v1=X[ 9]; t13=v0+v1; t15=v0-v1; v0=X[ 5]; v1=X[13]; t14=v0+v1; v0=imul(v0-v1); t12=t13+t14; t13=     t13-t14 ; t14=t15+v0; t15-=v0; }

  { value_type v0,v1; v0=t00+t08;           v1=t12+t04; Y.first0(   v0+v1); Y.second( 8,v0-v1); }
  { value_type v0,v1; v0=imul(t04-t12);     v1=t00-t08; Y.first ( 4,v1+v0); Y.second(12,v1-v0); }
  { value_type v0,v1; v0=sub_imul(t13,t13); v1=add_imul(t05,t05); v1=SIN16*(v0+v1); v0=t01+t09; Y.first( 2,v0+v1); Y.second(10,v0-v1); }
  { value_type v0,v1; v0=add_imul(t13,t13); v1=imul(t05)-t05;     v1=SIN16*(v1-v0); v0=t01-t09; Y.first( 6,v0+v1); Y.second(14,v0-v1); }
  { value_type v0,v1,v2,v3; v0=real_type( SIN08)*t14-imul(SIN24,t14); v1=real_type( SIN08)*t06+imul(SIN24,t06); v2=v1+v0; v3=imul(v1-v0); { value_type v4; v4=SIN16*(t11-t10); v0=t02+v4; v1=t02-v4; } Y.first (3 ,v0+v2); Y.second(11,v0-v2); Y.first(7,v1+v3); Y.second(15,v1-v3); }
  { value_type v0,v1,v2,v3; v0=SIN24*t15-imul(real_type( SIN08),t15); v1=SIN24*t07+imul(real_type( SIN08),t07); v2=v1+v0; v3=imul(v1-v0); { value_type v4; v4=SIN16*(t11+t10); v0=t03+v4; v1=t03-v4; } Y.first (1 ,v0+v2); Y.second(9 ,v0-v2); Y.first(5 ,v1+v3);Y.second(13,v1-v3); }
}

template<class G,class S>
inline void complex_fft_32(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33;

  { value_type v0,v1,v2,v3;    { value_type v4,v5; v4=X[ 0]; v5=X[16]; v0=v4+v5; v1=v4-v5; v4=X[ 8]; v5=X[24]; v2=v4+v5; v3=imul(v4-v5); } t00=v0+v2; t01=v0-v2; t02=v1+v3; t03=v1-v3; }
  { value_type v0,v1,v2,v3;    { value_type v4,v5; v4=X[ 4]; v5=X[20]; v0=v4+v5; v1=v4-v5; v4=X[28]; v5=X[12]; v2=v4+v5; v3=     v4-v5 ; } t04=v0+v2; t05=imul(v0-v2); v0=v3+v1; v2=imul(v3-v1); t06=SIN16*(v2+v0); t07=SIN16*(v2-v0); }
  { value_type v0,v1,v2,v3;    { value_type v4,v5; v4=X[30]; v5=X[14]; v0=v4+v5; v1=v4-v5; v4=X[ 6]; v5=X[22]; v2=v4+v5; v3=imul(v4-v5); } t12=v0+v2; v0=v0-v2; t15=add_imul(v0,v0); v0=v1+v3; v2=v1-v3; t13=SIN08*v0+imul(SIN24,v0); t14=SIN24*v2+imul(SIN08,v2); }
  { value_type v0,v1,v2,v3;    { value_type v4,v5; v4=X[ 2]; v5=X[18]; v0=v4+v5; v1=v4-v5; v4=X[10]; v5=X[26]; v2=v4+v5; v3=imul(v4-v5); } t08=v0+v2; v0=v0-v2; t11=sub_imul(v0,v0); v0=v1+v3; v2=v1-v3; t09=SIN08*v0-imul(SIN24,v0); t10=SIN24*v2-imul(SIN08,v2); }
  { value_type v0,v1,v2,v3,v6; { value_type v4,v5; v4=X[ 1]; v5=X[17]; v0=v4+v5; v1=v4-v5; v4=X[ 9]; v5=X[25]; v2=v4+v5; v3=imul(v4-v5); } t25=v1+v3; t26=v1-v3; t27=v0+v2; v6=v0-v2; { value_type v4,v5; v4=X[ 5]; v5=X[21]; v0=v4+v5; v1=v4-v5; v4=X[29]; v5=X[13]; v2=v4+v5; v3=     v4-v5 ; } t28=v0+v2; t29=t27-t28; v2=imul(v2-v0); t30=v6+v2; t31=v6-v2; v6=v3+v1; v2=imul(v3-v1); t32=SIN16*(v2+v6); t33=SIN16*(v2-v6); }
  { value_type v0,v1,v2,v3,v6; { value_type v4,v5; v4=X[31]; v5=X[15]; v0=v4+v5; v1=v4-v5; v4=X[ 7]; v5=X[23]; v2=v4+v5; v3=imul(v4-v5); } t16=v1+v3; t17=v1-v3; t18=v0+v2; v6=v0-v2; { value_type v4,v5; v4=X[ 3]; v5=X[19]; v0=v4+v5; v1=v4-v5; v4=X[27]; v5=X[11]; v2=v4+v5; v3=     v4-v5 ; } t19=v0+v2; v0=imul(v0-v2); t20=t18-t19; t21=v6+v0; t22=v6-v0; v6=v3+v1; v2=imul(v3-v1); t23=SIN16*(v2+v6); t24=SIN16*(v2-v6); }

  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t00+t04; v5=                          t12+t08 ; v0=v4+v5; v1=v4-v5; v4=t27+t28;                                                              v5=t18+t19;                                                             v2=                     v5+v4 ; v3=imul(v5-v4);                 } 
  Y.first0(  v0+v2); 
  Y.second(16,v0-v2); 
  Y.first( 8,v1+v3); 
  Y.second(24,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t00-t04; v5=imul(                     t12-t08); v0=v4+v5; v1=v4-v5; v4=t20+t29;                                                              v5=imul(t20-t29);                                                       v2=SIN16*(v4+v5); v3=SIN16*(v5-v4); } Y.first( 4,v0+v2); Y.second(20,v0-v2); Y.first(12,v1+v3); Y.second(28,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t01-t05; v5=     SIN16*(t11+t15); v0=v4+v5; v1=v4-v5;             v4=SIN24*t30-imul(SIN08,t30); v5=SIN24*t22+imul(SIN08,t22);              v2=                     v4+v5 ; v3=imul(v5-v4);                 } Y.first( 2,v0+v2); Y.second(18,v0-v2); Y.first(10,v1+v3); Y.second(26,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t01+t05; v5=imul(SIN16, t11-t15); v1=v4+v5; v0=v4-v5;             v4=SIN08*t31-imul(SIN24,t31); v5=SIN08*t21+imul(SIN24,t21);              v2=                     v4+v5 ; v3=imul(v5-v4);                 } Y.first( 6,v0+v2); Y.second(22,v0-v2); Y.first(14,v1+v3); Y.second(30,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t03+t06; v5=                          t14+t10 ; v0=v4+v5; v1=v4-v5; v4=t26+t32; v4=SIN28*v4-imul(SIN04,v4) ; v5=t17+t23; v5=SIN28*v5+imul(SIN04,v5); v2=                     v5+v4 ; v3=imul(v5-v4);                 } Y.first( 1,v0+v2); Y.second(17,v0-v2); Y.first( 9,v1+v3); Y.second(25,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t03-t06; v5=imul(                     t14-t10); v0=v4+v5; v1=v4-v5; v4=t26-t32; v4=SIN12*v4-imul(SIN20,v4) ; v5=t17-t23; v5=SIN12*v5+imul(SIN20,v5); v2=                     v4+v5 ; v3=imul(v5-v4);                 } Y.first( 5,v0+v2); Y.second(21,v0-v2); Y.first(13,v1+v3); Y.second(29,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t02+t07; v5=                          t13+t09 ; v0=v4+v5; v1=v4-v5; v4=t25+t33; v4=SIN20*v4-imul(SIN12,v4) ; v5=t16+t24; v5=SIN20*v5+imul(SIN12,v5); v2=                     v5+v4 ; v3=imul(v5-v4);                 } Y.first( 3,v0+v2); Y.second(19,v0-v2); Y.first(11,v1+v3); Y.second(27,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t02-t07; v5=imul(                     t13-t09); v0=v4+v5; v1=v4-v5; v4=t25-t33; v4=SIN04*v4-imul(SIN28,v4) ; v5=t16-t24; v5=SIN04*v5+imul(SIN28,v5); v2=                     v4+v5 ; v3=imul(v5-v4);                 } Y.first( 7,v0+v2); Y.second(23,v0-v2); Y.first(15,v1+v3); Y.second(31,v1-v3); }
}

template<class G,class S>
inline void complex_fft_64(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63;

  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 0]; u1=X[32]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[16]; u1=X[48]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[ 8]; u1=X[40]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[24]; u1=X[56]; v6=u0+u1; v7=u1-u0; } { t00=add_imul(v1,v3); t01=sub_imul(v1,v3); } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t02=u0+u1; t03=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v6-v4); t04=u0-u1; t05=u0+u1; } { value_type u0,u1; u0=add_imul(v7,v7); u1=sub_imul(v5,v5); t06=imul(SIN16,u0-u1); t07=SIN16*(u1+u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 2]; u1=X[34]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[18]; u1=X[50]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[26]; u1=X[58]; v4=u0+u1; { value_type w0; w0=u1-u0; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[10]; u1=X[42]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t16=u0+u1; t17=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t18=u0-u1; t19=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t20=u0-u1; t21=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t22=u0-u1; t23=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 4]; u1=X[36]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[20]; u1=X[52]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[28]; u1=X[60]; v4=u0+u1; v5=u1-u0; } { value_type u0,u1; u0=X[12]; u1=X[44]; v6=u0+u1; v7=u0-u1; } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t08=u0+u1; t09=imul(u1-u0); } { value_type u0; u0=v0-v2; t10=sub_imul(u0,u0); } { value_type u0; u0=v4-v6; t11=add_imul(u0,u0); } { value_type u0; u0=add_imul(v1,v3); t12=SIN08*u0-imul(SIN24,u0); } { value_type u0; u0=sub_imul(v1,v3); t13=SIN24*u0-imul(SIN08,u0); } { value_type u0; u0=add_imul(v5,v7); t14=SIN08*u0+imul(SIN24,u0); } { value_type u0; u0=sub_imul(v5,v7); t15=SIN24*u0+imul(SIN08,u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[30]; u1=X[62]; v0=u0+u1; v1=u1-u0; } { value_type u0,u1; u0=X[14]; u1=X[46]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[22]; u1=X[54]; v4=u0+u1; { value_type w0; w0=u1-u0; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[ 6]; u1=X[38]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t32=u0+u1; t33=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t34=u0-u1; t35=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t36=u0-u1; t37=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t38=u0-u1; t39=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 1]; u1=X[33]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[17]; u1=X[49]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[25]; u1=X[57]; v4=u0+u1; { value_type w0; w0=u1-u0; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[ 9]; u1=X[41]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t40=u0+u1; t41=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t42=u0-u1; t43=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t44=u0-u1; t45=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t46=u0-u1; t47=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[31]; u1=X[63]; v0=u0+u1; v1=u1-u0; } { value_type u0,u1; u0=X[15]; u1=X[47]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[23]; u1=X[55]; v4=u0+u1; { value_type w0; w0=u1-u0; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[ 7]; u1=X[39]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t24=u0+u1; t25=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t26=u0-u1; t27=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t28=u0-u1; t29=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t30=u0-u1; t31=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; { value_type w0,w1; w0=X[ 5]; w1=X[37]; v0=w0+w1; u0=w0-w1; } { value_type w0,w1; w0=X[21]; w1=X[53]; v1=w0+w1; u1=w0-w1; } v2=add_imul(u0,u1); v3=sub_imul(u0,u1); } { value_type u0,u1; { value_type w0,w1; w0=X[29]; w1=X[61]; v4=w0+w1; u0=w1-w0; } { value_type w0,w1; w0=X[13]; w1=X[45]; v5=w0+w1; u1=w0-w1; } v6=add_imul(u0,u1); v7=sub_imul(u0,u1); } { value_type u0,u1; u0=v0+v1; u1=v4+v5; t48=u0+u1; t49=imul(u1-u0); } { value_type u0,u1; u0=v0-v1; u0=sub_imul(u0,u0); u1=v4-v5; u1=add_imul(u1,u1); t50=imul(SIN16,u1-u0); t51=SIN16*(u0+u1); } { value_type u0,u1; u0=SIN08*v6+imul(SIN24,v6); u1=SIN08*v2-imul(SIN24,v2); t53=u0+u1; t52=imul(u0-u1); } { value_type u0,u1; u0=SIN24*v3-imul(SIN08,v3); u1=SIN24*v7+imul(SIN08,v7); t55=u1+u0; t54=imul(u1-u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; { value_type w0,w1; w0=X[ 3]; w1=X[35]; v0=w0+w1; u0=w0-w1; } { value_type w0,w1; w0=X[19]; w1=X[51]; v1=w0+w1; u1=w0-w1; } v2=add_imul(u0,u1); v3=sub_imul(u0,u1); } { value_type u0,u1; { value_type w0,w1; w0=X[27]; w1=X[59]; v4=w0+w1; u0=w1-w0; } { value_type w0,w1; w0=X[11]; w1=X[43]; v5=w0+w1; u1=w0-w1; } v6=add_imul(u0,u1); v7=sub_imul(u0,u1); } { value_type u0,u1; u0=v0+v1; u1=v4+v5; t56=u0+u1; t57=imul(u1-u0); } { value_type u0,u1; u0=v0-v1; u0=sub_imul(u0,u0); u1=v4-v5; u1=add_imul(u1,u1); t58=imul(SIN16,u1-u0); t59=SIN16*(u1+u0); } { value_type u0,u1; u0=SIN08*v2-imul(SIN24,v2); u1=SIN08*v6+imul(SIN24,v6); t61=u0+u1; t60=imul(u1-u0); } { value_type u0,u1; u0=SIN24*v7+imul(SIN08,v7); u1=SIN24*v3-imul(SIN08,v3); t63=u0+u1; t62=imul(u0-u1); } }

  { value_type v0,v1,v2,v3; { value_type u0,u1; u0=t02+t08; u1=     t16+t32;  v0=u0+u1; v1=u0-u1; } { value_type u0,u1;u0=t40+t48; u1=t24+t56; v2=imul(u1-u0); v3=u0+u1; } Y.first0(v0+v3); Y.second(32,v0-v3);  Y.first(16,v1+v2); Y.second(48,v1-v2); }
  { value_type v0,v1,v2,v3; { value_type u0,u1; u0=t02-t08; u1=imul(t32-t16); v0=u0+u1; v1=u0-u1; } { value_type u0; u0=t40-t48; v2=sub_imul(u0,u0); } { value_type u0; u0=t24-t56; v3=add_imul(u0,u0); } { value_type u0; u0=SIN16*(v2+v3); Y.first( 8,v0+u0); Y.second(40,v0-u0); } { value_type u0; u0=imul(SIN16,v3-v2); Y.first(24,v1+u0); Y.second(56,v1-u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; v0=t03-t09; v1=t03+t09;                                           { value_type u0; u0=t41+t49; v2=SIN24*u0-imul(SIN08,u0); } { value_type u0; u0=t25+t57; v3=SIN24*u0+imul(SIN08,u0); } { value_type u0,u1; u0=sub_imul(t17,t17); u1=add_imul(t33,t33); v4=imul(SIN16,u1-u0); v5=SIN16*(u0+u1); }   { value_type u0; u0=t41-t49; v6=SIN08*u0-imul(SIN24,u0); } { value_type u0; u0=t25-t57; v7=SIN08*u0+imul(SIN24,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first(12,u0+u1); Y.second(44,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(28,u0+u1); Y.second(60,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 4,u0+u1); Y.second(36,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(20,u0+u1); Y.second(52,u0-u1); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0; u0=     SIN16*(t10+t11); v0=t05-u0; v1=t05+u0; } { value_type u0; u0=t43+t51; v2=SIN28*u0-imul(SIN04,u0); } { value_type u0; u0=t27+t59; v3=SIN28*u0+imul(SIN04,u0); } { value_type u0,u1; u0=SIN24*t19-imul(SIN08,t19); u1=SIN24*t35+imul(SIN08,t35); v4=imul(u1-u0); v5=u0+u1; } { value_type u0; u0=t43-t51; v6=SIN12*u0-imul(SIN20,u0); } { value_type u0; u0=t27-t59; v7=SIN12*u0+imul(SIN20,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first(10,u0+u1); Y.second(42,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(26,u0+u1); Y.second(58,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 2,u0+u1); Y.second(34,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(18,u0+u1); Y.second(50,u0-u1); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0; u0=imul(SIN16 ,t11-t10); v0=t04-u0; v1=t04+u0; } { value_type u0; u0=t42+t50; v2=SIN20*u0-imul(SIN12,u0); } { value_type u0; u0=t26+t58; v3=SIN20*u0+imul(SIN12,u0); } { value_type u0,u1; u0=SIN08*t18-imul(SIN24,t18); u1=SIN08*t34+imul(SIN24,t34); v4=imul(u1-u0); v5=u0+u1; } { value_type u0; u0=t42-t50; v6=SIN04*u0-imul(SIN28,u0); } { value_type u0; u0=t26-t58; v7=SIN04*u0+imul(SIN28,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first(14,u0+u1); Y.second(46,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(30,u0+u1); Y.second(62,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 6,u0+u1); Y.second(38,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(22,u0+u1); Y.second(54,u0-u1); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t00+t06; u1=     t12+t14 ; v0=u0-u1; v1=u0+u1; } { value_type u0; u0=t45+t53; v2=SIN26*u0-imul(SIN06,u0); } { value_type u0; u0=t29+t61; v3=SIN26*u0+imul(SIN06,u0); } { value_type u0,u1; u0=SIN20*t21-imul(SIN12,t21); u1=SIN20*t37+imul(SIN12,t37); v4=imul(u1-u0); v5=u1+u0; } { value_type u0; u0=t45-t53; v6=SIN10*u0-imul(SIN22,u0); } { value_type u0; u0=t29-t61; v7=SIN10*u0+imul(SIN22,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first(11,u0+u1); Y.second(43,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(27,u0+u1); Y.second(59,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 3,u0+u1); Y.second(35,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(19,u0+u1); Y.second(51,u0-u1); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t01+t07; u1=     t13+t15 ; v0=u0-u1; v1=u0+u1; } { value_type u0; u0=t47+t55; v2=SIN30*u0-imul(SIN02,u0); } { value_type u0; u0=t31+t63; v3=SIN30*u0+imul(SIN02,u0); } { value_type u0,u1; u0=SIN28*t23-imul(SIN04,t23); u1=SIN28*t39+imul(SIN04,t39); v4=imul(u1-u0); v5=u1+u0; } { value_type u0; u0=t47-t55; v6=SIN14*u0-imul(SIN18,u0); } { value_type u0; u0=t31-t63; v7=SIN14*u0+imul(SIN18,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first( 9,u0+u1); Y.second(41,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(25,u0+u1); Y.second(57,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 1,u0+u1); Y.second(33,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(17,u0+u1); Y.second(49,u0-u1); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t00-t06; u1=imul(t14-t12); v0=u0-u1; v1=u0+u1; } { value_type u0; u0=t44+t52; v2=SIN18*u0-imul(SIN14,u0); } { value_type u0; u0=t28+t60; v3=SIN18*u0+imul(SIN14,u0); } { value_type u0,u1; u0=SIN04*t20-imul(SIN28,t20); u1=SIN04*t36+imul(SIN28,t36); v4=imul(u1-u0); v5=u1+u0; } { value_type u0; u0=t44-t52; v6=SIN02*u0-imul(SIN30,u0); } { value_type u0; u0=t28-t60; v7=SIN02*u0+imul(SIN30,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first(15,u0+u1); Y.second(47,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(31,u0+u1); Y.second(63,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 7,u0+u1); Y.second(39,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(23,u0+u1); Y.second(55,u0-u1); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t01-t07; u1=imul(t15-t13); v0=u0-u1; v1=u0+u1; } { value_type u0; u0=t46+t54; v2=SIN22*u0-imul(SIN10,u0); } { value_type u0; u0=t30+t62; v3=SIN22*u0+imul(SIN10,u0); } { value_type u0,u1; u0=SIN12*t22-imul(SIN20,t22); u1=SIN12*t38+imul(SIN20,t38); v4=imul(u1-u0); v5=u1+u0; } { value_type u0; u0=t46-t54; v6=SIN06*u0-imul(SIN26,u0); } { value_type u0; u0=t30-t62; v7=SIN06*u0+imul(SIN26,u0); } { value_type u0,u1; u0=v0+v4; u1=v6+v7; Y.first(13,u0+u1); Y.second(45,u0-u1); } { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); Y.first(29,u0+u1); Y.second(61,u0-u1); } { value_type u0,u1; u0=v1+v5; u1=v2+v3; Y.first( 5,u0+u1); Y.second(37,u0-u1); } { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); Y.first(21,u0+u1); Y.second(53,u0-u1); } }
}

template<int M> struct complex_fft_base_function {};
template<> struct complex_fft_base_function<1>  { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_1 (X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<2>  
{
  template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_2 (X,f); } 
  template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } 
};
template<> struct complex_fft_base_function<3>  { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_3 (X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<4>  { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_4 (X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<5>  { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_5 (X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<8>  { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_8 (X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<16> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_16(X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<32> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_32(X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };
template<> struct complex_fft_base_function<64> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_64(X,f); } template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W,const F &f) { (*this)(X*W,f); } };

template<int M> struct complex_fft_function
{
  template<class G1,         class G3        > inline void operator()(const Vector<G1> &X,                      Vector<G3> &Y                            ) { (*this)(X  ,Y    ,make_fft_adaptor()); }
  template<class G1,class G2,class G3        > inline void operator()(const Vector<G1> &X, const Vector<G2> &W, Vector<G3> &Y                            ) { (*this)(X,W,Y    ,make_fft_adaptor()); }
  template<class G1,class G2,class G3        > inline void operator()(const Vector<G1> &X, const Vector<G2> &W, Vector<G3> &Y0, Vector<G3> &Y1           ) { (*this)(X,W,Y0,Y1,make_fft_adaptor()); }

  template<class G1,         class G3,class F> inline void operator()(const Vector<G1> &X,                      Vector<G3> &Y                 ,const F &f) { complex_fft_base_function<M>()(X  ,typename F::template rebind<M>::other(f).make0(Y    )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W, Vector<G3> &Y                 ,const F &f) { complex_fft_base_function<M>()(X,W,typename F::template rebind<M>::other(f).make (Y    )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W, Vector<G3> &Y0, Vector<G3> &Y1,const F &f) { complex_fft_base_function<M>()(X,W,typename F::template rebind<M>::other(f).make (Y0,Y1)); }
};



#ifdef SSE

template<class V> inline void cross(const Vector<simd_vector_generator<2,complex<V> > > &Xa, const Vector<simd_vector_generator<2,complex<V> > > &Xb, Vector<simd_vector_generator<2,complex<V> > > &Ya, Vector<simd_vector_generator<2,complex<V> > > &Yb) { Ya=movelh(Xa,Xb); Yb=movehl(Xb,Xa); }

template<class V>
inline void complex_crossed_fft_0_2_eval(const Vector<simd_vector_generator<2,complex<V> > > &X0, const Vector<simd_vector_generator<2,complex<V> > > &X1, Vector<simd_vector_generator<2,complex<V> > > &Y0, Vector<simd_vector_generator<2,complex<V> > > &Y1)
{
  Vector<simd_vector_generator<2,complex<V> > > a,b;
  a = movelh(X0,X1);
  b = movehl(X1,X0);
  Y0=a+b; Y1=a-b;
}

template<class V>
inline void complex_crossed_fft_0_2_eval(const Vector<simd_vector_generator<2,complex<V> > > &Xa0, const Vector<simd_vector_generator<2,complex<V> > > &Xa1, Vector<simd_vector_generator<2,complex<V> > > &Ya0, Vector<simd_vector_generator<2,complex<V> > > &Ya1,
                                         const Vector<simd_vector_generator<2,complex<V> > > &Xb0, const Vector<simd_vector_generator<2,complex<V> > > &Xb1, Vector<simd_vector_generator<2,complex<V> > > &Yb0, Vector<simd_vector_generator<2,complex<V> > > &Yb1)
{
  Vector<simd_vector_generator<2,complex<V> > > a,b;
  a = movelh(Xa0,Xb0); b = movehl(Xb0,Xa0);
  Ya0=a+b; Yb0=a-b;
  a = movelh(Xa1,Xb1); b = imul(movehl(Xb1,Xa1));
  Ya1=a-b; Yb1=a+b;
}

template<class V,class T2>
inline void complex_crossed_fft_0_2_eval(const Vector<simd_vector_generator<2,complex<V> > > &Xa0, const Vector<simd_vector_generator<2,complex<V> > > &Xa1, const T2 &Wa, Vector<simd_vector_generator<2,complex<V> > > &Ya0, Vector<simd_vector_generator<2,complex<V> > > &Ya1,
                                         const Vector<simd_vector_generator<2,complex<V> > > &Xb0, const Vector<simd_vector_generator<2,complex<V> > > &Xb1, const T2 &Wb, Vector<simd_vector_generator<2,complex<V> > > &Yb0, Vector<simd_vector_generator<2,complex<V> > > &Yb1)
{
  Vector<simd_vector_generator<2,complex<V> > > a,b;
  a=movelh(Xa0,Xb0); b=Wa*movehl(Xb0,Xa0);
  Ya0=a+b; Yb0=a-b;
  a=movelh(Xa1,Xb1); b=Wb*movehl(Xb1,Xa1);
  Ya1=a+b; Yb1=a-b;
}

template<class V>
inline void complex_crossed_fft_0_2(const Vector<simd_vector_generator<2,complex<V> > > &X0, const Vector<simd_vector_generator<2,complex<V> > > &X1, Vector<simd_vector_generator<2,complex<V> > > &Y0, Vector<simd_vector_generator<2,complex<V> > > &Y1)
{
  Vector<simd_vector_generator<2,complex<V> > > a,b;
  a = movelh(X0,X1);
  b = mi1mul(movehl(X1,X0));
  Y0=a+b; Y1=a-b;
}

template<class V>
inline void complex_crossed_fft_0_2(const Vector<simd_vector_generator<2,complex<V> > > &X0, const Vector<simd_vector_generator<2,complex<V> > > &X1, const Vector<simd_vector_generator<2,complex<V> > > &W, Vector<simd_vector_generator<2,complex<V> > > &Y0, Vector<simd_vector_generator<2,complex<V> > > &Y1)
{
  typename TinyVector<2,complex<V> >::self a,b;
  a =   movelh(X0,X1);
  b = W*movehl(X1,X0);
  Y0=a+b; Y1=a-b;
}

template<class V,class T2>
inline void complex_crossed_fft_0_2(const Vector<simd_vector_generator<2,complex<V> > > &Xa0, const Vector<simd_vector_generator<2,complex<V> > > &Xa1              , Vector<simd_vector_generator<2,complex<V> > > &Ya0, Vector<simd_vector_generator<2,complex<V> > > &Ya1,
                                    const Vector<simd_vector_generator<2,complex<V> > > &Xb0, const Vector<simd_vector_generator<2,complex<V> > > &Xb1, const T2 &Wb, Vector<simd_vector_generator<2,complex<V> > > &Yb0, Vector<simd_vector_generator<2,complex<V> > > &Yb1)
{
  Vector<simd_vector_generator<2,complex<V> > > Ya2,Ya3;
  { Vector<simd_vector_generator<2,complex<V> > > a,b; a=movelh(Xa0,Xa1); b=mi1mul(movehl(Xa1,Xa0)); Ya2=a+b; Ya3=a-b; }
  { Vector<simd_vector_generator<2,complex<V> > > a,b; a=movelh(Xb0,Xb1); b=Wb    *movehl(Xb1,Xb0 ); Yb0=a+b; Yb1=a-b; }
  Ya0=movelh(Ya2,Yb0); Ya1=movehl(Yb0,Ya2);
  Yb0=movelh(Ya3,Yb1); Yb1=movehl(Yb1,Ya3);
}

template<class V,class T2>
inline void complex_crossed_fft_0_2(const Vector<simd_vector_generator<2,complex<V> > > &Xa0, const Vector<simd_vector_generator<2,complex<V> > > &Xa1, const T2 &Wa, Vector<simd_vector_generator<2,complex<V> > > &Ya0, Vector<simd_vector_generator<2,complex<V> > > &Ya1,
                                    const Vector<simd_vector_generator<2,complex<V> > > &Xb0, const Vector<simd_vector_generator<2,complex<V> > > &Xb1, const T2 &Wb, Vector<simd_vector_generator<2,complex<V> > > &Yb0, Vector<simd_vector_generator<2,complex<V> > > &Yb1)
{
  Vector<simd_vector_generator<2,complex<V> > > Ya2,Ya3;
  { typename TinyVector<2,complex<V> >::self a,b; a=movelh(Xa0,Xa1); b=Wa*movehl(Xa1,Xa0); Ya2=a+b; Ya3=a-b; }
  { typename TinyVector<2,complex<V> >::self a,b; a=movelh(Xb0,Xb1); b=Wb*movehl(Xb1,Xb0); Yb0=a+b; Yb1=a-b; }
  Ya0=movelh(Ya2,Yb0); Ya1=movehl(Yb0,Ya2);
  Yb0=movelh(Ya3,Yb1); Yb1=movehl(Yb1,Ya3);
}
#endif

template<int M,int N> struct complex_crossed_fft_base_function {};
template<> struct complex_crossed_fft_base_function<0,2>
{
  template<class G1,         class G3> inline void eval      (const Vector<G1> &X0 , const Vector<G1> &X1               , Vector<G3> &Y0 , Vector<G3> &Y1 ) { complex_crossed_fft_0_2_eval(X0,X1,  Y0,Y1); }
  template<class G1,         class G3> inline void eval      (const Vector<G1> &Xa0, const Vector<G1> &Xa1              , Vector<G3> &Ya0, Vector<G3> &Ya1, const Vector<G1> &Xb0, const Vector<G1> &Xb1,               Vector<G3> &Yb0, Vector<G3> &Yb1) { complex_crossed_fft_0_2_eval(Xa0,Xa1   ,Ya0,Ya1, Xb0,Xb1   ,Yb0,Yb1); }
  template<class G1,class T2,class G3> inline void eval      (const Vector<G1> &Xa0, const Vector<G1> &Xa1, const T2 &Wa, Vector<G3> &Ya0, Vector<G3> &Ya1, const Vector<G1> &Xb0, const Vector<G1> &Xb1, const T2 &Wb, Vector<G3> &Yb0, Vector<G3> &Yb1) { complex_crossed_fft_0_2_eval(Xa0,Xa1,Wa,Ya0,Ya1, Xb0,Xb1,Wb,Yb0,Yb1); }

  template<class G1,         class G3> inline void operator()(const Vector<G1> &X0 , const Vector<G1> &X1 ,               Vector<G3> &Y0 , Vector<G3> &Y1 ) { complex_crossed_fft_0_2(X0,X1,  Y0,Y1); }
  template<class G1,class T2,class G3> inline void operator()(const Vector<G1> &X0 , const Vector<G1> &X1 , const T2 &W , Vector<G3> &Y0 , Vector<G3> &Y1 ) { complex_crossed_fft_0_2(X0,X1,W,Y0,Y1); }

  template<class G1,class T2,class G3> inline void operator()(const Vector<G1> &Xa0, const Vector<G1> &Xa1              , Vector<G3> &Ya0, Vector<G3> &Ya1, const Vector<G1> &Xb0, const Vector<G1> &Xb1, const T2 &Wb, Vector<G3> &Yb0, Vector<G3> &Yb1) { complex_crossed_fft_0_2(Xa0,Xa1   ,Ya0,Ya1, Xb0,Xb1,Wb,Yb0,Yb1); }
  template<class G1,class T2,class G3> inline void operator()(const Vector<G1> &Xa0, const Vector<G1> &Xa1, const T2 &Wa, Vector<G3> &Ya0, Vector<G3> &Ya1, const Vector<G1> &Xb0, const Vector<G1> &Xb1, const T2 &Wb, Vector<G3> &Yb0, Vector<G3> &Yb1) { complex_crossed_fft_0_2(Xa0,Xa1,Wa,Ya0,Ya1, Xb0,Xb1,Wb,Yb0,Yb1); }
};

template<class G1,class S>
inline void complex_fft_1_2(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;

  complex_fft_2(X[0],Y);
}

template<class G1,class S>
inline void complex_fft_1_2(const Vector<G1> &Xa, const Vector<G1> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;

  value_type v0,v1;
  cross(Xa[0],Xb[0],v0,v1);
  Ya.first0(  v0+v1); Ya.second(1,v0-v1);
}


template<class G1,class S>
inline void complex_fft_2_2(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;

  value_type z0(X[0]), z1(X[1]);
  value_type t1,t2, a0,a2;
  t1=z0+z1; t2=z0-z1;
  complex_crossed_fft_base_function<0,2>()(t1,t2,a0,a2); Y.first0(  a0); Y.second(2,a2);
}

template<class G1,class S>
inline void complex_fft_2_2(const Vector<G1> &Xa, const Vector<G1> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3;
  {
    { value_type v0,v1,v2,v3; cross(Xa[0],Xb[0],v0,v1); cross(Xa[1],Xb[1],v2,v3); t0=v0+v2; t1=v0-v2; t2=v1+v3; t3=v1-v3; t3=imul(t3); }
  }
  Ya.first0(  t0+t2); Ya.second(2,t0-t2);Ya.first (1,t1-t3); Ya.second(3,t1+t3);
}

template<class G1,class S>
inline void complex_fft_3_2(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3;
  {
    value_type v0, v1;
    t0=X[0]; v0=X[1]; v1=X[2]; t1=v0+v1; t2=t0-real_type(0.5)*t1; t3=imul(SIN01_03,v1-v0);
  }

  value_type a0,a1;
  complex_fft_function<2>()(t0+t1,a0);
  Y.first0(  a0[0]);
  Y.second(3,a0[1]);

  complex_crossed_fft_base_function<0,2>()(t2+t3,t2-t3,value_type(complex_type(0.5,-SIN01_03),complex_type(-0.5,-SIN01_03)),a0,a1);
  Y.first (1,a0[0]); Y.first (2,a0[1]);
  Y.second(4,a1[0]); Y.second(5,a1[1]);
}

template<class G,class S>
inline void complex_fft_3_2(const Vector<G> &Xa, const Vector<G> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type ta0,ta1,ta2,ta3, a0,a1;
  value_type tb0,tb1,tb2,tb3, b0,b1;

  {
    value_type v0, v1;
    ta0=Xa[0]; v0=Xa[1]; v1=Xa[2]; ta1=v0+v1; ta2=ta0-real_type(0.5)*ta1; ta3=imul(SIN01_03,v1-v0);
    tb0=Xb[0]; v0=Xb[1]; v1=Xb[2]; tb1=v0+v1; tb2=tb0-real_type(0.5)*tb1; tb3=imul(SIN01_03,v1-v0);
  }

  a0=ta0+ta1;
  b0=tb0+tb1;
  complex_crossed_fft_base_function<0,2>().eval(a0,b0,a0,b0);
  Ya.first0(  a0);
  Ya.second(3,b0);

  a0=ta2+ta3; a1=ta2-ta3;
  b0=tb2+tb3; b1=tb2-tb3;
  complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type( 0.5,-SIN01_03),a0,a1, b0,b1,complex_type(-0.5,-SIN01_03),b0,b1);
  Ya.first(1,a0); Ya.first(2,a1); Ya.second(4,b0); Ya.second(5,b1);
}

template<class G1,class S>
inline void complex_fft_4_2(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3;
  { value_type a,b; a=X[0]; b=X[2]; t0=a+b; t1=a-b; a=X[1]; b=X[3]; t2=a+b; t3=a-b; t3=imul(t3); }
  { value_type a0,a1,b0,b1; complex_crossed_fft_base_function<0,2>()(t0+t2,t0-t2,a0,a1, t1-t3,t1+t3,value_type(complex_type(SIN16,-SIN16),complex_type(-SIN16,-SIN16)),b0,b1); Y.first0(  a0); Y.first (2,a1); Y.second(4,b0); Y.second(6,b1); }
}

template<class G1,class S>
inline void complex_fft_4_2(const Vector<G1> &Xa, const Vector<G1> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3,t4,t5,t6,t7;
  {
    { value_type v0,v1,v2,v3; cross(Xa[0],Xb[0],v0,v2); cross(Xa[2],Xb[2],v1,v3); t0=v0+v1; t1=     v0-v1;  t6=v2+v3; v2=v2-v3; t7=sub_imul(v2,v2); }
    { value_type v0,v1,v2,v3; cross(Xa[1],Xb[1],v0,v3); cross(Xa[3],Xb[3],v1,v2); t2=v0+v1; t3=imul(v0-v1); t4=v2+v3; v2=v2-v3; t5=add_imul(v2,v2); }
  }

  { value_type v0,v1; v0=t0+t2; v1=            t4+t6 ; Ya.first0(  v0+v1); Ya.second(4,v0-v1); }
  { value_type v0,v1; v0=t0-t2; v1=imul(       t4-t6); Ya.first (2,v0+v1); Ya.second(6,v0-v1); }
  { value_type v0,v1; v0=t1-t3; v1=     SIN16*(t5+t7); Ya.first (1,v0+v1); Ya.second(5,v0-v1); }
  { value_type v0,v1; v0=t1+t3; v1=imul(SIN16, t5-t7); Ya.first (3,v0+v1); Ya.second(7,v0-v1); }
}

template<class G1,class S>
inline void complex_fft_5_2(const Vector<G1> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3,t4;

  {
    value_type v0,v1;
    t0=X[0];
    { value_type v2,v3; v2=X[1]; v3=X[4]; v0=v2+v3; t2=imul(v2-v3); v2=X[2]; v3=X[3]; v1=v2+v3; t4=imul(v2-v3); }
    t3=real_type(SQRT_5_16)*(v0-v1); v0=v0+v1;
    t1=t0-real_type(0.25)*v0;
    t0=t0+v0;
  }

  {
    value_type a0;
    complex_fft_function<2>()(t0,a0);
    Y.first0(  a0[0]);
    Y.second(5,a0[1]);
  }
  {
    value_type v0,v1, a0,a1;
    v0=SIN02_05*t2+SIN01_05*t4; v1=t1+t3;
    complex_crossed_fft_base_function<0,2>()(v1-v0,v1+v0,value_type(complex_type(SIN03_10,-SIN01_05),complex_type(-SIN03_10,-SIN01_05)),a0,a1);
    Y.first (1,a0[0]); Y.first (4,a0[1]);
    Y.second(6,a1[0]); Y.second(9,a1[1]);
  }
  {
    value_type v0,v1, a0,a1;
    v0=SIN01_05*t2-SIN02_05*t4; v1=t1-t3;
    complex_crossed_fft_base_function<0,2>()(v1-v0,v1+v0,value_type(complex_type(SIN01_10,-SIN02_05),complex_type(-SIN01_10,-SIN02_05)),a0,a1);
    Y.first (2,a0[0]); Y.first (3,a0[1]);
    Y.second(7,a1[0]); Y.second(8,a1[1]);
  }
}

template<class G,class S>
inline void complex_fft_5_2(const Vector<G> &Xa, const Vector<G> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type ta0,ta1,ta2,ta3,ta4;
  value_type tb0,tb1,tb2,tb3,tb4;

  {
    value_type v0,v1;
    ta0=Xa[0]; { value_type v2,v3; v2=Xa[1]; v3=Xa[4]; v0=v2+v3; ta2=imul(v2-v3); v2=Xa[2]; v3=Xa[3]; v1=v2+v3; ta4=imul(v2-v3); } ta1=v0+v1; ta3=real_type(SQRT_5_16)*(v0-v1);
    tb0=Xb[0]; { value_type v2,v3; v2=Xb[1]; v3=Xb[4]; v0=v2+v3; tb2=imul(v2-v3); v2=Xb[2]; v3=Xb[3]; v1=v2+v3; tb4=imul(v2-v3); } tb1=v0+v1; tb3=real_type(SQRT_5_16)*(v0-v1);
  }

  {
    value_type a,b;
    a=ta0+ta1;
    b=tb0+tb1;
    complex_crossed_fft_base_function<0,2>().eval(a,b,a,b);
    Ya.first0(  a);
    Ya.second(5,b);
  }
  ta0-=real_type(0.25)*ta1;
  tb0-=real_type(0.25)*tb1;

  {
    value_type v0,v1, a0,a1,b0,b1;
    v0=SIN02_05*ta2+SIN01_05*ta4; v1=ta0+ta3; a0=v1-v0; a1=v1+v0;
    v0=SIN02_05*tb2+SIN01_05*tb4; v1=tb0+tb3; b0=v1-v0; b1=v1+v0;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN03_10,-SIN01_05),a0,a1, b0,b1,complex_type(-SIN03_10,-SIN01_05),b0,b1);
    Ya.first(1,a0); Ya.first(4,a1); Ya.second(6,b0); Ya.second(9,b1);
  }
  {
    value_type v0,v1, a0,a1,b0,b1;
    v0=SIN01_05*ta2-SIN02_05*ta4; v1=ta0-ta3; a0=v1-v0; a1=v1+v0;
    v0=SIN01_05*tb2-SIN02_05*tb4; v1=tb0-tb3; b0=v1-v0; b1=v1+v0;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN01_10,-SIN02_05),a0,a1, b0,b1,complex_type(-SIN01_10,-SIN02_05),b0,b1);
    Ya.first(2,a0); Ya.first(3,a1); Ya.second(7,b0); Ya.second(8,b1);
  }
}


template<class G,class S>
inline void complex_fft_8_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t1,t2,t3,t4,t5,t6,t7,t8;
  {
    value_type v0,v1;
    v0=X[0],v1=X[4]; t1=v0+v1; t2=v0-v1;
    v0=X[2],v1=X[6]; t3=v0+v1; t4=v0-v1; t4=imul(t4);
    v0=X[7],v1=X[3]; t5=v0+v1; v0=v0-v1; t6=add_imul(v0,v0);
    v0=X[1],v1=X[5]; t7=v0+v1; v0=v0-v1; t8=sub_imul(v0,v0);
  }

  {
    value_type a0,a1,b0,b1;
    { value_type v0,v1; v0=t1+t3; v1=t7+t5;                        a0=v0+v1; a1=v0-v1; }
    { value_type v0,v1; v0=t2-t4; v1=SIN16*(t8+t6);  b0=v0+v1; b1=v0-v1; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,a0,a1, b0,b1,value_type(complex_type(SIN24,-SIN08),complex_type(-SIN08,-SIN24)),b0,b1);
    Y.first0(   a0); Y.second( 8,b0);
    Y.first ( 4,a1); Y.second(12,b1);
  }
  {
    value_type a0,a1,b0,b1;
    { value_type v0,v1; v0=imul(t5-t7); v1=t1-t3;                     a0=v1+v0; a1=v1-v0; }
    { value_type v0,v1; v0=t2+t4; v1=imul(SIN16,t6-t8); b0=v0+v1; b1=v0-v1; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN16,-SIN16),complex_type(-SIN16,-SIN16)),a0,a1,
  	                                         b0,b1,value_type(complex_type(SIN08,-SIN24),complex_type(-SIN24,-SIN08)),b0,b1);
    Y.first ( 2,a0); Y.second(10,b0);
    Y.first ( 6,a1); Y.second(14,b1);
  }
}


template<class G1,class S>
inline void complex_fft_8_2(const Vector<G1> &Xa, const Vector<G1> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15;

  {
    { value_type v0,v1,v2,v3; cross(Xa[0],Xb[0],v0,v2); cross(Xa[4],Xb[4],v1,v3); t01=v0+v1; t03=v0-v1; t13=v2+v3; t15=v2-v3; }
    { value_type v0,v1,v2,v3; cross(Xa[2],Xb[2],v0,v2); cross(Xa[6],Xb[6],v1,v3); t02=v0+v1; v0=imul(v0-v1); t00=t01+t02; t01-=t02; t02=t03+v0; t03-=v0; t14=v2+v3; v2=imul(v2-v3);  t12=t13+t14; t13-=t14; t14=t15+v2; t15-=v2; }
    { value_type v0,v1,v2,v3; cross(Xa[3],Xb[3],v0,v3); cross(Xa[7],Xb[7],v1,v2); t09=v0+v1; t11=v1-v0; t05=v2+v3; t07=v2-v3; }
    { value_type v0,v1,v2,v3; cross(Xa[5],Xb[5],v0,v3); cross(Xa[1],Xb[1],v1,v2); t10=v0+v1; v0-=v1; t08=t09+t10; t09=imul(t09-t10); t10=t11-v0; t11=imul(t11+v0); t06=v2+v3; v2=imul(v2-v3); t04=t05+t06; t05-=t06; t06=t07+v2; t07-=v2; }
  }

  { value_type v0,v1; v0=t00 +t08; v1=t12+t04; Ya.first0(  v0+v1); Ya.second(8,v0-v1); }
  { value_type v0,v1; v0=imul(t04-t12); v1=t00-t08; Ya.first ( 4,v1+v0); Ya.second(12,v1-v0); }
  { value_type v0,v1; v0=sub_imul(t13,t13); v1=add_imul(t05,t05); v1=SIN16*(v0+v1); v0=t01+t09; Ya.first ( 2,v0+v1); Ya.second(10,v0-v1); }
  { value_type v0,v1; v0=add_imul(t13,t13); v1=imul(t05)-t05;     v1=SIN16*(v1-v0); v0=t01-t09; Ya.first ( 6,v0+v1); Ya.second(14,v0-v1); }
  { value_type v0,v1,v2,v3; v0=SIN08*t14-imul(SIN24,t14); v1=SIN08*t06 +imul(SIN24,t06 ); v2=v1+v0; v3=imul(v1-v0); { value_type v4; v4=SIN16*(t11-t10); v0=t02+v4; v1=t02-v4; } Ya.first ( 3,v0+v2); Ya.second(11,v0-v2); Ya.first ( 7,v1+v3); Ya.second(15,v1-v3); }
  { value_type v0,v1,v2,v3; v0=SIN24*t15-imul(SIN08,t15); v1=SIN24*t07 +imul(SIN08,t07 ); v2=v1+v0; v3=imul(v1-v0); { value_type v4; v4=SIN16*(t11+t10); v0=t03+v4; v1=t03-v4; } Ya.first ( 1,v0+v2); Ya.second( 9,v0-v2); Ya.first ( 5,v1+v3); Ya.second(13,v1-v3); }
}



template<class G,class S>
inline void complex_fft_16_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;

  {
    value_type v0,v1;
    v0=X[ 0]; v1=X[ 8]; t1 =v0+v1; t3 =v0-v1; v0=X[ 4]; v1=X[12]; t2 =v0+v1; v0-=v1; v0=imul(v0); t0=t1+t2; t1-=t2; t2=t3+v0; t3-=v0;
    v0=X[15]; v1=X[ 7]; t5 =v0+v1; t7 =v0-v1; v0=X[ 3]; v1=X[11]; t6 =v0+v1; v0-=v1; v0=imul(v0); t4=t5+t6; t5-=t6; t6=t7+v0; t7-=v0;
    v0=X[ 6]; v1=X[14]; t9 =v0+v1; t11=v1-v0; v0=X[10]; v1=X[ 2]; t10=v0+v1; v0-=v1; t8=t9+t10; t9-=t10; t9=imul(t9); t10=t11-v0; t11+=v0; t11=imul(t11);
    v0=X[ 1]; v1=X[ 9]; t13=v0+v1; t15=v0-v1; v0=X[ 5]; v1=X[13]; t14=v0+v1; v0-=v1; v0=imul(v0);  t12=t13+t14; t13-=t14; t14=t15+v0; t15-=v0;
  }

  {
    value_type a0,a1,b0,b1, b2,b3;
    { value_type v0,v1; v0=t0+t8; v1=t12+t4; a0=v0+v1; a1=v0-v1; }
    { value_type v0,v1,v2,v3,v4; v1=SIN24*t15-imul(SIN08,t15); v2=SIN24*t7 +imul(SIN08,t7 ); v3=v2+v1; v4=v2-v1; v4=imul(v4); v0=SIN16*(t11+t10); v1=t3+v0; v2=t3-v0; b0=v1+v3; b1=v1-v3; b2=v2+v4; b3=v2-v4; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,a0,a1, b0,b1,value_type(complex_type(SIN28,-SIN04),complex_type(-SIN04,-SIN28)),b0,b1); Y.first0(   a0); Y.second(16,b0); Y.first ( 8,a1); Y.second(24,b1);
    { value_type v0,v1; v0=imul(t4-t12); v1=t0-t8; a0=v1+v0; a1=v1-v0; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN16,-SIN16),complex_type(-SIN16,-SIN16)),a0,a1, b2,b3,value_type(complex_type(SIN12,-SIN20),complex_type(-SIN20,-SIN12)),b2,b3); Y.first ( 4,a0); Y.second(20,b2); Y.first (12,a1); Y.second(28,b3);
  }
  {
    value_type a0,a1,b0,b1, b2,b3;
    { value_type v0,v1; v0=sub_imul(t13,t13); v1=add_imul(t5,t5); v1+=v0; v1*=SIN16; v0=t1+t9; a0=v0+v1; a1=v0-v1; }
    { value_type v0,v1,v2,v3,v4; v1=SIN08*t14-imul(SIN24,t14); v2=SIN08*t6 +imul(SIN24,t6 ); v3=v2+v1; v4=v2-v1; v4=imul(v4); v0=SIN16*(t11-t10);  v1=t2+v0; v2=t2-v0; b0=v1+v3; b1=v1-v3; b2=v2+v4; b3=v2-v4; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN24,-SIN08),complex_type(-SIN08,-SIN24)),a0,a1, b0,b1,value_type(complex_type(SIN20,-SIN12),complex_type(-SIN12,-SIN20)),b0,b1); Y.first ( 2,a0); Y.second(18,b0); Y.first (10,a1); Y.second(26,b1);
    { value_type v0,v1; v0=add_imul(t13,t13); v1=imul(t5)-t5; v1-=v0; v1*=SIN16; v0=t1-t9; a0=v0+v1; a1=v0-v1; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN08,-SIN24),complex_type(-SIN24,-SIN08)),a0,a1, b2,b3,value_type(complex_type(SIN04,-SIN28),complex_type(-SIN28,-SIN04)),b2,b3); Y.first ( 6,a0); Y.second(22,b2); Y.first (14,a1); Y.second(30,b3);
  }
}

template<class G1,class S>
inline void complex_fft_16_2(const Vector<G1> &Xa, const Vector<G1> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33;

  {
    {
      value_type v0,v1,v2,v3, w0,w1,w2,w3, u6;
      { value_type r0,r1,s0,s1; cross(Xa[ 0],Xb[ 0],r0,s0); cross(Xa[ 8],Xb[ 8],r1,s1); v0=r0+r1; v1=r0-r1; w0=s0+s1; w1=s0-s1; cross(Xa[ 4],Xb[ 4],r0,s0); cross(Xa[12],Xb[12],r1,s1); v2=r0+r1; v3=imul(r0-r1); w2=s0+s1; w3=imul(s0-s1); } t00=v0+v2; t01=v0-v2; t02=v1+v3; t03=v1-v3; t25=w1+w3; t26=w1-w3; t27=w0+w2; u6=w0-w2;
      { value_type r0,r1,s0,s1; cross(Xa[ 2],Xb[ 2],r0,s0); cross(Xa[10],Xb[10],r1,s1); v0=r0+r1; v1=r0-r1; w0=s0+s1; w1=s0-s1; cross(Xa[ 6],Xb[ 6],r0,s0); cross(Xa[14],Xb[14],r1,s1); v2=r0+r1; v3=     r1-r0 ; w2=s0+s1; w3=     s1-s0 ; } t28=w0+w2; t29=t27-t28; w2=imul(w2-w0); t30=u6+w2; t31=u6-w2; u6=w3+w1; w2=imul(w3-w1); t32=SIN16*(w2+u6); t33=SIN16*(w2-u6); t04=v0+v2; t05=imul(v0-v2); v0=v3+v1; v2=imul(v3-v1); t06=SIN16*(v2+v0); t07=SIN16*(v2-v0);
    }
    {
      value_type v0,v1,v2,v3, w0,w1,w2,w3, u6;
      { value_type r0,r1,s0,s1; cross(Xa[15],Xb[15],r0,s0); cross(Xa[ 7],Xb[ 7],r1,s1); v0=r0+r1; v1=r0-r1; w0=s0+s1; w1=s0-s1; cross(Xa[ 3],Xb[ 3],r0,s0); cross(Xa[11],Xb[11],r1,s1); v2=r0+r1; v3=imul(r0-r1); w2=s0+s1; w3=imul(s0-s1); } t12=v0+v2; v0=v0-v2; t15=add_imul(v0,v0); v0=v1+v3; v2=v1-v3; t13=SIN08*v0+imul(SIN24,v0); t14=SIN24*v2+imul(SIN08,v2); t16=w1+w3; t17=w1-w3; t18=w0+w2; u6=w0-w2;
      { value_type r0,r1,s0,s1; cross(Xa[ 1],Xb[ 1],r0,s0); cross(Xa[ 9],Xb[ 9],r1,s1); v0=r0+r1; v1=r0-r1; w0=s0+s1; w1=s0-s1; cross(Xa[ 5],Xb[ 5],r0,s0); cross(Xa[13],Xb[13],r1,s1); v2=r0+r1; v3=imul(r0-r1); w2=s0+s1; w3=     s1-s0 ; } t19=w0+w2; w0=imul(w0-w2); t20=t18-t19; t21=u6+w0; t22=u6-w0; u6=w3+w1; w2=imul(w3-w1); t23=SIN16*(w2+u6); t24=SIN16*(w2-u6); t08=v0+v2; v0=v0-v2; t11=sub_imul(v0,v0)    ; v0=v1+v3; v2=v1-v3; t09 =SIN08*v0-imul(SIN24,v0); t10=SIN24*v2-imul(SIN08,v2);
    }
  }

  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t00+t04; v5=t12+t08; v0=v4+v5; v1=v4-v5; v4=t27+t28; v5=t18+t19; v2=v5+v4 ; v3=imul(v5-v4); } Ya.first0(   v0+v2); Ya.second(16,v0-v2); Ya.first ( 8,v1+v3); Ya.second(24,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t00-t04; v5=imul(t12-t08);       v0=v4+v5; v1=v4-v5; v4=t20+t29; v5=imul(t20-t29); v2=SIN16*(v4+v5); v3=SIN16*(v5-v4); } Ya.first ( 4,v0+v2); Ya.second(20,v0-v2); Ya.first (12,v1+v3); Ya.second(28,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t01-t05; v5=SIN16*(t11+t15);     v0=v4+v5; v1=v4-v5; v4=SIN24*t30-imul(SIN08,t30); v5=SIN24*t22+imul(SIN08,t22); v2=v4+v5 ; v3=imul(v5-v4); } Ya.first ( 2,v0+v2); Ya.second(18,v0-v2); Ya.first (10,v1+v3); Ya.second(26,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t01+t05; v5=imul(SIN16,t11-t15); v1=v4+v5; v0=v4-v5;v4=SIN08*t31-imul(SIN24,t31);v5=SIN08*t21+imul(SIN24,t21); v2=v4+v5 ; v3=imul(v5-v4); } Ya.first ( 6,v0+v2); Ya.second(22,v0-v2); Ya.first (14,v1+v3); Ya.second(30,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t03+t06; v5=t14+t10; v0=v4+v5;   v1=v4-v5; v4=t26+t32; v4=SIN28*v4-imul(SIN04,v4) ; v5=t17+t23; v5=SIN28*v5+imul(SIN04,v5); v2=v5+v4 ; v3=imul(v5-v4); } Ya.first ( 1,v0+v2); Ya.second(17,v0-v2); Ya.first ( 9,v1+v3); Ya.second(25,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t03-t06; v5=imul(t14-t10);       v0=v4+v5; v1=v4-v5; v4=t26-t32; v4=SIN12*v4-imul(SIN20,v4) ; v5=t17-t23; v5=SIN12*v5+imul(SIN20,v5); v2=v4+v5 ; v3=imul(v5-v4); } Ya.first ( 5,v0+v2); Ya.second(21,v0-v2); Ya.first (13,v1+v3); Ya.second(29,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t02+t07; v5=t13+t09;             v0=v4+v5; v1=v4-v5; v4=t25+t33; v4=SIN20*v4-imul(SIN12,v4) ; v5=t16+t24; v5=SIN20*v5+imul(SIN12,v5); v2=v5+v4 ; v3=imul(v5-v4); } Ya.first ( 3,v0+v2); Ya.second(19,v0-v2); Ya.first (11,v1+v3); Ya.second(27,v1-v3); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t02-t07; v5=imul(t13-t09);       v0=v4+v5; v1=v4-v5; v4=t25-t33; v4=SIN04*v4-imul(SIN28,v4) ; v5=t16-t24; v5=SIN04*v5+imul(SIN28,v5); v2=v4+v5 ; v3=imul(v5-v4); } Ya.first ( 7,v0+v2); Ya.second(23,v0-v2); Ya.first (15,v1+v3); Ya.second(31,v1-v3); }
}



template<class G,class S>
inline void complex_fft_32_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33;

  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[0]; v5=X[16]; v0=v4+v5; v1=v4-v5; v4=X[8]; v5=X[24]; v2=v4+v5; v3=imul(v4-v5); } t0=v0+v2; t1=v0-v2; t2=v1+v3; t3=v1-v3; }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[ 4]; v5=X[20]; v0=v4+v5; v1=v4-v5; v4=X[28]; v5=X[12]; v2=v4+v5;  v3=v4-v5; } t4=v0+v2; t5=imul(v0-v2); { value_type v4,v5; v4=v3+v1; v5=imul(v3-v1); t6=SIN16*(v5+v4); t7=SIN16*(v5-v4); } }
  { value_type v0,v1,v2,v3,v4; { value_type v5,v6; v5=X[ 2]; v6=X[18]; v0=v5+v6; v1=v5-v6; v5=X[10]; v6=X[26]; v2=v5+v6; v3=imul(v5-v6); } v4=v1+v3; t8=v0+v2; t9=SIN08*v4-imul(SIN24,v4); v4=v1-v3; t10=SIN24*v4-imul(SIN08,v4); v4=v0-v2; t11=v4-imul(v4); }
  { value_type v0,v1,v2,v3,v4; { value_type v5,v6; v5=X[30]; v6=X[14]; v0=v5+v6; v1=v5-v6; v5=X[ 6]; v6=X[22]; v2=v5+v6; v3=imul(v5-v6); } t12=v0+v2; v4=v1+v3; t13=SIN08*v4+imul(SIN24,v4); v4=v1-v3; t14=SIN24*v4+imul(SIN08,v4); v4=v0-v2; t15=add_imul(v4,v4); }
  { value_type v0,v1,v2,v3,v4,v5; { value_type v6,v7; v6=X[31]; v7=X[15]; v0=v6+v7; v1=v6-v7; v6=X[7]; v7=X[23]; v2=v6+v7; v3=imul(v6-v7); } t16=v1+v3; t17=v1-v3; t18=v0+v2; v0-=v2; { value_type v6,v7; v6=X[ 3]; v7=X[19]; v4=v6+v7; v5=v6-v7; v6=X[27]; v7=X[11]; v1=v6+v7; v3=v6-v7; } t19=v4+v1; v4-=v1; t20=t18-t19; { v4=imul(v4); t21=v0+v4; t22=v0-v4; v0=v3+v5; v2=imul(v3-v5); t23=SIN16*(v2+v0); t24=SIN16*(v2-v0); } }
  { value_type v0,v1,v2,v3,v4,v5; { value_type v6,v7; v6=X[ 1]; v7=X[17]; v0=v6+v7; v1=v6-v7; v6=X[9]; v7=X[25]; v2=v6+v7; v3=imul(v6-v7); } t25=v1+v3; t26=v1-v3; t27=v0+v2; v0-=v2; { value_type v6,v7; v6=X[ 5]; v7=X[21]; v4=v6+v7; v5=v6-v7; v6=X[29]; v7=X[13]; v1=v6+v7; v3=v6-v7; } t28=v4+v1; t29=t27-t28; { v1-=v4; v1=imul(v1); t30=v0+v1; t31=v0-v1; v0=v3+v5; v2=imul(v3-v5); t32=SIN16*(v2+v0); t33=SIN16*(v2-v0); } }

  {
    value_type a0,a1,a2,a3, b0,b1,b2,b3;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t0+t4; v5=t8 +t12; v0=v4+v5; v1=v4-v5; v4=t27+t28; v5=t18+t19; v2=v5+v4; v3=imul(v5-v4); } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t3+t6; v5=t10+t14; v0=v4+v5; v1=v4-v5; v4=t26+t32; v4=SIN28*v4-imul(SIN04,v4); v5=t17+t23; v5=SIN28*v5+imul(SIN04,v5); v2=v5+v4; v3=imul(v5-v4); } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,                                                                                           a0,a1, b0,b1,value_type(complex_type(SIN30,-SIN02),complex_type(-SIN02,-SIN30)),b0,b1); Y.first0(   a0); Y.second(32,b0); Y.first (16,a1); Y.second(48,b1);
    complex_crossed_fft_base_function<0,2>()(a2,a3,value_type(complex_type(SIN16,-SIN16),complex_type(-SIN16,-SIN16)),a2,a3, b2,b3,value_type(complex_type(SIN14,-SIN18),complex_type(-SIN18,-SIN14)),b2,b3); Y.first ( 8,a2); Y.second(40,b2); Y.first (24,a3); Y.second(56,b3);
  }
  {
    value_type a0,a1,a2,a3, b0,b1,b2,b3;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t0-t4; v5=imul(t12-t8); v0=v4+v5; v1=v4-v5; v4=t20+t29; v5=imul(t20-t29); v2=SIN16*(v4+v5); v3=SIN16*(v5-v4); } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t3-t6; v5=imul(t14-t10); v0=v4+v5; v1=v4-v5; v4=t26-t32; v4=SIN12*v4-imul(SIN20,v4); v5=t17-t23; v5=SIN12*v5+imul(SIN20,v5); v2=v4+v5; v3=imul(v5-v4); } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN24,-SIN08),complex_type(-SIN08,-SIN24)),a0,a1, b0,b1,value_type(complex_type(SIN22,-SIN10),complex_type(-SIN10,-SIN22)),b0,b1); Y.first ( 4,a0); Y.second(36,b0); Y.first (20,a1); Y.second(52,b1);
    complex_crossed_fft_base_function<0,2>()(a2,a3,value_type(complex_type(SIN08,-SIN24),complex_type(-SIN24,-SIN08)),a2,a3, b2,b3,value_type(complex_type(SIN06,-SIN26),complex_type(-SIN26,-SIN06)),b2,b3); Y.first (12,a2); Y.second(44,b2); Y.first (28,a3); Y.second(60,b3); }
  {
    value_type a0,a1,a2,a3, b0,b1,b2,b3;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t1-t5; v5=SIN16*(t11+t15); v0=v4+v5; v1=v4-v5; v4=SIN24*t30-imul(SIN08,t30); v5=SIN24*t22+imul(SIN08,t22); v2=v4+v5; v3=imul(v5-v4); } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t2+t7; v5=t9+t13; v0=v4+v5; v1=v4-v5; v4=t25+t33; v4=SIN20*v4-imul(SIN12,v4); v5=t16+t24; v5=SIN20*v5+imul(SIN12,v5); v2= v5+v4; v3=imul(v5-v4); } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN28,-SIN04),complex_type(-SIN04,-SIN28)),a0,a1, b0,b1,value_type(complex_type(SIN26,-SIN06),complex_type(-SIN06,-SIN26)),b0,b1); Y.first ( 2,a0); Y.second(34,b0); Y.first (18,a1); Y.second(50,b1);
    complex_crossed_fft_base_function<0,2>()(a2,a3,value_type(complex_type(SIN12,-SIN20),complex_type(-SIN20,-SIN12)),a2,a3, b2,b3,value_type(complex_type(SIN10,-SIN22),complex_type(-SIN22,-SIN10)),b2,b3); Y.first (10,a2); Y.second(42,b2); Y.first (26,a3); Y.second(58,b3);
  }
  {
    value_type a0,a1,a2,a3, b0,b1,b2,b3;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=t1+t5; v5=imul(SIN16,t11-t15); v1=v4+v5; v0=v4-v5; v4=SIN08*t31-imul(SIN24,t31); v5=SIN08*t21+imul(SIN24,t21); v2=v4+v5; v3=imul(v5-v4); } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5;  v4=t2-t7; v5=imul(t13-t9); v0=v4+v5; v1=v4-v5; v4=t25-t33; v4=SIN04*v4-imul(SIN28,v4); v5=t16-t24; v5=SIN04*v5+imul(SIN28,v5); v2=v4+v5; v3=imul(v5-v4); } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3; }
    complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN20,-SIN12),complex_type(-SIN12,-SIN20)),a0,a1, b0,b1,value_type(complex_type(SIN18,-SIN14),complex_type(-SIN14,-SIN18)),b0,b1); Y.first ( 6,a0); Y.second(38,b0); Y.first (22,a1); Y.second(54,b1);
    complex_crossed_fft_base_function<0,2>()(a2,a3,value_type(complex_type(SIN04,-SIN28),complex_type(-SIN28,-SIN04)),a2,a3, b2,b3,value_type(complex_type(SIN02,-SIN30),complex_type(-SIN30,-SIN02)),b2,b3); Y.first (14,a2); Y.second(46,b2); Y.first (30,a3); Y.second(62,b3);
  }
}

template<class G1,class S>
inline void complex_fft_32_2(const Vector<G1> &Xa, const Vector<G1> &Xb, const S &Ya)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type ta00,ta01,ta02,ta03,ta04,ta05,ta06,ta07,ta08,ta09,ta10,ta11,ta12,ta13,ta14,ta15,ta16,ta17,ta18,ta19,ta20,ta21,ta22,ta23,ta24,ta25,ta26,ta27,ta28,ta29,ta30,ta31,ta32,ta33;
  value_type tb00,tb01,tb02,tb03,tb04,tb05,tb06,tb07,tb08,tb09,tb10,tb11,tb12,tb13,tb14,tb15,tb16,tb17,tb18,tb19,tb20,tb21,tb22,tb23,tb24,tb25,tb26,tb27,tb28,tb29,tb30,tb31,tb32,tb33;

  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xa[ 0]; v5=Xa[16]; v0=v4+v5; v1=v4-v5; v4=Xa[ 8]; v5=Xa[24]; v2=v4+v5; v3=imul(v4-v5); } ta00 =v0+v2; ta01=v0-v2; ta02=v1+v3; ta03=v1-v3; { value_type v4,v5; v4=Xa[ 0]; v5=Xa[16]; v0=v4+v5; v1=v4-v5; v4=Xa[ 8]; v5=Xa[24]; v2=v4+v5; v3=imul(v4-v5); } ta00 =v0+v2; ta01=v0-v2; ta02=v1+v3; ta03=v1-v3; }
  { value_type v0,v1,v2,v3,v6; { value_type v4,v5; v4=Xa[ 1]; v5=Xa[17]; v0=v4+v5; v1=v4-v5; v4=Xa[ 9]; v5=Xa[25]; v2=v4+v5; v3=imul(v4-v5); } ta25=v1+v3; ta26=v1-v3; ta27=v0+v2; v6=v0-v2; { value_type v4,v5; v4=Xa[ 5]; v5=Xa[21]; v0=v4+v5; v1=v4-v5; v4=Xa[29]; v5=Xa[13]; v2=v4+v5; v3=     v4-v5 ; } ta28=v0+v2; ta29=ta27-ta28; v2=imul(v2-v0); ta30=v6+v2; ta31=v6-v2; v6=v3+v1; v2=imul(v3-v1); ta32=SIN16*(v2+v6); ta33=SIN16*(v2-v6); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xa[ 4]; v5=Xa[20]; v0=v4+v5; v1=v4-v5; v4=Xa[28]; v5=Xa[12]; v2=v4+v5; v3=     v4-v5 ; } ta04 =v0+v2; ta05=imul(v0-v2); v0=v3+v1; v2=imul(v3-v1); ta06=SIN16*(v2+v0); ta07=SIN16*(v2-v0); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xa[30]; v5=Xa[14]; v0=v4+v5; v1=v4-v5; v4=Xa[ 6]; v5=Xa[22]; v2=v4+v5; v3=imul(v4-v5); } ta12=v0+v2; v0=v0-v2; ta15=add_imul(v0,v0); v0=v1+v3; v2=v1-v3; ta13=SIN08*v0+imul(SIN24,v0); ta14=SIN24*v2+imul(SIN08,v2); }
  { value_type v0,v1,v2,v3,v6; { value_type v4,v5; v4=Xa[31]; v5=Xa[15]; v0=v4+v5; v1=v4-v5; v4=Xa[ 7]; v5=Xa[23]; v2=v4+v5; v3=imul(v4-v5); } ta16=v1+v3; ta17=v1-v3; ta18=v0+v2; v6=v0-v2; { value_type v4,v5; v4=Xa[ 3]; v5=Xa[19]; v0=v4+v5; v1=v4-v5; v4=Xa[27]; v5=Xa[11]; v2=v4+v5; v3=     v4-v5 ; } ta19=v0+v2; v0=imul(v0-v2); ta20=ta18-ta19; ta21=v6+v0; ta22=v6-v0; v6=v3+v1; v2=imul(v3-v1); ta23=SIN16*(v2+v6); ta24=SIN16*(v2-v6); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xa[ 2]; v5=Xa[18]; v0=v4+v5; v1=v4-v5; v4=Xa[10]; v5=Xa[26]; v2=v4+v5; v3=imul(v4-v5); } ta08 =v0+v2; v0=v0-v2; ta11=v0-imul(v0)    ; v0=v1+v3; v2=v1-v3; ta09 =SIN08*v0-imul(SIN24,v0); ta10=SIN24*v2-imul(SIN08,v2); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xb[ 0]; v5=Xb[16]; v0=v4+v5; v1=v4-v5; v4=Xb[ 8]; v5=Xb[24]; v2=v4+v5; v3=imul(v4-v5); } tb00 =v0+v2; tb01=v0-v2; tb02=v1+v3; tb03=v1-v3; { value_type v4,v5; v4=Xb[ 0]; v5=Xb[16]; v0=v4+v5; v1=v4-v5; v4=Xb[ 8]; v5=Xb[24]; v2=v4+v5; v3=imul(v4-v5); } tb00 =v0+v2; tb01=v0-v2; tb02=v1+v3; tb03=v1-v3; }
  { value_type v0,v1,v2,v3,v6; { value_type v4,v5; v4=Xb[ 1]; v5=Xb[17]; v0=v4+v5; v1=v4-v5; v4=Xb[ 9]; v5=Xb[25]; v2=v4+v5; v3=imul(v4-v5); } tb25=v1+v3; tb26=v1-v3; tb27=v0+v2; v6=v0-v2; { value_type v4,v5; v4=Xb[ 5]; v5=Xb[21]; v0=v4+v5; v1=v4-v5; v4=Xb[29]; v5=Xb[13]; v2=v4+v5; v3=     v4-v5 ; } tb28=v0+v2; tb29=tb27-tb28; v2=imul(v2-v0); tb30=v6+v2; tb31=v6-v2; v6=v3+v1; v2=imul(v3-v1); tb32=SIN16*(v2+v6); tb33=SIN16*(v2-v6); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xb[ 4]; v5=Xb[20]; v0=v4+v5; v1=v4-v5; v4=Xb[28]; v5=Xb[12]; v2=v4+v5; v3=     v4-v5 ; } tb04 =v0+v2; tb05=imul(v0-v2); v0=v3+v1; v2=imul(v3-v1); tb06=SIN16*(v2+v0); tb07=SIN16*(v2-v0); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xb[30]; v5=Xb[14]; v0=v4+v5; v1=v4-v5; v4=Xb[ 6]; v5=Xb[22]; v2=v4+v5; v3=imul(v4-v5); } tb12=v0+v2; v0=v0-v2; tb15=add_imul(v0,v0); v0=v1+v3; v2=v1-v3; tb13=SIN08*v0+imul(SIN24,v0); tb14=SIN24*v2+imul(SIN08,v2); }
  { value_type v0,v1,v2,v3,v6; { value_type v4,v5; v4=Xb[31]; v5=Xb[15]; v0=v4+v5; v1=v4-v5; v4=Xb[ 7]; v5=Xb[23]; v2=v4+v5; v3=imul(v4-v5); } tb16=v1+v3; tb17=v1-v3; tb18=v0+v2; v6=v0-v2; { value_type v4,v5; v4=Xb[ 3]; v5=Xb[19]; v0=v4+v5; v1=v4-v5; v4=Xb[27]; v5=Xb[11]; v2=v4+v5; v3=     v4-v5 ; } tb19=v0+v2; v0=imul(v0-v2); tb20=tb18-tb19; tb21=v6+v0; tb22=v6-v0; v6=v3+v1; v2=imul(v3-v1); tb23=SIN16*(v2+v6); tb24=SIN16*(v2-v6); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=Xb[ 2]; v5=Xb[18]; v0=v4+v5; v1=v4-v5; v4=Xb[10]; v5=Xb[26]; v2=v4+v5; v3=imul(v4-v5); } tb08 =v0+v2; v0=v0-v2; tb11=v0-imul(v0)    ; v0=v1+v3; v2=v1-v3; tb09 =SIN08*v0-imul(SIN24,v0); tb10=SIN24*v2-imul(SIN08,v2); }

  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta00+ta04; v5=                          ta12+ta08  ; v0=v4+v5; v1=v4-v5; v4=ta27+ta28;                                                              v5=ta18+ta19;                                                             v2=                     v5+v4 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb00+tb04; v5=                          tb12+tb08  ; v0=v4+v5; v1=v4-v5; v4=tb27+tb28;                                                              v5=tb18+tb19;                                                             v2=                     v5+v4 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,                                 a0,a1,b0,b1,                                  b0,b1);
    Ya.first0 (a0); Ya.first(16,a1); Ya.second(32,b0); Ya.second(48,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN16,-SIN16),a2,a3,b2,b3,complex_type(-SIN16,-SIN16),b2,b3);
    Ya.first(8,a2); Ya.first(24,a3); Ya.second(40,b2); Ya.second(56,b3);
  }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta00-ta04; v5=imul(                     ta12-ta08 ); v0=v4+v5; v1=v4-v5; v4=ta20+ta29;                                                              v5=imul(ta20-ta29);                                                       v2=SIN16*(v4+v5); v3=SIN16*(v5-v4); } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb00-tb04; v5=imul(                     tb12-tb08 ); v0=v4+v5; v1=v4-v5; v4=tb20+tb29;                                                              v5=imul(tb20-tb29);                                                       v2=SIN16*(v4+v5); v3=SIN16*(v5-v4); } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN24,-SIN08),a0,a1,b0,b1,complex_type(-SIN08,-SIN24),b0,b1);
    Ya.first( 4,a0); Ya.first(20,a1); Ya.second(36,b0); Ya.second(52,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN08,-SIN24),a2,a3,b2,b3,complex_type(-SIN24,-SIN08),b2,b3);
    Ya.first(12,a2); Ya.first(28,a3); Ya.second(44,b2); Ya.second(60,b3);
   }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta01-ta05; v5=     SIN16*(ta11+ta15); v0=v4+v5; v1=v4-v5;             v4=SIN24*ta30-imul(SIN08,ta30); v5=SIN24*ta22+imul(SIN08,ta22);              v2=                     v4+v5 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb01-tb05; v5=     SIN16*(tb11+tb15); v0=v4+v5; v1=v4-v5;             v4=SIN24*tb30-imul(SIN08,tb30); v5=SIN24*tb22+imul(SIN08,tb22);              v2=                     v4+v5 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN28,-SIN04),a0,a1,b0,b1,complex_type(-SIN04,-SIN28),b0,b1);
    Ya.first( 2,a0); Ya.first(18,a1); Ya.second(34,b0); Ya.second(50,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN12,-SIN20),a2,a3,b2,b3,complex_type(-SIN20,-SIN12),b2,b3);
    Ya.first(10,a2); Ya.first(26,a3); Ya.second(42,b2); Ya.second(58,b3);
  }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta01+ta05; v5=imul(SIN16, ta11-ta15); v1=v4+v5; v0=v4-v5;             v4=SIN08*ta31-imul(SIN24,ta31); v5=SIN08*ta21+imul(SIN24,ta21);              v2=                     v4+v5 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb01+tb05; v5=imul(SIN16, tb11-tb15); v1=v4+v5; v0=v4-v5;             v4=SIN08*tb31-imul(SIN24,tb31); v5=SIN08*tb21+imul(SIN24,tb21);              v2=                     v4+v5 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN20,-SIN12),a0,a1,b0,b1,complex_type(-SIN12,-SIN20),b0,b1);
    Ya.first( 6,a0); Ya.first(22,a1); Ya.second(38,b0); Ya.second(54,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN04,-SIN28),a2,a3,b2,b3,complex_type(-SIN28,-SIN04),b2,b3);
    Ya.first(14,a2); Ya.first(30,a3); Ya.second(46,b2); Ya.second(62,b3);
  }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta03+ta06; v5=                          ta14+ta10 ; v0=v4+v5; v1=v4-v5; v4=ta26+ta32; v4=SIN28*v4-imul(SIN04,v4) ; v5=ta17+ta23; v5=SIN28*v5+imul(SIN04,v5); v2=                     v5+v4 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb03+tb06; v5=                          tb14+tb10 ; v0=v4+v5; v1=v4-v5; v4=tb26+tb32; v4=SIN28*v4-imul(SIN04,v4) ; v5=tb17+tb23; v5=SIN28*v5+imul(SIN04,v5); v2=                     v5+v4 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN30,-SIN02),a0,a1,b0,b1,complex_type(-SIN02,-SIN30),b0,b1);
    Ya.first( 1,a0); Ya.first(17,a1); Ya.second(33,b0); Ya.second(49,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN14,-SIN18),a2,a3,b2,b3,complex_type(-SIN18,-SIN14),b2,b3);
    Ya.first( 9,a2); Ya.first(25,a3); Ya.second(41,b2); Ya.second(57,b3);
  }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta02+ta07; v5=                          ta13+ta09  ; v0=v4+v5; v1=v4-v5; v4=ta25+ta33; v4=SIN20*v4-imul(SIN12,v4) ; v5=ta16+ta24; v5=SIN20*v5+imul(SIN12,v5); v2=                     v5+v4 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb02+tb07; v5=                          tb13+tb09  ; v0=v4+v5; v1=v4-v5; v4=tb25+tb33; v4=SIN20*v4-imul(SIN12,v4) ; v5=tb16+tb24; v5=SIN20*v5+imul(SIN12,v5); v2=                     v5+v4 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN26,-SIN06),a0,a1,b0,b1,complex_type(-SIN06,-SIN26),b0,b1);
    Ya.first( 3,a0); Ya.first(19,a1); Ya.second(35,b0); Ya.second(51,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN10,-SIN22),a2,a3,b2,b3,complex_type(-SIN22,-SIN10),b2,b3);
    Ya.first(11,a2); Ya.first(27,a3); Ya.second(43,b2); Ya.second(59,b3);
  }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta03-ta06; v5=imul(                     ta14-ta10); v0=v4+v5; v1=v4-v5; v4=ta26-ta32; v4=SIN12*v4-imul(SIN20,v4) ; v5=ta17-ta23; v5=SIN12*v5+imul(SIN20,v5); v2=                     v4+v5 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb03-tb06; v5=imul(                     tb14-tb10); v0=v4+v5; v1=v4-v5; v4=tb26-tb32; v4=SIN12*v4-imul(SIN20,v4) ; v5=tb17-tb23; v5=SIN12*v5+imul(SIN20,v5); v2=                     v4+v5 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN22,-SIN10),a0,a1,b0,b1,complex_type(-SIN10,-SIN22),b0,b1);
    Ya.first( 5,a0); Ya.first(21,a1); Ya.second(37,b0); Ya.second(53,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN06,-SIN26),a2,a3,b2,b3,complex_type(-SIN26,-SIN06),b2,b3);
    Ya.first(13,a2); Ya.first(29,a3); Ya.second(45,b2); Ya.second(61,b3);
  }
  {
    value_type v0,v1,v2,v3, a0,a1,a2,a3,b0,b1,b2,b3;
    { value_type v4,v5; v4=ta02-ta07; v5=imul(                     ta13-ta09 ); v0=v4+v5; v1=v4-v5; v4=ta25-ta33; v4=SIN04*v4-imul(SIN28,v4) ; v5=ta16-ta24; v5=SIN04*v5+imul(SIN28,v5); v2=                     v4+v5 ; v3=imul(v5-v4);                 } a0=v0+v2; a1=v0-v2; a2=v1+v3; a3=v1-v3;
    { value_type v4,v5; v4=tb02-tb07; v5=imul(                     tb13-tb09 ); v0=v4+v5; v1=v4-v5; v4=tb25-tb33; v4=SIN04*v4-imul(SIN28,v4) ; v5=tb16-tb24; v5=SIN04*v5+imul(SIN28,v5); v2=                     v4+v5 ; v3=imul(v5-v4);                 } b0=v0+v2; b1=v0-v2; b2=v1+v3; b3=v1-v3;
    complex_crossed_fft_base_function<0,2>().eval(a0,a1,complex_type(SIN18,-SIN14),a0,a1,b0,b1,complex_type(-SIN14,-SIN18),b0,b1);
    Ya.first( 7,a0); Ya.first(23,a1); Ya.second(39,b0); Ya.second(55,b1);
    complex_crossed_fft_base_function<0,2>().eval(a2,a3,complex_type(SIN02,-SIN30),a2,a3,b2,b3,complex_type(-SIN30,-SIN02),b2,b3);
    Ya.first(15,a2); Ya.first(31,a3); Ya.second(47,b2); Ya.second(63,b3);
  }
}

template<class G,class S>
inline void complex_fft_64_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63;

  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 0]; u1=X[32]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[16]; u1=X[48]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[ 8]; u1=X[40]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[56]; u1=X[24]; v6=u0+u1; v7=u0-u1; } { t00=add_imul(v1,v3); t01=sub_imul(v1,v3); } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t02=u0+u1; t03=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v6-v4); t04=u0-u1; t05=u0+u1; } { value_type u0,u1; u0=add_imul(v7,v7); u1=sub_imul(v5,v5); t06=imul(SIN16,u0-u1); t07=SIN16*(u1+u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 4]; u1=X[36]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[20]; u1=X[52]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[60]; u1=X[28]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[12]; u1=X[44]; v6=u0+u1; v7=u0-u1; } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t08=u0+u1; t09=imul(u1-u0); } { value_type u0; u0=v0-v2; t10=sub_imul(u0,u0); } { value_type u0; u0=v4-v6; t11=add_imul(u0,u0); } { value_type u0; u0=add_imul(v1,v3); t12=SIN08*u0-imul(SIN24,u0); } { value_type u0; u0=sub_imul(v1,v3); t13=SIN24*u0-imul(SIN08,u0); } { value_type u0; u0=add_imul(v5,v7); t14=SIN08*u0+imul(SIN24,u0); } { value_type u0; u0=sub_imul(v5,v7); t15=SIN24*u0+imul(SIN08,u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 2]; u1=X[34]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[18]; u1=X[50]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[58]; u1=X[26]; v4=u0+u1; { value_type w0; w0=u0-u1; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[10]; u1=X[42]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t16=u0+u1; t17=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t18=u0-u1; t19=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t20=u0-u1; t21=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t22=u0-u1; t23=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[63]; u1=X[31]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[15]; u1=X[47]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[55]; u1=X[23]; v4=u0+u1; { value_type w0; w0=u0-u1; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[ 7]; u1=X[39]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t24=u0+u1; t25=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t26=u0-u1; t27=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t28=u0-u1; t29=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t30=u0-u1; t31=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[62]; u1=X[30]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[14]; u1=X[46]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[54]; u1=X[22]; v4=u0+u1; { value_type w0; w0=u0-u1; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[ 6]; u1=X[38]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t32=u0+u1; t33=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t34=u0-u1; t35=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t36=u0-u1; t37=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t38=u0-u1; t39=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 1]; u1=X[33]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[17]; u1=X[49]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[57]; u1=X[25]; v4=u0+u1; { value_type w0; w0=u0-u1; v5=add_imul(w0,w0); } } { value_type u0,u1; u0=X[ 9]; u1=X[41]; v6=u0+u1; { value_type w0; w0=u0-u1; v7=sub_imul(w0,w0); } } { value_type u0,u1; u0=v0+v2; u1=v4+v6; t40=u0+u1; t41=u0-u1; } { value_type u0,u1; u0=v0-v2; u1=imul(v4-v6); t42=u0-u1; t43=u0+u1; } { value_type u0,u1; u0=add_imul(v1,v3); u1=imul(SIN16,v5-v7); t44=u0-u1; t45=u0+u1; } { value_type u0,u1; u0=sub_imul(v1,v3); u1=SIN16*(v5+v7); t46=u0-u1; t47=u0+u1; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; { value_type w0,w1; w0=X[ 5]; w1=X[37]; v0=w0+w1; u0=w0-w1; } { value_type w0,w1; w0=X[21]; w1=X[53]; v1=w0+w1; u1=w0-w1; } v2=add_imul(u0,u1); v3=sub_imul(u0,u1); } { value_type u0,u1; { value_type w0,w1; w0=X[61]; w1=X[29]; v4=w0+w1; u0=w0-w1; } { value_type w0,w1; w0=X[13]; w1=X[45]; v5=w0+w1; u1=w0-w1; } v6=add_imul(u0,u1); v7=sub_imul(u0,u1); } { value_type u0,u1; u0=v0+v1; u1=v4+v5; t48=u0+u1; t49=imul(u1-u0); } { value_type u0,u1; u0=v0-v1; u0=sub_imul(u0,u0); u1=v4-v5; u1=add_imul(u1,u1); t50=imul(SIN16,u1-u0); t51=SIN16*(u0+u1); } { value_type u0,u1; u0=SIN08*v6+imul(SIN24,v6); u1=SIN08*v2-imul(SIN24,v2); t53=u0+u1; t52=imul(u0-u1); } { value_type u0,u1; u0=SIN24*v3-imul(SIN08,v3); u1=SIN24*v7+imul(SIN08,v7); t55=u1+u0; t54=imul(u1-u0); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; { value_type w0,w1; w0=X[ 3]; w1=X[35]; v0=w0+w1; u0=w0-w1; } { value_type w0,w1; w0=X[19]; w1=X[51]; v1=w0+w1; u1=w0-w1; } v2=add_imul(u0,u1); v3=sub_imul(u0,u1); } { value_type u0,u1; { value_type w0,w1; w0=X[59]; w1=X[27]; v4=w0+w1; u0=w0-w1; } { value_type w0,w1; w0=X[11]; w1=X[43]; v5=w0+w1; u1=w0-w1; } v6=add_imul(u0,u1); v7=sub_imul(u0,u1); } { value_type u0,u1; u0=v0+v1; u1=v4+v5; t56=u0+u1; t57=imul(u1-u0); } { value_type u0,u1; u0=v0-v1; u0=sub_imul(u0,u0); u1=v4-v5; u1=add_imul(u1,u1); t58=imul(SIN16,u1-u0); t59=SIN16*(u1+u0); } { value_type u0,u1; u0=SIN08*v2-imul(SIN24,v2); u1=SIN08*v6+imul(SIN24,v6); t61=u0+u1; t60=imul(u1-u0); } { value_type u0,u1; u0=SIN24*v7+imul(SIN08,v7); u1=SIN24*v3-imul(SIN08,v3); t63=u0+u1; t62=imul(u0-u1); } }
  {
    value_type v0,v1,v2,v3; { value_type u0,u1; u0=t02+t08; u1=     t16+t32;  v0=u0+u1; v1=u0-u1; } { value_type u0,u1;u0=t40+t48; u1=t24+t56; v2=imul(u1-u0); v3=u0+u1; }
    value_type w0,w1,w2,w3,w4,w5,w6,w7; { value_type u0,u1; u0=t01+t07; u1=     t13+t15 ; w0=u0-u1; w1=u0+u1; } { value_type u0; u0=t47+t55; w2=SIN30*u0-imul(SIN02,u0); } { value_type u0; u0=t31+t63; w3=SIN30*u0+imul(SIN02,u0); } { value_type u0,u1; u0=SIN28*t23-imul(SIN04,t23); u1=SIN28*t39+imul(SIN04,t39); w4=imul(u1-u0); w5=u1+u0; } { value_type u0; u0=t47-t55; w6=SIN14*u0-imul(SIN18,u0); } { value_type u0; u0=t31-t63; w7=SIN14*u0+imul(SIN18,u0); }
    { value_type a0,a1,b0,b1; a0=v0+v3; a1=v0-v3; { value_type u0,u1; u0=w1+w5; u1=     w2+w3;  b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,                                                                                           a0,a1, b0,b1,value_type(complex_type(SIN31,-SIN01),complex_type(-SIN01,-SIN31)),b0,b1); Y.first0(   a0); Y.second(64,b0); Y.first (32,a1); Y.second( 96,b1); }
    { value_type a0,a1,b0,b1; a0=v1+v2; a1=v1-v2; { value_type u0,u1; u0=w1-w5; u1=imul(w3-w2); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN16,-SIN16),complex_type(-SIN16,-SIN16)),a0,a1, b0,b1,value_type(complex_type(SIN15,-SIN17),complex_type(-SIN17,-SIN15)),b0,b1); Y.first (16,a0); Y.second(80,b0); Y.first (48,a1); Y.second(112,b1); }
    value_type v4,v5,v6,v7; { value_type u0,u1; u0=t02-t08; u1=imul(t32-t16); v4=u0+u1; v5=u0-u1; } { value_type u0; u0=t40-t48; v6=sub_imul(u0,u0); } { value_type u0; u0=t24-t56; v7=add_imul(u0,u0); }
    { value_type a0,a1,b0,b1; { value_type u0; u0=    SIN16*(v6+v7); a0=v4+u0; a1=v4-u0; } { value_type u0,u1; u0=w0+w4; u1=     w6+w7;  b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN24,-SIN08),complex_type(-SIN08,-SIN24)),a0,a1, b0,b1,value_type(complex_type(SIN23,-SIN09),complex_type(-SIN09,-SIN23)),b0,b1); Y.first ( 8,a0); Y.second(72,b0); Y.first (40,a1); Y.second(104,b1); }
    { value_type a0,a1,b0,b1; { value_type u0; u0=imul(SIN16,v7-v6); a0=v5+u0; a1=v5-u0; } { value_type u0,u1; u0=w0-w4; u1=imul(w7-w6); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN08,-SIN24),complex_type(-SIN24,-SIN08)),a0,a1, b0,b1,value_type(complex_type(SIN07,-SIN25),complex_type(-SIN25,-SIN07)),b0,b1); Y.first (24,a0); Y.second(88,b0); Y.first (56,a1); Y.second(120,b1); }
  }
  {
    value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0; u0=     SIN16*(t10+t11); v0=t05-u0; v1=t05+u0; } { value_type u0; u0=t43+t51; v2=SIN28*u0-imul(SIN04,u0); } { value_type u0; u0=t27+t59; v3=SIN28*u0+imul(SIN04,u0); } { value_type u0,u1; u0=SIN24*t19-imul(SIN08,t19); u1=SIN24*t35+imul(SIN08,t35); v4=imul(u1-u0); v5=u0+u1; } { value_type u0; u0=t43-t51; v6=SIN12*u0-imul(SIN20,u0); } { value_type u0; u0=t27-t59; v7=SIN12*u0+imul(SIN20,u0); }
    value_type w0,w1,w2,w3,w4,w5,w6,w7; { value_type u0,u1; u0=t00+t06; u1=     t12+t14 ; w0=u0-u1; w1=u0+u1; } { value_type u0; u0=t45+t53; w2=SIN26*u0-imul(SIN06,u0); } { value_type u0; u0=t29+t61; w3=SIN26*u0+imul(SIN06,u0); } { value_type u0,u1; u0=SIN20*t21-imul(SIN12,t21); u1=SIN20*t37+imul(SIN12,t37); w4=imul(u1-u0); w5=u1+u0; } { value_type u0; u0=t45-t53; w6=SIN10*u0-imul(SIN22,u0); } { value_type u0; u0=t29-t61; w7=SIN10*u0+imul(SIN22,u0); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v1+v5; u1=     v2+v3;  a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w1+w5; u1=     w2+w3;  b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN30,-SIN02),complex_type(-SIN02,-SIN30)),a0,a1, b0,b1,value_type(complex_type(SIN29,-SIN03),complex_type(-SIN03,-SIN29)),b0,b1); Y.first ( 2,a0); Y.second(66,b0); Y.first (34,a1); Y.second( 98,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v0+v4; u1=     v6+v7;  a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w0+w4; u1=     w6+w7;  b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN22,-SIN10),complex_type(-SIN10,-SIN22)),a0,a1, b0,b1,value_type(complex_type(SIN21,-SIN11),complex_type(-SIN11,-SIN21)),b0,b1); Y.first (10,a0); Y.second(74,b0); Y.first (42,a1); Y.second(106,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w1-w5; u1=imul(w3-w2); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN14,-SIN18),complex_type(-SIN18,-SIN14)),a0,a1, b0,b1,value_type(complex_type(SIN13,-SIN19),complex_type(-SIN19,-SIN13)),b0,b1); Y.first (18,a0); Y.second(82,b0); Y.first (50,a1); Y.second(114,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w0-w4; u1=imul(w7-w6); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN06,-SIN26),complex_type(-SIN26,-SIN06)),a0,a1, b0,b1,value_type(complex_type(SIN05,-SIN27),complex_type(-SIN27,-SIN05)),b0,b1); Y.first (26,a0); Y.second(90,b0); Y.first (58,a1); Y.second(122,b1); }
  }
  {
    value_type v0,v1,v2,v3,v4,v5,v6,v7; v0=t03-t09; v1=t03+t09;                                                 { value_type u0; u0=t41+t49; v2=SIN24*u0-imul(SIN08,u0); } { value_type u0; u0=t25+t57; v3=SIN24*u0+imul(SIN08,u0); } { value_type u0,u1; u0=sub_imul(t17,t17); u1=add_imul(t33,t33); v4=imul(SIN16,u1-u0); v5=SIN16*(u0+u1); }               { value_type u0; u0=t41-t49; v6=SIN08*u0-imul(SIN24,u0); } { value_type u0; u0=t25-t57; v7=SIN08*u0+imul(SIN24,u0); }
    value_type w0,w1,w2,w3,w4,w5,w6,w7; { value_type u0,u1; u0=t01-t07; u1=imul(t15-t13); w0=u0-u1; w1=u0+u1; } { value_type u0; u0=t46+t54; w2=SIN22*u0-imul(SIN10,u0); } { value_type u0; u0=t30+t62; w3=SIN22*u0+imul(SIN10,u0); } { value_type u0,u1; u0=SIN12*t22-imul(SIN20,t22); u1=SIN12*t38+imul(SIN20,t38); w4=imul(u1-u0); w5=u1+u0; } { value_type u0; u0=t46-t54; w6=SIN06*u0-imul(SIN26,u0); } { value_type u0; u0=t30-t62; w7=SIN06*u0+imul(SIN26,u0); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v1+v5; u1=     v2+v3;  a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w1+w5; u1=     w2+w3;  b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN28,-SIN04),complex_type(-SIN04,-SIN28)),a0,a1, b0,b1,value_type(complex_type(SIN27,-SIN05),complex_type(-SIN05,-SIN27)),b0,b1); Y.first ( 4,a0); Y.second(68,b0); Y.first (36,a1); Y.second(100,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v0+v4; u1=     v6+v7;  a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w0+w4; u1=     w6+w7;  b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN20,-SIN12),complex_type(-SIN12,-SIN20)),a0,a1, b0,b1,value_type(complex_type(SIN19,-SIN13),complex_type(-SIN13,-SIN19)),b0,b1); Y.first (12,a0); Y.second(76,b0); Y.first (44,a1); Y.second(108,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w1-w5; u1=imul(w3-w2); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN12,-SIN20),complex_type(-SIN20,-SIN12)),a0,a1, b0,b1,value_type(complex_type(SIN11,-SIN21),complex_type(-SIN21,-SIN11)),b0,b1); Y.first (20,a0); Y.second(84,b0); Y.first (52,a1); Y.second(116,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w0-w4; u1=imul(w7-w6); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN04,-SIN28),complex_type(-SIN28,-SIN04)),a0,a1, b0,b1,value_type(complex_type(SIN03,-SIN29),complex_type(-SIN29,-SIN03)),b0,b1); Y.first (28,a0); Y.second(92,b0); Y.first (60,a1); Y.second(124,b1); }
  }
  {
    value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0; u0=imul(SIN16 ,t11-t10); v0=t04-u0; v1=t04+u0; } { value_type u0; u0=t42+t50; v2=SIN20*u0-imul(SIN12,u0); } { value_type u0; u0=t26+t58; v3=SIN20*u0+imul(SIN12,u0); } { value_type u0,u1; u0=SIN08*t18-imul(SIN24,t18); u1=SIN08*t34+imul(SIN24,t34); v4=imul(u1-u0); v5=u0+u1; } { value_type u0; u0=t42-t50; v6=SIN04*u0-imul(SIN28,u0); } { value_type u0; u0=t26-t58; v7=SIN04*u0+imul(SIN28,u0); }
    value_type w0,w1,w2,w3,w4,w5,w6,w7; { value_type u0,u1; u0=t00-t06; u1=imul(t14-t12); w0=u0-u1; w1=u0+u1; } { value_type u0; u0=t44+t52; w2=SIN18*u0-imul(SIN14,u0); } { value_type u0; u0=t28+t60; w3=SIN18*u0+imul(SIN14,u0); } { value_type u0,u1; u0=SIN04*t20-imul(SIN28,t20); u1=SIN04*t36+imul(SIN28,t36); w4=imul(u1-u0); w5=u1+u0; } { value_type u0; u0=t44-t52; w6=SIN02*u0-imul(SIN30,u0); } { value_type u0; u0=t28-t60; w7=SIN02*u0+imul(SIN30,u0); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v1+v5; u1=     v2+v3;  a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w1+w5; u1=w2+w3; b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN26,-SIN06),complex_type(-SIN06,-SIN26)),a0,a1, b0,b1,value_type(complex_type(SIN25,-SIN07),complex_type(-SIN07,-SIN25)),b0,b1); Y.first ( 6,a0); Y.second(70,b0); Y.first (38,a1); Y.second(102,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v0+v4; u1=     v6+v7;  a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w0+w4; u1=w6+w7; b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN18,-SIN14),complex_type(-SIN14,-SIN18)),a0,a1, b0,b1,value_type(complex_type(SIN17,-SIN15),complex_type(-SIN15,-SIN17)),b0,b1); Y.first (14,a0); Y.second(78,b0); Y.first (46,a1); Y.second(110,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v1-v5; u1=imul(v3-v2); a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w1-w5; u1=imul(w3-w2); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN10,-SIN22),complex_type(-SIN22,-SIN10)),a0,a1, b0,b1,value_type(complex_type(SIN09,-SIN23),complex_type(-SIN23,-SIN09)),b0,b1); Y.first (22,a0); Y.second(86,b0); Y.first (54,a1); Y.second(118,b1); }
    { value_type a0,a1,b0,b1; { value_type u0,u1; u0=v0-v4; u1=imul(v7-v6); a0=u0+u1; a1=u0-u1; } { value_type u0,u1; u0=w0-w4; u1=imul(w7-w6); b0=u0+u1; b1=u0-u1; } complex_crossed_fft_base_function<0,2>()(a0,a1,value_type(complex_type(SIN02,-SIN30),complex_type(-SIN30,-SIN02)),a0,a1, b0,b1,value_type(complex_type(SIN01,-SIN31),complex_type(-SIN31,-SIN01)),b0,b1); Y.first (30,a0); Y.second(94,b0); Y.first (62,a1); Y.second(126,b1); }
  }
}


template<> struct complex_crossed_fft_base_function< 1,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_1_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_1_2 (Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function< 2,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_2_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_2_2 (Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function< 3,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_3_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_3_2 (Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function< 4,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_4_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_4_2 (Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function< 5,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_5_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_5_2 (Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function< 8,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_8_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_8_2 (Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function<16,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_16_2(X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_16_2(Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function<32,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_32_2(X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_32_2(Xa,Xb, f); } };
template<> struct complex_crossed_fft_base_function<64,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X , const F &f) { complex_fft_64_2(X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_64_2(Xa,Xb, f); } };

template<int M,int N> struct complex_crossed_fft_function
{
  template<class G1,         class G3,class F> inline void operator()(const Vector<G1> &X ,                       Vector<G3> &Y                  , const F &f) { complex_crossed_fft_base_function<M,N>()(X  ,typename F::template rebind<M*N>::other(f).make0(Y)); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X , const Vector<G2> &W , Vector<G3> &Y                  , const F &f) { complex_crossed_fft_base_function<M,N>()(X*W,typename F::template rebind<M*N>::other(f).make(Y)); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X , const Vector<G2> &W , Vector<G3> &Y0 , Vector<G3> &Y1, const F &f) { complex_crossed_fft_base_function<M,N>()(X*W,typename F::template rebind<M*N>::other(f).make(Y0,Y1)); }

  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G2> &Wa, Vector<G3> &Ya ,                  const Vector<G1> &Xb, const Vector<G2> &Wb, Vector<G3> &Yb                  , const F &f) { complex_crossed_fft_base_function<M,N>()(Xa*Wa, Xb*Wb, typename F::template rebind<M*N>::other(f).simd_make(Ya     )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G2> &Wa, Vector<G3> &Ya0, Vector<G3> &Ya1, const Vector<G1> &Xb, const Vector<G2> &Wb, Vector<G3> &Yb0, Vector<G3> &Yb1, const F &f) { complex_crossed_fft_base_function<M,N>()(Xa*Wa, Xb*Wb, typename F::template rebind<M*N>::other(f).simd_make(Ya0,Ya1)); }
};



template<class G,class V>
void fft_twiddles(Matrix<G> &W, const V &x, int m, int n, int l)
{
  typedef V real_type;
  typedef complex<real_type> complex_type;

  complex_type w0=complex_type(x,0);

#ifdef FFT_QUICK_TWIDDLES
  complex_type w1=x*polar<real_type>(1,-2*PI/l);
  complex_type wj=w0;
  for (int i=0; i<m; ++i) W(i,0)=w0;
  for (int j=1; j<n; ++j)
  {
    wj*=w1;
    complex_type w=wj;
    W(0,j)=w0;
    W(1,j)=w;
    for (int i=2; i<m; ++i)
      W(i,j)=(w*=wj);
  }
#else
  for (int i=0; i<m; ++i) W(i,0)=w0;
  for (int j=1; j<n; ++j) W(0,j)=w0;

  //for (int i=1; i<m; ++i)
  //  for (int j=1; j<n; ++j)
  //    W(i,j)=polar<real_type>(x,-(2*real_type(PI)*i*j)/l);

  int mn=min(m,n);
  for (int j=1; j<mn; ++j)
    for (int i=j; i<mn; ++i)
      W(i,j)=W(j,i)=polar<real_type>(x,-(2*real_type(PI)*i*j)/l);

  for (int j=mn; j<n; ++j)
    for (int i=1; i<m; ++i)
      W(i,j)=polar<real_type>(x,-(2*real_type(PI)*i*j)/l);

  for (int i=mn; i<m; ++i)
    for (int j=1; j<n; ++j)
      W(i,j)=polar<real_type>(x,-(2*real_type(PI)*i*j)/l);
#endif
}

template<class G,class V>
void fft_twiddles2(Matrix<G> &W, const V &x,int m, int n, int l, int N)
{
  typedef PROMOTE2(complex<float>,typename Vector<G>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;

  typename DataMatrix<complex_type>::self W2(W.nrows(),N*W.ncols(),(complex_type *)&W(0,0));
  fft_twiddles(W2,x,m,N*n,N*l);
}


template<int p,class V,class G1,class G2,class F>
inline void complex_mixed_fft_step1(const V &x, const Vector<G1> &X, Vector<G2> &Y,int m,int m2,const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef V real_type;

  int n=X.size()/m;
  int q=m/p;

  static __thread struct tls_twiddle { int m,n; real_type x; void *z; } twiddles[3*FFT_TWIDDLES]={ {0,0,0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.m==m && w.n==n && w.x==x) goto found; for (int i=1; i<3*FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.m==m && w.n==n && w.x==x) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.m=m; w.n=n; w.x=x; w.z=aligned_malloc(q*p*sizeof(value_type)); twiddles[0]=w; typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z); fft_twiddles(W,x,q,p,m); }
  found: typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z);

  int l=m/p;
  q=n*l;
#if defined(NDEBUG)
  complex_fft_function<p> fft;
  for (int i=0; i<n; ++i)
  {
    int mi=m*i, li=l*i;
    for (int j=0; j<l; ++j)
    {
      typename StrideArray<Vector<G2> >::self Ya=stride(Y,l,mi+j);
      typename StrideArray<Vector<G2> >::self Yb=stride(Y,l,q-mi-j);
      fft(stride(X,q,li+j),row(W,j),Ya,Yb,f);
    }
  }
#else
  complex_fft_function<p> fft;
  for (int i=0; i<n; ++i)
  {
    int mi=m*i, li=l*i;
    typename StrideArray<Vector<G2> >::self Ya=stride(Y,l,mi);
    fft(stride(X,q,li),row(W,0),Ya,f);
    for (int j=1; j<l; ++j)
    {
      typename StrideArray<Vector<G2> >::self Ya=stride(Y,l,mi+j);
      typename StrideArray<Vector<G2> >::self Yb=stride(Y,l,q-mi-j);
      fft(stride(X,q,li+j),row(W,j),Ya,Yb,f);
    }
  }
#endif
}

template<int p,class G1,class G2,class F>
void complex_mixed_fft_step0(const Vector<G1> &X, Vector<G2> &Y,int m,int m2,const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  int n = X.size()/p;
  int q = n;

  complex_fft_function<p> fft;
  for (int i=0; i<n; ++i)
  {
    typename SubArray<Vector<G2> >::self Ya=sub(Y,p*i,p);
    fft(stride(X,q,i),Ya,f);
  }
}


template<int N,int p,class V,class G1,class G2,class F>
void complex_mixed_fft_step2(const V &x,const Vector<G1> &X, Vector<G2> &Y,int m, int m2,const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef V real_type;

  int n=X.size()/m;
  int q=m/p;

  static __thread struct tls_twiddle { int m,n; real_type x; void *z; } twiddles[FFT_TWIDDLES]={ {0,0,0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.m==m && w.n==n && w.x==x) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.m==m && w.n==n && w.x==x) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.m=m; w.n=n; w.x=x; w.z=aligned_malloc(q*p*sizeof(value_type)); twiddles[0]=w; typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z); fft_twiddles2(W,x,q,p,m,N); }
  found: typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z);

#if defined(NDEBUG)
  complex_crossed_fft_function<p,N> fft;
  for (int i=0; i<q; ++i)
  {
    typename StrideArray<Vector<G2> >::self Y0 = stride(Y,q,i);
    typename StrideArray<Vector<G2> >::self Y1 = stride(Y,q,q-i);
    fft(stride(X,q,i),row(W,i),Y0,Y1,f);
  }
#else
  complex_crossed_fft_function<p,N> fft;
  typename StrideArray<Vector<G2> >::self Y0 = stride(Y,q,0);
  fft(stride(X,q,0),row(W,0),Y0,f);
  for (int i=1; i<q; ++i)
  {
    typename StrideArray<Vector<G2> >::self Y0 = stride(Y,q,i);
    typename StrideArray<Vector<G2> >::self Y1 = stride(Y,q,q-i);
    fft(stride(X,q,i),row(W,i),Y0,Y1,f);
  }
#endif
}


template<int N,int p,class V,class G1,class G2,class F>
void complex_mixed_fft_dense_step2(const V &x, const Vector<G1> &X, Vector<G2> &Y,int m,int m2,const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef V real_type;

  int n=X.size()/m;
  int q=m/p;

  static __thread struct tls_twiddle { int m,n; real_type x; void *z; } twiddles[FFT_TWIDDLES]={ {0,0,0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.m==m && w.n==n && w.x==x) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.m==m && w.n==n && w.x==x) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.m=m; w.n=n; w.x=x; w.z=aligned_malloc(q*p*sizeof(value_type)); twiddles[0]=w; typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z); fft_twiddles2(W,x,q,p,m,N); }
  found: typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z);

#if defined(NDEBUG)
  complex_crossed_fft_function<p,N> fft;
  for (int i=0; i<q; i+=2)
  {
    typename StrideArray<Vector<G2> >::self Ya0 = stride(Y,q,i);
    typename StrideArray<Vector<G2> >::self Ya1 = stride(Y,q,q-i);
    typename StrideArray<Vector<G2> >::self Yb0 = stride(Y,q,i+1);
    typename StrideArray<Vector<G2> >::self Yb1 = stride(Y,q,q-i-1);
    fft(stride(X,q,i),row(W,i),Ya0,Ya1,stride(X,q,i+1),row(W,i+1),Yb0,Yb1,f);
  }
#else
  complex_crossed_fft_function<p,N> fft;
  typename StrideArray<Vector<G2> >::self Y0 = stride(Y,q,0);
  fft(stride(X,q,0),row(W,0),Y0,f);
  for (int i=1; i<q; ++i)
  {
    typename StrideArray<Vector<G2> >::self Y0 = stride(Y,q,i);
    typename StrideArray<Vector<G2> >::self Y1 = stride(Y,q,q-i);
    fft(stride(X,q,i),row(W,i),Y0,Y1,f);
  }
#endif
}



struct complex_mixed_fft_step0_function
{
  template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p,int m,int m2,const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    typedef typename DenseVector<value_type>::self array_type;
    switch (p)
    {
      case  1: break;
#if FFT_LEVEL>=2
      case  2: complex_mixed_fft_step0< 2>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: complex_mixed_fft_step0< 4>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: complex_mixed_fft_step0< 8>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=16
      case 16: complex_mixed_fft_step0<16>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=32
      case 32: complex_mixed_fft_step0<32>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: complex_mixed_fft_step0< 3>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: complex_mixed_fft_step0< 5>(X,Y,m,m2,f); break;
#endif
      default:
        cerr << "fft error: factor unknown" << endl;
    }
  }

  template<class G1,class G2> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p,int m, int m2)
  {
    (*this)(X,Y,p,m,m2,make_fft_adaptor());
  }
};

struct complex_mixed_fft_step1_function
{
  template<class V,class G1,class G2,class F> inline void operator()(const V &x, const Vector<G1> &X, Vector<G2> &Y, int p, int m,int m2, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    typedef typename DenseVector<value_type>::self array_type;
    switch (p)
    {
      case  1: break;
#if FFT_LEVEL>=2
      case  2: complex_mixed_fft_step1< 2>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: complex_mixed_fft_step1< 4>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: complex_mixed_fft_step1< 8>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=16
      case 16: complex_mixed_fft_step1<16>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=32
      case 32: complex_mixed_fft_step1<32>(x,X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: complex_mixed_fft_step1< 3>(x,X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: complex_mixed_fft_step1< 5>(x,X,Y,m,m2,f); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }
  template<class V,class G1,class G2> inline void operator()(const V &x, const Vector<G1> &X, Vector<G2> &Y, int p, int m,int m2)
  {
    (*this)(x,X,Y,p,m,m2,make_fft_adaptor());
  }
};


template<int N>
struct complex_mixed_fft_step2_function
{
  template<class V,class G1,class G2,class F> inline void operator()(const V &x,const Vector<G1> &X, Vector<G2> &Y, int p,int m,int m2,const F &f)
  {
    switch (p)
    {
      case  1: complex_mixed_fft_step2<N, 1>(x,X,Y,m,m2,f); break;
#if FFT_LEVEL>=2
      case  2: complex_mixed_fft_step2<N, 2>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: complex_mixed_fft_step2<N, 4>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: complex_mixed_fft_step2<N, 8>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=16
      case 16: complex_mixed_fft_step2<N,16>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=32
      case 32: complex_mixed_fft_step2<N,32>(x,X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: complex_mixed_fft_step2<N, 3>(x,X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: complex_mixed_fft_step2<N, 5>(x,X,Y,m,m2,f); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }

#ifdef NDEBUG
  template<class V,class G1,class V2,class F> inline void operator()(const V &x,const Vector<G1> &X, Vector<data_vector_generator<V2> > &Y, int p,int m,int m2,const F &f)
  {
    switch (p)
    {
      case  1: complex_mixed_fft_dense_step2<N, 1>(x,X,Y,m,m2,f); break;
#if FFT_LEVEL>=2
      case  2: complex_mixed_fft_dense_step2<N, 2>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: complex_mixed_fft_dense_step2<N, 4>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: complex_mixed_fft_dense_step2<N, 8>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=16
      case 16: complex_mixed_fft_dense_step2<N,16>(x,X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=32
      case 32: complex_mixed_fft_dense_step2<N,32>(x,X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: complex_mixed_fft_dense_step2<N, 3>(x,X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: complex_mixed_fft_dense_step2<N, 5>(x,X,Y,m,m2,f); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }
#endif
};
template<> struct complex_mixed_fft_step2_function<0> { template<class V,class G1,class G2,class F> inline void operator()(const V &x,const Vector<G1> &X, Vector<G2> &Y, int p,int m,int m2,const F &f) { complex_mixed_fft_step1_function()(x,X,Y,p,m,m2,f); } };
template<> struct complex_mixed_fft_step2_function<1> { template<class V,class G1,class G2,class F> inline void operator()(const V &x,const Vector<G1> &X, Vector<G2> &Y, int p,int m,int m2,const F &f) { complex_mixed_fft_step1_function()(x,X,Y,p,m,m2,f); } };


inline void complex_fft_factors(int n, int *X)
{
#if FFT_LEVEL>=16
  switch (n)
  {
    case   64: X[0]=2; X[1]=16; X[2]= 4; return;
    case  256: X[0]=2; X[1]=16; X[2]=16; return;
    case 2048: X[0]=3; X[1]=16; X[2]=16; X[3]=8; return;
  };
#endif

  int k=0;

#if FFT_LEVEL>=32
  for (; n%32==0; n/=32) X[++k]=32;
#endif
#if FFT_LEVEL>=16
  for (; n%16==0; n/=16) X[++k]=16;
#endif
#if FFT_LEVEL>=8
  for (; n%8==0; n/=8) X[++k]=8;
#endif
#if FFT_LEVEL>=4
  for (; n%4==0; n/=4) X[++k]=4;
#endif
#if FFT_LEVEL>=2
  for (; n%2==0; n/=2) X[++k]=2;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
  for (; n%5==0; n/=5) X[++k]=5;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
  for (; n%3==0; n/=3) X[++k]=3;
#endif
  if (n!=1) { cerr << "fft factor unknown" << endl; }
  X[0]=k;
}

template<int N,class V,class G1,class G2,class G3,class F>
void complex_mixed_fft(const V &x,const Vector<G1> &X, Vector<G2> &T, Vector<G2> &U, Vector<G3> &Y, const F &f)
{
  int n=X.size();
  static __thread struct { int n; int factors[50]; } decomp={0};
  if (n==decomp.n); else { decomp.n=n; complex_fft_factors(n, decomp.factors); }

  int o=decomp.factors[0];

  int p=decomp.factors[1];
  int m=p, m2=p;

  if (o%2==1) swap(U,T);
  complex_mixed_fft_step0_function()(X,T,p,m,m2);
  for (int k=2; k<o; ++k)
  {
    p=decomp.factors[k];
    complex_mixed_fft_step1_function()(V(1),T,U,p,m*=p,m2*=p);
    swap(T,U);
  }
  p=decomp.factors[o];
  complex_mixed_fft_step2_function<N>()(x,T,Y,p,m*p,m2*p,f);
}


static __thread struct { int n; void *p; } fft_buf0 = {0,NULL};
static __thread struct { int n; void *p; } fft_buf1 = {0,NULL};
static __thread struct { int n; void *p; } fft_buf2 = {0,NULL};

template<int N>
struct complex_mixed_fft_function
{
  template<class V,class G1,class V2,class G2,class F>
  inline void operator()(const V &x,const Vector<G1> &X, Vector<data_vector_generator<V2> > &T, Vector<data_vector_generator<V2> > &U, Vector<G2> &Y, const F &f)
  {
    typedef V2 value_type;
    typedef typename TinyVector<N,value_type>::self elem_type;
    typename DenseVector<value_type>::self Z=X;
    typename DataVector<elem_type>::self X2(Z.size()/N, (elem_type *)&*Z.begin());
    typename DataVector<elem_type>::self T2(T.size()/N, (elem_type *)&*T.begin());
    typename DataVector<elem_type>::self U2(U.size()/N, (elem_type *)&*U.begin());
    complex_mixed_fft<N>(x,X2,T2,U2,Y,f);
  }

  template<class V,class V2,class G2,class F>
  inline void operator()(const V &x,const Vector<data_vector_generator<V2> > &X, Vector<data_vector_generator<V2> > &T, Vector<data_vector_generator<V2> > &U, Vector<G2> &Y, const F &f)
  {
    typedef V2 value_type;
    typedef typename TinyVector<2,value_type>::self elem_type;
    typename DataVector<elem_type>::self X2(X.size()/N, (elem_type *)&*X.begin());
    typename DataVector<elem_type>::self T2(T.size()/N, (elem_type *)&*T.begin());
    typename DataVector<elem_type>::self U2(U.size()/N, (elem_type *)&*U.begin());
    complex_mixed_fft<N>(x,X2,T2,U2,Y,f);
  }


  template<class V,class G1,class G2,class F>
  void operator()(const V &x,const Vector<G1> &X, Vector<G2> &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
    (*this)(x,X,T2,U2,Y,f);
  }

  template<class V,class V2,class F>
  void operator()(const V &x,const Vector<data_vector_generator<V2> > &X, Vector<data_vector_generator<V2> > &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,V2) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    typename DataVector<value_type>::self U2(Y.size(), (value_type *)&*Y.begin());
    typename DataVector<value_type>::self Y2(Y.size(), (value_type *)&*Y.begin());
    (*this)(x,X,T2,U2,Y2,f);
  }
};

template<>
struct complex_mixed_fft_function<0>
{
  template<class V,class G1,class V2,class G2,class F>
  inline void operator()(const V &x,const Vector<G1> &X, Vector<data_vector_generator<V2> > &T, Vector<data_vector_generator<V2> > &U, Vector<G2> &Y, const F &f)
  {
    complex_mixed_fft<0>(x,X,T,U,Y,f);
  }

  template<class V,class G1,class G2,class F>
  void operator()(const V &x,const Vector<G1> &X, Vector<G2> &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
    (*this)(x,X,T2,U2,Y,f);
  }

  template<class V,class V2,class F>
  void operator()(const V &x,const Vector<data_vector_generator<V2> > &X, Vector<data_vector_generator<V2> > &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,V2) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    typename DataVector<value_type>::self U2(X.size(), (value_type *)&*Y.begin());
    (*this)(x,X,T2,U2,Y,f);
  }
};

template<>
struct complex_mixed_fft_function<1>
{
  template<class V,class G1,class V2,class G2,class F>
  inline void operator()(const V &x,const Vector<G1> &X, Vector<data_vector_generator<V2> > &T, Vector<data_vector_generator<V2> > &U, Vector<G2> &Y, const F &f)
  {
    typedef V2 value_type;
    typedef typename TinyVector<1,value_type>::self elem_type;
    typename tinyBlockVector<const Vector<G1>,1>::self X2 = block<1>(X);
    typename DataVector<elem_type>::self T2(T.size()/1, (elem_type *)&*T.begin());
    typename DataVector<elem_type>::self U2(U.size()/1, (elem_type *)&*U.begin());
    complex_mixed_fft<1>(x,X2,T2,U2,Y,f);
  }

  template<class V,class V2,class G2,class F>
  inline void operator()(const V &x,const Vector<data_vector_generator<V2> > &X, Vector<data_vector_generator<V2> > &T, Vector<data_vector_generator<V2> > &U, Vector<G2> &Y, const F &f)
  {
    typedef V2 value_type;
    typedef typename TinyVector<1,value_type>::self elem_type;
    typename DataVector<elem_type>::self X2(X.size()/1, (elem_type *)&*X.begin());
    typename DataVector<elem_type>::self T2(T.size()/1, (elem_type *)&*T.begin());
    typename DataVector<elem_type>::self U2(U.size()/1, (elem_type *)&*U.begin());
    complex_mixed_fft<1>(x,X2,T2,U2,Y,f);
  }


  template<class V,class G1,class G2,class F>
  void operator()(const V &x,const Vector<G1> &X, Vector<G2> &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
    (*this)(x,X,T2,U2,Y,f);
  }

  template<class V,class V2,class F>
  void operator()(const V &x,const Vector<data_vector_generator<V2> > &X, Vector<data_vector_generator<V2> > &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,V2) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    typename DataVector<value_type>::self U2(Y.size(), (value_type *)&*Y.begin());
    typename DataVector<value_type>::self Y2(Y.size(), (value_type *)&*Y.begin());
    (*this)(x,X,T2,U2,Y2,f);
  }
};

template<      class G1,class G2,class F> void complex_mixed_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;
  complex_mixed_fft_function<0>()(1/real_type(f.normalization(X)),X,Y,f);
}
template<int N,class G1,class G2,class F> void complex_mixed_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;
  complex_mixed_fft_function<N>()(1/real_type(f.normalization(X)),X,Y,f);
}


template<class G1,class G2,class F,class V>
void aux_complex_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, const V &)
{
  switch (X.size())
  {
    case  1: complex_fft_function< 1>()(X,Y,f); break;
#if FFT_LEVEL>=2
    case  2: complex_fft_function< 2>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=4
    case  4: complex_fft_function< 4>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=8
    case  8: complex_fft_function< 8>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=16
    case 16: complex_fft_function<16>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=32
    case 32: complex_fft_function<32>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=64
    case 64: complex_fft_function<64>()(X,Y,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
    case  3: complex_fft_function< 3>()(X,Y,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
    case  5: complex_fft_function< 5>()(X,Y,f); break;
#endif
    default: complex_mixed_fft(X,Y,f);
  }
}

#ifdef SSE
template<class G1,class G2,class F>
void aux_complex_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, complex<float>)
{
  int n = X.size();
  switch (n)
  {
    case  1: complex_fft_function<1>()(X,Y,f); break;
#if FFT_LEVEL>=2
    case  2: complex_fft_function<2>()(X,Y,f); break;
    case  4: complex_crossed_fft_function< 2,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=4
    case  8: complex_crossed_fft_function< 4,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=8
    case 16: complex_crossed_fft_function< 8,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=16
    case 32: complex_crossed_fft_function<16,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=32
    case 64: complex_crossed_fft_function<32,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=64 && (!defined(__INTEL_COMPILER) || (__INTEL_COMPILER<900))
    case 128: complex_crossed_fft_function<64,2>()(block<2>(X),Y,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
    case  3: complex_fft_function<3>()(X,Y,f); break;
    case  6: complex_crossed_fft_function< 3,2>()(block<2>(X),Y,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
    case  5: complex_fft_function<5>()(X,Y,f); break;
    case 10: complex_crossed_fft_function< 5,2>()(block<2>(X),Y,f); break;
#endif
    default:
    {
      if (n%4==0) complex_mixed_fft<2>(X,Y,f);
#if !defined(FFT_ONLY_POW2)
      else  complex_mixed_fft(X,Y,f);
#else
      else cerr << "fft error: factor unknown" << endl;
#endif
    }
  }
}
#endif

#ifdef SSE2
template<class G1,class G2,class F>
void aux_complex_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, complex<double>)
{
  int n = X.size();
  switch (n)
  {
    case  1: complex_fft_function< 1>()(block<1>(X),Y,f); break;
#if FFT_LEVEL>=2
    case  2: complex_fft_function< 2>()(block<1>(X),Y,f); break;
#endif
#if FFT_LEVEL>=4
    case  4: complex_fft_function< 4>()(block<1>(X),Y,f); break;
#endif
#if FFT_LEVEL>=8
    case  8: complex_fft_function< 8>()(block<1>(X),Y,f); break;
#endif
#if FFT_LEVEL>=16
    case 16: complex_fft_function<16>()(block<1>(X),Y,f); break;
#endif
#if FFT_LEVEL>=32
    case 32: complex_fft_function<32>()(block<1>(X),Y,f); break;
#endif
#if FFT_LEVEL>=64
    case 64: complex_fft_function<64>()(block<1>(X),Y,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
    case  3: complex_fft_function< 3>()(block<1>(X),Y,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
    case  5: complex_fft_function< 5>()(block<1>(X),Y,f); break;
#endif
    default:
      complex_mixed_fft<1>(X,Y,f);
  }
}
#endif

//template<class G1,class G2,class F>
//void shift_dense_complex_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f)
//{
//  typedef Vector<G1> array_type;
//  typedef typename array_type::value_type value_type;
//  typedef typename value_type::value_type real_type;
//
//  int i0 = X.lower_bound();
//  int n = X.size();
//
//  Y.resize(n); Y.set_lower_bound(i0);
//
//  if (i0==0)
//  {
//    typename DenseVector<value_type>::self U;
//    typename DenseVector<value_type>::self V;
//
//    U.swap(const_cast<array_type &>(X));
//    Y.swap(V);
//    complex_fft(U,V,f);
//    Y.swap(V);
//    U.swap(const_cast<array_type &>(X));
//  }
//  else
//  {
//    static array_type W(0);
//    if ( W.size()!=n || W.lower_bound()!=i0)
//    {
//      W.resize(n); W.set_lower_bound(i0);
//      value_type w0=polar(real_type(1),(-2*real_type(PI)*i0)/n);
//      value_type w=w0;
//      W[i0]=value_type(1);
//      W[i0+1]=w;
//      for (int i=i0+2, imax=W.upper_bound()+1; i<imax; ++i)
//        W[i]=(w*=w0);
//    }
//    Y = W*X;
//
//    typename DenseVector<value_type>::self V;
//
//    Y.swap(V);
//    complex_fft_in_place(V,f);
//    Y.swap(V);
//  }
//}
//
////{secret}
//template<class G1,int C1,class G2,int C2,class F> void complex_fft(const Vector<shift_array_generator<Vector<G1>,C1> > &X, Vector<shift_array_generator<Vector<G2>,C2> > &Y, const F &f) { return shift_dense_complex_fft(X,Y,f); }
////{secret}
//template<class G1,int C1,class G2,int C2,class F> void complex_fft(const Vector<shift_dense_array_generator<Vector<G1>,C1> > &X, Vector<shift_dense_array_generator<Vector<G2>,C2> > &Y, const F &f) { return shift_dense_complex_fft(X,Y,f); }

template<class G1,class G2,class F> inline void aux_complex_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f) { aux_complex_fft(X,Y,f,typename Vector<G1>::value_type()); }

//{secret}
template<class G1,class G2,class F>
inline void complex_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f)
{
  //int i0 = X.lower_bound();
  //if (i0!=0)
  //  cerr << "error: fft of shifted vectors not impemented" << endl;
  //else
  {
    typename ArrayData<Vector<G2> >::self Y2 = data(Y);
    aux_complex_fft(data(X),Y2,f);
  }
}


template<class G1,class G2> inline void complex_fft(const Vector<G1> &X, Vector<G2> &Y) { complex_fft(X,Y,make_fft_adaptor()); }
template<class G> inline PROMOTE2(complex<float>,Vector<G>) complex_fft(const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) Y(X.size()); complex_fft(X,Y); return Y; }

template<class G1,class G2> inline void complex_norm_fft(const Vector<G1> &X, Vector<G2> &Y) { complex_fft(X,Y,make_norm_fft_adaptor()); }
template<class G1,class G2> inline void complex_abs_fft (const Vector<G1> &X, Vector<G2> &Y) { complex_fft(X,Y,make_abs_fft_adaptor ()); }
template<class G> inline typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self complex_norm_fft(const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size()); complex_norm_fft(X,Y); return Y; }
template<class G> inline typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self complex_abs_fft (const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size()); complex_abs_fft (X,Y); return Y; }

template<class G1,class G2> inline void complex_ifft    (const Vector<G1> &X, Vector<G2> &Y) { complex_fft(X,Y,make_ifft_adaptor(Y)); }
template<class G> inline PROMOTE2(complex<float>,Vector<G>) complex_ifft(const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) Y(X.size()); complex_ifft(X,Y); return Y; }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class G,class S>
inline void real_fft_1(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==1);

  Y.first0(X[0]);
}

template<class G,class S>
inline void real_fft_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==2);

  value_type t0,t1;

  t0=X[0]; t1=X[1];

  Y.first0(  t0+t1,value_type(0));
  Y.second(1,t0-t1,value_type(0));
}

template<class G,class S>
inline void real_fft_3(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==3);

  value_type z0=X[0],z1=X[1],z2=X[2];

  value_type t1;
  t1=z1+z2;
	Y.first0(  z0+t1,value_type(0));
	Y.first (1,z0-real_type(0.5)*t1,SIN01_03*(z2-z1));
}

template<class G,class S>
inline void real_fft_4(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==4);

  value_type t0,t1,t2,t3;

  { value_type v0,v1; v0=X[0]; v1=X[2]; t0=v0+v1; t1=v0-v1; v0=X[3]; v1=X[1]; t2=v0+v1; t3=v0-v1; }
  Y.first0(  t0+t2,value_type(0));
  Y.second(2,t0-t2,value_type(0));
  Y.first (1,t1,t3);
}

template<class G,class S>
inline void real_fft_5(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==5);

  value_type t0, t1, t2, t3, t4, t5;

	t0 = X[0];
	{
	  value_type v0, v1;
	  v0=X[1]; v1=X[4]; t1=v0+v1; t2=v1-v0;
	  v0=X[2]; v1=X[3]; t3=v0+v1; t4=v0-v1;
	}

	t5 = t1+t3;
	Y.first0(  t0+t5,value_type(0));

	{
	  value_type v0, v1;
	  v0=real_type(SQRT_5_16)*(t1-t3);
	  v1=t0-real_type(0.25)*t5;
	  Y.first(1,v1+v0,SIN02_05*t2-SIN01_05*t4);
	  Y.first(2,v1-v0,SIN02_05*t4+SIN01_05*t2);
	}
}

template<class G,class S>
inline void real_fft_8(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==8);

  value_type t0,t1,t2,t3,t4,t5,t6,t7;
  { value_type v0,v1,v2,v3; v0=X[0]; v1=X[4]; v2=v0+v1; t4=v0-v1; v0=X[2]; v1=X[6]; v3=v0+v1; t7=v0-v1; t2=v2+v3; t0=v2-v3; }
	{ value_type v0,v1,v2,v3,v4,v5; v0=X[7]; v1=X[3]; v2=v0+v1; v4=v0-v1; v0=X[1]; v1=X[5]; v3=v0+v1; v5=v0-v1; t5=SIN16*(v4+v5); t6=SIN16*(v4-v5); t3=v2+v3; t1=v2-v3; }

	Y.first(2,t0,t1);
	Y.first0(t2+t3,value_type(0)); Y.second(4,t2-t3,value_type(0));
	Y.first(1,t4+t5,t6-t7); Y.first(3,t4-t5,t6+t7);
}

template<class G,class S>
inline void real_fft_16(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==16);

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t15,t14;

  { value_type v0,v1; v0=X[ 0]; v1=X[ 8]; t00=v0+v1; t01=v0-v1; v0=X[ 4]; v1=X[12]; t02=v0-v1; v0+=v1; t03=t00+v0 ; t00-=v0; }
  { value_type v0,v1; v0=X[ 2]; v1=X[10]; t06=v0+v1; t07=v0-v1; v0=X[14]; v1=X[ 6]; t04=v0+v1; v0-=v1; t05=t04+t06; t04=t04-t06; t06=SIN16*(v0+t07); t07=SIN16*(v0-t07); }
  {
    value_type v6,v7;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[15];v5=X[ 7]; v0=v4+v5; v1=v4-v5; v4=X[3]; v5=X[11]; v2=v4+v5; v3=v4-v5; } t10=v0+v2; v6=v0-v2; t08=SIN08 *v1-SIN24*v3; t09=SIN24*v1+SIN08 *v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[1]; v5=X[ 9]; v0=v4+v5; v1=v4-v5; v4=X[5]; v5=X[13]; v2=v4+v5; v3=v4-v5; } t14=v0+v2; v7=v0-v2; t12=SIN08 *v1+SIN24*v3; t13=SIN24*v1-SIN08 *v3; }
    t15=SIN16; t11=t15*(v6+v7); t15=t15*(v6-v7);
  }

  Y.first(4,t03-t05,t10-t14);
  { Y.first(2,t00+t11,t15+t04); Y.first (6,t00-t11,t15-t04); t11=t03+t05; t15=t10+t14; Y.first0(t11+t15,value_type(0)); Y.second(8,t11-t15,value_type(0)); }
  { value_type v0,v1,v2,v3; v0=t01+t06; v1=t13+t09; v2=t07-t02; v3=t08-t12; Y.first(1,v0+v1,v2+v3); Y.first(7,v0-v1,v3-v2); }
  { value_type v0,v1,v2,v3; v0=t01-t06; v1=t12+t08; v2=t02+t07; v3=t09-t13; Y.first(3,v0+v1,v2+v3); Y.first(5,v0-v1,v3-v2); }
}

template<class G,class S>
inline void real_fft_32(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==32);

  value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;

  { value_type v0,   v2;          { value_type v4,v5; v4=X[ 0]; v5=X[16]; v0=v4+v5; t02=v4-v5; v4=X[ 8]; v5=X[24]; v2=v4+v5; t03=v4-v5; } t00=v0+v2; t01=v0-v2; }
  { value_type v0,v1,v2,v3;       { value_type v4,v5; v4=X[ 4]; v5=X[20]; v0=v4+v5; v1=v4-v5;  v4=X[28]; v5=X[12]; v2=v4+v5; v3=v4-v5; } t04=v2+v0; t05=v2-v0; t06=SIN16*(v3+v1);     t07=SIN16*(v3-v1); }
  { value_type v0,v1,v2,v3;       { value_type v4,v5; v4=X[30]; v5=X[14]; v0=v4+v5; v1=v4-v5;  v4=X[ 6]; v5=X[22]; v2=v4+v5; v3=v4-v5; } t08=v0+v2; t09=v0-v2; t10=SIN24*v1+SIN08*v3; t11=SIN08*v1-SIN24*v3; }
  { value_type v0,v1,v2,v3;       { value_type v4,v5; v4=X[ 2]; v5=X[18]; v0=v4+v5; v1=v4-v5;  v4=X[10]; v5=X[26]; v2=v4+v5; v3=v4-v5; } t12=v0+v2; t13=v0-v2; t14=SIN24*v1-SIN08*v3; t15=SIN08*v1+SIN24*v3; }
  { value_type v0,v1,v2,v3,v6,v7; { value_type v4,v5; v4=X[31]; v5=X[15]; v0=v4-v5; v1=v4+v5;  v4=X[ 7]; v5=X[23]; v2=v4-v5; v3=v4+v5; } t16=v1+v3; t19=v1-v3; { value_type v4,v5; v4=X[ 3]; v5=X[19]; v6=v4-v5; v7=v4+v5; v4=X[27]; v5=X[11]; v1=v4-v5; v3=v4+v5; } t17=v3+v7; t18=v3-v7; v3=SIN16*(v1-v6); v7=SIN16*(v1+v6); t20=v3-v2; t21=v3+v2;  t22=v0+v7; t23=v0-v7; }
  { value_type v0,v1,v2,v3,v6,v7; { value_type v4,v5; v4=X[ 1]; v5=X[17]; v0=v4-v5; v1=v4+v5;  v4=X[ 9]; v5=X[25]; v2=v4-v5; v3=v4+v5; } t24=v1+v3; t27=v1-v3; { value_type v4,v5; v4=X[ 5]; v5=X[21]; v6=v4-v5; v7=v4+v5; v4=X[29]; v5=X[13]; v1=v4-v5; v3=v4+v5; } t25=v3+v7; t26=v3-v7; v3=SIN16*(v1-v6); v7=SIN16*(v1+v6); t28=v3-v2; t29=v3+v2; t30=v0+v7; t31=v0-v7; }

  { value_type v0,v1,v2,v3,v4,v5;       v0=t00+t04; v1=t12+t08; v2=v0+v1; v3=t24+t25; v4=t16+t17; v5=v3+v4; Y.first0(v2+v5,value_type(0)); Y.second(16,v2-v5,value_type(0)); Y.first(8,v0-v1,v4-v3); }
  { value_type v0,v1,v2,v3;             v0=t24-t25; v1=t16-t17; v2=SIN16*(v1+v0); v3=SIN16*(v1-v0); v0=t00-t04; v1=t08-t12; Y.first(4 ,v0+v2,v1+v3); Y.first(12,v0-v2,v3-v1); }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; v1=SIN24*t27+SIN08*t26; v5=SIN24*t19-SIN08*t18; v4=v5+v1; v5=v5-v1; v1=SIN24*t26-SIN08*t27; v7=SIN08*t19+SIN24*t18; v6=v7+v1; v7=v7-v1; v1=SIN16*(t09+t13); v0=t01+v1; v1=t01-v1; v3=SIN16*(t09-t13); v2=v3-t05; v3=v3+t05; Y.first( 2,v0+v4,v6+v3); Y.first(14,v0-v4,v6-v3); Y.first( 6,v1+v7,v5+v2); Y.first(10,v1-v7,v5-v2); }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type v8,v9; v8=t02+t06; v9=t14+t10; v0=v8+v9; v1=v8-v9; v8=SIN28*t28-SIN04*t30; v9=SIN04*t22+SIN28*t20; v2=v9+v8; v3=v9-v8; } { value_type v8,v9; v8=SIN28*t30+SIN04*t28; v9=SIN28*t22-SIN04*t20; v4=v9+v8; v5=v9-v8; v8=t11-t15; v9=t07-t03; v6=v8-v9; v7=v8+v9; } Y.first( 1,v0+v4,v2+v7); Y.first(15,v0-v4,v2-v7); Y.first( 7,v1+v3,v5+v6); Y.first( 9,v1-v3,v5-v6); }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type v8,v9; v8=t02-t06; v9=t15+t11; v0=v8+v9; v1=v8-v9; v8=SIN20*t29-SIN12*t31; v9=SIN12*t23+SIN20*t21; v2=v9+v8; v3=v9-v8; } { value_type v8,v9; v8=SIN20*t31+SIN12*t29; v9=SIN20*t23-SIN12*t21; v4=v9+v8; v5=v9-v8; v8=t10-t14; v9=t03+t07; v6=v8-v9; v7=v8+v9; } Y.first( 3,v0+v4,v2+v7); Y.first(13,v0-v4,v2-v7); Y.first( 5,v1+v3,v5+v6); Y.first(11,v1-v3,v5-v6); }
}

template<class G,class S>
inline void real_fft_64(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type) value_type;
  typedef VALUE_TYPE(value_type) real_type;
  assert(X.size()==64);

	value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65;

	{ value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 0]; u1=X[32]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[56]; u1=X[24]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[16]; u1=X[48]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[ 8]; u1=X[40]; v6=u0+u1; v7=u0-u1; } t00=v0-v4; t01=v2-v6; t06=v0+v4; t07=v6+v2; { value_type u0; u0=SIN16*(v3+v7); t02=v1+u0; t03=v1-u0; } { value_type u0; u0=SIN16*(v3-v7); t04=u0-v5; t05=v5+u0; } t08=t06+t07; }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[57]; u1=X[25]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[ 1]; u1=X[33]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[ 9]; u1=X[41]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[17]; u1=X[49]; v6=u0-u1; v7=u0+u1; } t09=v3-v7; t10=v1-v5; { value_type u0,u1; u0=v3+v7; u1=v5+v1; t15=u0+u1; t16=u0-u1; } { value_type u0; u0=SIN16*(v0-v4); t11=u0-v6; t12=u0+v6; } { value_type u0; u0=SIN16*(v4+v0); t13=v2+u0; t14=v2-u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[63]; u1=X[31]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[15]; u1=X[47]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[ 7]; u1=X[39]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[55]; u1=X[23]; v6=u0-u1; v7=u0+u1; } t17=v1-v3; t18=v7-v5; { value_type u0,u1; u0=v1+v3; u1=v5+v7; t23=u0+u1; t24=u0-u1; } { value_type u0; u0=SIN16*(v4+v6); t19=v0+u0; t20=v0-u0; } { value_type u0; u0=SIN16*(v6-v4); t21=u0-v2; t22=u0+v2; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[62]; u1=X[30]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[54]; u1=X[22]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[14]; u1=X[46]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[ 6]; u1=X[38]; v6=u0+u1; v7=u0-u1; } { value_type u0; u0=SIN16*(v7+v3); t25=v1+u0; t26=v1-u0; } { value_type u0; u0=SIN16*(v3-v7); t27=u0-v5; t28=v5+u0; } { value_type u0,u1; u0=v0+v4; u1=v6+v2; t29=u0+u1; t30=u0-u1; } { value_type u0,u1; u0=v0-v4; u1=v2-v6; t31=SIN24*u0-SIN08*u1; t32=SIN24*u1+SIN08*u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 2]; u1=X[34]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[58]; u1=X[26]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[18]; u1=X[50]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[10]; u1=X[42]; v6=u0+u1; v7=u0-u1; } { value_type u0; u0=SIN16*(v3-v7); t33=u0-v5; t34=v5+u0; } { value_type u0; u0=SIN16*(v7+v3); t35=v1+u0; t36=v1-u0; } { value_type u0,u1; u0=v0+v4; u1=v6+v2; t37=u0+u1; t38=u0-u1; } { value_type u0,u1; u0=v0-v4; u1=v2-v6; t39=SIN24*u0+SIN08*u1; t40=SIN24*u1-SIN08*u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 4]; u1=X[36]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[12]; u1=X[44]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[20]; u1=X[52]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[60]; u1=X[28]; v6=u0+u1; v7=u0-u1; } t41=v0+v4; t42=v6+v2; t43=t41+t42;                         { value_type u0,u1; u0=v0-v4; u1=v6-v2; t44=SIN16*(u0+u1); t45=SIN16*(u1-u0); } { value_type u0,u1; u0=SIN24*v1-SIN08*v5; u1=SIN24*v7+SIN08*v3; t46=u0+u1; t47=u1-u0; } { value_type u0,u1; u0=SIN08*v7-SIN24*v3; u1=SIN08*v1+SIN24*v5; t48=u0-u1; t49=u1+u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[61]; u1=X[29]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[21]; u1=X[53]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[13]; u1=X[45]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[ 5]; u1=X[37]; v6=u0-u1; v7=u0+u1; } { value_type u0,u1; u0=v7+v3; u1=v1+v5; t52=u0+u1; t53=u1-u0; } { value_type u0,u1; u0=v7-v3; u1=v1-v5; t56=SIN16*(u0+u1); t57=SIN16*(u1-u0); } { value_type u0,u1; u0=SIN08*v0-SIN24*v4; u1=SIN08*v6+SIN24*v2; t50=u0-u1; t51=u1+u0; } { value_type u0,u1; u0=SIN24*v6-SIN08*v2; u1=SIN24*v0+SIN08*v4; t54=u0+u1; t55=u1-u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 3]; u1=X[35]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[11]; u1=X[43]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[19]; u1=X[51]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[59]; u1=X[27]; v6=u0-u1; v7=u0+u1; } { value_type u0,u1; u0=v1+v5; u1=v7+v3; t60=u0+u1; t61=u1-u0; } { value_type u0,u1; u0=v1-v5; u1=v7-v3; t64=SIN16*(u0+u1); t65=SIN16*(u1-u0); } { value_type u0,u1; u0=SIN24*v0-SIN08*v4; u1=SIN24*v6+SIN08*v2; t58=u0+u1; t59=u1-u0; } { value_type u0,u1; u0=SIN08*v6-SIN24*v2; u1=SIN08*v0+SIN24*v4; t62=u0-u1; t63=u1+u0; } }

  { value_type v0,v1,v2,v3,v4,v5; v0=t08+t43; v1=t37+t29; v2=v0+v1; v3=t15+t52; v4=t23+t60; v5=v3+v4; Y.first0(v2+v5,value_type(0)); Y.first(16,v0-v1,v4-v3); Y.second(32,v2-v5,value_type(0)); }
  { value_type v0,v1,v2,v3; v0=t08-t43; v1=t29-t37; { value_type u0,u1; u1=t15-t52; u0=t23-t60; v2=SIN16*(u0+u1); v3=SIN16*(u0-u1); } Y.first( 8,v0+v2,v1+v3); Y.first(24,v0-v2,v3-v1); }
  { value_type v0,v1,v2,v3,v4,v5,v7,v6; { value_type u0,u1; u0=t06-t07; u1=SIN16*(t38+t30); v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=SIN24*t53-SIN08*t16; u1=SIN24*t61+SIN08*t24; v2=u0+u1; v3=u1-u0; } { value_type u0,u1; u0=SIN24*t16+SIN08*t53; u1=SIN24*t24-SIN08*t61; v4=u0+u1; v5=u1-u0; } { value_type u0,u1; u0=SIN16*(t30-t38); u1=t42-t41; v6=u0+u1; v7=u0-u1;  } Y.first( 4,v0+v4,v2+v6); Y.first(28,v0-v4,v2-v6); Y.first(12,v1+v3,v5+v7); Y.first(20,v1-v3,v5-v7); }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t31-t39; u1=t45-t01; v2=u0-u1; v3=u1+u0; } { value_type u0,u1; u0=t00-t44; u1=t32-t40; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=t09-t56; u1=t57-t10; v4=SIN20*u0+SIN12*u1; v5=SIN20*u1-SIN12*u0; } { value_type u0,u1; u0=t17-t64; u1=t65-t18; v6=SIN20*u0-SIN12*u1; v7=SIN20*u1+SIN12*u0; } { value_type u0,u1; u0=v4+v6; u1=v5+v7; Y.first( 6,v0+u0,v3+u1); Y.first(26,v0-u0,u1-v3); } { value_type u0,u1; u0=v7-v5; u1=v6-v4; Y.first(10,v1+u0,v2+u1); Y.first(22,v1-u0,u1-v2); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t00+t44; u1=t39+t31; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=t40+t32; u1=t01+t45; v2=u0-u1; v3=u1+u0; } { value_type u0,u1; u0=t09+t56; u1=t10+t57; v4=SIN28*u0+SIN04*u1; v5=SIN28*u1-SIN04*u0; } { value_type u0,u1; u0=t17+t64; u1=t18+t65; v6=SIN28*u0-SIN04*u1; v7=SIN28*u1+SIN04*u0; } { value_type u0,u1; u0=v4+v6; u1=v5+v7; Y.first( 2,v0+u0,v3+u1); Y.first(30,v0-u0,u1-v3); } { value_type u0,u1; u0=v6-v4; u1=v7-v5; Y.first(14,v1+u1,v2+u0); Y.first(18,v1-u1,u0-v2); } }
  { value_type v00,v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12,v13,v14,v15; v00=t02-t46; v01=t02+t46; v02=t04+t48; v03=t48-t04; { value_type u0,u1; u0=SIN28*t27+SIN04*t25; u1=SIN28*t33-SIN04*t35; v04=u0-u1; v05=u1+u0; } { value_type u0,u1; u0=t19+t58; u1=t21+t62; v06=SIN30*u0-SIN02*u1; v07=SIN30*u1+SIN02*u0; } { value_type u0,u1; u0=SIN28*t25-SIN04*t27; u1=SIN28*t35+SIN04*t33; v08=u0-u1; v09=u1+u0; } { value_type u0,u1; u0=t50-t11; u1=t13-t54; v10=SIN18*u1+SIN14*u0; v11=SIN18*u0-SIN14*u1; } { value_type u0,u1; u0=t11+t50; u1=t13+t54; v12=SIN30*u1+SIN02*u0; v13=SIN30*u0-SIN02*u1; } { value_type u0,u1; u0=t19-t58; u1=t62-t21; v14=SIN18*u0-SIN14*u1; v15=SIN18*u1+SIN14*u0; } { value_type u0,u1,u2,u3; u0=v01+v09; u1=v12+v06; u2=v02+v05; u3=v13+v07; Y.first( 1,u0+u1,u2+u3); Y.first(31,u0-u1,u3-u2); } { value_type u0,u1,u2,u3; u0=v01-v09; u1=v07-v13; u2=v05-v02; u3=v06-v12; Y.first(15,u0+u1,u2+u3); Y.first(17,u0-u1,u3-u2); } { value_type u0,u1,u2,u3; u0=v00+v04; u1=v10+v14; u2=v03+v08; u3=v11+v15; Y.first( 7,u0+u1,u2+u3); Y.first(25,u0-u1,u3-u2); } { value_type u0,u1,u2,u3; u0=v00-v04; u1=v15-v11; u2=v08-v03; u3=v14-v10; Y.first( 9,u0+u1,u2+u3); Y.first(23,u0-u1,u3-u2); } }
  { value_type v00,v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12,v13,v14,v15; v00=t03+t49; v01=t03-t49; v02=t47-t05; v03=t05+t47; { value_type u0,u1; u0=SIN20*t36+SIN12*t34; u1=SIN20*t26-SIN12*t28; v04=u0+u1; v05=u1-u0; } { value_type u0,u1; u0=t20-t63; u1=t59-t22; v06=SIN22*u0-SIN10*u1; v07=SIN22*u1+SIN10*u0; } { value_type u0,u1; u0=SIN20*t34-SIN12*t36; u1=SIN20*t28+SIN12*t26; v08=u0+u1; v09=u1-u0; } { value_type u0,u1; u0=t14+t51; u1=t12+t55; v10=SIN26*u0+SIN06*u1; v11=SIN26*u1-SIN06*u0; } { value_type u0,u1; u0=t14-t51; u1=t55-t12; v12=SIN22*u0+SIN10*u1; v13=SIN22*u1-SIN10*u0; } { value_type u0,u1; u0=t20+t63; u1=t22+t59; v14=SIN26*u0-SIN06*u1; v15=SIN26*u1+SIN06*u0; } { value_type u0,u1,u2,u3; u0=v00+v04; u1=v10+v14; u2=v03+v08; u3=v11+v15; Y.first( 3,u0+u1,u2+u3); Y.first(29,u0-u1,u3-u2); } { value_type u0,u1,u2,u3; u0=v00-v04; u1=v15-v11; u2=v08-v03; u3=v14-v10; Y.first(13,u0+u1,u2+u3); Y.first(19,u0-u1,u3-u2); } { value_type u0,u1,u2,u3; u0=v01+v09; u1=v12+v06; u2=v02+v05; u3=v13+v07; Y.first( 5,u0+u1,u2+u3); Y.first(27,u0-u1,u3-u2); } { value_type u0,u1,u2,u3; u0=v01-v09; u1=v07-v13; u2=v05-v02; u3=v06-v12; Y.first(11,u0+u1,u2+u3); Y.first(21,u0-u1,u3-u2); } }
}


template<int M> struct real_fft_base_function {};
template<> struct real_fft_base_function< 1> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_1 (X,f); } };
template<> struct real_fft_base_function< 2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_2 (X,f); } };
template<> struct real_fft_base_function< 3> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_3 (X,f); } };
template<> struct real_fft_base_function< 4> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_4 (X,f); } };
template<> struct real_fft_base_function< 5> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_5 (X,f); } };
template<> struct real_fft_base_function< 8> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_8 (X,f); } };
template<> struct real_fft_base_function<16> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_16(X,f); } };
template<> struct real_fft_base_function<32> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_32(X,f); } };
template<> struct real_fft_base_function<64> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { real_fft_64(X,f); } };

template<int M> struct real_fft_function
{
  template<class G1,class G2        > inline void operator()(const Vector<G1> &X, Vector<G2> &Y            ) { (*this)(X,Y,make_fft_adaptor()); }
  template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, const F &f) { real_fft_base_function<M>()(X,typename F::template rebind<M>::other(f).real_make(Y)); }
};


#ifdef SSE

template<class G1,class G2,class G3>
inline void real_crossed_fft_0_4(const Vector<G1> &a0, const Vector<G2> &Wa0,const Vector<G2> &Wb0, Vector<G3> &y0)
{
  TinyVector<2,complex<float> >::self t0,t1;
  {
    TinyVector<2,complex<float> >::self v0,v1;
    {
      TinyVector<4,float>::self r0,r1;
      r0=a0*Wa0; r1=a0*Wb0;
      v0=complexlo(r0,r1); v1=complexhi(r0,r1);
    }
    t0=v0+v1; t1=mi1mul(v0-v1);
 }
  y0=movelh(t0,t1)+movehl(t1,t0);
}

template<class G1,class G3>
inline void real_crossed_fft_0_4(const Vector<G1> &a0, Vector<G3> &y0,Vector<G3> &z0, const Vector<G1> &a1, Vector<G3> &y1, Vector<G3> &z1)
{
  typedef float real_type;
  TinyVector<2,complex<float> >::self za0,zb0;
  TinyVector<2,complex<float> >::self za1,zb1;

  {
    TinyVector<2,complex<float> >::self t0,t1;
    {
      TinyVector<2,complex<float> >::self v0,v1;
      {
        TinyVector<4,float>::self r0,r1;
        r0=a0; r1=0;
        v0=complexlo(r0,r1); v1=complexhi(r0,r1);
      }
      t0=v0+v1; t1=i1mul(v0-v1);
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      v0=movelh(t0,t1);
      v1=movehl(t1,t0);
      za0=v0+v1; zb0=v0-v1;
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      {
        TinyVector<4,float>::self r0,r1;
        r0=a1*TinyVector<4,float>::self(SIN32, SIN16, SIN00,-SIN16);
        r1=a1*TinyVector<4,float>::self(SIN00,-SIN16,-SIN32,-SIN16);
        v0=complexlo(r0,r1); v1=complexhi(r0,r1);
      }
      t0=v0+v1; t1=i1mul(v0-v1);
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      v0=movelh(t0,t1);
      v1=movehl(t1,t0);
      za1=v0+v1; zb1=v0-v1;
    }
  }
  y0=movelh(za0,za1); z0=movehl(za0,za1);
  y1=movehl(zb1,zb0); z1=movelh(zb1,zb0);
}


template<class G1,class G2,class G3>
inline void real_crossed_fft_0_4(const Vector<G1> &a0, Vector<G3> &y0, Vector<G3> &z0, const Vector<G1> &a1, const Vector<G1> &b1, const Vector<G2> &Wa1,const Vector<G2> &Wb1, Vector<G3> &y1, Vector<G3> &z1)
{
  TinyVector<2,complex<float> >::self za0,zb0;
  TinyVector<2,complex<float> >::self za1,zb1;
  {
    TinyVector<2,complex<float> >::self t0,t1;
    {
      TinyVector<2,complex<float> >::self v0,v1;
      {
        TinyVector<4,float>::self r0,r1;
        r0=a0; r1=0;
        v0=complexlo(r0,r1); v1=complexhi(r0,r1);
      }
      t0=v0+v1; t1=i1mul(v0-v1);
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      v0=movelh(t0,t1);
      v1=movehl(t1,t0);
      za0=v0+v1; zb0=v0-v1;
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      {
        TinyVector<4,float>::self r0,r1;
        r0=a1*Wa1-b1*Wb1; r1=b1*Wa1+a1*Wb1;
        v0=complexlo(r0,r1); v1=complexhi(r0,r1);
      }
      t0=v0+v1; t1=i1mul(v0-v1);
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      v0=movelh(t0,t1);
      v1=movehl(t1,t0);
      za1=v0+v1; zb1=v0-v1;
    }
  }
  y0=movelh(za0,za1); z0=movehl(za0,za1);
  y1=movehl(zb1,zb0); z1=movelh(zb1,zb0);
}


template<class G1,class G2,class G3>
inline void real_crossed_fft_0_4(const Vector<G1> &a0, const Vector<G1> &b0, const Vector<G2> &Wa0,const Vector<G2> &Wb0, Vector<G3> &y0, Vector<G3> &z0, const Vector<G1> &a1, const Vector<G1> &b1, const Vector<G2> &Wa1,const Vector<G2> &Wb1, Vector<G3> &y1, Vector<G3> &z1)
{
  TinyVector<2,complex<float> >::self za0,zb0;
  TinyVector<2,complex<float> >::self za1,zb1;
  {
    TinyVector<2,complex<float> >::self t0,t1;
    {
      TinyVector<2,complex<float> >::self v0,v1;
      {
        TinyVector<4,float>::self r0,r1;
        r0=a0*Wa0-b0*Wb0; r1=b0*Wa0+a0*Wb0;
        v0=complexlo(r0,r1); v1=complexhi(r0,r1);
      }
      t0=v0+v1; t1=i1mul(v0-v1);
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      v0=movelh(t0,t1);
      v1=movehl(t1,t0);
      za0=v0+v1; zb0=v0-v1;
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      {
        TinyVector<4,float>::self r0,r1;
        r0=a1*Wa1-b1*Wb1; r1=b1*Wa1+a1*Wb1;
        v0=complexlo(r0,r1); v1=complexhi(r0,r1);
      }
      t0=v0+v1; t1=i1mul(v0-v1);
    }
    {
      TinyVector<2,complex<float> >::self v0,v1;
      v0=movelh(t0,t1);
      v1=movehl(t1,t0);
      za1=v0+v1; zb1=v0-v1;
    }
  }
  y0=movelh(za0,za1); z0=movehl(za0,za1);
  y1=movehl(zb1,zb0); z1=movelh(zb1,zb0);
}
#endif //SSE


#ifdef SSE2

template<class G1,class G2,class G3>
inline void real_crossed_fft_0_2(const Vector<G1> &a, const Vector<G2> &Wa,const Vector<G2> &Wb, Vector<G3> &y, Vector<G3> &z)
{
  TinyVector<1,complex<double> >::self v0,v1;
  {
    TinyVector<2,double>::self r0,r1;
    r0=a*Wa; r1=a*Wb;
    v0=complexlo(r0,r1); v1=complexhi(r0,r1);
  }
  y=v0+v1; z=v0-v1;
}

template<class G1,class G3>
inline void real_crossed_fft_0_2(const Vector<G1> &a, Vector<G3> &y, Vector<G3> &z)
{
  TinyVector<1,complex<double> >::self v0,v1;
  {
    TinyVector<2,double>::self r0,r1;
    r0=a; r1=0;
    v0=complexlo(r0,r1); v1=complexhi(r0,r1);
  }
  y=v0+v1; z=v0-v1;
}

template<class G1,class G2,class G3>
inline void real_crossed_fft_0_2(const Vector<G1> &a, const Vector<G1> &b, const Vector<G2> &Wa,const Vector<G2> &Wb, Vector<G3> &y, Vector<G3> &z)
{
  TinyVector<1,complex<double> >::self v0,v1;
  {
    TinyVector<2,double>::self r0,r1;
    r0=a*Wa-b*Wb; r1=b*Wa+a*Wb;
    v0=complexlo(r0,r1); v1=complexhi(r0,r1);
  }
  y=v0+v1; z=v0-v1;
}
#endif //SSE2

template<int M,int N> struct real_crossed_fft_function {};

template<> struct real_crossed_fft_function<0,2>
{
  template<class G1,class G2,class G3> inline void operator()(const Vector<G1> &a                     , const Vector<G2> &Wa,const Vector<G2> &Wb, Vector<G3> &y, Vector<G3> &z) { real_crossed_fft_0_2(a  ,Wa,Wb,y,z); }
  template<class G1,         class G3> inline void operator()(const Vector<G1> &a                                                                , Vector<G3> &y, Vector<G3> &z) { real_crossed_fft_0_2(a        ,y,z); }
  template<class G1,class G2,class G3> inline void operator()(const Vector<G1> &a, const Vector<G1> &b, const Vector<G2> &Wa,const Vector<G2> &Wb, Vector<G3> &y, Vector<G3> &z) { real_crossed_fft_0_2(a,b,Wa,Wb,y,z); }
};

template<> struct real_crossed_fft_function<0,4>
{
  template<class G1,class G2,class G3> inline void operator()(const Vector<G1> &a0, const Vector<G2> &Wa0,const Vector<G2> &Wb0, Vector<G3> &y0) { real_crossed_fft_0_4(a0,Wa0,Wb0,y0); }
  template<class G1,         class G3> inline void operator()(const Vector<G1> &a0, Vector<G3> &y0, Vector<G3> &z0, const Vector<G1> &a1, Vector<G3> &y1, Vector<G3> &z1) { real_crossed_fft_0_4(a0,y0,z0,a1,y1,z1); }
  template<class G1,class G2,class G3> inline void operator()(const Vector<G1> &a0, Vector<G3> &y0, Vector<G3> &z0, const Vector<G1> &a1, const Vector<G1> &b1, const Vector<G2> &Wa1,const Vector<G2> &Wb1, Vector<G3> &y1, Vector<G3> &z1) { real_crossed_fft_0_4(a0,y0,z0, a1,b1,Wa1,Wb1,y1,z1); }
  template<class G1,class G2,class G3> inline void operator()(const Vector<G1> &a0, const Vector<G1> &b0, const Vector<G2> &Wa0,const Vector<G2> &Wb0, Vector<G3> &y0, Vector<G3> &z0, const Vector<G1> &a1, const Vector<G1> &b1, const Vector<G2> &Wa1,const Vector<G2> &Wb1, Vector<G3> &y1, Vector<G3> &z1) { real_crossed_fft_0_4(a0,b0,Wa0,Wb0,y0,z0, a1,b1,Wa1,Wb1,y1,z1); }
};

template<class G,class S>
inline void real_fft_2_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<2,real_type>::self value_type;
  typedef typename TinyVector<1,complex_type>::self complex_value_type;
  assert(X.size()==2);

  value_type t0,t1;

  t0=X[0]; t1=X[1];

  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t0+t1                                 ,y,z); Y.first0(  y); Y.firstN(z); }
  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t0-t1,value_type(1,0),value_type(0,-1),y,z); Y.first (1,y); }
}

template<class G,class S>
inline void real_fft_2_4(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<4,real_type>::self value_type;
  typedef typename SimdVector<2,complex_type>::self complex_value_type;
  assert(X.size()==2);

  value_type t0,t1;

  t0=X[0]; t1=X[1];

  {
	  complex_value_type y0,z0,y1,z1;
    real_crossed_fft_function<0,4>()(t0+t1,y0,z0,t0-t1,y1,z1);
    Y.first0(y0);
    Y.firstN(z1[1]);
    Y.first(2,y1);
	}
}

template<class G,class S>
inline void real_fft_4_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<2,real_type>::self value_type;
  typedef typename TinyVector<1,complex_type>::self complex_value_type;
  assert(X.size()==4);

  value_type t0,t1,t2,t3;

  { value_type v0,v1; v0=X[0]; v1=X[2]; t0=v0+v1; t1=v0-v1; v0=X[3]; v1=X[1]; t2=v0+v1; t3=v0-v1; }

  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t0+t2                                                  ,y,z); Y.first0(  y); Y.firstN(  z); }
  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t1   ,t3,value_type(1,SIN16),value_type(0,-SIN16),y,z); Y.first (1,y); Y.second(5,z); }
  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t0-t2   ,value_type(1,0       ),value_type(0,-1       ),y,z); Y.first (2,y);                }
}

template<class G,class S>
inline void real_fft_4_4(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<4,real_type>::self value_type;
  typedef typename TinyVector<2,complex_type>::self complex_value_type;
  assert(X.size()==4);

  value_type t0,t1,t2,t3;

  { value_type v0,v1; v0=X[0]; v1=X[2]; t0=v0+v1; t1=v0-v1; v0=X[3]; v1=X[1]; t2=v0+v1; t3=v0-v1; }

  {
	  complex_value_type y0,z0,y1,z1;
    real_crossed_fft_function<0,4>()(t0+t2,y0,z0,t1,t3,value_type(SIN32,SIN24,SIN16,SIN08),value_type(SIN00,-SIN08,-SIN16,-SIN24),y1,z1);
    Y.first0(y0);
    Y.secondN(z1);
    Y.second (12,z0);
    Y.first (5,y1[1]);
	}
  {
	  complex_value_type y0;
    real_crossed_fft_function<0,4>()(t0-t2,value_type(SIN32, SIN16, SIN00, -SIN16),value_type(SIN00,-SIN16,-SIN32,-SIN16),y0);
    Y.first(2,y0[0]);
    Y.first(6,y0[1]);
	}
}

template<class G,class S>
inline void real_fft_8_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<2,real_type>::self value_type;
  typedef typename TinyVector<1,complex_type>::self complex_value_type;
  assert(X.size()==8);

  value_type t0,t1,t2,t3,t4,t5,t6,t7;
  { value_type v0,v1,v2,v3; v0=X[0]; v1=X[4]; v2=v0+v1; t4=v0-v1; v0=X[2]; v1=X[6]; v3=v0+v1; t7=v0-v1; t2=v2+v3; t0=v2-v3; }
	{ value_type v0,v1,v2,v3,v4,v5; v0=X[7]; v1=X[3]; v2=v0+v1; v4=v0-v1; v0=X[1]; v1=X[5]; v3=v0+v1; v5=v0-v1; t5=SIN16*(v4+v5); t6=SIN16*(v4-v5); t3=v2+v3; t1=v2-v3; }

  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t2+t3                                                       ,y,z); Y.first0(  y); Y.firstN(   z); }
  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t4+t5,t6-t7,value_type(SIN32, SIN24),value_type(SIN00,- SIN08),y,z); Y.first (1,y); Y.second( 9,z); }
	{ complex_value_type y,z; real_crossed_fft_function<0,2>()(t0   ,t1   ,value_type(SIN32, SIN16),value_type(0,- SIN16),y,z); Y.first (2,y); Y.second(10,z); }
  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t4-t5,t6+t7,value_type(SIN32, SIN08),value_type(0,-SIN24),y,z); Y.first (3,y); Y.second(11,z); }
  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t2-t3      ,value_type(SIN32, SIN00),value_type(0,-        1),y,z); Y.first (4,y);                 }
}

template<class G,class S>
inline void real_fft_8_4(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<4,real_type>::self value_type;
  typedef typename TinyVector<2,complex_type>::self complex_value_type;
  assert(X.size()==8);

  value_type t0,t1,t2,t3,t4,t5,t6,t7;
  { value_type v0,v1,v2,v3; v0=X[0]; v1=X[4]; v2=v0+v1; t4=v0-v1; v0=X[2]; v1=X[6]; v3=v0+v1; t7=v0-v1; t2=v2+v3; t0=v2-v3; }
	{ value_type v0,v1,v2,v3,v4,v5; v0=X[7]; v1=X[3]; v2=v0+v1; v4=v0-v1; v0=X[1]; v1=X[5]; v3=v0+v1; v5=v0-v1; t5=SIN16*(v4+v5); t6=SIN16*(v4-v5); t3=v2+v3; t1=v2-v3; }

  {
	  complex_value_type y0,z0,y1,z1;
    real_crossed_fft_function<0,4>()(t2+t3,y0,z0,t4+t5,t6-t7,value_type(SIN32, SIN28, SIN24, SIN20),value_type(SIN00,-SIN04,-SIN08,-SIN12),y1,z1);
    Y.first0(y0);      Y.secondN(   z1);
    Y.first (9,y1[1]); Y.second (24,z0);
	}
	{
	  complex_value_type y0,z0,y1,z1;
    real_crossed_fft_function<0,4>()(t0   ,t1   ,value_type(SIN32, SIN24, SIN16, SIN08),value_type(SIN00,-SIN08,-SIN16,-SIN24),y0,z0,
                                          t4-t5,t6+t7,value_type(SIN32, SIN20, SIN08,-SIN04),value_type(SIN00,-SIN12,-SIN24,-SIN28),y1,z1);
    Y.first( 2,y0); Y.second(26,z0);
    Y.first(10,y1); Y.second(18,z1);
	}
  {
	  complex_value_type y0;
    real_crossed_fft_function<0,4>()(t2-t3,value_type(SIN32, SIN16, SIN00, -SIN16),value_type(SIN00,-SIN16,-SIN32,-SIN16),y0);
    Y.first(4,y0[0]); Y.first(12,y0[1]);
	}
}

template<class G,class S>
inline void real_fft_16_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<2,real_type>::self value_type;
  typedef typename TinyVector<1,complex_type>::self complex_value_type;
  assert(X.size()==16);

  value_type t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t15,t14;

  { value_type v0,v1; v0=X[ 0]; v1=X[ 8]; t0=v0+v1; t1=v0-v1; v0=X[ 4]; v1=X[12]; t2=v0-v1; v0+=v1; t3=t0+v0; t0-=v0; }
  { value_type v0,v1; v0=X[ 2]; v1=X[10]; t6=v0+v1; t7=v0-v1; v0=X[14]; v1=X[ 6]; t4=v0+v1; v0=v0-v1; t5=t4+t6; t4=t4-t6; t6=SIN16*(v0+t7); t7=SIN16*(v0-t7); }
  {
    value_type v6,v7;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[15];v5=X[ 7]; v0=v4+v5; v1=v4-v5; v4=X[3]; v5=X[11]; v2=v4+v5; v3=v4-v5; } t10=v0+v2; v6=v0-v2; t8 =SIN08 *v1-SIN24*v3; t9 =SIN24*v1+SIN08 *v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[1]; v5=X[ 9]; v0=v4+v5; v1=v4-v5; v4=X[5]; v5=X[13]; v2=v4+v5; v3=v4-v5; } t14=v0+v2; v7=v0-v2; t12=SIN08 *v1+SIN24*v3; t13=SIN24*v1-SIN08 *v3; }
    t15=SIN16; t11=t15*(v6+v7); t15=t15*(v6-v7);
  }

  {
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(t0+t11,t15+t4,value_type(1,SIN24),value_type(0,-SIN08),y,z); Y.first( 2,y); Y.second(18,z); }
	  { complex_value_type y,z; real_crossed_fft_function<0,2>()(t0-t11,t15-t4,value_type(1,SIN08),value_type(0,-SIN24),y,z); Y.first( 6,y); Y.second(22,z); }
	}
  {
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(t3-t5,t10-t14,value_type(1,SIN16),value_type(0,-SIN16),y,z); Y.first( 4,y); Y.second(20,z); }
    value_type v0, v1; v0=t3+t5; v1=t10+t14;
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v1,y,z); Y.first0(y); Y.firstN(z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v1,value_type(1, 0),value_type(0,-1),y,z); Y.first( 8,y); }
  }
  {
    value_type v0,v1,v2,v3; v0=t1+t6; v1=t13+t9; v2=t7-t2; v3=t8-t12;
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v1,v2+v3,value_type(1,SIN28),value_type(0,-SIN04),y,z); Y.first (1,y); Y.second(17,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v1,v3-v2 ,value_type(1,SIN04),value_type(0,-SIN28),y,z); Y.first(7,y); Y.second(23,z); }
	}
	{
    value_type v0,v1,v2,v3; v0=t1-t6; v1=t12+t8; v2=t2+t7; v3=t9-t13;
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v1 ,v2+v3 ,value_type(1,SIN20),value_type(0,-SIN12),y,z); Y.first(3,y); Y.second(19,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v1 ,v3-v2 ,value_type(1,SIN12),value_type(0,-SIN20),y,z); Y.first(5,y); Y.second(21,z); }
  }
}

template<class G,class S>
inline void real_fft_32_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<2,real_type>::self value_type;
  typedef typename TinyVector<1,complex_type>::self complex_value_type;
  assert(X.size()==32);

  value_type t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;

  { value_type v0,   v2;    { value_type v4,v5; v4=X[ 0]; v5=X[16]; v0=v4+v5; t2=v4-v5; v4=X[ 8]; v5=X[24]; v2=v4+v5; t3=v4-v5; } t0=v0+v2; t1=v0-v2; }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[ 4]; v5=X[20]; v0=v4+v5; v1=v4-v5; v4=X[28]; v5=X[12]; v2=v4+v5; v3=v4-v5; } t4=v2+v0; t5=v2-v0; t6=SIN16*(v3+v1); t7=SIN16*(v3-v1); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[30]; v5=X[14]; v0=v4+v5; v1=v4-v5; v4=X[ 6]; v5=X[22]; v2=v4+v5; v3=v4-v5; } t8=v0+v2; t9=v0-v2; t10=SIN24*v1+SIN08*v3; t11=SIN08*v1-SIN24*v3; }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[ 2]; v5=X[18]; v0=v4+v5; v1=v4-v5; v4=X[10]; v5=X[26]; v2=v4+v5; v3=v4-v5; } t12=v0+v2; t13=v0-v2; t14=SIN24*v1-SIN08*v3; t15=SIN08*v1+SIN24*v3; }
  { value_type v0,v1,v2,v3,v6,v7; { value_type v4,v5; v4=X[31]; v5=X[15]; v0=v4-v5; v1=v4+v5; v4=X[ 7]; v5=X[23]; v2=v4-v5; v3=v4+v5; } t16=v1+v3; t19=v1-v3; { value_type v4,v5; v4=X[ 3]; v5=X[19]; v6=v4-v5; v7=v4+v5; v4=X[27]; v5=X[11]; v1=v4-v5; v3=v4+v5; } t17=v3+v7; t18=v3-v7; v3=SIN16*(v1-v6); v7=SIN16*(v1+v6); t20=v3-v2; t21=v3+v2; t22=v0+v7; t23=v0-v7; }
  { value_type v0,v1,v2,v3,v6,v7; { value_type v4,v5; v4=X[ 1]; v5=X[17]; v0=v4-v5; v1=v4+v5; v4=X[ 9]; v5=X[25]; v2=v4-v5; v3=v4+v5; } t24=v1+v3; t27=v1-v3; { value_type v4,v5; v4=X[ 5]; v5=X[21]; v6=v4-v5; v7=v4+v5; v4=X[29]; v5=X[13]; v1=v4-v5; v3=v4+v5; } t25=v3+v7; t26=v3-v7; v3=SIN16*(v1-v6); v7=SIN16*(v1+v6); t28=v3-v2; t29=v3+v2; t30=v0+v7; t31=v0-v7; }

  {
    value_type v0,v1,v2,v3,v4,v5;
    v0=t0+t4; v1=t12+t8; v2=v0+v1; v3=t24+t25; v4=t16+t17; v5=v3+v4;
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v1,v4-v3,value_type(1,SIN16),value_type(0,-SIN16),y,z); Y.first( 8,y); Y.second(40,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v2+v5,y,z); Y.first0(y); Y.firstN(z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v2-v5,value_type(1, 0),value_type(0,-1),y,z); Y.first(16,y); }
  }
  {
    value_type v0,v1,v2,v3;
    v0=t24-t25; v1=t16-t17; v2=SIN16*(v1+v0); v3=SIN16*(v1-v0); v0=t0-t4; v1=t8-t12;
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v2,v1+v3,value_type(1,SIN24),value_type(0,-SIN08),y,z); Y.first( 4,y); Y.second(36,z); }
	  { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v2,v3-v1,value_type(1,SIN08),value_type(0,-SIN24),y,z); Y.first(12,y); Y.second(44,z); }
  }
  {
    value_type v0,v1,v2,v3,v4,v5,v6,v7;
    v1=SIN24*t27+SIN08*t26; v5=SIN24*t19-SIN08*t18; v4=v5+v1; v5=v5-v1;
    v1=SIN24*t26-SIN08*t27; v7=SIN08*t19+SIN24*t18; v6=v7+v1; v7=v7-v1;
    v1=SIN16*(t9+t13); v0=t1+v1; v1=t1-v1;
    v3=SIN16*(t9-t13); v2=v3-t5; v3=v3+t5;
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v4,v6+v3,value_type(1,SIN28),value_type(0,-SIN04),y,z); Y.first( 2,y); Y.second(34,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v4,v6-v3,value_type(1,SIN04),value_type(0,-SIN28),y,z); Y.first(14,y); Y.second(46,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1+v7,v5+v2,value_type(1,SIN20),value_type(0,-SIN12),y,z); Y.first( 6,y); Y.second(38,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1-v7,v5-v2,value_type(1,SIN12),value_type(0,-SIN20),y,z); Y.first(10,y); Y.second(42,z); }
  }
  {
    value_type v0,v1,v2,v3,v4,v5,v6,v7;
    { value_type v8,v9; v8=t2+t6; v9=t14+t10; v0=v8+v9; v1=v8-v9; v8=SIN28*t28-SIN04*t30; v9=SIN04*t22+SIN28*t20; v2=v9+v8;  v3=v9-v8; }
    { value_type v8,v9; v8=SIN28*t30+SIN04*t28; v9=SIN28*t22-SIN04*t20; v4=v9+v8; v5=v9-v8; v8=t11-t15; v9=t7-t3; v6=v8-v9; v7=v8+v9; }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v4,v2+v7,value_type(1,SIN30),value_type(0,-SIN02),y,z); Y.first( 1,y); Y.second(33,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v4,v2-v7,value_type(1,SIN02),value_type(0,-SIN30),y,z); Y.first(15,y); Y.second(47,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1+v3,v5+v6,value_type(1,SIN18),value_type(0,-SIN14),y,z); Y.first( 7,y); Y.second(39,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1-v3,v5-v6,value_type(1,SIN14),value_type(0,-SIN18),y,z); Y.first( 9,y); Y.second(41,z); }
  }
  {
    value_type v0,v1,v2,v3,v4,v5,v6,v7;
    { value_type v8,v9; v8=t2-t6; v9=t15+t11; v0=v8+v9; v1=v8-v9; v8=SIN20*t29-SIN12*t31; v9=SIN12*t23+SIN20*t21; v2=v9+v8; v3=v9-v8; }
    { value_type v8,v9; v8=SIN20*t31+SIN12*t29; v9=SIN20*t23-SIN12*t21; v4=v9+v8; v5=v9-v8; v8=t10-t14; v9=t3+t7; v6=v8-v9; v7=v8+v9; }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v4,v2+v7,value_type(1,SIN26),value_type(0,-SIN06),y,z); Y.first( 3,y); Y.second(35,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v4,v2-v7,value_type(1,SIN06),value_type(0,-SIN26),y,z); Y.first(13,y); Y.second(45,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1+v3,v5+v6,value_type(1,SIN22),value_type(0,-SIN10),y,z); Y.first( 5,y); Y.second(37,z); }
    { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1-v3,v5-v6,value_type(1,SIN10),value_type(0,-SIN22),y,z); Y.first(11,y); Y.second(43,z); }
  }
}

template<class G,class S>
inline void real_fft_64_2(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<2,real_type>::self value_type;
  typedef typename TinyVector<1,complex_type>::self complex_value_type;
  assert(X.size()==64);

	value_type t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65;

	{ value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 0]; u1=X[32]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[56]; u1=X[24]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[16]; u1=X[48]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[ 8]; u1=X[40]; v6=u0+u1; v7=u0-u1; } t00=v0-v4; t01=v2-v6; t06=v0+v4; t07=v6+v2; { value_type u0; u0=SIN16*(v3+v7); t02=v1+u0; t03=v1-u0; } { value_type u0; u0=SIN16*(v3-v7); t04=u0-v5; t05=v5+u0; } t08=t06+t07; }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[57]; u1=X[25]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[ 1]; u1=X[33]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[ 9]; u1=X[41]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[17]; u1=X[49]; v6=u0-u1; v7=u0+u1; } t09=v3-v7; t10=v1-v5; { value_type u0,u1; u0=v3+v7; u1=v5+v1; t15=u0+u1; t16=u0-u1; } { value_type u0; u0=SIN16*(v0-v4); t11=u0-v6; t12=u0+v6; } { value_type u0; u0=SIN16*(v4+v0); t13=v2+u0; t14=v2-u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[63]; u1=X[31]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[15]; u1=X[47]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[ 7]; u1=X[39]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[55]; u1=X[23]; v6=u0-u1; v7=u0+u1; } t17=v1-v3; t18=v7-v5; { value_type u0,u1; u0=v1+v3; u1=v5+v7; t23=u0+u1; t24=u0-u1; } { value_type u0; u0=SIN16*(v4+v6); t19=v0+u0; t20=v0-u0; } { value_type u0; u0=SIN16*(v6-v4); t21=u0-v2; t22=u0+v2; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[62]; u1=X[30]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[54]; u1=X[22]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[14]; u1=X[46]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[ 6]; u1=X[38]; v6=u0+u1; v7=u0-u1; } { value_type u0; u0=SIN16*(v7+v3); t25=v1+u0; t26=v1-u0; } { value_type u0; u0=SIN16*(v3-v7); t27=u0-v5; t28=v5+u0; } { value_type u0,u1; u0=v0+v4; u1=v6+v2; t29=u0+u1; t30=u0-u1; } { value_type u0,u1; u0=v0-v4; u1=v2-v6; t31=SIN24*u0-SIN08*u1; t32=SIN24*u1+SIN08*u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 2]; u1=X[34]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[58]; u1=X[26]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[18]; u1=X[50]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[10]; u1=X[42]; v6=u0+u1; v7=u0-u1; } { value_type u0; u0=SIN16*(v3-v7); t33=u0-v5; t34=v5+u0; } { value_type u0; u0=SIN16*(v7+v3); t35=v1+u0; t36=v1-u0; } { value_type u0,u1; u0=v0+v4; u1=v6+v2; t37=u0+u1; t38=u0-u1; } { value_type u0,u1; u0=v0-v4; u1=v2-v6; t39=SIN24*u0+SIN08*u1; t40=SIN24*u1-SIN08*u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 4]; u1=X[36]; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=X[12]; u1=X[44]; v2=u0+u1; v3=u0-u1; } { value_type u0,u1; u0=X[20]; u1=X[52]; v4=u0+u1; v5=u0-u1; } { value_type u0,u1; u0=X[60]; u1=X[28]; v6=u0+u1; v7=u0-u1; } t41=v0+v4; t42=v6+v2; t43=t41+t42;                              { value_type u0,u1; u0=v0-v4; u1=v6-v2; t44=SIN16*(u0+u1); t45=SIN16*(u1-u0); } { value_type u0,u1; u0=SIN24*v1-SIN08*v5; u1=SIN24*v7+SIN08*v3; t46=u0+u1; t47=u1-u0; } { value_type u0,u1; u0=SIN08*v7-SIN24*v3; u1=SIN08*v1+SIN24*v5; t48=u0-u1; t49=u1+u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[61]; u1=X[29]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[21]; u1=X[53]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[13]; u1=X[45]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[ 5]; u1=X[37]; v6=u0-u1; v7=u0+u1; } { value_type u0,u1; u0=v7+v3; u1=v1+v5; t52=u0+u1; t53=u1-u0; } { value_type u0,u1; u0=v7-v3; u1=v1-v5; t56=SIN16*(u0+u1); t57=SIN16*(u1-u0); } { value_type u0,u1; u0=SIN08*v0-SIN24*v4; u1=SIN08*v6+SIN24*v2; t50=u0-u1; t51=u1+u0; } { value_type u0,u1; u0=SIN24*v6-SIN08*v2; u1=SIN24*v0+SIN08*v4; t54=u0+u1; t55=u1-u0; } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=X[ 3]; u1=X[35]; v0=u0-u1; v1=u0+u1; } { value_type u0,u1; u0=X[11]; u1=X[43]; v2=u0-u1; v3=u0+u1; } { value_type u0,u1; u0=X[19]; u1=X[51]; v4=u0-u1; v5=u0+u1; } { value_type u0,u1; u0=X[59]; u1=X[27]; v6=u0-u1; v7=u0+u1; } { value_type u0,u1; u0=v1+v5; u1=v7+v3; t60=u0+u1; t61=u1-u0; } { value_type u0,u1; u0=v1-v5; u1=v7-v3; t64=SIN16*(u0+u1); t65=SIN16*(u1-u0); } { value_type u0,u1; u0=SIN24*v0-SIN08*v4; u1=SIN24*v6+SIN08*v2; t58=u0+u1; t59=u1-u0; } { value_type u0,u1; u0=SIN08*v6-SIN24*v2; u1=SIN08*v0+SIN24*v4; t62=u0-u1; t63=u1+u0; } }

  { value_type v0,v1,v2,v3,v4,v5; v0=t08+t43; v1=t37+t29; v2=v0+v1; v3=t15+t52; v4=t23+t60; v5=v3+v4; { complex_value_type y,z; real_crossed_fft_function<0,2>()(v2+v5,y,z); Y.first0(y); Y.firstN(z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v1,v4-v3,value_type(1,SIN16),value_type(0,-SIN16),y,z); Y.first(16,y); Y.second(80,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v2-v5,value_type(1, 0),value_type(0,-1),y,z); Y.first(32,y); } }
  { value_type v0,v1,v2,v3; v0=t08-t43; v1=t29-t37; { value_type u0,u1; u1=t15-t52; u0=t23-t60; v2=SIN16*(u0+u1); v3=SIN16*(u0-u1); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v2,v1+v3,value_type(1,SIN24),value_type(0,-SIN08),y,z); Y.first( 8,y); Y.second(72,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v2,v3-v1,value_type(1,SIN08),value_type(0,-SIN24),y,z); Y.first(24,y); Y.second(88,z); } }
  { value_type v0,v1,v2,v3,v4,v5,v7,v6; { value_type u0,u1; u0=t06-t07; u1=SIN16*(t38+t30); v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=SIN24*t53-SIN08*t16; u1=SIN24*t61+SIN08*t24; v2=u0+u1; v3=u1-u0; } { value_type u0,u1; u0=SIN24*t16+SIN08*t53; u1=SIN24*t24-SIN08*t61; v4=u0+u1; v5=u1-u0; } { value_type u0,u1; u0=SIN16*(t30-t38); u1=t42-t41; v6=u0+u1; v7=u0-u1;  } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+v4,v2+v6,value_type(1,SIN28),value_type(0,-SIN04),y,z); Y.first( 4,y); Y.second(68,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-v4,v2-v6,value_type(1,SIN04),value_type(0,-SIN28),y,z); Y.first(28,y); Y.second(92,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1+v3,v5+v7,value_type(1,SIN20),value_type(0,-SIN12),y,z); Y.first(12,y); Y.second(76,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1-v3,v5-v7,value_type(1,SIN12),value_type(0,-SIN20),y,z); Y.first(20,y); Y.second(84,z); } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t31-t39; u1=t45-t01; v2=u0-u1; v3=u1+u0; } { value_type u0,u1; u0=t00-t44; u1=t32-t40; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=t09-t56; u1=t57-t10; v4=SIN20*u0+SIN12*u1; v5=SIN20*u1-SIN12*u0; } { value_type u0,u1; u0=t17-t64; u1=t65-t18; v6=SIN20*u0-SIN12*u1; v7=SIN20*u1+SIN12*u0; } { value_type u0,u1; u0=v4+v6; u1=v5+v7; { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+u0,v3+u1,value_type(1,SIN26),value_type(0,-SIN06),y,z); Y.first( 6,y); Y.second(70,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-u0,u1-v3,value_type(1,SIN06),value_type(0,-SIN26),y,z); Y.first(26,y); Y.second(90,z); } } { value_type u0,u1; u0=v7-v5; u1=v6-v4; { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1+u0,v2+u1,value_type(1,SIN22),value_type(0,-SIN10),y,z); Y.first(10,y); Y.second(74,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1-u0,u1-v2,value_type(1,SIN10),value_type(0,-SIN22),y,z); Y.first(22,y); Y.second(86,z); } } }
  { value_type v0,v1,v2,v3,v4,v5,v6,v7; { value_type u0,u1; u0=t00+t44; u1=t39+t31; v0=u0+u1; v1=u0-u1; } { value_type u0,u1; u0=t40+t32; u1=t01+t45; v2=u0-u1; v3=u1+u0; } { value_type u0,u1; u0=t09+t56; u1=t10+t57; v4=SIN28*u0+SIN04*u1; v5=SIN28*u1-SIN04*u0; } { value_type u0,u1; u0=t17+t64; u1=t18+t65; v6=SIN28*u0-SIN04*u1; v7=SIN28*u1+SIN04*u0; } { value_type u0,u1; u0=v4+v6; u1=v5+v7; { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0+u0,v3+u1,value_type(1,SIN30),value_type(0,-SIN02),y,z); Y.first( 2,y); Y.second(66,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v0-u0,u1-v3,value_type(1,SIN02),value_type(0,-SIN30),y,z); Y.first(30,y); Y.second(94,z); } } { value_type u0,u1; u0=v6-v4; u1=v7-v5; { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1+u1,v2+u0,value_type(1,SIN18),value_type(0,-SIN14),y,z); Y.first(14,y); Y.second(78,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(v1-u1,u0-v2,value_type(1,SIN14),value_type(0,-SIN18),y,z); Y.first(18,y); Y.second(82,z); } } }
  {
    value_type v00,v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12,v13,v14,v15; v00=t02-t46; v01=t02+t46; v02=t04+t48; v03=t48-t04; { value_type u0,u1; u0=SIN28*t27+SIN04*t25; u1=SIN28*t33-SIN04*t35; v04=u0-u1; v05=u1+u0; } { value_type u0,u1; u0=t19+t58; u1=t21+t62; v06=SIN30*u0-SIN02*u1; v07=SIN30*u1+SIN02*u0; } { value_type u0,u1; u0=SIN28*t25-SIN04*t27; u1=SIN28*t35+SIN04*t33; v08=u0-u1; v09=u1+u0; } { value_type u0,u1; u0=t50-t11; u1=t13-t54; v10=SIN18*u1+SIN14*u0; v11=SIN18*u0-SIN14*u1; } { value_type u0,u1; u0=t11+t50; u1=t13+t54; v12=SIN30*u1+SIN02*u0; v13=SIN30*u0-SIN02*u1; } { value_type u0,u1; u0=t19-t58; u1=t62-t21; v14=SIN18*u0-SIN14*u1; v15=SIN18*u1+SIN14*u0; }
    { value_type u0,u1,u2,u3; u0=v01+v09; u1=v12+v06; u2=v02+v05; u3=v13+v07; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN31),value_type(0,-SIN01),y,z); Y.first( 1,y); Y.second(65,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN01),value_type(0,-SIN31),y,z); Y.first(31,y); Y.second(95,z); } }
    { value_type u0,u1,u2,u3; u0=v01-v09; u1=v07-v13; u2=v05-v02; u3=v06-v12; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN17),value_type(0,-SIN15),y,z); Y.first(15,y); Y.second(79,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN15),value_type(0,-SIN17),y,z); Y.first(17,y); Y.second(81,z); } }
    { value_type u0,u1,u2,u3; u0=v00+v04; u1=v10+v14; u2=v03+v08; u3=v11+v15; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN25),value_type(0,-SIN07),y,z); Y.first( 7,y); Y.second(71,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN07),value_type(0,-SIN25),y,z); Y.first(25,y); Y.second(89,z); } }
    { value_type u0,u1,u2,u3; u0=v00-v04; u1=v15-v11; u2=v08-v03; u3=v14-v10; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN23),value_type(0,-SIN09),y,z); Y.first( 9,y); Y.second(73,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN09),value_type(0,-SIN23),y,z); Y.first(23,y); Y.second(87,z); } }
  }
  {
    value_type v00,v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12,v13,v14,v15; v00=t03+t49; v01=t03-t49; v02=t47-t05; v03=t05+t47; { value_type u0,u1; u0=SIN20*t36+SIN12*t34; u1=SIN20*t26-SIN12*t28; v04=u0+u1; v05=u1-u0; } { value_type u0,u1; u0=t20-t63; u1=t59-t22; v06=SIN22*u0-SIN10*u1; v07=SIN22*u1+SIN10*u0; } { value_type u0,u1; u0=SIN20*t34-SIN12*t36; u1=SIN20*t28+SIN12*t26; v08=u0+u1; v09=u1-u0; } { value_type u0,u1; u0=t14+t51; u1=t12+t55; v10=SIN26*u0+SIN06*u1; v11=SIN26*u1-SIN06*u0; } { value_type u0,u1; u0=t14-t51; u1=t55-t12; v12=SIN22*u0+SIN10*u1; v13=SIN22*u1-SIN10*u0; } { value_type u0,u1; u0=t20+t63; u1=t22+t59; v14=SIN26*u0-SIN06*u1; v15=SIN26*u1+SIN06*u0; }
    { value_type u0,u1,u2,u3; u0=v00+v04; u1=v10+v14; u2=v03+v08; u3=v11+v15; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN29),value_type(0,-SIN03),y,z); Y.first( 3,y); Y.second(67,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN03),value_type(0,-SIN29),y,z); Y.first(29,y); Y.second(93,z); } }
    { value_type u0,u1,u2,u3; u0=v00-v04; u1=v15-v11; u2=v08-v03; u3=v14-v10; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN19),value_type(0,-SIN13),y,z); Y.first(13,y); Y.second(77,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN13),value_type(0,-SIN19),y,z); Y.first(19,y); Y.second(83,z); } }
    { value_type u0,u1,u2,u3; u0=v01+v09; u1=v12+v06; u2=v02+v05; u3=v13+v07; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN27),value_type(0,-SIN05),y,z); Y.first( 5,y); Y.second(69,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN05),value_type(0,-SIN27),y,z); Y.first(27,y); Y.second(91,z); } }
    { value_type u0,u1,u2,u3; u0=v01-v09; u1=v07-v13; u2=v05-v02; u3=v06-v12; { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0+u1,u2+u3,value_type(1,SIN21),value_type(0,-SIN11),y,z); Y.first(11,y); Y.second(75,z); } { complex_value_type y,z; real_crossed_fft_function<0,2>()(u0-u1,u3-u2,value_type(1,SIN11),value_type(0,-SIN21),y,z); Y.first(21,y); Y.second(85,z); } }
  }
}

template<class G,class S>
inline void real_fft_16_4(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<4,real_type>::self value_type;
  typedef typename TinyVector<2,complex_type>::self complex_value_type;
  assert(X.size()==16);

  value_type t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t15,t14;

  { value_type v0,v1; v0=X[ 0]; v1=X[ 8]; t0=v0+v1; t1=v0-v1; v0=X[ 4]; v1=X[12]; t2=v0-v1; v0+=v1; t3=t0+v0; t0-=v0; }
  { value_type v0,v1; v0=X[ 2]; v1=X[10]; t6=v0+v1; t7=v0-v1; v0=X[14]; v1=X[ 6]; t4=v0+v1; v0-=v1; t5=t4+t6; t4=t4-t6; t6=SIN16*(v0+t7); t7=SIN16*(v0-t7); }
  {
    value_type v6,v7;
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[15];v5=X[ 7]; v0=v4+v5; v1=v4-v5; v4=X[3]; v5=X[11]; v2=v4+v5; v3=v4-v5; } t10=v0+v2; v6=v0-v2; t8 =SIN08 *v1-SIN24*v3; t9 =SIN24*v1+SIN08 *v3; }
    { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[1]; v5=X[ 9]; v0=v4+v5; v1=v4-v5; v4=X[5]; v5=X[13]; v2=v4+v5; v3=v4-v5; } t14=v0+v2; v7=v0-v2; t12=SIN08 *v1+SIN24*v3; t13=SIN24*v1-SIN08 *v3; }
    t15=SIN16; t11=t15*(v6+v7); t15=t15*(v6-v7);
  }

  {
    value_type v0,v1,v2,v3, v4, v5; v0=t1+t6; v1=t13+t9; v2=t7-t2; v3=t8-t12; v4=t3+t5; v5=t10+t14;
    { complex_value_type y0;          real_crossed_fft_function<0,4>()(v4-v5,        value_type(SIN32,SIN16,SIN00,-SIN16),value_type(SIN00,-SIN16,-SIN32,-SIN16),y0   ); Y.first( 8,y0[0]); Y.first(24,y0[1]); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v4+v5,                                                                                    y0,z0,  v0+v1,v2+v3,value_type(SIN32,SIN30,SIN28, SIN26),value_type(SIN00,-SIN02,-SIN04,-SIN06),y1,z1); Y.first0(   y0); Y.secondN(   z1); Y.first(17,y1[1]); Y.second(48,z0); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(t0-t11,t15-t4,value_type(SIN32,SIN20,SIN08,-SIN04),value_type(SIN00,-SIN12,-SIN24,-SIN28),y0,z0,  v0-v1,v3-v2,value_type(SIN32,SIN18,SIN04,-SIN10),value_type(SIN00,-SIN14,-SIN28,-SIN22),y1,z1); Y.first ( 6,y0); Y.second (54,z0); Y.first(22,y1   ); Y.second(38,z1); }
	}
	{
    value_type v0,v1,v2,v3; v0=t1-t6; v1=t12+t8; v2=t2+t7; v3=t9-t13;
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(t0+t11,t15+t4,value_type(SIN32,SIN28,SIN24,SIN20),value_type(SIN00,-SIN04,-SIN08,-SIN12),y0,z0,  v0+v1 ,v2+v3 ,value_type(SIN32,SIN26,SIN20,SIN14),value_type(SIN00,-SIN06,-SIN12,-SIN18),y1,z1); Y.first( 2,y0); Y.second(50,z0); Y.first(18,y1); Y.second(34,z1); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(t3-t5,t10-t14,value_type(SIN32,SIN24,SIN16,SIN08),value_type(SIN00,-SIN08,-SIN16,-SIN24),y0,z0,  v0-v1 ,v3-v2 ,value_type(SIN32,SIN22,SIN12,SIN02),value_type(SIN00,-SIN10,-SIN20,-SIN30),y1,z1); Y.first( 4,y0); Y.second(52,z0); Y.first(20,y1); Y.second(36,z1); }
  }
}

template<class G,class S>
inline void real_fft_32_4(const Vector<G> &X, const S &Y)
{
  typedef PROMOTE2(float,typename Vector<G>::value_type::value_type) real_type;
  typedef complex<real_type> complex_type;
  typedef typename TinyVector<4,real_type>::self value_type;
  typedef typename TinyVector<2,complex_type>::self complex_value_type;
  assert(X.size()==32);

  value_type t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;

  { value_type v0,   v2;    { value_type v4,v5; v4=X[ 0]; v5=X[16]; v0=v4+v5; t2=v4-v5; v4=X[ 8]; v5=X[24]; v2=v4+v5; t3=v4-v5; } t0=v0+v2; t1=v0-v2; }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[ 4]; v5=X[20]; v0=v4+v5; v1=v4-v5; v4=X[28]; v5=X[12]; v2=v4+v5; v3=v4-v5; } t4=v2+v0; t5=v2-v0; t6=SIN16*(v3+v1); t7=SIN16*(v3-v1); }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[30]; v5=X[14]; v0=v4+v5; v1=v4-v5; v4=X[ 6]; v5=X[22]; v2=v4+v5; v3=v4-v5; } t8=v0+v2; t9=v0-v2; t10=SIN24*v1+SIN08*v3; t11=SIN08*v1-SIN24*v3; }
  { value_type v0,v1,v2,v3; { value_type v4,v5; v4=X[ 2]; v5=X[18]; v0=v4+v5; v1=v4-v5; v4=X[10]; v5=X[26]; v2=v4+v5; v3=v4-v5; } t12=v0+v2; t13=v0-v2; t14=SIN24*v1-SIN08*v3; t15=SIN08*v1+SIN24*v3; }
  { value_type v0,v1,v2,v3,v4,v5; { value_type u0,u1; u0=X[31]; u1=X[15]; v0=u0-u1; v1=u0+u1; u0=X[ 7]; u1=X[23]; v2=u0-u1; v3=u0+u1; } t16=v1+v3; t19=v1-v3; { value_type u0,u1; u0=X[ 3]; u1=X[19]; v4=u0-u1; v5=u0+u1; u0=X[27]; u1=X[11]; v1=u0-u1; v3=u0+u1; } t17=v3+v5; t18=v3-v5; v3=SIN16*(v1-v4); v5=SIN16*(v1+v4); t20=v3-v2; t21=v3+v2;  t22=v0+v5; t23=v0-v5; }
  { value_type v0,v1,v2,v3,v4,v5; { value_type u0,u1; u0=X[ 1]; u1=X[17]; v0=u0-u1; v1=u0+u1; u0=X[ 9]; u1=X[25]; v2=u0-u1; v3=u0+u1; } t24=v1+v3; t27=v1-v3; { value_type u0,u1; u0=X[ 5]; u1=X[21]; v4=u0-u1; v5=u0+u1; u0=X[29]; u1=X[13]; v1=u0-u1; v3=u0+u1; } t25=v3+v5; t26=v3-v5; v3=SIN16*(v1-v4); v5=SIN16*(v1+v4); t28=v3-v2; t29=v3+v2; t30=v0+v5; t31=v0-v5; }

  {
    value_type v16,v17,v18,v19; v16=t24-t25; v17=t16-t17; v18=SIN16*(v17+v16); v19=SIN16*(v17-v16); v16=t0-t4; v17=t8-t12;
    value_type v8,v9,v10,v11,v12,v13,v14,v15; { value_type u0,u1; u0=t2-t6; u1=t15+t11; v8=u0+u1; v9=u0-u1; u0=SIN20*t29-SIN12*t31; u1=SIN12*t23+SIN20*t21; v10=u1+u0; v11=u1-u0; } { value_type u0,u1; u0=SIN20*t31+SIN12*t29; u1=SIN20*t23-SIN12*t21; v12=u1+u0; v13=u1-u0; u0=t10-t14; u1=t3+t7; v14=u0-u1; v15=u0+u1; }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v16+v18,v19+v17,value_type(SIN32,SIN28,SIN24, SIN20),value_type(SIN00,-SIN04,-SIN08,-SIN12),y0,z0,  v9+v11,v13+v14,value_type(SIN32,SIN27,SIN22, SIN17),value_type(SIN00,-SIN05,-SIN10,-SIN15),y1,z1); Y.first( 4,y0); Y.second(100,z0); Y.first(36,y1); Y.second( 68,z1); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v16-v18,v19-v17,value_type(SIN32,SIN20,SIN08,-SIN04),value_type(SIN00,-SIN12,-SIN24,-SIN28),y0,z0,  v8-v12,v10-v15,value_type(SIN32,SIN19,SIN06,-SIN07),value_type(SIN00,-SIN13,-SIN26,-SIN25),y1,z1); Y.first(12,y0); Y.second(108,z0); Y.first(44,y1); Y.second( 76,z1); }
    value_type v0,v1,v2,v3,v4,v5,v6,v7; v1=SIN24*t27+SIN08*t26; v5=SIN24*t19-SIN08*t18; v4=v5+v1; v5=v5-v1; v1=SIN24*t26-SIN08*t27; v7=SIN08*t19+SIN24*t18; v6=v7+v1; v7=v7-v1; v1=SIN16*(t9+t13); v0=t1+v1; v1=t1-v1; v3=SIN16*(t9-t13); v2=v3-t5; v3=v3+t5;
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v0 +v4 ,v6 +v3 ,value_type(SIN32,SIN30,SIN28,SIN26),value_type(SIN00,-SIN02,-SIN04,-SIN06),y0,z0,  v8+v12,v10+v15,value_type(SIN32,SIN29,SIN26, SIN23),value_type(SIN00,-SIN03,-SIN06,-SIN09),y1,z1); Y.first( 2,y0); Y.second( 98,z0); Y.first(34,y1); Y.second( 66,z1); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v1 -v7 ,v5 -v2 ,value_type(SIN32,SIN22,SIN12, SIN02),value_type(SIN00,-SIN10,-SIN20,-SIN30),y0,z0,  v9-v11,v13-v14,value_type(SIN32,SIN21,SIN10,-SIN01),value_type(SIN00,-SIN11,-SIN22,-SIN31),y1,z1); Y.first(10,y0); Y.second(106,z0); Y.first(42,y1); Y.second( 74,z1); }
    value_type v26,v27,v28,v29,v30,v31,v32,v33; { value_type u0,u1; u0=t2+t6; u1=t14+t10; v26=u0+u1; v27=u0-u1; u0=SIN28*t28-SIN04*t30; u1=SIN04*t22+SIN28*t20; v28=u1+u0;  v29=u1-u0; } { value_type u0,u1; u0=SIN28*t30+SIN04*t28; u1=SIN28*t22-SIN04*t20; v30=u1+u0; v31=u1-u0; u0=t11-t15; u1=t7-t3; v32=u0-u1; v33=u0+u1; }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v1 +v7 ,v5 +v2 ,value_type(SIN32,SIN26,SIN20,SIN14),value_type(SIN00,-SIN06,-SIN12,-SIN18),y0,z0,  v27+v29,v31+v32,value_type(SIN32,SIN25,SIN18,SIN11),value_type(SIN00,-SIN07,-SIN14,-SIN21),y1,z1); Y.first( 6,y0); Y.second(102,z0); Y.first(38,y1); Y.second( 70,z1); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v0 -v4 ,v6 -v3 ,value_type(SIN32,SIN18,SIN04,-SIN10),value_type(SIN00,-SIN14,-SIN28,-SIN22),y0,z0,  v26-v30,v28-v33,value_type(SIN32,SIN17,SIN02,-SIN13),value_type(SIN00,-SIN15,-SIN30,-SIN19),y1,z1); Y.first(14,y0); Y.second(110,z0); Y.first(46,y1); Y.second( 78,z1); }
    value_type v20,v21,v22,v23,v24,v25; v20=t0+t4; v21=t12+t8; v22=v20+v21; v23=t24+t25; v24=t16+t17; v25=v23+v24;
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v22+v25,                                                                                                  y0,z0,  v26+v30,v28+v33,value_type(SIN32,SIN31,SIN30,SIN29),value_type(SIN00,-SIN01,-SIN02,-SIN03),y1,z1); Y.first0(  y0); Y.secondN(   z1); Y.first(33,y1[1]); Y.second( 96,z0); }
    { complex_value_type y0,z0,y1,z1; real_crossed_fft_function<0,4>()(v20-v21,v24-v23,value_type(SIN32,SIN24,SIN16,SIN08),value_type(SIN00,-SIN08,-SIN16,-SIN24),y0,z0,  v27-v29,v31-v32,value_type(SIN32,SIN23,SIN14,SIN05),value_type(SIN00,-SIN09,-SIN18,-SIN27),y1,z1); Y.first( 8,y0); Y.second(104,z0); Y.first(40,y1   ); Y.second( 72,z1); }
    { complex_value_type y0; real_crossed_fft_function<0,4>()(v22-v25,value_type(SIN32, SIN16, SIN00, -SIN16),value_type(SIN00,-SIN16,-SIN32,-SIN16),y0); Y.first( 16,y0[0]); Y.first(48,y0[1]); }
  }
}

template<> struct real_crossed_fft_function< 2,2> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_2_2 (X,typename F::template rebind<  4>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function< 4,2> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_4_2 (X,typename F::template rebind<  8>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function< 8,2> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_8_2 (X,typename F::template rebind< 16>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function<16,2> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_16_2(X,typename F::template rebind< 32>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function<32,2> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_32_2(X,typename F::template rebind< 64>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function<64,2> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_64_2(X,typename F::template rebind<128>::other(f).real_make(Y)); } };

template<> struct real_crossed_fft_function< 2,4> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_2_4 (X,typename F::template rebind<  8>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function< 4,4> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_4_4 (X,typename F::template rebind< 16>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function< 8,4> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_8_4 (X,typename F::template rebind< 32>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function<16,4> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_16_4(X,typename F::template rebind< 64>::other(f).real_make(Y)); } };
template<> struct real_crossed_fft_function<32,4> { template<class G1, class G3,class F> inline void operator()(const Vector<G1> &X, Vector<G3> &Y, const F &f) { real_fft_32_4(X,typename F::template rebind<128>::other(f).real_make(Y)); } };

template<int M> struct real_complex_fft_function
{
  template<class G1,         class G3,class F> inline void operator()(const Vector<G1> &X,                      Vector<G3> &Y                 ,const F &f) { complex_fft_base_function<M>()(X  ,typename F::template rebind<M>::other(f).     make(Y    )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W, Vector<G3> &Y                 ,const F &f) { complex_fft_base_function<M>()(X,W,typename F::template rebind<M>::other(f).     make(Y    )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X, const Vector<G2> &W, Vector<G3> &Y0, Vector<G3> &Y1,const F &f) { complex_fft_base_function<M>()(X,W,typename F::template rebind<M>::other(f).pair_make(Y0,Y1)); }
};


template<int M,int N> struct real_complex_crossed_fft_base_function {};
template<> struct real_complex_crossed_fft_base_function< 1,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_1_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_1_2 (Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function< 2,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_2_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_2_2 (Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function< 3,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_3_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_3_2 (Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function< 4,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_4_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_4_2 (Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function< 5,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_5_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_5_2 (Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function< 8,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_8_2 (X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_8_2 (Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function<16,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_16_2(X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_16_2(Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function<32,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_32_2(X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_32_2(Xa,Xb,f); } };
template<> struct real_complex_crossed_fft_base_function<64,2> { template<class G1,class F> inline void operator()(const Vector<G1> &X, const F &f) { complex_fft_64_2(X,f); } template<class G1,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G1> &Xb, const F &f) { complex_fft_64_2(Xa,Xb,f); } };

template<int M,int N> struct real_complex_crossed_fft_function
{
  template<class G1,         class G3,class F> inline void operator()(const Vector<G1> &X ,                       Vector<G3> &Y                  , const F &f) { real_complex_crossed_fft_base_function<M,N>()(X  ,typename F::template rebind<M*N>::other(f).make     (Y    )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X , const Vector<G2> &W , Vector<G3> &Y                  , const F &f) { real_complex_crossed_fft_base_function<M,N>()(X*W,typename F::template rebind<M*N>::other(f).make     (Y    )); }
  template<class G1,         class G3,class F> inline void operator()(const Vector<G1> &X ,                       Vector<G3> &Y0 , Vector<G3> &Y1, const F &f) { real_complex_crossed_fft_base_function<M,N>()(X  ,typename F::template rebind<M*N>::other(f).pair_make(Y0,Y1)); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &X , const Vector<G2> &W , Vector<G3> &Y0 , Vector<G3> &Y1, const F &f) { real_complex_crossed_fft_base_function<M,N>()(X*W,typename F::template rebind<M*N>::other(f).pair_make(Y0,Y1)); }

  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G2> &Wa, Vector<G3> &Ya ,                  const Vector<G1> &Xb, const Vector<G2> &Wb, Vector<G3> &Yb                  , const F &f) { real_complex_crossed_fft_base_function<M,N>()(Xa*Wa, Xb*Wb, typename F::template rebind<M*N>::other(f).simd_make     (Ya     )); }
  template<class G1,class G2,class G3,class F> inline void operator()(const Vector<G1> &Xa, const Vector<G2> &Wa, Vector<G3> &Ya0, Vector<G3> &Ya1, const Vector<G1> &Xb, const Vector<G2> &Wb, Vector<G3> &Yb0, Vector<G3> &Yb1, const F &f) { real_complex_crossed_fft_base_function<M,N>()(Xa*Wa, Xb*Wb, typename F::template rebind<M*N>::other(f).simd_pair_make(Ya0,Ya1)); }
};



template<int p, class G1,class G2>
void real_mixed_fft_step0(const Vector<G1> &X, Vector<G2> &Y, int m, int m2, int i0=0, int is=1)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  int q = X.size()/p;
  int n = X.size()/m;

  real_fft_function<p> fft;

  for (int i=i0; i<n; i+=is)
  {
    typename SubArray<Vector<G2> >::self Ya=sub(Y,i*m2,m2);
    fft(stride(X,q,i),Ya);
  }
}

template<int p,class G1,class G,class G2,class F>
void real_mixed_fft_step1(const Vector<G1> &X, const Matrix<G> &W, Vector<G2> &Y, int m, int m2, const F &f, int i0, int is, int imax)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  int n = X.size()/m;
  int l=m/p, l2=m2/p;
  int q=l*n, q2=l2*n;

  real_complex_fft_function<p> fft;

  for (int i=i0; i<imax; i+=is)
  {
    int li=l2*i, mi=m2*i;
    for (int j=l/2; j>=0; --j)
    {
      typename StrideArray<Vector<G2> >::self Y0=stride(Y,l,mi+j);
      typename StrideArray<Vector<G2> >::self Y1=stride(Y,l,mi+(l-j));
      fft(stride(X,q2,li+j),row(W,j),Y0,Y1,f);
    };
  }
}


template<int p,class G1,class G2,class F>
inline void real_mixed_fft_step1(const Vector<G1> &X, Vector<G2> &Y,int m, int m2, const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  int n=X.size()/m;
  int q=m/p;
  real_type x=1;

  static __thread struct tls_twiddle { int m,n; real_type x; void *z; } twiddles[3*FFT_TWIDDLES]={ {0,0,0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.m==m && w.n==n && w.x==x) goto found; for (int i=1; i<3*FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.m==m && w.n==n && w.x==x) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.m=m; w.n=n; w.x=x; w.z=aligned_malloc((q/2+1)*p*sizeof(value_type)); twiddles[0]=w; typename DataMatrix<value_type>::self W(q/2+1,p, (value_type *)w.z); fft_twiddles(W,x,q/2+1,p,m); }
  found: typename DataMatrix<value_type>::self W(q/2+1,p, (value_type *)w.z);

  real_mixed_fft_step1<p>(X,W,Y,m,m2,f, 0,1,n);
}


template<int N,int p,class G1,class G2,class F>
void real_mixed_fft_step2(const Vector<G1> &X, Vector<G2> &Y, int m, int m2,const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  int n=X.size()/m;
  int q=m/p, q2=m2/p;
  real_type x=1;

  static __thread struct tls_twiddle { int m,n; real_type x; void *z; } twiddles[FFT_TWIDDLES]={ {0,0,0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.m==m && w.n==n && w.x==x) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.m==m && w.n==n && w.x==x) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.m=m; w.n=n; w.x=x; w.z=aligned_malloc(q*p*sizeof(value_type)); twiddles[0]=w; typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z); fft_twiddles2(W,x,q,p,m,N); }
  found: typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z);

  real_complex_crossed_fft_function<p,N> fft;

  for (int i=q/2; i>=0; --i)
  {
    typename StrideArray<Vector<G2> >::self Y0 = stride(Y,q,i);
    typename StrideArray<Vector<G2> >::self Y1 = stride(Y,q,q-i);
    fft(stride(X,q2,i),row(W,i),Y0,Y1,f);
  }
}

template<int N,int p,class G1,class G2,class F>
void real_mixed_fft_dense_step2(const Vector<G1> &X, Vector<G2> &Y,int m, int m2, const F &f)
{
  typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  int n=X.size()/m;
  int q=m/p, q2=m2/p;
  real_type x=1;

  static __thread struct tls_twiddle { int m,n; real_type x; void *z; } twiddles[FFT_TWIDDLES]={ {0,0,0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.m==m && w.n==n && w.x==x) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.m==m && w.n==n && w.x==x) { twiddles[0]=w; goto found; } }
  { aligned_free(w.z); w.m=m; w.n=n; w.x=x; w.z=aligned_malloc(q*p*sizeof(value_type)); twiddles[0]=w; typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z); fft_twiddles2(W,x,q,p,m,N); }
  found: typename DataMatrix<value_type>::self W(q,p, (value_type *)w.z);

  real_complex_crossed_fft_function<p,N> fft;

  for (int i=q/2; i>=0; i-=2)
  {
    typename StrideArray<Vector<G2> >::self Ya0 = stride(Y,q,i);
    typename StrideArray<Vector<G2> >::self Ya1 = stride(Y,q,q-i);
    typename StrideArray<Vector<G2> >::self Yb0 = stride(Y,q,i+1);
    typename StrideArray<Vector<G2> >::self Yb1 = stride(Y,q,q-i-1);
    fft(stride(X,q2,i),row(W,i),Ya0,Ya1, stride(X,q2,i+1),row(W,i+1),Yb0,Yb1, f);
  }
}


template<int N>
struct real_mixed_fft_step0_function
{
  template<class G1,class G2> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p, int m, int m2)
  {
    switch (p)
    {
#if FFT_LEVEL>=2
      case  2: real_mixed_fft_step0< 2>(X,Y,m,m2); break;
#endif
#if FFT_LEVEL>=4
      case  4: real_mixed_fft_step0< 4>(X,Y,m,m2); break;
#endif
#if FFT_LEVEL>=8
      case  8: real_mixed_fft_step0< 8>(X,Y,m,m2); break;
#endif
#if FFT_LEVEL>=16
      case 16: real_mixed_fft_step0<16>(X,Y,m,m2); break;
#endif
#if FFT_LEVEL>=32
      case 32: real_mixed_fft_step0<32>(X,Y,m,m2); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: real_mixed_fft_step0< 3>(X,Y,m,m2); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: real_mixed_fft_step0< 5>(X,Y,m,m2); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }
};

struct real_mixed_fft_step1_function
{
  template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p, int m, int m2, const F &f)
  {
    switch (p)
    {
#if FFT_LEVEL>=2
      case  2: real_mixed_fft_step1< 2>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: real_mixed_fft_step1< 4>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: real_mixed_fft_step1< 8>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL >= 16
      case 16: real_mixed_fft_step1<16>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL >= 32
      case 32: real_mixed_fft_step1<32>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: real_mixed_fft_step1< 3>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: real_mixed_fft_step1< 5>(X,Y,m,m2,f); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }
  template<class G1,class G2> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p, int m,int m2)
  {
    (*this)(X,Y,p,m,m2,make_fft_adaptor());
  }
};


template<int N>
struct real_mixed_fft_step2_function
{
  template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p,int m, int m2,const F &f)
  {
    switch (p)
    {
      case  1: real_mixed_fft_step2<N, 1>(X,Y,m,m2,f); break;
#if FFT_LEVEL>=2
      case  2: real_mixed_fft_step2<N, 2>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: real_mixed_fft_step2<N, 4>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: real_mixed_fft_step2<N, 8>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=16
      case 16: real_mixed_fft_step2<N,16>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=32
      case 32: real_mixed_fft_step2<N,32>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: real_mixed_fft_step2<N, 3>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: real_mixed_fft_step2<N, 5>(X,Y,m,m2,f); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }

#ifdef NDEBUG
  template<class G1,class V,class F> inline void operator()(const Vector<G1> &X, Vector<data_vector_generator<V> > &Y, int p,int m,int m2,const F &f)
  {
    switch (p)
    {
      case  1: real_mixed_fft_dense_step2<N, 1>(X,Y,m,m2,f); break;
#if FFT_LEVEL>=2
      case  2: real_mixed_fft_dense_step2<N, 2>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=4
      case  4: real_mixed_fft_dense_step2<N, 4>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=8
      case  8: real_mixed_fft_dense_step2<N, 8>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=16
      case 16: real_mixed_fft_dense_step2<N,16>(X,Y,m,m2,f); break;
#endif
#if FFT_LEVEL>=32
      case 32: real_mixed_fft_dense_step2<N,32>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
      case  3: real_mixed_fft_dense_step2<N, 3>(X,Y,m,m2,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
      case  5: real_mixed_fft_dense_step2<N, 5>(X,Y,m,m2,f); break;
#endif
      default: cerr << "fft error: factor unknown" << endl;
    }
  }
#endif
};

template<> struct real_mixed_fft_step2_function<0> { template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p,int m,int m2,const F &f) { real_mixed_fft_step1_function()(X,Y,p,m,m2,f); } };
template<> struct real_mixed_fft_step2_function<1> { template<class G1,class G2,class F> inline void operator()(const Vector<G1> &X, Vector<G2> &Y, int p,int m,int m2,const F &f) { real_mixed_fft_step1_function()(X,Y,p,m,m2,f); } };


inline void real_fft_factors(int n, int n2, int *X)
{
#if FFT_LEVEL>=16
  switch (n2)
  {
    case   128: X[0]=2; X[1]=16; X[2]= 8; return;
    case   256: X[0]=2; X[1]=16; X[2]=16; return;
    case  2048: X[0]=3; X[1]=16; X[2]=16; X[3]= 8; return;
    case  4096: X[0]=3; X[1]=16; X[2]=16; X[3]=16; return;
    case 32768: X[0]=4; X[1]=16; X[2]=16; X[3]=16; X[4]= 8; return;
    case 65536: X[0]=4; X[1]=16; X[2]=16; X[3]=16; X[4]=16; return;
  };
#endif

#if (FFT_LEVEL>=32)
  if (n%32==0) { X[1]=32; n=n2/32; } else
#endif
#if (FFT_LEVEL>=16)
  if (n%16==0) { X[1]=16; n=n2/16; } else
#endif
#if (FFT_LEVEL>= 8)
  if (n% 8==0) { X[1]= 8; n=n2/ 8; } else
#endif
#if (FFT_LEVEL>= 4)
  if (n% 4==0) { X[1]= 4; n=n2/ 4; } else
#endif
#if (FFT_LEVEL>= 2)
  if (n% 2==0) { X[1]= 2; n=n2/ 2; } else
#endif
#if (FFT_LEVEL>= 5) && !defined(FFT_ONLY_POW2)
  if (n% 5==0) { X[1]= 5; n=n2/ 5; } else
#endif
#if (FFT_LEVEL>= 3) && !defined(FFT_ONLY_POW2)
  if (n% 3==0) { X[1]= 3; n=n2/ 3; } else
#endif
  cerr << "fft factor unknown" << endl;

  int k=1;
#if (FFT_LEVEL>=32)
  for (; n%32==0; n/=32) X[++k]=32;
#endif
#if (FFT_LEVEL>=16)
  for (; n%16==0; n/=16) X[++k]=16;
#endif
#if (FFT_LEVEL>= 8)
  for (; n% 8==0; n/= 8) X[++k]= 8;
#endif
#if (FFT_LEVEL>= 4)
  for (; n% 4==0; n/= 4) X[++k]= 4;
#endif
#if (FFT_LEVEL>= 2)
  for (; n% 2==0; n/= 2) X[++k]= 2;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
  for (; n% 5==0; n/= 5) X[++k]= 5;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
  for (; n% 3==0; n/= 3) X[++k]= 3;
#endif
  if (n!=1) { cerr << "fft factor unknown" << endl; }
  X[0]=k;
}

template<int N,class G1,class G2,class G3,class F>
void real_mixed_fft(const Vector<G1> &X, Vector<G2> &T, Vector<G2> &U, Vector<G3> &Y, const F &f)
{
  int n=T.size();
  static __thread struct { int n; int factors[50]; } decomp = {0,0};
  if (n!=decomp.n) { decomp.n=n; real_fft_factors(X.size(),n, decomp.factors); }

  int o=decomp.factors[0];

  int p=decomp.factors[1];
  int M=N+(N==0);
  int m=p, m2=M*(p/(2*M)+1);

  if (o%2==1) swap(U,T);
  real_mixed_fft_step0_function<N>()(X,T,p,m,((N!=0)+1)*m2);
  for (int k=2; k<o; ++k)
  {
    p=decomp.factors[k];
    real_mixed_fft_step1_function()(T,U,p,m*=p,m2*=p);
    swap(T,U);
  }
  p=decomp.factors[o];
  real_mixed_fft_step2_function<N>()(T,Y,p,m*=p,m2*=p,f);
}


template<int N>
struct real_mixed_fft_function
{
  template<class G1,class V,class G2,class F>
  inline void operator()(const Vector<G1> &X, Vector<data_vector_generator<complex<V> > > &T, Vector<data_vector_generator<complex<V> > > &U, Vector<G2> &Y, const F &f)
  {
    typedef complex<V> value_type;
    typedef typename TinyVector<N,value_type>::self elem_type;
    typename tinyBlockVector<const Vector<G1>,2*N>::self X2 = block<2*N>(X);
    typename DataVector<elem_type>::self T2(T.size()/N, (elem_type *)&*T.begin());
    typename DataVector<elem_type>::self U2(U.size()/N, (elem_type *)&*U.begin());
    real_mixed_fft<N>(X2,T2,U2,Y,f);
  }

  template<class V,class G2,class F>
  inline void operator()(const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<complex<V> > > &T, Vector<data_vector_generator<complex<V> > > &U, Vector<G2> &Y, const F &f)
  {
    typedef complex<V> value_type;
    typedef typename TinyVector<2*N,V>::self real_elem_type;
    typedef typename TinyVector<N,value_type>::self elem_type;
    typename DataVector<real_elem_type>::self X2(X.size()/(2*N), (real_elem_type *)&*X.begin());
    typename DataVector<elem_type>::self T2(T.size()/N, (elem_type *)&*T.begin());
    typename DataVector<elem_type>::self U2(U.size()/N, (elem_type *)&*U.begin());
    real_mixed_fft<N>(X2,T2,U2,Y,f);
  }


  template<class G1,class G2,class F>
  inline void operator()(const Vector<G1> &X, Vector<G2> &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
    (*this)(X,T2,U2,Y,f);
  }

  template<class V,class F>
  inline void operator()(const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<complex<V> > > &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,V) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);

    if (Y.size()>=X.size())
    {
      typename DataVector<value_type>::self U2(X.size(), (value_type*)&*Y.begin());
      (*this)(X,T2,U2,Y,f);
    }
    else
    {
      //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
      if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
      (*this)(X,T2,U2,Y,f);
    }

  }
};

template<>
struct real_mixed_fft_function<0>
{
  template<class G1,class V,class G2,class F>
  inline void operator()(const Vector<G1> &X, Vector<data_vector_generator<complex<V> > > &T, Vector<data_vector_generator<complex<V> > > &U, Vector<G2> &Y, const F &f)
  {
    real_mixed_fft<0>(X,T,U,Y,f);
  }


  template<class G1,class G2,class F>
  inline void operator()(const Vector<G1> &X, Vector<G2> &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Vector<G1>::value_type) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
    (*this)(X,T2,U2,Y,f);
  }

  template<class V,class F>
  inline void operator()(const Vector<data_vector_generator<V> > &X, Vector<data_vector_generator<complex<V> > > &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,V) value_type;
    //typename DenseVector<value_type>::self T(X.size()); typename DataVector<value_type>::self T2(X.size(), (value_type*)&*T.begin());
    if (fft_buf0.n<int(sizeof(value_type)*X.size())) { fft_buf0.n=sizeof(value_type)*X.size(); aligned_free(fft_buf0.p); fft_buf0.p=aligned_malloc(fft_buf0.n); } typename DataVector<value_type>::self T2(X.size(), (value_type *)fft_buf0.p);
    if (Y.size()>=X.size())
    {
      typename DataVector<value_type>::self U2(X.size(), (value_type*)&*Y.begin());
      (*this)(X,T2,U2,Y,f);
    }
    else
    {
      //typename DenseVector<value_type>::self U(X.size()); typename DataVector<value_type>::self U2(X.size(), (value_type*)&*U.begin());
      if (fft_buf1.n<int(sizeof(value_type)*X.size())) { fft_buf1.n=sizeof(value_type)*X.size(); aligned_free(fft_buf1.p); fft_buf1.p=aligned_malloc(fft_buf1.n); } typename DataVector<value_type>::self U2(X.size(), (value_type *)fft_buf1.p);
      (*this)(X,T2,U2,Y,f);
    }
  }
};


template<class G1,class G2,class F> inline void real_mixed_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f) { real_mixed_fft_function<0>()(X,Y,f); }
template<class G1,class G2> inline void real_mixed_fft(const Vector<G1> &X, Vector<G2> &Y) { return real_mixed_fft(X,Y,make_fft_adaptor()); }
template<class G> inline PROMOTE2(complex<float>,Vector<G>) real_mixed_fft     (const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) T(X.size()); real_mixed_fft(X,T); return T; }

template<int N,class G1,class G2,class F> inline void real_mixed_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f) { real_mixed_fft_function<N>()(X,Y,f); }
template<int N,class G1,class G2> inline void real_mixed_fft(const Vector<G1> &X, Vector<G2> &Y) { return real_mixed_fft<N>(X,Y,make_fft_adaptor()); }
template<int N,class G> inline PROMOTE2(complex<float>,Vector<G>) real_mixed_fft(const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) T(X.size()); real_mixed_fft<N>(X,T); return T; }


template<class G1,class G2,class F,class V>
void aux_real_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, const V &)
{
  switch (X.size())
  {
    case  1: real_fft_function< 1>()(X,Y,f); break;
#if FFT_LEVEL>=2
    case  2: real_fft_function< 2>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=4
    case  4: real_fft_function< 4>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=8
    case  8: real_fft_function< 8>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=16
    case 16: real_fft_function<16>()(X,Y,f); break;
#endif
#if FFT_LEVEL>=32
    case 32: real_fft_function<32>()(X,Y,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
    case  3: real_fft_function< 3>()(X,Y,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
    case  5: real_fft_function< 5>()(X,Y,f); break;
#endif
    default: real_mixed_fft(X,Y,f);
  }
}

#ifdef SSE
template<class G1,class G2,class F>
void aux_real_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, float)
{
  int n = X.size();
  switch (n)
  {
    case   1: real_fft_function< 1>()(X,Y,f); break;
#if FFT_LEVEL>=2
    case   2: real_fft_function< 2>()(X,Y,f); break;
    case   4: real_fft_function< 4>()(X,Y,f); break;
    case   8: real_crossed_fft_function< 2,4>()(block<4>(X),Y,f); break;
#endif
#if FFT_LEVEL>=4
    case  16: real_crossed_fft_function< 4,4>()(block<4>(X),Y,f); break;
#endif
#if FFT_LEVEL>=8
    case  32: real_crossed_fft_function< 8,4>()(block<4>(X),Y,f); break;
#endif
#if FFT_LEVEL>=16
    case  64: real_crossed_fft_function<16,4>()(block<4>(X),Y,f); break;
#endif
#if FFT_LEVEL>=32
    case 128: real_crossed_fft_function<32,4>()(block<4>(X),Y,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
    case   3: real_fft_function< 3>()(X,Y,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
    case   5: real_fft_function< 5>()(X,Y,f); break;
#endif
    default:
      if (n%8==0) real_mixed_fft<2>(X,Y,f);
#if !defined(FFT_ONLY_POW2)
      else real_mixed_fft(X,Y,f);
#else
      else cerr << "fft error: factor unknown" << endl;
#endif
  }
}
#endif

#ifdef SSE2
template<class G1,class G2,class F>
void aux_real_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, double)
{
  int n = X.size();
  switch (n)
  {
#if FFT_LEVEL>=2
    case   2: real_fft_function< 2>()(X,Y,f); break;
    case   4: real_crossed_fft_function< 2,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=4
    case   8: real_crossed_fft_function< 4,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=8
    case  16: real_crossed_fft_function< 8,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=16
    case  32: real_crossed_fft_function<16,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=32
    case  64: real_crossed_fft_function<32,2>()(block<2>(X),Y,f); break;
#endif
#if FFT_LEVEL>=64
    case 128: real_crossed_fft_function<64,2>()(block<2>(X),Y,f); break;
#endif
#if (FFT_LEVEL>=3) && !defined(FFT_ONLY_POW2)
    case   3: real_fft_function< 3>()(X,Y,f); break;
#endif
#if (FFT_LEVEL>=5) && !defined(FFT_ONLY_POW2)
    case   5: real_fft_function< 5>()(X,Y,f); break;
#endif
    default:
      if (n%2==0)  real_mixed_fft<1>(X,Y,f);
#if !defined(FFT_ONLY_POW2)
      else real_mixed_fft(X,Y,f);
#else
      else cerr << "fft error: factor unknown" << endl;
#endif
  }
}
#endif

//template<class G1,class G2,class F>
//void shift_dense_real_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f)
//{
//  typedef Vector<G1> first_array_type;
//  typedef Vector<G2> second_array_type;
//  typedef typename first_array_type::value_type first_value_type;
//  typedef typename second_array_type::value_type second_value_type;
//  typedef PROMOTE2(complex<float>,first_array_type) array_type;
//  typedef typename array_type::value_type value_type;
//
//  int i0 = X.lower_bound();
//  int n = X.size();
//
//  Y.resize(n); Y.set_lower_bound(i0);
//
//  if (i0==0)
//  {
//    typename DenseVector<first_value_type>::self U;
//    typename DenseVector<second_value_type>::self V;
//
//    U.swap(const_cast<first_array_type &>(X));
//    Y.swap(V);
//    real_fft(U,V,f);
//    Y.swap(V);
//    U.swap(const_cast<first_array_type &>(X));
//  }
//  else
//  {
//    static array_type W(0);
//    if ( W.size()!=n || W.lower_bound()!=i0)
//    {
//      W.resize(n); W.set_lower_bound(i0);
//      value_type w0=polar(first_value_type(1),(-2*first_value_type(PI)*i0)/n);
//      value_type w=w0;
//      W[i0]=value_type(1);
//      W[i0+1]=w;
//      for (int i=i0+2, imax=W.upper_bound()+1; i<imax; ++i)
//        W[i]=(w*=w0);
//    }
//    Y = W*X;
//
//    typename DenseVector<value_type>::self V;
//
//    Y.swap(V);
//    complex_fft_in_place(V,f);
//    Y.swap(V);
//  }
//}


template<class G1,class G2,class F> void aux_real_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f) { aux_real_fft(X,Y,f,typename Vector<G1>::value_type()); }

//{secret}
template<class G1,class G2,class F>
inline void real_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f)
{
  //int i0 = X.lower_bound();
  //if (i0!=0)
  //  cerr << "error: fft of shifted vectors not impemented" << endl;
  //else
  {
    typename ArrayData<Vector<G2> >::self Y2 = data(Y);
    aux_real_fft(data(X),Y2,f);
  }
}


//{unsecret}
//Summary: FFT of a real signal
//Arguments: X - A real signal
//Output: Y - The half complex FFT
//Remarks:
//  The FFT of a signal with N real values is complex and has the following symmetry property:
//
//      Y[N-i] = conj(Y[i])
//
//  For some speed consideration, only the first half part of /Y/ is filled with correct values.
//  /Y/ should have at least N/2+1 complex elements.
//  Nevertheless /Y/ should best have N elements, /Y/ is then used as internal memory buffer resulting
//  in overwritten values in the second half part and a somewhat quicker execution speed.
//
//  In case of a 2D signal with MxN real values, only the left part of the /Y/ matrix is filled with correct values.
//  /Y/ should usualy have Mx(N/2+1) complex elements.
//  Nevertheless for an optimal execution on single precision signals with SSE instructions,
//  every line of /Y/ has to be aligned in memory, that is /Y/ should have Mx2*((N/2)/2+1) complex elements
//
//  If needed, the second half part of /Y/ can be reconstructed by the use of its symmetry property,
//  but the optimized /fft/ function that entirely fills /Y/ should be preferred.
//Example:
//  DenseVector<float>::self X(4, 1);
//  DenseVector<complex<float> >::self Y(X.size());
//  half_real_fft(X,Y);
//
//  DenseMatrix<float>::self X(4,4, 1);
//  DenseMatrix<complex<float> >::self Y(X.nrows(),2*((X.ncols()/2)/2+1);
//  half_real_fft(X,Y);
template<class G1,class G2> inline void half_real_fft(const Vector<G1> &X, Vector<G2> &Y) { real_fft(X,Y,make_fft_adaptor     ( )); }
template<class G1,class G2> inline void      real_fft(const Vector<G1> &X, Vector<G2> &Y) { real_fft(X,Y,make_full_fft_adaptor(Y)); }
template<class G> inline PROMOTE2(complex<float>,Vector<G>) real_fft(const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) Y(X.size()); real_fft     (X,Y); return Y; }


template<class G1,class G2> inline void      real_norm_fft(const Vector<G1> &X, Vector<G2> &Y) { real_fft(X,Y,make_norm_full_fft_adaptor(Y)); }
template<class G1,class G2> inline void half_real_norm_fft(const Vector<G1> &X, Vector<G2> &Y) { real_fft(X,Y,make_norm_fft_adaptor     ( )); }
template<class G> inline  typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self      real_norm_fft(const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size());      real_norm_fft(X,Y); return Y; }
template<class G> inline  typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self half_real_norm_fft(const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size()); half_real_norm_fft(X,Y); return Y; }

template<class G1,class G2> inline void      real_abs_fft(const Vector<G1> &X, Vector<G2> &Y) { real_fft(X,Y,make_abs_full_fft_adaptor(Y)); }
template<class G1,class G2> inline void half_real_abs_fft(const Vector<G1> &X, Vector<G2> &Y) { real_fft(X,Y,make_abs_fft_adaptor     ( )); }
template<class G> inline  typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self      real_abs_fft (const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size());      real_abs_fft (X,Y); return Y; }
template<class G> inline  typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self half_real_abs_fft (const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size()); half_real_abs_fft (X,Y); return Y; }



template<class G1,class G2,class V> void aux2_real_ifft(const Vector<G1> &X, Vector<G2> &T, Vector<data_vector_generator<V> > &Y)
{
  typedef typename Vector<G1>::value_type value_type;
  typedef typename value_type::value_type real_value_type;
  typedef typename value_traits<real_value_type   >::value_type real_type;
  typedef complex<real_type> complex_type;

  if (is_odd(Y.size())) cerr << "real_ifft error: algorithm only works with even lengths" << endl;
  int n=Y.size()/2, imax=(n+1)/2;
  real_type x=1/real_type(2*n);
  const int di=SimdBlockLength<value_type>::RET;

  static __thread struct tls_twiddle { int n; void *z; } twiddles[FFT_TWIDDLES]={ {0,NULL} };
  tls_twiddle w=twiddles[0]; if (w.n==n) goto found; for (int i=1; i<FFT_TWIDDLES; ++i) { tls_twiddle wi=twiddles[i]; twiddles[i]=w; w=wi; if (w.n==n) { twiddles[0]=w; goto found; } }
  {
    aligned_free(w.z); w.n=n; 
    w.z=aligned_malloc((n/2+1)*sizeof(complex_type)); 
    twiddles[0]=w; typename DataVector<complex_type>::self W(n/2+1, (complex_type *)w.z); 
    W=x*generate_vector(n/2+1,scale(polar_function<real_type>(),PI/n)); 
    }
  found: typename DataVector<complex_type>::self W(n/2+1, (value_type *)w.z);

  {
    real_value_type a=real(X[0]), b=real(X[n]);
    T[0]=x*value_type(a+b,a-b);
  }
  if (imax<4*di)
    for (int i=1; i<imax; ++i)
    {
      value_type x0=X[i], x1=X[n-i];
      value_type a=x*addsub(x0,x1), b=subadd(x0,x1);
      b=imul(b,W[i]);
      T[i]=conj(a-b);
      T[n-i]=a+b;
    }
  else
  {
    for (int i=1, dimax=__min(di,imax); i<dimax; ++i)
    {
      value_type x0=X[i], x1=X[n-i];
      value_type a=x*addsub(x0,x1), b=subadd(x0,x1);
      b=imul(b,W[i]);
      T[i]=conj(a-b);
      T[n-i]=a+b;
    }
    {
      typename SimdBlock<typename DataVector<complex_type>::self>::self W2=simd_block(W);
      typename SimdBlock<const Vector<G1> >::self X2=simd_block(X); typename SubArray <const Vector<G1> >::self X3=sub(X,1,X.size()-1); typename SimdUnalignedBlock<typename SubArray<const Vector<G1> >::self>::self X4=simd_unaligned_block(X3);
      typename SimdBlock<      Vector<G2> >::self T2=simd_block(T); typename SubArray <      Vector<G2> >::self T3=sub(T,1,T.size()-1); typename SimdUnalignedBlock<typename SubArray<      Vector<G2> >::self>::self T4=simd_unaligned_block(T3);
      int m=n/di;
      for (int i=1, imax=m/2; i<imax; ++i)
      {
        typename SimdBlock<const Vector<G1> >::value_type x0,x1,a,b;
        x0=X2[i]; x1=X4[m-1-i]; x1=flip(x1);
        a=x*addsub(x0,x1); b=subadd(x0,x1);
        b=imul(b,W2[i]);
        T2[i]=conj(a-b);
        T4[m-1-i]=flip(a+b);
      }
    }
    for (int i=(imax/di)*di; i<imax; ++i)
    {
      value_type x0=X[i], x1=X[n-i];
      value_type a=x*addsub(x0,x1), b=subadd(x0,x1);
      b=imul(b,W[i]);
      T[i]=conj(a-b);
      T[n-i]=a+b;
    }
  }
  if (is_even(n)) T[n/2]=x*real_type(2)*conj(X[n/2]);

  typename DataVector <value_type>::self Y2(n,(value_type *)&Y.front());
  complex_fft(T,Y2);
}

template<class G1,class G2,class G3> inline void aux2_real_ifft(const Vector<G1> &X, Vector<G2> &T, Vector<G3> &Y)
{
  typename DenseVector<typename Vector<G3>::value_type>::self U(Y.size());
  real_ifft(X,T,U);
  Y=U;
}

template<class G1,class G2,class G3> inline void real_ifft(const Vector<G1> &X, Vector<G2> &T, Vector<G3> &Y)
{
  typename ArrayData<Vector<G2> >::self T2 = data(T);
  typename ArrayData<Vector<G3> >::self Y2 = data(Y);
  aux2_real_ifft(data(X),T2,Y2);
}

template<class G1,class G2,class V> inline void aux_real_ifft(Vector<G1> &X, Vector<G2> &Y,const V &)
{
  typename SubArray<Vector<G1> >::self T=sub(X,0,Y.size()/2);
  real_ifft(X,T,Y);
}

template<class G1,class G2> inline void aux_real_ifft(Vector<G1> &X, Vector<G2> &Y)
{
  aux_real_ifft(X,Y,typename Vector<G1>::value_type());
}

//{unsecret}
//Summary: IFFT of a complex-symmetric signal
//Arguments: X - A complex or half-complex signal
//Output: Y - The real IFFT
//Remarks:
//  The complex input signal is supposed to correspond to the FFT of a real signal,
//  and therefore to respect the following symmetry property:
//
//          Y[N-i] = conj(Y[i])
//
//  Only the first values of the input signal (indexes 0 to N/2) are therefore considered.
//  For any signal of real values, the following statments equal X to within roundoff error:
//
//          real_ifft(fft(X))
//          real_ifft(half_real_fft(X))
//
//  The present algorithm has 2 limitations:
//
//    - It only works with input signal of even length.
//      For odd lengths, use the slowler /ifft/ function with a complex-symmetric input.
//
//          Y=real(ifft(X));
//
//    - In order to achieve better performances, it overwrites the vectotial input with arbitrary data.
//      (Matrix inputs are not overwritten).
//      Copy the input signal in a temporary array if you need to preserve it.
//
//          T=X;
//          Y=real_ifft(T);
//See: ^ifft^, ^fft^, ^half_real_fft^
template<class G1,class G2> inline void real_ifft(Vector<G1> &X, Vector<G2> &Y)
{
  typename ArrayData<Vector<G1> >::self X2 = data(X);
  typename ArrayData<Vector<G2> >::self Y2 = data(Y);
  aux_real_ifft(X2,Y2);
}
//{unsecret}
template<class G> inline typename DenseVector<typename Vector<G>::value_type::value_type >::self real_ifft(Vector<G> &X)
{
  typename DenseVector<typename Vector<G>::value_type::value_type>::self Y(X.size());
  real_ifft(X,Y);
  return Y;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

template<class G1,class G2,class F,class  V> inline void aux_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, const V          &) { real_fft    (X,Y,f); }
template<class G1,class G2,class F,class  V> inline void aux_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, const complex<V> &) { complex_fft (X,Y,f); }
template<class G1,class G2,class F,class G3> inline void aux_fft(const Vector<G1> &X, Vector<G2> &Y, const F &f, const Vector<G3> &) { aux_fft (X,Y,f,typename Vector<G3>::value_type()); }
//{secret}
template<class G1,class G2,class F> inline void fft(const Vector<G1> &X, Vector<G2> &Y, const F &f) { aux_fft (X,Y,f,typename Vector<G1>::value_type()); }
//{secret}
template<class G1,class G2        > inline void fft(Vector<G1> &X, Vector<G2> &Y, const make_ifft_adaptor &f) { ifft(X,Y); }


template<class G1,class G2,class  V> inline void aux_fft (const Vector<G1> &X, Vector<G2> &Y, const V          &) { real_fft    (X,Y); }
template<class G1,class G2,class  V> inline void aux_fft (const Vector<G1> &X, Vector<G2> &Y, const complex<V> &) { complex_fft (X,Y); }
template<class G1,class G2,class  G> inline void aux_fft (const Vector<G1> &X, Vector<G2> &Y, const Vector<G>  &) { aux_fft (X,Y,typename Vector<G>::value_type()); }
template<class G1,class G2,class  V> inline void aux_ifft(      Vector<G1> &X, Vector<G2> &Y, const V          &) { real_ifft   (X,Y); }
template<class G1,class G2,class  V> inline void aux_ifft(      Vector<G1> &X, Vector<G2> &Y, const complex<V> &) { complex_ifft(X,Y); }
template<class G1,class G2,class  G> inline void aux_ifft(      Vector<G1> &X, Vector<G2> &Y, const Vector<G>  &) { aux_ifft(X,Y,typename Vector<G>::value_type()); }


//{unsecret}
//Summary: FFT of a signal
//Arguments: X - The signal
//Output: Y - The complex FFT
//Return: Dense array of complex elements
//See: ^half_real_fft^
//Example:
//  DenseVector<complex<float> >::self X(4, 1);
//  DenseVector<complex<float> >::self Y(X.size());
//  fft(X,Y);
//
//  DenseVector<complex<float> >::self X(4, 1);
//  fft(X,X); // in place
//
//  DenseVector<double>::self X(4, 1);
//  DenseVector<complex<double> >::self Y;
//  Y=fft(X);
template<class G1,class G2> inline void fft(const Vector<G1> &X, Vector<G2> &Y) { aux_fft(X,Y,typename Vector<G1>::value_type()); }

//{unsecret}
//Summary: IFFT of a signal
//Arguments: X - The signal
//Output: Y - The complex IFFT
//Return: Dense array of complex elements
//Remarks:
//  If /Y/ contains real values, X is then supposed to be complex symetric and /real_ifft/ is performed (with its constraints).
//Example:
//  DenseVector<complex<double> >::self X(4, 1);
//  DenseVector<complex<double> >::self Y;
//  Y=ifft(fft(X));
//
//  DenseVector<float>::self X(4, 1);
//  DenseVector<complex<float> >::self Y(X.size()/2+1);
//  half_real_fft(X,Y); // quicker than fft(X,Y)
//  ifft(Y,X);
//See: ^fft^, ^real_ifft^
template<class G1,class G2> inline void ifft(Vector<G1> &X, Vector<G2> &Y) { aux_ifft(X,Y,typename Vector<G2>::value_type()); }

//{unsecret}
template<class G> inline PROMOTE2(complex<float>,Vector<G>) fft (const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) Y(X.size()); fft (X,Y); return Y; }
//{unsecret}
template<class G> inline PROMOTE2(complex<float>,Vector<G>) ifft(const Vector<G> &X) { PROMOTE2(complex<float>,Vector<G>) Y(X.size()); ifft(X,Y); return Y; }


template<class G1,class G2,class  V> inline void aux_norm_fft(const Vector<G1> &X, Vector<G2> &Y, const V          &) { real_norm_fft   (X,Y); }
template<class G1,class G2,class  V> inline void aux_norm_fft(const Vector<G1> &X, Vector<G2> &Y, const complex<V> &) { complex_norm_fft(X,Y); }
//{unsecret}
//Summary: FFT squared magnitude of a signal, energy spectrum
//Arguments: X - The signal
//Output:    Y - The FFT squared magnitude
//Return: Dense array of real elements
//See: ^fft^, ^abs_fft^
//Example:
//  DenseVector<complex<float> >::self X(4, 1);
//  DenseVector<float>::self Y;
//  Y=norm_fft(X);
template<class G1,class G2> inline void norm_fft(const Vector<G1> &X, Vector<G2> &Y) { aux_norm_fft(X,Y,typename Vector<G1>::value_type()); }
//{unsecret}
template<class G> inline typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self norm_fft(const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size()); norm_fft(X,Y); return Y; }


template<class G1,class G2,class  V> inline void aux_abs_fft (const Vector<G1> &X, Vector<G2> &Y, const V          &) { real_abs_fft    (X,Y); }
template<class G1,class G2,class  V> inline void aux_abs_fft (const Vector<G1> &X, Vector<G2> &Y, const complex<V> &) { complex_abs_fft (X,Y); }
//{unsecret}
//Summary: FFT magnitude of a signal, magnitude spectrum
//Arguments: X - The signal
//Output: Y - The FFT magnitude
//Return: Dense array of real elements
//See: ^fft^, ^norm_fft^
//Example:
//  DenseVector<float>::self X(4, 1);
//  DenseVector<float>::self Y(X.size());
//  abs_fft(X,Y);
template<class G1,class G2> inline void abs_fft (const Vector<G1> &X, Vector<G2> &Y) { aux_abs_fft(X,Y,typename Vector<G1>::value_type()); }
//{unsecret}
template<class G> inline typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self abs_fft (const Vector<G> &X) { typename DenseVector<typename complex_traits<typename Vector<G>::value_type>::value_type>::self Y(X.size()); abs_fft (X,Y); return Y; }


template<class G,class  V> inline void aux_fft_in_place(Vector<G> &X, const V          &) { real_fft(X,X,make_in_place_fft_adaptor()); }
template<class G,class  V> inline void aux_fft_in_place(Vector<G> &X, const complex<V> &) { complex_fft(X,X); }
//{unsecret}
//Summary: FFT in place of a signal
//Arguments: X - A signal
//Output: X - The FFT
//Remarks:
//  For complex signal, this function is equivalent to fft(X,X).
//
//  The FFT of a real signal with N terms is complex and has the following symmetry property:
//
//      Y[N-i] = conj(Y[i])
//
//  The arrangement of the compact half-complex terms uses the following scheme:
//
//    @untitled table
//    k = 0     X[k] = real(Y[k]), imaginary part is zero and not stored
//    k < N/2   X[k] = real(Y[k]), X[N-k] = imag(Y[k])
//    k = N/2   X[k] = real(Y[k]), imaginary part is zero and not stored
//    k > N/2   Terms not stored but can be reconstructed using the symmetry property.
//Example:
//  DenseVector<float>::self X(4, 1);
//  fft_in_place(X);
template<class G> inline void fft_in_place(Vector<G> &X) { aux_fft_in_place(X,typename Vector<G>::value_type()); }






//////////////////////////////////////////////////////////////////////////////////////////////////

template<class T1,class T2, class F>
struct row_fft_function
{
  T1 &X; T2 &Y; F &f;
  inline row_fft_function(T1 &x, T2 &y, F &f0) : X(x),Y(y), f(f0) {}
  inline void operator()(int i,int imax) const
  {
    for (; i<imax; ++i)
    {
      typename MatrixRow<T1>::self X2 = row(X,i);
      typename MatrixRow<T2>::self Y2 = row(Y,i);
      fft(X2,Y2,f);
    }
  }
};

template<class G1,class G2,class F,class V>
void aux_row_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f, const V&)
{
#if defined(NDEBUG) && defined(FFT_THREADING)
  row_fft_function<const Matrix<G1>,Matrix<G2>,const F> func(X,Y,f);
  int grain; if (4096*num_threads()<=X.nelms()) grain=(X.nrows()-1)/((X.nelms()-1)/4096+1)+1; else if (16384*num_threads()<=X.nelms()) grain=(X.nrows()-1)/num_threads()+1; else grain=16384/X.ncols();
  parallel_for(X.row_lower_bound(),X.row_upper_bound()+1,func,grain);
#else
  int imax=X.row_upper_bound();
  for (int i=X.row_lower_bound(); i<=imax; ++i)
  {
    typename MatrixRow<const Matrix<G1> >::self X2 = row(X,i);
    typename MatrixRow<      Matrix<G2> >::self Y2 = row(Y,i);
    fft(X2,Y2,f);
  }
#endif
}

template<class G1,class G2,class F,class V>
void aux_row_fft(      Matrix<G1> &X, Matrix<G2> &Y, const F &f, const V&)
{
#if defined(NDEBUG) && defined(FFT_THREADING)
  row_fft_function<Matrix<G1>,Matrix<G2>,const F> func(X,Y,f);
  int grain; if (4096*num_threads()<=X.nelms()) grain=(X.nrows()-1)/((X.nelms()-1)/4096+1)+1; else if (16384*num_threads()<=X.nelms()) grain=(X.nrows()-1)/num_threads()+1; else grain=16384/X.ncols();
  parallel_for(X.row_lower_bound(),X.row_upper_bound()+1,func,grain);
#else
  for (int i=X.row_lower_bound(); i<=X.row_upper_bound(); ++i)
  {
    typename MatrixRow<Matrix<G1> >::self X2 = row(X,i);
    typename MatrixRow<Matrix<G2> >::self Y2 = row(Y,i);
    fft(X2,Y2,f);
  }
#endif
}

template<class G1,class G2,class F>
inline void aux_row_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  aux_row_fft(X,Y,f,typename Matrix<G1>::value_type());
}

template<class G1,class G2,class F>
inline void aux_row_fft(      Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  aux_row_fft(X,Y,f,typename Matrix<G1>::value_type());
}

//{secret}
template<class G1,class G2,class F>
inline void row_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  typename ArrayData<const Matrix<G1> >::self X2 = data(X);
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_row_fft(X2,Y2,f);
}

//{secret}
template<class G1,class G2,class F>
inline void row_fft(      Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  typename ArrayData<Matrix<G1> >::self X2 = data(X);
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_row_fft(X2,Y2,f);
}


//{unsecret}
//Summary: FFT of every row of a real signal
//Arguments: X - The real signal
//Output: Y - The half complex FFT of every row
//Example:
//  DenseMatrix<float>::self X(4,4, 1);
//  DenseMatrix<complex<float> >::self Y(X.size());
//  row_fft(X,Y);
template<class G1,class G2> inline void half_real_row_fft(const Matrix<G1> &X, Matrix<G2> &Y) { row_fft(X,Y,make_fft_adaptor     ()); }


template<class G1,class G2,class V> inline void val_aux_row_fft(const Matrix<G1> &X, Matrix<G2> &Y, const complex<V> &) { row_fft(X,Y,make_fft_adaptor     ( )); }
template<class G1,class G2,class V> inline void val_aux_row_fft(const Matrix<G1> &X, Matrix<G2> &Y, const V          &) { row_fft(X,Y,make_full_fft_adaptor(Y)); }

//{unsecret}
//Summary: FFT of every row of a signal
//Arguments: X - A signal
//Output: Y - The complex FFT of every row
//Example:
//  DenseMatrix<float>::self X(4,4, 1);
//  DenseMatrix<complex<float> >:: self Y(X.size());
//  row_fft(X,Y);
//
//  DenseMatrix<complex<double> >::self X(4,4, 1);
//  DenseMatrix<complex<double> >:: self Y(X.size());
//  row_fft(X,Y);
template<class G1,class G2> inline void row_fft(const Matrix<G1> &X, Matrix<G2> &Y) { val_aux_row_fft(X,Y,typename Matrix<G1>::value_type()); }


//{unsecret}
//Summary: IFFT of every row of a signal
//Arguments: X - A signal
//Output: Y - The complex IFFT of every row
//Remarks:
//  If /Y/ contains real values, then X is supposed to be complex symetric and /real_ifft/ is performed on every row.
//See: ^ifft^, ^real_ifft^
//Example:
//  DenseMatrix<complex<float> >::self X(4,4, 1);
//  DenseMatrix<complex<float> >::self Y(X.size());
//  row_ifft(X,Y);
//
//  DenseMatrix<complex<double> >::self X(4,4, 1);
//  DenseMatrix<double>::self Y(X.size());
//  row_ifft(X,Y);
template<class G1,class G2> inline void row_ifft(      Matrix<G1> &X, Matrix<G2> &Y) { row_fft(X,Y,make_ifft_adaptor(Y)); }
template<class G1,class G2> inline void row_ifft(const Matrix<G1> &X, Matrix<G2> &Y) { row_fft(X,Y,make_ifft_adaptor(Y)); }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class T1,class T2, class F>
struct col_fft_function
{
  T1 &X; T2 &Y; F &f;
  inline col_fft_function(T1 &x, T2 &y, F &f0) : X(x),Y(y), f(f0) {}
  inline void operator()(int j,int jmax) const
  {
    for (; j<jmax; ++j)
    {
      typename MatrixCol<T1>::self X2 = col(X,j);
      typename MatrixCol<T2>::self Y2 = col(Y,j);
      fft(X2,Y2,f);
    }
  }
};

//{secret}
template<class G1,class G2,class F>
void matrix_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
#if defined(NDEBUG) && defined(FFT_THREADING)
  col_fft_function<const Matrix<G1>,Matrix<G2>,const F> func(X,Y,f);
  int grain; if (4096*num_threads()<=X.nelms()) grain=(X.ncols()-1)/((X.nelms()-1)/4096+1)+1; else if (16384*num_threads()<=X.nelms()) grain=(X.ncols()-1)/num_threads()+1; else grain=16384/X.nrows();
  parallel_for(X.col_lower_bound(),X.col_upper_bound()+1,func,grain);
#else
  int jmax=X.col_upper_bound();
  for (int j=X.col_lower_bound(); j<=jmax; ++j)
  {
    typename MatrixCol<const Matrix<G1> >::self X2 = col(X,j);
    typename MatrixCol<      Matrix<G2> >::self Y2 = col(Y,j);
    fft(X2,Y2,f);
  }
#endif
}

template<class G1,class G2,class F,class V>
inline void aux_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f, const V &)
{
  matrix_col_fft(X,Y,f);
}

#ifdef SSE
template<class G1,class G2,class F>
inline void aux_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f, complex<float>)
{
  if (X.ncols()%2==0 && Y.ncols()%2==0)
  {
    typename SimdBlock<Matrix<G2> >::self Y2=simd_block(Y);
    matrix_col_fft(simd_block(X),Y2,f);
  }
  else
    matrix_col_fft(X,Y,f);
}
#endif

#ifdef SSE2
template<class G1,class G2,class F>
inline void aux_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f, complex<double>)
{
  typename SimdBlock<Matrix<G2> >::self Y2=simd_block(Y);
  matrix_col_fft(simd_block(X),Y2,f);
}
#endif


//{secret}
template<class G1,class G2,class F>
void col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_col_fft(data(X),Y2,f,typename Vector<G1>::value_type());
}

template<class G1,class G2> inline void complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y) { col_fft(X,Y,make_fft_adaptor     ( )); }
template<class G1,class G2> inline void real_col_fft   (const Matrix<G1> &X, Matrix<G2> &Y) { col_fft(X,Y,make_full_fft_adaptor(Y)); }


template<class G1,class G2,class V> inline void aux_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const V          &) { col_fft(X,Y,make_full_fft_adaptor(Y)); }
template<class G1,class G2,class V> inline void aux_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, const complex<V> &) { col_fft(X,Y,make_fft_adaptor     ( )); }

//{secret}
//Summary: FFT of every column of a signal
//Arguments: X - A signal
//Output: Y - The complex FFT of every column
//Example:
//  DenseMatrix<float>:: self X(4,4, 1);
//  DenseMatrix<complex<float> >:: self Y(X.size());
//  col_fft(X,Y);
template<class G1,class G2> inline void col_fft(const Matrix<G1> &X, Matrix<G2> &Y) { aux_col_fft(X,Y,typename Matrix<G1>::value_type()); }

//{secret}
//Summary: FFT of every column of a real signal
//Arguments: X - The real signal
//Output: Y - The half complex FFT of every column
template<class G1,class G2> inline void half_real_col_fft(const Matrix<G1> &X, Matrix<G2> &Y) { col_fft(X,Y,make_fft_adaptor()); }



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class G1,class G2,class G3,class F>
inline void complex_matrix_fft(const Matrix<G1> &X, Matrix<G2> &T, Matrix<G3> &Y,const F &f)
{
  row_fft(X,T,f(T));
  col_fft(T,Y,f);
}

struct complex_matrix_fft_function
{
  template<class G1,class G2,class F>
  inline void operator()(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
  {
    typedef PROMOTE2(complex<float>,typename Matrix<G1>::value_type) value_type;
    //typename DenseMatrix<value_type>::self T(X.size());
    typename DenseVector<char>::self &T=fft_buf2; if (T.size()<int(sizeof(value_type))*X.nelms()) T.resize(sizeof(value_type)*X.nelms());

    typename DataMatrix<value_type>::self T2(X.size(),(value_type*)&*T.begin());
    complex_matrix_fft(X,T2,Y,f);
  }

  template<class V1,class V2,class F>
  inline void operator()(const Matrix<data_matrix_generator<V1> > &X, Matrix<data_matrix_generator<complex<V2> > > &Y, const F &f)
  {
    typedef complex<V2> value_type;
    typename DataMatrix<value_type>::self T2(Y.size(),&*Y.begin());
    complex_matrix_fft(X,T2,Y,f);
  }
};

template<class G1,class G2,class F> void aux_complex_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f) { complex_matrix_fft_function()(X,Y,f); }

template<class G1,class G2,class F> inline void complex_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_complex_fft(data(X),Y2,f);
}

template<class G1,class G2> inline void complex_fft     (const Matrix<G1> &X, Matrix<G2> &Y) { complex_fft(X,Y,make_fft_adaptor     ()); }
template<class G1,class G2> inline void complex_ifft    (const Matrix<G1> &X, Matrix<G2> &Y) { complex_fft(X,Y,make_ifft_adaptor    (Y)); }
template<class G1,class G2> inline void complex_norm_fft(const Matrix<G1> &X, Matrix<G2> &Y) { complex_fft(X,Y,make_norm_fft_adaptor()); }
template<class G1,class G2> inline void complex_abs_fft (const Matrix<G1> &X, Matrix<G2> &Y) { complex_fft(X,Y,make_abs_fft_adaptor ()); }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T1,class T2, class F>
struct half_matrix_half_complex_col_fft_function
{
  T1 &X; T2 &Y; F &f;
  inline half_matrix_half_complex_col_fft_function(T1 &x, T2 &y, F &f0) : X(x),Y(y), f(f0) {}
  inline void operator()(int j,int jmax) const
  {
    for (; j<jmax; ++j)
    {
      typename MatrixCol<T1>::self X2 = col(X,j);
      typename MatrixCol<T2>::self Y2 = col(Y,j);
      complex_fft(X2,Y2,f(Y2));
    }
  }
};

template<class G1,class G2,class F>
void half_matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f)
{
#if defined(NDEBUG) && defined(FFT_THREADING)
  half_matrix_half_complex_col_fft_function<const Matrix<G1>,Matrix<G2>,const F> func(X,Y,f);
  int grain, m=X.nrows(); n=n/2+1; if (4096*num_threads()<=m*n) grain=(n-1)/((m*n-1)/4096+1)+1; else if (16384*num_threads()<=m*n) grain=(n-1)/num_threads()+1; else grain=16384/m;
  parallel_for(0,n,func,grain);
#else
  for (int j=0, jmax=n/2; j<=jmax; ++j)
  {
    typename MatrixCol<Matrix<G2> >::self Y2=col(Y,j);
    complex_fft(col(X,j),Y2,f(Y2));
  }
#endif
}

template<int N,class G1,class G2,class F>
void half_matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f)
{
  typename SimdBlock<const Matrix<G1> >::self X2 = simd_block(X);
  typename SimdBlock<      Matrix<G2> >::self Y2 = simd_block(Y);

#if defined(NDEBUG) && defined(FFT_THREADING)
  half_matrix_half_complex_col_fft_function<typename SimdBlock<const Matrix<G1> >::self,typename SimdBlock<Matrix<G2> >::self,const F> func(X2,Y2,f);
  int grain, m=X.nrows(); n=n/(2*N)+1; if (4096*num_threads()<=m*n) grain=(n-1)/((m*n-1)/4096+1)+1; else if (16384*num_threads()<=m*n) grain=(n-1)/num_threads()+1; else grain=16384/m;
  parallel_for(0,n,func,grain);
#else
  for (int j=0, jmax=n/(2*N); j<=jmax; ++j)
  {
    typename MatrixCol<typename SimdBlock<Matrix<G2> >::self>::self T=col(Y2,j);
    complex_fft(col(X2,j),T,f(T));
  }
#endif
}


template<class T1,class T2, class F>
struct full_matrix_half_complex_col_fft_function
{
  T1 &X; T2 &Y; F &f; int n;
  inline full_matrix_half_complex_col_fft_function(T1 &x, T2 &y, F &f0, int n0) : X(x),Y(y), f(f0), n(n0) {}
  inline void operator()(int j,int jmax) const
  {
    for (; j<jmax; ++j)
    {
      typename MatrixCol<T1>::self X1=col(X,j);
      typename MatrixCol<T2>::self Y1=col(Y,j),Y2=col(Y,n-j);
      complex_fft(X1,Y1,f(Y1,Y2));
    }
  }
};

template<class G1,class G2,class F>
void full_matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f)
{
  {
    typename MatrixCol<Matrix<G2> >::self T=col(Y,0);
    complex_fft(col(X,0),T,f(T));
  }
#if defined(NDEBUG) && defined(FFT_THREADING)
  full_matrix_half_complex_col_fft_function<const Matrix<G1>,Matrix<G2>,const F> func(X,Y,f,n);
  int grain, m=X.nrows(), l=n/2; if (4096*num_threads()<=m*l) grain=(l-1)/((m*l-1)/4096+1)+1; else if (16384*num_threads()<=m*l) grain=(l-1)/num_threads()+1; else grain=16384/m;
  parallel_for(1,l+1,func,grain);
#else
  for (int j=1, jmax=n/2; j<=jmax; ++j)
  {
    typename MatrixCol<Matrix<G2> >::self T1=col(Y,j), T2=col(Y,n-j);
    complex_fft(col(X,j),T1,f(T1,T2));
  }
#endif
}

template<int N,class G1,class G2,class F>
void full_matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f)
{
  {
    typename MatrixCol<Matrix<G2> >::self T=col(Y,0);
    complex_fft(col(X,0),T,f(T));
  }
  for (int j=1; j<N; ++j)
  {
    typename MatrixCol<Matrix<G2> >::self T1=col(Y,j), T2=col(Y,n-j);
    complex_fft(col(X,j),T1,f(T1,T2));
  }

  typename SimdBlock<const Matrix<G1> >::self X2 = simd_block(X);
  typename SimdBlock<      Matrix<G2> >::self Y2 = simd_block(Y);
#if defined(NDEBUG) && defined(FFT_THREADING)
  full_matrix_half_complex_col_fft_function<typename SimdBlock<const Matrix<G1> >::self,typename SimdBlock<Matrix<G2> >::self,const F> func(X2,Y2,f,n/N);
  int grain, m=X.nrows(), l=(n/2+1)/N-1; if (4096*num_threads()<=m*l) grain=(l-1)/((m*l-1)/4096+1)+1; else if (16384*num_threads()<=m*l) grain=(l-1)/num_threads()+1; else grain=16384/m;
  parallel_for(1,1+l,func,grain);
#else
  for (int j=1, jmax=(n/2+1)/N; j<jmax; ++j)
  {
    typename MatrixCol<typename SimdBlock<Matrix<G2> >::self>::self T1=col(Y2,j), T2=col(Y2,n/N-j);
    complex_fft(col(X2,j),T1,f(T1,T2));
  }
#endif

  for (int j=((n/2+1)/N)*N, jmax=n/2; j<=jmax; ++j)
  {
    typename MatrixCol<Matrix<G2> >::self T1=col(Y,j), T2=col(Y,n-j);
    complex_fft(col(X,j),T1,f(T1,T2));
  }
}

template<      class G1,class G2        > inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_fft_adaptor                  &f) { half_matrix_half_complex_col_fft   (X,Y,n,f); }
template<int N,class G1,class G2        > inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_fft_adaptor                  &f) { half_matrix_half_complex_col_fft<N>(X,Y,n,f); }
template<      class G1,class G2        > inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_full_fft_adaptor             &f) { full_matrix_half_complex_col_fft   (X,Y,n,f); }
template<int N,class G1,class G2        > inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_full_fft_adaptor             &f) { full_matrix_half_complex_col_fft<N>(X,Y,n,f); }
template<      class G1,class G2        > inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_ifft_adaptor                 &f) { half_matrix_half_complex_col_fft   (X,Y,n,f); }
template<int N,class G1,class G2        > inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_ifft_adaptor                 &f) { half_matrix_half_complex_col_fft<N>(X,Y,n,f); }
template<      class G1,class G2,class F> inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_function_fft_adaptor     <F> &f) { half_matrix_half_complex_col_fft   (X,Y,n,f); }
template<int N,class G1,class G2,class F> inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_function_fft_adaptor     <F> &f) { half_matrix_half_complex_col_fft<N>(X,Y,n,f); }
template<      class G1,class G2,class F> inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_function_full_fft_adaptor<F> &f) { full_matrix_half_complex_col_fft   (X,Y,n,f); }
template<int N,class G1,class G2,class F> inline void matrix_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const make_function_full_fft_adaptor<F> &f) { full_matrix_half_complex_col_fft<N>(X,Y,n,f); }


template<class G1,class G2,class F,class V>
inline void aux_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f, const V &)
{
  matrix_half_complex_col_fft(X,Y,n,f);
}

#ifdef SSE
template<class G1,class G2,class F>
inline void aux_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f, complex<float>)
{
  if (X.ncols()%2==0 && Y.ncols()%2==0)
    matrix_half_complex_col_fft<2>(X,Y,n,f);
  else
    matrix_half_complex_col_fft(X,Y,n,f);
}
#endif

#ifdef SSE2
template<class G1,class G2,class F>
inline void aux_half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f, complex<double>)
{
  typename SimdBlock<Matrix<G2> >::self Y2=simd_block(Y);
  matrix_half_complex_col_fft(simd_block(X),Y2,n,f);
}
#endif


template<class G1,class G2,class F>
inline void half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y, int n, const F &f)
{
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_half_complex_col_fft(data(X),Y2,n,f,typename Vector<G1>::value_type());
}

template<class G1,class G2> void half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y,int n) { half_complex_col_fft(X,Y,n,make_fft_adaptor()); }
template<class G1,class G2> void half_complex_col_fft(const Matrix<G1> &X, Matrix<G2> &Y) { half_complex_col_fft(X,Y,X.ncols()); }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class G1,class G2,class G3,class F>
inline void real_matrix_fft(const Matrix<G1> &X, Matrix<G2> &T, Matrix<G3> &Y,const F &f)
{
  row_fft(X,T,f(T));
  half_complex_col_fft(T,Y,X.ncols(),f);
}

template<class G1,class G2,class G3,class F>
inline void real_matrix_ifft(const Matrix<G1> &X, Matrix<G2> &T, Matrix<G3> &Y,const F &f)
{
  half_complex_col_fft(X,T,Y.ncols(),f(T));
  row_fft(T,Y,f);
}

struct real_matrix_fft_function
{
  template<class G1,class G2,class F>
  inline void operator()(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
  {
    int m=X.nrows();
    int n=X.ncols();

    typedef PROMOTE2(complex<float>,typename Matrix<G1>::value_type) value_type;
    //typename DenseMatrix<value_type>::self T(m,n); typename DataMatrix<value_type>::self T2(m,n, (value_type*)&*T.begin());
    if (fft_buf2.n<int(sizeof(value_type)*(m*n))) { fft_buf2.n=sizeof(value_type)*(m*n); aligned_free(fft_buf2.p); fft_buf2.p=aligned_malloc(fft_buf2.n); } typename DataMatrix<value_type>::self T2(m,n, (value_type *)fft_buf2.p);
    real_matrix_fft(X,T2,Y,f);
  }

  template<class V1,class V2,class F>
  inline void operator()(const Matrix<data_matrix_generator<V1> > &X, Matrix<data_matrix_generator<complex<V2> > > &Y, const F &f)
  {
    typedef complex<V2> value_type;
    typename DataMatrix<value_type>::self T2(Y.size(),&*Y.begin());
    real_matrix_fft(X,T2,Y,f);
  }
};

struct real_matrix_ifft_function
{
  template<class G1,class G2,class F>
  inline void operator()(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
  {
    typedef typename Matrix<G1>::value_type value_type;
    int m=Y.nrows();
    int n=((Y.ncols()/2+1)/SimdBlockLength<value_type>::RET+1)*SimdBlockLength<value_type>::RET;
    //typename DenseMatrix<value_type>::self T(m,n); typename DataMatrix<value_type>::self T2(m,n, (value_type*)&*T.begin());
    if (fft_buf2.n<int(sizeof(value_type)*(m*n))) { fft_buf2.n=sizeof(value_type)*(m*n); aligned_free(fft_buf2.p); fft_buf2.p=aligned_malloc(fft_buf2.n); } typename DataMatrix<value_type>::self T2(m,n, (value_type *)fft_buf2.p);
    real_matrix_ifft(X,T2,Y,f);
  }
};

template<class G1,class G2,class F> void aux_real_fft(const Matrix<G1> &X, Matrix<G2> &Y, const                 F &f) { real_matrix_fft_function ()(X,Y,f); }
template<class G1,class G2        > void aux_real_fft(const Matrix<G1> &X, Matrix<G2> &Y, const make_ifft_adaptor &f) { real_matrix_ifft_function()(X,Y,f); }

template<class G1,class G2,class F> inline void real_fft(const Matrix<G1> &X, Matrix<G2> &Y, const F &f)
{
  typename ArrayData<Matrix<G2> >::self Y2 = data(Y);
  aux_real_fft(data(X),Y2,f);
}

//{unsecret}
template<class G1,class G2> inline void half_real_fft(const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_fft_adaptor          ( )); }
template<class G1,class G2> inline void real_fft          (const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_full_fft_adaptor     (Y)); }
template<class G1,class G2> inline void half_real_norm_fft(const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_norm_fft_adaptor     ( )); }
template<class G1,class G2> inline void real_norm_fft     (const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_norm_full_fft_adaptor( )); }
template<class G1,class G2> inline void half_real_abs_fft (const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_abs_fft_adaptor      ( )); }
template<class G1,class G2> inline void real_abs_fft      (const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_abs_full_fft_adaptor ( )); }
//{unsecret}
template<class G1,class G2> inline void real_ifft(const Matrix<G1> &X, Matrix<G2> &Y) { real_fft(X,Y,make_ifft_adaptor(Y)); }

//{unsecret}
template<class G> inline typename DenseMatrix<typename Matrix<G>::value_type::value_type >::self real_ifft(const Matrix<G> &X)
{
  typename DenseMatrix<typename Vector<G>::value_type::value_type>::self Y(X.size());
  real_ifft(X,Y);
  return Y;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class G1,class G2,class V> inline void aux_fft (const Matrix<G1> &X, Matrix<G2> &Y, const V          &) { real_fft    (X,Y); }
template<class G1,class G2,class V> inline void aux_fft (const Matrix<G1> &X, Matrix<G2> &Y, const complex<V> &) { complex_fft (X,Y); }
template<class G1,class G2,class V> inline void aux_ifft(const Matrix<G1> &X, Matrix<G2> &Y, const V          &) { real_ifft   (X,Y); }
template<class G1,class G2,class V> inline void aux_ifft(const Matrix<G1> &X, Matrix<G2> &Y, const complex<V> &) { complex_ifft(X,Y); }

//{unsecret}
template<class G1,class G2> inline void fft (const Matrix<G1> &X, Matrix<G2> &Y) { aux_fft (X,Y,typename Matrix<G1>::value_type()); }
//{unsecret}
template<class G1,class G2> inline void ifft(const Matrix<G1> &X, Matrix<G2> &Y) { aux_ifft(X,Y,typename Matrix<G2>::value_type()); }
//{unsecret}
template<class G> inline PROMOTE2(complex<float>,Matrix<G>) fft (const Matrix<G> &X) { PROMOTE2(complex<float>,Matrix<G>) Y(X.size()); fft (X,Y); return Y; }
//{unsecret}
template<class G> inline PROMOTE2(complex<float>,Matrix<G>) ifft(const Matrix<G> &X) { PROMOTE2(complex<float>,Matrix<G>) Y(X.size()); ifft(X,Y); return Y; }



//template<class G1,class G2,class V> inline void aux_norm_fft(const Matrix<G1> &X, Matrix<G2> &Y, const V          &) { return real_norm_fft   (X,Y); }
//template<class G1,class G2,class V> inline void aux_norm_fft(const Matrix<G1> &X, Matrix<G2> &Y, const complex<V> &) { return complex_norm_fft(X,Y); }
////{unsecret}
//template<class G1,class G2> inline void norm_fft(const Matrix<G1> &X, Matrix<G2> &Y) { return aux_norm_fft(X,Y,Matrix<G1>::value_type()); }
////{unsecret}
//template<class G1> inline PROMOTE2(complex<float>,Matrix<G1>) norm_fft(const Matrix<G1> &X) { PROMOTE2(complex<float>,Matrix<G1>) Y(X.size()); norm_fft(X,Y); return Y; }
//
//
//template<class G1,class G2,class V> inline void aux_abs_fft(const Matrix<G1> &X, Matrix<G2> &Y, const V          &) { return real_abs_fft   (X,Y); }
//template<class G1,class G2,class V> inline void aux_abs_fft(const Matrix<G1> &X, Matrix<G2> &Y, const complex<V> &) { return complex_abs_fft(X,Y); }
////{unsecret}
//template<class G1,class G2> inline void abs_fft(const Matrix<G1> &X, Matrix<G2> &Y) { return aux_abs_fft(X,Y,Matrix<G1>::value_type()); }
////{unsecret}
//template<class G1> inline PROMOTE2(complex<float>,Matrix<G1>) abs_fft(const Matrix<G1> &X) { PROMOTE2(complex<float>,Matrix<G1>) Y(X.size()); abs_fft(X,Y); return Y; }



template<class G> Vector<G> &bitreverse_order(Vector<G> &X)
{
  size_t i;
  size_t j = 0;

  size_t n = X.size();

  for (i = 0; i < n - 1; ++i)
  {
    size_t k=n/2;
    if (i<j) swap(X[i],X[j]);
    while (k<=j)
    {
      j=j-k; k=k/2;
    }
    j+=k;
  }
  return X;
}



template<class G>
Vector<G> &complex_radix2_fft_in_place(Vector<G> &X)
{
  typedef Vector<G> array_type;
  typedef typename array_type::size_type size_type;
  typedef PROMOTE2(complex<float>,typename array_type::value_type) value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  size_type n = X.size();

  assert(is_pow2(n));

  bitreverse_order(X);

  for (int dual=1; dual<n; dual*=2)
  {
    real_type theta = -PI/dual;
    real_type s  = sin(theta);
    real_type s2 = real_type(2)*sqr(sin(real_type(0.5)*theta));

    complex_type w(1,0);

    value_type wd;
    for (size_t a=0; a<dual; ++a, w+=s*complex_type(-w.imag(),w.real())-s2*w)
      for (size_t b=0; b<n; b+=2*dual)
      {
        value_type &Ti=X[b+a];
        value_type &Tj=X[b+a+dual];
        wd=w*Tj;
        Tj=Ti-wd;
        Ti+=wd;
      }
  }
  return X;
}


template<class G>
void real_radix2_fft_in_place(Vector<G> &X)
{
  typedef Vector<G> array_type;
  typedef typename array_type::size_type size_type;
  typedef PROMOTE2(float,typename array_type::value_type) value_type;
  typedef PROMOTE2(complex<float>,typename array_type::value_type) complex_value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  size_type n = X.size();
  assert(is_pow2(n));

  bitreverse_order(X);

  for (size_t p_1=1, p=2, q=n/2; p_1<n; p_1=p, p*=2, q/=2)
  {
    complex_type w(1,0);
    real_type theta = -2*PI/p;
    real_type s  = sin (theta);
    real_type s2 = 2*sqr(sin(theta/2));

    {
      value_type t0;
      for (size_type b=0; b<q; b++)
      {
        typename array_type::reference r0 = X[b*p];
        typename array_type::reference r1 = X[b*p+p_1];

        t0 = r0+r1;
        r1 = r0-r1;
        r0 = t0;
	    }
    }

    complex_value_type z0, z1, z2;
    for (size_t a=1; a<p_1/2; ++a)
    {
      w+=s*complex_type(-w.imag(),w.real())-s2*w;
      for (size_t b=0; b<q; ++b)
      {
        typename array_type::reference r0 = X[b*p+a];
        typename array_type::reference r1 = X[b*p+p_1-a];
        typename array_type::reference r2 = X[b*p+p_1+a];
        typename array_type::reference r3 = X[b*p+p-a];

        real(z0)=r2; imag(z0)=r3;
        real(z1)=r0; imag(z1)=r1;

        z0 *= w;

        z2 = z1+z0;
        z1 -= z0;

        r0 = real(z2);
        r3 = imag(z2);

        r1 = real(z1);
        r2 = -imag(z1);
      }
    }
    if (p_1>1) for (size_t b=0; b<q; b++) X[b*p+p-p_1/2]*=double(-1);
  }
}

template<class G>
void real_radix2_ifft_in_place(Vector<G> &X)
{
  typedef Vector<G> array_type;
  typedef typename array_type::size_type size_type;
  typedef PROMOTE2(float,typename array_type::value_type) value_type;
  typedef PROMOTE2(complex<float>,typename array_type::value_type) complex_value_type;
  typedef PROMOTE2(complex<float>,VALUE_TYPE(value_type)) complex_type;
  typedef typename complex_type::value_type real_type;

  size_type n = X.size();
  assert(is_pow2(n));

  for (size_t p_1=n/2, p=n, q=1; q<n; p=p_1, p_1/=2, q*=2)
  {
    complex_type w(1,0);
    real_type theta = 2*PI/p;
    real_type s  = sin (theta);
    real_type s2 = 2*sqr(sin(theta/2));

    {
      value_type t0;
      for (size_t b=0; b<q; b++)
      {
        typename array_type::reference r0 = X[b*p];
        typename array_type::reference r1 = X[b*p+p_1];

        t0 = r0+r1;
        r1 = r0-r1;
        r0 = t0;
	    }
	  }

    complex_value_type z0;
    for (size_t a=1; a<p_1/2; ++a)
    {
      w+=s*complex_type(-w.imag(),w.real())-s2*w;
      for (size_t b=0; b<q; ++b)
      {
        typename array_type::reference r0 = X[b*p+a];
        typename array_type::reference r1 = X[b*p+p-a];
        typename array_type::reference r2 = X[b*p+p_1-a];
        typename array_type::reference r3 = X[b*p+p_1+a];

        real(z0)=r0-r2; imag(z0)=r1+r3;
        z0 *= w;

        r0 += r2;
        r2 = r1-r3;

        r3 = real(z0);
        r1 = imag(z0);
      }
    }
    if (p_1>1) for (size_t b=0; b<q; b++) { X[b*p+p_1/2]*=2; X[b*p+p_1+p_1/2]*=-2; }

  }
  bitreverse_order(X);
  X/=n;
}


#ifdef FFT_PRECOMPILE

#define DECLARE_FFT(V)\
GENIAL_API void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_fft_adaptor           &f);\
GENIAL_API void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_fft_adaptor           &f);\
GENIAL_API void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_full_fft_adaptor      &f);
#define DEFINE_FFT(V)\
void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_fft_adaptor           &f) { aux_complex_fft(X,Y,f,complex<V>()); }\
void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_fft_adaptor           &f) { aux_real_fft   (X,Y,f,V()); }\
void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_full_fft_adaptor      &f) { aux_real_fft   (X,Y,f,V()); }

#define DECLARE_IFFT(V)\
GENIAL_API void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_ifft_adaptor &f);\
GENIAL_API void aux_real_ifft  (      Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<        V  > > &Y                            );
#define DEFINE_IFFT(V)\
void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<complex<V> > > &Y, const make_ifft_adaptor &f) { aux_complex_fft(X,Y,f,complex<V>()); }\
void aux_real_ifft  (      Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<        V  > > &Y                            ) { aux_real_ifft  (X,Y,  complex<V>()); }

#define DECLARE_ROW_FFT(V)\
GENIAL_API void aux_row_fft    (const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor           &f);\
GENIAL_API void aux_row_fft    (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor           &f);\
GENIAL_API void aux_row_fft    (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_full_fft_adaptor      &f);
#define DEFINE_ROW_FFT(V)\
void aux_row_fft    (const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor           &f) { aux_row_fft(X,Y,f,complex<V>()); }\
void aux_row_fft    (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor           &f) { aux_row_fft(X,Y,f,V         ()); }\
void aux_row_fft    (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_full_fft_adaptor      &f) { aux_row_fft(X,Y,f,V         ()); }

#define DECLARE_ROW_IFFT(V)\
GENIAL_API void aux_row_fft    (const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_ifft_adaptor &f);\
GENIAL_API void aux_row_fft    (      Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<        V  > > &Y, const make_ifft_adaptor &f);
#define DEFINE_ROW_IFFT(V)\
void aux_row_fft    (const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_ifft_adaptor &f) { aux_row_fft(X,Y,f,complex<V>()); }\
void aux_row_fft    (      Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<        V  > > &Y, const make_ifft_adaptor &f) { aux_row_fft(X,Y,f,V         ()); }

#define DECLARE_FFT2(V)\
GENIAL_API void aux_complex_fft(const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor      &f);\
GENIAL_API void aux_real_fft   (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor      &f);\
GENIAL_API void aux_real_fft   (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_full_fft_adaptor &f);
#define DEFINE_FFT2(V)\
void aux_complex_fft(const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor      &f) { complex_matrix_fft_function()(X,Y,f); }\
void aux_real_fft   (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_fft_adaptor      &f) { real_matrix_fft_function   ()(X,Y,f); }\
void aux_real_fft   (const Matrix<data_matrix_generator<        V  > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_full_fft_adaptor &f) { real_matrix_fft_function   ()(X,Y,f); }

#define DECLARE_IFFT2(V)\
GENIAL_API void aux_complex_fft(const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_ifft_adaptor     &f);\
GENIAL_API void aux_real_fft   (const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<        V  > > &Y, const make_ifft_adaptor     &f);
#define DEFINE_IFFT2(V)\
void aux_complex_fft(const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<complex<V> > > &Y, const make_ifft_adaptor     &f) { complex_matrix_fft_function()(X,Y,f); }\
void aux_real_fft   (const Matrix<data_matrix_generator<complex<V> > > &X, Matrix<data_matrix_generator<        V  > > &Y, const make_ifft_adaptor     &f) { real_matrix_ifft_function  ()(X,Y,f); }

#define DECLARE_NORM_FFT(V)\
GENIAL_API void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<        V  > > &Y, const make_norm_fft_adaptor      &f);\
GENIAL_API void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_norm_fft_adaptor      &f);\
GENIAL_API void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_norm_full_fft_adaptor &f);
#define DEFINE_NORM_FFT(V)\
void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<        V  > > &Y, const make_norm_fft_adaptor      &f) { aux_complex_fft(X,Y,f,complex<V>()); }\
void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_norm_fft_adaptor      &f) { aux_real_fft   (X,Y,f,V()); }\
void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_norm_full_fft_adaptor &f) { aux_real_fft   (X,Y,f,V()); }

#define DECLARE_ABS_FFT(V)\
GENIAL_API void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<        V  > > &Y, const make_abs_fft_adaptor       &f);\
GENIAL_API void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_abs_fft_adaptor       &f);\
GENIAL_API void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_abs_full_fft_adaptor  &f);
#define DEFINE_ABS_FFT(V)\
void aux_complex_fft(const Vector<data_vector_generator<complex<V> > > &X, Vector<data_vector_generator<        V  > > &Y, const make_abs_fft_adaptor       &f) { aux_complex_fft(X,Y,f,complex<V>()); }\
void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_abs_fft_adaptor       &f) { aux_real_fft   (X,Y,f,V()); }\
void aux_real_fft   (const Vector<data_vector_generator<        V  > > &X, Vector<data_vector_generator<        V  > > &Y, const make_abs_full_fft_adaptor  &f) { aux_real_fft   (X,Y,f,V()); }

DECLARE_FFT(float)
DECLARE_FFT(double)
DECLARE_IFFT(float)
DECLARE_IFFT(double)

DECLARE_ROW_FFT(float)
DECLARE_ROW_FFT(double)
DECLARE_ROW_IFFT(float)
DECLARE_ROW_IFFT(double)

DECLARE_FFT2(float)
DECLARE_FFT2(double)
DECLARE_IFFT2(float)
DECLARE_IFFT2(double)

DECLARE_NORM_FFT(float)
DECLARE_NORM_FFT(double)

DECLARE_ABS_FFT(float)
DECLARE_ABS_FFT(double)

#endif
#endif //__cplusplus


#ifdef __cplusplus
extern "C" {
#endif

#ifdef FFT_PRECOMPILE

//Group = FFT C Interface


//{unsecret}
//{group:FFT C Interface}
GENIAL_API void shfft(int nx, int ny, const float *px, float *py);
//{unsecret}
GENIAL_API void cfft(int n, const float *px, float *py);
//{unsecret}
GENIAL_API void sfft(int n, const float *px, float *py);
//{unsecret}
GENIAL_API void zfft(int n, const double *px, double *py);
//{unsecret}
GENIAL_API void dhfft(int nx, int ny, const double *px, double *py);
//{unsecret}
GENIAL_API void dfft(int n, const double *px, double *py);

//{unsecret}
GENIAL_API void cfft1(int m, int n, const float *px, float *py);
//{unsecret}
GENIAL_API void shfft1(int m, int nx, int ny, const float *px, float *py);
//{unsecret}
GENIAL_API void sfft1(int m, int n, const float *px, float *py);
//{unsecret}
GENIAL_API void zfft1(int m, int n, const double *px, double *py);
//{unsecret}
GENIAL_API void dhfft1(int m, int nx, int ny, const double *px, double *py);
//{unsecret}
GENIAL_API void dfft1(int m, int n, const double *px, double *py);

//{unsecret}
GENIAL_API void cfft2(int m, int n, const float *px, float *py);
//{unsecret}
GENIAL_API void shfft2(int m, int nx, int ny, const float *px, float *py);
//{unsecret}
GENIAL_API void sfft2(int m, int n, const float *px, float *py);
//{unsecret}
GENIAL_API void zfft2(int m, int n, const double *px, double *py);
//{unsecret}
GENIAL_API void dhfft2(int m, int nx, int ny, const double *px, double *py);
//{unsecret}
GENIAL_API void dfft2(int m, int n, const double *px, double *py);

//{unsecret}
GENIAL_API void cifft(int n, const float  *p, float  *q);
//{unsecret}
GENIAL_API void sifft(int n, float  *p, float  *q);
//{unsecret}
GENIAL_API void zifft(int n, const double *p, double *q);
//{unsecret}
GENIAL_API void difft(int n, double *p, double *q);

//{unsecret}
GENIAL_API void cifft1(int m, int n ,         const float  *p, float  *q);
//{unsecret}
GENIAL_API void sifft1(int m, int nx, int ny, const float  *p, float  *q);
//{unsecret}
GENIAL_API void zifft1(int m, int n ,         const double *p, double *q);
//{unsecret}
GENIAL_API void difft1(int m, int nx, int ny, const double *p, double *q);

//{unsecret}
GENIAL_API void cifft2(int m, int n ,         const float  *p, float  *q);
//{unsecret}
GENIAL_API void sifft2(int m, int nx, int ny, const float  *p, float  *q);
//{unsecret}
GENIAL_API void zifft2(int m, int n ,         const double *p, double *q);
//{unsecret}
GENIAL_API void difft2(int m, int nx, int ny, const double *p, double *q);

//{unsecret}
GENIAL_API void cnfft(int n, const float *px, float *py);
//{unsecret}
GENIAL_API void shnfft(int nx, int ny, const float *px, float *py);
//{unsecret}
GENIAL_API void snfft(int n, const float *px, float *py);
//{unsecret}
GENIAL_API void znfft(int n, const double *px, double *py);
//{unsecret}
GENIAL_API void dhnfft(int nx, int ny, const double *px, double *py);
//{unsecret}
GENIAL_API void dnfft(int n, const double *px, double *py);

//{unsecret}
GENIAL_API void cafft(int n, const float *px, float *py);
//{unsecret}
GENIAL_API void shafft(int nx, int ny, const float *px, float *py);
//{unsecret}
GENIAL_API void safft(int n, const float *px, float *py);
//{unsecret}
GENIAL_API void zafft(int n, const double *px, double *py);
//{unsecret}
GENIAL_API void dhafft(int nx, int ny, const double *px, double *py);
//{unsecret}
GENIAL_API void dafft(int n, const double *px, double *py);

#endif // FFT_PRECOMPILE

#ifdef __cplusplus
}
#endif



#endif //FFT_H


