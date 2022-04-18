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

#ifndef COMPLEX_H
#define COMPLEX_H

#include "cmath.h"
#include <complex>

namespace std
{

template<class T>
struct complex_traits
{
  typedef T value_type;
};

template<class T> struct complex_traits<      complex<T> > { typedef T value_type; };
template<class T> struct complex_traits<const complex<T> > { typedef T value_type; };


#ifdef _STLPORT_VERSION
template <class T> inline T &real(complex<T> &z) { return z._M_re; }
template <class T> inline T &imag(complex<T> &z) { return z._M_im; }
#endif

inline float        &conj(float        &x) { return x; }
inline const float  &conj(const float  &x) { return x; }
inline double       &conj(double       &x) { return x; }
inline const double &conj(const double &x) { return x; }

inline float        &real(float        &x) { return x; }
inline const float  &real(const float  &x) { return x; }
inline double       &real(double       &x) { return x; }
inline const double &real(const double &x) { return x; }

//#ifndef SSE
//inline float  abs(const complex<float > &z) { return sqrt(norm(z)); }
//#else
//inline float  abs(const complex<float > &z) { return abs(m128cf(z))[0]; }
//#endif
//#ifndef SSE2
//inline double abs(const complex<double> &z) { return sqrt(norm(z)); }
//#else
//inline double abs(const complex<double> &z) { return abs(m128cd(z))[0]; }
//#endif

inline float  abs(float  x, float  y) { return abs(complex<float >(x,y)); }
inline double abs(double x, double y) { return abs(complex<double>(x,y)); }

inline float  norm(float  x, float  y) { return norm(complex<float >(x,y)); }
inline double norm(double x, double y) { return norm(complex<double>(x,y)); }

template<class T> inline T first (const complex<T>& z) { return real(z); }
template<class T> inline T second(const complex<T>& z) { return imag(z); }

template<class V> inline complex<V> round(const complex<V> &x) { return complex<V>(::round(x.real()),::round(x.imag())); }

template<class V> inline complex<V> flip_ri(complex<V> z) { return complex<V>(z.imag(),z.real()); }

template<class V> inline complex<V> addsub(const complex<V> &x,const complex<V> &y) { return complex<V>(x.real()+y.real(), x.imag()-y.imag()); }
template<class V> inline complex<V> subadd(const complex<V> &x,const complex<V> &y) { return complex<V>(x.real()-y.real(), x.imag()+y.imag()); }

template<class V> inline complex<V> imul(    const complex<V> &z) { return complex<V>(  -z.imag(),  z.real()); }
template<class V> inline complex<V> imul(V a,const complex<V> &z) { return complex<V>(-a*z.imag(),a*z.real()); }

template<class V> inline complex<V> add_imul (const complex<V> &x,const complex<V> &y) { return complex<V>(x.real()-y.imag(), x.imag()+y.real()); }
template<class V> inline complex<V> sub_imul (const complex<V> &x,const complex<V> &y) { return complex<V>(x.real()+y.imag(), x.imag()-y.real()); }

template<class V> inline complex<V> mimul(    const complex<V> &z) { return complex<V>(   z.imag(),  z.real()); }
template<class V> inline complex<V> mimul(V a,const complex<V> &z) { return complex<V>( a*z.imag(),a*z.real()); }

template<class V> inline V cmul (const V &x,const V &y) { return conj(x)*y; }
template<class V> inline complex<V> cmul (const complex<V> &x,const complex<V> &y) { return complex<V>( real(x)*real(y)+imag(x)*imag(y),real(x)*imag(y)-imag(x)*real(y)); }
template<class V> inline complex<V> mimul(const complex<V> &x,const complex<V> &y) { return complex<V>( imag(x)*real(y)+real(x)*imag(y),imag(x)*imag(y)-real(x)*real(y)); }
template<class V> inline complex<V> icmul(const complex<V> &x,const complex<V> &y) { return complex<V>( imag(x)*real(y)-real(x)*imag(y),real(x)*real(y)+imag(x)*imag(y)); }
template<class V> inline complex<V> imul (const complex<V> &x,const complex<V> &y) { return complex<V>(-real(x)*imag(y)-imag(x)*real(y),real(x)*real(y)-imag(x)*imag(y)); }

template<class V> inline void load  (complex<V> &x, complex<float> *p) { x=*p; }
template<class V> inline void loadu (complex<V> &x, complex<float> *p) { x=*p; }


}

#endif

