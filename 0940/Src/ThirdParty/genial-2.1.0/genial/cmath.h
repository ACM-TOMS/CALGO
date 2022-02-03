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

#ifndef CMATH_H
#define CMATH_H

#include <cmath>

#ifndef PI
#define PI 3.141592653589793238462643383279 
#endif 

template<class T> inline T positive_mod(T x, T y) { T m = x%y; return (m>=0)?m:m+y; }
template<class T> inline T mod(T x, T y) { return x%y; }

template<class T> inline T sqr(T x) { return x*x; }
template<class T> inline T inv(T x) { return T(1)/x; }

template<class T> inline T sign(T x) { return (x<0)?-1:1; }

//#if (defined(__ICL) || defined(_MSC_VER))
//inline float  abs(const float  &x) { return (x<0)?-x:x;}
//inline double abs(const double &x) { return (x<0)?-x:x;}
//#endif

inline float  norm(float  x) { return sqr(x); }
inline double norm(double x) { return sqr(x); }

inline float  round(float  x) { return (float)floor(double(x)+double(0.5)); }
inline double round(double x) { return        floor(double(x)+double(0.5)); }

inline float  exp10(float  x) { return exp(x*float (2.30258509299)); }
inline double exp10(double x) { return exp(x*double(2.30258509299)); }

template<class T> inline const T &sum (const T &x) { return x; }
template<class T> inline const T &flip(const T &x) { return x; }

inline int rand(int n) { return rand()%n; }

#endif
