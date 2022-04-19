//-*- C++ -*-

#ifndef MATHFUNC_H
#define MATHFUNC_H

#include "cmplxtype.h"
#include <math.h>

inline double Abs(double x) { return ( fabs(x) ); }
inline double Abs(Complex x) { return ( abs(x) ); }

inline double Norm(double x) { return ( x * x ); }
inline double Norm(Complex x) { return norm(x); }

inline double Conj(double x) { return x; }
inline Complex Conj(Complex x) { return conj(x); }

inline double Real(double x) { return x; }
inline double Real(Complex x) { return real(x); }

#endif // MATHFUNC_H
