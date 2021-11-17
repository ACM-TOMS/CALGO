//-*- C++ -*-
#ifndef COMPLEX_TYPE_H_
#define COMPLEX_TYPE_H_

#include <complex>

#if defined(__GNUC__) || defined(__sgi)
typedef std::complex<double> Complex;
#elif defined(__DECCXX) || defined(__xlC__)
typedef complex Complex;
#else
typedef std::complex<double> Complex;
#endif

#endif // COMPLEX_TYPE_H_
