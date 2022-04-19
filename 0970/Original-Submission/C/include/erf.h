#ifndef ERF_H
#define	ERF_H
#include <math.h>

#if _MSC_VER <= 1700
/*
 * erf() and erfc() are not included into MSVC math library while they are part of GNU libm.
 * To compile STS using MSVC I had to add implementation of these functions below. Both
 * are taken from Numerical recipes in C.
 */

double erf(double x);
double erfc(double x);
#endif
#endif