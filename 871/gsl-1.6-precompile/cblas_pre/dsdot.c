#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

MpIeee cblas_dsdot(const int N, const MpIeee *X, const int incX, const MpIeee *Y,
             const int incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  double
#define BASE float
#include "source_dot_r.h"
#undef ACC_TYPE
#undef BASE
#undef INIT_VAL
}
