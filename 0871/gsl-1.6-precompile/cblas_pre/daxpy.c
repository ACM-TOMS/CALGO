#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_daxpy (const int N, const MpIeee alpha, const MpIeee *X, const int incX,
             MpIeee *Y, const int incY)
{
#define BASE double
#include "source_axpy_r.h"
#undef BASE
}
