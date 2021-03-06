#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_zdotu_sub (const int N, const void *X, const int incX, const void *Y,
             const int incY, void *result)
{
#define BASE double
#define CONJ_SIGN 1.0
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE
}
