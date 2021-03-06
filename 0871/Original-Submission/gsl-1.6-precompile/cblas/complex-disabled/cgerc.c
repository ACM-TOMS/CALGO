#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_cgerc (const enum CBLAS_ORDER order, const int M, const int N,
             const void *alpha, const void *X, const int incX, const void *Y,
             const int incY, void *A, const int lda)
{
#define BASE MpIeee
#include "source_gerc.h"
#undef BASE
}
