#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
             const int M, const int N, const MpIeee alpha, const MpIeee *A,
             const int lda, const MpIeee *X, const int incX,
             const MpIeee beta, MpIeee *Y, const int incY)
{
#define BASE MpIeee
#include "source_gemv_r.h"
#undef BASE
}
