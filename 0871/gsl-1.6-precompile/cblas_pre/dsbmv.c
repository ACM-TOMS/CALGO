#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dsbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const int N, const int K, const MpIeee alpha, const MpIeee *A,
             const int lda, const MpIeee *X, const int incX,
             const MpIeee beta, MpIeee *Y, const int incY)
{
#define BASE double
#include "source_sbmv.h"
#undef BASE
}
