#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_stbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
             const int N, const int K, const MpIeee *A, const int lda,
             MpIeee *X, const int incX)
{
#define BASE float
#include "source_tbmv_r.h"
#undef BASE
}
