#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_zhemm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
             const enum CBLAS_UPLO Uplo, const int M, const int N,
             const void *alpha, const void *A, const int lda, const void *B,
             const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE MpIeee
#include "source_hemm.h"
#undef BASE
}
