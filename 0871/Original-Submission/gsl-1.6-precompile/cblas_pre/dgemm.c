#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
             const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
             const int K, const MpIeee alpha, const MpIeee *A, const int lda,
             const MpIeee *B, const int ldb, const MpIeee beta, MpIeee *C,
             const int ldc)
{
#define BASE double
#include "source_gemm_r.h"
#undef BASE
}
