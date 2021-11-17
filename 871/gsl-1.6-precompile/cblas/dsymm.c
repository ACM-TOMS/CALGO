#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dsymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
             const enum CBLAS_UPLO Uplo, const int M, const int N,
             const MpIeee alpha, const MpIeee *A, const int lda,
             const MpIeee *B, const int ldb, const MpIeee beta, MpIeee *C,
             const int ldc)
{
#define BASE MpIeee
#include "source_symm_r.h"
#undef BASE
}
