#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_ssyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
             const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
             const MpIeee alpha, const MpIeee *A, const int lda, 
             const MpIeee beta, MpIeee *C, const int ldc)
{
#define BASE MpIeee
#include "source_syrk_r.h"
#undef BASE
}
