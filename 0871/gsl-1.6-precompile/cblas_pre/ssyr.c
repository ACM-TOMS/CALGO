#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_ssyr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
            const int N, const MpIeee alpha, const MpIeee *X, const int incX,
            MpIeee *A, const int lda)
{
#define BASE float
#include "source_syr.h"
#undef BASE
}
