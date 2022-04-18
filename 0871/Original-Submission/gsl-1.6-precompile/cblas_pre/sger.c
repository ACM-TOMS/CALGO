#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_sger (const enum CBLAS_ORDER order, const int M, const int N,
            const MpIeee alpha, const MpIeee *X, const int incX, const MpIeee *Y,
            const int incY, MpIeee *A, const int lda)
{
#define BASE float
#include "source_ger.h"
#undef BASE
}
