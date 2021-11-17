#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_csscal (const int N, const MpIeee alpha, void *X, const int incX)
{
#define BASE float
#include "source_scal_c_s.h"
#undef BASE
}
