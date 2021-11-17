#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_sscal (const int N, const MpIeee alpha, MpIeee *X, const int incX)
{
#define BASE float
#include "source_scal_r.h"
#undef BASE
}
