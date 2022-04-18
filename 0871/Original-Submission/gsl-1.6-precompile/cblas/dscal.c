#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dscal (const int N, const MpIeee alpha, MpIeee *X, const int incX)
{
#define BASE MpIeee
#include "source_scal_r.h"
#undef BASE
}
