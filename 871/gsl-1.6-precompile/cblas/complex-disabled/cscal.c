#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_cscal (const int N, const void *alpha, void *X, const int incX)
{
#define BASE MpIeee
#include "source_scal_c.h"
#undef BASE
}
