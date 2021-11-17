#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

CBLAS_INDEX
cblas_icamax (const int N, const void *X, const int incX)
{
#define BASE MpIeee
#include "source_iamax_c.h"
#undef BASE
}
