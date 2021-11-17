#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

CBLAS_INDEX
cblas_isamax (const int N, const MpIeee *X, const int incX)
{
#define BASE MpIeee
#include "source_iamax_r.h"
#undef BASE
}
