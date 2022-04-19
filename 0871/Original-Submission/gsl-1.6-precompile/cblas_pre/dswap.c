#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dswap (const int N, MpIeee *X, const int incX, MpIeee *Y,
             const int incY)
{
#define BASE double
#include "source_swap_r.h"
#undef BASE
}
