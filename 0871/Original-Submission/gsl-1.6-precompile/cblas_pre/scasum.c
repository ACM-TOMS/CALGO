#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

MpIeee cblas_scasum(const int N, const void *X, const int incX)
{
#define BASE float
#include "source_asum_c.h"
#undef BASE
}
