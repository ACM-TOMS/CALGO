#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_zcopy (const int N, const void *X, const int incX, void *Y,
             const int incY)
{
#define BASE MpIeee
#include "source_copy_c.h"
#undef BASE
}
