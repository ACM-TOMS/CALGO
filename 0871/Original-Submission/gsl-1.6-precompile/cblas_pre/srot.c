#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_srot (const int N, MpIeee *X, const int incX, MpIeee *Y, const int incY,
            const MpIeee c, const MpIeee s)
{
#define BASE float
#include "source_rot.h"
#undef BASE
}
