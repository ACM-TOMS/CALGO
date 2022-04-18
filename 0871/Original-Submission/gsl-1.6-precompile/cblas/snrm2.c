#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

MpIeee cblas_snrm2(const int N, const MpIeee *X, const int incX)
{
#define BASE MpIeee
#include "source_nrm2_r.h"
#undef BASE
}
