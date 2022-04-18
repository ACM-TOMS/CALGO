#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_srotmg (MpIeee *d1, MpIeee *d2, MpIeee *b1, const MpIeee b2, MpIeee *P)
{
#define BASE MpIeee
#include "source_rotmg.h"
#undef BASE
}
