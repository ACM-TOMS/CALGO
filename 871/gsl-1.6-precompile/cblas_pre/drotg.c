#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_drotg (MpIeee *a, MpIeee *b, MpIeee *c, MpIeee *s)
{
#define BASE double
#include "source_rotg.h"
#undef BASE
}
