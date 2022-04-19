#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* Compare the integral of f(x) with the integral of |f(x)|
   to determine if f(x) covers both positive and negative values */

static inline int
 test_positivity(MpIeee result, MpIeee resabs);

static inline int
 test_positivity(MpIeee result, MpIeee resabs)
{
  int  status=  (fabs (result) >= (1 - 50 * GSL_DBL_EPSILON) * resabs);

  return status;
}
