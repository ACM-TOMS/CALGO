#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

static MpIeee rat_eval(const MpIeee a[], const size_t na,
          const MpIeee b[], const size_t nb, const MpIeee x)
{
  size_t i, j;
  MpIeee u;MpIeee  v;MpIeee  r;

  u = a[na - 1];

  for (i = na - 1; i > 0; i--)
    {
      u = x * u + a[i - 1];
    }

  v = b[nb - 1];

  for (j = nb - 1; j > 0; j--)
    {
      v = x * v + b[j - 1];
    }

  r = u / v;

  return r;
}
