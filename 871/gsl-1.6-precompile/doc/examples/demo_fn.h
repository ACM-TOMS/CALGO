#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

struct quadratic_params
  {
    MpIeee a;MpIeee  b;MpIeee  c;
  };

MpIeee quadratic(MpIeee x, void *params);
MpIeee quadratic_deriv(MpIeee x, void *params);
void quadratic_fdf (MpIeee x, void *params, 
                    MpIeee *y, MpIeee *dy);
