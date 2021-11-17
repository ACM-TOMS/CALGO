#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

MpIeee quadratic(MpIeee x, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  MpIeee a=  p->a;
  MpIeee b=  p->b;
  MpIeee c=  p->c;

  return (a * x + b) * x + c;
}

MpIeee quadratic_deriv(MpIeee x, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  MpIeee a=  p->a;
  MpIeee b=  p->b;
  MpIeee c=  p->c;

  return MpIeee( "2.0" ) * a * x + b;
}

void
quadratic_fdf (MpIeee x, void *params, 
               MpIeee *y, MpIeee *dy)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  MpIeee a=  p->a;
  MpIeee b=  p->b;
  MpIeee c=  p->c;

  *y = (a * x + b) * x + c;
  *dy = MpIeee( "2.0" ) * a * x + b;
}
