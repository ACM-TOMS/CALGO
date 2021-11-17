#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

static inline int
 cheb_eval_mode_e(const cheb_series * cs,
                 const MpIeee x,
                 gsl_mode_t mode,
                 gsl_sf_result * result)
{
  int  j;
  MpIeee d=  MpIeee( "0.0" );
  MpIeee dd=  MpIeee( "0.0" );

  MpIeee y=  (MpIeee( "2." )*x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  int  eval_order;

  if(GSL_MODE_PREC(mode) == GSL_PREC_DOUBLE)
    eval_order = cs->order;
  else
    eval_order = cs->order_sp;

  for(j = eval_order; j>=1; j--) {
    MpIeee temp=  d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  result->val = y*d - dd + 0.5 * cs->c[0];
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->c[eval_order]);
  return GSL_SUCCESS;
}
