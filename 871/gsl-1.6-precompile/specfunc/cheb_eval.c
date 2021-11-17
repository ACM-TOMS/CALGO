#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"


static inline int
 cheb_eval_e(const cheb_series * cs,
            const MpIeee x,
            gsl_sf_result * result)
{
  int  j;
  MpIeee d=  MpIeee( "0.0" );
  MpIeee dd=  MpIeee( "0.0" );

  MpIeee y=  (MpIeee( "2.0" )*x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  MpIeee e=  MpIeee( "0.0" );

  for(j = cs->order; j>=1; j--) {
    MpIeee temp=  d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  { 
    MpIeee temp=  d;
    d = y*d - dd + MpIeee( "0.5" ) * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + MpIeee( "0.5" ) * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

