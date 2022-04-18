#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/atanint.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_expint.h>

#include "chebyshev.h"
#include "cheb_eval.c"


static MpIeee atanint_data[21] =  {
  MpIeee( "1.91040361296235937512" ),
 -MpIeee( "0.4176351437656746940e-01" ),
  MpIeee( "0.275392550786367434e-02" ),
 -MpIeee( "0.25051809526248881e-03" ),
  MpIeee( "0.2666981285121171e-04" ),
 -MpIeee( "0.311890514107001e-05" ),
  MpIeee( "0.38833853132249e-06" ),
 -MpIeee( "0.5057274584964e-07" ),
  MpIeee( "0.681225282949e-08" ),
 -MpIeee( "0.94212561654e-09" ),
  MpIeee( "0.13307878816e-09" ),
 -MpIeee( "0.1912678075e-10" ),
  MpIeee( "0.278912620e-11" ),
 -MpIeee( "0.41174820e-12" ),
  MpIeee( "0.6142987e-13" ),
 -MpIeee( "0.924929e-14" ),
  MpIeee( "0.140387e-14" ),
 -MpIeee( "0.21460e-15" ),
  MpIeee( "0.3301e-16" ),
 -MpIeee( "0.511e-17" ),
  MpIeee( "0.79e-18" ),
};
static cheb_series atanint_cs = {
  atanint_data,
  20,
  -1, 1,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_atanint_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee ax=  fabs(x);
  const MpIeee sgn=  GSL_SIGN(x);

  /* CHECK_POINTER(result) */

  if(ax == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax < 0.5*GSL_SQRT_DBL_EPSILON) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax <= 1.0) {
    const MpIeee t=  2.0 * (x*x - 0.5);
    gsl_sf_result result_c;
    cheb_eval_e(&atanint_cs, t, &result_c);
    result->val  = x * result_c.val;
    result->err  = x * result_c.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(ax < 1.0/GSL_SQRT_DBL_EPSILON) {
    const MpIeee t=  2.0 * (1.0/(x*x) - 0.5);
    gsl_sf_result result_c;
    cheb_eval_e(&atanint_cs, t, &result_c);
    result->val  = sgn * (0.5*M_PI*log(ax) + result_c.val/ax);
    result->err  = result_c.err/ax + fabs(result->val*GSL_DBL_EPSILON);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = sgn * 0.5*M_PI*log(ax);
    result->err = 2.0 * fabs(result->val * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_atanint(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_atanint_e(x, &result));
}
