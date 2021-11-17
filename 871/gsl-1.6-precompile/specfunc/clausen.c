#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/clausen.c
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
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_clausen.h>

#include "chebyshev.h"
#include "cheb_eval.c"

static MpIeee aclaus_data[15] =  {
  MpIeee( "2.142694363766688447e+00" ),
  MpIeee( "0.723324281221257925e-01" ),
  MpIeee( "0.101642475021151164e-02" ),
  MpIeee( "0.3245250328531645e-04" ),
  MpIeee( "0.133315187571472e-05" ),
  MpIeee( "0.6213240591653e-07" ),
  MpIeee( "0.313004135337e-08" ),
  MpIeee( "0.16635723056e-09" ),
  MpIeee( "0.919659293e-11" ),
  MpIeee( "0.52400462e-12" ),
  MpIeee( "0.3058040e-13" ),
  MpIeee( "0.18197e-14" ),
  MpIeee( "0.1100e-15" ),
  MpIeee( "0.68e-17" ),
  MpIeee( "0.4e-18" )
};
static cheb_series aclaus_cs = {
  aclaus_data,
  14,
  -1, 1,
  8  /* FIXME:  this is a guess, correct value needed here BJG */
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_clausen_e(MpIeee x, gsl_sf_result *result)
{
  const MpIeee x_cut=  M_PI * GSL_SQRT_DBL_EPSILON;

  MpIeee sgn=  MpIeee( "1.0" );
  int  status_red;

  if(x < MpIeee( "0.0" )) {
    x   = -x;
    sgn = -MpIeee( "1.0" );
  }

  /* Argument reduction to [0, 2pi) */
  status_red = gsl_sf_angle_restrict_pos_e(&x);

  /* Further reduction to [0,pi) */
  if(x > M_PI) {
    /* simulated extra precision: 2PI = p0 + p1 */
    const MpIeee p0=  6.28125;
    const MpIeee p1=  0.19353071795864769253e-02;
    x = (p0 - x) + p1;
    sgn = -sgn;
  }

  if(x == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
  }
  else if(x < x_cut) {
    result->val = x * (1.0 - log(x));
    result->err = x * GSL_DBL_EPSILON;
  }
  else {
    const MpIeee t=  2.0*(x*x / (M_PI*M_PI) - 0.5);
    gsl_sf_result result_c;
    cheb_eval_e(&aclaus_cs, t, &result_c);
    result->val = x * (result_c.val - log(x));
    result->err = x * (result_c.err + GSL_DBL_EPSILON);
  }

  result->val *= sgn;

  return status_red;
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_clausen(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_clausen_e(x, &result));
}
