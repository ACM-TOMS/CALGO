#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/lambert.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_lambert.h>

/* Started with code donated by K. Briggs; added
 * error estimates, GSL foo, and minor tweaks.
 * Some Lambert-ology from
 *  [Corless, Gonnet, Hare, and Jeffrey, "On Lambert's W Function".]
 */


/* Halley iteration (eqn. 5.12, Corless et al) */
static int
 halley_iteration(
  MpIeee x,
  MpIeee w_initial,
  unsigned int  max_iters,
  gsl_sf_result * result
  )
{
  MpIeee w=  w_initial;
  unsigned int  i;

  for(i=0; i<max_iters; i++) {
    MpIeee tol;
    const MpIeee e=  exp(w);
    const MpIeee p=  w + 1.0;
    MpIeee t=  w*e - x;
    /* printf("FOO: %20.16g  %20.16g\n", w, t); */

    if (w > MpIeee( "0" )) {
      t = (t/p)/e;  /* Newton iteration */
    } else {
      t /= e*p - MpIeee( "0.5" )*(p + MpIeee( "1.0" ))*t/p;  /* Halley iteration */
    };

    w -= t;

    tol = GSL_DBL_EPSILON * GSL_MAX_DBL(fabs(w), MpIeee( "1.0" )/(fabs(p)*e));

    if(fabs(t) < tol)
    {
      result->val = w;
      result->err = 2.0*tol;
      return GSL_SUCCESS;
    }
  }

  /* should never get here */
  result->val = w;
  result->err = fabs(w);
  return GSL_EMAXITER;
}


/* series which appears for q near zero;
 * only the argument is different for the different branches
 */
static MpIeee series_eval(MpIeee r)
{
  static const MpIeee c[12] =  {
    -1.0,
     2.331643981597124203363536062168,
    -1.812187885639363490240191647568,
     1.936631114492359755363277457668,
    -2.353551201881614516821543561516,
     3.066858901050631912893148922704,
    -4.175335600258177138854984177460,
     5.858023729874774148815053846119,
    -8.401032217523977370984161688514,
     12.250753501314460424,
    -18.100697012472442755,
     27.029044799010561650
  };
  const MpIeee t_8=  c[8] + r*(c[9] + r*(c[10] + r*c[11]));
  const MpIeee t_5=  c[5] + r*(c[6] + r*(c[7]  + r*t_8));
  const MpIeee t_1=  c[1] + r*(c[2] + r*(c[3]  + r*(c[4] + r*t_5)));
  return c[0] + r*t_1;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_lambert_W0_e(MpIeee x, gsl_sf_result * result)
{
  const MpIeee one_over_E=  1.0/M_E;
  const MpIeee q=  x + one_over_E;

  if(x == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(q < 0.0) {
    /* Strictly speaking this is an error. But because of the
     * arithmetic operation connecting x and q, I am a little
     * lenient in case of some epsilon overshoot. The following
     * answer is quite accurate in that case. Anyway, we have
     * to return GSL_EDOM.
     */
    result->val = -1.0;
    result->err =  sqrt(-q);
    return GSL_EDOM;
  }
  else if(q == 0.0) {
    result->val = -1.0;
    result->err =  GSL_DBL_EPSILON; /* cannot error is zero, maybe q == 0 by "accident" */
    return GSL_SUCCESS;
  }
  else if(q < 1.0e-03) {
    /* series near -1/E in sqrt(q) */
    const MpIeee r=  sqrt(q);
    result->val = series_eval(r);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    static const unsigned int  MAX_ITERS=  10;
    MpIeee w;

    if (x < MpIeee( "1.0" )) {
      /* obtain initial approximation from series near x=0;
       * no need for extra care, since the Halley iteration
       * converges nicely on this branch
       */
      const MpIeee p=  sqrt(2.0 * M_E * q);
      w = -MpIeee( "1.0" ) + p*(MpIeee( "1.0" ) + p*(-MpIeee( "1.0" )/MpIeee( "3.0" ) + p*MpIeee( "11.0" )/MpIeee( "72.0" ))); 
    }
    else {
      /* obtain initial approximation from rough asymptotic */
      w = log(x);
      if(x > MpIeee( "3.0" )) w -= log(w);
    }

    return halley_iteration(x, w, MAX_ITERS, result);
  }
}


int
 gsl_sf_lambert_Wm1_e(MpIeee x, gsl_sf_result * result)
{
  if(x > MpIeee( "0.0" )) {
    return gsl_sf_lambert_W0_e(x, result);
  }
  else if(x == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    static const unsigned int  MAX_ITERS=  32;
    const MpIeee one_over_E=  1.0/M_E;
    const MpIeee q=  x + one_over_E;
    MpIeee w;

    if (q < 0.0) {
      /* As in the W0 branch above, return some reasonable answer anyway. */
      result->val = -1.0; 
      result->err =  sqrt(-q);
      return GSL_EDOM;
    }

    if(x < -MpIeee( "1.0e-6" )) {
      /* Obtain initial approximation from series about q = 0,
       * as long as we're not very close to x = 0.
       * Use full series and try to bail out if q is too small,
       * since the Halley iteration has bad convergence properties
       * in finite arithmetic for q very small, because the
       * increment alternates and p is near zero.
       */
      const MpIeee r=  -sqrt(q);
      w = series_eval(r);
      if(q < 3.0e-3) {
        /* this approximation is good enough */
        result->val = w;
        result->err = 5.0 * GSL_DBL_EPSILON * fabs(w);
        return GSL_SUCCESS;
      }
    }
    else {
      /* Obtain initial approximation from asymptotic near zero. */
      const MpIeee L_1=  log(-x);
      const MpIeee L_2=  log(-L_1);
      w = L_1 - L_2 + L_2/L_1;
    }

    return halley_iteration(x, w, MAX_ITERS, result);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_lambert_W0(MpIeee x)
{
  EVAL_RESULT(gsl_sf_lambert_W0_e(x, &result));
}

MpIeee gsl_sf_lambert_Wm1(MpIeee x)
{
  EVAL_RESULT(gsl_sf_lambert_Wm1_e(x, &result));
}
