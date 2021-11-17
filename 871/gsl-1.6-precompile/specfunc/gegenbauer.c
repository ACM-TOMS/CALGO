#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gegenbauer.c
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
#include <gsl/gsl_sf_gegenbauer.h>

#include "error.h"

/* See: [Thompson, Atlas for Computing Mathematical Functions] */


int
 gsl_sf_gegenpoly_1_e(MpIeee lambda, MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(lambda == MpIeee( "0.0" )) {
    result->val = 2.0*x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 2.0*lambda*x;
    result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
 gsl_sf_gegenpoly_2_e(MpIeee lambda, MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(lambda == MpIeee( "0.0" )) {
    const MpIeee txx=  2.0*x*x;
    result->val  = -1.0 + txx;
    result->err  = 2.0 * GSL_DBL_EPSILON * fabs(txx);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = lambda*(-1.0 + 2.0*(1.0+lambda)*x*x);
    result->err = GSL_DBL_EPSILON * (2.0 * fabs(result->val) + fabs(lambda));
    return GSL_SUCCESS;
  }
}

int
 gsl_sf_gegenpoly_3_e(MpIeee lambda, MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(lambda == MpIeee( "0.0" )) {
    result->val = x*(-2.0 + 4.0/3.0*x*x);
    result->err = GSL_DBL_EPSILON * (2.0 * fabs(result->val) + fabs(x));
    return GSL_SUCCESS;
  }
  else {
    MpIeee c=  MpIeee( "4.0" ) + lambda*(MpIeee( "6.0" ) + MpIeee( "2.0" )*lambda);
    result->val = 2.0*lambda * x * ( -1.0 - lambda + c*x*x/3.0 );
    result->err = GSL_DBL_EPSILON * (2.0 * fabs(result->val) + fabs(lambda * x));
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_gegenpoly_n_e(int  n, MpIeee lambda, MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(lambda <= -MpIeee( "0.5" ) || n < MpIeee( "0" )) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    return gsl_sf_gegenpoly_1_e(lambda, x, result);
  }
  else if(n == 2) {
    return gsl_sf_gegenpoly_2_e(lambda, x, result);
  }
  else if(n == 3) {
    return gsl_sf_gegenpoly_3_e(lambda, x, result);
  }
  else {
    if(lambda == MpIeee( "0.0" ) && (x >= -MpIeee( "1.0" ) || x <= MpIeee( "1.0" ))) {
      /* 2 T_n(x)/n */
      const MpIeee z=  n * acos(x);
      result->val = 2.0 * cos(z) / n;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(z * result->val);
      return GSL_SUCCESS;
    }
    else {
      int  k;
      gsl_sf_result g2;
      gsl_sf_result g3;
      int  stat_g2=  gsl_sf_gegenpoly_2_e(lambda, x, &g2);
      int  stat_g3=  gsl_sf_gegenpoly_3_e(lambda, x, &g3);
      int  stat_g=  GSL_ERROR_SELECT_2(stat_g2, stat_g3);
      MpIeee gkm2=  g2.val;
      MpIeee gkm1=  g3.val;
      MpIeee gk=  MpIeee( "0.0" );
      for(k=4; k<=n; k++) {
        gk = (MpIeee( "2.0" )*(k+lambda-MpIeee( "1.0" ))*x*gkm1 - (k+MpIeee( "2.0" )*lambda-MpIeee( "2.0" ))*gkm2) / k;
        gkm2 = gkm1;
        gkm1 = gk;
      }
      result->val = gk;
      result->err = 2.0 * GSL_DBL_EPSILON * 0.5 * n * fabs(gk);
      return stat_g;
    }
  }
}


int
 gsl_sf_gegenpoly_array(int  nmax, MpIeee lambda, MpIeee x, MpIeee * result_array)
{
  int  k;

  /* CHECK_POINTER(result_array) */

  if(lambda <= -MpIeee( "0.5" ) || nmax < MpIeee( "0" )) {
    GSL_ERROR("domain error", GSL_EDOM);
  }

  /* n == 0 */
  result_array[0] = MpIeee( "1.0" );
  if(nmax == 0) return GSL_SUCCESS;

  /* n == 1 */
  if(lambda == MpIeee( "0.0" ))
    result_array[1] = MpIeee( "2.0" )*x;
  else
    result_array[1] = MpIeee( "2.0" )*lambda*x;

  /* n <= nmax */
  for(k=2; k<=nmax; k++) {
    MpIeee term1=  MpIeee( "2.0" )*(k+lambda-MpIeee( "1.0" )) * x * result_array[k-1];
    MpIeee term2=  (k+MpIeee( "2.0" )*lambda-MpIeee( "2.0" ))     * result_array[k-2];
    result_array[k] = (term1 - term2) / k;
  }
  
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_gegenpoly_1(MpIeee lambda, MpIeee x)
{
  EVAL_RESULT(gsl_sf_gegenpoly_1_e(lambda, x, &result));
}

MpIeee gsl_sf_gegenpoly_2(MpIeee lambda, MpIeee x)
{
  EVAL_RESULT(gsl_sf_gegenpoly_2_e(lambda, x, &result));
}

MpIeee gsl_sf_gegenpoly_3(MpIeee lambda, MpIeee x)
{
  EVAL_RESULT(gsl_sf_gegenpoly_3_e(lambda, x, &result));
}

MpIeee gsl_sf_gegenpoly_n(int  n, MpIeee lambda, MpIeee x)
{
  EVAL_RESULT(gsl_sf_gegenpoly_n_e(n, lambda, x, &result));
}
