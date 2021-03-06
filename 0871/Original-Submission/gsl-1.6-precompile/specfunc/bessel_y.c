#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_y.c
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
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "bessel.h"
#include "bessel_olver.h"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* [Abramowitz+Stegun, 10.1.3]
 * with lmax=15, precision ~ 15D for x < 3
 *
 * checked OK [GJ] Wed May 13 15:41:25 MDT 1998 
 */
static int  bessel_yl_small_x(int  l, const MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result num_fact;
  MpIeee den=  gsl_sf_pow_int(x, l+1);
  int  stat_df=  gsl_sf_doublefact_e(2*l-1, &num_fact);

  if(stat_df != GSL_SUCCESS || den == 0.0) {
    OVERFLOW_ERROR(result);
  }
  else {
    const int lmax = 200;
    MpIeee t=  -MpIeee( "0.5" )*x*x;
    MpIeee sum=  MpIeee( "1.0" );
    MpIeee t_coeff=  MpIeee( "1.0" );
    MpIeee t_power=  MpIeee( "1.0" );
    MpIeee delta;
    int  i;
    for(i=1; i<=lmax; i++) {
      t_coeff /= i*(MpIeee( "2" )*(i-l) - MpIeee( "1" ));
      t_power *= t;
      delta = t_power*t_coeff;
      sum += delta;
      if(fabs(delta/sum) < 0.5*GSL_DBL_EPSILON) break;
    }
    result->val = -num_fact.val/den * sum;
    result->err = GSL_DBL_EPSILON * fabs(result->val);

    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int  gsl_sf_bessel_y0_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(1.0/GSL_DBL_MAX > 0.0 && x < 1.0/GSL_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else {
    gsl_sf_result cos_result;
    const int stat = gsl_sf_cos_e(x, &cos_result);
    result->val  = -cos_result.val/x;
    result->err  = fabs(cos_result.err/x);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
}


int  gsl_sf_bessel_y1_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 1.0/GSL_SQRT_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else if(x < 0.25) {
    const MpIeee y=  x*x;
    const MpIeee c1=   1.0/2.0;
    const MpIeee c2=  -1.0/8.0;
    const MpIeee c3=   1.0/144.0;
    const MpIeee c4=  -1.0/5760.0;
    const MpIeee c5=   1.0/403200.0;
    const MpIeee c6=  -1.0/43545600.0;
    const MpIeee sum=  1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*c6)))));
    result->val = -sum/y;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result cos_result;
    gsl_sf_result sin_result;
    const int stat_cos = gsl_sf_cos_e(x, &cos_result);
    const int stat_sin = gsl_sf_sin_e(x, &sin_result);
    const MpIeee cx=  cos_result.val;
    const MpIeee sx=  sin_result.val;
    result->val  = -(cx/x + sx)/x;
    result->err  = (fabs(cos_result.err/x) + sin_result.err)/fabs(x);
    result->err += GSL_DBL_EPSILON * (fabs(sx/x) + fabs(cx/(x*x)));
    return GSL_ERROR_SELECT_2(stat_cos, stat_sin);
  }
}


int  gsl_sf_bessel_y2_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 1.0/GSL_ROOT3_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else if(x < 0.5) {
    const MpIeee y=  x*x;
    const MpIeee c1=   1.0/6.0;
    const MpIeee c2=   1.0/24.0;
    const MpIeee c3=  -1.0/144.0;
    const MpIeee c4=   1.0/3456.0;
    const MpIeee c5=  -1.0/172800.0;
    const MpIeee c6=   1.0/14515200.0;
    const MpIeee c7=  -1.0/1828915200.0;
    const MpIeee sum=  1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
    result->val = -3.0/(x*x*x) * sum;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result cos_result;
    gsl_sf_result sin_result;
    const int stat_cos = gsl_sf_cos_e(x, &cos_result);
    const int stat_sin = gsl_sf_sin_e(x, &sin_result);
    const MpIeee sx=  sin_result.val;
    const MpIeee cx=  cos_result.val;
    const MpIeee a=  3.0/(x*x);
    result->val  = (1.0 - a)/x * cx - a * sx;
    result->err  = cos_result.err * fabs((1.0 - a)/x) + sin_result.err * fabs(a);
    result->err += GSL_DBL_EPSILON * (fabs(cx/x) + fabs(sx/(x*x)));
    return GSL_ERROR_SELECT_2(stat_cos, stat_sin);
  }
}


int  gsl_sf_bessel_yl_e(int  l, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(l < 0 || x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(l == 0) {
    return gsl_sf_bessel_y0_e(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_y1_e(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_y2_e(x, result);
  }
  else if(x < 3.0) {
    return bessel_yl_small_x(l, x, result);
  }
  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1.0)) {
    int  status=  gsl_sf_bessel_Ynu_asympx_e(l + 0.5, x, result);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val *= pre;
    result->err *= pre;
    return status;
  }
  else if(l > 40) {
    int  status=  gsl_sf_bessel_Ynu_asymp_Olver_e(l + 0.5, x, result);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val *= pre;
    result->err *= pre;
    return status;
  }
  else {
    /* recurse upward */
    gsl_sf_result r_by;
    gsl_sf_result r_bym;
    int  stat_1=  gsl_sf_bessel_y1_e(x, &r_by);
    int  stat_0=  gsl_sf_bessel_y0_e(x, &r_bym);
    MpIeee bym=  r_bym.val;
    MpIeee by=  r_by.val;
    MpIeee byp;
    int  j;
    for(j=1; j<l; j++) { 
      byp = (MpIeee( "2" )*j+MpIeee( "1" ))/x*by - bym;
      bym = by;
      by  = byp;
    }
    result->val = by;
    result->err = fabs(result->val) * (GSL_DBL_EPSILON + fabs(r_by.err/r_by.val) + fabs(r_bym.err/r_bym.val));

    return GSL_ERROR_SELECT_2(stat_1, stat_0);
  }
}


int  gsl_sf_bessel_yl_array(const int lmax, const MpIeee x, MpIeee * result_array)
{
  /* CHECK_POINTER(result_array) */

  if(lmax < 1 || x <= 0.0) {
    GSL_ERROR ("error", GSL_EDOM);
  }
  else {
    gsl_sf_result r_yell;
    gsl_sf_result r_yellm1;
    int  stat_1=  gsl_sf_bessel_y1_e(x, &r_yell);
    int  stat_0=  gsl_sf_bessel_y0_e(x, &r_yellm1);
    MpIeee yellp1;
    MpIeee yell=  r_yell.val;
    MpIeee yellm1=  r_yellm1.val;
    int  ell;

    result_array[0] = yellm1;
    result_array[1] = yell;

    for(ell = 1; ell < lmax; ell++) {
      yellp1 = (MpIeee( "2" )*ell+MpIeee( "1" ))/x * yell - yellm1;
      result_array[ell+1] = yellp1;
      yellm1 = yell;
      yell   = yellp1;
    }

    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_y0(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_y0_e(x, &result));
}

MpIeee gsl_sf_bessel_y1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_y1_e(x, &result));
}

MpIeee gsl_sf_bessel_y2(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_y2_e(x, &result));
}

MpIeee gsl_sf_bessel_yl(const int l, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_yl_e(l, x, &result));
}


