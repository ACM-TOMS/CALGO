#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_k.c
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
#include <gsl/gsl_sf_bessel.h>

#include "error.h"
#include "check.h"

#include "bessel.h"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* [Abramowitz+Stegun, 10.2.4 + 10.2.6]
 * with lmax=15, precision ~ 15D for x < 3
 *
 * assumes l >= 1
 */
static int  bessel_kl_scaled_small_x(int  l, const MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result num_fact;
  MpIeee den=  gsl_sf_pow_int(x, l+1);
  int  stat_df=  gsl_sf_doublefact_e((unsigned int) (2*l-1), &num_fact);

  if(stat_df != GSL_SUCCESS || den == 0.0) {
    OVERFLOW_ERROR(result);
  }
  else {
    const int lmax = 50;
    gsl_sf_result ipos_term;
    MpIeee ineg_term;
    MpIeee sgn=  (GSL_IS_ODD(l) ? -MpIeee( "1.0" ) : MpIeee( "1.0" ));
    MpIeee ex=  exp(x);
    MpIeee t=  MpIeee( "0.5" )*x*x;
    MpIeee sum=  MpIeee( "1.0" );
    MpIeee t_coeff=  MpIeee( "1.0" );
    MpIeee t_power=  MpIeee( "1.0" );
    MpIeee delta;
    int  stat_il;
    int  i;

    for(i=1; i<lmax; i++) {
      t_coeff /= i*(MpIeee( "2" )*(i-l) - MpIeee( "1" ));
      t_power *= t;
      delta = t_power*t_coeff;
      sum += delta;
      if(fabs(delta/sum) < GSL_DBL_EPSILON) break;
    }

    stat_il = gsl_sf_bessel_il_scaled_e(l, x, &ipos_term);
    ineg_term =  sgn * num_fact.val/den * sum;
    result->val = -sgn * 0.5*M_PI * (ex*ipos_term.val - ineg_term);
    result->val *= ex;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_il;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_bessel_k0_scaled_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = M_PI/(2.0*x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_bessel_k1_scaled_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < (M_SQRTPI+1.0)/(M_SQRT2*GSL_SQRT_DBL_MAX)) {
    OVERFLOW_ERROR(result);
  }
  else {
    result->val = M_PI/(2.0*x) * (1.0 + 1.0/x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_bessel_k2_scaled_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0/GSL_ROOT3_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else {
    result->val = M_PI/(2.0*x) * (1.0 + 3.0/x * (1.0 + 1.0/x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_bessel_kl_scaled_e(int  l, const MpIeee x, gsl_sf_result * result)
{
  if(l < 0 || x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(l == 0) {
    return gsl_sf_bessel_k0_scaled_e(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_k1_scaled_e(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_k2_scaled_e(x, result);
  }
  else if(x < 3.0) {
    return bessel_kl_scaled_small_x(l, x, result);
  }
  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1)) {
    int  status=  gsl_sf_bessel_Knu_scaled_asympx_e(l + 0.5, x, result);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val *= pre;
    result->err *= pre;
    return status;
  }
  else if(GSL_MIN_DBL(MpIeee("0.29")/(l*l+1.0), MpIeee("0.5")/(l*l+1.0+x*x)) < GSL_ROOT3_DBL_EPSILON) {
    int  status=  gsl_sf_bessel_Knu_scaled_asymp_unif_e(l + 0.5, x, result);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val *= pre;
    result->err *= pre;
    return status;
  }
  else {
    /* recurse upward */
    gsl_sf_result r_bk;
    gsl_sf_result r_bkm;
    int  stat_1=  gsl_sf_bessel_k1_scaled_e(x, &r_bk);
    int  stat_0=  gsl_sf_bessel_k0_scaled_e(x, &r_bkm);
    MpIeee bkp;
    MpIeee bk=  r_bk.val;
    MpIeee bkm=  r_bkm.val;
    int  j;
    for(j=1; j<l; j++) { 
      bkp = (MpIeee( "2" )*j+MpIeee( "1" ))/x*bk + bkm;
      bkm = bk;
      bk  = bkp;
    }
    result->val  = bk;
    result->err  = fabs(bk) * (fabs(r_bk.err/r_bk.val) + fabs(r_bkm.err/r_bkm.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_1, stat_0);
  }
}

int  gsl_sf_bessel_kl_scaled_array(const int lmax, const MpIeee x, MpIeee * result_array)
{
  if(lmax < 1 || x <= 0.0) {
    GSL_ERROR("domain error", GSL_EDOM);
  }
  else {
    int  ell;
    MpIeee kellp1;MpIeee  kell;MpIeee  kellm1;
    gsl_sf_result r_kell;
    gsl_sf_result r_kellm1;
    gsl_sf_bessel_k1_scaled_e(x, &r_kell);
    gsl_sf_bessel_k0_scaled_e(x, &r_kellm1);
    kell   = r_kell.val;
    kellm1 = r_kellm1.val;
    result_array[0] = kellm1;
    result_array[1] = kell;
    for(ell = 1; ell < lmax; ell++) {
      kellp1 = (MpIeee( "2" )*ell+MpIeee( "1" ))/x * kell + kellm1;
      result_array[ell+1] = kellp1;
      kellm1 = kell;
      kell   = kellp1;
    }
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_k0_scaled(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_k0_scaled_e(x, &result));
}

MpIeee gsl_sf_bessel_k1_scaled(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_k1_scaled_e(x, &result));
}

MpIeee gsl_sf_bessel_k2_scaled(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_k2_scaled_e(x, &result));
}

MpIeee gsl_sf_bessel_kl_scaled(const int l, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_kl_scaled_e(l, x, &result));
}


