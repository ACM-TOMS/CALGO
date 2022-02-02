#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/laguerre.c
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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>

#include "error.h"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on the large 2b-4a asymptotic for 1F1
 * [Abramowitz+Stegun, 13.5.21]
 */
static
int
 laguerre_large_n(const int n, const MpIeee alpha, const MpIeee x,
                 gsl_sf_result * result)
{
  const MpIeee a=  -n;
  const MpIeee b=  alpha + 1.0;
  const MpIeee eta=  2.0*b - 4.0*a;
  const MpIeee cos2th=  x/eta;
  const MpIeee sin2th=  1.0 - cos2th;
  const MpIeee th=  acos(sqrt(cos2th));
  const MpIeee pre_h=  0.25*M_PI*M_PI*eta*eta*cos2th*sin2th;
  gsl_sf_result lg_b;
  gsl_sf_result lnfact;
  int  stat_lg=  gsl_sf_lngamma_e(b+n, &lg_b);
  int  stat_lf=  gsl_sf_lnfact_e(n, &lnfact);
  MpIeee pre_term1=  MpIeee( "0.5" )*(MpIeee( "1.0" )-b)*log(MpIeee( "0.25" )*x*eta);
  MpIeee pre_term2=  MpIeee( "0.25" )*log(pre_h);
  MpIeee lnpre_val=  lg_b.val - lnfact.val + MpIeee( "0.5" )*x + pre_term1 - pre_term2;
  MpIeee lnpre_err=  lg_b.err + lnfact.err + GSL_DBL_EPSILON * (fabs(pre_term1)+fabs(pre_term2));
  MpIeee ser_term1=  sin(a*M_PI);
  MpIeee ser_term2=  sin(MpIeee( "0.25" )*eta*(MpIeee( "2.0" )*th - sin(MpIeee( "2.0" )*th)) + MpIeee( "0.25" )*M_PI);
  MpIeee ser_val=  ser_term1 + ser_term2;
  MpIeee ser_err=  GSL_DBL_EPSILON * (fabs(ser_term1) + fabs(ser_term2));
  int  stat_e=  gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err, ser_val, ser_err, result);
  result->err += 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
  return GSL_ERROR_SELECT_3(stat_e, stat_lf, stat_lg);
}


/* Evaluate polynomial based on confluent hypergeometric representation.
 *
 * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
 *
 * assumes n > 0 and a != negative integer greater than -n
 */
static
int
 laguerre_n_cp(const int n, const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result lnfact;
  gsl_sf_result lg1;
  gsl_sf_result lg2;
  MpIeee s1;MpIeee  s2;
  int  stat_f=  gsl_sf_lnfact_e(n, &lnfact);
  int  stat_g1=  gsl_sf_lngamma_sgn_e(a+1.0+n, &lg1, &s1);
  int  stat_g2=  gsl_sf_lngamma_sgn_e(a+1.0, &lg2, &s2);
  MpIeee poly_1F1_val=  MpIeee( "1.0" );
  MpIeee poly_1F1_err=  MpIeee( "0.0" );
  int  stat_e;
  int  k;

  MpIeee lnpre_val=  (lg1.val - lg2.val) - lnfact.val;
  MpIeee lnpre_err=  lg1.err + lg2.err + lnfact.err + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(lnpre_val);

  for(k=n-1; k>=0; k--) {
    MpIeee t=  (-n+k)/(a+MpIeee( "1.0" )+k) * (x/(k+MpIeee( "1" )));
    MpIeee r=  t + MpIeee( "1.0" )/poly_1F1_val;
    if(r > MpIeee( "0.9" )*GSL_DBL_MAX/poly_1F1_val) {
      /* internal error only, don't call the error handler */
      INTERNAL_OVERFLOW_ERROR(result);
    }
    else {
      /* Collect the Horner terms. */
      poly_1F1_val  = MpIeee( "1.0" ) + t * poly_1F1_val;
      poly_1F1_err += GSL_DBL_EPSILON + fabs(t) * poly_1F1_err;
    }
  }

  stat_e = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err,
                                    poly_1F1_val, poly_1F1_err,
                                    result);

  return GSL_ERROR_SELECT_4(stat_e, stat_f, stat_g1, stat_g2);
}


/* Evaluate the polynomial based on the confluent hypergeometric
 * function in a safe way, with no restriction on the arguments.
 *
 * assumes x != 0
 */
static
int
 laguerre_n_poly_safe(const int n, const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  const MpIeee b=  a + 1.0;
  const MpIeee mx=  -x;
  const MpIeee tc_sgn=  (x < 0.0 ? 1.0 : (GSL_IS_ODD(n) ? -1.0 : 1.0));
  gsl_sf_result tc;
  int  stat_tc=  gsl_sf_taylorcoeff_e(n, fabs(x), &tc);

  if(stat_tc == GSL_SUCCESS) {
    MpIeee term=  tc.val * tc_sgn;
    MpIeee sum_val=  term;
    MpIeee sum_err=  tc.err;
    int  k;
    for(k=n-1; k>=0; k--) {
      term *= ((b+k)/(n-k))*(k+MpIeee( "1.0" ))/mx;
      sum_val += term;
      sum_err += MpIeee( "4.0" ) * GSL_DBL_EPSILON * fabs(term);
    }
    result->val = sum_val;
    result->err = sum_err + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(stat_tc == GSL_EOVRFLW) {
    result->val = 0.0; /* FIXME: should be Inf */
    result->err = 0.0;
    return stat_tc;
  }
  else {
    result->val = 0.0;
    result->err = 0.0;
    return stat_tc;
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_laguerre_1_e(const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = 1.0 + a - x;
    result->err = 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(a) + fabs(x));
    return GSL_SUCCESS;
  }
}

int
 gsl_sf_laguerre_2_e(const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(a == -2.0) {
    result->val = 0.5*x*x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee c0=  MpIeee( "0.5" ) * (MpIeee( "2.0" )+a)*(MpIeee( "1.0" )+a);
    MpIeee c1=  -(MpIeee( "2.0" )+a);
    MpIeee c2=  -MpIeee( "0.5" )/(MpIeee( "2.0" )+a);
    result->val  = c0 + c1*x*(1.0 + c2*x);
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(c0) + 2.0 * fabs(c1*x) * (1.0 + 2.0 * fabs(c2*x)));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
 gsl_sf_laguerre_3_e(const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(a == -2.0) {
    MpIeee x2_6=  x*x/MpIeee( "6.0" );
    result->val  = x2_6 * (3.0 - x);
    result->err  = x2_6 * (3.0 + fabs(x)) * 2.0 * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(a == -3.0) {
    result->val = -x*x/6.0;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee c0=  (MpIeee( "3.0" )+a)*(MpIeee( "2.0" )+a)*(MpIeee( "1.0" )+a) / MpIeee( "6.0" );
    MpIeee c1=  -c0 * MpIeee( "3.0" ) / (MpIeee( "1.0" )+a);
    MpIeee c2=  -MpIeee( "1.0" )/(MpIeee( "2.0" )+a);
    MpIeee c3=  -MpIeee( "1.0" )/(MpIeee( "3.0" )*(MpIeee( "3.0" )+a));
    result->val  = c0 + c1*x*(1.0 + c2*x*(1.0 + c3*x));
    result->err  = 1.0 + 2.0 * fabs(c3*x);
    result->err  = 1.0 + 2.0 * fabs(c2*x) * result->err;
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(c0) + 2.0 * fabs(c1*x) * result->err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_laguerre_n_e(const int n, const MpIeee a, const MpIeee x,
                           gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = 1.0 + a - x;
    result->err = 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(a) + fabs(x));
    return GSL_SUCCESS;
  }
  else if(x == 0.0) {
    MpIeee product=  a + MpIeee( "1.0" );
    int  k;
    for(k=2; k<=n; k++) {
      product *= (a + k)/k;
    }
    result->val = product;
    result->err = 2.0 * (n + 1.0) * GSL_DBL_EPSILON * fabs(product) + GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 0.0 && a > -1.0) {
    /* In this case all the terms in the polynomial
     * are of the same sign. Note that this also
     * catches overflows correctly.
     */
    return laguerre_n_cp(n, a, x, result);
  }
  else if(n < 5 || (x > 0.0 && a < -n-1)) {
    /* Either the polynomial will not lose too much accuracy
     * or all the terms are negative. In any case,
     * the error estimate here is good. We try both
     * explicit summation methods, as they have different
     * characteristics. One may underflow/overflow while the
     * other does not.
     */
    if(laguerre_n_cp(n, a, x, result) == GSL_SUCCESS)
      return GSL_SUCCESS;
    else
      return laguerre_n_poly_safe(n, a, x, result);
  }
  else if(n > 1.0e+07 && x > 0.0 && a > -1.0 && x < 2.0*(a+1.0)+4.0*n) {
    return laguerre_large_n(n, a, x, result);
  }
  else if(a > 0.0 || (x > 0.0 && a < -n-1)) {
    gsl_sf_result lg2;
    int  stat_lg2=  gsl_sf_laguerre_2_e(a, x, &lg2);
    MpIeee Lkm1=  MpIeee( "1.0" ) + a - x;
    MpIeee Lk=  lg2.val;
    MpIeee Lkp1;
    int  k;

    for(k=2; k<n; k++) {
      Lkp1 = (-(k+a)*Lkm1 + (MpIeee( "2.0" )*k+a+MpIeee( "1.0" )-x)*Lk)/(k+MpIeee( "1.0" ));
      Lkm1 = Lk;
      Lk   = Lkp1;
    }
    result->val = Lk;
    result->err = (fabs(lg2.err/lg2.val) + GSL_DBL_EPSILON) * n * fabs(Lk);
    return stat_lg2;
  }
  else {
    /* Despair... or magic? */
    return laguerre_n_poly_safe(n, a, x, result);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_laguerre_1(MpIeee a, MpIeee x)
{
  EVAL_RESULT(gsl_sf_laguerre_1_e(a, x, &result));
}

MpIeee gsl_sf_laguerre_2(MpIeee a, MpIeee x)
{
  EVAL_RESULT(gsl_sf_laguerre_2_e(a, x, &result));
}

MpIeee gsl_sf_laguerre_3(MpIeee a, MpIeee x)
{
  EVAL_RESULT(gsl_sf_laguerre_3_e(a, x, &result));
}

MpIeee gsl_sf_laguerre_n(int  n, MpIeee a, MpIeee x)
{
  EVAL_RESULT(gsl_sf_laguerre_n_e(n, a, x, &result));
}
