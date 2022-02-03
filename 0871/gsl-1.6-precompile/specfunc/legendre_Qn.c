#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/legendre_Qn.c
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
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_legendre.h>

#include "error.h"

/* Evaluate f_{ell+1}/f_ell
 * f_ell := Q^{b}_{a+ell}(x)
 * x > 1
 */
static
int
 legendreQ_CF1_xgt1(int  ell, MpIeee a, MpIeee b, MpIeee x, MpIeee * result)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  ell + MpIeee( "1.0" ) + a + b;
  MpIeee b1=  (MpIeee( "2.0" )*(ell+MpIeee( "1.0" )+a) + MpIeee( "1.0" )) * x;
  MpIeee An=  b1*Anm1 + a1*Anm2;
  MpIeee Bn=  b1*Bnm1 + a1*Bnm2;
  MpIeee an;MpIeee  bn;
  MpIeee fn=  An/Bn;

  while(n < maxiter) {
    MpIeee old_fn;
    MpIeee del;
    MpIeee lna;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    lna = ell + n + a;
    an = b*b - lna*lna;
    bn = (MpIeee( "2.0" )*lna + MpIeee( "1.0" )) * x;
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;

    if(fabs(del - 1.0) < 4.0*GSL_DBL_EPSILON) break;
  }

  *result = fn;

  if(n == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS; 
}


/* Uniform asymptotic for Q_l(x).
 * Assumes x > -1.0 and x != 1.0.
 * Discards second order and higher terms.
 */
static
int
 legendre_Ql_asymp_unif(const MpIeee ell, const MpIeee x, gsl_sf_result * result)
{
  if(x < 1.0) {
    MpIeee u=  ell + MpIeee( "0.5" );
    MpIeee th=  acos(x);
    gsl_sf_result Y0, Y1;
    int  stat_Y0;int   stat_Y1;
    int  stat_m;
    MpIeee pre;
    MpIeee B00;
    MpIeee sum;

    /* B00 = 1/8 (1 - th cot(th) / th^2
     * pre = sqrt(th/sin(th))
     */
    if(th < GSL_ROOT4_DBL_EPSILON) {
      B00 = (MpIeee( "1.0" ) + th*th/MpIeee( "15.0" ))/MpIeee( "24.0" );
      pre = MpIeee( "1.0" ) + th*th/MpIeee( "12.0" );
    }
    else {
      MpIeee sin_th=  sqrt(MpIeee( "1.0" ) - x*x);
      MpIeee cot_th=  x / sin_th;
      B00 = MpIeee( "1.0" )/MpIeee( "8.0" ) * (MpIeee( "1.0" ) - th * cot_th) / (th*th);
      pre = sqrt(th/sin_th);
    }

    stat_Y0 = gsl_sf_bessel_Y0_e(u*th, &Y0);
    stat_Y1 = gsl_sf_bessel_Y1_e(u*th, &Y1);

    sum = -MpIeee( "0.5" )*M_PI * (Y0.val + th/u * Y1.val * B00);

    stat_m = gsl_sf_multiply_e(pre, sum, result);
    result->err += 0.5*M_PI * fabs(pre) * (Y0.err + fabs(th/u*B00)*Y1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_3(stat_m, stat_Y0, stat_Y1);
  }
  else {
    MpIeee u=  ell + MpIeee( "0.5" );
    MpIeee xi=  acosh(x);
    gsl_sf_result K0_scaled, K1_scaled;
    int  stat_K0;int   stat_K1;
    int  stat_e;
    MpIeee pre;
    MpIeee B00;
    MpIeee sum;

    /* B00 = -1/8 (1 - xi coth(xi) / xi^2
     * pre = sqrt(xi/sinh(xi))
     */
    if(xi < GSL_ROOT4_DBL_EPSILON) {
      B00 = (MpIeee( "1.0" )-xi*xi/MpIeee( "15.0" ))/MpIeee( "24.0" );
      pre = MpIeee( "1.0" ) - xi*xi/MpIeee( "12.0" );
    }
    else {
      MpIeee sinh_xi=  sqrt(x*x - MpIeee( "1.0" ));
      MpIeee coth_xi=  x / sinh_xi;
      B00 = -MpIeee( "1.0" )/MpIeee( "8.0" ) * (MpIeee( "1.0" ) - xi * coth_xi) / (xi*xi);
      pre = sqrt(xi/sinh_xi);
    }

    stat_K0 = gsl_sf_bessel_K0_scaled_e(u*xi, &K0_scaled);
    stat_K1 = gsl_sf_bessel_K1_scaled_e(u*xi, &K1_scaled);

    sum = K0_scaled.val - xi/u * K1_scaled.val * B00;

    stat_e = gsl_sf_exp_mult_e(-u*xi, pre * sum, result);
    result->err  = GSL_DBL_EPSILON * fabs(result->val) * fabs(u*xi);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_3(stat_e, stat_K0, stat_K1);
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_legendre_Q0_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0 || x == 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(x*x < GSL_ROOT6_DBL_EPSILON) { /* |x| <~ 0.05 */
    const MpIeee c3=  1.0/3.0;
    const MpIeee c5=  1.0/5.0;
    const MpIeee c7=  1.0/7.0;
    const MpIeee c9=  1.0/9.0;
    const MpIeee c11=  1.0/11.0;
    const MpIeee y=  x * x;
    const MpIeee series=  1.0 + y*(c3 + y*(c5 + y*(c7 + y*(c9 + y*c11))));
    result->val = x * series;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(x);
    return GSL_SUCCESS;
  }
  else if(x < 1.0) {
    result->val = 0.5 * log((1.0+x)/(1.0-x));
    result->err  = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 10.0) {
    result->val = 0.5 * log((x+1.0)/(x-1.0));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x*GSL_DBL_MIN < 2.0) {
    const MpIeee y=  1.0/(x*x);
    const MpIeee c1=  1.0/3.0;
    const MpIeee c2=  1.0/5.0;
    const MpIeee c3=  1.0/7.0;
    const MpIeee c4=  1.0/9.0;
    const MpIeee c5=  1.0/11.0;
    const MpIeee c6=  1.0/13.0;
    const MpIeee c7=  1.0/15.0;
    result->val = (1.0/x) * (1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7)))))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    UNDERFLOW_ERROR(result);
  }
}


int
 gsl_sf_legendre_Q1_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0 || x == 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(x*x < GSL_ROOT6_DBL_EPSILON) { /* |x| <~ 0.05 */
    const MpIeee c3=  1.0/3.0;
    const MpIeee c5=  1.0/5.0;
    const MpIeee c7=  1.0/7.0;
    const MpIeee c9=  1.0/9.0;
    const MpIeee c11=  1.0/11.0;
    const MpIeee y=  x * x;
    const MpIeee series=  1.0 + y*(c3 + y*(c5 + y*(c7 + y*(c9 + y*c11))));
    result->val = x * x * series - 1.0;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 1.0){
    result->val = 0.5 * x * (log((1.0+x)/(1.0-x))) - 1.0;
    result->err  = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 6.0) {
    result->val = 0.5 * x * log((x+1.0)/(x-1.0)) - 1.0;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x*GSL_SQRT_DBL_MIN < 0.99/M_SQRT3) {
    const MpIeee y=  1/(x*x);
    const MpIeee c1=  3.0/5.0;
    const MpIeee c2=  3.0/7.0;
    const MpIeee c3=  3.0/9.0;
    const MpIeee c4=  3.0/11.0;
    const MpIeee c5=  3.0/13.0;
    const MpIeee c6=  3.0/15.0;
    const MpIeee c7=  3.0/17.0;
    const MpIeee c8=  3.0/19.0;
    const MpIeee sum=  1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*c8)))))));
    result->val = sum / (3.0*x*x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    UNDERFLOW_ERROR(result);
  }
}


int
 gsl_sf_legendre_Ql_e(const int l, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0 || x == 1.0 || l < 0) {
    DOMAIN_ERROR(result);
  }
  else if(l == 0) {
    return gsl_sf_legendre_Q0_e(x, result);
  }
  else if(l == 1) {
    return gsl_sf_legendre_Q1_e(x, result);
  }
  else if(l > 100000) {
    return legendre_Ql_asymp_unif(l, x, result);
  }
  else if(x < 1.0){
    /* Forward recurrence.
     */
    gsl_sf_result Q0, Q1;
    int  stat_Q0=  gsl_sf_legendre_Q0_e(x, &Q0);
    int  stat_Q1=  gsl_sf_legendre_Q1_e(x, &Q1);
    MpIeee Qellm1=  Q0.val;
    MpIeee Qell=  Q1.val;
    MpIeee Qellp1;
    int  ell;
    for(ell=1; ell<l; ell++) {
      Qellp1 = (x*(MpIeee( "2.0" )*ell + MpIeee( "1.0" )) * Qell - ell * Qellm1) / (ell + MpIeee( "1.0" ));
      Qellm1 = Qell;
      Qell   = Qellp1;
    }
    result->val = Qell;
    result->err = GSL_DBL_EPSILON * l * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_Q0, stat_Q1);
  }
  else {
    /* x > 1.0 */

    MpIeee rat;
    int  stat_CF1=  legendreQ_CF1_xgt1(l, 0.0, 0.0, x, &rat);
    int  stat_Q;
    MpIeee Qellp1=  rat * GSL_SQRT_DBL_MIN;
    MpIeee Qell=  GSL_SQRT_DBL_MIN;
    MpIeee Qellm1;
    int  ell;
    for(ell=l; ell>0; ell--) {
      Qellm1 = (x * (MpIeee( "2.0" )*ell + MpIeee( "1.0" )) * Qell - (ell+MpIeee( "1.0" )) * Qellp1) / ell;
      Qellp1 = Qell;
      Qell   = Qellm1;
    }

    if(fabs(Qell) > fabs(Qellp1)) {
      gsl_sf_result Q0;
      stat_Q = gsl_sf_legendre_Q0_e(x, &Q0);
      result->val = GSL_SQRT_DBL_MIN * Q0.val / Qell;
      result->err = l * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
      gsl_sf_result Q1;
      stat_Q = gsl_sf_legendre_Q1_e(x, &Q1);
      result->val = GSL_SQRT_DBL_MIN * Q1.val / Qellp1;
      result->err = l * GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_ERROR_SELECT_2(stat_Q, stat_CF1);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_legendre_Q0(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_Q0_e(x, &result));
}

MpIeee gsl_sf_legendre_Q1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_Q1_e(x, &result));
}

MpIeee gsl_sf_legendre_Ql(const int l, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_Ql_e(l, x, &result));
}
