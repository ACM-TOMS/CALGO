#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/coulomb_bound.c
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
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_coulomb.h>

#include "error.h"
#include "check.h"

/* normalization for hydrogenic wave functions */
static
int
 R_norm(const int n, const int l, const MpIeee Z, gsl_sf_result * result)
{
  MpIeee A=  MpIeee( "2.0" )*Z/n;
  MpIeee pre=  sqrt(A*A*A /(MpIeee( "2.0" )*n));
  gsl_sf_result ln_a, ln_b;
  gsl_sf_result ex;
  int  stat_a=  gsl_sf_lnfact_e(n+l, &ln_a);
  int  stat_b=  gsl_sf_lnfact_e(n-l-1, &ln_b);
  MpIeee diff_val=  MpIeee( "0.5" )*(ln_b.val - ln_a.val);
  MpIeee diff_err=  MpIeee( "0.5" )*(ln_b.err + ln_a.err) + GSL_DBL_EPSILON * fabs(diff_val);
  int  stat_e=  gsl_sf_exp_err_e(diff_val, diff_err, &ex);
  result->val  = pre * ex.val;
  result->err  = pre * ex.err;
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_ERROR_SELECT_3(stat_e, stat_a, stat_b);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_hydrogenicR_1_e(const MpIeee Z, const MpIeee r, gsl_sf_result * result)
{
  if(Z > 0.0 && r >= 0.0) {
    MpIeee A=  MpIeee( "2.0" )*Z;
    MpIeee norm=  A*sqrt(Z);
    MpIeee ea=  exp(-Z*r);
    result->val = norm*ea;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) * fabs(Z*r);
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR(result);
  }
}


int
 gsl_sf_hydrogenicR_e(const int n, const int l,
                        const MpIeee Z, const MpIeee r,
                        gsl_sf_result * result)
{
  if(n < 1 || l > n-1 || Z <= 0.0 || r < 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    MpIeee A=  MpIeee( "2.0" )*Z/n;
    gsl_sf_result norm;
    int  stat_norm=  R_norm(n, l, Z, &norm);
    MpIeee rho=  A*r;
    MpIeee ea=  exp(-MpIeee( "0.5" )*rho);
    MpIeee pp=  gsl_sf_pow_int(rho, l);
    gsl_sf_result lag;
    int  stat_lag=  gsl_sf_laguerre_n_e(n-l-1, 2*l+1, rho, &lag);
    MpIeee W_val=  norm.val * ea * pp;
    MpIeee W_err=  norm.err * ea * pp;
    W_err += norm.val * ((MpIeee( "0.5" )*rho + MpIeee( "1.0" )) * GSL_DBL_EPSILON) * ea * pp;
    W_err += norm.val * ea * ((l+MpIeee( "1.0" )) * GSL_DBL_EPSILON) * pp;
    result->val  = W_val * lag.val;
    result->err  = W_val * lag.err + W_err * fabs(lag.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    if (stat_lag == GSL_SUCCESS && stat_norm == GSL_SUCCESS) {
      CHECK_UNDERFLOW(result);
    };
    return GSL_ERROR_SELECT_2(stat_lag, stat_norm);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_hydrogenicR_1(const MpIeee Z, const MpIeee r)
{
  EVAL_RESULT(gsl_sf_hydrogenicR_1_e(Z, r, &result));
}


MpIeee gsl_sf_hydrogenicR(const int n, const int l, const MpIeee Z, const MpIeee r)
{
  EVAL_RESULT(gsl_sf_hydrogenicR_e(n, l, Z, r, &result));
}
