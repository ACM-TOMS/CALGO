#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/beta_inc.c
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
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>

#include "error.h"
#include "check.h"

static
int
 beta_cont_frac(
  const MpIeee a,
  const MpIeee b,
  const MpIeee x,
  gsl_sf_result * result
  )
{
  const unsigned int  max_iter=  512;        /* control iterations      */
  const MpIeee cutoff=  2.0 * GSL_DBL_MIN;  /* control the zero cutoff */
  unsigned int  iter_count=  0;
  MpIeee cf;

  /* standard initialization for continued fraction */
  MpIeee num_term=  MpIeee( "1.0" );
  MpIeee den_term=  MpIeee( "1.0" ) - (a+b)*x/(a+MpIeee( "1.0" ));
  if (fabs(den_term) < cutoff) den_term = cutoff;
  den_term = MpIeee( "1.0" )/den_term;
  cf = den_term;

  while(iter_count < max_iter) {
    const int k  = iter_count + 1;
    MpIeee coeff=  k*(b-k)*x/(((a-MpIeee( "1.0" ))+MpIeee( "2" )*k)*(a+MpIeee( "2" )*k));
    MpIeee delta_frac;

    /* first step */
    den_term = MpIeee( "1.0" ) + coeff*den_term;
    num_term = MpIeee( "1.0" ) + coeff/num_term;
    if(fabs(den_term) < cutoff) den_term = cutoff;
    if(fabs(num_term) < cutoff) num_term = cutoff;
    den_term  = MpIeee( "1.0" )/den_term;

    delta_frac = den_term * num_term;
    cf *= delta_frac;

    coeff = -(a+k)*(a+b+k)*x/((a+MpIeee( "2" )*k)*(a+MpIeee( "2" )*k+MpIeee( "1.0" )));

    /* second step */
    den_term = MpIeee( "1.0" ) + coeff*den_term;
    num_term = MpIeee( "1.0" ) + coeff/num_term;
    if(fabs(den_term) < cutoff) den_term = cutoff;
    if(fabs(num_term) < cutoff) num_term = cutoff;
    den_term = MpIeee( "1.0" )/den_term;

    delta_frac = den_term*num_term;
    cf *= delta_frac;

    if(fabs(delta_frac-1.0) < 2.0*GSL_DBL_EPSILON) break;

    ++iter_count;
  }

  result->val = cf;
  result->err = iter_count * 4.0 * GSL_DBL_EPSILON * fabs(cf);

  if(iter_count >= max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_beta_inc_e(
  const MpIeee a,
  const MpIeee b,
  const MpIeee x,
  gsl_sf_result * result
  )
{
  if(a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result ln_beta;
    gsl_sf_result ln_x;
    gsl_sf_result ln_1mx;
    gsl_sf_result prefactor;
    const int stat_ln_beta = gsl_sf_lnbeta_e(a, b, &ln_beta);
    const int stat_ln_1mx = gsl_sf_log_1plusx_e(-x, &ln_1mx);
    const int stat_ln_x = gsl_sf_log_e(x, &ln_x);
    const int stat_ln = GSL_ERROR_SELECT_3(stat_ln_beta, stat_ln_1mx, stat_ln_x);

    const MpIeee ln_pre_val=  -ln_beta.val + a * ln_x.val + b * ln_1mx.val;
    const MpIeee ln_pre_err=   ln_beta.err + fabs(a*ln_x.err) + fabs(b*ln_1mx.err);
    const int stat_exp = gsl_sf_exp_err_e(ln_pre_val, ln_pre_err, &prefactor);

    if(stat_ln != GSL_SUCCESS) {
      result->val = 0.0;
      result->err = 0.0;
      GSL_ERROR ("error", GSL_ESANITY);
    }

    if(x < (a + 1.0)/(a+b+2.0)) {
      /* Apply continued fraction directly. */
      gsl_sf_result cf;
      const int stat_cf = beta_cont_frac(a, b, x, &cf);
      int  stat;
      result->val = prefactor.val * cf.val / a;
      result->err = (fabs(prefactor.err * cf.val) + fabs(prefactor.val * cf.err))/a;

      stat = GSL_ERROR_SELECT_2(stat_exp, stat_cf);
      if(stat == GSL_SUCCESS) {
        CHECK_UNDERFLOW(result);
      }
      return stat;
    }
    else {
      /* Apply continued fraction after hypergeometric transformation. */
      gsl_sf_result cf;
      const int stat_cf = beta_cont_frac(b, a, 1.0-x, &cf);
      int  stat;
      const MpIeee term=  prefactor.val * cf.val / b;
      result->val  = 1.0 - term;
      result->err  = fabs(prefactor.err * cf.val)/b;
      result->err += fabs(prefactor.val * cf.err)/b;
      result->err += 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(term));
      stat = GSL_ERROR_SELECT_2(stat_exp, stat_cf);
      if(stat == GSL_SUCCESS) {
        CHECK_UNDERFLOW(result);
      }
      return stat;
    }
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_beta_inc(const MpIeee a, const MpIeee b, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_beta_inc_e(a, b, x, &result));
}
