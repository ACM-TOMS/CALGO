#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_Knu.c
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
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "bessel.h"
#include "bessel_temme.h"

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_bessel_Knu_scaled_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0 || nu < 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    int  N=  MpIeee(nu + 0.5).toInt();
    MpIeee mu=  nu - N;      /* -1/2 <= mu <= 1/2 */
    MpIeee K_mu;MpIeee  K_mup1;MpIeee  Kp_mu;
    MpIeee K_nu;MpIeee  K_nup1;MpIeee  K_num1;
    int  n;

    if(x < 2.0) {
      gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    else {
      gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }

    /* recurse forward to obtain K_num1, K_nu */
    K_nu   = K_mu;
    K_nup1 = K_mup1;

    for(n=0; n<N; n++) {
      K_num1 = K_nu;
      K_nu   = K_nup1;
      K_nup1 = MpIeee( "2.0" )*(mu+n+MpIeee( "1" ))/x * K_nu + K_num1;
    }

    result->val = K_nu;
    result->err = 2.0 * GSL_DBL_EPSILON * (N + 4.0) * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_bessel_Knu_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result b;
  int  stat_K=  gsl_sf_bessel_Knu_scaled_e(nu, x, &b);
  int  stat_e=  gsl_sf_exp_mult_err_e(-x, 0.0, b.val, b.err, result);
  return GSL_ERROR_SELECT_2(stat_e, stat_K);
}


int
 gsl_sf_bessel_lnKnu_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0 || nu < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(nu == 0.0) {
    gsl_sf_result K_scaled;
    /* This cannot underflow, and
     * it will not throw GSL_EDOM
     * since that is already checked.
     */
    gsl_sf_bessel_K0_scaled_e(x, &K_scaled);
    result->val  = -x + log(fabs(K_scaled.val));
    result->err  = GSL_DBL_EPSILON * fabs(x) + fabs(K_scaled.err/K_scaled.val);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 2.0 && nu > 1.0) {
    /* Make use of the inequality
     * Knu(x) <= 1/2 (2/x)^nu Gamma(nu),
     * which follows from the integral representation
     * [Abramowitz+Stegun, 9.6.23 (2)]. With this
     * we decide whether or not there is an overflow
     * problem because x is small.
     */
    MpIeee ln_bound;
    gsl_sf_result lg_nu;
    gsl_sf_lngamma_e(nu, &lg_nu);
    ln_bound = -M_LN2 - nu*log(MpIeee( "0.5" )*x) + lg_nu.val;
    if(ln_bound > GSL_LOG_DBL_MAX - MpIeee( "20.0" )) {
      /* x must be very small or nu very large (or both).
       */
      MpIeee xi=  MpIeee( "0.25" )*x*x;
      MpIeee sum=  MpIeee( "1.0" ) - xi/(nu-MpIeee( "1.0" ));
      if(nu > 2.0) sum +=  (xi/(nu-MpIeee( "1.0" ))) * (xi/(nu-MpIeee( "2.0" )));
      result->val  = ln_bound + log(sum);
      result->err  = lg_nu.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    /* can drop-through here */
  }


  {
    /* We passed the above tests, so no problem.
     * Evaluate as usual. Note the possible drop-through
     * in the above code!
     */
    gsl_sf_result K_scaled;
    gsl_sf_bessel_Knu_scaled_e(nu, x, &K_scaled);
    result->val  = -x + log(fabs(K_scaled.val));
    result->err  = GSL_DBL_EPSILON * fabs(x) + fabs(K_scaled.err/K_scaled.val);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_Knu_scaled(const MpIeee nu, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_Knu_scaled_e(nu, x, &result));
}

MpIeee gsl_sf_bessel_Knu(const MpIeee nu, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_Knu_e(nu, x, &result));
}

MpIeee gsl_sf_bessel_lnKnu(const MpIeee nu, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_lnKnu_e(nu, x, &result));
}
