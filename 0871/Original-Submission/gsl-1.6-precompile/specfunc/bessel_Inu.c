#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_Inu.c
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
 gsl_sf_bessel_Inu_scaled_e(MpIeee nu, MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < MpIeee( "0.0" ) || nu < MpIeee( "0.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(x*x < MpIeee( "10.0" )*(nu+MpIeee( "1.0" ))) {
    gsl_sf_result b;
    MpIeee ex=  exp(-x);
    int  stat=  gsl_sf_bessel_IJ_taylor_e(nu, x, 1, 100, GSL_DBL_EPSILON, &b);
    result->val  = b.val * ex;
    result->err  = b.err * ex;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(0.5/(nu*nu + x*x) < GSL_ROOT3_DBL_EPSILON) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_e(nu, x, result);
  }
  else {
    int  N=  MpIeee(nu + 0.5).toInt();
    MpIeee mu=  nu - N;      /* -1/2 <= mu <= 1/2 */ 
    MpIeee K_mu;MpIeee  K_mup1;MpIeee  Kp_mu;
    MpIeee K_nu;MpIeee  K_nup1;MpIeee  K_num1;
    MpIeee I_nu_ratio;
    int  stat_Irat;
    int  stat_Kmu;
    int  n;

    /* obtain K_mu, K_mup1 */
    if(x < MpIeee( "2.0" )) {
      stat_Kmu = gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    else {
      stat_Kmu = gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }

    /* recurse forward to obtain K_num1, K_nu */
    K_nu   = K_mu;
    K_nup1 = K_mup1;

    for(n=0; n<N; n++) {
      K_num1 = K_nu;
      K_nu   = K_nup1;
      K_nup1 = MpIeee( "2.0" )*(mu+n+MpIeee( "1" ))/x * K_nu + K_num1;
    }

    /* calculate I_{nu+1}/I_nu */
    stat_Irat = gsl_sf_bessel_I_CF1_ser(nu, x, &I_nu_ratio);

    /* solve for I_nu */
    result->val = 1.0/(x * (K_nup1 + I_nu_ratio * K_nu));
    result->err = GSL_DBL_EPSILON * (0.5*N + 2.0) * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_Kmu, stat_Irat);
  }
}


int
 gsl_sf_bessel_Inu_e(MpIeee nu, MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result b;
  int  stat_I=  gsl_sf_bessel_Inu_scaled_e(nu, x, &b);
  int  stat_e=  gsl_sf_exp_mult_err_e(x, fabs(x*GSL_DBL_EPSILON),
                                        b.val, b.err,
                                        result);
  return GSL_ERROR_SELECT_2(stat_e, stat_I);
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_Inu_scaled(MpIeee nu, MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_Inu_scaled_e(nu, x, &result));
}


MpIeee gsl_sf_bessel_Inu(MpIeee nu, MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_Inu_e(nu, x, &result));
}
