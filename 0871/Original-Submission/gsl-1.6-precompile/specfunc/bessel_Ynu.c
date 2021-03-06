#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_Ynu.c
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

#include "error.h"

#include "bessel.h"
#include "bessel_olver.h"
#include "bessel_temme.h"

/* Perform forward recurrence for Y_nu(x) and Y'_nu(x)
 *
 *        Y_{nu+1} =  nu/x Y_nu - Y'_nu
 *       Y'_{nu+1} = -(nu+1)/x Y_{nu+1} + Y_nu
 */
#if 0
static
int
 bessel_Y_recur(const MpIeee nu_min, const MpIeee x, const int kmax,
               const MpIeee Y_start, const MpIeee Yp_start,
               MpIeee * Y_end, MpIeee * Yp_end)
{
  MpIeee x_inv=  MpIeee( "1.0" )/x;
  MpIeee nu=  nu_min;
  MpIeee Y_nu=  Y_start;
  MpIeee Yp_nu=  Yp_start;
  int  k;

  for(k=1; k<=kmax; k++) {
    MpIeee nuox=  nu*x_inv;
    MpIeee Y_nu_save=  Y_nu;
    Y_nu  = -Yp_nu + nuox * Y_nu;
    Yp_nu = Y_nu_save - (nuox+x_inv) * Y_nu;
    nu += MpIeee( "1.0" );
  }
  *Y_end  = Y_nu;
  *Yp_end = Yp_nu;
  return GSL_SUCCESS;
}
#endif


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_bessel_Ynu_e(MpIeee nu, MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= MpIeee( "0.0" ) || nu < MpIeee( "0.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(nu > MpIeee( "50.0" )) {
    return gsl_sf_bessel_Ynu_asymp_Olver_e(nu, x, result);
  }
  else {
    /* -1/2 <= mu <= 1/2 */
    int  N=  MpIeee(nu + 0.5).toInt();
    MpIeee mu=  nu - N;

    gsl_sf_result Y_mu, Y_mup1;
    int  stat_mu;
    MpIeee Ynm1;
    MpIeee Yn;
    MpIeee Ynp1;
    int  n;

    if(x < MpIeee( "2.0" )) {
      /* Determine Ymu, Ymup1 directly. This is really
       * an optimization since this case could as well
       * be handled by a call to gsl_sf_bessel_JY_mu_restricted(),
       * as below.
       */
      stat_mu = gsl_sf_bessel_Y_temme(mu, x, &Y_mu, &Y_mup1);
    }
    else {
      /* Determine Ymu, Ymup1 and Jmu, Jmup1.
       */
      gsl_sf_result J_mu, J_mup1;
      stat_mu = gsl_sf_bessel_JY_mu_restricted(mu, x, &J_mu, &J_mup1, &Y_mu, &Y_mup1);
    }

    /* Forward recursion to get Ynu, Ynup1.
     */
    Ynm1 = Y_mu.val;
    Yn   = Y_mup1.val;
    for(n=1; n<=N; n++) {
      Ynp1 = MpIeee( "2.0" )*(mu+n)/x * Yn - Ynm1;
      Ynm1 = Yn;
      Yn   = Ynp1;
    }

    result->val  = Ynm1; /* Y_nu */
    result->err  = (N + 1.0) * fabs(Ynm1) * (fabs(Y_mu.err/Y_mu.val) + fabs(Y_mup1.err/Y_mup1.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(Ynm1);

    return stat_mu;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_Ynu(const MpIeee nu, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_Ynu_e(nu, x, &result));
}
