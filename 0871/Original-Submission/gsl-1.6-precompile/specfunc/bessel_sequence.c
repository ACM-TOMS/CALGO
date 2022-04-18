#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_sequence.c
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


#define DYDX_p(p,u,x) (-(p)/(x) + (((nu)*(nu))/((x)*(x))-1.0)*(u))
#define DYDX_u(p,u,x) (p)

static int
 rk_step(MpIeee nu, MpIeee x, MpIeee dx, MpIeee * Jp, MpIeee * J)
{
  MpIeee p_0=  *Jp;
  MpIeee u_0=  *J;

  MpIeee p_1=  dx * DYDX_p(p_0, u_0, x);
  MpIeee u_1=  dx * DYDX_u(p_0, u_0, x);

  MpIeee p_2=  dx * DYDX_p(p_0 + MpIeee( "0.5" )*p_1, u_0 + MpIeee( "0.5" )*u_1, x + MpIeee( "0.5" )*dx);
  MpIeee u_2=  dx * DYDX_u(p_0 + MpIeee( "0.5" )*p_1, u_0 + MpIeee( "0.5" )*u_1, x + MpIeee( "0.5" )*dx);

  MpIeee p_3=  dx * DYDX_p(p_0 + MpIeee( "0.5" )*p_2, u_0 + MpIeee( "0.5" )*u_2, x + MpIeee( "0.5" )*dx);
  MpIeee u_3=  dx * DYDX_u(p_0 + MpIeee( "0.5" )*p_2, u_0 + MpIeee( "0.5" )*u_2, x + MpIeee( "0.5" )*dx);

  MpIeee p_4=  dx * DYDX_p(p_0 + p_3, u_0 + u_3, x + dx);
  MpIeee u_4=  dx * DYDX_u(p_0 + p_3, u_0 + u_3, x + dx);

  *Jp = p_0 + p_1/MpIeee( "6.0" ) + p_2/MpIeee( "3.0" ) + p_3/MpIeee( "3.0" ) + p_4/MpIeee( "6.0" );
  *J  = u_0 + u_1/MpIeee( "6.0" ) + u_2/MpIeee( "3.0" ) + u_3/MpIeee( "3.0" ) + u_4/MpIeee( "6.0" );

  return GSL_SUCCESS;
}


int
 gsl_sf_bessel_sequence_Jnu_e(MpIeee nu, gsl_mode_t mode, size_t size, MpIeee * v)
{
  /* CHECK_POINTER(v) */

  if(nu < MpIeee( "0.0" )) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(size == 0) {
    GSL_ERROR ("error", GSL_EINVAL);
  }
  else {
    const gsl_prec_t goal   = GSL_MODE_PREC(mode);
    const MpIeee dx_array[] =  { 0.001, 0.03, 0.1 }; /* double, single, approx */
    const MpIeee dx_nominal=  dx_array[goal];

    const int cnu = (int) MpIeee(ceil(nu)).toInt();
    const MpIeee nu13=  pow(nu,1.0/3.0);
    const MpIeee smalls[] =  { 0.01, 0.02, 0.4, 0.7, 1.3, 2.0, 2.5, 3.2, 3.5, 4.5, 6.0 };
    const MpIeee x_small=  ( nu >= 10.0 ? nu - nu13 : smalls[cnu] );

    gsl_sf_result J0, J1;
    MpIeee Jp;MpIeee  J;
    MpIeee x;
    size_t i = 0;

    /* Calculate the first point. */
    x = v[0];
    gsl_sf_bessel_Jnu_e(nu, x, &J0);
    v[0] = J0.val;
    ++i;

    /* Step over the idiot case where the
     * first point was actually zero.
     */
    if(x == MpIeee( "0.0" )) {
      if(v[1] <= x) {
        /* Strict ordering failure. */
        GSL_ERROR ("error", GSL_EFAILED);
      }
      x = v[1];
      gsl_sf_bessel_Jnu_e(nu, x, &J0);
      v[1] = J0.val;
      ++i;
    }

    /* Calculate directly as long as the argument
     * is small. This is necessary because the
     * integration is not very good there.
     */
    while(v[i] < x_small && i < size) {
      if(v[i] <= x) {
        /* Strict ordering failure. */
        GSL_ERROR ("error", GSL_EFAILED);
      }
      x = v[i];
      gsl_sf_bessel_Jnu_e(nu, x, &J0);
      v[i] = J0.val;
      ++i;
    }

    /* At this point we are ready to integrate.
     * The value of x is the last calculated
     * point, which has the value J0; v[i] is
     * the next point we need to calculate. We
     * calculate nu+1 at x as well to get
     * the derivative, then we go forward.
     */
    gsl_sf_bessel_Jnu_e(nu+1.0, x, &J1);
    J  = J0.val;
    Jp = -J1.val + nu/x * J0.val;

    while(i < size) {
      const MpIeee dv=  v[i] - x;
      const int Nd    = (int) MpIeee(ceil(dv/dx_nominal)).toInt();
      const MpIeee dx=  dv / Nd;
      MpIeee xj;
      int  j;

      if(v[i] <= x) {
        /* Strict ordering failure. */
        GSL_ERROR ("error", GSL_EFAILED);
      }

      /* Integrate over interval up to next sample point.
       */
      for(j=0, xj=x; j<Nd; j++, xj += dx) {
        rk_step(nu, xj, dx, &Jp, &J);
      }

      /* Go to next interval. */
      x = v[i];
      v[i] = J;
      ++i;
    }

    return GSL_SUCCESS;
  }
}
