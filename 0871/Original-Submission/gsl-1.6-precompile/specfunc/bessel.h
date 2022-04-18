#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel.h
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

#ifndef _BESSEL_H_
#define _BESSEL_H_

#include <gsl/gsl_sf_result.h>


/* Taylor expansion for J_nu(x) or I_nu(x)
 *   sign = -1  ==> Jnu
 *   sign = +1  ==> Inu
 */
int  gsl_sf_bessel_IJ_taylor_e(const MpIeee nu, const MpIeee x,
                                 const int sign,
                                 const int kmax,
                                 const MpIeee threshold,
                                 gsl_sf_result * result
                                 );

int  gsl_sf_bessel_Jnu_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result);
int  gsl_sf_bessel_Ynu_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result);

int  gsl_sf_bessel_Inu_scaled_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result);
int  gsl_sf_bessel_Knu_scaled_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result);

int  gsl_sf_bessel_Inu_scaled_asymp_unif_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result);
int  gsl_sf_bessel_Knu_scaled_asymp_unif_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result);


/* ratio = J_{nu+1}(x) / J_nu(x)
 * sgn   = sgn(J_nu(x))
 */
int
 gsl_sf_bessel_J_CF1(const MpIeee nu, const MpIeee x, MpIeee * ratio, MpIeee * sgn);


/* ratio = I_{nu+1}(x) / I_nu(x)
 */
int
 gsl_sf_bessel_I_CF1_ser(const MpIeee nu, const MpIeee x, MpIeee * ratio);


/* Evaluate the Steed method continued fraction CF2 for
 *
 * (J' + i Y')/(J + i Y) := P + i Q
 */
int
 gsl_sf_bessel_JY_steed_CF2(const MpIeee nu, const MpIeee x,
                           MpIeee * P, MpIeee * Q);


int
 gsl_sf_bessel_JY_mu_restricted(const MpIeee mu, const MpIeee x,
                               gsl_sf_result * Jmu, gsl_sf_result * Jmup1,
                               gsl_sf_result * Ymu, gsl_sf_result * Ymup1);


int
 gsl_sf_bessel_K_scaled_steed_temme_CF2(const MpIeee nu, const MpIeee x,
                                       MpIeee * K_nu, MpIeee * K_nup1,
                                       MpIeee * Kp_nu);


/* These are of use in calculating the oscillating
 * Bessel functions.
 *   cos(y - pi/4 + eps)
 *   sin(y - pi/4 + eps)
 */
int  gsl_sf_bessel_cos_pi4_e(MpIeee y, MpIeee eps, gsl_sf_result * result);
int  gsl_sf_bessel_sin_pi4_e(MpIeee y, MpIeee eps, gsl_sf_result * result);


#endif /* !_BESSEL_H_ */
