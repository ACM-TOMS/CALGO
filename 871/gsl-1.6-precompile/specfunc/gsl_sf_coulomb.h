#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_coulomb.h
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

#ifndef __GSL_SF_COULOMB_H__
#define __GSL_SF_COULOMB_H__

#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Normalized hydrogenic bound states, radial dependence. */

/* R_1 := 2Z sqrt(Z) exp(-Z r)
 */
int  gsl_sf_hydrogenicR_1_e(const MpIeee Z, const MpIeee r, gsl_sf_result * result);
MpIeee gsl_sf_hydrogenicR_1(const MpIeee Z, const MpIeee r);

/* R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
 *
 * normalization such that psi(n,l,r) = R_n Y_{lm}
 */
int  gsl_sf_hydrogenicR_e(const int n, const int l, const MpIeee Z, const MpIeee r, gsl_sf_result * result);
MpIeee gsl_sf_hydrogenicR(const int n, const int l, const MpIeee Z, const MpIeee r);


/* Coulomb wave functions F_{lam_F}(eta,x), G_{lam_G}(eta,x)
 * and their derivatives; lam_G := lam_F - k_lam_G
 *
 * lam_F, lam_G > -0.5
 * x > 0.0
 *
 * Conventions of Abramowitz+Stegun.
 *
 * Because there can be a large dynamic range of values,
 * overflows are handled gracefully. If an overflow occurs,
 * GSL_EOVRFLW is signalled and exponent(s) are returned
 * through exp_F, exp_G. These are such that
 *
 *   F_L(eta,x)  =  fc[k_L] * exp(exp_F)
 *   G_L(eta,x)  =  gc[k_L] * exp(exp_G)
 *   F_L'(eta,x) = fcp[k_L] * exp(exp_F)
 *   G_L'(eta,x) = gcp[k_L] * exp(exp_G)
 */
int
 gsl_sf_coulomb_wave_FG_e(const MpIeee eta, const MpIeee x,
                            const MpIeee lam_F,
                            const int  k_lam_G,
                            gsl_sf_result * F, gsl_sf_result * Fp,
                            gsl_sf_result * G, gsl_sf_result * Gp,
                            MpIeee * exp_F, MpIeee * exp_G);


/* F_L(eta,x) as array */
int  gsl_sf_coulomb_wave_F_array(
  MpIeee lam_min, int  kmax,
  MpIeee eta, MpIeee x,
  MpIeee * fc_array,
  MpIeee * F_exponent);

/* F_L(eta,x), G_L(eta,x) as arrays */
int  gsl_sf_coulomb_wave_FG_array(MpIeee lam_min, int  kmax,
                                MpIeee eta, MpIeee x,
                                MpIeee * fc_array, MpIeee * gc_array,
                                MpIeee * F_exponent,
                                MpIeee * G_exponent);

/* F_L(eta,x), G_L(eta,x), F'_L(eta,x), G'_L(eta,x) as arrays */
int  gsl_sf_coulomb_wave_FGp_array(MpIeee lam_min, int  kmax,
                                MpIeee eta, MpIeee x,
                                MpIeee * fc_array, MpIeee * fcp_array,
                                MpIeee * gc_array, MpIeee * gcp_array,
                                MpIeee * F_exponent,
                                MpIeee * G_exponent);

/* Coulomb wave function divided by the argument,
 * F(eta, x)/x. This is the function which reduces to
 * spherical Bessel functions in the limit eta->0.
 */
int  gsl_sf_coulomb_wave_sphF_array(MpIeee lam_min, int  kmax,
                                        MpIeee eta, MpIeee x,
                                        MpIeee * fc_array,
                                        MpIeee * F_exponent);


/* Coulomb wave function normalization constant.
 * [Abramowitz+Stegun 14.1.8, 14.1.9]
 */
int  gsl_sf_coulomb_CL_e(MpIeee L, MpIeee eta, gsl_sf_result * result);
int  gsl_sf_coulomb_CL_array(MpIeee Lmin, int  kmax, MpIeee eta, MpIeee * cl);


__END_DECLS

#endif /* __GSL_SF_COULOMB_H__ */
