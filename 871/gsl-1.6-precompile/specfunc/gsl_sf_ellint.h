#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_ellint.h
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

/* Author: G. Jungman */

#ifndef __GSL_SF_ELLINT_H__
#define __GSL_SF_ELLINT_H__

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


/* Legendre form of complete elliptic integrals
 *
 * K(k) = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
 * E(k) = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
 *
 * exceptions: GSL_EDOM
 */
int  gsl_sf_ellint_Kcomp_e(MpIeee k, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_Kcomp(MpIeee k, gsl_mode_t mode);

int  gsl_sf_ellint_Ecomp_e(MpIeee k, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_Ecomp(MpIeee k, gsl_mode_t mode);


/* Legendre form of incomplete elliptic integrals
 *
 * F(phi,k)   = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * E(phi,k)   = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * P(phi,k,n) = Integral[(1 + n Sin[t]^2)^(-1)/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * D(phi,k,n) = R_D(1-Sin[phi]^2, 1-k^2 Sin[phi]^2, 1.0)
 *
 * F: [Carlson, Numerische Mathematik 33 (1979) 1, (4.1)]
 * E: [Carlson, ", (4.2)]
 * P: [Carlson, ", (4.3)]
 * D: [Carlson, ", (4.4)]
 *
 * exceptions: GSL_EDOM
 */
int  gsl_sf_ellint_F_e(MpIeee phi, MpIeee k, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_F(MpIeee phi, MpIeee k, gsl_mode_t mode);

int  gsl_sf_ellint_E_e(MpIeee phi, MpIeee k, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_E(MpIeee phi, MpIeee k, gsl_mode_t mode);

int  gsl_sf_ellint_P_e(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_P(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode);

int  gsl_sf_ellint_D_e(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_D(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode);


/* Carlson's symmetric basis of functions
 *
 * RC(x,y)   = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1)], {t,0,Inf}]
 * RD(x,y,z) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2), {t,0,Inf}]
 * RF(x,y,z) = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2), {t,0,Inf}]
 * RJ(x,y,z,p) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1), {t,0,Inf}]
 *
 * exceptions: GSL_EDOM
 */
int  gsl_sf_ellint_RC_e(MpIeee x, MpIeee y, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_RC(MpIeee x, MpIeee y, gsl_mode_t mode);

int  gsl_sf_ellint_RD_e(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_RD(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode);

int  gsl_sf_ellint_RF_e(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_RF(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode);

int  gsl_sf_ellint_RJ_e(MpIeee x, MpIeee y, MpIeee z, MpIeee p, gsl_mode_t mode, gsl_sf_result * result);
MpIeee gsl_sf_ellint_RJ(MpIeee x, MpIeee y, MpIeee z, MpIeee p, gsl_mode_t mode);


__END_DECLS

#endif /* __GSL_SF_ELLINT_H__ */
