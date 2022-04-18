#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_trig.h
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

#ifndef __GSL_SF_TRIG_H__
#define __GSL_SF_TRIG_H__

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


/* Sin(x) with GSL semantics. This is actually important
 * because we want to control the error estimate, and trying
 * to guess the error for the standard library implementation
 * every time it is used would be a little goofy.
 */
int  gsl_sf_sin_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_sin(const MpIeee x);


/* Cos(x) with GSL semantics.
 */
int  gsl_sf_cos_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_cos(const MpIeee x);


/* Hypot(x,y) with GSL semantics.
 */
int  gsl_sf_hypot_e(const MpIeee x, const MpIeee y, gsl_sf_result * result);
MpIeee gsl_sf_hypot(const MpIeee x, const MpIeee y);


/* Sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int  gsl_sf_complex_sin_e(const MpIeee zr, const MpIeee zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* Cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int  gsl_sf_complex_cos_e(const MpIeee zr, const MpIeee zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* Log(Sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int  gsl_sf_complex_logsin_e(const MpIeee zr, const MpIeee zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* Sinc(x) = sin(pi x) / (pi x)
 *
 * exceptions: none
 */
int  gsl_sf_sinc_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_sinc(const MpIeee x);


/* Log(Sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int  gsl_sf_lnsinh_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_lnsinh(const MpIeee x);


/* Log(Cosh(x))
 *
 * exceptions: none
 */
int  gsl_sf_lncosh_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_lncosh(const MpIeee x);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int  gsl_sf_polar_to_rect(const MpIeee r, const MpIeee theta, gsl_sf_result * x, gsl_sf_result * y);

/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int  gsl_sf_rect_to_polar(const MpIeee x, const MpIeee y, gsl_sf_result * r, gsl_sf_result * theta);

/* Sin(x) for quantity with an associated error.
 */
int  gsl_sf_sin_err_e(const MpIeee x, const MpIeee dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int  gsl_sf_cos_err_e(const MpIeee x, const MpIeee dx, gsl_sf_result * result);


/* Force an angle to lie in the range (-pi,pi].
 *
 * exceptions: GSL_ELOSS
 */
int  gsl_sf_angle_restrict_symm_e(MpIeee * theta);
MpIeee gsl_sf_angle_restrict_symm(const MpIeee theta);


/* Force an angle to lie in the range [0, 2pi)
 *
 * exceptions: GSL_ELOSS
 */
int  gsl_sf_angle_restrict_pos_e(MpIeee * theta);
MpIeee gsl_sf_angle_restrict_pos(const MpIeee theta);


int  gsl_sf_angle_restrict_symm_err_e(const MpIeee theta, gsl_sf_result * result);

int  gsl_sf_angle_restrict_pos_err_e(const MpIeee theta, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_TRIG_H__ */
