#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_expint.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
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

#ifndef __GSL_SF_EXPINT_H__
#define __GSL_SF_EXPINT_H__

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


/* E_1(x) := Re[ Integrate[ Exp[-xt]/t, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_expint_E1_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_E1(const MpIeee x);


/* E_2(x) := Re[ Integrate[ Exp[-xt]/t^2, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_expint_E2_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_E2(const MpIeee x);


/* E_1_scaled(x) := exp(x) E_1(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_expint_E1_scaled_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_E1_scaled(const MpIeee x);


/* E_2_scaled(x) := exp(x) E_2(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_expint_E2_scaled_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_E2_scaled(const MpIeee x);


/* Ei(x) := - PV Integrate[ Exp[-t]/t, {t,-x,Infinity}]
 *       :=   PV Integrate[ Exp[t]/t, {t,-Infinity,x}]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_expint_Ei_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_Ei(const MpIeee x);


/* Ei_scaled(x) := exp(-x) Ei(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_expint_Ei_scaled_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_Ei_scaled(const MpIeee x);


/* Shi(x) := Integrate[ Sinh[t]/t, {t,0,x}]
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_Shi_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_Shi(const MpIeee x);


/* Chi(x) := Re[ M_EULER + log(x) + Integrate[(Cosh[t]-1)/t, {t,0,x}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int      gsl_sf_Chi_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_Chi(const MpIeee x);


/* Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
 *
 * x >= 0.0
 * exceptions: GSL_EDOM
 */
int      gsl_sf_expint_3_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_expint_3(MpIeee x);


/* Si(x) := Integrate[ Sin[t]/t, {t,0,x}]
 *
 * exceptions: none
 */
int      gsl_sf_Si_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_Si(const MpIeee x);


/* Ci(x) := -Integrate[ Cos[t]/t, {t,x,Infinity}]
 *
 * x > 0.0
 * exceptions: GSL_EDOM 
 */
int      gsl_sf_Ci_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_Ci(const MpIeee x);


/* AtanInt(x) := Integral[ Arctan[t]/t, {t,0,x}]
 *
 *
 * exceptions:
 */
int      gsl_sf_atanint_e(const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_atanint(const MpIeee x);


__END_DECLS

#endif /* __GSL_SF_EXPINT_H__ */
