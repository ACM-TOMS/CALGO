#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_erf.h
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

#ifndef __GSL_SF_ERF_H__
#define __GSL_SF_ERF_H__

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


/* Complementary Error Function
 * erfc(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,x,Infinity}]
 *
 * exceptions: none
 */
int  gsl_sf_erfc_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_erfc(MpIeee x);


/* Log Complementary Error Function
 *
 * exceptions: none
 */
int  gsl_sf_log_erfc_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_log_erfc(MpIeee x);


/* Error Function
 * erf(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,0,x}]
 *
 * exceptions: none
 */
int  gsl_sf_erf_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_erf(MpIeee x);


/* Probability functions:
 * Z(x) :  Abramowitz+Stegun 26.2.1
 * Q(x) :  Abramowitz+Stegun 26.2.3
 *
 * exceptions: none
 */
int  gsl_sf_erf_Z_e(MpIeee x, gsl_sf_result * result);
int  gsl_sf_erf_Q_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_erf_Z(MpIeee x);
MpIeee gsl_sf_erf_Q(MpIeee x);


/* Hazard function, also known as the inverse Mill's ratio.
 *
 *   H(x) := Z(x)/Q(x)
 *         = Sqrt[2/Pi] Exp[-x^2 / 2] / Erfc[x/Sqrt[2]]
 *
 * exceptions: GSL_EUNDRFLW
 */
int  gsl_sf_hazard_e(MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hazard(MpIeee x);


__END_DECLS

#endif /* __GSL_SF_ERF_H__ */
