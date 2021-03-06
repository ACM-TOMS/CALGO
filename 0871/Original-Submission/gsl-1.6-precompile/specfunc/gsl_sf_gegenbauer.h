#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_gegenbauer.h
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

#ifndef __GSL_SF_GEGENBAUER_H__
#define __GSL_SF_GEGENBAUER_H__

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


/* Evaluate Gegenbauer polynomials
 * using explicit representations.
 *
 * exceptions: none
 */
int  gsl_sf_gegenpoly_1_e(MpIeee lambda, MpIeee x, gsl_sf_result * result);
int  gsl_sf_gegenpoly_2_e(MpIeee lambda, MpIeee x, gsl_sf_result * result);
int  gsl_sf_gegenpoly_3_e(MpIeee lambda, MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_gegenpoly_1(MpIeee lambda, MpIeee x);
MpIeee gsl_sf_gegenpoly_2(MpIeee lambda, MpIeee x);
MpIeee gsl_sf_gegenpoly_3(MpIeee lambda, MpIeee x);


/* Evaluate Gegenbauer polynomials.
 *
 * lambda > -1/2, n >= 0
 * exceptions: GSL_EDOM
 */
int  gsl_sf_gegenpoly_n_e(int  n, MpIeee lambda, MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_gegenpoly_n(int  n, MpIeee lambda, MpIeee x);


/* Calculate array of Gegenbauer polynomials
 * for n = (0, 1, 2, ... nmax)
 *
 * lambda > -1/2, nmax >= 0
 * exceptions: GSL_EDOM
 */
int  gsl_sf_gegenpoly_array(int  nmax, MpIeee lambda, MpIeee x, MpIeee * result_array);


__END_DECLS

#endif /* __GSL_SF_GEGENBAUER_H__ */
