#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_laguerre.h
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

#ifndef __GSL_SF_LAGUERRE_H__
#define __GSL_SF_LAGUERRE_H__

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


/* L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x) */


/* Evaluate generalized Laguerre polynomials
 * using explicit representations.
 *
 * exceptions: none
 */
int  gsl_sf_laguerre_1_e(const MpIeee a, const MpIeee x, gsl_sf_result * result);
int  gsl_sf_laguerre_2_e(const MpIeee a, const MpIeee x, gsl_sf_result * result);
int  gsl_sf_laguerre_3_e(const MpIeee a, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_laguerre_1(MpIeee a, MpIeee x);
MpIeee gsl_sf_laguerre_2(MpIeee a, MpIeee x);
MpIeee gsl_sf_laguerre_3(MpIeee a, MpIeee x);


/* Evaluate generalized Laguerre polynomials.
 *
 * a > -1.0
 * n >= 0
 * exceptions: GSL_EDOM
 */
int      gsl_sf_laguerre_n_e(const int n, const MpIeee a, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_laguerre_n(int  n, MpIeee a, MpIeee x);


__END_DECLS

#endif /* __GSL_SF_LAGUERRE_H__ */
