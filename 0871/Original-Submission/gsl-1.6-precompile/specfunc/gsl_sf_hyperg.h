#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/gsl_sf_hyperg.h
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

#ifndef __GSL_SF_HYPERG_H__
#define __GSL_SF_HYPERG_H__

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


/* Hypergeometric function related to Bessel functions
 * 0F1[c,x] =
 *            Gamma[c]    x^(1/2(1-c)) I_{c-1}(2 Sqrt[x])
 *            Gamma[c] (-x)^(1/2(1-c)) J_{c-1}(2 Sqrt[-x])
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int  gsl_sf_hyperg_0F1_e(MpIeee c, MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_0F1(const MpIeee c, const MpIeee x);


/* Confluent hypergeometric function  for integer parameters.
 * 1F1[m,n,x] = M(m,n,x)
 *
 * exceptions: 
 */
int  gsl_sf_hyperg_1F1_int_e(const int m, const int n, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_1F1_int(const int m, const int n, MpIeee x);


/* Confluent hypergeometric function.
 * 1F1[a,b,x] = M(a,b,x)
 *
 * exceptions:
 */
int  gsl_sf_hyperg_1F1_e(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_1F1(MpIeee a, MpIeee b, MpIeee x);


/* Confluent hypergeometric function for integer parameters.
 * U(m,n,x)
 *
 * exceptions:
 */
int  gsl_sf_hyperg_U_int_e(const int m, const int n, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_U_int(const int m, const int n, const MpIeee x);


/* Confluent hypergeometric function for integer parameters.
 * U(m,n,x)
 *
 * exceptions:
 */
int  gsl_sf_hyperg_U_int_e10_e(const int m, const int n, const MpIeee x, gsl_sf_result_e10 * result);


/* Confluent hypergeometric function.
 * U(a,b,x)
 *
 * exceptions:
 */
int  gsl_sf_hyperg_U_e(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_U(const MpIeee a, const MpIeee b, const MpIeee x);


/* Confluent hypergeometric function.
 * U(a,b,x)
 *
 * exceptions:
 */
int  gsl_sf_hyperg_U_e10_e(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result_e10 * result);


/* Gauss hypergeometric function 2F1[a,b,c,x]
 * |x| < 1
 *
 * exceptions:
 */
int  gsl_sf_hyperg_2F1_e(MpIeee a, MpIeee b, const MpIeee c, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_2F1(MpIeee a, MpIeee b, MpIeee c, MpIeee x);


/* Gauss hypergeometric function
 * 2F1[aR + I aI, aR - I aI, c, x]
 * |x| < 1
 *
 * exceptions:
 */
int  gsl_sf_hyperg_2F1_conj_e(const MpIeee aR, const MpIeee aI, const MpIeee c, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_2F1_conj(MpIeee aR, MpIeee aI, MpIeee c, MpIeee x);


/* Renormalized Gauss hypergeometric function
 * 2F1[a,b,c,x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int  gsl_sf_hyperg_2F1_renorm_e(const MpIeee a, const MpIeee b, const MpIeee c, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_2F1_renorm(MpIeee a, MpIeee b, MpIeee c, MpIeee x);


/* Renormalized Gauss hypergeometric function
 * 2F1[aR + I aI, aR - I aI, c, x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int  gsl_sf_hyperg_2F1_conj_renorm_e(const MpIeee aR, const MpIeee aI, const MpIeee c, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_2F1_conj_renorm(MpIeee aR, MpIeee aI, MpIeee c, MpIeee x);


/* Mysterious hypergeometric function. The series representation
 * is a divergent hypergeometric series. However, for x < 0 we
 * have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
 *
 * exceptions: GSL_EDOM
 */
int      gsl_sf_hyperg_2F0_e(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result);
MpIeee gsl_sf_hyperg_2F0(const MpIeee a, const MpIeee b, const MpIeee x);


__END_DECLS

#endif /* __GSL_SF_HYPERG_H__ */
