#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/legendre.h
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

/* Declare private but non-local support functions
 * used in various Legendre function evaluations.
 */

#include <gsl/gsl_sf_result.h>


/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu -> Inf
 * |x| < 1
 */
int
 gsl_sf_conicalP_xlt1_large_neg_mu_e(MpIeee mu, MpIeee tau, MpIeee x,
                                       gsl_sf_result * result, MpIeee * ln_multiplier);


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}, tau -> Inf
 * 1 < x
 */
int
 gsl_sf_conicalP_xgt1_neg_mu_largetau_e(const MpIeee mu, const MpIeee tau,
                                          const MpIeee x, MpIeee acosh_x,
                                          gsl_sf_result * result, MpIeee * ln_multiplier);


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}, tau -> Inf 
 * -1 < x < 1
 */
int
 gsl_sf_conicalP_xlt1_neg_mu_largetau_e(const MpIeee mu, const MpIeee tau,
                                          const MpIeee x, const MpIeee acos_x,
                                          gsl_sf_result * result, MpIeee * ln_multiplier);


/* P^{mu}_{-1/2 + I tau}
 * x->Inf
 *
 *  * This is effective to precision EPS for
 *
 *    (mu^2 + tau^2)/((1 + tau^2)^(1/2) x^2) < EPS^{1/3}
 *
 * since it goes only to a fixed order, based on the
 * representation in terms of hypegeometric functions
 * of argument 1/x^2.
 * [Zhurina+Karmazina, (3.8)]
 */
int
 gsl_sf_conicalP_large_x_e(const MpIeee mu, const MpIeee tau, const MpIeee x,
                             gsl_sf_result * result, MpIeee * ln_multiplier);
