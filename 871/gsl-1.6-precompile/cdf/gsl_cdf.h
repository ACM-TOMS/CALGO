#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/gsl_cdf.h
 * 
 * Copyright (C) 2002 Jason H. Stover.
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA.
 */

/* Author:  J. Stover */

#ifndef __GSL_CDF_H__
#define __GSL_CDF_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS 

MpIeee gsl_cdf_ugaussian_P (const MpIeee x);
MpIeee gsl_cdf_ugaussian_Q(const MpIeee x);

MpIeee gsl_cdf_ugaussian_Pinv(const MpIeee P);
MpIeee gsl_cdf_ugaussian_Qinv(const MpIeee Q);

MpIeee gsl_cdf_gaussian_P(const MpIeee x, const MpIeee sigma);
MpIeee gsl_cdf_gaussian_Q(const MpIeee x, const MpIeee sigma);

MpIeee gsl_cdf_gaussian_Pinv(const MpIeee P, const MpIeee sigma);
MpIeee gsl_cdf_gaussian_Qinv(const MpIeee Q, const MpIeee sigma);

MpIeee gsl_cdf_gamma_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_gamma_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_gamma_Pinv(const MpIeee P, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_gamma_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_cauchy_P(const MpIeee x, const MpIeee a);
MpIeee gsl_cdf_cauchy_Q(const MpIeee x, const MpIeee a);

MpIeee gsl_cdf_cauchy_Pinv(const MpIeee P, const MpIeee a);
MpIeee gsl_cdf_cauchy_Qinv(const MpIeee Q, const MpIeee a);

MpIeee gsl_cdf_laplace_P(const MpIeee x, const MpIeee a);
MpIeee gsl_cdf_laplace_Q(const MpIeee x, const MpIeee a);

MpIeee gsl_cdf_laplace_Pinv(const MpIeee P, const MpIeee a);
MpIeee gsl_cdf_laplace_Qinv(const MpIeee Q, const MpIeee a);

MpIeee gsl_cdf_rayleigh_P(const MpIeee x, const MpIeee sigma);
MpIeee gsl_cdf_rayleigh_Q(const MpIeee x, const MpIeee sigma);

MpIeee gsl_cdf_rayleigh_Pinv(const MpIeee P, const MpIeee sigma);
MpIeee gsl_cdf_rayleigh_Qinv(const MpIeee Q, const MpIeee sigma);

MpIeee gsl_cdf_chisq_P(const MpIeee x, const MpIeee nu);
MpIeee gsl_cdf_chisq_Q(const MpIeee x, const MpIeee nu);

MpIeee gsl_cdf_chisq_Pinv(const MpIeee P, const MpIeee nu);
MpIeee gsl_cdf_chisq_Qinv(const MpIeee Q, const MpIeee nu);

MpIeee gsl_cdf_exponential_P(const MpIeee x, const MpIeee mu);
MpIeee gsl_cdf_exponential_Q(const MpIeee x, const MpIeee mu);

MpIeee gsl_cdf_exponential_Pinv(const MpIeee P, const MpIeee mu);
MpIeee gsl_cdf_exponential_Qinv(const MpIeee Q, const MpIeee mu);

MpIeee gsl_cdf_exppow_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_exppow_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_tdist_P(const MpIeee x, const MpIeee nu);
MpIeee gsl_cdf_tdist_Q(const MpIeee x, const MpIeee nu);

MpIeee gsl_cdf_tdist_Pinv(const MpIeee P, const MpIeee nu);
MpIeee gsl_cdf_tdist_Qinv(const MpIeee Q, const MpIeee nu);

MpIeee gsl_cdf_fdist_P(const MpIeee x, const MpIeee nu1, const MpIeee nu2);
MpIeee gsl_cdf_fdist_Q(const MpIeee x, const MpIeee nu1, const MpIeee nu2);

MpIeee gsl_cdf_beta_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_beta_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_flat_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_flat_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_flat_Pinv(const MpIeee P, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_flat_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_lognormal_P(const MpIeee x, const MpIeee zeta, const MpIeee sigma);
MpIeee gsl_cdf_lognormal_Q(const MpIeee x, const MpIeee zeta, const MpIeee sigma);

MpIeee gsl_cdf_lognormal_Pinv(const MpIeee P, const MpIeee zeta, const MpIeee sigma);
MpIeee gsl_cdf_lognormal_Qinv(const MpIeee Q, const MpIeee zeta, const MpIeee sigma);

MpIeee gsl_cdf_gumbel1_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_gumbel1_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_gumbel1_Pinv(const MpIeee P, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_gumbel1_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_gumbel2_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_gumbel2_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_gumbel2_Pinv(const MpIeee P, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_gumbel2_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_weibull_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_weibull_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_weibull_Pinv(const MpIeee P, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_weibull_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_pareto_P(const MpIeee x, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_pareto_Q(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_pareto_Pinv(const MpIeee P, const MpIeee a, const MpIeee b);
MpIeee gsl_cdf_pareto_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b);

MpIeee gsl_cdf_logistic_P(const MpIeee x, const MpIeee a);
MpIeee gsl_cdf_logistic_Q(const MpIeee x, const MpIeee a);

MpIeee gsl_cdf_logistic_Pinv(const MpIeee P, const MpIeee a);
MpIeee gsl_cdf_logistic_Qinv(const MpIeee Q, const MpIeee a);

__END_DECLS

#endif /* __GSL_CDF_H__ */
