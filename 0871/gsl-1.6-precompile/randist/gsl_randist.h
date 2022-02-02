#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/gsl_randist.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
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

#ifndef __GSL_RANDIST_H__
#define __GSL_RANDIST_H__
#include <gsl/gsl_rng.h>

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

unsigned int  gsl_ran_bernoulli(const gsl_rng * r, MpIeee p);
MpIeee gsl_ran_bernoulli_pdf(const unsigned int  k, MpIeee p);

MpIeee gsl_ran_beta(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_beta_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

unsigned int  gsl_ran_binomial(const gsl_rng * r, MpIeee p, unsigned int  n);
unsigned int  gsl_ran_binomial_tpe(const gsl_rng * r, MpIeee pp, unsigned int  n);
MpIeee gsl_ran_binomial_pdf(const unsigned int  k, const MpIeee p, const unsigned int  n);

MpIeee gsl_ran_exponential(const gsl_rng * r, const MpIeee mu);
MpIeee gsl_ran_exponential_pdf(const MpIeee x, const MpIeee mu);

MpIeee gsl_ran_exppow(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_exppow_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_ran_cauchy(const gsl_rng * r, const MpIeee a);
MpIeee gsl_ran_cauchy_pdf(const MpIeee x, const MpIeee a);

MpIeee gsl_ran_chisq(const gsl_rng * r, const MpIeee nu);
MpIeee gsl_ran_chisq_pdf(const MpIeee x, const MpIeee nu);

void gsl_ran_dirichlet (const gsl_rng * r, const size_t K, const MpIeee alpha[], MpIeee theta[]);
MpIeee gsl_ran_dirichlet_pdf(const size_t K, const MpIeee alpha[], const MpIeee theta[]);
MpIeee gsl_ran_dirichlet_lnpdf(const size_t K, const MpIeee alpha[], const MpIeee theta[]);

MpIeee gsl_ran_erlang(const gsl_rng * r, const MpIeee a, const MpIeee n);
MpIeee gsl_ran_erlang_pdf(const MpIeee x, const MpIeee a, const MpIeee n);

MpIeee gsl_ran_fdist(const gsl_rng * r, const MpIeee nu1, const MpIeee nu2);
MpIeee gsl_ran_fdist_pdf(const MpIeee x, const MpIeee nu1, const MpIeee nu2);

MpIeee gsl_ran_flat(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_flat_pdf(MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_ran_gamma(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_gamma_int(const gsl_rng * r, const unsigned int  a);
MpIeee gsl_ran_gamma_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_ran_gaussian(const gsl_rng * r, const MpIeee sigma);
MpIeee gsl_ran_gaussian_ratio_method(const gsl_rng * r, const MpIeee sigma);
MpIeee gsl_ran_gaussian_pdf(const MpIeee x, const MpIeee sigma);

MpIeee gsl_ran_ugaussian(const gsl_rng * r);
MpIeee gsl_ran_ugaussian_ratio_method(const gsl_rng * r);
MpIeee gsl_ran_ugaussian_pdf(const MpIeee x);

MpIeee gsl_ran_gaussian_tail(const gsl_rng * r, const MpIeee a, const MpIeee sigma);
MpIeee gsl_ran_gaussian_tail_pdf(const MpIeee x, const MpIeee a, const MpIeee sigma);

MpIeee gsl_ran_ugaussian_tail(const gsl_rng * r, const MpIeee a);
MpIeee gsl_ran_ugaussian_tail_pdf(const MpIeee x, const MpIeee a);

void gsl_ran_bivariate_gaussian (const gsl_rng * r, MpIeee sigma_x, MpIeee sigma_y, MpIeee rho, MpIeee *x, MpIeee *y);
MpIeee gsl_ran_bivariate_gaussian_pdf(const MpIeee x, const MpIeee y, const MpIeee sigma_x, const MpIeee sigma_y, const MpIeee rho);

MpIeee gsl_ran_landau(const gsl_rng * r);
MpIeee gsl_ran_landau_pdf(const MpIeee x);

unsigned int  gsl_ran_geometric(const gsl_rng * r, const MpIeee p);
MpIeee gsl_ran_geometric_pdf(const unsigned int  k, const MpIeee p);

unsigned int  gsl_ran_hypergeometric(const gsl_rng * r, unsigned int  n1, unsigned int  n2, unsigned int  t);
MpIeee gsl_ran_hypergeometric_pdf(const unsigned int  k, const unsigned int  n1, const unsigned int  n2, unsigned int  t);

MpIeee gsl_ran_gumbel1(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_gumbel1_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_ran_gumbel2(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_gumbel2_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

MpIeee gsl_ran_logistic(const gsl_rng * r, const MpIeee a);
MpIeee gsl_ran_logistic_pdf(const MpIeee x, const MpIeee a);

MpIeee gsl_ran_lognormal(const gsl_rng * r, const MpIeee zeta, const MpIeee sigma);
MpIeee gsl_ran_lognormal_pdf(const MpIeee x, const MpIeee zeta, const MpIeee sigma);

unsigned int  gsl_ran_logarithmic(const gsl_rng * r, const MpIeee p);
MpIeee gsl_ran_logarithmic_pdf(const unsigned int  k, const MpIeee p);

void gsl_ran_multinomial (const gsl_rng * r, const size_t K,
                          const unsigned int  N, const MpIeee p[],
                          unsigned int  n[] );
MpIeee gsl_ran_multinomial_pdf(const size_t K,
                                const MpIeee p[], const unsigned int  n[] );
MpIeee gsl_ran_multinomial_lnpdf(const size_t K,
                           const MpIeee p[], const unsigned int  n[] );


unsigned int  gsl_ran_negative_binomial(const gsl_rng * r, MpIeee p, MpIeee n);
MpIeee gsl_ran_negative_binomial_pdf(const unsigned int  k, const MpIeee p, MpIeee n);

unsigned int  gsl_ran_pascal(const gsl_rng * r, MpIeee p, unsigned int  n);
MpIeee gsl_ran_pascal_pdf(const unsigned int  k, const MpIeee p, unsigned int  n);

MpIeee gsl_ran_pareto(const gsl_rng * r, MpIeee a, const MpIeee b);
MpIeee gsl_ran_pareto_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

unsigned int  gsl_ran_poisson(const gsl_rng * r, MpIeee mu);
void gsl_ran_poisson_array (const gsl_rng * r, size_t n, unsigned int  array[],
                            MpIeee mu);
MpIeee gsl_ran_poisson_pdf(const unsigned int  k, const MpIeee mu);

MpIeee gsl_ran_rayleigh(const gsl_rng * r, const MpIeee sigma);
MpIeee gsl_ran_rayleigh_pdf(const MpIeee x, const MpIeee sigma);

MpIeee gsl_ran_rayleigh_tail(const gsl_rng * r, const MpIeee a, const MpIeee sigma);
MpIeee gsl_ran_rayleigh_tail_pdf(const MpIeee x, const MpIeee a, const MpIeee sigma);

MpIeee gsl_ran_tdist(const gsl_rng * r, const MpIeee nu);
MpIeee gsl_ran_tdist_pdf(const MpIeee x, const MpIeee nu);

MpIeee gsl_ran_laplace(const gsl_rng * r, const MpIeee a);
MpIeee gsl_ran_laplace_pdf(const MpIeee x, const MpIeee a);

MpIeee gsl_ran_levy(const gsl_rng * r, const MpIeee c, const MpIeee alpha);
MpIeee gsl_ran_levy_skew(const gsl_rng * r, const MpIeee c, const MpIeee alpha, const MpIeee beta);

MpIeee gsl_ran_weibull(const gsl_rng * r, const MpIeee a, const MpIeee b);
MpIeee gsl_ran_weibull_pdf(const MpIeee x, const MpIeee a, const MpIeee b);

void gsl_ran_dir_2d (const gsl_rng * r, MpIeee * x, MpIeee * y);
void gsl_ran_dir_2d_trig_method (const gsl_rng * r, MpIeee * x, MpIeee * y);
void gsl_ran_dir_3d (const gsl_rng * r, MpIeee * x, MpIeee * y, MpIeee * z);
void gsl_ran_dir_nd (const gsl_rng * r, size_t n, MpIeee * x);

void gsl_ran_shuffle (const gsl_rng * r, void * base, size_t nmembm, size_t size);
int  gsl_ran_choose(const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size) ;
void gsl_ran_sample (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size) ;


typedef struct {                /* struct for Walker algorithm */
    size_t K;
    size_t *A;
    MpIeee *F;
} gsl_ran_discrete_t;

gsl_ran_discrete_t * gsl_ran_discrete_preproc (size_t K, const MpIeee *P);
void gsl_ran_discrete_free(gsl_ran_discrete_t *g);
size_t gsl_ran_discrete (const gsl_rng *r, const gsl_ran_discrete_t *g);
MpIeee gsl_ran_discrete_pdf(size_t k, const gsl_ran_discrete_t *g);


__END_DECLS

#endif /* __GSL_RANDIST_H__ */
