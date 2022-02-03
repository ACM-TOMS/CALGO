#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/test.c
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

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#define N 100000

/* Convient test dimension for multivariant distributions */
#define MULTI_DIM 10


void testMoments (MpIeee(*f) (void), const char *name,
                  MpIeee a, MpIeee b, MpIeee p);
void testPDF (MpIeee(*f) (void), MpIeee(*pdf) (MpIeee), const char *name);
void testDiscretePDF (MpIeee(*f) (void), MpIeee(*pdf) (unsigned int),
                      const char *name);

void test_shuffle (void);
void test_choose (void);
MpIeee test_beta(void);
MpIeee test_beta_pdf(MpIeee x);
MpIeee test_bernoulli(void);
MpIeee test_bernoulli_pdf(unsigned int  n);

MpIeee test_binomial(void);
MpIeee test_binomial_pdf(unsigned int  n);
MpIeee test_binomial_large(void);
MpIeee test_binomial_large_pdf(unsigned int  n);
MpIeee test_binomial_huge(void);
MpIeee test_binomial_huge_pdf(unsigned int  n);

MpIeee test_binomial_tpe(void);
MpIeee test_binomial_tpe_pdf(unsigned int  n);
MpIeee test_binomial_large_tpe(void);
MpIeee test_binomial_large_tpe_pdf(unsigned int  n);
MpIeee test_binomial_huge_tpe(void);
MpIeee test_binomial_huge_tpe_pdf(unsigned int  n);

MpIeee test_cauchy(void);
MpIeee test_cauchy_pdf(MpIeee x);
MpIeee test_chisq(void);
MpIeee test_chisq_pdf(MpIeee x);
MpIeee test_dirichlet(void);
MpIeee test_dirichlet_pdf(MpIeee x);
void test_dirichlet_moments (void);
MpIeee test_discrete1(void);
MpIeee test_discrete1_pdf(unsigned int  n);
MpIeee test_discrete2(void);
MpIeee test_discrete2_pdf(unsigned int  n);
MpIeee test_discrete3(void);
MpIeee test_discrete3_pdf(unsigned int  n);
MpIeee test_erlang(void);
MpIeee test_erlang_pdf(MpIeee x);
MpIeee test_exponential(void);
MpIeee test_exponential_pdf(MpIeee x);
MpIeee test_exppow0(void);
MpIeee test_exppow0_pdf(MpIeee x);
MpIeee test_exppow1(void);
MpIeee test_exppow1_pdf(MpIeee x);
MpIeee test_exppow1a(void);
MpIeee test_exppow1a_pdf(MpIeee x);
MpIeee test_exppow2(void);
MpIeee test_exppow2_pdf(MpIeee x);
MpIeee test_exppow2a(void);
MpIeee test_exppow2a_pdf(MpIeee x);
MpIeee test_fdist(void);
MpIeee test_fdist_pdf(MpIeee x);
MpIeee test_flat(void);
MpIeee test_flat_pdf(MpIeee x);
MpIeee test_gamma(void);
MpIeee test_gamma_pdf(MpIeee x);
MpIeee test_gamma1(void);
MpIeee test_gamma1_pdf(MpIeee x);
MpIeee test_gamma_int(void);
MpIeee test_gamma_int_pdf(MpIeee x);
MpIeee test_gamma_large(void);
MpIeee test_gamma_large_pdf(MpIeee x);
MpIeee test_gaussian(void);
MpIeee test_gaussian_pdf(MpIeee x);
MpIeee test_gaussian_ratio_method(void);
MpIeee test_gaussian_ratio_method_pdf(MpIeee x);
MpIeee test_gaussian_tail(void);
MpIeee test_gaussian_tail_pdf(MpIeee x);
MpIeee test_gaussian_tail1(void);
MpIeee test_gaussian_tail1_pdf(MpIeee x);
MpIeee test_gaussian_tail2(void);
MpIeee test_gaussian_tail2_pdf(MpIeee x);
MpIeee test_ugaussian(void);
MpIeee test_ugaussian_pdf(MpIeee x);
MpIeee test_ugaussian_ratio_method(void);
MpIeee test_ugaussian_ratio_method_pdf(MpIeee x);
MpIeee test_ugaussian_tail(void);
MpIeee test_ugaussian_tail_pdf(MpIeee x);
MpIeee test_bivariate_gaussian1(void);
MpIeee test_bivariate_gaussian1_pdf(MpIeee x);
MpIeee test_bivariate_gaussian2(void);
MpIeee test_bivariate_gaussian2_pdf(MpIeee x);
MpIeee test_bivariate_gaussian3(void);
MpIeee test_bivariate_gaussian3_pdf(MpIeee x);
MpIeee test_bivariate_gaussian4(void);
MpIeee test_bivariate_gaussian4_pdf(MpIeee x);
MpIeee test_gumbel1(void);
MpIeee test_gumbel1_pdf(MpIeee x);
MpIeee test_gumbel2(void);
MpIeee test_gumbel2_pdf(MpIeee x);
MpIeee test_geometric(void);
MpIeee test_geometric_pdf(unsigned int  x);
MpIeee test_geometric1(void);
MpIeee test_geometric1_pdf(unsigned int  x);
MpIeee test_hypergeometric1(void);
MpIeee test_hypergeometric1_pdf(unsigned int  x);
MpIeee test_hypergeometric2(void);
MpIeee test_hypergeometric2_pdf(unsigned int  x);
MpIeee test_hypergeometric3(void);
MpIeee test_hypergeometric3_pdf(unsigned int  x);
MpIeee test_hypergeometric4(void);
MpIeee test_hypergeometric4_pdf(unsigned int  x);
MpIeee test_hypergeometric5(void);
MpIeee test_hypergeometric5_pdf(unsigned int  x);
MpIeee test_hypergeometric6(void);
MpIeee test_hypergeometric6_pdf(unsigned int  x);
MpIeee test_landau(void);
MpIeee test_landau_pdf(MpIeee x);
MpIeee test_levy1(void);
MpIeee test_levy1_pdf(MpIeee x);
MpIeee test_levy2(void);
MpIeee test_levy2_pdf(MpIeee x);
MpIeee test_levy1a(void);
MpIeee test_levy1a_pdf(MpIeee x);
MpIeee test_levy2a(void);
MpIeee test_levy2a_pdf(MpIeee x);
MpIeee test_levy_skew1(void);
MpIeee test_levy_skew1_pdf(MpIeee x);
MpIeee test_levy_skew2(void);
MpIeee test_levy_skew2_pdf(MpIeee x);
MpIeee test_levy_skew1a(void);
MpIeee test_levy_skew1a_pdf(MpIeee x);
MpIeee test_levy_skew2a(void);
MpIeee test_levy_skew2a_pdf(MpIeee x);
MpIeee test_levy_skew1b(void);
MpIeee test_levy_skew1b_pdf(MpIeee x);
MpIeee test_levy_skew2b(void);
MpIeee test_levy_skew2b_pdf(MpIeee x);
MpIeee test_logistic(void);
MpIeee test_logistic_pdf(MpIeee x);
MpIeee test_lognormal(void);
MpIeee test_lognormal_pdf(MpIeee x);
MpIeee test_logarithmic(void);
MpIeee test_logarithmic_pdf(unsigned int  n);
MpIeee test_multinomial(void);
MpIeee test_multinomial_pdf(unsigned int  n);
MpIeee test_multinomial_large(void);
MpIeee test_multinomial_large_pdf(unsigned int  n);
void test_multinomial_moments (void);
MpIeee test_negative_binomial(void);
MpIeee test_negative_binomial_pdf(unsigned int  n);
MpIeee test_pascal(void);
MpIeee test_pascal_pdf(unsigned int  n);
MpIeee test_pareto(void);
MpIeee test_pareto_pdf(MpIeee x);
MpIeee test_poisson(void);
MpIeee test_poisson_pdf(unsigned int  x);
MpIeee test_poisson_large(void);
MpIeee test_poisson_large_pdf(unsigned int  x);
MpIeee test_dir2d(void);
MpIeee test_dir2d_pdf(MpIeee x);
MpIeee test_dir2d_trig_method(void);
MpIeee test_dir2d_trig_method_pdf(MpIeee x);
MpIeee test_dir3dxy(void);
MpIeee test_dir3dxy_pdf(MpIeee x);
MpIeee test_dir3dyz(void);
MpIeee test_dir3dyz_pdf(MpIeee x);
MpIeee test_dir3dzx(void);
MpIeee test_dir3dzx_pdf(MpIeee x);
MpIeee test_rayleigh(void);
MpIeee test_rayleigh_pdf(MpIeee x);
MpIeee test_rayleigh_tail(void);
MpIeee test_rayleigh_tail_pdf(MpIeee x);
MpIeee test_tdist1(void);
MpIeee test_tdist1_pdf(MpIeee x);
MpIeee test_tdist2(void);
MpIeee test_tdist2_pdf(MpIeee x);
MpIeee test_laplace(void);
MpIeee test_laplace_pdf(MpIeee x);
MpIeee test_weibull(void);
MpIeee test_weibull_pdf(MpIeee x);
MpIeee test_weibull1(void);
MpIeee test_weibull1_pdf(MpIeee x);

gsl_rng *r_global;

int
main (void)
{
  gsl_ieee_env_setup ();

  gsl_rng_env_setup ();
  r_global = gsl_rng_alloc (gsl_rng_default);

#define FUNC(x)  test_ ## x,                     "test gsl_ran_" #x
#define FUNC2(x) test_ ## x, test_ ## x ## _pdf, "test gsl_ran_" #x

  test_shuffle ();
  test_choose ();

  testMoments (FUNC (ugaussian), 0.0, 100.0, 0.5);
  testMoments (FUNC (ugaussian), -1.0, 1.0, 0.6826895);
  testMoments (FUNC (ugaussian), 3.0, 3.5, 0.0011172689);
  testMoments (FUNC (ugaussian_tail), 3.0, 3.5, 0.0011172689 / 0.0013498981);
  testMoments (FUNC (exponential), 0.0, 1.0, 1 - exp (-0.5));
  testMoments (FUNC (cauchy), 0.0, 10000.0, 0.5);

  testMoments (FUNC (discrete1), -0.5, 0.5, 0.59);
  testMoments (FUNC (discrete1), 0.5, 1.5, 0.40);
  testMoments (FUNC (discrete1), 1.5, 3.5, 0.01);

  testMoments (FUNC (discrete2), -0.5,  0.5, 1.0/45.0 );
  testMoments (FUNC (discrete2),  8.5,  9.5, 0 );
  
  testMoments (FUNC (discrete3), -0.5, 0.5, 0.05 );
  testMoments (FUNC (discrete3),  0.5, 1.5, 0.05 );
  testMoments (FUNC (discrete3), -0.5, 9.5, 0.5 );

  test_dirichlet_moments ();
  test_multinomial_moments ();

  testPDF (FUNC2 (beta));
  testPDF (FUNC2 (cauchy));
  testPDF (FUNC2 (chisq));
  testPDF (FUNC2 (dirichlet));
  testPDF (FUNC2 (erlang));
  testPDF (FUNC2 (exponential));

  testPDF (FUNC2 (exppow0));
  testPDF (FUNC2 (exppow1));
  testPDF (FUNC2 (exppow1a));
  testPDF (FUNC2 (exppow2));
  testPDF (FUNC2 (exppow2a));

  testPDF (FUNC2 (fdist));
  testPDF (FUNC2 (flat));
  testPDF (FUNC2 (gamma));
  testPDF (FUNC2 (gamma1));
  testPDF (FUNC2 (gamma_int));
  testPDF (FUNC2 (gamma_large));
  testPDF (FUNC2 (gaussian));
  testPDF (FUNC2 (gaussian_ratio_method));
  testPDF (FUNC2 (ugaussian));
  testPDF (FUNC2 (ugaussian_ratio_method));
  testPDF (FUNC2 (gaussian_tail));
  testPDF (FUNC2 (gaussian_tail1));
  testPDF (FUNC2 (gaussian_tail2));
  testPDF (FUNC2 (ugaussian_tail));

  testPDF (FUNC2 (bivariate_gaussian1));
  testPDF (FUNC2 (bivariate_gaussian2));
  testPDF (FUNC2 (bivariate_gaussian3));
  testPDF (FUNC2 (bivariate_gaussian4));

  testPDF (FUNC2 (gumbel1));
  testPDF (FUNC2 (gumbel2));
  testPDF (FUNC2 (landau));
  testPDF (FUNC2 (levy1));
  testPDF (FUNC2 (levy2));
  testPDF (FUNC2 (levy1a));
  testPDF (FUNC2 (levy2a));
  testPDF (FUNC2 (levy_skew1));
  testPDF (FUNC2 (levy_skew2));
  testPDF (FUNC2 (levy_skew1a));
  testPDF (FUNC2 (levy_skew2a));
  testPDF (FUNC2 (levy_skew1b));
  testPDF (FUNC2 (levy_skew2b));
  testPDF (FUNC2 (logistic));
  testPDF (FUNC2 (lognormal));
  testPDF (FUNC2 (pareto));
  testPDF (FUNC2 (rayleigh));
  testPDF (FUNC2 (rayleigh_tail));
  testPDF (FUNC2 (tdist1));
  testPDF (FUNC2 (tdist2));
  testPDF (FUNC2 (laplace));
  testPDF (FUNC2 (weibull));
  testPDF (FUNC2 (weibull1));

  testPDF (FUNC2 (dir2d));
  testPDF (FUNC2 (dir2d_trig_method));
  testPDF (FUNC2 (dir3dxy));
  testPDF (FUNC2 (dir3dyz));
  testPDF (FUNC2 (dir3dzx));

  testDiscretePDF (FUNC2 (discrete1));
  testDiscretePDF (FUNC2 (discrete2));
  testDiscretePDF (FUNC2 (discrete3));
  testDiscretePDF (FUNC2 (poisson));
  testDiscretePDF (FUNC2 (poisson_large));
  testDiscretePDF (FUNC2 (bernoulli));
  testDiscretePDF (FUNC2 (binomial));
  testDiscretePDF (FUNC2 (binomial_tpe));
  testDiscretePDF (FUNC2 (binomial_large));
  testDiscretePDF (FUNC2 (binomial_large_tpe));
  testDiscretePDF (FUNC2 (binomial_huge));
  testDiscretePDF (FUNC2 (binomial_huge_tpe));
  testDiscretePDF (FUNC2 (geometric));
  testDiscretePDF (FUNC2 (geometric1));
  testDiscretePDF (FUNC2 (hypergeometric1));
  testDiscretePDF (FUNC2 (hypergeometric2));
  testDiscretePDF (FUNC2 (hypergeometric3));
  testDiscretePDF (FUNC2 (hypergeometric4));
  testDiscretePDF (FUNC2 (hypergeometric5));
  testDiscretePDF (FUNC2 (hypergeometric6));
  testDiscretePDF (FUNC2 (logarithmic));
  testDiscretePDF (FUNC2 (multinomial));
  testDiscretePDF (FUNC2 (multinomial_large));
  testDiscretePDF (FUNC2 (negative_binomial));
  testDiscretePDF (FUNC2 (pascal));

  exit (gsl_test_summary ());
}

void
test_shuffle (void)
{
  MpIeee count[10][10];
  int  x[10] =  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  int  i;int   j;int   status=  0;

  for (i = 0; i < 10; i++)
    {
      for (j = 0; j < 10; j++)
        {
          count[i][j] = MpIeee( "0" );
        }
    }

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < 10; j++)
        x[j] = j;

      gsl_ran_shuffle (r_global, x, 10, sizeof (int));

      for (j = 0; j < 10; j++)
        count[x[j]][j]++;
    }

  for (i = 0; i < 10; i++)
    {
      for (j = 0; j < 10; j++)
        {
          MpIeee expected=  N / MpIeee( "10.0" );
          MpIeee d=  fabs (count[i][j] - expected);
          MpIeee sigma=  d / sqrt (expected);
          if (sigma > MpIeee( "5" ) && d > MpIeee( "1" ))
            {
              status = 1;
              gsl_test (status,
                        "gsl_ran_shuffle %d,%d (%g observed vs %g expected)",
                        i, j, count[i][j] / N, 0.1);
            }
        }
    }

  gsl_test (status, "gsl_ran_shuffle on {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}");

}

void
test_choose (void)
{
  MpIeee count[10];
  int  x[10] =  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  int  y[3] =  { 0, 1, 2 };
  int  i;int   j;int   status=  0;

  for (i = 0; i < 10; i++)
    {
      count[i] = MpIeee( "0" );
    }

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < 10; j++)
        x[j] = j;

      gsl_ran_choose (r_global, y, 3, x, 10, sizeof (int));

      for (j = 0; j < 3; j++)
        count[y[j]]++;
    }

  for (i = 0; i < 10; i++)
    {
      MpIeee expected=  MpIeee( "3.0" ) * N / MpIeee( "10.0" );
      MpIeee d=  fabs (count[i] - expected);
      MpIeee sigma=  d / sqrt (expected);
      if (sigma > MpIeee( "5" ) && d > MpIeee( "1" ))
        {
          status = 1;
          gsl_test (status,
                    "gsl_ran_choose %d (%g observed vs %g expected)",
                    i, count[i] / N, 0.1);
        }
    }

  gsl_test (status, "gsl_ran_choose (3) on {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}");

}




void
testMoments (MpIeee(*f) (void), const char *name,
             MpIeee a, MpIeee b, MpIeee p)
{
  int  i;
  MpIeee count=  MpIeee( "0" );MpIeee  expected;MpIeee  sigma;
  int  status;

  for (i = 0; i < N; i++)
    {
      MpIeee r=  f ();
      if (r < b && r > a)
        count++;
    }

  expected = p * N;
  sigma = fabs (count - expected) / sqrt (expected);

  status = (sigma > 3);

  gsl_test (status, "%s [%g,%g] (%g observed vs %g expected)",
            name, a, b, count / N, p);
}

#define BINS 100

void
testPDF (MpIeee(*f) (void), MpIeee(*pdf) (MpIeee), const char *name)
{
  MpIeee count[BINS];MpIeee  p[BINS];
  MpIeee a=  -MpIeee( "5.0" );MpIeee  b=  +MpIeee( "5.0" );
  MpIeee dx=  (b - a) / BINS;
  int  i;int   j;int   status=  0;int   status_i=  0;

  for (i = 0; i < BINS; i++)
    count[i] = MpIeee( "0" );

  for (i = 0; i < N; i++)
    {
      MpIeee r=  f ();
      if (r < b && r > a)
        {
          j = (int) ((r - a) / dx);
          count[j]++;
        }
    }

  for (i = 0; i < BINS; i++)
    {
      /* Compute an approximation to the integral of p(x) from x to
         x+dx using Simpson's rule */

      MpIeee x=  a + i * dx;
#define STEPS 100
      MpIeee sum=  MpIeee( "0" );

      if (fabs (x) < 1e-10)     /* hit the origin exactly */
        x = MpIeee( "0.0" );

      for (j = 1; j < STEPS; j++)
        sum += pdf (x + j * dx / STEPS);

      p[i] = MpIeee( "0.5" ) * (pdf (x) + MpIeee( "2" ) * sum + pdf (x + dx - MpIeee( "1" )e-MpIeee( "7" ))) * dx / STEPS;
    }

  for (i = 0; i < BINS; i++)
    {
      MpIeee x=  a + i * dx;
      MpIeee d=  fabs (count[i] - N * p[i]);
      if (p[i] != MpIeee( "0" ))
        {
          MpIeee s=  d / sqrt (N * p[i]);
          status_i = (s > 5) && (d > 1);
        }
      else
        {
          status_i = (count[i] != 0);
        }
      status |= status_i;
      if (status_i)
        gsl_test (status_i, "%s [%g,%g) (%g/%d=%g observed vs %g expected)",
                  name, x, x + dx, count[i], N, count[i] / N, p[i]);
    }

  if (status == 0)
    gsl_test (status, "%s, sampling against pdf over range [%g,%g) ",
              name, a, b);
}

void
testDiscretePDF (MpIeee(*f) (void), MpIeee(*pdf) (unsigned int),
                 const char *name)
{
  MpIeee count[BINS];MpIeee  p[BINS];
  unsigned int  i;
  int  status=  0;int   status_i=  0;

  for (i = 0; i < BINS; i++)
    count[i] = MpIeee( "0" );

  for (i = 0; i < N; i++)
    {
      int  r=  (int) (f ());
      if (r >= 0 && r < BINS)
        count[r]++;
    }

  for (i = 0; i < BINS; i++)
    p[i] = pdf (i);

  for (i = 0; i < BINS; i++)
    {
      MpIeee d=  fabs (count[i] - N * p[i]);
      if (p[i] != MpIeee( "0" ))
        {
          MpIeee s=  d / sqrt (N * p[i]);
          status_i = (s > 5) && (d > 1);
        }
      else
        {
          status_i = (count[i] != 0);
        }
      status |= status_i;
      if (status_i)
        gsl_test (status_i, "%s i=%d (%g observed vs %g expected)",
                  name, i, count[i] / N, p[i]);
    }

  if (status == 0)
    gsl_test (status, "%s, sampling against pdf over range [%d,%d) ",
              name, 0, BINS);
}



MpIeee test_beta(void)
{
  return gsl_ran_beta (r_global, MpIeee( "2.0" ), MpIeee( "3.0" ));
}

MpIeee test_beta_pdf(MpIeee x)
{
  return gsl_ran_beta_pdf (x, MpIeee( "2.0" ), MpIeee( "3.0" ));
}

MpIeee test_bernoulli(void)
{
  return gsl_ran_bernoulli (r_global, MpIeee( "0.3" ));
}

MpIeee test_bernoulli_pdf(unsigned int  n)
{
  return gsl_ran_bernoulli_pdf (n, MpIeee( "0.3" ));
}

MpIeee test_binomial(void)
{
  return gsl_ran_binomial (r_global, MpIeee( "0.3" ), MpIeee( "5" ));
}

MpIeee test_binomial_pdf(unsigned int  n)
{
  return gsl_ran_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "5" ));
}

MpIeee test_binomial_tpe(void)
{
  return gsl_ran_binomial_tpe (r_global, MpIeee( "0.3" ), MpIeee( "5" ));
}

MpIeee test_binomial_tpe_pdf(unsigned int  n)
{
  return gsl_ran_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "5" ));
}


MpIeee test_binomial_large(void)
{
  return gsl_ran_binomial (r_global, MpIeee( "0.3" ), MpIeee( "55" ));
}

MpIeee test_binomial_large_pdf(unsigned int  n)
{
  return gsl_ran_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "55" ));
}

MpIeee test_binomial_large_tpe(void)
{
  return gsl_ran_binomial_tpe (r_global, MpIeee( "0.3" ), MpIeee( "55" ));
}

MpIeee test_binomial_large_tpe_pdf(unsigned int  n)
{
  return gsl_ran_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "55" ));
}


MpIeee test_binomial_huge(void)
{
  return gsl_ran_binomial (r_global, MpIeee( "0.3" ), MpIeee( "5500" ));
}

MpIeee test_binomial_huge_pdf(unsigned int  n)
{
  return gsl_ran_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "5500" ));
}

MpIeee test_binomial_huge_tpe(void)
{
  return gsl_ran_binomial_tpe (r_global, MpIeee( "0.3" ), MpIeee( "5500" ));
}

MpIeee test_binomial_huge_tpe_pdf(unsigned int  n)
{
  return gsl_ran_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "5500" ));
}

MpIeee test_cauchy(void)
{
  return gsl_ran_cauchy (r_global, MpIeee( "2.0" ));
}

MpIeee test_cauchy_pdf(MpIeee x)
{
  return gsl_ran_cauchy_pdf (x, MpIeee( "2.0" ));
}

MpIeee test_chisq(void)
{
  return gsl_ran_chisq (r_global, MpIeee( "13.0" ));
}

MpIeee test_chisq_pdf(MpIeee x)
{
  return gsl_ran_chisq_pdf (x, MpIeee( "13.0" ));
}

MpIeee test_dir2d(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );MpIeee  theta;
  gsl_ran_dir_2d (r_global, &x, &y);
  theta = atan2 (x, y);
  return theta;
}

MpIeee test_dir2d_pdf(MpIeee x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return MpIeee( "1" ) / (MpIeee( "2" ) * M_PI);
    }
  else
    {
      return MpIeee( "0" );
    }
}

MpIeee test_dir2d_trig_method(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );MpIeee  theta;
  gsl_ran_dir_2d_trig_method (r_global, &x, &y);
  theta = atan2 (x, y);
  return theta;
}

MpIeee test_dir2d_trig_method_pdf(MpIeee x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return MpIeee( "1" ) / (MpIeee( "2" ) * M_PI);
    }
  else
    {
      return MpIeee( "0" );
    }
}

MpIeee test_dir3dxy(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );MpIeee  z=  MpIeee( "0" );MpIeee  theta;
  gsl_ran_dir_3d (r_global, &x, &y, &z);
  theta = atan2 (x, y);
  return theta;
}

MpIeee test_dir3dxy_pdf(MpIeee x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return MpIeee( "1" ) / (MpIeee( "2" ) * M_PI);
    }
  else
    {
      return MpIeee( "0" );
    }
}

MpIeee test_dir3dyz(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );MpIeee  z=  MpIeee( "0" );MpIeee  theta;
  gsl_ran_dir_3d (r_global, &x, &y, &z);
  theta = atan2 (y, z);
  return theta;
}

MpIeee test_dir3dyz_pdf(MpIeee x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return MpIeee( "1" ) / (MpIeee( "2" ) * M_PI);
    }
  else
    {
      return MpIeee( "0" );
    }
}

MpIeee test_dir3dzx(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );MpIeee  z=  MpIeee( "0" );MpIeee  theta;
  gsl_ran_dir_3d (r_global, &x, &y, &z);
  theta = atan2 (z, x);
  return theta;
}

MpIeee test_dir3dzx_pdf(MpIeee x)
{
  if (x > -M_PI && x <= M_PI)
    {
      return MpIeee( "1" ) / (MpIeee( "2" ) * M_PI);
    }
  else
    {
      return MpIeee( "0" );
    }
}

MpIeee test_dirichlet(void)
{
  /* This is a bit of a lame test, since when K=2, the Dirichlet distribution
     becomes a beta distribution */
  size_t K = 2;
  MpIeee alpha[2] =  { MpIeee( "2.5" ), MpIeee( "5.0" ) };
  MpIeee theta[2] =  { MpIeee( "0.0" ), MpIeee( "0.0" ) };

  gsl_ran_dirichlet (r_global, K, alpha, theta);

  return theta[0];
}

MpIeee test_dirichlet_pdf(MpIeee x)
{
  size_t K = 2;
  MpIeee alpha[2] =  { MpIeee( "2.5" ), MpIeee( "5.0" ) };
  MpIeee theta[2];

  if (x <= MpIeee( "0.0" ) || x >= MpIeee( "1.0" ))
    return MpIeee( "0.0" );                 /* Out of range */

  theta[0] = x;
  theta[1] = MpIeee( "1.0" ) - x;

  return gsl_ran_dirichlet_pdf (K, alpha, theta);
}


/* Check that the observed means of the Dirichlet variables are
   within reasonable statistical errors of their correct values. */

#define DIRICHLET_K 10

void
test_dirichlet_moments (void)
{
  MpIeee alpha[DIRICHLET_K];
  MpIeee theta[DIRICHLET_K];
  MpIeee theta_sum[DIRICHLET_K];

  MpIeee alpha_sum=  MpIeee( "0.0" );
  MpIeee mean;MpIeee  obs_mean;MpIeee  sd;MpIeee  sigma;
  int  status;int   k;int   n;

  for (k = 0; k < DIRICHLET_K; k++)
    {
      alpha[k] = gsl_ran_exponential (r_global, MpIeee( "0.1" ));
      alpha_sum += alpha[k];
      theta_sum[k] = MpIeee( "0.0" );
    }

  for (n = 0; n < N; n++)
    {
      gsl_ran_dirichlet (r_global, DIRICHLET_K, alpha, theta);
      for (k = 0; k < DIRICHLET_K; k++)
        theta_sum[k] += theta[k];
    }

  for (k = 0; k < DIRICHLET_K; k++)
    {
      mean = alpha[k] / alpha_sum;
      sd =
        sqrt ((alpha[k] * (MpIeee( "1." ) - alpha[k] / alpha_sum)) /
              (alpha_sum * (alpha_sum + MpIeee( "1." ))));
      obs_mean = theta_sum[k] / N;
      sigma = sqrt ((MpIeee) N) * fabs (mean - obs_mean) / sd;

      status = (sigma > 3.0);

      gsl_test (status,
                "test gsl_ran_dirichlet: mean (%g observed vs %g expected)",
                obs_mean, mean);
    }
}


/* Check that the observed means of the multinomial variables are
   within reasonable statistical errors of their correct values. */

void
test_multinomial_moments (void)
{
  const unsigned int  sum_n=  100;

  const MpIeee p[MULTI_DIM] = { 0.2, 0.20, 0.17, 0.14, 0.12,
                               0.07, 0.05, 0.02, 0.02, 0.01 };

  unsigned int   x[MULTI_DIM];
  MpIeee x_sum[MULTI_DIM];

  MpIeee mean;MpIeee  obs_mean;MpIeee  sd;MpIeee  sigma;
  int  status;int   k;int   n;

  for (k = 0; k < MULTI_DIM; k++)
    x_sum[k] =MpIeee( "0.0" );

  for (n = 0; n < N; n++)
    {
      gsl_ran_multinomial (r_global, MULTI_DIM, sum_n, p, x);
      for (k = 0; k < MULTI_DIM; k++)
        x_sum[k] += x[k];
    }

  for (k = 0; k < MULTI_DIM; k++)
    {
      mean = p[k] * sum_n;
      sd = p[k] * (MpIeee( "1." )-p[k]) * sum_n;

      obs_mean = x_sum[k] / N;
      sigma = sqrt ((MpIeee) N) * fabs (mean - obs_mean) / sd;

      status = (sigma > 3.0);

      gsl_test (status,
                "test gsl_ran_multinomial: mean (%g observed vs %g expected)",
                obs_mean, mean);
    }
}


static gsl_ran_discrete_t *g1 = NULL;
static gsl_ran_discrete_t *g2 = NULL;
static gsl_ran_discrete_t *g3 = NULL;

MpIeee test_discrete1(void)
{
  static MpIeee P[3] =  { MpIeee( "0.59" ), MpIeee( "0.4" ), MpIeee( "0.01" ) };
  if (g1 == NULL)
    {
      g1 = gsl_ran_discrete_preproc (3, P);
    }
  return gsl_ran_discrete (r_global, g1);
}

MpIeee test_discrete1_pdf(unsigned int  n)
{
  return gsl_ran_discrete_pdf ((size_t) n, g1);
}

MpIeee test_discrete2(void)
{
  static MpIeee P[10] =  { MpIeee( "1" ), MpIeee( "9" ), MpIeee( "3" ), MpIeee( "4" ), MpIeee( "5" ), MpIeee( "8" ), MpIeee( "6" ), MpIeee( "7" ), MpIeee( "2" ), MpIeee( "0" ) };
  if (g2 == NULL)
    {
      g2 = gsl_ran_discrete_preproc (10, P);
    }
  return gsl_ran_discrete (r_global, g2);
}

MpIeee test_discrete2_pdf(unsigned int  n)
{
  return gsl_ran_discrete_pdf ((size_t) n, g2);
}
MpIeee test_discrete3(void)
{
  static MpIeee P[20];
  if (g3 == NULL)
    { int  i;
      for (i=0; i<20; ++i) P[i]=MpIeee( "1.0" )/MpIeee( "20" );
      g3 = gsl_ran_discrete_preproc (20, P);
    }
  return gsl_ran_discrete (r_global, g3);
}

MpIeee test_discrete3_pdf(unsigned int  n)
{
  return gsl_ran_discrete_pdf ((size_t) n, g3);
}


MpIeee test_erlang(void)
{
  return gsl_ran_erlang (r_global, MpIeee( "3.0" ), MpIeee( "4.0" ));
}

MpIeee test_erlang_pdf(MpIeee x)
{
  return gsl_ran_erlang_pdf (x, MpIeee( "3.0" ), MpIeee( "4.0" ));
}

MpIeee test_exponential(void)
{
  return gsl_ran_exponential (r_global, MpIeee( "2.0" ));
}

MpIeee test_exponential_pdf(MpIeee x)
{
  return gsl_ran_exponential_pdf (x, MpIeee( "2.0" ));
}

MpIeee test_exppow0(void)
{
  return gsl_ran_exppow (r_global, MpIeee( "3.7" ), MpIeee( "0.3" ));
}

MpIeee test_exppow0_pdf(MpIeee x)
{
  return gsl_ran_exppow_pdf (x, MpIeee( "3.7" ), MpIeee( "0.3" ));
}

MpIeee test_exppow1(void)
{
  return gsl_ran_exppow (r_global, MpIeee( "3.7" ), MpIeee( "1.0" ));
}

MpIeee test_exppow1_pdf(MpIeee x)
{
  return gsl_ran_exppow_pdf (x, MpIeee( "3.7" ), MpIeee( "1.0" ));
}

MpIeee test_exppow1a(void)
{
  return gsl_ran_exppow (r_global, MpIeee( "3.7" ), MpIeee( "1.9" ));
}

MpIeee test_exppow1a_pdf(MpIeee x)
{
  return gsl_ran_exppow_pdf (x, MpIeee( "3.7" ), MpIeee( "1.9" ));
}

MpIeee test_exppow2(void)
{
  return gsl_ran_exppow (r_global, MpIeee( "3.7" ), MpIeee( "2.0" ));
}

MpIeee test_exppow2_pdf(MpIeee x)
{
  return gsl_ran_exppow_pdf (x, MpIeee( "3.7" ), MpIeee( "2.0" ));
}


MpIeee test_exppow2a(void)
{
  return gsl_ran_exppow (r_global, MpIeee( "3.7" ), MpIeee( "7.5" ));
}

MpIeee test_exppow2a_pdf(MpIeee x)
{
  return gsl_ran_exppow_pdf (x, MpIeee( "3.7" ), MpIeee( "7.5" ));
}

MpIeee test_fdist(void)
{
  return gsl_ran_fdist (r_global, MpIeee( "3.0" ), MpIeee( "4.0" ));
}

MpIeee test_fdist_pdf(MpIeee x)
{
  return gsl_ran_fdist_pdf (x, MpIeee( "3.0" ), MpIeee( "4.0" ));
}

MpIeee test_flat(void)
{
  return gsl_ran_flat (r_global, MpIeee( "3.0" ), MpIeee( "4.0" ));
}

MpIeee test_flat_pdf(MpIeee x)
{
  return gsl_ran_flat_pdf (x, MpIeee( "3.0" ), MpIeee( "4.0" ));
}

MpIeee test_gamma(void)
{
  return gsl_ran_gamma (r_global, MpIeee( "2.5" ), MpIeee( "2.17" ));
}

MpIeee test_gamma_pdf(MpIeee x)
{
  return gsl_ran_gamma_pdf (x, MpIeee( "2.5" ), MpIeee( "2.17" ));
}

MpIeee test_gamma1(void)
{
  return gsl_ran_gamma (r_global, MpIeee( "1.0" ), MpIeee( "2.17" ));
}

MpIeee test_gamma1_pdf(MpIeee x)
{
  return gsl_ran_gamma_pdf (x, MpIeee( "1.0" ), MpIeee( "2.17" ));
}


MpIeee test_gamma_int(void)
{
  return gsl_ran_gamma (r_global, MpIeee( "10.0" ), MpIeee( "2.17" ));
}

MpIeee test_gamma_int_pdf(MpIeee x)
{
  return gsl_ran_gamma_pdf (x, MpIeee( "10.0" ), MpIeee( "2.17" ));
}


MpIeee test_gamma_large(void)
{
  return gsl_ran_gamma (r_global, MpIeee( "20.0" ), MpIeee( "2.17" ));
}

MpIeee test_gamma_large_pdf(MpIeee x)
{
  return gsl_ran_gamma_pdf (x, MpIeee( "20.0" ), MpIeee( "2.17" ));
}


MpIeee test_gaussian(void)
{
  return gsl_ran_gaussian (r_global, MpIeee( "3.0" ));
}

MpIeee test_gaussian_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, MpIeee( "3.0" ));
}

MpIeee test_gaussian_ratio_method(void)
{
  return gsl_ran_gaussian_ratio_method (r_global, MpIeee( "3.0" ));
}

MpIeee test_gaussian_ratio_method_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, MpIeee( "3.0" ));
}

MpIeee test_gaussian_tail(void)
{
  return gsl_ran_gaussian_tail (r_global, MpIeee( "1.7" ), MpIeee( "0.25" ));
}

MpIeee test_gaussian_tail_pdf(MpIeee x)
{
  return gsl_ran_gaussian_tail_pdf (x, MpIeee( "1.7" ), MpIeee( "0.25" ));
}

MpIeee test_gaussian_tail1(void)
{
  return gsl_ran_gaussian_tail (r_global, -MpIeee( "1.7" ), MpIeee( "5.0" ));
}

MpIeee test_gaussian_tail1_pdf(MpIeee x)
{
  return gsl_ran_gaussian_tail_pdf (x, -MpIeee( "1.7" ), MpIeee( "5.0" ));
}

MpIeee test_gaussian_tail2(void)
{
  return gsl_ran_gaussian_tail (r_global, MpIeee( "0.1" ), MpIeee( "2.0" ));
}

MpIeee test_gaussian_tail2_pdf(MpIeee x)
{
  return gsl_ran_gaussian_tail_pdf (x, MpIeee( "0.1" ), MpIeee( "2.0" ));
}


MpIeee test_ugaussian(void)
{
  return gsl_ran_ugaussian (r_global);
}

MpIeee test_ugaussian_pdf(MpIeee x)
{
  return gsl_ran_ugaussian_pdf (x);
}

MpIeee test_ugaussian_ratio_method(void)
{
  return gsl_ran_ugaussian_ratio_method (r_global);
}

MpIeee test_ugaussian_ratio_method_pdf(MpIeee x)
{
  return gsl_ran_ugaussian_pdf (x);
}

MpIeee test_ugaussian_tail(void)
{
  return gsl_ran_ugaussian_tail (r_global, MpIeee( "3.0" ));
}

MpIeee test_ugaussian_tail_pdf(MpIeee x)
{
  return gsl_ran_ugaussian_tail_pdf (x, MpIeee( "3.0" ));
}

MpIeee test_bivariate_gaussian1(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, 0.3, &x, &y);
  return x;
}

MpIeee test_bivariate_gaussian1_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, MpIeee( "3.0" ));
}

MpIeee test_bivariate_gaussian2(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, 0.3, &x, &y);
  return y;
}

MpIeee test_bivariate_gaussian2_pdf(MpIeee y)
{
  int  i;int   n=  10;
  MpIeee sum=  MpIeee( "0" );
  MpIeee a=  -MpIeee( "10" );MpIeee  b=  MpIeee( "10" );MpIeee  dx=  (b - a) / n;
  for (i = 0; i < n; i++)
    {
      MpIeee x=  a + i * dx;
      sum += gsl_ran_bivariate_gaussian_pdf (x, y, MpIeee( "3.0" ), MpIeee( "2.0" ), MpIeee( "0.3" )) * dx;
    }
  return sum;
}


MpIeee test_bivariate_gaussian3(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, 0.3, &x, &y);
  return x + y;
}

MpIeee test_bivariate_gaussian3_pdf(MpIeee x)
{
  MpIeee sx=  MpIeee( "3.0" );MpIeee  sy=  MpIeee( "2.0" );MpIeee  r=  MpIeee( "0.3" );
  MpIeee su=  (sx + r * sy);
  MpIeee sv=  sy * sqrt (MpIeee( "1" ) - r * r);
  MpIeee sigma=  sqrt (su * su + sv * sv);

  return gsl_ran_gaussian_pdf (x, sigma);
}

MpIeee test_bivariate_gaussian4(void)
{
  MpIeee x=  MpIeee( "0" );MpIeee  y=  MpIeee( "0" );
  gsl_ran_bivariate_gaussian (r_global, 3.0, 2.0, -0.5, &x, &y);
  return x + y;
}

MpIeee test_bivariate_gaussian4_pdf(MpIeee x)
{
  MpIeee sx=  MpIeee( "3.0" );MpIeee  sy=  MpIeee( "2.0" );MpIeee  r=  -MpIeee( "0.5" );
  MpIeee su=  (sx + r * sy);
  MpIeee sv=  sy * sqrt (MpIeee( "1" ) - r * r);
  MpIeee sigma=  sqrt (su * su + sv * sv);

  return gsl_ran_gaussian_pdf (x, sigma);
}


MpIeee test_geometric(void)
{
  return gsl_ran_geometric (r_global, MpIeee( "0.5" ));
}

MpIeee test_geometric_pdf(unsigned int  n)
{
  return gsl_ran_geometric_pdf (n, MpIeee( "0.5" ));
}

MpIeee test_geometric1(void)
{
  return gsl_ran_geometric (r_global, MpIeee( "1.0" ));
}

MpIeee test_geometric1_pdf(unsigned int  n)
{
  return gsl_ran_geometric_pdf (n, MpIeee( "1.0" ));
}

MpIeee test_hypergeometric1(void)
{
  return gsl_ran_hypergeometric (r_global, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "4" ));
}

MpIeee test_hypergeometric1_pdf(unsigned int  n)
{
  return gsl_ran_hypergeometric_pdf (n, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "4" ));
}


MpIeee test_hypergeometric2(void)
{
  return gsl_ran_hypergeometric (r_global, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "11" ));
}

MpIeee test_hypergeometric2_pdf(unsigned int  n)
{
  return gsl_ran_hypergeometric_pdf (n, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "11" ));
}

MpIeee test_hypergeometric3(void)
{
  return gsl_ran_hypergeometric (r_global, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "1" ));
}

MpIeee test_hypergeometric3_pdf(unsigned int  n)
{
  return gsl_ran_hypergeometric_pdf (n, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "1" ));
}

MpIeee test_hypergeometric4(void)
{
  return gsl_ran_hypergeometric (r_global, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "20" ));
}

MpIeee test_hypergeometric4_pdf(unsigned int  n)
{
  return gsl_ran_hypergeometric_pdf (n, MpIeee( "5" ), MpIeee( "7" ), MpIeee( "20" ));
}

MpIeee test_hypergeometric5(void)
{
  return gsl_ran_hypergeometric (r_global, MpIeee( "2" ), MpIeee( "7" ), MpIeee( "5" ));
}

MpIeee test_hypergeometric5_pdf(unsigned int  n)
{
  return gsl_ran_hypergeometric_pdf (n, MpIeee( "2" ), MpIeee( "7" ), MpIeee( "5" ));
}


MpIeee test_hypergeometric6(void)
{
  return gsl_ran_hypergeometric (r_global, MpIeee( "2" ), MpIeee( "10" ), MpIeee( "3" ));
}

MpIeee test_hypergeometric6_pdf(unsigned int  n)
{
  return gsl_ran_hypergeometric_pdf (n, MpIeee( "2" ), MpIeee( "10" ), MpIeee( "3" ));
}




MpIeee test_gumbel1(void)
{
  return gsl_ran_gumbel1 (r_global, MpIeee( "3.12" ), MpIeee( "4.56" ));
}

MpIeee test_gumbel1_pdf(MpIeee x)
{
  return gsl_ran_gumbel1_pdf (x, MpIeee( "3.12" ), MpIeee( "4.56" ));
}

MpIeee test_gumbel2(void)
{
  return gsl_ran_gumbel2 (r_global, MpIeee( "3.12" ), MpIeee( "4.56" ));
}

MpIeee test_gumbel2_pdf(MpIeee x)
{
  return gsl_ran_gumbel2_pdf (x, MpIeee( "3.12" ), MpIeee( "4.56" ));
}

MpIeee test_landau(void)
{
  return gsl_ran_landau (r_global);
}

MpIeee test_landau_pdf(MpIeee x)
{
  return gsl_ran_landau_pdf (x);
}

MpIeee test_levy1(void)
{
  return gsl_ran_levy (r_global, MpIeee( "5.0" ), MpIeee( "1.0" ));
}

MpIeee test_levy1_pdf(MpIeee x)
{
  return gsl_ran_cauchy_pdf (x, MpIeee( "5.0" ));
}

MpIeee test_levy2(void)
{
  return gsl_ran_levy (r_global, MpIeee( "5.0" ), MpIeee( "2.0" ));
}

MpIeee test_levy2_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, sqrt (MpIeee( "2.0" )) * MpIeee( "5.0" ));
}

MpIeee test_levy1a(void)
{
  return gsl_ran_levy (r_global, MpIeee( "5.0" ), MpIeee( "1.01" ));
}

MpIeee test_levy1a_pdf(MpIeee x)
{
  return gsl_ran_cauchy_pdf (x, MpIeee( "5.0" ));
}

MpIeee test_levy2a(void)
{
  return gsl_ran_levy (r_global, MpIeee( "5.0" ), MpIeee( "1.99" ));
}

MpIeee test_levy2a_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, sqrt (MpIeee( "2.0" )) * MpIeee( "5.0" ));
}


MpIeee test_levy_skew1(void)
{
  return gsl_ran_levy_skew (r_global, MpIeee( "5.0" ), MpIeee( "1.0" ), MpIeee( "0.0" ));
}

MpIeee test_levy_skew1_pdf(MpIeee x)
{
  return gsl_ran_cauchy_pdf (x, MpIeee( "5.0" ));
}

MpIeee test_levy_skew2(void)
{
  return gsl_ran_levy_skew (r_global, MpIeee( "5.0" ), MpIeee( "2.0" ), MpIeee( "0.0" ));
}

MpIeee test_levy_skew2_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, sqrt (MpIeee( "2.0" )) * MpIeee( "5.0" ));
}

MpIeee test_levy_skew1a(void)
{
  return gsl_ran_levy_skew (r_global, MpIeee( "5.0" ), MpIeee( "1.01" ), MpIeee( "0.0" ));
}

MpIeee test_levy_skew1a_pdf(MpIeee x)
{
  return gsl_ran_cauchy_pdf (x, MpIeee( "5.0" ));
}

MpIeee test_levy_skew2a(void)
{
  return gsl_ran_levy_skew (r_global, MpIeee( "5.0" ), MpIeee( "1.99" ), MpIeee( "0.0" ));
}

MpIeee test_levy_skew2a_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, sqrt (MpIeee( "2.0" )) * MpIeee( "5.0" ));
}

MpIeee test_levy_skew1b(void)
{
  return gsl_ran_levy_skew (r_global, MpIeee( "5.0" ), MpIeee( "1.01" ), MpIeee( "0.001" ));
}

MpIeee test_levy_skew1b_pdf(MpIeee x)
{
  return gsl_ran_cauchy_pdf (x, MpIeee( "5.0" ));
}

MpIeee test_levy_skew2b(void)
{
  return gsl_ran_levy_skew (r_global, MpIeee( "5.0" ), MpIeee( "1.99" ), MpIeee( "0.001" ));
}

MpIeee test_levy_skew2b_pdf(MpIeee x)
{
  return gsl_ran_gaussian_pdf (x, sqrt (MpIeee( "2.0" )) * MpIeee( "5.0" ));
}


MpIeee test_logistic(void)
{
  return gsl_ran_logistic (r_global, MpIeee( "3.1" ));
}

MpIeee test_logistic_pdf(MpIeee x)
{
  return gsl_ran_logistic_pdf (x, MpIeee( "3.1" ));
}

MpIeee test_logarithmic(void)
{
  return gsl_ran_logarithmic (r_global, MpIeee( "0.4" ));
}

MpIeee test_logarithmic_pdf(unsigned int  n)
{
  return gsl_ran_logarithmic_pdf (n, MpIeee( "0.4" ));
}


MpIeee test_lognormal(void)
{
  return gsl_ran_lognormal (r_global, MpIeee( "2.7" ), MpIeee( "1.3" ));
}

MpIeee test_lognormal_pdf(MpIeee x)
{
  return gsl_ran_lognormal_pdf (x, MpIeee( "2.7" ), MpIeee( "1.3" ));
}

MpIeee test_multinomial(void)
{
  const size_t K = 3;
  const unsigned int  sum_n=  BINS;
  unsigned int  n[3];
  /* Test use of weights instead of probabilities. */
  const MpIeee p[] =  { 2., 7., 1.};

  gsl_ran_multinomial ( r_global, K, sum_n, p, n);

  return n[0];
}

MpIeee test_multinomial_pdf(unsigned int  n_0)
{
  /* The margional distribution of just 1 variate  is binomial. */
  size_t K = 2;
  /* Test use of weights instead of probabilities */
  MpIeee p[] =  { MpIeee( "0.4" ), MpIeee( "1.6" )};
  const unsigned int  sum_n=  BINS;
  unsigned int  n[2];

  n[0] = n_0;
  n[1] =sum_n - n_0;

  return gsl_ran_multinomial_pdf (K, p, n);
}


MpIeee test_multinomial_large(void)
{
  const unsigned int  sum_n=  BINS;
  unsigned int  n[MULTI_DIM];
  const MpIeee p[MULTI_DIM] =  { 0.2, 0.20, 0.17, 0.14, 0.12,
                                0.07, 0.05, 0.04, 0.01, 0.00  };

  gsl_ran_multinomial ( r_global, MULTI_DIM, sum_n, p, n);

  return n[0];
}

MpIeee test_multinomial_large_pdf(unsigned int  n_0)
{
  return test_multinomial_pdf(n_0);
}

MpIeee test_negative_binomial(void)
{
  return gsl_ran_negative_binomial (r_global, MpIeee( "0.3" ), MpIeee( "20.0" ));
}

MpIeee test_negative_binomial_pdf(unsigned int  n)
{
  return gsl_ran_negative_binomial_pdf (n, MpIeee( "0.3" ), MpIeee( "20.0" ));
}

MpIeee test_pascal(void)
{
  return gsl_ran_pascal (r_global, MpIeee( "0.8" ), MpIeee( "3" ));
}

MpIeee test_pascal_pdf(unsigned int  n)
{
  return gsl_ran_pascal_pdf (n, MpIeee( "0.8" ), MpIeee( "3" ));
}


MpIeee test_pareto(void)
{
  return gsl_ran_pareto (r_global, MpIeee( "1.9" ), MpIeee( "2.75" ));
}

MpIeee test_pareto_pdf(MpIeee x)
{
  return gsl_ran_pareto_pdf (x, MpIeee( "1.9" ), MpIeee( "2.75" ));
}

MpIeee test_rayleigh(void)
{
  return gsl_ran_rayleigh (r_global, MpIeee( "1.9" ));
}

MpIeee test_rayleigh_pdf(MpIeee x)
{
  return gsl_ran_rayleigh_pdf (x, MpIeee( "1.9" ));
}

MpIeee test_rayleigh_tail(void)
{
  return gsl_ran_rayleigh_tail (r_global, MpIeee( "2.7" ), MpIeee( "1.9" ));
}

MpIeee test_rayleigh_tail_pdf(MpIeee x)
{
  return gsl_ran_rayleigh_tail_pdf (x, MpIeee( "2.7" ), MpIeee( "1.9" ));
}


MpIeee test_poisson(void)
{
  return gsl_ran_poisson (r_global, MpIeee( "5.0" ));
}

MpIeee test_poisson_pdf(unsigned int  n)
{
  return gsl_ran_poisson_pdf (n, MpIeee( "5.0" ));
}

MpIeee test_poisson_large(void)
{
  return gsl_ran_poisson (r_global, MpIeee( "30.0" ));
}

MpIeee test_poisson_large_pdf(unsigned int  n)
{
  return gsl_ran_poisson_pdf (n, MpIeee( "30.0" ));
}


MpIeee test_tdist1(void)
{
  return gsl_ran_tdist (r_global, MpIeee( "1.75" ));
}

MpIeee test_tdist1_pdf(MpIeee x)
{
  return gsl_ran_tdist_pdf (x, MpIeee( "1.75" ));
}

MpIeee test_tdist2(void)
{
  return gsl_ran_tdist (r_global, MpIeee( "12.75" ));
}

MpIeee test_tdist2_pdf(MpIeee x)
{
  return gsl_ran_tdist_pdf (x, MpIeee( "12.75" ));
}


MpIeee test_laplace(void)
{
  return gsl_ran_laplace (r_global, MpIeee( "2.75" ));
}

MpIeee test_laplace_pdf(MpIeee x)
{
  return gsl_ran_laplace_pdf (x, MpIeee( "2.75" ));
}

MpIeee test_weibull(void)
{
  return gsl_ran_weibull (r_global, MpIeee( "3.14" ), MpIeee( "2.75" ));
}

MpIeee test_weibull_pdf(MpIeee x)
{
  return gsl_ran_weibull_pdf (x, MpIeee( "3.14" ), MpIeee( "2.75" ));
}


MpIeee test_weibull1(void)
{
  return gsl_ran_weibull (r_global, MpIeee( "2.97" ), MpIeee( "1.0" ));
}

MpIeee test_weibull1_pdf(MpIeee x)
{
  return gsl_ran_weibull_pdf (x, MpIeee( "2.97" ), MpIeee( "1.0" ));
}
