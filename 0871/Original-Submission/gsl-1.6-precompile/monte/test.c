#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* monte/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
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
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#define CONSTANT
#define PRODUCT
#define GAUSSIAN
#define DBLGAUSSIAN
#define TSUDA

#define PLAIN
#define MISER
#define VEGAS

MpIeee xl[11]  =  { MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ) };
MpIeee xu[11]  =  { MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ), MpIeee( "1" ) };
MpIeee xu2[11] =  { MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ), MpIeee( "2" ) };
MpIeee xu3[2]  =  { GSL_DBL_MAX, GSL_DBL_MAX };

MpIeee fconst(MpIeee x[], size_t d, void *params);
MpIeee f0(MpIeee x[], size_t d, void *params);
MpIeee f1(MpIeee x[], size_t d, void *params);
MpIeee f2(MpIeee x[], size_t d, void *params);
MpIeee f3(MpIeee x[], size_t d, void *params);

void my_error_handler (const char *reason, const char *file,
                       int  line, int  err);

struct problem {
  gsl_monte_function * f;
  MpIeee * xl;
  MpIeee * xu;
  size_t dim;
  size_t calls;
  MpIeee expected_result;
  MpIeee expected_error;
  char * description;
} ;

gsl_monte_function 
make_function (MpIeee(*f)(MpIeee*, size_t, void *), size_t d, void * p);

gsl_monte_function 
make_function (MpIeee(*f)(MpIeee*, size_t, void *), size_t d, void * p)
{
  gsl_monte_function f_new;

  f_new.f = f;
  f_new.dim = d;
  f_new.params = p;

  return f_new;
}


void 
add (struct problem * problems, int  * n, 
     gsl_monte_function * f, MpIeee xl[], MpIeee xu[], size_t dim, size_t calls,
     MpIeee result, MpIeee err, char * description);

void 
add (struct problem * problems, int  * n, 
     gsl_monte_function * f, MpIeee xl[], MpIeee xu[], size_t dim, size_t calls,
     MpIeee result, MpIeee err, char * description)
{
  int  i=  *n;
  problems[i].f = f;
  problems[i].xl = xl;
  problems[i].xu = xu;
  problems[i].dim = dim;
  problems[i].calls = calls;
  problems[i].expected_result = result;
  problems[i].expected_error = err;
  problems[i].description = description;
  (*n)++;
}

#define TRIALS 10

int
main (void)
{
  MpIeee result[TRIALS];MpIeee  error[TRIALS];
  MpIeee a=  MpIeee( "0.1" );
  MpIeee c=  (MpIeee( "1.0" ) + sqrt (MpIeee( "10.0" ))) / MpIeee( "9.0" );

  gsl_monte_function Fc = make_function(&fconst, 0, 0);
  gsl_monte_function F0 = make_function(&f0, 0, &a);
  gsl_monte_function F1 = make_function(&f1, 0, &a);
  gsl_monte_function F2 = make_function(&f2, 0, &a);
  gsl_monte_function F3 = make_function(&f3, 0, &c);

  /* The relationship between the variance of the function itself, the
     error on the integral and the number of calls is,

     sigma = sqrt(variance/N)

     where the variance is the <(f - <f>)^2> where <.> denotes the
     volume average (integral over the integration region divided by
     the volume) */

  int  n=  0;
  struct problem * I;
  struct problem problems[256];

#ifdef CONSTANT
    /* variance(Fc) = 0 */

    add(problems,&n, &Fc, xl, xu,  1, 1000, 1.0, 0.0, "constant, 1d");
    add(problems,&n, &Fc, xl, xu,  2, 1000, 1.0, 0.0, "constant, 2d");
    add(problems,&n, &Fc, xl, xu,  3, 1000, 1.0, 0.0, "constant, 3d");
    add(problems,&n, &Fc, xl, xu,  4, 1000, 1.0, 0.0, "constant, 4d");
    add(problems,&n, &Fc, xl, xu,  5, 1000, 1.0, 0.0, "constant, 5d");
    add(problems,&n, &Fc, xl, xu,  6, 1000, 1.0, 0.0, "constant, 6d");
    add(problems,&n, &Fc, xl, xu,  7, 1000, 1.0, 0.0, "constant, 7d");
    add(problems,&n, &Fc, xl, xu,  8, 1000, 1.0, 0.0, "constant, 8d");
    add(problems,&n, &Fc, xl, xu,  9, 1000, 1.0, 0.0, "constant, 9d");
    add(problems,&n, &Fc, xl, xu, 10, 1000, 1.0, 0.0, "constant, 10d");
#endif

#ifdef PRODUCT
    /* variance(F0) = (4/3)^d - 1 */

    add(problems,&n, &F0, xl, xu,  1, 3333,   1.0, 0.01, "product, 1d" );
    add(problems,&n, &F0, xl, xu,  2, 7777,   1.0, 0.01, "product, 2d" );
    add(problems,&n, &F0, xl, xu,  3, 13703,  1.0, 0.01, "product, 3d" );
    add(problems,&n, &F0, xl, xu,  4, 21604,  1.0, 0.01, "product, 4d" );
    add(problems,&n, &F0, xl, xu,  5, 32139,  1.0, 0.01, "product, 5d" );
    add(problems,&n, &F0, xl, xu,  6, 46186,  1.0, 0.01, "product, 6d" );
    add(problems,&n, &F0, xl, xu,  7, 64915,  1.0, 0.01, "product, 7d" );
    add(problems,&n, &F0, xl, xu,  8, 89887,  1.0, 0.01, "product, 8d" );
    add(problems,&n, &F0, xl, xu,  9, 123182, 1.0, 0.01, "product, 9d" );
    add(problems,&n, &F0, xl, xu, 10, 167577, 1.0, 0.01, "product, 10d" );
#endif

#ifdef GAUSSIAN
    /* variance(F1) = (1/(a sqrt(2 pi)))^d - 1 */

    add(problems,&n, &F1, xl, xu,  1, 298,      1.0, 0.1, "gaussian, 1d" );
    add(problems,&n, &F1, xl, xu,  2, 1492,     1.0, 0.1, "gaussian, 2d" );
    add(problems,&n, &F1, xl, xu,  3, 6249,     1.0, 0.1, "gaussian, 3d" );
    add(problems,&n, &F1, xl, xu,  4, 25230,    1.0, 0.1, "gaussian, 4d" );
    add(problems,&n, &F1, xl, xu,  5, 100953,   1.0, 0.1, "gaussian, 5d" );
    add(problems,&n, &F1, xl, xu,  6, 44782,    1.0, 0.3, "gaussian, 6d" );
    add(problems,&n, &F1, xl, xu,  7, 178690,   1.0, 0.3, "gaussian, 7d" );
    add(problems,&n, &F1, xl, xu,  8, 712904,   1.0, 0.3, "gaussian, 8d" );
    add(problems,&n, &F1, xl, xu,  9, 2844109,  1.0, 0.3, "gaussian, 9d" );
    add(problems,&n, &F1, xl, xu, 10, 11346390, 1.0, 0.3, "gaussian, 10d" );
#endif

#ifdef DBLGAUSSIAN
    /* variance(F2) = 0.5 * (1/(a sqrt(2 pi)))^d - 1 */

    add(problems,&n, &F2, xl, xu,  1, 9947,  1.0, 0.01, "double gaussian, 1d" );
    add(problems,&n, &F2, xl, xu,  2, 695,   1.0, 0.1, "double gaussian, 2d" );
    add(problems,&n, &F2, xl, xu,  3, 3074,  1.0, 0.1, "double gaussian, 3d" );
    add(problems,&n, &F2, xl, xu,  4, 12565, 1.0, 0.1, "double gaussian, 4d" );
    add(problems,&n, &F2, xl, xu,  5, 50426, 1.0, 0.1, "double gaussian, 5d" );
    add(problems,&n, &F2, xl, xu,  6, 201472, 1.0, 0.1, "double gaussian, 6d" );
    add(problems,&n, &F2, xl, xu,  7, 804056, 1.0, 0.1, "double gaussian, 7d" );
    add(problems,&n, &F2, xl, xu,  8, 356446, 1.0, 0.3, "double gaussian, 8d" );
    add(problems,&n, &F2, xl, xu,  9, 1422049, 1.0, 0.3, "double gaussian, 9d" );
    add(problems,&n, &F2, xl, xu, 10, 5673189, 1.0, 0.3, "double gaussian, 10d" );
#endif

#ifdef TSUDA
    /* variance(F3) = ((c^2 + c + 1/3)/(c(c+1)))^d - 1 */

    add(problems,&n, &F3, xl, xu,  1, 4928,   1.0, 0.01, "tsuda function, 1d" );
    add(problems,&n, &F3, xl, xu,  2, 12285,  1.0, 0.01, "tsuda function, 2d" );
    add(problems,&n, &F3, xl, xu,  3, 23268,  1.0, 0.01, "tsuda function, 3d" );
    add(problems,&n, &F3, xl, xu,  4, 39664,  1.0, 0.01, "tsuda function, 4d" );
    add(problems,&n, &F3, xl, xu,  5, 64141,  1.0, 0.01, "tsuda function, 5d" );
    add(problems,&n, &F3, xl, xu,  6, 100680, 1.0, 0.01, "tsuda function, 6d" );
    add(problems,&n, &F3, xl, xu,  7, 155227, 1.0, 0.01, "tsuda function, 7d" );
    add(problems,&n, &F3, xl, xu,  8, 236657, 1.0, 0.01, "tsuda function, 8d" );
    add(problems,&n, &F3, xl, xu,  9, 358219, 1.0, 0.01, "tsuda function, 9d" );
    add(problems,&n, &F3, xl, xu, 10, 539690, 1.0, 0.01, "tsuda function, 10d" );
#endif

    add(problems,&n,   0,  0,  0, 0,    0,   0,   0, 0  );


  /* gsl_set_error_handler (&my_error_handler); */
  gsl_ieee_env_setup ();
  gsl_rng_env_setup ();

#ifdef A
  {cout<<"testing allocation/input checks\n";}

  status = gsl_monte_plain_validate (s, xl, xu3, 1, 1);
  gsl_test (status != 0, "error if limits too large");
  status = gsl_monte_plain_validate (s, xl, xu, 0, 10);
  gsl_test (status != 0, "error if num_dim = 0");
  status = gsl_monte_plain_validate (s, xl, xu, 1, 0);
  gsl_test (status != 0, "error if calls = 0");
  status = gsl_monte_plain_validate (s, xu, xl, 1, 10);
  gsl_test (status != 0, "error if xu < xl");
#endif

#ifdef PLAIN
#define NAME "plain"
#define MONTE_STATE gsl_monte_plain_state
#define MONTE_ALLOC gsl_monte_plain_alloc
#define MONTE_INTEGRATE gsl_monte_plain_integrate
#define MONTE_FREE gsl_monte_plain_free
#define MONTE_SPEEDUP 1
#define MONTE_ERROR_TEST(err,expected) gsl_test_factor(err,expected, 5.0, NAME ", %s, abserr[%d]", I->description, i)
#include "test_main.c"
#undef NAME
#undef MONTE_STATE
#undef MONTE_ALLOC
#undef MONTE_INTEGRATE
#undef MONTE_FREE
#undef MONTE_ERROR_TEST
#undef MONTE_SPEEDUP
#endif

#ifdef MISER
#define NAME "miser"
#define MONTE_STATE gsl_monte_miser_state
#define MONTE_ALLOC gsl_monte_miser_alloc
#define MONTE_INTEGRATE gsl_monte_miser_integrate
#define MONTE_FREE gsl_monte_miser_free
#define MONTE_SPEEDUP 2
#define MONTE_ERROR_TEST(err,expected) gsl_test(err > 5.0 * expected, NAME ", %s, abserr[%d] (obs %g vs plain %g)", I->description, i, err, expected)
#include "test_main.c"
#undef NAME
#undef MONTE_STATE
#undef MONTE_ALLOC
#undef MONTE_INTEGRATE
#undef MONTE_FREE
#undef MONTE_ERROR_TEST
#undef MONTE_SPEEDUP
#endif

#ifdef VEGAS
#define NAME "vegas"
#define MONTE_STATE gsl_monte_vegas_state
#define MONTE_ALLOC gsl_monte_vegas_alloc
#define MONTE_INTEGRATE(f,xl,xu,dim,calls,r,s,res,err) { gsl_monte_vegas_integrate(f,xl,xu,dim,calls,r,s,res,err) ; if (s->chisq < 0.5 || s->chisq > 2) gsl_monte_vegas_integrate(f,xl,xu,dim,calls,r,s,res,err); }
#define MONTE_FREE gsl_monte_vegas_free
#define MONTE_SPEEDUP 3
#define MONTE_ERROR_TEST(err,expected) gsl_test(err > 3.0 * (expected == 0 ? 1.0/(I->calls/MONTE_SPEEDUP) : expected), NAME ", %s, abserr[%d] (obs %g vs exp %g)", I->description, i, err, expected)
#include "test_main.c"
#undef NAME
#undef MONTE_STATE
#undef MONTE_ALLOC
#undef MONTE_INTEGRATE
#undef MONTE_FREE
#undef MONTE_ERROR_TEST
#undef MONTE_SPEEDUP
#endif
      
  exit (gsl_test_summary ());
}

/* Simple constant function */
MpIeee fconst(MpIeee x[], size_t num_dim, void *params)
{
  return MpIeee( "1" );
}

/* Simple product function */
MpIeee f0(MpIeee x[], size_t num_dim, void *params)
{
  MpIeee prod=  MpIeee( "1.0" );
  unsigned int  i;

  for (i = 0; i < num_dim; ++i)
    {
      prod *= MpIeee( "2.0" ) * x[i];
    }

  return prod;
}

/* Gaussian centered at 1/2. */

MpIeee f1(MpIeee x[], size_t num_dim, void *params)
{
  MpIeee a=  *(MpIeee *)params;
  MpIeee sum=  MpIeee( "0." );

  unsigned int  i;
  for (i = 0; i < num_dim; i++)
    {
      MpIeee dx=  x[i] - MpIeee( "0.5" );
      sum += dx * dx;
    }
  return (pow (M_2_SQRTPI / (MpIeee( "2." ) * a), (MpIeee) num_dim) *
          exp (-sum / (a * a)));
}

/* double gaussian */
MpIeee f2(MpIeee x[], size_t num_dim, void *params)
{
  MpIeee a=  *(MpIeee *)params;
  MpIeee sum1=  MpIeee( "0." );
  MpIeee sum2=  MpIeee( "0." );

  unsigned int  i;
  for (i = 0; i < num_dim; i++)
    {
      MpIeee dx1=  x[i] - MpIeee( "1." ) / MpIeee( "3." );
      MpIeee dx2=  x[i] - MpIeee( "2." ) / MpIeee( "3." );
      sum1 += dx1 * dx1;
      sum2 += dx2 * dx2;
    }
  return MpIeee( "0.5" ) * pow (M_2_SQRTPI / (MpIeee( "2." ) * a), num_dim) 
    * (exp (-sum1 / (a * a)) + exp (-sum2 / (a * a)));
}

/* Tsuda's example */
MpIeee f3(MpIeee x[], size_t num_dim, void *params)
{
  MpIeee c=  *(MpIeee *)params;

  MpIeee prod=  MpIeee( "1." );

  unsigned int  i;

  for (i = 0; i < num_dim; i++)
    {
      prod *= c / (c + MpIeee( "1" )) * pow((c + MpIeee( "1" )) / (c + x[i]), MpIeee( "2.0" ));
    }

  return prod;
}


void
my_error_handler (const char *reason, const char *file, int  line, int  err)
{
  if (0)
    {cout<<"(caught ["<< file<<":"<< line<<": "<< reason<<" ("<< err<<")])\n";}
}
