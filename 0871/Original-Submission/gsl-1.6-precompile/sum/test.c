#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sum/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sum.h>

#include <gsl/gsl_ieee_utils.h>

#define N 50

void check_trunc (MpIeee * t, MpIeee expected, const char * desc);
void check_full (MpIeee * t, MpIeee expected, const char * desc);

int
main (void)
{
  gsl_ieee_env_setup ();

  {
    MpIeee t[N];
    int  n;

    const MpIeee zeta_2=  M_PI * M_PI / 6.0;

    /* terms for zeta(2) */

    for (n = 0; n < N; n++)
      {
        MpIeee np1=  n + MpIeee( "1.0" );
        t[n] = MpIeee( "1.0" ) / (np1 * np1);
      }

    check_trunc (t, zeta_2, "zeta(2)");
    check_full (t, zeta_2, "zeta(2)");
  }

  {
    MpIeee t[N];
    MpIeee x;MpIeee  y;
    int  n;

    /* terms for exp(10.0) */
    x = MpIeee( "10.0" );
    y = exp(x);

    t[0] = MpIeee( "1.0" );
    for (n = 1; n < N; n++)
      {
        t[n] = t[n - 1] * (x / n);
      }

    check_trunc (t, y, "exp(10)");
    check_full (t, y, "exp(10)");
  }

  {
    MpIeee t[N];
    MpIeee x;MpIeee  y;
    int  n;

    /* terms for exp(-10.0) */
    x = -MpIeee( "10.0" );
    y = exp(x);

    t[0] = MpIeee( "1.0" );
    for (n = 1; n < N; n++)
      {
        t[n] = t[n - 1] * (x / n);
      }

    check_trunc (t, y, "exp(-10)");
    check_full (t, y, "exp(-10)");
  }

  {
    MpIeee t[N];
    MpIeee x;MpIeee  y;
    int  n;

    /* terms for -log(1-x) */
    x = MpIeee( "0.5" );
    y = -log(MpIeee( "1" )-x);
    t[0] = x;
    for (n = 1; n < N; n++)
      {
        t[n] = t[n - 1] * (x * n) / (n + MpIeee( "1.0" ));
      }

    check_trunc (t, y, "-log(1/2)");
    check_full (t, y, "-log(1/2)");
  }

  {
    MpIeee t[N];
    MpIeee x;MpIeee  y;
    int  n;

    /* terms for -log(1-x) */
    x = -MpIeee( "1.0" );
    y = -log(MpIeee( "1" )-x);
    t[0] = x;
    for (n = 1; n < N; n++)
      {
        t[n] = t[n - 1] * (x * n) / (n + MpIeee( "1.0" ));
      }

    check_trunc (t, y, "-log(2)");
    check_full (t, y, "-log(2)");
  }

  {
    MpIeee t[N];
    int  n;

    MpIeee result=  MpIeee( "0.192594048773" );

    /* terms for an alternating asymptotic series */

    t[0] = MpIeee( "3.0" ) / (M_PI * M_PI);

    for (n = 1; n < N; n++)
      {
        t[n] = -t[n - 1] * (MpIeee( "4.0" ) * (n + MpIeee( "1.0" )) - MpIeee( "1.0" )) / (M_PI * M_PI);
      }

    check_trunc (t, result, "asymptotic series");
    check_full (t, result, "asymptotic series");
  }

  {
    MpIeee t[N];
    int  n;

    /* Euler's gamma from GNU Calc (precision = 32) */

    MpIeee result=  MpIeee( "0.5772156649015328606065120900824" ); 

    /* terms for Euler's gamma */

    t[0] = MpIeee( "1.0" );

    for (n = 1; n < N; n++)
      {
        t[n] = MpIeee( "1" )/(n+MpIeee( "1.0" )) + log(n/(n+MpIeee( "1.0" )));
      }

    check_trunc (t, result, "Euler's constant");
    check_full (t, result, "Euler's constant");
  }

  {
    MpIeee t[N];
    int  n;

    /* eta(1/2) = sum_{k=1}^{\infty} (-1)^(k+1) / sqrt(k)

       From Levin, Intern. J. Computer Math. B3:371--388, 1973.

       I=(1-sqrt(2))zeta(1/2)
        =(2/sqrt(pi))*integ(1/(exp(x^2)+1),x,0,inf) */

    MpIeee result=  MpIeee( "0.6048986434216305" );  /* approx */

    /* terms for eta(1/2) */

    for (n = 0; n < N; n++)
      {
        t[n] = (n%MpIeee( "2" ) ? -MpIeee( "1" ) : MpIeee( "1" )) * MpIeee( "1.0" ) /sqrt(n + MpIeee( "1.0" ));
      }

    check_trunc (t, result, "eta(1/2)");
    check_full (t, result, "eta(1/2)");
  }

  exit (gsl_test_summary ());
}

void
check_trunc (MpIeee * t, MpIeee expected, const char * desc)
{
  MpIeee sum_accel;MpIeee  prec;

  gsl_sum_levin_utrunc_workspace * w = gsl_sum_levin_utrunc_alloc (N);
  
  gsl_sum_levin_utrunc_accel (t, N, w, &sum_accel, &prec);
  gsl_test_rel (sum_accel, expected, 1e-8, "trunc result, %s", desc);

  /* No need to check precision for truncated result since this is not
     a meaningful number */

  gsl_sum_levin_utrunc_free (w);
}

void
check_full (MpIeee * t, MpIeee expected, const char * desc)
{
  MpIeee sum_accel;MpIeee  err_est;MpIeee  sd_actual;MpIeee  sd_est;
  
  gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc (N);

  gsl_sum_levin_u_accel (t, N, w, &sum_accel, &err_est);
  gsl_test_rel (sum_accel, expected, 1e-8, "full result, %s", desc);
  
  sd_est = -log10 (err_est/fabs(sum_accel));
  sd_actual = -log10 (DBL_EPSILON + fabs ((sum_accel - expected)/expected));

  /* Allow one digit of slop */

  gsl_test (sd_est > sd_actual + MpIeee( "1.0" ), "full significant digits, %s (%g vs %g)", desc, sd_est, sd_actual);

  gsl_sum_levin_u_free (w);
}
