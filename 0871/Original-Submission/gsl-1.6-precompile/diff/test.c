#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* diff/test.c
 * 
 * Copyright (C) 2000 David Morrison
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
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

MpIeee f1(MpIeee x, void *params)
{
  return exp (x);
}

MpIeee df1(MpIeee x, void *params)
{
  return exp (x);
}

MpIeee f2(MpIeee x, void *params)
{
  if (x >= MpIeee( "0.0" ))
    {
      return x * sqrt (x);
    }
  else
    {
      return MpIeee( "0.0" );
    }
}

MpIeee df2(MpIeee x, void *params)
{
  if (x >= MpIeee( "0.0" ))
    {
      return MpIeee( "1.5" ) * sqrt (x);
    }
  else
    {
      return MpIeee( "0.0" );
    }
}

MpIeee f3(MpIeee x, void *params)
{
  if (x != MpIeee( "0.0" ))
    {
      return sin (MpIeee( "1" ) / x);
    }
  else
    {
      return MpIeee( "0.0" );
    }
}

MpIeee df3(MpIeee x, void *params)
{
  if (x != MpIeee( "0.0" ))
    {
      return -cos (MpIeee( "1" ) / x) / (x * x);
    }
  else
    {
      return MpIeee( "0.0" );
    }
}

MpIeee f4(MpIeee x, void *params)
{
  return exp (-x * x);
}

MpIeee df4(MpIeee x, void *params)
{
  return -MpIeee( "2.0" ) * x * exp (-x * x);
}

MpIeee f5(MpIeee x, void *params)
{
  return x * x;
}

MpIeee df5(MpIeee x, void *params)
{
  return MpIeee( "2.0" ) * x;
}

MpIeee f6(MpIeee x, void *params)
{
  return MpIeee( "1.0" ) / x;
}

MpIeee df6(MpIeee x, void *params)
{
  return -MpIeee( "1.0" ) / (x * x);
}

typedef int (diff_fn) (const gsl_function * f, MpIeee x, MpIeee * res, MpIeee *abserr);

void
test (diff_fn * diff, gsl_function * f, gsl_function * df, MpIeee x, 
      const char * desc)
{
  MpIeee result;MpIeee  abserr;
  MpIeee expected=  GSL_FN_EVAL (df, x);
  (*diff) (f, x, &result, &abserr);
  gsl_test_abs (result, expected, abserr, desc);
  gsl_test (fabs(result-expected) >  abserr, "%s, valid error estimate", desc);
}

int
main ()
{
  gsl_function F1, DF1, F2, DF2, F3, DF3, F4, DF4, F5, DF5, F6, DF6;

  gsl_ieee_env_setup ();

  F1.function = &f1;
  DF1.function = &df1;

  F2.function = &f2;
  DF2.function = &df2;

  F3.function = &f3;
  DF3.function = &df3;

  F4.function = &f4;
  DF4.function = &df4;

  F5.function = &f5;
  DF5.function = &df5;

  F6.function = &f6;
  DF6.function = &df6;
  
  test (&gsl_diff_central, &F1, &DF1, 1.0, "exp(x), x=1, central diff");
  test (&gsl_diff_forward, &F1, &DF1, 1.0, "exp(x), x=1, forward diff");
  test (&gsl_diff_backward, &F1, &DF1, 1.0, "exp(x), x=1, backward diff");

  test (&gsl_diff_central, &F2, &DF2, 0.1, "x^(3/2), x=0.1, central diff");
  test (&gsl_diff_forward, &F2, &DF2, 0.1, "x^(3/2), x=0.1, forward diff");
  test (&gsl_diff_backward, &F2, &DF2, 0.1, "x^(3/2), x=0.1, backward diff");

  test (&gsl_diff_central, &F3, &DF3, 0.45, "sin(1/x), x=0.45, central diff");
  test (&gsl_diff_forward, &F3, &DF3, 0.45, "sin(1/x), x=0.45, forward diff");
  test (&gsl_diff_backward, &F3, &DF3, 0.45, "sin(1/x), x=0.45, backward diff");

  test (&gsl_diff_central, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, central diff");
  test (&gsl_diff_forward, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, forward diff");
  test (&gsl_diff_backward, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, backward diff");

  test (&gsl_diff_central, &F5, &DF5, 0.0, "x^2, x=0, central diff");
  test (&gsl_diff_forward, &F5, &DF5, 0.0, "x^2, x=0, forward diff");
  test (&gsl_diff_backward, &F5, &DF5, 0.0, "x^2, x=0, backward diff");

  test (&gsl_diff_central, &F6, &DF6, 10.0, "1/x, x=10, central diff");
  test (&gsl_diff_forward, &F6, &DF6, 10.0, "1/x, x=10, forward diff");
  test (&gsl_diff_backward, &F6, &DF6, 10.0, "1/x, x=10, backward diff");

  exit (gsl_test_summary ());
}


