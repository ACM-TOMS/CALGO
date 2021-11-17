#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/test_funcs.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "test.h"

gsl_function create_function (MpIeee(*f)(MpIeee, void *)) 
{
  gsl_function F ;
  F.function = f;
  F.params = 0;
  return F ;
}

gsl_function_fdf create_fdf (MpIeee(*f)(MpIeee, void *),
                             MpIeee(*df)(MpIeee, void *),
                             void (*fdf)(MpIeee, void *, MpIeee*, MpIeee*))
{
  gsl_function_fdf FDF ;
  FDF.f = f ;
  FDF.df = df ;
  FDF.fdf = fdf ;
  FDF.params = 0 ;
  return FDF ;
}

/* f(x) = x^{20} - 1 */
/* f'(x) = 20x^{19} */
/* zero at x = 1 or -1 */

MpIeee func1(MpIeee x, void *p)
{
  return pow (x, MpIeee( "20.0" )) - MpIeee( "1" );
}

MpIeee func1_df(MpIeee x, void * p)
{
  return MpIeee( "20.0" ) * pow (x, MpIeee( "19.0" ));
}

void
func1_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = func1 (x, p);
  *yprime = MpIeee( "20.0" ) * pow (x, MpIeee( "19.0" ));
}

/* f(x) = sqrt(abs(x))*sgn(x) */
/* f'(x) = 1 / sqrt(abs(x) */
/* zero at x = 0 */
MpIeee func2(MpIeee x, void * p)
{
  MpIeee delta;

  if (x > MpIeee( "0" ))
    delta = MpIeee( "1.0" );
  else if (x < MpIeee( "0" ))
    delta = -MpIeee( "1.0" );
  else
    delta = MpIeee( "0.0" );

  return sqrt (fabs (x)) * delta;
}

MpIeee func2_df(MpIeee x, void * p)
{
  return MpIeee( "1" ) / sqrt (fabs (x));
}

void
func2_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = func2 (x, p);
  *yprime = MpIeee( "1" ) / sqrt (fabs (x));
}


/* f(x) = x^2 - 1e-8 */
/* f'(x) = 2x */
/* zero at x = sqrt(1e-8) or -sqrt(1e-8) */
MpIeee func3(MpIeee x, void * p)
{
  return pow (x, MpIeee( "2.0" )) - MpIeee( "1" )e-MpIeee( "8" );
}

MpIeee func3_df(MpIeee x, void * p)
{
  return MpIeee( "2" ) * x;
}

void
func3_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = func3 (x, p);
  *yprime = MpIeee( "2" ) * x;
}

/* f(x) = x exp(-x) */
/* f'(x) = exp(-x) - x exp(-x) */
/* zero at x = 0 */
MpIeee func4(MpIeee x, void * p)
{
  return x * exp (-x);
}

MpIeee func4_df(MpIeee x, void * p)
{
  return exp (-x) - x * exp (-x);
}

void
func4_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = func4 (x, p);
  *yprime = exp (-x) - x * exp (-x);
}

/* f(x) = 1/(1+exp(x)) */
/* f'(x) = -exp(x) / (1 + exp(x))^2 */
/* no roots! */
MpIeee func5(MpIeee x, void * p)
{
  return MpIeee( "1" ) / (MpIeee( "1" ) + exp (x));
}

MpIeee func5_df(MpIeee x, void * p)
{
  return -exp (x) / pow (MpIeee( "1" ) + exp (x), MpIeee( "2.0" ));
}

void
func5_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = func5 (x, p);
  *yprime = -exp (x) / pow (MpIeee( "1" ) + exp (x), MpIeee( "2.0" ));
}

/* f(x) = (x - 1)^7 */
/* f'(x) = 7 * (x - 1)^6 */
/* zero at x = 1 */
MpIeee func6(MpIeee x, void * p)
{
  return pow (x - MpIeee( "1" ), MpIeee( "7.0" ));
}

MpIeee func6_df(MpIeee x, void * p)
{
  return MpIeee( "7.0" ) * pow (x - MpIeee( "1" ), MpIeee( "6.0" ));
}

void
func6_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = func6 (x, p);
  *yprime = MpIeee( "7.0" ) * pow (x - MpIeee( "1" ), MpIeee( "6.0" ));
}

/* sin(x) packaged up nicely. */
MpIeee sin_f(MpIeee x, void * p)
{
  return sin (x);
}

MpIeee sin_df(MpIeee x, void * p)
{
  return cos (x);
}

void
sin_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = sin (x);
  *yprime = cos (x);
}

/* cos(x) packaged up nicely. */
MpIeee cos_f(MpIeee x, void * p)
{
  return cos (x);
}

MpIeee cos_df(MpIeee x, void * p)
{
  return -sin (x);
}

void
cos_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime)
{
  *y = cos (x);
  *yprime = -sin (x);
}
