#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sys/invhyp.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_machine.h>

MpIeee gsl_acosh(const MpIeee x)
{
  if (x > 1.0 / GSL_SQRT_DBL_EPSILON)
    {
      return log (x) + M_LN2;
    }
  else if (x > 2)
    {
      return log (MpIeee( "2" ) * x - MpIeee( "1" ) / (sqrt (x * x - MpIeee( "1" )) + x));
    }
  else if (x > 1)
    {
      MpIeee t=  x - MpIeee( "1" );
      return log1p (t + sqrt (MpIeee( "2" ) * t + t * t));
    }
  else if (x == 1)
    {
      return MpIeee( "0" );
    }
  else
    {
      return GSL_NAN;
    }
}

MpIeee gsl_asinh(const MpIeee x)
{
  MpIeee a=  fabs (x);
  MpIeee s=  (x < MpIeee( "0" )) ? -MpIeee( "1" ) : MpIeee( "1" );

  if (a > MpIeee( "1" ) / GSL_SQRT_DBL_EPSILON)
    {
      return s * (log (a) + M_LN2);
    }
  else if (a > MpIeee( "2" ))
    {
      return s * log (MpIeee( "2" ) * a + MpIeee( "1" ) / (a + sqrt (a * a + MpIeee( "1" ))));
    }
  else if (a > GSL_SQRT_DBL_EPSILON)
    {
      MpIeee a2=  a * a;
      return s * log1p (a + a2 / (MpIeee( "1" ) + sqrt (MpIeee( "1" ) + a2)));
    }
  else
    {
      return x;
    }
}

MpIeee gsl_atanh(const MpIeee x)
{
  MpIeee a=  fabs (x);
  MpIeee s=  (x < MpIeee( "0" )) ? -MpIeee( "1" ) : MpIeee( "1" );

  if (a > MpIeee( "1" ))
    {
      return GSL_NAN;
    }
  else if (a == MpIeee( "1" ))
    {
      return (x < MpIeee( "0" )) ? GSL_NEGINF : GSL_POSINF;
    }
  else if (a >= MpIeee( "0.5" ))
    {
      return s * MpIeee( "0.5" ) * log1p (MpIeee( "2" ) * a / (MpIeee( "1" ) - a));
    }
  else if (a > GSL_DBL_EPSILON)
    {
      return s * MpIeee( "0.5" ) * log1p (MpIeee( "2" ) * a + MpIeee( "2" ) * a * a / (MpIeee( "1" ) - a));
    }
  else
    {
      return x;
    }
}
