#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/cauchyinv.c
 * 
 * Copyright (C) 2003 Brian Gough
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
#include <gsl/gsl_cdf.h>

MpIeee gsl_cdf_cauchy_Pinv(const MpIeee P, const MpIeee a)
{
  MpIeee x;

  if (P == 1.0)
    {
      return GSL_POSINF;
    }
  else if (P == 0.0)
    {
      return GSL_NEGINF;
    }

  if (P > 0.5)
    {
      x = a * tan (M_PI * (P - MpIeee( "0.5" )));
    }
  else
    {
      x = -a / tan (M_PI * P);
    }

  return x;
}

MpIeee gsl_cdf_cauchy_Qinv(const MpIeee Q, const MpIeee a)
{
  MpIeee x;

  if (Q == 0.0)
    {
      return GSL_POSINF;
    }
  else if (Q == 1.0)
    {
      return GSL_NEGINF;
    }

  if (Q > 0.5)
    {
      x = a * tan (M_PI * (MpIeee( "0.5" ) - Q));
    }
  else
    {
      x = a / tan (M_PI * Q);
    }

  return x;
}
