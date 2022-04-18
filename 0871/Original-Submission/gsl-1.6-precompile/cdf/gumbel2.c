#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/gumbel2.c
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

MpIeee gsl_cdf_gumbel2_P(const MpIeee x, const MpIeee a, const MpIeee b)
{
  MpIeee P;

  if (x == 0)
    {
      P = MpIeee( "0" );
    }
  else
    {
      MpIeee u=  pow (x, a);
      P = exp (-b / u);
    }

  return P;
}

MpIeee gsl_cdf_gumbel2_Q(const MpIeee x, const MpIeee a, const MpIeee b)
{
  MpIeee Q;

  if (x == 0)
    {
      Q = MpIeee( "1" );
    }
  else
    {
      MpIeee u=  pow (x, a);
      Q = -expm1 (-b / u);
    }

  return Q;
}
