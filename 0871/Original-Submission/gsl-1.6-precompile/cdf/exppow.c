#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/exppow.c
 * 
 * Copyright (C) 2004 Giulio Bottazzi
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>

/* The exponential power density is parametrized according to

   p(x) dx = (1/(2 a Gamma(1 + 1/b))) * exp(-|x/a|^b) dx

   so that the distribution reads

           / x<0   0.5 - Gamma_inc_P(1/b,|x/a|^b)
   P(x) = |  x=0   0.5
           \ x>0   0.5 + Gamma_inc_P(1/b,|x/a|^b)


   for x in (-infty,+infty) */

MpIeee gsl_cdf_exppow_P(const MpIeee x, const MpIeee a, const MpIeee b)
{
  const MpIeee u=  x / a;

  if (u < 0)
    {
      MpIeee P=  MpIeee( "0.5" ) * gsl_sf_gamma_inc_Q (MpIeee( "1.0" ) / b, pow (-u, b));
      return P;
    }
  else
    {
      MpIeee P=  MpIeee( "0.5" ) * (MpIeee( "1.0" ) + gsl_sf_gamma_inc_P (MpIeee( "1.0" ) / b, pow (u, b)));
      return P;
    }
}

MpIeee gsl_cdf_exppow_Q(const MpIeee x, const MpIeee a, const MpIeee b)
{
  const MpIeee u=  x / a;

  if (u < 0)
    {
      MpIeee Q=  MpIeee( "0.5" ) * (MpIeee( "1.0" ) + gsl_sf_gamma_inc_P (MpIeee( "1.0" ) / b, pow (-u, b)));
      return Q;
    }
  else
    {
      MpIeee Q=  MpIeee( "0.5" ) * gsl_sf_gamma_inc_Q (MpIeee( "1.0" ) / b, pow (u, b));
      return Q;
    }
}
