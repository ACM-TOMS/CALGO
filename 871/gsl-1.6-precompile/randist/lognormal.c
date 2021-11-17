#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/lognormal.c
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The lognormal distribution has the form 

   p(x) dx = 1/(x * sqrt(2 pi sigma^2)) exp(-(ln(x) - zeta)^2/2 sigma^2) dx

   for x > 0. Lognormal random numbers are the exponentials of
   gaussian random numbers */

MpIeee gsl_ran_lognormal(const gsl_rng * r, const MpIeee zeta, const MpIeee sigma)
{
  MpIeee u;MpIeee  v;MpIeee  r2;MpIeee  normal;MpIeee  z;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -MpIeee( "1" ) + MpIeee( "2" ) * gsl_rng_uniform (r);
      v = -MpIeee( "1" ) + MpIeee( "2" ) * gsl_rng_uniform (r);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > MpIeee( "1.0" ) || r2 == MpIeee( "0" ));

  normal = u * sqrt (-MpIeee( "2.0" ) * log (r2) / r2);

  z =  exp (sigma * normal + zeta);

  return z;
}

MpIeee gsl_ran_lognormal_pdf(const MpIeee x, const MpIeee zeta, const MpIeee sigma)
{
  if (x <= 0)
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee u=  (log (x) - zeta)/sigma;
      MpIeee p=  MpIeee( "1" ) / (x * fabs(sigma) * sqrt (MpIeee( "2" ) * M_PI)) * exp (-(u * u) /MpIeee( "2" ));
      return p;
    }
}
