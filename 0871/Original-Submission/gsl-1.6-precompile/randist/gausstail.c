#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/gausstail.c
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
#include <gsl/gsl_sf_erf.h>

MpIeee gsl_ran_gaussian_tail(const gsl_rng * r, const MpIeee a, const MpIeee sigma)
{
  /* Returns a gaussian random variable larger than a
   * This implementation does one-sided upper-tailed deviates.
   */

  MpIeee s=  a / sigma;

  if (s < MpIeee( "1" ))
    {
      /* For small s, use a direct rejection method. The limit s < 1
         can be adjusted to optimise the overall efficiency */

      MpIeee x;

      do
        {
          x = gsl_ran_gaussian (r, MpIeee( "1.0" ));
        }
      while (x < s);
      return x * sigma;
    }
  else
    {
      /* Use the "supertail" deviates from the last two steps
       * of Marsaglia's rectangle-wedge-tail method, as described
       * in Knuth, v2, 3rd ed, pp 123-128.  (See also exercise 11, p139,
       * and the solution, p586.)
       */

      MpIeee u;MpIeee  v;MpIeee  x;

      do
        {
          u = gsl_rng_uniform (r);
          do
            {
              v = gsl_rng_uniform (r);
            }
          while (v == MpIeee( "0.0" ));
          x = sqrt (s * s - MpIeee( "2" ) * log (v));
        }
      while (x * u > s);
      return x * sigma;
    }
}

MpIeee gsl_ran_gaussian_tail_pdf(const MpIeee x, const MpIeee a, const MpIeee sigma)
{
  if (x < a)
    {
      return MpIeee( "0" );
    }
  else
    {
      MpIeee N;MpIeee  p;
      MpIeee u=  x / sigma ;

      MpIeee f=  gsl_sf_erfc (a / (sqrt (MpIeee( "2.0" )) * sigma));

      N = MpIeee( "0.5" ) * f;

      p = (MpIeee( "1" ) / (N * sqrt (MpIeee( "2" ) * M_PI) * sigma)) * exp (-u * u / MpIeee( "2" ));

      return p;
    }
}

MpIeee gsl_ran_ugaussian_tail(const gsl_rng * r, const MpIeee a)
{
  return gsl_ran_gaussian_tail (r, a, MpIeee( "1.0" )) ;
}

MpIeee gsl_ran_ugaussian_tail_pdf(const MpIeee x, const MpIeee a)
{
  return gsl_ran_gaussian_tail_pdf (x, a, MpIeee( "1.0" )) ;
}
