#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/gamma.c
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static MpIeee gamma_large(const gsl_rng * r, const MpIeee a);
static MpIeee gamma_frac(const gsl_rng * r, const MpIeee a);

/* The Gamma distribution of order a>0 is defined by:

   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a1 and a2 with the same scale parameter b, then
   X+Y has gamma distribution of order a1+a2.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

MpIeee gsl_ran_gamma(const gsl_rng * r, const MpIeee a, const MpIeee b)
{
  /* assume a > 0 */
  unsigned int  na=  floor (a).toUnsignedInt();

  if (a == na)
    {
      return b * gsl_ran_gamma_int (r, na);
    }
  else if (na == 0)
    {
      return b * gamma_frac (r, a);
    }
  else
    {
      return b * (gsl_ran_gamma_int (r, na) + gamma_frac (r, a - na)) ;
    }
}

MpIeee gsl_ran_gamma_int(const gsl_rng * r, const unsigned int  a)
{
  if (a < 12)
    {
      unsigned int  i;
      MpIeee prod=  MpIeee( "1" );

      for (i = 0; i < a; i++)
        {
          prod *= gsl_rng_uniform_pos (r);
        }

      /* Note: for 12 iterations we are safe against underflow, since
         the smallest positive random number is O(2^-32). This means
         the smallest possible product is 2^(-12*32) = 10^-116 which
         is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gamma_large (r, (MpIeee) a);
    }
}

static MpIeee gamma_large(const gsl_rng * r, const MpIeee a)
{
  /* Works only if a > 1, and is most efficient if a is large

     This algorithm, reported in Knuth, is attributed to Ahrens.  A
     faster one, we are told, can be found in: J. H. Ahrens and
     U. Dieter, Computing 12 (1974) 223-246.  */

  MpIeee sqa;MpIeee  x;MpIeee  y;MpIeee  v;
  sqa = sqrt (MpIeee( "2" ) * a - MpIeee( "1" ));
  do
    {
      do
        {
          y = tan (M_PI * gsl_rng_uniform (r));
          x = sqa * y + a - MpIeee( "1" );
        }
      while (x <= MpIeee( "0" ));
      v = gsl_rng_uniform (r);
    }
  while (v > (MpIeee( "1" ) + y * y) * exp ((a - MpIeee( "1" )) * log (x / (a - MpIeee( "1" ))) - sqa * y));

  return x;
}

static MpIeee gamma_frac(const gsl_rng * r, const MpIeee a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  MpIeee p;MpIeee  q;MpIeee  x;MpIeee  u;MpIeee  v;
  p = M_E / (a + M_E);
  do
    {
      u = gsl_rng_uniform (r);
      v = gsl_rng_uniform_pos (r);

      if (u < p)
        {
          x = exp ((MpIeee( "1" ) / a) * log (v));
          q = exp (-x);
        }
      else
        {
          x = MpIeee( "1" ) - log (v);
          q = exp ((a - MpIeee( "1" )) * log (x));
        }
    }
  while (gsl_rng_uniform (r) >= q);

  return x;
}

MpIeee gsl_ran_gamma_pdf(const MpIeee x, const MpIeee a, const MpIeee b)
{
  if (x < 0)
    {
      return MpIeee( "0" ) ;
    }
  else if (x == 0)
    {
      if (a == 1)
        return MpIeee( "1" )/b ;
      else
        return MpIeee( "0" ) ;
    }
  else if (a == 1)
    {
      return exp(-x/b)/b ;
    }
  else 
    {
      MpIeee p;
      MpIeee lngamma=  gsl_sf_lngamma (a);
      p = exp ((a - MpIeee( "1" )) * log (x/b) - x/b - lngamma)/b;
      return p;
    }
}
