#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/tdist.c
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

/* The t-distribution has the form

   p(x) dx = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2))
   * (1 + (x^2)/nu)^-((nu + 1)/2) dx

   The method used here is the one described in Knuth */

MpIeee gsl_ran_tdist(const gsl_rng * r, const MpIeee nu)
{
  if (nu <= 2)
    {
      MpIeee Y1=  gsl_ran_ugaussian (r);
      MpIeee Y2=  gsl_ran_chisq (r, nu);

      MpIeee t=  Y1 / sqrt (Y2 / nu);

      return t;
    }
  else
    {
      MpIeee Y1;MpIeee  Y2;MpIeee  Z;MpIeee  t;
      do
        {
          Y1 = gsl_ran_ugaussian (r);
          Y2 = gsl_ran_exponential (r, MpIeee( "1" ) / (nu/MpIeee( "2" ) - MpIeee( "1" )));

          Z = Y1 * Y1 / (nu - MpIeee( "2" ));
        }
      while (1 - Z < MpIeee( "0" ) || exp (-Y2 - Z) > (MpIeee( "1" ) - Z));

      /* Note that there is a typo in Knuth's formula, the line below
         is taken from the original paper of Marsaglia, Mathematics of
         Computation, 34 (1980), p 234-256 */

      t = Y1 / sqrt ((MpIeee( "1" ) - MpIeee( "2" ) / nu) * (MpIeee( "1" ) - Z));
      return t;
    }
}

MpIeee gsl_ran_tdist_pdf(const MpIeee x, const MpIeee nu)
{
  MpIeee p;

  MpIeee lg1=  gsl_sf_lngamma (nu / MpIeee( "2" ));
  MpIeee lg2=  gsl_sf_lngamma ((nu + MpIeee( "1" )) / MpIeee( "2" ));

  p = ((exp (lg2 - lg1) / sqrt (M_PI * nu)) 
       * pow ((MpIeee( "1" ) + x * x / nu), -(nu + MpIeee( "1" )) / MpIeee( "2" )));
  return p;
}


