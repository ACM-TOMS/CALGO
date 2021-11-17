#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/beta.c
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

/* The beta distribution has the form

   p(x) dx = (Gamma(a + b)/(Gamma(a) Gamma(b))) x^(a-1) (1-x)^(b-1) dx

   The method used here is the one described in Knuth */

MpIeee gsl_ran_beta(const gsl_rng * r, const MpIeee a, const MpIeee b)
{
  MpIeee x1=  gsl_ran_gamma (r, a, MpIeee( "1.0" ));
  MpIeee x2=  gsl_ran_gamma (r, b, MpIeee( "1.0" ));

  return x1 / (x1 + x2);
}

MpIeee gsl_ran_beta_pdf(const MpIeee x, const MpIeee a, const MpIeee b)
{
  if (x < 0 || x > 1)
    {
      return MpIeee( "0" ) ;
    }
  else 
    {
      MpIeee p;

      MpIeee gab=  gsl_sf_lngamma (a + b);
      MpIeee ga=  gsl_sf_lngamma (a);
      MpIeee gb=  gsl_sf_lngamma (b);

      p = exp (gab - ga - gb) * pow (x, a - MpIeee( "1" )) * pow (MpIeee( "1" ) - x, b - MpIeee( "1" ));

      return p;
    }
}
