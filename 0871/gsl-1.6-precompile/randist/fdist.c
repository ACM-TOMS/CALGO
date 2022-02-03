#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/fdist.c
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The F distribution has the form

   p(x) dx = (nu1^(nu1/2) nu2^(nu2/2) Gamma((nu1 + nu2)/2) /
   Gamma(nu1/2) Gamma(nu2/2)) *
   x^(nu1/2 - 1) (nu2 + nu1 * x)^(-nu1/2 -nu2/2) dx

   The method used here is the one described in Knuth */

MpIeee gsl_ran_fdist(const gsl_rng * r, const MpIeee nu1, const MpIeee nu2)
{

  MpIeee Y1=   gsl_ran_gamma (r, nu1 / MpIeee( "2" ), MpIeee( "2.0" ));
  MpIeee Y2=   gsl_ran_gamma (r, nu2 / MpIeee( "2" ), MpIeee( "2.0" ));

  MpIeee f=  (Y1 * nu2) / (Y2 * nu1);

  return f;
}

MpIeee gsl_ran_fdist_pdf(const MpIeee x, const MpIeee nu1, const MpIeee nu2)
{
  if (x < 0)
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee p;
      MpIeee lglg=  (nu1 / MpIeee( "2" )) * log (nu1) + (nu2 / MpIeee( "2" )) * log (nu2) ;

      MpIeee lg12=  gsl_sf_lngamma ((nu1 + nu2) / MpIeee( "2" ));
      MpIeee lg1=  gsl_sf_lngamma (nu1 / MpIeee( "2" ));
      MpIeee lg2=  gsl_sf_lngamma (nu2 / MpIeee( "2" ));
      
      p = exp (lglg + lg12 - lg1 - lg2)
        * pow (x, nu1 / MpIeee( "2" ) - MpIeee( "1" )) * pow (nu2 + nu1 * x, -nu1 / MpIeee( "2" ) - nu2 / MpIeee( "2" ));
      
      return p;
    }
}
