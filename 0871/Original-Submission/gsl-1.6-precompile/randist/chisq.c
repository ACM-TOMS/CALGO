#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/chisq.c
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

/* The chisq distribution has the form

   p(x) dx = (1/(2*Gamma(nu/2))) (x/2)^(nu/2 - 1) exp(-x/2) dx

   for x = 0 ... +infty */

MpIeee gsl_ran_chisq(const gsl_rng * r, const MpIeee nu)
{
  MpIeee chisq=  MpIeee( "2" ) * gsl_ran_gamma (r, nu / MpIeee( "2" ), MpIeee( "1.0" ));
  return chisq;
}

MpIeee gsl_ran_chisq_pdf(const MpIeee x, const MpIeee nu)
{
  if (x <= 0)
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee p;
      MpIeee lngamma=  gsl_sf_lngamma (nu / MpIeee( "2" ));
      
      p = exp ((nu / MpIeee( "2" ) - MpIeee( "1" )) * log (x/MpIeee( "2" )) - x/MpIeee( "2" ) - lngamma) / MpIeee( "2" );
      return p;
    }
}
