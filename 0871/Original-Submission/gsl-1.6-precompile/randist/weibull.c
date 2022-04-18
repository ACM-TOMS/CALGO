#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/weibull.c
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

/* The Weibull distribution has the form,

   p(x) dx = (b/a) (x/a)^(b-1) exp(-(x/a)^b) dx

 */

MpIeee gsl_ran_weibull(const gsl_rng * r, const MpIeee a, const MpIeee b)
{
  MpIeee x=  gsl_rng_uniform_pos (r);

  MpIeee z=  pow (-log (x), MpIeee( "1" ) / b);

  return a * z;
}

MpIeee gsl_ran_weibull_pdf(const MpIeee x, const MpIeee a, const MpIeee b)
{
  if (x < 0)
    {
      return MpIeee( "0" ) ;
    }
  else if (x == 0)
    {
      if (b == 1)
        return MpIeee( "1" )/a ;
      else
        return MpIeee( "0" ) ;
    }
  else if (b == 1)
    {
      return exp(-x/a)/a ;
    }
  else
    {
      MpIeee p=  (b/a) * exp (-pow (x/a, b) + (b - MpIeee( "1" )) * log (x/a));
      return p;
    }
}
