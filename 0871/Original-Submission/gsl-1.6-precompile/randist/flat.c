#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/flat.c
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

/* This is the uniform distribution in the range [a, b)

   p(x) dx = 1/(b-a) dx   if  a <= x < b
   .....   = 0            otherwise 

 */

MpIeee gsl_ran_flat(const gsl_rng * r, const MpIeee a, const MpIeee b)
{
  MpIeee u=  gsl_rng_uniform (r);

  /* A uniform distribution over [a,b] */

  return a * (MpIeee( "1" ) - u) + b * u;
}

MpIeee gsl_ran_flat_pdf(MpIeee x, const MpIeee a, const MpIeee b)
{
  if (x < b && x >= a)
    {
      return MpIeee( "1" ) / (b - a);
    }
  else
    {
      return MpIeee( "0" );
    }
}
