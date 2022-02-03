#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/cauchy.c
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

/* The Cauchy probability distribution is 

   p(x) dx = (1/(pi a)) (1 + (x/a)^2)^(-1) dx

   It is also known as the Lorentzian probability distribution */

MpIeee gsl_ran_cauchy(const gsl_rng * r, const MpIeee a)
{
  MpIeee u;
  do
    {
      u = gsl_rng_uniform (r);
    }
  while (u == MpIeee( "0.5" ));

  return a * tan (M_PI * u);
}

MpIeee gsl_ran_cauchy_pdf(const MpIeee x, const MpIeee a)
{
  MpIeee u=  x / a;
  MpIeee p=  (MpIeee( "1" ) / (M_PI * a)) / (MpIeee( "1" ) + u * u);
  return p;
}
