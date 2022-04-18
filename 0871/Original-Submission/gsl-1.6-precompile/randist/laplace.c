#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/laplace.c
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

/* The two-sided exponential probability distribution is  

   p(x) dx = (1/(2 a)) * exp(-|x/a|) dx

   for -infty < x < infty. It is also known as the Laplace distribution.  */

MpIeee gsl_ran_laplace(const gsl_rng * r, const MpIeee a)
{
  MpIeee u;
  do
    {
      u = MpIeee( "2" ) * gsl_rng_uniform (r) - MpIeee( "1.0" );
    }
  while (u == MpIeee( "0.0" ));

  if (u < MpIeee( "0" ))
    {
      return a * log (-u);
    }
  else
    {
      return -a * log (u);
    }
}

MpIeee gsl_ran_laplace_pdf(const MpIeee x, const MpIeee a)
{
  MpIeee p=  (MpIeee( "1" )/(MpIeee( "2" )*a)) * exp (-fabs (x)/a);
  return p;
}

