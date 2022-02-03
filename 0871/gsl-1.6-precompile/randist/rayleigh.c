#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/rayleigh.c
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

/* The Rayleigh distribution has the form

   p(x) dx = (x / sigma^2) exp(-x^2/(2 sigma^2)) dx

   for x = 0 ... +infty */

MpIeee gsl_ran_rayleigh(const gsl_rng * r, const MpIeee sigma)
{
  MpIeee u=  gsl_rng_uniform_pos (r);

  return sigma * sqrt(-MpIeee( "2.0" ) * log (u));
}

MpIeee gsl_ran_rayleigh_pdf(const MpIeee x, const MpIeee sigma)
{
  if (x < 0)
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee u=  x / sigma ;
      MpIeee p=  (u / sigma) * exp(-u * u / MpIeee( "2.0" )) ;
      
      return p;
    }
}

/* The Rayleigh tail distribution has the form

   p(x) dx = (x / sigma^2) exp((a^2 - x^2)/(2 sigma^2)) dx

   for x = a ... +infty */

MpIeee gsl_ran_rayleigh_tail(const gsl_rng * r, const MpIeee a, const MpIeee sigma)
{
  MpIeee u=  gsl_rng_uniform_pos (r);

  return sqrt(a * a - MpIeee( "2.0" ) * sigma * sigma * log (u));
}

MpIeee gsl_ran_rayleigh_tail_pdf(const MpIeee x, const MpIeee a, const MpIeee sigma)
{
  if (x < a)
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee u=  x / sigma ;
      MpIeee v=  a / sigma ;

      MpIeee p=  (u / sigma) * exp((v + u) * (v - u) / MpIeee( "2.0" )) ;
      
      return p;
    }
}
