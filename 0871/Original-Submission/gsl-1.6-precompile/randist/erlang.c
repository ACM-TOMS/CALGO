#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/erlang.c
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

/* The sum of N samples from an exponential distribution gives an
   Erlang distribution

   p(x) dx = x^(n-1) exp (-x/a) / ((n-1)!a^n) dx

   for x = 0 ... +infty */

MpIeee gsl_ran_erlang(const gsl_rng * r, const MpIeee a, const MpIeee n)
{
  return gsl_ran_gamma (r, n, a);
}

MpIeee gsl_ran_erlang_pdf(const MpIeee x, const MpIeee a, const MpIeee n)
{
  if (x <= 0) 
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee p;
      MpIeee lngamma=  gsl_sf_lngamma (n);

      p = exp ((n - MpIeee( "1" )) * log (x/a) - x/a - lngamma) / a;
      return p;
    }
}
