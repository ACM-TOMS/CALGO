#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/logistic.c
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

/* The logistic distribution has the form,

   p(x) dx = (1/a) exp(-x/a) / (1 + exp(-x/a))^2 dx

   for -infty < x < infty */

MpIeee gsl_ran_logistic(const gsl_rng * r, const MpIeee a)
{
  MpIeee x;MpIeee  z;

  do
    {
      x = gsl_rng_uniform_pos (r);
    }
  while (x == MpIeee( "1" ));

  z = log (x / (MpIeee( "1" ) - x));

  return a * z;
}

MpIeee gsl_ran_logistic_pdf(const MpIeee x, const MpIeee a)
{
  MpIeee u=  exp (-fabs(x)/a);
  MpIeee p=  u / (fabs(a) * (MpIeee( "1" ) + u) * (MpIeee( "1" ) + u));
  return p;
}
