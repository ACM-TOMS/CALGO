#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/gumbel.c
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

/* The Type I Gumbel distribution has the form,

   p(x) dx = a b exp(-(b exp(-ax) + ax)) dx

   and the Type II Gumbel distribution has the form,

   p(x) dx = b a x^-(a+1) exp(-b x^-a)) dx

 */

MpIeee gsl_ran_gumbel1(const gsl_rng * r, const MpIeee a, const MpIeee b)
{
  MpIeee x=  gsl_rng_uniform_pos (r);

  MpIeee z=  (log(b) - log(-log(x))) / a;

  return z;
}

MpIeee gsl_ran_gumbel1_pdf(const MpIeee x, const MpIeee a, const MpIeee b)
{
  MpIeee p=  a * b *  exp (-(b * exp(-a * x) + a * x));
  return p;
}

MpIeee gsl_ran_gumbel2(const gsl_rng * r, const MpIeee a, const MpIeee b)
{
  MpIeee x=  gsl_rng_uniform_pos (r);

  MpIeee z=  pow(-b / log(x), MpIeee( "1" )/a);

  return z;
}

MpIeee gsl_ran_gumbel2_pdf(const MpIeee x, const MpIeee a, const MpIeee b)
{
  if (x <= 0)
    {
      return MpIeee( "0" ) ;
    }
  else
    {
      MpIeee p=  b * a *  pow(x,-(a+MpIeee( "1" ))) * exp (-b * pow(x, -a));
      return p;
    }
}




