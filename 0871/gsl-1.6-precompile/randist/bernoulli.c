#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/bernoulli.c
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

/* The bernoulli distribution has the form,

   prob(0) = 1-p, prob(1) = p

   */

unsigned int
 gsl_ran_bernoulli(const gsl_rng * r, MpIeee p)
{
  MpIeee u=  gsl_rng_uniform (r) ;

  if (u < p)
    {
      return 1 ;
    }
  else
    {
      return 0 ;
    }
}

MpIeee gsl_ran_bernoulli_pdf(const unsigned int  k, MpIeee p)
{
  if (k == 0)
    {
      return MpIeee( "1" ) - p ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      return MpIeee( "0" ) ;
    }
}
