#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/geometric.c
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

/* Geometric distribution (bernoulli trial with probability p) 

   prob(k) =  p (1 - p)^(k-1) for n = 1, 2, 3, ...

   It gives the distribution of "waiting times" for an event that
   occurs with probability p. */

unsigned int
 gsl_ran_geometric(const gsl_rng * r, const MpIeee p)
{
  MpIeee u=  gsl_rng_uniform_pos (r);

  unsigned int  k;

  if (p == 1)
    {
      k = 1;
    }
  else
    {
      k = MpIeee(log (u) / log (1 - p)).toUnsignedInt() + 1;
    }

  return k;
}

MpIeee gsl_ran_geometric_pdf(const unsigned int  k, const MpIeee p)
{
  if (k == 0)
    {
      return MpIeee( "0" ) ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      MpIeee P=  p * pow (MpIeee( "1" ) - p, k - MpIeee( "1.0" ));
      return P;
    }
}
