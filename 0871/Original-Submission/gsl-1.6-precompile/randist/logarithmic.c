#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/logarithmic.c
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

/* Logarithmic distribution 

   prob(n) =   p^n / (n log(1/(1-p)) for n = 1, 2, 3, ...

   We use Kemp's second accelerated generator, from Luc Devroye's book
   on "Non-Uniform Random Variate Generation", Springer */

unsigned int
 gsl_ran_logarithmic(const gsl_rng * r, const MpIeee p)
{
  MpIeee c=  log (MpIeee( "1" )-p) ;

  MpIeee v=  gsl_rng_uniform_pos (r);
  
  if (v >= p)
    {
      return 1 ;
    }
  else
    {
      MpIeee u=  gsl_rng_uniform_pos (r);      
      MpIeee q=  MpIeee( "1" ) - exp (c * u);

      if (v <= q*q)
        {
          MpIeee x=  MpIeee( "1" ) + log(v)/log(q) ;
          return x.toUnsignedInt() ;
        }
      else if (v <= q)
        {
          return 2;
        }
      else
        {
          return 1 ;
        }
    }
}

MpIeee gsl_ran_logarithmic_pdf(const unsigned int  k, const MpIeee p)
{
  if (k == 0)
    {
      return MpIeee( "0" ) ;
    }
  else 
    {
      MpIeee P=  pow(p, (MpIeee)k) / (MpIeee) k / log(MpIeee( "1" )/(MpIeee( "1" )-p)) ;
      return P;
    }
}
