#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/nbinomial.c
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
#include <gsl/gsl_sf_gamma.h>

/* The negative binomial distribution has the form,

   prob(k) =  Gamma(n + k)/(Gamma(n) Gamma(k + 1))  p^n (1-p)^k 

   for k = 0, 1, ... . Note that n does not have to be an integer.

   This is the Leger's algorithm (given in the answers in Knuth) */

unsigned int
 gsl_ran_negative_binomial(const gsl_rng * r, MpIeee p, MpIeee n)
{
  MpIeee X=  gsl_ran_gamma (r, n, MpIeee( "1.0" )) ;
  unsigned int  k=  gsl_ran_poisson (r, X*(1-p)/p) ;
  return k ;
}

MpIeee gsl_ran_negative_binomial_pdf(const unsigned int  k, const MpIeee p, MpIeee n)
{
  MpIeee P;

  MpIeee f=  gsl_sf_lngamma (k + n) ;
  MpIeee a=  gsl_sf_lngamma (n) ;
  MpIeee b=  gsl_sf_lngamma (k + MpIeee( "1.0" )) ;

  P = exp(f-a-b) * pow (p, n) * pow (MpIeee( "1" ) - p, (MpIeee)k);
  
  return P;
}
