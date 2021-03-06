#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/pascal.c
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

/* The Pascal distribution is a negative binomial with valued integer n

   prob(k) =  (n - 1 + k)!/(n!(k - 1)!) *  p^n (1-p)^k for k = 0, 1, ..., n

   */

unsigned int
 gsl_ran_pascal(const gsl_rng * r, MpIeee p, unsigned int  n)
{
  /* This is a separate interface for the pascal distribution so that
     it can be optimized differently from the negative binomial in
     future.
     
     e.g. if n < 10 it might be faster to generate the Pascal
     distributions as the sum of geometric variates directly.  */
  
  unsigned int  k=  gsl_ran_negative_binomial (r, p, (double) n);
  return k;
}

MpIeee gsl_ran_pascal_pdf(const unsigned int  k, const MpIeee p, unsigned int  n)
{
  MpIeee P=  gsl_ran_negative_binomial_pdf (k, p, (MpIeee) n);
  return P;
}
