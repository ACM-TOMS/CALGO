#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/falsepos.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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

/* falsepos.c -- falsepos root finding algorithm 

   The false position algorithm uses bracketing by linear interpolation.

   If a linear interpolation step would decrease the size of the
   bracket by less than a bisection step would then the algorithm
   takes a bisection step instead.
   
   The last linear interpolation estimate of the root is used. If a
   bisection step causes it to fall outside the brackets then it is
   replaced by the bisection estimate (x_upper + x_lower)/2.

*/

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "roots.h"

typedef struct
  {
    MpIeee f_lower;MpIeee  f_upper;
  }
falsepos_state_t;

static int  falsepos_init(void * vstate, gsl_function * f, MpIeee * root, MpIeee x_lower, MpIeee x_upper);
static int  falsepos_iterate(void * vstate, gsl_function * f, MpIeee * root, MpIeee * x_lower, MpIeee * x_upper);

static int
 falsepos_init(void * vstate, gsl_function * f, MpIeee * root, MpIeee x_lower, MpIeee x_upper)
{
  falsepos_state_t * state = (falsepos_state_t *) vstate;

  MpIeee f_lower;MpIeee  f_upper;

  *root = MpIeee( "0.5" ) * (x_lower + x_upper);

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  
  state->f_lower = f_lower;
  state->f_upper = f_upper;

  if ((f_lower < MpIeee( "0.0" ) && f_upper < MpIeee( "0.0" )) || (f_lower > MpIeee( "0.0" ) && f_upper > MpIeee( "0.0" )))
    {
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

static int
 falsepos_iterate(void * vstate, gsl_function * f, MpIeee * root, MpIeee * x_lower, MpIeee * x_upper)
{
  falsepos_state_t * state = (falsepos_state_t *) vstate;

  MpIeee x_linear;MpIeee  f_linear;
  MpIeee x_bisect;MpIeee  f_bisect;

  MpIeee x_left=  *x_lower ;
  MpIeee x_right=  *x_upper ;

  MpIeee f_lower=  state->f_lower; 
  MpIeee f_upper=  state->f_upper;

  MpIeee w;

  if (f_lower == MpIeee( "0.0" ))
    {
      *root = x_left ;
      *x_upper = x_left;
      return GSL_SUCCESS;
    }
  
  if (f_upper == MpIeee( "0.0" ))
    {
      *root = x_right ;
      *x_lower = x_right;
      return GSL_SUCCESS;
    }
      
  /* Draw a line between f(*lower_bound) and f(*upper_bound) and
     note where it crosses the X axis; that's where we will
     split the interval. */
  
  x_linear = x_right - (f_upper * (x_left - x_right) / (f_lower - f_upper));

  SAFE_FUNC_CALL (f, x_linear, &f_linear);
      
  if (f_linear == MpIeee( "0.0" ))
    {
      *root = x_linear;
      *x_lower = x_linear;
      *x_upper = x_linear;
      return GSL_SUCCESS;
    }
      
  /* Discard the half of the interval which doesn't contain the root. */
  
  if ((f_lower > MpIeee( "0.0" ) && f_linear < MpIeee( "0.0" )) || (f_lower < MpIeee( "0.0" ) && f_linear > MpIeee( "0.0" )))
    {
      *root = x_linear ;
      *x_upper = x_linear;
      state->f_upper = f_linear;
      w = x_linear - x_left ;
    }
  else
    {
      *root = x_linear ;
      *x_lower = x_linear;
      state->f_lower = f_linear;
      w = x_right - x_linear;
    }

  if (w < MpIeee( "0.5" ) * (x_right - x_left)) 
    {
      return GSL_SUCCESS ;
    }

  x_bisect = MpIeee( "0.5" ) * (x_left + x_right);

  SAFE_FUNC_CALL (f, x_bisect, &f_bisect);

  if ((f_lower > MpIeee( "0.0" ) && f_bisect < MpIeee( "0.0" )) || (f_lower < MpIeee( "0.0" ) && f_bisect > MpIeee( "0.0" )))
    {
      *x_upper = x_bisect;
      state->f_upper = f_bisect;
      if (*root > x_bisect)
        *root = MpIeee( "0.5" ) * (x_left + x_bisect) ;
    }
  else
    {
      *x_lower = x_bisect;
      state->f_lower = f_bisect;
      if (*root < x_bisect)
        *root = MpIeee( "0.5" ) * (x_bisect + x_right) ;
    }

  return GSL_SUCCESS;
}


static const gsl_root_fsolver_type falsepos_type =
{"falsepos",                            /* name */
 sizeof (falsepos_state_t),
 &falsepos_init,
 &falsepos_iterate};

const gsl_root_fsolver_type  * gsl_root_fsolver_falsepos = &falsepos_type;
