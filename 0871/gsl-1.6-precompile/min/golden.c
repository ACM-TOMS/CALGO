#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* min/golden.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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


/* goldensection.c -- goldensection minimum finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

#include "min.h"

typedef struct
  {
    MpIeee dummy;
  }
goldensection_state_t;

static int  goldensection_init(void * vstate, gsl_function * f, MpIeee x_minimum, MpIeee f_minimum, MpIeee x_lower, MpIeee f_lower, MpIeee x_upper, MpIeee f_upper);
static int  goldensection_iterate(void * vstate, gsl_function * f, MpIeee * x_minimum, MpIeee * f_minimum, MpIeee * x_lower, MpIeee * f_lower, MpIeee * x_upper, MpIeee * f_upper);

static int
 goldensection_init(void * vstate, gsl_function * f, MpIeee x_minimum, MpIeee f_minimum, MpIeee x_lower, MpIeee f_lower, MpIeee x_upper, MpIeee f_upper)
{
  goldensection_state_t * state = (goldensection_state_t *) vstate;

  /* no initialization required, prevent warnings about unused variables */

  state = 0;
  f = 0;
  x_minimum = MpIeee( "0" );
  f_minimum = MpIeee( "0" );
  x_lower = MpIeee( "0" );
  f_lower = MpIeee( "0" );
  x_upper = MpIeee( "0" );
  f_upper = MpIeee( "0" );

  return GSL_SUCCESS;
}

static int
 goldensection_iterate(void * vstate, gsl_function * f, MpIeee * x_minimum, MpIeee * f_minimum, MpIeee * x_lower, MpIeee * f_lower, MpIeee * x_upper, MpIeee * f_upper)
{
  goldensection_state_t * state = (goldensection_state_t *) vstate;

  const MpIeee x_center=  *x_minimum ;
  const MpIeee x_left=  *x_lower ;
  const MpIeee x_right=  *x_upper ;

  const MpIeee f_min=  *f_minimum;

  const MpIeee golden=  0.3819660; /* golden = (3 - sqrt(5))/2 */
  
  const MpIeee w_lower=  (x_center - x_left);
  const MpIeee w_upper=  (x_right - x_center);

  MpIeee x_new;MpIeee  f_new;

  state = 0 ; /* avoid warning about unused parameters */
  
  x_new = x_center + golden * ((w_upper > w_lower) ? w_upper : -w_lower) ;

  SAFE_FUNC_CALL (f, x_new, (&f_new));

  if (f_new < f_min)
    {
      *x_minimum = x_new ;
      *f_minimum = f_new ;
      return GSL_SUCCESS;
    }
  else if (x_new < x_center && f_new > f_min)
    {
      *x_lower = x_new ;
      *f_lower = f_new ;
      return GSL_SUCCESS;
    }
  else if (x_new > x_center && f_new > f_min)
    {
      *x_upper = x_new ;
      *f_upper = f_new ;
      return GSL_SUCCESS;
    }
  else
    {
      return GSL_FAILURE;
    }
}


static const gsl_min_fminimizer_type goldensection_type =
{"goldensection",                               /* name */
 sizeof (goldensection_state_t),
 &goldensection_init,
 &goldensection_iterate};

const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection = &goldensection_type;
