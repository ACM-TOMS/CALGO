#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/steffenson.c
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

/* steffenson.c -- steffenson root finding algorithm 

   This is Newton's method with an Aitken "delta-squared"
   acceleration of the iterates. This can improve the convergence on
   multiple roots where the ordinary Newton algorithm is slow.

   x[i+1] = x[i] - f(x[i]) / f'(x[i])

   x_accelerated[i] = x[i] - (x[i+1] - x[i])**2 / (x[i+2] - 2*x[i+1] + x[i])

   We can only use the accelerated estimate after three iterations,
   and use the unaccelerated value until then.

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
    MpIeee f;MpIeee  df;
    MpIeee x;
    MpIeee x_1;
    MpIeee x_2;
    int  count;
  }
steffenson_state_t;

static int  steffenson_init(void * vstate, gsl_function_fdf * fdf, MpIeee * root);
static int  steffenson_iterate(void * vstate, gsl_function_fdf * fdf, MpIeee * root);

static int
 steffenson_init(void * vstate, gsl_function_fdf * fdf, MpIeee * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;

  const MpIeee x=  *root ;

  state->f = GSL_FN_FDF_EVAL_F (fdf, x);
  state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;

  state->x = x;
  state->x_1 = 0.0;
  state->x_2 = 0.0;

  state->count = 1;

  return GSL_SUCCESS;

}

static int
 steffenson_iterate(void * vstate, gsl_function_fdf * fdf, MpIeee * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;
  
  MpIeee x_new;MpIeee  f_new;MpIeee  df_new;

  MpIeee x_1=  state->x_1 ;
  MpIeee x=  state->x ;

  if (state->df == 0.0)
    {
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    }

  x_new = x - (state->f / state->df);
  
  GSL_FN_FDF_EVAL_F_DF(fdf, x_new, &f_new, &df_new);

  state->x_2 = x_1 ;
  state->x_1 = x ;
  state->x = x_new;

  state->f = f_new ;
  state->df = df_new ;

/*  if (!finite (f_new))
    {
      GSL_ERROR ("function not continuous", GSL_EBADFUNC);
    }
*/
  if (state->count < 3)
    {
      *root = x_new ;
      state->count++ ;
    }
  else 
    {
      MpIeee u=  (x - x_1) ;
      MpIeee v=  (x_new - MpIeee( "2" ) * x + x_1);

      if (v == MpIeee( "0" ))
        *root = x_new;  /* avoid division by zero */
      else
        *root = x_1 - u * u / v ;  /* accelerated value */
    }

/*  if (!finite (df_new))
    {
      GSL_ERROR ("function not differentiable", GSL_EBADFUNC);
    }
 */
      
  return GSL_SUCCESS;
}


static const gsl_root_fdfsolver_type steffenson_type =
{"steffenson",                          /* name */
 sizeof (steffenson_state_t),
 &steffenson_init,
 &steffenson_iterate};

const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson = &steffenson_type;
