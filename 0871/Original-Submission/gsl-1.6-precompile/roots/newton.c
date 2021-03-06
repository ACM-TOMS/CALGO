#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/newton.c
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

/* newton.c -- newton root finding algorithm 

   This is the classical Newton-Raphson iteration.

   x[i+1] = x[i] - f(x[i])/f'(x[i])

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
  }
newton_state_t;

static int  newton_init(void * vstate, gsl_function_fdf * fdf, MpIeee * root);
static int  newton_iterate(void * vstate, gsl_function_fdf * fdf, MpIeee * root);

static int
 newton_init(void * vstate, gsl_function_fdf * fdf, MpIeee * root)
{
  newton_state_t * state = (newton_state_t *) vstate;

  const MpIeee x=  *root ;

  state->f = GSL_FN_FDF_EVAL_F (fdf, x);
  state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;

  return GSL_SUCCESS;

}

static int
 newton_iterate(void * vstate, gsl_function_fdf * fdf, MpIeee * root)
{
  newton_state_t * state = (newton_state_t *) vstate;
  
  MpIeee root_new;MpIeee  f_new;MpIeee  df_new;

  if (state->df == 0.0)
    {
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    }

  root_new = *root - (state->f / state->df);

  *root = root_new ;
  
  GSL_FN_FDF_EVAL_F_DF(fdf, root_new, &f_new, &df_new);

  state->f = f_new ;
  state->df = df_new ;

/*  if (!finite(f_new))
    {
      GSL_ERROR ("function not continuous", GSL_EBADFUNC);
    }

  if (!finite (df_new))
    {
      GSL_ERROR ("function not differentiable", GSL_EBADFUNC);
    }
 */     
  return GSL_SUCCESS;
}


static const gsl_root_fdfsolver_type newton_type =
{"newton",                              /* name */
 sizeof (newton_state_t),
 &newton_init,
 &newton_iterate};

const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_newton = &newton_type;
