#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multimin/steepest_descent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
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

/* steepest_descent.c -- the steepest descent algorithm */

/* Modified by Brian Gough to use single iteration structure */

#include <config.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  MpIeee step;
  MpIeee max_step;
  MpIeee tol;
  gsl_vector *x1;
  gsl_vector *g1;
}
steepest_descent_state_t;

static int
 steepest_descent_alloc(void *vstate, size_t n)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  state->x1 = gsl_vector_alloc (n);

  if (state->x1 == NULL)
    {
      GSL_ERROR ("failed to allocate space for x1", GSL_ENOMEM);
    }

  state->g1 = gsl_vector_alloc (n);

  if (state->g1 == NULL)
    {
      gsl_vector_free (state->x1);
      GSL_ERROR ("failed to allocate space for g1", GSL_ENOMEM);
    }

  return GSL_SUCCESS;
}

static int
 steepest_descent_set(void *vstate, gsl_multimin_function_fdf * fdf,
                      const gsl_vector * x, MpIeee *f,
                      gsl_vector * gradient, MpIeee step_size, MpIeee tol)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, x, f, gradient);

  state->step = step_size;
  state->max_step = step_size;
  state->tol = tol;

  return GSL_SUCCESS;
}


static void
steepest_descent_free (void *vstate)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  gsl_vector_free (state->x1);
  gsl_vector_free (state->g1);
}

static int
 steepest_descent_restart(void *vstate)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  state->step = state->max_step;

  return GSL_SUCCESS;
}

static int
 steepest_descent_iterate(void *vstate, gsl_multimin_function_fdf * fdf,
                          gsl_vector * x, MpIeee *f,
                          gsl_vector * gradient, gsl_vector * dx)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  gsl_vector *x1 = state->x1;
  gsl_vector *g1 = state->g1;

  MpIeee f0=  *f;
  MpIeee f1;
  MpIeee step=  state->step;MpIeee  tol=  state->tol;

  int  failed=  0;

  /* compute new trial point at x1= x - step * dir, where dir is the
     normalized gradient */

  MpIeee gnorm=  gsl_blas_dnrm2 (gradient);

  if (gnorm == MpIeee( "0.0" ))
    {
      gsl_vector_set_zero (dx);
      return GSL_ENOPROG;
    }

trial:
  gsl_vector_set_zero (dx);
  gsl_blas_daxpy (-step / gnorm, gradient, dx);

  gsl_vector_memcpy (x1, x);
  gsl_blas_daxpy (1.0, dx, x1);

  /* evaluate function and gradient at new point x1 */

  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, x1, &f1, g1);

  if (f1 > f0)
    {
      /* downhill step failed, reduce step-size and try again */

      failed = 1;
      step *= tol;
      goto trial;
    }

  if (failed)
    step *= tol;
  else
    step *= MpIeee( "2.0" );

  state->step = step;

  gsl_vector_memcpy (x, x1);
  gsl_vector_memcpy (gradient, g1);

  *f = f1;

  return GSL_SUCCESS;
}

static const gsl_multimin_fdfminimizer_type steepest_descent_type =
  { "steepest_descent",         /* name */
  sizeof (steepest_descent_state_t),
  &steepest_descent_alloc,
  &steepest_descent_set,
  &steepest_descent_iterate,
  &steepest_descent_restart,
  &steepest_descent_free
};

const gsl_multimin_fdfminimizer_type
  * gsl_multimin_fdfminimizer_steepest_descent = &steepest_descent_type;
