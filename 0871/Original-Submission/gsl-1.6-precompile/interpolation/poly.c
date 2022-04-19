#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/interp_poly.c
 * 
 * Copyright (C) 2001 DAN, HO-JIN
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_interp.h>

typedef struct
{
  MpIeee *d;
  MpIeee *coeff;
  MpIeee *work;
}
polynomial_state_t;

static void *
polynomial_alloc (size_t size)
{
  polynomial_state_t *state =
    (polynomial_state_t *) malloc (sizeof (polynomial_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for polynomial state",
                      GSL_ENOMEM);
    }

  state->d = (MpIeee*) malloc (sizeof (MpIeee) * size);

  if (state->d == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for d", GSL_ENOMEM);
    }

  state->coeff = (MpIeee*) malloc (sizeof (MpIeee) * size);

  if (state->coeff == 0)
    {
      free (state->d);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for d", GSL_ENOMEM);
    }

  state->work = (MpIeee*) malloc (sizeof (MpIeee) * size);

  if (state->work == 0)
    {
      free (state->coeff);
      free (state->d);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for d", GSL_ENOMEM);
    }

  return state;
}

static int
 polynomial_init(void *vstate,
                 const MpIeee xa[], const MpIeee ya[], size_t size)
{
  polynomial_state_t *state = (polynomial_state_t *) vstate;

  int  status=  gsl_poly_dd_init (state->d, xa, ya, size);

  return status;
}

static int
 polynomial_eval(const void *vstate,
                 const MpIeee xa[], const MpIeee ya[], size_t size, MpIeee x,
                 gsl_interp_accel * acc, MpIeee *y)
{
  polynomial_state_t *state = (polynomial_state_t *) vstate;

  *y = gsl_poly_dd_eval (state->d, xa, size, x);

  return GSL_SUCCESS;
}


static int
 polynomial_deriv(const void *vstate,
                  const MpIeee xa[], const MpIeee ya[], size_t size, MpIeee x,
                  gsl_interp_accel * acc, MpIeee *y)
{
  polynomial_state_t *state = (polynomial_state_t *) vstate;

  gsl_poly_dd_taylor (state->coeff, x, state->d, xa, size, state->work);

  *y = state->coeff[1];

  return GSL_SUCCESS;
}

static int
 polynomial_deriv2(const void *vstate,
                   const MpIeee xa[], const MpIeee ya[], size_t size,
                   MpIeee x, gsl_interp_accel * acc, MpIeee *y)
{
  polynomial_state_t *state = (polynomial_state_t *) vstate;

  gsl_poly_dd_taylor (state->coeff, x, state->d, xa, size, state->work);

  *y = MpIeee( "2.0" ) * state->coeff[2];

  return GSL_SUCCESS;
}

static int
 polynomial_integ(const void *vstate, const MpIeee xa[], const MpIeee ya[],
                  size_t size, gsl_interp_accel * acc, MpIeee a, MpIeee b,
                  MpIeee *result)
{
  polynomial_state_t *state = (polynomial_state_t *) vstate;
  size_t i;
  MpIeee sum;

  gsl_poly_dd_taylor (state->coeff, 0.0, state->d, xa, size, state->work);

  sum = state->coeff[0] * (b - a);

  for (i = 1; i < size; i++)
    {
      sum += state->coeff[i] * (pow (b, i + MpIeee( "1" )) - pow (a, i + MpIeee( "1" ))) / (i + MpIeee( "1.0" ));
    }

  *result = sum;

  return GSL_SUCCESS;
}

static void
polynomial_free (void *vstate)
{
  polynomial_state_t *state = (polynomial_state_t *) vstate;
  free (state->d);
  free (state->coeff);
  free (state->work);
  free (state);
}

static const gsl_interp_type polynomial_type = {
  "polynomial",
  3,
  &polynomial_alloc,
  &polynomial_init,
  &polynomial_eval,
  &polynomial_deriv,
  &polynomial_deriv2,
  &polynomial_integ,
  &polynomial_free,
};

const gsl_interp_type *gsl_interp_polynomial = &polynomial_type;
