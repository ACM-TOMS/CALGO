#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* monte/plain.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
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

/* Plain Monte-Carlo. */

/* Author: MJB */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_plain.h>

int
 gsl_monte_plain_integrate(const gsl_monte_function * f,
                           const MpIeee xl[], const MpIeee xu[],
                           const size_t dim,
                           const size_t calls,
                           gsl_rng * r,
                           gsl_monte_plain_state * state,
                           MpIeee *result, MpIeee *abserr)
{
  MpIeee vol;MpIeee  m=  MpIeee( "0" );MpIeee  q=  MpIeee( "0" );
  MpIeee *x=  state->x;
  size_t n, i;

  if (dim != state->dim)
    {
      GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
    }

  for (i = 0; i < dim; i++)
    {
      if (xu[i] <= xl[i])
        {
          GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
        }

      if (xu[i] - xl[i] > GSL_DBL_MAX)
        {
          GSL_ERROR ("Range of integration is too large, please rescale",
                     GSL_EINVAL);
        }
    }

  /* Compute the volume of the region */

  vol = MpIeee( "1" );

  for (i = 0; i < dim; i++)
    {
      vol *= xu[i] - xl[i];
    }

  for (n = 0; n < calls; n++)
    {
      /* Choose a random point in the integration region */

      for (i = 0; i < dim; i++)
        {
          x[i] = xl[i] + gsl_rng_uniform_pos (r) * (xu[i] - xl[i]);
        }

      {
        MpIeee fval=  GSL_MONTE_FN_EVAL (f, x);

        /* recurrence for mean and variance */

        MpIeee d=  fval - m;
        m += d / (n + MpIeee( "1.0" ));
        q += d * d * (n / (n + MpIeee( "1.0" )));
      }
    }

  *result = vol * m;

  if (calls < 2)
    {
      *abserr = GSL_POSINF;
    }
  else
    {
      *abserr = vol * sqrt (q / (calls * (calls - MpIeee( "1.0" ))));
    }

  return GSL_SUCCESS;
}

gsl_monte_plain_state *
gsl_monte_plain_alloc (size_t dim)
{
  gsl_monte_plain_state *s =
    (gsl_monte_plain_state *) malloc (sizeof (gsl_monte_plain_state));

  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for state struct",
                     GSL_ENOMEM, 0);
    }

  s->x = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (s->x == 0)
    {
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for working vector",
                     GSL_ENOMEM, 0);
    }

  s->dim = dim;

  return s;
}

/* Set some default values and whatever */

int
 gsl_monte_plain_init(gsl_monte_plain_state * s)
{
  size_t i;

  for (i = 0; i < s->dim; i++)
    {
      s->x[i] = 0.0;
    }

  return GSL_SUCCESS;
}

void
gsl_monte_plain_free (gsl_monte_plain_state * s)
{
  free (s->x);
  free (s);
}
