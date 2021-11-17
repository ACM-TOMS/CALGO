#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ode-initval/evolve.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman
 */
#include <config.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

gsl_odeiv_evolve *
gsl_odeiv_evolve_alloc (size_t dim)
{
  gsl_odeiv_evolve *e =
    (gsl_odeiv_evolve *) malloc (sizeof (gsl_odeiv_evolve));

  if (e == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for evolve struct",
                      GSL_ENOMEM);
    }

  e->y0 = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (e->y0 == 0)
    {
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for y0", GSL_ENOMEM);
    }

  e->yerr = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (e->yerr == 0)
    {
      free (e->y0);
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for yerr", GSL_ENOMEM);
    }

  e->dydt_in = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (e->dydt_in == 0)
    {
      free (e->yerr);
      free (e->y0);
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for dydt_in", GSL_ENOMEM);
    }

  e->dydt_out = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (e->dydt_out == 0)
    {
      free (e->dydt_in);
      free (e->yerr);
      free (e->y0);
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for dydt_out", GSL_ENOMEM);
    }

  e->dimension = dim;
  e->count = 0;
  e->failed_steps = 0;
  e->last_step = 0.0;

  return e;
}

int
 gsl_odeiv_evolve_reset(gsl_odeiv_evolve * e)
{
  e->count = 0;
  e->failed_steps = 0;
  e->last_step = 0.0;
  return GSL_SUCCESS;
}

void
gsl_odeiv_evolve_free (gsl_odeiv_evolve * e)
{
  free (e->dydt_out);
  free (e->dydt_in);
  free (e->yerr);
  free (e->y0);
  free (e);
}

/* Evolution framework method.
 *
 * Uses an adaptive step control object
 */
int
 gsl_odeiv_evolve_apply(gsl_odeiv_evolve * e,
                        gsl_odeiv_control * con,
                        gsl_odeiv_step * step,
                        const gsl_odeiv_system * dydt,
                        MpIeee *t, MpIeee t1, MpIeee *h, MpIeee y[])
{
  const MpIeee t0=  *t;
  MpIeee h0=  *h;
  int  step_status;
  int  final_step=  0;
  MpIeee dt=  t1 - t0;  /* remaining time, possibly less than h */

  if (e->dimension != step->dimension)
    {
      GSL_ERROR ("step dimension must match evolution size", GSL_EINVAL);
    }

  if ((dt < MpIeee( "0.0" ) && h0 > MpIeee( "0.0" )) || (dt > MpIeee( "0.0" ) && h0 < MpIeee( "0.0" )))
    {
      GSL_ERROR ("step direction must match interval direction", GSL_EINVAL);
    }

  /* No need to copy if we cannot control the step size. */

  if (con != NULL)
    {
      DBL_MEMCPY (e->y0, y, e->dimension);
    }

  /* Calculate initial dydt once if the method can benefit. */

  if (step->type->can_use_dydt_in)
    {
      int  status=  GSL_ODEIV_FN_EVAL (dydt, t0, y, e->dydt_in);

      if (status) 
        {
          return GSL_EBADFUNC;
        }
    }

try_step:
    
  if ((dt >= MpIeee( "0.0" ) && h0 > dt) || (dt < MpIeee( "0.0" ) && h0 < dt))
    {
      h0 = dt;
      final_step = 1;
    }
  else
    {
      final_step = 0;
    }

  if (step->type->can_use_dydt_in)
    {
      step_status =
        gsl_odeiv_step_apply (step, t0, h0, y, e->yerr, e->dydt_in,
                              e->dydt_out, dydt);
    }
  else
    {
      step_status =
        gsl_odeiv_step_apply (step, t0, h0, y, e->yerr, NULL, e->dydt_out,
                              dydt);
    }

  /* Check for stepper internal failure */

  if (step_status != GSL_SUCCESS) 
    {
      return step_status;
    }

  e->count++;
  e->last_step = h0;

  if (final_step)
    {
      *t = t1;
    }
  else
    {
      *t = t0 + h0;
    }

  if (con != NULL)
    {
      /* Check error and attempt to adjust the step. */
      const int hadjust_status 
        = gsl_odeiv_control_hadjust (con, step, y, e->yerr, e->dydt_out, &h0);

      if (hadjust_status == GSL_ODEIV_HADJ_DEC)
        {
          /* Step was decreased. Undo and go back to try again. */
          DBL_MEMCPY (y, e->y0, dydt->dimension);
          e->failed_steps++;
          goto try_step;
        }
    }

  *h = h0;  /* suggest step size for next time-step */

  return step_status;
}
