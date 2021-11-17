#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/interp.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#define DISCARD_STATUS(s) if ((s) != GSL_SUCCESS) { GSL_ERROR_VAL("interpolation error", (s),  GSL_NAN); }

gsl_interp *
gsl_interp_alloc (const gsl_interp_type * T, size_t size)
{
  gsl_interp * interp;

  if (size < T->min_size)
    {
      GSL_ERROR_NULL ("insufficient number of points for interpolation type",
                      GSL_EINVAL);
    }

  interp = (gsl_interp *) malloc (sizeof(gsl_interp));
  
  if (interp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for interp struct", 
                      GSL_ENOMEM);
    }
  
  interp->type = T;
  interp->size = size;

  if (interp->type->alloc == NULL)
    {
      interp->state = NULL;
      return interp;
    }

  interp->state = interp->type->alloc(size);
  
  if (interp->state == NULL)
    {
      free (interp);          
      GSL_ERROR_NULL ("failed to allocate space for interp state", GSL_ENOMEM);
    };
    
  return interp;
}

int
 gsl_interp_init(gsl_interp * interp, const MpIeee x_array[], const MpIeee y_array[], size_t size)
{
  if (size != interp->size)
    {
      GSL_ERROR ("data must match size of interpolation object", GSL_EINVAL);
    }
  
  interp->xmin = x_array[0];
  interp->xmax = x_array[size - 1];

  {
    int  status=  interp->type->init(interp->state, x_array, y_array, size);
    return status;
  }
}

const char *
gsl_interp_name(const gsl_interp * interp)
{
  return interp->type->name;
}

unsigned int
 gsl_interp_min_size(const gsl_interp * interp)
{
  return interp->type->min_size;
}

void
gsl_interp_free (gsl_interp * interp)
{
  if (interp->type->free)
    interp->type->free (interp->state);
  free (interp);
}



int
 gsl_interp_eval_e(const gsl_interp * interp,
                   const MpIeee xa[], const MpIeee ya[], MpIeee x,
                   gsl_interp_accel * a, MpIeee *y)
{
  if (x < interp->xmin)
    {
      *y = ya[0];
      return GSL_EDOM;
    }
  else if (x > interp->xmax)
    {
      *y = ya[interp->size - 1];
      return GSL_EDOM;
    }

  return interp->type->eval (interp->state, xa, ya, interp->size, x, a, y);
}

MpIeee gsl_interp_eval(const gsl_interp * interp,
                 const MpIeee xa[], const MpIeee ya[], MpIeee x,
                 gsl_interp_accel * a)
{
  MpIeee y;
  int  status=  interp->type->eval (interp->state, xa, ya, interp->size, x, a, &y);

  DISCARD_STATUS(status);

  return y;
}


int
 gsl_interp_eval_deriv_e(const gsl_interp * interp,
                         const MpIeee xa[], const MpIeee ya[], MpIeee x,
                         gsl_interp_accel * a,
                         MpIeee *dydx)
{
  if (x < interp->xmin)
    {
      *dydx = MpIeee( "0.0" );
      return GSL_EDOM;
    }
  else if (x > interp->xmax)
    {
      *dydx = MpIeee( "0.0" );
      return GSL_EDOM;
    }

  return interp->type->eval_deriv (interp->state, xa, ya, interp->size, x, a, dydx);
}

MpIeee gsl_interp_eval_deriv(const gsl_interp * interp,
                       const MpIeee xa[], const MpIeee ya[], MpIeee x,
                       gsl_interp_accel * a)
{
  MpIeee dydx;
  int  status=  interp->type->eval_deriv (interp->state, xa, ya, interp->size, x, a, &dydx);

  DISCARD_STATUS(status);

  return dydx;
}


int
 gsl_interp_eval_deriv2_e(const gsl_interp * interp,
                          const MpIeee xa[], const MpIeee ya[], MpIeee x,
                          gsl_interp_accel * a,
                          MpIeee * d2)
{
  if (x < interp->xmin)
    {
      *d2 = MpIeee( "0.0" );
      return GSL_EDOM;
    }
  else if (x > interp->xmax)
    {
      *d2 = MpIeee( "0.0" );
      return GSL_EDOM;
    }

  return interp->type->eval_deriv2 (interp->state, xa, ya, interp->size, x, a, d2);
}

MpIeee gsl_interp_eval_deriv2(const gsl_interp * interp,
                        const MpIeee xa[], const MpIeee ya[], MpIeee x,
                        gsl_interp_accel * a)
{
  MpIeee d2;
  int  status=  interp->type->eval_deriv2 (interp->state, xa, ya, interp->size, x, a, &d2);

  DISCARD_STATUS(status);

  return d2;
}


int
 gsl_interp_eval_integ_e(const gsl_interp * interp,
                         const MpIeee xa[], const MpIeee ya[],
                         MpIeee a, MpIeee b,
                         gsl_interp_accel * acc,
                         MpIeee * result)
{
  if (a > b || a < interp->xmin || b > interp->xmax)
    {
      *result = MpIeee( "0.0" );
      return GSL_EDOM;
    }
  else if(a == b)
    {
      *result = MpIeee( "0.0" );
      return GSL_SUCCESS;
    }

  return interp->type->eval_integ (interp->state, xa, ya, interp->size, acc, a, b, result);
}


MpIeee gsl_interp_eval_integ(const gsl_interp * interp,
                       const MpIeee xa[], const MpIeee ya[],
                       MpIeee a, MpIeee b,
                       gsl_interp_accel * acc)
{
  MpIeee result;
  int  status=  interp->type->eval_integ (interp->state, xa, ya, interp->size, acc, a, b, &result);

  DISCARD_STATUS(status);

  return result;
}


