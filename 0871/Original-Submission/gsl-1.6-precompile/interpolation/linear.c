#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/linear.c
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
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

static int
 linear_init(void * vstate,
             const MpIeee x_array[],
             const MpIeee y_array[],
             size_t size)
{
  return GSL_SUCCESS;
}

static
int
 linear_eval(const void * vstate,
             const MpIeee x_array[], const MpIeee y_array[], size_t size,
             MpIeee x,
             gsl_interp_accel * a,
             MpIeee *y)
{
  MpIeee x_lo;MpIeee  x_hi;
  MpIeee y_lo;MpIeee  y_hi;
  MpIeee dx;
  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  x_lo = x_array[index];
  x_hi = x_array[index + 1];
  y_lo = y_array[index];
  y_hi = y_array[index + 1];
  dx = x_hi - x_lo;
  if (dx > MpIeee( "0.0" ))
    {
      *y = y_lo + (x - x_lo) / dx * (y_hi - y_lo);
      return GSL_SUCCESS;
    }
  else
    {
      *y = MpIeee( "0.0" );
      return GSL_EINVAL;
    }
}


static
int
 linear_eval_deriv(const void * vstate,
                   const MpIeee x_array[], const MpIeee y_array[], size_t size,
                   MpIeee x,
                   gsl_interp_accel * a,
                   MpIeee *dydx)
{
  MpIeee x_lo;MpIeee  x_hi;
  MpIeee y_lo;MpIeee  y_hi;
  MpIeee dx;
  MpIeee dy;
  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  x_lo = x_array[index];
  x_hi = x_array[index + 1];
  y_lo = y_array[index];
  y_hi = y_array[index + 1];
  dx = x_hi - x_lo;
  dy = y_hi - y_lo;
  if (dx > MpIeee( "0.0" ))
    {
      *dydx = dy / dx;;
      return GSL_SUCCESS;
    }
  else
    {
      *dydx = MpIeee( "0.0" );
      return GSL_EINVAL;
    }
}


static
int
 linear_eval_deriv2(const void * vstate,
                    const MpIeee x_array[], const MpIeee y_array[], size_t size,
                    MpIeee x,
                    gsl_interp_accel * a,
                    MpIeee *y_pp)
{
  *y_pp = MpIeee( "0.0" );

  return GSL_SUCCESS;
}


static
int
 linear_eval_integ(const void * vstate,
                   const MpIeee x_array[], const MpIeee y_array[], size_t size,
                   gsl_interp_accel * acc,
                   MpIeee a, MpIeee b,
                   MpIeee * result)
{
  size_t i, index_a, index_b;
  
  if (acc != 0)
    {
      index_a = gsl_interp_accel_find (acc, x_array, size, a);
      index_b = gsl_interp_accel_find (acc, x_array, size, b);
    }
  else
    {
      index_a = gsl_interp_bsearch (x_array, a, 0, size - 1);
      index_b = gsl_interp_bsearch (x_array, b, 0, size - 1);
    }
  
    /* endpoints span more than one interval */

  *result = MpIeee( "0.0" );
  
  /* interior intervals */
  for(i=index_a; i<=index_b; i++) {
    const MpIeee x_hi=  x_array[i + 1];
    const MpIeee x_lo=  x_array[i];
    const MpIeee y_lo=  y_array[i];
    const MpIeee y_hi=  y_array[i + 1];
    const MpIeee dx=  x_hi - x_lo;

    if(dx != 0.0) {
      if (i == index_a || i == index_b)
        {
          MpIeee x1=  (i == index_a) ? a : x_lo;
          MpIeee x2=  (i == index_b) ? b : x_hi;
          const MpIeee D=  (y_hi-y_lo)/dx;
          *result += (x2-x1) * (y_lo + MpIeee( "0.5" )*D*((x2-x_lo)+(x1-x_lo)));
        }
      else
        {
          *result += MpIeee( "0.5" ) * dx * (y_lo + y_hi);
        }
    }
  }
    
  return GSL_SUCCESS;
}

static const gsl_interp_type linear_type = 
{
  "linear", 
  2,
  NULL, /* alloc, not applicable */
  &linear_init,
  &linear_eval,
  &linear_eval_deriv,
  &linear_eval_deriv2,
  &linear_eval_integ,
  NULL, /* free, not applicable */
};

const gsl_interp_type * gsl_interp_linear = &linear_type;
