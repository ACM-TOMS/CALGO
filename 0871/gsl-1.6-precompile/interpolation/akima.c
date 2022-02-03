#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/akima.c
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
#include <math.h>
#include <gsl/gsl_errno.h>
#include "integ_eval.h"
#include <gsl/gsl_interp.h>

typedef struct
{
  MpIeee * b;
  MpIeee * c;
  MpIeee * d;
  MpIeee * _m;
} akima_state_t;


/* common creation */
static void *
akima_alloc (size_t size)
{
  akima_state_t *state = (akima_state_t *) malloc (sizeof (akima_state_t));
  
  if (state == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for state", GSL_ENOMEM);
    }
  
  state->b = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->b == NULL)
    {
      free (state);
      GSL_ERROR_NULL("failed to allocate space for b", GSL_ENOMEM);
    }
  
  state->c = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->c == NULL)
    {
      free (state->b);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for c", GSL_ENOMEM);
    }
  
  state->d = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->d == NULL)
    {
      free (state->c);
      free (state->b);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for d", GSL_ENOMEM);
    }

  state->_m = (MpIeee*) malloc ((size + 4) * sizeof (MpIeee));

  if (state->_m == NULL)
    {
      free (state->d);
      free (state->c);
      free (state->b);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for _m", GSL_ENOMEM);
    }
  
  return state;
}


/* common calculation */
static void
akima_calc (const MpIeee x_array[], MpIeee b[],  MpIeee c[],  MpIeee d[], size_t size, MpIeee m[])
{
  size_t i;

  for (i = 0; i < (size - 1); i++)
    {
      const MpIeee NE=  fabs (m[i + 1] - m[i]) + fabs (m[i - 1] - m[i - 2]);
      if (NE == 0.0)
        {
          b[i] = m[i];
          c[i] = MpIeee( "0.0" );
          d[i] = MpIeee( "0.0" );
        }
      else
        {
          const MpIeee h_i=  x_array[i + 1] - x_array[i];
          const MpIeee NE_next=  fabs (m[i + 2] - m[i + 1]) + fabs (m[i] - m[i - 1]);
          const MpIeee alpha_i=  fabs (m[i - 1] - m[i - 2]) / NE;
          MpIeee alpha_ip1;
          MpIeee tL_ip1;
          if (NE_next == 0.0)
            {
              tL_ip1 = m[i];
            }
          else
            {
              alpha_ip1 = fabs (m[i] - m[i - 1]) / NE_next;
              tL_ip1 = (MpIeee( "1.0" ) - alpha_ip1) * m[i] + alpha_ip1 * m[i + 1];
            }
          b[i] = (MpIeee( "1.0" ) - alpha_i) * m[i - 1] + alpha_i * m[i];
          c[i] = (MpIeee( "3.0" ) * m[i] - MpIeee( "2.0" ) * b[i] - tL_ip1) / h_i;
          d[i] = (b[i] + tL_ip1 - MpIeee( "2.0" ) * m[i]) / (h_i * h_i);
        }
    }
}


static int
 akima_init(void * vstate, const MpIeee x_array[], const MpIeee y_array[],
            size_t size)
{
  akima_state_t *state = (akima_state_t *) vstate;

  MpIeee * m=  &state->_m[2];  /* //state->_m + 2;  offset so we can address the -1,-2
                                 components */

  size_t i;
  for (i = 0; i <= size - 2; i++)
    {
      m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);
    }
  
  /* non-periodic boundary conditions */
  m[-2] = MpIeee( "3.0" ) * m[0] - MpIeee( "2.0" ) * m[1];
  m[-1] = MpIeee( "2.0" ) * m[0] - m[1];
  m[size - 1] = MpIeee( "2.0" ) * m[size - 2] - m[size - 3];
  m[size] = MpIeee( "3.0" ) * m[size - 2] - MpIeee( "2.0" ) * m[size - 3];
  
  akima_calc (x_array, state->b, state->c, state->d, size, m);
  
  return GSL_SUCCESS;
}


static int
 akima_init_periodic(void * vstate,
                     const MpIeee x_array[],
                     const MpIeee y_array[],
                     size_t size)
{
  akima_state_t *state = (akima_state_t *) vstate;
  
  MpIeee * m=  &state->_m[2];/*state->_m + MpIeee( "2" ); /* offset so we can address the -1,-2
                                 components */

  size_t i;
  for (i = 0; i <= size - 2; i++)
    {
      m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);
    }
  
  /* periodic boundary conditions */
  m[-2] = m[size - 1 - 2];
  m[-1] = m[size - 1 - 1];
  m[size - 1] = m[0];
  m[size] = m[1];
  
  akima_calc (x_array, state->b, state->c, state->d, size, m);

  return GSL_SUCCESS;
}

static void
akima_free (void * vstate)
{
  akima_state_t *state = (akima_state_t *) vstate;

  free (state->b);
  free (state->c);
  free (state->d);
  free (state->_m);
  free (state);
}


static
int
 akima_eval(const void * vstate,
            const MpIeee x_array[], const MpIeee y_array[], size_t size,
            MpIeee x,
            gsl_interp_accel * a,
            MpIeee *y)
{
  const akima_state_t *state = (const akima_state_t *) vstate;

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
  {
    const MpIeee x_lo=  x_array[index];
    const MpIeee delx=  x - x_lo;
    const MpIeee b=  state->b[index];
    const MpIeee c=  state->c[index];
    const MpIeee d=  state->d[index];
    *y = y_array[index] + delx * (b + delx * (c + d * delx));
    return GSL_SUCCESS;
  }
}


static int
 akima_eval_deriv(const void * vstate,
                  const MpIeee x_array[], const MpIeee y_array[], size_t size,
                  MpIeee x,
                  gsl_interp_accel * a,
                  MpIeee *dydx)
{
  const akima_state_t *state = (const akima_state_t *) vstate;

  size_t index;

  DISCARD_POINTER(y_array); /* prevent warning about unused parameter */
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  {
    MpIeee x_lo=  x_array[index];
    MpIeee delx=  x - x_lo;
    MpIeee b=  state->b[index];
    MpIeee c=  state->c[index];
    MpIeee d=  state->d[index];
    *dydx = b + delx * (MpIeee( "2.0" ) * c + MpIeee( "3.0" ) * d * delx);
    return GSL_SUCCESS;
  }
}


static
int
 akima_eval_deriv2(const void * vstate,
                   const MpIeee x_array[], const MpIeee y_array[], size_t size,
                   MpIeee x,
                   gsl_interp_accel * a,
                   MpIeee *y_pp)
{
  const akima_state_t *state = (const akima_state_t *) vstate;

  size_t index;

  DISCARD_POINTER(y_array); /* prevent warning about unused parameter */

  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  {
    const MpIeee x_lo=  x_array[index];
    const MpIeee delx=  x - x_lo;
    const MpIeee c=  state->c[index];
    const MpIeee d=  state->d[index];
    *y_pp = MpIeee( "2.0" ) * c + MpIeee( "6.0" ) * d * delx;
    return GSL_SUCCESS;
  }
}


static
int
 akima_eval_integ(const void * vstate,
                  const MpIeee x_array[], const MpIeee y_array[], size_t size,
                  gsl_interp_accel * acc,
                  MpIeee a, MpIeee b,
                  MpIeee * result)
{
  const akima_state_t *state = (const akima_state_t *) vstate;

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
  
  *result = MpIeee( "0.0" );

  /* interior intervals */
  
  for(i=index_a; i<=index_b; i++) {
    const MpIeee x_hi=  x_array[i + 1];
    const MpIeee x_lo=  x_array[i];
    const MpIeee y_lo=  y_array[i];
    const MpIeee dx=  x_hi - x_lo;
    if(dx != 0.0) {

      if (i == index_a || i == index_b)
        {
          MpIeee x1=  (i == index_a) ? a : x_lo;
          MpIeee x2=  (i == index_b) ? b : x_hi;
          *result += integ_eval (y_lo, state->b[i], state->c[i], state->d[i],
                                 x_lo, x1, x2);
        }
      else
        {
          *result += dx * (y_lo 
                           + dx*(MpIeee( "0.5" )*state->b[i] 
                                 + dx*(state->c[i]/MpIeee( "3.0" ) 
                                       + MpIeee( "0.25" )*state->d[i]*dx)));
        }
    }
    else {
      *result = MpIeee( "0.0" );
      return GSL_FAILURE;
    }
  }
  
  return GSL_SUCCESS;
}


static const gsl_interp_type akima_type = 
{
  "akima", 
  5,
  &akima_alloc,
  &akima_init,
  &akima_eval,
  &akima_eval_deriv,
  &akima_eval_deriv2,
  &akima_eval_integ,
  &akima_free
};

const gsl_interp_type * gsl_interp_akima = &akima_type;

static const gsl_interp_type akima_periodic_type = 
{
  "akima-periodic", 
  5,
  &akima_alloc,
  &akima_init_periodic,
  &akima_eval,
  &akima_eval_deriv,
  &akima_eval_deriv2,
  &akima_eval_integ,
  &akima_free
};

const gsl_interp_type * gsl_interp_akima_periodic = &akima_periodic_type;
