#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ode-initval/rk4imp.c
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

/* Runge-Kutta 4, Gaussian implicit */

/* Author:  G. Jungman
*/

/* Error estimation by step doubling, see eg. Ascher, U.M., Petzold,
   L.R., Computer methods for ordinary differential and
   differential-algebraic equations, SIAM, Philadelphia, 1998.
   Method coefficients can also be found in it.
*/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

typedef struct
{
  MpIeee *k1nu;
  MpIeee *k2nu;
  MpIeee *ytmp1;
  MpIeee *ytmp2;
  MpIeee *y0;
  MpIeee *y0_orig;
  MpIeee *y_onestep;
}
rk4imp_state_t;

static void *
rk4imp_alloc (size_t dim)
{
  rk4imp_state_t *state = (rk4imp_state_t *) malloc (sizeof (rk4imp_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for rk4imp_state",
                      GSL_ENOMEM);
    }

  state->k1nu = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->k1nu == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for k1nu", GSL_ENOMEM);
    }

  state->k2nu = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->k2nu == 0)
    {
      free (state->k1nu);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for k2nu", GSL_ENOMEM);
    }

  state->ytmp1 = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->ytmp1 == 0)
    {
      free (state->k2nu);
      free (state->k1nu);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ytmp1", GSL_ENOMEM);
    }

  state->ytmp2 = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->ytmp2 == 0)
    {
      free (state->ytmp1);
      free (state->k2nu);
      free (state->k1nu);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ytmp2", GSL_ENOMEM);
    }

  state->y0 = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->y0 == 0)
    {
      free (state->ytmp2);
      free (state->ytmp1);
      free (state->k2nu);
      free (state->k1nu);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for y0", GSL_ENOMEM);
    }

  state->y0_orig = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->y0_orig == 0)
    {
      free (state->y0);
      free (state->ytmp2);
      free (state->ytmp1);
      free (state->k2nu);
      free (state->k1nu);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for y0_orig", GSL_ENOMEM);
    }

  state->y_onestep = (MpIeee*) malloc (dim * sizeof (MpIeee));

  if (state->y_onestep == 0)
    {
      free (state->y0_orig);
      free (state->y0);
      free (state->ytmp2);
      free (state->ytmp1);
      free (state->k2nu);
      free (state->k1nu);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for y_onestep", GSL_ENOMEM);
    }

  return state;
}

static int
 rk4imp_step(MpIeee *y, rk4imp_state_t *state, 
	     const MpIeee h, const MpIeee t, 
	     const size_t dim, const gsl_odeiv_system *sys)
{
  /* Makes a Runge-Kutta 4th order implicit advance with step size h.
     y0 is initial values of variables y. 

     The implicit matrix equations to solve are:

     Y1 = y0 + h * a11 * f(t + h * c1, Y1) + h * a12 * f(t + h * c2, Y2) 
     Y2 = y0 + h * a21 * f(t + h * c1, Y1) + h * a22 * f(t + h * c2, Y2) 

     y = y0 + h * b1 * f(t + h * c1, Y1) + h * b2 * f(t + h * c2, Y2)

     with constant coefficients a, b and c. For this method
     they are: b=[0.5 0.5] c=[(3-sqrt(3))/6 (3+sqrt(3))/6]
     a11=1/4, a12=(3-2*sqrt(3))/12, a21=(3+2*sqrt(3))/12 and a22=1/4
  */

  const MpIeee ir3=  1.0 / M_SQRT3;
  const int iter_steps = 3;
  int  nu;
  size_t i;

  MpIeee *const k1nu=  state->k1nu;
  MpIeee *const k2nu=  state->k2nu;
  MpIeee *const ytmp1=  state->ytmp1;
  MpIeee *const ytmp2=  state->ytmp2;

  /* iterative solution of Y1 and Y2.

     Note: This method does not check for convergence of the
     iterative solution! 
  */

  for (nu = 0; nu < iter_steps; nu++)
    {
      for (i = 0; i < dim; i++)
        {
          ytmp1[i] =
            y[i] + h * (MpIeee( "0.25" ) * k1nu[i] + MpIeee( "0.5" ) * (MpIeee( "0.5" ) - ir3) * k2nu[i]);
          ytmp2[i] =
            y[i] + h * (MpIeee( "0.25" ) * k2nu[i] + MpIeee( "0.5" ) * (MpIeee( "0.5" ) + ir3) * k1nu[i]);
        }
      {
        int  s= 
	  GSL_ODEIV_FN_EVAL (sys, t + 0.5 * h * (1.0 - ir3), ytmp1, k1nu);
	
	if (s != GSL_SUCCESS)
	  {
	    return GSL_EBADFUNC;
	  }    
      }
      {
        int  s= 
	  GSL_ODEIV_FN_EVAL (sys, t + 0.5 * h * (1.0 + ir3), ytmp2, k2nu);
	
	if (s != GSL_SUCCESS)
	  {
	    return GSL_EBADFUNC;
	  }    
      }
    }

  /* assignment */
  
  for (i = 0; i < dim; i++)
    {
      const MpIeee d_i=  0.5 * (k1nu[i] + k2nu[i]);
      y[i] += h * d_i;
    }

  return GSL_SUCCESS;
}

static int
 rk4imp_apply(void *vstate,
              size_t dim,
              MpIeee t,
              MpIeee h,
              MpIeee y[],
              MpIeee yerr[],
              const MpIeee dydt_in[],
              MpIeee dydt_out[], 
              const gsl_odeiv_system * sys)
{
  rk4imp_state_t *state = (rk4imp_state_t *) vstate;

  size_t i;

  MpIeee *y0=  state->y0;
  MpIeee *y0_orig=  state->y0_orig;
  MpIeee *y_onestep=  state->y_onestep;
  MpIeee *k1nu=  state->k1nu;
  MpIeee *k2nu=  state->k2nu;

  /* Initialization step */
  DBL_MEMCPY (y0, y, dim);

  /* Save initial values in case of failure */
  DBL_MEMCPY (y0_orig, y, dim);

  if (dydt_in != 0)
    {
      DBL_MEMCPY (k1nu, dydt_in, dim);
    }
  else
    {
      int  s=  GSL_ODEIV_FN_EVAL (sys, t, y, k1nu);

      if (s != GSL_SUCCESS)
        {
          return GSL_EBADFUNC;
        }
    }

  DBL_MEMCPY (k2nu, k1nu, dim);

  /* First traverse h with one step (save to y_onestep) */

  DBL_MEMCPY (y_onestep, y, dim);

  {
    int  s=  rk4imp_step (y_onestep, state, h, t, dim, sys);

    if (s != GSL_SUCCESS) 
      {
	return GSL_EBADFUNC;
      }
  }
  
 /* Then with two steps with half step length (save to y) */ 
  
  {
    int  s=  rk4imp_step (y, state, h/2.0, t, dim, sys);

    if (s != GSL_SUCCESS) 
      {
	/* Restore original y vector */
	DBL_MEMCPY (y, y0_orig, dim);
	return GSL_EBADFUNC;
      }
  }

  DBL_MEMCPY (y0, y, dim);

  {
    int  s=  GSL_ODEIV_FN_EVAL (sys, t + h/2.0, y, k1nu);

    if (s != GSL_SUCCESS)
      {
	/* Restore original y vector */
	DBL_MEMCPY (y, y0_orig, dim);
	return GSL_EBADFUNC;
      }
  }

  DBL_MEMCPY (k2nu, k1nu, dim);
  
  {
    int  s=  rk4imp_step (y, state, h/2.0, t + h/2.0, dim, sys);

    if (s != GSL_SUCCESS) 
      {
	/* Restore original y vector */
	DBL_MEMCPY (y, y0_orig, dim);
	return GSL_EBADFUNC;
      }
  }
  
  /* Derivatives at output */
  
  if (dydt_out != NULL) 
    {
      int  s=  GSL_ODEIV_FN_EVAL (sys, t + h, y, dydt_out);
      
      if (s != GSL_SUCCESS) {
	/* Restore original y vector */
	DBL_MEMCPY (y, y0_orig, dim);
	return GSL_EBADFUNC;
      } 
    }
    
  /* Error estimation */

  /* Denominator in step doubling error equation 
   *  yerr = 0.5 * | y(onestep) - y(twosteps) | / (2^order - 1)  
   */

  for (i = 0; i < dim; i++) 
    {
      yerr[i] = MpIeee( "8.0" ) * MpIeee( "0.5" ) * (y[i] - y_onestep[i]) / MpIeee( "15.0" );
    }
  
  return GSL_SUCCESS;
}

static int
 rk4imp_reset(void *vstate, size_t dim)
{
  rk4imp_state_t *state = (rk4imp_state_t *) vstate;

  DBL_ZERO_MEMSET (state->y_onestep, dim);
  DBL_ZERO_MEMSET (state->y0_orig, dim);
  DBL_ZERO_MEMSET (state->y0, dim);
  DBL_ZERO_MEMSET (state->k1nu, dim);
  DBL_ZERO_MEMSET (state->k2nu, dim);
  DBL_ZERO_MEMSET (state->ytmp1, dim);
  DBL_ZERO_MEMSET (state->ytmp2, dim);

  return GSL_SUCCESS;
}

static unsigned int
 rk4imp_order(void *vstate)
{
  rk4imp_state_t *state = (rk4imp_state_t *) vstate;
  state = 0; /* prevent warnings about unused parameters */
  return 4;
}

static void
rk4imp_free (void *vstate)
{
  rk4imp_state_t *state = (rk4imp_state_t *) vstate;

  free (state->y_onestep);
  free (state->y0_orig);
  free (state->y0);
  free (state->k1nu);
  free (state->k2nu);
  free (state->ytmp1);
  free (state->ytmp2);
  free (state);
}

static const gsl_odeiv_step_type rk4imp_type = { "rk4imp",      /* name */
  1,                             /* can use dydt_in? */
  1,                             /* gives exact dydt_out? */
  &rk4imp_alloc,
  &rk4imp_apply,
  &rk4imp_reset,
  &rk4imp_order,
  &rk4imp_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_rk4imp = &rk4imp_type;
