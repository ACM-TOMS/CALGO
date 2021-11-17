#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/cspline.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2004 Gerard Jungman
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include "integ_eval.h"
#include <gsl/gsl_interp.h>

typedef struct
{
  MpIeee * c;
  MpIeee * g;
  MpIeee * diag;
  MpIeee * offdiag;
} cspline_state_t;


/* common initialization */
static void *
cspline_alloc (size_t size)
{
  cspline_state_t * state = (cspline_state_t *) malloc (sizeof (cspline_state_t));

  if (state == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for state", GSL_ENOMEM);
    }
  
  state->c = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->c == NULL)
    {
      free (state);
      GSL_ERROR_NULL("failed to allocate space for c", GSL_ENOMEM);
    }

  state->g = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->g == NULL)
    {
      free (state->c);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for g", GSL_ENOMEM);
    }

  state->diag = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->diag == NULL)
    {
      free (state->g);
      free (state->c);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for diag", GSL_ENOMEM);
    }

  state->offdiag = (MpIeee*) malloc (size * sizeof (MpIeee));
  
  if (state->offdiag == NULL)
    {
      free (state->diag);
      free (state->g);
      free (state->c);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for offdiag", GSL_ENOMEM);
    }

  return state;
}


/* natural spline calculation
 * see [Engeln-Mullges + Uhlig, p. 254]
 */
static int
 cspline_init(void * vstate, const MpIeee xa[], const MpIeee ya[],
              size_t size)
{
  cspline_state_t *state = (cspline_state_t *) vstate;

  size_t i;
  size_t num_points = size;
  size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index - 1;    /* linear system is sys_size x sys_size */

  state->c[0] = 0.0;
  state->c[max_index] = 0.0;

  for (i = 0; i < sys_size; i++)
    {
      const MpIeee h_i=  xa[i + 1] - xa[i];
      const MpIeee h_ip1=  xa[i + 2] - xa[i + 1];
      const MpIeee ydiff_i=  ya[i + 1] - ya[i];
      const MpIeee ydiff_ip1=  ya[i + 2] - ya[i + 1];
      
      MpIeee g_i; //=  (h_i != 0.0) ? 1.0 / h_i : 0.0;
      if( h_i != MpIeee("0") ) g_i=MpIeee("1")/h_i;
      else g_i=MpIeee("0");
      
      MpIeee g_ip1;  //  =(h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      if( h_ip1 != MpIeee("0") ) g_ip1 = MpIeee("1") / h_ip1;
      else g_ip1 = MpIeee("0" );

      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i);
    }

  if (sys_size == 1)
    {
      state->c[1] = state->g[0] / state->diag[0];
      return GSL_SUCCESS;
    }
  else
    {
      gsl_vector_view g_vec = gsl_vector_view_array(state->g, sys_size);
      gsl_vector_view diag_vec = gsl_vector_view_array(state->diag, sys_size);
      gsl_vector_view offdiag_vec = gsl_vector_view_array(state->offdiag, sys_size - 1);
      gsl_vector_view solution_vec = gsl_vector_view_array ((state->c) + 1, sys_size);
      
      int  status=  gsl_linalg_solve_symm_tridiag(&diag_vec.vector, 
                                                 &offdiag_vec.vector, 
                                                 &g_vec.vector, 
                                                 &solution_vec.vector);
      return status;
    }
}


/* periodic spline calculation
 * see [Engeln-Mullges + Uhlig, p. 256]
 */
static int
 cspline_init_periodic(void * vstate, const MpIeee xa[], const MpIeee ya[],
                       size_t size)
{
  cspline_state_t *state = (cspline_state_t *) vstate;

  size_t i;
  size_t num_points = size;
  size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index;    /* linear system is sys_size x sys_size */

  if (sys_size == 2) {
    /* solve 2x2 system */
    
    const MpIeee h0=  xa[1] - xa[0];
    const MpIeee h1=  xa[2] - xa[1];
    const MpIeee h2=  xa[3] - xa[2];
    const MpIeee A=  2.0*(h0 + h1);
    const MpIeee B=  h0 + h1;
    MpIeee g[2];
    MpIeee det;
    
    g[0] = MpIeee( "3.0" ) * ((ya[2] - ya[1]) / h1 - (ya[1] - ya[0]) / h0);
    g[1] = MpIeee( "3.0" ) * ((ya[1] - ya[2]) / h2 - (ya[2] - ya[1]) / h1);
    
    det = MpIeee( "3.0" ) * (h0 + h1) * (h0 + h1);
    state->c[1] = ( A * g[0] - B * g[1])/det;
    state->c[2] = (-B * g[0] + A * g[1])/det;
    state->c[0] = state->c[2];
    
    return GSL_SUCCESS;
  } else {
    
    for (i = 0; i < sys_size-1; i++) {
      const MpIeee h_i=  xa[i + 1] - xa[i];
      const MpIeee h_ip1=  xa[i + 2] - xa[i + 1];
      const MpIeee ydiff_i=  ya[i + 1] - ya[i];
      const MpIeee ydiff_ip1=  ya[i + 2] - ya[i + 1];
      
      //const MpIeee g_i=  (h_i != 0.0) ? 1.0 / h_i : 0.0;
      //const MpIeee g_ip1=  (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;

      MpIeee g_i; //=  (h_i != 0.0) ? 1.0 / h_i : 0.0;
      if( h_i != MpIeee("0") ) g_i=MpIeee("1")/h_i;
      else g_i=MpIeee("0");
      
      MpIeee g_ip1;  //  =(h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      if( h_ip1 != MpIeee("0") ) g_ip1 = MpIeee("1") / h_ip1;
      else g_ip1 = MpIeee("0" );
      
      
      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i);
    }

    i = sys_size - 1;

    {
      const MpIeee h_i=  xa[i + 1] - xa[i];
      const MpIeee h_ip1=  xa[1] - xa[0];
      const MpIeee ydiff_i=  ya[i + 1] - ya[i];
      const MpIeee ydiff_ip1=  ya[1] - ya[0];
//      const MpIeee g_i=  (h_i != 0.0) ? 1.0 / h_i : 0.0;
//      const MpIeee g_ip1=  (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;

      MpIeee g_i; //=  (h_i != 0.0) ? 1.0 / h_i : 0.0;
      if( h_i != MpIeee("0") ) g_i=MpIeee("1")/h_i;
      else g_i=MpIeee("0");
      
      MpIeee g_ip1;  //  =(h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      if( h_ip1 != MpIeee("0") ) g_ip1 = MpIeee("1") / h_ip1;
      else g_ip1 = MpIeee("0" );


      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i);
    }
    
    {
      gsl_vector_view g_vec = gsl_vector_view_array(state->g, sys_size);
      gsl_vector_view diag_vec = gsl_vector_view_array(state->diag, sys_size);
      gsl_vector_view offdiag_vec = gsl_vector_view_array(state->offdiag, sys_size);
      gsl_vector_view solution_vec = gsl_vector_view_array ((state->c) + 1, sys_size);
      
      int  status=  gsl_linalg_solve_symm_cyc_tridiag(&diag_vec.vector, 
                                                     &offdiag_vec.vector, 
                                                     &g_vec.vector, 
                                                     &solution_vec.vector);
      state->c[0] = state->c[max_index];
      
      return status;
    }
  }
}


static
void
cspline_free (void * vstate)
{
  cspline_state_t *state = (cspline_state_t *) vstate;
  
  free (state->c);
  free (state->g);
  free (state->diag);
  free (state->offdiag);
  free (state);
}

/* function for common coefficient determination
 */
static inline void
coeff_calc (const MpIeee c_array[], MpIeee dy, MpIeee dx, size_t index,  
            MpIeee * b, MpIeee * c, MpIeee * d)
{
  const MpIeee c_i=  c_array[index];
  const MpIeee c_ip1=  c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + MpIeee( "2.0" ) * c_i) / MpIeee( "3.0" );
  *c = c_i;
  *d = (c_ip1 - c_i) / (MpIeee( "3.0" ) * dx);
}


static
int
 cspline_eval(const void * vstate,
              const MpIeee x_array[], const MpIeee y_array[], size_t size,
              MpIeee x,
              gsl_interp_accel * a,
              MpIeee *y)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  MpIeee x_lo;MpIeee  x_hi;
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
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > MpIeee( "0.0" ))
    {
      const MpIeee y_lo=  y_array[index];
      const MpIeee y_hi=  y_array[index + 1];
      const MpIeee dy=  y_hi - y_lo;
      MpIeee delx=  x - x_lo;
      MpIeee b_i;MpIeee  c_i;MpIeee  d_i; 
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *y = y_lo + delx * (b_i + delx * (c_i + delx * d_i));
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
 cspline_eval_deriv(const void * vstate,
                    const MpIeee x_array[], const MpIeee y_array[], size_t size,
                    MpIeee x,
                    gsl_interp_accel * a,
                    MpIeee *dydx)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  MpIeee x_lo;MpIeee  x_hi;
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
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > MpIeee( "0.0" ))
    {
      const MpIeee y_lo=  y_array[index];
      const MpIeee y_hi=  y_array[index + 1];
      const MpIeee dy=  y_hi - y_lo;
      MpIeee delx=  x - x_lo;
      MpIeee b_i;MpIeee  c_i;MpIeee  d_i; 
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *dydx = b_i + delx * (MpIeee( "2.0" ) * c_i + MpIeee( "3.0" ) * d_i * delx);
      return GSL_SUCCESS;
    }
  else
    {
      *dydx = MpIeee( "0.0" );
      return GSL_FAILURE;
    }
}


static
int
 cspline_eval_deriv2(const void * vstate,
                     const MpIeee x_array[], const MpIeee y_array[], size_t size,
                     MpIeee x,
                     gsl_interp_accel * a,
                     MpIeee * y_pp)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  MpIeee x_lo;MpIeee  x_hi;
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
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > MpIeee( "0.0" ))
    {
      const MpIeee y_lo=  y_array[index];
      const MpIeee y_hi=  y_array[index + 1];
      const MpIeee dy=  y_hi - y_lo;
      MpIeee delx=  x - x_lo;
      MpIeee b_i;MpIeee  c_i;MpIeee  d_i;
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *y_pp = MpIeee( "2.0" ) * c_i + MpIeee( "6.0" ) * d_i * delx;
      return GSL_SUCCESS;
    }
  else
    {
      *y_pp = MpIeee( "0.0" );
      return GSL_FAILURE;
    }
}


static
int
 cspline_eval_integ(const void * vstate,
                    const MpIeee x_array[], const MpIeee y_array[], size_t size,
                    gsl_interp_accel * acc,
                    MpIeee a, MpIeee b,
                    MpIeee * result)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

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
    const MpIeee y_hi=  y_array[i + 1];
    const MpIeee dx=  x_hi - x_lo;
    const MpIeee dy=  y_hi - y_lo;
    if(dx != 0.0) {
      MpIeee b_i;MpIeee  c_i;MpIeee  d_i; 
      coeff_calc(state->c, dy, dx, i,  &b_i, &c_i, &d_i);
      
      if (i == index_a || i == index_b)
        {
          MpIeee x1=  (i == index_a) ? a : x_lo;
          MpIeee x2=  (i == index_b) ? b : x_hi;
          *result += integ_eval(y_lo, b_i, c_i, d_i, x_lo, x1, x2);
        }
      else
        {
          *result += dx * (y_lo + dx*(MpIeee( "0.5" )*b_i + dx*(c_i/MpIeee( "3.0" ) + MpIeee( "0.25" )*d_i*dx)));
        }
    }
    else {
      *result = MpIeee( "0.0" );
      return GSL_FAILURE;
    }
  }
  
  return GSL_SUCCESS;
}

static const gsl_interp_type cspline_type = 
{
  "cspline", 
  3,
  &cspline_alloc,
  &cspline_init,
  &cspline_eval,
  &cspline_eval_deriv,
  &cspline_eval_deriv2,
  &cspline_eval_integ,
  &cspline_free
};

const gsl_interp_type * gsl_interp_cspline = &cspline_type;

static const gsl_interp_type cspline_periodic_type = 
{
  "cspline-periodic", 
  2,
  &cspline_alloc,
  &cspline_init_periodic,
  &cspline_eval,
  &cspline_eval_deriv,
  &cspline_eval_deriv2,
  &cspline_eval_integ,
  &cspline_free
};

const gsl_interp_type * gsl_interp_cspline_periodic = &cspline_periodic_type;


