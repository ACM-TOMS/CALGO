#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multimin/simplex.c
 * 
 * Copyright (C) 2002 Tuomo Keskitalo, Ivo Alxneit
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

/*
   - Originally written by Tuomo Keskitalo <tuomo.keskitalo@iki.fi>
   - Corrections to nmsimplex_iterate and other functions 
     by Ivo Alxneit <ivo.alxneit@psi.ch>
   - Additional help by Brian Gough <bjg@network-theory.co.uk>
*/

/* The Simplex method of Nelder and Mead,
   also known as the polytope search alogorithm. Ref:
   Nelder, J.A., Mead, R., Computer Journal 7 (1965) pp. 308-313.

   This implementation uses n+1 corner points in the simplex.
*/

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

typedef struct
{
  gsl_matrix *x1;               /* simplex corner points */
  gsl_vector *y1;               /* function value at corner points */
  gsl_vector *ws1;              /* workspace 1 for algorithm */
  gsl_vector *ws2;              /* workspace 2 for algorithm */
}
nmsimplex_state_t;

static MpIeee nmsimplex_move_corner(const MpIeee coeff, const nmsimplex_state_t * state,
                       size_t corner, gsl_vector * xc,
                       const gsl_multimin_function * f)
{
  /* moves a simplex corner scaled by coeff (negative value represents 
     mirroring by the middle point of the "other" corner points)
     and gives new corner in xc and function value at xc as a 
     return value 
   */

  gsl_matrix *x1 = state->x1;

  size_t i, j;
  MpIeee newval;MpIeee  mp;

  if (x1->size1 < 2)
    {
      GSL_ERROR ("simplex cannot have less than two corners!", GSL_EFAILED);
    }

  for (j = 0; j < x1->size2; j++)
    {
      mp = MpIeee( "0.0" );
      for (i = 0; i < x1->size1; i++)
        {
          if (i != corner)
            {
              mp += (gsl_matrix_get (x1, i, j));
            }
        }
      mp /= (MpIeee) (x1->size1 - MpIeee( "1" ));
      newval = mp - coeff * (mp - gsl_matrix_get (x1, corner, j));
      gsl_vector_set (xc, j, newval);
    }

  newval = GSL_MULTIMIN_FN_EVAL (f, xc);

  return newval;
}

static int
 nmsimplex_contract_by_best(nmsimplex_state_t * state, size_t best,
                            gsl_vector * xc, gsl_multimin_function * f)
{

  /* Function contracts the simplex in respect to 
     best valued corner. That is, all corners besides the 
     best corner are moved. */

  /* the xc vector is simply work space here */

  gsl_matrix *x1 = state->x1;
  gsl_vector *y1 = state->y1;

  size_t i, j;
  MpIeee newval;

  for (i = 0; i < x1->size1; i++)
    {
      if (i != best)
        {
          for (j = 0; j < x1->size2; j++)
            {
              newval = MpIeee( "0.5" ) * (gsl_matrix_get (x1, i, j)
                              + gsl_matrix_get (x1, best, j));
              gsl_matrix_set (x1, i, j, newval);
            }

          /* evaluate function in the new point */

          gsl_matrix_get_row (xc, x1, i);
          newval = GSL_MULTIMIN_FN_EVAL (f, xc);
          gsl_vector_set (y1, i, newval);
        }
    }

  return GSL_SUCCESS;
}

static int
 nmsimplex_calc_center(const nmsimplex_state_t * state, gsl_vector * mp)
{
  /* calculates the center of the simplex to mp */

  gsl_matrix *x1 = state->x1;

  size_t i, j;
  MpIeee val;

  for (j = 0; j < x1->size2; j++)
    {
      val = MpIeee( "0.0" );
      for (i = 0; i < x1->size1; i++)
        {
          val += gsl_matrix_get (x1, i, j);
        }
      val /= x1->size1;
      gsl_vector_set (mp, j, val);
    }

  return GSL_SUCCESS;
}

static MpIeee nmsimplex_size(nmsimplex_state_t * state)
{
  /* calculates simplex size as average sum of length of vectors 
     from simplex center to corner points:     

     ( sum ( || y - y_middlepoint || ) ) / n 
   */

  gsl_vector *s = state->ws1;
  gsl_vector *mp = state->ws2;

  gsl_matrix *x1 = state->x1;
  size_t i;

  MpIeee ss=  MpIeee( "0.0" );

  /* Calculate middle point */
  nmsimplex_calc_center (state, mp);

  for (i = 0; i < x1->size1; i++)
    {
      gsl_matrix_get_row (s, x1, i);
      gsl_blas_daxpy (-1.0, mp, s);
      ss += gsl_blas_dnrm2 (s);
    }

  return ss / (MpIeee) (x1->size1);
}

static int
 nmsimplex_alloc(void *vstate, size_t n)
{
  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  state->x1 = gsl_matrix_alloc (n + 1, n);

  if (state->x1 == NULL)
    {
      GSL_ERROR ("failed to allocate space for x1", GSL_ENOMEM);
    }

  state->y1 = gsl_vector_alloc (n + 1);

  if (state->y1 == NULL)
    {
      GSL_ERROR ("failed to allocate space for y", GSL_ENOMEM);
    }

  state->ws1 = gsl_vector_alloc (n);

  if (state->ws1 == NULL)
    {
      GSL_ERROR ("failed to allocate space for ws1", GSL_ENOMEM);
    }

  state->ws2 = gsl_vector_alloc (n);

  if (state->ws2 == NULL)
    {
      GSL_ERROR ("failed to allocate space for ws2", GSL_ENOMEM);
    }

  return GSL_SUCCESS;
}

static int
 nmsimplex_set(void *vstate, gsl_multimin_function * f,
               const gsl_vector * x,
               MpIeee *size, const gsl_vector * step_size)
{
  int  status;
  size_t i;
  MpIeee val;

  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  gsl_vector *xtemp = state->ws1;

  /* first point is the original x0 */

  val = GSL_MULTIMIN_FN_EVAL (f, x);
  gsl_matrix_set_row (state->x1, 0, x);
  gsl_vector_set (state->y1, 0, val);

  /* following points are initialized to x0 + step_size */

  for (i = 0; i < x->size; i++)
    {
      status = gsl_vector_memcpy (xtemp, x);

      if (status != 0)
        {
          GSL_ERROR ("vector memcopy failed", GSL_EFAILED);
        }

      val = gsl_vector_get (xtemp, i) + gsl_vector_get (step_size, i);
      gsl_vector_set (xtemp, i, val);
      val = GSL_MULTIMIN_FN_EVAL (f, xtemp);
      gsl_matrix_set_row (state->x1, i + 1, xtemp);
      gsl_vector_set (state->y1, i + 1, val);
    }

  /* Initialize simplex size */

  *size = nmsimplex_size (state);

  return GSL_SUCCESS;
}

static void
nmsimplex_free (void *vstate)
{
  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  gsl_matrix_free (state->x1);
  gsl_vector_free (state->y1);
  gsl_vector_free (state->ws1);
  gsl_vector_free (state->ws2);
}

static int
 nmsimplex_iterate(void *vstate, gsl_multimin_function * f,
                   gsl_vector * x, MpIeee *size, MpIeee *fval)
{

  /* Simplex iteration tries to minimize function f value */
  /* Includes corrections from Ivo Alxneit <ivo.alxneit@psi.ch> */

  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  /* xc and xc2 vectors store tried corner point coordinates */

  gsl_vector *xc = state->ws1;
  gsl_vector *xc2 = state->ws2;
  gsl_vector *y1 = state->y1;
  gsl_matrix *x1 = state->x1;

  size_t n = y1->size;
  size_t i;
  size_t hi = 0, s_hi = 0, lo = 0;
  MpIeee dhi;MpIeee  ds_hi;MpIeee  dlo;
  int  status;
  MpIeee val;MpIeee  val2;

  /* get index of highest, second highest and lowest point */

  dhi = ds_hi = dlo = gsl_vector_get (y1, 0);

  for (i = 1; i < n; i++)
    {
      val = (gsl_vector_get (y1, i));
      if (val < dlo)
        {
          dlo = val;
          lo = i;
        }
      else if (val > dhi)
        {
          ds_hi = dhi;
          s_hi = hi;
          dhi = val;
          hi = i;
        }
      else if (val > ds_hi)
        {
          ds_hi = val;
          s_hi = i;
        }
    }

  /* reflect the highest value */

  val = nmsimplex_move_corner (-MpIeee( "1.0" ), state, hi, xc, f);

  if (val < gsl_vector_get (y1, lo))
    {

      /* reflected point becomes lowest point, try expansion */

      val2 = nmsimplex_move_corner (-MpIeee( "2.0" ), state, hi, xc2, f);

      if (val2 < gsl_vector_get (y1, lo))
        {
          gsl_matrix_set_row (x1, hi, xc2);
          gsl_vector_set (y1, hi, val2);
        }
      else
        {
          gsl_matrix_set_row (x1, hi, xc);
          gsl_vector_set (y1, hi, val);
        }
    }

  /* reflection does not improve things enough */

  else if (val > gsl_vector_get (y1, s_hi))
    {
      if (val <= gsl_vector_get (y1, hi))
        {

          /* if trial point is better than highest point, replace 
             highest point */

          gsl_matrix_set_row (x1, hi, xc);
          gsl_vector_set (y1, hi, val);
        }

      /* try one dimensional contraction */

      val2 = nmsimplex_move_corner (MpIeee( "0.5" ), state, hi, xc2, f);

      if (val2 <= gsl_vector_get (y1, hi))
        {
          gsl_matrix_set_row (state->x1, hi, xc2);
          gsl_vector_set (y1, hi, val2);
        }

      else
        {

          /* contract the whole simplex in respect to the best point */

          status = nmsimplex_contract_by_best (state, lo, xc, f);
          if (status != 0)
            {
              GSL_ERROR ("nmsimplex_contract_by_best failed", GSL_EFAILED);
            }
        }
    }
  else
    {

      /* trial point is better than second highest point. 
         Replace highest point by it */

      gsl_matrix_set_row (x1, hi, xc);
      gsl_vector_set (y1, hi, val);
    }

  /* return lowest point of simplex as x */

  lo = gsl_vector_min_index (y1);
  gsl_matrix_get_row (x, x1, lo);
  *fval = gsl_vector_get (y1, lo);

  /* Update simplex size */

  *size = nmsimplex_size (state);

  return GSL_SUCCESS;
}

static const gsl_multimin_fminimizer_type nmsimplex_type = 
{ "nmsimplex",  /* name */
  sizeof (nmsimplex_state_t),
  &nmsimplex_alloc,
  &nmsimplex_set,
  &nmsimplex_iterate,
  &nmsimplex_free
};

const gsl_multimin_fminimizer_type
  * gsl_multimin_fminimizer_nmsimplex = &nmsimplex_type;
