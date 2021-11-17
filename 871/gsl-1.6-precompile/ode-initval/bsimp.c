#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ode-initval/bsimp.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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

/* Bulirsch-Stoer Implicit */

/* Author:  G. Jungman
 */
/* Bader-Deuflhard implicit extrapolative stepper.
 * [Numer. Math., 41, 373 (1983)]
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

#define SEQUENCE_COUNT 8
#define SEQUENCE_MAX   7

/* Bader-Deuflhard extrapolation sequence */
static const int bd_sequence[SEQUENCE_COUNT] =
  { 2, 6, 10, 14, 22, 34, 50, 70 };

typedef struct
{
  gsl_matrix *d;                /* workspace for extrapolation         */
  gsl_matrix *a_mat;            /* workspace for linear system matrix  */
  gsl_permutation *p_vec;       /* workspace for LU permutation        */

  MpIeee x[SEQUENCE_MAX];       /* workspace for extrapolation */

  /* state info */
  size_t k_current;
  size_t k_choice;
  MpIeee h_next;
  MpIeee eps;

  /* workspace for extrapolation step */
  MpIeee *yp;
  MpIeee *y_save;
  MpIeee *yerr_save;
  MpIeee *y_extrap_save;
  MpIeee *y_extrap_sequence;
  MpIeee *extrap_work;
  MpIeee *dfdt;
  MpIeee *y_temp;
  MpIeee *delta_temp;
  MpIeee *weight;
  gsl_matrix *dfdy;

  /* workspace for the basic stepper */
  MpIeee *rhs_temp;
  MpIeee *delta;

  /* order of last step */
  size_t order;
}
bsimp_state_t;

/* Compute weighting factor */

static void
compute_weights (const MpIeee y[], MpIeee w[], size_t dim)
{
  size_t i;
  for (i = 0; i < dim; i++)
    {
      MpIeee u=  fabs(y[i]);
      w[i] = (u > MpIeee( "0.0" )) ? u : MpIeee( "1.0" );
    }
}

/* Calculate a choice for the "order" of the method, using the
 * Deuflhard criteria.  
 */

static size_t
bsimp_deuf_kchoice (MpIeee eps, size_t dimension)
{
  const MpIeee safety_f=  0.25;
  const MpIeee small_eps=  safety_f * eps;

  MpIeee a_work[SEQUENCE_COUNT];
  MpIeee alpha[SEQUENCE_MAX][SEQUENCE_MAX];

  int  i;int   k;

  a_work[0] = bd_sequence[0] + MpIeee( "1.0" );

  for (k = 0; k < SEQUENCE_MAX; k++)
    {
      a_work[k + 1] = a_work[k] + bd_sequence[k + 1];
    }

  for (i = 0; i < SEQUENCE_MAX; i++)
    {
      alpha[i][i] = MpIeee( "1.0" );
      for (k = 0; k < i; k++)
        {
          const MpIeee tmp1=  a_work[k + 1] - a_work[i + 1];
          const MpIeee tmp2=  (a_work[i + 1] - a_work[0] + 1.0) * (2 * k + 1);
          alpha[k][i] = pow (small_eps, tmp1 / tmp2);
        }
    }

  a_work[0] += dimension;

  for (k = 0; k < SEQUENCE_MAX; k++)
    {
      a_work[k + 1] = a_work[k] + bd_sequence[k + 1];
    }

  for (k = 0; k < SEQUENCE_MAX - 1; k++)
    {
      if (a_work[k + 2] > a_work[k + 1] * alpha[k][k + 1])
        break;
    }

  return k;
}

static void
poly_extrap (gsl_matrix * d,
             const MpIeee x[],
             const unsigned int  i_step,
             const MpIeee x_i,
             const MpIeee y_i[],
             MpIeee y_0[], MpIeee y_0_err[], MpIeee work[], const size_t dim)
{
  size_t j, k;

  DBL_MEMCPY (y_0_err, y_i, dim);
  DBL_MEMCPY (y_0, y_i, dim);

  if (i_step == 0)
    {
      for (j = 0; j < dim; j++)
        {
          gsl_matrix_set (d, 0, j, y_i[j]);
        }
    }
  else
    {
      DBL_MEMCPY (work, y_i, dim);

      for (k = 0; k < i_step; k++)
        {
          MpIeee delta=  MpIeee( "1.0" ) / (x[i_step - k - 1] - x_i);
          const MpIeee f1=  delta * x_i;
          const MpIeee f2=  delta * x[i_step - k - 1];

          for (j = 0; j < dim; j++)
            {
              const MpIeee q_kj=  gsl_matrix_get (d, k, j);
              gsl_matrix_set (d, k, j, y_0_err[j]);
              delta = work[j] - q_kj;
              y_0_err[j] = f1 * delta;
              work[j] = f2 * delta;
              y_0[j] += y_0_err[j];
            }
        }

      for (j = 0; j < dim; j++)
        {
          gsl_matrix_set (d, i_step, j, y_0_err[j]);
        }
    }
}

/* Basic implicit Bulirsch-Stoer step.  Divide the step h_total into
 * n_step smaller steps and do the Bader-Deuflhard semi-implicit
 * iteration.  */

static int
 bsimp_step_local(void *vstate,
                  size_t dim,
                  const MpIeee t0,
                  const MpIeee h_total,
                  const unsigned int  n_step,
                  const MpIeee y[],
                  const MpIeee yp[],
                  const MpIeee dfdt[],
                  const gsl_matrix * dfdy,
                  MpIeee y_out[], 
                  const gsl_odeiv_system * sys)
{
  bsimp_state_t *state = (bsimp_state_t *) vstate;

  gsl_matrix *const a_mat = state->a_mat;
  gsl_permutation *const p_vec = state->p_vec;

  MpIeee *const delta=  state->delta;
  MpIeee *const y_temp=  state->y_temp;
  MpIeee *const delta_temp=  state->delta_temp;
  MpIeee *const rhs_temp=  state->rhs_temp;
  MpIeee *const w=  state->weight;

  gsl_vector_view y_temp_vec = gsl_vector_view_array (y_temp, dim);
  gsl_vector_view delta_temp_vec = gsl_vector_view_array (delta_temp, dim);
  gsl_vector_view rhs_temp_vec = gsl_vector_view_array (rhs_temp, dim);

  const MpIeee h=  h_total / n_step;
  MpIeee t=  t0 + h;

  MpIeee sum;

  /* This is the factor sigma referred to in equation 3.4 of the
     paper.  A relative change in y exceeding sigma indicates a
     runaway behavior. According to the authors suitable values for
     sigma are >>1.  I have chosen a value of 100*dim. BJG */

  const MpIeee max_sum=  100.0 * dim;

  int  signum;int   status;
  size_t i, j;
  size_t n_inter;

  /* Calculate the matrix for the linear system. */
  for (i = 0; i < dim; i++)
    {
      for (j = 0; j < dim; j++)
        {
          gsl_matrix_set (a_mat, i, j, -h * gsl_matrix_get (dfdy, i, j));
        }
      gsl_matrix_set (a_mat, i, i, gsl_matrix_get (a_mat, i, i) + 1.0);
    }

  /* LU decomposition for the linear system. */

  gsl_linalg_LU_decomp (a_mat, p_vec, &signum);

  /* Compute weighting factors */

  compute_weights (y, w, dim);

  /* Initial step. */

  for (i = 0; i < dim; i++)
    {
      y_temp[i] = h * (yp[i] + h * dfdt[i]);
    }

  gsl_linalg_LU_solve (a_mat, p_vec, &y_temp_vec.vector, &delta_temp_vec.vector);

  sum = MpIeee( "0.0" );

  for (i = 0; i < dim; i++)
    {
      const MpIeee di=  delta_temp[i];
      delta[i] = di;
      y_temp[i] = y[i] + di;
      sum += fabs(di) / w[i];
    }

  if (sum > max_sum) 
    {
      return GSL_EFAILED ;
    }

  /* Intermediate steps. */

  status = GSL_ODEIV_FN_EVAL (sys, t, y_temp, y_out);

  if (status)
    {
      return GSL_EBADFUNC;
    }

  for (n_inter = 1; n_inter < n_step; n_inter++)
    {
      for (i = 0; i < dim; i++)
        {
          rhs_temp[i] = h * y_out[i] - delta[i];
        }

      gsl_linalg_LU_solve (a_mat, p_vec, &rhs_temp_vec.vector, &delta_temp_vec.vector);

      sum = MpIeee( "0.0" );

      for (i = 0; i < dim; i++)
        {
          delta[i] += MpIeee( "2.0" ) * delta_temp[i];
          y_temp[i] += delta[i];
          sum += fabs(delta[i]) / w[i];
        }

      if (sum > max_sum) 
        {
          return GSL_EFAILED ;
        }

      t += h;

      status = GSL_ODEIV_FN_EVAL (sys, t, y_temp, y_out);

      if (status)
        {
          return GSL_EBADFUNC;
        }
    }


  /* Final step. */

  for (i = 0; i < dim; i++)
    {
      rhs_temp[i] = h * y_out[i] - delta[i];
    }

  gsl_linalg_LU_solve (a_mat, p_vec, &rhs_temp_vec.vector, &delta_temp_vec.vector);

  sum = MpIeee( "0.0" );

  for (i = 0; i < dim; i++)
    {
      y_out[i] = y_temp[i] + delta_temp[i];
      sum += fabs(delta_temp[i]) / w[i];
    }

  if (sum > max_sum) 
    {
      return GSL_EFAILED ;
    }

  return GSL_SUCCESS;
}


static void *
bsimp_alloc (size_t dim)
{
  bsimp_state_t *state = (bsimp_state_t *) malloc (sizeof (bsimp_state_t));

  state->d = gsl_matrix_alloc (SEQUENCE_MAX, dim);
  state->a_mat = gsl_matrix_alloc (dim, dim);
  state->p_vec = gsl_permutation_alloc (dim);

  state->yp = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->y_save = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->yerr_save = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->y_extrap_save = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->y_extrap_sequence = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->extrap_work = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->dfdt = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->y_temp = (MpIeee*) malloc (dim * sizeof (MpIeee));
  state->delta_temp = (MpIeee*) malloc (dim * sizeof(MpIeee));
  state->weight = (MpIeee*) malloc (dim * sizeof(MpIeee));

  state->dfdy = gsl_matrix_alloc (dim, dim);

  state->rhs_temp = (MpIeee*) malloc (dim * sizeof(MpIeee));
  state->delta = (MpIeee*) malloc (dim * sizeof (MpIeee));

  {
    size_t k_choice = bsimp_deuf_kchoice (GSL_SQRT_DBL_EPSILON, dim); /*FIXME: choice of epsilon? */
    state->k_choice = k_choice;
    state->k_current = k_choice;
    state->order = 2 * k_choice;
  }

  state->h_next = -GSL_SQRT_DBL_MAX;

  return state;
}

/* Perform the basic semi-implicit extrapolation
 * step, of size h, at a Deuflhard determined order.
 */
static int
 bsimp_apply(void *vstate,
             size_t dim,
             MpIeee t,
             MpIeee h,
             MpIeee y[],
             MpIeee yerr[],
             const MpIeee dydt_in[],
             MpIeee dydt_out[], 
             const gsl_odeiv_system * sys)
{
  bsimp_state_t *state = (bsimp_state_t *) vstate;

  MpIeee *const x=  state->x;
  MpIeee *const yp=  state->yp;
  MpIeee *const y_save=  state->y_save;
  MpIeee *const yerr_save=  state->yerr_save;
  MpIeee *const y_extrap_sequence=  state->y_extrap_sequence;
  MpIeee *const y_extrap_save=  state->y_extrap_save;
  MpIeee *const extrap_work=  state->extrap_work;
  MpIeee *const dfdt=  state->dfdt;
  gsl_matrix *d = state->d;
  gsl_matrix *dfdy = state->dfdy;

  const MpIeee t_local=  t;
  size_t i, k;

  if (h + t_local == t_local)
    {
      return GSL_EUNDRFLW;      /* FIXME: error condition */
    }

  DBL_MEMCPY (y_extrap_save, y, dim);

  /* Save inputs */
  DBL_MEMCPY (y_save, y, dim);
  DBL_MEMCPY (yerr_save, yerr, dim);
  
  /* Evaluate the derivative. */
  if (dydt_in != NULL)
    {
      DBL_MEMCPY (yp, dydt_in, dim);
    }
  else
    {
      int  s=  GSL_ODEIV_FN_EVAL (sys, t_local, y, yp);

      if (s != GSL_SUCCESS)
	{
          return GSL_EBADFUNC;
	}
    }

  /* Evaluate the Jacobian for the system. */
  {
    int  s=  GSL_ODEIV_JA_EVAL (sys, t_local, y, dfdy->data, dfdt);
  
    if (s != GSL_SUCCESS)
      {
        return GSL_EBADFUNC;
      }
  }

  /* Make a series of refined extrapolations,
   * up to the specified maximum order, which
   * was calculated based on the Deuflhard
   * criterion upon state initialization.  */
  for (k = 0; k <= state->k_current; k++)
    {
      const unsigned int  N=  bd_sequence[k];
      const MpIeee r=  (h / N);
      const MpIeee x_k=  r * r;

      int  status=  bsimp_step_local (state,
                                     dim, t_local, h, N,
                                     y_extrap_save, yp,
                                     dfdt, dfdy,
                                     y_extrap_sequence, 
                                     sys);

      if (status == GSL_EBADFUNC) 
        {
          return GSL_EBADFUNC;
        }

      if (status == GSL_EFAILED)
        {
          /* If the local step fails, set the error to infinity in
             order to force a reduction in the step size */

          for (i = 0; i < dim; i++)
            {
              yerr[i] = GSL_POSINF;
            }

          break;
        }

      x[k] = x_k;

      poly_extrap (d, x, k, x_k, y_extrap_sequence, y, yerr, extrap_work, dim);
    }

  /* Evaluate dydt_out[]. */

  if (dydt_out != NULL)
    {
      int  s=  GSL_ODEIV_FN_EVAL (sys, t + h, y, dydt_out);

      if (s != GSL_SUCCESS)
        {
          DBL_MEMCPY (y, y_save, dim);
          DBL_MEMCPY (yerr, yerr_save, dim);
          return GSL_EBADFUNC;
        }
    }

  return GSL_SUCCESS;
}

static unsigned int
 bsimp_order(void *vstate)
{
  bsimp_state_t *state = (bsimp_state_t *) vstate;
  return state->order;
}

static int
 bsimp_reset(void *vstate, size_t dim)
{
  bsimp_state_t *state = (bsimp_state_t *) vstate;

  state->h_next = 0;

  DBL_ZERO_MEMSET (state->yp, dim);

  return GSL_SUCCESS;
}


static void
bsimp_free (void * vstate)
{
  bsimp_state_t *state = (bsimp_state_t *) vstate;

  free (state->delta);
  free (state->rhs_temp);

  gsl_matrix_free (state->dfdy);

  free (state->weight);
  free (state->delta_temp);
  free (state->y_temp);
  free (state->dfdt);
  free (state->extrap_work);
  free (state->y_extrap_sequence);
  free (state->y_extrap_save);
  free (state->y_save);
  free (state->yerr_save);
  free (state->yp);

  gsl_permutation_free (state->p_vec);
  gsl_matrix_free (state->a_mat);
  gsl_matrix_free (state->d);
  free (state);
}

static const gsl_odeiv_step_type bsimp_type = { 
  "bsimp",                      /* name */
  1,                            /* can use dydt_in */
  1,                            /* gives exact dydt_out */
  &bsimp_alloc,
  &bsimp_apply,
  &bsimp_reset,
  &bsimp_order,
  &bsimp_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_bsimp = &bsimp_type;
