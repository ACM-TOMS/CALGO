#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* linalg/tridiag.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2002, 2004 Gerard Jungman,
 * Brian Gough, David Necas
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

/* Author: G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "tridiag.h"
#include <gsl/gsl_linalg.h>

/* for description of method see [Engeln-Mullges + Uhlig, p. 92]
 *
 *     diag[0]  offdiag[0]             0   .....
 *  offdiag[0]     diag[1]    offdiag[1]   .....
 *           0  offdiag[1]       diag[2]
 *           0           0    offdiag[2]   .....
 */
static
int 
 solve_tridiag(
  const MpIeee diag[], size_t d_stride,
  const MpIeee offdiag[], size_t o_stride,
  const MpIeee b[], size_t b_stride,
  MpIeee x[], size_t x_stride,
  size_t N)
{
  int  status;
  MpIeee *gamma= new MpIeee[N];  // (MpIeee *) malloc (N * sizeof (MpIeee)); similar to ../block/init_source.c where we allocate matrices with new, remember to also replace free by delete!
  MpIeee *alpha= new MpIeee[N];  // (MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee *c=     new MpIeee[N];  // (MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee *z=     new MpIeee[N];  // (MpIeee *) malloc (N * sizeof (MpIeee));

  if (gamma == 0 || alpha == 0 || c == 0 || z == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      size_t i, j;

      /* Cholesky decomposition
         A = L.D.L^t
         lower_diag(L) = gamma
         diag(D) = alpha
       */
      alpha[0] = diag[0];
      gamma[0] = offdiag[0] / alpha[0];

      for (i = 1; i < N - 1; i++)
        {
          alpha[i] = diag[d_stride * i] - offdiag[o_stride*(i - 1)] * gamma[i - 1];
          gamma[i] = offdiag[o_stride * i] / alpha[i];
        }

      if (N > 1) 
        {
          alpha[N - 1] = diag[d_stride * (N - 1)] - offdiag[o_stride*(N - 2)] * gamma[N - 2];
        }

      /* update RHS */
      z[0] = b[0];
      for (i = 1; i < N; i++)
        {
          z[i] = b[b_stride * i] - gamma[i - 1] * z[i - 1];
        }
      for (i = 0; i < N; i++)
        {
          c[i] = z[i] / alpha[i];
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = c[N - 1];
      if (N >= 2)
        {
          for (i = N - 2, j = 0; j <= N - 2; j++, i--)
            {
              x[x_stride * i] = c[i] - gamma[i] * x[x_stride * (i + 1)];
            }
        }

      status = GSL_SUCCESS;
    }

  if (z != 0)
    delete [] (z);
  if (c != 0)
    delete [] (c);
  if (alpha != 0)
    delete [] (alpha);
  if (gamma != 0)
    delete [] (gamma);

  return status;
}

/* plain gauss elimination, only not bothering with the zeroes
 *
 *       diag[0]  abovediag[0]             0   .....
 *  belowdiag[0]       diag[1]  abovediag[1]   .....
 *             0  belowdiag[1]       diag[2]
 *             0             0  belowdiag[2]   .....
 */
static
int 
 solve_tridiag_nonsym(
  const MpIeee diag[], size_t d_stride,
  const MpIeee abovediag[], size_t a_stride,
  const MpIeee belowdiag[], size_t b_stride,
  const MpIeee rhs[], size_t r_stride,
  MpIeee x[], size_t x_stride,
  size_t N)
{
  cout << "ok parameter passing went fine, we are in solve_tridiag_nonsym" << endl;
  int  status;
  MpIeee *alpha= new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee *z=     new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  cout << "allocated alpha and z"<<endl;
  
  if (alpha == 0 || z == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      size_t i, j;

      /* Bidiagonalization (eliminating belowdiag)
         & rhs update
         diag' = alpha
         rhs' = z
       */
      alpha[0] = diag[0];
      z[0] = rhs[0];

      for (i = 1; i < N; i++)
        {
          const MpIeee t=  belowdiag[b_stride*(i - 1)]/alpha[i-1];
          alpha[i] = diag[d_stride*i] - t*abovediag[a_stride*(i - 1)];
          z[i] = rhs[r_stride*i] - t*z[i-1];
          /* FIXME!!! */
          if (alpha[i] == MpIeee( "0" )) {
            status = GSL_EZERODIV;
            goto solve_tridiag_nonsym_END;
          }
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = z[N - 1]/alpha[N - 1];
      if (N >= 2)
        {
          for (i = N - 2, j = 0; j <= N - 2; j++, i--)
            {
              x[x_stride * i] = (z[i] - abovediag[a_stride*i] * x[x_stride * (i + 1)])/alpha[i];
            }
        }

      status = GSL_SUCCESS;
    }

solve_tridiag_nonsym_END:
  cout << "cleaning up" << endl;
  if (z != 0)
    delete [] (z);
  if (alpha != 0)
    delete [] (alpha);

  return status;
}

/* for description of method see [Engeln-Mullges + Uhlig, p. 96]
 *
 *      diag[0]  offdiag[0]             0   .....  offdiag[N-1]
 *   offdiag[0]     diag[1]    offdiag[1]   .....
 *            0  offdiag[1]       diag[2]
 *            0           0    offdiag[2]   .....
 *          ...         ...
 * offdiag[N-1]         ...
 *
 */
static
int 
 solve_cyc_tridiag(
  const MpIeee diag[], size_t d_stride,
  const MpIeee offdiag[], size_t o_stride,
  const MpIeee b[], size_t b_stride,
  MpIeee x[], size_t x_stride,
  size_t N)
{
  int  status;
  MpIeee * delta=  new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee * gamma=  new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee * alpha=  new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee * c=  	   new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee * z=      new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));

  if (delta == 0 || gamma == 0 || alpha == 0 || c == 0 || z == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      size_t i, j;
      MpIeee sum=  MpIeee( "0.0" );

      /* factor */

      if (N == 1) 
        {
          x[0] = b[0] / diag[0];
          return GSL_SUCCESS;
        }

      alpha[0] = diag[0];
      gamma[0] = offdiag[0] / alpha[0];
      delta[0] = offdiag[o_stride * (N-1)] / alpha[0];

      for (i = 1; i < N - 2; i++)
        {
          alpha[i] = diag[d_stride * i] - offdiag[o_stride * (i-1)] * gamma[i - 1];
          gamma[i] = offdiag[o_stride * i] / alpha[i];
          delta[i] = -delta[i - 1] * offdiag[o_stride * (i-1)] / alpha[i];
        }

      for (i = 0; i < N - 2; i++)
        {
          sum += alpha[i] * delta[i] * delta[i];
        }

      alpha[N - 2] = diag[d_stride * (N - 2)] - offdiag[o_stride * (N - 3)] * gamma[N - 3];

      gamma[N - 2] = (offdiag[o_stride * (N - 2)] - offdiag[o_stride * (N - 3)] * delta[N - 3]) / alpha[N - 2];

      alpha[N - 1] = diag[d_stride * (N - 1)] - sum - alpha[(N - 2)] * gamma[N - 2] * gamma[N - 2];

      /* update */
      z[0] = b[0];
      for (i = 1; i < N - 1; i++)
        {
          z[i] = b[b_stride * i] - z[i - 1] * gamma[i - 1];
        }
      sum = MpIeee( "0.0" );
      for (i = 0; i < N - 2; i++)
        {
          sum += delta[i] * z[i];
        }
      z[N - 1] = b[b_stride * (N - 1)] - sum - gamma[N - 2] * z[N - 2];
      for (i = 0; i < N; i++)
        {
          c[i] = z[i] / alpha[i];
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = c[N - 1];
      x[x_stride * (N - 2)] = c[N - 2] - gamma[N - 2] * x[x_stride * (N - 1)];
      if (N >= 3)
        {
          for (i = N - 3, j = 0; j <= N - 3; j++, i--)
            {
              x[x_stride * i] = c[i] - gamma[i] * x[x_stride * (i + 1)] - delta[i] * x[x_stride * (N - 1)];
            }
        }

      status = GSL_SUCCESS;
    }

  if (z != 0)
    delete [] (z);
  if (c != 0)
    delete [] (c);
  if (alpha != 0)
    delete [] (alpha);
  if (gamma != 0)
    delete [] (gamma);
  if (delta != 0)
    delete [] (delta);

  return status;
}

/* solve following system w/o the corner elements and then use
 * Sherman-Morrison formula to compensate for them
 *
 *        diag[0]  abovediag[0]             0   .....  belowdiag[N-1]
 *   belowdiag[0]       diag[1]  abovediag[1]   .....
 *              0  belowdiag[1]       diag[2]
 *              0             0  belowdiag[2]   .....
 *            ...           ...
 * abovediag[N-1]           ...
 */
static
int  solve_cyc_tridiag_nonsym(
  const MpIeee diag[], size_t d_stride,
  const MpIeee abovediag[], size_t a_stride,
  const MpIeee belowdiag[], size_t b_stride,
  const MpIeee rhs[], size_t r_stride,
  MpIeee x[], size_t x_stride,
  size_t N)
{
  int  status;
  MpIeee *alpha=  new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee *zb=     new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee *zu=     new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee *w=      new MpIeee[N];  //(MpIeee *) malloc (N * sizeof (MpIeee));
  MpIeee beta;

  if (alpha == 0 || zb == 0 || zu == 0 || w == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      /* Bidiagonalization (eliminating belowdiag)
         & rhs update
         diag' = alpha
         rhs' = zb
         rhs' for Aq=u is zu
       */
      zb[0] = rhs[0];
      if (diag[0] != 0) beta = -diag[0]; else beta = MpIeee( "1" );
      {
        const MpIeee q=  1 - abovediag[0]*belowdiag[0]/(diag[0]*diag[d_stride]);
        if (fabs(q/beta) > 0.5 && fabs(q/beta) < 2) {
          beta *= (fabs(q/beta) < MpIeee( "1" )) ? MpIeee( "0.5" ) : MpIeee( "2" );
        }
      }
      zu[0] = beta;
      alpha[0] = diag[0] - beta;


      { 
        size_t i;
        for (i = 1; i+1 < N; i++)
        {
          const MpIeee t=  belowdiag[b_stride*(i - 1)]/alpha[i-1];
          alpha[i] = diag[d_stride*i] - t*abovediag[a_stride*(i - 1)];
          zb[i] = rhs[r_stride*i] - t*zb[i-1];
          zu[i] = -t*zu[i-1];
          /* FIXME!!! */
          if (alpha[i] == MpIeee( "0" )) {
            status = GSL_EZERODIV;
            goto solve_cyc_tridiag_nonsym_END;
          }
        }
      }

      {
        const size_t i = N-1;
        const MpIeee t=  belowdiag[b_stride*(i - 1)]/alpha[i-1];
        alpha[i] = diag[d_stride*i]
                   - abovediag[a_stride*i]*belowdiag[b_stride*i]/beta
                   - t*abovediag[a_stride*(i - 1)];
        zb[i] = rhs[r_stride*i] - t*zb[i-1];
        zu[i] = abovediag[a_stride*i] - t*zu[i-1];
        /* FIXME!!! */
        if (alpha[i] == MpIeee( "0" )) {
          status = GSL_EZERODIV;
          goto solve_cyc_tridiag_nonsym_END;
        }
      }


      /* backsubstitution */
      {
        size_t i, j;
        w[N-1] = zu[N-1]/alpha[N-1];
        x[N-1] = zb[N-1]/alpha[N-1];
        for (i = N - 2, j = 0; j <= N - 2; j++, i--)
          {
            w[i] = (zu[i] - abovediag[a_stride*i] * w[i+1])/alpha[i];
            x[i*x_stride] = (zb[i] - abovediag[a_stride*i] * x[x_stride*(i + 1)])/alpha[i];
          }
      }
      
      /* Sherman-Morrison */
      {
        const MpIeee vw=  w[0] + belowdiag[b_stride*(N - 1)]/beta * w[N-1];
        const MpIeee vx=  x[0] + belowdiag[b_stride*(N - 1)]/beta * x[x_stride*(N - 1)];
        /* FIXME!!! */
        if (vw + 1 == 0) {
          status = GSL_EZERODIV;
          goto solve_cyc_tridiag_nonsym_END;
        }

        {
          size_t i;
          for (i = 0; i < N; i++)
            x[i] -= vx/(MpIeee( "1" ) + vw)*w[i];
        }
      }

      status = GSL_SUCCESS;
    }

solve_cyc_tridiag_nonsym_END:
  if (zb != 0)
    delete [] (zb);
  if (zu != 0)
    delete [] (zu);
  if (w != 0)
    delete [] (w);
  if (alpha != 0)
    delete [] (alpha);

  return status;
}

int
 gsl_linalg_solve_symm_tridiag(
  const gsl_vector * diag,
  const gsl_vector * offdiag,
  const gsl_vector * rhs,
  gsl_vector * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (offdiag->size != rhs->size-1)
    {
      GSL_ERROR ("size of offdiag must match rhs-1", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else 
    {
      return solve_tridiag(diag->data, diag->stride,
                           offdiag->data, offdiag->stride,
                           rhs->data, rhs->stride,
                           solution->data, solution->stride,
                           diag->size);
    }
}

int
 gsl_linalg_solve_tridiag(
  const gsl_vector * diag,
  const gsl_vector * abovediag,
  const gsl_vector * belowdiag,
  const gsl_vector * rhs,
  gsl_vector * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (abovediag->size != rhs->size-1)
    {
      GSL_ERROR ("size of abovediag must match rhs-1", GSL_EBADLEN);
    }
  else if (belowdiag->size != rhs->size-1)
    {
      GSL_ERROR ("size of belowdiag must match rhs-1", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else 
    {
      return solve_tridiag_nonsym(diag->data, diag->stride,
                                  abovediag->data, abovediag->stride,
                                  belowdiag->data, belowdiag->stride,
                                  rhs->data, rhs->stride,
                                  solution->data, solution->stride,
                                  diag->size);
    }
}


int
 gsl_linalg_solve_symm_cyc_tridiag(
  const gsl_vector * diag,
  const gsl_vector * offdiag,
  const gsl_vector * rhs,
  gsl_vector * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (offdiag->size != rhs->size)
    {
      GSL_ERROR ("size of offdiag must match rhs", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else if (diag->size < 3)
    {
      GSL_ERROR ("size of cyclic system must be 3 or more", GSL_EBADLEN);
    }
  else 
    {
      return solve_cyc_tridiag(diag->data, diag->stride,
                               offdiag->data, offdiag->stride,
                               rhs->data, rhs->stride,
                               solution->data, solution->stride,
                               diag->size);
    }
}

int
 gsl_linalg_solve_cyc_tridiag(
  const gsl_vector * diag,
  const gsl_vector * abovediag,
  const gsl_vector * belowdiag,
  const gsl_vector * rhs,
  gsl_vector * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (abovediag->size != rhs->size)
    {
      GSL_ERROR ("size of abovediag must match rhs", GSL_EBADLEN);
    }
  else if (belowdiag->size != rhs->size)
    {
      GSL_ERROR ("size of belowdiag must match rhs", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else if (diag->size < 3)
    {
      GSL_ERROR ("size of cyclic system must be 3 or more", GSL_EBADLEN);
    }
  else 
    {
      return solve_cyc_tridiag_nonsym(diag->data, diag->stride,
                                      abovediag->data, abovediag->stride,
                                      belowdiag->data, belowdiag->stride,
                                      rhs->data, rhs->stride,
                                      solution->data, solution->stride,
                                      diag->size);
    }
}
