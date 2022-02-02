/* multifit/multilinear.c
 * 
 * Copyright (C) 2000 Brian Gough
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

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* Fit
 *
 *  y = X c
 *
 *  where X is an M x N matrix of M observations for N variables.
 *
 */

int
gsl_multifit_linear (const gsl_matrix * X,
                     const gsl_vector * y,
                     gsl_vector * c,
                     gsl_matrix * cov,
                     double *chisq, gsl_multifit_linear_workspace * work)
{
  if (X->size1 != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match rows of matrix X",
         GSL_EBADLEN);
    }
  else if (X->size2 != c->size)
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1)
    {
      GSL_ERROR
        ("number of parameters does not match size of covariance matrix",
         GSL_EBADLEN);
    }
  else if (X->size1 != work->n || X->size2 != work->p)
    {
      GSL_ERROR
        ("size of workspace does not match size of observation matrix",
         GSL_EBADLEN);
    }
  else
    {
      const size_t n = X->size1;
      const size_t p = X->size2;

      size_t i, j;

      gsl_matrix *A = work->A;
      gsl_matrix *Q = work->Q;
      gsl_matrix *QSI = work->QSI;
      gsl_vector *S = work->S;
      gsl_vector *xt = work->xt;
      gsl_vector *D = work->D;

      /* Copy X to workspace,  A <= X */

      gsl_matrix_memcpy (A, X);

      /* Balance the columns of the matrix A */

      gsl_linalg_balance_columns (A, D);

      /* Decompose A into U S Q^T */

      gsl_linalg_SV_decomp_mod (A, QSI, Q, S, xt);

      /* Solve y = A c for c */

      gsl_blas_dgemv (CblasTrans, 1.0, A, y, 0.0, xt);

      /* Scale the matrix Q,  Q' = Q S^-1 */

      gsl_matrix_memcpy (QSI, Q);

      for (j = 0; j < p; j++)
        {
          gsl_vector_view column = gsl_matrix_column (QSI, j);
          double alpha = gsl_vector_get (S, j);
          if (alpha != 0)
            alpha = 1.0 / alpha;
          gsl_vector_scale (&column.vector, alpha);
        }

      gsl_vector_set_zero (c);

      gsl_blas_dgemv (CblasNoTrans, 1.0, QSI, xt, 0.0, c);

      /* Unscale the balancing factors */

      gsl_vector_div (c, D);

      /* Compute chisq, from residual r = y - X c */

      {
        double s2 = 0, r2 = 0;

        for (i = 0; i < n; i++)
          {
            double yi = gsl_vector_get (y, i);
            gsl_vector_const_view row = gsl_matrix_const_row (X, i);
            double y_est, ri;
            gsl_blas_ddot (&row.vector, c, &y_est);
            ri = yi - y_est;
            r2 += ri * ri;
          }

        s2 = r2 / (n - p);

        *chisq = r2;

        /* Form variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */

        for (i = 0; i < p; i++)
          {
            gsl_vector_view row_i = gsl_matrix_row (QSI, i);
            double d_i = gsl_vector_get (D, i);

            for (j = i; j < p; j++)
              {
                gsl_vector_view row_j = gsl_matrix_row (QSI, j);
                double d_j = gsl_vector_get (D, j);
                double s;

                gsl_blas_ddot (&row_i.vector, &row_j.vector, &s);

                gsl_matrix_set (cov, i, j, s * s2 / (d_i * d_j));
                gsl_matrix_set (cov, j, i, s * s2 / (d_i * d_j));
              }
          }
      }

      return GSL_SUCCESS;
    }
}

int
gsl_multifit_wlinear (const gsl_matrix * X,
                      const gsl_vector * w,
                      const gsl_vector * y,
                      gsl_vector * c,
                      gsl_matrix * cov,
                      double *chisq, gsl_multifit_linear_workspace * work)
{
  if (X->size1 != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match rows of matrix X",
         GSL_EBADLEN);
    }
  else if (X->size2 != c->size)
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (w->size != y->size)
    {
      GSL_ERROR ("number of weights does not match number of observations",
                 GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1)
    {
      GSL_ERROR
        ("number of parameters does not match size of covariance matrix",
         GSL_EBADLEN);
    }
  else if (X->size1 != work->n || X->size2 != work->p)
    {
      GSL_ERROR
        ("size of workspace does not match size of observation matrix",
         GSL_EBADLEN);
    }
  else
    {
      const size_t n = X->size1;
      const size_t p = X->size2;

      size_t i, j;

      gsl_matrix *A = work->A;
      gsl_matrix *Q = work->Q;
      gsl_matrix *QSI = work->QSI;
      gsl_vector *S = work->S;
      gsl_vector *t = work->t;
      gsl_vector *xt = work->xt;
      gsl_vector *D = work->D;

      /* Scale X,  A = sqrt(w) X */

      gsl_matrix_memcpy (A, X);

      for (i = 0; i < n; i++)
        {
          double wi = gsl_vector_get (w, i);

          if (wi < 0)
            wi = 0;

          {
            gsl_vector_view row = gsl_matrix_row (A, i);
            gsl_vector_scale (&row.vector, sqrt (wi));
          }
        }

      /* Balance the columns of the matrix A */

      gsl_linalg_balance_columns (A, D);

      /* Decompose A into U S Q^T */

      gsl_linalg_SV_decomp_mod (A, QSI, Q, S, xt);

      /* Solve sqrt(w) y = A c for c, by first computing t = sqrt(w) y */

      for (i = 0; i < n; i++)
        {
          double wi = gsl_vector_get (w, i);
          double yi = gsl_vector_get (y, i);
          if (wi < 0)
            wi = 0;
          gsl_vector_set (t, i, sqrt (wi) * yi);
        }

      gsl_blas_dgemv (CblasTrans, 1.0, A, t, 0.0, xt);

      /* Scale the matrix Q,  Q' = Q S^-1 */

      gsl_matrix_memcpy (QSI, Q);

      for (j = 0; j < p; j++)
        {
          gsl_vector_view column = gsl_matrix_column (QSI, j);
          double alpha = gsl_vector_get (S, j);
          if (alpha != 0)
            alpha = 1.0 / alpha;
          gsl_vector_scale (&column.vector, alpha);
        }

      gsl_vector_set_zero (c);

      /* Solution */

      gsl_blas_dgemv (CblasNoTrans, 1.0, QSI, xt, 0.0, c);

      /* Unscale the balancing factors */

      gsl_vector_div (c, D);

      /* Form covariance matrix cov = (Q S^-1) (Q S^-1)^T */

      for (i = 0; i < p; i++)
        {
          gsl_vector_view row_i = gsl_matrix_row (QSI, i);
          double d_i = gsl_vector_get (D, i);

          for (j = i; j < p; j++)
            {
              gsl_vector_view row_j = gsl_matrix_row (QSI, j);
              double d_j = gsl_vector_get (D, j);
              double s;

              gsl_blas_ddot (&row_i.vector, &row_j.vector, &s);

              gsl_matrix_set (cov, i, j, s / (d_i * d_j));
              gsl_matrix_set (cov, j, i, s / (d_i * d_j));
            }
        }

      /* Compute chisq, from residual r = y - X c */

      {
        double r2 = 0;

        for (i = 0; i < n; i++)
          {
            double yi = gsl_vector_get (y, i);
            double wi = gsl_vector_get (w, i);
            gsl_vector_const_view row = gsl_matrix_const_row (X, i);
            double y_est, ri;
            gsl_blas_ddot (&row.vector, c, &y_est);
            ri = yi - y_est;
            r2 += wi * ri * ri;
          }

        *chisq = r2;
      }

      return GSL_SUCCESS;
    }
}
