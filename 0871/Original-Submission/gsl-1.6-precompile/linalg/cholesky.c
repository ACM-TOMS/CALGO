#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* Cholesky Decomposition
 *
 * Copyright (C) 2000  Thomas Walter
 *
 * 3 May 2000: Modified for GSL by Brian Gough
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

/*
 * Cholesky decomposition of a symmetrix positive definite matrix.
 * This is useful to solve the matrix arising in
 *    periodic cubic splines
 *    approximating splines
 *
 * This algorthm does:
 *   A = L * L'
 * with
 *   L  := lower left triangle matrix
 *   L' := the transposed form of L.
 *
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int
 gsl_linalg_cholesky_decomp(gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t i,j,k;
      int  status=  0;

      /* Do the first 2 rows explicitly.  It is simple, and faster.  And
       * one can return if the matrix has only 1 or 2 rows.  
       */

      MpIeee A_00=  gsl_matrix_get (A, 0,0);
      
      MpIeee L_00=  sqrt(A_00);
      
      if (A_00 <= MpIeee( "0" ))
        {
          status = GSL_EDOM ;
        }

      gsl_matrix_set (A, 0, 0, L_00);
  
      if (M > 1)
        {
          MpIeee A_10=  gsl_matrix_get (A, 1,0);
          MpIeee A_11=  gsl_matrix_get (A, 1,1);
          
          MpIeee L_10=  A_10 / L_00;
          MpIeee diag=  A_11 - L_10 * L_10;
          MpIeee L_11=  sqrt(diag);
          
          if (diag <= MpIeee( "0" ))
            {
              status = GSL_EDOM;
            }

          gsl_matrix_set (A, 1, 0, L_10);        
          gsl_matrix_set (A, 1, 1, L_11);
        }
      
      for (k = 2; k < M; k++)
        {
          MpIeee A_kk=  gsl_matrix_get (A, k, k);
          
          for (i = 0; i < k; i++)
            {
              MpIeee sum=  MpIeee( "0" );

              MpIeee A_ki=  gsl_matrix_get (A, k, i);
              MpIeee A_ii=  gsl_matrix_get (A, i, i);

              gsl_vector_view ci = gsl_matrix_row (A, i);
              gsl_vector_view ck = gsl_matrix_row (A, k);

              if (i > 0) {
                gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
                gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);
                
                gsl_blas_ddot (&di.vector, &dk.vector, &sum);
              }

              A_ki = (A_ki - sum) / A_ii;
              gsl_matrix_set (A, k, i, A_ki);
            } 

          {
            gsl_vector_view ck = gsl_matrix_row (A, k);
            gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);
            
            MpIeee sum=  gsl_blas_dnrm2 (&dk.vector);
            MpIeee diag=  A_kk - sum * sum;

            MpIeee L_kk=  sqrt(diag);
            
            if (diag <= MpIeee( "0" ))
              {
                status = GSL_EDOM;
              }
            
            gsl_matrix_set (A, k, k, L_kk);
          }
        }

      /* Now copy the transposed lower triangle to the upper triangle,
       * the diagonal is common.  
       */
      
      for (i = 1; i < M; i++)
        {
          for (j = 0; j < i; j++)
            {
              MpIeee A_ij=  gsl_matrix_get (A, i, j);
              gsl_matrix_set (A, j, i, A_ij);
            }
        } 
      
      if (status == GSL_EDOM)
        {
          GSL_ERROR ("matrix must be positive definite", GSL_EDOM);
        }
      
      return GSL_SUCCESS;
    }
}


int
 gsl_linalg_cholesky_solve(const gsl_matrix * LLT,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);


      return GSL_SUCCESS;
    }
}

int
 gsl_linalg_cholesky_svx(const gsl_matrix * LLT,
                         gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);

      return GSL_SUCCESS;
    }
}



