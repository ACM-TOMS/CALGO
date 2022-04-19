#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* eigen/jacobi.c
 * 
 * Copyright (C) 2004 Brian Gough, Gerard Jungman
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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

/* Algorithm 8.4.3 - Cyclic Jacobi.  Golub & Van Loan, Matrix Computations */

static inline MpIeee symschur2(gsl_matrix * A, size_t p, size_t q, MpIeee *c, MpIeee *s)
{
  MpIeee Apq=  gsl_matrix_get (A, p, q);

  if (Apq != MpIeee( "0.0" ))
    {
      MpIeee App=  gsl_matrix_get (A, p, p);
      MpIeee Aqq=  gsl_matrix_get (A, q, q);
      MpIeee tau=  (Aqq - App) / (MpIeee( "2.0" ) * Apq);
      MpIeee t;MpIeee  c1;

      if (tau >= MpIeee( "0.0" ))
        {
          t = MpIeee( "1.0" ) / (tau + hypot (MpIeee( "1.0" ), tau));
        }
      else
        {
          t = -MpIeee( "1.0" ) / (-tau + hypot (MpIeee( "1.0" ), tau));
        }

      c1 = MpIeee( "1.0" ) / hypot (MpIeee( "1.0" ), t);

      *c = c1;
      *s = t * c1;
    }
  else
    {
      *c = MpIeee( "1.0" );
      *s = MpIeee( "0.0" );
    }

  /* reduction in off(A) is 2*(A_pq)^2 */

  return fabs (Apq);
}

inline static void
apply_jacobi_L (gsl_matrix * A, size_t p, size_t q, MpIeee c, MpIeee s)
{
  size_t j;
  const size_t N = A->size2;

  /* Apply rotation to matrix A,  A' = J^T A */

  for (j = 0; j < N; j++)
    {
      MpIeee Apj=  gsl_matrix_get (A, p, j);
      MpIeee Aqj=  gsl_matrix_get (A, q, j);
      gsl_matrix_set (A, p, j, Apj * c - Aqj * s);
      gsl_matrix_set (A, q, j, Apj * s + Aqj * c);
    }
}

inline static void
apply_jacobi_R (gsl_matrix * A, size_t p, size_t q, MpIeee c, MpIeee s)
{
  size_t i;
  const size_t M = A->size1;

  /* Apply rotation to matrix A,  A' = A J */

  for (i = 0; i < M; i++)
    {
      MpIeee Aip=  gsl_matrix_get (A, i, p);
      MpIeee Aiq=  gsl_matrix_get (A, i, q);
      gsl_matrix_set (A, i, p, Aip * c - Aiq * s);
      gsl_matrix_set (A, i, q, Aip * s + Aiq * c);
    }
}

inline static MpIeee norm(gsl_matrix * A)
{
  size_t i, j, M = A->size1, N = A->size2;
  MpIeee sum=  MpIeee( "0.0" );MpIeee  scale=  MpIeee( "0.0" );MpIeee  ssq=  MpIeee( "1.0" );

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          MpIeee Aij=  gsl_matrix_get (A, i, j);

          if (Aij != MpIeee( "0.0" ))
            {
              MpIeee ax=  fabs (Aij);

              if (scale < ax)
                {
                  ssq = MpIeee( "1.0" ) + ssq * (scale / ax) * (scale / ax);
                  scale = ax;
                }
              else
                {
                  ssq += (ax / scale) * (ax / scale);
                }
            }

        }
    }

  sum = scale * sqrt (ssq);

  return sum;
}

int
 gsl_eigen_jacobi(gsl_matrix * a,
                  gsl_vector * eval,
                  gsl_matrix * evec, unsigned int  max_rot, unsigned int  *nrot)
{
  size_t i, p, q;
  const size_t M = a->size1, N = a->size2;
  MpIeee red;MpIeee  redsum=  MpIeee( "0.0" );

  if (M != N)
    {
      GSL_ERROR ("eigenproblem requires square matrix", GSL_ENOTSQR);
    }
  else if (M != evec->size1 || M != evec->size2)
    {
      GSL_ERROR ("eigenvector matrix must match input matrix", GSL_EBADLEN);
    }
  else if (M != eval->size)
    {
      GSL_ERROR ("eigenvalue vector must match input matrix", GSL_EBADLEN);
    }

  gsl_vector_set_zero (eval);
  gsl_matrix_set_identity (evec);

  for (i = 0; i < max_rot; i++)
    {
      MpIeee nrm=  norm (a);

      if (nrm == MpIeee( "0.0" ))
        break;

      for (p = 0; p < N; p++)
        {
          for (q = p + 1; q < N; q++)
            {
              MpIeee c;MpIeee  s;

              red = symschur2 (a, p, q, &c, &s);
              redsum += red;

              /* Compute A <- J^T A J */
              apply_jacobi_L (a, p, q, c, s);
              apply_jacobi_R (a, p, q, c, s);

              /* Compute V <- V J */
              apply_jacobi_R (evec, p, q, c, s);
            }
        }
    }

  *nrot = i;

  for (p = 0; p < N; p++)
    {
      MpIeee ep=  gsl_matrix_get (a, p, p);
      gsl_vector_set (eval, p, ep);
    }

  if (i == max_rot)
    {
      return GSL_EMAXITER;
    }

  return GSL_SUCCESS;
}

int
 gsl_eigen_invert_jacobi(const gsl_matrix * a,
                         gsl_matrix * ainv, unsigned int  max_rot)
{
  if (a->size1 != a->size2 || ainv->size1 != ainv->size2)
    {
      GSL_ERROR("jacobi method requires square matrix", GSL_ENOTSQR);
    }
  else if (a->size1 != ainv->size2)
    {
     GSL_ERROR ("inverse matrix must match input matrix", GSL_EBADLEN);
    }
  
  {
    const size_t n = a->size2;
    size_t i,j,k;
    size_t nrot = 0;
    int  status;

    gsl_vector * eval = gsl_vector_alloc(n);
    gsl_matrix * evec = gsl_matrix_alloc(n, n);
    gsl_matrix * tmp = gsl_matrix_alloc(n, n);

    gsl_matrix_memcpy (tmp, a);

    status = gsl_eigen_jacobi(tmp, eval, evec, max_rot, &nrot);
      
    for(i=0; i<n; i++) 
      {
        for(j=0; j<n; j++) 
          {
            MpIeee ainv_ij=  MpIeee( "0.0" );
            
            for(k = 0; k<n; k++)
              {
                MpIeee f=  MpIeee( "1.0" ) / gsl_vector_get(eval, k);
                MpIeee vik=  gsl_matrix_get (evec, i, k);
                MpIeee vjk=  gsl_matrix_get (evec, j, k);
                ainv_ij += vik * vjk * f;
              }
            gsl_matrix_set (ainv, i, j, ainv_ij);
          }
      }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(tmp);

    if (status)
      {
        return status;
      }
    else
      {
        return GSL_SUCCESS;
      }
  }
}
