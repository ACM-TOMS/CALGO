#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* eigen/symm.c
 * 
 * Copyright (C) 2001 Brian Gough
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

/* Compute eigenvalues/eigenvectors of real symmetric matrix using
   reduction to tridiagonal form, followed by QR iteration with
   implicit shifts.

   See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 8.3
   */

#include "qrstep.c"

gsl_eigen_symm_workspace *
gsl_eigen_symm_alloc (const size_t n)
{
  gsl_eigen_symm_workspace *w;

  if (n == 0)
    {
      GSL_ERROR_NULL ("matrix dimension must be positive integer",
                      GSL_EINVAL);
    }

  w = ((gsl_eigen_symm_workspace *)
       malloc (sizeof (gsl_eigen_symm_workspace)));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->d = (MpIeee*) malloc (n * sizeof (MpIeee));

  if (w->d == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for diagonal", GSL_ENOMEM);
    }

  w->sd = (MpIeee*) malloc (n * sizeof (MpIeee));

  if (w->sd == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for subdiagonal", GSL_ENOMEM);
    }

  w->size = n;

  return w;
}

void
gsl_eigen_symm_free (gsl_eigen_symm_workspace * w)
{
  free (w->sd);
  free (w->d);
  free (w);
}


int
 gsl_eigen_symm(gsl_matrix * A, gsl_vector * eval,
                     gsl_eigen_symm_workspace * w)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("matrix must be square to compute eigenvalues", GSL_ENOTSQR);
    }
  else if (eval->size != A->size1)
    {
      GSL_ERROR ("eigenvalue vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      MpIeee *const d=  w->d;
      MpIeee *const sd=  w->sd;

      size_t a, b;

      /* handle special case */

      if (N == 1)
        {
          MpIeee A00=  gsl_matrix_get (A, 0,0);
          gsl_vector_set (eval, 0, A00);
          return GSL_SUCCESS;
        }

      /* use sd as the temporary workspace for the decomposition,
         since we can discard the tau result immediately if we are not
         computing eigenvectors */

      {
        gsl_vector_view d_vec = gsl_vector_view_array (d, N);
        gsl_vector_view sd_vec = gsl_vector_view_array (sd, N - 1);
        gsl_vector_view tau = gsl_vector_view_array (sd, N - 1);
        gsl_linalg_symmtd_decomp (A, &tau.vector);
        gsl_linalg_symmtd_unpack_T (A, &d_vec.vector, &sd_vec.vector);
      }
      
      /* Make an initial pass through the tridiagonal decomposition
         to remove off-diagonal elements which are effectively zero */
      
      chop_small_elements (N, d, sd);
      
      /* Progressively reduce the matrix until it is diagonal */
      
      b = N - 1;
      
      while (b > 0)
        {
          if (sd[b - 1] == MpIeee( "0.0" ) || (sd[b - 1]).isNan())
            {
              b--;
              continue;
            }
          
          /* Find the largest unreduced block (a,b) starting from b
             and working backwards */
          
          a = b - 1;
          
          while (a > 0)
            {
              if (sd[a - 1] == MpIeee( "0.0" ))
                {
                  break;
                }
              a--;
            }
          
          {
            const size_t n_block = b - a + 1;
            MpIeee *d_block=  d + a;
            MpIeee *sd_block=  sd + a;
            
            /* apply QR reduction with implicit deflation to the
               unreduced block */
            
            qrstep (n_block, d_block, sd_block, NULL, NULL);
            
            /* remove any small off-diagonal elements */
            
            chop_small_elements (n_block, d_block, sd_block);
          }
        }
      
      {
        gsl_vector_view d_vec = gsl_vector_view_array (d, N);
        gsl_vector_memcpy (eval, &d_vec.vector);
      }

      return GSL_SUCCESS;
    }
}
