#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fft/c_init.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

TYPE(gsl_fft_complex_wavetable) * 
FUNCTION(gsl_fft_complex_wavetable,alloc) (size_t n)
{
  int  status;
  size_t i;
  size_t n_factors;
  size_t t, product, product_1, q;
  MpIeee d_theta;

  TYPE(gsl_fft_complex_wavetable) * wavetable ;

  if (n == 0)
    {
      GSL_ERROR_VAL ("length n must be positive integer", GSL_EDOM, 0);
    }

  wavetable = (TYPE(gsl_fft_complex_wavetable) *) 
    malloc(sizeof(TYPE(gsl_fft_complex_wavetable)));

  if (wavetable == NULL)
    {
      GSL_ERROR_VAL ("failed to allocate struct", GSL_ENOMEM, 0);
    }

  wavetable->trig = (TYPE(gsl_complex) *) malloc (n * sizeof (TYPE(gsl_complex)));

  if (wavetable->trig == NULL)
    {
      free(wavetable) ; /* error in constructor, prevent memory leak */

      GSL_ERROR_VAL ("failed to allocate trigonometric lookup table", 
                        GSL_ENOMEM, 0);
    }

  wavetable->n = n ;

  status = fft_complex_factorize (n, &n_factors, wavetable->factor);

  if (status)
    {
      /* exception in constructor, avoid memory leak */

      free (wavetable->trig);
      free (wavetable);         

      GSL_ERROR_VAL ("factorization failed", GSL_EFACTOR, 0);
    };

  wavetable->nf = n_factors;

  d_theta = -MpIeee( "2.0" ) * M_PI / ((MpIeee) n);

  t = 0;
  product = 1;
  for (i = 0; i < n_factors; i++)
    {
      size_t j;
      const size_t factor = wavetable->factor[i];
      wavetable->twiddle[i] = wavetable->trig + t;
      product_1 = product;      /* product_1 = p_(i-1) */
      product *= factor;
      q = n / product;

      for (j = 1; j < factor; j++)
        {
          size_t k;
          size_t m = 0;
          for (k = 1; k <= q; k++)
            {
              MpIeee theta;
              m = m + j * product_1;
              m = m % n;
              theta = d_theta * m;      /*  d_theta*j*k*p_(i-1) */
              GSL_REAL(wavetable->trig[t]) = cos (theta);
              GSL_IMAG(wavetable->trig[t]) = sin (theta);

              t++;
            }
        }
    }

  if (t > n)
    {
      /* exception in constructor, avoid memory leak */

      free (wavetable->trig);
      free (wavetable);

      GSL_ERROR_VAL ("overflowed trigonometric lookup table", 
                        GSL_ESANITY, 0);
    }

  return wavetable;
}


TYPE(gsl_fft_complex_workspace) * 
FUNCTION(gsl_fft_complex_workspace,alloc) (size_t n)
{
  TYPE(gsl_fft_complex_workspace) * workspace ;

  if (n == 0)
    {
      GSL_ERROR_VAL ("length n must be positive integer", GSL_EDOM, 0);
    }

  workspace = (TYPE(gsl_fft_complex_workspace) *) 
    malloc(sizeof(TYPE(gsl_fft_complex_workspace)));

  if (workspace == NULL)
    {
      GSL_ERROR_VAL ("failed to allocate struct", GSL_ENOMEM, 0);
    }

  workspace->n = n ;

  workspace->scratch = (MpIeee*) malloc (2 * n * sizeof (MpIeee));

  if (workspace->scratch == NULL)
    {
      free(workspace) ; /* error in constructor, prevent memory leak */

      GSL_ERROR_VAL ("failed to allocate scratch space", GSL_ENOMEM, 0);
    }
  
  return workspace;
}


void
FUNCTION(gsl_fft_complex_wavetable,free) (TYPE(gsl_fft_complex_wavetable) * wavetable)
{

  /* release trigonometric lookup tables */

  free (wavetable->trig);
  wavetable->trig = NULL;

  free (wavetable) ;
}

void
FUNCTION(gsl_fft_complex_workspace,free) (TYPE(gsl_fft_complex_workspace) * workspace)
{
  /* release scratch space */

  free (workspace->scratch);
  workspace->scratch = NULL;
  free (workspace) ;
}


int
 FUNCTION(gsl_fft_complex,memcpy) (TYPE(gsl_fft_complex_wavetable) * dest,
                                  TYPE(gsl_fft_complex_wavetable) * src)
{
  int  i;int   n;int   nf;

  if (dest->n != src->n) 
    {
      GSL_ERROR ("length of src and dest do not match", GSL_EINVAL);
    } 
  
  n = dest->n ;
  nf = dest->nf ;

  memcpy(dest->trig, src->trig, n * sizeof (MpIeee)) ;
  
  for (i = 0 ; i < nf ; i++)
    {
      dest->twiddle[i] = dest->trig + (src->twiddle[i] - src->trig) ;
    }

  return 0 ;
}
