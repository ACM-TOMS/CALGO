#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fft/test_real.c
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

#include "bitreverse.h"
#include "signals.h"
#include "compare.h"

void FUNCTION(test_real,func) (size_t stride, size_t n);
void FUNCTION(test_real,bitreverse_order) (size_t stride, size_t n);
void FUNCTION(test_real,radix2) (size_t stride, size_t n);

void FUNCTION(test_real,func) (size_t stride, size_t n) 
{
  size_t i ;
  int  status;

  TYPE(gsl_fft_real_wavetable) * rw ;
  TYPE(gsl_fft_halfcomplex_wavetable) * hcw ;
  TYPE(gsl_fft_real_workspace) * rwork ;

  MpIeee * real_data=  (MpIeee *) malloc (n * stride * sizeof (MpIeee));
  MpIeee * complex_data=  (MpIeee *) malloc (MpIeee( "2" ) * n * stride * sizeof (MpIeee));
  MpIeee * complex_tmp=  (MpIeee *) malloc (MpIeee( "2" ) * n * stride * sizeof (MpIeee));
  MpIeee * fft_complex_data=  (MpIeee *) malloc (MpIeee( "2" ) * n * stride * sizeof (MpIeee));

  for (i = 0 ; i < n * stride ; i++)
    {
      real_data[i] = (MpIeee)i ;
    }

  for (i = 0 ; i < 2 * n * stride ; i++)
    {
      complex_data[i] = (MpIeee)(i + MpIeee( "1000.0" )) ;
      complex_tmp[i] = (MpIeee)(i + MpIeee( "2000.0" )) ;
      fft_complex_data[i] = (MpIeee)(i + MpIeee( "3000.0" )) ;
    }

  gsl_set_error_handler (NULL); /* abort on any errors */
  
  /* mixed radix real fft */
  
  rw = FUNCTION(gsl_fft_real_wavetable,alloc) (n);
  gsl_test (rw == 0, NAME(gsl_fft_real_wavetable) 
            "_alloc, n = %d, stride = %d", n, stride);

  rwork = FUNCTION(gsl_fft_real_workspace,alloc) (n);
  gsl_test (rwork == 0, NAME(gsl_fft_real_workspace) 
            "_alloc, n = %d", n);
    
  FUNCTION(fft_signal,real_noise) (n, stride, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, 2 * n * stride * sizeof (MpIeee));

  for (i = 0; i < n; i++)
    {
      real_data[i*stride] = REAL(complex_data,stride,i);
    }
  
  FUNCTION(gsl_fft_real,transform) (real_data, stride, n, rw, rwork);
  FUNCTION(gsl_fft_halfcomplex,unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("dft", fft_complex_data,
                                              "fft of noise", complex_data,
                                              stride, n, 1e6);
  gsl_test (status, NAME(gsl_fft_real) 
            " with signal_real_noise, n = %d, stride = %d", n, stride);
  
  /* compute the inverse fft */

  hcw = FUNCTION(gsl_fft_halfcomplex_wavetable,alloc) (n);
  gsl_test (hcw == 0, NAME(gsl_fft_halfcomplex_wavetable) 
            "_alloc, n = %d, stride = %d", n, stride);

  status = FUNCTION(gsl_fft_halfcomplex,transform) (real_data, stride, n, hcw, rwork);
  
  for (i = 0; i < n; i++)
    {
      real_data[i*stride] /= n;
    }
  
  FUNCTION(gsl_fft_real,unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("orig", complex_tmp,
                                              "fft inverse", complex_data,
                                              stride, n, 1e6);

  gsl_test (status, NAME(gsl_fft_halfcomplex) 
            " with data from signal_noise, n = %d, stride = %d", n, stride);

  FUNCTION(gsl_fft_real_workspace,free) (rwork);
  FUNCTION(gsl_fft_real_wavetable,free) (rw);
  FUNCTION(gsl_fft_halfcomplex_wavetable,free) (hcw);

  free(real_data) ;
  free(complex_data) ;
  free(complex_tmp) ;
  free(fft_complex_data) ;
}


void 
FUNCTION(test_real,bitreverse_order) (size_t stride, size_t n) 
{
  int  status;
  size_t logn, i ;

  MpIeee * tmp=  (MpIeee *) malloc (n * stride * sizeof (MpIeee));
  MpIeee * data=  (MpIeee *) malloc (n * stride * sizeof (MpIeee));
  MpIeee * reversed_data=  (MpIeee *) malloc (n * stride * sizeof (MpIeee));
  
  for (i = 0; i <  stride * n; i++) 
    {
      data[i] = (MpIeee)i ;
    }

  memcpy (tmp, data, n * stride * sizeof(MpIeee)) ;

  logn = 0 ;  while (n > (1U <<logn)) {logn++;} ;

  /* do a naive bit reversal as a baseline for testing the other routines */

  for (i = 0; i < n; i++) 
    {
      size_t i_tmp = i ;
      size_t j = 0 ;
      size_t bit ;

      for (bit = 0; bit < logn; bit++)
        {
          j <<= 1;              /* reverse shift i into j */
          j |= i_tmp & 1;
          i_tmp >>= 1;
        }

      reversed_data[j*stride] = data[i*stride] ;
    }

  FUNCTION(fft_real,bitreverse_order) (data, stride, n, logn);

  status = FUNCTION(compare_real,results) ("naive bit reverse", 
                                           reversed_data,
                                           "fft_complex_bitreverse_order", 
                                           data,
                                           stride, n, 1e6);

  gsl_test (status, NAME(fft_real) "_bitreverse_order, n = %d", n);

  free (reversed_data) ;
  free (data) ;
  free (tmp) ;
}


void FUNCTION(test_real,radix2) (size_t stride, size_t n) 
{
  size_t i ;
  int  status;

  MpIeee * real_data=  (MpIeee *) malloc (n * stride * sizeof (MpIeee));
  MpIeee * complex_data=  (MpIeee *) malloc (MpIeee( "2" ) * n * stride * sizeof (MpIeee));
  MpIeee * complex_tmp=  (MpIeee *) malloc (MpIeee( "2" ) * n * stride * sizeof (MpIeee));
  MpIeee * fft_complex_data=  (MpIeee *) malloc (MpIeee( "2" ) * n * stride * sizeof (MpIeee));

  for (i = 0 ; i < n * stride ; i++)
    {
      real_data[i] = (MpIeee)i ;
    }

  for (i = 0 ; i < 2 * n * stride ; i++)
    {
      complex_data[i] = (MpIeee)(i + MpIeee( "1000.0" )) ;
      complex_tmp[i] = (MpIeee)(i + MpIeee( "2000.0" )) ;
      fft_complex_data[i] = (MpIeee)(i + MpIeee( "3000.0" )) ;
    }

  gsl_set_error_handler (NULL); /* abort on any errors */
  
  /* radix-2 real fft */
  
  FUNCTION(fft_signal,real_noise) (n, stride, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, 2 * n * stride * sizeof (MpIeee));

  for (i = 0; i < n; i++)
    {
      real_data[i*stride] = REAL(complex_data,stride,i);
    }
  
  FUNCTION(gsl_fft_real,radix2_transform) (real_data, stride, n);
  FUNCTION(gsl_fft_halfcomplex,radix2_unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("dft", fft_complex_data,
                                              "fft of noise", complex_data,
                                              stride, n, 1e6);
  gsl_test (status, NAME(gsl_fft_real) 
            "_radix2 with signal_real_noise, n = %d, stride = %d", n, stride);
  
  /* compute the inverse fft */
  
  status = FUNCTION(gsl_fft_halfcomplex,radix2_transform) (real_data, stride, n);
  
  for (i = 0; i < n; i++)
    {
      real_data[i*stride] /= n;
    }
  
  FUNCTION(gsl_fft_real,unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("orig", complex_tmp,
                                              "fft inverse", complex_data,
                                              stride, n, 1e6);

  gsl_test (status, NAME(gsl_fft_halfcomplex) 
            "_radix2 with data from signal_noise, n = %d, stride = %d", n, stride);


  free(real_data) ;
  free(complex_data) ;
  free(complex_tmp) ;
  free(fft_complex_data) ;
}
