#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sum/levin_u.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sum.h>

int
 gsl_sum_levin_u_accel(const MpIeee *array, const size_t array_size,
                       gsl_sum_levin_u_workspace * w, 
                       MpIeee *sum_accel, MpIeee *abserr)
{
  return gsl_sum_levin_u_minmax (array, array_size,
                                 0, array_size - 1, w, sum_accel, abserr);
}

int
 gsl_sum_levin_u_minmax(const MpIeee *array, const size_t array_size,
                        const size_t min_terms, const size_t max_terms,
                        gsl_sum_levin_u_workspace * w,
                        MpIeee *sum_accel, MpIeee *abserr)
{
  if (array_size == 0)
    {
      *sum_accel = MpIeee( "0.0" );
      *abserr = MpIeee( "0.0" );
      w->sum_plain = 0.0;
      w->terms_used = 0;
      return GSL_SUCCESS;
    }
  else if (array_size == 1)
    {
      *sum_accel = array[0];
      *abserr = MpIeee( "0.0" );
      w->sum_plain = array[0];
      w->terms_used = 1;
      return GSL_SUCCESS;
    }
  else
    {
      const MpIeee SMALL=  0.01;
      const size_t nmax = GSL_MAX (max_terms, array_size) - 1;
      MpIeee noise_n=  MpIeee( "0.0" );MpIeee  noise_nm1=  MpIeee( "0.0" );
      MpIeee trunc_n=  MpIeee( "0.0" );MpIeee  trunc_nm1=  MpIeee( "0.0" );
      MpIeee actual_trunc_n=  MpIeee( "0.0" );MpIeee  actual_trunc_nm1=  MpIeee( "0.0" );
      MpIeee result_n=  MpIeee( "0.0" );MpIeee  result_nm1=  MpIeee( "0.0" );
      MpIeee variance=  MpIeee( "0" );
      size_t n;
      unsigned int  i;
      int  better=  0;
      int  before=  0;
      int  converging=  0;
      MpIeee least_trunc=  GSL_DBL_MAX;
      MpIeee least_trunc_noise=  GSL_DBL_MAX;
      MpIeee least_trunc_result;

      /* Calculate specified minimum number of terms.  No convergence
         tests are made, and no truncation information is stored.  */

      for (n = 0; n < min_terms; n++)
        {
          const MpIeee t=  array[n];
          result_nm1 = result_n;
          gsl_sum_levin_u_step (t, n, nmax, w, &result_n);
        }

      least_trunc_result = result_n;

      variance = MpIeee( "0" );
      for (i = 0; i < n; i++)
        {
          MpIeee dn=  w->dsum[i] * GSL_MACH_EPS * array[i];
          variance += dn * dn;
        }
      noise_n = sqrt (variance);

      /* Calculate up to maximum number of terms.  Check truncation
         condition.  */

      for (; n <= nmax; n++)
        {
          const MpIeee t=  array[n];

          result_nm1 = result_n;
          gsl_sum_levin_u_step (t, n, nmax, w, &result_n);

          /* Compute the truncation error directly */

          actual_trunc_nm1 = actual_trunc_n;
          actual_trunc_n = fabs (result_n - result_nm1);

          /* Average results to make a more reliable estimate of the
             real truncation error */

          trunc_nm1 = trunc_n;
          trunc_n = MpIeee( "0.5" ) * (actual_trunc_n + actual_trunc_nm1);

          noise_nm1 = noise_n;
          variance = MpIeee( "0" );

          for (i = 0; i <= n; i++)
            {
              MpIeee dn=  w->dsum[i] * GSL_MACH_EPS * array[i];
              variance += dn * dn;
            }

          noise_n = sqrt (variance);

          /* Determine if we are in the convergence region.  */

          better = (trunc_n < trunc_nm1 || trunc_n < SMALL * fabs (result_n));
          converging = converging || (better && before);
          before = better;

          if (converging)
            {
              if (trunc_n < least_trunc)
                {
                  /* Found a low truncation point in the convergence
                     region. Save it. */

                  least_trunc_result = result_n;
                  least_trunc = trunc_n;
                  least_trunc_noise = noise_n;
                }

              if (noise_n > trunc_n / MpIeee( "3.0" ))
                break;

              if (trunc_n < MpIeee( "10.0" ) * GSL_MACH_EPS * fabs (result_n))
                break;
            }

        }

      if (converging)
        {
          /* Stopped in the convergence region.  Return result and
             error estimate.  */

          *sum_accel = least_trunc_result;
          *abserr = GSL_MAX_DBL (least_trunc, least_trunc_noise);
          w->terms_used = n;
          return GSL_SUCCESS;
        }
      else
        {
          /* Never reached the convergence region.  Use the last
             calculated values.  */

          *sum_accel = result_n;
          *abserr = GSL_MAX_DBL (trunc_n, noise_n);
          w->terms_used = n;
          return GSL_SUCCESS;
        }
    }
}


int
 gsl_sum_levin_u_step(const MpIeee term, const size_t n, const size_t nmax,
                      gsl_sum_levin_u_workspace * w, MpIeee *sum_accel)
{

#define I(i,j) ((i)*(nmax+1) + (j))

  if (n == 0)
    {
      *sum_accel = term;
      w->sum_plain = term;

      w->q_den[0] = 1.0 / term;
      w->q_num[0] = 1.0;

      w->dq_den[I (0, 0)] = -1.0 / (term * term);
      w->dq_num[I (0, 0)] = 0.0;

      w->dsum[0] = 1.0;

      return GSL_SUCCESS;
    }
  else
    {
      MpIeee result;
      MpIeee factor=  MpIeee( "1.0" );
      MpIeee ratio=  (MpIeee) n / (n + MpIeee( "1.0" ));
      unsigned int  i;
      int  j;

      w->sum_plain += term;

      w->q_den[n] = 1.0 / (term * (n + 1.0) * (n + 1.0));
      w->q_num[n] = w->sum_plain * w->q_den[n];

      for (i = 0; i < n; i++)
        {
          w->dq_den[I (i, n)] = 0;
          w->dq_num[I (i, n)] = w->q_den[n];
        }

      w->dq_den[I (n, n)] = -w->q_den[n] / term;
      w->dq_num[I (n, n)] =
        w->q_den[n] + w->sum_plain * (w->dq_den[I (n, n)]);

      for (j = n - 1; j >= 0; j--)
        {
          MpIeee c=  factor * (j + MpIeee( "1" )) / (n + MpIeee( "1" ));
          factor *= ratio;
          w->q_den[j] = w->q_den[j + 1] - c * w->q_den[j];
          w->q_num[j] = w->q_num[j + 1] - c * w->q_num[j];

          for (i = 0; i < n; i++)
            {
              w->dq_den[I (i, j)] =
                w->dq_den[I (i, j + 1)] - c * w->dq_den[I (i, j)];
              w->dq_num[I (i, j)] =
                w->dq_num[I (i, j + 1)] - c * w->dq_num[I (i, j)];
            }

          w->dq_den[I (n, j)] = w->dq_den[I (n, j + 1)];
          w->dq_num[I (n, j)] = w->dq_num[I (n, j + 1)];
        }

      result = w->q_num[0] / w->q_den[0];

      *sum_accel = result;

      for (i = 0; i <= n; i++)
        {
          w->dsum[i] =
            (w->dq_num[I (i, 0)] -
             result * w->dq_den[I (i, 0)]) / w->q_den[0];
        }

      return GSL_SUCCESS;
    }
}
