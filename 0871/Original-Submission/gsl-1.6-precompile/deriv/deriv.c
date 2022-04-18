#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* deriv/deriv.c
 * 
 * Copyright (C) 2004 Brian Gough
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>

static void
central_deriv (const gsl_function * f, MpIeee x, MpIeee h,
               MpIeee *result, MpIeee *abserr_round, MpIeee *abserr_trunc)
{
  /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
     x+h/2, x+h). Note that the central point is not used.  

     Compute the error using the difference between the 5-point and
     the 3-point rule (x-h,x,x+h). Again the central point is not
     used. */

  MpIeee fm1=  GSL_FN_EVAL (f, x - h);
  MpIeee fp1=  GSL_FN_EVAL (f, x + h);

  MpIeee fmh=  GSL_FN_EVAL (f, x - h / MpIeee( "2" ));
  MpIeee fph=  GSL_FN_EVAL (f, x + h / MpIeee( "2" ));

  MpIeee r3=  MpIeee( "0.5" ) * (fp1 - fm1);
  MpIeee r5=  (MpIeee( "4.0" ) / MpIeee( "3.0" )) * (fph - fmh) - (MpIeee( "1.0" ) / MpIeee( "3.0" )) * r3;

  MpIeee e3=  (fabs (fp1) + fabs (fm1)) * GSL_DBL_EPSILON;
  MpIeee e5=  MpIeee( "2.0" ) * (fabs (fph) + fabs (fmh)) * GSL_DBL_EPSILON + e3;

  MpIeee dy=  GSL_MAX (fabs (r3), fabs (r5)) * fabs (x) * GSL_DBL_EPSILON;

  /* The truncation error in the r5 approximation itself is O(h^4).
     However, for safety, we estimate the error from r5-r3, which is
     O(h^2).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r5. */

  *result = r5 / h;
  *abserr_trunc = fabs ((r5 - r3) / h); /* Estimated truncation error O(h^2) */
  *abserr_round = fabs (e5 / h) + dy;   /* Rounding error (cancellations) */
}

int
 gsl_deriv_central(const gsl_function * f, MpIeee x, MpIeee h,
                   MpIeee *result, MpIeee *abserr)
{
  MpIeee r_0;MpIeee  round;MpIeee  trunc;MpIeee  error;
  central_deriv (f, x, h, &r_0, &round, &trunc);
  error = round + trunc;

  if (round < trunc && (round > MpIeee( "0" ) && trunc > MpIeee( "0" )))
    {
      MpIeee r_opt;MpIeee  round_opt;MpIeee  trunc_opt;MpIeee  error_opt;

      /* Compute an optimised stepsize to minimize the total error,
         using the scaling of the truncation error (O(h^2)) and
         rounding error (O(1/h)). */

      MpIeee h_opt=  h * pow (round / (MpIeee( "2.0" ) * trunc), MpIeee( "1.0" ) / MpIeee( "3.0" ));
      central_deriv (f, x, h_opt, &r_opt, &round_opt, &trunc_opt);
      error_opt = round_opt + trunc_opt;

      /* Check that the new error is smaller, and that the new derivative 
         is consistent with the error bounds of the original estimate. */

      if (error_opt < error && fabs (r_opt - r_0) < MpIeee( "4.0" ) * error)
        {
          r_0 = r_opt;
          error = error_opt;
        }
    }

  *result = r_0;
  *abserr = error;

  return GSL_SUCCESS;
}


static void
forward_deriv (const gsl_function * f, MpIeee x, MpIeee h,
               MpIeee *result, MpIeee *abserr_round, MpIeee *abserr_trunc)
{
  /* Compute the derivative using the 4-point rule (x+h/4, x+h/2,
     x+3h/4, x+h).

     Compute the error using the difference between the 4-point and
     the 2-point rule (x+h/2,x+h).  */

  MpIeee f1=  GSL_FN_EVAL (f, x + h / MpIeee( "4.0" ));
  MpIeee f2=  GSL_FN_EVAL (f, x + h / MpIeee( "2.0" ));
  MpIeee f3=  GSL_FN_EVAL (f, x + (MpIeee( "3.0" ) / MpIeee( "4.0" )) * h);
  MpIeee f4=  GSL_FN_EVAL (f, x + h);

  MpIeee r2=  MpIeee( "2.0" )*(f4 - f2);
  MpIeee r4=  (MpIeee( "22.0" ) / MpIeee( "3.0" )) * (f4 - f3) - (MpIeee( "62.0" ) / MpIeee( "3.0" )) * (f3 - f2) +
    (MpIeee( "52.0" ) / MpIeee( "3.0" )) * (f2 - f1);

  /* Estimate the rounding error for r4 */

  MpIeee e4=  MpIeee( "2" ) * MpIeee( "20.67" ) * (fabs (f4) + fabs (f3) + fabs (f2) + fabs (f1)) * GSL_DBL_EPSILON;

  MpIeee dy=  GSL_MAX (fabs (r2), fabs (r4)) * fabs (x) * GSL_DBL_EPSILON;

  /* The truncation error in the r4 approximation itself is O(h^3).
     However, for safety, we estimate the error from r4-r2, which is
     O(h).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r4. */

  *result = r4 / h;
  *abserr_trunc = fabs ((r4 - r2) / h); /* Estimated truncation error O(h) */
  *abserr_round = fabs (e4 / h) + dy;
}

int
 gsl_deriv_forward(const gsl_function * f, MpIeee x, MpIeee h,
                   MpIeee *result, MpIeee *abserr)
{
  MpIeee r_0;MpIeee  round;MpIeee  trunc;MpIeee  error;
  forward_deriv (f, x, h, &r_0, &round, &trunc);
  error = round + trunc;

  if (round < trunc && (round > MpIeee( "0" ) && trunc > MpIeee( "0" )))
    {
      MpIeee r_opt;MpIeee  round_opt;MpIeee  trunc_opt;MpIeee  error_opt;

      /* Compute an optimised stepsize to minimize the total error,
         using the scaling of the estimated truncation error (O(h)) and
         rounding error (O(1/h)). */

      MpIeee h_opt=  h * pow (round / (trunc), MpIeee( "1.0" ) / MpIeee( "2.0" ));
      forward_deriv (f, x, h_opt, &r_opt, &round_opt, &trunc_opt);
      error_opt = round_opt + trunc_opt;

      /* Check that the new error is smaller, and that the new derivative 
         is consistent with the error bounds of the original estimate. */

      if (error_opt < error && fabs (r_opt - r_0) < MpIeee( "4.0" ) * error)
        {
          r_0 = r_opt;
          error = error_opt;
        }
    }

  *result = r_0;
  *abserr = error;

  return GSL_SUCCESS;
}

int
 gsl_deriv_backward(const gsl_function * f, MpIeee x, MpIeee h,
                    MpIeee *result, MpIeee *abserr)
{
  return gsl_deriv_forward (f, x, -h, result, abserr);
}
