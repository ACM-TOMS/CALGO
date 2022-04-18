#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/tdistinv.c
 *
 * Copyright (C) 2002 Jason H. Stover.
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <stdio.h>

static MpIeee inv_cornish_fisher(MpIeee z, MpIeee nu)
{
  MpIeee a=  MpIeee( "1" ) / (nu - MpIeee( "0.5" ));
  MpIeee b=  MpIeee( "48.0" ) / (a * a);

  MpIeee cf1=  z * (MpIeee( "3" ) + z * z);
  MpIeee cf2=  z * (MpIeee( "945" ) + z * z * (MpIeee( "360" ) + z * z * (MpIeee( "63" ) + z * z * MpIeee( "4" ))));

  MpIeee y=  z - cf1 / b + cf2 / (MpIeee( "10" ) * b * b);

  MpIeee t=  GSL_SIGN (z) * sqrt (nu * expm1 (a * y * y));

  return t;
}


MpIeee gsl_cdf_tdist_Pinv(const MpIeee P, const MpIeee nu)
{
  MpIeee x;MpIeee  ptail;

  if (P == 1.0)
    {
      return GSL_POSINF;
    }
  else if (P == 0.0)
    {
      return GSL_NEGINF;
    }

  if (nu == 1.0)
    {
      x = tan (M_PI * (P - MpIeee( "0.5" )));
    }
  else if (nu == 2.0)
    {
      MpIeee a=  MpIeee( "2" ) * P - MpIeee( "1" );
      x = a / sqrt (MpIeee( "2" ) * (MpIeee( "1" ) - a * a));
    }

  ptail = (P < MpIeee( "0.5" )) ? P : MpIeee( "1" ) - P;

  if (sqrt (M_PI * nu / 2) * ptail > pow (MpIeee( "0.05" ), nu / MpIeee( "2" )))
    {
      MpIeee xg=  gsl_cdf_ugaussian_Pinv (P);
      x = inv_cornish_fisher (xg, nu);
    }
  else
    {
      /* Use an asymptotic expansion of the tail of integral */

      MpIeee beta=  gsl_sf_beta (MpIeee( "0.5" ), nu / MpIeee( "2" ));

      if (P < 0.5)
        {
          x = -sqrt (nu) * pow (beta * nu * P, -MpIeee( "1.0" ) / nu);
        }
      else
        {
          x = sqrt (nu) * pow (beta * nu * (MpIeee( "1" ) - P), -MpIeee( "1.0" ) / nu);
        }

      /* Correct nu -> nu/(1+nu/x^2) in the leading term to account
         for higher order terms. This avoids overestimating x, which
         makes the iteration unstable due to the rapidly decreasing
         tails of the distribution. */

      x /= sqrt (MpIeee( "1" ) + nu / (x * x));
    }

  {
    MpIeee dP;MpIeee  phi;

  start:
    dP = P - gsl_cdf_tdist_P (x, nu);
    phi = gsl_ran_tdist_pdf (x, nu);

    if (dP == MpIeee( "0.0" ))
      goto end;

    {
      MpIeee lambda=  dP / phi;
      MpIeee step0=  lambda;
      MpIeee step1=  ((nu + MpIeee( "1" )) * x / (x * x + nu)) * (lambda * lambda / MpIeee( "4.0" ));

      MpIeee step=  step0;

      if (fabs (step1) < fabs (step0))
        {
          step += step1;
        }

      if (P > 0.5 && x + step < 0)
        x /= MpIeee( "2" );
      else if (P < 0.5 && x + step > 0)
        x /= MpIeee( "2" );
      else
        x += step;

      if (fabs (step) > 1e-10 * fabs (x))
        goto start;
    }
  }

end:

  return x;
}

MpIeee gsl_cdf_tdist_Qinv(const MpIeee Q, const MpIeee nu)
{
  MpIeee x;MpIeee  qtail;

  if (Q == 0.0)
    {
      return GSL_POSINF;
    }
  else if (Q == 1.0)
    {
      return GSL_NEGINF;
    }

  if (nu == 1.0)
    {
      x = tan (M_PI * (MpIeee( "0.5" ) - Q));
    }
  else if (nu == 2.0)
    {
      MpIeee a=  MpIeee( "2" ) * (MpIeee( "1" ) - Q) - MpIeee( "1" );
      x = a / sqrt (MpIeee( "2" ) * (MpIeee( "1" ) - a * a));
    }

  qtail = (Q < MpIeee( "0.5" )) ? Q : MpIeee( "1" ) - Q;

  if (sqrt (M_PI * nu / 2) * qtail > pow (MpIeee( "0.05" ), nu / MpIeee( "2" )))
    {
      MpIeee xg=  gsl_cdf_ugaussian_Qinv (Q);
      x = inv_cornish_fisher (xg, nu);
    }
  else
    {
      /* Use an asymptotic expansion of the tail of integral */

      MpIeee beta=  gsl_sf_beta (MpIeee( "0.5" ), nu / MpIeee( "2" ));

      if (Q < 0.5)
        {
          x = sqrt (nu) * pow (beta * nu * Q, -MpIeee( "1.0" ) / nu);
        }
      else
        {
          x = -sqrt (nu) * pow (beta * nu * (MpIeee( "1" ) - Q), -MpIeee( "1.0" ) / nu);
        }

      /* Correct nu -> nu/(1+nu/x^2) in the leading term to account
         for higher order terms. This avoids overestimating x, which
         makes the iteration unstable due to the rapidly decreasing
         tails of the distribution. */

      x /= sqrt (MpIeee( "1" ) + nu / (x * x));
    }

  {
    MpIeee dQ;MpIeee  phi;

  start:
    dQ = Q - gsl_cdf_tdist_Q (x, nu);
    phi = gsl_ran_tdist_pdf (x, nu);

    if (dQ == MpIeee( "0.0" ))
      goto end;

    {
      MpIeee lambda=  - dQ / phi;
      MpIeee step0=  lambda;
      MpIeee step1=  ((nu + MpIeee( "1" )) * x / (x * x + nu)) * (lambda * lambda / MpIeee( "4.0" ));

      MpIeee step=  step0;

      if (fabs (step1) < fabs (step0))
        {
          step += step1;
        }

      if (Q < 0.5 && x + step < 0)
        x /= MpIeee( "2" );
      else if (Q > 0.5 && x + step > 0)
        x /= MpIeee( "2" );
      else
        x += step;

      if (fabs (step) > 1e-10 * fabs (x))
        goto start;
    }
  }

end:

  return x;
}
