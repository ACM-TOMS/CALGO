#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/gammainv.c
 * 
 * Copyright (C) 2003 Brian Gough
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

MpIeee gsl_cdf_gamma_Pinv(const MpIeee P, const MpIeee a, const MpIeee b)
{
  MpIeee x;

  if (P == 1.0)
    {
      return GSL_POSINF;
    }
  else if (P == 0.0)
    {
      return MpIeee( "0.0" );
    }

  /* Consider, small, large and intermediate cases separately.  The
     boundaries at 0.05 and 0.95 have not been optimised, but seem ok
     for an initial approximation. */

  if (P < 0.05)
    {
      MpIeee x0=  exp ((gsl_sf_lngamma (a) + log (P)) / a);
      x = x0;
    }
  else if (P > 0.95)
    {
      MpIeee x0=  -log1p (-P) + gsl_sf_lngamma (a);
      x = x0;
    }
  else
    {
      MpIeee xg=  gsl_cdf_ugaussian_Pinv (P);
      MpIeee x0=  (xg < -sqrt (a)) ? a : sqrt (a) * xg + a;
      x = x0;
    }

  /* Use Lagrange's interpolation for E(x)/phi(x0) to work backwards
     to an improved value of x (Abramowitz & Stegun, 3.6.6) 

     where E(x)=P-integ(phi(u),u,x0,x) and phi(u) is the pdf.
   */

  {
    MpIeee lambda;MpIeee  dP;MpIeee  phi;

  start:
    dP = P - gsl_cdf_gamma_P (x, a, MpIeee( "1.0" ));
    phi = gsl_ran_gamma_pdf (x, a, MpIeee( "1.0" ));

    if (dP == MpIeee( "0.0" ))
      goto end;

    lambda = dP / GSL_MAX (MpIeee( "2" ) * fabs (dP / x), phi);

    {
      MpIeee step0=  lambda;
      MpIeee step1=  -((a - MpIeee( "1" )) / x - MpIeee( "1" )) * lambda * lambda / MpIeee( "4.0" );

      MpIeee step=  step0;
      if (fabs (step1) < fabs (step0))
        step += step1;

      if (x + step > MpIeee( "0" ))
        x += step;
      else
        {
          x /= MpIeee( "2.0" );
        }

      if (fabs (step0) > 1e-10 * x)
        goto start;
    }

  }

end:
  return b * x;
}

MpIeee gsl_cdf_gamma_Qinv(const MpIeee Q, const MpIeee a, const MpIeee b)
{
  MpIeee x;

  if (Q == 1.0)
    {
      return MpIeee( "0.0" );
    }
  else if (Q == 0.0)
    {
      return GSL_POSINF;
    }

  /* Consider, small, large and intermediate cases separately.  The
     boundaries at 0.05 and 0.95 have not been optimised, but seem ok
     for an initial approximation. */

  if (Q < 0.05)
    {
      MpIeee x0=  -log (Q) + gsl_sf_lngamma (a);
      x = x0;
    }
  else if (Q > 0.95)
    {
      MpIeee x0=  exp ((gsl_sf_lngamma (a) + log1p (-Q)) / a);
      x = x0;
    }
  else
    {
      MpIeee xg=  gsl_cdf_ugaussian_Qinv (Q);
      MpIeee x0=  (xg < -sqrt (a)) ? a : sqrt (a) * xg + a;
      x = x0;
    }

  /* Use Lagrange's interpolation for E(x)/phi(x0) to work backwards
     to an improved value of x (Abramowitz & Stegun, 3.6.6) 

     where E(x)=P-integ(phi(u),u,x0,x) and phi(u) is the pdf.
   */

  {
    MpIeee lambda;MpIeee  dQ;MpIeee  phi;

  start:
    dQ = Q - gsl_cdf_gamma_Q (x, a, MpIeee( "1.0" ));
    phi = gsl_ran_gamma_pdf (x, a, MpIeee( "1.0" ));

    if (dQ == MpIeee( "0.0" ))
      goto end;

    lambda = -dQ / GSL_MAX (MpIeee( "2" ) * fabs (dQ / x), phi);

    {
      MpIeee step0=  lambda;
      MpIeee step1=  -((a - MpIeee( "1" )) / x - MpIeee( "1" )) * lambda * lambda / MpIeee( "4.0" );

      MpIeee step=  step0;
      if (fabs (step1) < fabs (step0))
        step += step1;

      if (x + step > MpIeee( "0" ))
        x += step;
      else
        {
          x /= MpIeee( "2.0" );
        }

      if (fabs (step0) > 1e-10 * x)
        goto start;
    }

  }

end:
  return b * x;
}
