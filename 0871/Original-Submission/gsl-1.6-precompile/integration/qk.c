#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/qk.c
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

#include <config.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "err.c"

void
gsl_integration_qk (const int n, 
                    const MpIeee xgk[], const MpIeee wg[], const MpIeee wgk[],
                    MpIeee fv1[], MpIeee fv2[],
                    const gsl_function * f, MpIeee a, MpIeee b,
                    MpIeee *result, MpIeee *abserr,
                    MpIeee *resabs, MpIeee *resasc)
{

  const MpIeee center=  0.5 * (a + b);
  const MpIeee half_length=  0.5 * (b - a);
  const MpIeee abs_half_length=  fabs (half_length);
  const MpIeee f_center=  GSL_FN_EVAL (f, center);

  MpIeee result_gauss=  MpIeee( "0" );
  MpIeee result_kronrod=  f_center * wgk[n - 1];

  MpIeee result_abs=  fabs (result_kronrod);
  MpIeee result_asc=  MpIeee( "0" );
  MpIeee mean=  MpIeee( "0" );MpIeee  err=  MpIeee( "0" );

  int  j;

  if (n % 2 == 0)
    {
      result_gauss = f_center * wg[n / 2 - 1];
    }

  for (j = 0; j < (n - 1) / 2; j++)
    {
      const int jtw = j * 2 + 1;        /* j=1,2,3 jtw=2,4,6 */
      const MpIeee abscissa=  half_length * xgk[jtw];
      const MpIeee fval1=  GSL_FN_EVAL (f, center - abscissa);
      const MpIeee fval2=  GSL_FN_EVAL (f, center + abscissa);
      const MpIeee fsum=  fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg[j] * fsum;
      result_kronrod += wgk[jtw] * fsum;
      result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
    }

  for (j = 0; j < n / 2; j++)
    {
      int  jtwm1=  j * 2;
      const MpIeee abscissa=  half_length * xgk[jtwm1];
      const MpIeee fval1=  GSL_FN_EVAL (f, center - abscissa);
      const MpIeee fval2=  GSL_FN_EVAL (f, center + abscissa);
      fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk[jtwm1] * (fval1 + fval2);
      result_abs += wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    };

  mean = result_kronrod * MpIeee( "0.5" );

  result_asc = wgk[n - 1] * fabs (f_center - mean);

  for (j = 0; j < n - 1; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

  /* scale by the width of the integration region */

  err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  *result = result_kronrod;
  *resabs = result_abs;
  *resasc = result_asc;
  *abserr = rescale_error (err, result_abs, result_asc);

}
