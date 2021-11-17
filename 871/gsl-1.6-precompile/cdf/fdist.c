#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/f.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "beta_inc.c"

/*
 * Lower tail.
 */
MpIeee gsl_cdf_fdist_P(const MpIeee x, const MpIeee nu1, const MpIeee nu2)
{
  MpIeee P;
  MpIeee r=  nu2 / nu1;

  if (x < r)
    {
      MpIeee u=  x / (r + x);

      P = beta_inc_AXPY (MpIeee( "1.0" ), MpIeee( "0.0" ), nu1 / MpIeee( "2.0" ), nu2 / MpIeee( "2.0" ), u);
    }
  else
    {
      MpIeee u=  r / (r + x);

      P = beta_inc_AXPY (-MpIeee( "1.0" ), MpIeee( "1.0" ), nu2 / MpIeee( "2.0" ), nu1 / MpIeee( "2.0" ), u);
    }

  return P;
}

/*
 * Upper tail.
 */
MpIeee gsl_cdf_fdist_Q(const MpIeee x, const MpIeee nu1, const MpIeee nu2)
{
  MpIeee P;
  MpIeee r=  nu2 / nu1;

  if (x < r)
    {
      MpIeee u=  x / (r + x);

      P = beta_inc_AXPY (-MpIeee( "1.0" ), MpIeee( "1.0" ), nu1 / MpIeee( "2.0" ), nu2 / MpIeee( "2.0" ), u);
    }
  else
    {
      MpIeee u=  r / (r + x);

      P = beta_inc_AXPY (MpIeee( "1.0" ), MpIeee( "0.0" ), nu2 / MpIeee( "2.0" ), nu1 / MpIeee( "2.0" ), u);
    }

  return P;
}
