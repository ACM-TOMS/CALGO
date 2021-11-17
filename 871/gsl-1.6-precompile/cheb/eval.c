#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cheb/eval.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
#include <gsl/gsl_chebyshev.h>

/* For efficiency there are separate implementations of each of these
   functions */

MpIeee gsl_cheb_eval(const gsl_cheb_series * cs, const MpIeee x)
{
  size_t i;
  MpIeee d1=  MpIeee( "0.0" );
  MpIeee d2=  MpIeee( "0.0" );

  MpIeee y=  (MpIeee( "2.0" ) * x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  for (i = cs->order; i >= 1; i--)
    {
      MpIeee temp=  d1;
      d1 = y2 * d1 - d2 + cs->c[i];
      d2 = temp;
    }

  return y * d1 - d2 + MpIeee( "0.5" ) * cs->c[0];
}

MpIeee gsl_cheb_eval_n(const gsl_cheb_series * cs, const size_t n, const MpIeee x)
{
  size_t i;
  MpIeee d1=  MpIeee( "0.0" );
  MpIeee d2=  MpIeee( "0.0" );

  size_t eval_order = GSL_MIN (n, cs->order);

  MpIeee y=  (MpIeee( "2.0" ) * x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  for (i = eval_order; i >= 1; i--)
    {
      MpIeee temp=  d1;
      d1 = y2 * d1 - d2 + cs->c[i];
      d2 = temp;
    }

  return y * d1 - d2 + MpIeee( "0.5" ) * cs->c[0];
}


int
 gsl_cheb_eval_err(const gsl_cheb_series * cs, const MpIeee x,
                   MpIeee *result, MpIeee *abserr)
{
  size_t i;
  MpIeee d1=  MpIeee( "0.0" );
  MpIeee d2=  MpIeee( "0.0" );

  MpIeee y=  (MpIeee( "2." ) * x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  MpIeee absc=  MpIeee( "0.0" );

  for (i = cs->order; i >= 1; i--)
    {
      MpIeee temp=  d1;
      d1 = y2 * d1 - d2 + cs->c[i];
      d2 = temp;
    }

  *result = y * d1 - d2 + MpIeee( "0.5" ) * cs->c[0];

  /* Estimate cumulative numerical error */

  for (i = 0; i <= cs->order; i++)
    {
      absc += fabs(cs->c[i]);
    }

  /* Combine truncation error and numerical error */

  *abserr = fabs (cs->c[cs->order]) + absc * GSL_DBL_EPSILON;

  return GSL_SUCCESS;
}

int
 gsl_cheb_eval_n_err(const gsl_cheb_series * cs,
                     const size_t n, const MpIeee x,
                     MpIeee *result, MpIeee *abserr)
{
  size_t i;
  MpIeee d1=  MpIeee( "0.0" );
  MpIeee d2=  MpIeee( "0.0" );

  MpIeee y=  (MpIeee( "2." ) * x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  MpIeee absc=  MpIeee( "0.0" );

  size_t eval_order = GSL_MIN (n, cs->order);

  for (i = eval_order; i >= 1; i--)
    {
      MpIeee temp=  d1;
      d1 = y2 * d1 - d2 + cs->c[i];
      d2 = temp;
    }

  *result = y * d1 - d2 + MpIeee( "0.5" ) * cs->c[0];

  /* Estimate cumulative numerical error */

  for (i = 0; i <= eval_order; i++)
    {
      absc += fabs(cs->c[i]);
    }

  /* Combine truncation error and numerical error */

  *abserr = fabs (cs->c[eval_order]) + absc * GSL_DBL_EPSILON;

  return GSL_SUCCESS;
}

int
 gsl_cheb_eval_mode_e(const gsl_cheb_series * cs,
                      const MpIeee x, gsl_mode_t mode,
                      MpIeee *result, MpIeee *abserr)
{
  size_t i;
  MpIeee d1=  MpIeee( "0.0" );
  MpIeee d2=  MpIeee( "0.0" );

  MpIeee y=  (MpIeee( "2." ) * x - cs->a - cs->b) / (cs->b - cs->a);
  MpIeee y2=  MpIeee( "2.0" ) * y;

  size_t eval_order;

  if (GSL_MODE_PREC (mode) == GSL_PREC_DOUBLE)
    eval_order = cs->order;
  else
    eval_order = cs->order_sp;

  for (i = eval_order; i >= 1; i--)
    {
      MpIeee temp=  d1;
      d1 = y2 * d1 - d2 + cs->c[i];
      d2 = temp;
    }

  *result = y * d1 - d2 + MpIeee( "0.5" ) * cs->c[0];
  *abserr = fabs (*result) * GSL_DBL_EPSILON + fabs (cs->c[eval_order]);

  return GSL_SUCCESS;
}

MpIeee gsl_cheb_eval_mode(const gsl_cheb_series * cs,
                    const MpIeee x, gsl_mode_t mode)
{
  MpIeee result;MpIeee  abserr;
  int  status=  gsl_cheb_eval_mode_e (cs, x, mode, &result, &abserr);

  if (status != GSL_SUCCESS) 
    {
      GSL_ERROR_VAL("gsl_cheb_eval_mode", status, result);
    };

  return result;
}


