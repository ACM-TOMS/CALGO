#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/interp_poly.c
 * 
 * Copyright (C) 2001 DAN, HO-JIN
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

/* Modified for standalone use in polynomial directory, B.Gough 2001 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

int
 gsl_poly_dd_init(MpIeee dd[], const MpIeee xa[], const MpIeee ya[],
                  size_t size)
{
  size_t i, j;

  /* Newton's divided differences */

  dd[0] = ya[0];

  for (j = size - 1; j >= 1; j--)
    {
      dd[j] = (ya[j] - ya[j - 1]) / (xa[j] - xa[j - 1]);
    }

  for (i = 2; i < size; i++)
    {
      for (j = size - 1; j >= i; j--)
        {
          dd[j] = (dd[j] - dd[j - 1]) / (xa[j] - xa[j - i]);
        }
    }

  return GSL_SUCCESS;
}

#ifndef HIDE_INLINE_STATIC
MpIeee gsl_poly_dd_eval(const MpIeee dd[], const MpIeee xa[], const size_t size, const MpIeee x)
{
  size_t i;
  MpIeee y=  dd[size - 1];

  for (i = size - 1; i--;)
    {
      y = dd[i] + (x - xa[i]) * y;
    }

  return y;
}
#endif

int
 gsl_poly_dd_taylor(MpIeee c[], MpIeee xp, 
                    const MpIeee dd[], const MpIeee xa[], size_t size,
                    MpIeee w[])
{
  size_t i, j;

  for (i = 0; i < size; i++)
    {
      c[i] = MpIeee( "0.0" );
      w[i] = MpIeee( "0.0" );
    }

  w[size - 1] = MpIeee( "1.0" );

  c[0] = dd[0];

  for (i = size - 1; i > 0 && i--;)
    {
      w[i] = -w[i + 1] * (xa[size - 2 - i] - xp);

      for (j = i + 1; j < size - 1; j++)
        {
          w[j] = w[j] - w[j + 1] * (xa[size - 2 - i] - xp);
        }

      for (j = i; j < size; j++)
        {
          c[j - i] += w[j] * dd[size - i - 1];
        }
    }

  return GSL_SUCCESS;
}  
