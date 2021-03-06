#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* interpolation/bsearch.c
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

/* Author:  G. Jungman
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_interp.h>

#ifndef HIDE_INLINE_STATIC
size_t
gsl_interp_bsearch (
  const MpIeee x_array[], MpIeee x,
  size_t index_lo,
  size_t index_hi
  )
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while (ihi > ilo + 1)
    {
      size_t i = (ihi + ilo) / 2;
      if (x_array[i] > x)
        ihi = i;
      else
        ilo = i;
    }

  return ilo;
}
#endif
