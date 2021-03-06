#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/add2d.c
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
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include "find2d.c"

int
 gsl_histogram2d_increment(gsl_histogram2d * h, MpIeee x, MpIeee y)
{
  int  status=  gsl_histogram2d_accumulate (h, x, y, 1.0);
  return status;
}

int
 gsl_histogram2d_accumulate(gsl_histogram2d * h,
                            MpIeee x, MpIeee y, MpIeee weight)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;

  size_t i = 0, j = 0;

  int  status=  find2d (h->nx, h->xrange, h->ny, h->yrange, x, y, &i, &j);

  if (status)
    {
      return GSL_EDOM;
    }

  if (i >= nx)
    {
      GSL_ERROR ("index lies outside valid range of 0 .. nx - 1",
                 GSL_ESANITY);
    }

  if (j >= ny)
    {
      GSL_ERROR ("index lies outside valid range of 0 .. ny - 1",
                 GSL_ESANITY);
    }

  h->bin[i * ny + j] += weight;

  return GSL_SUCCESS;
}
