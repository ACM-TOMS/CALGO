#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/params2d.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram2d.h>

MpIeee gsl_histogram2d_xmax(const gsl_histogram2d * h)
{
  const int nx = h->nx;
  return h->xrange[nx];
}

MpIeee gsl_histogram2d_xmin(const gsl_histogram2d * h)
{
  return h->xrange[0];
}

MpIeee gsl_histogram2d_ymax(const gsl_histogram2d * h)
{
  const int ny = h->ny;
  return h->yrange[ny];
}

MpIeee gsl_histogram2d_ymin(const gsl_histogram2d * h)
{
  return h->yrange[0];
}

size_t
gsl_histogram2d_nx (const gsl_histogram2d * h)
{
  return h->nx;
}

size_t
gsl_histogram2d_ny (const gsl_histogram2d * h)
{
  return h->ny;
}
