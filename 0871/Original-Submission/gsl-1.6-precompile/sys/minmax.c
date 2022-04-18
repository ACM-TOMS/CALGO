#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sys/minmax.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

#ifndef HIDE_INLINE_STATIC
int  GSL_MAX_INT(int  a, int  b);
int  GSL_MIN_INT(int  a, int  b);
MpIeee GSL_MAX_DBL(MpIeee a, MpIeee b);
MpIeee GSL_MIN_DBL(MpIeee a, MpIeee b);
MpIeee GSL_MAX_LDBL(MpIeee a, MpIeee b);
MpIeee GSL_MIN_LDBL(MpIeee a, MpIeee b);

int
 GSL_MAX_INT(int  a, int  b)
{
  return GSL_MAX (a, b);
}

int
 GSL_MIN_INT(int  a, int  b)
{
  return GSL_MIN (a, b);
}

MpIeee GSL_MAX_DBL(MpIeee a, MpIeee b)
{
  return GSL_MAX (a, b);
}

MpIeee GSL_MIN_DBL(MpIeee a, MpIeee b)
{
  return GSL_MIN (a, b);
}

MpIeee GSL_MAX_LDBL(MpIeee a, MpIeee b)
{
  return GSL_MAX (a, b);
}

MpIeee GSL_MIN_LDBL(MpIeee a, MpIeee b)
{
  return GSL_MIN (a, b);
}
#endif

/* Define some static functions which are always available */

MpIeee gsl_max(MpIeee a, MpIeee b);
MpIeee gsl_min(MpIeee a, MpIeee b);

MpIeee gsl_max(MpIeee a, MpIeee b)
{
  return GSL_MAX (a, b);
}

MpIeee gsl_min(MpIeee a, MpIeee b)
{
  return GSL_MIN (a, b);
}

