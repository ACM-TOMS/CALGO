#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* vector/minmax_source.c
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

MpIeee FUNCTION(gsl_vector,max) (const TYPE(gsl_vector) * v)
{
  /* finds the largest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  MpIeee max=  v->data[0 * stride];
  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee x=  v->data[i*stride];
      if (x > max)
        max = x;
    }

  return max;
}

MpIeee FUNCTION(gsl_vector,min) (const TYPE(gsl_vector) * v)
{
  /* finds the smallest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  MpIeee min=  v->data[0 * stride];
  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee x=  v->data[i*stride];
      if (x < min)
        min = x;
    }

  return min;
}

void
FUNCTION(gsl_vector,minmax) (const TYPE(gsl_vector) * v,
                             MpIeee * min_out, 
                             MpIeee * max_out)
{
  /* finds the smallest and largest elements of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  MpIeee max=  v->data[0 * stride];
  MpIeee min=  v->data[0 * stride];

  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee x=  v->data[i*stride];
      if (x < min)
        {
          min = x;
        }
      if (x > max)
        {
          max = x;
        }
    }

  *min_out = min;
  *max_out = max;
}


size_t 
FUNCTION(gsl_vector,max_index) (const TYPE(gsl_vector) * v)
{
  /* finds the largest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  MpIeee max=  v->data[0 * stride];
  size_t imax = 0;
  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee x=  v->data[i*stride];
      if (x > max)
        {
          max = x;
          imax = i;
        }
    }

  return imax;
}

size_t 
FUNCTION(gsl_vector,min_index) (const TYPE(gsl_vector) * v)
{
  /* finds the smallest element of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  MpIeee min=  v->data[0 * stride];
  size_t imin = 0;
  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee x=  v->data[i*stride];
      if (x < min)
        {
          min = x;
          imin = i;
        }
    }

  return imin;
}


void
FUNCTION(gsl_vector,minmax_index) (const TYPE(gsl_vector) * v,
                                   size_t * imin_out, 
                                   size_t * imax_out)
{
  /* finds the smallest and largest elements of a vector */

  const size_t N = v->size ;
  const size_t stride = v->stride ;

  size_t imin = 0, imax = 0;
  MpIeee max=  v->data[0 * stride];
  MpIeee min=  v->data[0 * stride];

  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee x=  v->data[i*stride];
      if (x < min)
        {
          min = x;
          imin = i;
        }
      if (x > max)
        {
          max = x;
          imax = i;
        }
    }

  *imin_out = imin;
  *imax_out = imax;
}


