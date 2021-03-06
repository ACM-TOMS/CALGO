#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* vector/vector_source.c
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

#ifndef HIDE_INLINE_STATIC
MpIeee FUNCTION(gsl_vector, get) (const TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          MpIeee zero=  ZERO;
          GSL_ERROR_VAL ("index out of range", GSL_EINVAL, zero);
        }
    }

  /* The following line is a generalization of return v->data[i] */

  return *(MpIeee *) (v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION (gsl_vector, set) (TYPE (gsl_vector) * v, const size_t i, MpIeee x)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          GSL_ERROR_VOID ("index out of range", GSL_EINVAL);
        }
    }

  /* The following line is a generalization of v->data[i] = x */

  *(MpIeee*) (v->data + MULTIPLICITY * i * v->stride) = x;
}

MpIeee *
FUNCTION(gsl_vector, ptr) (TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          GSL_ERROR_NULL ("index out of range", GSL_EINVAL);
        }
    }

  return (MpIeee *) (v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION (gsl_vector, const_ptr) (const TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          GSL_ERROR_NULL ("index out of range", GSL_EINVAL);
        }
    }

  /* The following line is a generalization of return v->data[i] */

  return (const BASE *) (v->data + MULTIPLICITY * i * v->stride);
}
#endif
