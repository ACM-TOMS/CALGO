#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* vector/reim_source.c
 * 
 * Copyright (C) 2001 Brian Gough
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

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
  REAL_TYPE(gsl_vector) s = NULL_VECTOR;

  s.data = v->data;
  s.size = v->size;
  s.stride = MULTIPLICITY * v->stride;
  s.block = 0;  /* FIXME: should be v->block, but cannot point to
                   block of different type */
  s.owner = 0;

  {
    QUALIFIED_REAL_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
    ((REAL_VIEW(gsl_vector,view) *)(&view))->vector = s;
    return view;
  }
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
  REAL_TYPE(gsl_vector) s = NULL_VECTOR;

  s.data = v->data + 1;
  s.size = v->size;
  s.stride = MULTIPLICITY * v->stride;
  s.block = 0;  /* FIXME: cannot point to block of different type */
  s.owner = 0;

  {
    QUALIFIED_REAL_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
    ((REAL_VIEW(gsl_vector,view) *)(&view))->vector = s;
    return view;
  }
}

