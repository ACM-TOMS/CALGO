#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* matrix/rowcol_source.c
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

QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
  
  if (i >= m->size1)
    {
      GSL_ERROR_VAL ("row index is out of range", GSL_EINVAL, view);
    }
  
  {
    TYPE(gsl_vector) v = NULL_VECTOR;
    
    v.data = m->data + i * MULTIPLICITY * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;
    
    ((VIEW(_gsl_vector, view) *)&view)->vector = v;
    return view;
  }
}

QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
  
  if (j >= m->size2)
    {
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, view);
    }

  {
    TYPE(gsl_vector) v = NULL_VECTOR;
    
    v.data = m->data + j * MULTIPLICITY;
    v.size = m->size1;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    ((VIEW(_gsl_vector, view) *)&view)->vector = v;
    return view;
  }
}

QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;

  TYPE(gsl_vector) v = NULL_VECTOR;
  v.data = m->data;
  v.size = GSL_MIN(m->size1,m->size2);
  v.stride = m->tda + 1;
  v.block = m->block;
  v.owner = 0;

  ((VIEW(_gsl_vector, view) *)&view)->vector = v;
  return view;
}

QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
                                    const size_t k)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
  
  if (k >= m->size1)
    {
      GSL_ERROR_VAL ("subdiagonal index is out of range", GSL_EINVAL, view);
    }

  {
    TYPE(gsl_vector) v = NULL_VECTOR;
    
    v.data = m->data + k * MULTIPLICITY * m->tda;
    v.size = GSL_MIN(m->size1 - k, m->size2);
    v.stride = m->tda + 1;
    v.block = m->block;
    v.owner = 0;
    
    ((VIEW(_gsl_vector, view) *)&view)->vector = v;
    return view;
  }
}

QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
                                      const size_t k)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;


  if (k >= m->size2)
    {
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, view);
    }

  {
    TYPE(gsl_vector) v = NULL_VECTOR;
    
    v.data = m->data + k * MULTIPLICITY;
    v.size = GSL_MIN(m->size1, m->size2 - k);
    v.stride = m->tda + 1;
    v.block = m->block;
    v.owner = 0;

    ((VIEW(_gsl_vector, view) *)&view)->vector = v;
    return view;
  }
}

