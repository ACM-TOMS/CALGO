#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* matrix/matrix_source.c
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
MpIeee FUNCTION(gsl_matrix, get) (const TYPE (gsl_matrix) * m,
                            const size_t i, const size_t j)
{
  MpIeee zero=  ZERO;

  if (gsl_check_range)
    {
      if (i >= m->size1)        /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_VAL ("first index out of range", GSL_EINVAL, zero);
        }
      else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_VAL ("second index out of range", GSL_EINVAL, zero);
        }
    }
  
  //cout << " |get i="<<i<<" j="<<j<<" "<<flush;
  //cout << " get index="<<MULTIPLICITY * (i * m->tda + j)<<"| " <<flush;
  return *(MpIeee *) (m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION (gsl_matrix, set) (TYPE (gsl_matrix) * m,
                            const size_t i, const size_t j,
                            const BASE x)
{
  if (gsl_check_range)
    {
      if (i >= m->size1)        /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_VOID ("first index out of range", GSL_EINVAL);
          //cerr << "gslmatrix:first index out of range"<<endl;
        }
      else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_VOID ("second index out of range", GSL_EINVAL);
          //cerr << "gslmatrix:second index out of range"<<endl;
        }
    }
  *(MpIeee*) (m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


MpIeee *
FUNCTION(gsl_matrix, ptr) (TYPE (gsl_matrix) * m,
                            const size_t i, const size_t j)
{
  if (gsl_check_range)
    {
      if (i >= m->size1)        /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_NULL ("first index out of range", GSL_EINVAL);
          //cerr << "gslmatrix:first index out of range"<<endl;
        }
      else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_NULL ("second index out of range", GSL_EINVAL);
          //cerr << "gslmatrix:second index out of range"<<endl;
        }
    }
  return (MpIeee *) (m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION (gsl_matrix, const_ptr) (const TYPE (gsl_matrix) * m,
                                  const size_t i, const size_t j)
{
  if (gsl_check_range)
    {
      if (i >= m->size1)        /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_NULL ("first index out of range", GSL_EINVAL);
          //cerr << "gslmatrix:first index out of range"<<endl;

        }
      else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
        {
          //GSL_ERROR_NULL ("second index out of range", GSL_EINVAL);
          //cerr << "gslmatrix:first index out of range"<<endl;
        }
    }
  return (const BASE *) (m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

