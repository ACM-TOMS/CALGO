#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multimin/convergence.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
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
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

int
 gsl_multimin_test_gradient(const gsl_vector *g, MpIeee epsabs)
{
  MpIeee norm;

  if (epsabs < MpIeee( "0.0" ))
    {
      GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
    }

  norm = gsl_blas_dnrm2(g);
  
  if (norm < epsabs)
    {
      return GSL_SUCCESS;
    }

  return GSL_CONTINUE;
}

int
 gsl_multimin_test_size(const MpIeee size, MpIeee epsabs)
{
  if (epsabs < MpIeee( "0.0" ))
    {
      GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
    }
  
  if (size < epsabs)
    {
      return GSL_SUCCESS;
    }

  return GSL_CONTINUE;
}
