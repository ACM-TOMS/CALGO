#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* min/convergence.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

int
 gsl_min_test_interval(MpIeee x_lower, MpIeee x_upper, MpIeee epsabs, MpIeee epsrel)
{
  const MpIeee lower=  x_lower;
  const MpIeee upper=  x_upper;

  const MpIeee abs_lower=  fabs(lower) ;
  const MpIeee abs_upper=  fabs(upper) ;

  MpIeee min_abs;MpIeee  tolerance;

  if (epsrel < MpIeee( "0.0" ))
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);
  
  if (epsabs < MpIeee( "0.0" ))
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  if (lower > upper)
    GSL_ERROR ("lower bound larger than upper_bound", GSL_EINVAL);

  if ((lower > 0 && upper > 0) || (lower < 0 && upper < 0)) 
    {
      min_abs = GSL_MIN_DBL(abs_lower, abs_upper) ;
    }
  else
    {
      min_abs = MpIeee( "0" );
    }

  tolerance = epsabs + epsrel * min_abs  ;
  
  if (fabs(upper - lower) < tolerance)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

