#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/convergence.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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
#include <gsl/gsl_roots.h>

int
 gsl_root_test_interval(MpIeee x_lower, MpIeee x_upper, MpIeee epsabs, MpIeee epsrel)
{
  const MpIeee abs_lower=  fabs(x_lower) ;
  const MpIeee abs_upper=  fabs(x_upper) ;

  MpIeee min_abs;MpIeee  tolerance;

  if (epsrel < MpIeee( "0.0" ))
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);
  
  if (epsabs < MpIeee( "0.0" ))
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  if (x_lower > x_upper)
    GSL_ERROR ("lower bound larger than upper bound", GSL_EINVAL);

  if ((x_lower > MpIeee( "0.0" ) && x_upper > MpIeee( "0.0" )) || (x_lower < MpIeee( "0.0" ) && x_upper < MpIeee( "0.0" ))) 
    {
      min_abs = GSL_MIN_DBL(abs_lower, abs_upper) ;
    }
  else
    {
      min_abs = MpIeee( "0" );
    }

  tolerance = epsabs + epsrel * min_abs  ;
  
  if (fabs(x_upper - x_lower) < tolerance)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

int
 gsl_root_test_delta(MpIeee x1, MpIeee x0, MpIeee epsabs, MpIeee epsrel)
{
  const MpIeee tolerance=  epsabs + epsrel * fabs(x1)  ;

  if (epsrel < MpIeee( "0.0" ))
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);
  
  if (epsabs < MpIeee( "0.0" ))
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
  
  if (fabs(x1 - x0) < tolerance)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

int
 gsl_root_test_residual(MpIeee f, MpIeee epsabs)
{
  if (epsabs < MpIeee( "0.0" ))
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
 
  if (fabs(f) < epsabs)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

