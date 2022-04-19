#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* min/test_funcs.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include "test.h"

gsl_function create_function (MpIeee(*f)(MpIeee, void *)) 
{
  gsl_function F ;
  F.function = f ;
  F.params = 0 ;
  return F ;
}

MpIeee f_cos(MpIeee x, void * p)
{
  p = 0;  /* avoid warning about unused parameter */
  return cos(x);
}

/* f(x) = x^4 - 1 */
/* minimum at x = 0 */

MpIeee func1(MpIeee x, void * p)
{
  p = 0;  /* avoid warning about unused parameter */
  return pow (x, MpIeee( "4.0" )) - MpIeee( "1" );
}

/* f(x) = sqrt(|x|) */
/* minimum at x = 0 */

MpIeee func2(MpIeee x, void * p)
{
  p = 0;  /* avoid warning about unused parameter */
  return sqrt(fabs(x));
}


/* f(x) = 1 for x < 1 and -exp(-x) for x >= 1 */
/* minimum at x = 1 */

MpIeee func3(MpIeee x, void * p)
{
  p = 0;  /* avoid warning about unused parameter */

  if (x < MpIeee( "1" ))
    return MpIeee( "1" ) ;
  else
    return - exp(-x) ;
}

/* f(x) = x - 30/(1+1e5*(x-0.8)**2) */
/* minimum near x = 0.8 */

MpIeee func4(MpIeee x, void * p)
{
  p = 0;  /* avoid warning about unused parameter */

  return x - MpIeee( "30.0" ) / (MpIeee( "1.0" ) + MpIeee( "1" )e5 * pow(x-MpIeee( "0.8" ), MpIeee( "2.0" )));
}

