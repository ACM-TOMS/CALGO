#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/pow_int.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_pow_int.h>


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error handling *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_pow_int_e(MpIeee x, int  n, gsl_sf_result * result)
{
  MpIeee value=  MpIeee( "1.0" );
  int  count=  0;

  /* CHECK_POINTER(result) */


  if(n < 0) {
    n = -n;

    if(x == MpIeee( "0.0" )) {
      MpIeee u=  MpIeee( "1.0" ) / x;
      result->val = (n % 2) ? u : (u * u) ;  /* correct sign of infinity */
      result->err = GSL_POSINF;
      GSL_ERROR ("overflow", GSL_EOVRFLW);
    }

    x = MpIeee( "1.0" )/x;
  }

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(GSL_IS_ODD(n)) value *= x;
     n >>= 1;
     x *= x;
     ++count;
  } while (n);

  result->val = value;
  result->err = 2.0 * GSL_DBL_EPSILON * (count + 1.0) * fabs(value); 

  return GSL_SUCCESS;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_pow_int(const MpIeee x, const int n)
{
  EVAL_RESULT(gsl_sf_pow_int_e(x, n, &result));
}
