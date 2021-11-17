#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_I1.c
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
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#define ROOT_EIGHT (2.0*M_SQRT2)


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besi1(), besi1e() */

/* chebyshev expansions

 series for bi1        on the interval  0.          to  9.00000d+00
                                        with weighted error   2.40e-17
                                         log weighted error  16.62
                               significant figures required  16.23
                                    decimal places required  17.14

 series for ai1        on the interval  1.25000d-01 to  3.33333d-01
                                        with weighted error   6.98e-17
                                         log weighted error  16.16
                               significant figures required  14.53
                                    decimal places required  16.82

 series for ai12       on the interval  0.          to  1.25000d-01
                                       with weighted error   3.55e-17
                                        log weighted error  16.45
                              significant figures required  14.69
                                   decimal places required  17.12
*/

static MpIeee bi1_data[11] =  {
  -MpIeee( "0.001971713261099859" ),
   MpIeee( "0.407348876675464810" ),
   MpIeee( "0.034838994299959456" ),
   MpIeee( "0.001545394556300123" ),
   MpIeee( "0.000041888521098377" ),
   MpIeee( "0.000000764902676483" ),
   MpIeee( "0.000000010042493924" ),
   MpIeee( "0.000000000099322077" ),
   MpIeee( "0.000000000000766380" ),
   MpIeee( "0.000000000000004741" ),
   MpIeee( "0.000000000000000024" )
};
static cheb_series bi1_cs = {
  bi1_data,
  10,
  -1, 1,
  10
};

static MpIeee ai1_data[21] =  {
  -MpIeee( "0.02846744181881479" ),
  -MpIeee( "0.01922953231443221" ),
  -MpIeee( "0.00061151858579437" ),
  -MpIeee( "0.00002069971253350" ),
   MpIeee( "0.00000858561914581" ),
   MpIeee( "0.00000104949824671" ),
  -MpIeee( "0.00000029183389184" ),
  -MpIeee( "0.00000001559378146" ),
   MpIeee( "0.00000001318012367" ),
  -MpIeee( "0.00000000144842341" ),
  -MpIeee( "0.00000000029085122" ),
   MpIeee( "0.00000000012663889" ),
  -MpIeee( "0.00000000001664947" ),
  -MpIeee( "0.00000000000166665" ),
   MpIeee( "0.00000000000124260" ),
  -MpIeee( "0.00000000000027315" ),
   MpIeee( "0.00000000000002023" ),
   MpIeee( "0.00000000000000730" ),
  -MpIeee( "0.00000000000000333" ),
   MpIeee( "0.00000000000000071" ),
  -MpIeee( "0.00000000000000006" )
};
static cheb_series ai1_cs = {
  ai1_data,
  20,
  -1, 1,
  11
};

static MpIeee ai12_data[22] =  {
   MpIeee( "0.02857623501828014" ),
  -MpIeee( "0.00976109749136147" ),
  -MpIeee( "0.00011058893876263" ),
  -MpIeee( "0.00000388256480887" ),
  -MpIeee( "0.00000025122362377" ),
  -MpIeee( "0.00000002631468847" ),
  -MpIeee( "0.00000000383538039" ),
  -MpIeee( "0.00000000055897433" ),
  -MpIeee( "0.00000000001897495" ),
   MpIeee( "0.00000000003252602" ),
   MpIeee( "0.00000000001412580" ),
   MpIeee( "0.00000000000203564" ),
  -MpIeee( "0.00000000000071985" ),
  -MpIeee( "0.00000000000040836" ),
  -MpIeee( "0.00000000000002101" ),
   MpIeee( "0.00000000000004273" ),
   MpIeee( "0.00000000000001041" ),
  -MpIeee( "0.00000000000000382" ),
  -MpIeee( "0.00000000000000186" ),
   MpIeee( "0.00000000000000033" ),
   MpIeee( "0.00000000000000028" ),
  -MpIeee( "0.00000000000000003" )
};
static cheb_series ai12_cs = {
  ai12_data,
  21,
  -1, 1,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_bessel_I1_scaled_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee xmin=  2.0 * GSL_DBL_MIN;
  const MpIeee x_small=  ROOT_EIGHT * GSL_SQRT_DBL_EPSILON;
  const MpIeee y=  fabs(x);

  /* CHECK_POINTER(result) */

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(y < xmin) {
    UNDERFLOW_ERROR(result);
  }
  else if(y < x_small) {
    result->val = 0.5*x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    const MpIeee ey=  exp(-y);
    gsl_sf_result c;
    cheb_eval_e(&bi1_cs, y*y/4.5-1.0, &c);
    result->val  = x * ey * (0.875 + c.val);
    result->err  = ey * c.err + y * GSL_DBL_EPSILON * fabs(result->val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(y <= 8.0) {
    const MpIeee sy=  sqrt(y);
    gsl_sf_result c;
    MpIeee b;
    MpIeee s;
    cheb_eval_e(&ai1_cs, (48.0/y-11.0)/5.0, &c);
    b = (MpIeee( "0.375" ) + c.val) / sy;
    s = (x > MpIeee( "0.0" ) ? MpIeee( "1.0" ) : -MpIeee( "1.0" ));
    result->val  = s * b;
    result->err  = c.err / sy;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee sy=  sqrt(y);
    gsl_sf_result c;
    MpIeee b;
    MpIeee s;
    cheb_eval_e(&ai12_cs, 16.0/y-1.0, &c);
    b = (MpIeee( "0.375" ) + c.val) / sy;
    s = (x > MpIeee( "0.0" ) ? MpIeee( "1.0" ) : -MpIeee( "1.0" ));
    result->val  = s * b;
    result->err  = c.err / sy;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_bessel_I1_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee xmin=  2.0 * GSL_DBL_MIN;
  const MpIeee x_small=  ROOT_EIGHT * GSL_SQRT_DBL_EPSILON;
  const MpIeee y=  fabs(x);

  /* CHECK_POINTER(result) */

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(y < xmin) {
    UNDERFLOW_ERROR(result);
  }
  else if(y < x_small) {
    result->val = 0.5*x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    gsl_sf_result c;
    cheb_eval_e(&bi1_cs, y*y/4.5-1.0, &c);
    result->val  = x * (0.875 + c.val);
    result->err  = y * c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(y < GSL_LOG_DBL_MAX) {
    const MpIeee ey=  exp(y);
    gsl_sf_result I1_scaled;
    gsl_sf_bessel_I1_scaled_e(x, &I1_scaled);
    result->val  = ey * I1_scaled.val;
    result->err  = ey * I1_scaled.err + y * GSL_DBL_EPSILON * fabs(result->val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_I1_scaled(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_I1_scaled_e(x, &result));
}

MpIeee gsl_sf_bessel_I1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_I1_e(x, &result));
}
