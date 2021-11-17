#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_I0.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on SLATEC besi0 */

/* chebyshev expansions

 series for bi0        on the interval  0.          to  9.00000d+00
                                        with weighted error   2.46e-18
                                         log weighted error  17.61
                               significant figures required  17.90
                                    decimal places required  18.15

 series for ai0        on the interval  1.25000d-01 to  3.33333d-01
                                        with weighted error   7.87e-17
                                         log weighted error  16.10
                               significant figures required  14.69
                                    decimal places required  16.76


 series for ai02       on the interval  0.          to  1.25000d-01
                                        with weighted error   3.79e-17
                                         log weighted error  16.42
                               significant figures required  14.86
                                    decimal places required  17.09
*/

static MpIeee bi0_data[12] =  {
  -MpIeee( ".07660547252839144951" ),
  MpIeee( "1.92733795399380827000" ),
   MpIeee( ".22826445869203013390" ), 
   MpIeee( ".01304891466707290428" ),
   MpIeee( ".00043442709008164874" ),
   MpIeee( ".00000942265768600193" ),
   MpIeee( ".00000014340062895106" ),
   MpIeee( ".00000000161384906966" ),
   MpIeee( ".00000000001396650044" ),
   MpIeee( ".00000000000009579451" ),
   MpIeee( ".00000000000000053339" ),
   MpIeee( ".00000000000000000245" )
};
static cheb_series bi0_cs = {
  bi0_data,
  11,
  -1, 1,
  11
};

static MpIeee ai0_data[21] =  {
   MpIeee( ".07575994494023796" ), 
   MpIeee( ".00759138081082334" ),
   MpIeee( ".00041531313389237" ),
   MpIeee( ".00001070076463439" ),
  -MpIeee( ".00000790117997921" ),
  -MpIeee( ".00000078261435014" ),
   MpIeee( ".00000027838499429" ),
   MpIeee( ".00000000825247260" ),
  -MpIeee( ".00000001204463945" ),
   MpIeee( ".00000000155964859" ),
   MpIeee( ".00000000022925563" ),
  -MpIeee( ".00000000011916228" ),
   MpIeee( ".00000000001757854" ),
   MpIeee( ".00000000000112822" ),
  -MpIeee( ".00000000000114684" ),
   MpIeee( ".00000000000027155" ),
  -MpIeee( ".00000000000002415" ),
  -MpIeee( ".00000000000000608" ),
   MpIeee( ".00000000000000314" ),
  -MpIeee( ".00000000000000071" ),
   MpIeee( ".00000000000000007" )
};
static cheb_series ai0_cs = {
  ai0_data,
  20,
  -1, 1,
  13
};

static MpIeee ai02_data[22] =  {
   MpIeee( ".05449041101410882" ),
   MpIeee( ".00336911647825569" ),
   MpIeee( ".00006889758346918" ),
   MpIeee( ".00000289137052082" ),
   MpIeee( ".00000020489185893" ),
   MpIeee( ".00000002266668991" ),
   MpIeee( ".00000000339623203" ),
   MpIeee( ".00000000049406022" ),
   MpIeee( ".00000000001188914" ),
  -MpIeee( ".00000000003149915" ),
  -MpIeee( ".00000000001321580" ),
  -MpIeee( ".00000000000179419" ),
   MpIeee( ".00000000000071801" ),
   MpIeee( ".00000000000038529" ),
   MpIeee( ".00000000000001539" ),
  -MpIeee( ".00000000000004151" ),
  -MpIeee( ".00000000000000954" ),
   MpIeee( ".00000000000000382" ),
   MpIeee( ".00000000000000176" ),
  -MpIeee( ".00000000000000034" ),
  -MpIeee( ".00000000000000027" ),
   MpIeee( ".00000000000000003" )
};
static cheb_series ai02_cs = {
  ai02_data,
  21,
  -1, 1,
  11
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_bessel_I0_scaled_e(const MpIeee x, gsl_sf_result * result)
{
  MpIeee y=  fabs(x);

  /* CHECK_POINTER(result) */

  if(y < MpIeee( "2.0" ) * GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - y;
    result->err = 0.5*y*y;
    return GSL_SUCCESS;
  }
  else if(y <= MpIeee( "3.0" )) {
    const MpIeee ey=  exp(-y);
    gsl_sf_result c;
    cheb_eval_e(&bi0_cs, y*y/4.5-1.0, &c);
    result->val = ey * (2.75 + c.val);
    result->err = GSL_DBL_EPSILON * fabs(result->val) + ey * c.err;
    return GSL_SUCCESS;
  }
  else if(y <= MpIeee( "8.0" )) {
    const MpIeee sy=  sqrt(y);
    gsl_sf_result c;
    cheb_eval_e(&ai0_cs, (48.0/y-11.0)/5.0, &c);
    result->val  = (0.375 + c.val) / sy;
    result->err  = 2.0 * GSL_DBL_EPSILON * (0.375 + fabs(c.val)) / sy;
    result->err += c.err / sy;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee sy=  sqrt(y);
    gsl_sf_result c;
    cheb_eval_e(&ai02_cs, 16.0/y-1.0, &c);
    result->val = (0.375 + c.val) / sy;
    result->err  = 2.0 * GSL_DBL_EPSILON * (0.375 + fabs(c.val)) / sy;
    result->err += c.err / sy;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_bessel_I0_e(const MpIeee x, gsl_sf_result * result)
{
  MpIeee y=  fabs(x);

  /* CHECK_POINTER(result) */

  if(y < MpIeee( "2.0" ) * GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0;
    result->err = 0.5*y*y;
    return GSL_SUCCESS;
  }
  else if(y <= MpIeee( "3.0" )) {
    gsl_sf_result c;
    cheb_eval_e(&bi0_cs, y*y/4.5-1.0, &c);
    result->val  = 2.75 + c.val;
    result->err  = GSL_DBL_EPSILON * (2.75 + fabs(c.val));
    result->err += c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(y < GSL_LOG_DBL_MAX - MpIeee( "1.0" )) {
    const MpIeee ey=  exp(y);
    gsl_sf_result b_scaled;
    gsl_sf_bessel_I0_scaled_e(x, &b_scaled);
    result->val  = ey * b_scaled.val;
    result->err  = ey * b_scaled.err + y*GSL_DBL_EPSILON*fabs(result->val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_I0_scaled(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_I0_scaled_e(x, &result); ) 
}

MpIeee gsl_sf_bessel_I0(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_I0_e(x, &result); ) 
}
