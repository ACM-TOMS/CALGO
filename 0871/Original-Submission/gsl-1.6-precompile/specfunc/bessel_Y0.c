#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_Y0.c
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
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "bessel.h"
#include "bessel_amp_phase.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besy0, 1980 version, w. fullerton */

/* chebyshev expansions

 series for by0        on the interval  0.          to  1.60000d+01
                                        with weighted error   1.20e-17
                                         log weighted error  16.92
                               significant figures required  16.15
                                    decimal places required  17.48
*/

static MpIeee by0_data[13] =  {
  -MpIeee( "0.011277839392865573" ),
  -MpIeee( "0.128345237560420350" ),
  -MpIeee( "0.104378847997942490" ),
   MpIeee( "0.023662749183969695" ),
  -MpIeee( "0.002090391647700486" ),
   MpIeee( "0.000103975453939057" ),
  -MpIeee( "0.000003369747162423" ),
   MpIeee( "0.000000077293842676" ),
  -MpIeee( "0.000000001324976772" ),
   MpIeee( "0.000000000017648232" ),
  -MpIeee( "0.000000000000188105" ),
   MpIeee( "0.000000000000001641" ),
  -MpIeee( "0.000000000000000011" )
};
static cheb_series by0_cs = {
  by0_data,
  12,
  -1, 1,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_bessel_Y0_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee two_over_pi=  2.0/M_PI;
  const MpIeee xmax=  1.0/GSL_DBL_EPSILON;

  /* CHECK_POINTER(result) */

  if (x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 4.0) {
    gsl_sf_result J0;
    gsl_sf_result c;
    int  stat_J0=  gsl_sf_bessel_J0_e(x, &J0);
    cheb_eval_e(&by0_cs, 0.125*x*x-1.0, &c);
    result->val = two_over_pi*(-M_LN2 + log(x))*J0.val + 0.375 + c.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + c.err;
    return stat_J0;
  }
  else if(x < xmax) {
    /* Leading behaviour of phase is x, which is exact,
     * so the error is bounded.
     */
    const MpIeee z=  32.0/(x*x) - 1.0;
    gsl_sf_result c1;
    gsl_sf_result c2;
    gsl_sf_result sp;
    const int stat_c1 = cheb_eval_e(&_gsl_sf_bessel_amp_phase_bm0_cs,  z, &c1);
    const int stat_c2 = cheb_eval_e(&_gsl_sf_bessel_amp_phase_bth0_cs, z, &c2);
    const int stat_sp = gsl_sf_bessel_sin_pi4_e(x, c2.val/x, &sp);
    const MpIeee sqrtx=  sqrt(x);
    const MpIeee ampl=  (0.75 + c1.val) / sqrtx;
    result->val  = ampl * sp.val;
    result->err  = fabs(sp.val) * c1.err/sqrtx + fabs(ampl) * sp.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_3(stat_sp, stat_c1, stat_c2);
  }
  else {
    UNDERFLOW_ERROR(result);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_Y0(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_Y0_e(x, &result));
}
