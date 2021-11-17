#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_K1.c
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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besk1(), besk1e() */

/* chebyshev expansions 

 series for bk1        on the interval  0.          to  4.00000d+00
                                        with weighted error   7.02e-18
                                         log weighted error  17.15
                               significant figures required  16.73
                                    decimal places required  17.67

 series for ak1        on the interval  1.25000d-01 to  5.00000d-01
                                        with weighted error   6.06e-17
                                         log weighted error  16.22
                               significant figures required  15.41
                                    decimal places required  16.83

 series for ak12       on the interval  0.          to  1.25000d-01
                                        with weighted error   2.58e-17
                                         log weighted error  16.59
                               significant figures required  15.22
                                    decimal places required  17.16
*/

static MpIeee bk1_data[11] =  {
   MpIeee( "0.0253002273389477705" ),
  -MpIeee( "0.3531559607765448760" ), 
  -MpIeee( "0.1226111808226571480" ), 
  -MpIeee( "0.0069757238596398643" ),
  -MpIeee( "0.0001730288957513052" ),
  -MpIeee( "0.0000024334061415659" ),
  -MpIeee( "0.0000000221338763073" ),
  -MpIeee( "0.0000000001411488392" ),
  -MpIeee( "0.0000000000006666901" ),
  -MpIeee( "0.0000000000000024274" ),
  -MpIeee( "0.0000000000000000070" )
};

static cheb_series bk1_cs = {
  bk1_data,
  10,
  -1, 1,
  8
};

static MpIeee ak1_data[17] =  {
   MpIeee( "0.27443134069738830" ), 
   MpIeee( "0.07571989953199368" ),
  -MpIeee( "0.00144105155647540" ),
   MpIeee( "0.00006650116955125" ),
  -MpIeee( "0.00000436998470952" ),
   MpIeee( "0.00000035402774997" ),
  -MpIeee( "0.00000003311163779" ),
   MpIeee( "0.00000000344597758" ),
  -MpIeee( "0.00000000038989323" ),
   MpIeee( "0.00000000004720819" ),
  -MpIeee( "0.00000000000604783" ),
   MpIeee( "0.00000000000081284" ),
  -MpIeee( "0.00000000000011386" ),
   MpIeee( "0.00000000000001654" ),
  -MpIeee( "0.00000000000000248" ),
   MpIeee( "0.00000000000000038" ),
  -MpIeee( "0.00000000000000006" )
};
static cheb_series ak1_cs = {
  ak1_data,
  16,
  -1, 1,
  9
};

static MpIeee ak12_data[14] =  {
   MpIeee( "0.06379308343739001" ),
   MpIeee( "0.02832887813049721" ),
  -MpIeee( "0.00024753706739052" ),
   MpIeee( "0.00000577197245160" ),
  -MpIeee( "0.00000020689392195" ),
   MpIeee( "0.00000000973998344" ),
  -MpIeee( "0.00000000055853361" ),
   MpIeee( "0.00000000003732996" ),
  -MpIeee( "0.00000000000282505" ),
   MpIeee( "0.00000000000023720" ),
  -MpIeee( "0.00000000000002176" ),
   MpIeee( "0.00000000000000215" ),
  -MpIeee( "0.00000000000000022" ),
   MpIeee( "0.00000000000000002" )
};
static cheb_series ak12_cs = {
  ak12_data,
  13,
  -1, 1,
  7
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_bessel_K1_scaled_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    OVERFLOW_ERROR(result);
  }
  else if(x <= 2.0) {
    const MpIeee lx=  log(x);
    const MpIeee ex=  exp(x);
    int  stat_I1;
    gsl_sf_result I1;
    gsl_sf_result c;
    cheb_eval_e(&bk1_cs, 0.5*x*x-1.0, &c);
    stat_I1 = gsl_sf_bessel_I1_e(x, &I1);
    result->val  = ex * ((lx-M_LN2)*I1.val + (0.75 + c.val)/x);
    result->err  = ex * (c.err/x + fabs(lx)*I1.err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I1;
  }
  else if(x <= 8.0) {
    const MpIeee sx=  sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak1_cs, (16.0/x-5.0)/3.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee sx=  sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak12_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_bessel_K1_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    OVERFLOW_ERROR(result);
  }
  else if(x <= 2.0) {
    const MpIeee lx=  log(x);
    int  stat_I1;
    gsl_sf_result I1;
    gsl_sf_result c;
    cheb_eval_e(&bk1_cs, 0.5*x*x-1.0, &c);
    stat_I1 = gsl_sf_bessel_I1_e(x, &I1);
    result->val  = (lx-M_LN2)*I1.val + (0.75 + c.val)/x;
    result->err  = c.err/x + fabs(lx)*I1.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I1;
  }
  else {
    gsl_sf_result K1_scaled;
    int  stat_K1=  gsl_sf_bessel_K1_scaled_e(x, &K1_scaled);
    int  stat_e=  gsl_sf_exp_mult_err_e(-x, 0.0,
                                           K1_scaled.val, K1_scaled.err,
                                           result);
    result->err = fabs(result->val) * (GSL_DBL_EPSILON*fabs(x) + K1_scaled.err/K1_scaled.val);
    return GSL_ERROR_SELECT_2(stat_e, stat_K1);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_K1_scaled(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_K1_scaled_e(x, &result));
}

MpIeee gsl_sf_bessel_K1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_K1_e(x, &result));
}
