#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/log.c
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
#include <gsl/gsl_sf_log.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev expansion for log(1 + x(t))/x(t)
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static MpIeee lopx_data[21] =  {
  MpIeee( "2.16647910664395270521272590407" ),
 -MpIeee( "0.28565398551049742084877469679" ),
  MpIeee( "0.01517767255690553732382488171" ),
 -MpIeee( "0.00200215904941415466274422081" ),
  MpIeee( "0.00019211375164056698287947962" ),
 -MpIeee( "0.00002553258886105542567601400" ),
  MpIeee( "2.9004512660400621301999384544e-06" ),
 -MpIeee( "3.8873813517057343800270917900e-07" ),
  MpIeee( "4.7743678729400456026672697926e-08" ),
 -MpIeee( "6.4501969776090319441714445454e-09" ),
  MpIeee( "8.2751976628812389601561347296e-10" ),
 -MpIeee( "1.1260499376492049411710290413e-10" ),
  MpIeee( "1.4844576692270934446023686322e-11" ),
 -MpIeee( "2.0328515972462118942821556033e-12" ),
  MpIeee( "2.7291231220549214896095654769e-13" ),
 -MpIeee( "3.7581977830387938294437434651e-14" ),
  MpIeee( "5.1107345870861673561462339876e-15" ),
 -MpIeee( "7.0722150011433276578323272272e-16" ),
  MpIeee( "9.7089758328248469219003866867e-17" ),
 -MpIeee( "1.3492637457521938883731579510e-17" ),
  MpIeee( "1.8657327910677296608121390705e-18" )
};
static cheb_series lopx_cs = {
  lopx_data,
  20,
  -1, 1,
  10
};

/* Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static MpIeee lopxmx_data[20] =  {
 -MpIeee( "1.12100231323744103373737274541" ),
  MpIeee( "0.19553462773379386241549597019" ),
 -MpIeee( "0.01467470453808083971825344956" ),
  MpIeee( "0.00166678250474365477643629067" ),
 -MpIeee( "0.00018543356147700369785746902" ),
  MpIeee( "0.00002280154021771635036301071" ),
 -MpIeee( "2.8031253116633521699214134172e-06" ),
  MpIeee( "3.5936568872522162983669541401e-07" ),
 -MpIeee( "4.6241857041062060284381167925e-08" ),
  MpIeee( "6.0822637459403991012451054971e-09" ),
 -MpIeee( "8.0339824424815790302621320732e-10" ),
  MpIeee( "1.0751718277499375044851551587e-10" ),
 -MpIeee( "1.4445310914224613448759230882e-11" ),
  MpIeee( "1.9573912180610336168921438426e-12" ),
 -MpIeee( "2.6614436796793061741564104510e-13" ),
  MpIeee( "3.6402634315269586532158344584e-14" ),
 -MpIeee( "4.9937495922755006545809120531e-15" ),
  MpIeee( "6.8802890218846809524646902703e-16" ),
 -MpIeee( "9.5034129794804273611403251480e-17" ),
  MpIeee( "1.3170135013050997157326965813e-17" )
};
static cheb_series lopxmx_cs = {
  lopxmx_data,
  19,
  -1, 1,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

#ifndef HIDE_INLINE_STATIC
int
 gsl_sf_log_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_log_abs_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}
#endif

int
 gsl_sf_complex_log_e(const MpIeee zr, const MpIeee zi, gsl_sf_result * lnr, gsl_sf_result * theta)
{
  /* CHECK_POINTER(lnr) */
  /* CHECK_POINTER(theta) */

  if(zr != 0.0 || zi != 0.0) {
    const MpIeee ax=  fabs(zr);
    const MpIeee ay=  fabs(zi);
    const MpIeee min=  GSL_MIN(ax, ay);
    const MpIeee max=  GSL_MAX(ax, ay);
    lnr->val = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
    lnr->err = 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
    theta->val = atan2(zi, zr);
    theta->err = GSL_DBL_EPSILON * fabs(lnr->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR_2(lnr, theta);
  }
}


int
 gsl_sf_log_1plusx_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < GSL_ROOT6_DBL_EPSILON) {
    const MpIeee c1=  -0.5;
    const MpIeee c2=   1.0/3.0;
    const MpIeee c3=  -1.0/4.0;
    const MpIeee c4=   1.0/5.0;
    const MpIeee c5=  -1.0/6.0;
    const MpIeee c6=   1.0/7.0;
    const MpIeee c7=  -1.0/8.0;
    const MpIeee c8=   1.0/9.0;
    const MpIeee c9=  -1.0/10.0;
    const MpIeee t=   c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
    result->val = x * (1.0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*t)))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    MpIeee t=  MpIeee( "0.5" )*(MpIeee( "8.0" )*x + MpIeee( "1.0" ))/(x+MpIeee( "2.0" ));
    gsl_sf_result c;
    cheb_eval_e(&lopx_cs, t, &c);
    result->val = x * c.val;
    result->err = fabs(x * c.err);
    return GSL_SUCCESS;
  }
  else {
    result->val = log(1.0 + x);
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_log_1plusx_mx_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < GSL_ROOT5_DBL_EPSILON) {
    const MpIeee c1=  -0.5;
    const MpIeee c2=   1.0/3.0;
    const MpIeee c3=  -1.0/4.0;
    const MpIeee c4=   1.0/5.0;
    const MpIeee c5=  -1.0/6.0;
    const MpIeee c6=   1.0/7.0;
    const MpIeee c7=  -1.0/8.0;
    const MpIeee c8=   1.0/9.0;
    const MpIeee c9=  -1.0/10.0;
    const MpIeee t=   c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
    result->val = x*x * (c1 + x*(c2 + x*(c3 + x*(c4 + x*t))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    MpIeee t=  MpIeee( "0.5" )*(MpIeee( "8.0" )*x + MpIeee( "1.0" ))/(x+MpIeee( "2.0" ));
    gsl_sf_result c;
    cheb_eval_e(&lopxmx_cs, t, &c);
    result->val = x*x * c.val;
    result->err = x*x * c.err;
    return GSL_SUCCESS;
  }
  else {
    const MpIeee lterm=  log(1.0 + x);
    result->val = lterm - x;
    result->err = GSL_DBL_EPSILON * (fabs(lterm) + fabs(x));
    return GSL_SUCCESS;
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_log(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_log_e(x, &result));
}

MpIeee gsl_sf_log_abs(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_log_abs_e(x, &result));
}

MpIeee gsl_sf_log_1plusx(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_e(x, &result));
}

MpIeee gsl_sf_log_1plusx_mx(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_mx_e(x, &result));
}
