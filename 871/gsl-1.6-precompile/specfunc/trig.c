#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/trig.c
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
#include <gsl/gsl_sf_trig.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/* sinh(x) series
 * double-precision for |x| < 1.0
 */
inline
static
int
 sinh_series(const MpIeee x, MpIeee * result)
{
  const MpIeee y=  x*x;
  const MpIeee c0=  1.0/6.0;
  const MpIeee c1=  1.0/120.0;
  const MpIeee c2=  1.0/5040.0;
  const MpIeee c3=  1.0/362880.0;
  const MpIeee c4=  1.0/39916800.0;
  const MpIeee c5=  1.0/6227020800.0;
  const MpIeee c6=  1.0/1307674368000.0;
  const MpIeee c7=  1.0/355687428096000.0;
  *result = x*(MpIeee( "1.0" ) + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))));
  return GSL_SUCCESS;
}


/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
inline
static
int
 cosh_m1_series(const MpIeee x, MpIeee * result)
{
  const MpIeee y=  x*x;
  const MpIeee c0=  0.5;
  const MpIeee c1=  1.0/24.0;
  const MpIeee c2=  1.0/720.0;
  const MpIeee c3=  1.0/40320.0;
  const MpIeee c4=  1.0/3628800.0;
  const MpIeee c5=  1.0/479001600.0;
  const MpIeee c6=  1.0/87178291200.0;
  const MpIeee c7=  1.0/20922789888000.0;
  const MpIeee c8=  1.0/6402373705728000.0;
  *result = y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))));
  return GSL_SUCCESS;
}


/* Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
 */
static MpIeee sinc_data[17] =  {
  MpIeee( "1.133648177811747875422" ),
 -MpIeee( "0.532677564732557348781" ),
 -MpIeee( "0.068293048346633177859" ),
  MpIeee( "0.033403684226353715020" ),
  MpIeee( "0.001485679893925747818" ),
 -MpIeee( "0.000734421305768455295" ),
 -MpIeee( "0.000016837282388837229" ),
  MpIeee( "0.000008359950146618018" ),
  MpIeee( "0.000000117382095601192" ),
 -MpIeee( "0.000000058413665922724" ),
 -MpIeee( "0.000000000554763755743" ),
  MpIeee( "0.000000000276434190426" ),
  MpIeee( "0.000000000001895374892" ),
 -MpIeee( "0.000000000000945237101" ),
 -MpIeee( "0.000000000000004900690" ),
  MpIeee( "0.000000000000002445383" ),
  MpIeee( "0.000000000000000009925" )
};
static cheb_series sinc_cs = {
  sinc_data,
  16,
  -1, 1,
  10
};


/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (sin(x)/x - 1)/(x*x)
 */
static MpIeee sin_data[12] =  {
  -MpIeee( "0.3295190160663511504173" ),
   MpIeee( "0.0025374284671667991990" ),
   MpIeee( "0.0006261928782647355874" ),
  -MpIeee( "4.6495547521854042157541e-06" ),
  -MpIeee( "5.6917531549379706526677e-07" ),
   MpIeee( "3.7283335140973803627866e-09" ),
   MpIeee( "3.0267376484747473727186e-10" ),
  -MpIeee( "1.7400875016436622322022e-12" ),
  -MpIeee( "1.0554678305790849834462e-13" ),
   MpIeee( "5.3701981409132410797062e-16" ),
   MpIeee( "2.5984137983099020336115e-17" ),
  -MpIeee( "1.1821555255364833468288e-19" )
};
static cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  11
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
 */
static MpIeee cos_data[11] =  {
  MpIeee( "0.165391825637921473505668118136" ),
 -MpIeee( "0.00084852883845000173671196530195" ),
 -MpIeee( "0.000210086507222940730213625768083" ),
  MpIeee( "1.16582269619760204299639757584e-6" ),
  MpIeee( "1.43319375856259870334412701165e-7" ),
 -MpIeee( "7.4770883429007141617951330184e-10" ),
 -MpIeee( "6.0969994944584252706997438007e-11" ),
  MpIeee( "2.90748249201909353949854872638e-13" ),
  MpIeee( "1.77126739876261435667156490461e-14" ),
 -MpIeee( "7.6896421502815579078577263149e-17" ),
 -MpIeee( "3.7363121133079412079201377318e-18" )
};
static cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* I would have prefered just using the library sin() function.
 * But after some experimentation I decided that there was
 * no good way to understand the error; library sin() is just a black box.
 * So we have to roll our own.
 */
int
 gsl_sf_sin_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const MpIeee P1=  7.85398125648498535156e-1;
    const MpIeee P2=  3.77489470793079817668e-8;
    const MpIeee P3=  2.69515142907905952645e-15;

    const MpIeee sgn_x=  GSL_SIGN(x);
    const MpIeee abs_x=  fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const MpIeee x2=  x*x;
      result->val = x * (1.0 - x2/6.0);
      result->err = fabs(x*x2*x2 / 100.0);
      return GSL_SUCCESS;
    }
    else {
      MpIeee sgn_result=  sgn_x;
      MpIeee y=  floor(abs_x/(MpIeee( "0.25" )*M_PI));
      int  octant=  y - ldexp(floor(ldexp(y,-3)),3);
      int  stat_cs;
      MpIeee z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += MpIeee( "1.0" );
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }
      
      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result sin_cs_result;
        const MpIeee t=  8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result cos_cs_result;
        const MpIeee t=  8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
 gsl_sf_cos_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const MpIeee P1=  7.85398125648498535156e-1;
    const MpIeee P2=  3.77489470793079817668e-8;
    const MpIeee P3=  2.69515142907905952645e-15;

    const MpIeee abs_x=  fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const MpIeee x2=  x*x;
      result->val = 1.0 - 0.5*x2;
      result->err = fabs(x2*x2/12.0);
      return GSL_SUCCESS;
    }
    else {
      MpIeee sgn_result=  MpIeee( "1.0" );
      MpIeee y=  floor(abs_x/(MpIeee( "0.25" )*M_PI));
      int  octant=  y - ldexp(floor(ldexp(y,-3)),3);
      int  stat_cs;
      MpIeee z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += MpIeee( "1.0" );
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }

      if(octant > 1) {
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result cos_cs_result;
        const MpIeee t=  8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result sin_cs_result;
        const MpIeee t=  8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
 gsl_sf_hypot_e(const MpIeee x, const MpIeee y, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0 && y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const MpIeee a=  fabs(x);
    const MpIeee b=  fabs(y);
    const MpIeee min=  GSL_MIN_DBL(a,b);
    const MpIeee max=  GSL_MAX_DBL(a,b);
    const MpIeee rat=  min/max;
    const MpIeee root_term=  sqrt(1.0 + rat*rat);

    if(max < GSL_DBL_MAX/root_term) {
      result->val = max * root_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
}


int
 gsl_sf_complex_sin_e(const MpIeee zr, const MpIeee zi,
                        gsl_sf_result * szr, gsl_sf_result * szi)
{
  /* CHECK_POINTER(szr) */
  /* CHECK_POINTER(szi) */

  if(fabs(zi) < 1.0) {
    MpIeee ch_m1;MpIeee  sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    szr->val = sin(zr)*(ch_m1 + 1.0);
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    MpIeee ex=  exp(zi);
    MpIeee ch=  MpIeee( "0.5" )*(ex+MpIeee( "1.0" )/ex);
    MpIeee sh=  MpIeee( "0.5" )*(ex-MpIeee( "1.0" )/ex);
    szr->val = sin(zr)*ch;
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(szr, szi);
  }
}


int
 gsl_sf_complex_cos_e(const MpIeee zr, const MpIeee zi,
                        gsl_sf_result * czr, gsl_sf_result * czi)
{
  /* CHECK_POINTER(czr) */
  /* CHECK_POINTER(czi) */

  if(fabs(zi) < 1.0) {
    MpIeee ch_m1;MpIeee  sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    czr->val =  cos(zr)*(ch_m1 + 1.0);
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    MpIeee ex=  exp(zi);
    MpIeee ch=  MpIeee( "0.5" )*(ex+MpIeee( "1.0" )/ex);
    MpIeee sh=  MpIeee( "0.5" )*(ex-MpIeee( "1.0" )/ex);
    czr->val =  cos(zr)*ch;
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(czr,czi);
  }
}


int
 gsl_sf_complex_logsin_e(const MpIeee zr, const MpIeee zi,
                           gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  /* CHECK_POINTER(lszr) */
  /* CHECK_POINTER(lszi) */

  if(zi > 60.0) {
    lszr->val = -M_LN2 + zi;
    lszi->val =  0.5*M_PI - zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else {
    gsl_sf_result sin_r, sin_i;
    int  status;
    gsl_sf_complex_sin_e(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_e(sin_r.val, sin_i.val, lszr, lszi);
    if(status == GSL_EDOM) {
      DOMAIN_ERROR_2(lszr, lszi);
    }
  }
  return gsl_sf_angle_restrict_symm_e(&(lszi->val));
}


int
 gsl_sf_lnsinh_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 1.0) {
    MpIeee eps;
    sinh_series(x, &eps);
    result->val = log(eps);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_lncosh_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    MpIeee eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log_1plusx_e(eps, result);
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 + exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*
inline int gsl_sf_sincos_e(const double theta, double * s, double * c)
{
  double tan_half = tan(0.5 * theta);
  double den = 1. + tan_half*tan_half;
  double cos_theta = (1.0 - tan_half*tan_half) / den;
  double sin_theta = 2.0 * tan_half / den;
}
*/

int
 gsl_sf_polar_to_rect(const MpIeee r, const MpIeee theta,
                          gsl_sf_result * x, gsl_sf_result * y)
{
  MpIeee t=  theta;
  int  status=  gsl_sf_angle_restrict_symm_e(&t);
  MpIeee c=  cos(t);
  MpIeee s=  sin(t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err  = r * fabs(s * GSL_DBL_EPSILON * t);
  x->err += 2.0 * GSL_DBL_EPSILON * fabs(x->val);
  y->err  = r * fabs(c * GSL_DBL_EPSILON * t);
  y->err += 2.0 * GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
 gsl_sf_rect_to_polar(const MpIeee x, const MpIeee y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  int  stat_h=  gsl_sf_hypot_e(x, y, r);
  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = 2.0 * GSL_DBL_EPSILON * fabs(theta->val);
    return stat_h;
  }
  else {
    DOMAIN_ERROR(theta);
  }
}


int  gsl_sf_angle_restrict_symm_err_e(const MpIeee theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const MpIeee P1=  4 * 7.8539812564849853515625e-01;
  const MpIeee P2=  4 * 3.7748947079307981766760e-08;
  const MpIeee P3=  4 * 2.6951514290790594840552e-15;
  const MpIeee TwoPi=  2*(P1 + P2 + P3);

  const MpIeee y=  2*floor(theta/TwoPi);
  MpIeee r=  ((theta - y*P1) - y*P2) - y*P3;

  if(r >  M_PI) r -= TwoPi;
  result->val = r;

  if(theta > 0.0625/GSL_DBL_EPSILON) {
    result->err = fabs(result->val);
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(theta > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = GSL_SQRT_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_angle_restrict_pos_err_e(const MpIeee theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const MpIeee P1=  4 * 7.85398125648498535156e-01;
  const MpIeee P2=  4 * 3.77489470793079817668e-08;
  const MpIeee P3=  4 * 2.69515142907905952645e-15;
  const MpIeee TwoPi=  2*(P1 + P2 + P3);

  const MpIeee y=  2*floor(theta/TwoPi);

  result->val = ((theta - y*P1) - y*P2) - y*P3;

  if(theta > 0.0625/GSL_DBL_EPSILON) {
    result->err = fabs(result->val);
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(theta > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = GSL_SQRT_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_angle_restrict_symm_e(MpIeee * theta)
{
  gsl_sf_result r;
  int  stat=  gsl_sf_angle_restrict_symm_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int  gsl_sf_angle_restrict_pos_e(MpIeee * theta)
{
  gsl_sf_result r;
  int  stat=  gsl_sf_angle_restrict_pos_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int  gsl_sf_sin_err_e(const MpIeee x, const MpIeee dx, gsl_sf_result * result)
{
  int  stat_s=  gsl_sf_sin_e(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}


int  gsl_sf_cos_err_e(const MpIeee x, const MpIeee dx, gsl_sf_result * result)
{
  int  stat_c=  gsl_sf_cos_e(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}


#if 0
int
 gsl_sf_sin_pi_x_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(-100.0 < x && x < 100.0) {
    result->val = sin(M_PI * x) / (M_PI * x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee N=  floor(x + 0.5);
    const MpIeee f=  x - N;

    if(N < INT_MAX && N > INT_MIN) {
      /* Make it an integer if we can. Saves another
       * call to floor().
       */
      const int intN    = (int)N;
      const MpIeee sign=  ( GSL_IS_ODD(intN) ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI * f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }
    else if(N > 2.0/GSL_DBL_EPSILON || N < -2.0/GSL_DBL_EPSILON) {
      /* All integer-valued floating point numbers
       * bigger than 2/eps=2^53 are actually even.
       */
      result->val = 0.0;
      result->err = 0.0;
    }
    else {
      const MpIeee resN=  N - 2.0*floor(0.5*N); /* 0 for even N, 1 for odd N */
      const MpIeee sign=  ( fabs(resN) > 0.5 ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI*f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_SUCCESS;
  }
}
#endif


int  gsl_sf_sinc_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const MpIeee ax=  fabs(x);

    if(ax < 0.8) {
      /* Do not go to the limit of the fit since
       * there is a zero there and the Chebyshev
       * accuracy will go to zero.
       */
      return cheb_eval_e(&sinc_cs, 2.0*ax-1.0, result);
    }
    else if(ax < 100.0) {
      /* Small arguments are no problem.
       * We trust the library sin() to
       * roughly machine precision.
       */
      result->val = sin(M_PI * ax)/(M_PI * ax);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Large arguments must be handled separately.
       */
      const MpIeee r=  M_PI*ax;
      gsl_sf_result s;
      int  stat_s=  gsl_sf_sin_e(r, &s);
      result->val = s.val/r;
      result->err = s.err/r + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_s;
    }
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_sin(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_sin_e(x, &result));
}

MpIeee gsl_sf_cos(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_cos_e(x, &result));
}

MpIeee gsl_sf_hypot(const MpIeee x, const MpIeee y)
{
  EVAL_RESULT(gsl_sf_hypot_e(x, y, &result));
}

MpIeee gsl_sf_lnsinh(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_lnsinh_e(x, &result));
}

MpIeee gsl_sf_lncosh(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_lncosh_e(x, &result));
}

MpIeee gsl_sf_angle_restrict_symm(const MpIeee theta)
{
  MpIeee result=  theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_symm_e(&result));
}

MpIeee gsl_sf_angle_restrict_pos(const MpIeee theta)
{
  MpIeee result=  theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_pos_e(&result));
}

#if 0
MpIeee gsl_sf_sin_pi_x(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_sin_pi_x_e(x, &result));
}
#endif

MpIeee gsl_sf_sinc(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_sinc_e(x, &result));
}
