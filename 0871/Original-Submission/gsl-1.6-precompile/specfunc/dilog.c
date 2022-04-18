#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/dilog.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_dilog.h>


/* Evaluate series for real dilog(x)
 * Sum[ x^k / k^2, {k,1,Infinity}]
 *
 * Converges rapidly for |x| < 1/2.
 */
static
int
 dilog_series_1(const MpIeee x, gsl_sf_result * result)
{
  const int kmax = 1000;
  MpIeee sum=  x;
  MpIeee term=  x;
  int  k;
  for(k=2; k<kmax; k++) {
    const MpIeee rk=  (k-1.0)/k;
    term *= x;
    term *= rk*rk;
    sum += term;
    if(fabs(term/sum) < GSL_DBL_EPSILON) break;
  }

  result->val  = sum;
  result->err  = 2.0 * fabs(term);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  if(k == kmax)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Compute the associated series
 *
 *   sum_{k=1}{infty} r^k / (k^2 (k+1))
 *
 * This is a series which appears in the one-step accelerated
 * method, which splits out one elementary function from the
 * full definition of Li_2(x). See below.
 */
static int
 series_2(MpIeee r, gsl_sf_result * result)
{
  static const int kmax = 100;
  MpIeee rk=  r;
  MpIeee sum=  MpIeee( "0.5" ) * r;
  int  k;
  for(k=2; k<10; k++)
  {
    MpIeee ds;
    rk *= r;
    ds = rk/(k*k*(k+MpIeee( "1.0" )));
    sum += ds;
  }
  for(; k<kmax; k++)
  {
    MpIeee ds;
    rk *= r;
    ds = rk/(k*k*(k+MpIeee( "1.0" )));
    sum += ds;
    if(fabs(ds/sum) < 0.5*GSL_DBL_EPSILON) break;
  }

  result->val = sum;
  result->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(sum);

  return GSL_SUCCESS;
}


/* Compute Li_2(x) using the accelerated series representation.
 *
 * Li_2(x) = 1 + (1-x)ln(1-x)/x + series_2(x)
 *
 * assumes: -1 < x < 1
 */
static int
 dilog_series_2(MpIeee x, gsl_sf_result * result)
{
  const int stat_s3 = series_2(x, result);
  MpIeee t;
  if(x > MpIeee( "0.01" ))
    t = (MpIeee( "1.0" ) - x) * log(MpIeee( "1.0" )-x) / x;
  else
  {
    static const MpIeee c3=  1.0/3.0;
    static const MpIeee c4=  1.0/4.0;
    static const MpIeee c5=  1.0/5.0;
    static const MpIeee c6=  1.0/6.0;
    static const MpIeee c7=  1.0/7.0;
    static const MpIeee c8=  1.0/8.0;
    const MpIeee t68=  c6 + x*(c7 + x*c8);
    const MpIeee t38=  c3 + x *(c4 + x *(c5 + x * t68));
    t = (x - MpIeee( "1.0" )) * (MpIeee( "1.0" ) + x*(MpIeee( "0.5" ) + x*t38));
  }
  result->val += 1.0 + t;
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(t);
  return stat_s3;
}


/* Calculates Li_2(x) for real x. Assumes x >= 0.0.
 */
static
int
 dilog_xge0(const MpIeee x, gsl_sf_result * result)
{
  if(x > 2.0) {
    gsl_sf_result ser;
    const int stat_ser = dilog_series_2(1.0/x, &ser);
    const MpIeee log_x=  log(x);
    const MpIeee t1=  M_PI*M_PI/3.0;
    const MpIeee t2=  ser.val;
    const MpIeee t3=  0.5*log_x*log_x;
    result->val  = t1 - t2 - t3;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
    result->val += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 1.01) {
    gsl_sf_result ser;
    const int stat_ser = dilog_series_2(1.0 - 1.0/x, &ser);
    const MpIeee log_x=  log(x);
    const MpIeee log_term=  log_x * (log(1.0-1.0/x) + 0.5*log_x);
    const MpIeee t1=  M_PI*M_PI/6.0;
    const MpIeee t2=  ser.val;
    const MpIeee t3=  log_term;
    result->val  = t1 + t2 - t3;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 1.0) {
    /* series around x = 1.0 */
    const MpIeee eps=  x - 1.0;
    const MpIeee lne=  log(eps);
    const MpIeee c0=  M_PI*M_PI/6.0;
    const MpIeee c1=    1.0 - lne;
    const MpIeee c2=  -(1.0 - 2.0*lne)/4.0;
    const MpIeee c3=   (1.0 - 3.0*lne)/9.0;
    const MpIeee c4=  -(1.0 - 4.0*lne)/16.0;
    const MpIeee c5=   (1.0 - 5.0*lne)/25.0;
    const MpIeee c6=  -(1.0 - 6.0*lne)/36.0;
    const MpIeee c7=   (1.0 - 7.0*lne)/49.0;
    const MpIeee c8=  -(1.0 - 8.0*lne)/64.0;
    result->val = c0+eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = M_PI*M_PI/6.0;
    result->err = 2.0 * GSL_DBL_EPSILON * M_PI*M_PI/6.0;
    return GSL_SUCCESS;
  }
  else if(x > 0.5) {
    gsl_sf_result ser;
    const int stat_ser = dilog_series_2(1.0-x, &ser);
    const MpIeee log_x=  log(x);
    const MpIeee t1=  M_PI*M_PI/6.0;
    const MpIeee t2=  ser.val;
    const MpIeee t3=  log_x*log(1.0-x);
    result->val  = t1 - t2 - t3;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 0.25) {
    return dilog_series_2(x, result);
  }
  else if(x > 0.0) {
    return dilog_series_1(x, result);
  }
  else {
    /* x == 0.0 */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
}


/* Evaluate the series representation for Li2(z):
 *
 *   Li2(z) = Sum[ |z|^k / k^2 Exp[i k arg(z)], {k,1,Infinity}]
 *   |z|    = r
 *   arg(z) = theta
 *   
 * Assumes 0 < r < 1.
 * It is used only for small r.
 */
static
int
 dilogc_series_1(
  const MpIeee r,
  const MpIeee x,
  const MpIeee y,
  gsl_sf_result * real_result,
  gsl_sf_result * imag_result
  )
{
  const MpIeee cos_theta=  x/r;
  const MpIeee sin_theta=  y/r;
  const MpIeee alpha=  1.0 - cos_theta;
  const MpIeee beta=  sin_theta;
  MpIeee ck=  cos_theta;
  MpIeee sk=  sin_theta;
  MpIeee rk=  r;
  MpIeee real_sum=  r*ck;
  MpIeee imag_sum=  r*sk;
  const int kmax = 50 + (int)(22.0/(-log(r))); /* tuned for double-precision */
  int  k;
  for(k=2; k<kmax; k++) {
    MpIeee dr;MpIeee  di;
    MpIeee ck_tmp=  ck;
    ck = ck - (alpha*ck + beta*sk);
    sk = sk - (alpha*sk - beta*ck_tmp);
    rk *= r;
    dr = rk/((MpIeee)k*k) * ck;
    di = rk/((MpIeee)k*k) * sk;
    real_sum += dr;
    imag_sum += di;
    if(fabs((dr*dr + di*di)/(real_sum*real_sum + imag_sum*imag_sum)) < GSL_DBL_EPSILON*GSL_DBL_EPSILON) break;
  }

  real_result->val = real_sum;
  real_result->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(real_sum);
  imag_result->val = imag_sum;
  imag_result->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(imag_sum);

  return GSL_SUCCESS;
}


/* Compute
 *
 *   sum_{k=1}{infty} z^k / (k^2 (k+1))
 *
 * This is a series which appears in the one-step accelerated
 * method, which splits out one elementary function from the
 * full definition of Li_2.
 */
static int
 series_2_c(
  MpIeee r,
  MpIeee x,
  MpIeee y,
  gsl_sf_result * sum_re,
  gsl_sf_result * sum_im
  )
{
  const MpIeee cos_theta=  x/r;
  const MpIeee sin_theta=  y/r;
  const MpIeee alpha=  1.0 - cos_theta;
  const MpIeee beta=  sin_theta;
  MpIeee ck=  cos_theta;
  MpIeee sk=  sin_theta;
  MpIeee rk=  r;
  MpIeee real_sum=  MpIeee( "0.5" ) * r*ck;
  MpIeee imag_sum=  MpIeee( "0.5" ) * r*sk;
  const int kmax = 30 + (int)(18.0/(-log(r))); /* tuned for double-precision */
  int  k;
  for(k=2; k<kmax; k++)
  {
    MpIeee dr;MpIeee  di;
    const MpIeee ck_tmp=  ck;
    ck = ck - (alpha*ck + beta*sk);
    sk = sk - (alpha*sk - beta*ck_tmp);
    rk *= r;
    dr = rk/((MpIeee)k*k*(k+MpIeee( "1.0" ))) * ck;
    di = rk/((MpIeee)k*k*(k+MpIeee( "1.0" ))) * sk;
    real_sum += dr;
    imag_sum += di;
    if(fabs((dr*dr + di*di)/(real_sum*real_sum + imag_sum*imag_sum)) < GSL_DBL_EPSILON*GSL_DBL_EPSILON) break;
  }

  sum_re->val = real_sum;
  sum_re->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(real_sum);
  sum_im->val = imag_sum;
  sum_im->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(imag_sum);

  return GSL_SUCCESS;
}


/* Compute Li_2(z) using the one-step accelerated series.
 *
 * Li_2(z) = 1 + (1-z)ln(1-z)/z + series_2_c(z)
 *
 * z = r exp(i theta)
 * assumes: r < 1
 * assumes: r > epsilon, so that we take no special care with log(1-z)
 */
static
int
 dilogc_series_2(
  const MpIeee r,
  const MpIeee x,
  const MpIeee y,
  gsl_sf_result * real_dl,
  gsl_sf_result * imag_dl
  )
{
  if(r == 0.0)
  {
    real_dl->val = 0.0;
    imag_dl->val = 0.0;
    real_dl->err = 0.0;
    imag_dl->err = 0.0;
    return GSL_SUCCESS;
  }
  else
  {
    gsl_sf_result sum_re;
    gsl_sf_result sum_im;
    const int stat_s3 = series_2_c(r, x, y, &sum_re, &sum_im);

    /* t = ln(1-z)/z */
    gsl_sf_result ln_omz_r;
    gsl_sf_result ln_omz_theta;
    const int stat_log = gsl_sf_complex_log_e(1.0-x, -y, &ln_omz_r, &ln_omz_theta);
    const MpIeee t_x=  ( ln_omz_r.val * x + ln_omz_theta.val * y)/(r*r);
    const MpIeee t_y=  (-ln_omz_r.val * y + ln_omz_theta.val * x)/(r*r);

    /* r = (1-z) ln(1-z)/z */
    const MpIeee r_x=  (1.0 - x) * t_x + y * t_y;
    const MpIeee r_y=  (1.0 - x) * t_y - y * t_x;

    real_dl->val = sum_re.val + r_x + 1.0;
    imag_dl->val = sum_im.val + r_y;
    real_dl->err = sum_re.err + 2.0*GSL_DBL_EPSILON*(fabs(real_dl->val) + fabs(r_x));
    imag_dl->err = sum_im.err + 2.0*GSL_DBL_EPSILON*(fabs(imag_dl->val) + fabs(r_y));
    return GSL_ERROR_SELECT_2(stat_s3, stat_log);
  }
}


/* Evaluate a series for Li_2(z) when |z| is near 1.
 * This is uniformly good away from z=1.
 *
 *   Li_2(z) = Sum[ a^n/n! H_n(theta), {n, 0, Infinity}]
 *
 * where
 *   H_n(theta) = Sum[ e^(i m theta) m^n / m^2, {m, 1, Infinity}]
 *   a = ln(r)
 *
 *  H_0(t) = Gl_2(t) + i Cl_2(t)
 *  H_1(t) = 1/2 ln(2(1-c)) + I atan2(-s, 1-c)
 *  H_2(t) = -1/2 + I/2 s/(1-c)
 *  H_3(t) = -1/2 /(1-c)
 *  H_4(t) = -I/2 s/(1-c)^2
 *  H_5(t) = 1/2 (2 + c)/(1-c)^2
 *  H_6(t) = I/2 s/(1-c)^5 (8(1-c) - s^2 (3 + c))
 */
static
int
 dilogc_series_3(
  const MpIeee r,
  const MpIeee x,
  const MpIeee y,
  gsl_sf_result * real_result,
  gsl_sf_result * imag_result
  )
{
  const MpIeee theta=  atan2(y, x);
  const MpIeee cos_theta=  x/r;
  const MpIeee sin_theta=  y/r;
  const MpIeee a=  log(r);
  const MpIeee omc=  1.0 - cos_theta;
  const MpIeee omc2=  omc*omc;
  MpIeee H_re[7];
  MpIeee H_im[7];
  MpIeee an;MpIeee  nfact;
  MpIeee sum_re;MpIeee  sum_im;
  gsl_sf_result Him0;
  int  n;

  H_re[0] = M_PI*M_PI/MpIeee( "6.0" ) + MpIeee( "0.25" )*(theta*theta - MpIeee( "2.0" )*M_PI*fabs(theta));
  gsl_sf_clausen_e(theta, &Him0);
  H_im[0] = Him0.val;

  H_re[1] = -MpIeee( "0.5" )*log(MpIeee( "2.0" )*omc);
  H_im[1] = -atan2(-sin_theta, omc);

  H_re[2] = -MpIeee( "0.5" );
  H_im[2] = MpIeee( "0.5" ) * sin_theta/omc;

  H_re[3] = -MpIeee( "0.5" )/omc;
  H_im[3] = MpIeee( "0.0" );

  H_re[4] = MpIeee( "0.0" );
  H_im[4] = -MpIeee( "0.5" )*sin_theta/omc2;

  H_re[5] = MpIeee( "0.5" ) * (MpIeee( "2.0" ) + cos_theta)/omc2;
  H_im[5] = MpIeee( "0.0" );

  H_re[6] = MpIeee( "0.0" );
  H_im[6] = MpIeee( "0.5" ) * sin_theta/(omc2*omc2*omc) * (MpIeee( "8.0" )*omc - sin_theta*sin_theta*(MpIeee( "3.0" ) + cos_theta));

  sum_re = H_re[0];
  sum_im = H_im[0];
  an = MpIeee( "1.0" );
  nfact = MpIeee( "1.0" );
  for(n=1; n<=6; n++) {
    MpIeee t;
    an *= a;
    nfact *= n;
    t = an/nfact;
    sum_re += t * H_re[n];
    sum_im += t * H_im[n];
  }

  real_result->val = sum_re;
  real_result->err = 2.0 * 6.0 * GSL_DBL_EPSILON * fabs(sum_re) + fabs(an/nfact);
  imag_result->val = sum_im;
  imag_result->err = 2.0 * 6.0 * GSL_DBL_EPSILON * fabs(sum_im) + Him0.err + fabs(an/nfact);

  return GSL_SUCCESS;
}


/* Calculate complex dilogarithm Li_2(z) in the fundamental region,
 * which we take to be the intersection of the unit disk with the
 * half-space x < MAGIC_SPLIT_VALUE. It turns out that 0.732 is a
 * nice choice for MAGIC_SPLIT_VALUE since then points mapped out
 * of the x > MAGIC_SPLIT_VALUE region and into another part of the
 * unit disk are bounded in radius by MAGIC_SPLIT_VALUE itself.
 *
 * If |z| < 0.98 we use a direct series summation. Otherwise z is very
 * near the unit circle, and the series_2 expansion is used; see above.
 * Because the fundamental region is bounded away from z = 1, this
 * works well.
 */
static
int
 dilogc_fundamental(MpIeee r, MpIeee x, MpIeee y, gsl_sf_result * real_dl, gsl_sf_result * imag_dl)
{
  if(r > MpIeee( "0.98" ))  
    return dilogc_series_3(r, x, y, real_dl, imag_dl);
  else if(r > MpIeee( "0.25" ))
    return dilogc_series_2(r, x, y, real_dl, imag_dl);
  else
    return dilogc_series_1(r, x, y, real_dl, imag_dl);
}


/* Compute Li_2(z) for z in the unit disk, |z| < 1. If z is outside
 * the fundamental region, which means that it is too close to z = 1,
 * then it is reflected into the fundamental region using the identity
 *
 *   Li2(z) = -Li2(1-z) + zeta(2) - ln(z) ln(1-z).
 */
static
int
 dilogc_unitdisk(MpIeee x, MpIeee y, gsl_sf_result * real_dl, gsl_sf_result * imag_dl)
{
  static const MpIeee MAGIC_SPLIT_VALUE=  0.732;
  static const MpIeee zeta2=  M_PI*M_PI/6.0;
  const MpIeee r=  sqrt(x*x + y*y);

  if(x > MAGIC_SPLIT_VALUE)
  {
    /* Reflect away from z = 1 if we are too close. The magic value
     * insures that the reflected value of the radius satisfies the
     * related inequality r_tmp < MAGIC_SPLIT_VALUE.
     */
    const MpIeee x_tmp=  1.0 - x;
    const MpIeee y_tmp=      - y;
    const MpIeee r_tmp=  sqrt(x_tmp*x_tmp + y_tmp*y_tmp);
    /* const double cos_theta_tmp = x_tmp/r_tmp; */
    /* const double sin_theta_tmp = y_tmp/r_tmp; */

    gsl_sf_result result_re_tmp;
    gsl_sf_result result_im_tmp;

    const int stat_dilog = dilogc_fundamental(r_tmp, x_tmp, y_tmp, &result_re_tmp, &result_im_tmp);

    const MpIeee lnz=   log(r);               /*  log(|z|)   */
    const MpIeee lnomz=   log(r_tmp);           /*  log(|1-z|) */
    const MpIeee argz=   atan2(y, x);          /*  arg(z) assuming principal branch */
    const MpIeee argomz=   atan2(y_tmp, x_tmp);  /*  arg(1-z)   */
    real_dl->val  = -result_re_tmp.val + zeta2 - lnz*lnomz + argz*argomz;
    real_dl->err  =  result_re_tmp.err;
    real_dl->err +=  2.0 * GSL_DBL_EPSILON * (zeta2 + fabs(lnz*lnomz) + fabs(argz*argomz));
    imag_dl->val  = -result_im_tmp.val - argz*lnomz - argomz*lnz;
    imag_dl->err  =  result_im_tmp.err;
    imag_dl->err +=  2.0 * GSL_DBL_EPSILON * (fabs(argz*lnomz) + fabs(argomz*lnz));

    return stat_dilog;
  }
  else
  {
    return dilogc_fundamental(r, x, y, real_dl, imag_dl);
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int
 gsl_sf_dilog_e(const MpIeee x, gsl_sf_result * result)
{
  if(x >= 0.0) {
    return dilog_xge0(x, result);
  }
  else {
    gsl_sf_result d1, d2;
    int  stat_d1=  dilog_xge0( -x, &d1);
    int  stat_d2=  dilog_xge0(x*x, &d2);
    result->val  = -d1.val + 0.5 * d2.val;
    result->err  =  d1.err + 0.5 * d2.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_d1, stat_d2);
  }
}


int
 gsl_sf_complex_dilog_xy_e(
  const MpIeee x,
  const MpIeee y,
  gsl_sf_result * real_dl,
  gsl_sf_result * imag_dl
  )
{
  const MpIeee zeta2=  M_PI*M_PI/6.0;
  const MpIeee r2=  x*x + y*y;

  if(y == 0.0)
  {
    if(x >= 1.0)
    {
      imag_dl->val = -M_PI * log(x);
      imag_dl->err = 2.0 * GSL_DBL_EPSILON * fabs(imag_dl->val);
    }
    else
    {
      imag_dl->val = 0.0;
      imag_dl->err = 0.0;
    }
    return gsl_sf_dilog_e(x, real_dl);
  }
  else if(fabs(r2 - 1.0) < GSL_DBL_EPSILON)
  {
    /* Lewin A.2.4.1 and A.2.4.2 */

    const MpIeee theta=  atan2(y, x);
    const MpIeee term1=  theta*theta/4.0;
    const MpIeee term2=  M_PI*fabs(theta)/2.0;
    real_dl->val = zeta2 + term1 - term2;
    real_dl->err = 2.0 * GSL_DBL_EPSILON * (zeta2 + term1 + term2);
    return gsl_sf_clausen_e(theta, imag_dl);
  }
  else if(r2 < 1.0)
  {
    return dilogc_unitdisk(x, y, real_dl, imag_dl);
  }
  else
  {
    /* Reduce argument to unit disk. */
    const MpIeee r=  sqrt(r2);
    const MpIeee x_tmp=   x/r2;
    const MpIeee y_tmp=  -y/r2;
    /* const double r_tmp = 1.0/r; */
    gsl_sf_result result_re_tmp, result_im_tmp;

    const int stat_dilog =
      dilogc_unitdisk(x_tmp, y_tmp, &result_re_tmp, &result_im_tmp);

    /* Unwind the inversion.
     *
     *  Li_2(z) + Li_2(1/z) = -zeta(2) - 1/2 ln(-z)^2
     */
    const MpIeee theta=  atan2(y, x);
    const MpIeee theta_abs=  fabs(theta);
    const MpIeee theta_sgn=  ( theta < 0.0 ? -1.0 : 1.0 );
    const MpIeee ln_minusz_re=  log(r);
    const MpIeee ln_minusz_im=  theta_sgn * (theta_abs - M_PI);
    const MpIeee lmz2_re=  ln_minusz_re*ln_minusz_re - ln_minusz_im*ln_minusz_im;
    const MpIeee lmz2_im=  2.0*ln_minusz_re*ln_minusz_im;
    real_dl->val = -result_re_tmp.val - 0.5 * lmz2_re - zeta2;
    real_dl->err =  result_re_tmp.err + 2.0*GSL_DBL_EPSILON*(0.5 * fabs(lmz2_re) + zeta2);
    imag_dl->val = -result_im_tmp.val - 0.5 * lmz2_im;
    imag_dl->err =  result_im_tmp.err + 2.0*GSL_DBL_EPSILON*fabs(lmz2_im);
    return stat_dilog;
  }
}


int
 gsl_sf_complex_dilog_e(
  const MpIeee r,
  const MpIeee theta,
  gsl_sf_result * real_dl,
  gsl_sf_result * imag_dl
  )
{
  const MpIeee cos_theta=  cos(theta);
  const MpIeee sin_theta=  sin(theta);
  const MpIeee x=  r * cos_theta;
  const MpIeee y=  r * sin_theta;
  return gsl_sf_complex_dilog_xy_e(x, y, real_dl, imag_dl);
}


int
 gsl_sf_complex_spence_xy_e(
  const MpIeee x,
  const MpIeee y,
  gsl_sf_result * real_sp,
  gsl_sf_result * imag_sp
  )
{
  const MpIeee oms_x=  1.0 - x;
  const MpIeee oms_y=      - y;
  return gsl_sf_complex_dilog_xy_e(oms_x, oms_y, real_sp, imag_sp);
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_dilog(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_dilog_e(x, &result));
}
