#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/hyperg_1F1.c
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
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_hyperg.h>

#include "error.h"
#include "hyperg.h"

#define _1F1_INT_THRESHOLD (100.0*GSL_DBL_EPSILON)


/* Asymptotic result for 1F1(a, b, x)  x -> -Infinity.
 * Assumes b-a != neg integer and b != neg integer.
 */
static
int
 hyperg_1F1_asymp_negx(const MpIeee a, const MpIeee b, const MpIeee x,
                     gsl_sf_result * result)
{
  gsl_sf_result lg_b;
  gsl_sf_result lg_bma;
  MpIeee sgn_b;
  MpIeee sgn_bma;

  int  stat_b=  gsl_sf_lngamma_sgn_e(b,   &lg_b,   &sgn_b);
  int  stat_bma=  gsl_sf_lngamma_sgn_e(b-a, &lg_bma, &sgn_bma);

  if(stat_b == GSL_SUCCESS && stat_bma == GSL_SUCCESS) {
    gsl_sf_result F;
    int  stat_F=  gsl_sf_hyperg_2F0_series_e(a, 1.0+a-b, -1.0/x, -1, &F);
    if(F.val != 0) {
      MpIeee ln_term_val=  a*log(-x);
      MpIeee ln_term_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(a) + fabs(ln_term_val));
      MpIeee ln_pre_val=  lg_b.val - lg_bma.val - ln_term_val;
      MpIeee ln_pre_err=  lg_b.err + lg_bma.err + ln_term_err;
      int  stat_e=  gsl_sf_exp_mult_err_e(ln_pre_val, ln_pre_err,
                                            sgn_bma*sgn_b*F.val, F.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_F);
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
      return stat_F;
    }
  }
  else {
    DOMAIN_ERROR(result);
  }
}


/* Asymptotic result for 1F1(a, b, x)  x -> +Infinity
 * Assumes b != neg integer and a != neg integer
 */
static
int
 hyperg_1F1_asymp_posx(const MpIeee a, const MpIeee b, const MpIeee x,
                      gsl_sf_result * result)
{
  gsl_sf_result lg_b;
  gsl_sf_result lg_a;
  MpIeee sgn_b;
  MpIeee sgn_a;

  int  stat_b=  gsl_sf_lngamma_sgn_e(b, &lg_b, &sgn_b);
  int  stat_a=  gsl_sf_lngamma_sgn_e(a, &lg_a, &sgn_a);

  if(stat_a == GSL_SUCCESS && stat_b == GSL_SUCCESS) {
    gsl_sf_result F;
    int  stat_F=  gsl_sf_hyperg_2F0_series_e(b-a, 1.0-a, 1.0/x, -1, &F);
    if(stat_F == GSL_SUCCESS && F.val != 0) {
      MpIeee lnx=  log(x);
      MpIeee ln_term_val=  (a-b)*lnx;
      MpIeee ln_term_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(a) + fabs(b)) * fabs(lnx)
                         + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(a-b);
      MpIeee ln_pre_val=  lg_b.val - lg_a.val + ln_term_val + x;
      MpIeee ln_pre_err=  lg_b.err + lg_a.err + ln_term_err + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(x);
      int  stat_e=  gsl_sf_exp_mult_err_e(ln_pre_val, ln_pre_err,
                                            sgn_a*sgn_b*F.val, F.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_F);
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
      return stat_F;
    }
  }
  else {
    DOMAIN_ERROR(result);
  }
}


/* Asymptotic result for x < 2b-4a, 2b-4a large.
 * [Abramowitz+Stegun, 13.5.21]
 *
 * assumes 0 <= x/(2b-4a) <= 1
 */
static
int
 hyperg_1F1_large2bm4a(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  MpIeee eta=  MpIeee( "2.0" )*b - MpIeee( "4.0" )*a;
  MpIeee cos2th=  x/eta;
  MpIeee sin2th=  MpIeee( "1.0" ) - cos2th;
  MpIeee th=  acos(sqrt(cos2th));
  MpIeee pre_h=  MpIeee( "0.25" )*M_PI*M_PI*eta*eta*cos2th*sin2th;
  gsl_sf_result lg_b;
  int  stat_lg=  gsl_sf_lngamma_e(b, &lg_b);
  MpIeee t1=  MpIeee( "0.5" )*(MpIeee( "1.0" )-b)*log(MpIeee( "0.25" )*x*eta);
  MpIeee t2=  MpIeee( "0.25" )*log(pre_h);
  MpIeee lnpre_val=  lg_b.val + MpIeee( "0.5" )*x + t1 - t2;
  MpIeee lnpre_err=  lg_b.err + MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(MpIeee( "0.5" )*x) + fabs(t1) + fabs(t2));
  MpIeee s1=  sin(a*M_PI);
  MpIeee s2=  sin(MpIeee( "0.25" )*eta*(MpIeee( "2.0" )*th - sin(MpIeee( "2.0" )*th)) + MpIeee( "0.25" )*M_PI);
  MpIeee ser_val=  s1 + s2;
  MpIeee ser_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(s1) + fabs(s2));
  int  stat_e=  gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err,
                                        ser_val, ser_err,
                                        result);
  return GSL_ERROR_SELECT_2(stat_e, stat_lg);
}


/* Luke's rational approximation.
 * See [Luke, Algorithms for the Computation of Mathematical Functions, p.182]
 *
 * Like the case of the 2F1 rational approximations, these are
 * probably guaranteed to converge for x < 0, barring gross
 * numerical instability in the pre-asymptotic regime.
 */
static
int
 hyperg_1F1_luke(const MpIeee a, const MpIeee c, const MpIeee xin,
                gsl_sf_result * result)
{
  const MpIeee RECUR_BIG=  1.0e+50;
  const int nmax = 5000;
  int  n=  3;
  const MpIeee x=  -xin;
  const MpIeee x3=  x*x*x;
  const MpIeee t0=  a/c;
  const MpIeee t1=  (a+1.0)/(2.0*c);
  const MpIeee t2=  (a+2.0)/(2.0*(c+1.0));
  MpIeee F=  MpIeee( "1.0" );
  MpIeee prec;

  MpIeee Bnm3=  MpIeee( "1.0" );                                  /* B0 */
  MpIeee Bnm2=  MpIeee( "1.0" ) + t1 * x;                         /* B1 */
  MpIeee Bnm1=  MpIeee( "1.0" ) + t2 * x * (MpIeee( "1.0" ) + t1/MpIeee( "3.0" ) * x);    /* B2 */
 
  MpIeee Anm3=  MpIeee( "1.0" );                                                      /* A0 */
  MpIeee Anm2=  Bnm2 - t0 * x;                                            /* A1 */
  MpIeee Anm1=  Bnm1 - t0*(MpIeee( "1.0" ) + t2*x)*x + t0 * t1 * (c/(c+MpIeee( "1.0" ))) * x*x;   /* A2 */

  while(1) {
    MpIeee npam1=  n + a - MpIeee( "1" );
    MpIeee npcm1=  n + c - MpIeee( "1" );
    MpIeee npam2=  n + a - MpIeee( "2" );
    MpIeee npcm2=  n + c - MpIeee( "2" );
    MpIeee tnm1=  MpIeee( "2" )*n - MpIeee( "1" );
    MpIeee tnm3=  MpIeee( "2" )*n - MpIeee( "3" );
    MpIeee tnm5=  MpIeee( "2" )*n - MpIeee( "5" );
    MpIeee F1=   (n-a-MpIeee( "2" )) / (MpIeee( "2" )*tnm3*npcm1);
    MpIeee F2=   (n+a)*npam1 / (MpIeee( "4" )*tnm1*tnm3*npcm2*npcm1);
    MpIeee F3=  -npam2*npam1*(n-a-MpIeee( "2" )) / (MpIeee( "8" )*tnm3*tnm3*tnm5*(n+c-MpIeee( "3" ))*npcm2*npcm1);
    MpIeee E=  -npam1*(n-c-MpIeee( "1" )) / (MpIeee( "2" )*tnm3*npcm2*npcm1);

    MpIeee An=  (MpIeee( "1.0" )+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    MpIeee Bn=  (MpIeee( "1.0" )+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    MpIeee r=  An/Bn;

    prec = fabs((F - r)/F);
    F = r;

    if(prec < GSL_DBL_EPSILON || n > nmax) break;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An   /= RECUR_BIG;
      Bn   /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
      Anm3 /= RECUR_BIG;
      Bnm3 /= RECUR_BIG;
    }
    else if(fabs(An) < 1.0/RECUR_BIG || fabs(Bn) < 1.0/RECUR_BIG) {
      An   *= RECUR_BIG;
      Bn   *= RECUR_BIG;
      Anm1 *= RECUR_BIG;
      Bnm1 *= RECUR_BIG;
      Anm2 *= RECUR_BIG;
      Bnm2 *= RECUR_BIG;
      Anm3 *= RECUR_BIG;
      Bnm3 *= RECUR_BIG;
    }

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }

  result->val  = F;
  result->err  = 2.0 * fabs(F * prec);
  result->err += 2.0 * GSL_DBL_EPSILON * (n-1.0) * fabs(F);

  return GSL_SUCCESS;
}


/* Series for 1F1(1,b,x)
 * b > 0
 */
static
int
 hyperg_1F1_1_series(const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  MpIeee sum_val=  MpIeee( "1.0" );
  MpIeee sum_err=  MpIeee( "0.0" );
  MpIeee term=  MpIeee( "1.0" );
  MpIeee n=  MpIeee( "1.0" );
  while(fabs(term/sum_val) > 2.0*GSL_DBL_EPSILON) {
    term *= x/(b+n-MpIeee( "1" ));
    sum_val += term;
    sum_err += MpIeee( "2.0" ) * MpIeee( "4.0" ) * GSL_DBL_EPSILON * fabs(term);
    n += MpIeee( "1.0" );
  }
  result->val  = sum_val;
  result->err  = sum_err;
  result->err += 2.0 * fabs(term);
  return GSL_SUCCESS;
}


/* 1F1(1,b,x)
 * b >= 1, b integer
 */
static
int
 hyperg_1F1_1_int(const int b, const MpIeee x, gsl_sf_result * result)
{
  if(b < 1) {
    DOMAIN_ERROR(result);
  }
  else if(b == 1) {
    return gsl_sf_exp_e(x, result);
  }
  else if(b == 2) {
    return gsl_sf_exprel_e(x, result);
  }
  else if(b == 3) {
    return gsl_sf_exprel_2_e(x, result);
  }
  else {
    return gsl_sf_exprel_n_e(b-1, x, result);
  }
}


/* 1F1(1,b,x)
 * b >=1, b real
 *
 * checked OK: [GJ] Thu Oct  1 16:46:35 MDT 1998
 */
static
int
 hyperg_1F1_1(const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  MpIeee ax=  fabs(x);
  MpIeee ib=  floor(b + MpIeee( "0.1" ));

  if(b < 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(b == 1.0) {
    return gsl_sf_exp_e(x, result);
  }
  else if(b >= 1.4*ax) {
    return hyperg_1F1_1_series(b, x, result);
  }
  else if(fabs(b - ib) < _1F1_INT_THRESHOLD && ib < INT_MAX) {
    return hyperg_1F1_1_int((int)ib, x, result);
  }
  else if(x > 0.0) {
    if(x > 100.0 && b < 0.75*x) {
      return hyperg_1F1_asymp_posx(1.0, b, x, result);
    }
    else if(b < 1.0e+05) {
      /* Recurse backward on b, from a
       * chosen offset point. For x > 0,
       * which holds here, this should
       * be a stable direction.
       */
      const MpIeee off=  ceil(1.4*x-b) + 1.0;
      MpIeee bp=  b + off;
      gsl_sf_result M;
      int  stat_s=  hyperg_1F1_1_series(bp, x, &M);
      const MpIeee err_rat=  M.err / fabs(M.val);
      while(bp > b+MpIeee( "0.1" )) {
        /* M(1,b-1) = x/(b-1) M(1,b) + 1 */
        bp -= MpIeee( "1.0" );
        M.val  = 1.0 + x/bp * M.val;
      }
      result->val  = M.val;
      result->err  = err_rat * fabs(M.val);
      result->err += 2.0 * GSL_DBL_EPSILON * (fabs(off)+1.0) * fabs(M.val);
      return stat_s;
    }
    else {
      return hyperg_1F1_large2bm4a(1.0, b, x, result);
    }
  }
  else {
    /* x <= 0 and b not large compared to |x|
     */
    if(ax < MpIeee( "10.0" ) && b < MpIeee( "10.0" )) {
      return hyperg_1F1_1_series(b, x, result);
    }
    else if(ax >= MpIeee( "100.0" ) && GSL_MAX_DBL(fabs(MpIeee( "2.0" )-b),MpIeee( "1.0" )) < MpIeee( "0.99" )*ax) {
      return hyperg_1F1_asymp_negx(1.0, b, x, result);
    }
    else {
      return hyperg_1F1_luke(1.0, b, x, result);
    }
  }
}


/* 1F1(a,b,x)/Gamma(b) for b->0
 * [limit of Abramowitz+Stegun 13.3.7]
 */
static
int
 hyperg_1F1_renorm_b0(const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  MpIeee eta=  a*x;
  if(eta > MpIeee( "0.0" )) {
    MpIeee root_eta=  sqrt(eta);
    gsl_sf_result I1_scaled;
    int  stat_I=  gsl_sf_bessel_I1_scaled_e(2.0*root_eta, &I1_scaled);
    if(I1_scaled.val <= 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_ERROR_SELECT_2(stat_I, GSL_EDOM);
    }
    else {
      const MpIeee lnr_val=  0.5*x + 0.5*log(eta) + fabs(x) + log(I1_scaled.val);
      const MpIeee lnr_err=  GSL_DBL_EPSILON * (1.5*fabs(x) + 1.0) + fabs(I1_scaled.err/I1_scaled.val);
      return gsl_sf_exp_err_e(lnr_val, lnr_err, result);
    }
  }
  else if(eta == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* eta < 0 */
    MpIeee root_eta=  sqrt(-eta);
    gsl_sf_result J1;
    int  stat_J=  gsl_sf_bessel_J1_e(2.0*root_eta, &J1);
    if(J1.val <= 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_ERROR_SELECT_2(stat_J, GSL_EDOM);
    }
    else {
      const MpIeee t1=  0.5*x;
      const MpIeee t2=  0.5*log(-eta);
      const MpIeee t3=  fabs(x);
      const MpIeee t4=  log(J1.val);
      const MpIeee lnr_val=  t1 + t2 + t3 + t4;
      const MpIeee lnr_err=  GSL_DBL_EPSILON * (1.5*fabs(x) + 1.0) + fabs(J1.err/J1.val);
      gsl_sf_result ex;
      int  stat_e=  gsl_sf_exp_err_e(lnr_val, lnr_err, &ex);
      result->val = -ex.val;
      result->err =  ex.err;
      return stat_e;
    }
  }
  
}


/* 1F1'(a,b,x)/1F1(a,b,x)
 * Uses Gautschi's version of the CF.
 * [Gautschi, Math. Comp. 31, 994 (1977)]
 *
 * Supposedly this suffers from the "anomalous convergence"
 * problem when b < x. I have seen anomalous convergence
 * in several of the continued fractions associated with
 * 1F1(a,b,x). This particular CF formulation seems stable
 * for b > x. However, it does display a painful artifact
 * of the anomalous convergence; the convergence plateaus
 * unless b >>> x. For example, even for b=1000, x=1, this
 * method locks onto a ratio which is only good to about
 * 4 digits. Apparently the rest of the digits are hiding
 * way out on the plateau, but finite-precision lossage
 * means you will never get them.
 */
#if 0
static
int
 hyperg_1F1_CF1_p(const MpIeee a, const MpIeee b, const MpIeee x, MpIeee * result)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  MpIeee( "1.0" );
  MpIeee b1=  MpIeee( "1.0" );
  MpIeee An=  b1*Anm1 + a1*Anm2;
  MpIeee Bn=  b1*Bnm1 + a1*Bnm2;
  MpIeee an;MpIeee  bn;
  MpIeee fn=  An/Bn;

  while(n < maxiter) {
    MpIeee old_fn;
    MpIeee del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = (a+n)*x/((b-x+n-MpIeee( "1" ))*(b-x+n));
    bn = MpIeee( "1.0" );
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;
    
    if(fabs(del - 1.0) < 10.0*GSL_DBL_EPSILON) break;
  }

  *result = a/(b-x) * fn;

  if(n == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}
#endif /* 0 */


/* 1F1'(a,b,x)/1F1(a,b,x)
 * Uses Gautschi's series transformation of the
 * continued fraction. This is apparently the best
 * method for getting this ratio in the stable region.
 * The convergence is monotone and supergeometric
 * when b > x.
 * Assumes a >= -1.
 */
static
int
 hyperg_1F1_CF1_p_ser(const MpIeee a, const MpIeee b, const MpIeee x, MpIeee * result)
{
  if(a == 0.0) {
    *result = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else {
    const int maxiter = 5000;
    MpIeee sum=  MpIeee( "1.0" );
    MpIeee pk=  MpIeee( "1.0" );
    MpIeee rhok=  MpIeee( "0.0" );
    int  k;
    for(k=1; k<maxiter; k++) {
      MpIeee ak=  (a + k)*x/((b-x+k-MpIeee( "1.0" ))*(b-x+k));
      rhok = -ak*(MpIeee( "1.0" ) + rhok)/(MpIeee( "1.0" ) + ak*(MpIeee( "1.0" )+rhok));
      pk  *= rhok;
      sum += pk;
      if(fabs(pk/sum) < 2.0*GSL_DBL_EPSILON) break;
    }
    *result = a/(b-x) * sum;
    if(k == maxiter)
      GSL_ERROR ("error", GSL_EMAXITER);
    else
      return GSL_SUCCESS;
  }
}


/* 1F1(a+1,b,x)/1F1(a,b,x)
 *
 * I think this suffers from typical "anomalous convergence".
 * I could not find a region where it was truly useful.
 */
#if 0
static
int
 hyperg_1F1_CF1(const MpIeee a, const MpIeee b, const MpIeee x, MpIeee * result)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  b - a - MpIeee( "1.0" );
  MpIeee b1=  b - x - MpIeee( "2.0" )*(a+MpIeee( "1.0" ));
  MpIeee An=  b1*Anm1 + a1*Anm2;
  MpIeee Bn=  b1*Bnm1 + a1*Bnm2;
  MpIeee an;MpIeee  bn;
  MpIeee fn=  An/Bn;

  while(n < maxiter) {
    MpIeee old_fn;
    MpIeee del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = (a + n - MpIeee( "1.0" )) * (b - a - n);
    bn = b - x - MpIeee( "2.0" )*(a+n);
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;
    
    if(fabs(del - 1.0) < 10.0*GSL_DBL_EPSILON) break;
  }

  *result = fn;
  if(n == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}
#endif /* 0 */


/* 1F1(a,b+1,x)/1F1(a,b,x)
 *
 * This seemed to suffer from "anomalous convergence".
 * However, I have no theory for this recurrence.
 */
#if 0
static
int
 hyperg_1F1_CF1_b(const MpIeee a, const MpIeee b, const MpIeee x, MpIeee * result)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  b + MpIeee( "1.0" );
  MpIeee b1=  (b + MpIeee( "1.0" )) * (b - x);
  MpIeee An=  b1*Anm1 + a1*Anm2;
  MpIeee Bn=  b1*Bnm1 + a1*Bnm2;
  MpIeee an;MpIeee  bn;
  MpIeee fn=  An/Bn;

  while(n < maxiter) {
    MpIeee old_fn;
    MpIeee del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = (b + n) * (b + n - MpIeee( "1.0" ) - a) * x;
    bn = (b + n) * (b + n - MpIeee( "1.0" ) - x);
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;
    
    if(fabs(del - 1.0) < 10.0*GSL_DBL_EPSILON) break;
  }

  *result = fn;
  if(n == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}
#endif /* 0 */


/* 1F1(a,b,x)
 * |a| <= 1, b > 0
 */
static
int
 hyperg_1F1_small_a_bgt0(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  const MpIeee bma=  b-a;
  const MpIeee oma=  1.0-a;
  const MpIeee ap1mb=  1.0+a-b;
  const MpIeee abs_bma=  fabs(bma);
  const MpIeee abs_oma=  fabs(oma);
  const MpIeee abs_ap1mb=  fabs(ap1mb);

  const MpIeee ax=  fabs(x);

  if(a == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == 1.0 && b >= 1.0) {
    return hyperg_1F1_1(b, x, result);
  }
  else if(a == -1.0) {
    result->val  = 1.0 + a/b * x;
    result->err  = GSL_DBL_EPSILON * (1.0 + fabs(a/b * x));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(b >= 1.4*ax) {
    return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
  }
  else if(x > 0.0) {
    if(x > 100.0 && abs_bma*abs_oma < 0.5*x) {
      return hyperg_1F1_asymp_posx(a, b, x, result);
    }
    else if(b < 5.0e+06) {
      /* Recurse backward on b from
       * a suitably high point.
       */
      const MpIeee b_del=  ceil(1.4*x-b) + 1.0;
      MpIeee bp=  b + b_del;
      gsl_sf_result r_Mbp1;
      gsl_sf_result r_Mb;
      MpIeee Mbp1;
      MpIeee Mb;
      MpIeee Mbm1;
      int  stat_0=  gsl_sf_hyperg_1F1_series_e(a, bp+1.0, x, &r_Mbp1);
      int  stat_1=  gsl_sf_hyperg_1F1_series_e(a, bp,     x, &r_Mb);
      const MpIeee err_rat=  fabs(r_Mbp1.err/r_Mbp1.val) + fabs(r_Mb.err/r_Mb.val);
      Mbp1 = r_Mbp1.val;
      Mb   = r_Mb.val;
      while(bp > b+MpIeee( "0.1" )) {
        /* Do backward recursion. */
        Mbm1 = ((x+bp-MpIeee( "1.0" ))*Mb - x*(bp-a)/bp*Mbp1)/(bp-MpIeee( "1.0" ));
        bp -= MpIeee( "1.0" );
        Mbp1 = Mb;
        Mb   = Mbm1;
      }
      result->val  = Mb;
      result->err  = err_rat * (fabs(b_del)+1.0) * fabs(Mb);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(Mb);
      return GSL_ERROR_SELECT_2(stat_0, stat_1);
    }
    else {
      return hyperg_1F1_large2bm4a(a, b, x, result);
    }
  }
  else {
    /* x < 0 and b not large compared to |x|
     */
    if(ax < 10.0 && b < 10.0) {
      return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
    }
    else if(ax >= 100.0 && GSL_MAX(abs_ap1mb,1.0) < 0.99*ax) {
      return hyperg_1F1_asymp_negx(a, b, x, result);
    }
    else {
      return hyperg_1F1_luke(a, b, x, result);
    }
  }
}


/* 1F1(b+eps,b,x)
 * |eps|<=1, b > 0
 */
static
int
 hyperg_1F1_beps_bgt0(const MpIeee eps, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  if(b > fabs(x) && fabs(eps) < GSL_SQRT_DBL_EPSILON) {
    /* If b-a is very small and x/b is not too large we can
     * use this explicit approximation.
     *
     * 1F1(b+eps,b,x) = exp(ax/b) (1 - eps x^2 (v2 + v3 x + ...) + ...)
     *
     *   v2 = a/(2b^2(b+1))
     *   v3 = a(b-2a)/(3b^3(b+1)(b+2))
     *   ...
     *
     * See [Luke, Mathematical Functions and Their Approximations, p.292]
     *
     * This cannot be used for b near a negative integer or zero.
     * Also, if x/b is large the deviation from exp(x) behaviour grows.
     */
    MpIeee a=  b + eps;
    gsl_sf_result exab;
    int  stat_e=  gsl_sf_exp_e(a*x/b, &exab);
    MpIeee v2=  a/(MpIeee( "2.0" )*b*b*(b+MpIeee( "1.0" )));
    MpIeee v3=  a*(b-MpIeee( "2.0" )*a)/(MpIeee( "3.0" )*b*b*b*(b+MpIeee( "1.0" ))*(b+MpIeee( "2.0" )));
    MpIeee v=  v2 + v3 * x;
    MpIeee f=  (MpIeee( "1.0" ) - eps*x*x*v);
    result->val  = exab.val * f;
    result->err  = exab.err * fabs(f);
    result->err += fabs(exab.val) * GSL_DBL_EPSILON * (1.0 + fabs(eps*x*x*v));
    result->err += 4.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_e;
  }
  else {
    /* Otherwise use a Kummer transformation to reduce
     * it to the small a case.
     */
    gsl_sf_result Kummer_1F1;
    int  stat_K=  hyperg_1F1_small_a_bgt0(-eps, b, -x, &Kummer_1F1);
    if(Kummer_1F1.val != 0.0) {
      int  stat_e=  gsl_sf_exp_mult_err_e(x, 2.0*GSL_DBL_EPSILON*fabs(x),
                                            Kummer_1F1.val, Kummer_1F1.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_K);
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
      return stat_K;
    }
  }
}


/* 1F1(a,2a,x) = Gamma(a + 1/2) E(x) (|x|/4)^(-a+1/2) scaled_I(a-1/2,|x|/2)
 *
 * E(x) = exp(x) x > 0
 *      = 1      x < 0
 *
 * a >= 1/2
 */
static
int
 hyperg_1F1_beq2a_pos(const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result I;
    int  stat_I=  gsl_sf_bessel_Inu_scaled_e(a-0.5, 0.5*fabs(x), &I);
    gsl_sf_result lg;
    int  stat_g=  gsl_sf_lngamma_e(a + 0.5, &lg);
    MpIeee ln_term=  (MpIeee( "0.5" )-a)*log(MpIeee( "0.25" )*fabs(x));
    MpIeee lnpre_val=  lg.val + GSL_MAX_DBL(x,MpIeee( "0.0" )) + ln_term;
    MpIeee lnpre_err=  lg.err + GSL_DBL_EPSILON * (fabs(ln_term) + fabs(x));
    int  stat_e=  gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err,
                                          I.val, I.err,
                                          result);
    return GSL_ERROR_SELECT_3(stat_e, stat_g, stat_I);
  }
}


/* Determine middle parts of diagonal recursion along b=2a
 * from two endpoints, i.e.
 *
 * given:  M(a,b)      and  M(a+1,b+2)
 * get:    M(a+1,b+1)  and  M(a,b+1)
 */
#if 0
inline
static
int
 hyperg_1F1_diag_step(const MpIeee a, const MpIeee b, const MpIeee x,
                     const MpIeee Mab, const MpIeee Map1bp2,
                     MpIeee * Map1bp1, MpIeee * Mabp1)
{
  if(a == b) {
    *Map1bp1 = Mab;
    *Mabp1   = Mab - x/(b+MpIeee( "1.0" )) * Map1bp2;
  }
  else {
    *Map1bp1 = Mab - x * (a-b)/(b*(b+MpIeee( "1.0" ))) * Map1bp2;
    *Mabp1   = (a * *Map1bp1 - b * Mab)/(a-b);
  }
  return GSL_SUCCESS;
}
#endif /* 0 */


/* Determine endpoint of diagonal recursion.
 *
 * given:  M(a,b)    and  M(a+1,b+2)
 * get:    M(a+1,b)  and  M(a+1,b+1)
 */
#if 0
inline
static
int
 hyperg_1F1_diag_end_step(const MpIeee a, const MpIeee b, const MpIeee x,
                         const MpIeee Mab, const MpIeee Map1bp2,
                         MpIeee * Map1b, MpIeee * Map1bp1)
{
  *Map1bp1 = Mab - x * (a-b)/(b*(b+MpIeee( "1.0" ))) * Map1bp2;
  *Map1b   = Mab + x/b * *Map1bp1;
  return GSL_SUCCESS;
}
#endif /* 0 */


/* Handle the case of a and b both positive integers.
 * Assumes a > 0 and b > 0.
 */
static
int
 hyperg_1F1_ab_posint(const int a, const int b, const MpIeee x, gsl_sf_result * result)
{
  MpIeee ax=  fabs(x);

  if(a == b) {
    return gsl_sf_exp_e(x, result);             /* 1F1(a,a,x) */
  }
  else if(a == 1) {
    return gsl_sf_exprel_n_e(b-1, x, result);   /* 1F1(1,b,x) */
  }
  else if(b == a + 1) {
    gsl_sf_result K;
    int  stat_K=  gsl_sf_exprel_n_e(a, -x, &K);  /* 1F1(1,1+a,-x) */
    int  stat_e=  gsl_sf_exp_mult_err_e(x, 2.0 * GSL_DBL_EPSILON * fabs(x),
                                          K.val, K.err,
                                          result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K);
  }
  else if(a == b + 1) {
    gsl_sf_result ex;
    int  stat_e=  gsl_sf_exp_e(x, &ex);
    result->val  = ex.val * (1.0 + x/b);
    result->err  = ex.err * (1.0 + x/b);
    result->err += ex.val * GSL_DBL_EPSILON * (1.0 + fabs(x/b));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_e;
  }
  else if(a == b + 2) {
    gsl_sf_result ex;
    int  stat_e=  gsl_sf_exp_e(x, &ex);
    MpIeee poly=  (MpIeee( "1.0" ) + x/b*(MpIeee( "2.0" ) + x/(b+MpIeee( "1.0" ))));
    result->val  = ex.val * poly;
    result->err  = ex.err * fabs(poly);
    result->err += ex.val * GSL_DBL_EPSILON * (1.0 + fabs(x/b) * (2.0 + fabs(x/(b+1.0))));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_e;
  }
  else if(b == 2*a) {
    return hyperg_1F1_beq2a_pos(a, x, result);  /* 1F1(a,2a,x) */
  }
  else if(   ( b < 10 && a < 10 && ax < MpIeee( "5.0" ) )
          || ( b > a*ax )
          || ( b > a && ax < MpIeee( "5.0" ) )
    ) {
    return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
  }
  else if(b > a && b >= 2*a + x) {
    /* Use the Gautschi CF series, then
     * recurse backward to a=0 for normalization.
     * This will work for either sign of x.
     */
    MpIeee rap;
    int  stat_CF1=  hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    MpIeee ra=  MpIeee( "1.0" ) + x/a * rap;
    MpIeee Ma=  GSL_SQRT_DBL_MIN;
    MpIeee Map1=  ra * Ma;
    MpIeee Mnp1=  Map1;
    MpIeee Mn=  Ma;
    MpIeee Mnm1;
    int  n;
    for(n=a; n>0; n--) {
      Mnm1 = (n * Mnp1 - (MpIeee( "2" )*n-b+x) * Mn) / (b-n);
      Mnp1 = Mn;
      Mn   = Mnm1;
    }
    result->val = Ma/Mn;
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(a) + 1.0) * fabs(Ma/Mn);
    return stat_CF1;
  }
  else if(b > a && b < 2*a + x && b > x) {
    /* Use the Gautschi series representation of
     * the continued fraction. Then recurse forward
     * to the a=b line for normalization. This will
     * work for either sign of x, although we do need
     * to check for b > x, for when x is positive.
     */
    MpIeee rap;
    int  stat_CF1=  hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    MpIeee ra=  MpIeee( "1.0" ) + x/a * rap;
    gsl_sf_result ex;
    int  stat_ex;

    MpIeee Ma=  GSL_SQRT_DBL_MIN;
    MpIeee Map1=  ra * Ma;
    MpIeee Mnm1=  Ma;
    MpIeee Mn=  Map1;
    MpIeee Mnp1;
    int  n;
    for(n=a+1; n<b; n++) {
      Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
      Mnm1 = Mn;
      Mn   = Mnp1;
    }

    stat_ex = gsl_sf_exp_e(x, &ex);  /* 1F1(b,b,x) */
    result->val  = ex.val * Ma/Mn;
    result->err  = ex.err * fabs(Ma/Mn);
    result->err += 4.0 * GSL_DBL_EPSILON * (fabs(b-a)+1.0) * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_ex, stat_CF1);
  }
  else if(x >= 0.0) {

    if(b < a) {
      /* The point b,b is below the b=2a+x line.
       * Forward recursion on a from b,b+1 is possible.
       * Note that a > b + 1 as well, since we already tried a = b + 1.
       */
      if(x + log(fabs(x/b)) < GSL_LOG_DBL_MAX-2.0) {
        MpIeee ex=  exp(x);
        int  n;
        MpIeee Mnm1=  ex;                 /* 1F1(b,b,x)   */
        MpIeee Mn=  ex * (MpIeee( "1.0" ) + x/b);   /* 1F1(b+1,b,x) */
        MpIeee Mnp1;
        for(n=b+1; n<a; n++) {
          Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
          Mnm1 = Mn;
          Mn   = Mnp1;
        }
        result->val  = Mn;
        result->err  = (x + 1.0) * GSL_DBL_EPSILON * fabs(Mn);
        result->err *= fabs(a-b)+1.0;
        return GSL_SUCCESS;
      }
      else {
        OVERFLOW_ERROR(result);
      }
    }
    else {
      /* b > a
       * b < 2a + x 
       * b <= x (otherwise we would have finished above)
       *
       * Gautschi anomalous convergence region. However, we can
       * recurse forward all the way from a=0,1 because we are
       * always underneath the b=2a+x line.
       */
      gsl_sf_result r_Mn;
      MpIeee Mnm1=  MpIeee( "1.0" );    /* 1F1(0,b,x) */
      MpIeee Mn;            /* 1F1(1,b,x)  */
      MpIeee Mnp1;
      int  n;
      gsl_sf_exprel_n_e(b-1, x, &r_Mn);
      Mn = r_Mn.val;
      for(n=1; n<a; n++) {
        Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
        Mnm1 = Mn;
        Mn   = Mnp1;
      }
      result->val  = Mn;
      result->err  = fabs(Mn) * (1.0 + fabs(a)) * fabs(r_Mn.err / r_Mn.val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(Mn);
      return GSL_SUCCESS;
    }
  }
  else {
    /* x < 0
     * b < a (otherwise we would have tripped one of the above)
     */

    if(a <= 0.5*(b-x) || a >= -x) {
      /* Gautschi continued fraction is in the anomalous region,
       * so we must find another way. We recurse down in b,
       * from the a=b line.
       */
      MpIeee ex=  exp(x);
      MpIeee Manp1=  ex;
      MpIeee Man=  ex * (MpIeee( "1.0" ) + x/(a-MpIeee( "1.0" )));
      MpIeee Manm1;
      int  n;
      for(n=a-1; n>b; n--) {
        Manm1 = (-n*(MpIeee( "1" )-n-x)*Man - x*(n-a)*Manp1)/(n*(n-MpIeee( "1.0" )));
        Manp1 = Man;
        Man = Manm1;
      }
      result->val  = Man;
      result->err  = (fabs(x) + 1.0) * GSL_DBL_EPSILON * fabs(Man);
      result->err *= fabs(b-a)+1.0;
      return GSL_SUCCESS;
    }
    else {
      /* Pick a0 such that b ~= 2a0 + x, then
       * recurse down in b from a0,a0 to determine
       * the values near the line b=2a+x. Then recurse
       * forward on a from a0.
       */
      int  a0=  ceil(0.5*(b-x));
      MpIeee Ma0b;    /* M(a0,b)   */
      MpIeee Ma0bp1;  /* M(a0,b+1) */
      MpIeee Ma0p1b;  /* M(a0+1,b) */
      MpIeee Mnm1;
      MpIeee Mn;
      MpIeee Mnp1;
      int  n;
      {
        MpIeee ex=  exp(x);
        MpIeee Ma0np1=  ex;
        MpIeee Ma0n=  ex * (MpIeee( "1.0" ) + x/(a0-MpIeee( "1.0" )));
        MpIeee Ma0nm1;
        for(n=a0-1; n>b; n--) {
          Ma0nm1 = (-n*(MpIeee( "1" )-n-x)*Ma0n - x*(n-a0)*Ma0np1)/(n*(n-MpIeee( "1.0" )));
          Ma0np1 = Ma0n;
          Ma0n = Ma0nm1;
        }
        Ma0bp1 = Ma0np1;
        Ma0b   = Ma0n;
        Ma0p1b = (b*(a0+x)*Ma0b + x*(a0-b)*Ma0bp1)/(a0*b);
      }

      /* Initialise the recurrence correctly BJG */

      if (a0 >= a)
        { 
          Mn = Ma0b;
        }
      else if (a0 + 1>= a)
        {
          Mn = Ma0p1b;
        }
      else
        {
          Mnm1 = Ma0b;
          Mn   = Ma0p1b;

          for(n=a0+1; n<a; n++) {
            Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
            Mnm1 = Mn;
            Mn   = Mnp1;
          }
        }

      result->val  = Mn;
      result->err  = (fabs(x) + 1.0) * GSL_DBL_EPSILON *  fabs(Mn);
      result->err *= fabs(b-a)+1.0;
      return GSL_SUCCESS;
    }
  }
}


/* Evaluate a <= 0, a integer, cases directly. (Polynomial; Horner)
 * When the terms are all positive, this
 * must work. We will assume this here.
 */
static
int
 hyperg_1F1_a_negint_poly(const int a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  if(a == 0) {
    result->val = 1.0;
    result->err = 1.0;
    return GSL_SUCCESS;
  }
  else {
    int  N=  -a;
    MpIeee poly=  MpIeee( "1.0" );
    int  k;
    for(k=N-1; k>=0; k--) {
      MpIeee t=  (a+k)/(b+k) * (x/(k+MpIeee( "1" )));
      MpIeee r=  t + MpIeee( "1.0" )/poly;
      if(r > MpIeee( "0.9" )*GSL_DBL_MAX/poly) {
        OVERFLOW_ERROR(result);
      }
      else {
        poly *= r;  /* P_n = 1 + t_n P_{n-1} */
      }
    }
    result->val = poly;
    result->err = 2.0 * (sqrt(N) + 1.0) * GSL_DBL_EPSILON * fabs(poly);
    return GSL_SUCCESS;
  }
}


/* Evaluate negative integer a case by relation
 * to Laguerre polynomials. This is more general than
 * the direct polynomial evaluation, but is safe
 * for all values of x.
 *
 * 1F1(-n,b,x) = n!/(b)_n Laguerre[n,b-1,x]
 *             = n B(b,n) Laguerre[n,b-1,x]
 *
 * assumes b is not a negative integer
 */
static
int
 hyperg_1F1_a_negint_lag(const int a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  const int n = -a;

  gsl_sf_result lag;
  const int stat_l = gsl_sf_laguerre_n_e(n, b-1.0, x, &lag);
  if(b < 0.0) {
    gsl_sf_result lnfact;
    gsl_sf_result lng1;
    gsl_sf_result lng2;
    MpIeee s1;MpIeee  s2;
    const int stat_f  = gsl_sf_lnfact_e(n, &lnfact);
    const int stat_g1 = gsl_sf_lngamma_sgn_e(b + n, &lng1, &s1);
    const int stat_g2 = gsl_sf_lngamma_sgn_e(b, &lng2, &s2);
    const MpIeee lnpre_val=  lnfact.val - (lng1.val - lng2.val);
    const MpIeee lnpre_err=  lnfact.err + lng1.err + lng2.err
      + 2.0 * GSL_DBL_EPSILON * fabs(lnpre_val);
    const int stat_e = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err,
                                                s1*s2*lag.val, lag.err,
                                                result);
    return GSL_ERROR_SELECT_5(stat_e, stat_l, stat_g1, stat_g2, stat_f);
  }
  else {
    gsl_sf_result lnbeta;
    gsl_sf_lnbeta_e(b, n, &lnbeta);
    if(fabs(lnbeta.val) < 0.1) {
      /* As we have noted, when B(x,y) is near 1,
       * evaluating log(B(x,y)) is not accurate.
       * Instead we evaluate B(x,y) directly.
       */
      const MpIeee ln_term_val=  log(1.25*n);
      const MpIeee ln_term_err=  2.0 * GSL_DBL_EPSILON * ln_term_val;
      gsl_sf_result beta;
      int  stat_b=  gsl_sf_beta_e(b, n, &beta);
      int  stat_e=  gsl_sf_exp_mult_err_e(ln_term_val, ln_term_err,
                                            lag.val, lag.err,
                                            result);
      result->val *= beta.val/1.25;
      result->err *= beta.val/1.25;
      return GSL_ERROR_SELECT_3(stat_e, stat_l, stat_b);
    }
    else {
      /* B(x,y) was not near 1, so it is safe to use
       * the logarithmic values.
       */
      const MpIeee ln_n=  log(n);
      const MpIeee ln_term_val=  lnbeta.val + ln_n;
      const MpIeee ln_term_err=  lnbeta.err + 2.0 * GSL_DBL_EPSILON * fabs(ln_n);
      int  stat_e=  gsl_sf_exp_mult_err_e(ln_term_val, ln_term_err,
                                            lag.val, lag.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_l);
    }
  }
}


/* Handle negative integer a case for x > 0 and
 * generic b.
 *
 * Combine [Abramowitz+Stegun, 13.6.9 + 13.6.27]
 * M(-n,b,x) = (-1)^n / (b)_n U(-n,b,x) = n! / (b)_n Laguerre^(b-1)_n(x)
 */
#if 0
static
int
 hyperg_1F1_a_negint_U(const int a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  const int n = -a;
  const MpIeee sgn=  ( GSL_IS_ODD(n) ? -1.0 : 1.0 );
  MpIeee sgpoch;
  gsl_sf_result lnpoch;
  gsl_sf_result U;
  const int stat_p = gsl_sf_lnpoch_sgn_e(b, n, &lnpoch, &sgpoch);
  const int stat_U = gsl_sf_hyperg_U_e(-n, b, x, &U);
  const int stat_e = gsl_sf_exp_mult_err_e(-lnpoch.val, lnpoch.err,
                                              sgn * sgpoch * U.val, U.err,
                                              result);
  return GSL_ERROR_SELECT_3(stat_e, stat_U, stat_p);
}
#endif


/* Assumes a <= -1,  b <= -1, and b <= a.
 */
static
int
 hyperg_1F1_ab_negint(const int a, const int b, const MpIeee x, gsl_sf_result * result)
{
  if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x > 0.0) {
    return hyperg_1F1_a_negint_poly(a, b, x, result);
  }
  else {
    /* Apply a Kummer transformation to make x > 0 so
     * we can evaluate the polynomial safely. Of course,
     * this assumes b <= a, which must be true for
     * a<0 and b<0, since otherwise the thing is undefined.
     */
    gsl_sf_result K;
    int  stat_K=  hyperg_1F1_a_negint_poly(b-a, b, -x, &K);
    int  stat_e=  gsl_sf_exp_mult_err_e(x, 2.0 * GSL_DBL_EPSILON * fabs(x),
                                          K.val, K.err,
                                          result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K);
  }
}


/* [Abramowitz+Stegun, 13.1.3]
 *
 * M(a,b,x) = Gamma(1+a-b)/Gamma(2-b) x^(1-b) *
 *            { Gamma(b)/Gamma(a) M(1+a-b,2-b,x) - (b-1) U(1+a-b,2-b,x) }
 *
 * b not an integer >= 2
 * a-b not a negative integer
 */
static
int
 hyperg_1F1_U(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  const MpIeee bp=  2.0 - b;
  const MpIeee ap=  a - b + 1.0;

  gsl_sf_result lg_ap, lg_bp;
  MpIeee sg_ap;
  int  stat_lg0=  gsl_sf_lngamma_sgn_e(ap, &lg_ap, &sg_ap);
  int  stat_lg1=  gsl_sf_lngamma_e(bp, &lg_bp);
  int  stat_lg2=  GSL_ERROR_SELECT_2(stat_lg0, stat_lg1);
  MpIeee t1=  (bp-MpIeee( "1.0" )) * log(x);
  MpIeee lnpre_val=  lg_ap.val - lg_bp.val + t1;
  MpIeee lnpre_err=  lg_ap.err + lg_bp.err + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(t1);

  gsl_sf_result lg_2mbp, lg_1papmbp;
  MpIeee sg_2mbp;MpIeee  sg_1papmbp;
  int  stat_lg3=  gsl_sf_lngamma_sgn_e(2.0-bp,    &lg_2mbp,    &sg_2mbp);
  int  stat_lg4=  gsl_sf_lngamma_sgn_e(1.0+ap-bp, &lg_1papmbp, &sg_1papmbp);
  int  stat_lg5=  GSL_ERROR_SELECT_2(stat_lg3, stat_lg4);
  MpIeee lnc1_val=  lg_2mbp.val - lg_1papmbp.val;
  MpIeee lnc1_err=  lg_2mbp.err + lg_1papmbp.err
                    + GSL_DBL_EPSILON * (fabs(lg_2mbp.val) + fabs(lg_1papmbp.val));

  gsl_sf_result M;
  gsl_sf_result_e10 U;
  int  stat_F=  gsl_sf_hyperg_1F1_e(ap, bp, x, &M);
  int  stat_U=  gsl_sf_hyperg_U_e10_e(ap, bp, x, &U);
  int  stat_FU=  GSL_ERROR_SELECT_2(stat_F, stat_U);

  gsl_sf_result_e10 term_M;
  int  stat_e0=  gsl_sf_exp_mult_err_e10_e(lnc1_val, lnc1_err,
                                             sg_2mbp*sg_1papmbp*M.val, M.err,
                                             &term_M);

  const MpIeee ombp=  1.0 - bp;
  const MpIeee Uee_val=  U.e10*M_LN10;
  const MpIeee Uee_err=  2.0 * GSL_DBL_EPSILON * fabs(Uee_val);
  const MpIeee Mee_val=  term_M.e10*M_LN10;
  const MpIeee Mee_err=  2.0 * GSL_DBL_EPSILON * fabs(Mee_val);
  int  stat_e1;

  /* Do a little dance with the exponential prefactors
   * to avoid overflows in intermediate results.
   */
  if(Uee_val > Mee_val) {
    const MpIeee factorM_val=  exp(Mee_val-Uee_val);
    const MpIeee factorM_err=  2.0 * GSL_DBL_EPSILON * (fabs(Mee_val-Uee_val)+1.0) * factorM_val;
    const MpIeee inner_val=  term_M.val*factorM_val - ombp*U.val;
    const MpIeee inner_err= 
        term_M.err*factorM_val + fabs(ombp) * U.err
      + fabs(term_M.val) * factorM_err
      + GSL_DBL_EPSILON * (fabs(term_M.val*factorM_val) + fabs(ombp*U.val));
    stat_e1 = gsl_sf_exp_mult_err_e(lnpre_val+Uee_val, lnpre_err+Uee_err,
                                       sg_ap*inner_val, inner_err,
                                       result);
  }
  else {
    const MpIeee factorU_val=  exp(Uee_val - Mee_val);
    const MpIeee factorU_err=  2.0 * GSL_DBL_EPSILON * (fabs(Mee_val-Uee_val)+1.0) * factorU_val;
    const MpIeee inner_val=  term_M.val - ombp*factorU_val*U.val;
    const MpIeee inner_err= 
        term_M.err + fabs(ombp*factorU_val*U.err)
      + fabs(ombp*factorU_err*U.val)
      + GSL_DBL_EPSILON * (fabs(term_M.val) + fabs(ombp*factorU_val*U.val));
    stat_e1 = gsl_sf_exp_mult_err_e(lnpre_val+Mee_val, lnpre_err+Mee_err,
                                       sg_ap*inner_val, inner_err,
                                       result);
  }

  return GSL_ERROR_SELECT_5(stat_e1, stat_e0, stat_FU, stat_lg5, stat_lg2);
}


/* Handle case of generic positive a, b.
 * Assumes b-a is not a negative integer.
 */
static
int
 hyperg_1F1_ab_pos(const MpIeee a, const MpIeee b,
                  const MpIeee x,
                  gsl_sf_result * result)
{
  const MpIeee ax=  fabs(x);

  if(   ( b < 10.0 && a < 10.0 && ax < 5.0 )
     || ( b > a*ax )
     || ( b > a && ax < 5.0 )
    ) {
    return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
  }
  else if(   x < -100.0
          && GSL_MAX_DBL(fabs(a),1.0)*GSL_MAX_DBL(fabs(1.0+a-b),1.0) < 0.7*fabs(x)
    ) {
    /* Large negative x asymptotic.
     */
    return hyperg_1F1_asymp_negx(a, b, x, result);
  }
  else if(   x > 100.0
          && GSL_MAX_DBL(fabs(b-a),1.0)*GSL_MAX_DBL(fabs(1.0-a),1.0) < 0.7*fabs(x)
    ) {
    /* Large positive x asymptotic.
     */
    return hyperg_1F1_asymp_posx(a, b, x, result);
  }
  else if(fabs(b-a) <= 1.0) {
    /* Directly handle b near a.
     */
    return hyperg_1F1_beps_bgt0(a-b, b, x, result);  /* a = b + eps */
  }

  else if(b > a && b >= 2*a + x) {
    /* Use the Gautschi CF series, then
     * recurse backward to a near 0 for normalization.
     * This will work for either sign of x.
     */ 
    MpIeee rap;
    int  stat_CF1=  hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    MpIeee ra=  MpIeee( "1.0" ) + x/a * rap;

    MpIeee Ma=  GSL_SQRT_DBL_MIN;
    MpIeee Map1=  ra * Ma;
    MpIeee Mnp1=  Map1;
    MpIeee Mn=  Ma;
    MpIeee Mnm1;
    gsl_sf_result Mn_true;
    int  stat_Mt;
    MpIeee n;
    for(n=a; n>MpIeee( "0.5" ); n -= MpIeee( "1.0" )) {
      Mnm1 = (n * Mnp1 - (MpIeee( "2.0" )*n-b+x) * Mn) / (b-n);
      Mnp1 = Mn;
      Mn   = Mnm1;
    }

    stat_Mt = hyperg_1F1_small_a_bgt0(n, b, x, &Mn_true);

    result->val  = (Ma/Mn) * Mn_true.val;
    result->err  = fabs(Ma/Mn) * Mn_true.err;
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(a)+1.0) * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_Mt, stat_CF1);
  }
  else if(b > a && b < 2*a + x && b > x) {
    /* Use the Gautschi series representation of
     * the continued fraction. Then recurse forward
     * to near the a=b line for normalization. This will
     * work for either sign of x, although we do need
     * to check for b > x, which is relevant when x is positive.
     */
    gsl_sf_result Mn_true;
    int  stat_Mt;
    MpIeee rap;
    int  stat_CF1=  hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    MpIeee ra=  MpIeee( "1.0" ) + x/a * rap;
    MpIeee Ma=  GSL_SQRT_DBL_MIN;
    MpIeee Mnm1=  Ma;
    MpIeee Mn=  ra * Mnm1;
    MpIeee Mnp1;
    MpIeee n;
    for(n=a+MpIeee( "1.0" ); n<b-MpIeee( "0.5" ); n += MpIeee( "1.0" )) {
      Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
      Mnm1 = Mn;
      Mn   = Mnp1;
    }
    stat_Mt = hyperg_1F1_beps_bgt0(n-b, b, x, &Mn_true);
    result->val  = Ma/Mn * Mn_true.val;
    result->err  = fabs(Ma/Mn) * Mn_true.err;
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(b-a)+1.0) * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_Mt, stat_CF1);
  }
  else if(x >= 0.0) {

    if(b < a) {
      /* Forward recursion on a from a=b+eps-1,b+eps.
       */
      MpIeee N=  floor(a-b);
      MpIeee eps=  a - b - N;
      gsl_sf_result r_M0;
      gsl_sf_result r_M1;
      int  stat_0=  hyperg_1F1_beps_bgt0(eps-1.0, b, x, &r_M0);
      int  stat_1=  hyperg_1F1_beps_bgt0(eps,     b, x, &r_M1);
      MpIeee M0=  r_M0.val;
      MpIeee M1=  r_M1.val;

      MpIeee Mam1=  M0;
      MpIeee Ma=  M1;
      MpIeee Map1;
      MpIeee ap;
      MpIeee start_pair=  fabs(M0) + fabs(M1);
      MpIeee minim_pair=  GSL_DBL_MAX;
      MpIeee pair_ratio;
      MpIeee rat_0=  fabs(r_M0.err/r_M0.val);
      MpIeee rat_1=  fabs(r_M1.err/r_M1.val);
      for(ap=b+eps; ap<a-MpIeee( "0.1" ); ap += MpIeee( "1.0" )) {
        Map1 = ((b-ap)*Mam1 + (MpIeee( "2.0" )*ap-b+x)*Ma)/ap;
        Mam1 = Ma;
        Ma   = Map1;
        minim_pair = GSL_MIN_DBL(fabs(Mam1) + fabs(Ma), minim_pair);
      }
      pair_ratio = start_pair/minim_pair;
      result->val  = Ma;
      result->err  = 2.0 * (rat_0 + rat_1 + GSL_DBL_EPSILON) * (fabs(b-a)+1.0) * fabs(Ma);
      result->err += 2.0 * (rat_0 + rat_1) * pair_ratio*pair_ratio * fabs(Ma);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(Ma);
      return GSL_ERROR_SELECT_2(stat_0, stat_1);
    }
    else {
      /* b > a
       * b < 2a + x 
       * b <= x
       *
       * Recurse forward on a from a=eps,eps+1.
       */
      MpIeee eps=  a - floor(a);
      gsl_sf_result r_Mnm1;
      gsl_sf_result r_Mn;
      int  stat_0=  hyperg_1F1_small_a_bgt0(eps,     b, x, &r_Mnm1);
      int  stat_1=  hyperg_1F1_small_a_bgt0(eps+1.0, b, x, &r_Mn);
      MpIeee Mnm1=  r_Mnm1.val;
      MpIeee Mn=  r_Mn.val;
      MpIeee Mnp1;

      MpIeee n;
      MpIeee start_pair=  fabs(Mn) + fabs(Mnm1);
      MpIeee minim_pair=  GSL_DBL_MAX;
      MpIeee pair_ratio;
      MpIeee rat_0=  fabs(r_Mnm1.err/r_Mnm1.val);
      MpIeee rat_1=  fabs(r_Mn.err/r_Mn.val);
      for(n=eps+MpIeee( "1.0" ); n<a-MpIeee( "0.1" ); n++) {
        Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
        Mnm1 = Mn;
        Mn   = Mnp1;
        minim_pair = GSL_MIN_DBL(fabs(Mn) + fabs(Mnm1), minim_pair);
      }
      pair_ratio = start_pair/minim_pair;
      result->val  = Mn;
      result->err  = 2.0 * (rat_0 + rat_1 + GSL_DBL_EPSILON) * (fabs(a)+1.0) * fabs(Mn);
      result->err += 2.0 * (rat_0 + rat_1) * pair_ratio*pair_ratio * fabs(Mn);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(Mn);
      return GSL_ERROR_SELECT_2(stat_0, stat_1);
    }
  }
  else {
    /* x < 0
     * b < a
     */

    if(a <= 0.5*(b-x) || a >= -x) {
      /* Recurse down in b, from near the a=b line, b=a+eps,a+eps-1.
       */
      MpIeee N=  floor(a - b);
      MpIeee eps=  MpIeee( "1.0" ) + N - a + b;
      gsl_sf_result r_Manp1;
      gsl_sf_result r_Man;
      int  stat_0=  hyperg_1F1_beps_bgt0(-eps,    a+eps,     x, &r_Manp1);
      int  stat_1=  hyperg_1F1_beps_bgt0(1.0-eps, a+eps-1.0, x, &r_Man);
      MpIeee Manp1=  r_Manp1.val;
      MpIeee Man=  r_Man.val;
      MpIeee Manm1;

      MpIeee n;
      MpIeee start_pair=  fabs(Manp1) + fabs(Man);
      MpIeee minim_pair=  GSL_DBL_MAX;
      MpIeee pair_ratio;
      MpIeee rat_0=  fabs(r_Manp1.err/r_Manp1.val);
      MpIeee rat_1=  fabs(r_Man.err/r_Man.val);
      for(n=a+eps-MpIeee( "1.0" ); n>b+MpIeee( "0.1" ); n -= MpIeee( "1.0" )) {
        Manm1 = (-n*(MpIeee( "1" )-n-x)*Man - x*(n-a)*Manp1)/(n*(n-MpIeee( "1.0" )));
        Manp1 = Man;
        Man = Manm1;
        minim_pair = GSL_MIN_DBL(fabs(Manp1) + fabs(Man), minim_pair);
      }

      /* FIXME: this is a nasty little hack; there is some
         (transient?) instability in this recurrence for some
         values. I can tell when it happens, which is when
         this pair_ratio is large. But I do not know how to
         measure the error in terms of it. I guessed quadratic
         below, but it is probably worse than that.
         */
      pair_ratio = start_pair/minim_pair;
      result->val  = Man;
      result->err  = 2.0 * (rat_0 + rat_1 + GSL_DBL_EPSILON) * (fabs(b-a)+1.0) * fabs(Man);
      result->err *= pair_ratio*pair_ratio + 1.0;
      return GSL_ERROR_SELECT_2(stat_0, stat_1);
    }
    else {
      /* Pick a0 such that b ~= 2a0 + x, then
       * recurse down in b from a0,a0 to determine
       * the values near the line b=2a+x. Then recurse
       * forward on a from a0.
       */
      MpIeee epsa=  a - floor(a);
      MpIeee a0=  floor(MpIeee( "0.5" )*(b-x)) + epsa;
      MpIeee N=  floor(a0 - b);
      MpIeee epsb=  MpIeee( "1.0" ) + N - a0 + b;
      MpIeee Ma0b;
      MpIeee Ma0bp1;
      MpIeee Ma0p1b;
      int  stat_a0;
      MpIeee Mnm1;
      MpIeee Mn;
      MpIeee Mnp1;
      MpIeee n;
      MpIeee err_rat;
      {
        gsl_sf_result r_Ma0np1;
        gsl_sf_result r_Ma0n;
        int  stat_0=  hyperg_1F1_beps_bgt0(-epsb,    a0+epsb,     x, &r_Ma0np1);
        int  stat_1=  hyperg_1F1_beps_bgt0(1.0-epsb, a0+epsb-1.0, x, &r_Ma0n);
        MpIeee Ma0np1=  r_Ma0np1.val;
        MpIeee Ma0n=  r_Ma0n.val;
        MpIeee Ma0nm1;

        err_rat = fabs(r_Ma0np1.err/r_Ma0np1.val) + fabs(r_Ma0n.err/r_Ma0n.val);

        for(n=a0+epsb-MpIeee( "1.0" ); n>b+MpIeee( "0.1" ); n -= MpIeee( "1.0" )) {
          Ma0nm1 = (-n*(MpIeee( "1" )-n-x)*Ma0n - x*(n-a0)*Ma0np1)/(n*(n-MpIeee( "1.0" )));
          Ma0np1 = Ma0n;
          Ma0n = Ma0nm1;
        }
        Ma0bp1 = Ma0np1;
        Ma0b   = Ma0n;
        Ma0p1b = (b*(a0+x)*Ma0b+x*(a0-b)*Ma0bp1)/(a0*b); /* right-down hook */
        stat_a0 = GSL_ERROR_SELECT_2(stat_0, stat_1);
      }

          
      /* Initialise the recurrence correctly BJG */

      if (a0 >= a - MpIeee( "0.1" ))
        { 
          Mn = Ma0b;
        }
      else if (a0 + 1>= a - 0.1)
        {
          Mn = Ma0p1b;
        }
      else
        {
          Mnm1 = Ma0b;
          Mn   = Ma0p1b;

          for(n=a0+MpIeee( "1.0" ); n<a-MpIeee( "0.1" ); n += MpIeee( "1.0" )) {
            Mnp1 = ((b-n)*Mnm1 + (MpIeee( "2" )*n-b+x)*Mn)/n;
            Mnm1 = Mn;
            Mn   = Mnp1;
          }
        }

      result->val = Mn;
      result->err = (err_rat + GSL_DBL_EPSILON) * (fabs(b-a)+1.0) * fabs(Mn);
      return stat_a0;
    }
  }
}


/* Assumes b != integer
 * Assumes a != integer when x > 0
 * Assumes b-a != neg integer when x < 0
 */
static
int
 hyperg_1F1_ab_neg(const MpIeee a, const MpIeee b, const MpIeee x,
                  gsl_sf_result * result)
{
  const MpIeee bma=  b - a;
  const MpIeee abs_x=  fabs(x);
  const MpIeee abs_a=  fabs(a);
  const MpIeee abs_b=  fabs(b);
  const MpIeee size_a=  GSL_MAX(abs_a, 1.0);
  const MpIeee size_b=  GSL_MAX(abs_b, 1.0);
  const int bma_integer = ( bma - floor(bma+0.5) < _1F1_INT_THRESHOLD );

  if(   (abs_a < 10.0 && abs_b < 10.0 && abs_x < 5.0)
     || (b > 0.8*GSL_MAX(fabs(a),1.0)*fabs(x))
    ) {
    return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
  }
  else if(   x > 0.0
          && size_b > size_a
          && size_a*log(M_E*x/size_b) < GSL_LOG_DBL_EPSILON+7.0
    ) {
    /* Series terms are positive definite up until
     * there is a sign change. But by then the
     * terms are small due to the last condition.
     */
    return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
  }
  else if(   (abs_x < 5.0 && fabs(bma) < 10.0 && abs_b < 10.0)
          || (b > 0.8*GSL_MAX_DBL(fabs(bma),1.0)*abs_x)
    ) {
    /* Use Kummer transformation to render series safe.
     */
    gsl_sf_result Kummer_1F1;
    int  stat_K=  gsl_sf_hyperg_1F1_series_e(bma, b, -x, &Kummer_1F1);
    int  stat_e=  gsl_sf_exp_mult_err_e(x, GSL_DBL_EPSILON * fabs(x),
                                      Kummer_1F1.val, Kummer_1F1.err,
                                      result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K);
  }
  else if(   x < -30.0
          && GSL_MAX_DBL(fabs(a),1.0)*GSL_MAX_DBL(fabs(1.0+a-b),1.0) < 0.99*fabs(x)
    ) {
    /* Large negative x asymptotic.
     * Note that we do not check if b-a is a negative integer.
     */
    return hyperg_1F1_asymp_negx(a, b, x, result);
  }
  else if(   x > 100.0
          && GSL_MAX_DBL(fabs(bma),1.0)*GSL_MAX_DBL(fabs(1.0-a),1.0) < 0.99*fabs(x)
    ) {
    /* Large positive x asymptotic.
     * Note that we do not check if a is a negative integer.
     */
    return hyperg_1F1_asymp_posx(a, b, x, result);
  }
  else if(x > 0.0 && !(bma_integer && bma > 0.0)) {
    return hyperg_1F1_U(a, b, x, result);
  }
  else {
    /* FIXME:  if all else fails, try the series... BJG */
    if (x < 0.0) {
      /* Apply Kummer Transformation */
      int  status=  gsl_sf_hyperg_1F1_series_e(b-a, b, -x, result);
      MpIeee K_factor=  exp(x);
      result->val *= K_factor;
      result->err *= K_factor;
      return status;
    } else {
      int  status=  gsl_sf_hyperg_1F1_series_e(a, b, x, result);
      return status;
    }

    /* Sadness... */
    /* result->val = 0.0; */
    /* result->err = 0.0; */
    /* GSL_ERROR ("error", GSL_EUNIMPL); */
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_hyperg_1F1_int_e(const int a, const int b, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == b) {
    return gsl_sf_exp_e(x, result);
  }
  else if(b == 0) {
    DOMAIN_ERROR(result);
  }
  else if(a == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(b < 0 && (a < b || a > 0)) {
    /* Standard domain error due to singularity. */
    DOMAIN_ERROR(result);
  }
  else if(x > 100.0  && GSL_MAX_DBL(1.0,fabs(b-a))*GSL_MAX_DBL(1.0,fabs(1-a)) < 0.5 * x) {
    /* x -> +Inf asymptotic */
    return hyperg_1F1_asymp_posx(a, b, x, result);
  }
  else if(x < -100.0 && GSL_MAX_DBL(1.0,fabs(a))*GSL_MAX_DBL(1.0,fabs(1+a-b)) < 0.5 * fabs(x)) {
    /* x -> -Inf asymptotic */
    return hyperg_1F1_asymp_negx(a, b, x, result);
  }
  else if(a < 0 && b < 0) {
    return hyperg_1F1_ab_negint(a, b, x, result);
  }
  else if(a < 0 && b > 0) {
    /* Use Kummer to reduce it to the positive integer case.
     * Note that b > a, strictly, since we already trapped b = a.
     */
    gsl_sf_result Kummer_1F1;
    int  stat_K=  hyperg_1F1_ab_posint(b-a, b, -x, &Kummer_1F1);
    int  stat_e=  gsl_sf_exp_mult_err_e(x, GSL_DBL_EPSILON * fabs(x),
                                      Kummer_1F1.val, Kummer_1F1.err,
                                      result); 
    return GSL_ERROR_SELECT_2(stat_e, stat_K);
  }
  else {
    /* a > 0 and b > 0 */
    return hyperg_1F1_ab_posint(a, b, x, result);
  }
}


int
 gsl_sf_hyperg_1F1_e(const MpIeee a, const MpIeee b, const MpIeee x,
                       gsl_sf_result * result
                       )
{
  const MpIeee bma=  b - a;
  const MpIeee rinta=  floor(a + 0.5);
  const MpIeee rintb=  floor(b + 0.5);
  const MpIeee rintbma=  floor(bma + 0.5);
  const int a_integer   = ( fabs(a-rinta) < _1F1_INT_THRESHOLD && rinta > INT_MIN && rinta < INT_MAX );
  const int b_integer   = ( fabs(b-rintb) < _1F1_INT_THRESHOLD && rintb > INT_MIN && rintb < INT_MAX );
  const int bma_integer = ( fabs(bma-rintbma) < _1F1_INT_THRESHOLD && rintbma > INT_MIN && rintbma < INT_MAX );
  const int b_neg_integer   = ( b < -0.1 && b_integer );
  const int a_neg_integer   = ( a < -0.1 && a_integer );
  const int bma_neg_integer = ( bma < -0.1 &&  bma_integer );

  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    /* Testing for this before testing a and b
     * is somewhat arbitrary. The result is that
     * we have 1F1(a,0,0) = 1.
     */
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(b == 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(a == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == b) {
    /* case: a==b; exp(x)
     * It's good to test exact equality now.
     * We also test approximate equality later.
     */
    return gsl_sf_exp_e(x, result);
  }
  else if(fabs(b) < _1F1_INT_THRESHOLD) {
    /* Note that neither a nor b is zero, since
     * we eliminated that with the above tests.
     */
    if(fabs(a) < _1F1_INT_THRESHOLD) {
      /* a and b near zero: 1 + a/b (exp(x)-1)
       */
      gsl_sf_result exm1;
      int  stat_e=  gsl_sf_expm1_e(x, &exm1);
      MpIeee sa=  ( a > MpIeee( "0.0" ) ? MpIeee( "1.0" ) : -MpIeee( "1.0" ) );
      MpIeee sb=  ( b > MpIeee( "0.0" ) ? MpIeee( "1.0" ) : -MpIeee( "1.0" ) );
      MpIeee lnab=  log(fabs(a/b)); /* safe */
      gsl_sf_result hx;
      int  stat_hx=  gsl_sf_exp_mult_err_e(lnab, GSL_DBL_EPSILON * fabs(lnab),
                                             sa * sb * exm1.val, exm1.err,
                                             &hx);
      result->val = (hx.val == GSL_DBL_MAX ? hx.val : 1.0 + hx.val);  /* FIXME: excessive paranoia ? what is DBL_MAX+1 ?*/
      result->err = hx.err;
      return GSL_ERROR_SELECT_2(stat_hx, stat_e);
    }
    else {
      /* b near zero and a not near zero
       */
      const MpIeee m_arg=  1.0/(0.5*b);
      gsl_sf_result F_renorm;
      int  stat_F=  hyperg_1F1_renorm_b0(a, x, &F_renorm);
      int  stat_m=  gsl_sf_multiply_err_e(m_arg, 2.0 * GSL_DBL_EPSILON * m_arg,
                                            0.5*F_renorm.val, 0.5*F_renorm.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_m, stat_F);
    }
  }
  else if(a_integer && b_integer) {
    /* Check for reduction to the integer case.
     * Relies on the arbitrary "near an integer" test.
     */
    return gsl_sf_hyperg_1F1_int_e((int)rinta, (int)rintb, x, result);
  }
  else if(b_neg_integer && !(a_neg_integer && a > b)) {
    /* Standard domain error due to
     * uncancelled singularity.
     */
    DOMAIN_ERROR(result);
  }
  else if(a_neg_integer) {
    return hyperg_1F1_a_negint_lag((int)rinta, b, x, result);
  }
  else if(b > 0.0) {
    if(-1.0 <= a && a <= 1.0) {
      /* Handle small a explicitly.
       */
      return hyperg_1F1_small_a_bgt0(a, b, x, result);
    }
    else if(bma_neg_integer) {
      /* Catch this now, to avoid problems in the
       * generic evaluation code.
       */
      gsl_sf_result Kummer_1F1;
      int  stat_K=  hyperg_1F1_a_negint_lag((int)rintbma, b, -x, &Kummer_1F1);
      int  stat_e=  gsl_sf_exp_mult_err_e(x, GSL_DBL_EPSILON * fabs(x),
                                            Kummer_1F1.val, Kummer_1F1.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_K);
    }
    else if(a < 0.0) {
      /* Use Kummer to reduce it to the generic positive case.
       * Note that b > a, strictly, since we already trapped b = a.
       * Also b-(b-a)=a, and a is not a negative integer here,
       * so the generic evaluation is safe.
       */
      gsl_sf_result Kummer_1F1;
      int  stat_K=  hyperg_1F1_ab_pos(b-a, b, -x, &Kummer_1F1);
      int  stat_e=  gsl_sf_exp_mult_err_e(x, GSL_DBL_EPSILON * fabs(x),
                                            Kummer_1F1.val, Kummer_1F1.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_K);
    }
    else {
      /* a > 0.0 */
      return hyperg_1F1_ab_pos(a, b, x, result);
    }
  }
  else {
    /* b < 0.0 */

    if(bma_neg_integer && x < 0.0) {
      /* Handle this now to prevent problems
       * in the generic evaluation.
       */
      gsl_sf_result K;
      int  stat_K;
      int  stat_e;
      if(a < 0.0) {
        /* Kummer transformed version of safe polynomial.
         * The condition a < 0 is equivalent to b < b-a,
         * which is the condition required for the series
         * to be positive definite here.
         */
        stat_K = hyperg_1F1_a_negint_poly((int)rintbma, b, -x, &K);
      }
      else {
        /* Generic eval for negative integer a. */
        stat_K = hyperg_1F1_a_negint_lag((int)rintbma, b, -x, &K);
      }
      stat_e = gsl_sf_exp_mult_err_e(x, GSL_DBL_EPSILON * fabs(x),
                                        K.val, K.err,
                                        result);
      return GSL_ERROR_SELECT_2(stat_e, stat_K);
    }
    else if(a > 0.0) {
      /* Use Kummer to reduce it to the generic negative case.
       */
      gsl_sf_result K;
      int  stat_K=  hyperg_1F1_ab_neg(b-a, b, -x, &K);
      int  stat_e=  gsl_sf_exp_mult_err_e(x, GSL_DBL_EPSILON * fabs(x),
                                            K.val, K.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_K);
    }
    else {
      return hyperg_1F1_ab_neg(a, b, x, result);
    }
  }
}


  
#if 0  
    /* Luke in the canonical case.
   */
  if(x < 0.0 && !a_neg_integer && !bma_neg_integer) {
    MpIeee prec;
    return hyperg_1F1_luke(a, b, x, result, &prec);
  }


  /* Luke with Kummer transformation.
   */
  if(x > 0.0 && !a_neg_integer && !bma_neg_integer) {
    MpIeee prec;
    MpIeee Kummer_1F1;
    MpIeee ex;
    int  stat_F=  hyperg_1F1_luke(b-a, b, -x, &Kummer_1F1, &prec);
    int  stat_e=  gsl_sf_exp_e(x, &ex);
    if(stat_F == GSL_SUCCESS && stat_e == GSL_SUCCESS) {
      MpIeee lnr=  log(fabs(Kummer_1F1)) + x;
      if(lnr < GSL_LOG_DBL_MAX) {
        *result = ex * Kummer_1F1;
        return GSL_SUCCESS;
      }
      else {
        *result = GSL_POSINF; 
        GSL_ERROR ("overflow", GSL_EOVRFLW);
      }
    }
    else if(stat_F != GSL_SUCCESS) {
      *result = 0.0;
      return stat_F;
    }
    else {
      *result = 0.0;
      return stat_e;
    }
  }
#endif



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_hyperg_1F1_int(const int m, const int n, MpIeee x)
{
  EVAL_RESULT(gsl_sf_hyperg_1F1_int_e(m, n, x, &result));
}

MpIeee gsl_sf_hyperg_1F1(MpIeee a, MpIeee b, MpIeee x)
{
  EVAL_RESULT(gsl_sf_hyperg_1F1_e(a, b, x, &result));
}
