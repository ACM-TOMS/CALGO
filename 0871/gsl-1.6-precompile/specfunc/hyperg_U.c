#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/hyperg_U.c
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_hyperg.h>

#include "error.h"
#include "hyperg.h"

#define INT_THRESHOLD (1000.0*GSL_DBL_EPSILON)

#define SERIES_EVAL_OK(a,b,x) ((fabs(a) < 5 && b < 5 && x < 2.0) || (fabs(a) <  10 && b < 10 && x < 1.0))

#define ASYMP_EVAL_OK(a,b,x) (GSL_MAX_DBL(fabs(a),1.0)*GSL_MAX_DBL(fabs(1.0+a-b),1.0) < 0.99*fabs(x))


/* Log[U(a,2a,x)]
 * [Abramowitz+stegun, 13.6.21]
 * Assumes x > 0, a > 1/2.
 */
static
int
 hyperg_lnU_beq2a(const MpIeee a, const MpIeee x, gsl_sf_result * result)
{
  const MpIeee lx=  log(x);
  const MpIeee nu=  a - 0.5;
  const MpIeee lnpre=  0.5*(x - M_LNPI) - nu*lx;
  gsl_sf_result lnK;
  gsl_sf_bessel_lnKnu_e(nu, 0.5*x, &lnK);
  result->val  = lnpre + lnK.val;
  result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(0.5*x) + 0.5*M_LNPI + fabs(nu*lx));
  result->err += lnK.err;
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}


/* Evaluate u_{N+1}/u_N by Steed's continued fraction method.
 *
 * u_N := Gamma[a+N]/Gamma[a] U(a + N, b, x)
 *
 * u_{N+1}/u_N = (a+N) U(a+N+1,b,x)/U(a+N,b,x)
 */
static
int
 hyperg_U_CF1(const MpIeee a, const MpIeee b, const int N, const MpIeee x,
             MpIeee * result, int  * count)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 20000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  -(a + N);
  MpIeee b1=   (b - MpIeee( "2.0" )*a - x - MpIeee( "2.0" )*(N+MpIeee( "1" )));
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
    an = -(a + N + n - b)*(a + N + n - MpIeee( "1.0" ));
    bn =  (b - MpIeee( "2.0" )*a - x - MpIeee( "2.0" )*(N+n));
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
  *count  = n;

  if(n == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Large x asymptotic for  x^a U(a,b,x)
 * Based on SLATEC D9CHU() [W. Fullerton]
 *
 * Uses a rational approximation due to Luke.
 * See [Luke, Algorithms for the Computation of Special Functions, p. 252]
 *     [Luke, Utilitas Math. (1977)]
 *
 * z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
 *
 * This assumes that a is not a negative integer and
 * that 1+a-b is not a negative integer. If one of them
 * is, then the 2F0 actually terminates, the above
 * relation is an equality, and the sum should be
 * evaluated directly [see below].
 */
static
int
 d9chu(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  const MpIeee EPS=  8.0 * GSL_DBL_EPSILON;  /* EPS = 4.0D0*D1MACH(4)   */
  const int maxiter = 500;
  MpIeee aa[4];MpIeee  bb[4];
  int  i;

  MpIeee bp=  MpIeee( "1.0" ) + a - b;
  MpIeee ab=  a*bp;
  MpIeee ct2=  MpIeee( "2.0" ) * (x - ab);
  MpIeee sab=  a + bp;
  
  MpIeee ct3=  sab + MpIeee( "1.0" ) + ab;
  MpIeee anbn=  ct3 + sab + MpIeee( "3.0" );
  MpIeee ct1=  MpIeee( "1.0" ) + MpIeee( "2.0" )*x/anbn;

  bb[0] = MpIeee( "1.0" );
  aa[0] = MpIeee( "1.0" );

  bb[1] = MpIeee( "1.0" ) + MpIeee( "2.0" )*x/ct3;
  aa[1] = MpIeee( "1.0" ) + ct2/ct3;
  
  bb[2] = MpIeee( "1.0" ) + MpIeee( "6.0" )*ct1*x/ct3;
  aa[2] = MpIeee( "1.0" ) + MpIeee( "6.0" )*ab/anbn + MpIeee( "3.0" )*ct1*ct2/ct3;

  for(i=4; i<maxiter; i++) {
    int  j;
    MpIeee c2;
    MpIeee d1z;
    MpIeee g1;MpIeee  g2;MpIeee  g3;
    MpIeee x2i1=  MpIeee( "2" )*i - MpIeee( "3" );
    ct1   = x2i1/(x2i1-MpIeee( "2.0" ));
    anbn += x2i1 + sab;
    ct2   = (x2i1 - MpIeee( "1.0" ))/anbn;
    c2    = x2i1*ct2 - MpIeee( "1.0" );
    d1z   = MpIeee( "2.0" )*x2i1*x/anbn;
    
    ct3 = sab*ct2;
    g1  = d1z + ct1*(c2+ct3);
    g2  = d1z - c2;
    g3  = ct1*(MpIeee( "1.0" ) - ct3 - MpIeee( "2.0" )*ct2);
    
    bb[3] = g1*bb[2] + g2*bb[1] + g3*bb[0];
    aa[3] = g1*aa[2] + g2*aa[1] + g3*aa[0];
    
    if(fabs(aa[3]*bb[0]-aa[0]*bb[3]) < EPS*fabs(bb[3]*bb[0])) break;
    
    for(j=0; j<3; j++) {
      aa[j] = aa[j+1];
      bb[j] = bb[j+1];
    }
  }
  
  result->val = aa[3]/bb[3];
  result->err = 8.0 * GSL_DBL_EPSILON * fabs(result->val);
  
  if(i == maxiter) {
    GSL_ERROR ("error", GSL_EMAXITER);
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Evaluate asymptotic for z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
 * We check for termination of the 2F0 as a special case.
 * Assumes x > 0.
 * Also assumes a,b are not too large compared to x.
 */
static
int
 hyperg_zaU_asymp(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result *result)
{
  const MpIeee ap=  a;
  const MpIeee bp=  1.0 + a - b;
  const MpIeee rintap=  floor(ap + 0.5);
  const MpIeee rintbp=  floor(bp + 0.5);
  const int ap_neg_int = ( ap < 0.0 && fabs(ap - rintap) < INT_THRESHOLD );
  const int bp_neg_int = ( bp < 0.0 && fabs(bp - rintbp) < INT_THRESHOLD );

  if(ap_neg_int || bp_neg_int) {
    /* Evaluate 2F0 polynomial.
     */
    MpIeee mxi=  -MpIeee( "1.0" )/x;
    MpIeee nmax=  -(int)(GSL_MIN(ap,bp) - MpIeee( "0.1" ));
    MpIeee tn=  MpIeee( "1.0" );
    MpIeee sum=  MpIeee( "1.0" );
    MpIeee n=  MpIeee( "1.0" );
    MpIeee sum_err=  MpIeee( "0.0" );
    while(n <= nmax) {
      MpIeee apn=  (ap+n-MpIeee( "1.0" ));
      MpIeee bpn=  (bp+n-MpIeee( "1.0" ));
      tn  *= ((apn/n)*mxi)*bpn;
      sum += tn;
      sum_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(tn);
      n += MpIeee( "1.0" );
    }
    result->val  = sum;
    result->err  = sum_err;
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(nmax)+1.0) * fabs(sum);
    return GSL_SUCCESS;
  }
  else {
    return d9chu(a,b,x,result);
  }
}


/* Evaluate finite sum which appears below.
 */
static
int
 hyperg_U_finite_sum(int  N, MpIeee a, MpIeee b, MpIeee x, MpIeee xeps,
                    gsl_sf_result * result)
{
  int  i;
  MpIeee sum_val;
  MpIeee sum_err;

  if(N <= 0) {
    MpIeee t_val=  MpIeee( "1.0" );
    MpIeee t_err=  MpIeee( "0.0" );
    gsl_sf_result poch;
    int  stat_poch;

    sum_val = MpIeee( "1.0" );
    sum_err = MpIeee( "0.0" );
    for(i=1; i<= -N; i++) {
      const MpIeee xi1=  i - 1;
      const MpIeee mult=  (a+xi1)*x/((b+xi1)*(xi1+1.0));
      t_val *= mult;
      t_err += fabs(mult) * t_err + fabs(t_val) * MpIeee( "8.0" ) * MpIeee( "2.0" ) * GSL_DBL_EPSILON;
      sum_val += t_val;
      sum_err += t_err;
    }

    stat_poch = gsl_sf_poch_e(1.0+a-b, -a, &poch);

    result->val  = sum_val * poch.val;
    result->err  = fabs(sum_val) * poch.err + sum_err * fabs(poch.val);
    result->err += fabs(poch.val) * (fabs(N) + 2.0) * GSL_DBL_EPSILON * fabs(sum_val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    result->err *= 2.0; /* FIXME: fudge factor... why is the error estimate too small? */
    return stat_poch;
  }
  else {
    const int M = N - 2;
    if(M < 0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result gbm1;
      gsl_sf_result gamr;
      int  stat_gbm1;
      int  stat_gamr;
      MpIeee t_val=  MpIeee( "1.0" );
      MpIeee t_err=  MpIeee( "0.0" );

      sum_val = MpIeee( "1.0" );
      sum_err = MpIeee( "0.0" );
      for(i=1; i<=M; i++) {
        const MpIeee mult=  (a-b+i)*x/((1.0-b+i)*i);
        t_val *= mult;
        t_err += t_err * fabs(mult) + fabs(t_val) * MpIeee( "8.0" ) * MpIeee( "2.0" ) * GSL_DBL_EPSILON;
        sum_val += t_val;
        sum_err += t_err;
      }

      stat_gbm1 = gsl_sf_gamma_e(b-1.0, &gbm1);
      stat_gamr = gsl_sf_gammainv_e(a,  &gamr);

      if(stat_gbm1 == GSL_SUCCESS) {
        gsl_sf_result powx1N;
        int  stat_p=  gsl_sf_pow_int_e(x, 1-N, &powx1N);
        MpIeee pe_val=  powx1N.val * xeps;
        MpIeee pe_err=  powx1N.err * fabs(xeps) + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(pe_val);
        MpIeee coeff_val=  gbm1.val * gamr.val * pe_val;
        MpIeee coeff_err=  gbm1.err * fabs(gamr.val * pe_val)
                         + gamr.err * fabs(gbm1.val * pe_val)
                         + fabs(gbm1.val * gamr.val) * pe_err
                         + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(coeff_val);

        result->val  = sum_val * coeff_val;
        result->err  = fabs(sum_val) * coeff_err + sum_err * fabs(coeff_val);
        result->err += 2.0 * GSL_DBL_EPSILON * (M+2.0) * fabs(result->val);
        result->err *= 2.0; /* FIXME: fudge factor... why is the error estimate too small? */
        return stat_p;
      }
      else {
        result->val = 0.0;
        result->err = 0.0;
        return stat_gbm1;
      }
    }
  }
}


/* Based on SLATEC DCHU() [W. Fullerton]
 * Assumes x > 0.
 * This is just a series summation method, and
 * it is not good for large a.
 *
 * I patched up the window for 1+a-b near zero. [GJ]
 */
static
int
 hyperg_U_series(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  const MpIeee EPS=  2.0 * GSL_DBL_EPSILON;  /* EPS = D1MACH(3) */
  const MpIeee SQRT_EPS=  M_SQRT2 * GSL_SQRT_DBL_EPSILON;

  if(fabs(1.0 + a - b) < SQRT_EPS) {
    /* Original Comment: ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X
     */
    /* We can however do the following:
     * U(a,b,x) = U(a,a+1,x) when 1+a-b=0
     * and U(a,a+1,x) = x^(-a).
     */
    MpIeee lnr=  -a * log(x);
    int  stat_e=   gsl_sf_exp_e(lnr, result);
    result->err += 2.0 * SQRT_EPS * fabs(result->val);
    return stat_e;
  }
  else {
    MpIeee aintb=  ( b < MpIeee( "0.0" ) ? ceil(b-MpIeee( "0.5" )) : floor(b+MpIeee( "0.5" )) );
    MpIeee beps=  b - aintb;
    int  N=  aintb;
    
    MpIeee lnx=  log(x);
    MpIeee xeps=  exp(-beps*lnx);

    /* Evaluate finite sum.
     */
    gsl_sf_result sum;
    int  stat_sum=  hyperg_U_finite_sum(N, a, b, x, xeps, &sum);


    /* Evaluate infinite sum. */

    int  istrt=  ( N < 1 ? 1-N : 0 );
    MpIeee xi=  istrt;

    gsl_sf_result gamr;
    gsl_sf_result powx;
    int  stat_gamr=  gsl_sf_gammainv_e(1.0+a-b, &gamr);
    int  stat_powx=  gsl_sf_pow_int_e(x, istrt, &powx);
    MpIeee sarg=  beps*M_PI;
    MpIeee sfact=  ( sarg != MpIeee( "0.0" ) ? sarg/sin(sarg) : MpIeee( "1.0" ) );
    MpIeee factor_val=  sfact * ( GSL_IS_ODD(N) ? -MpIeee( "1.0" ) : MpIeee( "1.0" ) ) * gamr.val * powx.val;
    MpIeee factor_err=  fabs(gamr.val) * powx.err + fabs(powx.val) * gamr.err
                      + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(factor_val);

    gsl_sf_result pochai;
    gsl_sf_result gamri1;
    gsl_sf_result gamrni;
    int  stat_pochai=  gsl_sf_poch_e(a, xi, &pochai);
    int  stat_gamri1=  gsl_sf_gammainv_e(xi + 1.0, &gamri1);
    int  stat_gamrni=  gsl_sf_gammainv_e(aintb + xi, &gamrni);
    int  stat_gam123=  GSL_ERROR_SELECT_3(stat_gamr, stat_gamri1, stat_gamrni);
    int  stat_gamall=  GSL_ERROR_SELECT_4(stat_sum, stat_gam123, stat_pochai, stat_powx);

    gsl_sf_result pochaxibeps;
    gsl_sf_result gamrxi1beps;
    int  stat_pochaxibeps=  gsl_sf_poch_e(a, xi-beps, &pochaxibeps);
    int  stat_gamrxi1beps=  gsl_sf_gammainv_e(xi + 1.0 - beps, &gamrxi1beps);

    int  stat_all=  GSL_ERROR_SELECT_3(stat_gamall, stat_pochaxibeps, stat_gamrxi1beps);

    MpIeee b0_val=  factor_val * pochaxibeps.val * gamrni.val * gamrxi1beps.val;
    MpIeee b0_err=   fabs(factor_val * pochaxibeps.val * gamrni.val) * gamrxi1beps.err
                   + fabs(factor_val * pochaxibeps.val * gamrxi1beps.val) * gamrni.err
                   + fabs(factor_val * gamrni.val * gamrxi1beps.val) * pochaxibeps.err
                   + fabs(pochaxibeps.val * gamrni.val * gamrxi1beps.val) * factor_err
                   + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(b0_val);

    if(fabs(xeps-1.0) < 0.5) {
      /*
       C  X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE
       C  CAREFUL IN EVALUATING THE DIFFERENCES.
       */
      int  i;
      gsl_sf_result pch1ai;
      gsl_sf_result pch1i;
      gsl_sf_result poch1bxibeps;
      int  stat_pch1ai=  gsl_sf_pochrel_e(a + xi, -beps, &pch1ai);
      int  stat_pch1i=  gsl_sf_pochrel_e(xi + 1.0 - beps, beps, &pch1i);
      int  stat_poch1bxibeps=  gsl_sf_pochrel_e(b+xi, -beps, &poch1bxibeps);
      MpIeee c0_t1_val=  beps*pch1ai.val*pch1i.val;
      MpIeee c0_t1_err=  fabs(beps) * fabs(pch1ai.val) * pch1i.err
                       + fabs(beps) * fabs(pch1i.val)  * pch1ai.err
                       + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(c0_t1_val);
      MpIeee c0_t2_val=  -poch1bxibeps.val + pch1ai.val - pch1i.val + c0_t1_val;
      MpIeee c0_t2_err=   poch1bxibeps.err + pch1ai.err + pch1i.err + c0_t1_err
                       + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(c0_t2_val);
      MpIeee c0_val=  factor_val * pochai.val * gamrni.val * gamri1.val * c0_t2_val;
      MpIeee c0_err=   fabs(factor_val * pochai.val * gamrni.val * gamri1.val) * c0_t2_err
                     + fabs(factor_val * pochai.val * gamrni.val * c0_t2_val) * gamri1.err
                     + fabs(factor_val * pochai.val * gamri1.val * c0_t2_val) * gamrni.err
                     + fabs(factor_val * gamrni.val * gamri1.val * c0_t2_val) * pochai.err
                     + fabs(pochai.val * gamrni.val * gamri1.val * c0_t2_val) * factor_err
                     + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(c0_val);
      /*
       C  XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
       */
      gsl_sf_result dexprl;
      int  stat_dexprl=  gsl_sf_exprel_e(-beps*lnx, &dexprl);
      MpIeee xeps1_val=  lnx * dexprl.val;
      MpIeee xeps1_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * (MpIeee( "1.0" ) + fabs(beps*lnx)) * fabs(dexprl.val)
                       + fabs(lnx) * dexprl.err
                       + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(xeps1_val);
      MpIeee dchu_val=  sum.val + c0_val + xeps1_val*b0_val;
      MpIeee dchu_err=  sum.err + c0_err
                      + fabs(xeps1_val)*b0_err + xeps1_err * fabs(b0_val)
                      + fabs(b0_val*lnx)*dexprl.err
                      + MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(sum.val) + fabs(c0_val) + fabs(xeps1_val*b0_val));
      MpIeee xn=  N;
      MpIeee t_val;
      MpIeee t_err;

      stat_all = GSL_ERROR_SELECT_5(stat_all, stat_dexprl, stat_poch1bxibeps, stat_pch1i, stat_pch1ai);

      for(i=1; i<2000; i++) {
        const MpIeee xi=  istrt + i;
        const MpIeee xi1=  istrt + i - 1;
        const MpIeee tmp=  (a-1.0)*(xn+2.0*xi-1.0) + xi*(xi-beps);
        const MpIeee b0_multiplier=  (a+xi1-beps)*x/((xn+xi1)*(xi-beps));
        const MpIeee c0_multiplier_1=  (a+xi1)*x/((b+xi1)*xi);
        const MpIeee c0_multiplier_2=  tmp / (xi*(b+xi1)*(a+xi1-beps));
        b0_val *= b0_multiplier;
        b0_err += fabs(b0_multiplier) * b0_err + fabs(b0_val) * MpIeee( "8.0" ) * MpIeee( "2.0" ) * GSL_DBL_EPSILON;
        c0_val  = c0_multiplier_1 * c0_val - c0_multiplier_2 * b0_val;
        c0_err  =  fabs(c0_multiplier_1) * c0_err
                 + fabs(c0_multiplier_2) * b0_err
                 + fabs(c0_val) * MpIeee( "8.0" ) * MpIeee( "2.0" ) * GSL_DBL_EPSILON
                 + fabs(b0_val * c0_multiplier_2) * MpIeee( "16.0" ) * MpIeee( "2.0" ) * GSL_DBL_EPSILON;
        t_val  = c0_val + xeps1_val*b0_val;
        t_err  = c0_err + fabs(xeps1_val)*b0_err;
        t_err += fabs(b0_val*lnx) * dexprl.err;
        t_err += fabs(b0_val)*xeps1_err;
        dchu_val += t_val;
        dchu_err += t_err;
        if(fabs(t_val) < EPS*fabs(dchu_val)) break;
      }

      result->val  = dchu_val;
      result->err  = 2.0 * dchu_err;
      result->err += 2.0 * fabs(t_val);
      result->err += 4.0 * GSL_DBL_EPSILON * (i+2.0) * fabs(dchu_val);
      result->err *= 2.0; /* FIXME: fudge factor */

      if(i >= 2000) {
        GSL_ERROR ("error", GSL_EMAXITER);
      }
      else {
        return stat_all;
      }
    }
    else {
      /*
       C  X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE
       C  STRAIGHTFORWARD FORMULATION IS STABLE.
       */
      int  i;
      MpIeee dchu_val;
      MpIeee dchu_err;
      MpIeee t_val;
      MpIeee t_err;
      gsl_sf_result dgamrbxi;
      int  stat_dgamrbxi=  gsl_sf_gammainv_e(b+xi, &dgamrbxi);
      MpIeee a0_val=  factor_val * pochai.val * dgamrbxi.val * gamri1.val / beps;
      MpIeee a0_err=   fabs(factor_val * pochai.val * dgamrbxi.val / beps) * gamri1.err
                     + fabs(factor_val * pochai.val * gamri1.val / beps) * dgamrbxi.err
                     + fabs(factor_val * dgamrbxi.val * gamri1.val / beps) * pochai.err
                     + fabs(pochai.val * dgamrbxi.val * gamri1.val / beps) * factor_err
                     + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(a0_val);
      stat_all = GSL_ERROR_SELECT_2(stat_all, stat_dgamrbxi);

      b0_val = xeps * b0_val / beps;
      b0_err = fabs(xeps / beps) * b0_err + MpIeee( "4.0" ) * GSL_DBL_EPSILON * fabs(b0_val);
      dchu_val = sum.val + a0_val - b0_val;
      dchu_err = sum.err + a0_err + b0_err
        + MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(sum.val) + fabs(a0_val) + fabs(b0_val));

      for(i=1; i<2000; i++) {
        MpIeee xi=  istrt + i;
        MpIeee xi1=  istrt + i - MpIeee( "1" );
        MpIeee a0_multiplier=  (a+xi1)*x/((b+xi1)*xi);
        MpIeee b0_multiplier=  (a+xi1-beps)*x/((aintb+xi1)*(xi-beps));
        a0_val *= a0_multiplier;
        a0_err += fabs(a0_multiplier) * a0_err;
        b0_val *= b0_multiplier;
        b0_err += fabs(b0_multiplier) * b0_err;
        t_val = a0_val - b0_val;
        t_err = a0_err + b0_err;
        dchu_val += t_val;
        dchu_err += t_err;
        if(fabs(t_val) < EPS*fabs(dchu_val)) break;
      }

      result->val  = dchu_val;
      result->err  = 2.0 * dchu_err;
      result->err += 2.0 * fabs(t_val);
      result->err += 4.0 * GSL_DBL_EPSILON * (i+2.0) * fabs(dchu_val);
      result->err *= 2.0; /* FIXME: fudge factor */

      if(i >= 2000) {
        GSL_ERROR ("error", GSL_EMAXITER);
      }
      else {
        return stat_all;
      }
    }
  }
}


/* Assumes b > 0 and x > 0.
 */
static
int
 hyperg_U_small_ab(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  if(a == -1.0) {
    /* U(-1,c+1,x) = Laguerre[c,0,x] = -b + x
     */
    result->val  = -b + x;
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(b) + fabs(x));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(a == 0.0) {
    /* U(0,c+1,x) = Laguerre[c,0,x] = 1
     */
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ASYMP_EVAL_OK(a,b,x)) {
    MpIeee p=  pow(x, -a);
    gsl_sf_result asymp;
    int  stat_asymp=  hyperg_zaU_asymp(a, b, x, &asymp);
    result->val  = asymp.val * p;
    result->err  = asymp.err * p;
    result->err += fabs(asymp.val) * GSL_DBL_EPSILON * fabs(a) * p;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_asymp;
  }
  else {
    return hyperg_U_series(a, b, x, result);
  }
}


/* Assumes b > 0 and x > 0.
 */
static
int
 hyperg_U_small_a_bgt0(const MpIeee a, const MpIeee b, const MpIeee x,
                      gsl_sf_result * result,
                      MpIeee * ln_multiplier)
{
  if(a == 0.0) {
    result->val = 1.0;
    result->err = 1.0;
    *ln_multiplier = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else if(   (b > 5000.0 && x < 0.90 * fabs(b))
          || (b >  500.0 && x < 0.50 * fabs(b))
    ) {
    int  stat=  gsl_sf_hyperg_U_large_b_e(a, b, x, result, ln_multiplier);
    if(stat == GSL_EOVRFLW)
      return GSL_SUCCESS;
    else
      return stat;
  }
  else if(b > 15.0) {
    /* Recurse up from b near 1.
     */
    MpIeee eps=  b - floor(b);
    MpIeee b0=  MpIeee( "1.0" ) + eps;
    gsl_sf_result r_Ubm1;
    gsl_sf_result r_Ub;
    int  stat_0=  hyperg_U_small_ab(a, b0,     x, &r_Ubm1);
    int  stat_1=  hyperg_U_small_ab(a, b0+1.0, x, &r_Ub);
    MpIeee Ubm1=  r_Ubm1.val;
    MpIeee Ub=  r_Ub.val;
    MpIeee Ubp1;
    MpIeee bp;

    for(bp = b0+MpIeee( "1.0" ); bp<b-MpIeee( "0.1" ); bp += MpIeee( "1.0" )) {
      Ubp1 = ((MpIeee( "1.0" )+a-bp)*Ubm1 + (bp+x-MpIeee( "1.0" ))*Ub)/x;
      Ubm1 = Ub;
      Ub   = Ubp1;
    }
    result->val  = Ub;
    result->err  = (fabs(r_Ubm1.err/r_Ubm1.val) + fabs(r_Ub.err/r_Ub.val)) * fabs(Ub);
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(b-b0)+1.0) * fabs(Ub);
    *ln_multiplier = MpIeee( "0.0" );
    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
  else {
    *ln_multiplier = MpIeee( "0.0" );
    return hyperg_U_small_ab(a, b, x, result);
  }
}


/* We use this to keep track of large
 * dynamic ranges in the recursions.
 * This can be important because sometimes
 * we want to calculate a very large and
 * a very small number and the answer is
 * the product, of order 1. This happens,
 * for instance, when we apply a Kummer
 * transform to make b positive and
 * both x and b are large.
 */
#define RESCALE_2(u0,u1,factor,count)      \
do {                                       \
  MpIeee au0=  fabs(u0);                   \
  if(au0 > factor) {                       \
    u0 /= factor;                          \
    u1 /= factor;                          \
    count++;                               \
  }                                        \
  else if(au0 < MpIeee( "1.0" )/factor) {              \
    u0 *= factor;                          \
    u1 *= factor;                          \
    count--;                               \
  }                                        \
} while (0)


/* Specialization to b >= 1, for integer parameters.
 * Assumes x > 0.
 */
static
int
 hyperg_U_int_bge1(const int a, const int b, const MpIeee x,
                  gsl_sf_result_e10 * result)
{
  if(a == 0) {
    result->val = 1.0;
    result->err = 0.0;
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else if(a == -1) {
    result->val  = -b + x;
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(b) + fabs(x));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    result->e10  = 0;
    return GSL_SUCCESS;
  }
  else if(b == a + 1) {
    /* U(a,a+1,x) = x^(-a)
     */
    return gsl_sf_exp_e10_e(-a*log(x), result);
  }
  else if(ASYMP_EVAL_OK(a,b,x)) {
    const MpIeee ln_pre_val=  -a*log(x);
    const MpIeee ln_pre_err=  2.0 * GSL_DBL_EPSILON * fabs(ln_pre_val);
    gsl_sf_result asymp;
    int  stat_asymp=  hyperg_zaU_asymp(a, b, x, &asymp);
    int  stat_e=  gsl_sf_exp_mult_err_e10_e(ln_pre_val, ln_pre_err,
                                              asymp.val, asymp.err,
                                              result);
    return GSL_ERROR_SELECT_2(stat_e, stat_asymp);
  }
  else if(SERIES_EVAL_OK(a,b,x)) {
    gsl_sf_result ser;
    const int stat_ser = hyperg_U_series(a, b, x, &ser);
    result->val = ser.val;
    result->err = ser.err;
    result->e10 = 0;
    return stat_ser;
  }
  else if(a < 0) {
    /* Recurse backward from a = -1,0.
     */
    int  scale_count=  0;
    const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
    gsl_sf_result lnm;
    gsl_sf_result y;
    MpIeee lnscale;
    MpIeee Uap1=  MpIeee( "1.0" );     /* U(0,b,x)  */
    MpIeee Ua=  -b + x;  /* U(-1,b,x) */
    MpIeee Uam1;
    int  ap;

    for(ap=-1; ap>a; ap--) {
      Uam1 = ap*(b-ap-MpIeee( "1.0" ))*Uap1 + (x+MpIeee( "2.0" )*ap-b)*Ua;
      Uap1 = Ua;
      Ua   = Uam1;
      RESCALE_2(Ua,Uap1,scale_factor,scale_count);
    }

    lnscale = log(scale_factor);
    lnm.val = scale_count*lnscale;
    lnm.err = 2.0 * GSL_DBL_EPSILON * fabs(lnm.val);
    y.val = Ua;
    y.err = 4.0 * GSL_DBL_EPSILON * (fabs(a)+1.0) * fabs(Ua);
    return gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
  }
  else if(b >= 2.0*a + x) {
    /* Recurse forward from a = 0,1.
     */
    int  scale_count=  0;
    const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
    gsl_sf_result r_Ua;
    gsl_sf_result lnm;
    gsl_sf_result y;
    MpIeee lnscale;
    MpIeee lm;
    int  stat_1=  hyperg_U_small_a_bgt0(1.0, b, x, &r_Ua, &lm);  /* U(1,b,x) */
    int  stat_e;
    MpIeee Uam1=  MpIeee( "1.0" );  /* U(0,b,x) */
    MpIeee Ua=  r_Ua.val;
    MpIeee Uap1;
    int  ap;

    Uam1 *= exp(-lm);

    for(ap=1; ap<a; ap++) {
      Uap1 = -(Uam1 + (b-MpIeee( "2.0" )*ap-x)*Ua)/(ap*(MpIeee( "1.0" )+ap-b));
      Uam1 = Ua;
      Ua   = Uap1;
      RESCALE_2(Ua,Uam1,scale_factor,scale_count);
    }

    lnscale = log(scale_factor);
    lnm.val = lm + scale_count * lnscale;
    lnm.err = 2.0 * GSL_DBL_EPSILON * (fabs(lm) + fabs(scale_count*lnscale));
    y.val  = Ua;
    y.err  = fabs(r_Ua.err/r_Ua.val) * fabs(Ua);
    y.err += 2.0 * GSL_DBL_EPSILON * (fabs(a) + 1.0) * fabs(Ua);
    stat_e = gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
    return GSL_ERROR_SELECT_2(stat_e, stat_1);
  }
  else {
    if(b <= x) {
      /* Recurse backward either to the b=a+1 line
       * or to a=0, whichever we hit.
       */
      const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
      int  scale_count=  0;
      int  stat_CF1;
      MpIeee ru;
      int  CF1_count;
      int  a_target;
      MpIeee lnU_target;
      MpIeee Ua;
      MpIeee Uap1;
      MpIeee Uam1;
      int  ap;

      if(b < a + 1) {
        a_target = b-1;
        lnU_target = -a_target*log(x);
      }
      else {
        a_target = 0;
        lnU_target = MpIeee( "0.0" );
      }

      stat_CF1 = hyperg_U_CF1(a, b, 0, x, &ru, &CF1_count);

      Ua   = MpIeee( "1.0" );
      Uap1 = ru/a * Ua;
      for(ap=a; ap>a_target; ap--) {
        Uam1 = -((b-MpIeee( "2.0" )*ap-x)*Ua + ap*(MpIeee( "1.0" )+ap-b)*Uap1);
        Uap1 = Ua;
        Ua   = Uam1;
        RESCALE_2(Ua,Uap1,scale_factor,scale_count);
      }

      if(Ua == MpIeee( "0.0" )) {
        result->val = 0.0;
        result->err = 0.0;
        result->e10 = 0;
        GSL_ERROR ("error", GSL_EZERODIV);
      }
      else {
        MpIeee lnscl=  -scale_count*log(scale_factor);
        MpIeee lnpre_val=  lnU_target + lnscl;
        MpIeee lnpre_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(lnU_target) + fabs(lnscl));
        MpIeee oUa_err=  MpIeee( "2.0" ) * (fabs(a_target-a) + CF1_count + MpIeee( "1.0" )) * GSL_DBL_EPSILON * fabs(MpIeee( "1.0" )/Ua);
        int  stat_e=  gsl_sf_exp_mult_err_e10_e(lnpre_val, lnpre_err,
                                                  1.0/Ua, oUa_err,
                                                  result);
        return GSL_ERROR_SELECT_2(stat_e, stat_CF1);
      }
    }
    else {
      /* Recurse backward to near the b=2a+x line, then
       * determine normalization by either direct evaluation
       * or by a forward recursion. The direct evaluation
       * is needed when x is small (which is precisely
       * when it is easy to do).
       */
      const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
      int  scale_count_for=  0;
      int  scale_count_bck=  0;
      int  a0=  1;
      int  a1=  a0 + ceil(0.5*(b-x) - a0);
      MpIeee Ua1_bck_val;
      MpIeee Ua1_bck_err;
      MpIeee Ua1_for_val;
      MpIeee Ua1_for_err;
      int  stat_for;
      int  stat_bck;
      gsl_sf_result lm_for;

      {
        /* Recurse back to determine U(a1,b), sans normalization.
         */
        MpIeee ru;
        int  CF1_count;
        int  stat_CF1=  hyperg_U_CF1(a, b, 0, x, &ru, &CF1_count);
        MpIeee Ua=  MpIeee( "1.0" );
        MpIeee Uap1=  ru/a * Ua;
        MpIeee Uam1;
        int  ap;
        for(ap=a; ap>a1; ap--) {
          Uam1 = -((b-MpIeee( "2.0" )*ap-x)*Ua + ap*(MpIeee( "1.0" )+ap-b)*Uap1);
          Uap1 = Ua;
          Ua   = Uam1;
          RESCALE_2(Ua,Uap1,scale_factor,scale_count_bck);
        }
        Ua1_bck_val = Ua;
        Ua1_bck_err = MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(a1-a)+CF1_count+MpIeee( "1.0" )) * fabs(Ua);
        stat_bck = stat_CF1;
      }

      if(b == 2*a1 && a1 > 1) {
        /* This can happen when x is small, which is
         * precisely when we need to be careful with
         * this evaluation.
         */
        hyperg_lnU_beq2a((MpIeee)a1, x, &lm_for);
        Ua1_for_val = MpIeee( "1.0" );
        Ua1_for_err = MpIeee( "0.0" );
        stat_for = GSL_SUCCESS;
      }
      else if(b == 2*a1 - 1 && a1 > 1) {
        /* Similar to the above. Happens when x is small.
         * Use
         *   U(a,2a-1) = (x U(a,2a) - U(a-1,2(a-1))) / (2a - 2)
         */
        gsl_sf_result lnU00, lnU12;
        gsl_sf_result U00, U12;
        hyperg_lnU_beq2a(a1-1.0, x, &lnU00);
        hyperg_lnU_beq2a(a1,     x, &lnU12);
        if(lnU00.val > lnU12.val) {
          lm_for.val = lnU00.val;
          lm_for.err = lnU00.err;
          U00.val = 1.0;
          U00.err = 0.0;
          gsl_sf_exp_err_e(lnU12.val - lm_for.val, lnU12.err + lm_for.err, &U12);
        }
        else {
          lm_for.val = lnU12.val;
          lm_for.err = lnU12.err;
          U12.val = 1.0;
          U12.err = 0.0;
          gsl_sf_exp_err_e(lnU00.val - lm_for.val, lnU00.err + lm_for.err, &U00);
        }
        Ua1_for_val  = (x * U12.val - U00.val) / (MpIeee( "2.0" )*a1 - MpIeee( "2.0" ));
        Ua1_for_err  = (fabs(x)*U12.err + U00.err) / fabs(MpIeee( "2.0" )*a1 - MpIeee( "2.0" ));
        Ua1_for_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(Ua1_for_val);
        stat_for = GSL_SUCCESS;
      }
      else {
        /* Recurse forward to determine U(a1,b) with
         * absolute normalization.
         */
        gsl_sf_result r_Ua;
        MpIeee Uam1=  MpIeee( "1.0" );  /* U(a0-1,b,x) = U(0,b,x) */
        MpIeee Ua;
        MpIeee Uap1;
        int  ap;
        MpIeee lm_for_local;
        stat_for = hyperg_U_small_a_bgt0(a0, b, x, &r_Ua, &lm_for_local); /* U(1,b,x) */
        Ua = r_Ua.val;
        Uam1 *= exp(-lm_for_local);
        lm_for.val = lm_for_local;
        lm_for.err = 0.0;

        for(ap=a0; ap<a1; ap++) {
          Uap1 = -(Uam1 + (b-MpIeee( "2.0" )*ap-x)*Ua)/(ap*(MpIeee( "1.0" )+ap-b));
          Uam1 = Ua;
          Ua   = Uap1;
          RESCALE_2(Ua,Uam1,scale_factor,scale_count_for);
        }
        Ua1_for_val  = Ua;
        Ua1_for_err  = fabs(Ua) * fabs(r_Ua.err/r_Ua.val);
        Ua1_for_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(a1-a0)+MpIeee( "1.0" )) * fabs(Ua1_for_val);
      }

      /* Now do the matching to produce the final result.
       */
      if(Ua1_bck_val == MpIeee( "0.0" )) {
        result->val = 0.0;
        result->err = 0.0;
        result->e10 = 0;
        GSL_ERROR ("error", GSL_EZERODIV);
      }
      else if(Ua1_for_val == MpIeee( "0.0" )) {
        /* Should never happen. */
        UNDERFLOW_ERROR_E10(result);
      }
      else {
        MpIeee lns=  (scale_count_for - scale_count_bck)*log(scale_factor);
        MpIeee ln_for_val=  log(fabs(Ua1_for_val));
        MpIeee ln_for_err=  GSL_DBL_EPSILON + fabs(Ua1_for_err/Ua1_for_val);
        MpIeee ln_bck_val=  log(fabs(Ua1_bck_val));
        MpIeee ln_bck_err=  GSL_DBL_EPSILON + fabs(Ua1_bck_err/Ua1_bck_val);
        MpIeee lnr_val=  lm_for.val + ln_for_val - ln_bck_val + lns;
        MpIeee lnr_err=  lm_for.err + ln_for_err + ln_bck_err
          + MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(lm_for.val) + fabs(ln_for_val) + fabs(ln_bck_val) + fabs(lns));
        MpIeee sgn=  GSL_SIGN(Ua1_for_val) * GSL_SIGN(Ua1_bck_val);
        int  stat_e=  gsl_sf_exp_err_e10_e(lnr_val, lnr_err, result);
        result->val *= sgn;
        return GSL_ERROR_SELECT_3(stat_e, stat_bck, stat_for);
      }
    }
  }
}


/* Handle b >= 1 for generic a,b values.
 */
static
int
 hyperg_U_bge1(const MpIeee a, const MpIeee b, const MpIeee x,
              gsl_sf_result_e10 * result)
{
  const MpIeee rinta=  floor(a+0.5);
  const int a_neg_integer = (a < 0.0 && fabs(a - rinta) < INT_THRESHOLD);

  if(a == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else if(a_neg_integer && fabs(rinta) < INT_MAX) {
    /* U(-n,b,x) = (-1)^n n! Laguerre[n,b-1,x]
     */
    const int n = -(int)rinta;
    const MpIeee sgn=  (GSL_IS_ODD(n) ? -1.0 : 1.0);
    gsl_sf_result lnfact;
    gsl_sf_result L;
    const int stat_L = gsl_sf_laguerre_n_e(n, b-1.0, x, &L);
    gsl_sf_lnfact_e(n, &lnfact);
    {
      const int stat_e = gsl_sf_exp_mult_err_e10_e(lnfact.val, lnfact.err,
                                                      sgn*L.val, L.err,
                                                      result);
      return GSL_ERROR_SELECT_2(stat_e, stat_L);
    }
  }
  else if(ASYMP_EVAL_OK(a,b,x)) {
    const MpIeee ln_pre_val=  -a*log(x);
    const MpIeee ln_pre_err=  2.0 * GSL_DBL_EPSILON * fabs(ln_pre_val);
    gsl_sf_result asymp;
    int  stat_asymp=  hyperg_zaU_asymp(a, b, x, &asymp);
    int  stat_e=  gsl_sf_exp_mult_err_e10_e(ln_pre_val, ln_pre_err,
                                              asymp.val, asymp.err,
                                              result);
    return GSL_ERROR_SELECT_2(stat_e, stat_asymp);
  }
  else if(fabs(a) <= 1.0) {
    gsl_sf_result rU;
    MpIeee ln_multiplier;
    int  stat_U=  hyperg_U_small_a_bgt0(a, b, x, &rU, &ln_multiplier);
    int  stat_e=  gsl_sf_exp_mult_err_e10_e(ln_multiplier, 2.0*GSL_DBL_EPSILON*fabs(ln_multiplier), rU.val, rU.err, result);
    return GSL_ERROR_SELECT_2(stat_U, stat_e);
  }
  else if(SERIES_EVAL_OK(a,b,x)) {
    gsl_sf_result ser;
    const int stat_ser = hyperg_U_series(a, b, x, &ser);
    result->val = ser.val;
    result->err = ser.err;
    result->e10 = 0;
    return stat_ser;
  }
  else if(a < 0.0) {
    /* Recurse backward on a and then upward on b.
     */
    const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
    const MpIeee a0=  a - floor(a) - 1.0;
    const MpIeee b0=  b - floor(b) + 1.0;
    int  scale_count=  0;
    MpIeee lm_0;MpIeee  lm_1;
    MpIeee lm_max;
    gsl_sf_result r_Uap1;
    gsl_sf_result r_Ua;
    int  stat_0=  hyperg_U_small_a_bgt0(a0+1.0, b0, x, &r_Uap1, &lm_0);
    int  stat_1=  hyperg_U_small_a_bgt0(a0,     b0, x, &r_Ua,   &lm_1);
    int  stat_e;
    MpIeee Uap1=  r_Uap1.val;
    MpIeee Ua=  r_Ua.val;
    MpIeee Uam1;
    MpIeee ap;
    lm_max = GSL_MAX(lm_0, lm_1);
    Uap1 *= exp(lm_0-lm_max);
    Ua   *= exp(lm_1-lm_max);

    /* Downward recursion on a.
     */
    for(ap=a0; ap>a+MpIeee( "0.1" ); ap -= MpIeee( "1.0" )) {
      Uam1 = ap*(b0-ap-MpIeee( "1.0" ))*Uap1 + (x+MpIeee( "2.0" )*ap-b0)*Ua;
      Uap1 = Ua;
      Ua   = Uam1;
      RESCALE_2(Ua,Uap1,scale_factor,scale_count);
    }

    if(b < 2.0) {
      /* b == b0, so no recursion necessary
       */
      const MpIeee lnscale=  log(scale_factor);
      gsl_sf_result lnm;
      gsl_sf_result y;
      lnm.val = lm_max + scale_count * lnscale;
      lnm.err = 2.0 * GSL_DBL_EPSILON * (fabs(lm_max) + scale_count * fabs(lnscale));
      y.val  = Ua;
      y.err  = fabs(r_Uap1.err/r_Uap1.val) * fabs(Ua);
      y.err += fabs(r_Ua.err/r_Ua.val) * fabs(Ua);
      y.err += 2.0 * GSL_DBL_EPSILON * (fabs(a-a0) + 1.0) * fabs(Ua);
      y.err *= fabs(lm_0-lm_max) + 1.0;
      y.err *= fabs(lm_1-lm_max) + 1.0;
      stat_e = gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
    }
    else {
      /* Upward recursion on b.
       */
      const MpIeee err_mult=  fabs(b-b0) + fabs(a-a0) + 1.0;
      const MpIeee lnscale=  log(scale_factor);
      gsl_sf_result lnm;
      gsl_sf_result y;

      MpIeee Ubm1=  Ua;                                 /* U(a,b0)   */
      MpIeee Ub=  (a*(b0-a-MpIeee( "1.0" ))*Uap1 + (a+x)*Ua)/x;   /* U(a,b0+1) */
      MpIeee Ubp1;
      MpIeee bp;
      for(bp=b0+MpIeee( "1.0" ); bp<b-MpIeee( "0.1" ); bp += MpIeee( "1.0" )) {
        Ubp1 = ((MpIeee( "1.0" )+a-bp)*Ubm1 + (bp+x-MpIeee( "1.0" ))*Ub)/x;
        Ubm1 = Ub;
        Ub   = Ubp1;
        RESCALE_2(Ub,Ubm1,scale_factor,scale_count);
      }

      lnm.val = lm_max + scale_count * lnscale;
      lnm.err = 2.0 * GSL_DBL_EPSILON * (fabs(lm_max) + fabs(scale_count * lnscale));
      y.val = Ub;
      y.err  = 2.0 * err_mult * fabs(r_Uap1.err/r_Uap1.val) * fabs(Ub);
      y.err += 2.0 * err_mult * fabs(r_Ua.err/r_Ua.val) * fabs(Ub);
      y.err += 2.0 * GSL_DBL_EPSILON * err_mult * fabs(Ub);
      y.err *= fabs(lm_0-lm_max) + 1.0;
      y.err *= fabs(lm_1-lm_max) + 1.0;
      stat_e = gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
    }
    return GSL_ERROR_SELECT_3(stat_e, stat_0, stat_1);
  }
  else if(b >= 2*a + x) {
    /* Recurse forward from a near zero.
     * Note that we cannot cross the singularity at
     * the line b=a+1, because the only way we could
     * be in that little wedge is if a < 1. But we
     * have already dealt with the small a case.
     */
    int  scale_count=  0;
    const MpIeee a0=  a - floor(a);
    const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
    MpIeee lnscale;
    MpIeee lm_0;MpIeee  lm_1;MpIeee  lm_max;
    gsl_sf_result r_Uam1;
    gsl_sf_result r_Ua;
    int  stat_0=  hyperg_U_small_a_bgt0(a0-1.0, b, x, &r_Uam1, &lm_0);
    int  stat_1=  hyperg_U_small_a_bgt0(a0,     b, x, &r_Ua,   &lm_1);
    int  stat_e;
    gsl_sf_result lnm;
    gsl_sf_result y;
    MpIeee Uam1=  r_Uam1.val;
    MpIeee Ua=  r_Ua.val;
    MpIeee Uap1;
    MpIeee ap;
    lm_max = GSL_MAX(lm_0, lm_1);
    Uam1 *= exp(lm_0-lm_max);
    Ua   *= exp(lm_1-lm_max);

    for(ap=a0; ap<a-MpIeee( "0.1" ); ap += MpIeee( "1.0" )) {
      Uap1 = -(Uam1 + (b-MpIeee( "2.0" )*ap-x)*Ua)/(ap*(MpIeee( "1.0" )+ap-b));
      Uam1 = Ua;
      Ua   = Uap1;
      RESCALE_2(Ua,Uam1,scale_factor,scale_count);
    }

    lnscale = log(scale_factor);
    lnm.val = lm_max + scale_count * lnscale;
    lnm.err = 2.0 * GSL_DBL_EPSILON * (fabs(lm_max) + fabs(scale_count * lnscale));
    y.val  = Ua;
    y.err  = fabs(r_Uam1.err/r_Uam1.val) * fabs(Ua);
    y.err += fabs(r_Ua.err/r_Ua.val) * fabs(Ua);
    y.err += 2.0 * GSL_DBL_EPSILON * (fabs(a-a0) + 1.0) * fabs(Ua);
    y.err *= fabs(lm_0-lm_max) + 1.0;
    y.err *= fabs(lm_1-lm_max) + 1.0;
    stat_e = gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
    return GSL_ERROR_SELECT_3(stat_e, stat_0, stat_1);
  }
  else {
    if(b <= x) {
      /* Recurse backward to a near zero.
       */
      const MpIeee a0=  a - floor(a);
      const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
      int  scale_count=  0;
      gsl_sf_result lnm;
      gsl_sf_result y;
      MpIeee lnscale;
      MpIeee lm_0;
      MpIeee Uap1;
      MpIeee Ua;
      MpIeee Uam1;
      gsl_sf_result U0;
      MpIeee ap;
      MpIeee ru;
      MpIeee r;
      int  CF1_count;
      int  stat_CF1=  hyperg_U_CF1(a, b, 0, x, &ru, &CF1_count);
      int  stat_U0;
      int  stat_e;
      r = ru/a;
      Ua   = GSL_SQRT_DBL_MIN;
      Uap1 = r * Ua;
      for(ap=a; ap>a0+MpIeee( "0.1" ); ap -= MpIeee( "1.0" )) {
        Uam1 = -((b-MpIeee( "2.0" )*ap-x)*Ua + ap*(MpIeee( "1.0" )+ap-b)*Uap1);
        Uap1 = Ua;
        Ua   = Uam1;
        RESCALE_2(Ua,Uap1,scale_factor,scale_count);
      }

      stat_U0 = hyperg_U_small_a_bgt0(a0, b, x, &U0, &lm_0);

      lnscale = log(scale_factor);
      lnm.val = lm_0 - scale_count * lnscale;
      lnm.err = 2.0 * GSL_DBL_EPSILON * (fabs(lm_0) + fabs(scale_count * lnscale));
      y.val  = GSL_SQRT_DBL_MIN*(U0.val/Ua);
      y.err  = GSL_SQRT_DBL_MIN*(U0.err/fabs(Ua));
      y.err += 2.0 * GSL_DBL_EPSILON * (fabs(a0-a) + CF1_count + 1.0) * fabs(y.val);
      stat_e = gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
      return GSL_ERROR_SELECT_3(stat_e, stat_U0, stat_CF1);
    }
    else {
      /* Recurse backward to near the b=2a+x line, then
       * forward from a near zero to get the normalization.
       */
      int  scale_count_for=  0;
      int  scale_count_bck=  0;
      const MpIeee scale_factor=  GSL_SQRT_DBL_MAX;
      const MpIeee eps=  a - floor(a);
      const MpIeee a0=  ( eps == 0.0 ? 1.0 : eps );
      const MpIeee a1=  a0 + ceil(0.5*(b-x) - a0);
      gsl_sf_result lnm;
      gsl_sf_result y;
      MpIeee lm_for;
      MpIeee lnscale;
      MpIeee Ua1_bck;
      MpIeee Ua1_for;
      int  stat_for;
      int  stat_bck;
      int  stat_e;
      int  CF1_count;

      {
        /* Recurse back to determine U(a1,b), sans normalization.
         */
        MpIeee Uap1;
        MpIeee Ua;
        MpIeee Uam1;
        MpIeee ap;
        MpIeee ru;
        MpIeee r;
        int  stat_CF1=  hyperg_U_CF1(a, b, 0, x, &ru, &CF1_count);
        r = ru/a;
        Ua   = GSL_SQRT_DBL_MIN;
        Uap1 = r * Ua;
        for(ap=a; ap>a1+MpIeee( "0.1" ); ap -= MpIeee( "1.0" )) {
          Uam1 = -((b-MpIeee( "2.0" )*ap-x)*Ua + ap*(MpIeee( "1.0" )+ap-b)*Uap1);
          Uap1 = Ua;
          Ua   = Uam1;
          RESCALE_2(Ua,Uap1,scale_factor,scale_count_bck);
        }
        Ua1_bck = Ua;
        stat_bck = stat_CF1;
      }
      {
        /* Recurse forward to determine U(a1,b) with
         * absolute normalization.
         */
        gsl_sf_result r_Uam1;
        gsl_sf_result r_Ua;
        MpIeee lm_0;MpIeee  lm_1;
        int  stat_0=  hyperg_U_small_a_bgt0(a0-1.0, b, x, &r_Uam1, &lm_0);
        int  stat_1=  hyperg_U_small_a_bgt0(a0,     b, x, &r_Ua,   &lm_1);
        MpIeee Uam1=  r_Uam1.val;
        MpIeee Ua=  r_Ua.val;
        MpIeee Uap1;
        MpIeee ap;

        lm_for = GSL_MAX(lm_0, lm_1);
        Uam1 *= exp(lm_0 - lm_for);
        Ua   *= exp(lm_1 - lm_for);

        for(ap=a0; ap<a1-MpIeee( "0.1" ); ap += MpIeee( "1.0" )) {
          Uap1 = -(Uam1 + (b-MpIeee( "2.0" )*ap-x)*Ua)/(ap*(MpIeee( "1.0" )+ap-b));
          Uam1 = Ua;
          Ua   = Uap1;
          RESCALE_2(Ua,Uam1,scale_factor,scale_count_for);
        }
        Ua1_for = Ua;
        stat_for = GSL_ERROR_SELECT_2(stat_0, stat_1);
      }

      lnscale = log(scale_factor);
      lnm.val = lm_for + (scale_count_for - scale_count_bck)*lnscale;
      lnm.err = 2.0 * GSL_DBL_EPSILON * (fabs(lm_for) + fabs(scale_count_for - scale_count_bck)*fabs(lnscale));
      y.val = GSL_SQRT_DBL_MIN*Ua1_for/Ua1_bck;
      y.err = 2.0 * GSL_DBL_EPSILON * (fabs(a-a0) + CF1_count + 1.0) * fabs(y.val);
      stat_e = gsl_sf_exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result);
      return GSL_ERROR_SELECT_3(stat_e, stat_bck, stat_for);
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int
 gsl_sf_hyperg_U_int_e10_e(const int a, const int b, const MpIeee x,
                             gsl_sf_result_e10 * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR_E10(result);
  }
  else {
    if(b >= 1) {
      return hyperg_U_int_bge1(a, b, x, result);
    }
    else {
      /* Use the reflection formula
       * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
       */
      gsl_sf_result_e10 U;
      MpIeee ln_x=  log(x);
      int  ap=  1 + a - b;
      int  bp=  2 - b;
      int  stat_e;
      int  stat_U=  hyperg_U_int_bge1(ap, bp, x, &U);
      MpIeee ln_pre_val=  (MpIeee( "1.0" )-b)*ln_x;
      MpIeee ln_pre_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * (fabs(b)+MpIeee( "1.0" )) * fabs(ln_x);
      ln_pre_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(MpIeee( "1.0" )-b); /* error in log(x) */
      stat_e = gsl_sf_exp_mult_err_e10_e(ln_pre_val + U.e10*M_LN10, ln_pre_err,
                                            U.val, U.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_U);
    }
  }
}


int
 gsl_sf_hyperg_U_e10_e(const MpIeee a, const MpIeee b, const MpIeee x,
                         gsl_sf_result_e10 * result)
{
  const MpIeee rinta=  floor(a + 0.5);
  const MpIeee rintb=  floor(b + 0.5);
  const int a_integer = ( fabs(a - rinta) < INT_THRESHOLD );
  const int b_integer = ( fabs(b - rintb) < INT_THRESHOLD );

  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR_E10(result);
  }
  else if(a == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else if(a_integer && b_integer) {
    return gsl_sf_hyperg_U_int_e10_e(rinta, rintb, x, result);
  }
  else {
    if(b >= 1.0) {
      /* Use b >= 1 function.
       */
      return hyperg_U_bge1(a, b, x, result);
    }
    else {
      /* Use the reflection formula
       * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
       */
      const MpIeee lnx=  log(x);
      const MpIeee ln_pre_val=  (1.0-b)*lnx;
      const MpIeee ln_pre_err=  fabs(lnx) * 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(b));
      const MpIeee ap=  1.0 + a - b;
      const MpIeee bp=  2.0 - b;
      gsl_sf_result_e10 U;
      int  stat_U=  hyperg_U_bge1(ap, bp, x, &U);
      int  stat_e=  gsl_sf_exp_mult_err_e10_e(ln_pre_val + U.e10*M_LN10, ln_pre_err,
                                            U.val, U.err,
                                            result);
      return GSL_ERROR_SELECT_2(stat_e, stat_U);
    }
  }
}


int
 gsl_sf_hyperg_U_int_e(const int a, const int b, const MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result_e10 re;
  int  stat_U=  gsl_sf_hyperg_U_int_e10_e(a, b, x, &re);
  int  stat_c=  gsl_sf_result_smash_e(&re, result);
  return GSL_ERROR_SELECT_2(stat_c, stat_U);
}


int
 gsl_sf_hyperg_U_e(const MpIeee a, const MpIeee b, const MpIeee x, gsl_sf_result * result)
{
  gsl_sf_result_e10 re;
  int  stat_U=  gsl_sf_hyperg_U_e10_e(a, b, x, &re);
  int  stat_c=  gsl_sf_result_smash_e(&re, result);
  return GSL_ERROR_SELECT_2(stat_c, stat_U);
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_hyperg_U_int(const int a, const int b, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_hyperg_U_int_e(a, b, x, &result));
}

MpIeee gsl_sf_hyperg_U(const MpIeee a, const MpIeee b, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_hyperg_U_e(a, b, x, &result));
}
