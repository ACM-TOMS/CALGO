#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/legendre_H3d.c
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
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_legendre.h>

#include "error.h"

#include "legendre.h"

/* See [Abbott+Schaefer, Ap.J. 308, 546 (1986)] for
 * enough details to follow what is happening here.
 */


/* Logarithm of normalization factor, Log[N(ell,lambda)].
 * N(ell,lambda) = Product[ lambda^2 + n^2, {n,0,ell} ]
 *               = |Gamma(ell + 1 + I lambda)|^2  lambda sinh(Pi lambda) / Pi
 * Assumes ell >= 0.
 */
static
int
 legendre_H3d_lnnorm(const int ell, const MpIeee lambda, MpIeee * result)
{
  MpIeee abs_lam=  fabs(lambda);

  if(abs_lam == MpIeee( "0.0" )) {
    *result = MpIeee( "0.0" );
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(lambda > (ell + 1.0)/GSL_ROOT3_DBL_EPSILON) {
    /* There is a cancellation between the sinh(Pi lambda)
     * term and the log(gamma(ell + 1 + i lambda) in the
     * result below, so we show some care and save some digits.
     * Note that the above guarantees that lambda is large,
     * since ell >= 0. We use Stirling and a simple expansion
     * of sinh.
     */
    MpIeee rat=  (ell+MpIeee( "1.0" ))/lambda;
    MpIeee ln_lam2ell2=  MpIeee( "2.0" )*log(lambda) + log(MpIeee( "1.0" ) + rat*rat);
    MpIeee lg_corrected=  -MpIeee( "2.0" )*(ell+MpIeee( "1.0" )) + M_LNPI + (ell+MpIeee( "0.5" ))*ln_lam2ell2 + MpIeee( "1.0" )/(MpIeee( "288.0" )*lambda*lambda);
    MpIeee angle_terms=  lambda * MpIeee( "2.0" ) * rat * (MpIeee( "1.0" ) - rat*rat/MpIeee( "3.0" ));
    *result = log(abs_lam) + lg_corrected + angle_terms - M_LNPI;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result lg_r;
    gsl_sf_result lg_theta;
    gsl_sf_result ln_sinh;
    gsl_sf_lngamma_complex_e(ell+1.0, lambda, &lg_r, &lg_theta);
    gsl_sf_lnsinh_e(M_PI * abs_lam, &ln_sinh);
    *result = log(abs_lam) + ln_sinh.val + MpIeee( "2.0" )*lg_r.val - M_LNPI;
    return GSL_SUCCESS;
  }
}


/* Calculate series for small eta*lambda.
 * Assumes eta > 0, lambda != 0.
 *
 * This is just the defining hypergeometric for the Legendre function.
 *
 * P^{mu}_{-1/2 + I lam}(z) = 1/Gamma(l+3/2) ((z+1)/(z-1)^(mu/2)
 *                            2F1(1/2 - I lam, 1/2 + I lam; l+3/2; (1-z)/2)
 * We use
 *       z = cosh(eta)
 * (z-1)/2 = sinh^2(eta/2)
 *
 * And recall
 * H3d = sqrt(Pi Norm /(2 lam^2 sinh(eta))) P^{-l-1/2}_{-1/2 + I lam}(cosh(eta))
 */
static
int
 legendre_H3d_series(const int ell, const MpIeee lambda, const MpIeee eta,
                    gsl_sf_result * result)
{
  const int nmax = 5000;
  const MpIeee shheta=  sinh(0.5*eta);
  const MpIeee ln_zp1=  M_LN2 + log(1.0 + shheta*shheta);
  const MpIeee ln_zm1=  M_LN2 + 2.0*log(shheta);
  const MpIeee zeta=  -shheta*shheta;
  gsl_sf_result lg_lp32;
  MpIeee term=  MpIeee( "1.0" );
  MpIeee sum=  MpIeee( "1.0" );
  MpIeee sum_err=  MpIeee( "0.0" );
  gsl_sf_result lnsheta;
  MpIeee lnN;
  MpIeee lnpre_val;MpIeee  lnpre_err;MpIeee  lnprepow;
  int  stat_e;
  int  n;

  gsl_sf_lngamma_e(ell + 3.0/2.0, &lg_lp32);
  gsl_sf_lnsinh_e(eta, &lnsheta);
  legendre_H3d_lnnorm(ell, lambda, &lnN);
  lnprepow = MpIeee( "0.5" )*(ell + MpIeee( "0.5" )) * (ln_zm1 - ln_zp1);
  lnpre_val  = lnprepow + MpIeee( "0.5" )*(lnN + M_LNPI - M_LN2 - lnsheta.val) - lg_lp32.val - log(fabs(lambda));
  lnpre_err  = lnsheta.err + lg_lp32.err + GSL_DBL_EPSILON * fabs(lnpre_val);
  lnpre_err += MpIeee( "2.0" )*GSL_DBL_EPSILON * (fabs(lnN) + M_LNPI + M_LN2);
  lnpre_err += MpIeee( "2.0" )*GSL_DBL_EPSILON * (MpIeee( "0.5" )*(ell + MpIeee( "0.5" )) * (fabs(ln_zm1) + fabs(ln_zp1)));
  for(n=1; n<nmax; n++) {
    MpIeee aR=  n - MpIeee( "0.5" );
    term *= (aR*aR + lambda*lambda)*zeta/(ell + n + MpIeee( "0.5" ))/n;
    sum  += term;
    sum_err += MpIeee( "2.0" )*GSL_DBL_EPSILON*fabs(term);
    if(fabs(term/sum) < 2.0 * GSL_DBL_EPSILON) break;
  }

  stat_e = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err, sum, fabs(term)+sum_err, result);
  return GSL_ERROR_SELECT_2(stat_e, (n==nmax ? GSL_EMAXITER : GSL_SUCCESS));
}


/* Evaluate legendre_H3d(ell+1)/legendre_H3d(ell)
 * by continued fraction.
 */
#if 0
static
int
 legendre_H3d_CF1(const int ell, const MpIeee lambda, const MpIeee coth_eta,
                 gsl_sf_result * result)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  sqrt(lambda*lambda + (ell+MpIeee( "1.0" ))*(ell+MpIeee( "1.0" )));
  MpIeee b1=  (MpIeee( "2.0" )*ell + MpIeee( "3.0" )) * coth_eta;
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
    an = -(lambda*lambda + ((MpIeee)ell + n)*((MpIeee)ell + n));
    bn = (MpIeee( "2.0" )*ell + MpIeee( "2.0" )*n + MpIeee( "1.0" )) * coth_eta;
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
    
    if(fabs(del - 1.0) < 4.0*GSL_DBL_EPSILON) break;
  }

  result->val = fn;
  result->err = 2.0 * GSL_DBL_EPSILON * (sqrt(n)+1.0) * fabs(fn);

  if(n >= maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}
#endif /* 0 */


/* Evaluate legendre_H3d(ell+1)/legendre_H3d(ell)
 * by continued fraction. Use the Gautschi (Euler)
 * equivalent series.
 */
 /* FIXME: Maybe we have to worry about this. The a_k are
  * not positive and there can be a blow-up. It happened
  * for J_nu once or twice. Then we should probably use
  * the method above.
  */
static
int
 legendre_H3d_CF1_ser(const int ell, const MpIeee lambda, const MpIeee coth_eta,
                     gsl_sf_result * result)
{
  const MpIeee pre=  sqrt(lambda*lambda+(ell+1.0)*(ell+1.0))/((2.0*ell+3)*coth_eta);
  const int maxk = 20000;
  MpIeee tk=  MpIeee( "1.0" );
  MpIeee sum=  MpIeee( "1.0" );
  MpIeee rhok=  MpIeee( "0.0" );
  MpIeee sum_err=  MpIeee( "0.0" );
  int  k;
 
  for(k=1; k<maxk; k++) {
    MpIeee tlk=  (MpIeee( "2.0" )*ell + MpIeee( "1.0" ) + MpIeee( "2.0" )*k);
    MpIeee l1k=  (ell + MpIeee( "1.0" ) + k);
    MpIeee ak=  -(lambda*lambda + l1k*l1k)/(tlk*(tlk+MpIeee( "2.0" ))*coth_eta*coth_eta);
    rhok = -ak*(MpIeee( "1.0" ) + rhok)/(MpIeee( "1.0" ) + ak*(MpIeee( "1.0" ) + rhok));
    tk  *= rhok;
    sum += tk;
    sum_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * k * fabs(tk);
    if(fabs(tk/sum) < GSL_DBL_EPSILON) break;
  }

  result->val  = pre * sum;
  result->err  = fabs(pre * tk);
  result->err += fabs(pre * sum_err);
  result->err += 4.0 * GSL_DBL_EPSILON * fabs(result->val);

  if(k >= maxk)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_legendre_H3d_0_e(const MpIeee lambda, const MpIeee eta, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(eta < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(eta == 0.0 || lambda == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const MpIeee lam_eta=  lambda * eta;
    gsl_sf_result s;
    gsl_sf_sin_err_e(lam_eta, 2.0*GSL_DBL_EPSILON * fabs(lam_eta), &s);
    if(eta > -0.5*GSL_LOG_DBL_EPSILON) {
      MpIeee f=  MpIeee( "2.0" ) / lambda * exp(-eta);
      result->val  = f * s.val;
      result->err  = fabs(f * s.val) * (fabs(eta) + 1.0) * GSL_DBL_EPSILON;
      result->err += fabs(f) * s.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
      MpIeee f=  MpIeee( "1.0" )/(lambda*sinh(eta));
      result->val  = f * s.val;
      result->err  = fabs(f * s.val) * (fabs(eta) + 1.0) * GSL_DBL_EPSILON;
      result->err += fabs(f) * s.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_legendre_H3d_1_e(const MpIeee lambda, const MpIeee eta, gsl_sf_result * result)
{
  const MpIeee xi=  fabs(eta*lambda);
  const MpIeee lsq=  lambda*lambda;
  const MpIeee lsqp1=  lsq + 1.0;

  /* CHECK_POINTER(result) */

  if(eta < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(eta == 0.0 || lambda == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(xi < GSL_ROOT5_DBL_EPSILON && eta < GSL_ROOT5_DBL_EPSILON) {
    MpIeee etasq=  eta*eta;
    MpIeee xisq=  xi*xi;
    MpIeee term1=  (etasq + xisq)/MpIeee( "3.0" );
    MpIeee term2=  -(MpIeee( "2.0" )*etasq*etasq + MpIeee( "5.0" )*etasq*xisq + MpIeee( "3.0" )*xisq*xisq)/MpIeee( "90.0" );
    MpIeee sinh_term=  MpIeee( "1.0" ) - eta*eta/MpIeee( "6.0" ) * (MpIeee( "1.0" ) - MpIeee( "7.0" )/MpIeee( "60.0" )*eta*eta);
    MpIeee pre=  sinh_term/sqrt(lsqp1) / eta;
    result->val  = pre * (term1 + term2);
    result->err  = pre * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee sin_term;     /*  Sin(xi)/xi     */
    MpIeee cos_term;     /*  Cos(xi)        */
    MpIeee coth_term;    /*  eta/Tanh(eta)  */
    MpIeee sinh_term;    /*  eta/Sinh(eta)  */
    MpIeee sin_term_err;
    MpIeee cos_term_err;
    MpIeee t1;
    MpIeee pre_val;
    MpIeee pre_err;
    MpIeee term1;
    MpIeee term2;
    if(xi < GSL_ROOT5_DBL_EPSILON) {
      sin_term = MpIeee( "1.0" ) - xi*xi/MpIeee( "6.0" ) * (MpIeee( "1.0" ) - xi*xi/MpIeee( "20.0" ));
      cos_term = MpIeee( "1.0" ) - MpIeee( "0.5" )*xi*xi * (MpIeee( "1.0" ) - xi*xi/MpIeee( "12.0" ));
      sin_term_err = GSL_DBL_EPSILON;
      cos_term_err = GSL_DBL_EPSILON;
    }
    else {
      gsl_sf_result sin_xi_result;
      gsl_sf_result cos_xi_result;
      gsl_sf_sin_e(xi, &sin_xi_result);
      gsl_sf_cos_e(xi, &cos_xi_result);
      sin_term = sin_xi_result.val/xi;
      cos_term = cos_xi_result.val;
      sin_term_err = sin_xi_result.err/fabs(xi);
      cos_term_err = cos_xi_result.err;
    }
    if(eta < GSL_ROOT5_DBL_EPSILON) {
      coth_term = MpIeee( "1.0" ) + eta*eta/MpIeee( "3.0" ) * (MpIeee( "1.0" ) - eta*eta/MpIeee( "15.0" ));
      sinh_term = MpIeee( "1.0" ) - eta*eta/MpIeee( "6.0" ) * (MpIeee( "1.0" ) - MpIeee( "7.0" )/MpIeee( "60.0" )*eta*eta);
    }
    else {
      coth_term = eta/tanh(eta);
      sinh_term = eta/sinh(eta);
    }
    t1 = sqrt(lsqp1) * eta;
    pre_val = sinh_term/t1;
    pre_err = MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(pre_val);
    term1 = sin_term*coth_term;
    term2 = cos_term;
    result->val  = pre_val * (term1 - term2);
    result->err  = pre_err * fabs(term1 - term2);
    result->err += pre_val * (sin_term_err * coth_term + cos_term_err);
    result->err += pre_val * fabs(term1-term2) * (fabs(eta) + 1.0) * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_legendre_H3d_e(const int ell, const MpIeee lambda, const MpIeee eta,
                         gsl_sf_result * result)
{
  const MpIeee abs_lam=  fabs(lambda);
  const MpIeee lsq=  abs_lam*abs_lam;
  const MpIeee xi=  abs_lam * eta;
  const MpIeee cosh_eta=  cosh(eta);

  /* CHECK_POINTER(result) */

  if(eta < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(eta > GSL_LOG_DBL_MAX) {
    /* cosh(eta) is too big. */
    OVERFLOW_ERROR(result);
  }
  else if(ell == 0) {
    return gsl_sf_legendre_H3d_0_e(lambda, eta, result);
  }
  else if(ell == 1) {
    return gsl_sf_legendre_H3d_1_e(lambda, eta, result);
  }
  else if(eta == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(xi < 1.0) {
    return legendre_H3d_series(ell, lambda, eta, result);
  }
  else if((ell*ell+lsq)/sqrt(1.0+lsq)/(cosh_eta*cosh_eta) < 5.0*GSL_ROOT3_DBL_EPSILON) {
    /* Large argument.
     */
    gsl_sf_result P;
    MpIeee lm;
    int  stat_P=  gsl_sf_conicalP_large_x_e(-ell-0.5, lambda, cosh_eta, &P, &lm);
    if(P.val == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return stat_P;
    }
    else {
      MpIeee lnN;
      gsl_sf_result lnsh;
      MpIeee ln_abslam;
      MpIeee lnpre_val;MpIeee  lnpre_err;
      int  stat_e;
      gsl_sf_lnsinh_e(eta, &lnsh);
      legendre_H3d_lnnorm(ell, lambda, &lnN);
      ln_abslam = log(abs_lam);
      lnpre_val  = MpIeee( "0.5" )*(M_LNPI + lnN - M_LN2 - lnsh.val) - ln_abslam;
      lnpre_err  = lnsh.err;
      lnpre_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * (MpIeee( "0.5" )*(M_LNPI + M_LN2 + fabs(lnN)) + fabs(ln_abslam));
      lnpre_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(lnpre_val);
      stat_e = gsl_sf_exp_mult_err_e(lnpre_val + lm, lnpre_err, P.val, P.err, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_P);
    }
  }
  else if(abs_lam > 1000.0*ell*ell) {
    /* Large degree.
     */
    gsl_sf_result P;
    MpIeee lm;
    int  stat_P=  gsl_sf_conicalP_xgt1_neg_mu_largetau_e(ell+0.5,
                                                           lambda,
                                                           cosh_eta, eta,
                                                           &P, &lm);
    if(P.val == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return stat_P;
    }
    else {
      MpIeee lnN;
      gsl_sf_result lnsh;
      MpIeee ln_abslam;
      MpIeee lnpre_val;MpIeee  lnpre_err;
      int  stat_e;
      gsl_sf_lnsinh_e(eta, &lnsh);
      legendre_H3d_lnnorm(ell, lambda, &lnN);
      ln_abslam = log(abs_lam);
      lnpre_val  = MpIeee( "0.5" )*(M_LNPI + lnN - M_LN2 - lnsh.val) - ln_abslam;
      lnpre_err  = lnsh.err;
      lnpre_err += GSL_DBL_EPSILON * (MpIeee( "0.5" )*(M_LNPI + M_LN2 + fabs(lnN)) + fabs(ln_abslam));
      lnpre_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(lnpre_val);
      stat_e = gsl_sf_exp_mult_err_e(lnpre_val + lm, lnpre_err, P.val, P.err, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_P);
    }
  }
  else {
    /* Backward recurrence.
     */
    const MpIeee coth_eta=  1.0/tanh(eta);
    const MpIeee coth_err_mult=  fabs(eta) + 1.0;
    gsl_sf_result rH;
    int  stat_CF1=  legendre_H3d_CF1_ser(ell, lambda, coth_eta, &rH);
    MpIeee Hlm1;
    MpIeee Hl=  GSL_SQRT_DBL_MIN;
    MpIeee Hlp1=  rH.val * Hl;
    int  lp;
    for(lp=ell; lp>0; lp--) {
      MpIeee root_term_0=  sqrt(lambda*lambda + (MpIeee)lp*lp);
      MpIeee root_term_1=  sqrt(lambda*lambda + (lp+MpIeee( "1.0" ))*(lp+MpIeee( "1.0" )));
      Hlm1 = ((MpIeee( "2.0" )*lp + MpIeee( "1.0" ))*coth_eta*Hl - root_term_1 * Hlp1)/root_term_0;
      Hlp1 = Hl;
      Hl   = Hlm1;
    }

    if(fabs(Hl) > fabs(Hlp1)) {
      gsl_sf_result H0;
      int  stat_H0=  gsl_sf_legendre_H3d_0_e(lambda, eta, &H0);
      result->val  = GSL_SQRT_DBL_MIN/Hl * H0.val;
      result->err  = GSL_SQRT_DBL_MIN/fabs(Hl) * H0.err;
      result->err += fabs(rH.err/rH.val) * (ell+1.0) * coth_err_mult * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_H0, stat_CF1);
    }
    else {
      gsl_sf_result H1;
      int  stat_H1=  gsl_sf_legendre_H3d_1_e(lambda, eta, &H1);
      result->val  = GSL_SQRT_DBL_MIN/Hlp1 * H1.val;
      result->err  = GSL_SQRT_DBL_MIN/fabs(Hlp1) * H1.err;
      result->err += fabs(rH.err/rH.val) * (ell+1.0) * coth_err_mult * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_H1, stat_CF1);
    }
  }
}


int
 gsl_sf_legendre_H3d_array(const int lmax, const MpIeee lambda, const MpIeee eta, MpIeee * result_array)
{
  /* CHECK_POINTER(result_array) */

 if(eta < 0.0 || lmax < 0) {
    int  ell;
    for(ell=0; ell<=lmax; ell++) result_array[ell] = MpIeee( "0.0" );
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(eta > GSL_LOG_DBL_MAX) {
    /* cosh(eta) is too big. */
    int  ell;
    for(ell=0; ell<=lmax; ell++) result_array[ell] = MpIeee( "0.0" );
    GSL_ERROR ("overflow", GSL_EOVRFLW);
  }
  else if(lmax == 0) {
    gsl_sf_result H0;
    int  stat=  gsl_sf_legendre_H3d_e(0, lambda, eta, &H0);
    result_array[0] = H0.val;
    return stat;
  }
  else {
    /* Not the most efficient method. But what the hell... it's simple.
     */
    gsl_sf_result r_Hlp1;
    gsl_sf_result r_Hl;
    int  stat_lmax=  gsl_sf_legendre_H3d_e(lmax,   lambda, eta, &r_Hlp1);
    int  stat_lmaxm1=  gsl_sf_legendre_H3d_e(lmax-1, lambda, eta, &r_Hl);
    int  stat_max=  GSL_ERROR_SELECT_2(stat_lmax, stat_lmaxm1);

    const MpIeee coth_eta=  1.0/tanh(eta);
    int  stat_recursion=  GSL_SUCCESS;
    MpIeee Hlp1=  r_Hlp1.val;
    MpIeee Hl=  r_Hl.val;
    MpIeee Hlm1;
    int  ell;

    result_array[lmax]   = Hlp1;
    result_array[lmax-1] = Hl;

    for(ell=lmax-1; ell>0; ell--) {
      MpIeee root_term_0=  sqrt(lambda*lambda + (MpIeee)ell*ell);
      MpIeee root_term_1=  sqrt(lambda*lambda + (ell+MpIeee( "1.0" ))*(ell+MpIeee( "1.0" )));
      Hlm1 = ((MpIeee( "2.0" )*ell + MpIeee( "1.0" ))*coth_eta*Hl - root_term_1 * Hlp1)/root_term_0;
      result_array[ell-1] = Hlm1;
      if(!(Hlm1 < GSL_DBL_MAX)) stat_recursion = GSL_EOVRFLW;
      Hlp1 = Hl;
      Hl   = Hlm1;
    }

    return GSL_ERROR_SELECT_2(stat_recursion, stat_max);
  }
}
  

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_legendre_H3d_0(const MpIeee lambda, const MpIeee eta)
{
  EVAL_RESULT(gsl_sf_legendre_H3d_0_e(lambda, eta, &result));
}

MpIeee gsl_sf_legendre_H3d_1(const MpIeee lambda, const MpIeee eta)
{
  EVAL_RESULT(gsl_sf_legendre_H3d_1_e(lambda, eta, &result));
}

MpIeee gsl_sf_legendre_H3d(const int l, const MpIeee lambda, const MpIeee eta)
{
  EVAL_RESULT(gsl_sf_legendre_H3d_e(l, lambda, eta, &result));
}
