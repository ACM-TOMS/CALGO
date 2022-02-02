#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel.c
 * 
 * Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003 Gerard Jungman
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
/* Miscellaneous support functions for Bessel function evaluations.
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>

#include "error.h"

#include "bessel_amp_phase.h"
#include "bessel_temme.h"
#include "bessel.h"

#define CubeRoot2_  1.25992104989487316476721060728



/* Debye functions [Abramowitz+Stegun, 9.3.9-10] */

inline static MpIeee debye_u1(const MpIeee * tpow)
{
  return (MpIeee( "3.0" )*tpow[1] - MpIeee( "5.0" )*tpow[3])/MpIeee( "24.0" );
}

inline static MpIeee debye_u2(const MpIeee * tpow)
{
  return (MpIeee( "81.0" )*tpow[2] - MpIeee( "462.0" )*tpow[4] + MpIeee( "385.0" )*tpow[6])/MpIeee( "1152.0" );
}

inline
static MpIeee debye_u3(const MpIeee * tpow)
{
  return (MpIeee( "30375.0" )*tpow[3] - MpIeee( "369603.0" )*tpow[5] + MpIeee( "765765.0" )*tpow[7] - MpIeee( "425425.0" )*tpow[9])/MpIeee( "414720.0" );
}

inline
static MpIeee debye_u4(const MpIeee * tpow)
{
  return (MpIeee( "4465125.0" )*tpow[4] - MpIeee( "94121676.0" )*tpow[6] + MpIeee( "349922430.0" )*tpow[8] - 
          MpIeee( "446185740.0" )*tpow[10] + MpIeee( "185910725.0" )*tpow[12])/MpIeee( "39813120.0" );
}

inline
static MpIeee debye_u5(const MpIeee * tpow)
{
  return (MpIeee( "1519035525.0" )*tpow[5]     - MpIeee( "49286948607.0" )*tpow[7] + 
          MpIeee( "284499769554.0" )*tpow[9]   - MpIeee( "614135872350.0" )*tpow[11] + 
          MpIeee( "566098157625.0" )*tpow[13]  - MpIeee( "188699385875.0" )*tpow[15])/MpIeee( "6688604160.0" );
}

#if 0
inline
static MpIeee debye_u6(const MpIeee * tpow)
{
  return (MpIeee( "2757049477875.0" )*tpow[6] - MpIeee( "127577298354750.0" )*tpow[8] + 
          MpIeee( "1050760774457901.0" )*tpow[10] - MpIeee( "3369032068261860.0" )*tpow[12] + 
          MpIeee( "5104696716244125.0" )*tpow[14] - MpIeee( "3685299006138750.0" )*tpow[16] + 
          MpIeee( "1023694168371875.0" )*tpow[18])/MpIeee( "4815794995200.0" );
}
#endif


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_bessel_IJ_taylor_e(const MpIeee nu, const MpIeee x,
                             const int sign,
                             const int kmax,
                             const MpIeee threshold,
                             gsl_sf_result * result
                             )
{
  /* CHECK_POINTER(result) */

  if(nu < 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    if(nu == 0.0) {
      result->val = 1.0;
      result->err = 0.0;
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
    }
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result prefactor;   /* (x/2)^nu / Gamma(nu+1) */
    gsl_sf_result sum;

    int  stat_pre;
    int  stat_sum;
    int  stat_mul;

    if(nu == 0.0) {
      prefactor.val = 1.0;
      prefactor.err = 0.0;
      stat_pre = GSL_SUCCESS;
    }
    else if(nu < INT_MAX-1) {
      /* Separate the integer part and use
       * y^nu / Gamma(nu+1) = y^N /N! y^f / (N+1)_f,
       * to control the error.
       */
      const int    N = (int)floor(nu + 0.5).toInt();
      const MpIeee f=  nu - N;
      gsl_sf_result poch_factor;
      gsl_sf_result tc_factor;
      const int stat_poch = gsl_sf_poch_e(N+1.0, f, &poch_factor);
      const int stat_tc   = gsl_sf_taylorcoeff_e(N, 0.5*x, &tc_factor);
      const MpIeee p=  pow(0.5*x,f);
      prefactor.val  = tc_factor.val * p / poch_factor.val;
      prefactor.err  = tc_factor.err * p / poch_factor.val;
      prefactor.err += fabs(prefactor.val) / poch_factor.val * poch_factor.err;
      prefactor.err += 2.0 * GSL_DBL_EPSILON * fabs(prefactor.val);
      stat_pre = GSL_ERROR_SELECT_2(stat_tc, stat_poch);
    }
    else {
      gsl_sf_result lg;
      const int stat_lg = gsl_sf_lngamma_e(nu+1.0, &lg);
      const MpIeee term1=  nu*log(0.5*x);
      const MpIeee term2=  lg.val;
      const MpIeee ln_pre=  term1 - term2;
      const MpIeee ln_pre_err=  GSL_DBL_EPSILON * (fabs(term1)+fabs(term2)) + lg.err;
      const int stat_ex = gsl_sf_exp_err_e(ln_pre, ln_pre_err, &prefactor);
      stat_pre = GSL_ERROR_SELECT_2(stat_ex, stat_lg);
    }

    /* Evaluate the sum.
     * [Abramowitz+Stegun, 9.1.10]
     * [Abramowitz+Stegun, 9.6.7]
     */
    {
      const MpIeee y=  sign * 0.25 * x*x;
      MpIeee sumk=  MpIeee( "1.0" );
      MpIeee term=  MpIeee( "1.0" );
      int  k;

      for(k=1; k<=kmax; k++) {
        term *= y/((nu+k)*k);
        sumk += term;
        if(fabs(term/sumk) < threshold) break;
      }

      sum.val = sumk;
      sum.err = threshold * fabs(sumk);

      stat_sum = ( k >= kmax ? GSL_EMAXITER : GSL_SUCCESS );
    }

    stat_mul = gsl_sf_multiply_err_e(prefactor.val, prefactor.err,
                                        sum.val, sum.err,
                                        result);

    return GSL_ERROR_SELECT_3(stat_mul, stat_pre, stat_sum);
  }
}


/* x >> nu*nu+1
 * error ~ O( ((nu*nu+1)/x)^4 )
 *
 * empirical error analysis:
 *   choose  GSL_ROOT4_MACH_EPS * x > (nu*nu + 1)
 *
 * This is not especially useful. When the argument gets
 * large enough for this to apply, the cos() and sin()
 * start loosing digits. However, this seems inevitable
 * for this particular method.
 *
 * Wed Jun 25 14:39:38 MDT 2003 [GJ]
 * This function was inconsistent since the Q term did not
 * go to relative order eps^2. That's why the error estimate
 * originally given was screwy (it didn't make sense that the
 * "empirical" error was coming out O(eps^3)).
 * With Q to proper order, the error is O(eps^4).
 */
int
 gsl_sf_bessel_Jnu_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  MpIeee mu=  MpIeee( "4.0" )*nu*nu;
  MpIeee mum1=  mu-1.0;
  MpIeee mum9=  mu-MpIeee( "9.0" );
  MpIeee mum25=  mu-25.0;
  MpIeee chi=  x - (MpIeee( "0.5" )*nu + MpIeee( "0.25" ))*M_PI;
  MpIeee P=  MpIeee( "1.0" ) - mum1*mum9/(MpIeee( "128.0" )*x*x);
  MpIeee Q=  mum1/(MpIeee( "8.0" )*x) * (MpIeee( "1.0" ) - mum9*mum25/(MpIeee( "384.0" )*x*x));
  MpIeee pre=  sqrt(MpIeee( "2.0" )/(M_PI*x));
  MpIeee c=  cos(chi);
  MpIeee s=  sin(chi);
  MpIeee r=  mu/x;
  result->val  = pre * (c*P - s*Q);
  result->err  = pre * GSL_DBL_EPSILON * (1.0 + fabs(x)) * (fabs(c*P) + fabs(s*Q));
  result->err += pre * fabs(0.1*r*r*r*r);
  return GSL_SUCCESS;
}


/* x >> nu*nu+1
 */
int
 gsl_sf_bessel_Ynu_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  MpIeee ampl;
  MpIeee theta;
  MpIeee alpha=  x;
  MpIeee beta=  -MpIeee( "0.5" )*nu*M_PI;
  int  stat_a=  gsl_sf_bessel_asymp_Mnu_e(nu, x, &ampl);
  int  stat_t=  gsl_sf_bessel_asymp_thetanu_corr_e(nu, x, &theta);
  MpIeee sin_alpha=  sin(alpha);
  MpIeee cos_alpha=  cos(alpha);
  MpIeee sin_chi=  sin(beta + theta);
  MpIeee cos_chi=  cos(beta + theta);
  MpIeee sin_term=  sin_alpha * cos_chi + sin_chi * cos_alpha;
  MpIeee sin_term_mag=  fabs(sin_alpha * cos_chi) + fabs(sin_chi * cos_alpha);
  result->val  = ampl * sin_term;
  result->err  = fabs(ampl) * GSL_DBL_EPSILON * sin_term_mag;
  result->err += fabs(result->val) * 2.0 * GSL_DBL_EPSILON;

  if(fabs(alpha) > 1.0/GSL_DBL_EPSILON) {
    result->err *= 0.5 * fabs(alpha);
  }
  else if(fabs(alpha) > 1.0/GSL_SQRT_DBL_EPSILON) {
    result->err *= 256.0 * fabs(alpha) * GSL_SQRT_DBL_EPSILON;
  }

  return GSL_ERROR_SELECT_2(stat_t, stat_a);
}


/* x >> nu*nu+1
 */
int
 gsl_sf_bessel_Inu_scaled_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  MpIeee mu=  MpIeee( "4.0" )*nu*nu;
  MpIeee mum1=  mu-1.0;
  MpIeee mum9=  mu-MpIeee( "9.0" );
  MpIeee pre=  MpIeee( "1.0" )/sqrt(MpIeee( "2.0" )*M_PI*x);
  MpIeee r=  mu/x;
  result->val = pre * (1.0 - mum1/(8.0*x) + mum1*mum9/(128.0*x*x));
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + pre * fabs(0.1*r*r*r);
  return GSL_SUCCESS;
}

/* x >> nu*nu+1
 */
int
 gsl_sf_bessel_Knu_scaled_asympx_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  MpIeee mu=  MpIeee( "4.0" )*nu*nu;
  MpIeee mum1=  mu-MpIeee( "1.0" );
  MpIeee mum9=  mu-9.0;
  MpIeee pre=  sqrt(M_PI/(MpIeee( "2.0" )*x));
  MpIeee r=  nu/x;
  result->val = pre * (1.0 + mum1/(8.0*x) + mum1*mum9/(128.0*x*x));
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + pre * fabs(0.1*r*r*r);
  return GSL_SUCCESS;
}


/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
 *
 * error:
 *   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
 *   It is not hard to show that |u_N(t)| is small for such t.
 *   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
 *   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
 *   seen below as nu ~ 100. For general MACH_EPS it will be 
 *                     nu > 0.5 / MACH_EPS^(1/6)
 *   When t is small, the bound is even better because |u_N(t)| vanishes
 *   as t->0. In fact u_N(t) ~ C t^N as t->0, with C ~= 0.1.
 *   We write
 *                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
 *   therefore
 *                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
 *   and this is the general form.
 *
 * empirical error analysis, assuming 14 digit requirement:
 *   choose   x > 50.000 nu   ==>  nu >   3
 *   choose   x > 10.000 nu   ==>  nu >  15
 *   choose   x >  2.000 nu   ==>  nu >  50
 *   choose   x >  1.000 nu   ==>  nu >  75
 *   choose   x >  0.500 nu   ==>  nu >  80
 *   choose   x >  0.100 nu   ==>  nu >  83
 *
 * This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
 * since the polynomial term will be evaluated near t=1, so the bound
 * on nu will become constant for small x. Furthermore, increasing x with
 * nu fixed will decrease the error.
 */
int
 gsl_sf_bessel_Inu_scaled_asymp_unif_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  int  i;
  MpIeee z=  x/nu;
  MpIeee root_term=  sqrt(MpIeee( "1.0" ) + z*z);
  MpIeee pre=  MpIeee( "1.0" )/sqrt(MpIeee( "2.0" )*M_PI*nu * root_term);
  MpIeee eta=  root_term + log(z/(MpIeee( "1.0" )+root_term));
  MpIeee ex_arg=  ( z < MpIeee( "1.0" )/GSL_ROOT3_DBL_EPSILON ? nu*(-z + eta) : -MpIeee( "0.5" )*nu/z*(MpIeee( "1.0" ) - MpIeee( "1.0" )/(MpIeee( "12.0" )*z*z)) );
  gsl_sf_result ex_result;
  int  stat_ex=  gsl_sf_exp_e(ex_arg, &ex_result);
  if(stat_ex == GSL_SUCCESS) {
    MpIeee t=  MpIeee( "1.0" )/root_term;
    MpIeee sum;
    MpIeee tpow[16];
    tpow[0] = MpIeee( "1.0" );
    for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
    sum = MpIeee( "1.0" ) + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + debye_u3(tpow)/(nu*nu*nu)
          + debye_u4(tpow)/(nu*nu*nu*nu) + debye_u5(tpow)/(nu*nu*nu*nu*nu);
    result->val  = pre * ex_result.val * sum;
    result->err  = pre * ex_result.val / (nu*nu*nu*nu*nu*nu);
    result->err += pre * ex_result.err * fabs(sum);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 0.0;
    result->err = 0.0;
    return stat_ex;
  }
}


/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.8]
 *
 * error:
 *   identical to that above for Inu_scaled
 */
int
 gsl_sf_bessel_Knu_scaled_asymp_unif_e(const MpIeee nu, const MpIeee x, gsl_sf_result * result)
{
  int  i;
  MpIeee z=  x/nu;
  MpIeee root_term=  sqrt(MpIeee( "1.0" ) + z*z);
  MpIeee pre=  sqrt(M_PI/(MpIeee( "2.0" )*nu*root_term));
  MpIeee eta=  root_term + log(z/(MpIeee( "1.0" )+root_term));
  MpIeee ex_arg=  ( z < MpIeee( "1.0" )/GSL_ROOT3_DBL_EPSILON ? nu*(z - eta) : MpIeee( "0.5" )*nu/z*(MpIeee( "1.0" ) + MpIeee( "1.0" )/(MpIeee( "12.0" )*z*z)) );
  gsl_sf_result ex_result;
  int  stat_ex=  gsl_sf_exp_e(ex_arg, &ex_result);
  if(stat_ex == GSL_SUCCESS) {
    MpIeee t=  MpIeee( "1.0" )/root_term;
    MpIeee sum;
    MpIeee tpow[16];
    tpow[0] = MpIeee( "1.0" );
    for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
    sum = MpIeee( "1.0" ) - debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) - debye_u3(tpow)/(nu*nu*nu)
          + debye_u4(tpow)/(nu*nu*nu*nu) - debye_u5(tpow)/(nu*nu*nu*nu*nu);
    result->val  = pre * ex_result.val * sum;
    result->err  = pre * ex_result.err * fabs(sum);
    result->err += pre * ex_result.val / (nu*nu*nu*nu*nu*nu);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 0.0;
    result->err = 0.0;
    return stat_ex;
  }
}


/* Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
 */
int
 gsl_sf_bessel_JY_mu_restricted(const MpIeee mu, const MpIeee x,
                               gsl_sf_result * Jmu, gsl_sf_result * Jmup1,
                               gsl_sf_result * Ymu, gsl_sf_result * Ymup1)
{
  /* CHECK_POINTER(Jmu) */
  /* CHECK_POINTER(Jmup1) */
  /* CHECK_POINTER(Ymu) */
  /* CHECK_POINTER(Ymup1) */

  if(x < 0.0 || fabs(mu) > 0.5) {
    Jmu->val   = 0.0;
    Jmu->err   = 0.0;
    Jmup1->val = 0.0;
    Jmup1->err = 0.0;
    Ymu->val   = 0.0;
    Ymu->err   = 0.0;
    Ymup1->val = 0.0;
    Ymup1->err = 0.0;
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(x == 0.0) {
    if(mu == 0.0) {
      Jmu->val   = 1.0;
      Jmu->err   = 0.0;
    }
    else {
      Jmu->val   = 0.0;
      Jmu->err   = 0.0;
    }
    Jmup1->val = 0.0;
    Jmup1->err = 0.0;
    Ymu->val   = 0.0;
    Ymu->err   = 0.0;
    Ymup1->val = 0.0;
    Ymup1->err = 0.0;
    GSL_ERROR ("error", GSL_EDOM);
  }
  else {
    int  stat_Y;
    int  stat_J;

    if(x < 2.0) {
      /* Use Taylor series for J and the Temme series for Y.
       * The Taylor series for J requires nu > 0, so we shift
       * up one and use the recursion relation to get Jmu, in
       * case mu < 0.
       */
      gsl_sf_result Jmup2;
      int  stat_J1=  gsl_sf_bessel_IJ_taylor_e(mu+1.0, x, -1, 100, GSL_DBL_EPSILON,  Jmup1);
      int  stat_J2=  gsl_sf_bessel_IJ_taylor_e(mu+2.0, x, -1, 100, GSL_DBL_EPSILON, &Jmup2);
      MpIeee c=  MpIeee( "2.0" )*(mu+MpIeee( "1.0" ))/x;
      Jmu->val  = c * Jmup1->val - Jmup2.val;
      Jmu->err  = c * Jmup1->err + Jmup2.err;
      Jmu->err += 2.0 * GSL_DBL_EPSILON * fabs(Jmu->val);
      stat_J = GSL_ERROR_SELECT_2(stat_J1, stat_J2);
      stat_Y = gsl_sf_bessel_Y_temme(mu, x, Ymu, Ymup1);
      return GSL_ERROR_SELECT_2(stat_J, stat_Y);
    }
    else if(x < 1000.0) {
      MpIeee P;MpIeee  Q;
      MpIeee J_ratio;
      MpIeee J_sgn;
      const int stat_CF1 = gsl_sf_bessel_J_CF1(mu, x, &J_ratio, &J_sgn);
      const int stat_CF2 = gsl_sf_bessel_JY_steed_CF2(mu, x, &P, &Q);
      MpIeee Jprime_J_ratio=  mu/x - J_ratio;
      MpIeee gamma=  (P - Jprime_J_ratio)/Q;
      Jmu->val = J_sgn * sqrt(2.0/(M_PI*x) / (Q + gamma*(P-Jprime_J_ratio)));
      Jmu->err = 4.0 * GSL_DBL_EPSILON * fabs(Jmu->val);
      Jmup1->val = J_ratio * Jmu->val;
      Jmup1->err = fabs(J_ratio) * Jmu->err;
      Ymu->val = gamma * Jmu->val;
      Ymu->err = fabs(gamma) * Jmu->err;
      Ymup1->val = Ymu->val * (mu/x - P - Q/gamma);
      Ymup1->err = Ymu->err * fabs(mu/x - P - Q/gamma) + 4.0*GSL_DBL_EPSILON*fabs(Ymup1->val);
      return GSL_ERROR_SELECT_2(stat_CF1, stat_CF2);
    }
    else {
      /* Use asymptotics for large argument.
       */
      const int stat_J0 = gsl_sf_bessel_Jnu_asympx_e(mu,     x, Jmu);
      const int stat_J1 = gsl_sf_bessel_Jnu_asympx_e(mu+1.0, x, Jmup1);
      const int stat_Y0 = gsl_sf_bessel_Ynu_asympx_e(mu,     x, Ymu);
      const int stat_Y1 = gsl_sf_bessel_Ynu_asympx_e(mu+1.0, x, Ymup1);
      stat_J = GSL_ERROR_SELECT_2(stat_J0, stat_J1);
      stat_Y = GSL_ERROR_SELECT_2(stat_Y0, stat_Y1);
      return GSL_ERROR_SELECT_2(stat_J, stat_Y);
    }
  }
}


int
 gsl_sf_bessel_J_CF1(const MpIeee nu, const MpIeee x,
                    MpIeee * ratio, MpIeee * sgn)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 10000;
  int  n=  1;
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  x/(MpIeee( "2.0" )*(nu+MpIeee( "1.0" )));
  MpIeee An=  Anm1 + a1*Anm2;
  MpIeee Bn=  Bnm1 + a1*Bnm2;
  MpIeee an;
  MpIeee fn=  An/Bn;
  MpIeee dn=  a1;
  MpIeee s=  MpIeee( "1.0" );

  while(n < maxiter) {
    MpIeee old_fn;
    MpIeee del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = -x*x/(MpIeee( "4.0" )*(nu+n-MpIeee( "1.0" ))*(nu+n));
    An = Anm1 + an*Anm2;
    Bn = Bnm1 + an*Bnm2;

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

    dn = MpIeee( "1.0" ) / (MpIeee( "2.0" )*(nu+n)/x - dn);
    if(dn < MpIeee( "0.0" )) s = -s;

    if(fabs(del - 1.0) < 2.0*GSL_DBL_EPSILON) break;
  }

  *ratio = fn;
  *sgn   = s;

  if(n >= maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}



/* Evaluate the continued fraction CF1 for J_{nu+1}/J_nu
 * using Gautschi (Euler) equivalent series.
 * This exhibits an annoying problem because the
 * a_k are not positive definite (in fact they are all negative).
 * There are cases when rho_k blows up. Example: nu=1,x=4.
 */
#if 0
int
 gsl_sf_bessel_J_CF1_ser(const MpIeee nu, const MpIeee x,
                        MpIeee * ratio, MpIeee * sgn)
{
  const int maxk = 20000;
  MpIeee tk=  MpIeee( "1.0" );
  MpIeee sum=  MpIeee( "1.0" );
  MpIeee rhok=  MpIeee( "0.0" );
  MpIeee dk=  MpIeee( "0.0" );
  MpIeee s=  MpIeee( "1.0" );
  int  k;

  for(k=1; k<maxk; k++) {
    MpIeee ak=  -MpIeee( "0.25" ) * (x/(nu+k)) * x/(nu+k+MpIeee( "1.0" ));
    rhok = -ak*(MpIeee( "1.0" ) + rhok)/(MpIeee( "1.0" ) + ak*(MpIeee( "1.0" ) + rhok));
    tk  *= rhok;
    sum += tk;

    dk = MpIeee( "1.0" ) / (MpIeee( "2.0" )/x - (nu+k-MpIeee( "1.0" ))/(nu+k) * dk);
    if(dk < MpIeee( "0.0" )) s = -s;

    if(fabs(tk/sum) < GSL_DBL_EPSILON) break;
  }

  *ratio = x/(MpIeee( "2.0" )*(nu+MpIeee( "1.0" ))) * sum;
  *sgn   = s;

  if(k == maxk)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}
#endif


/* Evaluate the continued fraction CF1 for I_{nu+1}/I_nu
 * using Gautschi (Euler) equivalent series.
 */
int
 gsl_sf_bessel_I_CF1_ser(const MpIeee nu, const MpIeee x, MpIeee * ratio)
{
  const int maxk = 20000;
  MpIeee tk=  MpIeee( "1.0" );
  MpIeee sum=  MpIeee( "1.0" );
  MpIeee rhok=  MpIeee( "0.0" );
  int  k;

  for(k=1; k<maxk; k++) {
    MpIeee ak=  MpIeee( "0.25" ) * (x/(nu+k)) * x/(nu+k+MpIeee( "1.0" ));
    rhok = -ak*(MpIeee( "1.0" ) + rhok)/(MpIeee( "1.0" ) + ak*(MpIeee( "1.0" ) + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < GSL_DBL_EPSILON) break;
  }

  *ratio = x/(MpIeee( "2.0" )*(nu+MpIeee( "1.0" ))) * sum;

  if(k == maxk)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


int
 gsl_sf_bessel_JY_steed_CF2(const MpIeee nu, const MpIeee x,
                           MpIeee * P, MpIeee * Q)
{
  const int max_iter = 10000;
  const MpIeee SMALL=  1.0e-100;

  int  i=  1;

  MpIeee x_inv=  MpIeee( "1.0" )/x;
  MpIeee a=  MpIeee( "0.25" ) - nu*nu;
  MpIeee p=  -MpIeee( "0.5" )*x_inv;
  MpIeee q=  MpIeee( "1.0" );
  MpIeee br=  MpIeee( "2.0" )*x;
  MpIeee bi=  MpIeee( "2.0" );
  MpIeee fact=  a*x_inv/(p*p + q*q);
  MpIeee cr=  br + q*fact;
  MpIeee ci=  bi + p*fact;
  MpIeee den=  br*br + bi*bi;
  MpIeee dr=  br/den;
  MpIeee di=  -bi/den;
  MpIeee dlr=  cr*dr - ci*di;
  MpIeee dli=  cr*di + ci*dr;
  MpIeee temp=  p*dlr - q*dli;
  q = p*dli + q*dlr;
  p = temp;
  for (i=2; i<=max_iter; i++) {
    a  += MpIeee( "2" )*(i-MpIeee( "1" ));
    bi += MpIeee( "2.0" );
    dr = a*dr + br;
    di = a*di + bi;
    if(fabs(dr)+fabs(di) < SMALL) dr = SMALL;
    fact = a/(cr*cr+ci*ci);
    cr = br + cr*fact;
    ci = bi - ci*fact;
    if(fabs(cr)+fabs(ci) < SMALL) cr = SMALL;
    den = dr*dr + di*di;
    dr /= den;
    di /= -den;
    dlr = cr*dr - ci*di;
    dli = cr*di + ci*dr;
    temp = p*dlr - q*dli;
    q = p*dli + q*dlr;
    p = temp;
    if(fabs(dlr-1.0)+fabs(dli) < GSL_DBL_EPSILON) break;
  }

  *P = p;
  *Q = q;

  if(i == max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 *
 * This is unstable for small x; x > 2 is a good cutoff.
 * Also requires |nu| < 1/2.
 */
int
 gsl_sf_bessel_K_scaled_steed_temme_CF2(const MpIeee nu, const MpIeee x,
                                       MpIeee * K_nu, MpIeee * K_nup1,
                                       MpIeee * Kp_nu)
{
  const int maxiter = 10000;

  int  i=  1;
  MpIeee bi=  MpIeee( "2.0" )*(MpIeee( "1.0" ) + x);
  MpIeee di=  MpIeee( "1.0" )/bi;
  MpIeee delhi=  di;
  MpIeee hi=  di;

  MpIeee qi=  MpIeee( "0.0" );
  MpIeee qip1=  MpIeee( "1.0" );

  MpIeee ai=  -(MpIeee( "0.25" ) - nu*nu);
  MpIeee a1=  ai;
  MpIeee ci=  -ai;
  MpIeee Qi=  -ai;

  MpIeee s=  MpIeee( "1.0" ) + Qi*delhi;

  for(i=2; i<=maxiter; i++) {
    MpIeee dels;
    MpIeee tmp;
    ai -= MpIeee( "2.0" )*(i-MpIeee( "1" ));
    ci  = -ai*ci/i;
    tmp  = (qi - bi*qip1)/ai;
    qi   = qip1;
    qip1 = tmp;
    Qi += ci*qip1;
    bi += MpIeee( "2.0" );
    di  = MpIeee( "1.0" )/(bi + ai*di);
    delhi = (bi*di - MpIeee( "1.0" )) * delhi;
    hi += delhi;
    dels = Qi*delhi;
    s += dels;
    if(fabs(dels/s) < GSL_DBL_EPSILON) break;
  }
  
  hi *= -a1;
  
  *K_nu   = sqrt(M_PI/(MpIeee( "2.0" )*x)) / s;
  *K_nup1 = *K_nu * (nu + x + MpIeee( "0.5" ) - hi)/x;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;
  if(i == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


int  gsl_sf_bessel_cos_pi4_e(MpIeee y, MpIeee eps, gsl_sf_result * result)
{
  const MpIeee sy=  sin(y);
  const MpIeee cy=  cos(y);
  const MpIeee s=  sy + cy;
  const MpIeee d=  sy - cy;
  const MpIeee abs_sum=  fabs(cy) + fabs(sy);
  MpIeee seps;
  MpIeee ceps;
  if(fabs(eps) < GSL_ROOT5_DBL_EPSILON) {
    const MpIeee e2=  eps*eps;
    seps = eps * (MpIeee( "1.0" ) - e2/MpIeee( "6.0" ) * (MpIeee( "1.0" ) - e2/MpIeee( "20.0" )));
    ceps = MpIeee( "1.0" ) - e2/MpIeee( "2.0" ) * (MpIeee( "1.0" ) - e2/MpIeee( "12.0" ));
  }
  else {
    seps = sin(eps);
    ceps = cos(eps);
  }
  result->val = (ceps * s - seps * d)/ M_SQRT2;
  result->err = 2.0 * GSL_DBL_EPSILON * (fabs(ceps) + fabs(seps)) * abs_sum / M_SQRT2;

  /* Try to account for error in evaluation of sin(y), cos(y).
   * This is a little sticky because we don't really know
   * how the library routines are doing their argument reduction.
   * However, we will make a reasonable guess.
   * FIXME ?
   */
  if(y > MpIeee( "1.0" )/GSL_DBL_EPSILON) {
    result->err *= 0.5 * y;
  }
  else if(y > MpIeee( "1.0" )/GSL_SQRT_DBL_EPSILON) {
    result->err *= 256.0 * y * GSL_SQRT_DBL_EPSILON;
  }

  return GSL_SUCCESS;
}


int  gsl_sf_bessel_sin_pi4_e(MpIeee y, MpIeee eps, gsl_sf_result * result)
{
  const MpIeee sy=  sin(y);
  const MpIeee cy=  cos(y);
  const MpIeee s=  sy + cy;
  const MpIeee d=  sy - cy;
  const MpIeee abs_sum=  fabs(cy) + fabs(sy);
  MpIeee seps;
  MpIeee ceps;
  if(fabs(eps) < GSL_ROOT5_DBL_EPSILON) {
    const MpIeee e2=  eps*eps;
    seps = eps * (MpIeee( "1.0" ) - e2/MpIeee( "6.0" ) * (MpIeee( "1.0" ) - e2/MpIeee( "20.0" )));
    ceps = MpIeee( "1.0" ) - e2/MpIeee( "2.0" ) * (MpIeee( "1.0" ) - e2/MpIeee( "12.0" ));
  }
  else {
    seps = sin(eps);
    ceps = cos(eps);
  }
  result->val = (ceps * d + seps * s)/ M_SQRT2;
  result->err = 2.0 * GSL_DBL_EPSILON * (fabs(ceps) + fabs(seps)) * abs_sum / M_SQRT2;

  /* Try to account for error in evaluation of sin(y), cos(y).
   * See above.
   * FIXME ?
   */
  if(y > MpIeee( "1.0" )/GSL_DBL_EPSILON) {
    result->err *= 0.5 * y;
  }
  else if(y > MpIeee( "1.0" )/GSL_SQRT_DBL_EPSILON) {
    result->err *= 256.0 * y * GSL_SQRT_DBL_EPSILON;
  }

  return GSL_SUCCESS;
}


/************************************************************************
 *                                                                      *
  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from
  G.N.Watson, A Treatise on the Theory of Bessel Functions,
  2nd Edition (Cambridge University Press, 1944).
  Higher terms in expansion for x near l given by
  Airey in Phil. Mag. 31, 520 (1916).

  This approximation is accurate to near 0.1% at the boundaries
  between the asymptotic regions; well away from the boundaries
  the accuracy is better than 10^{-5}.
 *                                                                      *
 ************************************************************************/
#if 0
MpIeee besselJ_meissel(MpIeee nu, MpIeee x)
{
  MpIeee beta=  pow(nu, MpIeee( "0.325" ));
  MpIeee result;

  /* Fitted matching points.   */
  MpIeee llimit=  MpIeee( "1.1" ) * beta;
  MpIeee ulimit=  MpIeee( "1.3" ) * beta;

  MpIeee nu2=  nu * nu;

  if (nu < MpIeee( "5." ) && x < MpIeee( "1." ))
    {
      /* Small argument and order. Use a Taylor expansion. */
      int  k;
      MpIeee xo2=  MpIeee( "0.5" ) * x;
      MpIeee gamfactor=  pow(nu,nu) * exp(-nu) * sqrt(nu * MpIeee( "2." ) * M_PI)
        * (MpIeee( "1." ) + MpIeee( "1." )/(MpIeee( "12." )*nu) + MpIeee( "1." )/(MpIeee( "288." )*nu*nu));
      MpIeee prefactor=  pow(xo2, nu) / gamfactor;
      MpIeee C[5];

      C[0] = MpIeee( "1." );
      C[1] = -C[0] / (nu+MpIeee( "1." ));
      C[2] = -C[1] / (MpIeee( "2." )*(nu+MpIeee( "2." )));
      C[3] = -C[2] / (MpIeee( "3." )*(nu+MpIeee( "3." )));
      C[4] = -C[3] / (MpIeee( "4." )*(nu+MpIeee( "4." )));
      
      result = MpIeee( "0." );
      for(k=0; k<5; k++)
        result += C[k] * pow(xo2, MpIeee( "2." )*k);

      result *= prefactor;
    }
  else if(x < nu - llimit)
    {
      /* Small x region: x << l.    */
      MpIeee z=  x / nu;
      MpIeee z2=  z*z;
      MpIeee rtomz2=  sqrt(MpIeee( "1." )-z2);
      MpIeee omz2_2=  (MpIeee( "1." )-z2)*(MpIeee( "1." )-z2);

      /* Calculate Meissel exponent. */
      MpIeee term1=  MpIeee( "1." )/(MpIeee( "24." )*nu) * ((MpIeee( "2." )+MpIeee( "3." )*z2)/((MpIeee( "1." )-z2)*rtomz2) -MpIeee( "2." ));
      MpIeee term2=  - z2*(MpIeee( "4." ) + z2)/(MpIeee( "16." )*nu2*(MpIeee( "1." )-z2)*omz2_2);
      MpIeee V_nu=  term1 + term2;
      
      /* Calculate the harmless prefactor. */
      MpIeee sterlingsum=  MpIeee( "1." ) + MpIeee( "1." )/(MpIeee( "12." )*nu) + MpIeee( "1." )/(MpIeee( "288" )*nu2);
      MpIeee harmless=  MpIeee( "1." ) / (sqrt(rtomz2*MpIeee( "2." )*M_PI*nu) * sterlingsum);

      /* Calculate the logarithm of the nu dependent prefactor. */
      MpIeee ln_nupre=  rtomz2 + log(z) - log(MpIeee( "1." ) + rtomz2);

      result = harmless * exp(nu*ln_nupre - V_nu);
    } 
  else if(x < nu + ulimit)
    {         
      /* Intermediate region 1: x near nu. */
      MpIeee eps=  MpIeee( "1." )-nu/x;
      MpIeee eps_x=  eps * x;
      MpIeee eps_x_2=  eps_x * eps_x;
      MpIeee xo6=  x/MpIeee( "6." );
      MpIeee B[6];
      static MpIeee gam[6] =  {MpIeee( "2.67894" ), MpIeee( "1.35412" ), MpIeee( "1." ), MpIeee( "0.89298" ), MpIeee( "0.902745" ), MpIeee( "1." )};
      static MpIeee sf[6] =  {MpIeee( "0.866025" ), MpIeee( "0.866025" ), MpIeee( "0." ), -MpIeee( "0.866025" ), -MpIeee( "0.866025" ), MpIeee( "0." )};
      
      /* Some terms are identically zero, because sf[] can be zero.
       * Some terms do not appear in the result.
       */
      B[0] = MpIeee( "1." );
      B[1] = eps_x;
      /* B[2] = 0.5 * eps_x_2 - 1./20.; */
      B[3] = eps_x * (eps_x_2/MpIeee( "6." ) - MpIeee( "1." )/MpIeee( "15." ));
      B[4] = eps_x_2 * (eps_x_2 - MpIeee( "1." ))/MpIeee( "24." ) + MpIeee( "1." )/MpIeee( "280." );
      /* B[5] = eps_x * (eps_x_2*(0.5*eps_x_2 - 1.)/60. + 43./8400.); */

      result  = B[0] * gam[0] * sf[0] / pow(xo6, MpIeee( "1." )/MpIeee( "3." ));
      result += B[1] * gam[1] * sf[1] / pow(xo6, MpIeee( "2." )/MpIeee( "3." ));
      result += B[3] * gam[3] * sf[3] / pow(xo6, MpIeee( "4." )/MpIeee( "3." ));
      result += B[4] * gam[4] * sf[4] / pow(xo6, MpIeee( "5." )/MpIeee( "3." ));

      result /= (MpIeee( "3." )*M_PI);
    }
  else 
    {
      /* Region of very large argument. Use expansion
       * for x>>l, and we need not be very exacting.
       */
      MpIeee secb=  x/nu;
      MpIeee sec2b=  secb*secb;
      
      MpIeee cotb=  MpIeee( "1." )/sqrt(sec2b-MpIeee( "1." ));      /* cotb=cot(beta) */

      MpIeee beta=  acos(nu/x);
      MpIeee trigarg=  nu/cotb - nu*beta - MpIeee( "0.25" ) * M_PI;
      
      MpIeee cot3b=  cotb * cotb * cotb;
      MpIeee cot6b=  cot3b * cot3b;

      MpIeee sum1;MpIeee  sum2;MpIeee  expterm;MpIeee  prefactor;MpIeee  trigcos;

      sum1  = MpIeee( "2.0" ) + MpIeee( "3.0" ) * sec2b;
      trigarg -= sum1 * cot3b / (MpIeee( "24.0" ) * nu);

      trigcos = cos(trigarg);

      sum2 = MpIeee( "4.0" ) + sec2b;
      expterm = sum2 * sec2b * cot6b / (MpIeee( "16.0" ) * nu2);

      expterm = exp(-expterm);
      prefactor = sqrt(MpIeee( "2." ) * cotb / (nu * M_PI));
      
      result = prefactor * expterm * trigcos;
    }

  return  result;
}
#endif
