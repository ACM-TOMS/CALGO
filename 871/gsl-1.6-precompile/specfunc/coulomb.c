#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/coulomb.c
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

/* Evaluation of Coulomb wave functions F_L(eta, x), G_L(eta, x),
 * and their derivatives. A combination of Steed's method, asymptotic
 * results, and power series.
 *
 * Steed's method:
 *  [Barnett, CPC 21, 297 (1981)]
 * Power series and other methods:
 *  [Biedenharn et al., PR 97, 542 (1954)]
 *  [Bardin et al., CPC 3, 73 (1972)]
 *  [Abad+Sesma, CPC 71, 110 (1992)]
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coulomb.h>

#include "error.h"

/* the L=0 normalization constant
 * [Abramowitz+Stegun 14.1.8]
 */
static
MpIeee C0sq(MpIeee eta)
{
  MpIeee twopieta=  MpIeee( "2.0" )*M_PI*eta;

  if(fabs(eta) < GSL_DBL_EPSILON) {
    return MpIeee( "1.0" );
  }
  else if(twopieta > GSL_LOG_DBL_MAX) {
    return MpIeee( "0.0" );
  }
  else {
    gsl_sf_result scale;
    gsl_sf_expm1_e(twopieta, &scale);
    return twopieta/scale.val;
  }
}


/* the full definition of C_L(eta) for any valid L and eta
 * [Abramowitz and Stegun 14.1.7]
 * This depends on the complex gamma function. For large
 * arguments the phase of the complex gamma function is not
 * very accurately determined. However the modulus is, and that
 * is all that we need to calculate C_L.
 *
 * This is not valid for L <= -3/2  or  L = -1.
 */
static
int
 CLeta(MpIeee L, MpIeee eta, gsl_sf_result * result)
{
  gsl_sf_result ln1; /* log of numerator Gamma function */
  gsl_sf_result ln2; /* log of denominator Gamma function */
  MpIeee sgn=  MpIeee( "1.0" );
  MpIeee arg_val;MpIeee  arg_err;

  if(fabs(eta/(L+1.0)) < GSL_DBL_EPSILON) {
    gsl_sf_lngamma_e(L+1.0, &ln1);
  }
  else {
    gsl_sf_result p1;                 /* phase of numerator Gamma -- not used */
    gsl_sf_lngamma_complex_e(L+1.0, eta, &ln1, &p1); /* should be ok */
  }

  gsl_sf_lngamma_e(2.0*(L+1.0), &ln2);
  if(L < -MpIeee( "1.0" )) sgn = -sgn;

  arg_val  = L*M_LN2 - MpIeee( "0.5" )*eta*M_PI + ln1.val - ln2.val;
  arg_err  = ln1.err + ln2.err;
  arg_err += GSL_DBL_EPSILON * (fabs(L*M_LN2) + fabs(MpIeee( "0.5" )*eta*M_PI));
  return gsl_sf_exp_err_e(arg_val, arg_err, result);
}


int
 gsl_sf_coulomb_CL_e(MpIeee lam, MpIeee eta, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(lam <= -MpIeee( "1.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(lam) < GSL_DBL_EPSILON) {
    /* saves a calculation of complex_lngamma(), otherwise not necessary */
    result->val = sqrt(C0sq(eta));
    result->err = 2.0 * GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    return CLeta(lam, eta, result);
  }
}


/* cl[0] .. cl[kmax] = C_{lam_min}(eta) .. C_{lam_min+kmax}(eta)
 */
int
 gsl_sf_coulomb_CL_array(MpIeee lam_min, int  kmax, MpIeee eta, MpIeee * cl)
{
  int  k;
  gsl_sf_result cl_0;
  gsl_sf_coulomb_CL_e(lam_min, eta, &cl_0);
  cl[0] = cl_0.val;

  for(k=1; k<=kmax; k++) {
    MpIeee L=  lam_min + k;
    cl[k] = cl[k-1] * sqrt(L*L + eta*eta)/(L*(MpIeee( "2.0" )*L+MpIeee( "1.0" )));
  }

  return GSL_SUCCESS;
}


/* Evaluate the series for Phi_L(eta,x) and Phi_L*(eta,x)
 * [Abramowitz+Stegun 14.1.5]
 * [Abramowitz+Stegun 14.1.13]
 *
 * The sequence of coefficients A_k^L is
 * manifestly well-controlled for L >= -1/2
 * and eta < 10.
 *
 * This makes sense since this is the region
 * away from threshold, and you expect
 * the evaluation to become easier as you
 * get farther from threshold.
 *
 * Empirically, this is quite well-behaved for
 *   L >= -1/2
 *   eta < 10
 *   x   < 10
 */
#if 0
static
int
 coulomb_Phi_series(const MpIeee lam, const MpIeee eta, const MpIeee x,
                   MpIeee * result, MpIeee * result_star)
{
  int  kmin=    5;
  int  kmax=  200;
  int  k;
  MpIeee Akm2=  MpIeee( "1.0" );
  MpIeee Akm1=  eta/(lam+MpIeee( "1.0" ));
  MpIeee Ak;

  MpIeee xpow=  x;
  MpIeee sum=  Akm2 + Akm1*x;
  MpIeee sump=  (lam+MpIeee( "1.0" ))*Akm2 + (lam+MpIeee( "2.0" ))*Akm1*x;
  MpIeee prev_abs_del=  fabs(Akm1*x);
  MpIeee prev_abs_del_p=  (lam+MpIeee( "2.0" )) * prev_abs_del;

  for(k=2; k<kmax; k++) {
    MpIeee del;
    MpIeee del_p;
    MpIeee abs_del;
    MpIeee abs_del_p;

    Ak = (MpIeee( "2.0" )*eta*Akm1 - Akm2)/(k*(MpIeee( "2.0" )*lam + MpIeee( "1.0" ) + k));

    xpow *= x;
    del   = Ak*xpow;
    del_p = (k+lam+MpIeee( "1.0" ))*del;
    sum  += del;
    sump += del_p;

    abs_del   = fabs(del);
    abs_del_p = fabs(del_p);

    if(          abs_del/(fabs(sum)+abs_del)          < GSL_DBL_EPSILON
       &&   prev_abs_del/(fabs(sum)+prev_abs_del)     < GSL_DBL_EPSILON
       &&      abs_del_p/(fabs(sump)+abs_del_p)       < GSL_DBL_EPSILON
       && prev_abs_del_p/(fabs(sump)+prev_abs_del_p)  < GSL_DBL_EPSILON
       && k > kmin
       ) break;

    /* We need to keep track of the previous delta because when
     * eta is near zero the odd terms of the sum are very small
     * and this could lead to premature termination.
     */
    prev_abs_del   = abs_del;
    prev_abs_del_p = abs_del_p;

    Akm2 = Akm1;
    Akm1 = Ak;
  }

  *result      = sum;
  *result_star = sump;

  if(k==kmax) {
    GSL_ERROR ("error", GSL_EMAXITER);
  }
  else {
    return GSL_SUCCESS;
  }
}
#endif /* 0 */


/* Determine the connection phase, phi_lambda.
 * See coulomb_FG_series() below. We have
 * to be careful about sin(phi)->0. Note that
 * there is an underflow condition for large 
 * positive eta in any case.
 */
static
int
 coulomb_connection(const MpIeee lam, const MpIeee eta,
                   MpIeee * cos_phi, MpIeee * sin_phi)
{
  if(eta > -GSL_LOG_DBL_MIN/2.0*M_PI-1.0) {
    *cos_phi = MpIeee( "1.0" );
    *sin_phi = MpIeee( "0.0" );
    GSL_ERROR ("error", GSL_EUNDRFLW);
  }
  else if(eta > -GSL_LOG_DBL_EPSILON/(4.0*M_PI)) {
    const MpIeee eps=  2.0 * exp(-2.0*M_PI*eta);
    const MpIeee tpl=  tan(M_PI * lam);
    const MpIeee dth=  eps * tpl / (tpl*tpl + 1.0);
    *cos_phi = -MpIeee( "1.0" ) + MpIeee( "0.5" ) * dth*dth;
    *sin_phi = -dth;
    return GSL_SUCCESS;
  }
  else {
    MpIeee X=  tanh(M_PI * eta) / tan(M_PI * lam);
    MpIeee phi=  -atan(X) - (lam + MpIeee( "0.5" )) * M_PI;
    *cos_phi = cos(phi);
    *sin_phi = sin(phi);
    return GSL_SUCCESS;
  }
}


/* Evaluate the Frobenius series for F_lam(eta,x) and G_lam(eta,x).
 * Homegrown algebra. Evaluates the series for F_{lam} and
 * F_{-lam-1}, then uses
 *    G_{lam} = (F_{lam} cos(phi) - F_{-lam-1}) / sin(phi)
 * where
 *    phi = Arg[Gamma[1+lam+I eta]] - Arg[Gamma[-lam + I eta]] - (lam+1/2)Pi
 *        = Arg[Sin[Pi(-lam+I eta)] - (lam+1/2)Pi
 *        = atan2(-cos(lam Pi)sinh(eta Pi), -sin(lam Pi)cosh(eta Pi)) - (lam+1/2)Pi
 *
 *        = -atan(X) - (lam+1/2) Pi,  X = tanh(eta Pi)/tan(lam Pi)
 *
 * Not appropriate for lam <= -1/2, lam = 0, or lam >= 1/2.
 */
static
int
 coulomb_FG_series(const MpIeee lam, const MpIeee eta, const MpIeee x,
                  gsl_sf_result * F, gsl_sf_result * G)
{
  const int max_iter = 800;
  gsl_sf_result ClamA;
  gsl_sf_result ClamB;
  int  stat_A=  CLeta(lam, eta, &ClamA);
  int  stat_B=  CLeta(-lam-1.0, eta, &ClamB);
  const MpIeee tlp1=  2.0*lam + 1.0;
  const MpIeee pow_x=  pow(x, lam);
  MpIeee cos_phi_lam;
  MpIeee sin_phi_lam;

  MpIeee uA_mm2=  MpIeee( "1.0" );                  /* uA sum is for F_{lam} */
  MpIeee uA_mm1=  x*eta/(lam+MpIeee( "1.0" ));
  MpIeee uA_m;
  MpIeee uB_mm2=  MpIeee( "1.0" );                  /* uB sum is for F_{-lam-1} */
  MpIeee uB_mm1=  -x*eta/lam;
  MpIeee uB_m;
  MpIeee A_sum=  uA_mm2 + uA_mm1;
  MpIeee B_sum=  uB_mm2 + uB_mm1;
  MpIeee A_abs_del_prev=  fabs(A_sum);
  MpIeee B_abs_del_prev=  fabs(B_sum);
  gsl_sf_result FA, FB;
  int  m=  2;

  int  stat_conn=  coulomb_connection(lam, eta, &cos_phi_lam, &sin_phi_lam);

  if(stat_conn == GSL_EUNDRFLW) {
    F->val = 0.0;  /* FIXME: should this be set to Inf too like G? */
    F->err = 0.0;
    OVERFLOW_ERROR(G);
  }

  while(m < max_iter) {
    MpIeee abs_dA;
    MpIeee abs_dB;
    uA_m = x*(MpIeee( "2.0" )*eta*uA_mm1 - x*uA_mm2)/(m*(m+tlp1));
    uB_m = x*(MpIeee( "2.0" )*eta*uB_mm1 - x*uB_mm2)/(m*(m-tlp1));
    A_sum += uA_m;
    B_sum += uB_m;
    abs_dA = fabs(uA_m);
    abs_dB = fabs(uB_m);
    if(m > 15) {
      /* Don't bother checking until we have gone out a little ways;
       * a minor optimization. Also make sure to check both the
       * current and the previous increment because the odd and even
       * terms of the sum can have very different behaviour, depending
       * on the value of eta.
       */
      MpIeee max_abs_dA=  GSL_MAX(abs_dA, A_abs_del_prev);
      MpIeee max_abs_dB=  GSL_MAX(abs_dB, B_abs_del_prev);
      MpIeee abs_A=  fabs(A_sum);
      MpIeee abs_B=  fabs(B_sum);
      if(   max_abs_dA/(max_abs_dA + abs_A) < 4.0*GSL_DBL_EPSILON
         && max_abs_dB/(max_abs_dB + abs_B) < 4.0*GSL_DBL_EPSILON
         ) break;
    }
    A_abs_del_prev = abs_dA;
    B_abs_del_prev = abs_dB;
    uA_mm2 = uA_mm1;
    uA_mm1 = uA_m;
    uB_mm2 = uB_mm1;
    uB_mm1 = uB_m;
    m++;
  }

  FA.val = A_sum * ClamA.val * pow_x * x;
  FA.err = fabs(A_sum) * ClamA.err * pow_x * x + 2.0*GSL_DBL_EPSILON*fabs(FA.val);
  FB.val = B_sum * ClamB.val / pow_x;
  FB.err = fabs(B_sum) * ClamB.err / pow_x + 2.0*GSL_DBL_EPSILON*fabs(FB.val);

  F->val = FA.val;
  F->err = FA.err;

  G->val = (FA.val * cos_phi_lam - FB.val)/sin_phi_lam;
  G->err = (FA.err * fabs(cos_phi_lam) + FB.err)/fabs(sin_phi_lam);

  if(m >= max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_ERROR_SELECT_2(stat_A, stat_B);
}


/* Evaluate the Frobenius series for F_0(eta,x) and G_0(eta,x).
 * See [Bardin et al., CPC 3, 73 (1972), (14)-(17)];
 * note the misprint in (17): nu_0=1 is correct, not nu_0=0.
 */
static
int
 coulomb_FG0_series(const MpIeee eta, const MpIeee x,
                   gsl_sf_result * F, gsl_sf_result * G)
{
  const int max_iter = 800;
  const MpIeee x2=  x*x;
  const MpIeee tex=  2.0*eta*x;
  gsl_sf_result C0;
  int  stat_CL=  CLeta(0.0, eta, &C0);
  gsl_sf_result r1pie;
  int  psi_stat=  gsl_sf_psi_1piy_e(eta, &r1pie);
  MpIeee u_mm2=  MpIeee( "0.0" );  /* u_0 */
  MpIeee u_mm1=  x;    /* u_1 */
  MpIeee u_m;
  MpIeee v_mm2=  MpIeee( "1.0" );                               /* nu_0 */
  MpIeee v_mm1=  tex*(MpIeee( "2.0" )*M_EULER-MpIeee( "1.0" )+r1pie.val);   /* nu_1 */
  MpIeee v_m;
  MpIeee u_sum=  u_mm2 + u_mm1;
  MpIeee v_sum=  v_mm2 + v_mm1;
  MpIeee u_abs_del_prev=  fabs(u_sum);
  MpIeee v_abs_del_prev=  fabs(v_sum);
  int  m=  2;
  MpIeee u_sum_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(u_sum);
  MpIeee v_sum_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(v_sum);
  MpIeee ln2x=  log(MpIeee( "2.0" )*x);

  while(m < max_iter) {
    MpIeee abs_du;
    MpIeee abs_dv;
    MpIeee m_mm1=  m*(m-MpIeee( "1.0" ));
    u_m = (tex*u_mm1 - x2*u_mm2)/m_mm1;
    v_m = (tex*v_mm1 - x2*v_mm2 - MpIeee( "2.0" )*eta*(MpIeee( "2" )*m-MpIeee( "1" ))*u_m)/m_mm1;
    u_sum += u_m;
    v_sum += v_m;
    abs_du = fabs(u_m);
    abs_dv = fabs(v_m);
    u_sum_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * abs_du;
    v_sum_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * abs_dv;
    if(m > 15) {
      /* Don't bother checking until we have gone out a little ways;
       * a minor optimization. Also make sure to check both the
       * current and the previous increment because the odd and even
       * terms of the sum can have very different behaviour, depending
       * on the value of eta.
       */
      MpIeee max_abs_du=  GSL_MAX(abs_du, u_abs_del_prev);
      MpIeee max_abs_dv=  GSL_MAX(abs_dv, v_abs_del_prev);
      MpIeee abs_u=  fabs(u_sum);
      MpIeee abs_v=  fabs(v_sum);
      if(   max_abs_du/(max_abs_du + abs_u) < 40.0*GSL_DBL_EPSILON
         && max_abs_dv/(max_abs_dv + abs_v) < 40.0*GSL_DBL_EPSILON
         ) break;
    }
    u_abs_del_prev = abs_du;
    v_abs_del_prev = abs_dv;
    u_mm2 = u_mm1;
    u_mm1 = u_m;
    v_mm2 = v_mm1;
    v_mm1 = v_m;
    m++;
  }

  F->val  = C0.val * u_sum;
  F->err  = C0.err * fabs(u_sum);
  F->err += fabs(C0.val) * u_sum_err;
  F->err += 2.0 * GSL_DBL_EPSILON * fabs(F->val);

  G->val  = (v_sum + 2.0*eta*u_sum * ln2x) / C0.val;
  G->err  = (fabs(v_sum) + fabs(2.0*eta*u_sum * ln2x)) / fabs(C0.val) * fabs(C0.err/C0.val);
  G->err += (v_sum_err + fabs(2.0*eta*u_sum_err*ln2x)) / fabs(C0.val);
  G->err += 2.0 * GSL_DBL_EPSILON * fabs(G->val);

  if(m == max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_ERROR_SELECT_2(psi_stat, stat_CL);
}


/* Evaluate the Frobenius series for F_{-1/2}(eta,x) and G_{-1/2}(eta,x).
 * Homegrown algebra.
 */
static
int
 coulomb_FGmhalf_series(const MpIeee eta, const MpIeee x,
                       gsl_sf_result * F, gsl_sf_result * G)
{
  const int max_iter = 800;
  const MpIeee rx=  sqrt(x);
  const MpIeee x2=  x*x;
  const MpIeee tex=  2.0*eta*x;
  gsl_sf_result Cmhalf;
  int  stat_CL=  CLeta(-0.5, eta, &Cmhalf);
  MpIeee u_mm2=  MpIeee( "1.0" );                      /* u_0 */
  MpIeee u_mm1=  tex * u_mm2;              /* u_1 */
  MpIeee u_m;
  MpIeee v_mm2;MpIeee  v_mm1;MpIeee  v_m;
  MpIeee f_sum;MpIeee  g_sum;
  MpIeee tmp1;
  gsl_sf_result rpsi_1pe;
  gsl_sf_result rpsi_1p2e;
  int  m=  2;

  gsl_sf_psi_1piy_e(eta,     &rpsi_1pe);
  gsl_sf_psi_1piy_e(2.0*eta, &rpsi_1p2e);

  v_mm2 = MpIeee( "2.0" )*M_EULER - M_LN2 - rpsi_1pe.val + MpIeee( "2.0" )*rpsi_1p2e.val;
  v_mm1 = tex*(v_mm2 - MpIeee( "2.0" )*u_mm2);

  f_sum = u_mm2 + u_mm1;
  g_sum = v_mm2 + v_mm1;

  while(m < max_iter) {
    MpIeee m2=  m*m;
    u_m = (tex*u_mm1 - x2*u_mm2)/m2;
    v_m = (tex*v_mm1 - x2*v_mm2 - MpIeee( "2.0" )*m*u_m)/m2;
    f_sum += u_m;
    g_sum += v_m;
    if(   f_sum != MpIeee( "0.0" )
       && g_sum != MpIeee( "0.0" )
       && (fabs(u_m/f_sum) + fabs(v_m/g_sum) < MpIeee( "10.0" )*GSL_DBL_EPSILON)) break;
    u_mm2 = u_mm1;
    u_mm1 = u_m;
    v_mm2 = v_mm1;
    v_mm1 = v_m;
    m++;
  }
  
  F->val = Cmhalf.val * rx * f_sum;
  F->err = Cmhalf.err * fabs(rx * f_sum) + 2.0*GSL_DBL_EPSILON*fabs(F->val);

  tmp1 = f_sum*log(x);
  G->val = -rx*(tmp1 + g_sum)/Cmhalf.val;
  G->err = fabs(rx)*(fabs(tmp1) + fabs(g_sum))/fabs(Cmhalf.val) * fabs(Cmhalf.err/Cmhalf.val);

  if(m == max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return stat_CL;
}


/* Evolve the backwards recurrence for F,F'.
 *
 *    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
 *    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
 * where
 *    R_lam = sqrt(1 + (eta/lam)^2)
 *    S_lam = lam/x + eta/lam
 *
 */
static
int
 coulomb_F_recur(MpIeee lam_min, int  kmax,
                MpIeee eta, MpIeee x,
                MpIeee F_lam_max, MpIeee Fp_lam_max,
                MpIeee * F_lam_min, MpIeee * Fp_lam_min)
{
  MpIeee x_inv=  MpIeee( "1.0" )/x;
  MpIeee fcl=  F_lam_max;
  MpIeee fpl=  Fp_lam_max;
  MpIeee lam_max=  lam_min + kmax;
  MpIeee lam=  lam_max;
  int  k;

  for(k=kmax-1; k>=0; k--) {
    MpIeee el=  eta/lam;
    MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
    MpIeee sl=  el  + lam*x_inv;
    MpIeee fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fpl    =  fc_lm1*sl - fcl*rl;
    fcl    =  fc_lm1;
    lam -= MpIeee( "1.0" );
  }

  *F_lam_min  = fcl;
  *Fp_lam_min = fpl;  
  return GSL_SUCCESS;
}


/* Evolve the forward recurrence for G,G'.
 *
 *   G_{lam+1}  = (S_lam G_lam - G_lam')/R_lam
 *   G_{lam+1}' = R_{lam+1} G_lam - S_lam G_{lam+1}
 *
 * where S_lam and R_lam are as above in the F recursion.
 */
static
int
 coulomb_G_recur(const MpIeee lam_min, const int kmax,
                const MpIeee eta, const MpIeee x,
                const MpIeee G_lam_min, const MpIeee Gp_lam_min,
                MpIeee * G_lam_max, MpIeee * Gp_lam_max)
{
  MpIeee x_inv=  MpIeee( "1.0" )/x;
  MpIeee gcl=  G_lam_min;
  MpIeee gpl=  Gp_lam_min;
  MpIeee lam=  lam_min + MpIeee( "1.0" );
  int  k;

  for(k=1; k<=kmax; k++) {
    MpIeee el=  eta/lam;
    MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
    MpIeee sl=  el + lam*x_inv;
    MpIeee gcl1=  (sl*gcl - gpl)/rl;
    gpl   = rl*gcl - sl*gcl1;
    gcl   = gcl1;
    lam += MpIeee( "1.0" );
  }
  
  *G_lam_max  = gcl;
  *Gp_lam_max = gpl;
  return GSL_SUCCESS;
}


/* Evaluate the first continued fraction, giving
 * the ratio F'/F at the upper lambda value.
 * We also determine the sign of F at that point,
 * since it is the sign of the last denominator
 * in the continued fraction.
 */
static
int
 coulomb_CF1(MpIeee lambda,
            MpIeee eta, MpIeee x,
            MpIeee * fcl_sign,
            MpIeee * result,
            int  * count)
{
  const MpIeee CF1_small=  1.e-30;
  const MpIeee CF1_abort=  1.0e+05;
  const MpIeee CF1_acc=  2.0*GSL_DBL_EPSILON;
  const MpIeee x_inv=  1.0/x;
  const MpIeee px=  lambda + 1.0 + CF1_abort;

  MpIeee pk=  lambda + MpIeee( "1.0" );
  MpIeee F=  eta/pk + pk*x_inv;
  MpIeee D;MpIeee  C;
  MpIeee df;

  *fcl_sign = MpIeee( "1.0" );
  *count = 0;

  if(fabs(F) < CF1_small) F = CF1_small;
  D = MpIeee( "0.0" );
  C = F;

  do {
    MpIeee pk1=  pk + MpIeee( "1.0" );
    MpIeee ek=  eta / pk;
    MpIeee rk2=  MpIeee( "1.0" ) + ek*ek;
    MpIeee tk=  (pk + pk1)*(x_inv + ek/pk1);
    D   =  tk - rk2 * D;
    C   =  tk - rk2 / C;
    if(fabs(C) < CF1_small) C = CF1_small;
    if(fabs(D) < CF1_small) D = CF1_small;
    D = MpIeee( "1.0" )/D;
    df = D * C;
    F  = F * df;
    if(D < MpIeee( "0.0" )) {
      /* sign of result depends on sign of denominator */
      *fcl_sign = - *fcl_sign;
    }
    pk = pk1;
    if( pk > px ) {
      *result = F;
      GSL_ERROR ("error", GSL_ERUNAWAY);
    }
    ++(*count);
  }
  while(fabs(df-1.0) > CF1_acc);
  
  *result = F;
  return GSL_SUCCESS;
}


#if 0
static
int
 old_coulomb_CF1(const MpIeee lambda,
                MpIeee eta, MpIeee x,
                MpIeee * fcl_sign,
                MpIeee * result)
{
  const MpIeee CF1_abort=  1.e5;
  const MpIeee CF1_acc=  10.0*GSL_DBL_EPSILON;
  const MpIeee x_inv=  1.0/x;
  const MpIeee px=  lambda + 1.0 + CF1_abort;
  
  MpIeee pk=  lambda + MpIeee( "1.0" );
  
  MpIeee D;
  MpIeee df;

  MpIeee F;
  MpIeee p;
  MpIeee pk1;
  MpIeee ek;
  
  MpIeee fcl=  MpIeee( "1.0" );

  MpIeee tk;

  while(1) {
    ek = eta/pk;
    F = (ek + pk*x_inv)*fcl + (fcl - MpIeee( "1.0" ))*x_inv;
    pk1 = pk + MpIeee( "1.0" );
    if(fabs(eta*x + pk*pk1) > CF1_acc) break;
    fcl = (MpIeee( "1.0" ) + ek*ek)/(MpIeee( "1.0" ) + eta*eta/(pk1*pk1));
    pk = MpIeee( "2.0" ) + pk;
  }

  D  = MpIeee( "1.0" )/((pk + pk1)*(x_inv + ek/pk1));
  df = -fcl*(MpIeee( "1.0" ) + ek*ek)*D;
  
  if(fcl != MpIeee( "1.0" )) fcl = -MpIeee( "1.0" );
  if(D    < MpIeee( "0.0" )) fcl = -fcl;
  
  F = F + df;

  p = MpIeee( "1.0" );
  do {
    pk = pk1;
    pk1 = pk + MpIeee( "1.0" );
    ek  = eta / pk;
    tk  = (pk + pk1)*(x_inv + ek/pk1);
    D   =  tk - D*(MpIeee( "1.0" )+ek*ek);
    if(fabs(D) < sqrt(CF1_acc)) {
      p += MpIeee( "1.0" );
      if(p > MpIeee( "2.0" )) {
        {cout<<"HELP............\n";}
      }
    }
    D = MpIeee( "1.0" )/D;
    if(D < MpIeee( "0.0" )) {
      /* sign of result depends on sign of denominator */
      fcl = -fcl;
    }
    df = df*(D*tk - MpIeee( "1.0" ));
    F  = F + df;
    if( pk > px ) {
      GSL_ERROR ("error", GSL_ERUNAWAY);
    }
  }
  while(fabs(df) > fabs(F)*CF1_acc);
  
  *fcl_sign = fcl;
  *result = F;
  return GSL_SUCCESS;
}
#endif /* 0 */


/* Evaluate the second continued fraction to 
 * obtain the ratio
 *    (G' + i F')/(G + i F) := P + i Q
 * at the specified lambda value.
 */
static
int
 coulomb_CF2(const MpIeee lambda, const MpIeee eta, const MpIeee x,
            MpIeee * result_P, MpIeee * result_Q, int  * count)
{
  int  status=  GSL_SUCCESS;

  const MpIeee CF2_acc=  4.0*GSL_DBL_EPSILON;
  const MpIeee CF2_abort=  2.0e+05;

  const MpIeee wi=  2.0*eta;
  const MpIeee x_inv=  1.0/x;
  const MpIeee e2mm1=  eta*eta + lambda*(lambda + 1.0);
  
  MpIeee ar=  -e2mm1;
  MpIeee ai=   eta;

  MpIeee br=   MpIeee( "2.0" )*(x - eta);
  MpIeee bi=   MpIeee( "2.0" );

  MpIeee dr=   br/(br*br + bi*bi);
  MpIeee di=  -bi/(br*br + bi*bi);

  MpIeee dp=  -x_inv*(ar*di + ai*dr);
  MpIeee dq=   x_inv*(ar*dr - ai*di);

  MpIeee A;MpIeee  B;MpIeee  C;MpIeee  D;

  MpIeee pk=   MpIeee( "0.0" );
  MpIeee P=   MpIeee( "0.0" );
  MpIeee Q=   MpIeee( "1.0" ) - eta*x_inv;

  *count = 0;
 
  do {
    P += dp;
    Q += dq;
    pk += MpIeee( "2.0" );
    ar += pk;
    ai += wi;
    bi += MpIeee( "2.0" );
    D  = ar*dr - ai*di + br;
    di = ai*dr + ar*di + bi;
    C  = MpIeee( "1.0" )/(D*D + di*di);
    dr =  C*D;
    di = -C*di;
    A  = br*dr - bi*di - MpIeee( "1." );
    B  = bi*dr + br*di;
    C  = dp*A  - dq*B;
    dq = dp*B  + dq*A;
    dp = C;
    if(pk > CF2_abort) {
      status = GSL_ERUNAWAY;
      break;
    }
    ++(*count);
  }
  while(fabs(dp)+fabs(dq) > (fabs(P)+fabs(Q))*CF2_acc);

  if(Q < CF2_abort*GSL_DBL_EPSILON*fabs(P)) {
    status = GSL_ELOSS;
  }

  *result_P = P;
  *result_Q = Q;
  return status;
}


/* WKB evaluation of F, G. Assumes  0 < x < turning point.
 * Overflows are trapped, GSL_EOVRFLW is signalled,
 * and an exponent is returned such that:
 *
 *   result_F = fjwkb * exp(-exponent)
 *   result_G = gjwkb * exp( exponent)
 *
 * See [Biedenharn et al. Phys. Rev. 97, 542-554 (1955), Section IV]
 *
 * Unfortunately, this is not very accurate in general. The
 * test cases typically have 3-4 digits of precision. One could
 * argue that this is ok for general use because, for instance,
 * F is exponentially small in this region and so the absolute
 * accuracy is still roughly acceptable. But it would be better
 * to have a systematic method for improving the precision. See
 * the Abad+Sesma method discussion below.
 */
static
int
 coulomb_jwkb(const MpIeee lam, const MpIeee eta, const MpIeee x,
             gsl_sf_result * fjwkb, gsl_sf_result * gjwkb,
             MpIeee * exponent)
{
  const MpIeee llp1=  lam*(lam+1.0) + 6.0/35.0;
  const MpIeee llp1_eff=  GSL_MAX(llp1, 0.0);
  const MpIeee rho_ghalf=  sqrt(x*(2.0*eta - x) + llp1_eff);
  const MpIeee sinh_arg=  sqrt(llp1_eff/(eta*eta+llp1_eff)) * rho_ghalf / x;
  const MpIeee sinh_inv=  log(sinh_arg + sqrt(1.0 + sinh_arg*sinh_arg));

  const MpIeee phi=  fabs(rho_ghalf - eta*atan2(rho_ghalf,x-eta) - sqrt(llp1_eff) * sinh_inv);

  const MpIeee zeta_half=  pow(3.0*phi/2.0, 1.0/3.0);
  const MpIeee prefactor=  sqrt(M_PI*phi*x/(6.0 * rho_ghalf));
  
  MpIeee F=  prefactor * MpIeee( "3.0" )/zeta_half;
  MpIeee G=  prefactor * MpIeee( "3.0" )/zeta_half; /* Note the sqrt(3) from Bi normalization */
  MpIeee F_exp;
  MpIeee G_exp;
  
  const MpIeee airy_scale_exp=  phi;
  gsl_sf_result ai;
  gsl_sf_result bi;
  gsl_sf_airy_Ai_scaled_e(zeta_half*zeta_half, GSL_MODE_DEFAULT, &ai);
  gsl_sf_airy_Bi_scaled_e(zeta_half*zeta_half, GSL_MODE_DEFAULT, &bi);
  F *= ai.val;
  G *= bi.val;
  F_exp = log(F) - airy_scale_exp;
  G_exp = log(G) + airy_scale_exp;

  if(G_exp >= GSL_LOG_DBL_MAX) {
    fjwkb->val = F;
    gjwkb->val = G;
    fjwkb->err = 1.0e-3 * fabs(F); /* FIXME: real error here ... could be smaller */
    gjwkb->err = 1.0e-3 * fabs(G);
    *exponent = airy_scale_exp;
    GSL_ERROR ("error", GSL_EOVRFLW);
  }
  else {
    fjwkb->val = exp(F_exp);
    gjwkb->val = exp(G_exp);
    fjwkb->err = 1.0e-3 * fabs(fjwkb->val);
    gjwkb->err = 1.0e-3 * fabs(gjwkb->val);
    *exponent = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
}


/* Asymptotic evaluation of F and G below the minimal turning point.
 *
 * This is meant to be a drop-in replacement for coulomb_jwkb().
 * It uses the expressions in [Abad+Sesma]. This requires some
 * work because I am not sure where it is valid. They mumble
 * something about |x| < |lam|^(-1/2) or 8|eta x| > lam when |x| < 1.
 * This seems true, but I thought the result was based on a uniform
 * expansion and could be controlled by simply using more terms.
 */
#if 0
static
int
 coulomb_AS_xlt2eta(const MpIeee lam, const MpIeee eta, const MpIeee x,
                   gsl_sf_result * f_AS, gsl_sf_result * g_AS,
                   MpIeee * exponent)
{
  /* no time to do this now... */
}
#endif /* 0 */



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_coulomb_wave_FG_e(const MpIeee eta, const MpIeee x,
                            const MpIeee lam_F,
                            const int  k_lam_G,      /* lam_G = lam_F - k_lam_G */
                            gsl_sf_result * F, gsl_sf_result * Fp,
                            gsl_sf_result * G, gsl_sf_result * Gp,
                            MpIeee * exp_F, MpIeee * exp_G)
{
  const MpIeee lam_G=  lam_F - k_lam_G;

  if(x < 0.0 || lam_F <= -0.5 || lam_G <= -0.5) {
    GSL_SF_RESULT_SET(F,  0.0, 0.0);
    GSL_SF_RESULT_SET(Fp, 0.0, 0.0);
    GSL_SF_RESULT_SET(G,  0.0, 0.0);
    GSL_SF_RESULT_SET(Gp, 0.0, 0.0);
    *exp_F = MpIeee( "0.0" );
    *exp_G = MpIeee( "0.0" );
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(x == 0.0) {
    gsl_sf_result C0;
    CLeta(0.0, eta, &C0);
    GSL_SF_RESULT_SET(F,  0.0, 0.0);
    GSL_SF_RESULT_SET(Fp, 0.0, 0.0);
    GSL_SF_RESULT_SET(G,  0.0, 0.0); /* FIXME: should be Inf */
    GSL_SF_RESULT_SET(Gp, 0.0, 0.0); /* FIXME: should be Inf */
    *exp_F = MpIeee( "0.0" );
    *exp_G = MpIeee( "0.0" );
    if(lam_F == 0.0){
      GSL_SF_RESULT_SET(Fp, C0.val, C0.err);
    }
    if(lam_G == 0.0) {
      GSL_SF_RESULT_SET(Gp, 1.0/C0.val, fabs(C0.err/C0.val)/fabs(C0.val));
    }
    GSL_ERROR ("domain error", GSL_EDOM);
    /* After all, since we are asking for G, this is a domain error... */
  }
  else if(x < 1.2 && 2.0*M_PI*eta < 0.9*(-GSL_LOG_DBL_MIN) && fabs(eta*x) < 10.0) {
    /* Reduce to a small lambda value and use the series
     * representations for F and G. We cannot allow eta to
     * be large and positive because the connection formula
     * for G_lam is badly behaved due to an underflow in sin(phi_lam) 
     * [see coulomb_FG_series() and coulomb_connection() above].
     * Note that large negative eta is ok however.
     */
    const MpIeee SMALL=  GSL_SQRT_DBL_EPSILON;
    const int N    = (int)MpIeee(lam_F + MpIeee("0.5")).toInt();
    const int span = GSL_MAX(k_lam_G, N);
    const MpIeee lam_min=  lam_F - N;    /* -1/2 <= lam_min < 1/2 */
    MpIeee F_lam_F;MpIeee  Fp_lam_F;
    MpIeee G_lam_G;MpIeee  Gp_lam_G;
    MpIeee F_lam_F_err;MpIeee  Fp_lam_F_err;
    MpIeee Fp_over_F_lam_F;
    MpIeee F_sign_lam_F;
    MpIeee F_lam_min_unnorm;MpIeee  Fp_lam_min_unnorm;
    MpIeee Fp_over_F_lam_min;
    gsl_sf_result F_lam_min;
    gsl_sf_result G_lam_min, Gp_lam_min;
    MpIeee F_scale;
    MpIeee Gerr_frac;
    MpIeee F_scale_frac_err;
    MpIeee F_unnorm_frac_err;

    /* Determine F'/F at lam_F. */
    int  CF1_count;
    int  stat_CF1=  coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count);

    int  stat_ser;
    int  stat_Fr;
    int  stat_Gr;

    /* Recurse down with unnormalized F,F' values. */
    F_lam_F  = SMALL;
    Fp_lam_F = Fp_over_F_lam_F * F_lam_F;
    if(span != 0) {
      stat_Fr = coulomb_F_recur(lam_min, span, eta, x,
                                F_lam_F, Fp_lam_F,
                                &F_lam_min_unnorm, &Fp_lam_min_unnorm
                                );
    }
    else {
      F_lam_min_unnorm  =  F_lam_F;
      Fp_lam_min_unnorm = Fp_lam_F;
      stat_Fr = GSL_SUCCESS;
    }

    /* Determine F and G at lam_min. */
    if(lam_min == -0.5) {
      stat_ser = coulomb_FGmhalf_series(eta, x, &F_lam_min, &G_lam_min);
    }
    else if(lam_min == 0.0) {
      stat_ser = coulomb_FG0_series(eta, x, &F_lam_min, &G_lam_min);
    }
    else if(lam_min == 0.5) {
      /* This cannot happen. */
      F->val  = F_lam_F;
      F->err  = 2.0 * GSL_DBL_EPSILON * fabs(F->val);
      Fp->val = Fp_lam_F;
      Fp->err = 2.0 * GSL_DBL_EPSILON * fabs(Fp->val);
      G->val  = G_lam_G;
      G->err  = 2.0 * GSL_DBL_EPSILON * fabs(G->val);
      Gp->val = Gp_lam_G;
      Gp->err = 2.0 * GSL_DBL_EPSILON * fabs(Gp->val);
      *exp_F = MpIeee( "0.0" );
      *exp_G = MpIeee( "0.0" );
      GSL_ERROR ("error", GSL_ESANITY);
    }
    else {
      stat_ser = coulomb_FG_series(lam_min, eta, x, &F_lam_min, &G_lam_min);
    }

    /* Determine remaining quantities. */
    Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm;
    Gp_lam_min.val  = Fp_over_F_lam_min*G_lam_min.val - 1.0/F_lam_min.val;
    Gp_lam_min.err  = fabs(Fp_over_F_lam_min)*G_lam_min.err;
    Gp_lam_min.err += fabs(1.0/F_lam_min.val) * fabs(F_lam_min.err/F_lam_min.val);
    F_scale     = F_lam_min.val / F_lam_min_unnorm;

    /* Apply scale to the original F,F' values. */
    F_scale_frac_err  = fabs(F_lam_min.err/F_lam_min.val);
    F_unnorm_frac_err = MpIeee( "2.0" )*GSL_DBL_EPSILON*(CF1_count+span+MpIeee( "1" ));
    F_lam_F     *= F_scale;
    F_lam_F_err  = fabs(F_lam_F) * (F_unnorm_frac_err + F_scale_frac_err);
    Fp_lam_F    *= F_scale;
    Fp_lam_F_err = fabs(Fp_lam_F) * (F_unnorm_frac_err + F_scale_frac_err);

    /* Recurse up to get the required G,G' values. */
    stat_Gr = coulomb_G_recur(lam_min, GSL_MAX(N-k_lam_G,0), eta, x,
                              G_lam_min.val, Gp_lam_min.val,
                              &G_lam_G, &Gp_lam_G
                              );

    F->val  = F_lam_F;
    F->err  = F_lam_F_err;
    F->err += 2.0 * GSL_DBL_EPSILON * fabs(F_lam_F);

    Fp->val  = Fp_lam_F;
    Fp->err  = Fp_lam_F_err;
    Fp->err += 2.0 * GSL_DBL_EPSILON * fabs(Fp_lam_F);

    Gerr_frac = fabs(G_lam_min.err/G_lam_min.val) + fabs(Gp_lam_min.err/Gp_lam_min.val);

    G->val  = G_lam_G;
    G->err  = Gerr_frac * fabs(G_lam_G);
    G->err += 2.0 * (CF1_count+1) * GSL_DBL_EPSILON * fabs(G->val);

    Gp->val  = Gp_lam_G;
    Gp->err  = Gerr_frac * fabs(Gp->val);
    Gp->err += 2.0 * (CF1_count+1) * GSL_DBL_EPSILON * fabs(Gp->val);

    *exp_F = MpIeee( "0.0" );
    *exp_G = MpIeee( "0.0" );

    return GSL_ERROR_SELECT_4(stat_ser, stat_CF1, stat_Fr, stat_Gr);
  }
  else if(x < 2.0*eta) {
    /* Use WKB approximation to obtain F and G at the two
     * lambda values, and use the Wronskian and the
     * continued fractions for F'/F to obtain F' and G'.
     */
    gsl_sf_result F_lam_F, G_lam_F;
    gsl_sf_result F_lam_G, G_lam_G;
    MpIeee exp_lam_F;MpIeee  exp_lam_G;
    int  stat_lam_F;
    int  stat_lam_G;
    int  stat_CF1_lam_F;
    int  stat_CF1_lam_G;
    int  CF1_count;
    MpIeee Fp_over_F_lam_F;
    MpIeee Fp_over_F_lam_G;
    MpIeee F_sign_lam_F;
    MpIeee F_sign_lam_G;

    stat_lam_F = coulomb_jwkb(lam_F, eta, x, &F_lam_F, &G_lam_F, &exp_lam_F);
    if(k_lam_G == 0) {
      stat_lam_G = stat_lam_F;
      F_lam_G = F_lam_F;
      G_lam_G = G_lam_F;
      exp_lam_G = exp_lam_F;
    }
    else {
      stat_lam_G = coulomb_jwkb(lam_G, eta, x, &F_lam_G, &G_lam_G, &exp_lam_G);
    }

    stat_CF1_lam_F = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count);
    if(k_lam_G == 0) {
      stat_CF1_lam_G  = stat_CF1_lam_F;
      F_sign_lam_G    = F_sign_lam_F;
      Fp_over_F_lam_G = Fp_over_F_lam_F;
    }
    else {
      stat_CF1_lam_G = coulomb_CF1(lam_G, eta, x, &F_sign_lam_G, &Fp_over_F_lam_G, &CF1_count);
    }

    F->val = F_lam_F.val;
    F->err = F_lam_F.err;

    G->val = G_lam_G.val;
    G->err = G_lam_G.err;

    Fp->val  = Fp_over_F_lam_F * F_lam_F.val;
    Fp->err  = fabs(Fp_over_F_lam_F) * F_lam_F.err;
    Fp->err += 2.0*GSL_DBL_EPSILON*fabs(Fp->val);

    Gp->val  = Fp_over_F_lam_G * G_lam_G.val - 1.0/F_lam_G.val;
    Gp->err  = fabs(Fp_over_F_lam_G) * G_lam_G.err;
    Gp->err += fabs(1.0/F_lam_G.val) * fabs(F_lam_G.err/F_lam_G.val);

    *exp_F = exp_lam_F;
    *exp_G = exp_lam_G;

    if(stat_lam_F == GSL_EOVRFLW || stat_lam_G == GSL_EOVRFLW) {
      GSL_ERROR ("overflow", GSL_EOVRFLW);
    }
    else {
      return GSL_ERROR_SELECT_2(stat_lam_F, stat_lam_G);
    }
  }
  else {
    /* x > 2 eta, so we know that we can find a lambda value such
     * that x is above the turning point. We do this, evaluate
     * using Steed's method at that oscillatory point, then
     * use recursion on F and G to obtain the required values.
     *
     * lam_0   = a value of lambda such that x is below the turning point
     * lam_min = minimum of lam_0 and the requested lam_G, since
     *           we must go at least as low as lam_G
     */
    const MpIeee SMALL=  GSL_SQRT_DBL_EPSILON;
    const MpIeee C=  sqrt(1.0 + 4.0*x*(x-2.0*eta));
    const int N = ceil(lam_F - C + 0.5);
    const MpIeee lam_0=  lam_F - GSL_MAX(N, 0);
    const MpIeee lam_min=  GSL_MIN(lam_0, lam_G);
    MpIeee F_lam_F;MpIeee  Fp_lam_F;
    MpIeee G_lam_G;MpIeee  Gp_lam_G;
    MpIeee F_lam_min_unnorm;MpIeee  Fp_lam_min_unnorm;
    MpIeee F_lam_min;MpIeee  Fp_lam_min;
    MpIeee G_lam_min;MpIeee  Gp_lam_min;
    MpIeee Fp_over_F_lam_F;
    MpIeee Fp_over_F_lam_min;
    MpIeee F_sign_lam_F;
    MpIeee P_lam_min;MpIeee  Q_lam_min;
    MpIeee alpha;
    MpIeee gamma;
    MpIeee F_scale;

    int  CF1_count;
    int  CF2_count;
    int  stat_CF1=  coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count);
    int  stat_CF2;
    int  stat_Fr;
    int  stat_Gr;

    int  F_recur_count;
    int  G_recur_count;

    MpIeee err_amplify;

    F_lam_F  = SMALL;
    Fp_lam_F = Fp_over_F_lam_F * F_lam_F;

    /* Backward recurrence to get F,Fp at lam_min */
    F_recur_count = GSL_MAX(k_lam_G, N);
    stat_Fr = coulomb_F_recur(lam_min, F_recur_count, eta, x,
                              F_lam_F, Fp_lam_F,
                              &F_lam_min_unnorm, &Fp_lam_min_unnorm
                              );
    Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm;

    /* Steed evaluation to complete evaluation of F,Fp,G,Gp at lam_min */
    stat_CF2 = coulomb_CF2(lam_min, eta, x, &P_lam_min, &Q_lam_min, &CF2_count);
    alpha = Fp_over_F_lam_min - P_lam_min;
    gamma = alpha/Q_lam_min;
    F_lam_min  = F_sign_lam_F / sqrt(alpha*alpha/Q_lam_min + Q_lam_min);
    Fp_lam_min = Fp_over_F_lam_min * F_lam_min;
    G_lam_min  = gamma * F_lam_min;
    Gp_lam_min = (P_lam_min * gamma - Q_lam_min) * F_lam_min;

    /* Apply scale to values of F,Fp at lam_F (the top). */
    F_scale = F_lam_min / F_lam_min_unnorm;    
    F_lam_F  *= F_scale;
    Fp_lam_F *= F_scale;

    /* Forward recurrence to get G,Gp at lam_G (the top). */
    G_recur_count = GSL_MAX(N-k_lam_G,0);
    stat_Gr = coulomb_G_recur(lam_min, G_recur_count, eta, x,
                              G_lam_min, Gp_lam_min,
                              &G_lam_G, &Gp_lam_G
                              );

    err_amplify = CF1_count + CF2_count + F_recur_count + G_recur_count + MpIeee( "1" );

    F->val  = F_lam_F;
    F->err  = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(F->val);

    Fp->val = Fp_lam_F;
    Fp->err = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(Fp->val);

    G->val  = G_lam_G;
    G->err  = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(G->val);

    Gp->val = Gp_lam_G;
    Gp->err = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(Gp->val);

    *exp_F = MpIeee( "0.0" );
    *exp_G = MpIeee( "0.0" );

    return GSL_ERROR_SELECT_4(stat_CF1, stat_CF2, stat_Fr, stat_Gr);
  }
}


int
 gsl_sf_coulomb_wave_F_array(MpIeee lam_min, int  kmax,
                                 MpIeee eta, MpIeee x, 
                                 MpIeee * fc_array,
                                 MpIeee * F_exp)
{
  if(x == MpIeee( "0.0" )) {
    int  k;
    *F_exp = MpIeee( "0.0" );
    for(k=0; k<=kmax; k++) {
      fc_array[k] = MpIeee( "0.0" );
    }
    if(lam_min == MpIeee( "0.0" )){
      gsl_sf_result f_0;
      CLeta(0.0, eta, &f_0);
      fc_array[0] = f_0.val;
    }
    return GSL_SUCCESS;
  }
  else {
    const MpIeee x_inv=  1.0/x;
    const MpIeee lam_max=  lam_min + kmax;
    gsl_sf_result F, Fp;
    gsl_sf_result G, Gp;
    MpIeee G_exp;

    int  stat_FG=  gsl_sf_coulomb_wave_FG_e(eta, x, lam_max, 0,
                                              &F, &Fp, &G, &Gp, F_exp, &G_exp);

    MpIeee fcl=  F.val;
    MpIeee fpl=  Fp.val;
    MpIeee lam=  lam_max;
    int  k;

    fc_array[kmax] = F.val;

    for(k=kmax-1; k>=0; k--) {
      MpIeee el=  eta/lam;
      MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
      MpIeee sl=  el  + lam*x_inv;
      MpIeee fc_lm1=  (fcl*sl + fpl)/rl;
      fc_array[k]   = fc_lm1;
      fpl           =  fc_lm1*sl - fcl*rl;
      fcl           =  fc_lm1;
      lam -= MpIeee( "1.0" );
    }

    return stat_FG;
  }
}


int
 gsl_sf_coulomb_wave_FG_array(MpIeee lam_min, int  kmax,
                                  MpIeee eta, MpIeee x,
                                  MpIeee * fc_array, MpIeee * gc_array,
                                  MpIeee * F_exp, MpIeee * G_exp)
{
  const MpIeee x_inv=  1.0/x;
  const MpIeee lam_max=  lam_min + kmax;
  gsl_sf_result F, Fp;
  gsl_sf_result G, Gp;

  int  stat_FG=  gsl_sf_coulomb_wave_FG_e(eta, x, lam_max, kmax,
                                            &F, &Fp, &G, &Gp, F_exp, G_exp);

  MpIeee fcl=  F.val;
  MpIeee fpl=  Fp.val;
  MpIeee lam=  lam_max;
  int  k;

  MpIeee gcl;MpIeee  gpl;

  fc_array[kmax] = F.val;

  for(k=kmax-1; k>=0; k--) {
    MpIeee el=  eta/lam;
    MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
    MpIeee sl=  el  + lam*x_inv;
    MpIeee fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fc_array[k] = fc_lm1;
    fpl         =  fc_lm1*sl - fcl*rl;
    fcl         =  fc_lm1;
    lam -= MpIeee( "1.0" );
  }

  gcl = G.val;
  gpl = Gp.val;
  lam = lam_min + MpIeee( "1.0" );

  gc_array[0] = G.val;

  for(k=1; k<=kmax; k++) {
    MpIeee el=  eta/lam;
    MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
    MpIeee sl=  el + lam*x_inv;
    MpIeee gcl1=  (sl*gcl - gpl)/rl;
    gc_array[k] = gcl1;
    gpl         = rl*gcl - sl*gcl1;
    gcl         = gcl1;
    lam += MpIeee( "1.0" );
  }

  return stat_FG;
}


int
 gsl_sf_coulomb_wave_FGp_array(MpIeee lam_min, int  kmax,
                                   MpIeee eta, MpIeee x,
                                   MpIeee * fc_array, MpIeee * fcp_array,
                                   MpIeee * gc_array, MpIeee * gcp_array,
                                   MpIeee * F_exp, MpIeee * G_exp)

{
  const MpIeee x_inv=  1.0/x;
  const MpIeee lam_max=  lam_min + kmax;
  gsl_sf_result F, Fp;
  gsl_sf_result G, Gp;

  int  stat_FG=  gsl_sf_coulomb_wave_FG_e(eta, x, lam_max, kmax,
                                            &F, &Fp, &G, &Gp, F_exp, G_exp);

  MpIeee fcl=  F.val;
  MpIeee fpl=  Fp.val;
  MpIeee lam=  lam_max;
  int  k;

  MpIeee gcl;MpIeee  gpl;

  fc_array[kmax]  = F.val;
  fcp_array[kmax] = Fp.val;

  for(k=kmax-1; k>=0; k--) {
    MpIeee el=  eta/lam;
    MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
    MpIeee sl=  el  + lam*x_inv;
    MpIeee fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fc_array[k]  = fc_lm1;
    fpl          = fc_lm1*sl - fcl*rl;
    fcp_array[k] = fpl;
    fcl          =  fc_lm1;
    lam -= MpIeee( "1.0" );
  }

  gcl = G.val;
  gpl = Gp.val;
  lam = lam_min + MpIeee( "1.0" );

  gc_array[0]  = G.val;
  gcp_array[0] = Gp.val;

  for(k=1; k<=kmax; k++) {
    MpIeee el=  eta/lam;
    MpIeee rl=  sqrt(MpIeee( "1.0" ) + el*el);
    MpIeee sl=  el + lam*x_inv;
    MpIeee gcl1=  (sl*gcl - gpl)/rl;
    gc_array[k]  = gcl1;
    gpl          = rl*gcl - sl*gcl1;
    gcp_array[k] = gpl;
    gcl          = gcl1;
    lam += MpIeee( "1.0" );
  }

  return stat_FG;
}


int
 gsl_sf_coulomb_wave_sphF_array(MpIeee lam_min, int  kmax,
                                    MpIeee eta, MpIeee x,
                                    MpIeee * fc_array,
                                    MpIeee * F_exp)
{
  int  k;

  if(x < MpIeee( "0.0" ) || lam_min < -MpIeee( "0.5" )) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(x < MpIeee( "10.0" )/GSL_DBL_MAX) {
    for(k=0; k<=kmax; k++) {
      fc_array[k] = MpIeee( "0.0" );
    }
    if(lam_min == MpIeee( "0.0" )) {
      fc_array[0] = sqrt(C0sq(eta));
    }
    *F_exp = MpIeee( "0.0" );
    if(x == MpIeee( "0.0" ))
      return GSL_SUCCESS;
    else
      GSL_ERROR ("underflow", GSL_EUNDRFLW);
  }
  else {
    int  k;
    int  stat_F=  gsl_sf_coulomb_wave_F_array(lam_min, kmax,
                                                  eta, x, 
                                                  fc_array,
                                                  F_exp);

    for(k=0; k<=kmax; k++) {
      fc_array[k] = fc_array[k] / x;
    }
    return stat_F;
  }
}


