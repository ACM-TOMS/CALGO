#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/legendre_con.c
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
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_legendre.h>

#include "error.h"
#include "legendre.h"

#define Root_2OverPi_  0.797884560802865355879892
#define locEPS         (1000.0*GSL_DBL_EPSILON)


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


#define RECURSE_LARGE  (1.0e-5*GSL_DBL_MAX)
#define RECURSE_SMALL  (1.0e+5*GSL_DBL_MIN)


/* Continued fraction for f_{ell+1}/f_ell
 * f_ell := P^{-mu-ell}_{-1/2 + I tau}(x),  x < 1.0
 *
 * Uses standard CF method from Temme's book.
 */
static
int
 conicalP_negmu_xlt1_CF1(const MpIeee mu, const int ell, const MpIeee tau,
                        const MpIeee x, gsl_sf_result * result)
{
  const MpIeee RECUR_BIG=  GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int  n=  1;
  MpIeee xi=  x/(sqrt(MpIeee( "1.0" )-x)*sqrt(MpIeee( "1.0" )+x));
  MpIeee Anm2=  MpIeee( "1.0" );
  MpIeee Bnm2=  MpIeee( "0.0" );
  MpIeee Anm1=  MpIeee( "0.0" );
  MpIeee Bnm1=  MpIeee( "1.0" );
  MpIeee a1=  MpIeee( "1.0" );
  MpIeee b1=  MpIeee( "2.0" )*(mu + ell + MpIeee( "1.0" )) * xi;
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
    an = tau*tau + (mu - MpIeee( "0.5" ) + ell + n)*(mu - MpIeee( "0.5" ) + ell + n);
    bn = MpIeee( "2.0" )*(ell + mu + n) * xi;
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
    
    if(fabs(del - 1.0) < 2.0*GSL_DBL_EPSILON) break;
  }

  result->val = fn;
  result->err = 4.0 * GSL_DBL_EPSILON * (sqrt(n) + 1.0) * fabs(fn);

  if(n >= maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Continued fraction for f_{ell+1}/f_ell
 * f_ell := P^{-mu-ell}_{-1/2 + I tau}(x),  x >= 1.0
 *
 * Uses Gautschi (Euler) equivalent series.
 */
static
int
 conicalP_negmu_xgt1_CF1(const MpIeee mu, const int ell, const MpIeee tau,
                        const MpIeee x, gsl_sf_result * result)
{ 
  const int maxk = 20000;
  const MpIeee gamma=  1.0-1.0/(x*x);
  const MpIeee pre=  sqrt(x-1.0)*sqrt(x+1.0) / (x*(2.0*(ell+mu+1.0)));
  MpIeee tk=  MpIeee( "1.0" );
  MpIeee sum=  MpIeee( "1.0" );
  MpIeee rhok=  MpIeee( "0.0" );
  int  k;
 
  for(k=1; k<maxk; k++) {
    MpIeee tlk=  MpIeee( "2.0" )*(ell + mu + k);
    MpIeee l1k=  (ell + mu - MpIeee( "0.5" ) + MpIeee( "1.0" ) + k);
    MpIeee ak=  -(tau*tau + l1k*l1k)/(tlk*(tlk+MpIeee( "2.0" ))) * gamma;
    rhok = -ak*(MpIeee( "1.0" ) + rhok)/(MpIeee( "1.0" ) + ak*(MpIeee( "1.0" ) + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < GSL_DBL_EPSILON) break;
  }

  result->val  = pre * sum;
  result->err  = fabs(pre * tk);
  result->err += 2.0 * GSL_DBL_EPSILON * (sqrt(k) + 1.0) * fabs(pre*sum);

  if(k >= maxk)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Implementation of large negative mu asymptotic
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */

inline
static MpIeee olver_U1(MpIeee beta2, MpIeee p)
{
  return (p-MpIeee( "1.0" ))/(MpIeee( "24.0" )*(MpIeee( "1.0" )+beta2)) * (MpIeee( "3.0" ) + beta2*(MpIeee( "2.0" ) + MpIeee( "5.0" )*p*(MpIeee( "1.0" )+p)));
}

inline
static MpIeee olver_U2(MpIeee beta2, MpIeee p)
{
  MpIeee beta4=  beta2*beta2;
  MpIeee p2=  p*p;
  MpIeee poly1=   MpIeee( "4.0" )*beta4 + MpIeee( "84.0" )*beta2 - MpIeee( "63.0" );
  MpIeee poly2=  MpIeee( "16.0" )*beta4 + MpIeee( "90.0" )*beta2 - MpIeee( "81.0" );
  MpIeee poly3=  beta2*p2*(MpIeee( "97.0" )*beta2 - MpIeee( "432.0" ) + MpIeee( "77.0" )*p*(beta2-MpIeee( "6.0" )) - MpIeee( "385.0" )*beta2*p2*(MpIeee( "1.0" ) + p));
  return (MpIeee( "1.0" )-p)/(MpIeee( "1152.0" )*(MpIeee( "1.0" )+beta2)) * (poly1 + poly2 + poly3);
}

static const MpIeee U3c1[] =  {   -1307.0,   -1647.0,    3375.0,    3675.0 };
static const MpIeee U3c2[] =  {   29366.0,   35835.0, -252360.0, -272630.0,
                                276810.0,  290499.0 };
static const MpIeee U3c3[] =  {  -29748.0,   -8840.0, 1725295.0, 1767025.0,
                              -7313470.0, -754778.0, 6309875.0, 6480045.0 };
static const MpIeee U3c4[] =  {    2696.0,    -16740.0,   -524250.0,  -183975.0,
                              14670540.0,  14172939.0, -48206730.0, -48461985.0,
                              36756720.0,  37182145.0 };
static const MpIeee U3c5[] =  {       9136.0,      22480.0,     12760.0,
                                  -252480.0,   -4662165.0,   -1705341.0,
                                 92370135.0,   86244015.0, -263678415.0,
                               -260275015.0, 185910725.0,  185910725.0 };

#if 0
static MpIeee olver_U3(MpIeee beta2, MpIeee p)
{
  MpIeee beta4=  beta2*beta2;
  MpIeee beta6=  beta4*beta2;
  MpIeee opb2s=  (MpIeee( "1.0" )+beta2)*(MpIeee( "1.0" )+beta2);
  MpIeee den=  MpIeee( "39813120.0" ) * opb2s*opb2s;
  MpIeee poly1=  gsl_poly_eval(U3c1, MpIeee( "4" ), p);
  MpIeee poly2=  gsl_poly_eval(U3c2, MpIeee( "6" ), p);
  MpIeee poly3=  gsl_poly_eval(U3c3, MpIeee( "8" ), p);
  MpIeee poly4=  gsl_poly_eval(U3c4, MpIeee( "10" ), p);
  MpIeee poly5=  gsl_poly_eval(U3c5, MpIeee( "12" ), p);
  
  return (p-MpIeee( "1.0" ))*(     MpIeee( "1215.0" )*poly1 + MpIeee( "324.0" )*beta2*poly2
                 + MpIeee( "54.0" )*beta4*poly3 +  MpIeee( "12.0" )*beta6*poly4
                 + beta4*beta4*poly5
                 ) / den;
}
#endif /* 0 */


/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu -> Inf
 * |x| < 1
 *
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */
int
 gsl_sf_conicalP_xlt1_large_neg_mu_e(MpIeee mu, MpIeee tau, MpIeee x,
                                       gsl_sf_result * result, MpIeee * ln_multiplier)
{
  MpIeee beta=  tau/mu;
  MpIeee beta2=  beta*beta;
  MpIeee S=  beta * acos((MpIeee( "1.0" )-beta2)/(MpIeee( "1.0" )+beta2));
  MpIeee p=  x/sqrt(beta2*(MpIeee( "1.0" )-x*x) + MpIeee( "1.0" ));
  gsl_sf_result lg_mup1;
  int  lg_stat=  gsl_sf_lngamma_e(mu+1.0, &lg_mup1);
  MpIeee ln_pre_1=   MpIeee( "0.5" )*mu*(S - log(MpIeee( "1.0" )+beta2) + log((MpIeee( "1.0" )-p)/(MpIeee( "1.0" )+p))) - lg_mup1.val;
  MpIeee ln_pre_2=  -MpIeee( "0.25" ) * log(MpIeee( "1.0" ) + beta2*(MpIeee( "1.0" )-x));
  MpIeee ln_pre_3=  -tau * atan(p*beta);
  MpIeee ln_pre=  ln_pre_1 + ln_pre_2 + ln_pre_3;
  MpIeee sum=  MpIeee( "1.0" ) - olver_U1(beta2, p)/mu + olver_U2(beta2, p)/(mu*mu);

  if(sum == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
    *ln_multiplier = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else {
    int  stat_e=  gsl_sf_exp_mult_e(ln_pre, sum, result);
    if(stat_e != GSL_SUCCESS) {
      result->val = sum;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(sum);
      *ln_multiplier = ln_pre;
    }
    else {
      *ln_multiplier = MpIeee( "0.0" );
    }
    return lg_stat;
  }
}


/* Implementation of large tau asymptotic
 *
 * A_n^{-mu}, B_n^{-mu}  [Olver, p.465, 469]
 */

inline
static MpIeee olver_B0_xi(MpIeee mu, MpIeee xi)
{
  return (MpIeee( "1.0" ) - MpIeee( "4.0" )*mu*mu)/(MpIeee( "8.0" )*xi) * (MpIeee( "1.0" )/tanh(xi) - MpIeee( "1.0" )/xi);
}

static MpIeee olver_A1_xi(MpIeee mu, MpIeee xi, MpIeee x)
{
  MpIeee B=  olver_B0_xi(mu, xi);
  MpIeee psi;
  if(fabs(x - 1.0) < GSL_ROOT4_DBL_EPSILON) {
    MpIeee y=  x - MpIeee( "1.0" );
    MpIeee s=  -MpIeee( "1.0" )/MpIeee( "3.0" ) + y*(MpIeee( "2.0" )/MpIeee( "15.0" ) - y *(MpIeee( "61.0" )/MpIeee( "945.0" ) - MpIeee( "452.0" )/MpIeee( "14175.0" )*y));
    psi = (MpIeee( "4.0" )*mu*mu - MpIeee( "1.0" ))/MpIeee( "16.0" ) * s;
  }
  else {
    psi = (MpIeee( "4.0" )*mu*mu - MpIeee( "1.0" ))/MpIeee( "16.0" ) * (MpIeee( "1.0" )/(x*x-MpIeee( "1.0" )) - MpIeee( "1.0" )/(xi*xi));
  }
  return MpIeee( "0.5" )*xi*xi*B*B + (mu+MpIeee( "0.5" ))*B - psi + mu/MpIeee( "6.0" )*(MpIeee( "0.25" ) - mu*mu);
}

inline
static MpIeee olver_B0_th(MpIeee mu, MpIeee theta)
{
  return -(MpIeee( "1.0" ) - MpIeee( "4.0" )*mu*mu)/(MpIeee( "8.0" )*theta) * (MpIeee( "1.0" )/tan(theta) - MpIeee( "1.0" )/theta);
}

static MpIeee olver_A1_th(MpIeee mu, MpIeee theta, MpIeee x)
{
  MpIeee B=  olver_B0_th(mu, theta);
  MpIeee psi;
  if(fabs(x - 1.0) < GSL_ROOT4_DBL_EPSILON) {
    MpIeee y=  MpIeee( "1.0" ) - x;
    MpIeee s=  -MpIeee( "1.0" )/MpIeee( "3.0" ) + y*(MpIeee( "2.0" )/MpIeee( "15.0" ) - y *(MpIeee( "61.0" )/MpIeee( "945.0" ) - MpIeee( "452.0" )/MpIeee( "14175.0" )*y));
    psi = (MpIeee( "4.0" )*mu*mu - MpIeee( "1.0" ))/MpIeee( "16.0" ) * s;
  }
  else {
    psi = (MpIeee( "4.0" )*mu*mu - MpIeee( "1.0" ))/MpIeee( "16.0" ) * (MpIeee( "1.0" )/(x*x-MpIeee( "1.0" )) + MpIeee( "1.0" )/(theta*theta));
  }
  return -MpIeee( "0.5" )*theta*theta*B*B + (mu+MpIeee( "0.5" ))*B - psi + mu/MpIeee( "6.0" )*(MpIeee( "0.25" ) - mu*mu);
}


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * 1 < x
 * tau -> Inf 
 * [Olver, p. 469]
 */
int
 gsl_sf_conicalP_xgt1_neg_mu_largetau_e(const MpIeee mu, const MpIeee tau,
                                          const MpIeee x, MpIeee acosh_x,
                                          gsl_sf_result * result, MpIeee * ln_multiplier)
{
  MpIeee xi=  acosh_x;
  MpIeee ln_xi_pre;
  MpIeee ln_pre;
  MpIeee sumA;MpIeee  sumB;MpIeee  sum;
  MpIeee arg;
  gsl_sf_result J_mup1;
  gsl_sf_result J_mu;
  MpIeee J_mum1;

  if(xi < GSL_ROOT4_DBL_EPSILON) {
    ln_xi_pre = -xi*xi/MpIeee( "6.0" );           /* log(1.0 - xi*xi/6.0) */
  }
  else {
    gsl_sf_result lnshxi;
    gsl_sf_lnsinh_e(xi, &lnshxi);
    ln_xi_pre = log(xi) - lnshxi.val;     /* log(xi/sinh(xi) */
  }

  ln_pre = MpIeee( "0.5" )*ln_xi_pre - mu*log(tau);

  arg = tau*xi;

  gsl_sf_bessel_Jnu_e(mu + 1.0,   arg, &J_mup1);
  gsl_sf_bessel_Jnu_e(mu,         arg, &J_mu);
  J_mum1 = -J_mup1.val + MpIeee( "2.0" )*mu/arg*J_mu.val;      /* careful of mu < 1 */

  sumA = MpIeee( "1.0" ) - olver_A1_xi(-mu, xi, x)/(tau*tau);
  sumB = olver_B0_xi(-mu, xi);
  sum  = J_mu.val * sumA - xi/tau * J_mum1 * sumB;

  if(sum == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
    *ln_multiplier = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else {
    int  stat_e=  gsl_sf_exp_mult_e(ln_pre, sum, result);
    if(stat_e != GSL_SUCCESS) {
      result->val = sum;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(sum);
      *ln_multiplier = ln_pre;
    }
    else {
      *ln_multiplier = MpIeee( "0.0" );
    }
    return GSL_SUCCESS;
  }
}


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * -1 < x < 1
 * tau -> Inf 
 * [Olver, p. 473]
 */
int
 gsl_sf_conicalP_xlt1_neg_mu_largetau_e(const MpIeee mu, const MpIeee tau,
                                          const MpIeee x, const MpIeee acos_x,
                                          gsl_sf_result * result, MpIeee * ln_multiplier)
{
  MpIeee theta=  acos_x;
  MpIeee ln_th_pre;
  MpIeee ln_pre;
  MpIeee sumA;MpIeee  sumB;MpIeee  sum;MpIeee  sumerr;
  MpIeee arg;
  gsl_sf_result I_mup1, I_mu;
  MpIeee I_mum1;

  if(theta < GSL_ROOT4_DBL_EPSILON) {
    ln_th_pre = theta*theta/MpIeee( "6.0" );   /* log(1.0 + theta*theta/6.0) */
  }
  else {
    ln_th_pre = log(theta/sin(theta));
  }

  ln_pre = MpIeee( "0.5" ) * ln_th_pre - mu * log(tau);

  arg = tau*theta;
  gsl_sf_bessel_Inu_e(mu + 1.0,   arg, &I_mup1);
  gsl_sf_bessel_Inu_e(mu,         arg, &I_mu);
  I_mum1 = I_mup1.val + MpIeee( "2.0" )*mu/arg * I_mu.val; /* careful of mu < 1 */

  sumA = MpIeee( "1.0" ) - olver_A1_th(-mu, theta, x)/(tau*tau);
  sumB = olver_B0_th(-mu, theta);
  sum  = I_mu.val * sumA - theta/tau * I_mum1 * sumB;
  sumerr  = fabs(I_mu.err * sumA);
  sumerr += fabs(I_mup1.err * theta/tau * sumB);
  sumerr += fabs(I_mu.err   * theta/tau * sumB * MpIeee( "2.0" ) * mu/arg);

  if(sum == MpIeee( "0.0" )) {
    result->val = 0.0;
    result->err = 0.0;
    *ln_multiplier = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else {
    int  stat_e=  gsl_sf_exp_mult_e(ln_pre, sum, result);
    if(stat_e != GSL_SUCCESS) {
      result->val  = sum;
      result->err  = sumerr;
      result->err += GSL_DBL_EPSILON * fabs(sum);
      *ln_multiplier = ln_pre;
    }
    else {
      *ln_multiplier = MpIeee( "0.0" );
    }
    return GSL_SUCCESS;
  }
}


/* Hypergeometric function which appears in the
 * large x expansion below:
 *
 *   2F1(1/4 - mu/2 - I tau/2, 3/4 - mu/2 - I tau/2, 1 - I tau, y)
 *
 * Note that for the usage below y = 1/x^2;
 */
static
int
 conicalP_hyperg_large_x(const MpIeee mu, const MpIeee tau, const MpIeee y,
                        MpIeee * reF, MpIeee * imF)
{
  const int kmax = 1000;
  const MpIeee re_a=  0.25 - 0.5*mu;
  const MpIeee re_b=  0.75 - 0.5*mu;
  const MpIeee re_c=  1.0;
  const MpIeee im_a=  -0.5*tau;
  const MpIeee im_b=  -0.5*tau;
  const MpIeee im_c=  -tau;

  MpIeee re_sum=  MpIeee( "1.0" );
  MpIeee im_sum=  MpIeee( "0.0" );
  MpIeee re_term=  MpIeee( "1.0" );
  MpIeee im_term=  MpIeee( "0.0" );
  int  k;

  for(k=1; k<=kmax; k++) {
    MpIeee re_ak=  re_a + k - MpIeee( "1.0" );
    MpIeee re_bk=  re_b + k - MpIeee( "1.0" );
    MpIeee re_ck=  re_c + k - MpIeee( "1.0" );
    MpIeee im_ak=  im_a;
    MpIeee im_bk=  im_b;
    MpIeee im_ck=  im_c;
    MpIeee den=  re_ck*re_ck + im_ck*im_ck;
    MpIeee re_multiplier=  ((re_ak*re_bk - im_ak*im_bk)*re_ck + im_ck*(im_ak*re_bk + re_ak*im_bk)) / den;
    MpIeee im_multiplier=  ((im_ak*re_bk + re_ak*im_bk)*re_ck - im_ck*(re_ak*re_bk - im_ak*im_bk)) / den;
    MpIeee re_tmp=  re_multiplier*re_term - im_multiplier*im_term;
    MpIeee im_tmp=  im_multiplier*re_term + re_multiplier*im_term;
    MpIeee asum=  fabs(re_sum) + fabs(im_sum);
    re_term = y/k * re_tmp;
    im_term = y/k * im_tmp;
    if(fabs(re_term/asum) < GSL_DBL_EPSILON && fabs(im_term/asum) < GSL_DBL_EPSILON) break;
    re_sum += re_term;
    im_sum += im_term;
  }

  *reF = re_sum;
  *imF = im_sum;

  if(k == kmax)
    GSL_ERROR ("error", GSL_EMAXITER);
  else  
    return GSL_SUCCESS;
}


/* P^{mu}_{-1/2 + I tau}
 * x->Inf
 */
int
 gsl_sf_conicalP_large_x_e(const MpIeee mu, const MpIeee tau, const MpIeee x,
                             gsl_sf_result * result, MpIeee * ln_multiplier)
{
  /* 2F1 term
   */
  MpIeee y=  ( x < MpIeee( "0.5" )*GSL_SQRT_DBL_MAX ? MpIeee( "1.0" )/(x*x) : MpIeee( "0.0" ) );
  MpIeee reF;MpIeee  imF;
  int  stat_F=  conicalP_hyperg_large_x(mu, tau, y, &reF, &imF);

  /* f = Gamma(+i tau)/Gamma(1/2 - mu + i tau)
   * FIXME: shift so it's better for tau-> 0
   */
  gsl_sf_result lgr_num, lgth_num;
  gsl_sf_result lgr_den, lgth_den;
  int  stat_gn=  gsl_sf_lngamma_complex_e(0.0,tau,&lgr_num,&lgth_num);
  int  stat_gd=  gsl_sf_lngamma_complex_e(0.5-mu,tau,&lgr_den,&lgth_den);

  MpIeee angle=  lgth_num.val - lgth_den.val + atan2(imF,reF);

  MpIeee lnx=  log(x);
  MpIeee lnxp1=  log(x+MpIeee( "1.0" ));
  MpIeee lnxm1=  log(x-MpIeee( "1.0" ));
  MpIeee lnpre_const=  MpIeee( "0.5" )*M_LN2 - MpIeee( "0.5" )*M_LNPI;
  MpIeee lnpre_comm=  (mu-MpIeee( "0.5" ))*lnx - MpIeee( "0.5" )*mu*(lnxp1 + lnxm1);
  MpIeee lnpre_err=    GSL_DBL_EPSILON * (MpIeee( "0.5" )*M_LN2 + MpIeee( "0.5" )*M_LNPI)
                      + GSL_DBL_EPSILON * fabs((mu-MpIeee( "0.5" ))*lnx)
                      + GSL_DBL_EPSILON * fabs(MpIeee( "0.5" )*mu)*(fabs(lnxp1)+fabs(lnxm1));

  /*  result = pre*|F|*|f| * cos(angle - tau * (log(x)+M_LN2))
   */
  gsl_sf_result cos_result;
  int  stat_cos=  gsl_sf_cos_e(angle + tau*(log(x) + M_LN2), &cos_result);
  int  status=  GSL_ERROR_SELECT_4(stat_cos, stat_gd, stat_gn, stat_F);
  if(cos_result.val == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return status;
  }
  else {
    MpIeee lnFf_val=  MpIeee( "0.5" )*log(reF*reF+imF*imF) + lgr_num.val - lgr_den.val;
    MpIeee lnFf_err=  lgr_num.err + lgr_den.err + GSL_DBL_EPSILON * fabs(lnFf_val);
    MpIeee lnnoc_val=  lnpre_const + lnpre_comm + lnFf_val;
    MpIeee lnnoc_err=  lnpre_err + lnFf_err + GSL_DBL_EPSILON * fabs(lnnoc_val);
    int  stat_e=  gsl_sf_exp_mult_err_e(lnnoc_val, lnnoc_err,
                                          cos_result.val, cos_result.err,
                                          result);
    if(stat_e == GSL_SUCCESS) {
      *ln_multiplier = MpIeee( "0.0" );
    }
    else {
      result->val  = cos_result.val;
      result->err  = cos_result.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      *ln_multiplier = lnnoc_val;
    }
    return status;
  }
}


/* P^{mu}_{-1/2 + I tau}  first hypergeometric representation
 * -1 < x < 1
 * This is more effective for |x| small, however it will work w/o
 * reservation for any x < 0 because everything is positive
 * definite in that case.
 *
 * [Kolbig,   (3)] (note typo in args of gamma functions)
 * [Bateman, (22)] (correct form)
 */
static
int
 conicalP_xlt1_hyperg_A(MpIeee mu, MpIeee tau, MpIeee x, gsl_sf_result * result)
{
  MpIeee x2=  x*x;
  MpIeee err_amp=  MpIeee( "1.0" ) + MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-fabs(x)));
  MpIeee pre_val=  M_SQRTPI / pow(MpIeee( "0.5" )*sqrt(MpIeee( "1" )-x2), mu);
  MpIeee pre_err=  err_amp * GSL_DBL_EPSILON * (fabs(mu)+MpIeee( "1.0" )) * fabs(pre_val) ;
  gsl_sf_result ln_g1, ln_g2, arg_g1, arg_g2;
  gsl_sf_result F1, F2;
  gsl_sf_result pre1, pre2;
  MpIeee t1_val;MpIeee  t1_err;
  MpIeee t2_val;MpIeee  t2_err;

  int  stat_F1=  gsl_sf_hyperg_2F1_conj_e(0.25 - 0.5*mu, 0.5*tau, 0.5, x2, &F1);
  int  stat_F2=  gsl_sf_hyperg_2F1_conj_e(0.75 - 0.5*mu, 0.5*tau, 1.5, x2, &F2);
  int  status=  GSL_ERROR_SELECT_2(stat_F1, stat_F2);

  gsl_sf_lngamma_complex_e(0.75 - 0.5*mu, -0.5*tau, &ln_g1, &arg_g1);
  gsl_sf_lngamma_complex_e(0.25 - 0.5*mu, -0.5*tau, &ln_g2, &arg_g2);

  gsl_sf_exp_err_e(-2.0*ln_g1.val, 2.0*ln_g1.err, &pre1);
  gsl_sf_exp_err_e(-2.0*ln_g2.val, 2.0*ln_g2.err, &pre2);
  pre2.val *= -2.0*x;
  pre2.err *=  2.0*fabs(x);
  pre2.err +=  GSL_DBL_EPSILON * fabs(pre2.val);

  t1_val = pre1.val * F1.val;
  t1_err = fabs(pre1.val) * F1.err + pre1.err * fabs(F1.val);
  t2_val = pre2.val * F2.val;
  t2_err = fabs(pre2.val) * F2.err + pre2.err * fabs(F2.val);

  result->val  = pre_val * (t1_val + t2_val);
  result->err  = pre_val * (t1_err + t2_err);
  result->err += pre_err * fabs(t1_val + t2_val);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return status;
}


/* P^{mu}_{-1/2 + I tau}
 * defining hypergeometric representation
 * [Abramowitz+Stegun, 8.1.2]
 * 1 < x < 3
 * effective for x near 1
 *
 */
#if 0
static
int
 conicalP_def_hyperg(MpIeee mu, MpIeee tau, MpIeee x, MpIeee * result)
{
  MpIeee F;
  int  stat_F=  gsl_sf_hyperg_2F1_conj_renorm_e(0.5, tau, 1.0-mu, 0.5*(1.0-x), &F);
  *result = pow((x+MpIeee( "1.0" ))/(x-MpIeee( "1.0" )), MpIeee( "0.5" )*mu) * F;
  return stat_F;
}
#endif /* 0 */


/* P^{mu}_{-1/2 + I tau}  second hypergeometric representation
 * [Zhurina+Karmazina, (3.1)] 
 * -1 < x < 3
 * effective for x near 1
 *
 */
#if 0
static
int
 conicalP_xnear1_hyperg_C(MpIeee mu, MpIeee tau, MpIeee x, MpIeee * result)
{
  MpIeee ln_pre;MpIeee  arg_pre;
  MpIeee ln_g1;MpIeee  arg_g1;
  MpIeee ln_g2;MpIeee  arg_g2;
  MpIeee F;

  int  stat_F=  gsl_sf_hyperg_2F1_conj_renorm_e(0.5+mu, tau, 1.0+mu, 0.5*(1.0-x), &F);

  gsl_sf_lngamma_complex_e(0.5+mu, tau, &ln_g1, &arg_g1);
  gsl_sf_lngamma_complex_e(0.5-mu, tau, &ln_g2, &arg_g2);

  ln_pre  = mu*M_LN2 - MpIeee( "0.5" )*mu*log(fabs(x*x-MpIeee( "1.0" ))) + ln_g1 - ln_g2;
  arg_pre = arg_g1 - arg_g2;

  *result = exp(ln_pre) * F;
  return stat_F;
}
#endif /* 0 */


/* V0, V1 from Kolbig, m = 0
 */
static
int
 conicalP_0_V(const MpIeee t, const MpIeee f, const MpIeee tau, const MpIeee sgn,
             MpIeee * V0, MpIeee * V1)
{
  MpIeee C[8];
  MpIeee T[8];
  MpIeee H[8];
  MpIeee V[12];
  int  i;
  T[0] = MpIeee( "1.0" );
  H[0] = MpIeee( "1.0" );
  V[0] = MpIeee( "1.0" );
  for(i=1; i<=7; i++) {
    T[i] = T[i-1] * t;
    H[i] = H[i-1] * (t*f);
  }
  for(i=1; i<=11; i++) {
    V[i] = V[i-1] * tau;
  }

  C[0] = MpIeee( "1.0" );
  C[1] = (H[1]-MpIeee( "1.0" ))/(MpIeee( "8.0" )*T[1]);
  C[2] = (MpIeee( "9.0" )*H[2] + MpIeee( "6.0" )*H[1] - MpIeee( "15.0" ) - sgn*MpIeee( "8.0" )*T[2])/(MpIeee( "128.0" )*T[2]);
  C[3] = MpIeee( "5.0" )*(MpIeee( "15.0" )*H[3] + MpIeee( "27.0" )*H[2] + MpIeee( "21.0" )*H[1] - MpIeee( "63.0" ) - sgn*T[2]*(MpIeee( "16.0" )*H[1]+MpIeee( "24.0" )))/(MpIeee( "1024.0" )*T[3]);
  C[4] = MpIeee( "7.0" )*(MpIeee( "525.0" )*H[4] + MpIeee( "1500.0" )*H[3] + MpIeee( "2430.0" )*H[2] + MpIeee( "1980.0" )*H[1] - MpIeee( "6435.0" )
              + MpIeee( "192.0" )*T[4] - sgn*T[2]*(MpIeee( "720.0" )*H[2]+MpIeee( "1600.0" )*H[1]+MpIeee( "2160.0" ))
              ) / (MpIeee( "32768.0" )*T[4]);
  C[5] = MpIeee( "21.0" )*(MpIeee( "2835.0" )*H[5] + MpIeee( "11025.0" )*H[4] + MpIeee( "24750.0" )*H[3] + MpIeee( "38610.0" )*H[2]
               + MpIeee( "32175.0" )*H[1] - MpIeee( "109395.0" ) + T[4]*(MpIeee( "1984.0" )*H[1]+MpIeee( "4032.0" ))
               - sgn*T[2]*(MpIeee( "4800.0" )*H[3]+MpIeee( "15120.0" )*H[2]+MpIeee( "26400.0" )*H[1]+MpIeee( "34320.0" ))
               ) / (MpIeee( "262144.0" )*T[5]);
  C[6] = MpIeee( "11.0" )*(MpIeee( "218295.0" )*H[6] + MpIeee( "1071630.0" )*H[5] + MpIeee( "3009825.0" )*H[4] + MpIeee( "6142500.0" )*H[3]
               + MpIeee( "9398025.0" )*H[2] + MpIeee( "7936110.0" )*H[1] - MpIeee( "27776385.0" )
               + T[4]*(MpIeee( "254016.0" )*H[2]+MpIeee( "749952.0" )*H[1]+MpIeee( "1100736.0" ))
               - sgn*T[2]*(MpIeee( "441000.0" )*H[4] + MpIeee( "1814400.0" )*H[3] + MpIeee( "4127760.0" )*H[2]
                         + MpIeee( "6552000.0" )*H[1] + MpIeee( "8353800.0" ) + MpIeee( "31232.0" )*T[4]
                         )
               ) / (MpIeee( "4194304.0" )*T[6]);

  *V0 = C[0] + (-MpIeee( "4.0" )*C[3]/T[1]+C[4])/V[4]
             + (-MpIeee( "192.0" )*C[5]/T[3]+MpIeee( "144.0" )*C[6]/T[2])/V[8]
             + sgn * (-C[2]/V[2]
                      + (-MpIeee( "24.0" )*C[4]/T[2]+MpIeee( "12.0" )*C[5]/T[1]-C[6])/V[6] 
                      + (-MpIeee( "1920.0" )*C[6]/T[4])/V[10]
                      );
  *V1 = C[1]/V[1] + (MpIeee( "8.0" )*(C[3]/T[2]-C[4]/T[1])+C[5])/V[5]
                  + (MpIeee( "384.0" )*C[5]/T[4] - MpIeee( "768.0" )*C[6]/T[3])/V[9]
                  + sgn * ((MpIeee( "2.0" )*C[2]/T[1]-C[3])/V[3]
                           + (MpIeee( "48.0" )*C[4]/T[3]-MpIeee( "72.0" )*C[5]/T[2] + MpIeee( "18.0" )*C[6]/T[1])/V[7]
                           + (MpIeee( "3840.0" )*C[6]/T[5])/V[11]
                           );

  return GSL_SUCCESS;
}


/* V0, V1 from Kolbig, m = 1
 */
static
int
 conicalP_1_V(const MpIeee t, const MpIeee f, const MpIeee tau, const MpIeee sgn,
             MpIeee * V0, MpIeee * V1)
{
  MpIeee Cm1;
  MpIeee C[8];
  MpIeee T[8];
  MpIeee H[8];
  MpIeee V[12];
  int  i;
  T[0] = MpIeee( "1.0" );
  H[0] = MpIeee( "1.0" );
  V[0] = MpIeee( "1.0" );
  for(i=1; i<=7; i++) {
    T[i] = T[i-1] * t;
    H[i] = H[i-1] * (t*f);
  }
  for(i=1; i<=11; i++) {
    V[i] = V[i-1] * tau;
  }

  Cm1  = -MpIeee( "1.0" );
  C[0] = MpIeee( "3.0" )*(MpIeee( "1.0" )-H[1])/(MpIeee( "8.0" )*T[1]);
  C[1] = (-MpIeee( "15.0" )*H[2]+MpIeee( "6.0" )*H[1]+MpIeee( "9.0" )+sgn*MpIeee( "8.0" )*T[2])/(MpIeee( "128.0" )*T[2]);
  C[2] = MpIeee( "3.0" )*(-MpIeee( "35.0" )*H[3] - MpIeee( "15.0" )*H[2] + MpIeee( "15.0" )*H[1] + MpIeee( "35.0" ) + sgn*T[2]*(MpIeee( "32.0" )*H[1]+MpIeee( "8.0" )))/(MpIeee( "1024.0" )*T[3]);
  C[3] = (-MpIeee( "4725.0" )*H[4] - MpIeee( "6300.0" )*H[3] - MpIeee( "3150.0" )*H[2] + MpIeee( "3780.0" )*H[1] + MpIeee( "10395.0" )
          -MpIeee( "1216.0" )*T[4] + sgn*T[2]*(MpIeee( "6000.0" )*H[2]+MpIeee( "5760.0" )*H[1]+MpIeee( "1680.0" ))) / (MpIeee( "32768.0" )*T[4]);
  C[4] = MpIeee( "7.0" )*(-MpIeee( "10395.0" )*H[5] - MpIeee( "23625.0" )*H[4] - MpIeee( "28350.0" )*H[3] - MpIeee( "14850.0" )*H[2]
              +MpIeee( "19305.0" )*H[1] + MpIeee( "57915.0" ) - T[4]*(MpIeee( "6336.0" )*H[1]+MpIeee( "6080.0" ))
              + sgn*T[2]*(MpIeee( "16800.0" )*H[3] + MpIeee( "30000.0" )*H[2] + MpIeee( "25920.0" )*H[1] + MpIeee( "7920.0" ))
              ) / (MpIeee( "262144.0" )*T[5]);
  C[5] = (-MpIeee( "2837835.0" )*H[6] - MpIeee( "9168390.0" )*H[5] - MpIeee( "16372125.0" )*H[4] - MpIeee( "18918900" )*H[3]
          -MpIeee( "10135125.0" )*H[2] + MpIeee( "13783770.0" )*H[1] + MpIeee( "43648605.0" )
          -T[4]*(MpIeee( "3044160.0" )*H[2] + MpIeee( "5588352.0" )*H[1] + MpIeee( "4213440.0" ))
          +sgn*T[2]*(MpIeee( "5556600.0" )*H[4] + MpIeee( "14817600.0" )*H[3] + MpIeee( "20790000.0" )*H[2]
                     + MpIeee( "17297280.0" )*H[1] + MpIeee( "5405400.0" ) + MpIeee( "323072.0" )*T[4]
                     )
          ) / (MpIeee( "4194304.0" )*T[6]);
  C[6] = MpIeee( "0.0" );

  *V0 = C[0] + (-MpIeee( "4.0" )*C[3]/T[1]+C[4])/V[4]
             + (-MpIeee( "192.0" )*C[5]/T[3]+MpIeee( "144.0" )*C[6]/T[2])/V[8]
             + sgn * (-C[2]/V[2]
                      + (-MpIeee( "24.0" )*C[4]/T[2]+MpIeee( "12.0" )*C[5]/T[1]-C[6])/V[6] 
                      );
  *V1 = C[1]/V[1] + (MpIeee( "8.0" )*(C[3]/T[2]-C[4]/T[1])+C[5])/V[5]
                  + (MpIeee( "384.0" )*C[5]/T[4] - MpIeee( "768.0" )*C[6]/T[3])/V[9]
                  + sgn * (Cm1*V[1] + (MpIeee( "2.0" )*C[2]/T[1]-C[3])/V[3]
                           + (MpIeee( "48.0" )*C[4]/T[3]-MpIeee( "72.0" )*C[5]/T[2] + MpIeee( "18.0" )*C[6]/T[1])/V[7]
                           );

  return GSL_SUCCESS;
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* P^0_{-1/2 + I lambda}
 */
int
 gsl_sf_conicalP_0_e(const MpIeee lambda, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 1.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(lambda == 0.0) {
    gsl_sf_result K;
    int  stat_K;
    if(x < 1.0) {
      const MpIeee th=  acos(x);
      const MpIeee s=  sin(0.5*th);
      stat_K = gsl_sf_ellint_Kcomp_e(s, GSL_MODE_DEFAULT, &K);
      result->val  = 2.0/M_PI * K.val;
      result->err  = 2.0/M_PI * K.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_K;
    }
    else {
      const MpIeee xi=  acosh(x);
      const MpIeee c=  cosh(0.5*xi);
      const MpIeee t=  tanh(0.5*xi);
      stat_K = gsl_sf_ellint_Kcomp_e(t, GSL_MODE_DEFAULT, &K);
      result->val  = 2.0/M_PI / c * K.val;
      result->err  = 2.0/M_PI / c * K.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_K;
    }
  }
  else if(   (x <= 0.0 && lambda < 1000.0)
          || (x <  0.1 && lambda < 17.0)
          || (x <  0.2 && lambda < 5.0 )
    ) {
    return conicalP_xlt1_hyperg_A(0.0, lambda, x, result);
  }
  else if(   (x <= 0.2 && lambda < 17.0)
          || (x <= 1.5 && lambda < 20.0)
    ) {
    return gsl_sf_hyperg_2F1_conj_e(0.5, lambda, 1.0, (1.0-x)/2, result);
  }
  else if(1.5 < x && lambda < GSL_MAX(x,20.0)) {
    gsl_sf_result P;
    MpIeee lm;
    int  stat_P=  gsl_sf_conicalP_large_x_e(0.0, lambda, x,
                                              &P, &lm
                                              );
    int  stat_e=  gsl_sf_exp_mult_err_e(lm, 2.0*GSL_DBL_EPSILON * fabs(lm),
                                          P.val, P.err,
                                          result);
    return GSL_ERROR_SELECT_2(stat_e, stat_P);
  }
  else {
    MpIeee V0;MpIeee  V1;
    if(x < 1.0) {
      MpIeee th=  acos(x);
      MpIeee sth=  sqrt(MpIeee( "1.0" )-x*x);  /* sin(th) */
      gsl_sf_result I0, I1;
      int  stat_I0=  gsl_sf_bessel_I0_scaled_e(th * lambda, &I0);
      int  stat_I1=  gsl_sf_bessel_I1_scaled_e(th * lambda, &I1);
      int  stat_I=  GSL_ERROR_SELECT_2(stat_I0, stat_I1);
      int  stat_V=  conicalP_0_V(th, x/sth, lambda, -1.0, &V0, &V1);
      MpIeee bessterm=  V0 * I0.val + V1 * I1.val;
      MpIeee besserr=  fabs(V0) * I0.err + fabs(V1) * I1.err;
      MpIeee arg1=  th*lambda;
      MpIeee sqts=  sqrt(th/sth);
      int  stat_e=  gsl_sf_exp_mult_err_e(arg1, 4.0 * GSL_DBL_EPSILON * fabs(arg1),
                                            sqts * bessterm, sqts * besserr,
                                            result);
      return GSL_ERROR_SELECT_3(stat_e, stat_V, stat_I);
    }
    else {
      MpIeee sh=  sqrt(x-MpIeee( "1.0" ))*sqrt(x+MpIeee( "1.0" ));  /* sinh(xi)      */
      MpIeee xi=  log(x + sh);              /* xi = acosh(x) */
      gsl_sf_result J0, J1;
      int  stat_J0=  gsl_sf_bessel_J0_e(xi * lambda, &J0);
      int  stat_J1=  gsl_sf_bessel_J1_e(xi * lambda, &J1);
      int  stat_J=  GSL_ERROR_SELECT_2(stat_J0, stat_J1);
      int  stat_V=  conicalP_0_V(xi, x/sh, lambda, 1.0, &V0, &V1);
      MpIeee bessterm=  V0 * J0.val + V1 * J1.val;
      MpIeee besserr=  fabs(V0) * J0.err + fabs(V1) * J1.err;
      MpIeee pre_val=  sqrt(xi/sh);
      MpIeee pre_err=  MpIeee( "2.0" ) * fabs(pre_val);
      result->val  = pre_val * bessterm;
      result->err  = pre_val * besserr;
      result->err += pre_err * fabs(bessterm);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_V, stat_J);
    }
  }
}


/* P^1_{-1/2 + I lambda}
 */
int
 gsl_sf_conicalP_1_e(const MpIeee lambda, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(lambda == 0.0) {
    gsl_sf_result K, E;
    int  stat_K;int   stat_E;
    if(x == 1.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(x < 1.0) {
      if(1.0-x < GSL_SQRT_DBL_EPSILON) {
        MpIeee err_amp=  GSL_MAX_DBL(MpIeee( "1.0" ), MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-x)));
        result->val = 0.25/M_SQRT2 * sqrt(1.0-x) * (1.0 + 5.0/16.0 * (1.0-x));
        result->err = err_amp * 3.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
      else {
        const MpIeee th=  acos(x);
        const MpIeee s=  sin(0.5*th);
        const MpIeee c2=  1.0 - s*s;
        const MpIeee sth=  sin(th);
        const MpIeee pre=  2.0/(M_PI*sth);
        stat_K = gsl_sf_ellint_Kcomp_e(s, GSL_MODE_DEFAULT, &K);
        stat_E = gsl_sf_ellint_Ecomp_e(s, GSL_MODE_DEFAULT, &E);
        result->val  = pre * (E.val - c2 * K.val);
        result->err  = pre * (E.err + fabs(c2) * K.err);
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return stat_K;
      }
    }
    else {
      if(x-1.0 < GSL_SQRT_DBL_EPSILON) {
        MpIeee err_amp=  GSL_MAX_DBL(MpIeee( "1.0" ), MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-x)));
        result->val = -0.25/M_SQRT2 * sqrt(x-1.0) * (1.0 - 5.0/16.0 * (x-1.0));
        result->err = err_amp * 3.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
      else {
        const MpIeee xi=  acosh(x);
        const MpIeee c=  cosh(0.5*xi);
        const MpIeee t=  tanh(0.5*xi);
        const MpIeee sxi=  sinh(xi);
        const MpIeee pre=  2.0/(M_PI*sxi) * c;
        stat_K = gsl_sf_ellint_Kcomp_e(t, GSL_MODE_DEFAULT, &K);
        stat_E = gsl_sf_ellint_Ecomp_e(t, GSL_MODE_DEFAULT, &E);
        result->val  = pre * (E.val - K.val);
        result->err  = pre * (E.err + K.err);
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return stat_K;
      }
    }
  }
  else if(   (x <= 0.0 && lambda < 1000.0)
          || (x <  0.1 && lambda < 17.0)
          || (x <  0.2 && lambda < 5.0 )
    ) {
    return conicalP_xlt1_hyperg_A(1.0, lambda, x, result);
  }
  else if(   (x <= 0.2 && lambda < 17.0)
          || (x <  1.5 && lambda < 20.0)
    ) {
    const MpIeee arg=  fabs(x*x - 1.0);
    const MpIeee sgn=  GSL_SIGN(1.0 - x);
    const MpIeee pre=  0.5*(lambda*lambda + 0.25) * sgn * sqrt(arg);
    gsl_sf_result F;
    int  stat_F=  gsl_sf_hyperg_2F1_conj_e(1.5, lambda, 2.0, (1.0-x)/2, &F);
    result->val  = pre * F.val;
    result->err  = fabs(pre) * F.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_F;
  }
  else if(1.5 <= x && lambda < GSL_MAX(x,20.0)) {
    gsl_sf_result P;
    MpIeee lm;
    int  stat_P=  gsl_sf_conicalP_large_x_e(1.0, lambda, x,
                                              &P, &lm
                                              );
    int  stat_e=  gsl_sf_exp_mult_err_e(lm, 2.0 * GSL_DBL_EPSILON * fabs(lm),
                                          P.val, P.err,
                                          result);
    return GSL_ERROR_SELECT_2(stat_e, stat_P);
  }
  else {
    MpIeee V0;MpIeee  V1;
    if(x < 1.0) {
      const MpIeee sqrt_1mx=  sqrt(1.0 - x);
      const MpIeee sqrt_1px=  sqrt(1.0 + x);
      const MpIeee th=  acos(x);
      const MpIeee sth=  sqrt_1mx * sqrt_1px;  /* sin(th) */
      gsl_sf_result I0, I1;
      int  stat_I0=  gsl_sf_bessel_I0_scaled_e(th * lambda, &I0);
      int  stat_I1=  gsl_sf_bessel_I1_scaled_e(th * lambda, &I1);
      int  stat_I=  GSL_ERROR_SELECT_2(stat_I0, stat_I1);
      int  stat_V=  conicalP_1_V(th, x/sth, lambda, -1.0, &V0, &V1);
      MpIeee bessterm=  V0 * I0.val + V1 * I1.val;
      MpIeee besserr=   fabs(V0) * I0.err + fabs(V1) * I1.err
                       + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(V0 * I0.val)
                       + MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(V1 * I1.val);
      MpIeee arg1=  th * lambda;
      MpIeee sqts=  sqrt(th/sth);
      int  stat_e=  gsl_sf_exp_mult_err_e(arg1, 2.0 * GSL_DBL_EPSILON * fabs(arg1),
                                            sqts * bessterm, sqts * besserr,
                                            result);
      result->err *= 1.0/sqrt_1mx;
      return GSL_ERROR_SELECT_3(stat_e, stat_V, stat_I);
    }
    else {
      const MpIeee sqrt_xm1=  sqrt(x - 1.0);
      const MpIeee sqrt_xp1=  sqrt(x + 1.0);
      const MpIeee sh=  sqrt_xm1 * sqrt_xp1;  /* sinh(xi)      */
      const MpIeee xi=  log(x + sh);          /* xi = acosh(x) */
      const MpIeee xi_lam=  xi * lambda;
      gsl_sf_result J0, J1;
      const int stat_J0 = gsl_sf_bessel_J0_e(xi_lam, &J0);
      const int stat_J1 = gsl_sf_bessel_J1_e(xi_lam, &J1);
      const int stat_J  = GSL_ERROR_SELECT_2(stat_J0, stat_J1);
      const int stat_V  = conicalP_1_V(xi, x/sh, lambda, 1.0, &V0, &V1);
      const MpIeee bessterm=  V0 * J0.val + V1 * J1.val;
      const MpIeee besserr=  fabs(V0) * J0.err + fabs(V1) * J1.err
                       + 512.0 * 2.0 * GSL_DBL_EPSILON * fabs(V0 * J0.val)
                       + 512.0 * 2.0 * GSL_DBL_EPSILON * fabs(V1 * J1.val)
                       + GSL_DBL_EPSILON * fabs(xi_lam * V0 * J1.val)
                       + GSL_DBL_EPSILON * fabs(xi_lam * V1 * J0.val);
      const MpIeee pre=  sqrt(xi/sh);
      result->val  = pre * bessterm;
      result->err  = pre * besserr * sqrt_xp1 / sqrt_xm1;
      result->err += 4.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_V, stat_J);
    }
  }
}


/* P^{1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.8, 8.6.12]
 * checked OK [GJ] Fri May  8 12:24:36 MDT 1998 
 */
int  gsl_sf_conicalP_half_e(const MpIeee lambda, const MpIeee x,
                              gsl_sf_result * result
                              )
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 1.0) {
    MpIeee err_amp=  MpIeee( "1.0" ) + MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-fabs(x)));
    MpIeee ac=  acos(x);
    MpIeee den=  sqrt(sqrt(MpIeee( "1.0" )-x)*sqrt(MpIeee( "1.0" )+x));
    result->val  = Root_2OverPi_ / den * cosh(ac * lambda);
    result->err  = err_amp * 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    result->err *= fabs(ac * lambda) + 1.0;
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1 */
    MpIeee err_amp=  MpIeee( "1.0" ) + MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-fabs(x)));
    MpIeee sq_term=  sqrt(x-MpIeee( "1.0" ))*sqrt(x+MpIeee( "1.0" ));
    MpIeee ln_term=  log(x + sq_term);
    MpIeee den=  sqrt(sq_term);
    MpIeee carg_val=  lambda * ln_term;
    MpIeee carg_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(carg_val);
    gsl_sf_result cos_result;
    int  stat_cos=  gsl_sf_cos_err_e(carg_val, carg_err, &cos_result);
    result->val  = Root_2OverPi_ / den * cos_result.val;
    result->err  = err_amp * Root_2OverPi_ / den * cos_result.err;
    result->err += 4.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_cos;
  }
}


/* P^{-1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.9, 8.6.14]
 * checked OK [GJ] Fri May  8 12:24:43 MDT 1998 
 */
int  gsl_sf_conicalP_mhalf_e(const MpIeee lambda, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 1.0) {
    MpIeee ac=  acos(x);
    MpIeee den=  sqrt(sqrt(MpIeee( "1.0" )-x)*sqrt(MpIeee( "1.0" )+x));
    MpIeee arg=  ac * lambda;
    MpIeee err_amp=  MpIeee( "1.0" ) + MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-fabs(x)));
    if(fabs(arg) < GSL_SQRT_DBL_EPSILON) {
      result->val  = Root_2OverPi_ / den * ac;
      result->err  = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      result->err *= err_amp;
    }
    else {
      result->val  = Root_2OverPi_ / (den*lambda) * sinh(arg);
      result->err  = GSL_DBL_EPSILON * (fabs(arg)+1.0) * fabs(result->val);
      result->err *= err_amp;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1 */
    MpIeee sq_term=  sqrt(x-MpIeee( "1.0" ))*sqrt(x+MpIeee( "1.0" ));
    MpIeee ln_term=  log(x + sq_term);
    MpIeee den=  sqrt(sq_term);
    MpIeee arg_val=  lambda * ln_term;
    MpIeee arg_err=  MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(arg_val);
    if(arg_val < GSL_SQRT_DBL_EPSILON) {
      result->val = Root_2OverPi_ / den * ln_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result sin_result;
      int  stat_sin=  gsl_sf_sin_err_e(arg_val, arg_err, &sin_result);
      result->val  = Root_2OverPi_ / (den*lambda) * sin_result.val;
      result->err  = Root_2OverPi_ / fabs(den*lambda) * sin_result.err;
      result->err += 3.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_sin;
    }
  }
}


int  gsl_sf_conicalP_sph_reg_e(const int l, const MpIeee lambda,
                                 const MpIeee x,
                                 gsl_sf_result * result
                                 )
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0 || l < -1) {
    DOMAIN_ERROR(result);
  }
  else if(l == -1) {
    return gsl_sf_conicalP_half_e(lambda, x, result);
  }
  else if(l == 0) {
    return gsl_sf_conicalP_mhalf_e(lambda, x, result);
  }
  else if(x == 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 0.0) {
    MpIeee c=  MpIeee( "1.0" )/sqrt(MpIeee( "1.0" )-x*x);
    gsl_sf_result r_Pellm1;
    gsl_sf_result r_Pell;
    int  stat_0=  gsl_sf_conicalP_half_e(lambda, x, &r_Pellm1);  /* P^( 1/2) */
    int  stat_1=  gsl_sf_conicalP_mhalf_e(lambda, x, &r_Pell);   /* P^(-1/2) */
    int  stat_P=  GSL_ERROR_SELECT_2(stat_0, stat_1);
    MpIeee Pellm1=  r_Pellm1.val;
    MpIeee Pell=  r_Pell.val;
    MpIeee Pellp1;
    int  ell;

    for(ell=0; ell<l; ell++) {
      MpIeee d=  (ell+MpIeee( "1.0" ))*(ell+MpIeee( "1.0" )) + lambda*lambda;
      Pellp1 = (Pellm1 - (MpIeee( "2.0" )*ell+MpIeee( "1.0" ))*c*x * Pell) / d;
      Pellm1 = Pell;
      Pell   = Pellp1;
    }

    result->val  = Pell;
    result->err  = (0.5*l + 1.0) * GSL_DBL_EPSILON * fabs(Pell);
    result->err += GSL_DBL_EPSILON * l * fabs(result->val);
    return stat_P;
  }
  else if(x < 1.0) {
    const MpIeee xi=  x/(sqrt(1.0-x)*sqrt(1.0+x));
    gsl_sf_result rat;
    gsl_sf_result Phf;
    int  stat_CF1=  conicalP_negmu_xlt1_CF1(0.5, l, lambda, x, &rat);
    int  stat_Phf=  gsl_sf_conicalP_half_e(lambda, x, &Phf);
    MpIeee Pellp1=  rat.val * GSL_SQRT_DBL_MIN;
    MpIeee Pell=  GSL_SQRT_DBL_MIN;
    MpIeee Pellm1;
    int  ell;

    for(ell=l; ell>=0; ell--) {
      MpIeee d=  (ell+MpIeee( "1.0" ))*(ell+MpIeee( "1.0" )) + lambda*lambda;
      Pellm1 = (MpIeee( "2.0" )*ell+MpIeee( "1.0" ))*xi * Pell + d * Pellp1;
      Pellp1 = Pell;
      Pell   = Pellm1;
    }

    result->val  = GSL_SQRT_DBL_MIN * Phf.val / Pell;
    result->err  = GSL_SQRT_DBL_MIN * Phf.err / fabs(Pell);
    result->err += fabs(rat.err/rat.val) * (l + 1.0) * fabs(result->val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_Phf, stat_CF1);
  }
  else if(x == 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1.0 */

    const MpIeee xi=  x/sqrt((x-1.0)*(x+1.0));
    gsl_sf_result rat;
    int  stat_CF1=  conicalP_negmu_xgt1_CF1(0.5, l, lambda, x, &rat);
    int  stat_P;
    MpIeee Pellp1=  rat.val * GSL_SQRT_DBL_MIN;
    MpIeee Pell=  GSL_SQRT_DBL_MIN;
    MpIeee Pellm1;
    int  ell;

    for(ell=l; ell>=0; ell--) {
      MpIeee d=  (ell+MpIeee( "1.0" ))*(ell+MpIeee( "1.0" )) + lambda*lambda;
      Pellm1 = (MpIeee( "2.0" )*ell+MpIeee( "1.0" ))*xi * Pell - d * Pellp1;
      Pellp1 = Pell;
      Pell   = Pellm1;
    }

    if(fabs(Pell) > fabs(Pellp1)){
      gsl_sf_result Phf;
      stat_P = gsl_sf_conicalP_half_e(lambda, x, &Phf);
      result->val  =       GSL_SQRT_DBL_MIN * Phf.val / Pell;
      result->err  = 2.0 * GSL_SQRT_DBL_MIN * Phf.err / fabs(Pell);
      result->err += 2.0 * fabs(rat.err/rat.val) * (l + 1.0) * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
      gsl_sf_result Pmhf;
      stat_P = gsl_sf_conicalP_mhalf_e(lambda, x, &Pmhf);
      result->val  =       GSL_SQRT_DBL_MIN * Pmhf.val / Pellp1;
      result->err  = 2.0 * GSL_SQRT_DBL_MIN * Pmhf.err / fabs(Pellp1);
      result->err += 2.0 * fabs(rat.err/rat.val) * (l + 1.0) * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_ERROR_SELECT_2(stat_P, stat_CF1);
  }
}


int  gsl_sf_conicalP_cyl_reg_e(const int m, const MpIeee lambda,
                                 const MpIeee x,
                                 gsl_sf_result * result
                                 )
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0 || m < -1) {
    DOMAIN_ERROR(result);
  }
  else if(m == -1) {
    return gsl_sf_conicalP_1_e(lambda, x, result);
  }
  else if(m == 0) {
    return gsl_sf_conicalP_0_e(lambda, x, result);
  }
  else if(x == 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 0.0) {
    MpIeee c=  MpIeee( "1.0" )/sqrt(MpIeee( "1.0" )-x*x);
    gsl_sf_result r_Pkm1;
    gsl_sf_result r_Pk;
    int  stat_0=  gsl_sf_conicalP_1_e(lambda, x, &r_Pkm1);  /* P^1 */
    int  stat_1=  gsl_sf_conicalP_0_e(lambda, x, &r_Pk);    /* P^0 */
    int  stat_P=  GSL_ERROR_SELECT_2(stat_0, stat_1);
    MpIeee Pkm1=  r_Pkm1.val;
    MpIeee Pk=  r_Pk.val;
    MpIeee Pkp1;
    int  k;

    for(k=0; k<m; k++) {
      MpIeee d=  (k+MpIeee( "0.5" ))*(k+MpIeee( "0.5" )) + lambda*lambda;
      Pkp1 = (Pkm1 - MpIeee( "2.0" )*k*c*x * Pk) / d;
      Pkm1 = Pk;
      Pk   = Pkp1;
    }

    result->val  = Pk;
    result->err  = (m + 2.0) * GSL_DBL_EPSILON * fabs(Pk);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return stat_P;
  }
  else if(x < 1.0) {
    const MpIeee xi=  x/(sqrt(1.0-x)*sqrt(1.0+x));
    gsl_sf_result rat;
    gsl_sf_result P0;
    int  stat_CF1=  conicalP_negmu_xlt1_CF1(0.0, m, lambda, x, &rat);
    int  stat_P0=  gsl_sf_conicalP_0_e(lambda, x, &P0);
    MpIeee Pkp1=  rat.val * GSL_SQRT_DBL_MIN;
    MpIeee Pk=  GSL_SQRT_DBL_MIN;
    MpIeee Pkm1;
    int  k;

    for(k=m; k>0; k--) {
      MpIeee d=  (k+MpIeee( "0.5" ))*(k+MpIeee( "0.5" )) + lambda*lambda;
      Pkm1 = MpIeee( "2.0" )*k*xi * Pk + d * Pkp1;
      Pkp1 = Pk;
      Pk   = Pkm1;
    }

    result->val  = GSL_SQRT_DBL_MIN * P0.val / Pk;
    result->err  = 2.0 * GSL_SQRT_DBL_MIN * P0.err / fabs(Pk);
    result->err += 2.0 * fabs(rat.err/rat.val) * (m + 1.0) * fabs(result->val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_P0, stat_CF1);
  }
  else if(x == 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1.0 */

    const MpIeee xi=  x/sqrt((x-1.0)*(x+1.0));
    gsl_sf_result rat;
    int  stat_CF1=  conicalP_negmu_xgt1_CF1(0.0, m, lambda, x, &rat);
    int  stat_P;
    MpIeee Pkp1=  rat.val * GSL_SQRT_DBL_MIN;
    MpIeee Pk=  GSL_SQRT_DBL_MIN;
    MpIeee Pkm1;
    int  k;

    for(k=m; k>-1; k--) {
      MpIeee d=  (k+MpIeee( "0.5" ))*(k+MpIeee( "0.5" )) + lambda*lambda;
      Pkm1 = MpIeee( "2.0" )*k*xi * Pk - d * Pkp1;
      Pkp1 = Pk;
      Pk   = Pkm1;
    }

    if(fabs(Pk) > fabs(Pkp1)){
      gsl_sf_result P1;
      stat_P = gsl_sf_conicalP_1_e(lambda, x, &P1);
      result->val  = GSL_SQRT_DBL_MIN * P1.val / Pk;
      result->err  = 2.0 * GSL_SQRT_DBL_MIN * P1.err / fabs(Pk);
      result->err += 2.0 * fabs(rat.err/rat.val) * (m+2.0) * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
      gsl_sf_result P0;
      stat_P = gsl_sf_conicalP_0_e(lambda, x, &P0);
      result->val  = GSL_SQRT_DBL_MIN * P0.val / Pkp1;
      result->err  = 2.0 * GSL_SQRT_DBL_MIN * P0.err / fabs(Pkp1);
      result->err += 2.0 * fabs(rat.err/rat.val) * (m+2.0) * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_ERROR_SELECT_2(stat_P, stat_CF1);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_conicalP_0(const MpIeee lambda, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_conicalP_0_e(lambda, x, &result));
}

MpIeee gsl_sf_conicalP_1(const MpIeee lambda, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_conicalP_1_e(lambda, x, &result));
}

MpIeee gsl_sf_conicalP_half(const MpIeee lambda, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_conicalP_half_e(lambda, x, &result));
}

MpIeee gsl_sf_conicalP_mhalf(const MpIeee lambda, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_conicalP_mhalf_e(lambda, x, &result));
}

MpIeee gsl_sf_conicalP_sph_reg(const int l, const MpIeee lambda, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_conicalP_sph_reg_e(l, lambda, x, &result));
}

MpIeee gsl_sf_conicalP_cyl_reg(const int m, const MpIeee lambda, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_conicalP_cyl_reg_e(m, lambda, x, &result));
}
