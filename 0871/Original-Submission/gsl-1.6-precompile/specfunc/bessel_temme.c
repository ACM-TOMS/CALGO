#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_temme.c
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

/* Calculate series for Y_nu and K_nu for small x and nu.
 * This is applicable for x < 2 and |nu|<=1/2.
 * These functions assume x > 0.
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_mode.h>
#include "bessel_temme.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/* nu = (x+1)/4, -1<x<1, 1/(2nu)(1/Gamma[1-nu]-1/Gamma[1+nu]) */
static MpIeee g1_dat[14] =  {
  -MpIeee( "1.14516408366268311786898152867" ),
   MpIeee( "0.00636085311347084238122955495" ),
   MpIeee( "0.00186245193007206848934643657" ),
   MpIeee( "0.000152833085873453507081227824" ),
   MpIeee( "0.000017017464011802038795324732" ),
  -MpIeee( "6.4597502923347254354668326451e-07" ),
  -MpIeee( "5.1819848432519380894104312968e-08" ),
   MpIeee( "4.5189092894858183051123180797e-10" ),
   MpIeee( "3.2433227371020873043666259180e-11" ),
   MpIeee( "6.8309434024947522875432400828e-13" ),
   MpIeee( "2.8353502755172101513119628130e-14" ),
  -MpIeee( "7.9883905769323592875638087541e-16" ),
  -MpIeee( "3.3726677300771949833341213457e-17" ),
  -MpIeee( "3.6586334809210520744054437104e-20" )
};
static cheb_series g1_cs = {
  g1_dat,
  13,
  -1, 1,
  7
};

/* nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu]) */
static MpIeee g2_dat[15] =  
{
  MpIeee( "1.882645524949671835019616975350" ),
 -MpIeee( "0.077490658396167518329547945212" ),  
 -MpIeee( "0.018256714847324929419579340950" ),
  MpIeee( "0.0006338030209074895795923971731" ),
  MpIeee( "0.0000762290543508729021194461175" ),
 -MpIeee( "9.5501647561720443519853993526e-07" ),
 -MpIeee( "8.8927268107886351912431512955e-08" ),
 -MpIeee( "1.9521334772319613740511880132e-09" ),
 -MpIeee( "9.4003052735885162111769579771e-11" ),
  MpIeee( "4.6875133849532393179290879101e-12" ),
  MpIeee( "2.2658535746925759582447545145e-13" ),
 -MpIeee( "1.1725509698488015111878735251e-15" ),
 -MpIeee( "7.0441338200245222530843155877e-17" ),
 -MpIeee( "2.4377878310107693650659740228e-18" ),
 -MpIeee( "7.5225243218253901727164675011e-20" )
};
static cheb_series g2_cs = {
  g2_dat,
  14,
  -1, 1,
  8
};


static
int
 gsl_sf_temme_gamma(const MpIeee nu, MpIeee * g_1pnu, MpIeee * g_1mnu, MpIeee * g1, MpIeee * g2)
{
  const MpIeee anu=  fabs(nu);    /* functions are even */
  const MpIeee x=  4.0*anu - 1.0;
  gsl_sf_result r_g1;
  gsl_sf_result r_g2;
  cheb_eval_e(&g1_cs, x, &r_g1);
  cheb_eval_e(&g2_cs, x, &r_g2);
  *g1 = r_g1.val;
  *g2 = r_g2.val;
  *g_1mnu = MpIeee( "1.0" )/(r_g2.val + nu * r_g1.val);
  *g_1pnu = MpIeee( "1.0" )/(r_g2.val - nu * r_g1.val);
  return GSL_SUCCESS;
}


int
 gsl_sf_bessel_Y_temme(const MpIeee nu, const MpIeee x,
                      gsl_sf_result * Ynu,
                      gsl_sf_result * Ynup1)
{
  const int max_iter = 15000;
  
  const MpIeee half_x=  0.5 * x;
  const MpIeee ln_half_x=  log(half_x);
  const MpIeee half_x_nu=  exp(nu*ln_half_x);
  const MpIeee pi_nu=  M_PI * nu;
  const MpIeee alpha=  pi_nu / 2.0;
  const MpIeee sigma=  -nu * ln_half_x;
  const MpIeee sinrat=  (fabs(pi_nu) < GSL_DBL_EPSILON ?  MpIeee("1.0") : pi_nu/sin(pi_nu));
  const MpIeee sinhrat=  (fabs(sigma) < GSL_DBL_EPSILON ? MpIeee("1.0") : sinh(sigma)/sigma);
  const MpIeee sinhalf=  (fabs(alpha) < GSL_DBL_EPSILON ? MpIeee("1.0") : sin(alpha)/alpha);
  const MpIeee sin_sqr=  nu*M_PI*M_PI*0.5 * sinhalf*sinhalf;
  
  MpIeee sum0;MpIeee  sum1;
  MpIeee fk;MpIeee  pk;MpIeee  qk;MpIeee  hk;MpIeee  ck;
  int  k=  0;
  int  stat_iter;

  MpIeee g_1pnu;MpIeee  g_1mnu;MpIeee  g1;MpIeee  g2;
  int  stat_g=  gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);

  fk = MpIeee( "2.0" )/M_PI * sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = MpIeee( "1.0" )/M_PI /half_x_nu * g_1pnu;
  qk = MpIeee( "1.0" )/M_PI *half_x_nu * g_1mnu;
  hk = pk;
  ck = MpIeee( "1.0" );

  sum0 = fk + sin_sqr * qk;
  sum1 = pk;

  while(k < max_iter) {
    MpIeee del0;
    MpIeee del1;
    MpIeee gk;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= -half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    gk  = fk + sin_sqr * qk;
    hk  = -k*gk + pk; 
    del0 = ck * gk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if(fabs(del0) < 0.5*(1.0 + fabs(sum0))*GSL_DBL_EPSILON) break;
  }

  Ynu->val   = -sum0;
  Ynu->err   = (2.0 + 0.5*k) * GSL_DBL_EPSILON * fabs(Ynu->val);
  Ynup1->val = -sum1 * 2.0/x;
  Ynup1->err = (2.0 + 0.5*k) * GSL_DBL_EPSILON * fabs(Ynup1->val);

  stat_iter = ( k >= max_iter ? GSL_EMAXITER : GSL_SUCCESS );
  return GSL_ERROR_SELECT_2(stat_iter, stat_g);
}


int
 gsl_sf_bessel_K_scaled_temme(const MpIeee nu, const MpIeee x,
                             MpIeee * K_nu, MpIeee * K_nup1, MpIeee * Kp_nu)
{
  const int max_iter = 15000;

  const MpIeee half_x=  0.5 * x;
  const MpIeee ln_half_x=  log(half_x);
  const MpIeee half_x_nu=  exp(nu*ln_half_x);
  const MpIeee pi_nu=  M_PI * nu;
  const MpIeee sigma=  -nu * ln_half_x;
  const MpIeee sinrat=  (fabs(pi_nu) < GSL_DBL_EPSILON ? MpIeee("1.0") : pi_nu/sin(pi_nu));
  const MpIeee sinhrat=  (fabs(sigma) < GSL_DBL_EPSILON ? MpIeee("1.0") : sinh(sigma)/sigma);
  const MpIeee ex=  exp(x);

  MpIeee sum0;MpIeee  sum1;
  MpIeee fk;MpIeee  pk;MpIeee  qk;MpIeee  hk;MpIeee  ck;
  int  k=  0;
  int  stat_iter;

  MpIeee g_1pnu;MpIeee  g_1mnu;MpIeee  g1;MpIeee  g2;
  int  stat_g=  gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);

  fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = MpIeee( "0.5" )/half_x_nu * g_1pnu;
  qk = MpIeee( "0.5" )*half_x_nu * g_1mnu;
  hk = pk;
  ck = MpIeee( "1.0" );
  sum0 = fk;
  sum1 = hk;
  while(k < max_iter) {
    MpIeee del0;
    MpIeee del1;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    hk  = -k*fk + pk;
    del0 = ck * fk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if(fabs(del0) < 0.5*fabs(sum0)*GSL_DBL_EPSILON) break;
  }
  
  *K_nu   = sum0 * ex;
  *K_nup1 = sum1 * MpIeee( "2.0" )/x * ex;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;

  stat_iter = ( k == max_iter ? GSL_EMAXITER : GSL_SUCCESS );
  return GSL_ERROR_SELECT_2(stat_iter, stat_g);
}
