#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/legendre_poly.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
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
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_legendre.h>

#include "error.h"



/* Calculate P_m^m(x) from the analytic result:
 *   P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
 *            = 1 , m = 0
 */
static MpIeee legendre_Pmm(int  m, MpIeee x)
{
  if(m == 0)
  {
    return MpIeee( "1.0" );
  }
  else
  {
    MpIeee p_mm=  MpIeee( "1.0" );
    MpIeee root_factor=  sqrt(MpIeee( "1.0" )-x)*sqrt(MpIeee( "1.0" )+x);
    MpIeee fact_coeff=  MpIeee( "1.0" );
    int  i;
    for(i=1; i<=m; i++)
    {
      p_mm *= -fact_coeff * root_factor;
      fact_coeff += MpIeee( "2.0" );
    }
    return p_mm;
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_legendre_P1_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
}

int
 gsl_sf_legendre_P2_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = 0.5*(3.0*x*x - 1.0);
    result->err = GSL_DBL_EPSILON * (fabs(3.0*x*x) + 1.0);
    return GSL_SUCCESS;
  }
}

int
 gsl_sf_legendre_P3_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = 0.5*x*(5.0*x*x - 3.0);
    result->err = GSL_DBL_EPSILON * (fabs(result->val) + 0.5 * fabs(x) * (fabs(5.0*x*x) + 3.0));
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_legendre_Pl_e(const int l, const MpIeee x, gsl_sf_result * result)
{ 
  /* CHECK_POINTER(result) */

  if(l < 0 || x < -1.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(l == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 1) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 2) {
    result->val = 0.5 * (3.0*x*x - 1.0);
    result->err = GSL_DBL_EPSILON * (fabs(3.0*x*x) + 1.0);
    /*result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
      removed this old bogus estimate [GJ]
      */
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == -1.0) {
    result->val = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l < 100000) {
    /* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */

    MpIeee p_ellm2=  MpIeee( "1.0" );    /* P_0(x) */
    MpIeee p_ellm1=  x;      /* P_1(x) */
    MpIeee p_ell=  p_ellm1;
    int  ell;

    for(ell=2; ell <= l; ell++){
      p_ell = (x*(MpIeee( "2" )*ell-MpIeee( "1" ))*p_ellm1 - (ell-MpIeee( "1" ))*p_ellm2) / ell;
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;
    }

    result->val = p_ell;
    result->err = (0.5 * ell + 1.0) * GSL_DBL_EPSILON * fabs(p_ell);
    return GSL_SUCCESS;
  }
  else {
    /* Asymptotic expansion.
     * [Olver, p. 473]
     */
    MpIeee u=  l + MpIeee( "0.5" );
    MpIeee th=  acos(x);
    gsl_sf_result J0;
    gsl_sf_result Jm1;
    int  stat_J0=  gsl_sf_bessel_J0_e(u*th, &J0);
    int  stat_Jm1=  gsl_sf_bessel_Jn_e(-1, u*th, &Jm1);
    MpIeee pre;
    MpIeee B00;
    MpIeee c1;

    /* B00 = 1/8 (1 - th cot(th) / th^2
     * pre = sqrt(th/sin(th))
     */
    if(th < GSL_ROOT4_DBL_EPSILON) {
      B00 = (MpIeee( "1.0" ) + th*th/MpIeee( "15.0" ))/MpIeee( "24.0" );
      pre = MpIeee( "1.0" ) + th*th/MpIeee( "12.0" );
    }
    else {
      MpIeee sin_th=  sqrt(MpIeee( "1.0" ) - x*x);
      MpIeee cot_th=  x / sin_th;
      B00 = MpIeee( "1.0" )/MpIeee( "8.0" ) * (MpIeee( "1.0" ) - th * cot_th) / (th*th);
      pre = sqrt(th/sin_th);
    }

    c1 = th/u * B00;

    result->val  = pre * (J0.val + c1 * Jm1.val);
    result->err  = pre * (J0.err + fabs(c1) * Jm1.err);
    result->err += GSL_SQRT_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_J0, stat_Jm1);
  }
}


int
 gsl_sf_legendre_Pl_array(const int lmax, const MpIeee x, MpIeee * result_array)
{
  /* CHECK_POINTER(result_array) */

  if(lmax < 0 || x < -1.0 || x > 1.0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(lmax == 0) {
    result_array[0] = MpIeee( "1.0" );
    return GSL_SUCCESS;
  }
  else if(lmax == 1) {
    result_array[0] = MpIeee( "1.0" );
    result_array[1] = x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */

    MpIeee p_ellm2=  MpIeee( "1.0" );    /* P_0(x) */
    MpIeee p_ellm1=  x;    /* P_1(x) */
    MpIeee p_ell=  p_ellm1;
    int  ell;

    result_array[0] = MpIeee( "1.0" );
    result_array[1] = x;

    for(ell=2; ell <= lmax; ell++){
      p_ell = (x*(MpIeee( "2" )*ell-MpIeee( "1" ))*p_ellm1 - (ell-MpIeee( "1" ))*p_ellm2) / ell;
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;
      result_array[ell] = p_ell;
    }

    return GSL_SUCCESS;
  }
}


int
 gsl_sf_legendre_Pl_deriv_array(const int lmax, const MpIeee x, MpIeee * result_array, MpIeee * result_deriv_array)
{
  int  stat_array=  gsl_sf_legendre_Pl_array(lmax, x, result_array);

  if(lmax >= 0) result_deriv_array[0] = MpIeee( "0.0" );
  if(lmax >= 1) result_deriv_array[1] = MpIeee( "1.0" );

  if(stat_array == GSL_SUCCESS)
  {
    int  ell;

    if(fabs(x - 1.0)*(lmax+1.0)*(lmax+1.0) <  GSL_SQRT_DBL_EPSILON)
    {
      /* x is near 1 */
      for(ell = 2; ell <= lmax; ell++)
      {
        const MpIeee pre=  0.5 * ell * (ell+1.0);
        result_deriv_array[ell] = pre * (MpIeee( "1.0" ) - MpIeee( "0.25" ) * (MpIeee( "1.0" )-x) * (ell+MpIeee( "2.0" ))*(ell-MpIeee( "1.0" )));
      }
    }
    else if(fabs(x + 1.0)*(lmax+1.0)*(lmax+1.0) <  GSL_SQRT_DBL_EPSILON)
    {
      /* x is near -1 */
      for(ell = 2; ell <= lmax; ell++)
      {
        const MpIeee sgn=  ( GSL_IS_ODD(ell) ? 1.0 : -1.0 ); /* derivative is odd in x for even ell */
        const MpIeee pre=  sgn * 0.5 * ell * (ell+1.0);
        result_deriv_array[ell] = pre * (MpIeee( "1.0" ) - MpIeee( "0.25" ) * (MpIeee( "1.0" )+x) * (ell+MpIeee( "2.0" ))*(ell-MpIeee( "1.0" )));
      }
    }
    else
    {
      const MpIeee diff_a=  1.0 + x;
      const MpIeee diff_b=  1.0 - x;
      for(ell = 2; ell <= lmax; ell++)
      {
        result_deriv_array[ell] = - ell * (x * result_array[ell] - result_array[ell-1]) / (diff_a * diff_b);
      }
    }

    return GSL_SUCCESS;
  }
  else
  {
    return stat_array;
  }
}


int
 gsl_sf_legendre_Plm_e(const int l, const int m, const MpIeee x, gsl_sf_result * result)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  const MpIeee dif=  l-m;
  const MpIeee sum=  l+m;
  const MpIeee t_d=  ( dif == 0.0 ? 0.0 : 0.5 * dif * (log(dif)-1.0) );
  const MpIeee t_s=  ( dif == 0.0 ? 0.0 : 0.5 * sum * (log(sum)-1.0) );
  const MpIeee exp_check=  0.5 * log(2.0*l+1.0) + t_d - t_s;

  /* CHECK_POINTER(result) */

  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out. */
    OVERFLOW_ERROR(result);
  }
  else {
    /* Account for the error due to the
     * representation of 1-x.
     */
    const MpIeee err_amp=  1.0 / (GSL_DBL_EPSILON + fabs(1.0-fabs(x)));

    /* P_m^m(x) and P_{m+1}^m(x) */
    MpIeee p_mm=  legendre_Pmm(m, x);
    MpIeee p_mmp1=  x * (MpIeee( "2" )*m + MpIeee( "1" )) * p_mm;

    if(l == m){
      result->val = p_mm;
      result->err = err_amp * 2.0 * GSL_DBL_EPSILON * fabs(p_mm);
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      result->val = p_mmp1;
      result->err = err_amp * 2.0 * GSL_DBL_EPSILON * fabs(p_mmp1);
      return GSL_SUCCESS;
    }
    else{
      /* upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
       * start at P(m,m), P(m+1,m)
       */

      MpIeee p_ellm2=  p_mm;
      MpIeee p_ellm1=  p_mmp1;
      MpIeee p_ell=  MpIeee( "0.0" );
      int  ell;

      for(ell=m+2; ell <= l; ell++){
        p_ell = (x*(MpIeee( "2" )*ell-MpIeee( "1" ))*p_ellm1 - (ell+m-MpIeee( "1" ))*p_ellm2) / (ell-m);
        p_ellm2 = p_ellm1;
        p_ellm1 = p_ell;
      }

      result->val = p_ell;
      result->err = err_amp * (0.5*(l-m) + 1.0) * GSL_DBL_EPSILON * fabs(p_ell);

      return GSL_SUCCESS;
    }
  }
}


int
 gsl_sf_legendre_Plm_array(const int lmax, const int m, const MpIeee x, MpIeee * result_array)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  const MpIeee dif=  lmax-m;
  const MpIeee sum=  lmax+m;
  const MpIeee t_d=  ( dif == 0.0 ? 0.0 : 0.5 * dif * (log(dif)-1.0) );
  const MpIeee t_s=  ( dif == 0.0 ? 0.0 : 0.5 * sum * (log(sum)-1.0) );
  const MpIeee exp_check=  0.5 * log(2.0*lmax+1.0) + t_d - t_s;

  /* CHECK_POINTER(result_array) */

  if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(m > 0 && (x == 1.0 || x == -1.0)) {
    int  ell;
    for(ell=m; ell<=lmax; ell++) result_array[ell-m] = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out. */
    GSL_ERROR ("overflow", GSL_EOVRFLW);
  }
  else {
    MpIeee p_mm=  legendre_Pmm(m, x);
    MpIeee p_mmp1=  x * (MpIeee( "2.0" )*m + MpIeee( "1.0" )) * p_mm;

    if(lmax == m){
      result_array[0] = p_mm;
      return GSL_SUCCESS;
    }
    else if(lmax == m + 1) {
      result_array[0] = p_mm;
      result_array[1] = p_mmp1;
      return GSL_SUCCESS;
    }
    else {
      MpIeee p_ellm2=  p_mm;
      MpIeee p_ellm1=  p_mmp1;
      MpIeee p_ell=  MpIeee( "0.0" );
      int  ell;

      result_array[0] = p_mm;
      result_array[1] = p_mmp1;

      for(ell=m+2; ell <= lmax; ell++){
        p_ell = (x*(MpIeee( "2.0" )*ell-MpIeee( "1.0" ))*p_ellm1 - (ell+m-MpIeee( "1" ))*p_ellm2) / (ell-m);
        p_ellm2 = p_ellm1;
        p_ellm1 = p_ell;
        result_array[ell-m] = p_ell;
      }

      return GSL_SUCCESS;
    }
  }
}


int
 gsl_sf_legendre_Plm_deriv_array(
  const int lmax, const int m, const MpIeee x,
  MpIeee * result_array,
  MpIeee * result_deriv_array)
{
  if(m < 0 || m > lmax)
  {
    GSL_ERROR("m < 0 or m > lmax", GSL_EDOM);
  }
  else if(m == 0)
  {
    /* It is better to do m=0 this way, so we can more easily
     * trap the divergent case which can occur when m == 1.
     */
    return gsl_sf_legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array);
  }
  else
  {
    int  stat_array=  gsl_sf_legendre_Plm_array(lmax, m, x, result_array);

    if(stat_array == GSL_SUCCESS)
    {
      int  ell;

      if(m == 1 && (1.0 - fabs(x) < GSL_DBL_EPSILON))
      {
        /* This divergence is real and comes from the cusp-like
         * behaviour for m = 1. For example, P[1,1] = - Sqrt[1-x^2].
         */
        GSL_ERROR("divergence near |x| = 1.0 since m = 1", GSL_EOVRFLW);
      }
      else if(m == 2 && (1.0 - fabs(x) < GSL_DBL_EPSILON))
      {
        /* m = 2 gives a finite nonzero result for |x| near 1 */
        if(fabs(x - 1.0) < GSL_DBL_EPSILON)
        {
          for(ell = m; ell <= lmax; ell++) result_deriv_array[ell-m] = -MpIeee( "0.25" ) * x * (ell - MpIeee( "1.0" ))*ell*(ell+MpIeee( "1.0" ))*(ell+MpIeee( "2.0" ));
        }
        else if(fabs(x + 1.0) < GSL_DBL_EPSILON)
        {
          for(ell = m; ell <= lmax; ell++)
          {
            const MpIeee sgn=  ( GSL_IS_ODD(ell) ? 1.0 : -1.0 );
            result_deriv_array[ell-m] = -MpIeee( "0.25" ) * sgn * x * (ell - MpIeee( "1.0" ))*ell*(ell+MpIeee( "1.0" ))*(ell+MpIeee( "2.0" ));
          }
        }
        return GSL_SUCCESS;
      }
      else 
      {
        /* m > 2 is easier to deal with since the endpoints always vanish */
        if(1.0 - fabs(x) < GSL_DBL_EPSILON)
        {
          for(ell = m; ell <= lmax; ell++) result_deriv_array[ell-m] = MpIeee( "0.0" );
          return GSL_SUCCESS;
        }
        else
        {
          const MpIeee diff_a=  1.0 + x;
          const MpIeee diff_b=  1.0 - x;
          result_deriv_array[0] = - m * x / (diff_a * diff_b) * result_array[0];
          if(lmax-m >= 1) result_deriv_array[1] = (MpIeee( "2.0" ) * m + MpIeee( "1.0" )) * (x * result_deriv_array[0] + result_array[0]);
          for(ell = m+2; ell <= lmax; ell++)
          {
            result_deriv_array[ell-m] = - (ell * x * result_array[ell-m] - (ell+m) * result_array[ell-1-m]) / (diff_a * diff_b);
          }
          return GSL_SUCCESS;
        }
      }
    }
    else
    {
      return stat_array;
    }
  }
}


int
 gsl_sf_legendre_sphPlm_e(const int l, int  m, const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(m == 0) {
    gsl_sf_result P;
    int  stat_P=  gsl_sf_legendre_Pl_e(l, x, &P);
    MpIeee pre=  sqrt((MpIeee( "2.0" )*l + MpIeee( "1.0" ))/(MpIeee( "4.0" )*M_PI));
    result->val  = pre * P.val;
    result->err  = pre * P.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_P;
  }
  else if(x == 1.0 || x == -1.0) {
    /* m > 0 here */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* m > 0 and |x| < 1 here */

    /* Starting value for recursion.
     * Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
     */
    gsl_sf_result lncirc;
    gsl_sf_result lnpoch;
    MpIeee lnpre_val;
    MpIeee lnpre_err;
    gsl_sf_result ex_pre;
    MpIeee sr;
    const MpIeee sgn=  ( GSL_IS_ODD(m) ? -1.0 : 1.0);
    const MpIeee y_mmp1_factor=  x * sqrt(2.0*m + 3.0);
    MpIeee y_mm;MpIeee  y_mm_err;
    MpIeee y_mmp1;
    gsl_sf_log_1plusx_e(-x*x, &lncirc);
    gsl_sf_lnpoch_e(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
    lnpre_val = -MpIeee( "0.25" )*M_LNPI + MpIeee( "0.5" ) * (lnpoch.val + m*lncirc.val);
    lnpre_err = MpIeee( "0.25" )*M_LNPI*GSL_DBL_EPSILON + MpIeee( "0.5" ) * (lnpoch.err + fabs(m)*lncirc.err);
    gsl_sf_exp_err_e(lnpre_val, lnpre_err, &ex_pre);
    sr     = sqrt((MpIeee( "2.0" )+MpIeee( "1.0" )/m)/(MpIeee( "4.0" )*M_PI));
    y_mm   = sgn * sr * ex_pre.val;
    y_mmp1 = y_mmp1_factor * y_mm;
    y_mm_err  = MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(y_mm) + sr * ex_pre.err;
    y_mm_err *= MpIeee( "1.0" ) + MpIeee( "1.0" )/(GSL_DBL_EPSILON + fabs(MpIeee( "1.0" )-x));

    if(l == m){
      result->val  = y_mm;
      result->err  = y_mm_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(y_mm);
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      result->val  = y_mmp1;
      result->err  = fabs(y_mmp1_factor) * y_mm_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(y_mmp1);
      return GSL_SUCCESS;
    }
    else{
      MpIeee y_ell=  MpIeee( "0.0" );
      int  ell;

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      for(ell=m+2; ell <= l; ell++){
        const MpIeee rat1=  (const MpIeee)(ell-m)/(const MpIeee)(ell+m);
        const MpIeee rat2=  (ell-m-1.0)/(ell+m-1.0);
        const MpIeee factor1=  sqrt(rat1*(2*ell+1)*(2*ell-1));
        const MpIeee factor2=  sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*y_mmp1*factor1 - (ell+m-MpIeee( "1" ))*y_mm*factor2) / (ell-m);
        y_mm   = y_mmp1;
        y_mmp1 = y_ell;
      }

      result->val  = y_ell;
      result->err  = (0.5*(l-m) + 1.0) * GSL_DBL_EPSILON * fabs(y_ell);
      result->err += fabs(y_mm_err/y_mm) * fabs(y_ell);

      return GSL_SUCCESS;
    }
  }
}


int
 gsl_sf_legendre_sphPlm_array(const int lmax, int  m, const MpIeee x, MpIeee * result_array)
{
  /* CHECK_POINTER(result_array) */

  if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(m > 0 && (x == 1.0 || x == -1.0)) {
    int  ell;
    for(ell=m; ell<=lmax; ell++) result_array[ell-m] = MpIeee( "0.0" );
    return GSL_SUCCESS;
  }
  else {
    MpIeee y_mm;
    MpIeee y_mmp1;

    if(m == 0) {
      y_mm   = MpIeee( "0.5" )/M_SQRTPI;          /* Y00 = 1/sqrt(4pi) */
      y_mmp1 = x * M_SQRT3 * y_mm;
    }
    else {
      /* |x| < 1 here */

      gsl_sf_result lncirc;
      gsl_sf_result lnpoch;
      MpIeee lnpre;
      const MpIeee sgn=  ( GSL_IS_ODD(m) ? -1.0 : 1.0);
      gsl_sf_log_1plusx_e(-x*x, &lncirc);
      gsl_sf_lnpoch_e(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
      lnpre = -MpIeee( "0.25" )*M_LNPI + MpIeee( "0.5" ) * (lnpoch.val + m*lncirc.val);
      y_mm   = sqrt((MpIeee( "2.0" )+MpIeee( "1.0" )/m)/(MpIeee( "4.0" )*M_PI)) * sgn * exp(lnpre);
      y_mmp1 = x * sqrt(MpIeee( "2.0" )*m + MpIeee( "3.0" )) * y_mm;
    }

    if(lmax == m){
      result_array[0] = y_mm;
      return GSL_SUCCESS;
    }
    else if(lmax == m + 1) {
      result_array[0] = y_mm;
      result_array[1] = y_mmp1;
      return GSL_SUCCESS;
    }
    else{
      MpIeee y_ell;
      int  ell;

      result_array[0] = y_mm;
      result_array[1] = y_mmp1;

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      for(ell=m+2; ell <= lmax; ell++){
        const MpIeee rat1=  (const MpIeee)(ell-m)/(const MpIeee)(ell+m);
        const MpIeee rat2=  (ell-m-1.0)/(ell+m-1.0);
        const MpIeee factor1=  sqrt(rat1*(2*ell+1)*(2*ell-1));
        const MpIeee factor2=  sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*y_mmp1*factor1 - (ell+m-MpIeee( "1" ))*y_mm*factor2) / (ell-m);
        y_mm   = y_mmp1;
        y_mmp1 = y_ell;
        result_array[ell-m] = y_ell;
      }
    }

    return GSL_SUCCESS;
  }
}


int
 gsl_sf_legendre_sphPlm_deriv_array(
  const int lmax, const int m, const MpIeee x,
  MpIeee * result_array,
  MpIeee * result_deriv_array)
{
  if(m < 0 || lmax < m || x < -1.0 || x > 1.0)
  {
    GSL_ERROR ("domain", GSL_EDOM);
  }
  else if(m == 0)
  {
    /* m = 0 is easy to trap */
    const int stat_array = gsl_sf_legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array);
    int  ell;
    for(ell = 0; ell <= lmax; ell++)
    {
      const MpIeee prefactor=  sqrt((2.0 * ell + 1.0)/(4.0*M_PI));
      result_array[ell] *= prefactor;
      result_deriv_array[ell] *= prefactor;
    }
    return stat_array;
  }
  else if(m == 1)
  {
    /* Trapping m = 1 is necessary because of the possible divergence.
     * Recall that this divergence is handled properly in ..._Plm_deriv_array(),
     * and the scaling factor is not large for small m, so we just scale.
     */
    const int stat_array = gsl_sf_legendre_Plm_deriv_array(lmax, m, x, result_array, result_deriv_array);
    int  ell;
    for(ell = 1; ell <= lmax; ell++)
    {
      const MpIeee prefactor=  sqrt((2.0 * ell + 1.0)/(ell + 1.0) / (4.0*M_PI*ell));
      result_array[ell-1] *= prefactor;
      result_deriv_array[ell-1] *= prefactor;
    }
    return stat_array;
  }
  else
  {
    /* as for the derivative of P_lm, everything is regular for m >= 2 */

    int  stat_array=  gsl_sf_legendre_sphPlm_array(lmax, m, x, result_array);

    if(stat_array == GSL_SUCCESS)
    {
      int  ell;

      if(1.0 - fabs(x) < GSL_DBL_EPSILON)
      {
        for(ell = m; ell <= lmax; ell++) result_deriv_array[ell-m] = MpIeee( "0.0" );
        return GSL_SUCCESS;
      }
      else
      {
        const MpIeee diff_a=  1.0 + x;
        const MpIeee diff_b=  1.0 - x;
        result_deriv_array[0] = - m * x / (diff_a * diff_b) * result_array[0];
        if(lmax-m >= 1) result_deriv_array[1] = sqrt(MpIeee( "2.0" ) * m + MpIeee( "3.0" )) * (x * result_deriv_array[0] + result_array[0]);
        for(ell = m+2; ell <= lmax; ell++)
        {
          const MpIeee c1=  sqrt(((2.0*ell+1.0)/(2.0*ell-1.0)) * ((const MpIeee)(ell-m)/(const MpIeee)(ell+m)));
          result_deriv_array[ell-m] = - (ell * x * result_array[ell-m] - c1 * (ell+m) * result_array[ell-1-m]) / (diff_a * diff_b);
        }
        return GSL_SUCCESS;
      }
    }
    else
    {
      return stat_array;
    }
  }
}


#ifndef HIDE_INLINE_STATIC
int
 gsl_sf_legendre_array_size(const int lmax, const int m)
{
  return lmax-m+1;
}
#endif


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_legendre_P1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_P1_e(x, &result));
}

MpIeee gsl_sf_legendre_P2(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_P2_e(x, &result));
}

MpIeee gsl_sf_legendre_P3(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_P3_e(x, &result));
}

MpIeee gsl_sf_legendre_Pl(const int l, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_Pl_e(l, x, &result));
}

MpIeee gsl_sf_legendre_Plm(const int l, const int m, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_Plm_e(l, m, x, &result));
}

MpIeee gsl_sf_legendre_sphPlm(const int l, const int m, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_legendre_sphPlm_e(l, m, x, &result));
}

