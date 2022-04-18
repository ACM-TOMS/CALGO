#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_j.c
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "bessel.h"
#include "bessel_olver.h"

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_bessel_j0_e(const MpIeee x, gsl_sf_result * result)
{
  MpIeee ax=  fabs(x);

  /* CHECK_POINTER(result) */

  if(ax < MpIeee( "0.5" )) {
    const MpIeee y=  x*x;
    const MpIeee c1=  -1.0/6.0;
    const MpIeee c2=   1.0/120.0;
    const MpIeee c3=  -1.0/5040.0;
    const MpIeee c4=   1.0/362880.0;
    const MpIeee c5=  -1.0/39916800.0;
    const MpIeee c6=   1.0/6227020800.0;
    result->val = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*c6)))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result sin_result;
    const int stat = gsl_sf_sin_e(x, &sin_result);
    result->val  = sin_result.val/x;
    result->err  = fabs(sin_result.err/x);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
}


int  gsl_sf_bessel_j1_e(const MpIeee x, gsl_sf_result * result)
{
  MpIeee ax=  fabs(x);

  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax < MpIeee( "3.1" )*GSL_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(ax < MpIeee( "0.25" )) {
    const MpIeee y=  x*x;
    const MpIeee c1=  -1.0/10.0;
    const MpIeee c2=   1.0/280.0;
    const MpIeee c3=  -1.0/15120.0;
    const MpIeee c4=   1.0/1330560.0;
    const MpIeee c5=  -1.0/172972800.0;
    const MpIeee sum=  1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*c5))));
    result->val = x/3.0 * sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result cos_result;
    gsl_sf_result sin_result;
    const int stat_cos = gsl_sf_cos_e(x, &cos_result);
    const int stat_sin = gsl_sf_sin_e(x, &sin_result);
    const MpIeee cos_x=  cos_result.val;
    const MpIeee sin_x=  sin_result.val;
    result->val  = (sin_x/x - cos_x)/x;
    result->err  = (fabs(sin_result.err/x) + fabs(cos_result.err))/fabs(x);
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(sin_x/(x*x)) + fabs(cos_x/x));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_cos, stat_sin);
  }
}


int  gsl_sf_bessel_j2_e(const MpIeee x, gsl_sf_result * result)
{
  MpIeee ax=  fabs(x);

  /* CHECK_POINTER(result) */
  
  if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax < MpIeee( "4.0" )*GSL_SQRT_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(ax < MpIeee( "1.3" )) {
    const MpIeee y=  x*x;
    const MpIeee c1=  -1.0/14.0;
    const MpIeee c2=   1.0/504.0;
    const MpIeee c3=  -1.0/33264.0;
    const MpIeee c4=   1.0/3459456.0;
    const MpIeee c5=  -1.0/518918400;
    const MpIeee c6=   1.0/105859353600.0;
    const MpIeee c7=  -1.0/28158588057600.0;
    const MpIeee c8=   1.0/9461285587353600.0;
    const MpIeee c9=  -1.0/3916972233164390400.0;
    const MpIeee sum=  1.0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*c9))))))));
    result->val = y/15.0 * sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result cos_result;
    gsl_sf_result sin_result;
    const int stat_cos = gsl_sf_cos_e(x, &cos_result);
    const int stat_sin = gsl_sf_sin_e(x, &sin_result);
    const MpIeee cos_x=  cos_result.val;
    const MpIeee sin_x=  sin_result.val;
    const MpIeee f=  (3.0/(x*x) - 1.0);
    result->val  = (f * sin_x - 3.0*cos_x/x)/x;
    result->err  = fabs(f * sin_result.err/x) + fabs((3.0*cos_result.err/x)/x);
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(f*sin_x/x) + 3.0*fabs(cos_x/(x*x)));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_cos, stat_sin);
  }
}


int
 gsl_sf_bessel_jl_e(const int l, const MpIeee x, gsl_sf_result * result)
{
  if(l < 0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = ( l > 0 ? 0.0 : 1.0 );
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 0) {
    return gsl_sf_bessel_j0_e(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_j1_e(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_j2_e(x, result);
  }
  else if(x*x < 10.0*(l+0.5)/M_E) {
    gsl_sf_result b;
    int  status=  gsl_sf_bessel_IJ_taylor_e(l+0.5, x, -1, 50, GSL_DBL_EPSILON, &b);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val  = pre * b.val;
    result->err  = pre * b.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return status;
  }
  else if(GSL_ROOT4_DBL_EPSILON * x > (l*l + l + 1.0)) {
    gsl_sf_result b;
    int  status=  gsl_sf_bessel_Jnu_asympx_e(l + 0.5, x, &b);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val = pre * b.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + pre * b.err;
    return status;
  }
  else if(l > 1.0/GSL_ROOT6_DBL_EPSILON) {
    gsl_sf_result b;
    int  status=  gsl_sf_bessel_Jnu_asymp_Olver_e(l + 0.5, x, &b);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val = pre * b.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + pre * b.err;
    return status;
  }
  else if(x > 1000.0)
  {
    /* We need this to avoid feeding large x to CF1; note that
     * due to the above check, we know that n <= 50.
     */
    gsl_sf_result b;
    int  status=  gsl_sf_bessel_Jnu_asympx_e(l + 0.5, x, &b);
    MpIeee pre=  sqrt((MpIeee( "0.5" )*M_PI)/x);
    result->val = pre * b.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + pre * b.err;
    return status;  
  }
  else {
    MpIeee sgn;
    MpIeee ratio;
    int  stat_CF1=  gsl_sf_bessel_J_CF1(l+0.5, x, &ratio, &sgn);
    MpIeee jellp1=  GSL_SQRT_DBL_EPSILON * ratio;
    MpIeee jell=  GSL_SQRT_DBL_EPSILON;
    MpIeee jellm1;
    int  ell;
    for(ell = l; ell > 0; ell--) {
      jellm1 = -jellp1 + (MpIeee( "2" )*ell + MpIeee( "1" ))/x * jell;
      jellp1 = jell;
      jell   = jellm1;
    }

    if(fabs(jell) > fabs(jellp1)) {
      gsl_sf_result j0_result;
      int  stat_j0=  gsl_sf_bessel_j0_e(x, &j0_result);
      MpIeee pre=  GSL_SQRT_DBL_EPSILON / jell;
      result->val  = j0_result.val * pre;
      result->err  = j0_result.err * fabs(pre);
      result->err += 2.0 * GSL_DBL_EPSILON * (0.5*l + 1.0) * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_j0, stat_CF1);
    }
    else {
      gsl_sf_result j1_result;
      int  stat_j1=  gsl_sf_bessel_j1_e(x, &j1_result);
      MpIeee pre=  GSL_SQRT_DBL_EPSILON / jellp1;
      result->val  = j1_result.val * pre;
      result->err  = j1_result.err * fabs(pre);
      result->err += 2.0 * GSL_DBL_EPSILON * (0.5*l + 1.0) * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_j1, stat_CF1);
    }
  }
}


int
 gsl_sf_bessel_jl_array(const int lmax, const MpIeee x, MpIeee * result_array)
{
  /* CHECK_POINTER(result_array) */

  if(lmax < 0 || x < 0.0) {
    int  j;
    for(j=0; j<=lmax; j++) result_array[j] = MpIeee( "0.0" );
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(x == 0.0) {
    int  j;
    for(j=1; j<=lmax; j++) result_array[j] = MpIeee( "0.0" );
    result_array[0] = MpIeee( "1.0" );
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result r_jellp1;
    gsl_sf_result r_jell;
    int  stat_0=  gsl_sf_bessel_jl_e(lmax+1, x, &r_jellp1);
    int  stat_1=  gsl_sf_bessel_jl_e(lmax,   x, &r_jell);
    MpIeee jellp1=  r_jellp1.val;
    MpIeee jell=  r_jell.val;
    MpIeee jellm1;
    int  ell;

    result_array[lmax] = jell;
    for(ell = lmax; ell >= 1; ell--) {
      jellm1 = -jellp1 + (MpIeee( "2" )*ell + MpIeee( "1" ))/x * jell;
      jellp1 = jell;
      jell   = jellm1;
      result_array[ell-1] = jellm1;
    }

    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
}


int  gsl_sf_bessel_jl_steed_array(const int lmax, const MpIeee x, MpIeee * jl_x)
{
  /* CHECK_POINTER(jl_x) */

  if(lmax < 0 || x < 0.0) {
    int  j;
    for(j=0; j<=lmax; j++) jl_x[j] = MpIeee( "0.0" );
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(x == 0.0) {
    int  j;
    for(j=1; j<=lmax; j++) jl_x[j] = MpIeee( "0.0" );
    jl_x[0] = MpIeee( "1.0" );
    return GSL_SUCCESS;
  }
  else if(x < 2.0*GSL_ROOT4_DBL_EPSILON) {
    /* first two terms of Taylor series */
    MpIeee inv_fact=  MpIeee( "1.0" );  /* 1/(1 3 5 ... (2l+1)) */
    MpIeee x_l=  MpIeee( "1.0" );  /* x^l */
    int  l;
    for(l=0; l<=lmax; l++) {
      jl_x[l]  = x_l * inv_fact;
      jl_x[l] *= MpIeee( "1.0" ) - MpIeee( "0.5" )*x*x/(MpIeee( "2.0" )*l+MpIeee( "3.0" ));
      inv_fact /= MpIeee( "2.0" )*l+MpIeee( "3.0" );
      x_l      *= x;
    }
    return GSL_SUCCESS;
  }
  else {
    /* Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)] */
    MpIeee x_inv=  MpIeee( "1.0" )/x;
    MpIeee W=  MpIeee( "2.0" )*x_inv;
    MpIeee F=  MpIeee( "1.0" );
    MpIeee FP=  (lmax+MpIeee( "1.0" )) * x_inv;
    MpIeee B=  MpIeee( "2.0" )*FP + x_inv;
    MpIeee end=  B + MpIeee( "20000.0" )*W;
    MpIeee D=  MpIeee( "1.0" )/B;
    MpIeee del=  -D;
    
    FP += del;
    
    /* continued fraction */
    do {
      B += W;
      D = MpIeee( "1.0" )/(B-D);
      del *= (B*D - MpIeee( "1." ));
      FP += del;
      if(D < MpIeee( "0.0" )) F = -F;
      if(B > end) {
        GSL_ERROR ("error", GSL_EMAXITER);
      }
    }
    while(fabs(del) >= fabs(FP) * GSL_DBL_EPSILON);
    
    FP *= F;
    
    if(lmax > 0) {
      /* downward recursion */
      MpIeee XP2=  FP;
      MpIeee PL=  lmax * x_inv;
      int  L=  lmax;
      int  LP;
      jl_x[lmax] = F;
      for(LP = 1; LP<=lmax; LP++) {
        jl_x[L-1] = PL * jl_x[L] + XP2;
        FP = PL*jl_x[L-1] - jl_x[L];
        XP2 = FP;
        PL -= x_inv;
        --L;
      }
      F = jl_x[0];
    }
    
    /* normalization */
    W = x_inv / sqrt(FP*FP + F*F);
    jl_x[0] = W*F;
    if(lmax > 0) {
      int  L;
      for(L=1; L<=lmax; L++) {
        jl_x[L] *= W;
      }
    }

    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_bessel_j0(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_j0_e(x, &result));
}

MpIeee gsl_sf_bessel_j1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_j1_e(x, &result));
}

MpIeee gsl_sf_bessel_j2(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_j2_e(x, &result));
}

MpIeee gsl_sf_bessel_jl(const int l, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_bessel_jl_e(l, x, &result));
}

