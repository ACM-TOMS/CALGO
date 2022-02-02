#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/ellint.c
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

/* Author: G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_sf_ellint.h>

#include "error.h"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

inline
static MpIeee locMAX3(MpIeee x, MpIeee y, MpIeee z)
{
  MpIeee xy=  GSL_MAX(x, y);
  return GSL_MAX(xy, z);
}

inline
static MpIeee locMAX4(MpIeee x, MpIeee y, MpIeee z, MpIeee w)
{
  MpIeee xy=  GSL_MAX(x,  y);
  MpIeee xyz=  GSL_MAX(xy, z);
  return GSL_MAX(xyz, w);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


/* based on Carlson's algorithms:
   [B. C. Carlson Numer. Math. 33, 1 (1979)]
   
   see also:
   [B.C. Carlson, Special Functions of Applied Mathematics (1977)]
 */

/* According to Carlson's algorithm, the errtol parameter
   typically effects the relative error in the following way:

   relative error < 16 errtol^6 / (1 - 2 errtol)

     errtol     precision
     ------     ----------
     0.001       1.0e-17
     0.003       2.0e-14 
     0.01        2.0e-11
     0.03        2.0e-8
     0.1         2.0e-5
*/


int
 gsl_sf_ellint_RC_e(MpIeee x, MpIeee y, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee lolim=  5.0 * GSL_DBL_MIN;
  const MpIeee uplim=  0.2 * GSL_DBL_MAX;
  const gsl_prec_t goal = GSL_MODE_PREC(mode);
  const MpIeee errtol=  ( goal == GSL_PREC_DOUBLE ? 0.001 : 0.03 );
  const MpIeee prec=  gsl_prec_eps[goal];

  if(x < MpIeee( "0.0" ) || y < MpIeee( "0.0" ) || x + y < lolim) {
    DOMAIN_ERROR(result);
  }
  else if(GSL_MAX(x, y) < uplim) { 
    const MpIeee c1=  1.0 / 7.0;
    const MpIeee c2=  9.0 / 22.0;
    MpIeee xn=  x;
    MpIeee yn=  y;
    MpIeee mu;MpIeee  sn;MpIeee  lamda;MpIeee  s;
    while(1) {
      mu = (xn + yn + yn) / MpIeee( "3.0" );
      sn = (yn + mu) / mu - MpIeee( "2.0" );
      if (fabs(sn) < errtol) break;
      lamda = MpIeee( "2.0" ) * sqrt(xn) * sqrt(yn) + yn;
      xn = (xn + lamda) * MpIeee( "0.25" );
      yn = (yn + lamda) * MpIeee( "0.25" );
    }
    s = sn * sn * (MpIeee( "0.3" ) + sn * (c1 + sn * (MpIeee( "0.375" ) + sn * c2)));
    result->val = (1.0 + s) / sqrt(mu);
    result->err = prec * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR(result);
  }
}


int
 gsl_sf_ellint_RD_e(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode, gsl_sf_result * result)
{
  const gsl_prec_t goal = GSL_MODE_PREC(mode);
  const MpIeee errtol=  ( goal == GSL_PREC_DOUBLE ? 0.001 : 0.03 );
  const MpIeee prec=  gsl_prec_eps[goal];
  const MpIeee lolim=  2.0/pow(GSL_DBL_MAX, 2.0/3.0);
  const MpIeee uplim=  pow(0.1*errtol/GSL_DBL_MIN, 2.0/3.0);

  if(GSL_MIN(x,y) < 0.0 || GSL_MIN(x+y,z) < lolim) {
    DOMAIN_ERROR(result);
  }
  else if(locMAX3(x,y,z) < uplim) {
    const MpIeee c1=  3.0 / 14.0;
    const MpIeee c2=  1.0 /  6.0;
    const MpIeee c3=  9.0 / 22.0;
    const MpIeee c4=  3.0 / 26.0;
    MpIeee xn=  x;
    MpIeee yn=  y;
    MpIeee zn=  z;
    MpIeee sigma=  MpIeee( "0.0" );
    MpIeee power4=  MpIeee( "1.0" );
    MpIeee ea;MpIeee  eb;MpIeee  ec;MpIeee  ed;MpIeee  ef;MpIeee  s1;MpIeee  s2;
    MpIeee mu;MpIeee  xndev;MpIeee  yndev;MpIeee  zndev;
    while(1) {
      MpIeee xnroot;MpIeee  ynroot;MpIeee  znroot;MpIeee  lamda;
      MpIeee epslon;
      mu = (xn + yn + MpIeee( "3.0" ) * zn) * MpIeee( "0.2" );
      xndev = (mu - xn) / mu;
      yndev = (mu - yn) / mu;
      zndev = (mu - zn) / mu;
      epslon = locMAX3(fabs(xndev), fabs(yndev), fabs(zndev));
      if (epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      sigma  += power4 / (znroot * (zn + lamda));
      power4 *= MpIeee( "0.25" );
      xn = (xn + lamda) * MpIeee( "0.25" );
      yn = (yn + lamda) * MpIeee( "0.25" );
      zn = (zn + lamda) * MpIeee( "0.25" );
    }
    ea = xndev * yndev;
    eb = zndev * zndev;
    ec = ea - eb;
    ed = ea - MpIeee( "6.0" ) * eb;
    ef = ed + ec + ec;
    s1 = ed * (- c1 + MpIeee( "0.25" ) * c3 * ed - MpIeee( "1.5" ) * c4 * zndev * ef);
    s2 = zndev * (c2 * ef + zndev * (- c3 * ec + zndev * c4 * ea));
    result->val = 3.0 * sigma + power4 * (1.0 + s1 + s2) / (mu * sqrt(mu));
    result->err = prec * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR(result);
  }
}


int
 gsl_sf_ellint_RF_e(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee lolim=  5.0 * GSL_DBL_MIN;
  const MpIeee uplim=  0.2 * GSL_DBL_MAX;
  const gsl_prec_t goal = GSL_MODE_PREC(mode);
  const MpIeee errtol=  ( goal == GSL_PREC_DOUBLE ? 0.001 : 0.03 );
  const MpIeee prec=  gsl_prec_eps[goal];

  if(x < MpIeee( "0.0" ) || y < MpIeee( "0.0" ) || z < MpIeee( "0.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(x+y < lolim || x+z < lolim || y+z < lolim) {
    DOMAIN_ERROR(result);
  }
  else if(locMAX3(x,y,z) < uplim) { 
    const MpIeee c1=  1.0 / 24.0;
    const MpIeee c2=  3.0 / 44.0;
    const MpIeee c3=  1.0 / 14.0;
    MpIeee xn=  x;
    MpIeee yn=  y;
    MpIeee zn=  z;
    MpIeee mu;MpIeee  xndev;MpIeee  yndev;MpIeee  zndev;MpIeee  e2;MpIeee  e3;MpIeee  s;
    while(1) {
      MpIeee epslon;MpIeee  lamda;
      MpIeee xnroot;MpIeee  ynroot;MpIeee  znroot;
      mu = (xn + yn + zn) / MpIeee( "3.0" );
      xndev = MpIeee( "2.0" ) - (mu + xn) / mu;
      yndev = MpIeee( "2.0" ) - (mu + yn) / mu;
      zndev = MpIeee( "2.0" ) - (mu + zn) / mu;
      epslon = locMAX3(fabs(xndev), fabs(yndev), fabs(zndev));
      if (epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      xn = (xn + lamda) * MpIeee( "0.25" );
      yn = (yn + lamda) * MpIeee( "0.25" );
      zn = (zn + lamda) * MpIeee( "0.25" );
    }
    e2 = xndev * yndev - zndev * zndev;
    e3 = xndev * yndev * zndev;
    s = MpIeee( "1.0" ) + (c1 * e2 - MpIeee( "0.1" ) - c2 * e3) * e2 + c3 * e3;
    result->val = s / sqrt(mu);
    result->err = prec * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR(result);
  }
}


int
 gsl_sf_ellint_RJ_e(MpIeee x, MpIeee y, MpIeee z, MpIeee p, gsl_mode_t mode, gsl_sf_result * result)
{
  const gsl_prec_t goal = GSL_MODE_PREC(mode);
  const MpIeee errtol=  ( goal == GSL_PREC_DOUBLE ? 0.001 : 0.03 );
  const MpIeee prec=  gsl_prec_eps[goal];
  const MpIeee lolim=        pow(5.0 * GSL_DBL_MIN, 1.0/3.0);
  const MpIeee uplim=  0.3 * pow(0.2 * GSL_DBL_MAX, 1.0/3.0);

  if(x < MpIeee( "0.0" ) || y < MpIeee( "0.0" ) || z < MpIeee( "0.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(x + y < lolim || x + z < lolim || y + z < lolim || p < lolim) {
    DOMAIN_ERROR(result);
  }
  else if(locMAX4(x,y,z,p) < uplim) {
    const MpIeee c1=  3.0 / 14.0;
    const MpIeee c2=  1.0 /  3.0;
    const MpIeee c3=  3.0 / 22.0;
    const MpIeee c4=  3.0 / 26.0;
    MpIeee xn=  x;
    MpIeee yn=  y;
    MpIeee zn=  z;
    MpIeee pn=  p;
    MpIeee sigma=  MpIeee( "0.0" );
    MpIeee power4=  MpIeee( "1.0" );
    MpIeee mu;MpIeee  xndev;MpIeee  yndev;MpIeee  zndev;MpIeee  pndev;
    MpIeee ea;MpIeee  eb;MpIeee  ec;MpIeee  e2;MpIeee  e3;MpIeee  s1;MpIeee  s2;MpIeee  s3;
    while(1) {
      MpIeee xnroot;MpIeee  ynroot;MpIeee  znroot;
      MpIeee lamda;MpIeee  alfa;MpIeee  beta;
      MpIeee epslon;
      gsl_sf_result rcresult;
      int  rcstatus;
      mu = (xn + yn + zn + pn + pn) * MpIeee( "0.2" );
      xndev = (mu - xn) / mu;
      yndev = (mu - yn) / mu;
      zndev = (mu - zn) / mu;
      pndev = (mu - pn) / mu;
      epslon = locMAX4(fabs(xndev), fabs(yndev), fabs(zndev), fabs(pndev));
      if(epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      alfa = pn * (xnroot + ynroot + znroot) + xnroot * ynroot * znroot;
      alfa = alfa * alfa;
      beta = pn * (pn + lamda) * (pn + lamda);
      rcstatus = gsl_sf_ellint_RC_e(alfa, beta, mode, &rcresult);
      if(rcstatus != GSL_SUCCESS) {
        result->val = 0.0;
        result->err = 0.0;
        return rcstatus;
      }
      sigma  += power4 * rcresult.val;
      power4 *= MpIeee( "0.25" );
      xn = (xn + lamda) * MpIeee( "0.25" );
      yn = (yn + lamda) * MpIeee( "0.25" );
      zn = (zn + lamda) * MpIeee( "0.25" );
      pn = (pn + lamda) * MpIeee( "0.25" );
    }
    
    ea = xndev * (yndev + zndev) + yndev * zndev;
    eb = xndev * yndev * zndev;
    ec = pndev * pndev;
    e2 = ea - MpIeee( "3.0" ) * ec;
    e3 = eb + MpIeee( "2.0" ) * pndev * (ea - ec);
    s1 = MpIeee( "1.0" ) + e2 * (- c1 + MpIeee( "0.75" ) * c3 * e2 - MpIeee( "1.5" ) * c4 * e3);
    s2 = eb * (MpIeee( "0.5" ) * c2 + pndev * (- c3 - c3 + pndev * c4));
    s3 = pndev * ea * (c2 - pndev * c3) - c2 * pndev * ec;
    result->val = 3.0 * sigma + power4 * (s1 + s2 + s3) / (mu * sqrt(mu));
    result->err = prec * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR(result);
  }
}


/* [Carlson, Numer. Math. 33 (1979) 1, (4.1)] */
int
 gsl_sf_ellint_F_e(MpIeee phi, MpIeee k, gsl_mode_t mode, gsl_sf_result * result)
{
  MpIeee sin_phi=  sin(phi);
  MpIeee sin2_phi=  sin_phi*sin_phi;
  MpIeee x=  MpIeee( "1.0" ) - sin2_phi;
  MpIeee y=  MpIeee( "1.0" ) - k*k*sin2_phi;
  gsl_sf_result rf;
  int  status=  gsl_sf_ellint_RF_e(x, y, 1.0, mode, &rf);
  result->val = sin_phi * rf.val;
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(sin_phi*rf.err);
  return status;
}


/* [Carlson, Numer. Math. 33 (1979) 1, (4.2)] */
int
 gsl_sf_ellint_E_e(MpIeee phi, MpIeee k, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee sin_phi=  sin(phi);
  const MpIeee sin2_phi=  sin_phi  * sin_phi;
  const MpIeee x=  1.0 - sin2_phi;
  const MpIeee y=  1.0 - k*k*sin2_phi;
  if(x < GSL_DBL_EPSILON) {
    return gsl_sf_ellint_Ecomp_e(k, mode, result);
  }
  else {
    gsl_sf_result rf;
    gsl_sf_result rd;
    const MpIeee sin3_phi=  sin2_phi * sin_phi;
    const int rfstatus = gsl_sf_ellint_RF_e(x, y, 1.0, mode, &rf);
    const int rdstatus = gsl_sf_ellint_RD_e(x, y, 1.0, mode, &rd);
    result->val  = sin_phi * rf.val - k*k/3.0 * sin3_phi * rd.val;
    result->err  = GSL_DBL_EPSILON * fabs(sin_phi * rf.val);
    result->err += fabs(sin_phi*rf.err);
    result->err += k*k/3.0 * GSL_DBL_EPSILON * fabs(sin3_phi * rd.val);
    result->err += k*k/3.0 * fabs(sin3_phi*rd.err);
    return GSL_ERROR_SELECT_2(rfstatus, rdstatus);
  }
}


/* [Carlson, Numer. Math. 33 (1979) 1, (4.3)] */
int
 gsl_sf_ellint_P_e(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee sin_phi=  sin(phi);
  const MpIeee sin2_phi=  sin_phi  * sin_phi;
  const MpIeee sin3_phi=  sin2_phi * sin_phi;
  const MpIeee x=  1.0 - sin2_phi;
  const MpIeee y=  1.0 - k*k*sin2_phi;
  gsl_sf_result rf;
  gsl_sf_result rj;
  const int rfstatus = gsl_sf_ellint_RF_e(x, y, 1.0, mode, &rf);
  const int rjstatus = gsl_sf_ellint_RJ_e(x, y, 1.0, 1.0 + n*sin2_phi, mode, &rj);
  result->val  = sin_phi * rf.val - n/3.0*sin3_phi * rj.val;
  result->err  = GSL_DBL_EPSILON * fabs(sin_phi * rf.val);
  result->err += n/3.0 * fabs(sin3_phi*rj.err);
  return GSL_ERROR_SELECT_2(rfstatus, rjstatus);
}


/* [Carlson, Numer. Math. 33 (1979) 1, (4.4)] */
int
 gsl_sf_ellint_D_e(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee sin_phi=  sin(phi);
  const MpIeee sin2_phi=  sin_phi  * sin_phi;
  const MpIeee sin3_phi=  sin2_phi * sin_phi;
  const MpIeee x=  1.0 - sin2_phi;
  const MpIeee y=  1.0 - k*k*sin2_phi;
  gsl_sf_result rd;
  const int status = gsl_sf_ellint_RD_e(x, y, 1.0, mode, &rd);
  result->val = sin3_phi/3.0 * rd.val;
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(sin3_phi/3.0 * rd.err);
  return status;
}


/* [Carlson, Numer. Math. 33 (1979) 1, (4.5)] */
int
 gsl_sf_ellint_Kcomp_e(MpIeee k, gsl_mode_t mode, gsl_sf_result * result)
{
  if(k*k >= MpIeee( "1.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(k*k >= MpIeee( "1.0" ) - GSL_SQRT_DBL_EPSILON) {
    /* [Abramowitz+Stegun, 17.3.33] */
    const MpIeee y=  1.0 - k*k;
    const MpIeee a[] =  { 1.38629436112, 0.09666344259, 0.03590092383 };
    const MpIeee b[] =  { 0.5, 0.12498593597, 0.06880248576 };
    const MpIeee ta=  a[0] + y*(a[1] + y*a[2]);
    const MpIeee tb=  -log(y) * (b[0] * y*(b[1] + y*b[2]));
    result->val = ta + tb;
    result->err = 2.0 * GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    /* This was previously computed as,

         return gsl_sf_ellint_RF_e(0.0, 1.0 - k*k, 1.0, mode, result);

       but this underestimated the total error for small k, since the 
       argument y=1-k^2 is not exact (there is an absolute error of
       GSL_DBL_EPSILON near y=0 due to cancellation in the subtraction).
       Taking the singular behavior of -log(y) above gives an error
       of 0.5*epsilon/y near y=0. (BJG) */

    MpIeee y=  MpIeee( "1.0" ) - k*k;
    int  status=  gsl_sf_ellint_RF_e(0.0, y, 1.0, mode, result);
    result->err += 0.5 * GSL_DBL_EPSILON / y;
    return status ;
  }
}


/* [Carlson, Numer. Math. 33 (1979) 1, (4.6)] */
int
 gsl_sf_ellint_Ecomp_e(MpIeee k, gsl_mode_t mode, gsl_sf_result * result)
{
  if(k*k >= MpIeee( "1.0" )) {
    DOMAIN_ERROR(result);
  }
  else if(k*k >= MpIeee( "1.0" ) - GSL_SQRT_DBL_EPSILON) {
    /* [Abramowitz+Stegun, 17.3.36] */
    const MpIeee y=  1.0 - k*k;
    const MpIeee a[] =  { 0.44325141463, 0.06260601220, 0.04757383546 };
    const MpIeee b[] =  { 0.24998368310, 0.09200180037, 0.04069697526 };
    const MpIeee ta=  1.0 + y*(a[0] + y*(a[1] + a[2]*y));
    const MpIeee tb=  -y*log(y) * (b[0] + y*(b[1] + b[2]*y));
    result->val = ta + tb;
    result->err = 2.0 * GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result rf;
    gsl_sf_result rd;
    const MpIeee y=  1.0 - k*k;
    const int rfstatus = gsl_sf_ellint_RF_e(0.0, y, 1.0, mode, &rf);
    const int rdstatus = gsl_sf_ellint_RD_e(0.0, y, 1.0, mode, &rd);
    result->val = rf.val - k*k/3.0 * rd.val;
    result->err = rf.err + k*k/3.0 * rd.err;
    return GSL_ERROR_SELECT_2(rfstatus, rdstatus);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_ellint_Kcomp(MpIeee k, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_Kcomp_e(k, mode, &result));
}

MpIeee gsl_sf_ellint_Ecomp(MpIeee k, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_Ecomp_e(k, mode, &result));
}

MpIeee gsl_sf_ellint_F(MpIeee phi, MpIeee k, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_F_e(phi, k, mode, &result));
}

MpIeee gsl_sf_ellint_E(MpIeee phi, MpIeee k, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_E_e(phi, k, mode, &result));
}

MpIeee gsl_sf_ellint_P(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_P_e(phi, k, n, mode, &result));
}

MpIeee gsl_sf_ellint_D(MpIeee phi, MpIeee k, MpIeee n, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_D_e(phi, k, n, mode, &result));
}

MpIeee gsl_sf_ellint_RC(MpIeee x, MpIeee y, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_RC_e(x, y, mode, &result));
}

MpIeee gsl_sf_ellint_RD(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_RD_e(x, y, z, mode, &result));
}

MpIeee gsl_sf_ellint_RF(MpIeee x, MpIeee y, MpIeee z, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_RF_e(x, y, z, mode, &result));
}

MpIeee gsl_sf_ellint_RJ(MpIeee x, MpIeee y, MpIeee z, MpIeee p, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_ellint_RJ_e(x, y, z, p, mode, &result));
}
