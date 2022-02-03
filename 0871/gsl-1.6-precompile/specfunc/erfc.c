#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/erfc.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 Gerard Jungman
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

/* Author:  J. Theiler (modifications by G. Jungman) */

/*
 * See Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)
 * (This applies only to the erfc8 stuff, which is the part
 *  of the original code that survives. I have replaced much of
 *  the other stuff with Chebyshev fits. These are simpler and
 *  more precise than the original approximations. [GJ])
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_erf.h>

#include "check.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#define LogRootPi_  0.57236494292470008706


static MpIeee erfc8_sum(MpIeee x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static MpIeee P[] =  {
      MpIeee( "2.97886562639399288862" ),
      MpIeee( "7.409740605964741794425" ),
      MpIeee( "6.1602098531096305440906" ),
      MpIeee( "5.019049726784267463450058" ),
      MpIeee( "1.275366644729965952479585264" ),
      MpIeee( "0.5641895835477550741253201704" )
  };
  static MpIeee Q[] =  {
      MpIeee( "3.3690752069827527677" ),
      MpIeee( "9.608965327192787870698" ),
      MpIeee( "17.08144074746600431571095" ),
      MpIeee( "12.0489519278551290360340491" ),
      MpIeee( "9.396034016235054150430579648" ),
      MpIeee( "2.260528520767326969591866945" ),
      MpIeee( "1.0" )
  };
  MpIeee num= MpIeee( "0.0" );MpIeee  den= MpIeee( "0.0" );
  int  i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}

inline
static MpIeee erfc8(MpIeee x)
{
  MpIeee e;
  e = erfc8_sum(x);
  e *= exp(-x*x);
  return e;
}

inline
static MpIeee log_erfc8(MpIeee x)
{
  MpIeee e;
  e = erfc8_sum(x);
  e = log(e) - x*x;
  return e;
}

#if 0
/* Abramowitz+Stegun, 7.2.14 */
static MpIeee erfcasympsum(MpIeee x)
{
  int  i;
  MpIeee e=  MpIeee( "1." );
  MpIeee coef=  MpIeee( "1." );
  for (i=1; i<5; ++i) {
    /* coef *= -(2*i-1)/(2*x*x); ??? [GJ] */
    coef *= -(MpIeee( "2" )*i+MpIeee( "1" ))/(i*(MpIeee( "4" )*x*x*x*x));
    e += coef;
    /*
    if (fabs(coef) < 1.0e-15) break;
    if (fabs(coef) > 1.0e10) break;
    
    [GJ]: These tests are not useful. This function is only
    used below. Took them out; they gum up the pipeline.
    */
  }
  return e;
}
#endif /* 0 */


/* Abramowitz+Stegun, 7.1.5 */
static int  erfseries(MpIeee x, gsl_sf_result * result)
{
  MpIeee coef=  x;
  MpIeee e=  coef;
  MpIeee del;
  int  k;
  for (k=1; k<30; ++k) {
    coef *= -x*x/k;
    del   = coef/(MpIeee( "2.0" )*k+MpIeee( "1.0" ));
    e += del;
  }
  result->val = 2.0 / M_SQRTPI * e;
  result->err = 2.0 / M_SQRTPI * (fabs(del) + GSL_DBL_EPSILON);
  return GSL_SUCCESS;
}


/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */
static MpIeee erfc_xlt1_data[20] =  {
  MpIeee( "1.06073416421769980345174155056" ),
 -MpIeee( "0.42582445804381043569204735291" ),
  MpIeee( "0.04955262679620434040357683080" ),
  MpIeee( "0.00449293488768382749558001242" ),
 -MpIeee( "0.00129194104658496953494224761" ),
 -MpIeee( "0.00001836389292149396270416979" ),
  MpIeee( "0.00002211114704099526291538556" ),
 -MpIeee( "5.23337485234257134673693179020e-7" ),
 -MpIeee( "2.78184788833537885382530989578e-7" ),
  MpIeee( "1.41158092748813114560316684249e-8" ),
  MpIeee( "2.72571296330561699984539141865e-9" ),
 -MpIeee( "2.06343904872070629406401492476e-10" ),
 -MpIeee( "2.14273991996785367924201401812e-11" ),
  MpIeee( "2.22990255539358204580285098119e-12" ),
  MpIeee( "1.36250074650698280575807934155e-13" ),
 -MpIeee( "1.95144010922293091898995913038e-14" ),
 -MpIeee( "6.85627169231704599442806370690e-16" ),
  MpIeee( "1.44506492869699938239521607493e-16" ),
  MpIeee( "2.45935306460536488037576200030e-18" ),
 -MpIeee( "9.29599561220523396007359328540e-19" )
};
static cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  12
};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static MpIeee erfc_x15_data[25] =  {
  MpIeee( "0.44045832024338111077637466616" ),
 -MpIeee( "0.143958836762168335790826895326" ),
  MpIeee( "0.044786499817939267247056666937" ),
 -MpIeee( "0.013343124200271211203618353102" ),
  MpIeee( "0.003824682739750469767692372556" ),
 -MpIeee( "0.001058699227195126547306482530" ),
  MpIeee( "0.000283859419210073742736310108" ),
 -MpIeee( "0.000073906170662206760483959432" ),
  MpIeee( "0.000018725312521489179015872934" ),
 -MpIeee( "4.62530981164919445131297264430e-6" ),
  MpIeee( "1.11558657244432857487884006422e-6" ),
 -MpIeee( "2.63098662650834130067808832725e-7" ),
  MpIeee( "6.07462122724551777372119408710e-8" ),
 -MpIeee( "1.37460865539865444777251011793e-8" ),
  MpIeee( "3.05157051905475145520096717210e-9" ),
 -MpIeee( "6.65174789720310713757307724790e-10" ),
  MpIeee( "1.42483346273207784489792999706e-10" ),
 -MpIeee( "3.00141127395323902092018744545e-11" ),
  MpIeee( "6.22171792645348091472914001250e-12" ),
 -MpIeee( "1.26994639225668496876152836555e-12" ),
  MpIeee( "2.55385883033257575402681845385e-13" ),
 -MpIeee( "5.06258237507038698392265499770e-14" ),
  MpIeee( "9.89705409478327321641264227110e-15" ),
 -MpIeee( "1.90685978789192181051961024995e-15" ),
  MpIeee( "3.50826648032737849245113757340e-16" )
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static MpIeee erfc_x510_data[20] =  {
  MpIeee( "1.11684990123545698684297865808" ),
  MpIeee( "0.003736240359381998520654927536" ),
 -MpIeee( "0.000916623948045470238763619870" ),
  MpIeee( "0.000199094325044940833965078819" ),
 -MpIeee( "0.000040276384918650072591781859" ),
  MpIeee( "7.76515264697061049477127605790e-6" ),
 -MpIeee( "1.44464794206689070402099225301e-6" ),
  MpIeee( "2.61311930343463958393485241947e-7" ),
 -MpIeee( "4.61833026634844152345304095560e-8" ),
  MpIeee( "8.00253111512943601598732144340e-9" ),
 -MpIeee( "1.36291114862793031395712122089e-9" ),
  MpIeee( "2.28570483090160869607683087722e-10" ),
 -MpIeee( "3.78022521563251805044056974560e-11" ),
  MpIeee( "6.17253683874528285729910462130e-12" ),
 -MpIeee( "9.96019290955316888445830597430e-13" ),
  MpIeee( "1.58953143706980770269506726000e-13" ),
 -MpIeee( "2.51045971047162509999527428316e-14" ),
  MpIeee( "3.92607828989125810013581287560e-15" ),
 -MpIeee( "6.07970619384160374392535453420e-16" ),
  MpIeee( "9.12600607264794717315507477670e-17" )
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

#if 0
inline
static MpIeee erfc_asymptotic(MpIeee x)
{
  return exp(-x*x)/x * erfcasympsum(x) / M_SQRTPI;
}
inline
static MpIeee log_erfc_asymptotic(MpIeee x)
{
  return log(erfcasympsum(x)/x) - x*x - LogRootPi_;
}
#endif /* 0 */


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_erfc_e(MpIeee x, gsl_sf_result * result)
{
  const MpIeee ax=  fabs(x);
  MpIeee e_val;MpIeee  e_err;

  /* CHECK_POINTER(result) */

  if(ax <= 1.0) {
    MpIeee t=  MpIeee( "2.0" )*ax - MpIeee( "1.0" );
    gsl_sf_result c;
    cheb_eval_e(&erfc_xlt1_cs, t, &c);
    e_val = c.val;
    e_err = c.err;
  }
  else if(ax <= 5.0) {
    MpIeee ex2=  exp(-x*x);
    MpIeee t=  MpIeee( "0.5" )*(ax-MpIeee( "3.0" ));
    gsl_sf_result c;
    cheb_eval_e(&erfc_x15_cs, t, &c);
    e_val = ex2 * c.val;
    e_err = ex2 * (c.err + MpIeee( "2.0" )*fabs(x)*GSL_DBL_EPSILON);
  }
  else if(ax < 10.0) {
    MpIeee exterm=  exp(-x*x) / ax;
    MpIeee t=  (MpIeee( "2.0" )*ax - MpIeee( "15.0" ))/MpIeee( "5.0" );
    gsl_sf_result c;
    cheb_eval_e(&erfc_x510_cs, t, &c);
    e_val = exterm * c.val;
    e_err = exterm * (c.err + MpIeee( "2.0" )*fabs(x)*GSL_DBL_EPSILON + GSL_DBL_EPSILON);
  }
  else {
    e_val = erfc8(ax);
    e_err = (x*x + MpIeee( "1.0" )) * GSL_DBL_EPSILON * fabs(e_val);
  }

  if(x < MpIeee( "0.0" )) {
    result->val  = 2.0 - e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    result->val  = e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}


int  gsl_sf_log_erfc_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x*x < MpIeee( "10.0" )*GSL_ROOT6_DBL_EPSILON) {
    const MpIeee y=  x / M_SQRTPI;
    /* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
    const MpIeee c3=  (4.0 - M_PI)/3.0;
    const MpIeee c4=  2.0*(1.0 - M_PI/3.0);
    const MpIeee c5=  -0.001829764677455021;  /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
    const MpIeee c6=   0.02629651521057465;   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
    const MpIeee c7=  -0.01621575378835404;
    const MpIeee c8=   0.00125993961762116;
    const MpIeee c9=   0.00556964649138;
    const MpIeee c10=  -0.0045563339802;
    const MpIeee c11=   0.0009461589032;
    const MpIeee c12=   0.0013200243174;
    const MpIeee c13=  -0.00142906;
    const MpIeee c14=   0.00048204;
    MpIeee series=  c8 + y*(c9 + y*(c10 + y*(c11 + y*(c12 + y*(c13 + c14*y)))));
    series = y*(MpIeee( "1.0" ) + y*(MpIeee( "1.0" ) + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*series)))))));
    result->val = -2.0 * series;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  /*
  don't like use of log1p(); added above series stuff for small x instead, should be ok [GJ]
  else if (fabs(x) < 1.0) {
    gsl_sf_result result_erf;
    gsl_sf_erf_e(x, &result_erf);
    result->val  = log1p(-result_erf.val);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  */
  else if(x > MpIeee( "8.0" )) {
    result->val = log_erfc8(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result result_erfc;
    gsl_sf_erfc_e(x, &result_erfc);
    result->val  = log(result_erfc.val);
    result->err  = fabs(result_erfc.err / result_erfc.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_erf_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    return erfseries(x, result);
  }
  else {
    gsl_sf_result result_erfc;
    gsl_sf_erfc_e(x, &result_erfc);
    result->val  = 1.0 - result_erfc.val;
    result->err  = result_erfc.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_erf_Z_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const MpIeee ex2=  exp(-x*x/2.0);
    result->val  = ex2 / (M_SQRT2 * M_SQRTPI);
    result->err  = fabs(x * result->val) * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_erf_Q_e(MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    gsl_sf_result result_erfc;
    int  stat=  gsl_sf_erfc_e(x/M_SQRT2, &result_erfc);
    result->val  = 0.5 * result_erfc.val;
    result->err  = 0.5 * result_erfc.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
}


int  gsl_sf_hazard_e(MpIeee x, gsl_sf_result * result)
{
  if(x < MpIeee( "25.0" ))
  {
    gsl_sf_result result_ln_erfc;
    const int stat_l = gsl_sf_log_erfc_e(x/M_SQRT2, &result_ln_erfc);
    const MpIeee lnc=  -0.22579135264472743236; /* ln(sqrt(2/pi)) */
    const MpIeee arg=  lnc - 0.5*x*x - result_ln_erfc.val;
    const int stat_e = gsl_sf_exp_e(arg, result);
    result->err += 3.0 * (1.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
    result->err += fabs(result_ln_erfc.err * result->val);
    return GSL_ERROR_SELECT_2(stat_l, stat_e);
  }
  else
  {
    const MpIeee ix2=  1.0/(x*x);
    const MpIeee corrB=  1.0 - 9.0*ix2 * (1.0 - 11.0*ix2);
    const MpIeee corrM=  1.0 - 5.0*ix2 * (1.0 - 7.0*ix2 * corrB);
    const MpIeee corrT=  1.0 - ix2 * (1.0 - 3.0*ix2*corrM);
    result->val = x / corrT;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_erfc(MpIeee x)
{
  EVAL_RESULT(gsl_sf_erfc_e(x, &result));
}

MpIeee gsl_sf_log_erfc(MpIeee x)
{
  EVAL_RESULT(gsl_sf_log_erfc_e(x, &result));
}

MpIeee gsl_sf_erf(MpIeee x)
{
  EVAL_RESULT(gsl_sf_erf_e(x, &result));
}

MpIeee gsl_sf_erf_Z(MpIeee x)
{
  EVAL_RESULT(gsl_sf_erf_Z_e(x, &result));
}

MpIeee gsl_sf_erf_Q(MpIeee x)
{
  EVAL_RESULT(gsl_sf_erf_Q_e(x, &result));
}

MpIeee gsl_sf_hazard(MpIeee x)
{
  EVAL_RESULT(gsl_sf_hazard_e(x, &result));
}

