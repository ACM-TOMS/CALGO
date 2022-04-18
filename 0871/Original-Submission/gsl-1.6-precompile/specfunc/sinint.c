#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/sinint.c
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
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_expint.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC r9sifg.f, W. Fullerton */

/*
 series for f1   on the interval  2.00000e-02 to  6.25000e-02
                                        with weighted error   2.82e-17
                                         log weighted error  16.55
                               significant figures required  15.36
                                    decimal places required  17.20
*/
static MpIeee f1_data[20] =  {
   -MpIeee( "0.1191081969051363610" ),
   -MpIeee( "0.0247823144996236248" ),
    MpIeee( "0.0011910281453357821" ),
   -MpIeee( "0.0000927027714388562" ),
    MpIeee( "0.0000093373141568271" ),
   -MpIeee( "0.0000011058287820557" ),
    MpIeee( "0.0000001464772071460" ),
   -MpIeee( "0.0000000210694496288" ),
    MpIeee( "0.0000000032293492367" ),
   -MpIeee( "0.0000000005206529618" ),
    MpIeee( "0.0000000000874878885" ),
   -MpIeee( "0.0000000000152176187" ),
    MpIeee( "0.0000000000027257192" ),
   -MpIeee( "0.0000000000005007053" ),
    MpIeee( "0.0000000000000940241" ),
   -MpIeee( "0.0000000000000180014" ),
    MpIeee( "0.0000000000000035063" ),
   -MpIeee( "0.0000000000000006935" ),
    MpIeee( "0.0000000000000001391" ),
   -MpIeee( "0.0000000000000000282" )
};
static cheb_series f1_cs = {
  f1_data,
  19,
  -1, 1,
  10
};

/*

 series for f2   on the interval  0.00000e+00 to  2.00000e-02
                                        with weighted error   4.32e-17
                                         log weighted error  16.36
                               significant figures required  14.75
                                    decimal places required  17.10
*/
static MpIeee f2_data[29] =  {
   -MpIeee( "0.0348409253897013234" ),
   -MpIeee( "0.0166842205677959686" ),
    MpIeee( "0.0006752901241237738" ),
   -MpIeee( "0.0000535066622544701" ),
    MpIeee( "0.0000062693421779007" ),
   -MpIeee( "0.0000009526638801991" ),
    MpIeee( "0.0000001745629224251" ),
   -MpIeee( "0.0000000368795403065" ),
    MpIeee( "0.0000000087202677705" ),
   -MpIeee( "0.0000000022601970392" ),
    MpIeee( "0.0000000006324624977" ),
   -MpIeee( "0.0000000001888911889" ),
    MpIeee( "0.0000000000596774674" ),
   -MpIeee( "0.0000000000198044313" ),
    MpIeee( "0.0000000000068641396" ),
   -MpIeee( "0.0000000000024731020" ),
    MpIeee( "0.0000000000009226360" ),
   -MpIeee( "0.0000000000003552364" ),
    MpIeee( "0.0000000000001407606" ),
   -MpIeee( "0.0000000000000572623" ),
    MpIeee( "0.0000000000000238654" ),
   -MpIeee( "0.0000000000000101714" ),
    MpIeee( "0.0000000000000044259" ),
   -MpIeee( "0.0000000000000019634" ),
    MpIeee( "0.0000000000000008868" ),
   -MpIeee( "0.0000000000000004074" ),
    MpIeee( "0.0000000000000001901" ),
   -MpIeee( "0.0000000000000000900" ),
    MpIeee( "0.0000000000000000432" )
};
static cheb_series f2_cs = {
  f2_data,
  28,
  -1, 1,
  14
};

/*

 series for g1   on the interval  2.00000e-02 to  6.25000e-02
                                        with weighted error   5.48e-17
                                         log weighted error  16.26
                               significant figures required  15.47
                                    decimal places required  16.92
*/
static MpIeee g1_data[21] =  {
   -MpIeee( "0.3040578798253495954" ),
   -MpIeee( "0.0566890984597120588" ),
    MpIeee( "0.0039046158173275644" ),
   -MpIeee( "0.0003746075959202261" ),
    MpIeee( "0.0000435431556559844" ),
   -MpIeee( "0.0000057417294453025" ),
    MpIeee( "0.0000008282552104503" ),
   -MpIeee( "0.0000001278245892595" ),
    MpIeee( "0.0000000207978352949" ),
   -MpIeee( "0.0000000035313205922" ),
    MpIeee( "0.0000000006210824236" ),
   -MpIeee( "0.0000000001125215474" ),
    MpIeee( "0.0000000000209088918" ),
   -MpIeee( "0.0000000000039715832" ),
    MpIeee( "0.0000000000007690431" ),
   -MpIeee( "0.0000000000001514697" ),
    MpIeee( "0.0000000000000302892" ),
   -MpIeee( "0.0000000000000061400" ),
    MpIeee( "0.0000000000000012601" ),
   -MpIeee( "0.0000000000000002615" ),
    MpIeee( "0.0000000000000000548" )
};
static cheb_series g1_cs = {
  g1_data,
  20,
  -1, 1,
  13
};

/*

 series for g2   on the interval  0.00000e+00 to  2.00000e-02
                                        with weighted error   5.01e-17
                                         log weighted error  16.30
                               significant figures required  15.12
                                    decimal places required  17.07
*/
static MpIeee g2_data[34] =  {
   -MpIeee( "0.0967329367532432218" ),
   -MpIeee( "0.0452077907957459871" ),
    MpIeee( "0.0028190005352706523" ),
   -MpIeee( "0.0002899167740759160" ),
    MpIeee( "0.0000407444664601121" ),
   -MpIeee( "0.0000071056382192354" ),
    MpIeee( "0.0000014534723163019" ),
   -MpIeee( "0.0000003364116512503" ),
    MpIeee( "0.0000000859774367886" ),
   -MpIeee( "0.0000000238437656302" ),
    MpIeee( "0.0000000070831906340" ),
   -MpIeee( "0.0000000022318068154" ),
    MpIeee( "0.0000000007401087359" ),
   -MpIeee( "0.0000000002567171162" ),
    MpIeee( "0.0000000000926707021" ),
   -MpIeee( "0.0000000000346693311" ),
    MpIeee( "0.0000000000133950573" ),
   -MpIeee( "0.0000000000053290754" ),
    MpIeee( "0.0000000000021775312" ),
   -MpIeee( "0.0000000000009118621" ),
    MpIeee( "0.0000000000003905864" ),
   -MpIeee( "0.0000000000001708459" ),
    MpIeee( "0.0000000000000762015" ),
   -MpIeee( "0.0000000000000346151" ),
    MpIeee( "0.0000000000000159996" ),
   -MpIeee( "0.0000000000000075213" ),
    MpIeee( "0.0000000000000035970" ),
   -MpIeee( "0.0000000000000017530" ),
    MpIeee( "0.0000000000000008738" ),
   -MpIeee( "0.0000000000000004487" ),
    MpIeee( "0.0000000000000002397" ),
   -MpIeee( "0.0000000000000001347" ),
    MpIeee( "0.0000000000000000801" ),
   -MpIeee( "0.0000000000000000501" )
};
static cheb_series g2_cs = {
  g2_data,
  33,
  -1, 1,
  20
};


/* x >= 4.0 */
static void fg_asymp(const MpIeee x, gsl_sf_result * f, gsl_sf_result * g)
{
  /*
      xbig = sqrt (1.0/r1mach(3))
      xmaxf = exp (amin1(-alog(r1mach(1)), alog(r1mach(2))) - 0.01)
      xmaxg = 1.0/sqrt(r1mach(1))
      xbnd = sqrt(50.0)
  */
  const MpIeee xbig=  1.0/GSL_SQRT_DBL_EPSILON;
  const MpIeee xmaxf=  1.0/GSL_DBL_MIN;
  const MpIeee xmaxg=  1.0/GSL_SQRT_DBL_MIN;
  const MpIeee xbnd=  7.07106781187;

  const MpIeee x2=  x*x;

  if(x <= xbnd) {
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_e(&f1_cs, (1.0/x2-0.04125)/0.02125, &result_c1);
    cheb_eval_e(&g1_cs, (1.0/x2-0.04125)/0.02125, &result_c2);
    f->val = (1.0 + result_c1.val)/x;
    g->val = (1.0 + result_c2.val)/x2;
    f->err = result_c1.err/x  + 2.0 * GSL_DBL_EPSILON * fabs(f->val);
    g->err = result_c2.err/x2 + 2.0 * GSL_DBL_EPSILON * fabs(g->val);
  }
  else if(x <= xbig) {
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_e(&f2_cs, 100.0/x2-1.0, &result_c1);
    cheb_eval_e(&g2_cs, 100.0/x2-1.0, &result_c2);
    f->val = (1.0 + result_c1.val)/x;
    g->val = (1.0 + result_c2.val)/x2;
    f->err = result_c1.err/x  + 2.0 * GSL_DBL_EPSILON * fabs(f->val);
    g->err = result_c2.err/x2 + 2.0 * GSL_DBL_EPSILON * fabs(g->val);
  }
  else {
    f->val = (x < xmaxf ? 1.0/x  : 0.0);
    g->val = (x < xmaxg ? 1.0/x2 : 0.0);
    f->err = 2.0 * GSL_DBL_EPSILON * fabs(f->val);
    g->err = 2.0 * GSL_DBL_EPSILON * fabs(g->val);
  }

  return;
}


/* based on SLATEC si.f, W. Fullerton

 series for si   on the interval  0.00000e+00 to  1.60000e+01
                                        with weighted error   1.22e-17
                                         log weighted error  16.91
                               significant figures required  16.37
                                    decimal places required  17.45
*/

static MpIeee si_data[12] =  {
  -MpIeee( "0.1315646598184841929" ),
  -MpIeee( "0.2776578526973601892" ),
   MpIeee( "0.0354414054866659180" ),
  -MpIeee( "0.0025631631447933978" ),
   MpIeee( "0.0001162365390497009" ),
  -MpIeee( "0.0000035904327241606" ),
   MpIeee( "0.0000000802342123706" ),
  -MpIeee( "0.0000000013562997693" ),
   MpIeee( "0.0000000000179440722" ),
  -MpIeee( "0.0000000000001908387" ),
   MpIeee( "0.0000000000000016670" ),
  -MpIeee( "0.0000000000000000122" )
};

static cheb_series si_cs = {
  si_data,
  11,
  -1, 1,
  9
};

/*
 series for ci   on the interval  0.00000e+00 to  1.60000e+01
                                        with weighted error   1.94e-18
                                         log weighted error  17.71
                               significant figures required  17.74
                                    decimal places required  18.27
*/
static MpIeee ci_data[13] =  {
   -MpIeee( "0.34004281856055363156" ),
   -MpIeee( "1.03302166401177456807" ),
    MpIeee( "0.19388222659917082877" ),
   -MpIeee( "0.01918260436019865894" ),
    MpIeee( "0.00110789252584784967" ),
   -MpIeee( "0.00004157234558247209" ),
    MpIeee( "0.00000109278524300229" ),
   -MpIeee( "0.00000002123285954183" ),
    MpIeee( "0.00000000031733482164" ),
   -MpIeee( "0.00000000000376141548" ),
    MpIeee( "0.00000000000003622653" ),
   -MpIeee( "0.00000000000000028912" ),
    MpIeee( "0.00000000000000000194" )
};
static cheb_series ci_cs = {
  ci_data,
  12,
  -1, 1,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_Si_e(const MpIeee x, gsl_sf_result * result)
{
  MpIeee ax=  fabs(x);
  
  /* CHECK_POINTER(result) */

  if(ax < GSL_SQRT_DBL_EPSILON) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax <= MpIeee( "4.0" )) {
    gsl_sf_result result_c;
    cheb_eval_e(&si_cs, (x*x-8.0)*0.125, &result_c);
    result->val  =  x * (0.75 + result_c.val);
    result->err  = ax * result_c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* Note there is no loss of precision
     * here bcause of the leading constant.
     */
    gsl_sf_result f;
    gsl_sf_result g;
    fg_asymp(ax, &f, &g);
    result->val  = 0.5 * M_PI - f.val*cos(ax) - g.val*sin(ax);
    result->err  = f.err + g.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    if(x < 0.0) result->val = -result->val;
    return GSL_SUCCESS;
  }
}


int  gsl_sf_Ci_e(const MpIeee x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x <= 4.0) {
    const MpIeee lx=  log(x);
    const MpIeee y=  (x*x-8.0)*0.125;
    gsl_sf_result result_c;
    cheb_eval_e(&ci_cs, y, &result_c);
    result->val  = lx - 0.5 + result_c.val;
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(lx) + 0.5) + result_c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result sin_result;
    gsl_sf_result cos_result;
    int  stat_sin=  gsl_sf_sin_e(x, &sin_result);
    int  stat_cos=  gsl_sf_cos_e(x, &cos_result);
    gsl_sf_result f;
    gsl_sf_result g;
    fg_asymp(x, &f, &g);
    result->val  = f.val*sin_result.val - g.val*cos_result.val;
    result->err  = fabs(f.err*sin_result.val);
    result->err += fabs(g.err*cos_result.val);
    result->err += fabs(f.val*sin_result.err);
    result->err += fabs(g.val*cos_result.err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_sin, stat_cos);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_Si(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_Si_e(x, &result));
}

MpIeee gsl_sf_Ci(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_Ci_e(x, &result));
}
