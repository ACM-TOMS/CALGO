#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/debye.c
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
#include <gsl/gsl_sf_debye.h>

#include "error.h"
#include "check.h"

#include "chebyshev.h"
#include "cheb_eval.c"

static MpIeee adeb1_data[17] =  {
   MpIeee( "2.4006597190381410194" ),
   MpIeee( "0.1937213042189360089" ),
  -MpIeee( "0.62329124554895770e-02" ),
   MpIeee( "0.3511174770206480e-03" ),
  -MpIeee( "0.228222466701231e-04" ),
   MpIeee( "0.15805467875030e-05" ),
  -MpIeee( "0.1135378197072e-06" ),
   MpIeee( "0.83583361188e-08" ),
  -MpIeee( "0.6264424787e-09" ),
   MpIeee( "0.476033489e-10" ),
  -MpIeee( "0.36574154e-11" ),
   MpIeee( "0.2835431e-12" ),
  -MpIeee( "0.221473e-13" ),
   MpIeee( "0.17409e-14" ),
  -MpIeee( "0.1376e-15" ),
   MpIeee( "0.109e-16" ),
  -MpIeee( "0.9e-18" )
};
static cheb_series adeb1_cs = {
  adeb1_data,
  16,
  -1.0, 1.0,
  9
};

static MpIeee adeb2_data[18] =  {
   MpIeee( "2.5943810232570770282" ),
   MpIeee( "0.2863357204530719834" ),
  -MpIeee( "0.102062656158046713e-01" ),
   MpIeee( "0.6049109775346844e-03" ),
  -MpIeee( "0.405257658950210e-04" ),
   MpIeee( "0.28633826328811e-05" ),
  -MpIeee( "0.2086394303065e-06" ),
   MpIeee( "0.155237875826e-07" ),
  -MpIeee( "0.11731280087e-08" ),
   MpIeee( "0.897358589e-10" ),
  -MpIeee( "0.69317614e-11" ),
   MpIeee( "0.5398057e-12" ),
  -MpIeee( "0.423241e-13" ),
   MpIeee( "0.33378e-14" ),
  -MpIeee( "0.2645e-15" ),
   MpIeee( "0.211e-16" ),
  -MpIeee( "0.17e-17" ),
   MpIeee( "0.1e-18" )
};
static cheb_series adeb2_cs = {
  adeb2_data,
  17,
  -1.0, 1.0,
  10
};

static MpIeee adeb3_data[17] =  {
   MpIeee( "2.707737068327440945" ),
   MpIeee( "0.340068135211091751" ),
  -MpIeee( "0.12945150184440869e-01" ),
   MpIeee( "0.7963755380173816e-03" ),
  -MpIeee( "0.546360009590824e-04" ),
   MpIeee( "0.39243019598805e-05" ),
  -MpIeee( "0.2894032823539e-06" ),
   MpIeee( "0.217317613962e-07" ),
  -MpIeee( "0.16542099950e-08" ),
   MpIeee( "0.1272796189e-09" ),
  -MpIeee( "0.987963460e-11" ),
   MpIeee( "0.7725074e-12" ),
  -MpIeee( "0.607797e-13" ),
   MpIeee( "0.48076e-14" ),
  -MpIeee( "0.3820e-15" ),
   MpIeee( "0.305e-16" ),
  -MpIeee( "0.24e-17" )
};
static cheb_series adeb3_cs = {
  adeb3_data,
  16,
  -1.0, 1.0,
  10
};

static MpIeee adeb4_data[17] =  {
   MpIeee( "2.781869415020523460" ),
   MpIeee( "0.374976783526892863" ),
  -MpIeee( "0.14940907399031583e-01" ),
   MpIeee( "0.945679811437042e-03" ),
  -MpIeee( "0.66132916138933e-04" ),
   MpIeee( "0.4815632982144e-05" ),
  -MpIeee( "0.3588083958759e-06" ),
   MpIeee( "0.271601187416e-07" ),
  -MpIeee( "0.20807099122e-08" ),
   MpIeee( "0.1609383869e-09" ),
  -MpIeee( "0.125470979e-10" ),
   MpIeee( "0.9847265e-12" ),
  -MpIeee( "0.777237e-13" ),
   MpIeee( "0.61648e-14" ),
  -MpIeee( "0.4911e-15" ),
   MpIeee( "0.393e-16" ),
  -MpIeee( "0.32e-17" )
};
static cheb_series adeb4_cs = {
  adeb4_data,
  16,
  -1.0, 1.0,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_sf_debye_1_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee val_infinity=  1.64493406684822644;
  const MpIeee xcut=  -GSL_LOG_DBL_MIN;

  /* CHECK_POINTER(result) */

  if(x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - 0.25*x + x*x/36.0;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const MpIeee t=  x*x/8.0 - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&adeb1_cs, t, &c);
    result->val = c.val - 0.25 * x;
    result->err = c.err + 0.25 * x * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const MpIeee ex=  exp(-x);
    MpIeee sum=  MpIeee( "0.0" );
    MpIeee xk=  nexp * x;
    MpIeee rk=  nexp;
    int  i;
    for(i=nexp; i>=1; i--) {
      sum *= ex;
      sum += (MpIeee( "1.0" ) + MpIeee( "1.0" )/xk)/rk;
      rk -= MpIeee( "1.0" );
      xk -= x;
    }
    result->val = val_infinity/x - sum*ex;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    result->val = (val_infinity - exp(-x)*(x+1.0)) / x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = val_infinity/x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

    
int  gsl_sf_debye_2_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee val_infinity=  4.80822761263837714;
  const MpIeee xcut=  -GSL_LOG_DBL_MIN;

  /* CHECK_POINTER(result) */

  if(x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*M_SQRT2*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - x/3.0 + x*x/24.0;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const MpIeee t=  x*x/8.0 - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&adeb2_cs, t, &c);
    result->val = c.val - x/3.0;
    result->err = c.err + GSL_DBL_EPSILON * x/3.0;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const MpIeee ex=  exp(-x);
    MpIeee xk=  nexp * x;
    MpIeee rk=  nexp;
    MpIeee sum=  MpIeee( "0.0" );
    int  i;
    for(i=nexp; i>=1; i--) {
      sum *= ex;
      sum += (MpIeee( "1.0" ) + MpIeee( "2.0" )/xk + MpIeee( "2.0" )/(xk*xk)) / rk;
      rk -= MpIeee( "1.0" );
      xk -= x;
    }
    result->val = val_infinity/(x*x) - 2.0 * sum * ex;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    const MpIeee x2=  x*x;
    const MpIeee sum=  2.0 + 2.0*x + x2;
    result->val = (val_infinity - 2.0 * sum * exp(-x)) / x2;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = (val_infinity/x)/x;
    result->err = GSL_DBL_EPSILON * result->val;
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_debye_3_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee val_infinity=  19.4818182068004875;
  const MpIeee xcut=  -GSL_LOG_DBL_MIN;

  /* CHECK_POINTER(result) */

  if(x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*M_SQRT2*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - 3.0*x/8.0 + x*x/20.0;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const MpIeee t=  x*x/8.0 - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&adeb3_cs, t, &c);
    result->val = c.val - 0.375*x;
    result->err = c.err + GSL_DBL_EPSILON * 0.375*x;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const MpIeee ex=  exp(-x);
    MpIeee xk=  nexp * x;
    MpIeee rk=  nexp;
    MpIeee sum=  MpIeee( "0.0" );
    int  i;
    for(i=nexp; i>=1; i--) {
      MpIeee xk_inv=  MpIeee( "1.0" )/xk;
      sum *= ex;
      sum += (((MpIeee( "6.0" )*xk_inv + MpIeee( "6.0" ))*xk_inv + MpIeee( "3.0" ))*xk_inv + MpIeee( "1.0" )) / rk;
      rk -= MpIeee( "1.0" );
      xk -= x;
    }
    result->val = val_infinity/(x*x*x) - 3.0 * sum * ex;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    const MpIeee x3=  x*x*x;
    const MpIeee sum=  6.0 + 6.0*x + 3.0*x*x + x3;
    result->val = (val_infinity - 3.0 * sum * exp(-x)) / x3;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    result->val = ((val_infinity/x)/x)/x;
    result->err = GSL_DBL_EPSILON * result->val;
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


int  gsl_sf_debye_4_e(const MpIeee x, gsl_sf_result * result)
{
  const MpIeee val_infinity=  99.5450644937635129;
  const MpIeee xcut=  -GSL_LOG_DBL_MIN;

  /* CHECK_POINTER(result) */

  if(x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*M_SQRT2*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - 2.0*x/5.0 + x*x/18.0;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const MpIeee t=  x*x/8.0 - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&adeb4_cs, t, &c);
    result->val = c.val - 2.0*x/5.0;
    result->err = c.err + GSL_DBL_EPSILON * 2.0*x/5.0;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const MpIeee ex=  exp(-x);
    MpIeee xk=  nexp * x;
    MpIeee rk=  nexp;
    MpIeee sum=  MpIeee( "0.0" );
    int  i;
    for(i=nexp; i>=1; i--) {
      MpIeee xk_inv=  MpIeee( "1.0" )/xk;
      sum *= ex;
      sum += ((((MpIeee( "24.0" )*xk_inv + MpIeee( "24.0" ))*xk_inv + MpIeee( "12.0" ))*xk_inv + MpIeee( "4.0" ))*xk_inv + MpIeee( "1.0" )) / rk;
      rk -= MpIeee( "1.0" );
      xk -= x;
    }
    result->val = val_infinity/(x*x*x*x) - 4.0 * sum * ex;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    const MpIeee x2=  x*x;
    const MpIeee x4=  x2*x2;
    const MpIeee sum=  24.0 + 24.0*x + 12.0*x2 + 4.0*x2*x + x4;
    result->val = (val_infinity - 4.0 * sum * exp(-x)) / x4;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    result->val = (((val_infinity/x)/x)/x)/x;
    result->err = GSL_DBL_EPSILON * result->val;
    CHECK_UNDERFLOW(result);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_debye_1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_debye_1_e(x, &result));
}

MpIeee gsl_sf_debye_2(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_debye_2_e(x, &result));
}

MpIeee gsl_sf_debye_3(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_debye_3_e(x, &result));
}

MpIeee gsl_sf_debye_4(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_debye_4_e(x, &result));
}
