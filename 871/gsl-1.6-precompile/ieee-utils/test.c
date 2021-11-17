#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ieee-utils/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>

#if HAVE_IRIX_IEEE_INTERFACE
/* don't test denormals on IRIX */
#else
#if HAVE_IEEE_DENORMALS
#define TEST_DENORMAL 1
#endif
#endif

#ifndef FLT_MIN
#define FLT_MIN 1.17549435e-38f
#endif

#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38f
#endif

#ifndef DBL_MIN
#define DBL_MIN 2.2250738585072014e-308
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e+308
#endif

int
main (void)
{
  MpIeee zerof=  MpIeee( "0.0f" );MpIeee  minus_onef=  -MpIeee( "1.0f" ) ;
  MpIeee zero=  MpIeee( "0.0" );MpIeee  minus_one=  -MpIeee( "1.0" ) ;

  /* Check for +ZERO (float) */

  {
    MpIeee f=  MpIeee( "0.0f" );
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "float x = 0, sign is +");
    gsl_test_int (r.exponent, -127, "float x = 0, exponent is -127");
    gsl_test_str (r.mantissa, mantissa, "float x = 0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "float x = 0, type is ZERO");
  }

  /* Check for -ZERO (float) */

  {
    MpIeee f=  minus_onef;
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;
    
    while (f < MpIeee( "0" )) {
      f *= MpIeee( "0.1f" );
    }

    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 1, "float x = -1*0, sign is -");
    gsl_test_int (r.exponent, -127, "float x = -1*0, exponent is -127");
    gsl_test_str (r.mantissa, mantissa, "float x = -1*0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "float x = -1*0, type is ZERO");
  }

  /* Check for a positive NORMAL number (e.g. 2.1) (float) */

  {
    MpIeee f=  MpIeee( "2.1f" );
    const char mantissa[] = "00001100110011001100110";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "float x = 2.1, sign is +");
    gsl_test_int (r.exponent, 1, "float x = 2.1, exponent is 1");
    gsl_test_str (r.mantissa, mantissa, "float x = 2.1, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "float x = 2.1, type is NORMAL");
  }


  /* Check for a negative NORMAL number (e.g. -1.3304...) (float) */

  {
    MpIeee f=  -MpIeee( "1.3303577090924210f" ) ;
    const char mantissa[] = "01010100100100100101001";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 1, "float x = -1.3304..., sign is -");
    gsl_test_int (r.exponent, 0, "float x = -1.3304..., exponent is 0");
    gsl_test_str (r.mantissa, mantissa, "float x = -1.3304..., mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "float x = -1.3304..., type is NORMAL");
  }

  /* Check for a large positive NORMAL number (e.g. 3.37e31) (float) */

  {
    MpIeee f=  MpIeee( "3.37e31f" );
    const char mantissa[] = "10101001010110101001001";
    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "float x = 3.37e31, sign is +");
    gsl_test_int (r.exponent, 104, "float x = 3.37e31, exponent is 104");
    gsl_test_str (r.mantissa, mantissa, "float x = 3.37e31, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "float x = 3.37e31, type is NORMAL");
  }

  /* Check for a small positive NORMAL number (e.g. 3.37e-31) (float) */

  {
    MpIeee f=  MpIeee( "3.37e-31f" );
    const char mantissa[] = "10110101011100110111011";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "float x = 3.37e-31, sign is +");
    gsl_test_int (r.exponent, -102, "float x = 3.37e-31, exponent is -102");
    gsl_test_str (r.mantissa, mantissa, "float x = 3.37e-31, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "float x = 3.37e-31, type is NORMAL");
  }

  /* Check for FLT_MIN (smallest possible number that is not denormal) */

  {
    MpIeee f=  FLT_MIN;  
    const char mantissa[] = "00000000000000000000000";
    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "float x = FLT_MIN, sign is +");
    gsl_test_int (r.exponent, -126, "float x = FLT_MIN, exponent is -126");
    gsl_test_str (r.mantissa, mantissa, "float x = FLT_MIN, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "float x = FLT_MIN, type is NORMAL");
  }

  /* Check for FLT_MAX (largest possible number that is not Inf) */

  {
    MpIeee f=  FLT_MAX;
    const char mantissa[] = "11111111111111111111111";

    gsl_ieee_float_rep r;
    gsl_ieee_float_to_rep (&f, &r);

    gsl_test_int (r.sign, 0, "float x = FLT_MAX, sign is +");
    gsl_test_int (r.exponent, 127, "float x = FLT_MAX, exponent is 127");
    gsl_test_str (r.mantissa, mantissa, "float x = FLT_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "float x = FLT_MAX, type is NORMAL");
  }


  /* Check for DENORMAL numbers (e.g. FLT_MIN/2^n) */

#ifdef TEST_DENORMAL
  {
    MpIeee f=  FLT_MIN;  
    char mantissa[] = "10000000000000000000000";

    int  i;
    gsl_ieee_float_rep r;

    for (i = 0; i < 23; i++)
      {
        MpIeee x=  f / (MpIeee)pow (MpIeee( "2.0" ), MpIeee( "1" ) + (MpIeee) i);
        mantissa[i] = '1';
        gsl_ieee_float_to_rep (&x, &r);

        gsl_test_int (r.sign, 0, "float x = FLT_MIN/2^%d, sign is +", i + 1);
        gsl_test_int (r.exponent, -127,
                      "float x = FLT_MIN/2^%d, exponent is -127", i + 1);
        gsl_test_str (r.mantissa, mantissa,
                      "float x = FLT_MIN/2^%d, mantissa", i + 1);
        gsl_test_int (r.type, GSL_IEEE_TYPE_DENORMAL,
                      "float x = FLT_MIN/2^%d, type is DENORMAL", i + 1);
        mantissa[i] = '0';
      }
  }
#endif

  /* Check for positive INFINITY (e.g. 2*FLT_MAX) */

  {
    MpIeee f=  FLT_MAX;  
    const char mantissa[] = "00000000000000000000000";

    gsl_ieee_float_rep r;

    MpIeee x;
    x = MpIeee( "2" ) * f;
    gsl_ieee_float_to_rep (&x, &r);

    gsl_test_int (r.sign, 0, "float x = 2*FLT_MAX, sign is +");
    gsl_test_int (r.exponent, 128, "float x = 2*FLT_MAX, exponent is 128");
    gsl_test_str (r.mantissa, mantissa, "float x = 2*FLT_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF, "float x = -2*FLT_MAX, type is INF");
  }

  /* Check for negative INFINITY (e.g. -2*FLT_MAX) */

  {
    MpIeee f=  FLT_MAX;  
    const char mantissa[] = "00000000000000000000000";

    gsl_ieee_float_rep r;

    MpIeee x;
    x = -MpIeee( "2" ) * f;
    gsl_ieee_float_to_rep (&x, &r);

    gsl_test_int (r.sign, 1, "float x = -2*FLT_MAX, sign is -");
    gsl_test_int (r.exponent, 128, "float x = -2*FLT_MAX, exponent is 128");
    gsl_test_str (r.mantissa, mantissa, "float x = -2*FLT_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF, "float x = -2*FLT_MAX, type is INF");
  }

  /* Check for NAN (e.g. Inf - Inf) (float) */

  {
    gsl_ieee_float_rep r;
    MpIeee x=  MpIeee( "1.0f" );MpIeee  y=  MpIeee( "2.0f" );MpIeee  z=  zerof;

    x = x / z;
    y = y / z;
    z = y - x;

    gsl_ieee_float_to_rep (&z, &r);

    /* We don't check the sign and we don't check the mantissa because
       they could be anything for a NaN */

    gsl_test_int (r.exponent, 128, "float x = NaN, exponent is 128");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NAN, "float x = NaN, type is NAN");
  }


  /* Check for +ZERO */

  {
    MpIeee d=  MpIeee( "0.0" );
    const char mantissa[]
      = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "double x = 0, sign is +");
    gsl_test_int (r.exponent, -1023, "double x = 0, exponent is -1023");
    gsl_test_str (r.mantissa, mantissa, "double x = 0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "double x = 0, type is ZERO");
  }

  /* Check for -ZERO */

  {
    MpIeee d=   minus_one;
    const char mantissa[]
      = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;

    while (d < MpIeee( "0" )) {
      d *= MpIeee( "0.1" );
    }

    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 1, "double x = -1*0, sign is -");
    gsl_test_int (r.exponent, -1023, "double x = -1*0, exponent is -1023");
    gsl_test_str (r.mantissa, mantissa, "double x = -1*0, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_ZERO, "double x = -1*0, type is ZERO");
  }

  /* Check for a positive NORMAL number (e.g. 2.1) */

  {
    MpIeee d=  MpIeee( "2.1" );
    const char mantissa[]
      = "0000110011001100110011001100110011001100110011001101";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "double x = 2.1, sign is +");
    gsl_test_int (r.exponent, 1, "double x = 2.1, exponent is 1");
    gsl_test_str (r.mantissa, mantissa, "double x = 2.1, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL, "double x = 2.1, type is NORMAL");
  }


  /* Check for a negative NORMAL number (e.g. -1.3304...) */

  {
    MpIeee d=  -MpIeee( "1.3303577090924210146738460025517269968986511230468750" );
    const char mantissa[]
      = "0101010010010010010100101010010010001000100011101110";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 1, "double x = -1.3304..., sign is -");
    gsl_test_int (r.exponent, 0, "double x = -1.3304..., exponent is 0");
    gsl_test_str (r.mantissa, mantissa, "double x = -1.3304..., mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "double x = -1.3304..., type is NORMAL");
  }

  /* Check for a large positive NORMAL number (e.g. 3.37e297) */

  {
    MpIeee d=  MpIeee( "3.37e297" );
    const char mantissa[]
      = "0100100111001001100101111001100000100110011101000100";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "double x = 3.37e297, sign is +");
    gsl_test_int (r.exponent, 988, "double x = 3.37e297, exponent is 998");
    gsl_test_str (r.mantissa, mantissa, "double x = 3.37e297, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "double x = 3.37e297, type is NORMAL");
  }

  /* Check for a small positive NORMAL number (e.g. 3.37e-297) */

  {
    MpIeee d=  MpIeee( "3.37e-297" );
    const char mantissa[]
    = "0001101000011011101011100001110010100001001100110111";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "double x = 3.37e-297, sign is +");
    gsl_test_int (r.exponent, -985, "double x = 3.37e-297, exponent is -985");
    gsl_test_str (r.mantissa, mantissa, "double x = 3.37e-297, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "double x = 3.37e-297, type is NORMAL");
  }

  /* Check for DBL_MIN (smallest possible number that is not denormal) */

  {
    MpIeee d=  DBL_MIN;
    const char mantissa[]
      = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "double x = DBL_MIN, sign is +");
    gsl_test_int (r.exponent, -1022, "double x = DBL_MIN, exponent is -1022");
    gsl_test_str (r.mantissa, mantissa, "double x = DBL_MIN, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "double x = DBL_MIN, type is NORMAL");
  }

  /* Check for DBL_MAX (largest possible number that is not Inf) */

  {
    MpIeee d=  DBL_MAX;
    const char mantissa[]
    = "1111111111111111111111111111111111111111111111111111";
    gsl_ieee_double_rep r;
    gsl_ieee_double_to_rep (&d, &r);

    gsl_test_int (r.sign, 0, "double x = DBL_MAX, sign is +");
    gsl_test_int (r.exponent, 1023, "double x = DBL_MAX, exponent is 1023");
    gsl_test_str (r.mantissa, mantissa, "double x = DBL_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NORMAL,
                  "double x = DBL_MAX, type is NORMAL");
  }

  /* Check for DENORMAL numbers (e.g. DBL_MIN/2^n) */

#ifdef TEST_DENORMAL
  {
    MpIeee d=  DBL_MIN;
    char mantissa[]
      = "1000000000000000000000000000000000000000000000000000";
    int  i;
    gsl_ieee_double_rep r;

    for (i = 0; i < 52; i++)
      {
        MpIeee x=  d / pow (MpIeee( "2.0" ), MpIeee( "1" ) + (MpIeee) i);
        mantissa[i] = '1';
        gsl_ieee_double_to_rep (&x, &r);

        gsl_test_int (r.sign, 0, "double x = DBL_MIN/2^%d, sign is +", i + 1);
        gsl_test_int (r.exponent, -1023,
                      "double x = DBL_MIN/2^%d, exponent", i + 1);
        gsl_test_str (r.mantissa, mantissa,
                      "double x = DBL_MIN/2^%d, mantissa", i + 1);
        gsl_test_int (r.type, GSL_IEEE_TYPE_DENORMAL,
                      "double x = DBL_MIN/2^%d, type is DENORMAL", i + 1);
        mantissa[i] = '0';
      }
  }
#endif

  /* Check for positive INFINITY (e.g. 2*DBL_MAX) */

  {
    MpIeee d=  DBL_MAX;
    const char mantissa[]
      = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;

    MpIeee x;
    x = MpIeee( "2.0" ) * d;
    gsl_ieee_double_to_rep (&x, &r);

    gsl_test_int (r.sign, 0, "double x = 2*DBL_MAX, sign is +");
    gsl_test_int (r.exponent, 1024, "double x = 2*DBL_MAX, exponent is 1024");
    gsl_test_str (r.mantissa, mantissa, "double x = 2*DBL_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF, "double x = 2*DBL_MAX, type is INF");
  }

  /* Check for negative INFINITY (e.g. -2*DBL_MAX) */

  {
    MpIeee d=  DBL_MAX;
    const char mantissa[]
      = "0000000000000000000000000000000000000000000000000000";
    gsl_ieee_double_rep r;

    MpIeee x;
    x = -MpIeee( "2.0" ) * d;
    gsl_ieee_double_to_rep (&x, &r);

    gsl_test_int (r.sign, 1, "double x = -2*DBL_MAX, sign is -");
    gsl_test_int (r.exponent, 1024, "double x = -2*DBL_MAX, exponent is 1024");
    gsl_test_str (r.mantissa, mantissa, "double x = -2*DBL_MAX, mantissa");
    gsl_test_int (r.type, GSL_IEEE_TYPE_INF,"double x = -2*DBL_MAX, type is INF");
  }

  /* Check for NAN (e.g. Inf - Inf) */

  {
    gsl_ieee_double_rep r;
    MpIeee x=  MpIeee( "1.0" );MpIeee  y=  MpIeee( "2.0" );MpIeee  z=  zero;

    x = x / z;
    y = y / z;
    z = y - x;

    gsl_ieee_double_to_rep (&z, &r);

    /* We don't check the sign and we don't check the mantissa because
       they could be anything for a NaN */

    gsl_test_int (r.exponent, 1024, "double x = NaN, exponent is 1024");
    gsl_test_int (r.type, GSL_IEEE_TYPE_NAN, "double x = NaN, type is NAN");
  }

  exit (gsl_test_summary ());
}
