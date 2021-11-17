/****************************************************************************
 * RealPaver v. 0.4                                                         *
 *--------------------------------------------------------------------------*
 * Author: Laurent Granvilliers                                             *
 * Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes  *
 * Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique  *
 *--------------------------------------------------------------------------*
 * RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated       *
 * COPYRIGHT file for more details.                                         *
 *--------------------------------------------------------------------------*
 * default_fpu.h                                                            *
 ****************************************************************************/

#ifndef __default_fpu_h
#define __default_fpu_h

#include <math.h>
#include "config.h"

extern double ceil();
extern double floor();

extern double log2();
extern double log();
extern double exp();

extern double sqrt();

extern double cos();
extern double sin();
extern double tan();
extern double cosh();
extern double sinh();
extern double tanh();
extern double acos();
extern double asin();
extern double atan();
extern double acosh();
extern double asinh();
extern double atanh();

extern double fabs();
extern double nextafter();

extern double pow(double,double);


/****************************************************************************
 *                                  SPARC                                   *
 ****************************************************************************/
#if SYSTEM_SPARC

#  include <floatingpoint.h>
#  include <sys/ieeefp.h>

#  define RoundDownward fpsetround(3)
#  define RoundUpward   fpsetround(2)
#  define RoundNearest  fpsetround(0)

   static const double IBBasicInfinity = 1.0/0.0;
#  define IBBasicMaxDouble 1.797693134862315708e+308

/****************************************************************************
 *                             PC i386 & linux                              *
 ****************************************************************************/
#elif SYSTEM_LINUX_IX86

#  include <fenv.h>
#  include <values.h>

#  define RoundDownward    fesetround(FE_DOWNWARD)
#  define RoundUpward      fesetround(FE_UPWARD)
#  define RoundNearest     fesetround(FE_TONEAREST)

#  define IBBasicInfinity  HUGE_VAL
#  define IBBasicMaxDouble MAXDOUBLE

/****************************************************************************
 *                                MIPS SGI                                  *
 ****************************************************************************/
#elif SYSTEM_SGI

#  include <float.h>
#  include <ieeefp.h>

#  define RoundDownward fpsetround(3)
#  define RoundUpward   fpsetround(2)
#  define RoundNearest  fpsetround(0)

   static const double IBBasicInfinity = 1.0/0.0;
#  define IBBasicMaxDouble DBL_MAX

#endif
/****************************************************************************/

#define IBBasicConstPi       3.14159265358979323846
#define IBBasicConstHalfPi   1.57079632679489661923
#define IBBasicConstE        2.7182818284590452354
#define IBBasicConstLn2      0.69314718055994530942

#define IBBasicMinDouble     -IBBasicMaxDouble
#define IBBasicNegInfinity   -IBBasicInfinity
#define IBBasicPosInfinity   +IBBasicInfinity

#define IBBasicIsOdd(n)      (((n)%2)==1)
#define IBBasicIsEven(n)     (((n)%2)==0)

#define IBBasicMin(x,y)      (((x)<(y))? x : y)
#define IBBasicMax(x,y)      (((x)<(y))? y : x)


static inline double IBBasicBisectionPoint(double x, double y, int h, int n) {
/***************************************************************************
* Subdivision of interval [x,y] :
*      x   x+(1/h)(y-x)   x+(2/h)(y-x)   ...   x+(h/h)(y-x)==y  
*
*   Returns x + (n/h)(y-x)
*
*   Output guarantee: if [x,y] contains at least three floating-point numbers,
*   then the result is included in ]x,y[
*
*   Input Conditions: x <= y and not (x==y==-oo or x==y==+oo), 1 <= n <= h
*/
  double a = IBBasicMax(x,IBBasicMinDouble),
         b = IBBasicMin(IBBasicMaxDouble,y),
         c;

  if (a==b) {    /* [r,r] [-oo,MinReal], [MaxReal,+oo] */
    return a;
  }
  else if (nextafter(a,IBBasicPosInfinity)==b) {  /* [r,r+], [-oo,succ MinReal], [pred MaxReal,+oo] */
    if (x==IBBasicNegInfinity) {
      return IBBasicMinDouble;
    }
    else if (y==IBBasicPosInfinity) {
      return IBBasicMaxDouble;
    }
    else {
      return a;
    }
  }
  else {    /* [r,s], s> succ r */
    c = a + (((double)n)/((double)h))*(b-a);
    if ((c>a) && (c<b)) {
      return c;
    }
    else {
      return nextafter(a,IBBasicPosInfinity);
    }
  }
}


static double IBBasicPowReal(double x, int n)
/***************************************************************************
*  Returns x^n
*/
{
  double y = x, z = 1.0;
  for( ;; )
  {
    if( IBBasicIsOdd(n) )
    {
      z *= y;
      n >>= 1;
      if( n==0 ) return( z );
    }
    else n >>= 1;
    y *= y;
  }
}


static double IBBasicNthRoot(double x, int n, int round, double epsilon)
/***************************************************************************
*  Returns the n-th root of x
*  round=1 (downward) ; round=2 (upward)
*/
{
  double y;

  if( x==1.0 )
  {
    return( 1.0 );
  }
  else if( x==0.0 )
  {
    return( 0.0 );
  }
  else if( n==2 )
  {
    if (round==1) {
      RoundDownward;
    }
    else {
      RoundUpward;
    }
    return( sqrt(x) );
  }
  else
  {
    RoundNearest;
    if( x<0.0 )
    {
      y = -(exp((1.0/n)*log(-(x))));
    }
    else
    {
      y = exp((1.0/n)*log(x));
    }
    if (round==1) {
      RoundDownward;
      return( y-epsilon );
    }
    else {
      RoundUpward;
      return( y+epsilon );
    }
  }
}


/* specific bisection, next/prev-double, rounding functions */

#define IBBasicCenter(x,y)    IBBasicBisectionPoint(x,y,2,1)  /* Center of [x,y] */

/* let [x,y] = [x,a] union [a,b] union [b,y], in three "equal" parts */
#define IBBasicThird(x,y)     IBBasicBisectionPoint(x,y,3,1)  /* Returns a */
#define IBBasicTwoThirds(x,y) IBBasicBisectionPoint(x,y,3,2)  /* Returns b */


#define IBBasicNextDouble(x) nextafter(x,IBBasicPosInfinity)  /* Double after x */
#define IBBasicPrevDouble(x) nextafter(x,IBBasicNegInfinity)  /* Double before x */

#define IBBasicRoundDown()   RoundDownward      /* Rounding towards -oo */
#define IBBasicRoundUp()     RoundUpward        /* Rounding towards +oo */
#define IBBasicRoundNear()   RoundNearest       /* Rounding to nearest */


#endif
