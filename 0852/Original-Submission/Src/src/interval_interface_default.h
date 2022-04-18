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
 * interval_interface_default.h                                             *
 ****************************************************************************/

#ifndef __interval_interface_default_h
#define __interval_interface_default_h

#  include <stdio.h>
#  include "default_fpu.h"
#  include "default_interval.h"

/* Interval type */
#  define IBInterval                 IBBasicBounds

/* Left bound of interval i */
#  define IBLeftI(i)                 IBBasicLeftI(i)

/* Right bound of interval i */
#  define IBRightI(i)                IBBasicRightI(i)

/* +oo */
#  define IBPosInfinity              IBBasicPosInfinity

/* -oo */
#  define IBNegInfinity              IBBasicNegInfinity

/* Returns x+, the floating-point number following x */
#define IBNextReal(x)                IBBasicNextDouble(x)

/* Returns x-, the floating-point number followed by x */
#define IBPrevReal(x)                IBBasicPrevDouble(x)

/* Rounding modes */
#define IBRoundDown()                IBBasicRoundDown() 
#define IBRoundUp()                  IBBasicRoundUp()
#define IBRoundNear()                IBBasicRoundNear()

/* Returns true if interval i is empty */
#  define IBIsEmptyI(i)              _IBBasicEmptyI(i)

/* Returns true if intervals i1 and i2 are equal */
#  define IBIsEqualII(i1,i2)         _IBBasicIeqI(i1,i2)

/* Returns true if interval i1 is different from i2 */
#  define IBIsDifferentII(i1,i2)     _IBBasicIdiffI(i1,i2)

/* Returns true if double x belongs to interval i */
#  define IBIsDoubleInI(i,x)         _IBBasicDoubleInI(i,x)

/* Returns true if interval i is reduced to one real number */
#  define IBIsIntervalPoint(i)       _IBBasicIsDoubleI(i)

/* Returns true if interval i is reduced to one natural number */
#  define IBIsIntervalIntPoint(i)     _IBBasicIsIntegerI(i)

/* Returns true if interval i is reduced to 0 */
#  define IBIsReducedToZeroI(i)      _IBBasicZeroI(i)

/* Returns true if interval i1 is included in interval i2 */
#  define IBIsIncludedII(i1,i2)      _IBBasicIncludedII(i1,i2)

/* Returns true both intervals i1 and i2 are disjoint */
#  define IBIsDisjointII(i1,i2)      _IBBasicDisjointII(i1,i2)

/* Returns true if at least one bound of i is infinite */
#  define IBIsInfiniteI(i)           _IBBasicInfinite(i)

/* Returns the width of interval i rounded towards +oo */
#  define IBWidthOfI(i)              _IBBasicWidthI(i)

/* Returns the distance between i1 and i2 rounded towards +oo */
#  define IBDistanceBetweenII(i1,i2) _IBBasicDistanceII(i1,i2)

/* Returns the midpoint of i */
#  define IBMidpointOfI(i)           _IBBasicMidI(i)

/* Bisection of i=[a,b] in three equal parts :
   |-|-|-|
   a x   b
   Returns x */
#  define IBThirdOfI(i)              _IBBasicThirdI(i)

/* Bisection of i=[a,b] in three equal parts :
   |-|-|-|
   a   y b
   Returns y */
#  define IBTwoThirdsOfI(i)          _IBBasicTwoThirdsI(i)

/* Returns true if i contains at most two floating point numbers :
   i is either [a,a] or [a,a+] */
#  define IBIsCanonicalI(i)          _IBBasicCanonicalI(i)

/* i := empty interval */
#  define IBSetToEmptyI(i)           _IBBasicSetEmptyI(i)

/* i := [x1,x2] */
#  define IBSetBoundsOfI(i,x1,x2)    _IBBasicSetI(i,x1,x2)

/* i := source */
#  define IBCopyII(i,source)         _IBBasicCopyI(i,source)

/* Returns a pointer to a new interval in memory */
static inline IBInterval* IBCreateNewI() {
  return IBBasicNewI();
}

/* Returns a pointer to a new interval in memory which value is [-oo,+oo] */
static inline IBInterval* IBCreateNewRealDomainI() {
  return IBBasicNewLargestI();
}

/* Returns a pointer to a new interval in memory which value is a copy of *i */
static inline IBInterval* IBCreateAndCopyNewI(IBInterval* i) {
  return IBBasicNewCopyI(i);
}

/* Returns a pointer to a new interval in memory which value is [x1,x2] */
static inline IBInterval* IBCreateAndSetNewI(double x1, double x2) {
  return IBBasicSetNewI(x1,x2);
}

/* *i := [-oo,+oo] */
static inline void IBSetToRealDomain(IBInterval* i) {
  IBBasicToLargestI(i);
}

/* *i := [ceil(inf i), floor(sup i)] */
static inline void IBSetToIntegerDomain(IBInterval* i) {
  IBBasicToIntegerI(i);
}

/* j := intersection of i1 and i2 */
static inline void IBIntersectionII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicInterII(j,i1,i2);
}

/* Print the value of *i in output file out :
   - digits is the number of digits of bounds to be printed
   - The value of mode is:
     * IBPrintIntervalBounds:   *i is written '[a,b]'
     * IBPrintIntervalMidError: *i is written 'midpoint + [-e,e]'
*/
static inline void IBPrintI(FILE *out, IBInterval* i, int digits, int mode, int verbose) {
  IBBasicWriteI(out,i,digits,mode,verbose);
}

/* conversion of a string representing a float to an interval */
static inline void IBStringToI(char *s, IBInterval* i) {
  IBBasicStringToI(s,i);
}

/* j := i1 + i2 */
static inline void IBAdditionII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAddII(j,i1,i2);
}

/* j := i1 + i2, i1 being a pointer to an interval point [x,x] */
static inline void IBAdditionRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAddRI(j,i1,i2);
}

/* j := i1 - i2 */
static inline void IBSubstractionII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSubII(j,i1,i2);
}

/* j := i1 - i2, i1 being a pointer to an interval point [x,x] */
static inline void IBSubstractionRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSubRI(j,i1,i2);
}

/* j := i1 - i2, i2 being a pointer to an interval point [x,x] */
static inline void IBSubstractionIR(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSubIR(j,i1,i2);
}

/* j := -i */
static inline void IBNegationI(IBInterval* j, IBInterval* i, IBInterval* useless) {
  IBBasicNegI(j,i,useless);
}

/* j := i1 * i2 */
static inline void IBMultiplicationII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicMulII(j,i1,i2);
}

/* j := i1 * i2, i1 being a pointer to an interval point [x,x] */
static inline void IBMultiplicationRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicMulRI(j,i1,i2);
}

/* j := i1 * i2, i1 being a pointer to an interval point [x,x], x<=0 */
static inline void IBMultiplicationRnegI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicMulRnegI(j,i1,i2);
}

/* j := i1 * i2, i1 being a pointer to an interval point [x,x], x>=0 */
static inline void IBMultiplicationRposI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicMulRposI(j,i1,i2);
}

/* j := i1 / i2 */
static inline void IBDivisionII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivII(j,i1,i2);
}

/* j := i1 / i2, i2 being a pointer to an interval point [x,x] */
static inline void IBDivisionIR(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivIR(j,i1,i2);
}

/* j := i1 / i2, i1 being a pointer to an interval point [x,x] */
static inline void IBDivisionRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivRI(j,i1,i2);
}

/* j := i1 / i2, i2 being a pointer to an interval point [x,x], x <= 0 */
static inline void IBDivisionIRneg(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivIRneg(j,i1,i2);
}

/* j := i1 / i2, i2 being a pointer to an interval point [x,x], x>=0 */
static inline void IBDivisionIRpos(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivIRpos(j,i1,i2);
}

/* j := i1 / i2, i1 being a pointer to an interval point [x,x], x<=0 */
static inline void IBDivisionRnegI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivRnegI(j,i1,i2);
}

/* j := i1 / i2, i1 being a pointer to an interval point [x,x], x>=0 */
static inline void IBDivisionRposI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicDivRposI(j,i1,i2);
}

/* j := square(i1) */
static inline void IBSquareI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSqrI(j,i1,i2);
}

/* j := square_root(i1) */
static inline void IBSquareRoot(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSqrtI(j,i1,i2);
}

/* j := i1 power i2, i2 being a pointer to an interval point [n,n], n natural */
static inline void IBPowerIR(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicPowI(j,i1,i2);
}

/* j := exp(i1) */
static inline void IBExponentialI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicExpI(j,i1,i2);
}

/* j := log(i1) with base e */
static inline void IBLogarithmI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicLogI(j,i1,i2);
}

/* j := min(i1,i2) */
static inline void IBMinimumOfII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicMinimumII(j,i1,i2);
}

/* j := max(i1,i2) */
static inline void IBMaximumOfII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicMaximumII(j,i1,i2);
}

/* j := cos(i1) */
static inline void IBCosineI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicCosI(j,i1,i2);
}

/* j := sin(i1) */
static inline void IBSineI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSinI(j,i1,i2);
}

/* j := tan(i1) */
static inline void IBTangentI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicTanI(j,i1,i2);
}

/* j := cosh(i1) */
static inline void IBCosHypI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicCoshI(j,i1,i2);
}

/* j := sinh(i1) */
static inline void IBSinHypI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicSinhI(j,i1,i2);
}

/* j := tanh(i1) */
static inline void IBTanHypI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicTanhI(j,i1,i2);
}

/* j := acos(i1) */
static inline void IBArcCosI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAcosI(j,i1,i2);
}

/* j := asin(i1) */
static inline void IBArcSinI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAsinI(j,i1,i2);
}

/* j := atan(i1) */
static inline void IBArcTanI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAtanI(j,i1,i2);
}

/* j := acosh(i1) */
static inline void IBArcCosHypI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAcoshI(j,i1,i2);
}

/* j := asinh(i1) */
static inline void IBArcSinHypI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAsinhI(j,i1,i2);
}

/* j := atanh(i1) */
static inline void IBArcTanHypI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBBasicAtanhI(j,i1,i2);
}

/*-- Computes (i1/i2) using the extended division over intervals
     Returns:
           1 if j := num/den
           2 if (j union k) := num/den */
static inline int IBExtendedDivisionII(IBInterval* j, IBInterval* k,
                                            IBInterval* i1, IBInterval* i2) {
  return IBBasicExtendedDivisionII(j,k,i1,i2);
}

/* j := j intersection (i1 / i2), / is the extended division over intervals
   Returns 1 if j is not modified, 0 otherwise */
static inline int IBExtendedDivisionInterII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  return IBBasicExtendedDivisionInterII(j,i1,i2);
}

/* relational n-th root of i */
static inline int IBNthRootRelationalI(IBInterval* j, IBInterval* k,
                                       IBInterval* i, IBInterval* n) {
  return IBBasicNthRootRelI(j,k,i,n);
}

/* relational hyperbolic cosine of i */
static inline int IBCoshRelationalI(IBInterval* j, IBInterval* k, IBInterval* i) {
  return IBBasicCoshRelI(j,k,i);
}

/* relational sine of i */
static inline int IBSinRelationalI(IBInterval* j, IBInterval* k, IBInterval* l,
                                   IBInterval* i) {
  return IBBasicSinRelI(j,k,l,i);
}

/* j := x * i, x>=0 */
static inline void IBMulRealposI(IBInterval* j, double x, IBInterval* i) {
  IBBasicMulRposIinternal(j,x,i);
}

/* j := x * i */
static inline void IBMulRealI(IBInterval* j, double x, IBInterval* i) {
  IBBasicMulRIinternal(j,x,i);
}

/* j := x / i, x>=0 */
static inline void IBDivRealposI(IBInterval* j, double x, IBInterval* i) {
  IBBasicDivRposIinternal(j,x,i);
}

/* j := i^n */
static inline void IBPowerIN(IBInterval* j, IBInterval* i, int n) {
  IBBasicPowIinternal(j,i,n);
}

/* j := j intersection (m - e/d), 0 not in d
   Returns 0 if j is not modified, 1 otherwise */
static inline int IBNewtonStepNonzeroII(IBInterval* j, IBInterval* m, IBInterval* e, IBInterval* d) {
  return IBBasicNewtonNonzeroII(j,m,e,d);
}

/* j := j intersection (m - e/d), 0 in d
   Returns 0 if j is not modified, 1 otherwise */
static inline int IBNewtonStepZeroII(IBInterval* j, IBInterval* m, IBInterval* e, IBInterval* d) {
  return IBBasicNewtonZeroII(j,m,e,d);
}

/* i := hull({pi}) */
static inline void IBSetEnclosePi(IBInterval* i) {
  IBBasicSetToPi(i);
}

/* i := hull({pi/2}) */
static inline void IBSetEncloseHalfPi(IBInterval* i) {
  IBBasicSetToHalfPi(i);
}

/* i := hull({ln(2)}) */
static inline void IBSetEncloseLn2(IBInterval* i) {
  IBBasicSetToLn2(i);
}

/* i := hull({e}) */
static inline void IBSetEncloseE(IBInterval* i) {
  IBBasicSetToE(i);
}

/* Initialization for the interval arithmetic module */
static inline void IBInitInterval() {
  IBBasicIntervalInit();
}

/* Profiling information for the interval arithmetic module */
static inline void IBProfileInterval() {
}


#endif
