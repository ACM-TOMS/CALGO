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
 * default_interval.h                                                       *
 ****************************************************************************/

#ifndef __default_interval_h
#define __default_interval_h


#include "default_fpu.h"
#include "profile.h"
#include <stdio.h>


/*------ definitions ------*/

typedef struct
{      double left;
       double right;
} IBBasicBounds;
typedef IBBasicBounds IBBasicItv[1];  /* using IBBasicItv wrt. IBBasicBounds is optimal */

/* Macros for type IBBasicBounds */
#define IBBasicLeftI(i)  i.left
#define IBBasicRightI(i) i.right

/* Macros for type IBBasicItv */
#define IBBasicMinI(i)  IBBasicLeftI(i[0])
#define IBBasicMaxI(i)  IBBasicRightI(i[0])

/* Symbols used for printing */
#define IBBasicPrintIntervalBounds   1
#define IBBasicPrintIntervalMidError 2


static IBBasicItv IBBasicItvConstPi;
static IBBasicItv IBBasicItvConst_2_Pi;
static IBBasicItv IBBasicItvConst_4_Pi;
static IBBasicItv IBBasicItvConst_1_Pi_2;
static IBBasicItv IBBasicItvConst_3_Pi_2;
static IBBasicItv IBBasicItvConst_5_Pi_2;
static IBBasicItv IBBasicItvConst_7_Pi_2;


/*------ macros ------*/

/* _XXX is for type IBBasicBounds, XXX is for type IBBasicItv */

/* i empty ? */
#define _IBBasicEmptyI(i)         (!(IBBasicLeftI(i)<=IBBasicRightI(i)))
#define  IBBasicEmptyI(i)         _IBBasicEmptyI(i[0])

/* i1=i2 ? */
#define _IBBasicIeqI(i1,i2)       ((IBBasicLeftI(i1)==IBBasicLeftI(i2)) && \
                                  (IBBasicRightI(i1)==IBBasicRightI(i2)))
#define  IBBasicIeqI(i1,i2)       _IBBasicIeqI(i1[0],i2[0])

/* i1!=i2 ? */
#define _IBBasicIdiffI(i1,i2)     ((IBBasicLeftI(i1)!=IBBasicLeftI(i2)) || \
                                  (IBBasicRightI(i1)!=IBBasicRightI(i2)))
#define  IBBasicIdiffI(i1,i2)     _IBBasicIdiffI(i1[0],i2[0])


#define _IBBasicDoubleInI(i,x)    (((x<IBBasicLeftI(i)) || (x>IBBasicRightI(i))) ? 0 : 1)
#define  IBBasicDoubleInI(i,x)    _IBBasicDoubleInI(i[0],x)

/* i=[r,r] ? */
#define _IBBasicIsDoubleI(i)      (IBBasicLeftI(i)==IBBasicRightI(i))
#define  IBBasicIsDoubleI(i)      _IBBasicIsDoubleI(i[0])

/* i=[n,n] ? */
#define _IBBasicIsIntegerI(i)     ( (IBBasicLeftI(i)==IBBasicRightI(i)) && \
                                    (IBBasicLeftI(i)==((int)IBBasicLeftI(i))) )
#define  IBBasicIsIntegerI(i)     _IBBasicIsIntegerI(i[0])

/* i=[0,0] ? */
#define _IBBasicZeroI(i)          ((IBBasicLeftI(i)==0.0) && (IBBasicRightI(i)==0.0))
#define  IBBasicZeroI(i)          _IBBasicZeroI(i[0])

/* i1 included in i2 ? */
#define _IBBasicIncludedII(i1,i2) ((IBBasicLeftI(i1)>=IBBasicLeftI(i2)) && \
                                  (IBBasicRightI(i1)<=IBBasicRightI(i2)))
#define  IBBasicIncludedII(i1,i2) _IBBasicIncludedII(i1[0],i2[0])

/* i1 and i2 disjoint ? */
#define _IBBasicDisjointII(i1,i2) ((IBBasicRightI(i1)<IBBasicLeftI(i2)) || \
                                  (IBBasicRightI(i2)<IBBasicLeftI(i1)))
#define  IBBasicDisjointII(i1,i2) _IBBasicDisjointII(i1[0],i2[0])


/* i has an infinite bound ? */
#define _IBBasicInfinite(i)       ((IBBasicLeftI(i)==IBBasicNegInfinity) || \
                                  ((IBBasicRightI(i)==IBBasicPosInfinity)))
#define  IBBasicInfinite(i)       _IBBasicInfinite(i[0])

/* i has no infinite bound ? */
#define _IBBasicFinite(i)         ((IBBasicLeftI(i)!=IBBasicNegInfinity) && \
                                  ((IBBasicRightI(i)!=IBBasicPosInfinity)))
#define  IBBasicFinite(i)         _IBBasicFinite(i[0])

/* width of i */
#define _IBBasicWidthI(i)         (IBBasicRightI(i)-IBBasicLeftI(i))
#define  IBBasicWidthI(i)         _IBBasicWidthI(i[0])

/* distance between i and j */
#define _IBBasicDistanceII(i,j)   ((IBBasicLeftI(j)-IBBasicLeftI(i))+ \
                                   (IBBasicRightI(i)-IBBasicRightI(j)))
#define  IBBasicDistanceII(i,j)   _IBBasicDistanceII(i[0],j[0])

/* real number at the center of i */
#define _IBBasicMidI(i)           IBBasicCenter(IBBasicLeftI(i),IBBasicRightI(i))
#define  IBBasicMidI(i)           _IBBasicMidI(i[0])      

/* see IBBasicThird for explanation */
#define _IBBasicThirdI(i)         IBBasicThird(IBBasicLeftI(i),IBBasicRightI(i))
#define  IBBasicThirdI(i)         _IBBasicThirdI(i[0])

/* see IBBasicTwoThirds for explanation */
#define _IBBasicTwoThirdsI(i)     IBBasicTwoThirds(IBBasicLeftI(i),IBBasicRightI(i))
#define  IBBasicTwoThirdsI(i)     _IBBasicTwoThirdsI(i[0])

/* i contains at most two floating-point numbers */
#define _IBBasicCanonicalI(i) \
   ( ((IBBasicLeftI(i)==IBBasicNegInfinity) && (IBBasicRightI(i)==IBBasicMinDouble)) || \
     ((IBBasicLeftI(i)==IBBasicMaxDouble) && (IBBasicRightI(i)==IBBasicPosInfinity)) || \
     (IBBasicRightI(i)<=IBBasicNextDouble(IBBasicLeftI(i))) )
#define  IBBasicCanonicalI(i)     _IBBasicCanonicalI(i[0])

/* i := empty interval */
#define _IBBasicSetEmptyI(i)      IBBasicRightI(i) = IBBasicNegInfinity; \
                                  IBBasicLeftI(i) = IBBasicPosInfinity
#define  IBBasicSetEmptyI(i)      _IBBasicSetEmptyI(i[0])

/* i := [x1,x2] */
#define _IBBasicSetI(i,x1,x2)     IBBasicLeftI(i)=x1 ; IBBasicRightI(i)=x2
#define  IBBasicSetI(i,x1,x2)     _IBBasicSetI(i[0],x1,x2)

/* i := source */
#define _IBBasicCopyI(i,source)   IBBasicLeftI(i) = IBBasicLeftI(source); \
                                  IBBasicRightI(i) = IBBasicRightI(source)
#define  IBBasicCopyI(i,source)   _IBBasicCopyI(i[0],source[0])



/*------ interval functions ------*/

/* Initialization of interval module; must be called before use */
void IBBasicIntervalInit();

/* allocation */
IBBasicBounds *IBBasicNewI        ();
IBBasicBounds *IBBasicNewLargestI ();
IBBasicBounds *IBBasicNewCopyI    (IBBasicItv i);
IBBasicBounds *IBBasicSetNewI     (double x1, double x2);

/* operations */
void           IBBasicToLargestI  (IBBasicItv i);
void           IBBasicToIntegerI  (IBBasicItv i);
void           IBBasicAbsI        (IBBasicItv Result, IBBasicItv i);

/* output */
void           IBBasicWriteI      (FILE *out, IBBasicItv i, int digits, int mode, int verbose);

/* conversion */
void           IBBasicStringToI   (char* s, IBBasicItv i);


/* intersection */
void           IBBasicInterII     (IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);

/* arithmetic operations and elementary functions */
void IBBasicAddII    (IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);      /* i1+i2   */
void IBBasicAddRI    (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* x+i     */
void IBBasicSubII    (IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);      /* i1-i2   */
void IBBasicSubRI    (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* x-i     */
void IBBasicSubIR    (IBBasicItv Result, IBBasicItv i,  IBBasicItv x);       /* i-x     */
void IBBasicNegI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* -i      */
void IBBasicMulII    (IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);      /* i1*i2   */
void IBBasicMulRI    (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* x*i     */
void IBBasicMulRnegI (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* -r*i    */
void IBBasicMulRposI (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* +r*i    */
void IBBasicDivII    (IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);      /* i1/i2   */
void IBBasicDivIR    (IBBasicItv Result, IBBasicItv i,  IBBasicItv x);       /* i/x     */
void IBBasicDivRI    (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* x/i     */
void IBBasicDivIRneg (IBBasicItv Result, IBBasicItv i,  IBBasicItv x);       /* i/-r    */
void IBBasicDivIRpos (IBBasicItv Result, IBBasicItv i,  IBBasicItv x);       /* i/+r    */
void IBBasicDivRnegI (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* -r/i    */
void IBBasicDivRposI (IBBasicItv Result, IBBasicItv x,  IBBasicItv i);       /* +r/i    */
void IBBasicSqrI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* i**2    */
void IBBasicSqrtI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* sqrt(i) */
void IBBasicPowI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv n);       /* i**n    */
void IBBasicExpI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* exp(i)  */
void IBBasicLogI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* log(i)  */
void IBBasicMinimumII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);      /* min(i1,i2) */
void IBBasicMaximumII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2);      /* max(i1,i2) */
void IBBasicCosI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* cos(i) */
void IBBasicSinI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* sin(i) */
void IBBasicTanI     (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* tan(i) */
void IBBasicCoshI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* cosh(i) */
void IBBasicSinhI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* sinh(i) */
void IBBasicTanhI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* tanh(i) */
void IBBasicAcosI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* acos(i) */
void IBBasicAsinI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* asin(i) */
void IBBasicAtanI    (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* atan(i) */
void IBBasicAcoshI   (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* acosh(i) */
void IBBasicAsinhI   (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* asinh(i) */
void IBBasicAtanhI   (IBBasicItv Result, IBBasicItv i,  IBBasicItv useless); /* atanh(i) */


/*-- auxiliary functions --*/
void IBBasicMulRposIinternal (IBBasicItv Result, double x, IBBasicItv i);
void IBBasicMulRIinternal    (IBBasicItv Result, double x, IBBasicItv i);
void IBBasicDivRposIinternal (IBBasicItv Result, double x, IBBasicItv i);
void IBBasicPowIinternal     (IBBasicItv Result, IBBasicItv i, int n);
void IBBasicAddRIinternal    (IBBasicItv Result, double x, IBBasicItv i);
void IBBasicSubRIinternal    (IBBasicItv Result, double x, IBBasicItv i);
void IBBasicSubIRinternal    (IBBasicItv Result, IBBasicItv i, double x);


/*-- Result := Result inter (num/den) where num/den uses the extended
     division over intervals -> used in Gauss-Seidel iterations --*/
int IBBasicExtendedDivisionInterII (IBBasicItv Result, IBBasicItv num, IBBasicItv den);

/*-- Computes (num/den) using the extended division over intervals
     Returns:
           1 if num/den = Result1
           2 if num/den = Result1 union Result2 --*/
int IBBasicExtendedDivisionII(IBBasicItv Result1, IBBasicItv Result2,
                              IBBasicItv num, IBBasicItv den);

/*-- Computes the relational n-root of i
     Returns:
           0 if the result is empty (when i<0)
           1 if the n-root of i = Result1
           2 if the n-root of i = Result1 union Result2 --*/
int IBBasicNthRootRelI(IBBasicItv Result1, IBBasicItv Result2,
                       IBBasicItv i, IBBasicItv n);

/*-- Computes the relational hyperbolic cosine of i
     Returns:
           0 if the result is empty (when i<1.0)
           1 if the result = Result1
           2 if the result = Result1 union Result2 --*/
int IBBasicCoshRelI(IBBasicItv Result1, IBBasicItv Result2, IBBasicItv i);

/*-- Computes the relational sine of i => result in [-pi, +pi]
     Returns:
           0 if the result is empty (when i<1.0)
           1 if the result = Result1
           2 if the result = Result1 union Result2
           3 if the result = Result1 union Result2 union Result3 --*/
int IBBasicSinRelI(IBBasicItv Result1, IBBasicItv Result2, IBBasicItv Result3, IBBasicItv i);

/*-- i := smallest interval containing Pi */
void IBBasicSetToPi(IBBasicItv i);

/*-- i := smallest interval containing Pi/2 */
void IBBasicSetToHalfPi(IBBasicItv i);

/*-- i := smallest interval containing Ln(2) */
void IBBasicSetToLn2(IBBasicItv i);

/*-- i := smallest interval containing e */
void IBBasicSetToE(IBBasicItv i);


/*-- Result <- mid - eval/deriv   such that deriv does not contain 0.0
     Returns 0 if Result is not modified, 1 otherwise
     Used in Newton narrowing operator --*/
int IBBasicNewtonNonzeroII(IBBasicItv Result, IBBasicItv mid, IBBasicItv eval, IBBasicItv deriv);

/*-- Equivalent function for the case 0.0 in deriv */
int IBBasicNewtonZeroII(IBBasicItv Result, IBBasicItv mid, IBBasicItv eval, IBBasicItv deriv);

#endif
