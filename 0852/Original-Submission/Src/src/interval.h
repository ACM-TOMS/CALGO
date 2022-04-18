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
 * interval.h                                                               *
 ****************************************************************************/

#ifndef __interval_h
#define __interval_h

#include "profile.h"
#include "interval_interface.h"
#include "config.h"
#include <stdio.h>
#include <math.h>

/* variables used for profiling */
#if SOFTWARE_PROFILE
unsigned long
  IBNumberAdd,
  IBNumberSub,
  IBNumberMul,
  IBNumberDiv,
  IBNumberExtDiv,
  IBNumberNthRootRel,
  IBNumberSqr,
  IBNumberSqrt,
  IBNumberPow,
  IBNumberExp,
  IBNumberLog,
  IBNumberMin,
  IBNumberMax,
  IBNumberCos,
  IBNumberSin,
  IBNumberTan,
  IBNumberCosh,
  IBNumberSinh,
  IBNumberTanh,
  IBNumberAcos,
  IBNumberAsin,
  IBNumberAtan,
  IBNumberAcosh,
  IBNumberAsinh,
  IBNumberAtanh,
  IBNumberCoshRel,
  IBNumberSinRel;
#endif


/* mathematical functions */
extern double ceil();
extern double floor();
extern double log2();
extern double log();
extern double sqrt();
extern double cos();
extern double acos();
extern double sin();
extern double tan();
extern double exp();
extern double fabs();

/* operations on integers and real numbers */
#define IBAbs(x)   ((x<=0.0) ? -(x) : x)  /* absolute value of a real */
#define IBOdd(n)   (((n)%2)==1)           /* n is odd ? */
#define IBEven(n)  (((n)%2)==0)           /* n is even ? */
#define IBMin(x,y) (((x)<(y))? x : y)     /* minimum of x and y */
#define IBMax(x,y) (((x)<(y))? y : x)     /* maximum of x and y */

#define IBNextDouble(x) IBNextReal(x)     /* successor of real x */
#define IBPrevDouble(x) IBPrevReal(x)     /* predecessor of real x */

#define IBStrToInt(s)    atoi(s)          /* conversion of string to integer */
#define IBStrToDouble(s) atof(s)          /* conversion of string to double */

/* modes for interval printing */
#define IBPrintIntervalBounds       1
#define IBPrintIntervalMidError     2

/* interval type */
typedef IBInterval                  IBItv[1];

/* left and right bounds */
#define IBMinI(i)                   IBLeftI(i[0])
#define IBMaxI(i)                   IBRightI(i[0])

/* test functions */
#define IBEmptyI(i)                 IBIsEmptyI(i[0])             /* i empty ? */
#define IBIeqI(i1,i2)               IBIsEqualII(i1[0],i2[0])     /* i1==i2 ? */
#define IBIdiffI(i1,i2)             IBIsDifferentII(i1[0],i2[0]) /* i1!=i2 ? */
#define IBDoubleInI(i,x)            IBIsDoubleInI(i[0],x)        /* x in i ? */
#define IBIsDoubleI(i)              IBIsIntervalPoint(i[0])      /* i==[r,r] ? */
#define IBIsIntegerI(i)             IBIsIntervalIntPoint(i[0])   /* i==[n,n] ? */
#define IBIsZeroI(i)                IBIsReducedToZeroI(i[0])     /* i==[0,0] ? */
#define IBIncludedII(i1,i2)         IBIsIncludedII(i1[0],i2[0])  /* i1 included in i2 ? */
#define IBCanonicalI(i)             IBIsCanonicalI(i[0])         /* at most two floats in i */
#define IBDisjointII(i1,i2)         IBIsDisjointII(i1[0],i2[0])  /* i1 and i2 are disjoint */
#define IBInfiniteI(i)              IBIsInfiniteI(i[0])          /* at least one infinite bound? */

/* operations on intervals */
#define IBWidthI(i)                 IBWidthOfI(i[0])
#define IBDistanceII(i1,i2)         IBDistanceBetweenII(i1[0],i2[0])
#define IBMidI(i)                   IBMidpointOfI(i[0])
#define IBThirdI(i)                 IBThirdOfI(i[0])
#define IBTwoThirdsI(i)             IBTwoThirdsOfI(i[0])

/* modification */
#define IBSetEmptyI(i)              IBSetToEmptyI(i[0])        /* i := emptyset */
#define IBToLargestI(i)             IBSetToRealDomain(i)       /* i := [-oo,+oo] */
#define IBToIntegerI(i)             IBSetToIntegerDomain(i)    /* i := [ceil(inf i),floor(sup i)] */
#define IBSetI(i,x1,x2)             IBSetBoundsOfI(i[0],x1,x2) /* i := [x1,x2] */
#define IBCopyI(i,source)           IBCopyII(i[0],source[0])   /* i := source */

/* allocation in memory */
#define IBNewI()                    IBCreateNewI()
#define IBNewLargestI()             IBCreateNewRealDomainI()
#define IBNewCopyI(i)               IBCreateAndCopyNewI(i)
#define IBSetNewI(x1,x2)            IBCreateAndSetNewI(x1,x2)

/* enclosures of mathematical constants */
#define IBSetToPiI(i)               IBSetEnclosePi(i)
#define IBSetToHalfPiI(i)           IBSetEncloseHalfPi(i)
#define IBSetToLn2I(i)              IBSetEncloseLn2(i)
#define IBSetToEI(i)                IBSetEncloseE(i)

/* intersection */
#define IBInterII(j,i1,i2)          IBIntersectionII(j,i1,i2)

/* printing */
#define IBWriteI(f,i,d,p)           IBPrintI(f,i,d,p,0)
#define IBWriteIverb(f,i,d,p)       IBPrintI(f,i,d,p,1)



/* conversion of string to interval; the string represents a float
   The result is correctly rounded interval */
static inline void IBStrToI(char *s, IBInterval*i) {
  IBStringToI(s,i);
}

/* arithmetic operations and elementary functions:
   - IBDefxxx is the name of the C function
   - IBxxx is the call of IBDefxxx
*/


static inline void IBAddII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBAdditionII(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAdd++;
#endif
}
#define IBDefAddII IBAddII

static inline void IBAddRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBAdditionRI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAdd++;
#endif
}
#define IBDefAddRI IBAddRI

static inline void IBSubII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSubstractionII(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSub++;
#endif
}
#define IBDefSubII IBSubII

static inline void IBSubRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSubstractionRI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSub++;
#endif
}
#define IBDefSubRI IBSubRI

static inline void IBSubIR(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSubstractionIR(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSub++;
#endif
}
#define IBDefSubIR IBSubIR

static inline void IBNegI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBNegationI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSub++;
#endif
}
#define IBDefNegI IBNegI

static inline void IBMulII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBMultiplicationII(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberMul++;
#endif
}
#define IBDefMulII IBMulII

static inline void IBMulRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBMultiplicationRI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberMul++;
#endif
}
#define IBDefMulRI IBMulRI

static inline void IBMulRnegI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBMultiplicationRnegI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberMul++;
#endif
}
#define IBDefMulRnegI IBMulRnegI

static inline void IBMulRposI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBMultiplicationRposI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberMul++;
#endif
}
#define IBDefMulRposI IBMulRposI

static inline void IBDivII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionII(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivII IBDivII

static inline void IBDivIR(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionIR(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivIR IBDivIR

static inline void IBDivRI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionRI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivRI IBDivRI

static inline void IBDivIRneg(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionIRneg(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivIRneg IBDivIRneg

static inline void IBDivIRpos(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionIRpos(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivIRpos IBDivIRpos

static inline void IBDivRnegI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionRnegI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivRnegI IBDivRnegI

static inline void IBDivRposI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBDivisionRposI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}
#define IBDefDivRposI IBDivRposI

static inline void IBSqrI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSquareI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSqr++;
#endif
}
#define IBDefSqrI IBSqrI

static inline void IBSqrtI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSquareRoot(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSqrt++;
#endif
}
#define IBDefSqrtI IBSqrtI

static inline void IBPowI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBPowerIR(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberPow++;
#endif
}
#define IBDefPowI IBPowI

static inline void IBExpI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBExponentialI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberExp++;
#endif
}
#define IBDefExpI IBExpI

static inline void IBLogI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBLogarithmI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberLog++;
#endif
}
#define IBDefLogI IBLogI

static inline void IBMinimumII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBMinimumOfII(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberMin++;
#endif
}
#define IBDefMinimumII IBMinimumII

static inline void IBMaximumII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBMaximumOfII(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberMax++;
#endif
}
#define IBDefMaximumII IBMaximumII


static inline void IBCosI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBCosineI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberCos++;
#endif
}
#define IBDefCosI IBCosI

static inline void IBSinI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSineI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSin++;
#endif
}
#define IBDefSinI IBSinI

static inline void IBTanI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBTangentI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberTan++;
#endif
}
#define IBDefTanI IBTanI

static inline void IBCoshI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBCosHypI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberCosh++;
#endif
}
#define IBDefCoshI IBCoshI

static inline void IBSinhI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBSinHypI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberSinh++;
#endif
}
#define IBDefSinhI IBSinhI

static inline void IBTanhI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBTanHypI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberTanh++;
#endif
}
#define IBDefTanhI IBTanhI

static inline void IBAcosI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBArcCosI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAcos++;
#endif
}
#define IBDefAcosI IBAcosI

static inline void IBAsinI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBArcSinI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAsin++;
#endif
}
#define IBDefAsinI IBAsinI

static inline void IBAtanI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBArcTanI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAtan++;
#endif
}
#define IBDefAtanI IBAtanI

static inline void IBAcoshI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBArcCosHypI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAcosh++;
#endif
}
#define IBDefAcoshI IBAcoshI

static inline void IBAsinhI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBArcSinHypI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAsinh++;
#endif
}
#define IBDefAsinhI IBAsinhI

static inline void IBAtanhI(IBInterval* j, IBInterval* i1, IBInterval* i2) {
  IBArcTanHypI(j,i1,i2);
#if SOFTWARE_PROFILE
  IBNumberAtanh++;
#endif
}
#define IBDefAtanhI IBAtanhI

static inline int IBExtDivII(IBInterval* j, IBInterval* k, IBInterval* i1, IBInterval* i2) {
#if SOFTWARE_PROFILE
  IBNumberExtDiv++;
#endif
  return IBExtendedDivisionII(j,k,i1,i2);
}

static inline int IBExtDivInterII(IBInterval* j, IBInterval* i1, IBInterval* i2) {
#if SOFTWARE_PROFILE
  IBNumberExtDiv++;
#endif
  return IBExtendedDivisionInterII(j,i1,i2);
}

static inline int IBNthRootRelI(IBInterval* j, IBInterval* k,
                                IBInterval* i, IBInterval* n) {
#if SOFTWARE_PROFILE
  IBNumberNthRootRel++;
#endif
  return IBNthRootRelationalI(j,k,i,n);
}

static inline int IBCoshRelI(IBInterval* j, IBInterval* k, IBInterval* i) {
#if SOFTWARE_PROFILE
  IBNumberCoshRel++;
#endif
  return IBCoshRelationalI(j,k,i);
}

static inline int IBSinRelI(IBInterval* j, IBInterval* k, IBInterval* l,
                            IBInterval* i) {
#if SOFTWARE_PROFILE
  IBNumberSinRel++;
#endif
  return IBSinRelationalI(j,k,l,i);
}

static inline void IBMulRposIinternal(IBInterval* j, double x, IBInterval* i2) {
  IBMulRealposI(j,x,i2);
#if SOFTWARE_PROFILE
  IBNumberMul++;
#endif
}

static inline void IBMulRIinternal(IBInterval* j, double x, IBInterval* i) {
  IBMulRealI(j,x,i);
#if SOFTWARE_PROFILE
  IBNumberMul++;
#endif
}

static inline void IBDivRposIinternal(IBInterval* j, double x, IBInterval* i) {
  IBDivRealposI(j,x,i);
#if SOFTWARE_PROFILE
  IBNumberDiv++;
#endif
}

static inline void IBPowIinternal(IBInterval* j, IBInterval* i, int n) {
  IBPowerIN(j,i,n);
#if SOFTWARE_PROFILE
  IBNumberPow++;
#endif
}

#define IBNewtonNonzeroII(j,m,e,d)  IBNewtonStepNonzeroII(j,m,e,d)
#define IBNewtonZeroII(j,m,e,d)     IBNewtonStepZeroII(j,m,e,d)

/* generic type for binary interval operation */
typedef void (* IBEvalOpI)(IBItv, IBItv, IBItv);


/* initialization of the interval arithmetic module */
void IBInitIntervalConstants();

static inline void IBInitIA() {
#if SOFTWARE_PROFILE
  IBNumberAdd = 
  IBNumberSub = 
  IBNumberMul = 
  IBNumberDiv = 
  IBNumberDiv = 
  IBNumberSqr = 
  IBNumberSqrt = 
  IBNumberPow = 
  IBNumberExp = 
  IBNumberLog =
  IBNumberMin =
  IBNumberMax =
  IBNumberCos =
  IBNumberSin =
  IBNumberTan =
  IBNumberCosh =
  IBNumberSinh =
  IBNumberTanh =
  IBNumberAcos =
  IBNumberAsin =
  IBNumberAtan =
  IBNumberAcosh =
  IBNumberAsinh =
  IBNumberAtanh =
  IBNumberCoshRel =
  IBNumberSinRel
  = 0 ;
#endif

  IBInitIntervalConstants();
  IBInitInterval();
}


/* output of profiling information */
static inline int _IBNbDigits(unsigned long x) {
  int d = 1;
  while (x>=10) { ++d; x/=10;}
  return d;
}

static inline void _IBprintlong(char *s, unsigned long x, int esp) {
  int d = _IBNbDigits(x),
      nbcar = d+((d-1)/3),
      i = IBMax(nbcar,esp),
      j,
      k,
      l = 0;

  s[i--] = '\0';
  for (j=0; j<=i-nbcar; ++j) {
    s[j] = ' ';
  }

  for (k=i; k>=j; --k) {
    if (l==3) {
      s[k] = ',';
      l=0;
    }
    else {
      s[k] = '0' + (x % 10);
      x /= 10;
      ++l;
    }
  }
}

#if SOFTWARE_PROFILE
static inline void IBProfileIA() {
  char s[20];
  int esp = 0, n;
  unsigned long total;
  printf("  Number of interval operations:\n");

  total = IBNumberAdd+IBNumberSub+IBNumberMul+IBNumberDiv+IBNumberExtDiv+
          IBNumberNthRootRel+IBNumberSqr+IBNumberSqrt+IBNumberPow+IBNumberExp+
          IBNumberLog+IBNumberMin+IBNumberMax+
          IBNumberCos+IBNumberSin+IBNumberTan+
          IBNumberCosh+IBNumberSinh+IBNumberTanh+
          IBNumberAcos+IBNumberAsin+IBNumberAtan+
          IBNumberAcosh+IBNumberAsinh+IBNumberAtanh+
          IBNumberCoshRel+IBNumberSinRel;

  if (esp<(n=_IBNbDigits(IBNumberAdd)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberSub)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberMul)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberDiv)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberExtDiv))) esp=n;
  if (esp<(n=_IBNbDigits(IBNumberNthRootRel))) esp=n;
  if (esp<(n=_IBNbDigits(IBNumberSqr)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberSqrt)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberPow)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberExp)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberLog)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberMin)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberMax)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberCos)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberSin)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberTan)))    esp=n;
  if (esp<(n=_IBNbDigits(IBNumberCosh)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberSinh)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberTanh)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberAcos)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberAsin)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberAtan)))   esp=n;
  if (esp<(n=_IBNbDigits(IBNumberAcosh)))  esp=n;
  if (esp<(n=_IBNbDigits(IBNumberAsinh)))  esp=n;
  if (esp<(n=_IBNbDigits(IBNumberAtanh)))  esp=n;
  if (esp<(n=_IBNbDigits(IBNumberCoshRel))) esp=n;
  if (esp<(n=_IBNbDigits(IBNumberSinRel)))  esp=n;

  if (esp<(n=_IBNbDigits(total)))          esp=n;


  esp+=(esp-1)/3;

  if (IBNumberAdd>0) {
    _IBprintlong(s,IBNumberAdd,esp);
    printf("             add: %s\n",s);
  }

  if (IBNumberSub>0) {
    _IBprintlong(s,IBNumberSub,esp);
    printf("             sub: %s\n",s);
  }

  if (IBNumberMul>0) {
    _IBprintlong(s,IBNumberMul,esp);
    printf("             mul: %s\n",s);
  }

  if (IBNumberDiv>0) {
    _IBprintlong(s,IBNumberDiv+IBNumberExtDiv,esp);
    printf("             div: %s\n",s);
  }

  if (IBNumberExtDiv>0) {
    _IBprintlong(s,IBNumberExtDiv,esp);
    printf("            ediv: %s\n",s);
  }

  if (IBNumberSqr>0) {
    _IBprintlong(s,IBNumberSqr,esp);
    printf("             sqr: %s\n",s);
  }

  if (IBNumberSqrt>0) {
    _IBprintlong(s,IBNumberSqrt,esp);
    printf("            sqrt: %s\n",s);
  }

  if (IBNumberNthRootRel>0) {
    _IBprintlong(s,IBNumberNthRootRel,esp);
    printf("       n-th root: %s\n",s);
  }

  if (IBNumberPow>0) {
    _IBprintlong(s,IBNumberPow,esp);
    printf("             pow: %s\n",s);
  }

  if (IBNumberExp>0) {
    _IBprintlong(s,IBNumberExp,esp);
    printf("             exp: %s\n",s);
  }

  if (IBNumberLog>0) {
    _IBprintlong(s,IBNumberLog,esp);
    printf("             log: %s\n",s);
  }

  if (IBNumberMin>0) {
    _IBprintlong(s,IBNumberMin,esp);
    printf("             min: %s\n",s);
  }

  if (IBNumberMax>0) {
    _IBprintlong(s,IBNumberMax,esp);
    printf("             max: %s\n",s);
  }

  if (IBNumberCos>0) {
    _IBprintlong(s,IBNumberCos,esp);
    printf("             cos: %s\n",s);
  }

  if (IBNumberSin>0) {
    _IBprintlong(s,IBNumberSin,esp);
    printf("             sin: %s\n",s);
  }

  if (IBNumberTan>0) {
    _IBprintlong(s,IBNumberTan,esp);
    printf("             tan: %s\n",s);
  }

  if (IBNumberCosh>0) {
    _IBprintlong(s,IBNumberCosh,esp);
    printf("            cosh: %s\n",s);
  }

  if (IBNumberSinh>0) {
    _IBprintlong(s,IBNumberSinh,esp);
    printf("            sinh: %s\n",s);
  }

  if (IBNumberTanh>0) {
    _IBprintlong(s,IBNumberTanh,esp);
    printf("            tanh: %s\n",s);
  }

  if (IBNumberAcos>0) {
    _IBprintlong(s,IBNumberAcos,esp);
    printf("            acos: %s\n",s);
  }

  if (IBNumberAsin>0) {
    _IBprintlong(s,IBNumberAsin,esp);
    printf("            asin: %s\n",s);
  }

  if (IBNumberAtan>0) {
    _IBprintlong(s,IBNumberAtan,esp);
    printf("            atan: %s\n",s);
  }

  if (IBNumberAcosh>0) {
    _IBprintlong(s,IBNumberAcosh,esp);
    printf("           acosh: %s\n",s);
  }

  if (IBNumberAsinh>0) {
    _IBprintlong(s,IBNumberAsinh,esp);
    printf("           asinh: %s\n",s);
  }

  if (IBNumberAtanh>0) {
    _IBprintlong(s,IBNumberAtanh,esp);
    printf("           atanh: %s\n",s);
  }

  if (IBNumberCoshRel>0) {
    _IBprintlong(s,IBNumberCoshRel,esp);
    printf("         coshrel: %s\n",s);
  }

  if (IBNumberSinRel>0) {
    _IBprintlong(s,IBNumberSinRel,esp);
    printf("          sinrel: %s\n",s);
  }

  if (total>0) {
    printf("                  ");
    for( n=0; n<strlen(s); ++n )
    {
      printf("-");
    }
    _IBprintlong(s,total,esp);
    printf("\n           total: %s\n",s);
  }

  if( total==0 )
  {
    printf("     0\n");
  }

  IBProfileInterval();
}
#endif

#endif
