      REAL FUNCTION alog2(x)
************************************************************************
*     (Single-precision log_2(x))
*     Return the logarithm to the base 2 of x.
*     (01-May-2000)
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           alog
*
*     Built-in functions
*
      REAL                alog
*
*     Argument variables
*
      REAL                x
*
      INCLUDE 'algtwo.inc'
*
*     Systems with base-2 floating-point arithmetic (IEEE 754,
*     Compaq/DEC VAX, Cray, ...):
*
      alog2 = alog(x) * alg2in
*
*     Systems with base-16 floating-point arithmetic (IBM S/360): this
*     avoids the loss of three bits in the alg2in value: algtwo is
*     stored to full precision.  Alternatively, to avoid a slow
*     division, split algtwo into two parts, the larger exactly
*     representable:
*
*         alog2 = 0.6875e+00 +
*                 5.64718055994530941723212145817656807550013436026000e-03
*               = alg2hi + alg2lo
*
*     and then compute the result with two multiplies and one add,
*     where the parentheses are essential:
*
*         a = alog(x)
*         alog2 = (a*alg2hi) + (a*alg2lo)
*
*     Here is the slow way:
*
*     alog2 = alog(x) / algtwo
*
      END
