      DOUBLE PRECISION FUNCTION dlog2(x)
************************************************************************
*     (Double-precision log_2(x))
*     Return the logarithm to the base 2 of x.
*     (01-May-2000)
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dlog
*
*     Built-in functions
*
      DOUBLE PRECISION    dlog
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
      INCLUDE 'dlgtwo.inc'
*
*     Systems with base-2 floating-point arithmetic (IEEE 754,
*     Compaq/DEC VAX, Cray, ...):
*
      dlog2 = dlog(x) * dlg2in
*
*     Systems with base-16 floating-point arithmetic (IBM S/360): this
*     avoids the loss of three bits in the dlg2in value: dlgtwo is
*     stored to full precision.  Alternatively, to avoid a slow
*     division, split dlgtwo into two parts, the larger exactly
*     representable:
*
*         dlog2 = 0.693145751953125d+00 + 5.76999904754328571215d-08
*               = dlg2hi + dlg2lo
*
*     and then compute the result with two multiplies and one add,
*     where the parentheses are essential:
*
*         d = dlog(x)
*         dlog2 = (d*dlg2hi) + (d*dlg2lo)
*
*     Here is the slow way:
*
*     dlog2 = dlog(x) / dlgtwo
*
      END
