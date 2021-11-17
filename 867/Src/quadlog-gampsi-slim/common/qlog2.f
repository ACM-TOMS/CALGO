      REAL*16 FUNCTION qlog2(x)
************************************************************************
*     (Quadruple-precision log_2(x))
*     Return the logarithm to the base 2 of x.
*     (01-May-2000)
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           qlog
*
*     Built-in functions
*
      REAL*16             qlog
*
*     Argument variables
*
      REAL*16             x
*
      INCLUDE 'qlgtwo.inc'
*
*     Systems with base-2 floating-point arithmetic (IEEE 754,
*     Compaq/DEC VAX, Cray, ...):
*
      qlog2 = qlog(x) * qlg2in
*
*     Systems with base-16 floating-point arithmetic (IBM S/360): this
*     avoids the loss of three bits in the dlg2in value: dlgtwo is
*     stored to full precision.  Alternatively, to avoid a slow
*     division, split dlgtwo into two parts, the larger exactly
*     representable:
*
*         qlog2 = 0.693145751953125q+00 +
*                 5.76999904754328571214581765680755001e-08
*               = dlg2hi + dlg2lo
*
*     and then compute the result with two multiplies and one add,
*     where the parentheses are essential:
*
*         q = qlog(x)
*         qlog2 = (q*dlg2hi) + (q*dlg2lo)
*
*     Here is the slow way:
*
*     qlog2 = qlog(x) / dlgtwo
*
      END
