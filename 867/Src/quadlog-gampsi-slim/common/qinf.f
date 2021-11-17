      REAL*16 FUNCTION  qinf()
************************************************************************
*     (Quadruple-precision +Infinity)
*     Return quadruple-precision +Infinity, or else on non-IEEE 754
*     systems, the largest representable floating-point number.
*
*     For IEEE 754 systems, each call to this function intentionally
*     produces a trappable zero divide, rather than saving the
*     computed value on the first call, and then just returning the
*     saved value on subsequent calls.
*
*     This function exists because of at least one abberant software
*     implementation of quadruple-precision arithmetic (on IBM RS/6000
*     AIX 4.x), which produces NaN, instead of Infinity, for the square
*     of large numbers.  Fortunately, it correctly produces Infinity
*     for 1.0/0.0, so that is how we generate it here.
*
*     Relegating the computation of Infinity to a separate function
*     also provides a convenient single debugger breakpoint location.
*     [12-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            qstorf
*
      REAL*16             qstorf
*
*     Parameter variables
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
      qinf = one / qstorf(zero)
*
      END
