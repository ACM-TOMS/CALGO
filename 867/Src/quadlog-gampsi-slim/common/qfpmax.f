      REAL*16 FUNCTION  qfpmax()
************************************************************************
*     (Quadruple-precision floating-point maximum)
*     Return the largest finite representable quadruple-precision
*     floating-point number.
*     (03-Jun-2000)
************************************************************************
*
*     Parameter variables
*
      REAL*16             half
      PARAMETER           (half = 0.5q+00)
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
      REAL*16             two
      PARAMETER           (two = 2.0q+00)
*
      REAL*16             three
      PARAMETER           (three = 3.0q+00)
*
      REAL*16             four
      PARAMETER           (four = 4.0q+00)
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
*     Local variables
*
      REAL*16             arg,         fpmax,       w,           x
      REAL*16             x2,          x3,          y,           z
      REAL*16             zplusy
*
      INTEGER             ntimes,      pmax
      LOGICAL             first,       qisInf,      qisNaN
*
*     INCLUDE 'stdio.inc'
*
      SAVE                first,       fpmax
*
      DATA first /.TRUE./
*
      qisInf(arg) = (arg .NE. zero) .AND. (two*arg .EQ. arg)
      qisNaN(arg) = (arg .NE. arg)
*
      IF (first) THEN
           first = .FALSE.
*
*          Compute the largest finite representable power of two by
*          successive doubling.  This code CRUCIALLY DEPENDS on IEEE 754
*          nonstop computing: an overflow must not terminate execution!
*          The calls to qstore() are critical: they foil machines (e.g.,
*          Honeywell, Intel x86, Motorola 68K) in which the
*          floating-point registers have range and/or precision beyond
*          that of floating-point values in memory.
*
*          On well-designed floating-point systems, this code would be
*          relatively simple, and would generate only two overflows to
*          Infinity (in the computation of x2 and x3 at the start of the
*          loop 10 body), and never a NaN.
*
*          Unfortunately, IBM RS/6000 AIX and SGI IRIX quadruple
*          precision implementations are anomalous, and require more
*          complex coding.  Doubling by addition fails on the IBM
*          system, where near the overflow limit, the software
*          implementation produces NaNQ instead of Infinity: doubling by
*          (slower) multiplication eliminates this misbehavior.  Because
*          of the possibility of NaNs, the loop termination conditions
*          must also contain checks for them.
*
           x = one
           ntimes = 0
   10      x2 = two*x
           CALL qstore(x2)
           x3 = three*x
           CALL qstore(x3)
           IF ((x3 .NE. x2) .AND.
     X          .NOT.qisNaN(x) .AND.
     X          .NOT.qisNaN(x2) .AND.
     X          .NOT.qisNaN(x3)) THEN
                CALL qstore(x)
                x = two*x
                ntimes = ntimes + 1
                GO TO 10
           END IF
*
*          Compute the large finite representable number by successive
*          appending of a one bit to the significand.  We have to work
*          with (x/2) instead of x, so that the last addition does not
*          overflow.
*
           z = half * x
           y = half * z
           pmax = 1
   20      zplusy = z + y
           CALL qstore(zplusy)
           IF (.NOT.qisNaN(zplusy) .AND.
     X          .NOT.qisNaN(y) .AND.
     X          .NOT.qisNaN(z) .AND.
     X          (zplusy .NE. z) .AND.
     X          (zplusy .LT. x)) THEN
                CALL qstore(z)
                CALL qstore(y)
                z = zplusy
                y = half * y
                pmax = pmax + 1
*               WRITE (stdout,'(i5, 1p, e45.35e4, 2x, z32.32)') pmax,z,z
*
*              The tricky loop continuation test below handles the
*              anomalous IBM RS/6000 case.  At iteration 52, we have
*
*                  z = 8.988465674311579039686448570265e307
*                    = 7fdfffff ffffffff 00000000 00000000
*                  y = 4.989600773836799529140931782592e291
*                    = 7c800000 00000000 00000000 00000000
*
*              and then
*
*                  (z + y) = 8.988465674311579039686448570265e307
*                          = 7fe00000 00000000 fc800000 00000000
*
*                2*(z + y) = Infinity
*
*                (z + y) + (z + y) = NaNQ
*
*              If the loop were allowed to proceed one more iteration,
*              it would incorrectly produce
*
*                  z = 8.988465674311579039686448570265e307
*                    = 7fe00000 00000000 fc800000 00000000
*
*              thereby losing the string of consecutive one bits, and
*              producing a pair with opposite signs, sigh...  It does
*              not seem possible to produce the desired
*
*                  z = 7fefffff ffffffff 7c7fffff ffffffff
*
*              by any obvious computation that does not require
*              knowledge of p = 53.  For example, we could
*              generate
*
*                  z+(2**(-54))*z = 1.797693134862315807937289714053e308
*                                 = 7fefffff ffffffff 7c8fffff ffffffff
*
*              but I do not wish to have any machine-dependent magic
*              constants here.
*
*              On other systems, the qisNaN() tests at entry will
*              properly terminate this loop.
*
               IF (.NOT.qisNaN((z + y) + (z + y))) GO TO 20
           END IF
*
           fpmax = two * z
           IF (qisInf(fpmax)) THEN
*
*               On SGI IRIX 6.5, at iteration 52 above, we have
*
*                   z = 8.98846567431157854072637118658521800e+0307
*                     = 7FDFFFFF FFFFFFFF 00000000 00000000
*
*               and at iteration 53, we have
*
*                   z = 8.98846567431157903968644857026517100e+0307
*                     = 7FE00000 00000000 FC800000 00000000
*
*               This is even more anomalous, since the expected result
*               at iteration 53 is to have the first double word set to
*               dfpmax():
*
*                   z = 1.79769313486231570814527423731704360e+0308
*                     = 7FEFFFFF FFFFFFFF 00000000 00000000
*
*               The result is that the loop exits with
*
*                   z = 8.988465674311577542806216419225312e+0307
*                     = 7fe00000 00000000 f9300000 00000000
*
*               and doubling that value produces Infinity.  We
*               therefore recompute fpmax by creating a value
*               w whose binary representation is 1.11111...,
*               and then multiplying that by z.  The result is
*               again, erroneously, Infinity:
*
*                    x = four / three
*                    call qstore(x)
*                    y = x - one
*                    call qstore(y)
*                    w = three*y
*                    call qstore(w)
*                    w = one + w
*                    call qstore(w)
*                    fpmax = w * z
*
*               There seems to be no alternative to introducing an
*               ugly magic constant, whose value has been determined
*               by numerical experiments, sigh...
*
                w = (two - two**(-52))
                fpmax = w*z
*
           END IF
*          WRITE (stdout,'(''DONE:'', 1p, e45.35e4, 2x, z32.32)')
*    X         fpmax,fpmax
*
      END IF
*
      qfpmax = fpmax
*
      END
