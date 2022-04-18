      REAL*16 FUNCTION  qiran(x,y)
************************************************************************
*     (Quadruple-precision pseudo-random integer in (x..y))
*     Return a pseudo-random integer value, represented in quadruple
*     precision, in the range (x..y), excluding endpoint y, where x >= y
*     (a relation that is NOT checked).
*
*     The underlying pseudo-random number generator is dran(), which
*     produces about 58 random bits.
*
*     The range of representable integers is 0 .. (2**p - 1), where p
*     is the number of bits in the significand of a
*     quadruple-precision number.
*
*     In IEEE 754 quadruple-precision arithmetic, p = 113,
*     corresponding to the range
*     0 .. 10384593717069655257060992658440192 (about 0 .. 1.04e+34).
*     (30-Jun-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            dran,        qstorf
*
      DOUBLE PRECISION    dran
*
      REAL*16             qnint,       qstorf
*
*     Parameter variables
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
*     Argument variables
*
      REAL*16             x,           y
*
*     Local variables
*
      DOUBLE PRECISION    qtmp
*
*     NB: The qstorf() function calls are absolutely CRITICAL here: they
*     foil the use of fused-multiply-add instructions (Hal/Fujitsu
*     SPARC64-GP, HP PA-RISC, IBM Power and PowerPC, IBM S/390 G5 IEEE
*     754, Intel IA-64, National Semiconductor NS32381, SGI MIPS IV,
*     ...), and longer registers (Intel x86, Motorola 68K), which can
*     unacceptably lead to results that differ between architectures.
*     Test procedures require identical sequences of pseudo-random
*     numbers on all platforms.
*
      qtmp = dran()
      qiran = qnint(qstorf(x + qstorf(qstorf(qtmp)*qstorf(y - x))))
*
      IF (qiran .LT. x) THEN
          qiran = x
      ELSE IF (qiran .GE. (y - one)) THEN
          qiran = y - one
      END IF
*
      END
