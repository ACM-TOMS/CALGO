      DOUBLE PRECISION FUNCTION  diran(x,y)
************************************************************************
*     (Double-precision pseudo-random integer in (x..y))
*     Return a pseudo-random integer value, represented in double
*     precision, in the range (x..y), excluding endpoint y, where x >= y
*     (a relation that is NOT checked).
*
*     The range of representable integers is 0 .. (2**p - 1), where p is
*     the number of bits in the significand of a double-precision
*     number.
*
*     In IEEE 754 double-precision arithmetic, p = 53, corresponding to
*     the range 0 .. 9007199254740991.
*     (30-Jul-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            dran,        dstorf
*
      DOUBLE PRECISION    dnint,       dran,        dstorf
*
*     Parameter variables
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
*     Argument variables
*
      DOUBLE PRECISION    x,           y
*
*     NB: The dstorf() function calls are absolutely CRITICAL here: they
*     foil the use of fused-multiply-add instructions (Hal/Fujitsu
*     SPARC64-GP, HP PA-RISC, IBM Power and PowerPC, IBM S/390 G5 IEEE
*     754, Intel IA-64, National Semiconductor NS32381, SGI MIPS IV,
*     ...), and longer registers (Intel x86, Motorola 68K), which can
*     unacceptably lead to results that differ between architectures.
*     Test procedures require identical sequences of pseudo-random
*     numbers on all platforms.
*
      diran = dnint(dstorf(x + dstorf(dstorf(dran())*dstorf(y - x))))
*
      IF (diran .LT. x) THEN
          diran = x
      ELSE IF (diran .GE. dstorf(y - one)) THEN
          diran = y - one
      END IF
*
      END
