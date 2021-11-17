      REAL FUNCTION  airan(x,y)
************************************************************************
*     (Single-precision pseudo-random integer in (x..y))
*     Return a pseudo-random integer value, represented in single
*     precision, in the range (x..y), excluding endpoint y, where x >= y
*     (a relation that is NOT checked).
*
*     The range of representable integers is 0 .. (2**p - 1), where p
*     is the number of bits in the significand of a single-precision
*     number.
*
*     In IEEE 754 single-precision arithmetic, p = 23, corresponding
*     to the range 0 .. 8388608 (about 0 .. 8.28e+6).
*     (30-Jul-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            astorf,      ran
*
      REAL                anint,       astorf,      ran
*
*     Parameter variables
*
      REAL                one
      PARAMETER           (one = 1.0)
*
*     Argument variables
*
      REAL                x,           y
*
*     NB: The astorf() function calls are absolutely CRITICAL here: they
*     foil the use of fused-multiply-add instructions (Hal/Fujitsu
*     SPARC64-GP, HP PA-RISC, IBM Power and PowerPC, IBM S/390 G5 IEEE
*     754, Intel IA-64, National Semiconductor NS32381, SGI MIPS IV,
*     ...), and longer registers (Intel x86, Motorola 68K), which can
*     unacceptably lead to results that differ between architectures.
*     Test procedures require identical sequences of pseudo-random
*     numbers on all platforms.
*
      airan = anint(astorf(x + astorf(astorf(ran())*astorf(y - x))))
*
      IF (airan .LT. x) THEN
          airan = x
      ELSE IF (airan .GE. astorf(y - one)) THEN
          airan = y - one
      END IF
*
      END
