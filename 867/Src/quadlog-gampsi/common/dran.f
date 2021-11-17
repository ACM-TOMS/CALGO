      DOUBLE PRECISION FUNCTION dran()
************************************************************************
*     (Double-precision pseudo-random number)
*     Generate and return a double-precision pseudo-random number from
*     the interval (0.0, 1.0).  The significand of the returned value
*     should have about 58 pseudo-random bits.  This is sufficient for
*     most arithmetic systems: in IEEE 754 arithmetic, the fractional
*     part of the double-precision significand has only 52 bits.
*
*     The initial generator seed is the same on the first call to this
*     function after every program startup, so that the sequence of
*     pseudo-random number is reproducible.  This routine has no
*     provision for the user to alter the initial seed.
*
*     The algorithm is based on ``ACM Algorithm 266: Pseudo-Random
*     Numbers'', by M. C. Pike and I. D.  Hill, Communications of the
*     ACM, Vol. 8, No. 10, 605--606, October 1965, modified by
*     Hansson, and later used in the book ``Software Manual for the
*     Elementary Functions'', by W. J. Cody, Jr. and W. Waite,
*     Prentice-Hall (1980), ISBN 0-13-822064-6.
*     [14-Jul-2000]
************************************************************************
*
*     The single precision version of this subprogram is intended for
*     use on computers with fixed-point wordlength of at least 29
*     bits.  It is best if the floating point significand has at most
*     29 bits.
*
*     Following Cody and Waite's recommendation (p. 14), we produce a
*     pair of random numbers and use ran1 + 2**(-29)*ran2 in an
*     attempt to generate about 58 random bits.
*     [14-Jul-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dble
*
*     Built-in functions
*
      DOUBLE PRECISION    dble
*
*     External functions
*
      EXTERNAL            dfloat,      dstorf
*
      DOUBLE PRECISION    dfloat,      dstorf
*
*     Local variables
*
      INTEGER             iy
*
*     The variable iy must be saved to provide a memory across calls
*     to this function:
*
      SAVE                iy
*
      DATA                iy     /100001/
*
*     iy is always in the range (1..2796202) <= (1..2^22), so that it
*     can be represented exactly when converted to a single-precision
*     floating-point number in DEC VAX or IEEE 754 arithmetic.  That
*     conversion might lose 1 bit in IBM S/360 arithmetic, due to its
*     wobbling precision (because of hexadecimal normalization) of 21
*     to 24 fraction bits.
*
      iy = iy * 125
      iy = iy - (iy/2796203) * 2796203
*
*     NB: The dstorf() calls here are CRITICAL: they prevent compilers
*     from replacing the divide by a multiply with the reciprocal of the
*     constant, which results in a completely different sequence of
*     pseudo-random numbers, sigh...  At least Lahey/Fujitsu lf95 and
*     Portland Group, Inc. pgf77, pgf90, pghpf, have this misfeature.
*
*     In the last statement, the outermost dstorf() call also prevents
*     use of a floating-point multiply-add instruction, which has
*     different precision: several RISC architectures have such
*     instructions (Hal/Fujitsu SPARC64-GP, HP PA-RISC, IBM Power and
*     PowerPC, IBM S/390 G5 IEEE 754, Intel IA-64, National
*     Semiconductor NS32381, SGI MIPS IV, ...)
*
      dran = dfloat(iy) / dstorf(2796203.0d+00)
*
      iy = iy * 125
      iy = iy - (iy/2796203) * 2796203
      dran = dran + dstorf(dstorf(dfloat(iy) / dstorf(2796203.0d+00)) /
     X    dstorf(536870912.0d+00))
*     ---------- last card of dran ----------
      END
