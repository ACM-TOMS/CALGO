      LOGICAL FUNCTION isdnan(x)
************************************************************************
*     (Is double-precision x a NaN?)
*     Return .TRUE. if x is a NaN, and .FALSE. otherwise.
*
*     This function should be implementable as a simple inline test
*     for inequality of x with itself:
*
*         isdnan = (x .ne. x)
*
*     in ALL compilers for ALL programming languages on ALL systems
*     with IEEE 754 arithmetic.
*
*     Unfortunately, some compilers, even without optimization,
*     incorrectly reduce this test to .FALSE.  This happens with all
*     optimization levels on SGI IRIX 6.x f77 and f90 compilers.
*     Thus, we have to obfuscate the test by wrapping one operand in
*     a function call.  This successfully foiled the SGI compilers,
*     without requiring disassembly and examination of the bit
*     patterns of x.
*     [10-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dstorf
*
      DOUBLE PRECISION    dstorf
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
      isdnan = (dstorf(x) .NE. x)
*
      END
**
**     This more complex version is needed to work around PGI compiler
**     bugs: It also works on Compaq/DEC, GNU/Linux, HP, IBM, SGI, and
**     Sun systems, but does need the MILSTD/Fortran 95 bit intrinsics.
**
*      LOGICAL FUNCTION isdnan(x)
*      INTEGER             iand,        ishft
*      DOUBLE PRECISION    x
*      DOUBLE PRECISION    xx
*      INTEGER             hi,          ix(2),       lo
*      EQUIVALENCE         (ix(1), xx)
*      xx = 1.0d+00
*      IF (ix(1) .EQ. 0) THEN
*          hi = 2
*          lo = 1
*      ELSE
*          hi = 1
*          lo = 2
*      END IF
*      xx = x
*      isdnan = ((iand(2047,ishft(ix(hi), -20)) .EQ. 2047) .AND.
*     X         ((iand(ix(hi),1048575) .NE. 0) .OR. (ix(lo) .NE. 0)))
*      END
