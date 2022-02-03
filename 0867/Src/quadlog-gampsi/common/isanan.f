      LOGICAL FUNCTION isanan(x)
************************************************************************
*     (Is single-precision x a NaN?)
*     Return .TRUE. if x is a NaN, and .FALSE. otherwise.
*
*     This function should be implementable as a simple inline test
*     for inequality of x with itself:
*
*         isanan = (x .ne. x)
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
      EXTERNAL            astorf
*
      REAL                astorf
*
*     Argument variables
*
      REAL                x
*
      isanan = (astorf(x) .NE. x)
*
      END
