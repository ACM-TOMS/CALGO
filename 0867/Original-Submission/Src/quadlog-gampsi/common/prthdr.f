      SUBROUTINE prthdr(unit, header)
************************************************************************
*     (Print header)
*     Print a standard test output header.
*     (05-May-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            deps,        dlog2
*
      DOUBLE PRECISION    deps,        dlog2
*
      INTEGER             idnint
*
*     Parameter variables
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
*     Argument variables
*
      CHARACTER*(*)       header
*
      INTEGER             unit
*
*     Local variables
*
      DOUBLE PRECISION    ulp
*
      ulp = deps(one)
      WRITE (unit, '(2X, ''Test integral: '', A)') header
      WRITE (unit, '(/, 2X, ''1 ulp = '', 1P, E10.2)') ulp
      WRITE (unit, '(/, 2X, ''Floating-point precision appears '', ' //
     X     ' ''to be (1 - log2(ulp) =)'', I3, '' bits'')')
     X     (1 - idnint(dlog2(ulp)))
*
      END
