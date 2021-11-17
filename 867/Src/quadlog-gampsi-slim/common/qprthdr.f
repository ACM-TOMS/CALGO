      SUBROUTINE qprthd(unit, header)
************************************************************************
*     (Print header)
*     Print a standard test output header.
*     (05-May-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            qeps,        qlog2
*
      REAL*16             qeps,        qlog2
*
      INTEGER             iqnint
*
*     Parameter variables
*
      REAL*16             ONE
      PARAMETER           (ONE = 1.0q+00)
*
*     Argument variables
*
      CHARACTER*(*)       header
*
      INTEGER             unit
*
*     Local variables
*
      REAL*16             ulp
*
      ulp = qeps(ONE)
      WRITE (unit, '(2X, ''Test integral: '', A)') header
      WRITE (unit, '(/, 2X, ''1 ulp = '', 1P, E10.2)') ulp
      WRITE (unit, '(/, 2X, ''Floating-point precision appears '', ' //
     X     ' ''to be (1 - log2(ulp) =)'', I4, '' bits'')')
     X     (1 - iqnint(qlog2(ulp)))
*
      END
