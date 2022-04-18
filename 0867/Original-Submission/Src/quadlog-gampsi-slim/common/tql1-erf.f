      SUBROUTINE tql1(n,d,e,ierr)
*
*     Argument variables
*
      DOUBLE PRECISION    d(*),        e(*)
*
      INTEGER             ierr,        n
*
      CALL dsterf (n, d, e(2), ierr)
*
      END
