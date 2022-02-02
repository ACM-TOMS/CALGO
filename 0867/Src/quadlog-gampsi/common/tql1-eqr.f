      SUBROUTINE tql1(n,d,e,ierr)
*
*     Argument variables
*
      DOUBLE PRECISION    d(*),        e(*)
*
      INTEGER             ierr,        n
*
      DOUBLE PRECISION    z(1,1)
*
      CALL dsteqr('N', n, d, e(2), z, 1, z, ierr)
*
      END
