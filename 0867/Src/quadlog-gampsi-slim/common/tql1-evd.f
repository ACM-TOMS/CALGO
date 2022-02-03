      SUBROUTINE tql1(n,d,e,ierr)
*
*     Argument variables
*
      DOUBLE PRECISION    d(*),        e(*)
*
      INTEGER             ierr,        n
*
      INCLUDE 'maxpts.inc'
*
      DOUBLE PRECISION    z(1,1)
      DOUBLE PRECISION    work(MAXPTS)
      INTEGER             iwork(MAXPTS)
*
      CALL dstevd('N', n, d, e(2), z, 1, work, MAXPTS, iwork, MAXPTS, 
     X    ierr)
*
      END
