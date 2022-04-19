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
      DOUBLE PRECISION    dlamch
      EXTERNAL            dlamch
      INTEGER             ifail,       k,           m
      DOUBLE PRECISION    z(1,1)
      DOUBLE PRECISION    w(MAXPTS)
      DOUBLE PRECISION    work(5*MAXPTS)
      INTEGER             iwork(5*MAXPTS)
*
      CALL dstevx('N', 'A', n, d, e(2), 0.0d0, 0.0d0, 0, 0,
     X      dlamch('s'), m, w, z, 1, work, iwork, ifail, ierr)
*
      DO 10 k = 1,m
          d(k) = w(k)
   10 CONTINUE
*
      IF (m .NE. n) ierr = 99999
*
      END
