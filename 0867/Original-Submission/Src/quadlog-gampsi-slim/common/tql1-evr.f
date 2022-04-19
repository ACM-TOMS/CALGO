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
      INTEGER             k,           m
      DOUBLE PRECISION    z(1,1)
      DOUBLE PRECISION    w(MAXPTS)
      DOUBLE PRECISION    work(20*MAXPTS)
      INTEGER             isuppz(2*MAXPTS)
      INTEGER             iwork(10*MAXPTS)
*
      CALL dstevr('N', 'A', n, d, e(2), 0.0d0, 0.0d0, 0, 0,
     X      dlamch('s'), m, w, z, 1, isuppz, work, 20*MAXPTS, iwork,
     X      10*MAXPTS, ierr)
*
      DO 10 k = 1,m
          d(k) = w(k)
   10 CONTINUE
*
      IF (m .NE. n) ierr = 99999
*
      END
