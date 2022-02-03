      PROGRAM chkdgam
************************************************************************
*     (Check Gamma())
*     Check the rational Pad{\'e} polynomial computation for bit loss.
*
*     Worst cases detected in 2**16 (= 65,536) tests:
*
*     p(*): ratio = 1.000 ==> 999 bits lost at x = 1.000000000000000
*     n = 8 for a,b = -66456.143820240540663 66456.143820240540663
*
*     p(*): ratio = 1.000 ==> 19 bits lost at x = 1.851928710937500
*     n = 4 for a,b = -629.333403550504045 629.331155312818396
*
*       Count Bits lost
*           1 999
*           1 19
*           2 17
*           4 16
*           7 15
*          15 14
*          31 13
*          61 12
*         121 11
*         242 10
*         483  9
*         970  8
*        1935  7
*        3874  6
*        7749  5
*        8136  4
*       13563  3
*       26153  2
*       43306  1
*
*     [from awk '{print $6}' foo5.out | sort | uniq -c | sort -k 2nr]
*
*     There were 106654 bit-loss errors in 65536 tests.
*
*     [17-Jul-2000]
************************************************************************
*
*     Parameter variables
*
      INTEGER             npq
      PARAMETER           (npq = 8)
*
*     Local variables
*
      DOUBLE PRECISION    p(npq),      q(npq)
*
      DATA p/ -1.71618513886549492533811d+0,
     X    2.47656508055759199108314d+1, -3.79804256470945635097577d+2,
     X    6.29331155312818442661052d+2, 8.66966202790413211295064d+2,
     X    -3.14512729688483675254357d+4, -3.61444134186911729807069d+4,
     X    6.64561438202405440627855d+4/
      DATA q/ -3.08402300119738975254353d+1,
     X    3.15350626979604161529144d+2, -1.01515636749021914166146d+3,
     X    -3.10777167157231109440444d+3, 2.25381184209801510330112d+4,
     X    4.75584627752788110767815d+3, -1.34659959864969306392456d+5,
     X    -1.15132259675553483497211d+5/
*
      CALL check(p,q,npq,1.0d+00, 2.0d+00)
*
      END


      SUBROUTINE check(p,q,npq,xmin,xmax)
************************************************************************
*     (Check)
*     Sum the rational polynomial in Gamma(x) and check for bit loss for
*     x in the interval [xmin,xmax].
*     [17-Jul-2000]
************************************************************************
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
      INTEGER             nstep
      PARAMETER           (nstep = 2**16)
*
*     Argument variables
*
      DOUBLE PRECISION    p(npq),      q(npq),      xmax,        xmin
*
      INTEGER             npq
*
*     Local variables
*
      DOUBLE PRECISION    dx,          xden,        xnum,        z
*
      INTEGER             i,           n
*
      dx = nstep
      dx = (xmax - xmin)/dx
      DO 200 n = 0, nstep
          xnum = zero
          xden = one
          z = n
          z = xmin + dx*z
          DO 100 i = 1, npq
              CALL chkrat('p(*)', xnum, p(i), z, i)
              xnum = (xnum + p(i))*z
              CALL chkrat('q(*)', xden*z, q(i), z, i)
              xden = xden*z + q(i)
  100     CONTINUE
  200 CONTINUE
      END


      SUBROUTINE chkrat(name, a, b, x, n)
************************************************************************
*     (Check ratio)
*     Check whether forming a - b will suffer bit loss, and if so,
*     report how much.  name, x, and n identify the sum, coordinate, and
*     term where the loss occurs.
*     [17-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dlog2
*
      DOUBLE PRECISION    dlog2
*
      INTEGER             idint
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
      INCLUDE 'stdio.inc'
*
*     Argument variables
*
      CHARACTER*(*)       name
*
      DOUBLE PRECISION    a,           b,           x
*
      INTEGER             n
*
*     Local variables
*
      DOUBLE PRECISION    ratio
*
      INTEGER             lost
*
      IF (b .ne. zero) THEN
          ratio = -a/b
          IF ((0.5d+00 .LT. ratio) .AND. (ratio .LT. 2.0d+00)) THEN
              IF (ratio .LT. 1.0d+00) THEN
                  lost = idint(-dlog2(1.0d+00 - ratio))
              ELSE IF (ratio .GT. 1.0d+00) THEN
                  lost = idint(1.0d+00 - dlog2(ratio - 1.0d+00))
              ELSE
                  lost = 999
              END IF
              WRITE (stdout,10000) name, ratio, lost, x, n, a, b
          END IF
      END IF
10000 FORMAT (a, ': ratio = ', f5.3, ' ==> ', i3, ' bits lost at x = ',
     X     f20.15, ' n = ', i2, ' for a,b = ', 2f25.15)
      END
