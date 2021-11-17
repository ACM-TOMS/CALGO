      PROGRAM chkqpsi
************************************************************************
*     (Check psi())
*     Check the rational Pad{\'e} polynomial computation for bit loss.
*
*     Worst case detected in 2**16 (= 65,536) tests:
*
*     p(*): ratio = 1.000 ==> 21 bits lost at x = 1.426834106445312
*     n = 13 for a,b = -0.144660124871207 0.144660192675276
*
*       Count Bits lost
*           1 21
*           1 16
*           2 15
*           3 14
*           7 13
*          14 12
*          31 11
*          57 10
*         117  9
*         235  8
*         467  7
*         931  6
*        1849  5
*        3655  4
*        7154  3
*       13891  2
*       22781  1
*
*     [from awk '{print $6}' foo5.out | sort | uniq -c | sort -k 2nr]
*
*     There were 51196 bit-loss errors in 65536 tests.
*
*     [17-Jul-2000]
************************************************************************
*
*     Parameter variables
*
      INTEGER             npq
      PARAMETER           (npq = 14)
*
*     Local variables
*
      REAL*16             p(npq),      q(npq)
*
      DATA p /
     X      -0.6772663063963943146981268595814598872655q-12,
     X      -0.1904154711899938785500147210861454676340q-09,
     X      -0.1726192322043028966696316266272174162958q-07,
     X      -0.7333999829727342891789570524825086814001q-06,
     X      -0.00001710576674500151359323507634447548937917q+00,
     X      -0.0002371376370700527096496043978898572305626q+00,
     X      -0.002024948057319534411645576803950819640566q+00,
     X      -0.01064326319851109499497561483359956482945q+00,
     X      -0.03265025436904749915512141722032013854796q+00,
     X      -0.04719679609603306221330377971163682064420q+00,
     X      +0.01158988887302360942682627845186940086886q+00,
     X      +0.1313153310780367574374243137841517761616q+00,
     X      +0.1446601926752758648034366367256299249635q+00,
     X      +0.4318157621453326544170589138424553447123q-01/
*
      DATA q /
     X      +0.1971487489941648656045337002495236786607q-16,
     X      -0.9852098855954226555277489615001155249312q-13,
     X      -0.3664893595885658648806040064650999793322q-10,
     X      -0.4044855440099553317864295999011778607957q-08,
     X      -0.2042964222406585085540925078138880683266q-06,
     X      -0.5639905867723081880176591156954063084098q-05,
     X      -0.9336351590846479443215427621915763860448q-04,
     X      -0.9749440447927778183252890461834995868319q-03,
     X      -0.6578026344782533458110337514920679055906q-02,
     X      -0.2874659673327186360739751996881547675474q-01,
     X      -0.7973525176413495859494221340047492083962q-01,
     X      -0.1332331954635853538515193955365492687544q+00,
     X      -0.1197351104491078296112314625458178989835q+00,
     X      -0.4318157621453326544170572147905978566876q-1/
*
      CALL check(p,q,npq,1.0q+00, 2.0q+00)
*
      END


      SUBROUTINE check(p,q,npq,xmin,xmax)
************************************************************************
*     (Check)
*     Sum the rational polynomial in psi(x) and check for bit loss for
*     x in the interval [xmin,xmax].
*     [17-Jul-2000]
************************************************************************
*
*     Parameter variables
*
      INTEGER             nstep
      PARAMETER           (nstep = 2**16)
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
*     Argument variables
*
      INTEGER             npq
*
      REAL*16             p(npq),      q(npq),      xmax,        xmin
*
*     Local variables
*
      INTEGER             i,           n
*
      REAL*16             dx,          xden,        xnum,        z
*
      dx = nstep
      dx = (xmax - xmin)/dx
      DO 200 n = 0, nstep
          xnum = p(1)
          xden = q(1)
          z = n
          z = xmin + dx*z
          DO 100 i = 2, npq
              CALL chkrat('p(*)', xnum*z, p(i), z, i)
              xnum = xnum*z + p(i)
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
      EXTERNAL            qlog2
*
      INTEGER             iqint
*
      REAL*16             qlog2
*
*     Parameter variables
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
      INCLUDE 'stdio.inc'
*
*     Argument variables
*
      CHARACTER*(*)       name
*
      INTEGER             n
*
      REAL*16             a,           b,           x
*
*     Local variables
*
      INTEGER             lost
*
      REAL*16             ratio
*
      IF (b .ne. zero) THEN
          ratio = -a/b
          IF ((0.5q+00 .LT. ratio) .AND. (ratio .LT. 2.0q+00)) THEN
              IF (ratio .LT. 1.0q+00) THEN
                  lost = iqint(-qlog2(1.0q+00 - ratio))
              ELSE IF (ratio .GT. 1.0q+00) THEN
                  lost = iqint(1.0q+00 - qlog2(ratio - 1.0q+00))
              ELSE
                  lost = 999
              END IF
              WRITE (stdout,10000) name, ratio, lost, x, n, a, b
          END IF
      END IF
10000 FORMAT (a, ': ratio = ', f5.3, ' ==> ', i3, ' bits lost at x = ',
     X     f20.15, ' n = ', i2, ' for a,b = ', 2f25.15)
      END
