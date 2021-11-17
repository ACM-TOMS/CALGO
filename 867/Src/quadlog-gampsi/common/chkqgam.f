      PROGRAM chkqgam
************************************************************************
*     (Check Gamma())
*     Check the rational Pad{\'e} polynomial computation for bit loss.
*
*     Worst case detected in 2**16 (= 65,536) tests:
*
*         q(*): ratio = 1.083 ==> 4 bits lost at x = 1.000000000000000
*         n = 15 for a,b = 0.075864423806676 -0.070061270864165
*
*       Count Bits lost
*        2806 4
*        8553 3
*       26151 2
*       68542 1
*
*     [from awk '{print $6}' foo5.out | sort | uniq -c | sort -k 2nr]
*
*     There were 106052 bit-loss errors in 65536 tests.
*
*     [17-Jul-2000]
************************************************************************
*
*     Parameter variables
*
      INTEGER             npq
      PARAMETER           (npq = 17)
*
*     Local variables
*
      REAL*16             p(npq),      q(npq)
*
      DATA p /
     X      -5.77431901115732472295911874064153596940068136763167q-15,
     X      -1.56797490320050003213800346859245679901144331350448q-13,
     X      -3.78028631870137741397074534983799321203124415752340q-12,
     X      -6.45223850299882992346318581519945955137073079605824q-11,
     X      -9.57773260999682906148204876875358386860072064396121q-10,
     X      -1.18789946304711998961669755697047234308459820137236q-08,
     X      -1.29688059161157116936429224516158782474583552996390q-07,
     X      -1.22958541797707216613371499028802040668869394697011q-06,
     X      -1.02878289023582748093908396579316255215815995700520q-05,
     X      -7.50768338661706013484760924568582708783543038072889q-05,
     X      -4.78432978090750714710102993275260246273562014185706q-04,
     X      -2.61447025355310872578328515594831816304418138338835q-03,
     X      -1.21088373976211806683231520253652259470388300978605q-02,
     X      -4.58104373823350172120266058384144773790769448288121q-02,
     X      -1.36282200431838962638629833438264619264021627817802q-01,
     X      -2.85254037936504666809863185761056270762945836369187q-01,
     X      -3.52007543423276049518405042481564070179437479536879q-01 /
      DATA q /
     X      +6.86512547163641702960320902587557214739248632008009q-15,
     X      -7.59042505734871930619643757876547538111341894071448q-13,
     X      +3.56369935591054923181426610402720794774793212263523q-11,
     X      -9.43188492095640351711627001164301376867057746330490q-10,
     X      +1.52962101119166702455321211675458807858606912755968q-08,
     X      -1.48011763725540222994501936734900379008433539545861q-07,
     X      +6.32155894965970412510362786824597689209436354389813q-07,
     X      +2.83278123740864061993110584129497191676829717666322q-06,
     X      -5.12887804591908758756833611374530299268933223209224q-05,
     X      +1.97318185287111363511605005713815989285041401943088q-04,
     X      +6.06366157082808864433608154469176011937606857946674q-04,
     X      -6.48549852121071540402006206401631855319227236825582q-03,
     X      +4.19174110072042154435320268459010973870885020092510q-03,
     X      +7.74024543519734033480513779208663722839440995654525q-02,
     X      -7.00612708641652600486431341078527239938813932883259q-02,
     X      -4.88438306163926152366209271007899649217322649799316q-01,
     X      -3.52007543423276049518369071878409146334353800751657q-01 /
*
      CALL check(p,q,npq,1.0q+00, 2.0q+00)
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
          DO 100 i = 1, npq
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
