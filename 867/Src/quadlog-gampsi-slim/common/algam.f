      REAL FUNCTION algam(x)
************************************************************************
*     (Single-precision ln(abs(Gamma(x))))
*     Return ln(abs(Gamma(x))), where x is any representable value.
*     Unlike gamma(x), which overflows for even modest x, the return
*     of this function is finite and representable for all x > 0.
*     [03-Aug-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           abs,         alog
*
*     Built-in functions
*
      REAL                abs,         alog
*
*     External functions
*
      EXTERNAL            anan,        gamma,       isanan
*
      LOGICAL             isanan
*
      REAL                anan,        gamma
*
*     Parameter variables
*
      INTEGER             ncuts
      PARAMETER           (ncuts = 2)
*
      REAL                half
      PARAMETER           (half = 0.5)
*
      REAL                one
      PARAMETER           (one = 1.0)
*
      REAL                zero
      PARAMETER           (zero = 0.0)
*
      INCLUDE 'al2pi2.inc'
*
*     Argument variables
*
      REAL                x
*
*     Local variables
*
      INTEGER             i,           ncuse,       nterms(ncuts)
*
      REAL                cutoff(ncuts),            sum,         xsqinv
*
      INCLUDE 'acgam.inc'
*
*     We can use the asymptotic series with nterms(i) terms when
*     x < cutoff(i).  However, to avoid bit loss, we cannot have any
*     cutoff value below 12.491.  The tables must be arranged
*     according to increasing x, so that the search loop can bail out
*     as soon as a suitable cutoff is found.  Since the tables are
*     very small, smart compilers will unroll the search loop, and a
*     really good one will not even store unused elements in c(*)!
*
      DATA cutoff / 12.492, 377.0 /
      DATA nterms /  2,       1 /
*
      IF (isanan(x)) THEN
          algam = anan()
      ELSE IF (x .LE. zero) THEN
          algam = anan()
      ELSE IF (x .LT. cutoff(1)) THEN
          algam = alog(gamma(x))
      ELSE
*
*         Use the asymptotic series for ln(Gamma(x)), choosing the
*         minimal number of terms from the nterms(*) table, and
*         searching backwards, since large arguments are common, and
*         will thus cause early loop exit.  We initialize ncuse to the
*         maximal number of terms, to avoid optimizing compiler
*         use-before-set complaints.
*
          ncuse = nterms(1)
          DO 100 i = ncuts, 1, -1
              ncuse = nterms(i)
              IF (x .GE. cutoff(i)) GO TO 200
  100     CONTINUE
*
*         Sum the four parts in order from smallest to largest,
*         positive parts first.
*
  200     xsqinv = one/(x * x)
          sum = c(ncuse)
          DO 300 i = (ncuse - 1), 1, -1
               sum = sum*xsqinv + c(i)
  300     CONTINUE
          algam = ((x - half)*alog(x) + (ln2pi2 + sum/x)) - x
      END IF
*
      END
