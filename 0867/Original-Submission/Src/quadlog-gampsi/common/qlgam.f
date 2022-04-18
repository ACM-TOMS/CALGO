      REAL*16 FUNCTION qlgam(x)
************************************************************************
*     (Quadruple-precision ln(abs(Gamma(x))))
*     Return ln(abs(Gamma(x))), where x is any representable value.
*     Unlike gamma(x), which overflows for even modest x, the return
*     of this function is finite and representable for all x > 0.
*     [03-Aug-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           qabs,        qlog
*
*     Built-in functions
*
      REAL*16             qabs,        qlog
*
*     External functions
*
      EXTERNAL            isqnan,      qnan,       qgamma
*
      LOGICAL             isqnan
*
      REAL*16             qnan,        qgamma
*
*     Parameter variables
*
      INTEGER             ncuts
      PARAMETER           (ncuts = 12)
*
      REAL*16             half
      PARAMETER           (half = 0.5q+00)
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
      INCLUDE 'ql2pi2.inc'
*
*     Argument variables
*
      REAL*16             x
*
*     Local variables
*
      INTEGER             i,           ncuse,       nterms(ncuts)
*
      REAL*16             cutoff(ncuts),            sum,         xsqinv
*
      INCLUDE 'qcgam.inc'
*
*     We can use the asymptotic series with nterms(i) terms when
*     x < cutoff(i).  However, to avoid bit loss, we cannot have any
*     cutoff value below 12.491.  The tables must be arranged
*     according to increasing x, so that the search loop can bail out
*     as soon as a suitable cutoff is found.  Since the tables are
*     very small, smart compilers will unroll the search loop, and a
*     really good one will not even store unused elements in c(*)!
*
      DATA cutoff /  12.492q+00, 13.0q+00, 14.0q+00, 15.0q+00, 17.0q+00,
     X      19.0q+00, 23.0q+00, 31.0q+00, 47.0q+00, 95.0q+00, 336.0q+00,
     X      971.0q+00 /
      DATA nterms /   29, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 5 /
*
      IF (isqnan(x)) THEN
          qlgam = qnan()
      ELSE IF (x .LE. zero) THEN
          qlgam = qnan()
      ELSE IF (x .LT. cutoff(1)) THEN
          qlgam = qlog(qgamma(x))
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
          qlgam = ((x - half)*qlog(x) + (ln2pi2 + sum/x)) - x
      END IF
*
      END
