      DOUBLE PRECISION FUNCTION dlgam(x)
************************************************************************
*     (Double-precision ln(abs(Gamma(x))))
*     Return ln(abs(Gamma(x))), where x is any representable value.
*     Unlike gamma(x), which overflows for even modest x, the return
*     of this function is finite and representable for all x > 0.
*     [03-Aug-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dabs,        dlog
*
*     Built-in functions
*
      DOUBLE PRECISION    dabs,        dlog
*
*     External functions
*
      EXTERNAL            dgamma,      dnan,        isdnan
*
      DOUBLE PRECISION    dgamma,      dnan
*
      LOGICAL             isdnan
*
*     Parameter variables
*
      DOUBLE PRECISION    half
      PARAMETER           (half = 0.5d+00)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
      INCLUDE 'dl2pi2.inc'
*
      INTEGER             ncuts
      PARAMETER           (ncuts = 6)
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
*     Local variables
*
      DOUBLE PRECISION    cutoff(ncuts),            sum,         xsqinv
*
      INTEGER             i,           ncuse,       nterms(ncuts)
*
      INCLUDE 'dcgam.inc'
*
*     We can use the asymptotic series with nterms(i) terms when
*     x < cutoff(i).  However, to avoid bit loss, we cannot have any
*     cutoff value below 12.491.  The tables must be arranged
*     according to increasing x, so that the search loop can bail out
*     as soon as a suitable cutoff is found.  Since the tables are
*     very small, smart compilers will unroll the search loop, and a
*     really good one will not even store unused elements in c(*)!
*
      DATA cutoff /  12.492d+00, 13.0d+00, 19.0d+00, 35.0d+00,
     X              112.0d+00,   1415.0d+00 /
      DATA nterms /   7,           6,       5,        4,
     X                3,           2 /
*
      IF (isdnan(x)) THEN
          dlgam = dnan()
      ELSE IF (x .LE. zero) THEN
          dlgam = dnan()
      ELSE IF (x .LT. cutoff(1)) THEN
          dlgam = dlog(dgamma(x))
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
          dlgam = ((x - half)*dlog(x) + (ln2pi2 + sum/x)) - x
      END IF
*
      END
