      PROGRAM demgjq
*     (Double-precision x^p integration test)
*     Evaluate the integral given at the end of paper 1,
*     doing all quadrature computations in double precision.
*
*        int_{-1}^{+1} (1-x)^\alpha (1+x)^\beta \ln ((1+x)/2) (1-x)^p dx
*
*     [01-Sep-2003]
*
*     External functions
*
      EXTERNAL            deps,        dgamma,      dpsi,        dvsum
      EXTERNAL            pad
*     EXTERNAL            fphand
*
      CHARACTER*30        pad
*
      DOUBLE PRECISION    dble,        deps,        dgamma,      dlog
      DOUBLE PRECISION    dpsi,        dvsum
*
*     Statement functions
*
      DOUBLE PRECISION    f,           fprime
*
      INCLUDE 'dlgtwo.inc'
*
*     Parameter variables
*
      DOUBLE PRECISION    FOUR
      PARAMETER           (FOUR = 4.0D+00)
      DOUBLE PRECISION    HALF
      PARAMETER           (HALF = 0.5D+00)
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0D+00)
      DOUBLE PRECISION    ONEP5
      PARAMETER           (ONEP5 = 1.5D+00)
      DOUBLE PRECISION    TWO
      PARAMETER           (TWO = 2.0D+00)
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0D+00)
*
      INCLUDE 'maxpts.inc'
*
*     Local variables
*
      CHARACTER*30        c1,          c2
*
      DOUBLE PRECISION    a(0:MAXPTS), alpha,       b(0:MAXPTS), beta
      DOUBLE PRECISION    deltaw(MAXPTS),           deltax(MAXPTS)
      DOUBLE PRECISION    fexact,      fp,          parg,        rle1
      DOUBLE PRECISION    rle2,        rle3,        rle4
      DOUBLE PRECISION    s(0:MAXPTS), sum1,        sum2,        sum3
      DOUBLE PRECISION    sum4,        t(0:MAXPTS), t1(3,MAXPTS)
      DOUBLE PRECISION    t2(MAXPTS),  t3(MAXPTS),  t4(MAXPTS),  ulp
      DOUBLE PRECISION    ulps1,       ulps2,       ulps3,       ulps4
      DOUBLE PRECISION    w(MAXPTS),   ww(MAXPTS),  x(MAXPTS),   xarg
      DOUBLE PRECISION    xx(MAXPTS),  yy(MAXPTS),  zz(MAXPTS)
*
      INTEGER             i,           ierr,        n,           nquad
      INTEGER             p,           pmax
*
      f(parg,xarg) = (ONE - xarg)**parg
      fprime(parg,xarg) = -parg * (ONE - xarg)**(parg - ONE)
*
*     CALL ieee_handler( 'set', 'common', fphand)
*     CALL ieee_handler( 'set', 'underflow', fphand)
*
      ulp = deps(ONE)
*
      READ (5, *, END = 1200) alpha, beta, nquad, pmax
      IF (nquad .GT. MAXPTS) STOP '(nquad .GT. MAXPTS)'
      IF (nquad .LT. 1) STOP '(nquad .LT. 1)'
      IF (alpha .LE. -ONE) STOP '(alpha .LE. -ONE)'
      IF (beta .LE. -ONE) STOP '(beta .LE. -ONE)'
      WRITE (6, 10000) alpha, beta, nquad, pmax
*
      CALL gjqrc(a, b, s, t, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'gjqrc() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
*
      WRITE (6, 90000)
      WRITE (6, 96000) (i, a(i), b(i), i = 0, nquad)
      WRITE (6, 95000)
      WRITE (6, 96000) (i, t(i), s(i), i = 0, nquad)
*
      CALL gjqf(xx, ww, yy, zz, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'gjqf() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
      WRITE (6, 11000)
      WRITE (6, 13000)
      WRITE (6, 14000) (i, xx(i), ww(i), i = 1, nquad)
      WRITE (6, 15000)
      WRITE (6, 14000) (i, yy(i), zz(i), i = 1, nquad)
*
      CALL gjqfd(x, w, deltaw, deltax, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'gjqfd() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
      WRITE (6, 12000)
      WRITE (6, 16000)
      WRITE (6, 14000) (i, x(i), w(i), i = 1, nquad)
      WRITE (6, 17000)
      WRITE (6, 14000) (i, deltaw(i), deltax(i), i = 1, nquad)
*
*     The nodes and weights for the nonlogarithmic quadature should be
*     identical from both gjqf() and gjqfd(): check this.
*
      ierr = 0
      DO 100 i = 1, nquad
          IF (x(i) .NE. xx(i)) ierr = ierr + 1
          IF (w(i) .NE. ww(i)) ierr = ierr + 1
  100 CONTINUE
      IF (ierr .NE. 0) WRITE (6, *) 'WARNING: unexpected difference',
     X    'in (x,w) from gjqf() and gjqfd()'
*
      WRITE (6, 70000)
      WRITE (6, 30000)
      DO 300 p = 0, pmax
           fp = dble(p)
           fexact = TWO**(alpha + beta + ONE + fp) *
     x         ((dgamma(ONE + beta) * dgamma(alpha + ONE + fp)) /
     x          dgamma(alpha + beta + TWO + fp)) *
     x         (dpsi(ONE + beta) - dpsi(alpha + beta + TWO + fp))
           DO 200 i = 1, nquad
                t1(1,i) = deltaw(i) * f(fp,x(i))
                t1(2,i) = deltax(i) * fprime(fp,x(i))
                t1(3,i) = -DLGTWO * w(i) * f(fp,x(i))
                t2(i) = w(i) * dlog((ONE + x(i))/TWO) * f(fp,x(i))
  200      CONTINUE
           sum1 = dvsum(t1, 3 * nquad)
           sum2 = dvsum(t2, nquad)
           rle1 = (fexact - sum1)/fexact
           ulps1 = rle1 / ulp
           rle2 = (fexact - sum2)/fexact
           ulps2 = rle2 / ulp
           WRITE (6, 20000) p, fexact, sum1, rle1, ulps1, sum2, rle2,
     X         ulps2
  300 CONTINUE

      WRITE (6, 70000)
      WRITE (6, 40000)
      DO 500 p = 0, pmax
           fp = dble(p)
           fexact = TWO**(alpha + beta + ONE + fp) *
     x         ((dgamma(ONE + beta) * dgamma(alpha + ONE + fp)) /
     x          dgamma(alpha + beta + TWO + fp)) *
     x         (dpsi(ONE + beta) - dpsi(alpha + beta + TWO + fp))
           DO 400 i = 1, nquad
                t3(i) = -zz(i)*f(fp,yy(i))
                t4(i) = ww(i)*dlog((ONE + xx(i))/TWO)*f(fp,xx(i))
  400      CONTINUE
           sum3 = dvsum(t3, nquad)
           sum4 = dvsum(t4, nquad)
           rle3 = (fexact - sum3)/fexact
           ulps3 = rle3 / ulp
           rle4 = (fexact - sum4)/fexact
           ulps4 = rle4 / ulp
           WRITE (6, 20000) p, fexact, sum3, rle3, ulps3, sum4, rle4,
     X         ulps4
  500 CONTINUE

      WRITE (6, 70000)
      WRITE (6, 50000)
      DO 700 p = 0, pmax
           fp = dble(p)
           fexact = TWO**(alpha + beta + ONE + fp) *
     x         ((dgamma(ONE + beta) * dgamma(alpha + ONE + fp)) /
     x          dgamma(alpha + beta + TWO + fp)) *
     x         (dpsi(ONE + beta) - dpsi(alpha + beta + TWO + fp))
           DO 600 i = 1, nquad
                t1(1,i) = deltaw(i) * f(fp,x(i))
                t1(2,i) = deltax(i) * fprime(fp,x(i))
                t1(3,i) = -DLGTWO * w(i) * f(fp,x(i))
                t2(i) = w(i)*dlog((ONE + x(i))/TWO)*f(fp,x(i))
                t3(i) = -zz(i)*f(fp,yy(i))
                t4(i) = ww(i)*dlog((ONE + xx(i))/TWO)*f(fp,xx(i))
  600      CONTINUE
           sum1 = dvsum(t1, 3 * nquad)
           sum2 = dvsum(t2, nquad)
           sum3 = dvsum(t3, nquad)
           sum4 = dvsum(t4, nquad)
           rle1 = (fexact - sum1)/fexact
           rle2 = (fexact - sum2)/fexact
           rle3 = (fexact - sum3)/fexact
           rle4 = (fexact - sum4)/fexact
           WRITE (6, 60000) p, rle1, rle2, rle3, rle4
  700 CONTINUE
*
      WRITE (6,80000)
      WRITE (6,85000)
      DO 900 p = 0, pmax
           fp = dble(p)
           fexact = TWO**(alpha + fp + beta + ONE) *
     X          dgamma(alpha + fp + ONE) * dgamma(beta + ONE) /
     X          dgamma(alpha + fp + beta + TWO)
           DO 800 i = 1, nquad
                t2(i) = w(i) * f(fp,x(i))
                t4(i) = ww(i) * f(fp,x(i))
  800      CONTINUE
           sum2 = dvsum(t2, nquad)
           sum4 = dvsum(t4, nquad)
           rle2 = (fexact - sum2)/fexact
           ulps2 = rle2 / ulp
           rle4 = (fexact - sum4)/fexact
           ulps4 = rle4 / ulp
           WRITE (6, 20000) p, fexact, sum2, rle2, ulps2,
     X          sum4, rle4, ulps4
  900 CONTINUE
*
      IF ((alpha .EQ. ZERO) .AND.
     X    ( (BETA .EQ. -HALF) .OR.
     X      (BETA .EQ. ZERO) .OR.
     X      (BETA .EQ. +HALF) ) ) THEN
           WRITE (6, 97000)
           WRITE (6, 97500)
           WRITE (6, 98000)
           WRITE (6, 99000) 0, (ONE + b(0))/TWO, a(0)/TWO**(beta + ONE)
           WRITE (6, 99000) (i, (ONE + b(i))/TWO, a(i)/FOUR,
     X         i = 1, nquad)
      END IF
*
      IF ((alpha .EQ. beta) .AND.
     X     ((alpha .EQ. -HALF) .OR. (alpha .EQ. HALF) .OR.
     X      (alpha .EQ. ONE) .OR. (alpha .EQ. ONEP5))) THEN
          WRITE (6, 99100)
          DO 1100 n = 2, 20
              CALL gjqfd(x, w, deltaw, deltax, alpha, beta, n,ierr)
              IF (ierr .NE. 0) WRITE (6, *)
     X            'gjqfd() returns ierr =', ierr
              IF (ierr .NE. 0) STOP
              WRITE (6,99500) alpha
              WRITE (6,99200) n
              WRITE (6,99300)
              DO 1000 i = 1, n
                  c1 = pad(x(i))
                  c2 = pad(w(i))
                  WRITE (6,99400) i, c1, c2
 1000         CONTINUE
 1100     CONTINUE
      END IF
*
 1200 CONTINUE
*
10000 FORMAT (2x, 'alpha =', f10.6, '  beta = ', f10.6, '  nquad =', i5,
     x     '  pmax =', i4)
11000 FORMAT (//, 2x, 'Nodes and weights from gjqf()')
12000 FORMAT (//, 2x, 'Nodes and weights from gjqfd()')
13000 FORMAT (2x, 4x, 'n', 1x, 11x, 'x(n)', 10x, 1x, 11x, 'w(n)')
14000 FORMAT (2x, i5, 1x, 1pe25.16, 1x, 1pe25.16)
15000 FORMAT (/, 2x, 4x, 'n', 2x, 10x, 'y(n)', 11x, 2x, 9x, 'z(n)')
16000 FORMAT (2x, 4x, 'n', 1x, 11x, 'x(n)', 10x, 1x, 11x, 'w(n)')
17000 FORMAT (/, 2x, 4x, 'n', 1x, 10x, 'deltaw(n)', 6x,
     X    1x, 10x, 'deltax(n)')
20000 FORMAT (2x, i5, 1pe25.16, 2(1pe25.16, 1pe12.2, 1x, 0pf9.2))
30000 FORMAT (//, 2x, 4x, 'p', 7x, 'Exact Integral', 16x, 'gjqfd log',
     X    6x, 'RelErr', 10x, 'Ulps', 10x, 'gjqfd nolog', 6x, 'RelErr',
     X    10x, 'Ulps')
40000 FORMAT (//,2x, 4x, 'p', 7x, 'Exact Integral', 17x, 'gjqf log',
     X    6x, 'RelErr', 10x, 'Ulps', 11x, 'gjqf nolog', 6x, 'RelErr',
     X    10x, 'Ulps')
50000 FORMAT (//, 2x, 'Comparative relative errors', /,
     X    2x, 5x, '--------------------- gjqfd ----------------------',
     X    1x, ' -------------------- gjqf ----------------------', /,
     X    2x, 22x, '(1-x)**p',
     X    8x,  'log(x) * (1-x)**p',
     X    17x, '(1-x)**p',
     X    8x,  'log(x) * (1-x)**p', /,
     X    2x, 4x, 'p', 13x, 'log-Jacobi', 9x, 'Gauss-log-Jacobi',
     X    13x, 'Gauss-Jacobi', 13x, 'Gauss-Jacobi')
60000 FORMAT (2x, i5, 1p4e25.2)
70000 FORMAT (//, 2x, '\\int_{-1}^{+1} (1-x)^\\alpha (1+x)^\\beta',
     X    ' \\ln ((1+x)/2) (1-x)^p dx')
80000 FORMAT (//, 2x, '\\int_{-1}^{+1} (1-x)^\\alpha (1+x)^\\beta',
     X    ' (1-x)^p dx')
85000 FORMAT (//,2x, 4x, 'p', 7x, 'Exact Integral', 19x, 'gjqfd',
     X    7x, 'RelErr', 10x, 'Ulps', 14x, 'gjqf', 9x, 'RelErr',
     X    10x, 'Ulps')
90000 FORMAT (//, 2x, 'Monic polynomial recursion coefficients'/
     X       1x, 3x, 'n', 2x, 10x, 'a(n)', 11x, 2x, 10x, 'b(n)')
95000 FORMAT (//, 2x, 'Monic polynomial zeroth and first moments'/
     X       1x, 3x, 'n', 2x, 10x, 't(n)', 11x, 2x, 10x, 's(n)')
96000 FORMAT (1x, i4, 2x, 1pe25.16, 2x, e25.16)
97000 FORMAT (//,
     X       2x, 'For beta = sigma = -1/2, 0, +1/2, see the subset in ',
     X           'Table III, p. 33 in:', /,
     X      10x, 'Walter Gautschi', /,
     X      10x, 'Algorithm 726: ORTHPOL --- A Package of Routines ',
     X           'for Generating', /,
     X      10x, 'Orthogonal Polynomials and Gauss-Type Quadrature ',
     X           'Rules', /,
     X      10x, 'ACM Transactions on Mathematical Software 20(1) ',
     X           '21--62, March 1994', //,
     X       2x, 'For beta = sigma = -1/2, see the complete list in ',
     X           'Table 1 in:', /,
     X      10x, 'Walter Gautschi', /,
     X      10x, 'A class of slowly convergent series and their ',
     X           'summation by Gaussian quadrature', /,
     X      10x, 'Mathematics of Computation 57(195) 309--324, ',
     X           'July 1991')
97500 FORMAT (/,
     X       2x, 'Relationship of recursion coefficients:', /,
     X      10x, 'Gautschi alpha(k) = (1 + b(k))/2         this work, ',
     X           'k = 0, 1, 2, ...', /,
     X      10x, 'Gautschi  beta(0) = a(0)/2**(sigma + 1)  this work ',
     X           '(Gautschi sigma = our beta)', /,
     X      10x, 'Gautschi  beta(k) = a(k)/4               this work, ',
     X           'k = 1, 2, 3, ...')
98000 FORMAT (//, 5x, 'k', 2x, 8x, 'alpha(k)', 8x, 2x, 8x, 'beta(k)')
99000 FORMAT (2x, i4, 2x, 0pe24.17, 2x, 0pe24.17)
99100 FORMAT(/,
     X       2x, 'Nodes and weights for ordinary Gauss-Jacobi ',
     X           'quadrature: see', /,
     X       2x, 'Table 2 (pp. 157--169) in the book', //,
     X      10x, 'A. H. Stroud and Don Secrest', /,
     X      10x, 'Gaussian Quadrature Formulas', /,
     X      10x, 'Prentice-Hall 1966', //,
     X       2x, 'For alpha = beta, the weights and nodes are ',
     X           'symmetric about zero,', /,
     X       2x, 'so the book tabulates only half of them.  In this ',
     X           'program, symmetry', /,
     X       2x, 'is not enforced by the underlying method, so minor ',
     X           'deviations from', /,
     X       2x, 'strict symmetry are expected.')
99200 FORMAT(2x, 4x, 30x, 'N = ', i2)
99300 FORMAT(2x, 4x, 2x, 13x, 'x(i)', 13x, 2x, 13x, 'A(i)')
99400 FORMAT(2x, i4, 2x, a, 2x, a)
99500 FORMAT(///, 2x, 4x, 26x, 'alpha = ', F4.1)
      END

      CHARACTER*30 FUNCTION pad(x)
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
*     Local variables
*
      CHARACTER*24        s
*
      WRITE (s, '(0pe24.15e4)') x
*
      pad = s(1:8) // '_' // s(9:13) // '_' // s(14:18) // '_' //
     X    s(19:)
*
      END
