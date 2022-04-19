      PROGRAM dmqgjq
*     (Double-precision x^p integration test)
*     Evaluate the integral given at the end of paper 1,
*     doing all quadrature computations in double precision.
*
*        int_{-1}^{+1} (1-x)^\alpha (1+x)^\beta \ln ((1+x)/2) (1-x)^p dx
*
*     When alpha = 0, we report recursion coefficients which can be
*     compared with values published by Gautschi in 1991 and 1994.
*     [01-Sep-2003]
*
*     External functions
*
*     EXTERNAL            fphand
      EXTERNAL            pad,         qeps,        qgamma,      qpsi
      EXTERNAL            qvsum
*
      CHARACTER*45        pad
*
      REAL*16             qeps,        qext,        qgamma,      qlog
      REAL*16             qpsi,        qvsum
*
*     Statement functions
*
      REAL*16             f,           fprime
*
      INCLUDE 'qlgtwo.inc'
*
*     Parameter variables
*
      REAL*16             FOUR
      PARAMETER           (FOUR = 4.0Q+00)
      REAL*16             HALF
      PARAMETER           (HALF = 0.5Q+00)
      REAL*16             ONE
      PARAMETER           (ONE = 1.0Q+00)
      REAL*16             ONEP5
      PARAMETER           (ONEP5 = 1.5Q+00)
      REAL*16             TWO
      PARAMETER           (TWO = 2.0Q+00)
      REAL*16             ZERO
      PARAMETER           (ZERO = 0.0Q+00)
*
      INCLUDE 'maxpts.inc'
*
*     Local variables
*
      CHARACTER*45        c1,          c2
*
      INTEGER             i,           ierr,        n,           nquad
      INTEGER             p,           pmax
*
      REAL*16             a(0:MAXPTS), alpha,       b(0:MAXPTS), beta
      REAL*16             deltaw(MAXPTS),           deltax(MAXPTS)
      REAL*16             fexact,      fp,          parg,        rle1
      REAL*16             rle2,        rle3,        rle4
      REAL*16             s(0:MAXPTS), sum1,        sum2,        sum3
      REAL*16             sum4,        t(0:MAXPTS), t1(3,MAXPTS)
      REAL*16             t2(MAXPTS),  t3(MAXPTS),  t4(MAXPTS),  ulp
      REAL*16             ulps1,       ulps2,       ulps3,       ulps4
      REAL*16             w(MAXPTS),   ww(MAXPTS),  x(MAXPTS),   xarg
      REAL*16             xx(MAXPTS),  yy(MAXPTS),  zz(MAXPTS)
*
      f(parg,xarg) = (ONE - xarg)**parg
      fprime(parg,xarg) = -parg * (ONE - xarg)**(parg - ONE)
*
*     CALL ieee_handler( 'set', 'common', fphand)
*     CALL ieee_handler( 'set', 'underflow', fphand)
*
      ulp = qeps(ONE)
*
      READ (5, *, END = 1200) alpha, beta, nquad, pmax
      IF (nquad .GT. MAXPTS) STOP '(nquad .GT. MAXPTS)'
      IF (nquad .LT. 1) STOP '(nquad .LT. 1)'
      IF (alpha .LE. -ONE) STOP '(alpha .LE. -ONE)'
      IF (beta .LE. -ONE) STOP '(beta .LE. -ONE)'
      WRITE (6, 10000) alpha, beta, nquad, pmax
*
      CALL qgjqrc(a, b, s, t, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'qgjqrc() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
*
      WRITE (6, 90000)
      WRITE (6, 96000) (i, a(i), b(i), i = 0, nquad)
      WRITE (6, 95000)
      WRITE (6, 96000) (i, t(i), s(i), i = 0, nquad)
*
      CALL qgjqf(xx, ww, yy, zz, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'qgjqf() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
      WRITE (6, 11000)
      WRITE (6, 13000)
      WRITE (6, 14000) (i, xx(i), ww(i), i = 1, nquad)
      WRITE (6, 15000)
      WRITE (6, 14000) (i, yy(i), zz(i), i = 1, nquad)
*
      CALL qgjqfd(x, w, deltaw, deltax, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'qgjqfd() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
      WRITE (6, 12000)
      WRITE (6, 16000)
      WRITE (6, 14000) (i, x(i), w(i), i = 1, nquad)
      WRITE (6, 17000)
      WRITE (6, 14000) (i, deltaw(i), deltax(i), i = 1, nquad)
*
*     The nodes and weights for the nonlogarithmic quadature should be
*     identical from both qgjqf() and qgjqfd(): check this.
*
      ierr = 0
      DO 100 i = 1, nquad
          IF (x(i) .NE. xx(i)) ierr = ierr + 1
          IF (w(i) .NE. ww(i)) ierr = ierr + 1
  100 CONTINUE
      IF (ierr .NE. 0) WRITE (6, *) 'WARNING: unexpected difference',
     X    'in (x,w) from qgjqf() and qgjqfd()'
*
      WRITE (6, 70000)
      WRITE (6, 30000)
      DO 300 p = 0, pmax
           fp = qext(p)
           fexact = TWO**(alpha + beta + ONE + fp) *
     X         ((qgamma(ONE + beta) * qgamma(alpha + ONE + fp)) /
     X          qgamma(alpha + beta + TWO + fp)) *
     X         (qpsi(ONE + beta) - qpsi(alpha + beta + TWO + fp))
           DO 200 i = 1, nquad
                t1(1,i) = deltaw(i) * f(fp,x(i))
                t1(2,i) = deltax(i) * fprime(fp,x(i))
                t1(3,i) = -QLGTWO * w(i) * f(fp,x(i))
                t2(i) = w(i) * qlog((ONE + x(i))/TWO) * f(fp,x(i))
  200      CONTINUE
           sum1 = qvsum(t1, 3 * nquad)
           sum2 = qvsum(t2, nquad)
           rle1 = (fexact - sum1)/fexact
           ulps1 = rle1 / ulp
           rle2 = (fexact - sum2)/fexact
           ulps2 = rle2 / ulp
           WRITE (6, 20000) p, fexact, sum1, rle1, ulps1, sum2, rle2,
     X         ulps2
  300 CONTINUE
*
      WRITE (6, 70000)
      WRITE (6, 40000)
      DO 500 p = 0, pmax
           fp = qext(p)
           fexact = TWO**(alpha + beta + ONE + fp) *
     X         ((qgamma(ONE + beta) * qgamma(alpha + ONE + fp)) /
     X          qgamma(alpha + beta + TWO + fp)) *
     X         (qpsi(ONE + beta) - qpsi(alpha + beta + TWO + fp))
           DO 400 i = 1, nquad
                t3(i) = -zz(i)*f(fp,yy(i))
                t4(i) = ww(i)*qlog((ONE + xx(i))/TWO)*f(fp,xx(i))
  400      CONTINUE
           sum3 = qvsum(t3, nquad)
           sum4 = qvsum(t4, nquad)
           rle3 = (fexact - sum3)/fexact
           ulps3 = rle3 / ulp
           rle4 = (fexact - sum4)/fexact
           ulps4 = rle4 / ulp
           WRITE (6, 20000) p, fexact, sum3, rle3, ulps3, sum4, rle4,
     X         ulps4
  500 CONTINUE
*
      WRITE (6, 70000)
      WRITE (6, 50000)
      DO 700 p = 0, pmax
           fp = qext(p)
           fexact = TWO**(alpha + beta + ONE + fp) *
     X         ((qgamma(ONE + beta) * qgamma(alpha + ONE + fp)) /
     X          qgamma(alpha + beta + TWO + fp)) *
     X         (qpsi(ONE + beta) - qpsi(alpha + beta + TWO + fp))
           DO 600 i = 1, nquad
                t1(1,i) = deltaw(i) * f(fp,x(i))
                t1(2,i) = deltax(i) * fprime(fp,x(i))
                t1(3,i) = -QLGTWO * w(i) * f(fp,x(i))
                t2(i) = w(i)*qlog((ONE + x(i))/TWO)*f(fp,x(i))
                t3(i) = -zz(i)*f(fp,yy(i))
                t4(i) = ww(i)*qlog((ONE + xx(i))/TWO)*f(fp,xx(i))
  600      CONTINUE
           sum1 = qvsum(t1, 3 * nquad)
           sum2 = qvsum(t2, nquad)
           sum3 = qvsum(t3, nquad)
           sum4 = qvsum(t4, nquad)
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
           fp = qext(p)
           fexact = TWO**(alpha + fp + beta + ONE) *
     X          qgamma(alpha + fp + ONE) * qgamma(beta + ONE) /
     X          qgamma(alpha + fp + beta + TWO)
           DO 800 i = 1, nquad
                t2(i) = w(i) * f(fp,x(i))
                t4(i) = ww(i) * f(fp,x(i))
  800      CONTINUE
           sum2 = qvsum(t2, nquad)
           sum4 = qvsum(t4, nquad)
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
              CALL qgjqfd(x, w, deltaw, deltax, alpha, beta, n,ierr)
              IF (ierr .NE. 0) WRITE (6, *)
     X            'qgjqfd() returns ierr =', ierr
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
     X     '  pmax =', i4)
11000 FORMAT (//, 2x, 'Nodes and weights from qgjqf()')
12000 FORMAT (//, 2x, 'Nodes and weights from qgjqfd()')
13000 FORMAT (2x, 4x, 'n', 1x, 20x, 'x(n)', 19x, 1x, 20x, 'w(n)')
14000 FORMAT (2x, i5, 1x, 1pe43.34, 1x, 1pe43.34)
15000 FORMAT (/, 2x, 4x, 'n', 2x, 19x, 'y(n)', 20x, 2x, 18x, 'z(n)')
16000 FORMAT (2x, 4x, 'n', 1x, 20x, 'x(n)', 19x, 1x, 20x, 'w(n)')
17000 FORMAT (/, 2x, 4x, 'n', 1x, 19x, 'deltaw(n)', 15x,
     X    1x, 19x, 'deltax(n)')
20000 FORMAT (2x, i5, 1pe43.34, 2(1pe43.34, 1x, 1pe11.2, 1x, 0pf9.2))
30000 FORMAT (//, 2x, 4x, 'p', 16x, 'Exact Integral', 33x, 'qgjqfd log',
     X    15x, 'RelErr', 10x, 'Ulps', 18x, 'qgjqfd nolog',
     X    15x, 'RelErr', 10x, 'Ulps')
40000 FORMAT (//,2x, 4x, 'p', 16x, 'Exact Integral', 34x, 'qgjqf log',
     X    15x, 'RelErr', 10x, 'Ulps', 20x, 'qgjqf nolog', 14x, 'RelErr',
     X    10x, 'Ulps')
50000 FORMAT (//, 2x, 'Comparative relative errors', /,
     X    2x, 5x, '--------------------- qgjqfd ----------------------',
     X    1x, ' -------------------- qgjqf ----------------------', /,
     X    2x, 22x, '(1-x)**p',
     X    8x,  'log(x) * (1-x)**p',
     X    17x, '(1-x)**p',
     X    8x,  'log(x) * (1-x)**p', /,
     X    2x, 4x, 'p', 13x, 'log-Jacobi', 9x, 'Gauss-log-Jacobi',
     X    13x, 'Gauss-Jacobi', 13x, 'Gauss-Jacobi')
60000 FORMAT (2x, i5, 4(1x,1pe24.2))
70000 FORMAT (//, 2x, '\\int_{-1}^{+1} (1-x)^\\alpha (1+x)^\\beta',
     X    ' \\ln ((1+x)/2) (1-x)^p dx')
80000 FORMAT (//, 2x, '\\int_{-1}^{+1} (1-x)^\\alpha (1+x)^\\beta',
     X    ' (1-x)^p dx')
85000 FORMAT (//,2x, 4x, 'p', 16x, 'Exact Integral', 36x, 'qgjqfd',
     X    16x, 'RelErr', 10x, 'Ulps', 23x, 'qgjqf', 17x, 'RelErr',
     X    10x, 'Ulps')
90000 FORMAT (//, 2x, 'Monic polynomial recursion coefficients'/
     X       1x, 3x, 'n', 2x, 20x, 'a(n)', 21x, 2x, 20x, 'b(n)')
95000 FORMAT (//, 2x, 'Monic polynomial zeroth and first moments'/
     X       1x, 3x, 'n', 2x, 20x, 't(n)', 21x, 2x, 20x, 's(n)')
96000 FORMAT (1x, i4, 2x, 1pe45.34, 2x, e45.34)
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
98000 FORMAT (//, 5x, 'k', 2x, 12x, 'alpha(k)', 12x, 2x, 12x, 'beta(k)')
99000 FORMAT (2x, i4, 2x, 0pe32.25, 2x, 0pe32.25)
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
99200 FORMAT(2x, 4x, 45x, 'N = ', i2)
99300 FORMAT(2x, 4x, 2x, 20x, 'x(i)', 21x, 2x, 20x, 'A(i)')
99400 FORMAT(2x, i4, 2x, a, 2x, a)
99500 FORMAT(///, 2x, 4x, 41x, 'alpha = ', F4.1)
      END

      CHARACTER*45 FUNCTION pad(x)
*
*     Argument variables
*
      REAL*16             x
*
*     Local variables
*
      CHARACTER*39        s
*
      WRITE (s, '(0pe39.30e4)') x
*
      pad = s(1:8) // '_' // s(9:13) // '_' // s(14:18) // '_' //
     X    s(19:23) // '_' // s(24:28) // '_' // s(29:33) // '_' //
     X    s(34:)
*
      END
