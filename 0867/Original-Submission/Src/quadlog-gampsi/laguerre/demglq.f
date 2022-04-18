      PROGRAM demglq
*     (Double-precision x^p integration test)
*     Evaluate the difficult integral given at the end of paper 1,
*     doing all quadrature computations in double precision.
*
*         int_0^\infty x^\alpha \exp(-x) \ln x x^p dx
*
*     [01-Sep-2003]
*
*     External functions
*
      EXTERNAL            deps,        dgamma,      dlog2,       dpsi
      EXTERNAL            dvsum,       pad
*     EXTERNAL            fphand
*
      CHARACTER*30        pad
*
      DOUBLE PRECISION    dble,        deps,        dgamma,      dlog
      DOUBLE PRECISION    dlog2,       dpsi,        dvsum
*
*     Statement functions
*
      DOUBLE PRECISION    f,           fprime
*
*     Parameter variables
*
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0D+00)
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0D+00)
*
      INCLUDE 'maxpts.inc'
*
*     Local variables
*
      CHARACTER*30        c1,          c2
*
      DOUBLE PRECISION    a(0:MAXPTS), alpha,       b(0:MAXPTS)
      DOUBLE PRECISION    deltaw(MAXPTS),           deltax(MAXPTS)
      DOUBLE PRECISION    fexact,      fp,          parg,        rle1
      DOUBLE PRECISION    rle2,        rle3,        rle4
      DOUBLE PRECISION    s(0:MAXPTS), sum1,        sum2,        sum3
      DOUBLE PRECISION    sum4,        t(0:MAXPTS), t1(2,MAXPTS)
      DOUBLE PRECISION    t2(MAXPTS),  t3(2,MAXPTS)
      DOUBLE PRECISION    t4(MAXPTS),  ulp,         ulps1,       ulps2
      DOUBLE PRECISION    ulps3,       ulps4,       w(MAXPTS)
      DOUBLE PRECISION    ww(MAXPTS),  wwxm1(MAXPTS)
      DOUBLE PRECISION    x(MAXPTS),   xarg,        xx(MAXPTS)
      DOUBLE PRECISION    yy(MAXPTS),  zz(MAXPTS)
*
*     DOUBLE PRECISION    tmp1,        tmp2,        tmp3
*
      INTEGER             i,           ierr,        n,           nquad
      INTEGER             p,           pmax
*
      f(parg,xarg) = xarg**parg
      fprime(parg,xarg) = parg * xarg**(parg - ONE)
*
*     CALL ieee_handler( 'set', 'common', fphand)
*     CALL ieee_handler( 'set', 'underflow', fphand)
*
      ulp = deps(ONE)
*
      READ (5, *, END = 1200) alpha, nquad, pmax
      IF (nquad .GT. MAXPTS) STOP '(nquad .GT. MAXPTS)'
      IF (nquad .LT. 1) STOP '(nquad .LT. 1)'
      IF (alpha .LE. -ONE) STOP '(alpha .LE. -ONE)'
      WRITE (6, 10000) alpha, nquad, pmax
*
      CALL glqrc(a, b, s, t, alpha, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'glqrc() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
*
      WRITE (6, 90000)
      WRITE (6, 96000) (i, a(i), b(i), i = 0, nquad)
      WRITE (6, 95000)
      WRITE (6, 96000) (i, t(i), s(i), i = 0, nquad)
*
      CALL glqf(xx, ww, wwxm1, yy, zz, alpha, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'glqf() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
      WRITE (6, 11000)
      WRITE (6, 13000)
      WRITE (6, 14000) (i, xx(i), ww(i), wwxm1(i), i = 1, nquad)
      WRITE (6, 15000)
      WRITE (6, 16000) (i, yy(i), zz(i), i = 1, nquad)
*
      CALL glqfd(x, w, deltaw, deltax, alpha, nquad, ierr)
      IF (ierr .NE. 0) WRITE (6, *) 'glqfd() returns ierr =', ierr
      IF (ierr .NE. 0) STOP
      WRITE (6, 12000)
      WRITE (6, 17000)
      WRITE (6, 16000) (i, x(i), w(i), i = 1, nquad)
      WRITE (6, 18000)
      WRITE (6, 16000) (i, deltaw(i), deltax(i), i = 1, nquad)
*
*     The nodes and weights for the nonlogarithmic quadature should be
*     identical from both glqf() and glqfd(): check this.
*
      ierr = 0
      DO 100 i = 1, nquad
          IF (x(i) .NE. xx(i)) ierr = ierr + 1
          IF (w(i) .NE. ww(i)) ierr = ierr + 1
  100 CONTINUE
      IF (ierr .NE. 0) WRITE (6, *) 'WARNING: unexpected difference',
     X    'in (x,w) from glqf() and glqfd()'
*
      WRITE (6, 70000)
      WRITE (6, 30000)
      DO 300 p = 0, pmax
           fp = dble(p)
           fexact = dgamma(alpha + ONE + fp) * dpsi(alpha + ONE + fp)
           DO 200 i = 1, nquad
                t1(1,i) = deltaw(i)*f(fp,x(i))
                t1(2,i) = deltax(i)*fprime(fp,x(i))
                t2(i) = w(i)*dlog(x(i))*f(fp,x(i))
  200      CONTINUE
           sum1 = dvsum(t1, 2 * nquad)
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
           fexact = dgamma(alpha + ONE + fp) * dpsi(alpha + ONE + fp)
           DO 400 i = 1, nquad
*               tmp1 = wwxm1(i) * f(fp,xx(i))
*               tmp2 = zz(i) * f(fp,yy(i))
*               tmp3 = dabs((tmp1 - tmp2)/tmp1)
*               if (tmp3 .LE. 0.5d+00)
*    X              write (6, *) 'DEBUG: i =', i, ' bitloss = ',
*    X                  idint(dabs(dlog2(tmp3)) + 0.5d+00)
                t3(1,i) = wwxm1(i) * f(fp,xx(i))
                t3(2,i) = -zz(i) * f(fp,yy(i))
                t4(i) = ww(i)*dlog(xx(i))*f(fp,xx(i))
  400      CONTINUE
           sum3 = dvsum(t3, 2 * nquad)
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
           fexact = dgamma(alpha + ONE + fp) * dpsi(alpha + ONE + fp)
           DO 600 i = 1, nquad
                t1(1,i) = deltaw(i)*f(fp,x(i))
                t1(2,i) = deltax(i)*fprime(fp,x(i))
                t2(i) = w(i)*dlog(x(i))*f(fp,x(i))
                t3(1,i) = wwxm1(i) * f(fp,xx(i))
                t3(2,i) = -zz(i) * f(fp,yy(i))
                t4(i) = ww(i)*dlog(xx(i))*f(fp,xx(i))
  600      CONTINUE
           sum1 = dvsum(t1, 2 * nquad)
           sum2 = dvsum(t2, nquad)
           sum3 = dvsum(t3, 2 * nquad)
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
           fexact = dgamma(alpha + fp + ONE)
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
      IF (alpha .EQ. ZERO) THEN
          WRITE (6, 99100)
          DO 1100 n = 2, 68
              IF ((n .LE. 32) .OR. (mod(n,4) .EQ. 0)) THEN
                  CALL glqfd(x, w, deltaw, deltax, alpha, n, ierr)
                  IF (ierr .NE. 0) WRITE (6, *)
     X                'gjqfd() returns ierr =', ierr
                  IF (ierr .NE. 0) STOP
                  WRITE (6,99200) n
                  WRITE (6,99300)
                  DO 1000 i = 1, n
                      c1 = pad(x(i))
                      c2 = pad(w(i))
                      WRITE (6,99400) i, c1, c2
 1000             CONTINUE
              END IF
 1100     CONTINUE

      END IF
*
 1200 CONTINUE
*
10000 FORMAT (2x, 'alpha =', f10.6, '  nquad =', i5, '  pmax =', i4)
11000 FORMAT (//, 2x, 'Nodes and weights from glqf()')
12000 FORMAT (//, 2x, 'Nodes and weights from glqfd()')
13000 FORMAT (2x, 4x, 'n', 1x, 11x, 'x(n)', 10x, 1x, 11x, 'w(n)', 10x,
     X    1x, 11x, 'wxm1(n)')
14000 FORMAT (2x, i5, 1x, 1pe25.16, 1x, 1pe25.16, 1x, 1pe25.16)
15000 FORMAT (/, 2x, 4x, 'n', 2x, 10x, 'y(n)', 11x, 2x, 9x, 'z(n)')
16000 FORMAT (2x, i5, 1x, 1pe25.16, 1x, 1pe25.16)
17000 FORMAT (2x, 4x, 'n', 1x, 11x, 'x(n)', 10x, 1x, 11x, 'w(n)')
18000 FORMAT (/, 2x, 4x, 'n', 1x, 10x, 'deltaw(n)', 6x,
     X    1x, 10x, 'deltax(n)')
20000 FORMAT (2x, i5, 1pe25.16, 2(1pe25.16, 1x, 1pe11.2, 1x, 0pf9.2))
30000 FORMAT (//, 2x, 4x, 'p', 7x, 'Exact Integral', 16x, 'glqfd log',
     X    6x, 'RelErr', 10x, 'Ulps', 10x, 'glqfd nolog', 6x, 'RelErr',
     X    10x, 'Ulps')
40000 FORMAT (//,2x, 4x, 'p', 7x, 'Exact Integral', 17x, 'glqf log',
     X    6x, 'RelErr', 10x, 'Ulps', 11x, 'glqf nolog', 6x, 'RelErr',
     X    10x, 'Ulps')
50000 FORMAT (//, 2x, 'Comparative relative errors', /,
     X    2x, 5x, '--------------------- glqfd ----------------------',
     X    1x, ' -------------------- glqf ----------------------', /,
     X    2x, 22x, '(1-x)**p',
     X    8x,  'log(x) * (1-x)**p',
     X    17x, '(1-x)**p',
     X    8x,  'log(x) * (1-x)**p', /,
     X    2x, 4x, 'p', 13x, 'log-Laguerre', 7x, 'Gauss-log-Laguerre',
     X    11x, 'Gauss-Laguerre', 11x, 'Gauss-Laguerre')
60000 FORMAT (2x, i5, 1x, 1pe24.2, 1x, 1pe24.2, 1x, 1pe24.2,
     X    1x, 1pe24.2)
70000 FORMAT (//, 2x,
     X    '\\int_0^\\infty x^\\alpha \\exp(-x) \\ln x x^p dx')
80000 FORMAT (//, 2x, '\\int_0^\\infty x^\\alpha \\exp(-x) x^p dx')
85000 FORMAT (//,2x, 4x, 'p', 7x, 'Exact Integral', 19x, 'glqfd',
     X    7x, 'RelErr', 10x, 'Ulps', 14x, 'glqf', 9x, 'RelErr',
     X    10x, 'Ulps')
90000 FORMAT (//, 2x, 'Monic polynomial recursion coefficients'/
     X       1x, 3x, 'n', 2x, 10x, 'a(n)', 11x, 2x, 10x, 'b(n)')
95000 FORMAT (//, 2x, 'Monic polynomial zeroth and first moments'/
     X       1x, 3x, 'n', 2x, 10x, 't(n)', 11x, 2x, 10x, 's(n)')
96000 FORMAT (1x, i4, 2x, 1pe25.16, 2x, e25.16)
99100 FORMAT(/,
     X       2x, 'Weights and nodes for ordinary Gauss-Laguerre ',
     X           'quadrature: see', /,
     X       2x, 'Table 6 (pp. 253--274) in the book', //,
     X       10x, 'A. H. Stroud and Don Secrest', /,
     X       10x, 'Gaussian Quadrature Formulas', /,
     X       10x, 'Prentice-Hall 1966', /)
99200 FORMAT(//, 2x, 4x, 30x, 'N = ', i2)
99300 FORMAT(2x, 4x, 2x, 13x, 'x(i)', 13x, 2x, 13x, 'A(i)')
99400 FORMAT(2x, i4, 2x, a, 2x, a)
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
