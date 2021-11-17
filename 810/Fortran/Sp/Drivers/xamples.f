C  SUBROUTINE EXAMP

C     **********
C     MARCH 1, 2001; P.B. BAILEY, W.N. EVERITT AND A. ZETTL
C     VERSION 1.2
C     **********

      SUBROUTINE EXAMP
C
C
C     THIS SUBROUTINE CONTAINS A SELECTION OF COEFFICIENT FUNCTIONS
C     p,q,w (AND POSSIBLY SUITABLE FUNCTIONS u,v WITH WHICH TO 
C     DEFINE BOUNDARY CONDITIONS AT LIMIT CIRCLE ENDPOINTS) 
C     WHICH DEFINE SOME INTERESTING STURM-LIOUVILLE BOUNDARY
C     VALUE PROBLEMS. IT CAN BE CALLED BY THE MAIN PROGRAM, DRIVE.
C
C     .. Scalars in Common ..
      REAL A,ALPHA,B,BETA,C,D,E,GAMMA,H,K,L,NU,S
      INTEGER NUMBER
C     ..
C     .. Local Scalars ..
      REAL TMP
      INTEGER MUMBER
      CHARACTER ANS
C     ..
C     .. Common blocks ..
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA,S,A,B,C,D,E
C     ..
      WRITE (*,FMT=*)
     +  ' Here is a collection of 32 differential equations '
      WRITE (*,FMT=*) ' which can be used with SLEIGN2.  By typing an '
      WRITE (*,FMT=*)
     +  ' integer from 1 to 32, one of these differential '
      WRITE (*,FMT=*)
     +  ' equations is selected, whereupon its coefficient '
      WRITE (*,FMT=*) ' functions p,q,w will be displayed along with a '
      WRITE (*,FMT=*) ' brief description of its singular points.  The'
      WRITE (*,FMT=*) ' endpoints a, b of the interval over which the '
      WRITE (*,FMT=*)
     +  ' differential equation is integrated are specified '
      WRITE (*,FMT=*) ' later; any interval which does not contain '
      WRITE (*,FMT=*) ' singular points in its interior is acceptable. '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' DO YOU WISH TO CONTINUE ? (Y/N) '
      READ (*,FMT=*) ANS
      IF (.NOT. (ANS.EQ.'y'.OR.ANS.EQ.'Y')) STOP
   35 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) '  1 IS THE LEGENDRE EQUATION '
      WRITE (*,FMT=*) '  2 IS THE BESSEL EQUATION '
      WRITE (*,FMT=*) '  3 IS THE HALVORSEN EQUATION '
      WRITE (*,FMT=*) '  4 IS THE BOYD EQUATION '
      WRITE (*,FMT=*) '  5 IS THE REGULARIZED BOYD EQUATION '
      WRITE (*,FMT=*) '  6 IS THE SEARS-TITCHMARSH EQUATION '
      WRITE (*,FMT=*) '  7 IS THE BEZ EQUATION '
      WRITE (*,FMT=*) '  8 IS THE LAPLACE TIDAL WAVE EQUATION '
      WRITE (*,FMT=*) '  9 IS THE LATZKO EQUATION '
      WRITE (*,FMT=*) ' 10 IS A WEAKLY REGULAR EQUATION '
      WRITE (*,FMT=*) ' 11 IS THE PLUM EQUATION '
      WRITE (*,FMT=*) ' 12 IS THE MATHIEU PERIODIC EQUATION '
      WRITE (*,FMT=*) ' 13 IS THE HYDROGEN ATOM EQUATION '
      WRITE (*,FMT=*) ' 14 IS THE MARLETTA EQUATION '
      WRITE (*,FMT=*) ' 15 IS THE HARMONIC OSCILLATOR EQUATION '
      WRITE (*,FMT=*) ' 16 IS THE JACOBI EQUATION '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' Press any key to continue. '
      READ (*,FMT=9010) ANS
 9010 FORMAT (A1)
      WRITE (*,FMT=*) ' 17 IS THE ROTATION MORSE OSCILLATOR EQUATION '
      WRITE (*,FMT=*) ' 18 IS THE DUNSCH EQUATION '
      WRITE (*,FMT=*) ' 19 IS THE DONSCH EQUATION '
      WRITE (*,FMT=*) ' 20 IS THE KRALL EQUATION '
      WRITE (*,FMT=*) ' 21 IS THE FOURIER EQUATION '
      WRITE (*,FMT=*) ' 22 IS THE LAGUERRE EQUATION '
      WRITE (*,FMT=*) ' 23 IS THE LAGUERRE/LIOUVILLE FORM EQUATION '
      WRITE (*,FMT=*) ' 24 IS THE JACOBI/LIOUVILLE FORM EQUATION '
      WRITE (*,FMT=*) ' 25 IS THE MEISSNER EQUATION '
      WRITE (*,FMT=*) ' 26 IS THE LOHNER EQUATION '
      WRITE (*,FMT=*) ' 27 IS THE JOERGENS EQUATION '
      WRITE (*,FMT=*) ' 28 IS THE BEHNKE-GOERISCH EQUATION '
      WRITE (*,FMT=*) ' 29 IS THE WHITTAKER EQUATION '
      WRITE (*,FMT=*) ' 30 IS THE LITTLEWOOD-MCLEOD EQUATION '
      WRITE (*,FMT=*) ' 31 IS THE MORSE EQUATION '
      WRITE (*,FMT=*) ' 32 IS THE HEUN EQUATION '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ENTER THE NUMBER OF YOUR CHOICE: '
      READ (*,FMT=*) NUMBER
      IF (NUMBER.LT.1 .OR. NUMBER.GT.32) GO TO 35
      MUMBER = NUMBER
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26,27,28,29,30,31,32) NUMBER
C
    1 CONTINUE
      WRITE (*,FMT=*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1 - x*x, q = 1/4, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT CIRCLE, NON-OSCILLATORY AT -1. '
      WRITE (*,FMT=*) ' LIMIT CIRCLE, NON-OSCILLATORY AT +1. '
      GO TO 40
C
    2 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = (nu*nu-0.25)/x*x, w = 1 '
      WRITE (*,FMT=*) '            (nu a parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = 0: '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR -1.LT.nu.LT.1 BUT nu*nu.NE.0.25. '
      WRITE (*,FMT=*) '    REGULAR FOR nu*nu = 0.25. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR nu*nu.GE.1.0. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +INFINITY: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR ALL nu. '
      MUMBER = 101
      GO TO 40
C
    3 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 0, w = exp(-2/x)/x**4 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' WEAKLY REGULAR AT 0. '
      WRITE (*,FMT=*) ' LIMIT CIRCLE, NON-OSCILLATORY AT +INFINITY. '
      GO TO 40
C
    4 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,0) '
      WRITE (*,FMT=*)
     +  '                            AND on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = -1/x, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT CIRCLE, NON-OSCILLATORY AT 0+ AND 0-. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
    5 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,0) '
      WRITE (*,FMT=*)
     +  '                            AND on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = r*r, q = -r*r*(ln|x|)**2, w = r*r '
      WRITE (*,FMT=*) '      where r = exp(-(x*ln(|x|)-x)) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' WEAKLY REGULAR AT 0+ AND 0-. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
    6 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = x, q = -x, w = 1/x '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT 0. '
      WRITE (*,FMT=*) ' LIMIT CIRCLE, OSCILLATORY AT +INFINITY. '
      GO TO 40
C
    7 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,0) '
      WRITE (*,FMT=*)
     +  '                            AND on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = x, q = -1/x, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT CIRCLE, OSCILLATORY AT 0+ AND 0-. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
    8 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1/x, q = (k/x**2) + (k**2/x), w = 1 '
      WRITE (*,FMT=*) '         (k a non-zero parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT CIRCLE, NON-OSCILLATORY AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 102
      GO TO 40
C
    9 CONTINUE
      WRITE (*,FMT=*) '    -(p*y'')'' + q*y = lambda*w*y on (0,1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1 - x**7, q = 0, w = x**7 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' WEAKLY REGULAR AT 0. '
      WRITE (*,FMT=*) ' LIMIT CIRCLE, NON-OSCILLATORY AT +1. '
      GO TO 40
C
   10 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = sqrt(x), q = 0, w = 1/sqrt(x) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' WEAKLY REGULAR AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   11 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 100*cos(x)**2, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   12 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 2*k*cos(2x), w = 1 '
      WRITE (*,FMT=*) '       (k a non-zero parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 103
      GO TO 40
C
   13 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = k/x + h/x**2, w = 1 '
      WRITE (*,FMT=*) '        q = k/x + h/x**2 + 1  if h.lt.-0.25  '
      WRITE (*,FMT=*) '        (h,k parameters) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = 0: '
      WRITE (*,FMT=*) '    REGULAR for h = k = 0. '
      WRITE (*,FMT=*) '    LIMIT-CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR h = 0 AND ALL k.NE.0. '
      WRITE (*,FMT=*) '    LIMIT-CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*)
     +  '       FOR -0.25.LE.h.LT.0.75 BUT h.NE.0, AND ALL k. '
      WRITE (*,FMT=*) '    LIMIT-CIRCLE, OSCILLATORY '
      WRITE (*,FMT=*) '       FOR h.LT.-0.25 AND ALL k. '
      WRITE (*,FMT=*) '      (Here, 1 has been added to the usual q so '
      WRITE (*,FMT=*) '       at least some eigenvalues are positive.) '
      WRITE (*,FMT=*) '    LIMIT POINT FOR h.GE.0.75 AND ALL k. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +INFINITY: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR ALL h,k. '
      MUMBER = 104
      GO TO 40
C
   14 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' p = 1, q = 3.0*(X-31.0)/(4.0*(X+1.0)*(4.0+X)**2), '
      WRITE (*,FMT=*) ' w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' REGULAR AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   15 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = x*x, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   16 CONTINUE
      WRITE (*,FMT=*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = (1-x)**(alpha+1)*(1+x)**(beta+1), '
      WRITE (*,FMT=*) ' q = 0, w = (1-x)**alpha*(1+x)**beta '
      WRITE (*,FMT=*) '            (alpha, beta parameters) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = -1.0: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR beta.LE.-1. '
      WRITE (*,FMT=*) '    WEAKLY REGULAR FOR -1.LT.beta.LT.0. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR 0.LE.beta.LT.1. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR beta.GE.1. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +1.0: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.LE.-1. '
      WRITE (*,FMT=*) '    WEAKLY REGULAR FOR -1.LT.alpha.LT.0. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR 0.LE.alpha.LT.1. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.GE.1. '
      MUMBER = 105
      GO TO 40
C
   17 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 2/x**2 - 2000(2e-e*e), w = 1 '
      WRITE (*,FMT=*) '       where e = exp(-1.7(x-1.3)) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   18 CONTINUE
      WRITE (*,FMT=*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' p = 1 - x*x, q = 2*alpha**2/(1+x) + 2*beta**2/(1-x),'
      WRITE (*,FMT=*)
     +  ' w = 1        (alpha, beta non-negative parameters) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = -1.0: '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR 0.LE.alpha.LT.0.5. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.GE.0.5. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +1.0: '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR 0.LE.beta.LT.0.5. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR beta.GE.0.5. '
      MUMBER = 106
      GO TO 40
C
   19 CONTINUE
      WRITE (*,FMT=*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' p = 1 - x*x, q = -2*gamma**2/(1+x) +2*beta**2/(1-x),'
      WRITE (*,FMT=*) ' w = 1            (gamma, beta parameters) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = -1.0: '
      WRITE (*,FMT=*)
     +  '    LIMIT CIRCLE, NON-OSCILLATORY FOR gamma = 0. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, OSCILLATORY FOR gamma.GT.0. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +1.0: '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR 0.LE.beta.LT.0.5. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR beta.GE.0.5. '
      MUMBER = 107
      GO TO 40
C
   20 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 1 - (k**2+0.25)/x**2, w = 1 '
      WRITE (*,FMT=*) '        (k a positive parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT CIRCLE, OSCILLATORY AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 108
      GO TO 40
C
   21 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 0, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   22 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = x**(alpha+1)*exp(-x), w = x**alpha*exp(-x) '
      WRITE (*,FMT=*) ' q = 0        (alpha a parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = 0.: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.LE.-1. '
      WRITE (*,FMT=*) '    WEAKLY REGULAR FOR -1.LT.alpha.LT.0. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR 0.LE.alpha.LT.1. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.GE.1. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +INFINITY: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR ALL alpha. '
      MUMBER = 109
      GO TO 40
C
   23 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, w = 1 '
      WRITE (*,FMT=*)
     +  ' q = (alpha**2-0.25)/x**2 - (alpha+1)/2 + x**2/16'
      WRITE (*,FMT=*) '              (alpha a parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = 0.: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.LE.-1. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*)
     +  '       FOR -1.LT.alpha.LT.1 BUT alpha**2.NE.0.25. '
      WRITE (*,FMT=*) '    REGULAR FOR alpha**2 = 0.25. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.GE.1. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +INFINITY: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR ALL alpha. '
      MUMBER = 110
      GO TO 40
C
   24 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-pi/2,pi/2) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, w = 1 '
      WRITE (*,FMT=*) ' q = (beta**2-0.25)/(4*tan((x+pi/2)/2)**2)+ '
      WRITE (*,FMT=*) '     (alpha**2-0.25)/(4*tan((x-pi/2)/2)**2)- '
      WRITE (*,FMT=*) '     (4*alpha*beta+4*alpha+4*beta+3)/8 '
      WRITE (*,FMT=*) '            (alpha, beta parameters) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT X = -pi/2: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR beta.LE.-1. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '       FOR -1.LT.beta.LT.1 BUT beta**2.NE.0.25. '
      WRITE (*,FMT=*) '    REGULAR FOR beta**2 = 0.25. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR beta.GE.1. '
      WRITE (*,FMT=*) ' THE ENDPOINT X = +pi/2: '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.LE.-1. '
      WRITE (*,FMT=*) '    LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*)
     +  '       FOR -1.LT.alpha.LT.1 BUT alpha**2.NE.0.25. '
      WRITE (*,FMT=*) '    REGULAR FOR alpha**2 = 0.25. '
      WRITE (*,FMT=*) '    LIMIT POINT FOR alpha.GE.1. '
      MUMBER = 111
      GO TO 40
C
   25 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 0 '
      WRITE (*,FMT=*) ' w = 1 when x.le.0. '
      WRITE (*,FMT=*) '   = 9 when x.gt.0. '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   26 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = -1000*x, w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   27 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 0.25*exp(2*x) - k*exp(x), w = 1 '
      WRITE (*,FMT=*) '                (k a parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 112
      GO TO 40
C
   28 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE (*,FMT=*) '                                    +infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = k*cos(x)**2, w = 1 '
      WRITE (*,FMT=*) '          (k a parameter) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 113
      GO TO 40
C
   29 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 0.25 + (k**2-1)/(4*x**2)), w = 1/x '
      WRITE (*,FMT=*) '       (k a positive parameter .GE. 1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 114
      GO TO 40
C
   30 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = x*sin(x), w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' REGULAR AT 0. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   31 CONTINUE
      WRITE (*,FMT=*)
     +  '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,+infinity) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = 1, q = 9*exp(-2*x) - 18*exp(-x), w = 1 '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' LIMIT POINT AT -INFINITY. '
      WRITE (*,FMT=*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   32 CONTINUE
      WRITE (*,FMT=*) '    -(p*y'')'' + q*y = lambda*w*y on (0,1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' p = x**c*(1-x)**d*(x+s)**e '
      WRITE (*,FMT=*) ' q = a*b*x**c*(1-x)**(d-1)*(x+s)**(e-1) '
      WRITE (*,FMT=*) ' w = x**(c-1)*(1-x)**(d-1)*(x+s)**(e-1) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' where the parameters s,a,b,c,d,e satisfy '
      WRITE (*,FMT=*) ' the conditions '
      WRITE (*,FMT=*) '   s > 0, a.ge.b.ge.1, c.ge.1, d.ge.1, '
      WRITE (*,FMT=*) '   a+b+1-c-d-e = 0. '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' THE ENDPOINT 0: '
      WRITE (*,FMT=*) '   LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '     FOR 1.LE.c.LT.2; '
      WRITE (*,FMT=*) '   LIMIT POINT FOR 2.LE.c. '
      WRITE (*,FMT=*) ' THE ENDPOINT 1: '
      WRITE (*,FMT=*) '   LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE (*,FMT=*) '     FOR 1.LE.d.LT.2; '
      WRITE (*,FMT=*) '   LIMIT POINT FOR 2.LE.d. '
      MUMBER = 115
      GO TO 40
C
   40 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' IS THIS THE CORRECT DIFFERENTIAL EQUATION ? (Y/N) '
      READ (*,FMT=*) ANS
      IF (.NOT. (ANS.EQ.'y'.OR.ANS.EQ.'Y')) GO TO 35
      IF (MUMBER.EQ.NUMBER) RETURN
C
C     Now enter any parameters needed for these D.E.'s, or defaults.
C
      MUMBER = MUMBER - 100
      GO TO (41,42,43,44,45,46,47,48,49,50,51,52,53,54,55) MUMBER
C
   41 CONTINUE
      NU = 1.0D0
      WRITE (*,FMT=*) ' Choose real parameter nu, nu = '
      READ (*,FMT=*) NU
      RETURN
C
   42 CONTINUE
      K = 1.0D0
      WRITE (*,FMT=*) ' Choose real parameter k.ne.0., k = '
      READ (*,FMT=*) K
      RETURN
C
   43 CONTINUE
      K = 1.0D0
      WRITE (*,FMT=*) ' Choose real parameter k = '
      READ (*,FMT=*) K
      RETURN
C
   44 CONTINUE
      H = 1.0D0
      K = 1.0D0
      WRITE (*,FMT=*) ' Choose real parameters k,h = '
      READ (*,FMT=*) K,H
      RETURN
C
   45 CONTINUE
      BETA = 0.1D0
      ALPHA = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameters alpha,beta = '
      READ (*,FMT=*) ALPHA,BETA
      RETURN
C
   46 CONTINUE
      ALPHA = 0.1D0
      BETA = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameters alpha,beta = '
      READ (*,FMT=*) ALPHA,BETA
      RETURN
C
   47 CONTINUE
      GAMMA = 0.1D0
      BETA = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameters gamma,beta = '
      READ (*,FMT=*) GAMMA,BETA
      RETURN
C
   48 CONTINUE
      K = 1.0D0
      WRITE (*,FMT=*) ' Choose real parameter k.gt.0., k = '
      READ (*,FMT=*) K
      RETURN
C
   49 CONTINUE
      ALPHA = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameter alpha = '
      READ (*,FMT=*) ALPHA
      RETURN
C
   50 CONTINUE
      ALPHA = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameter alpha = '
      READ (*,FMT=*) ALPHA
      RETURN
C
   51 CONTINUE
      BETA = 0.1D0
      ALPHA = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameters alpha,beta = '
      READ (*,FMT=*) ALPHA,BETA
      RETURN
C
   52 CONTINUE
      K = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameter k = '
      READ (*,FMT=*) K
      RETURN
C
   53 CONTINUE
      K = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameter k = '
      READ (*,FMT=*) K
      RETURN
C
   54 CONTINUE
      K = 0.1D0
      WRITE (*,FMT=*) ' Choose real parameter k = '
      READ (*,FMT=*) K
      RETURN
C
   55 CONTINUE
      WRITE (*,FMT=*) ' Choose real parameter s > 0 '
      READ (*,FMT=*) S
      WRITE (*,FMT=*) ' Before entering the numerical values of the '
      WRITE (*,FMT=*) '   parameters a,b,c,d,e the user is advised to '
      WRITE (*,FMT=*) '   make a preliminary choice of these numbers '
      WRITE (*,FMT=*) '   that are consistent with the conditions '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) '    a.ge.b, b.ge.1 '
      WRITE (*,FMT=*) '    1.le.d .le. (a+b-1) '
      WRITE (*,FMT=*) '    1.le.c .le. (a+b-d) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' (Parameter e will be set equal to a+b+1-c-d.) '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' Choose parameters a, b  = '
      READ (*,FMT=*) A,B
      TMP = A + B - 1.0D0
      WRITE (*,FMT=*) ' Choose parameter d, 1.LE.d .LE. ',TMP
      READ (*,FMT=*) D
      TMP = A + B - D
      WRITE (*,FMT=*) ' Choose parameter c, 1.LE.c .LE. ',TMP
      READ (*,FMT=*) C
      E = A + B + 1.0D0 - C - D
      RETURN
      END
C
      REAL FUNCTION P(X)

C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Scalars in Common ..
      REAL A,ALPHA,B,BETA,C,D,E,GAMMA,H,K,L,NU,S
      INTEGER NUMBER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,LOG,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA,S,A,B,C,D,E
C     ..
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26,27,28,29,30,31,32) NUMBER
C
    1 CONTINUE
      P = 1.0D0 - X*X
      RETURN
C
    2 CONTINUE
      P = 1.0D0
      RETURN
C
    3 CONTINUE
      P = 1.0D0
      RETURN
C
    4 CONTINUE
      P = 1.0D0
      RETURN
C
    5 CONTINUE
      P = EXP(-2.0D0* (X*LOG(ABS(X))-X))
      RETURN
C
    6 CONTINUE
      P = X
      RETURN
C
    7 CONTINUE
      P = X
      RETURN
C
    8 CONTINUE
      P = 1.0D0/X
      RETURN
C
    9 CONTINUE
      P = 1.0D0 - X**7
      RETURN
C
   10 CONTINUE
      P = SQRT(X)
      RETURN
C
   11 CONTINUE
      P = 1.0D0
      RETURN
C
   12 CONTINUE
      P = 1.0D0
      RETURN
C
   13 CONTINUE
      P = 1.0D0
      RETURN
C
   14 CONTINUE
      P = 1.0D0
      RETURN
C
   15 CONTINUE
      P = 1.0D0
      RETURN
C
   16 CONTINUE
      IF (ABS(X).EQ.1.0D0) THEN
          P = 0.0D0
      ELSE IF (ALPHA.NE.-1.0D0 .AND. BETA.NE.-1.0D0) THEN
          P = (1.0D0-X)** (ALPHA+1.0D0)* (1.0D0+X)** (BETA+1.0D0)
      ELSE IF (ALPHA.NE.-1.0D0 .AND. BETA.EQ.-1.0D0) THEN
          P = (1.0D0-X)** (ALPHA+1.0D0)
      ELSE IF (ALPHA.EQ.-1.0D0 .AND. BETA.NE.-1.0D0) THEN
          P = (1.0D0+X)** (BETA+1.0D0)
      ELSE
          P = 1.0D0
      END IF
      RETURN
C
   17 CONTINUE
      P = 1.0D0
      RETURN
C
   18 CONTINUE
      P = 1.0D0 - X*X
      RETURN
C
   19 CONTINUE
      P = 1.0D0 - X*X
      RETURN
C
   20 CONTINUE
      P = 1.0D0
      RETURN
C
   21 CONTINUE
      P = 1.0D0
      RETURN
C
   22 CONTINUE
      P = EXP(-X)
      IF (ALPHA.NE.-1.0D0) P = P*X** (ALPHA+1.0D0)
      RETURN
C
   23 CONTINUE
      P = 1.0D0
      RETURN
C
   24 CONTINUE
      P = 1.0D0
      RETURN
C
   25 CONTINUE
      P = 1.0D0
      RETURN
C
   26 CONTINUE
      P = 1.0D0
      RETURN
C
   27 CONTINUE
      P = 1.0D0
      RETURN
C
   28 CONTINUE
      P = 1.0D0
      RETURN
C
   29 CONTINUE
      P = 1.0D0
      RETURN
C
   30 CONTINUE
      P = 1.0D0
      RETURN
C
   31 CONTINUE
      P = 1.0D0
      RETURN
C
   32 CONTINUE
      P = X**C* (1.0D0-X)**D* (X+S)**E
      RETURN
      END
C
      REAL FUNCTION Q(X)
c
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Scalars in Common ..
      REAL A,ALPHA,B,BETA,C,D,E,GAMMA,H,K,L,NU,S
      INTEGER NUMBER
C     ..
C     .. Local Scalars ..
      REAL EE,HPI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,EXP,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA,S,A,B,C,D,E
C     ..
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26,27,28,29,30,31,32) NUMBER
C
    1 CONTINUE
      Q = 0.25D0
      RETURN
C
    2 CONTINUE
      Q = 0.0D0
      IF (NU.NE.-0.5D0 .AND. NU.NE.0.5D0) Q = (NU*NU-0.25D0)/X**2
      RETURN
C
    3 CONTINUE
      Q = 0.0D0
      RETURN
C
    4 CONTINUE
      Q = -1.0D0/X
      RETURN
C
    5 CONTINUE
      Q = -EXP(-2.0D0* (X*LOG(ABS(X))-X))*LOG(ABS(X))**2
      RETURN
C
    6 CONTINUE
      Q = -X
      RETURN
C
    7 CONTINUE
      Q = -1.0D0/X
      RETURN
C
    8 CONTINUE
      Q = K/X**2 + K**2/X
      RETURN
C
    9 CONTINUE
      Q = 0.0D0
      RETURN
C
   10 CONTINUE
      Q = 0.0D0
      RETURN
C
   11 CONTINUE
      Q = 100.0D0*COS(X)**2
      RETURN
C
   12 CONTINUE
      Q = 2.0D0*K*COS(2.0D0*X)
      RETURN
C
   13 CONTINUE
      Q = K/X + H/ (X*X)
      IF (H.LT.-0.25D0) Q = Q + 1.0D0
      RETURN
C
   14 CONTINUE
      Q = 3.0D0* (X-31.0D0)/ (4.0D0* (X+1.0D0)* (X+4.0D0)**2)
      RETURN
C
   15 CONTINUE
      Q = X*X
      RETURN
C
   16 CONTINUE
      Q = 0.0D0
      RETURN
C
   17 CONTINUE
      L = 1.0D0
      EE = EXP(-1.7D0* (X-1.3D0))
      Q = L* (L+1.0D0)/X**2 - 2000.0D0*EE* (2.0D0-EE)
      RETURN
C
   18 CONTINUE
      IF (ALPHA.NE.0.0D0 .AND. BETA.NE.0.0D0) THEN
          Q = 2.0D0*ALPHA**2/ (1.0D0+X) + 2.0D0*BETA**2/ (1.0D0-X)
      ELSE IF (ALPHA.EQ.0.0D0 .AND. BETA.NE.0.0D0) THEN
          Q = 2.0D0*BETA**2/ (1.0D0-X)
      ELSE IF (ALPHA.NE.0.0D0 .AND. BETA.EQ.0.0D0) THEN
          Q = 2.0D0*ALPHA**2/ (1.0D0+X)
      ELSE
          Q = 0.0D0
      END IF
      RETURN
C
   19 CONTINUE
      Q = -2.0D0*GAMMA**2/ (1.0D0+X) + 2.0D0*BETA**2/ (1.0D0-X)
      RETURN
C
   20 CONTINUE
      Q = 1.0D0 - (K**2+0.25D0)/X**2
      RETURN
C
   21 CONTINUE
      Q = 0.0D0
      RETURN
C
   22 CONTINUE
      Q = 0.0D0
      RETURN
C
   23 CONTINUE
      Q = - (ALPHA+1.0D0)/2.0D0 + X**2/16.0D0
      IF (ALPHA.NE.0.5D0 .AND. ALPHA.NE.-0.5D0) Q = Q +
     +    (ALPHA**2-0.25D0)/X**2
      RETURN
C
   24 CONTINUE
      HPI = 2.0D0*ATAN(1.0D0)
      IF (BETA*BETA.NE.0.25D0 .AND. ALPHA*ALPHA.NE.0.25D0) THEN
          Q = (BETA**2-0.25D0)/ (4.0D0*TAN((X+HPI)/2.0D0)**2) +
     +        (ALPHA**2-0.25D0)/ (4.0D0*TAN((X-HPI)/2.0D0)**2) -
     +        (ALPHA*BETA+ALPHA+BETA+0.75D0)/2.0D0
      ELSE IF (BETA*BETA.EQ.0.25D0 .AND. ALPHA*ALPHA.NE.0.25D0) THEN
          Q = (ALPHA**2-0.25D0)/ (4.0D0*TAN((X-HPI)/2.0D0)**2) -
     +        (ALPHA*BETA+ALPHA+BETA+0.75D0)/2.0D0
      ELSE IF (BETA*BETA.NE.0.25D0 .AND. ALPHA*ALPHA.EQ.0.25D0) THEN
          Q = (BETA**2-0.25D0)/ (4.0D0*TAN((X+HPI)/2.0D0)**2) -
     +        (ALPHA*BETA+ALPHA+BETA+0.75D0)/2.0D0
      ELSE
          Q = - (ALPHA*BETA+ALPHA+BETA+0.75D0)/2.0D0
      END IF
      RETURN
C
   25 CONTINUE
      Q = 0.0D0
      RETURN
C
   26 CONTINUE
      Q = -1000.0D0*X
      RETURN
C
   27 CONTINUE
      EE = EXP(X)
      Q = EE* (0.25D0*EE-K)
      RETURN
C
   28 CONTINUE
      Q = K*COS(X)**2
      RETURN
C
   29 CONTINUE
      Q = 0.25D0 + (K**2-1.0D0)/ (4.0D0*X**2)
      RETURN
C
   30 CONTINUE
      Q = X*SIN(X)
      RETURN
C
   31 CONTINUE
      EE = EXP(-X)
      Q = 9.0D0*EE*EE - 18.0D0*EE
      RETURN
C
   32 CONTINUE
      Q = A*B*X**C* (1.0D0-X)** (D-1.0D0)* (X+S)** (E-1.0D0)
      RETURN
      END
C
      REAL FUNCTION W(X)

C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Scalars in Common ..
      REAL A,ALPHA,B,BETA,C,D,E,GAMMA,H,K,L,NU,S
      INTEGER NUMBER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,LOG,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA,S,A,B,C,D,E
C     ..
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26,27,28,29,30,31,32) NUMBER
C
    1 CONTINUE
      W = 1.0D0
      RETURN
C
    2 CONTINUE
      W = 1.0D0
      RETURN
C
    3 CONTINUE
      W = 0.0D0
      IF (X.NE.0.0D0) W = EXP(-2.0D0/X)/X**4
      RETURN
C
    4 CONTINUE
      W = 1.0D0
      RETURN
C
    5 CONTINUE
      W = EXP(-2.0D0* (X*LOG(ABS(X))-X))
      RETURN
C
    6 CONTINUE
      W = 1.0D0/X
      RETURN
C
    7 CONTINUE
      W = 1.0D0
      RETURN
C
    8 CONTINUE
      W = 1.0D0
      RETURN
C
    9 CONTINUE
      W = X**7
      RETURN
C
   10 CONTINUE
      W = 1.0D0/SQRT(X)
      RETURN
C
   11 CONTINUE
      W = 1.0D0
      RETURN
C
   12 CONTINUE
      W = 1.0D0
      RETURN
C
   13 CONTINUE
      W = 1.0D0
      RETURN
C
   14 CONTINUE
      W = 1.0D0
      RETURN
C
   15 CONTINUE
      W = 1.0D0
      RETURN
C
   16 CONTINUE
      IF (ALPHA.NE.0.0D0 .AND. BETA.NE.0.0D0) THEN
          W = (1.0D0-X)**ALPHA* (1.0D0+X)**BETA
      ELSE IF (ALPHA.NE.0.0D0 .AND. BETA.EQ.0.0D0) THEN
          W = (1.0D0-X)**ALPHA
      ELSE IF (ALPHA.EQ.0.0D0 .AND. BETA.NE.0.0D0) THEN
          W = (1.0D0+X)**BETA
      ELSE
          W = 1.0D0
      END IF
      RETURN
C
   17 CONTINUE
      W = 1.0D0
      RETURN
C
   18 CONTINUE
      W = 1.0D0
      RETURN
C
   19 CONTINUE
      W = 1.0D0
      RETURN
C
   20 CONTINUE
      W = 1.0D0
      RETURN
C
   21 CONTINUE
      W = 1.0D0
      RETURN
C
   22 CONTINUE
      W = EXP(-X)
      IF (ALPHA.NE.0.0D0) W = W* (X**ALPHA)
      RETURN
C
   23 CONTINUE
      W = 1.0D0
      RETURN
C
   24 CONTINUE
      W = 1.0D0
      RETURN
C
   25 CONTINUE
      W = 9.0D0
      IF (X.LE.0.0D0) W = 1.0D0
      RETURN
C
   26 CONTINUE
      W = 1.0D0
      RETURN
C
   27 CONTINUE
      W = 1.0D0
      RETURN
C
   28 CONTINUE
      W = 1.0D0
      RETURN
C
   29 CONTINUE
      W = 1.0D0/X
      RETURN
C
   30 CONTINUE
      W = 1.0D0
      RETURN
C
   31 CONTINUE
      W = 1.0D0
      RETURN
C
   32 CONTINUE
      W = X** (C-1.0D0)* (1.0D0-X)** (D-1.0D0)* (X+S)** (E-1.0D0)
      RETURN
      END
C
      SUBROUTINE UV(X,U,PUP,V,PVP,HU,HV)
C
C     HERE, HU MEANS -(pu')' + qu.
C
C     .. Scalar Arguments ..
      REAL HU,HV,PUP,PVP,U,V,X
C     ..
C     .. Scalars in Common ..
      REAL A,ALPHA,B,BETA,C,D,E,GAMMA,H,K,L,NU,S
      INTEGER NUMBER
C     ..
C     .. Local Scalars ..
      REAL CC,EE,HPI,L2,SQ,SS,TX
C     ..
C     .. External Functions ..
      REAL Q,W
      EXTERNAL Q,W
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,EXP,LOG,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA,S,A,B,C,D,E
C     ..
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26,27,28,29,30,31,32) NUMBER
C
    1 CONTINUE
      U = 1.0D0
      PUP = 0.0D0
      V = 0.5D0*LOG((1.0D0+X)/ (1.0D0-X))
      PVP = 1.0D0
      HU = 0.25D0*U
      HV = 0.25D0*V
      RETURN
C
    2 CONTINUE
      IF (NU.NE.-0.5D0 .AND. NU.NE.0.5D0 .AND. NU.NE.0.0D0) THEN
          U = X** (NU+0.5D0)
          PUP = (NU+0.5D0)*X** (NU-0.5D0)
          V = X** (-NU+0.5D0)
          PVP = (-NU+0.5D0)*X** (-NU-0.5D0)
      ELSE IF (NU.EQ.-0.5D0) THEN
          U = X
          PUP = 1.0D0
          V = 1.0D0
          PVP = 0.0D0
      ELSE IF (NU.EQ.0.5D0) THEN
          U = X
          PUP = 1.0D0
          V = -1.0D0
          PVP = 0.0D0
      ELSE IF (NU.EQ.0.0D0) THEN
          U = SQRT(X)
          PUP = 0.5D0/U
          V = U*LOG(X)
          PVP = (0.5D0*LOG(X)+1.0D0)/U
      END IF
      HU = 0.0D0
      HV = 0.0D0
      RETURN
C
    3 CONTINUE
      U = 1.0D0
      V = X
      PUP = 0.0D0
      PVP = 1.0D0
      HU = 0.0D0
      HV = 0.0D0
      RETURN
C
    4 CONTINUE
      TX = LOG(ABS(X))
      U = X
      PUP = 1.0D0
      V = 1.0D0 - X*TX
      PVP = -1.0D0 - TX
      HU = -1.0D0
      HV = TX
      RETURN
C
    5 CONTINUE
      RETURN
C
    6 CONTINUE
      U = (COS(X)+SIN(X))/SQRT(X)
      V = (COS(X)-SIN(X))/SQRT(X)
      PUP = -0.5D0*U + X*V
      PVP = -0.5D0*V - X*U
      HU = -0.25D0*U/X
      HV = -0.25D0*V/X
      RETURN
C
    7 CONTINUE
      TX = LOG(ABS(X))
      U = COS(TX)
      V = SIN(TX)
      PUP = -V
      PVP = U
      HU = 0.0D0
      HV = 0.0D0
      RETURN
C
    8 CONTINUE
      V = X - 1.0D0/K
      U = X*X
      PVP = 1.0D0/X
      PUP = 2.0D0
      HU = K + X*K**2
      HV = K**2
      RETURN
C
    9 CONTINUE
      U = 1.0D0
      V = -LOG(1.0D0-X)
      PUP = 0.0D0
      PVP = (((((X+1.0D0)*X+1.0D0)*X+1.0D0)*X+1.0D0)*X+1.0D0)*X + 1.0D0
      HU = 0.0D0
      HV = - (((((6.0D0*X+5.0D0)*X+4.0D0)*X+3.0D0)*X+2.0D0)*X+1.0D0)
      RETURN
C
   10 CONTINUE
      U = 2.0D0*SQRT(X)
      V = 1.0D0
      PUP = 1.0D0
      PVP = 0.0D0
      HU = 0.0D0
      HV = 0.0D0
      RETURN
C
   11 CONTINUE
      RETURN
C
   12 CONTINUE
      RETURN
C
   13 CONTINUE
      IF (H.GT.-0.25D0) THEN
          L = SQRT(H+0.25D0)
          IF (H.EQ.0.0D0) THEN
              U = X
              V = 1.0D0 + K*X*LOG(X)
              PUP = 1.0D0
              PVP = K* (1.0D0+LOG(X))
              HU = K
              HV = K*K*LOG(X)
          ELSE
              U = X** (0.5D0+L)
              V = X** (0.5D0-L) + (K/ (1.0D0-2.0D0*L))*X** (1.5D0-L)
              PUP = (0.5D0+L)*X** (L-0.5D0)
              PVP = (0.5D0-L)*X** (-L-0.5D0) +
     +              K* (1.5D0-L)/ (1.0D0-2.0D0*L)*X** (0.5D0-L)
              HU = K*X** (L-0.5D0)
              HV = K**2/ (1.0D0-2.0D0*L)*X** (0.5D0-L)
          END IF
      ELSE IF (H.LT.-0.25D0) THEN
          L2 = - (H+0.25D0)
          L = SQRT(L2)
          CC = COS(L*LOG(X))
          SS = SIN(L*LOG(X))
          SQ = SQRT(X)
          U = SQ* ((1.D0-0.25D0*K*X/H)*CC+0.5D0*K*L*X*SS)
          V = SQ* ((1.D0-0.25D0*K*X/H)*SS+0.5D0*K*L*X*CC)
          PUP = (0.5D0*CC-L*SS)/SQ - 0.5D0*K*SQ* ((0.5D0-H)*CC+L*SS)/H
          PVP = (0.5D0*SS+L*CC)/SQ + 0.5D0*K*SQ* ((H-0.5D0)*SS+L*CC)/H
          HU = 0.5D0*K*K*SQ* ((H+0.5D0)*CC+L*SS)/ (H*H) + U
          HV = 0.5D0*K*K*SQ* ((H+0.5D0)*SS-L*CC)/ (H*H) + V
      ELSE IF (H.EQ.-0.25D0) THEN
          SQ = SQRT(X)
          U = SQ + K*X*SQ
          V = 2.0D0*SQ + (SQ+K*X*SQ)*LOG(X)
          PUP = 0.5D0* (1.0D0/SQ+3.D0*K*SQ)
          PVP = 2.0D0/SQ + K*SQ + 0.5D0* (1.0D0/SQ+3.0D0*K*SQ)*LOG(X)
          HU = K*K*SQ
          HV = K*K*SQ*LOG(X)
      END IF
      RETURN
C
   14 CONTINUE
      RETURN
C
   15 CONTINUE
      RETURN
C
   16 CONTINUE
      IF (X.LT.0.0D0) THEN
          IF (BETA.GT.-1.0D0 .AND. BETA.LT.0.0D0) THEN
              U = (1.0D0+X)** (-BETA)
              V = 1.0D0
              IF (ALPHA.NE.-1.0D0) THEN
                  PUP = -BETA* (1.0D0-X)** (ALPHA+1.0D0)
                  HU = -BETA* (ALPHA+1.0D0)* (1.0D0-X)**ALPHA
              ELSE
                  PUP = -BETA
                  HU = 0.0D0
              END IF
              PVP = 0.0D0
              HV = 0.0D0
          ELSE IF (BETA.EQ.0.0D0) THEN
              U = 1.0D0
              V = LOG((1.0D0+X)/ (1.0D0-X))
              PUP = 0.0D0
              PVP = 2.0D0* (1.0D0-X)**ALPHA
              HU = 0.0D0
              HV = 2.0D0*ALPHA* (1.0D0-X)** (ALPHA-1.0D0)
          ELSE IF (BETA.GT.0.0D0 .AND. BETA.LT.1.0D0) THEN
              U = 1.0D0
              V = (1.0D0+X)** (-BETA)
              PUP = 0.0D0
              IF (ALPHA.NE.-1.0D0) THEN
                  PVP = -BETA* (1.0D0-X)** (ALPHA+1.0D0)
                  HV = -BETA* (ALPHA+1.0D0)* (1.0D0-X)**ALPHA
              ELSE
                  PVP = -BETA
                  HV = 0.0D0
              END IF
              HU = 0.0D0
          END IF
      ELSE IF (X.GE.0.0D0) THEN
          IF (ALPHA.GT.-1.0D0 .AND. ALPHA.LT.0.0D0) THEN
              U = (1.0D0-X)** (-ALPHA)
              V = 1.0D0
              IF (BETA.NE.-1.0D0) THEN
                  PUP = ALPHA* (1.0D0+X)** (BETA+1.0D0)
                  HU = -ALPHA* (BETA+1.0D0)* (1.0D0+X)**BETA
              ELSE
                  PUP = ALPHA
                  HU = 0.0D0
              END IF
              PVP = 0.0D0
              HV = 0.0D0
          ELSE IF (ALPHA.EQ.0.0D0) THEN
              U = 1.0D0
              V = LOG((1.0D0+X)/ (1.0D0-X))
              PUP = 0.0D0
              PVP = 2.0D0* (1.0D0+X)**BETA
              HU = 0.0D0
              HV = -2.0D0*BETA* (1.0D0+X)** (BETA-1.0D0)
          ELSE IF (ALPHA.GT.0.0D0 .AND. ALPHA.LT.1.0D0) THEN
              U = 1.0D0
              V = (1.0D0-X)** (-ALPHA)
              PUP = 0.0D0
              HU = 0.0D0
              IF (BETA.NE.-1.0D0) THEN
                  PVP = ALPHA* (1.0D0+X)** (BETA+1.0D0)
                  HV = -ALPHA* (BETA+1.0D0)* (1.0D0+X)**BETA
              ELSE
                  PVP = ALPHA
                  HV = 0.0D0
              END IF
          END IF
      END IF
      RETURN
C
   17 CONTINUE
      RETURN
C
   18 CONTINUE
      IF (X.LT.0.0D0) THEN
          IF (ALPHA.EQ.0.0D0) THEN
              U = 1.0D0
              V = 0.5D0*LOG((1.0D0+X)/ (1.0D0-X))
              PUP = 0.0D0
              PVP = 1.0D0
              HU = Q(X)
              HV = Q(X)*V
          ELSE IF (ALPHA.GT.0.0D0 .AND. ALPHA.LT.0.5D0) THEN
              U = (1.0D0+X)**ALPHA
              V = (1.0D0+X)** (-ALPHA)
              PUP = ALPHA* (1.0D0-X)*U
              PVP = -ALPHA* (1.0D0-X)*V
              HU = ALPHA* (ALPHA+1.0D0)*U + 2.0D0*BETA**2*U/ (1.0D0-X)
              HV = ALPHA* (ALPHA-1.0D0)*V + 2.0D0*BETA**2*V/ (1.0D0-X)
          ELSE
          END IF
      ELSE IF (X.GE.0.0D0) THEN
          IF (BETA.EQ.0.0D0) THEN
              U = 1.0D0
              V = 0.5D0*LOG((1.0D0+X)/ (1.0D0-X))
              PUP = 0.0D0
              PVP = 1.0D0
              HU = Q(X)
              HV = Q(X)*V
          ELSE IF (BETA.GT.0.0D0 .AND. BETA.LT.0.5D0) THEN
              U = (1.0D0-X)**BETA
              V = (1.0D0-X)** (-BETA)
              PUP = -BETA* (1.0D0+X)*U
              PVP = BETA* (1.0D0+X)*V
              HU = BETA* (BETA+1.0D0)*U + 2.0D0*ALPHA**2*U/ (1.0D0+X)
              HV = BETA* (BETA-1.0D0)*V + 2.0D0*ALPHA**2*V/ (1.0D0+X)
          ELSE
          END IF
      END IF
      RETURN
C
   19 CONTINUE
      IF (X.LT.0.0D0) THEN
          IF (GAMMA.EQ.0.0D0) THEN
              U = 1.0D0
              V = 0.5D0*LOG((1.0D0+X)/ (1.0D0-X))
              PUP = 0.0D0
              PVP = 1.0D0
              HU = Q(X)*U
              HV = Q(X)*V
          ELSE
              U = COS(GAMMA*LOG(1.0D0+X))
              V = SIN(GAMMA*LOG(1.0D0+X))
              PUP = -GAMMA* (1.0D0-X)*V
              PVP = GAMMA* (1.0D0-X)*U
              HU = -GAMMA**2*U - GAMMA*V + 2.0D0*BETA**2*U/ (1.0D0-X)
              HV = -GAMMA**2*V + GAMMA*U + 2.0D0*BETA**2*V/ (1.0D0-X)
          END IF
      ELSE IF (X.GE.0.0D0) THEN
          IF (BETA.EQ.0.0D0) THEN
              U = 1.0D0
              V = 0.5D0*LOG((1.0D0+X)/ (1.0D0-X))
              PUP = 0.0D0
              PVP = 1.0D0
              HU = Q(X)*U
              HV = Q(X)*V
          ELSE IF (BETA.GT.0.0D0 .AND. BETA.LT.0.5D0) THEN
              U = (1.0D0-X)**BETA
              V = (1.0D0-X)** (-BETA)
              PUP = -BETA* (1.0D0+X)*U
              PVP = BETA* (1.0D0+X)*V
              HU = BETA* (BETA+1.0D0)*U - 2.0D0*GAMMA**2*U/ (1.0D0+X)
              HV = BETA* (BETA-1.0D0)*V - 2.0D0*GAMMA**2*V/ (1.0D0+X)
          ELSE
          END IF
      END IF
      RETURN
C
   20 CONTINUE
      U = SQRT(X)*COS(K*LOG(X))
      V = SQRT(X)*SIN(K*LOG(X))
      PUP = 0.5D0*U/X - K*V/X
      PVP = 0.5D0*V/X + K*U/X
      HU = U
      HV = V
      RETURN
C
   21 CONTINUE
      RETURN
C
   22 CONTINUE
      EE = EXP(-X)
      IF (ALPHA.GT.-1.0D0 .AND. ALPHA.LT.0.0D0) THEN
          U = X** (-ALPHA)
          V = 1.0D0
          PUP = -ALPHA*EE
          PVP = 0.0D0
          HU = -ALPHA*EE
          HV = 0.0D0
      ELSE IF (ALPHA.EQ.0.0D0) THEN
          U = 1.0D0
          V = LOG(X)
          PUP = 0.0D0
          PVP = EE
          HU = 0.0D0
          HV = EE
      ELSE IF (ALPHA.GT.0.0D0 .AND. ALPHA.LT.1.0D0) THEN
          U = 1.0D0
          V = X** (-ALPHA)
          PUP = 0.0D0
          PVP = -ALPHA*EE
          HU = 0.0D0
          HV = -ALPHA*EE
      ELSE
      END IF
      RETURN
C
   23 CONTINUE
      IF (ALPHA.GT.-1.0D0 .AND. ALPHA.LT.0.0D0) THEN
          IF (ALPHA.NE.-0.5D0) THEN
              U = X** (0.5D0-ALPHA)
              V = X** (0.5D0+ALPHA)
              PUP = (0.5D0-ALPHA)*X** (-0.5D0-ALPHA)
              PVP = (0.5D0+ALPHA)*X** (-0.5D0+ALPHA)
          ELSE
              U = X
              V = 1.0D0
              PUP = 1.0D0
              PVP = 0.0D0
          END IF
          TX = X**2/16.0D0 - (ALPHA+1.0D0)/2.0D0
          HU = TX*U
          HV = TX*V
      ELSE IF (ALPHA.EQ.0.0D0) THEN
          U = SQRT(X)
          V = U*LOG(X)
          PUP = 0.5D0/U
          PVP = (1.0D0+0.5D0*LOG(X))/U
          HU = (X**2/16.0D0-0.5D0)*U
          HV = (X**2/16.0D0-0.5D0)*V
      ELSE IF (ALPHA.GT.0.0D0 .AND. ALPHA.LT.1.0D0) THEN
          IF (ALPHA.NE.0.5D0) THEN
              U = X** (0.5D0+ALPHA)
              V = X** (0.5D0-ALPHA)
              PUP = (0.5D0+ALPHA)*X** (-0.5D0+ALPHA)
              PVP = (0.5D0-ALPHA)*X** (-0.5D0-ALPHA)
          ELSE
              U = X
              V = 1.0D0
              PUP = 1.0D0
              PVP = 0.0D0
          END IF
          TX = X**2/16.0D0 - (ALPHA+1.0D0)/2.0D0
          HU = TX*U
          HV = TX*V
      ELSE
      END IF
      RETURN
C
   24 CONTINUE
      HPI = 2.0D0*ATAN(1.0D0)
      IF (X.GE.0.0D0) THEN
          IF (ALPHA.GT.-1.0D0 .AND. ALPHA.LT.0.0D0) THEN
              U = (HPI-X)** (0.5D0-ALPHA)
              V = (HPI-X)** (0.5D0+ALPHA)
              PUP = - (0.5D0-ALPHA)* (HPI-X)** (-0.5D0-ALPHA)
              PVP = - (0.5D0+ALPHA)* (HPI-X)** (-0.5D0+ALPHA)
              HU = (0.25D0-ALPHA**2)* (HPI-X)** (-1.5D0-ALPHA) + Q(X)*U
              HV = (0.25D0-ALPHA**2)* (HPI-X)** (-1.5D0+ALPHA) + Q(X)*V
          ELSE IF (ALPHA.EQ.0.0D0) THEN
              U = SQRT(HPI-X)
              V = U*LOG(HPI-X)
              PUP = -0.5D0/U
              PVP = - (1.0D0+0.5D0*LOG(HPI-X))/U
              HU = 0.25D0/ ((HPI-X)*U) + Q(X)*U
              HV = 0.25D0*LOG(HPI-X)/ ((HPI-X)*U) + Q(X)*V
          ELSE IF (ALPHA.GT.0.0D0 .AND. ALPHA.LT.1.0D0) THEN
              U = (HPI-X)** (0.5D0+ALPHA)
              V = (HPI-X)** (0.5D0-ALPHA)
              PUP = - (0.5D0+ALPHA)* (HPI-X)** (-0.5D0+ALPHA)
              PVP = - (0.5D0-ALPHA)* (HPI-X)** (-0.5D0-ALPHA)
              HU = (0.25D0-ALPHA**2)* (HPI-X)** (-1.5D0+ALPHA) + Q(X)*U
              HV = (0.25D0-ALPHA**2)* (HPI-X)** (-1.5D0-ALPHA) + Q(X)*V
          END IF
      ELSE
          IF (BETA.GT.-1.0D0 .AND. BETA.LT.0.0D0) THEN
              U = (HPI+X)** (0.5D0-BETA)
              V = (HPI+X)** (0.5D0+BETA)
              PUP = (0.5D0-BETA)* (HPI+X)** (-0.5D0-BETA)
              PVP = (0.5D0+BETA)* (HPI+X)** (-0.5D0+BETA)
              HU = (0.25D0-BETA**2)* (HPI+X)** (-1.5D0-BETA) + Q(X)*U
              HV = (0.25D0-BETA**2)* (HPI+X)** (-1.5D0+BETA) + Q(X)*V
          ELSE IF (BETA.EQ.0.0D0) THEN
              U = SQRT(HPI+X)
              V = U*LOG(HPI+X)
              PUP = 0.5D0/U
              PVP = (1.0D0+0.5D0*LOG(HPI+X))/U
              HU = 0.25D0/ ((HPI+X)*U) + Q(X)*U
              HV = 0.25D0*LOG(HPI+X)/ ((HPI+X)*U) + Q(X)*V
          ELSE IF (BETA.GT.0.0D0 .AND. BETA.LT.1.0D0) THEN
              U = (HPI+X)** (0.5D0+BETA)
              V = (HPI+X)** (0.5D0-BETA)
              PUP = (0.5D0+BETA)* (HPI+X)** (-0.5D0+BETA)
              PVP = (0.5D0-BETA)* (HPI+X)** (-0.5D0-BETA)
              HU = (0.25D0-BETA**2)* (HPI+X)** (-1.5D0+BETA) + Q(X)*U
              HV = (0.25D0-BETA**2)* (HPI+X)** (-1.5D0-BETA) + Q(X)*V
          END IF
      END IF
      RETURN
C
   25 CONTINUE
      RETURN
C
   26 CONTINUE
      RETURN
C
   27 CONTINUE
      RETURN
C
   28 CONTINUE
      RETURN
C
   29 CONTINUE
      RETURN
C
   30 CONTINUE
      RETURN
C
   31 CONTINUE
      RETURN
C
   32 CONTINUE
      U = 1.0D0
      PUP = 0.0D0
      HU = Q(X)
      IF (X.LE.0.5D0) THEN
          IF (C.EQ.1.0D0) THEN
              V = LOG(X)
              PVP = (1.0D0-X)* (X+S)*W(X)
              HV = D* (X+S)*W(X) - E* (1.0D0-X)*W(X) + V*Q(X)
          ELSE
              V = X** (1.0D0-C)
              PVP = (1.0D0-C)* (1.0D0-X)**D* (X+S)**E
              HV = (1.0D0-C)* (D* (1.0D0-X)** (D-1.0D0)* (X+S)**E-
     +             E* (1.0D0-X)**D* (X+S)** (E-1.0D0)) + V*Q(X)
          END IF
      ELSE
          IF (D.EQ.1.0D0) THEN
              V = LOG(1.0D0-X)
              PVP = -X* (X+S)*W(X)
              HV = C* (X+S)*W(X) + E*X*W(X) + V*Q(X)
          ELSE
              V = (1.0D0-X)** (1.0D0-D)
              PVP = - (1.0D0-D)*X**C* (X+S)**E
              HV = (1.0D0-D)* (C*X** (C-1.0D0)* (X+S)**E+
     +             E*X**C* (X+S)** (E-1.0D0)) + V*Q(X)
          END IF
      END IF
      RETURN
      END
