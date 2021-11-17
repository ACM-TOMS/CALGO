C     OCTOBER 15, 1995; P.B. BAILEY, W.N. EVERITT, B. GARBOW AND A. ZETTL
C
      SUBROUTINE EXAMP()
      CHARACTER ANS
      INTEGER MUMBER,NUMBER
      DOUBLE PRECISION NU,H,K,L,ALPHA,BETA,GAMMA
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA
C
C     THIS SUBROUTINE CONTAINS A SELECTION OF COEFFICIENT FUNCTIONS
C     p,q,w (AND POSSIBLY SUITABLE FUNCTIONS u,v TO DETERMINE SINGULAR
C     B.C.) WHICH DEFINE SOME INTERESTING STURM-LIOUVILLE BOUNDARY
C     VALUE PROBLEMS. IT CAN BE CALLED BY THE MAIN PROGRAM, DRIVE.
C
      WRITE(*,*) ' Here is a collection of 29 differential equations '
      WRITE(*,*) ' which can be used with SLEIGN2.  By typing an '
      WRITE(*,*) ' integer from 1 to 29, one of these differential '
      WRITE(*,*) ' equations is selected, whereupon its coefficient '
      WRITE(*,*) ' functions p,q,w will be displayed along with a '
      WRITE(*,*) ' brief description of its singular points.  The'
      WRITE(*,*) ' endpoints a, b of the interval over which the '
      WRITE(*,*) ' differential equation is integrated are specified '
      WRITE(*,*) ' later; any interval which does not contain '
      WRITE(*,*) ' singular points in its interior is acceptable. '
      WRITE(*,*)
      WRITE(*,*) ' DO YOU WISH TO CONTINUE ? (Y/N) '
      READ(*,9010) ANS
      IF (.NOT.(ANS.EQ.'y' .OR. ANS.EQ.'Y')) STOP
   30 CONTINUE
         WRITE(*,*)
         WRITE(*,*) '  1 IS THE LEGENDRE EQUATION '
         WRITE(*,*) '  2 IS THE BESSEL EQUATION '
         WRITE(*,*) '  3 IS THE HALVORSEN EQUATION '
         WRITE(*,*) '  4 IS THE BOYD EQUATION '
         WRITE(*,*) '  5 IS THE REGULARIZED BOYD EQUATION '
         WRITE(*,*) '  6 IS THE SEARS-TITCHMARSH EQUATION '
         WRITE(*,*) '  7 IS THE BEZ EQUATION '
         WRITE(*,*) '  8 IS THE LAPLACE TIDAL WAVE EQUATION '
         WRITE(*,*) '  9 IS THE LATZKO EQUATION '
         WRITE(*,*) ' 10 IS A WEAKLY REGULAR EQUATION '
         WRITE(*,*) ' 11 IS THE PLUM EQUATION '
         WRITE(*,*) ' 12 IS THE MATHIEU PERIODIC EQUATION '
         WRITE(*,*) ' 13 IS THE HYDROGEN ATOM EQUATION '
         WRITE(*,*) ' 14 IS THE MARLETTA EQUATION '
         WRITE(*,*) ' 15 IS THE HARMONIC OSCILLATOR EQUATION '
      WRITE(*,*)
      WRITE(*,*) ' Press any key to continue. '
      READ(*,9010) ANS
 9010 FORMAT(A1)
         WRITE(*,*) ' 16 IS THE JACOBI EQUATION '
         WRITE(*,*) ' 17 IS THE ROTATION MORSE OSCILLATOR EQUATION '
         WRITE(*,*) ' 18 IS THE DUNSCH EQUATION '
         WRITE(*,*) ' 19 IS THE DONSCH EQUATION '
         WRITE(*,*) ' 20 IS THE KRALL EQUATION '
         WRITE(*,*) ' 21 IS THE FOURIER EQUATION '
         WRITE(*,*) ' 22 IS THE LAGUERRE EQUATION '
         WRITE(*,*) ' 23 IS THE LAGUERRE/LIOUVILLE FORM EQUATION '
         WRITE(*,*) ' 24 IS THE JACOBI/LIOUVILLE FORM EQUATION '
         WRITE(*,*) ' 25 IS THE MEISSNER EQUATION '
         WRITE(*,*) ' 26 IS THE LOHNER EQUATION '
         WRITE(*,*) ' 27 IS THE JOERGENS EQUATION '
         WRITE(*,*) ' 28 IS THE BEHNKE-GOERISCH EQUATION '
         WRITE(*,*) ' 29 IS THE WHITTAKER EQUATION '
         WRITE(*,*)
         WRITE(*,*) ' ENTER THE NUMBER OF YOUR CHOICE: '
         READ(*,*) NUMBER
         IF (NUMBER.LT.1 .OR. NUMBER.GT.29) GO TO 30
      MUMBER = NUMBER
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     1       16,17,18,19,20,21,22,23,24,25,26,27,28,29), NUMBER
C
    1 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1 - x*x, q = 1/4, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT CIRCLE, NON-OSCILLATORY AT -1. '
      WRITE(*,*) ' LIMIT CIRCLE, NON-OSCILLATORY AT +1. '
      GO TO 40
C
    2 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = (nu*nu-0.25)/x*x, w = 1 '
      WRITE(*,*) '            (nu a parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT FOR ALL nu AT + INFINITY'
      WRITE(*,*) ' AT 0: LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '       FOR -1.LT.nu.LT.1 BUT nu*nu.NE.0.25. '
      WRITE(*,*) '       REGULAR FOR nu*nu = 0.25. '
      WRITE(*,*) '       LIMIT POINT FOR nu*nu.GE.1.0. '
      MUMBER = 101
      GO TO 40
C
    3 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 0, w = exp(-2/x)/x**4 '
      WRITE(*,*)
      WRITE(*,*) ' WEAKLY REGULAR AT 0. '
      WRITE(*,*) ' LIMIT CIRCLE, NON-OSCILLATORY AT +INFINITY. '
      GO TO 40
C
    4 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,0) '
      WRITE(*,*) '                            AND on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = -1/x, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT CIRCLE, NON-OSCILLATORY AT 0+ AND 0-. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
    5 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,0) '
      WRITE(*,*) '                            AND on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = r*r, q = -r*r*(ln|x|)**2, w = r*r '
      WRITE(*,*) '      where r = exp(-(x*ln(|x|)-x)) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' WEAKLY REGULAR AT 0+ AND 0-. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
    6 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = x, q = -x, w = 1/x '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT 0. '
      WRITE(*,*) ' LIMIT CIRCLE, OSCILLATORY AT +INFINITY. '
      GO TO 40
C
    7 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity,0) '
      WRITE(*,*) '                            AND on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = x, q = -1/x, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT CIRCLE, OSCILLATORY AT 0+ AND 0-. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
    8 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1/x, q = (k/x**2) + (k**2/x), w = 1 '
      WRITE(*,*) '         (k a non-zero parameter) '
      WRITE(*,*)
      WRITE(*,*) 'AT 0: LIMIT CIRCLE, NON-OSCILLATORY FOR ALL k. '
      WRITE(*,*) 'AT +INFINITY: LIMIT POINT AT FOR ALL k. '
      MUMBER = 102
      GO TO 40
C
    9 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,1) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1 - x**7, q = 0, w = x**7 '
      WRITE(*,*)
      WRITE(*,*) ' WEAKLY REGULAR AT 0. '
      WRITE(*,*) ' LIMIT CIRCLE, NON-OSCILLATORY AT +1. '
      GO TO 40
C
   10 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = sqrt(x), q = 0, w = 1/sqrt(x) '
      WRITE(*,*)
      WRITE(*,*) ' WEAKLY REGULAR AT 0. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   11 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 100*cos(x)**2, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   12 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 2*k*cos(2x), w = 1 '
      WRITE(*,*) '       (k a non-zero parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 103
      GO TO 40
C
   13 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = k/x + h/x**2, w = 1 '
      WRITE(*,*) '        q = k/x + h/x**2 + 1  if h.lt.-0.25  '
      WRITE(*,*) '        (h,k parameters) '
      WRITE(*,*)
      WRITE(*,*) 'LIMIT POINT FOR ALL h,k AT +INFINITY '
      WRITE(*,*) 'AT 0: '
      WRITE(*,*) '     REGULAR for h = k = 0. '
      WRITE(*,*) '     LIMIT-CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '       FOR h = 0 AND ALL k.NE.0. '
      WRITE(*,*) '     LIMIT-CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '       FOR -0.25.LE.h.LT.0.75 BUT h.NE.0, AND ALL k. '
      WRITE(*,*) '     LIMIT-CIRCLE, OSCILLATORY '
      WRITE(*,*) '       FOR h.LT.-0.25 AND ALL k. '
      WRITE(*,*) '      (Here, 1 has been added to q so that '
      WRITE(*,*) '       at least some eigenvalues are positive.) '
      WRITE(*,*) '     LIMIT POINT FOR h.GE.0.75 AND ALL k. '
      MUMBER = 104
      GO TO 40
C
   14 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 3.0*(X-31.0)/(4.0*(X+1.0)*(4.0+X)**2), '
      WRITE(*,*) ' w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' REGULAR AT 0. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   15 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = x*x, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   16 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE(*,*)
      WRITE(*,*) ' p = (1-x)**(alpha+1)*(1+x)**(beta+1), '
      WRITE(*,*) ' q = 0, w = (1-x)**alpha*(1+x)**beta '
      WRITE(*,*) '            (alpha, beta parameters) '
      WRITE(*,*)
      WRITE(*,*) ' AT -1.0: '
      WRITE(*,*) '         LIMIT POINT FOR beta.LE.-1. '
      WRITE(*,*) '         WEAKLY REGULAR FOR -1.LT.beta.LT.0. '
      WRITE(*,*) '         LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '         FOR 0.LE.beta.LT.1. '
      WRITE(*,*) '         LIMIT POINT FOR beta.GE.1. '
      WRITE(*,*) ' AT +1.0: '
      WRITE(*,*) '         LIMIT POINT FOR alpha.LE.-1. '
      WRITE(*,*) '         WEAKLY REGULAR FOR -1.LT.alpha.LT.0. '
      WRITE(*,*) '         LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '         FOR 0.LE.alpha.LT.1. '
      WRITE(*,*) '         LIMIT POINT FOR alpha.GE.1. '
      MUMBER = 105
      GO TO 40
C
   17 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 2/x**2 - 2000(2e-e*e), w = 1 '
      WRITE(*,*) '       where e = exp(-1.7(x-1.3)) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT 0. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   18 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1 - x*x, q = 2*alpha**2/(1+x) + 2*beta**2/(1-x),'
      WRITE(*,*) ' w = 1        (alpha, beta non-negative parameters) '
      WRITE(*,*)
      WRITE(*,*) 'AT -1.0: '
      WRITE(*,*) '        LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '        FOR 0.LE.alpha.LT.0.5. '
      WRITE(*,*) '        LIMIT POINT FOR alpha.GE.0.5. '
      WRITE(*,*) 'AT +1.0: '
      WRITE(*,*) '        LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '        FOR 0.LE.beta.LT.0.5. '
      WRITE(*,*) '        LIMIT POINT FOR beta.GE.0.5. '
      MUMBER = 106
      GO TO 40
C
   19 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-1,1) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1 - x*x, q = -2*gamma**2/(1+x) +2*beta**2/(1-x),'
      WRITE(*,*) ' w = 1            (gamma, beta parameters) '
      WRITE(*,*)
      WRITE(*,*) 'AT -1.0: '
      WRITE(*,*) '        LIMIT CIRCLE, NON-OSCILLATORY FOR gamma = 0.'
      WRITE(*,*) '        LIMIT CIRCLE, OSCILLATORY FOR gamma.GT.0. '
      WRITE(*,*) 'AT +1.0: '
      WRITE(*,*) '        LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '        FOR 0.LE.beta.LT.0.5. '
      WRITE(*,*) '        LIMIT POINT FOR beta.GE.0.5. '
      MUMBER = 107
      GO TO 40
C
   20 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 1 - (k**2+0.25)/x**2, w = 1 '
      WRITE(*,*) '        (k a positive parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT CIRCLE, OSCILLATORY AT 0. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 108
      GO TO 40
C
   21 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 0, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   22 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = x**(alpha+1)*exp(-x), w = x**alpha*exp(-x) '
      WRITE(*,*) ' q = 0        (alpha a parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT FOR ALL alpha AT +INFINITY'
      WRITE(*,*) 'AT 0: '
      WRITE(*,*) '     LIMIT POINT FOR alpha.LE.-1. '
      WRITE(*,*) '     WEAKLY REGULAR FOR -1.LT.alpha.LT.0. '
      WRITE(*,*) '     LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '       FOR 0.LE.alpha.LT.1. '
      WRITE(*,*) '     LIMIT POINT FOR alpha.GE.1. '
      MUMBER = 109
      GO TO 40
C
   23 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, w = 1 '
      WRITE(*,*) ' q = (alpha**2-0.25)/x**2 - (alpha+1)/2 + x**2/16'
      WRITE(*,*) '              (alpha a parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT FOR ALL alpha AT +INFINITY'
      WRITE(*,*) 'AT 0: '
      WRITE(*,*) '     LIMIT POINT FOR alpha.LE.-1. '
      WRITE(*,*) '     LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '       FOR -1.LT.alpha.LT.1 BUT alpha**2.NE.0.25. '
      WRITE(*,*) '    REGULAR FOR alpha**2 = 0.25. '
      WRITE(*,*) '    LIMIT POINT FOR alpha.GE.1. '
      MUMBER = 110
      GO TO 40
C
   24 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-pi/2,pi/2) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, w = 1 '
      WRITE(*,*) ' q = (beta**2-0.25)/(4*tan((x+pi/2)/2)**2)+ '
      WRITE(*,*) '     (alpha**2-0.25)/(4*tan((x-pi/2)/2)**2)- '
      WRITE(*,*) '     (4*alpha*beta+4*alpha+4*beta+3)/8 '
      WRITE(*,*) '            (alpha, beta parameters) '
      WRITE(*,*)
      WRITE(*,*) 'AT -pi/2: '
      WRITE(*,*) '         LIMIT POINT FOR beta.LE.-1. '
      WRITE(*,*) '         LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '         FOR -1.LT.beta.LT.1 BUT beta**2.NE.0.25.'
      WRITE(*,*) '         REGULAR FOR beta**2 = 0.25. '
      WRITE(*,*) '         LIMIT POINT FOR beta.GE.1. '
      WRITE(*,*) 'AT +pi/2: '
      WRITE(*,*) '         LIMIT POINT FOR alpha.LE.-1. '
      WRITE(*,*) '         LIMIT CIRCLE, NON-OSCILLATORY '
      WRITE(*,*) '         FOR -1.LT.alpha.LT.1 BUT alpha**2.NE.0.25.'
      WRITE(*,*) '         REGULAR FOR alpha**2 = 0.25. '
      WRITE(*,*) '         LIMIT POINT FOR alpha.GE.1. '
      MUMBER = 111
      GO TO 40
C
   25 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 0 '
      WRITE(*,*) ' w = 1 when x.le.0. '
      WRITE(*,*) '   = 9 when x.gt.0. '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   26 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = -1000*x, w = 1 '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      GO TO 40
C
   27 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 0.25*exp(2*x) - k*exp(x), w = 1 '
      WRITE(*,*) '                (k a parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 112
      GO TO 40
C
   28 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (-infinity, '
      WRITE(*,*) '                                    +infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = k*cos(x)**2, w = 1 '
      WRITE(*,*) '          (k a parameter) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT -INFINITY. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 113
      GO TO 40
C
   29 CONTINUE
      WRITE(*,*) '    -(p*y'')'' + q*y = lambda*w*y on (0,+infinity) '
      WRITE(*,*)
      WRITE(*,*) ' p = 1, q = 0.25 + (k**2-1)/(4*x**2)), w = 1/x '
      WRITE(*,*) '       (k a positive parameter .GE. 1) '
      WRITE(*,*)
      WRITE(*,*) ' LIMIT POINT AT 0. '
      WRITE(*,*) ' LIMIT POINT AT +INFINITY. '
      MUMBER = 114
      GO TO 40
C
   40 CONTINUE
      WRITE(*,*)
      WRITE(*,*) ' IS THIS THE CORRECT DIFFERENTIAL EQUATION ? (Y/N) '
      READ(*,9010) ANS
      IF (.NOT.(ANS.EQ.'y' .OR. ANS.EQ.'Y')) GO TO 30
      IF (MUMBER.EQ.NUMBER) RETURN
C
C     Now enter any parameters needed for these D.E.'s, or defaults.
C
      MUMBER = MUMBER - 100
      GO TO (41,42,43,44,45,46,47,48,49,50,51,52,53,54), MUMBER
C
   41 CONTINUE
      NU = 1.0
      WRITE(*,*) ' Choose real parameter nu, nu = '
      READ(*,*) NU
      RETURN
C
   42 CONTINUE
      K = 1.0
      WRITE(*,*) ' Choose real parameter k.ne.0., k = '
      READ(*,*) K
      RETURN
C
   43 CONTINUE
      K = 1.0
      WRITE(*,*) ' Choose real parameter k = '
      READ(*,*) K
      RETURN
C
   44 CONTINUE
      H = 1.0
      K = 1.0
      WRITE(*,*) ' Choose real parameters h,k = '
      READ(*,*) H,K
      RETURN
C
   45 CONTINUE
      BETA = 0.1
      ALPHA = 0.1
      WRITE(*,*) ' Choose real parameters beta,alpha = '
      READ(*,*) BETA,ALPHA
      RETURN
C
   46 CONTINUE
      ALPHA = 0.1
      BETA = 0.1
      WRITE(*,*) ' Choose real parameters alpha,beta = '
      READ(*,*) ALPHA,BETA
      RETURN
C
   47 CONTINUE
      GAMMA = 0.1
      BETA = 0.1
      WRITE(*,*) ' Choose real parameters gamma,beta = '
      READ(*,*) GAMMA,BETA
      RETURN
C
   48 CONTINUE
      K = 1.0
      WRITE(*,*) ' Choose real parameter k.gt.0., k = '
      READ(*,*) K
      RETURN
C
   49 CONTINUE
      ALPHA = 0.1
      WRITE(*,*) ' Choose real parameter alpha = '
      READ(*,*) ALPHA
      RETURN
C
   50 CONTINUE
      ALPHA = 0.1
      WRITE(*,*) ' Choose real parameter alpha = '
      READ(*,*) ALPHA
      RETURN
C
   51 CONTINUE
      BETA = 0.1
      ALPHA = 0.1
      WRITE(*,*) ' Choose real parameters beta,alpha = '
      READ(*,*) BETA,ALPHA
      RETURN
C
   52 CONTINUE
      K = 0.1
      WRITE(*,*) ' Choose real parameter k = '
      READ(*,*) K
      RETURN
C
   53 CONTINUE
      K = 0.1
      WRITE(*,*) ' Choose real parameter k = '
      READ(*,*) K
      RETURN
C
   54 CONTINUE
      K = 0.1
      WRITE(*,*) ' Choose real parameter k = '
      READ(*,*) K
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION P(X)
      DOUBLE PRECISION X
      INTEGER NUMBER
      DOUBLE PRECISION NU,H,K,L,ALPHA,BETA,GAMMA
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     1       16,17,18,19,20,21,22,23,24,25,26,27,28,29), NUMBER
C
    1 CONTINUE
      P = 1.0 - X*X
      RETURN
C
    2 CONTINUE
      P = 1.0
      RETURN
C
    3 CONTINUE
      P = 1.0
      RETURN
C
    4 CONTINUE
      P = 1.0
      RETURN
C
    5 CONTINUE
      P = EXP(-2.0*(X*LOG(ABS(X))-X))
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
      P = 1.0/X
      RETURN
C
    9 CONTINUE
      P = 1.0 - X**7
      RETURN
C
   10 CONTINUE
      P = SQRT(X)
      RETURN
C
   11 CONTINUE
      P = 1.0
      RETURN
C
   12 CONTINUE
      P = 1.0
      RETURN
C
   13 CONTINUE
      P = 1.0
      RETURN
C
   14 CONTINUE
      P = 1.0
      RETURN
C
   15 CONTINUE
      P = 1.0
      RETURN
C
   16 CONTINUE
      IF (ALPHA.NE.-1.0 .AND. BETA.NE.-1.0) THEN
         P = (1.0-X)**(ALPHA+1.0)*(1.0+X)**(BETA+1.0)
      ELSE IF (ALPHA.NE.-1.0 .AND. BETA.EQ.-1.0) THEN
         P = (1.0-X)**(ALPHA+1.0)
      ELSE IF (ALPHA.EQ.-1.0 .AND. BETA.NE.-1.0) THEN
         P = (1.0+X)**(BETA+1.0)
      ELSE
         P = 1.0
         END IF
      RETURN
C
   17 CONTINUE
      P = 1.0
      RETURN
C
   18 CONTINUE
      P = 1.0 - X*X
      RETURN
C
   19 CONTINUE
      P = 1.0 - X*X
      RETURN
C
   20 CONTINUE
      P = 1.0
      RETURN
C
   21 CONTINUE
      P = 1.0
      RETURN
C
   22 CONTINUE
      P = EXP(-X)
      IF (ALPHA.NE.-1.0) P = P*X**(ALPHA+1.0)
      RETURN
C
   23 CONTINUE
      P = 1.0
      RETURN
C
   24 CONTINUE
      P = 1.0
      RETURN
C
   25 CONTINUE
      P = 1.0
      RETURN
C
   26 CONTINUE
      P = 1.0
      RETURN
C
   27 CONTINUE
      P = 1.0
      RETURN
C
   28 CONTINUE
      P = 1.0
      RETURN
C
   29 CONTINUE
      P = 1.0
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION Q(X)
      DOUBLE PRECISION X
      INTEGER NUMBER
      DOUBLE PRECISION NU,H,K,L,ALPHA,BETA,GAMMA,E,HPI
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA
c
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     1       16,17,18,19,20,21,22,23,24,25,26,27,28,29), NUMBER
C
    1 CONTINUE
      Q = 0.25
      RETURN
C
    2 CONTINUE
      Q = 0.0
      IF (NU.NE.-0.5 .AND. NU.NE.0.5) Q = (NU*NU-0.25)/X**2
      RETURN
C
    3 CONTINUE
      Q = 0.0
      RETURN
C
    4 CONTINUE
      Q = -1.0/X
      RETURN
C
    5 CONTINUE
      Q = -EXP(-2.0*(X*LOG(ABS(X))-X))*LOG(ABS(X))**2
      RETURN
C
    6 CONTINUE
      Q = -X
      RETURN
C
    7 CONTINUE
      Q = -1.0/X
      RETURN
C
    8 CONTINUE
      Q = K/X**2 + K**2/X
      RETURN
C
    9 CONTINUE
      Q = 0.0
      RETURN
C
   10 CONTINUE
      Q = 0.0
      RETURN
C
   11 CONTINUE
      Q = 100.0*COS(X)**2
      RETURN
C
   12 CONTINUE
      Q = 2.0*K*COS(2.0*X)
      RETURN
C
   13 CONTINUE
      Q = K/X + H/(X*X)
      IF (H .LT. -0.25) Q = Q + 1.0
      RETURN
C
   14 CONTINUE
      Q = 3.0*(X-31.0)/(4.0*(X+1.0)*(X+4.0)**2)
      RETURN
C
   15 CONTINUE
      Q = X*X
      RETURN
C
   16 CONTINUE
      Q = 0.0
      RETURN
C
   17 CONTINUE
      L = 1.0
      E = EXP(-1.7*(X-1.3))
      Q = L*(L+1.0)/X**2 - 2000.0*E*(2.0-E)
      RETURN
C
   18 CONTINUE
      IF (ALPHA.NE.0.0 .AND. BETA.NE.0.0) THEN
         Q = 2.0*ALPHA**2/(1.0+X) + 2.0*BETA**2/(1.0-X)
      ELSE IF (ALPHA.EQ.0.0 .AND. BETA.NE.0.0) THEN
         Q = 2.0*BETA**2/(1.0-X)
      ELSE IF (ALPHA.NE.0.0 .AND. BETA.EQ.0.0) THEN
         Q = 2.0*ALPHA**2/(1.0+X)
      ELSE
         Q = 0.0
         END IF
      RETURN
C
   19 CONTINUE
      Q = -2.0*GAMMA**2/(1.0+X) + 2.0*BETA**2/(1.0-X)
      RETURN
C
   20 CONTINUE
      Q = 1.0 - (K**2+0.25)/X**2
      RETURN
C
   21 CONTINUE
      Q = 0.0
      RETURN
C
   22 CONTINUE
      Q = 0.0
      RETURN
C
   23 CONTINUE
      Q = -(ALPHA+1.0)/2.0 + X**2/16.0
      IF (ALPHA.NE.0.5 .AND. ALPHA.NE.-0.5) Q = Q + (ALPHA**2-0.25)/X**2
      RETURN
C
   24 CONTINUE
      HPI = 2.0*ATAN(1.0)
      IF (BETA*BETA.NE.0.25 .AND. ALPHA*ALPHA.NE.0.25) THEN
         Q = (BETA**2-0.25)/(4.0*TAN((X+HPI)/2.0)**2) +
     1       (ALPHA**2-0.25)/(4.0*TAN((X-HPI)/2.0)**2) -
     2       (ALPHA*BETA+ALPHA+BETA+0.75)/2.0
      ELSE IF (BETA*BETA.EQ.0.25 .AND. ALPHA*ALPHA.NE.0.25) THEN
         Q = (ALPHA**2-0.25)/(4.0*TAN((X-HPI)/2.0)**2) -
     1       (ALPHA*BETA+ALPHA+BETA+0.75)/2.0
      ELSE IF (BETA*BETA.NE.0.25 .AND. ALPHA*ALPHA.EQ.0.25) THEN
         Q = (BETA**2-0.25)/(4.0*TAN((X+HPI)/2.0)**2) -
     1       (ALPHA*BETA+ALPHA+BETA+0.75)/2.0
      ELSE
         Q = -(ALPHA*BETA+ALPHA+BETA+0.75)/2.0
         END IF
      RETURN
C
   25 CONTINUE
      Q = 0.0
      RETURN
C
   26 CONTINUE
      Q = -1000.0*X
      RETURN
C
   27 CONTINUE
      E = EXP(X)
      Q = E*(0.25*E-K)
      RETURN
C
   28 CONTINUE
      Q = K*COS(X)**2
      RETURN
C
   29 CONTINUE
      Q = 0.25 + (K**K-1.0)/(4.0*X**2)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION W(X)
      DOUBLE PRECISION X
      INTEGER NUMBER
      DOUBLE PRECISION NU,H,K,L,ALPHA,BETA,GAMMA
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     1       16,17,18,19,20,21,22,23,24,25,26,27,28,29), NUMBER
C
    1 CONTINUE
      W = 1.0
      RETURN
C
    2 CONTINUE
      W = 1.0
      RETURN
C
    3 CONTINUE
      W = 0.0
      IF (X.NE.0.0) W = EXP(-2.0/X)/X**4
      RETURN
C
    4 CONTINUE
      W = 1.0
      RETURN
C
    5 CONTINUE
      W = EXP(-2.0*(X*LOG(ABS(X))-X))
      RETURN
C
    6 CONTINUE
      W = 1.0/X
      RETURN
C
    7 CONTINUE
      W = 1.0
      RETURN
C
    8 CONTINUE
      W = 1.0
      RETURN
C
    9 CONTINUE
      W = X**7
      RETURN
C
   10 CONTINUE
      W = 1.0/SQRT(X)
      RETURN
C
   11 CONTINUE
      W = 1.0
      RETURN
C
   12 CONTINUE
      W = 1.0
      RETURN
C
   13 CONTINUE
      W = 1.0
      RETURN
C
   14 CONTINUE
      W = 1.0
      RETURN
C
   15 CONTINUE
      W = 1.0
      RETURN
C
   16 CONTINUE
      IF (ALPHA.NE.0.0 .AND. BETA.NE.0.0) THEN
         W = (1.0-X)**ALPHA*(1.0+X)**BETA
      ELSE IF (ALPHA.NE.0.0 .AND. BETA.EQ.0.0) THEN
         W = (1.0-X)**ALPHA
      ELSE IF (ALPHA.EQ.0.0 .AND. BETA.NE.0.0) THEN
         W = (1.0+X)**BETA
      ELSE
         W = 1.0
         END IF
      RETURN
C
   17 CONTINUE
      W = 1.0
      RETURN
C
   18 CONTINUE
      W = 1.0
      RETURN
C
   19 CONTINUE
      W = 1.0
      RETURN
C
   20 CONTINUE
      W = 1.0
      RETURN
C
   21 CONTINUE
      W = 1.0
      RETURN
C
   22 CONTINUE
      W = EXP(-X)
      IF (ALPHA.NE.0.0) W = W*(X**ALPHA)
      RETURN
C
   23 CONTINUE
      W = 1.0
      RETURN
C
   24 CONTINUE
      W = 1.0
      RETURN
C
   25 CONTINUE
      W = 9.0
      IF (X.LE.0.0) W = 1.0
      RETURN
C
   26 CONTINUE
      W = 1.0
      RETURN
C
   27 CONTINUE
      W = 1.0
      RETURN
C
   28 CONTINUE
      W = 1.0
      RETURN
C
   29 CONTINUE
      W = 1.0/X
      RETURN
      END
C
      SUBROUTINE UV(X,U,PUP,V,PVP,HU,HV)
      DOUBLE PRECISION X,U,PUP,V,PVP,HU,HV
      INTEGER NUMBER
      DOUBLE PRECISION NU,H,K,L,ALPHA,BETA,GAMMA,E,TX,HPI,L2,SQ,C,S
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,K,L,ALPHA,BETA,GAMMA
      DOUBLE PRECISION Q
      EXTERNAL Q
C
C     HERE, HU MEANS -(pu')' + qu.
C
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     1       16,17,18,19,20,21,22,23,24,25,26,27,28,29), NUMBER
C
    1 CONTINUE
      U = 1.0
      PUP = 0.0
      V = 0.5*LOG((1.0+X)/(1.0-X))
      PVP = 1.0
      HU = 0.25*U
      HV = 0.25*V
      RETURN
C
    2 CONTINUE
      IF (NU.NE.-0.5 .AND. NU.NE.0.5 .AND. NU.NE.0.0) THEN
         U = X**(NU+0.5)
         PUP = (NU+0.5)*X**(NU-0.5)
         V = X**(-NU+0.5)
         PVP = (-NU+0.5)*X**(-NU-0.5)
      ELSE IF (NU.EQ.-0.5) THEN
         U = X
         PUP = 1.0
         V = 1.0
         PVP = 0.0
      ELSE IF (NU.EQ.0.5) THEN
         U = X
         PUP = 1.0
         V = -1.0
         PVP = 0.0
      ELSE IF (NU.EQ.0.0) THEN
         U = SQRT(X)
         PUP = 0.5/U
         V = U*LOG(X)
         PVP = (0.5*LOG(X)+1.0)/U
         END IF
      HU = 0.0
      HV = 0.0
      RETURN
C
    3 CONTINUE
      U = 1.0
      V = X
      PUP = 0.0
      PVP = 1.0
      HU = 0.0
      HV = 0.0
      RETURN
C
    4 CONTINUE
      TX = LOG(ABS(X))
      U = X
      PUP = 1.0
      V = 1.0 - X*TX
      PVP = -1.0 - TX
      HU = -1.0
      HV = TX
      RETURN
C
    5 CONTINUE
      RETURN
C
    6 CONTINUE
      U = (COS(X)+SIN(X))/SQRT(X)
      V = (COS(X)-SIN(X))/SQRT(X)
      PUP = -0.5*U + X*V
      PVP = -0.5*V - X*U
      HU = -0.25*U/X
      HV = -0.25*V/X
      RETURN
C
    7 CONTINUE
      TX = LOG(ABS(X))
      U = COS(TX)
      V = SIN(TX)
      PUP = -V
      PVP = U
      HU = 0.0
      HV = 0.0
      RETURN
C
    8 CONTINUE
      V = X - 1.0/K
      U = X*X
      PVP = 1.0/X
      PUP = 2.0
      HU = K + X*K**2
      HV = K**2
      RETURN
C
    9 CONTINUE
      U = 1.0
      V = -LOG(1.0-X)
      PUP = 0.0
      PVP = (((((X+1.0)*X+1.0)*X+1.0)*X+1.0)*X+1.0)*X+1.0
      HU = 0.0
      HV = -(((((6.0*X+5.0)*X+4.0)*X+3.0)*X+2.0)*X+1.0)
      RETURN
C
   10 CONTINUE
      U = 2.0*SQRT(X)
      V = 1.0
      PUP = 1.0
      PVP = 0.0
      HU = 0.0
      HV = 0.0
      RETURN
C
   11 CONTINUE
      RETURN
C
   12 CONTINUE
      RETURN
C
   13 CONTINUE
      IF (H .GT. -0.25) THEN
        L = SQRT(H+0.25)
        IF (H.EQ.0.0) THEN
           U = X
           V = 1.0 + K*X*LOG(X)
           PUP = 1.0
           PVP = K*(1.0+LOG(X))
           HU = K
           HV = K*K*LOG(X)
        ELSE
           U = X**(0.5+L)
           V = X**(0.5-L) + (K/(1.0-2.0*L))*X**(1.5-L)
           PUP = (0.5+L)*X**(L-0.5)
           PVP = (0.5-L)*X**(-L-0.5) + K*(1.5-L)/(1.0-2.0*L)*X**(0.5-L)
           HU = K*X**(L-0.5)
           HV = K**2/(1.0-2.0*L)*X**(0.5-L)
           END IF
      ELSE IF (H .LT. -0.25) THEN
        L2 = -(H+0.25)
        L = SQRT(L2)
        C = COS(L*LOG(X))
        S = SIN(L*LOG(X))
        SQ = SQRT(X)
        U = SQ*((1.-0.25*K*X/H)*C + 0.5*K*L*X*S)
        V = SQ*((1.-0.25*K*X/H)*S + 0.5*K*L*X*C)
        PUP = (0.5*C - L*S)/SQ - 0.5*K*SQ*((0.5-H)*C + L*S)/H
        PVP = (0.5*S + L*C)/SQ + 0.5*K*SQ*((H-0.5)*S + L*C)/H
        HU = 0.5*K*K*SQ*((H+0.5)*C + L*S)/(H*H) + U
        HV = 0.5*K*K*SQ*((H+0.5)*S - L*C)/(H*H) + V
      ELSE IF (H .EQ. -0.25) THEN
        SQ = SQRT(X)
        U = SQ + K*X*SQ
        V = 2.0*SQ + (SQ + K*X*SQ)*LOG(X)
        PUP = 0.5*(1.0/SQ + 3.*K*SQ)
        PVP = 2.0/SQ + K*SQ + 0.5*(1.0/SQ + 3.0*K*SQ)*LOG(X)
        HU = K*K*SQ
        HV = K*K*SQ*LOG(X)
      ENDIF
      RETURN
C
   14 CONTINUE
      RETURN
C
   15 CONTINUE
      RETURN
C
   16 CONTINUE
      IF (X.LT.0.0) THEN
         IF (BETA.GT.-1.0 .AND. BETA.LT.0.0) THEN
            U = (1.0+X)**(-BETA)
            V = 1.0
            IF (ALPHA.NE.-1.0) THEN
               PUP = -BETA*(1.0-X)**(ALPHA+1.0)
               HU = -BETA*(ALPHA+1.0)*(1.0-X)**ALPHA
            ELSE
               PUP = -BETA
               HU = 0.0
               END IF
            PVP = 0.0
            HV = 0.0
         ELSE IF (BETA.EQ.0.0) THEN
            U = 1.0
            V =  LOG((1.0+X)/(1.0-X))
            PUP = 0.0
            PVP = 2.0*(1.0-X)**ALPHA
            HU = 0.0
            HV = 2.0*ALPHA*(1.0-X)**(ALPHA-1.0)
         ELSE IF (BETA.GT.0.0 .AND. BETA.LT.1.0) THEN
            U = 1.0
            V = (1.0+X)**(-BETA)
            PUP = 0.0
            IF (ALPHA.NE.-1.0) THEN
               PVP = -BETA*(1.0-X)**(ALPHA+1.0)
               HV = -BETA*(ALPHA+1.0)*(1.0-X)**ALPHA
            ELSE
               PVP = -BETA
               HV = 0.0
               END IF
            HU = 0.0
            END IF
      ELSE IF (X.GE.0.0) THEN
         IF (ALPHA.GT.-1.0 .AND. ALPHA.LT.0.0) THEN
            U = (1.0-X)**(-ALPHA)
            V = 1.0
            IF (BETA.NE.-1.0) THEN
               PUP = ALPHA*(1.0+X)**(BETA+1.0)
               HU = -ALPHA*(BETA+1.0)*(1.0+X)**BETA
            ELSE
               PUP = ALPHA
               HU = 0.0
               END IF
            PVP = 0.0
            HV = 0.0
         ELSE IF (ALPHA.EQ.0.0) THEN
            U = 1.0
            V = LOG((1.0+X)/(1.0-X))
            PUP = 0.0
            PVP = 2.0*(1.0+X)**BETA
            HU = 0.0
            HV = -2.0*BETA*(1.0+X)**(BETA-1.0)
         ELSE IF (ALPHA.GT.0.0 .AND. ALPHA.LT.1.0) THEN
            U = 1.0
            V = (1.0-X)**(-ALPHA)
            PUP = 0.0
            HU = 0.0
            IF (BETA.NE.-1.0) THEN
               PVP = ALPHA*(1.0+X)**(BETA+1.0)
               HV = -ALPHA*(BETA+1.0)*(1.0+X)**BETA
            ELSE
               PVP = ALPHA
               HV = 0.0
               END IF
            END IF
         END IF
      RETURN
C
   17 CONTINUE
      RETURN
C
   18 CONTINUE
      IF (X.LT.0.0) THEN
         IF (ALPHA.EQ.0.0) THEN
            U = 1.0
            V = 0.5*LOG((1.0+X)/(1.0-X))
            PUP = 0.0
            PVP = 1.0
            HU = Q(X)
            HV = Q(X)*V
         ELSE IF (ALPHA.GT.0.0 .AND. ALPHA.LT.0.5) THEN
            U = (1.0+X)**ALPHA
            V = (1.0+X)**(-ALPHA)
            PUP = ALPHA*(1.0-X)*U
            PVP = -ALPHA*(1.0-X)*V
            HU = ALPHA*(ALPHA+1.0)*U + 2.0*BETA**2*U/(1.0-X)
            HV = ALPHA*(ALPHA-1.0)*V + 2.0*BETA**2*V/(1.0-X)
         ELSE
            END IF
      ELSE IF (X.GE.0.0) THEN
         IF (BETA.EQ.0.0) THEN
            U = 1.0
            V = 0.5*LOG((1.0+X)/(1.0-X))
            PUP = 0.0
            PVP = 1.0
            HU = Q(X)
            HV = Q(X)*V
         ELSE IF (BETA.GT.0.0 .AND. BETA.LT.0.5) THEN
            U = (1.0-X)**BETA
            V = (1.0-X)**(-BETA)
            PUP = -BETA*(1.0+X)*U
            PVP = BETA*(1.0+X)*V
            HU = BETA*(BETA+1.0)*U + 2.0*ALPHA**2*U/(1.0+X)
            HV = BETA*(BETA-1.0)*V + 2.0*ALPHA**2*V/(1.0+X)
         ELSE
            END IF
         END IF
      RETURN
C
   19 CONTINUE
      IF (X.LT.0.0) THEN
         IF (GAMMA.EQ.0.0) THEN
            U = 1.0
            V = 0.5*LOG((1.0+X)/(1.0-X))
            PUP = 0.0
            PVP = 1.0
            HU = Q(X)*U
            HV = Q(X)*V
         ELSE
            U = COS(GAMMA*LOG(1.0+X))
            V = SIN(GAMMA*LOG(1.0+X))
            PUP = -GAMMA*(1.0-X)*V
            PVP =  GAMMA*(1.0-X)*U
            HU = -GAMMA**2*U - GAMMA*V + 2.0*BETA**2*U/(1.0-X)
            HV = -GAMMA**2*V + GAMMA*U + 2.0*BETA**2*V/(1.0-X)
            END IF
      ELSE IF (X.GE.0.0) THEN
         IF (BETA.EQ.0.0) THEN
            U = 1.0
            V = 0.5*LOG((1.0+X)/(1.0-X))
            PUP = 0.0
            PVP = 1.0
            HU = Q(X)*U
            HV = Q(X)*V
         ELSE IF (BETA.GT.0.0 .AND. BETA.LT.0.5) THEN
            U = (1.0-X)**BETA
            V = (1.0-X)**(-BETA)
            PUP = -BETA*(1.0+X)*U
            PVP = BETA*(1.0+X)*V
            HU = BETA*(BETA+1.0)*U - 2.0*GAMMA**2*U/(1.0+X)
            HV = BETA*(BETA-1.0)*V - 2.0*GAMMA**2*V/(1.0+X)
         ELSE
            END IF
         END IF
      RETURN
C
   20 CONTINUE
      U = SQRT(X)*COS(K*LOG(X))
      V = SQRT(X)*SIN(K*LOG(X))
      PUP = 0.5*U/X - K*V/X
      PVP = 0.5*V/X + K*U/X
      HU = U
      HV = V
      RETURN
C
   21 CONTINUE
      RETURN
C
   22 CONTINUE
      E = EXP(-X)
      IF (ALPHA.GT.-1.0 .AND. ALPHA.LT.0.0) THEN
         U = X**(-ALPHA)
         V = 1.0
         PUP = -ALPHA*E
         PVP = 0.0
         HU = -ALPHA*E
         HV = 0.0
      ELSE IF (ALPHA.EQ.0.0) THEN
         U = 1.0
         V = LOG(X)
         PUP = 0.0
         PVP = E
         HU = 0.0
         HV = E
      ELSE IF (ALPHA.GT.0.0 .AND. ALPHA.LT.1.0) THEN
         U = 1.0
         V = X**(-ALPHA)
         PUP = 0.0
         PVP = -ALPHA*E
         HU = 0.0
         HV = -ALPHA*E
      ELSE
         END IF
      RETURN
C
   23 CONTINUE
      IF (ALPHA.GT.-1.0 .AND. ALPHA.LT.0.0) THEN
         IF (ALPHA.NE.-0.5) THEN
            U = X**(0.5-ALPHA)
            V = X**(0.5+ALPHA)
            PUP = (0.5-ALPHA)*X**(-0.5-ALPHA)
            PVP = (0.5+ALPHA)*X**(-0.5+ALPHA)
         ELSE
            U = X
            V = 1.0
            PUP = 1.0
            PVP = 0.0
            END IF
         TX = X**2/16.0-(ALPHA+1.0)/2.0
         HU = TX*U
         HV = TX*V
      ELSE IF (ALPHA.EQ.0.0) THEN
         U = SQRT(X)
         V = U*LOG(X)
         PUP = 0.5/U
         PVP = (1.0+0.5*LOG(X))/U
         HU = (X**2/16.0-0.5)*U
         HV = (X**2/16.0-0.5)*V
      ELSE IF (ALPHA.GT.0.0 .AND. ALPHA.LT.1.0) THEN
         IF (ALPHA.NE.0.5) THEN
            U = X**(0.5+ALPHA)
            V = X**(0.5-ALPHA)
            PUP = (0.5+ALPHA)*X**(-0.5+ALPHA)
            PVP = (0.5-ALPHA)*X**(-0.5-ALPHA)
         ELSE
            U = X
            V = 1.0
            PUP = 1.0
            PVP = 0.0
            END IF
         TX = X**2/16.0-(ALPHA+1.0)/2.0
         HU = TX*U
         HV = TX*V
      ELSE
         END IF
      RETURN
C
   24 CONTINUE
      HPI = 2.0*ATAN(1.0)
      IF (X.GE.0.0) THEN
         IF (ALPHA.GT.-1.0 .AND. ALPHA.LT.0.0) THEN
            U = (HPI-X)**(0.5-ALPHA)
            V = (HPI-X)**(0.5+ALPHA)
            PUP = -(0.5-ALPHA)*(HPI-X)**(-0.5-ALPHA)
            PVP = -(0.5+ALPHA)*(HPI-X)**(-0.5+ALPHA)
            HU = (0.25-ALPHA**2)*(HPI-X)**(-1.5-ALPHA) + Q(X)*U
            HV = (0.25-ALPHA**2)*(HPI-X)**(-1.5+ALPHA) + Q(X)*V
         ELSE IF (ALPHA.EQ.0.0) THEN
            U = SQRT(HPI-X)
            V = U*LOG(HPI-X)
            PUP = -0.5/U
            PVP = -(1.0+0.5*LOG(HPI-X))/U
            HU = 0.25/((HPI-X)*U) + Q(X)*U
            HV = 0.25*LOG(HPI-X)/((HPI-X)*U) + Q(X)*V
         ELSE IF (ALPHA.GT.0.0 .AND. ALPHA.LT.1.0) THEN
            U = (HPI-X)**(0.5+ALPHA)
            V = (HPI-X)**(0.5-ALPHA)
            PUP = -(0.5+ALPHA)*(HPI-X)**(-0.5+ALPHA)
            PVP = -(0.5-ALPHA)*(HPI-X)**(-0.5-ALPHA)
            HU = (0.25-ALPHA**2)*(HPI-X)**(-1.5+ALPHA) + Q(X)*U
            HV = (0.25-ALPHA**2)*(HPI-X)**(-1.5-ALPHA) + Q(X)*V
            END IF
      ELSE
         IF (BETA.GT.-1.0 .AND. BETA.LT.0.0) THEN
            U = (HPI+X)**(0.5-BETA)
            V = (HPI+X)**(0.5+BETA)
            PUP = (0.5-BETA)*(HPI+X)**(-0.5-BETA)
            PVP = (0.5+BETA)*(HPI+X)**(-0.5+BETA)
            HU = (0.25-BETA**2)*(HPI+X)**(-1.5-BETA) + Q(X)*U
            HV = (0.25-BETA**2)*(HPI+X)**(-1.5+BETA) + Q(X)*V
         ELSE IF (BETA.EQ.0.0) THEN
            U = SQRT(HPI+X)
            V = U*LOG(HPI+X)
            PUP = 0.5/U
            PVP = (1.0+0.5*LOG(HPI+X))/U
            HU = 0.25/((HPI+X)*U) + Q(X)*U
            HV = 0.25*LOG(HPI+X)/((HPI+X)*U) + Q(X)*V
         ELSE IF (BETA.GT.0.0 .AND. BETA.LT.1.0) THEN
            U = (HPI+X)**(0.5+BETA)
            V = (HPI+X)**(0.5-BETA)
            PUP = (0.5+BETA)*(HPI+X)**(-0.5+BETA)
            PVP = (0.5-BETA)*(HPI+X)**(-0.5-BETA)
            HU = (0.25-BETA**2)*(HPI+X)**(-1.5+BETA) + Q(X)*U
            HV = (0.25-BETA**2)*(HPI+X)**(-1.5-BETA) + Q(X)*V
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
      END

