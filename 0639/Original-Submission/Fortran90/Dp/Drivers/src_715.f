C--**--CH3394--715--P:Gen--30:3:2000
C--**--CH2209--715--C:Gen--12:9:1999
CS    REAL FUNCTION ANORM(ARG)
      DOUBLE PRECISION FUNCTION ANORM(ARG)
C------------------------------------------------------------------
C
C This function evaluates the normal distribution function:
C
C                              / x    
C                     1       |       -t*t/2
C          P(x) = ----------- |      e       dt
C                 sqrt(2 pi)  |
C                             /-oo
C
C   The main computation evaluates near-minimax approximations
C   derived from those in "Rational Chebyshev approximations for
C   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
C   This transportable program uses rational functions that
C   theoretically approximate the normal distribution function to
C   at least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   XMIN  = the smallest positive floating-point number.
C
C Then the following machine-dependent constants must be declared 
C   in DATA statements.  IEEE values are provided as a default.
C
C   EPS   = argument below which anorm(x) may be represented by
C           0.5  and above which  x*x  will not underflow.
C           A conservative value is the largest machine number X
C           such that   1.0 + X = 1.0   to machine precision.
C   XLOW  = the most negative argument for which ANORM does not
C           vanish.  This is the negative of the solution to 
C                    W(x) * (1-1/x**2) = XMIN,
C           where W(x) = exp(-x*x/2)/[x*sqrt(2*pi)].
C   XUPPR = positive argument beyond which anorm = 1.0.  A
C           conservative value is the solution to the equation
C                    exp(-x*x/2) = EPS,
C           i.e., XUPPR = sqrt[-2 ln(eps)].
C
C   Approximate values for some important machines are:
C
C                          XMIN        EPS        XLOW    XUPPR
C
C  CDC 7600      (S.P.)  3.13E-294   7.11E-15   -36.641   8.072
C  CRAY-1        (S.P.)  4.58E-246   7.11E-157 -106.521  26.816
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  1.18E-38    5.96E-8    -12.949   5.768
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  2.23D-308   1.11D-16   -37.519   8.572
C  IBM 195       (D.P.)  5.40D-79    1.39D-17   -18.781   8.811
C  VAX D-Format  (D.P.)  2.94D-39    1.39D-17   -13.055   8.811
C  VAX G-Format  (D.P.)  5.56D-309   1.11D-16   -37.556   8.572
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ANORM = 0     for  ARG .LE. XLOW.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 15, 1992
C
C------------------------------------------------------------------
      INTEGER I
CS    REAL
      DOUBLE PRECISION
     1     A,ARG,B,C,D,DEL,EPS,HALF,P,ONE,Q,RESULT,SIXTEN,
     2     SQRPI,THRSH,ROOT32,X,XLOW,XDEN,XNUM,Y,XSQ,XUPPR,ZERO
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
C------------------------------------------------------------------
C  Mathematical constants
C
C  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
C  THRSH is the argument for which anorm = 0.75.
C------------------------------------------------------------------
CS    DATA ONE,HALF,ZERO,SIXTEN/1.0E0,0.5E0,0.0E0,1.60E1/,
CS   1     SQRPI/3.9894228040143267794E-1/,THRSH/0.66291E0/,
CS   2     ROOT32/5.656854248E0/
      DATA ONE,HALF,ZERO,SIXTEN/1.0D0,0.5D0,0.0D0,1.60D1/,
     1     SQRPI/3.9894228040143267794D-1/,THRSH/0.66291D0/,
     2     ROOT32/5.656854248D0/
C------------------------------------------------------------------
C  Machine-dependent constants
C------------------------------------------------------------------
CS    DATA EPS/5.96E-8/,XLOW/-12.949E0/,XUPPR/5.768E0/
      DATA EPS/1.11D-16/,XLOW/-37.519D0/,XUPPR/8.572D0/
C------------------------------------------------------------------
C  Coefficients for approximation in first interval
C------------------------------------------------------------------
CS    DATA A/2.2352520354606839287E00,1.6102823106855587881E02,
CS   1       1.0676894854603709582E03,1.8154981253343561249E04,
CS   2       6.5682337918207449113E-2/
CS    DATA B/4.7202581904688241870E01,9.7609855173777669322E02,
CS   1       1.0260932208618978205E04,4.5507789335026729956E04/
      DATA A/2.2352520354606839287D00,1.6102823106855587881D02,
     1       1.0676894854603709582D03,1.8154981253343561249D04,
     2       6.5682337918207449113D-2/
      DATA B/4.7202581904688241870D01,9.7609855173777669322D02,
     1       1.0260932208618978205D04,4.5507789335026729956D04/
C------------------------------------------------------------------
C  Coefficients for approximation in second interval
C------------------------------------------------------------------
CS    DATA C/3.9894151208813466764E-1,8.8831497943883759412E00,
CS   1       9.3506656132177855979E01,5.9727027639480026226E02,
CS   2       2.4945375852903726711E03,6.8481904505362823326E03,
CS   3       1.1602651437647350124E04,9.8427148383839780218E03,
CS   4       1.0765576773720192317E-8/
CS    DATA D/2.2266688044328115691E01,2.3538790178262499861E02,
CS   1       1.5193775994075548050E03,6.4855582982667607550E03,
CS   2       1.8615571640885098091E04,3.4900952721145977266E04,
CS   3       3.8912003286093271411E04,1.9685429676859990727E04/
      DATA C/3.9894151208813466764D-1,8.8831497943883759412D00,
     1       9.3506656132177855979D01,5.9727027639480026226D02,
     2       2.4945375852903726711D03,6.8481904505362823326D03,
     3       1.1602651437647350124D04,9.8427148383839780218D03,
     4       1.0765576773720192317D-8/
      DATA D/2.2266688044328115691D01,2.3538790178262499861D02,
     1       1.5193775994075548050D03,6.4855582982667607550D03,
     2       1.8615571640885098091D04,3.4900952721145977266D04,
     3       3.8912003286093271411D04,1.9685429676859990727D04/
C------------------------------------------------------------------
C  Coefficients for approximation in third interval
C------------------------------------------------------------------
CS    DATA P/2.1589853405795699E-1,1.274011611602473639E-1,
CS   1       2.2235277870649807E-2,1.421619193227893466E-3,
CS   2       2.9112874951168792E-5,2.307344176494017303E-2/
CS    DATA Q/1.28426009614491121E00,4.68238212480865118E-1,
CS   1       6.59881378689285515E-2,3.78239633202758244E-3,
CS   2       7.29751555083966205E-5/
      DATA P/2.1589853405795699D-1,1.274011611602473639D-1,
     1       2.2235277870649807D-2,1.421619193227893466D-3,
     2       2.9112874951168792D-5,2.307344176494017303D-2/
      DATA Q/1.28426009614491121D00,4.68238212480865118D-1,
     1       6.59881378689285515D-2,3.78239633202758244D-3,
     2       7.29751555083966205D-5/
C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y .LE. THRSH) THEN
C------------------------------------------------------------------
C  Evaluate  anorm  for  |X| <= 0.66291
C------------------------------------------------------------------
            XSQ = ZERO
            IF (Y .GT. EPS) XSQ = X * X
            XNUM = A(5)*XSQ
            XDEN = XSQ
            DO 20 I = 1, 3
               XNUM = (XNUM + A(I)) * XSQ
               XDEN = (XDEN + B(I)) * XSQ
   20       CONTINUE
            RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
            RESULT = HALF + RESULT
C------------------------------------------------------------------
C  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
C------------------------------------------------------------------
         ELSE IF (Y .LE. ROOT32) THEN
            XNUM = C(9)*Y
            XDEN = Y
            DO 120 I = 1, 7
               XNUM = (XNUM + C(I)) * Y
               XDEN = (XDEN + D(I)) * Y
  120       CONTINUE
            RESULT = (XNUM + C(8)) / (XDEN + D(8)) 
            XSQ = AINT(Y*SIXTEN)/SIXTEN
            DEL = (Y-XSQ)*(Y+XSQ)
            RESULT = EXP(-XSQ*XSQ*HALF)*EXP(-DEL*HALF)*RESULT
            IF (X .GT. ZERO) RESULT = ONE - RESULT
C------------------------------------------------------------------
C  Evaluate  anorm  for |X| > sqrt(32)
C------------------------------------------------------------------
         ELSE
            RESULT = ZERO
            IF ((X .GE. XLOW) .AND. (X .LT. XUPPR)) THEN
               XSQ = ONE / (X * X)
               XNUM = P(6)*XSQ
               XDEN = XSQ
               DO 240 I = 1, 4
                  XNUM = (XNUM + P(I)) * XSQ
                  XDEN = (XDEN + Q(I)) * XSQ
  240          CONTINUE
               RESULT = XSQ *(XNUM + P(5)) / (XDEN + Q(5))
               RESULT = (SQRPI -  RESULT) / Y
               XSQ = AINT(X*SIXTEN)/SIXTEN
               DEL = (X-XSQ)*(X+XSQ)
               RESULT = EXP(-XSQ*XSQ*HALF)*EXP(-DEL*HALF)*RESULT
            END IF
            IF (X .GT. ZERO) RESULT = ONE - RESULT
      END IF
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
      ANORM = RESULT
C---------- Last card of ANORM ----------
      END
      SUBROUTINE CALCI0(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functions of the first kind
C   and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
C   arguments X.  It contains two function type subprograms, BESI0
C   and BESEI0, and one subroutine type subprogram, CALCI0.
C   The calling statements for the primary entries are
C
C                   Y=BESI0(X)
C   and
C                   Y=BESEI0(X)
C
C   where the entry points correspond to the functions I0(X) and
C   EXP(-ABS(X))*I0(X), respectively.  The routine CALCI0 is
C   intended for internal packet use only, all computations within
C   the packet being concentrated in this routine.  The function
C   subprograms invoke CALCI0 with the statement
C          CALL CALCI0(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALCI0
C       Call              ARG                  RESULT          JINT
C
C     BESI0(ARG)    ABS(ARG) .LE. XMAX        I0(ARG)           1
C     BESEI0(ARG)    any real ARG        EXP(-ABS(ARG))*I0(ARG) 2
C
C   The main computation evaluates slightly modified forms of
C   minimax approximations generated by Blair and Edwards, Chalk
C   River (Atomic Energy of Canada Limited) Report AECL-4928,
C   October, 1974.  This transportable program is patterned after
C   the machine-dependent FUNPACK packet NATSI0, but cannot match
C   that version for efficiency or accuracy.  This version uses
C   rational functions that theoretically approximate I-SUB-0(X)
C   to at least 18 significant decimal digits.  The accuracy
C   achieved depends on the arithmetic system, the compiler, the
C   intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   maxexp = Smallest power of beta that overflows
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   XSMALL = Positive argument such that 1.0 - X = 1.0 to
C            machine precision for all ABS(X) .LE. XSMALL.
C   XINF =   Largest positive machine number; approximately
C            beta**maxexp
C   XMAX =   Largest argument acceptable to BESI0;  Solution to
C            equation:
C               W(X) * (1+1/(8*X)+9/(128*X**2) = beta**maxexp
C            where  W(X) = EXP(X)/SQRT(2*PI*X)
C
C
C     Approximate values for some important machines are:
C
C                          beta       maxexp       XSMALL
C
C CRAY-1        (S.P.)       2         8191       3.55E-15
C Cyber 180/855
C   under NOS   (S.P.)       2         1070       3.55E-15
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       2          128       2.98E-8
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       2         1024       5.55D-17
C IBM 3033      (D.P.)      16           63       6.95D-18
C VAX           (S.P.)       2          127       2.98E-8
C VAX D-Format  (D.P.)       2          127       6.95D-18
C VAX G-Format  (D.P.)       2         1023       5.55D-17
C
C
C                               XINF          XMAX
C
C CRAY-1        (S.P.)       5.45E+2465     5682.810
C Cyber 180/855
C   under NOS   (S.P.)       1.26E+322       745.893
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       3.40E+38         91.900
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       1.79D+308       713.986
C IBM 3033      (D.P.)       7.23D+75        178.182
C VAX           (S.P.)       1.70D+38         91.203
C VAX D-Format  (D.P.)       1.70D+38         91.203
C VAX G-Format  (D.P.)       8.98D+307       713.293
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns XINF for BESI0 for ABS(ARG) .GT. XMAX.
C
C
C  Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Authors: W. J. Cody and L. Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C  Latest modification: March 12, 1992
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1       A,ARG,B,EXP40,FORTY,ONE,ONE5,P,PP,Q,QQ,RESULT,
     2       REC15,SUMP,SUMQ,TWO25,X,XINF,XMAX,XSMALL,XX
      DIMENSION P(15),PP(8),Q(5),QQ(7)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ONE5/15.0E0/,EXP40/2.353852668370199854E17/,
CS   1     FORTY/40.0E0/,REC15/6.6666666666666666666E-2/,
CS   2     TWO25/225.0E0/
      DATA ONE/1.0D0/,ONE5/15.0D0/,EXP40/2.353852668370199854D17/,
     1     FORTY/40.0D0/,REC15/6.6666666666666666666D-2/,
     2     TWO25/225.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XSMALL/2.98E-8/,XINF/3.40E38/,XMAX/91.9E0/
      DATA XSMALL/5.55D-17/,XINF/1.79D308/,XMAX/713.986D0/
C--------------------------------------------------------------------
C  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
C--------------------------------------------------------------------
CS    DATA  P/-5.2487866627945699800E-18,-1.5982226675653184646E-14,
CS   1        -2.6843448573468483278E-11,-3.0517226450451067446E-08,
CS   2        -2.5172644670688975051E-05,-1.5453977791786851041E-02,
CS   3        -7.0935347449210549190E+00,-2.4125195876041896775E+03,
CS   4        -5.9545626019847898221E+05,-1.0313066708737980747E+08,
CS   5        -1.1912746104985237192E+10,-8.4925101247114157499E+11,
CS   6        -3.2940087627407749166E+13,-5.5050369673018427753E+14,
CS   7        -2.2335582639474375249E+15/
CS    DATA  Q/-3.7277560179962773046E+03, 6.5158506418655165707E+06,
CS   1        -6.5626560740833869295E+09, 3.7604188704092954661E+12,
CS   2        -9.7087946179594019126E+14/
      DATA  P/-5.2487866627945699800D-18,-1.5982226675653184646D-14,
     1        -2.6843448573468483278D-11,-3.0517226450451067446D-08,
     2        -2.5172644670688975051D-05,-1.5453977791786851041D-02,
     3        -7.0935347449210549190D+00,-2.4125195876041896775D+03,
     4        -5.9545626019847898221D+05,-1.0313066708737980747D+08,
     5        -1.1912746104985237192D+10,-8.4925101247114157499D+11,
     6        -3.2940087627407749166D+13,-5.5050369673018427753D+14,
     7        -2.2335582639474375249D+15/
      DATA  Q/-3.7277560179962773046D+03, 6.5158506418655165707D+06,
     1        -6.5626560740833869295D+09, 3.7604188704092954661D+12,
     2        -9.7087946179594019126D+14/
C--------------------------------------------------------------------
C  Coefficients for 15.0 .LE. ABS(ARG)
C--------------------------------------------------------------------
CS    DATA PP/-3.9843750000000000000E-01, 2.9205384596336793945E+00,
CS   1        -2.4708469169133954315E+00, 4.7914889422856814203E-01,
CS   2        -3.7384991926068969150E-03,-2.6801520353328635310E-03,
CS   3         9.9168777670983678974E-05,-2.1877128189032726730E-06/
CS    DATA QQ/-3.1446690275135491500E+01, 8.5539563258012929600E+01,
CS   1        -6.0228002066743340583E+01, 1.3982595353892851542E+01,
CS   2        -1.1151759188741312645E+00, 3.2547697594819615062E-02,
CS   3        -5.5194330231005480228E-04/
      DATA PP/-3.9843750000000000000D-01, 2.9205384596336793945D+00,
     1        -2.4708469169133954315D+00, 4.7914889422856814203D-01,
     2        -3.7384991926068969150D-03,-2.6801520353328635310D-03,
     3         9.9168777670983678974D-05,-2.1877128189032726730D-06/
      DATA QQ/-3.1446690275135491500D+01, 8.5539563258012929600D+01,
     1        -6.0228002066743340583D+01, 1.3982595353892851542D+01,
     2        -1.1151759188741312645D+00, 3.2547697594819615062D-02,
     3        -5.5194330231005480228D-04/
C--------------------------------------------------------------------
      X = ABS(ARG)
      IF (X .LT. XSMALL) THEN
            RESULT = ONE
         ELSE IF (X .LT. ONE5) THEN
C--------------------------------------------------------------------
C  XSMALL .LE.  ABS(ARG)  .LT. 15.0
C--------------------------------------------------------------------
            XX = X * X
            SUMP = P(1)
            DO 50 I = 2, 15
              SUMP = SUMP * XX + P(I)
   50       CONTINUE
            XX = XX - TWO25
            SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))*XX+Q(5)
            RESULT = SUMP / SUMQ
            IF (JINT .EQ. 2) RESULT = RESULT * EXP(-X)
         ELSE IF (X .GE. ONE5) THEN
            IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
                  RESULT = XINF
               ELSE
C--------------------------------------------------------------------
C  15.0  .LE.  ABS(ARG)
C--------------------------------------------------------------------
                  XX = ONE / X - REC15
                  SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+PP(4))*XX+
     1                   PP(5))*XX+PP(6))*XX+PP(7))*XX+PP(8)
                  SUMQ = ((((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+
     1                   QQ(4))*XX+QQ(5))*XX+QQ(6))*XX+QQ(7)
                  RESULT = SUMP / SUMQ
                  IF (JINT .EQ. 2) THEN
                        RESULT = (RESULT - PP(1)) / SQRT(X)
                     ELSE
C--------------------------------------------------------------------
C  Calculation reformulated to avoid premature overflow
C--------------------------------------------------------------------
                        IF (X .LE.(XMAX-ONE5)) THEN
                              A = EXP(X)
                              B = ONE
                           ELSE
                              A = EXP(X-FORTY)
                              B = EXP40
                        END IF
                        RESULT = ((RESULT*A-PP(1)*A)/SQRT(X))*B
                  END IF
            END IF
      END IF
C--------------------------------------------------------------------
C  Return for ABS(ARG) .LT. XSMALL
C--------------------------------------------------------------------
      RETURN
C----------- Last line of CALCI0 -----------
      END
CS    REAL FUNCTION BESI0(X)
      DOUBLE PRECISION FUNCTION BESI0(X)
C--------------------------------------------------------------------
C
C This long precision subprogram computes approximate values for
C   modified Bessel functions of the first kind of order zero for
C   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALCI0(X,RESULT,JINT)
      BESI0=RESULT
      RETURN
C---------- Last line of BESI0 ----------
      END
CS    REAL FUNCTION BESEI0(X)
      DOUBLE PRECISION FUNCTION BESEI0(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the first kind of order zero
C   multiplied by EXP(-ABS(X)), where EXP is the
C   exponential function, ABS is the absolute value, and X
C   is any argument.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=2
      CALL CALCI0(X,RESULT,JINT)
      BESEI0=RESULT
      RETURN
C---------- Last line of BESEI0 ----------
      END
      SUBROUTINE CALCI1(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functioons of the first kind
C   and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
C   arguments X.  It contains two function type subprograms, BESI1
C   and BESEI1, and one subroutine type subprogram, CALCI1.
C   The calling statements for the primary entries are
C
C                   Y=BESI1(X)
C   and
C                   Y=BESEI1(X)
C
C   where the entry points correspond to the functions I1(X) and
C   EXP(-ABS(X))*I1(X), respectively.  The routine CALCI1 is
C   intended for internal packet use only, all computations within
C   the packet being concentrated in this routine.  The function
C   subprograms invoke CALCI1 with the statement
C          CALL CALCI1(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALCI1
C       Call              ARG                  RESULT          JINT
C
C     BESI1(ARG)    ABS(ARG) .LE. XMAX        I1(ARG)           1
C     BESEI1(ARG)    any real ARG        EXP(-ABS(ARG))*I1(ARG) 2
C
C   The main computation evaluates slightly modified forms of
C   minimax approximations generated by Blair and Edwards, Chalk
C   River (Atomic Energy of Canada Limited) Report AECL-4928,
C   October, 1974.  This transportable program is patterned after
C   the machine-dependent FUNPACK packet NATSI1, but cannot match
C   that version for efficiency or accuracy.  This version uses
C   rational functions that theoretically approximate I-SUB-1(X)
C   to at least 18 significant decimal digits.  The accuracy
C   achieved depends on the arithmetic system, the compiler, the
C   intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   maxexp = Smallest power of beta that overflows
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   XSMALL = Positive argument such that 1.0 - X = 1.0 to
C            machine precision for all ABS(X) .LE. XSMALL.
C   XINF =   Largest positive machine number; approximately
C            beta**maxexp
C   XMAX =   Largest argument acceptable to BESI1;  Solution to
C            equation:
C               EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp
C
C
C     Approximate values for some important machines are:
C
C                          beta       maxexp       XSMALL
C
C CRAY-1        (S.P.)       2         8191       3.55E-15
C Cyber 180/855
C   under NOS   (S.P.)       2         1070       3.55E-15
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       2          128       2.98E-8
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       2         1024       5.55D-17
C IBM 3033      (D.P.)      16           63       6.95D-18
C VAX           (S.P.)       2          127       2.98E-8
C VAX D-Format  (D.P.)       2          127       6.95D-18
C VAX G-Format  (D.P.)       2         1023       5.55D-17
C
C
C                               XINF          XMAX
C
C CRAY-1        (S.P.)       5.45E+2465     5682.810
C Cyber 180/855
C   under NOS   (S.P.)       1.26E+322       745.894
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       3.40E+38         91.906
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       1.79D+308       713.987
C IBM 3033      (D.P.)       7.23D+75        178.185
C VAX           (S.P.)       1.70D+38         91.209
C VAX D-Format  (D.P.)       1.70D+38         91.209
C VAX G-Format  (D.P.)       8.98D+307       713.293
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for ABS(ARG) .GT. XMAX.
C
C
C Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Authors: W. J. Cody and L. Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL  60439
C
C  Latest modification: March 13, 1992
C
C--------------------------------------------------------------------
      INTEGER J,JINT
CS    REAL
      DOUBLE PRECISION
     1    A,ARG,B,EXP40,FORTY,HALF,ONE,ONE5,P,PBAR,PP,Q,QQ,REC15,
     2    RESULT,SUMP,SUMQ,TWO25,X,XINF,XMAX,XSMALL,XX,ZERO
      DIMENSION P(15),PP(8),Q(5),QQ(6)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ONE5/15.0E0/,EXP40/2.353852668370199854E17/,
CS   1     FORTY/40.0E0/,REC15/6.6666666666666666666E-2/,
CS   2     TWO25/225.0E0/,HALF/0.5E0/,ZERO/0.0E0/
      DATA ONE/1.0D0/,ONE5/15.0D0/,EXP40/2.353852668370199854D17/,
     1     FORTY/40.0D0/,REC15/6.6666666666666666666D-2/,
     2     TWO25/225.0D0/,HALF/0.5D0/,ZERO/0.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XSMALL/2.98E-8/,XINF/3.4E38/,XMAX/91.906E0/
      DATA XSMALL/5.55D-17/,XINF/1.79D308/,XMAX/713.987D0/
C--------------------------------------------------------------------
C  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
C--------------------------------------------------------------------
CS    DATA P/-1.9705291802535139930E-19,-6.5245515583151902910E-16,
CS   1       -1.1928788903603238754E-12,-1.4831904935994647675E-09,
CS   2       -1.3466829827635152875E-06,-9.1746443287817501309E-04,
CS   3       -4.7207090827310162436E-01,-1.8225946631657315931E+02,
CS   4       -5.1894091982308017540E+04,-1.0588550724769347106E+07,
CS   5       -1.4828267606612366099E+09,-1.3357437682275493024E+11,
CS   6       -6.9876779648010090070E+12,-1.7732037840791591320E+14,
CS   7       -1.4577180278143463643E+15/
      DATA P/-1.9705291802535139930D-19,-6.5245515583151902910D-16,
     1       -1.1928788903603238754D-12,-1.4831904935994647675D-09,
     2       -1.3466829827635152875D-06,-9.1746443287817501309D-04,
     3       -4.7207090827310162436D-01,-1.8225946631657315931D+02,
     4       -5.1894091982308017540D+04,-1.0588550724769347106D+07,
     5       -1.4828267606612366099D+09,-1.3357437682275493024D+11,
     6       -6.9876779648010090070D+12,-1.7732037840791591320D+14,
     7       -1.4577180278143463643D+15/
CS    DATA Q/-4.0076864679904189921E+03, 7.4810580356655069138E+06,
CS   1       -8.0059518998619764991E+09, 4.8544714258273622913E+12,
CS   2       -1.3218168307321442305E+15/
      DATA Q/-4.0076864679904189921D+03, 7.4810580356655069138D+06,
     1       -8.0059518998619764991D+09, 4.8544714258273622913D+12,
     2       -1.3218168307321442305D+15/
C--------------------------------------------------------------------
C  Coefficients for 15.0 .LE. ABS(ARG)
C--------------------------------------------------------------------
CS    DATA PP/-6.0437159056137600000E-02, 4.5748122901933459000E-01,
CS   1        -4.2843766903304806403E-01, 9.7356000150886612134E-02,
CS   2        -3.2457723974465568321E-03,-3.6395264712121795296E-04,
CS   3         1.6258661867440836395E-05,-3.6347578404608223492E-07/
      DATA PP/-6.0437159056137600000D-02, 4.5748122901933459000D-01,
     1        -4.2843766903304806403D-01, 9.7356000150886612134D-02,
     2        -3.2457723974465568321D-03,-3.6395264712121795296D-04,
     3         1.6258661867440836395D-05,-3.6347578404608223492D-07/
CS    DATA QQ/-3.8806586721556593450E+00, 3.2593714889036996297E+00,
CS   1        -8.5017476463217924408E-01, 7.4212010813186530069E-02,
CS   2        -2.2835624489492512649E-03, 3.7510433111922824643E-05/
      DATA QQ/-3.8806586721556593450D+00, 3.2593714889036996297D+00,
     1        -8.5017476463217924408D-01, 7.4212010813186530069D-02,
     2        -2.2835624489492512649D-03, 3.7510433111922824643D-05/
CS    DATA PBAR/3.98437500E-01/
      DATA PBAR/3.98437500D-01/
C--------------------------------------------------------------------
      X = ABS(ARG)
      IF (X .LT. XSMALL) THEN
C--------------------------------------------------------------------
C  Return for ABS(ARG) .LT. XSMALL
C--------------------------------------------------------------------
            RESULT = HALF * X
         ELSE IF (X .LT. ONE5) THEN
C--------------------------------------------------------------------
C  XSMALL .LE. ABS(ARG) .LT. 15.0
C--------------------------------------------------------------------
            XX = X * X
            SUMP = P(1)
            DO 50 J = 2, 15
               SUMP = SUMP * XX + P(J)
   50          CONTINUE
            XX = XX - TWO25
            SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))
     1           *XX+Q(5)
            RESULT = (SUMP / SUMQ) * X
            IF (JINT .EQ. 2) RESULT = RESULT * EXP(-X)
         ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
                  RESULT = XINF
         ELSE
C--------------------------------------------------------------------
C  15.0 .LE. ABS(ARG)
C--------------------------------------------------------------------
            XX = ONE / X - REC15
            SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+
     1           PP(4))*XX+PP(5))*XX+PP(6))*XX+PP(7))*XX+PP(8)
            SUMQ = (((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+
     1           QQ(4))*XX+QQ(5))*XX+QQ(6)
            RESULT = SUMP / SUMQ
            IF (JINT .NE. 1) THEN
                  RESULT = (RESULT + PBAR) / SQRT(X)
               ELSE
C--------------------------------------------------------------------
C  Calculation reformulated to avoid premature overflow
C--------------------------------------------------------------------
                  IF (X .GT. XMAX-ONE5) THEN
                        A = EXP(X-FORTY)
                        B = EXP40
                     ELSE
                        A = EXP(X)
                        B = ONE
                  END IF
                  RESULT = ((RESULT * A + PBAR * A) /
     1                  SQRT(X)) * B
C--------------------------------------------------------------------
C  Error return for ABS(ARG) .GT. XMAX
C--------------------------------------------------------------------
            END IF
      END IF
      IF (ARG .LT. ZERO) RESULT = -RESULT
      RETURN
C----------- Last line of CALCI1 -----------
      END
CS    REAL FUNCTION BESI1(X)
      DOUBLE PRECISION FUNCTION BESI1(X)
C--------------------------------------------------------------------
C
C This long precision subprogram computes approximate values for
C   modified Bessel functions of the first kind of order one for
C   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI1).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALCI1(X,RESULT,JINT)
      BESI1=RESULT
      RETURN
C---------- Last line of BESI1 ----------
      END
CS    REAL FUNCTION BESEI1(X)
      DOUBLE PRECISION FUNCTION BESEI1(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the first kind of order one
C   multiplied by EXP(-ABS(X)), where EXP is the
C   exponential function, ABS is the absolute value, and X
C   is any argument.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=2
      CALL CALCI1(X,RESULT,JINT)
      BESEI1=RESULT
      RETURN
C---------- Last line of BESEI1 ----------
      END
      SUBROUTINE CALJY0(ARG,RESULT,JINT)
C---------------------------------------------------------------------
C
C This packet computes zero-order Bessel functions of the first and
C   second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
C   for Y0, and |X| <= XMAX for J0.  It contains two function-type
C   subprograms,  BESJ0  and  BESY0,  and one subroutine-type
C   subprogram,  CALJY0.  The calling statements for the primary
C   entries are:
C
C           Y = BESJ0(X)
C   and
C           Y = BESY0(X),
C
C   where the entry points correspond to the functions J0(X) and Y0(X),
C   respectively.  The routine  CALJY0  is intended for internal packet
C   use only, all computations within the packet being concentrated in
C   this one routine.  The function subprograms invoke  CALJY0  with
C   the statement
C           CALL CALJY0(ARG,RESULT,JINT),
C   where the parameter usage is as follows:
C
C      Function                  Parameters for CALJY0
C       call              ARG             RESULT          JINT
C
C     BESJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
C     BESY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1
C
C   The main computation uses unpublished minimax rational
C   approximations for X .LE. 8.0, and an approximation from the
C   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
C   New York, 1968, for arguments larger than 8.0   Part of this
C   transportable packet is patterned after the machine-dependent
C   FUNPACK program BESJ0(X), but cannot match that version for
C   efficiency or accuracy.  This version uses rational functions
C   that are theoretically accurate to at least 18 significant decimal
C   digits for X <= 8, and at least 18 decimal places for X > 8.  The
C   accuracy achieved depends on the arithmetic system, the compiler,
C   the intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C
C The following machine-dependent constants must be declared in
C   DATA statements.  IEEE values are provided as a default.
C
C   XINF   = largest positive machine number
C   XMAX   = largest acceptable argument.  The functions AINT, SIN
C            and COS must perform properly for  ABS(X) .LE. XMAX.
C            We recommend that XMAX be a small integer multiple of
C            sqrt(1/eps), where eps is the smallest positive number
C            such that  1+eps > 1.
C   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
C            to machine precision for all  ABS(X) .LE. XSMALL.
C            We recommend that  XSMALL < sqrt(eps)/beta, where beta
C            is the floating-point radix (usually 2 or 16).
C
C     Approximate values for some important machines are
C
C                          eps      XMAX     XSMALL      XINF
C
C  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
C  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
C  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
C  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
C  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
C  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
C  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns the value zero for  X .GT. XMAX, and returns
C    -XINF when BESLY0 is called with a negative or zero argument.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, COS, LOG, SIN, SQRT
C
C
C  Latest modification: March 13, 1992
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1       ARG,AX,CONS,DOWN,EIGHT,FIVE5,FOUR,ONE,ONEOV8,PI2,PJ0,
     2       PJ1,PLG,PROD,PY0,PY1,PY2,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,
     3       QY2,Q0,Q1,RESJ,RESULT,R0,R1,SIXTY4,THREE,TWOPI,TWOPI1,
     4       TWOPI2,TWO56,UP,W,WSQ,XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,
     5       XJ1,XJ01,XJ02,XJ11,XJ12,XY,XY0,XY01,XY02,XY1,XY11,XY12,
     6       XY2,XY21,XY22,Z,ZERO,ZSQ
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(6),PY1(7),PY2(8),P0(6),P1(6),
     1          QJ0(5),QJ1(7),QLG(4),QY0(5),QY1(6),QY2(7),Q0(5),Q1(5)
C-------------------------------------------------------------------
C  Mathematical constants
C    CONS = ln(.5) + Euler's gamma
C-------------------------------------------------------------------
CS    DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0E0,1.0E0,3.0E0,4.0E0,8.0E0/,
CS   1     FIVE5,SIXTY4,ONEOV8,P17/5.5E0,64.0E0,0.125E0,1.716E-1/,
CS   2     TWO56,CONS/256.0E0,-1.1593151565841244881E-1/,
CS   3     PI2,TWOPI/6.3661977236758134308E-1,6.2831853071795864769E0/,
CS   4     TWOPI1,TWOPI2/6.28125E0,1.9353071795864769253E-3/
      DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0D0,1.0D0,3.0D0,4.0D0,8.0D0/,
     1     FIVE5,SIXTY4,ONEOV8,P17/5.5D0,64.0D0,0.125D0,1.716D-1/,
     2     TWO56,CONS/256.0D0,-1.1593151565841244881D-1/,
     3     PI2,TWOPI/6.3661977236758134308D-1,6.2831853071795864769D0/,
     4     TWOPI1,TWOPI2/6.28125D0,1.9353071795864769253D-3/
C-------------------------------------------------------------------
C  Machine-dependent constants
C-------------------------------------------------------------------
CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-04/,XINF/3.40E+38/
      DATA XMAX/2.68D+08/,XSMALL/3.72D-09/,XINF/1.79D+308/
C-------------------------------------------------------------------
C  Zeroes of Bessel functions
C-------------------------------------------------------------------
CS    DATA XJ0/2.4048255576957727686E+0/,XJ1/5.5200781102863106496E+0/,
CS   1     XY0/8.9357696627916752158E-1/,XY1/3.9576784193148578684E+0/,
CS   2     XY2/7.0860510603017726976E+0/,
CS   3     XJ01/ 616.0E+0/, XJ02/-1.4244423042272313784E-03/,
CS   4     XJ11/1413.0E+0/, XJ12/ 5.4686028631064959660E-04/,
CS   5     XY01/ 228.0E+0/, XY02/ 2.9519662791675215849E-03/,
CS   6     XY11/1013.0E+0/, XY12/ 6.4716931485786837568E-04/,
CS   7     XY21/1814.0E+0/, XY22/ 1.1356030177269762362E-04/
      DATA XJ0/2.4048255576957727686D+0/,XJ1/5.5200781102863106496D+0/,
     1     XY0/8.9357696627916752158D-1/,XY1/3.9576784193148578684D+0/,
     2     XY2/7.0860510603017726976D+0/,
     3     XJ01/ 616.0D+0/, XJ02/-1.4244423042272313784D-03/,
     4     XJ11/1413.0D+0/, XJ12/ 5.4686028631064959660D-04/,
     5     XY01/ 228.0D+0/, XY02/ 2.9519662791675215849D-03/,
     6     XY11/1013.0D+0/, XY12/ 6.4716931485786837568D-04/,
     7     XY21/1814.0D+0/, XY22/ 1.1356030177269762362D-04/
C-------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a)
C--------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PJ0/6.6302997904833794242E+06,-6.2140700423540120665E+08,
CS   1         2.7282507878605942706E+10,-4.1298668500990866786E+11,
CS   2        -1.2117036164593528341E-01, 1.0344222815443188943E+02,
CS   3        -3.6629814655107086448E+04/
CS    DATA QJ0/4.5612696224219938200E+05, 1.3985097372263433271E+08,
CS   1         2.6328198300859648632E+10, 2.3883787996332290397E+12,
CS   2         9.3614022392337710626E+02/
      DATA PJ0/6.6302997904833794242D+06,-6.2140700423540120665D+08,
     1         2.7282507878605942706D+10,-4.1298668500990866786D+11,
     2        -1.2117036164593528341D-01, 1.0344222815443188943D+02,
     3        -3.6629814655107086448D+04/
      DATA QJ0/4.5612696224219938200D+05, 1.3985097372263433271D+08,
     1         2.6328198300859648632D+10, 2.3883787996332290397D+12,
     2         9.3614022392337710626D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
C-------------------------------------------------------------------
CS    DATA PJ1/4.4176707025325087628E+03, 1.1725046279757103576E+04,
CS   1         1.0341910641583726701E+04,-7.2879702464464618998E+03,
CS   2        -1.2254078161378989535E+04,-1.8319397969392084011E+03,
CS   3         4.8591703355916499363E+01, 7.4321196680624245801E+02/
CS    DATA QJ1/3.3307310774649071172E+02,-2.9458766545509337327E+03,
CS   1         1.8680990008359188352E+04,-8.4055062591169562211E+04,
CS   2         2.4599102262586308984E+05,-3.5783478026152301072E+05,
CS   3        -2.5258076240801555057E+01/
      DATA PJ1/4.4176707025325087628D+03, 1.1725046279757103576D+04,
     1         1.0341910641583726701D+04,-7.2879702464464618998D+03,
     2        -1.2254078161378989535D+04,-1.8319397969392084011D+03,
     3         4.8591703355916499363D+01, 7.4321196680624245801D+02/
      DATA QJ1/3.3307310774649071172D+02,-2.9458766545509337327D+03,
     1         1.8680990008359188352D+04,-8.4055062591169562211D+04,
     2         2.4599102262586308984D+05,-3.5783478026152301072D+05,
     3        -2.5258076240801555057D+01/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
C        XSMALL  <  |X|  <=  3.0
C--------------------------------------------------------------------
CS    DATA PY0/1.0102532948020907590E+04,-2.1287548474401797963E+06,
CS   1         2.0422274357376619816E+08,-8.3716255451260504098E+09,
CS   2         1.0723538782003176831E+11,-1.8402381979244993524E+01/
CS    DATA QY0/6.6475986689240190091E+02, 2.3889393209447253406E+05,
CS   1         5.5662956624278251596E+07, 8.1617187777290363573E+09,
CS   2         5.8873865738997033405E+11/
      DATA PY0/1.0102532948020907590D+04,-2.1287548474401797963D+06,
     1         2.0422274357376619816D+08,-8.3716255451260504098D+09,
     2         1.0723538782003176831D+11,-1.8402381979244993524D+01/
      DATA QY0/6.6475986689240190091D+02, 2.3889393209447253406D+05,
     1         5.5662956624278251596D+07, 8.1617187777290363573D+09,
     2         5.8873865738997033405D+11/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
C        3.0  <  |X|  <=  5.5
C--------------------------------------------------------------------
CS    DATA PY1/-1.4566865832663635920E+04, 4.6905288611678631510E+06,
CS   1         -6.9590439394619619534E+08, 4.3600098638603061642E+10,
CS   2         -5.5107435206722644429E+11,-2.2213976967566192242E+13,
CS   3          1.7427031242901594547E+01/
CS    DATA QY1/ 8.3030857612070288823E+02, 4.0669982352539552018E+05,
CS   1          1.3960202770986831075E+08, 3.4015103849971240096E+10,
CS   2          5.4266824419412347550E+12, 4.3386146580707264428E+14/
      DATA PY1/-1.4566865832663635920D+04, 4.6905288611678631510D+06,
     1         -6.9590439394619619534D+08, 4.3600098638603061642D+10,
     2         -5.5107435206722644429D+11,-2.2213976967566192242D+13,
     3          1.7427031242901594547D+01/
      DATA QY1/ 8.3030857612070288823D+02, 4.0669982352539552018D+05,
     1          1.3960202770986831075D+08, 3.4015103849971240096D+10,
     2          5.4266824419412347550D+12, 4.3386146580707264428D+14/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
C        5.5  <  |X|  <=  8.0
C--------------------------------------------------------------------
CS    DATA PY2/ 2.1363534169313901632E+04,-1.0085539923498211426E+07,
CS   1          2.1958827170518100757E+09,-1.9363051266772083678E+11,
CS   2         -1.2829912364088687306E+11, 6.7016641869173237784E+14,
CS   3         -8.0728726905150210443E+15,-1.7439661319197499338E+01/
CS    DATA QY2/ 8.7903362168128450017E+02, 5.3924739209768057030E+05,
CS   1          2.4727219475672302327E+08, 8.6926121104209825246E+10,
CS   2          2.2598377924042897629E+13, 3.9272425569640309819E+15,
CS   3          3.4563724628846457519E+17/
      DATA PY2/ 2.1363534169313901632D+04,-1.0085539923498211426D+07,
     1          2.1958827170518100757D+09,-1.9363051266772083678D+11,
     2         -1.2829912364088687306D+11, 6.7016641869173237784D+14,
     3         -8.0728726905150210443D+15,-1.7439661319197499338D+01/
      DATA QY2/ 8.7903362168128450017D+02, 5.3924739209768057030D+05,
     1          2.4727219475672302327D+08, 8.6926121104209825246D+10,
     2          2.2598377924042897629D+13, 3.9272425569640309819D+15,
     3          3.4563724628846457519D+17/
C-------------------------------------------------------------------
C  Coefficients for Hart,s approximation,  |X| > 8.0
C-------------------------------------------------------------------
CS    DATA P0/3.4806486443249270347E+03, 2.1170523380864944322E+04,
CS   1        4.1345386639580765797E+04, 2.2779090197304684302E+04,
CS   2        8.8961548424210455236E-01, 1.5376201909008354296E+02/
CS    DATA Q0/3.5028735138235608207E+03, 2.1215350561880115730E+04,
CS   1        4.1370412495510416640E+04, 2.2779090197304684318E+04,
CS   2        1.5711159858080893649E+02/
CS    DATA P1/-2.2300261666214198472E+01,-1.1183429920482737611E+02,
CS   1        -1.8591953644342993800E+02,-8.9226600200800094098E+01,
CS   2        -8.8033303048680751817E-03,-1.2441026745835638459E+00/
CS    DATA Q1/1.4887231232283756582E+03, 7.2642780169211018836E+03,
CS   1        1.1951131543434613647E+04, 5.7105024128512061905E+03,
CS   2        9.0593769594993125859E+01/
      DATA P0/3.4806486443249270347D+03, 2.1170523380864944322D+04,
     1        4.1345386639580765797D+04, 2.2779090197304684302D+04,
     2        8.8961548424210455236D-01, 1.5376201909008354296D+02/
      DATA Q0/3.5028735138235608207D+03, 2.1215350561880115730D+04,
     1        4.1370412495510416640D+04, 2.2779090197304684318D+04,
     2        1.5711159858080893649D+02/
      DATA P1/-2.2300261666214198472D+01,-1.1183429920482737611D+02,
     1        -1.8591953644342993800D+02,-8.9226600200800094098D+01,
     2        -8.8033303048680751817D-03,-1.2441026745835638459D+00/
      DATA Q1/1.4887231232283756582D+03, 7.2642780169211018836D+03,
     1        1.1951131543434613647D+04, 5.7105024128512061905D+03,
     2        9.0593769594993125859D+01/
C-------------------------------------------------------------------
C  Check for error conditions
C-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. (ARG .LE. ZERO)) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) GO TO 800
      IF (AX .LE. XSMALL) THEN
         IF (JINT .EQ. 0) THEN
               RESULT = ONE
            ELSE
               RESULT = PI2 * (LOG(AX) + CONS)
         END IF
         GO TO 2000
      END IF
C-------------------------------------------------------------------
C  Calculate J0 for appropriate interval, preserving
C     accuracy near the zero of J0
C-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(5) * ZSQ + PJ0(6)) * ZSQ + PJ0(7)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            WSQ = ONE - ZSQ / SIXTY4
            XNUM = PJ1(7) * WSQ + PJ1(8)
            XDEN = WSQ + QJ1(7)
            DO 220 I = 1, 6
               XNUM = XNUM * WSQ + PJ1(I)
               XDEN = XDEN * WSQ + QJ1(I)
  220       CONTINUE
            PROD = (AX + XJ1) * ((AX - XJ11/TWO56) - XJ12)
      END IF
      RESULT = PROD * XNUM / XDEN
      IF (JINT .EQ. 0) GO TO 2000
C-------------------------------------------------------------------
C  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
C    where xn is a zero of Y0
C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE IF (AX .LE. FIVE5) THEN
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
         ELSE
            UP = (AX-XY21/TWO56)-XY22
            XY = XY2
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
C-------------------------------------------------------------------
C  Now calculate Y0 for appropriate interval, preserving
C     accuracy near the zero of Y0
C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            XNUM = PY0(6) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 5
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE IF (AX .LE. FIVE5) THEN
            XNUM = PY1(7) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 6
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
         ELSE
            XNUM = PY2(8) * ZSQ + PY2(1)
            XDEN = ZSQ + QY2(1)
            DO 380 I = 2, 7
               XNUM = XNUM * ZSQ + PY2(I)
               XDEN = XDEN * ZSQ + QY2(I)
  380       CONTINUE
      END IF
      RESULT = RESJ + UP * DOWN * XNUM / XDEN
      GO TO 2000
C-------------------------------------------------------------------
C  Calculate J0 or Y0 for |ARG|  >  8.0
C-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AX / TWOPI
      W = AINT(W) + ONEOV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(5) * ZSQ + P0(6)
      XDEN = ZSQ + Q0(5)
      UP = P1(5) * ZSQ + P1(6)
      DOWN = ZSQ + Q1(5)
      DO 850 I = 1, 4
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = SQRT(PI2/AX) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = SQRT(PI2/AX) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
 2000 RETURN
C---------- Last line of CALJY0 ----------
      END
CS    REAL FUNCTION BESJ0(X)
      DOUBLE PRECISION FUNCTION BESJ0(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the first kind of order zero for arguments  |X| <= XMAX
C   (see comments heading CALJY0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      JINT=0
      CALL CALJY0(X,RESULT,JINT)
      BESJ0 = RESULT
      RETURN
C---------- Last line of BESJ0 ----------
      END
CS    REAL FUNCTION BESY0(X)
      DOUBLE PRECISION FUNCTION BESY0(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the second kind of order zero for arguments 0 < X <= XMAX
C   (see comments heading CALJY0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALJY0(X,RESULT,JINT)
      BESY0 = RESULT
      RETURN
C---------- Last line of BESY0 ----------
      END
      SUBROUTINE CALJY1(ARG,RESULT,JINT)
C---------------------------------------------------------------------
C
C This packet computes first-order Bessel functions of the first and
C   second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
C   for Y1, and |X| <= XMAX for J1.  It contains two function-type
C   subprograms,  BESJ1  and  BESY1,  and one subroutine-type
C   subprogram,  CALJY1.  The calling statements for the primary
C   entries are:
C
C           Y = BESJ1(X)
C   and
C           Y = BESY1(X),
C
C   where the entry points correspond to the functions J1(X) and Y1(X),
C   respectively.  The routine  CALJY1  is intended for internal packet
C   use only, all computations within the packet being concentrated in
C   this one routine.  The function subprograms invoke  CALJY1  with
C   the statement
C           CALL CALJY1(ARG,RESULT,JINT),
C   where the parameter usage is as follows:
C
C      Function                  Parameters for CALJY1
C       call              ARG             RESULT          JINT
C
C     BESJ1(ARG)     |ARG| .LE. XMAX       J1(ARG)          0
C     BESY1(ARG)   0 .LT. ARG .LE. XMAX    Y1(ARG)          1
C
C   The main computation uses unpublished minimax rational
C   approximations for X .LE. 8.0, and an approximation from the
C   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
C   New York, 1968, for arguments larger than 8.0   Part of this
C   transportable packet is patterned after the machine-dependent
C   FUNPACK program BESJ1(X), but cannot match that version for
C   efficiency or accuracy.  This version uses rational functions
C   that are theoretically accurate to at least 18 significant decimal
C   digits for X <= 8, and at least 18 decimal places for X > 8.  The
C   accuracy achieved depends on the arithmetic system, the compiler,
C   the intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C
C The following machine-dependent constants must be declared in
C   DATA statements.  IEEE values are provided as a default.
C
C   XINF   = largest positive machine number
C   XMAX   = largest acceptable argument.  The functions AINT, SIN
C            and COS must perform properly for  ABS(X) .LE. XMAX.
C            We recommend that XMAX be a small integer multiple of
C            sqrt(1/eps), where eps is the smallest positive number
C            such that  1+eps > 1.
C   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
C            to machine precision for all  ABS(X) .LE. XSMALL.
C            We recommend that  XSMALL < sqrt(eps)/beta, where beta
C            is the floating-point radix (usually 2 or 16).
C
C     Approximate values for some important machines are
C
C                          eps      XMAX     XSMALL      XINF
C
C  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
C  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
C  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
C  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
C  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
C  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
C  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns the value zero for  X .GT. XMAX, and returns
C    -XINF when BESLY1 is called with a negative or zero argument.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, COS, LOG, SIN, SQRT
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 13, 1992
C
C--------------------------------------------------------------------
      INTEGER I,JINT
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(7),PY1(9),P0(6),P1(6),
     1          QJ0(5),QJ1(7),QLG(4),QY0(6),QY1(8),Q0(6),Q1(6)
CS    REAL
      DOUBLE PRECISION
     1   ARG,AX,DOWN,EIGHT,FOUR,HALF,PI2,PJ0,PJ1,PLG,PROD,PY0,
     2   PY1,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,Q0,Q1,RESJ,RESULT,
     3   RTPI2,R0,R1,THROV8,TWOPI,TWOPI1,TWOPI2,TWO56,UP,W,WSQ,
     4   XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,XJ1,XJ01,XJ02,XJ11,XJ12,
     5   XY,XY0,XY01,XY02,XY1,XY11,XY12,Z,ZERO,ZSQ
C-------------------------------------------------------------------
C  Mathematical constants
C-------------------------------------------------------------------
CS    DATA EIGHT/8.0E0/,
CS   1     FOUR/4.0E0/,HALF/0.5E0/,THROV8/0.375E0/,
CS   2     PI2/6.3661977236758134308E-1/,P17/1.716E-1/
CS   3     TWOPI/6.2831853071795864769E+0/,ZERO/0.0E0/,
CS   4     TWOPI1/6.28125E0/,TWOPI2/1.9353071795864769253E-03/
CS   5     TWO56/256.0E+0/,RTPI2/7.9788456080286535588E-1/
      DATA EIGHT/8.0D0/,
     1     FOUR/4.0D0/,HALF/0.5D0/,THROV8/0.375D0/,
     2     PI2/6.3661977236758134308D-1/,P17/1.716D-1/
     3     TWOPI/6.2831853071795864769D+0/,ZERO/0.0D0/,
     4     TWOPI1/6.28125D0/,TWOPI2/1.9353071795864769253D-03/
     5     TWO56/256.0D+0/,RTPI2/7.9788456080286535588D-1/
C-------------------------------------------------------------------
C  Machine-dependent constants
C-------------------------------------------------------------------
CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-04/,XINF/3.40E+38/
      DATA XMAX/2.68D+08/,XSMALL/3.72D-09/,XINF/1.79D+308/
C-------------------------------------------------------------------
C  Zeroes of Bessel functions
C-------------------------------------------------------------------
CS    DATA XJ0/3.8317059702075123156E+0/,XJ1/7.0155866698156187535E+0/,
CS   1     XY0/2.1971413260310170351E+0/,XY1/5.4296810407941351328E+0/,
CS   2     XJ01/ 981.0E+0/, XJ02/-3.2527979248768438556E-04/,
CS   3     XJ11/1796.0E+0/, XJ12/-3.8330184381246462950E-05/,
CS   4     XY01/ 562.0E+0/, XY02/ 1.8288260310170351490E-03/,
CS   5     XY11/1390.0E+0/, XY12/-6.4592058648672279948E-06/
      DATA XJ0/3.8317059702075123156D+0/,XJ1/7.0155866698156187535D+0/,
     1     XY0/2.1971413260310170351D+0/,XY1/5.4296810407941351328D+0/,
     2     XJ01/ 981.0D+0/, XJ02/-3.2527979248768438556D-04/,
     3     XJ11/1796.0D+0/, XJ12/-3.8330184381246462950D-05/,
     4     XY01/ 562.0D+0/, XY02/ 1.8288260310170351490D-03/,
     5     XY11/1390.0D+0/, XY12/-6.4592058648672279948D-06/
C-------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a)
C--------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PJ0/9.8062904098958257677E+05,-1.1548696764841276794E+08,
CS   1       6.6781041261492395835E+09,-1.4258509801366645672E+11,
CS   2      -4.4615792982775076130E+03, 1.0650724020080236441E+01,
CS   3      -1.0767857011487300348E-02/
CS    DATA QJ0/5.9117614494174794095E+05, 2.0228375140097033958E+08,
CS   1       4.2091902282580133541E+10, 4.1868604460820175290E+12,
CS   2       1.0742272239517380498E+03/
      DATA PJ0/9.8062904098958257677D+05,-1.1548696764841276794D+08,
     1       6.6781041261492395835D+09,-1.4258509801366645672D+11,
     2      -4.4615792982775076130D+03, 1.0650724020080236441D+01,
     3      -1.0767857011487300348D-02/
      DATA QJ0/5.9117614494174794095D+05, 2.0228375140097033958D+08,
     1       4.2091902282580133541D+10, 4.1868604460820175290D+12,
     2       1.0742272239517380498D+03/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
C-------------------------------------------------------------------
CS    DATA PJ1/4.6179191852758252280E+00,-7.1329006872560947377E+03,
CS   1       4.5039658105749078904E+06,-1.4437717718363239107E+09,
CS   2       2.3569285397217157313E+11,-1.6324168293282543629E+13,
CS   3       1.1357022719979468624E+14, 1.0051899717115285432E+15/
CS    DATA QJ1/1.1267125065029138050E+06, 6.4872502899596389593E+08,
CS   1       2.7622777286244082666E+11, 8.4899346165481429307E+13,
CS   2       1.7128800897135812012E+16, 1.7253905888447681194E+18,
CS   3       1.3886978985861357615E+03/
      DATA PJ1/4.6179191852758252280D+00,-7.1329006872560947377D+03,
     1       4.5039658105749078904D+06,-1.4437717718363239107D+09,
     2       2.3569285397217157313D+11,-1.6324168293282543629D+13,
     3       1.1357022719979468624D+14, 1.0051899717115285432D+15/
      DATA QJ1/1.1267125065029138050D+06, 6.4872502899596389593D+08,
     1       2.7622777286244082666D+11, 8.4899346165481429307D+13,
     2       1.7128800897135812012D+16, 1.7253905888447681194D+18,
     3       1.3886978985861357615D+03/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
C        XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PY0/2.2157953222280260820E+05,-5.9157479997408395984E+07,
CS   1         7.2144548214502560419E+09,-3.7595974497819597599E+11,
CS   2         5.4708611716525426053E+12, 4.0535726612579544093E+13,
CS   3        -3.1714424660046133456E+02/
CS    DATA QY0/8.2079908168393867438E+02, 3.8136470753052572164E+05,
CS   1         1.2250435122182963220E+08, 2.7800352738690585613E+10,
CS   2         4.1272286200406461981E+12, 3.0737873921079286084E+14/
      DATA PY0/2.2157953222280260820D+05,-5.9157479997408395984D+07,
     1         7.2144548214502560419D+09,-3.7595974497819597599D+11,
     2         5.4708611716525426053D+12, 4.0535726612579544093D+13,
     3        -3.1714424660046133456D+02/
      DATA QY0/8.2079908168393867438D+02, 3.8136470753052572164D+05,
     1         1.2250435122182963220D+08, 2.7800352738690585613D+10,
     2         4.1272286200406461981D+12, 3.0737873921079286084D+14/
C--------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
C        4.0  <  |X|  <=  8.0
C--------------------------------------------------------------------
CS    DATA PY1/ 1.9153806858264202986E+06,-1.1957961912070617006E+09,
CS   1          3.7453673962438488783E+11,-5.9530713129741981618E+13,
CS   2          4.0686275289804744814E+15,-2.3638408497043134724E+16,
CS   3         -5.6808094574724204577E+18, 1.1514276357909013326E+19,
CS   4         -1.2337180442012953128E+03/
CS    DATA QY1/ 1.2855164849321609336E+03, 1.0453748201934079734E+06,
CS   1          6.3550318087088919566E+08, 3.0221766852960403645E+11,
CS   2          1.1187010065856971027E+14, 3.0837179548112881950E+16,
CS   3          5.6968198822857178911E+18, 5.3321844313316185697E+20/
      DATA PY1/ 1.9153806858264202986D+06,-1.1957961912070617006D+09,
     1          3.7453673962438488783D+11,-5.9530713129741981618D+13,
     2          4.0686275289804744814D+15,-2.3638408497043134724D+16,
     3         -5.6808094574724204577D+18, 1.1514276357909013326D+19,
     4         -1.2337180442012953128D+03/
      DATA QY1/ 1.2855164849321609336D+03, 1.0453748201934079734D+06,
     1          6.3550318087088919566D+08, 3.0221766852960403645D+11,
     2          1.1187010065856971027D+14, 3.0837179548112881950D+16,
     3          5.6968198822857178911D+18, 5.3321844313316185697D+20/
C-------------------------------------------------------------------
C  Coefficients for Hart,s approximation,  |X| > 8.0
C-------------------------------------------------------------------
CS    DATA P0/-1.0982405543459346727E+05,-1.5235293511811373833E+06,
CS   1         -6.6033732483649391093E+06,-9.9422465050776411957E+06,
CS   2         -4.4357578167941278571E+06,-1.6116166443246101165E+03/
CS    DATA Q0/-1.0726385991103820119E+05,-1.5118095066341608816E+06,
CS   1         -6.5853394797230870728E+06,-9.9341243899345856590E+06,
CS   2         -4.4357578167941278568E+06,-1.4550094401904961825E+03/
CS    DATA P1/ 1.7063754290207680021E+03, 1.8494262873223866797E+04,
CS   1          6.6178836581270835179E+04, 8.5145160675335701966E+04,
CS   2          3.3220913409857223519E+04, 3.5265133846636032186E+01/
CS    DATA Q1/ 3.7890229745772202641E+04, 4.0029443582266975117E+05,
CS   1          1.4194606696037208929E+06, 1.8194580422439972989E+06,
CS   2          7.0871281941028743574E+05, 8.6383677696049909675E+02/
      DATA P0/-1.0982405543459346727D+05,-1.5235293511811373833D+06,
     1         -6.6033732483649391093D+06,-9.9422465050776411957D+06,
     2         -4.4357578167941278571D+06,-1.6116166443246101165D+03/
      DATA Q0/-1.0726385991103820119D+05,-1.5118095066341608816D+06,
     1         -6.5853394797230870728D+06,-9.9341243899345856590D+06,
     2         -4.4357578167941278568D+06,-1.4550094401904961825D+03/
      DATA P1/ 1.7063754290207680021D+03, 1.8494262873223866797D+04,
     1          6.6178836581270835179D+04, 8.5145160675335701966D+04,
     2          3.3220913409857223519D+04, 3.5265133846636032186D+01/
      DATA Q1/ 3.7890229745772202641D+04, 4.0029443582266975117D+05,
     1          1.4194606696037208929D+06, 1.8194580422439972989D+06,
     2          7.0871281941028743574D+05, 8.6383677696049909675D+02/
C-------------------------------------------------------------------
C  Check for error conditions
C-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. ((ARG .LE. ZERO) .OR.
     1   ((ARG .LT. HALF) .AND. (AX*XINF .LT. PI2)))) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) THEN
            GO TO 800
         ELSE IF (AX .LE. XSMALL) THEN
            IF (JINT .EQ. 0) THEN
                  RESULT = ARG * HALF
               ELSE
                  RESULT = -PI2 / AX
            END IF
            GO TO 2000
      END IF
C-------------------------------------------------------------------
C  Calculate J1 for appropriate interval, preserving
C     accuracy near the zero of J1
C-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(7) * ZSQ + PJ0(6)) * ZSQ + PJ0(5)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ARG * ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            XNUM = PJ1(1)
            XDEN = (ZSQ + QJ1(7)) * ZSQ + QJ1(1)
            DO 220 I = 2, 6
               XNUM = XNUM * ZSQ + PJ1(I)
               XDEN = XDEN * ZSQ + QJ1(I)
  220       CONTINUE
            XNUM = XNUM * (AX - EIGHT) * (AX + EIGHT) + PJ1(7)
            XNUM = XNUM * (AX - FOUR) * (AX + FOUR) + PJ1(8)
            PROD = ARG * ((AX - XJ11/TWO56) - XJ12) * (AX + XJ1)
      END IF
      RESULT = PROD * (XNUM / XDEN)
      IF (JINT .EQ. 0) GO TO 2000
C-------------------------------------------------------------------
C  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
C    where xn is a zero of Y1
C-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
C-------------------------------------------------------------------
C  Now calculate Y1 for appropriate interval, preserving
C     accuracy near the zero of Y1
C-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            XNUM = PY0(7) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 6
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE
            XNUM = PY1(9) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 8
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
      END IF
      RESULT = RESJ + (UP*DOWN/AX) * XNUM / XDEN
      GO TO 2000
C-------------------------------------------------------------------
C  Calculate J1 or Y1 for |ARG|  >  8.0
C-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AINT(AX/TWOPI) + THROV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(6)
      XDEN = ZSQ + Q0(6)
      UP = P1(6)
      DOWN = ZSQ + Q1(6)
      DO 850 I = 1, 5
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = (RTPI2/SQRT(AX)) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = (RTPI2/SQRT(AX)) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
      IF ((JINT .EQ. 0) .AND. (ARG .LT. ZERO)) RESULT = -RESULT
 2000 RETURN
C---------- Last card of CALJY1 ----------
      END
CS    REAL FUNCTION BESJ1(X)
      DOUBLE PRECISION FUNCTION BESJ1(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the first kind of order zero for arguments  |X| <= XMAX
C   (see comments heading CALJY1).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1   RESULT,X
C--------------------------------------------------------------------
      JINT=0
      CALL CALJY1(X,RESULT,JINT)
      BESJ1 = RESULT
      RETURN
C---------- Last card of BESJ1 ----------
      END
      SUBROUTINE CALCK1(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functions of the second kind
C   and order one,  K1(X)  and  EXP(X)*K1(X), for real arguments X.
C   It contains two function type subprograms, BESK1  and  BESEK1,
C   and one subroutine type subprogram, CALCK1.  The calling
C   statements for the primary entries are
C
C                   Y=BESK1(X)
C   and
C                   Y=BESEK1(X)
C
C   where the entry points correspond to the functions K1(X) and
C   EXP(X)*K1(X), respectively.  The routine CALCK1 is intended
C   for internal packet use only, all computations within the
C   packet being concentrated in this routine.  The function
C   subprograms invoke CALCK1 with the statement
C          CALL CALCK1(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                      Parameters for CALCK1
C        Call             ARG                  RESULT          JINT
C
C     BESK1(ARG)  XLEAST .LT. ARG .LT. XMAX    K1(ARG)          1
C     BESEK1(ARG)     XLEAST .LT. ARG       EXP(ARG)*K1(ARG)    2
C
C   The main computation evaluates slightly modified forms of near
C   minimax rational approximations generated by Russon and Blair,
C   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
C   1969.  This transportable program is patterned after the
C   machine-dependent FUNPACK packet NATSK1, but cannot match that
C   version for efficiency or accuracy.  This version uses rational
C   functions that theoretically approximate K-SUB-1(X) to at
C   least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   XLEAST = Smallest acceptable argument, i.e., smallest machine
C            number X such that 1/X is machine representable.
C   XSMALL = Argument below which BESK1(X) and BESEK1(X) may
C            each be represented by 1/X.  A safe value is the
C            largest X such that  1.0 + X = 1.0  to machine
C            precision.
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMAX   = Largest argument acceptable to BESK1;  Solution to
C            equation:
C               W(X) * (1+3/8X-15/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C
C
C     Approximate values for some important machines are:
C
C                           beta       minexp       maxexp
C
C  CRAY-1        (S.P.)       2        -8193         8191
C  Cyber 180/185
C    under NOS   (S.P.)       2         -975         1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2         -126          128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2        -1022         1024
C  IBM 3033      (D.P.)      16          -65           63
C  VAX D-Format  (D.P.)       2         -128          127
C  VAX G-Format  (D.P.)       2        -1024         1023
C
C
C                         XLEAST     XSMALL      XINF       XMAX
C
C CRAY-1                1.84E-2466  3.55E-15  5.45E+2465  5674.858
C Cyber 180/855
C   under NOS   (S.P.)  3.14E-294   1.77E-15  1.26E+322    672.789
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.18E-38    5.95E-8   3.40E+38      85.343
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  2.23D-308   1.11D-16  1.79D+308    705.343
C IBM 3033      (D.P.)  1.39D-76    1.11D-16  7.23D+75     177.855
C VAX D-Format  (D.P.)  5.88D-39    6.95D-18  1.70D+38      86.721
C VAX G-Format  (D.P.)  1.12D-308   5.55D-17  8.98D+307    706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for ARG .LE. 0.0 and the
C   BESK1 entry returns the value 0.0 for ARG .GT. XMAX.
C
C
C  Intrinsic functions required are:
C
C     LOG, SQRT, EXP
C
C
C  Authors: W. J. Cody and Laura Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C  Latest modification: March 14, 1992
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1    ARG,F,G,ONE,P,PP,Q,QQ,RESULT,SUMF,SUMG,
     2    SUMP,SUMQ,X,XINF,XMAX,XLEAST,XSMALL,XX,ZERO
      DIMENSION P(5),Q(3),PP(11),QQ(9),F(5),G(3)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ZERO/0.0E0/
      DATA ONE/1.0D0/,ZERO/0.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XLEAST/1.18E-38/,XSMALL/5.95E-8/,XINF/3.40E+38/,
CS   1     XMAX/85.343E+0/
      DATA XLEAST/2.23D-308/,XSMALL/1.11D-16/,XINF/1.79D+308/,
     1     XMAX/705.343D+0/
C--------------------------------------------------------------------
C  Coefficients for  XLEAST .LE.  ARG  .LE. 1.0
C--------------------------------------------------------------------
CS    DATA   P/ 4.8127070456878442310E-1, 9.9991373567429309922E+1,
CS   1          7.1885382604084798576E+3, 1.7733324035147015630E+5,
CS   2          7.1938920065420586101E+5/
CS    DATA   Q/-2.8143915754538725829E+2, 3.7264298672067697862E+4,
CS   1         -2.2149374878243304548E+6/
CS    DATA   F/-2.2795590826955002390E-1,-5.3103913335180275253E+1,
CS   1         -4.5051623763436087023E+3,-1.4758069205414222471E+5,
CS   2         -1.3531161492785421328E+6/
CS    DATA   G/-3.0507151578787595807E+2, 4.3117653211351080007E+4,
CS   2         -2.7062322985570842656E+6/
      DATA   P/ 4.8127070456878442310D-1, 9.9991373567429309922D+1,
     1          7.1885382604084798576D+3, 1.7733324035147015630D+5,
     2          7.1938920065420586101D+5/
      DATA   Q/-2.8143915754538725829D+2, 3.7264298672067697862D+4,
     1         -2.2149374878243304548D+6/
      DATA   F/-2.2795590826955002390D-1,-5.3103913335180275253D+1,
     1         -4.5051623763436087023D+3,-1.4758069205414222471D+5,
     2         -1.3531161492785421328D+6/
      DATA   G/-3.0507151578787595807D+2, 4.3117653211351080007D+4,
     2         -2.7062322985570842656D+6/
C--------------------------------------------------------------------
C  Coefficients for  1.0 .LT.  ARG
C--------------------------------------------------------------------
CS    DATA  PP/ 6.4257745859173138767E-2, 7.5584584631176030810E+0,
CS   1          1.3182609918569941308E+2, 8.1094256146537402173E+2,
CS   2          2.3123742209168871550E+3, 3.4540675585544584407E+3,
CS   3          2.8590657697910288226E+3, 1.3319486433183221990E+3,
CS   4          3.4122953486801312910E+2, 4.4137176114230414036E+1,
CS   5          2.2196792496874548962E+0/
CS    DATA  QQ/ 3.6001069306861518855E+1, 3.3031020088765390854E+2,
CS   1          1.2082692316002348638E+3, 2.1181000487171943810E+3,
CS   2          1.9448440788918006154E+3, 9.6929165726802648634E+2,
CS   3          2.5951223655579051357E+2, 3.4552228452758912848E+1,
CS   4          1.7710478032601086579E+0/
      DATA  PP/ 6.4257745859173138767D-2, 7.5584584631176030810D+0,
     1          1.3182609918569941308D+2, 8.1094256146537402173D+2,
     2          2.3123742209168871550D+3, 3.4540675585544584407D+3,
     3          2.8590657697910288226D+3, 1.3319486433183221990D+3,
     4          3.4122953486801312910D+2, 4.4137176114230414036D+1,
     5          2.2196792496874548962D+0/
      DATA  QQ/ 3.6001069306861518855D+1, 3.3031020088765390854D+2,
     1          1.2082692316002348638D+3, 2.1181000487171943810D+3,
     2          1.9448440788918006154D+3, 9.6929165726802648634D+2,
     3          2.5951223655579051357D+2, 3.4552228452758912848D+1,
     4          1.7710478032601086579D+0/
C--------------------------------------------------------------------
      X = ARG
      IF (X .LT. XLEAST) THEN
C--------------------------------------------------------------------
C  Error return for  ARG  .LT. XLEAST
C--------------------------------------------------------------------
            RESULT = XINF
         ELSE IF (X .LE. ONE) THEN
C--------------------------------------------------------------------
C  XLEAST .LE.  ARG  .LE. 1.0
C--------------------------------------------------------------------
            IF (X .LT. XSMALL) THEN
C--------------------------------------------------------------------
C  Return for small ARG
C--------------------------------------------------------------------
                  RESULT = ONE / X
               ELSE
                  XX = X * X
                  SUMP = ((((P(1)*XX + P(2))*XX + P(3))*XX + P(4))*XX
     1                   + P(5))*XX + Q(3)
                  SUMQ = ((XX + Q(1))*XX + Q(2))*XX + Q(3)
                  SUMF = (((F(1)*XX + F(2))*XX + F(3))*XX + F(4))*XX
     1                   + F(5)
                  SUMG = ((XX + G(1))*XX + G(2))*XX + G(3)
                  RESULT = (XX * LOG(X) * SUMF/SUMG + SUMP/SUMQ) / X
                  IF (JINT .EQ. 2) RESULT = RESULT * EXP(X)
            END IF
         ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
C--------------------------------------------------------------------
C  Error return for  ARG  .GT. XMAX
C--------------------------------------------------------------------
            RESULT = ZERO
         ELSE
C--------------------------------------------------------------------
C  1.0 .LT.  ARG
C--------------------------------------------------------------------
            XX = ONE / X
            SUMP = PP(1)
            DO 120 I = 2, 11
               SUMP = SUMP * XX + PP(I)
  120       CONTINUE
            SUMQ = XX
            DO 140 I = 1, 8
               SUMQ = (SUMQ + QQ(I)) * XX
  140       CONTINUE
            SUMQ = SUMQ + QQ(9)
            RESULT = SUMP / SUMQ / SQRT(X)
            IF (JINT .EQ. 1) RESULT = RESULT * EXP(-X)
      END IF
      RETURN
C---------- Last line of CALCK1 ----------
      END
CS    REAL FUNCTION BESK1(X)
      DOUBLE PRECISION FUNCTION BESK1(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order one
C   for arguments  XLEAST .LE. ARG .LE. XMAX.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL 
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 1
      CALL CALCK1(X,RESULT,JINT)
      BESK1 = RESULT
      RETURN
C---------- Last line of BESK1 ----------
      END
CS    REAL FUNCTION BESEK1(X)
      DOUBLE PRECISION FUNCTION BESEK1(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order one
C   multiplied by the exponential function, for arguments
C   XLEAST .LE. ARG .LE. XMAX.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 2
      CALL CALCK1(X,RESULT,JINT)
      BESEK1 = RESULT
      RETURN
C---------- Last line of BESEK1 ----------
      END
CS    REAL FUNCTION BESY1(X)
      DOUBLE PRECISION FUNCTION BESY1(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the second kind of order zero for arguments 0 < X <= XMAX
C   (see comments heading CALJY1).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1   RESULT,X
C--------------------------------------------------------------------
      JINT=1
      CALL CALJY1(X,RESULT,JINT)
      BESY1 = RESULT
      RETURN
C---------- Last card of BESY1 ----------
      END
      SUBROUTINE CALCK0(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functions of the second kind
C   and order zero, K0(X) and EXP(X)*K0(X), for real
C   arguments X.  It contains two function type subprograms, BESK0
C   and BESEK0, and one subroutine type subprogram, CALCK0.
C   the calling statements for the primary entries are
C
C                   Y=BESK0(X)
C   and
C                   Y=BESEK0(X)
C
C   where the entry points correspond to the functions K0(X) and
C   EXP(X)*K0(X), respectively.  The routine CALCK0 is
C   intended for internal packet use only, all computations within
C   the packet being concentrated in this routine.  The function
C   subprograms invoke CALCK0 with the statement
C          CALL CALCK0(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALCK0
C       Call              ARG                  RESULT          JINT
C
C     BESK0(ARG)   0 .LT. ARG .LE. XMAX       K0(ARG)           1
C     BESEK0(ARG)     0 .LT. ARG           EXP(ARG)*K0(ARG)     2
C
C   The main computation evaluates slightly modified forms of near
C   minimax rational approximations generated by Russon and Blair,
C   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
C   1969.  This transportable program is patterned after the
C   machine-dependent FUNPACK packet NATSK0, but cannot match that
C   version for efficiency or accuracy.  This version uses rational
C   functions that theoretically approximate K-SUB-0(X) to at
C   least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   XSMALL = Argument below which BESK0 and BESEK0 may
C            each be represented by a constant and a log.
C            largest X such that  1.0 + X = 1.0  to machine
C            precision.
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMAX   = Largest argument acceptable to BESK0;  Solution to
C            equation:
C               W(X) * (1-1/8X+9/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C
C
C     Approximate values for some important machines are:
C
C
C                           beta       minexp       maxexp
C
C  CRAY-1        (S.P.)       2        -8193         8191
C  Cyber 180/185
C    under NOS   (S.P.)       2         -975         1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2         -126          128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2        -1022         1024
C  IBM 3033      (D.P.)      16          -65           63
C  VAX D-Format  (D.P.)       2         -128          127
C  VAX G-Format  (D.P.)       2        -1024         1023
C
C
C                          XSMALL       XINF         XMAX
C
C CRAY-1        (S.P.)    3.55E-15   5.45E+2465    5674.858
C Cyber 180/855
C   under NOS   (S.P.)    1.77E-15   1.26E+322      672.788
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)    5.95E-8    3.40E+38        85.337
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)    1.11D-16   1.79D+308      705.342
C IBM 3033      (D.P.)    1.11D-16   7.23D+75       177.852
C VAX D-Format  (D.P.)    6.95D-18   1.70D+38        86.715
C VAX G-Format  (D.P.)    5.55D-17   8.98D+307      706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for ARG .LE. 0.0, and the
C  BESK0 entry returns the value 0.0 for ARG .GT. XMAX.
C
C
C  Intrinsic functions required are:
C
C     EXP, LOG, SQRT
C
C  Latest modification: March 14, 1992
C
C  Authors: W. J. Cody and Laura Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1    ARG,F,G,ONE,P,PP,Q,QQ,RESULT,SUMF,SUMG,SUMP,SUMQ,TEMP,
     2    X,XINF,XMAX,XSMALL,XX,ZERO
      DIMENSION P(6),Q(2),PP(10),QQ(10),F(4),G(3)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ZERO/0.0E0/
      DATA ONE/1.0D0/,ZERO/0.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XSMALL/5.95E-8/,XINF/3.40E+38/,XMAX/ 85.337E0/
      DATA XSMALL/1.11D-16/,XINF/1.79D+308/,XMAX/705.342D0/
C--------------------------------------------------------------------
C
C     Coefficients for XSMALL .LE.  ARG  .LE. 1.0
C
C--------------------------------------------------------------------
CS    DATA   P/ 5.8599221412826100000E-04, 1.3166052564989571850E-01,
CS   1          1.1999463724910714109E+01, 4.6850901201934832188E+02,
CS   2          5.9169059852270512312E+03, 2.4708152720399552679E+03/
CS    DATA   Q/-2.4994418972832303646E+02, 2.1312714303849120380E+04/
CS    DATA   F/-1.6414452837299064100E+00,-2.9601657892958843866E+02,
CS   1         -1.7733784684952985886E+04,-4.0320340761145482298E+05/
CS    DATA   G/-2.5064972445877992730E+02, 2.9865713163054025489E+04,
CS   1         -1.6128136304458193998E+06/
      DATA   P/ 5.8599221412826100000D-04, 1.3166052564989571850D-01,
     1          1.1999463724910714109D+01, 4.6850901201934832188D+02,
     2          5.9169059852270512312D+03, 2.4708152720399552679D+03/
      DATA   Q/-2.4994418972832303646D+02, 2.1312714303849120380D+04/
      DATA   F/-1.6414452837299064100D+00,-2.9601657892958843866D+02,
     1         -1.7733784684952985886D+04,-4.0320340761145482298D+05/
      DATA   G/-2.5064972445877992730D+02, 2.9865713163054025489D+04,
     1         -1.6128136304458193998D+06/
C--------------------------------------------------------------------
C
C     Coefficients for  1.0 .LT. ARG
C
C--------------------------------------------------------------------
CS    DATA  PP/ 1.1394980557384778174E+02, 3.6832589957340267940E+03,
CS   1          3.1075408980684392399E+04, 1.0577068948034021957E+05,
CS   2          1.7398867902565686251E+05, 1.5097646353289914539E+05,
CS   3          7.1557062783764037541E+04, 1.8321525870183537725E+04,
CS   4          2.3444738764199315021E+03, 1.1600249425076035558E+02/
CS    DATA  QQ/ 2.0013443064949242491E+02, 4.4329628889746408858E+03,
CS   1          3.1474655750295278825E+04, 9.7418829762268075784E+04,
CS   2          1.5144644673520157801E+05, 1.2689839587977598727E+05,
CS   3          5.8824616785857027752E+04, 1.4847228371802360957E+04,
CS   4          1.8821890840982713696E+03, 9.2556599177304839811E+01/
      DATA  PP/ 1.1394980557384778174D+02, 3.6832589957340267940D+03,
     1          3.1075408980684392399D+04, 1.0577068948034021957D+05,
     2          1.7398867902565686251D+05, 1.5097646353289914539D+05,
     3          7.1557062783764037541D+04, 1.8321525870183537725D+04,
     4          2.3444738764199315021D+03, 1.1600249425076035558D+02/
      DATA  QQ/ 2.0013443064949242491D+02, 4.4329628889746408858D+03,
     1          3.1474655750295278825D+04, 9.7418829762268075784D+04,
     2          1.5144644673520157801D+05, 1.2689839587977598727D+05,
     3          5.8824616785857027752D+04, 1.4847228371802360957D+04,
     4          1.8821890840982713696D+03, 9.2556599177304839811D+01/
C--------------------------------------------------------------------
      X = ARG
      IF (X .GT. ZERO) THEN
            IF (X .LE. ONE) THEN
C--------------------------------------------------------------------
C     0.0 .LT.  ARG  .LE. 1.0
C--------------------------------------------------------------------
                  TEMP = LOG(X)
                  IF (X .LT. XSMALL) THEN
C--------------------------------------------------------------------
C     Return for small ARG
C--------------------------------------------------------------------
                        RESULT = P(6)/Q(2) - TEMP
                     ELSE
                        XX = X * X
                        SUMP = ((((P(1)*XX + P(2))*XX + P(3))*XX +
     1                         P(4))*XX + P(5))*XX + P(6)
                        SUMQ = (XX + Q(1))*XX + Q(2)
                        SUMF = ((F(1)*XX + F(2))*XX + F(3))*XX + F(4)
                        SUMG = ((XX + G(1))*XX + G(2))*XX + G(3)
                        RESULT = SUMP/SUMQ - XX*SUMF*TEMP/SUMG - TEMP
                        IF (JINT .EQ. 2) RESULT = RESULT * EXP(X)
                  END IF
               ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
C--------------------------------------------------------------------
C     Error return for ARG .GT. XMAX
C--------------------------------------------------------------------
                  RESULT = ZERO
               ELSE
C--------------------------------------------------------------------
C     1.0 .LT. ARG
C--------------------------------------------------------------------
                  XX = ONE / X
                  SUMP = PP(1)
                  DO 120 I = 2, 10
                     SUMP = SUMP*XX + PP(I)
  120             CONTINUE
                  SUMQ = XX
                  DO 140 I = 1, 9
                     SUMQ = (SUMQ + QQ(I))*XX
  140             CONTINUE
                  SUMQ = SUMQ + QQ(10)
                  RESULT = SUMP / SUMQ / SQRT(X)
                  IF (JINT .EQ. 1) RESULT = RESULT * EXP(-X)
            END IF
         ELSE
C--------------------------------------------------------------------
C     Error return for ARG .LE. 0.0
C--------------------------------------------------------------------
            RESULT = XINF
      END IF
C--------------------------------------------------------------------
C     Update error counts, etc.
C--------------------------------------------------------------------
      RETURN
C---------- Last line of CALCK0 ----------
      END
CS    REAL FUNCTION BESK0(X)
      DOUBLE PRECISION FUNCTION BESK0(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order zero
C   for arguments 0.0 .LT. ARG .LE. XMAX (see comments heading
C   CALCK0).
C
C  Authors: W. J. Cody and Laura Stoltz
C
C  Latest Modification: March 14, 1992
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL 
      DOUBLE PRECISION
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 1
      CALL CALCK0(X,RESULT,JINT)
      BESK0 = RESULT
      RETURN
C---------- Last line of BESK0 ----------
      END
CS    REAL FUNCTION BESEK0(X)
      DOUBLE PRECISION FUNCTION BESEK0(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order zero
C   multiplied by the Exponential function, for arguments
C   0.0 .LT. ARG.
C
C  Authors: W. J. Cody and Laura Stoltz
C
C  Latest Modification: March 14, 1992
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL 
      DOUBLE PRECISION 
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 2
      CALL CALCK0(X,RESULT,JINT)
      BESEK0 = RESULT
      RETURN
C---------- Last line of BESEK0 ----------
      END
CS    REAL FUNCTION DAW(XX)
      DOUBLE PRECISION FUNCTION DAW(XX)
C----------------------------------------------------------------------
C
C This function program evaluates Dawson's integral, 
C
C                       2  / x   2
C                     -x   |    t
C             F(x) = e     |   e    dt
C                          |
C                          / 0
C
C   for a real argument x.
C
C   The calling sequence for this function is 
C
C                   Y=DAW(X)
C
C   The main computation uses rational Chebyshev approximations
C   published in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek
C   and Thacher.  This transportable program is patterned after the
C   machine-dependent FUNPACK program DDAW(X), but cannot match that
C   version for efficiency or accuracy.  This version uses rational
C   approximations that are theoretically accurate to about 19
C   significant decimal digits.  The accuracy achieved depends on the
C   arithmetic system, the compiler, the intrinsic functions, and
C   proper selection of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   XINF   = largest positive machine number
C   XMIN   = the smallest positive machine number.
C   EPS    = smallest positive number such that 1+eps > 1.
C            Approximately  beta**(-p), where beta is the machine
C            radix and p is the number of significant base-beta
C            digits in a floating-point number.
C
C Then the following machine-dependent constants must be declared 
C   in DATA statements.  IEEE values are provided as a default.
C
C   XMAX   = absolute argument beyond which DAW(X) underflows.
C            XMAX = min(0.5/xmin, xinf).
C   XSMALL = absolute argument below DAW(X)  may be represented
C            by X.  We recommend XSMALL = sqrt(eps).
C   XLARGE = argument beyond which DAW(X) may be represented by
C            1/(2x).  We recommend XLARGE = 1/sqrt(eps).
C
C     Approximate values for some important machines are
C
C                        beta  p     eps     xmin       xinf  
C
C  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322
C  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
C  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
C  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
C                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
C   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307
C
C                         XSMALL     XLARGE     XMAX    
C
C  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293
C  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307
C  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75
C  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38
C                (D.P.)  3.73E-09   2.68E+08  1.70E+38
C   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns 0.0 for |X| > XMAX.
C
C Intrinsic functions required are:
C
C     ABS
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division 
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 9, 1992
C
C----------------------------------------------------------------------
      INTEGER I
CS    REAL
      DOUBLE PRECISION
     1     FRAC,HALF,ONE,ONE225,P1,P2,P3,P4,Q1,Q2,Q3,Q4,SIX25,
     2     SUMP,SUMQ,TWO5,W2,X,XX,Y,XLARGE,XMAX,XSMALL,ZERO
      DIMENSION P1(10),P2(10),P3(10),P4(10),Q1(10),Q2(9),Q3(9),Q4(9)
C----------------------------------------------------------------------
C  Mathematical constants.
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE/0.0E0,0.5E0,1.0E0/,
CS   1     SIX25,ONE225,TWO5/6.25E0,12.25E0,25.0E0/
      DATA ZERO,HALF,ONE/0.0D0,0.5D0,1.0D0/,
     1     SIX25,ONE225,TWO5/6.25D0,12.25D0,25.0D0/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA XSMALL/2.44E-04/, XLARGE/4.10E+03/, XMAX/4.25E+37/
      DATA XSMALL/1.05D-08/, XLARGE/9.49D+07/, XMAX/2.24D+307/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation for  |x| < 2.5
C----------------------------------------------------------------------
CS    DATA P1/-2.69020398788704782410E-12, 4.18572065374337710778E-10,
CS   1        -1.34848304455939419963E-08, 9.28264872583444852976E-07,
CS   2        -1.23877783329049120592E-05, 4.07205792429155826266E-04,
CS   3        -2.84388121441008500446E-03, 4.70139022887204722217E-02,
CS   4        -1.38868086253931995101E-01, 1.00000000000000000004E+00/
CS    DATA Q1/ 1.71257170854690554214E-10, 1.19266846372297253797E-08,
CS   1         4.32287827678631772231E-07, 1.03867633767414421898E-05,
CS   2         1.78910965284246249340E-04, 2.26061077235076703171E-03,
CS   3         2.07422774641447644725E-02, 1.32212955897210128811E-01,
CS   4         5.27798580412734677256E-01, 1.00000000000000000000E+00/
      DATA P1/-2.69020398788704782410D-12, 4.18572065374337710778D-10,
     1        -1.34848304455939419963D-08, 9.28264872583444852976D-07,
     2        -1.23877783329049120592D-05, 4.07205792429155826266D-04,
     3        -2.84388121441008500446D-03, 4.70139022887204722217D-02,
     4        -1.38868086253931995101D-01, 1.00000000000000000004D+00/
      DATA Q1/ 1.71257170854690554214D-10, 1.19266846372297253797D-08,
     1         4.32287827678631772231D-07, 1.03867633767414421898D-05,
     2         1.78910965284246249340D-04, 2.26061077235076703171D-03,
     3         2.07422774641447644725D-02, 1.32212955897210128811D-01,
     4         5.27798580412734677256D-01, 1.00000000000000000000D+00/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation in J-fraction form
C     for  x in [2.5, 3.5)
C----------------------------------------------------------------------
CS    DATA P2/-1.70953804700855494930E+00,-3.79258977271042880786E+01,
CS   1         2.61935631268825992835E+01, 1.25808703738951251885E+01,
CS   2        -2.27571829525075891337E+01, 4.56604250725163310122E+00,
CS   3        -7.33080089896402870750E+00, 4.65842087940015295573E+01,
CS   4        -1.73717177843672791149E+01, 5.00260183622027967838E-01/
CS    DATA Q2/ 1.82180093313514478378E+00, 1.10067081034515532891E+03,
CS   1        -7.08465686676573000364E+00, 4.53642111102577727153E+02,
CS   2         4.06209742218935689922E+01, 3.02890110610122663923E+02,
CS   3         1.70641269745236227356E+02, 9.51190923960381458747E+02,
CS   4         2.06522691539642105009E-01/
      DATA P2/-1.70953804700855494930D+00,-3.79258977271042880786D+01,
     1         2.61935631268825992835D+01, 1.25808703738951251885D+01,
     2        -2.27571829525075891337D+01, 4.56604250725163310122D+00,
     3        -7.33080089896402870750D+00, 4.65842087940015295573D+01,
     4        -1.73717177843672791149D+01, 5.00260183622027967838D-01/
      DATA Q2/ 1.82180093313514478378D+00, 1.10067081034515532891D+03,
     1        -7.08465686676573000364D+00, 4.53642111102577727153D+02,
     2         4.06209742218935689922D+01, 3.02890110610122663923D+02,
     3         1.70641269745236227356D+02, 9.51190923960381458747D+02,
     4         2.06522691539642105009D-01/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation in J-fraction form
C     for  x in [3.5, 5.0]
C----------------------------------------------------------------------
CS    DATA P3/-4.55169503255094815112E+00,-1.86647123338493852582E+01,
CS   1        -7.36315669126830526754E+00,-6.68407240337696756838E+01,
CS   2         4.84507265081491452130E+01, 2.69790586735467649969E+01,
CS   3        -3.35044149820592449072E+01, 7.50964459838919612289E+00,
CS   4        -1.48432341823343965307E+00, 4.99999810924858824981E-01/
CS    DATA Q3/ 4.47820908025971749852E+01, 9.98607198039452081913E+01,
CS   1         1.40238373126149385228E+01, 3.48817758822286353588E+03,
CS   2        -9.18871385293215873406E+00, 1.24018500009917163023E+03,
CS   3        -6.88024952504512254535E+01,-2.31251575385145143070E+00,
CS   4         2.50041492369922381761E-01/
      DATA P3/-4.55169503255094815112D+00,-1.86647123338493852582D+01,
     1        -7.36315669126830526754D+00,-6.68407240337696756838D+01,
     2         4.84507265081491452130D+01, 2.69790586735467649969D+01,
     3        -3.35044149820592449072D+01, 7.50964459838919612289D+00,
     4        -1.48432341823343965307D+00, 4.99999810924858824981D-01/
      DATA Q3/ 4.47820908025971749852D+01, 9.98607198039452081913D+01,
     1         1.40238373126149385228D+01, 3.48817758822286353588D+03,
     2        -9.18871385293215873406D+00, 1.24018500009917163023D+03,
     3        -6.88024952504512254535D+01,-2.31251575385145143070D+00,
     4         2.50041492369922381761D-01/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation in J-fraction form
C     for  |x| > 5.0
C----------------------------------------------------------------------
CS    DATA P4/-8.11753647558432685797E+00,-3.84043882477454453430E+01,
CS   1        -2.23787669028751886675E+01,-2.88301992467056105854E+01,
CS   2        -5.99085540418222002197E+00,-1.13867365736066102577E+01,
CS   3        -6.52828727526980741590E+00,-4.50002293000355585708E+00,
CS   4        -2.50000000088955834952E+00, 5.00000000000000488400E-01/
CS    DATA Q4/ 2.69382300417238816428E+02, 5.04198958742465752861E+01,
CS   1         6.11539671480115846173E+01, 2.08210246935564547889E+02,
CS   2         1.97325365692316183531E+01,-1.22097010558934838708E+01,
CS   3        -6.99732735041547247161E+00,-2.49999970104184464568E+00,
CS   4         7.49999999999027092188E-01/
      DATA P4/-8.11753647558432685797D+00,-3.84043882477454453430D+01,
     1        -2.23787669028751886675D+01,-2.88301992467056105854D+01,
     2        -5.99085540418222002197D+00,-1.13867365736066102577D+01,
     3        -6.52828727526980741590D+00,-4.50002293000355585708D+00,
     4        -2.50000000088955834952D+00, 5.00000000000000488400D-01/
      DATA Q4/ 2.69382300417238816428D+02, 5.04198958742465752861D+01,
     1         6.11539671480115846173D+01, 2.08210246935564547889D+02,
     2         1.97325365692316183531D+01,-1.22097010558934838708D+01,
     3        -6.99732735041547247161D+00,-2.49999970104184464568D+00,
     4         7.49999999999027092188D-01/
C----------------------------------------------------------------------
      X = XX
      IF (ABS(X) .GT. XLARGE) THEN
            IF (ABS(X) .LE. XMAX) THEN
                  DAW = HALF / X
               ELSE
                  DAW = ZERO
            END IF
         ELSE IF (ABS(X) .LT. XSMALL) THEN
            DAW = X
         ELSE
            Y = X * X
            IF (Y .LT. SIX25) THEN
C----------------------------------------------------------------------
C  ABS(X) .LT. 2.5 
C----------------------------------------------------------------------
                  SUMP = P1(1)
                  SUMQ = Q1(1)
                  DO 100 I = 2, 10
                     SUMP = SUMP * Y + P1(I)
                     SUMQ = SUMQ * Y + Q1(I)
  100             CONTINUE
                  DAW = X * SUMP / SUMQ
               ELSE IF (Y .LT. ONE225) THEN
C----------------------------------------------------------------------
C  2.5 .LE. ABS(X) .LT. 3.5 
C----------------------------------------------------------------------
                  FRAC = ZERO
                  DO 200 I = 1, 9
  200                FRAC = Q2(I) / (P2(I) + Y + FRAC)
                  DAW = (P2(10) + FRAC) / X
               ELSE IF (Y .LT. TWO5) THEN
C----------------------------------------------------------------------
C  3.5 .LE. ABS(X) .LT. 5.0 
C---------------------------------------------------------------------
                  FRAC = ZERO
                  DO 300 I = 1, 9
  300                FRAC = Q3(I) / (P3(I) + Y + FRAC)
                  DAW = (P3(10) + FRAC) / X
               ELSE
C----------------------------------------------------------------------
C  5.0 .LE. ABS(X) .LE. XLARGE 
C------------------------------------------------------------------
                  W2 = ONE / X / X
                  FRAC = ZERO
                  DO 400 I = 1, 9
  400                FRAC = Q4(I) / (P4(I) + Y + FRAC)
                  FRAC = P4(10) + FRAC
                  DAW = (HALF + HALF * W2 * FRAC) / X
            END IF
      END IF
      RETURN
C---------- Last line of DAW ----------
      END
      SUBROUTINE CALCEI(ARG,RESULT,INT)
C----------------------------------------------------------------------
C
C This Fortran 77 packet computes the exponential integrals Ei(x),
C  E1(x), and  exp(-x)*Ei(x)  for real arguments  x  where
C
C           integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
C  Ei(x) =
C          -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
C
C  and where the first integral is a principal value integral.
C  The packet contains three function type subprograms: EI, EONE,
C  and EXPEI;  and one subroutine type subprogram: CALCEI.  The
C  calling statements for the primary entries are
C
C                 Y = EI(X),            where  X .NE. 0,
C
C                 Y = EONE(X),          where  X .GT. 0,
C  and
C                 Y = EXPEI(X),         where  X .NE. 0,
C
C  and where the entry points correspond to the functions Ei(x),
C  E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
C  is intended for internal packet use only, all computations within
C  the packet being concentrated in this routine.  The function
C  subprograms invoke CALCEI with the Fortran statement
C         CALL CALCEI(ARG,RESULT,INT)
C  where the parameter usage is as follows
C
C     Function                  Parameters for CALCEI
C       Call                 ARG             RESULT         INT
C
C      EI(X)              X .NE. 0          Ei(X)            1
C      EONE(X)            X .GT. 0         -Ei(-X)           2
C      EXPEI(X)           X .NE. 0          exp(-X)*Ei(X)    3
C
C  The main computation involves evaluation of rational Chebyshev
C  approximations published in Math. Comp. 22, 641-649 (1968), and
C  Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
C  transportable program is patterned after the machine-dependent
C  FUNPACK packet  NATSEI,  but cannot match that version for
C  efficiency or accuracy.  This version uses rational functions
C  that theoretically approximate the exponential integrals to
C  at least 18 significant decimal digits.  The accuracy achieved
C  depends on the arithmetic system, the compiler, the intrinsic
C  functions, and proper selection of the machine-dependent
C  constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta = radix for the floating-point system.
C   minexp = smallest representable power of beta.
C   maxexp = smallest power of beta that overflows.
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   XBIG = largest argument acceptable to EONE; solution to
C          equation:
C                     exp(-x)/x * (1 + 1/x) = beta ** minexp.
C   XINF = largest positive machine number; approximately
C                     beta ** maxexp
C   XMAX = largest argument acceptable to EI; solution to
C          equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
C
C     Approximate values for some important machines are:
C
C                           beta      minexp      maxexp
C
C  CRAY-1        (S.P.)       2       -8193        8191
C  Cyber 180/185
C    under NOS   (S.P.)       2        -975        1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2        -126         128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2       -1022        1024
C  IBM 3033      (D.P.)      16         -65          63
C  VAX D-Format  (D.P.)       2        -128         127
C  VAX G-Format  (D.P.)       2       -1024        1023
C
C                           XBIG       XINF       XMAX
C
C  CRAY-1        (S.P.)    5670.31  5.45E+2465   5686.21
C  Cyber 180/185
C    under NOS   (S.P.)     669.31  1.26E+322     748.28
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)      82.93  3.40E+38       93.24
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)     701.84  1.79D+308     716.35
C  IBM 3033      (D.P.)     175.05  7.23D+75      179.85
C  VAX D-Format  (D.P.)      84.30  1.70D+38       92.54
C  VAX G-Format  (D.P.)     703.22  8.98D+307     715.66
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The following table shows the types of error that may be
C  encountered in this routine and the function value supplied
C  in each case.
C
C       Error       Argument         Function values for
C                    Range         EI      EXPEI     EONE
C
C     UNDERFLOW  (-)X .GT. XBIG     0        -         0
C     OVERFLOW      X .GE. XMAX    XINF      -         -
C     ILLEGAL X       X = 0       -XINF    -XINF     XINF
C     ILLEGAL X      X .LT. 0       -        -     USE ABS(X)
C
C Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics abd Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 9, 1992
C
C----------------------------------------------------------------------
      INTEGER I,INT
CS    REAL
      DOUBLE PRECISION
     1       A,ARG,B,C,D,EXP40,E,EI,F,FOUR,FOURTY,FRAC,HALF,ONE,P,
     2       PLG,PX,P037,P1,P2,Q,QLG,QX,Q1,Q2,R,RESULT,S,SIX,SUMP,
     3       SUMQ,T,THREE,TWELVE,TWO,TWO4,W,X,XBIG,XINF,XMAX,XMX0,
     4       X0,X01,X02,X11,Y,YSQ,ZERO
      DIMENSION  A(7),B(6),C(9),D(9),E(10),F(10),P(10),Q(10),R(10),
     1   S(9),P1(10),Q1(9),P2(10),Q2(9),PLG(4),QLG(4),PX(10),QX(10)
C----------------------------------------------------------------------
C  Mathematical constants
C   EXP40 = exp(40)
C   X0 = zero of Ei
C   X01/X11 + X02 = zero of Ei to extra precision
C----------------------------------------------------------------------
CS    DATA ZERO,P037,HALF,ONE,TWO/0.0E0,0.037E0,0.5E0,1.0E0,2.0E0/,
CS   1     THREE,FOUR,SIX,TWELVE,TWO4/3.0E0,4.0E0,6.0E0,12.E0,24.0E0/,
CS   2     FOURTY,EXP40/40.0E0,2.3538526683701998541E17/,
CS   3     X01,X11,X02/381.5E0,1024.0E0,-5.1182968633365538008E-5/,
CS   4     X0/3.7250741078136663466E-1/
      DATA ZERO,P037,HALF,ONE,TWO/0.0D0,0.037D0,0.5D0,1.0D0,2.0D0/,
     1     THREE,FOUR,SIX,TWELVE,TWO4/3.0D0,4.0D0,6.0D0,12.D0,24.0D0/,
     2     FOURTY,EXP40/40.0D0,2.3538526683701998541D17/,
     3     X01,X11,X02/381.5D0,1024.0D0,-5.1182968633365538008D-5/,
     4     X0/3.7250741078136663466D-1/
C----------------------------------------------------------------------
C Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA XINF/3.40E+38/,XMAX/93.246E0/,XBIG/82.93E0/
      DATA XINF/1.79D+308/,XMAX/716.351D0/,XBIG/701.84D0/
C----------------------------------------------------------------------
C Coefficients  for -1.0 <= X < 0.0
C----------------------------------------------------------------------
CS    DATA A/1.1669552669734461083368E2, 2.1500672908092918123209E3,
CS   1       1.5924175980637303639884E4, 8.9904972007457256553251E4,
CS   2       1.5026059476436982420737E5,-1.4815102102575750838086E5,
CS   3       5.0196785185439843791020E0/
CS    DATA B/4.0205465640027706061433E1, 7.5043163907103936624165E2,
CS   1       8.1258035174768735759855E3, 5.2440529172056355429883E4,
CS   2       1.8434070063353677359298E5, 2.5666493484897117319268E5/
      DATA A/1.1669552669734461083368D2, 2.1500672908092918123209D3,
     1       1.5924175980637303639884D4, 8.9904972007457256553251D4,
     2       1.5026059476436982420737D5,-1.4815102102575750838086D5,
     3       5.0196785185439843791020D0/
      DATA B/4.0205465640027706061433D1, 7.5043163907103936624165D2,
     1       8.1258035174768735759855D3, 5.2440529172056355429883D4,
     2       1.8434070063353677359298D5, 2.5666493484897117319268D5/
C----------------------------------------------------------------------
C Coefficients for -4.0 <= X < -1.0
C----------------------------------------------------------------------
CS    DATA C/3.828573121022477169108E-1, 1.107326627786831743809E+1,
CS   1       7.246689782858597021199E+1, 1.700632978311516129328E+2,
CS   2       1.698106763764238382705E+2, 7.633628843705946890896E+1,
CS   3       1.487967702840464066613E+1, 9.999989642347613068437E-1,
CS   4       1.737331760720576030932E-8/
CS    DATA D/8.258160008564488034698E-2, 4.344836335509282083360E+0,
CS   1       4.662179610356861756812E+1, 1.775728186717289799677E+2,
CS   2       2.953136335677908517423E+2, 2.342573504717625153053E+2,
CS   3       9.021658450529372642314E+1, 1.587964570758947927903E+1,
CS   4       1.000000000000000000000E+0/
      DATA C/3.828573121022477169108D-1, 1.107326627786831743809D+1,
     1       7.246689782858597021199D+1, 1.700632978311516129328D+2,
     2       1.698106763764238382705D+2, 7.633628843705946890896D+1,
     3       1.487967702840464066613D+1, 9.999989642347613068437D-1,
     4       1.737331760720576030932D-8/
      DATA D/8.258160008564488034698D-2, 4.344836335509282083360D+0,
     1       4.662179610356861756812D+1, 1.775728186717289799677D+2,
     2       2.953136335677908517423D+2, 2.342573504717625153053D+2,
     3       9.021658450529372642314D+1, 1.587964570758947927903D+1,
     4       1.000000000000000000000D+0/
C----------------------------------------------------------------------
C Coefficients for X < -4.0
C----------------------------------------------------------------------
CS    DATA E/1.3276881505637444622987E+2,3.5846198743996904308695E+4,
CS   1       1.7283375773777593926828E+5,2.6181454937205639647381E+5,
CS   2       1.7503273087497081314708E+5,5.9346841538837119172356E+4,
CS   3       1.0816852399095915622498E+4,1.0611777263550331766871E03,
CS   4       5.2199632588522572481039E+1,9.9999999999999999087819E-1/
CS    DATA F/3.9147856245556345627078E+4,2.5989762083608489777411E+5,
CS   1       5.5903756210022864003380E+5,5.4616842050691155735758E+5,
CS   2       2.7858134710520842139357E+5,7.9231787945279043698718E+4,
CS   3       1.2842808586627297365998E+4,1.1635769915320848035459E+3,
CS   4       5.4199632588522559414924E+1,1.0E0/
      DATA E/1.3276881505637444622987D+2,3.5846198743996904308695D+4,
     1       1.7283375773777593926828D+5,2.6181454937205639647381D+5,
     2       1.7503273087497081314708D+5,5.9346841538837119172356D+4,
     3       1.0816852399095915622498D+4,1.0611777263550331766871D03,
     4       5.2199632588522572481039D+1,9.9999999999999999087819D-1/
      DATA F/3.9147856245556345627078D+4,2.5989762083608489777411D+5,
     1       5.5903756210022864003380D+5,5.4616842050691155735758D+5,
     2       2.7858134710520842139357D+5,7.9231787945279043698718D+4,
     3       1.2842808586627297365998D+4,1.1635769915320848035459D+3,
     4       5.4199632588522559414924D+1,1.0D0/
C----------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
C----------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C----------------------------------------------------------------------
C Coefficients for  0.0 < X < 6.0,
C  ratio of Chebyshev polynomials
C----------------------------------------------------------------------
CS    DATA P/-1.2963702602474830028590E01,-1.2831220659262000678155E03,
CS   1       -1.4287072500197005777376E04,-1.4299841572091610380064E06,
CS   2       -3.1398660864247265862050E05,-3.5377809694431133484800E08,
CS   3        3.1984354235237738511048E08,-2.5301823984599019348858E10,
CS   4        1.2177698136199594677580E10,-2.0829040666802497120940E11/
CS    DATA Q/ 7.6886718750000000000000E01,-5.5648470543369082846819E03,
CS   1        1.9418469440759880361415E05,-4.2648434812177161405483E06,
CS   2        6.4698830956576428587653E07,-7.0108568774215954065376E08,
CS   3        5.4229617984472955011862E09,-2.8986272696554495342658E10,
CS   4        9.8900934262481749439886E10,-8.9673749185755048616855E10/
      DATA P/-1.2963702602474830028590D01,-1.2831220659262000678155D03,
     1       -1.4287072500197005777376D04,-1.4299841572091610380064D06,
     2       -3.1398660864247265862050D05,-3.5377809694431133484800D08,
     3        3.1984354235237738511048D08,-2.5301823984599019348858D10,
     4        1.2177698136199594677580D10,-2.0829040666802497120940D11/
      DATA Q/ 7.6886718750000000000000D01,-5.5648470543369082846819D03,
     1        1.9418469440759880361415D05,-4.2648434812177161405483D06,
     2        6.4698830956576428587653D07,-7.0108568774215954065376D08,
     3        5.4229617984472955011862D09,-2.8986272696554495342658D10,
     4        9.8900934262481749439886D10,-8.9673749185755048616855D10/
C----------------------------------------------------------------------
C J-fraction coefficients for 6.0 <= X < 12.0
C----------------------------------------------------------------------
CS    DATA R/-2.645677793077147237806E00,-2.378372882815725244124E00,
CS   1       -2.421106956980653511550E01, 1.052976392459015155422E01,
CS   2        1.945603779539281810439E01,-3.015761863840593359165E01,
CS   3        1.120011024227297451523E01,-3.988850730390541057912E00,
CS   4        9.565134591978630774217E00, 9.981193787537396413219E-1/
CS    DATA S/ 1.598517957704779356479E-4, 4.644185932583286942650E00,
CS   1        3.697412299772985940785E02,-8.791401054875438925029E00,
CS   2        7.608194509086645763123E02, 2.852397548119248700147E01,
CS   3        4.731097187816050252967E02,-2.369210235636181001661E02,
CS   4        1.249884822712447891440E00/
      DATA R/-2.645677793077147237806D00,-2.378372882815725244124D00,
     1       -2.421106956980653511550D01, 1.052976392459015155422D01,
     2        1.945603779539281810439D01,-3.015761863840593359165D01,
     3        1.120011024227297451523D01,-3.988850730390541057912D00,
     4        9.565134591978630774217D00, 9.981193787537396413219D-1/
      DATA S/ 1.598517957704779356479D-4, 4.644185932583286942650D00,
     1        3.697412299772985940785D02,-8.791401054875438925029D00,
     2        7.608194509086645763123D02, 2.852397548119248700147D01,
     3        4.731097187816050252967D02,-2.369210235636181001661D02,
     4        1.249884822712447891440D00/
C----------------------------------------------------------------------
C J-fraction coefficients for 12.0 <= X < 24.0
C----------------------------------------------------------------------
CS    DATA P1/-1.647721172463463140042E00,-1.860092121726437582253E01,
CS   1        -1.000641913989284829961E01,-2.105740799548040450394E01,
CS   2        -9.134835699998742552432E-1,-3.323612579343962284333E01,
CS   3         2.495487730402059440626E01, 2.652575818452799819855E01,
CS   4        -1.845086232391278674524E00, 9.999933106160568739091E-1/
CS    DATA Q1/ 9.792403599217290296840E01, 6.403800405352415551324E01,
CS   1         5.994932325667407355255E01, 2.538819315630708031713E02,
CS   2         4.429413178337928401161E01, 1.192832423968601006985E03,
CS   3         1.991004470817742470726E02,-1.093556195391091143924E01,
CS   4         1.001533852045342697818E00/
      DATA P1/-1.647721172463463140042D00,-1.860092121726437582253D01,
     1        -1.000641913989284829961D01,-2.105740799548040450394D01,
     2        -9.134835699998742552432D-1,-3.323612579343962284333D01,
     3         2.495487730402059440626D01, 2.652575818452799819855D01,
     4        -1.845086232391278674524D00, 9.999933106160568739091D-1/
      DATA Q1/ 9.792403599217290296840D01, 6.403800405352415551324D01,
     1         5.994932325667407355255D01, 2.538819315630708031713D02,
     2         4.429413178337928401161D01, 1.192832423968601006985D03,
     3         1.991004470817742470726D02,-1.093556195391091143924D01,
     4         1.001533852045342697818D00/
C----------------------------------------------------------------------
C J-fraction coefficients for  X .GE. 24.0
C----------------------------------------------------------------------
CS    DATA P2/ 1.75338801265465972390E02,-2.23127670777632409550E02,
CS   1        -1.81949664929868906455E01,-2.79798528624305389340E01,
CS   2        -7.63147701620253630855E00,-1.52856623636929636839E01,
CS   3        -7.06810977895029358836E00,-5.00006640413131002475E00,
CS   4        -3.00000000320981265753E00, 1.00000000000000485503E00/
CS    DATA Q2/ 3.97845977167414720840E04, 3.97277109100414518365E00,
CS   1         1.37790390235747998793E02, 1.17179220502086455287E02,
CS   2         7.04831847180424675988E01,-1.20187763547154743238E01,
CS   3        -7.99243595776339741065E00,-2.99999894040324959612E00,
CS   4         1.99999999999048104167E00/
      DATA P2/ 1.75338801265465972390D02,-2.23127670777632409550D02,
     1        -1.81949664929868906455D01,-2.79798528624305389340D01,
     2        -7.63147701620253630855D00,-1.52856623636929636839D01,
     3        -7.06810977895029358836D00,-5.00006640413131002475D00,
     4        -3.00000000320981265753D00, 1.00000000000000485503D00/
      DATA Q2/ 3.97845977167414720840D04, 3.97277109100414518365D00,
     1         1.37790390235747998793D02, 1.17179220502086455287D02,
     2         7.04831847180424675988D01,-1.20187763547154743238D01,
     3        -7.99243595776339741065D00,-2.99999894040324959612D00,
     4         1.99999999999048104167D00/
C----------------------------------------------------------------------
      X = ARG
      IF (X .EQ. ZERO) THEN
            EI = -XINF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF ((X .LT. ZERO) .OR. (INT .EQ. 2)) THEN
C----------------------------------------------------------------------
C Calculate EI for negative argument or for E1.
C----------------------------------------------------------------------
            Y = ABS(X)
            IF (Y .LE. ONE) THEN
                  SUMP = A(7) * Y + A(1)
                  SUMQ = Y + B(1)
                  DO 110 I = 2, 6
                     SUMP = SUMP * Y + A(I)
                     SUMQ = SUMQ * Y + B(I)
  110             CONTINUE
                  EI = LOG(Y) - SUMP / SUMQ
                  IF (INT .EQ. 3) EI = EI * EXP(Y)
               ELSE IF (Y .LE. FOUR) THEN
                  W = ONE / Y
                  SUMP = C(1)
                  SUMQ = D(1)
                  DO 130 I = 2, 9
                     SUMP = SUMP * W + C(I)
                     SUMQ = SUMQ * W + D(I)
  130             CONTINUE
                  EI = - SUMP / SUMQ
                  IF (INT .NE. 3) EI = EI * EXP(-Y)
               ELSE
                  IF ((Y .GT. XBIG) .AND. (INT .LT. 3)) THEN
                        EI = ZERO
                     ELSE
                        W = ONE / Y
                        SUMP = E(1)
                        SUMQ = F(1)
                        DO 150 I = 2, 10
                           SUMP = SUMP * W + E(I)
                           SUMQ = SUMQ * W + F(I)
  150                   CONTINUE
                        EI = -W * (ONE - W * SUMP / SUMQ )
                        IF (INT .NE. 3) EI = EI * EXP(-Y)
                  END IF
            END IF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF (X .LT. SIX) THEN
C----------------------------------------------------------------------
C  To improve conditioning, rational approximations are expressed
C    in terms of Chebyshev polynomials for 0 <= X < 6, and in
C    continued fraction form for larger X.
C----------------------------------------------------------------------
            T = X + X
            T = T / THREE - TWO
            PX(1) = ZERO
            QX(1) = ZERO
            PX(2) = P(1)
            QX(2) = Q(1)
            DO 210 I = 2, 9
               PX(I+1) = T * PX(I) - PX(I-1) + P(I)
               QX(I+1) = T * QX(I) - QX(I-1) + Q(I)
  210       CONTINUE
            SUMP = HALF * T * PX(10) - PX(9) + P(10)
            SUMQ = HALF * T * QX(10) - QX(9) + Q(10)
            FRAC = SUMP / SUMQ
            XMX0 = (X - X01/X11) - X02
            IF (ABS(XMX0) .GE. P037) THEN
                  EI = LOG(X/X0) + XMX0 * FRAC
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
               ELSE
C----------------------------------------------------------------------
C Special approximation to  ln(X/X0)  for X close to X0
C----------------------------------------------------------------------
                  Y = XMX0 / (X + X0)
                  YSQ = Y*Y
                  SUMP = PLG(1)
                  SUMQ = YSQ + QLG(1)
                  DO 220 I = 2, 4
                     SUMP = SUMP*YSQ + PLG(I)
                     SUMQ = SUMQ*YSQ + QLG(I)
  220             CONTINUE
                  EI = (SUMP / (SUMQ*(X+X0)) + FRAC) * XMX0
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
            END IF
         ELSE IF (X .LT. TWELVE) THEN
            FRAC = ZERO
            DO 230 I = 1, 9
               FRAC = S(I) / (R(I) + X + FRAC)
  230       CONTINUE
            EI = (R(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE IF (X .LE. TWO4) THEN
            FRAC = ZERO
            DO 240 I = 1, 9
               FRAC = Q1(I) / (P1(I) + X + FRAC)
  240       CONTINUE
            EI = (P1(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE
            IF ((X .GE. XMAX) .AND. (INT .LT. 3)) THEN
                  EI = XINF
               ELSE
                  Y = ONE / X
                  FRAC = ZERO
                  DO 250 I = 1, 9
                     FRAC = Q2(I) / (P2(I) + X + FRAC)
  250             CONTINUE
                  FRAC = P2(10) + FRAC
                  EI = Y + Y * Y * FRAC
                  IF (INT .NE. 3) THEN
                        IF (X .LE. XMAX-TWO4) THEN
                              EI = EI * EXP(X)
                           ELSE
C----------------------------------------------------------------------
C Calculation reformulated to avoid premature overflow
C----------------------------------------------------------------------
                              EI = (EI * EXP(X-FOURTY)) * EXP40
                        END IF
                  END IF
            END IF
      END IF
      RESULT = EI
      RETURN
C---------- Last line of CALCEI ----------
      END
CS    REAL FUNCTION EI(X)
      DOUBLE PRECISION FUNCTION EI(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   exponential integral  Ei(x), where  x  is real.
C
C  Author: W. J. Cody
C
C  Latest modification: March 9, 1992
C
C--------------------------------------------------------------------
      INTEGER INT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      INT = 1
      CALL CALCEI(X,RESULT,INT)
      EI = RESULT
      RETURN
C---------- Last line of EI ----------
      END
CS    REAL FUNCTION EXPEI(X)
      DOUBLE PRECISION FUNCTION EXPEI(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
C   integral, and  x  is real.
C
C  Author: W. J. Cody
C
C  Latest modification: March 9, 1992
C
C--------------------------------------------------------------------
      INTEGER INT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      INT = 3
      CALL CALCEI(X,RESULT,INT)
      EXPEI = RESULT
      RETURN
C---------- Last line of EXPEI ----------
      END
CS    REAL FUNCTION EONE(X)
      DOUBLE PRECISION FUNCTION EONE(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   exponential integral E1(x), where  x  is real.
C
C  Author: W. J. Cody
C
C  Latest modification: March 9, 1992
C
C--------------------------------------------------------------------
      INTEGER INT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      INT = 2
      CALL CALCEI(X,RESULT,INT)
      EONE = RESULT
      RETURN
C---------- Last line of EONE ----------
      END
      SUBROUTINE CALERF(ARG,RESULT,JINT)
C------------------------------------------------------------------
C
C This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
C   for a real argument  x.  It contains three FUNCTION type
C   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
C   and one SUBROUTINE type subprogram, CALERF.  The calling
C   statements for the primary entries are:
C
C                   Y=ERF(X)     (or   Y=DERF(X)),
C
C                   Y=ERFC(X)    (or   Y=DERFC(X)),
C   and
C                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
C
C   The routine  CALERF  is intended for internal packet use only,
C   all computations within the packet being concentrated in this
C   routine.  The function subprograms invoke  CALERF  with the
C   statement
C
C          CALL CALERF(ARG,RESULT,JINT)
C
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALERF
C       call              ARG                  Result          JINT
C
C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
C     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
C     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
C
C   The main computation evaluates near-minimax approximations
C   from "Rational Chebyshev approximations for the error function"
C   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
C   transportable program uses rational functions that theoretically
C   approximate  erf(x)  and  erfc(x)  to at least 18 significant
C   decimal digits.  The accuracy achieved depends on the arithmetic
C   system, the compiler, the intrinsic functions, and proper
C   selection of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   XMIN   = the smallest positive floating-point number.
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   XINF   = the largest positive finite floating-point number.
C   XNEG   = the largest negative argument acceptable to ERFCX;
C            the negative of the solution to the equation
C            2*exp(x*x) = XINF.
C   XSMALL = argument below which erf(x) may be represented by
C            2*x/sqrt(pi)  and above which  x*x  will not underflow.
C            A conservative value is the largest machine number X
C            such that   1.0 + X = 1.0   to machine precision.
C   XBIG   = largest argument acceptable to ERFC;  solution to
C            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
C            W(x) = exp(-x*x)/[x*sqrt(pi)].
C   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
C            machine precision.  A conservative value is
C            1/[2*sqrt(XSMALL)]
C   XMAX   = largest acceptable argument to ERFCX; the minimum
C            of XINF and 1/[sqrt(pi)*XMIN].
C
C   Approximate values for some important machines are:
C
C                          XMIN       XINF        XNEG     XSMALL
C
C  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
C  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
C  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
C  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
C  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
C  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
C
C
C                          XBIG       XHUGE       XMAX
C
C  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
C  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
C  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
C  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
C  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
C  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ERFC = 0      for  ARG .GE. XBIG;
C
C                       ERFCX = XINF  for  ARG .LT. XNEG;
C      and
C                       ERFCX = 0     for  ARG .GE. XMAX.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 12, 1992
C
C------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1     A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI,
     2     TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL,
     3     Y,YSQ,ZERO
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
C------------------------------------------------------------------
C  Mathematical constants
C------------------------------------------------------------------
CS    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
CS   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
CS   2     SIXTEN/16.0E0/
      DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/,
     1     SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/,
     2     SIXTEN/16.0D0/
C------------------------------------------------------------------
C  Machine-dependent constants
C------------------------------------------------------------------
CS    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
CS   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
      DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/,
     1     XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
C------------------------------------------------------------------
C  Coefficients for approximation to  erf  in first interval
C------------------------------------------------------------------
CS    DATA A/3.16112374387056560E00,1.13864154151050156E02,
CS   1       3.77485237685302021E02,3.20937758913846947E03,
CS   2       1.85777706184603153E-1/
CS    DATA B/2.36012909523441209E01,2.44024637934444173E02,
CS   1       1.28261652607737228E03,2.84423683343917062E03/
      DATA A/3.16112374387056560D00,1.13864154151050156D02,
     1       3.77485237685302021D02,3.20937758913846947D03,
     2       1.85777706184603153D-1/
      DATA B/2.36012909523441209D01,2.44024637934444173D02,
     1       1.28261652607737228D03,2.84423683343917062D03/
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in second interval
C------------------------------------------------------------------
CS    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
CS   1       6.61191906371416295E01,2.98635138197400131E02,
CS   2       8.81952221241769090E02,1.71204761263407058E03,
CS   3       2.05107837782607147E03,1.23033935479799725E03,
CS   4       2.15311535474403846E-8/
CS    DATA D/1.57449261107098347E01,1.17693950891312499E02,
CS   1       5.37181101862009858E02,1.62138957456669019E03,
CS   2       3.29079923573345963E03,4.36261909014324716E03,
CS   3       3.43936767414372164E03,1.23033935480374942E03/
      DATA C/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/
      DATA D/1.57449261107098347D01,1.17693950891312499D02,
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in third interval
C------------------------------------------------------------------
CS    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
CS   1       1.25781726111229246E-1,1.60837851487422766E-2,
CS   2       6.58749161529837803E-4,1.63153871373020978E-2/
CS    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
CS   1       5.27905102951428412E-1,6.05183413124413191E-2,
CS   2       2.33520497626869185E-3/
      DATA P/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/
      DATA Q/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/
C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y .LE. THRESH) THEN
C------------------------------------------------------------------
C  Evaluate  erf  for  |X| <= 0.46875
C------------------------------------------------------------------
            YSQ = ZERO
            IF (Y .GT. XSMALL) YSQ = Y * Y
            XNUM = A(5)*YSQ
            XDEN = YSQ
            DO 20 I = 1, 3
               XNUM = (XNUM + A(I)) * YSQ
               XDEN = (XDEN + B(I)) * YSQ
   20       CONTINUE
            RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
            IF (JINT .NE. 0) RESULT = ONE - RESULT
            IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
            GO TO 800
C------------------------------------------------------------------
C  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
C------------------------------------------------------------------
         ELSE IF (Y .LE. FOUR) THEN
            XNUM = C(9)*Y
            XDEN = Y
            DO 120 I = 1, 7
               XNUM = (XNUM + C(I)) * Y
               XDEN = (XDEN + D(I)) * Y
  120       CONTINUE
            RESULT = (XNUM + C(8)) / (XDEN + D(8))
            IF (JINT .NE. 2) THEN
               YSQ = AINT(Y*SIXTEN)/SIXTEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
C------------------------------------------------------------------
C  Evaluate  erfc  for |X| > 4.0
C------------------------------------------------------------------
         ELSE
            RESULT = ZERO
            IF (Y .GE. XBIG) THEN
               IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
               IF (Y .GE. XHUGE) THEN
                  RESULT = SQRPI / Y
                  GO TO 300
               END IF
            END IF
            YSQ = ONE / (Y * Y)
            XNUM = P(6)*YSQ
            XDEN = YSQ
            DO 240 I = 1, 4
               XNUM = (XNUM + P(I)) * YSQ
               XDEN = (XDEN + Q(I)) * YSQ
  240       CONTINUE
            RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
            RESULT = (SQRPI -  RESULT) / Y
            IF (JINT .NE. 2) THEN
               YSQ = AINT(Y*SIXTEN)/SIXTEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
      END IF
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
  300 IF (JINT .EQ. 0) THEN
            RESULT = (HALF - RESULT) + HALF
            IF (X .LT. ZERO) RESULT = -RESULT
         ELSE IF (JINT .EQ. 1) THEN
            IF (X .LT. ZERO) RESULT = TWO - RESULT
         ELSE
            IF (X .LT. ZERO) THEN
               IF (X .LT. XNEG) THEN
                     RESULT = XINF
                  ELSE
                     YSQ = AINT(X*SIXTEN)/SIXTEN
                     DEL = (X-YSQ)*(X+YSQ)
                     Y = EXP(YSQ*YSQ) * EXP(DEL)
                     RESULT = (Y+Y) - RESULT
               END IF
            END IF
      END IF
  800 RETURN
C---------- Last card of CALERF ----------
      END
CS    REAL FUNCTION ERF(X)
      DOUBLE PRECISION FUNCTION DERF(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for erf(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 0
      CALL CALERF(X,RESULT,JINT)
CS    ERF = RESULT
      DERF = RESULT
      RETURN
C---------- Last card of DERF ----------
      END
CS    REAL FUNCTION ERFC(X)
      DOUBLE PRECISION FUNCTION DERFC(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 1
      CALL CALERF(X,RESULT,JINT)
CS    ERFC = RESULT
      DERFC = RESULT
      RETURN
C---------- Last card of DERFC ----------
      END
CS    REAL FUNCTION ERFCX(X)
      DOUBLE PRECISION FUNCTION DERFCX(X)
C------------------------------------------------------------------
C
C This subprogram computes approximate values for exp(x*x) * erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, March 30, 1987
C
C------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 2
      CALL CALERF(X,RESULT,JINT)
CS    ERFCX = RESULT
      DERFCX = RESULT
      RETURN
C---------- Last card of DERFCX ----------
      END
CS    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C
C Then the following machine-dependent constants must be declared 
C   in DATA statements.  IEEE values are provided as a default.
C
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: March 12, 1992
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
CS    REAL 
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
CS   1     SQRTPI/0.9189385332046727417803297E0/,
CS   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
CS   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
CS    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
CS   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
CS   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
CS   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
CS    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
CS   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
CS   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
CS   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
CS900 GAMMA = RES
  900 DGAMMA = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END
CS    REAL FUNCTION ALGAMA(X)
      DOUBLE PRECISION FUNCTION DLGAMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the LOG(GAMMA) function for a positive real
C   argument X.  Computation is based on an algorithm outlined in
C   references 1 and 2.  The program uses rational functions that
C   theoretically approximate LOG(GAMMA) to at least 18 significant
C   decimal digits.  The approximation for X > 12 is from reference
C   3, while approximations for X < 12.0 are similar to those in
C   reference 1, but are unpublished.  The accuracy achieved depends
C   on the arithmetic system, the compiler, the intrinsic functions,
C   and proper selection of the machine-dependent constants.
C
C
C*********************************************************************
C*********************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - largest argument for which LN(GAMMA(X)) is representable
C          in the machine, i.e., the solution to the equation
C                  LN(GAMMA(XBIG)) = beta**maxexp
C
C Then the following machine-dependent constants must be declared 
C   in DATA statements.  IEEE values are provided as a default.
C
C XINF   - largest machine representable floating-point number;
C          approximately beta**maxexp.
C EPS    - The smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C FRTBIG - Rough estimate of the fourth root of XBIG
C
C
C     Approximate values for some important machines are:
C
C                           beta      maxexp         XBIG
C
C CRAY-1        (S.P.)        2        8191       9.62E+2461
C Cyber 180/855
C   under NOS   (S.P.)        2        1070       1.72E+319
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)        2         128       4.08E+36
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)        2        1024       2.55D+305
C IBM 3033      (D.P.)       16          63       4.29D+73
C VAX D-Format  (D.P.)        2         127       2.05D+36
C VAX G-Format  (D.P.)        2        1023       1.28D+305
C
C
C                           XINF        EPS        FRTBIG
C
C CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
C Cyber 180/855
C   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
C IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
C VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
C VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
C
C**************************************************************
C**************************************************************
C
C Error returns
C
C  The program returns the value XINF for X .LE. 0.0 or when
C     overflow would occur.  The computation is believed to 
C     be free of underflow and overflow.
C
C
C Intrinsic functions required are:
C
C      LOG
C
C
C References:
C
C  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
C     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
C     1967, pp. 198-203.
C
C  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
C     1969.
C 
C  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
C     York, 1968.
C
C
C  Authors: W. J. Cody and L. Stoltz
C           Argonne National Laboratory
C
C  Latest modification: March 9, 1992
C
C----------------------------------------------------------------------
      INTEGER I
CS    REAL      
      DOUBLE PRECISION
     1    C,CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT68,P1,P2,P4,
     2    Q1,Q2,Q4,RES,SQRTPI,THRHAL,TWELVE,TWO,X,XBIG,XDEN,XINF,
     3    XM1,XM2,XM4,XNUM,Y,YSQ,ZERO
      DIMENSION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,ZERO/1.0E0,0.5E0,12.0E0,0.0E0/,
CS   1     FOUR,THRHAL,TWO,PNT68/4.0E0,1.5E0,2.0E0,0.6796875E0/,
CS   2     SQRTPI/0.9189385332046727417803297E0/
      DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/,
     1     FOUR,THRHAL,TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/,
     2     SQRTPI/0.9189385332046727417803297D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XINF,EPS,FRTBIG/4.08E36,3.401E38,1.19E-7,1.42E9/
      DATA XBIG,XINF,EPS,FRTBIG/2.55D305,1.79D308,2.22D-16,2.25D76/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (0.5,1.5).
C----------------------------------------------------------------------
CS    DATA D1/-5.772156649015328605195174E-1/
CS    DATA P1/4.945235359296727046734888E0,2.018112620856775083915565E2,
CS   1        2.290838373831346393026739E3,1.131967205903380828685045E4,
CS   2        2.855724635671635335736389E4,3.848496228443793359990269E4,
CS   3        2.637748787624195437963534E4,7.225813979700288197698961E3/
CS    DATA Q1/6.748212550303777196073036E1,1.113332393857199323513008E3,
CS   1        7.738757056935398733233834E3,2.763987074403340708898585E4,
CS   2        5.499310206226157329794414E4,6.161122180066002127833352E4,
CS   3        3.635127591501940507276287E4,8.785536302431013170870835E3/
      DATA D1/-5.772156649015328605195174D-1/
      DATA P1/4.945235359296727046734888D0,2.018112620856775083915565D2,
     1        2.290838373831346393026739D3,1.131967205903380828685045D4,
     2        2.855724635671635335736389D4,3.848496228443793359990269D4,
     3        2.637748787624195437963534D4,7.225813979700288197698961D3/
      DATA Q1/6.748212550303777196073036D1,1.113332393857199323513008D3,
     1        7.738757056935398733233834D3,2.763987074403340708898585D4,
     2        5.499310206226157329794414D4,6.161122180066002127833352D4,
     3        3.635127591501940507276287D4,8.785536302431013170870835D3/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     Approximation over (1.5,4.0).
C----------------------------------------------------------------------
CS    DATA D2/4.227843350984671393993777E-1/
CS    DATA P2/4.974607845568932035012064E0,5.424138599891070494101986E2,
CS   1        1.550693864978364947665077E4,1.847932904445632425417223E5,
CS   2        1.088204769468828767498470E6,3.338152967987029735917223E6,
CS   3        5.106661678927352456275255E6,3.074109054850539556250927E6/
CS    DATA Q2/1.830328399370592604055942E2,7.765049321445005871323047E3,
CS   1        1.331903827966074194402448E5,1.136705821321969608938755E6,
CS   2        5.267964117437946917577538E6,1.346701454311101692290052E7,
CS   3        1.782736530353274213975932E7,9.533095591844353613395747E6/
      DATA D2/4.227843350984671393993777D-1/
      DATA P2/4.974607845568932035012064D0,5.424138599891070494101986D2,
     1        1.550693864978364947665077D4,1.847932904445632425417223D5,
     2        1.088204769468828767498470D6,3.338152967987029735917223D6,
     3        5.106661678927352456275255D6,3.074109054850539556250927D6/
      DATA Q2/1.830328399370592604055942D2,7.765049321445005871323047D3,
     1        1.331903827966074194402448D5,1.136705821321969608938755D6,
     2        5.267964117437946917577538D6,1.346701454311101692290052D7,
     3        1.782736530353274213975932D7,9.533095591844353613395747D6/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     Approximation over (4.0,12.0).
C----------------------------------------------------------------------
CS    DATA D4/1.791759469228055000094023E0/
CS    DATA P4/1.474502166059939948905062E4,2.426813369486704502836312E6,
CS   1        1.214755574045093227939592E8,2.663432449630976949898078E9,
CS   2      2.940378956634553899906876E10,1.702665737765398868392998E11,
CS   3      4.926125793377430887588120E11,5.606251856223951465078242E11/
CS    DATA Q4/2.690530175870899333379843E3,6.393885654300092398984238E5,
CS   2        4.135599930241388052042842E7,1.120872109616147941376570E9,
CS   3      1.488613728678813811542398E10,1.016803586272438228077304E11,
CS   4      3.417476345507377132798597E11,4.463158187419713286462081E11/
      DATA D4/1.791759469228055000094023D0/
      DATA P4/1.474502166059939948905062D4,2.426813369486704502836312D6,
     1        1.214755574045093227939592D8,2.663432449630976949898078D9,
     2      2.940378956634553899906876D10,1.702665737765398868392998D11,
     3      4.926125793377430887588120D11,5.606251856223951465078242D11/
      DATA Q4/2.690530175870899333379843D3,6.393885654300092398984238D5,
     2        4.135599930241388052042842D7,1.120872109616147941376570D9,
     3      1.488613728678813811542398D10,1.016803586272438228077304D11,
     4      3.417476345507377132798597D11,4.463158187419713286462081D11/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
      Y = X
      IF ((Y .GT. ZERO) .AND. (Y .LE. XBIG)) THEN
            IF (Y .LE. EPS) THEN
                  RES = -LOG(Y)
               ELSE IF (Y .LE. THRHAL) THEN
C----------------------------------------------------------------------
C  EPS .LT. X .LE. 1.5
C----------------------------------------------------------------------
                  IF (Y .LT. PNT68) THEN
                        CORR = -LOG(Y)
                        XM1 = Y
                     ELSE
                        CORR = ZERO
                        XM1 = (Y - HALF) - HALF
                  END IF
                  IF ((Y .LE. HALF) .OR. (Y .GE. PNT68)) THEN
                        XDEN = ONE
                        XNUM = ZERO
                        DO 140 I = 1, 8
                           XNUM = XNUM*XM1 + P1(I)
                           XDEN = XDEN*XM1 + Q1(I)
  140                   CONTINUE
                        RES = CORR + (XM1 * (D1 + XM1*(XNUM/XDEN)))
                     ELSE
                        XM2 = (Y - HALF) - HALF
                        XDEN = ONE
                        XNUM = ZERO
                        DO 220 I = 1, 8
                           XNUM = XNUM*XM2 + P2(I)
                           XDEN = XDEN*XM2 + Q2(I)
  220                   CONTINUE
                        RES = CORR + XM2 * (D2 + XM2*(XNUM/XDEN))
                  END IF
               ELSE IF (Y .LE. FOUR) THEN
C----------------------------------------------------------------------
C  1.5 .LT. X .LE. 4.0
C----------------------------------------------------------------------
                  XM2 = Y - TWO
                  XDEN = ONE
                  XNUM = ZERO
                  DO 240 I = 1, 8
                     XNUM = XNUM*XM2 + P2(I)
                     XDEN = XDEN*XM2 + Q2(I)
  240             CONTINUE
                  RES = XM2 * (D2 + XM2*(XNUM/XDEN))
               ELSE IF (Y .LE. TWELVE) THEN
C----------------------------------------------------------------------
C  4.0 .LT. X .LE. 12.0
C----------------------------------------------------------------------
                  XM4 = Y - FOUR
                  XDEN = -ONE
                  XNUM = ZERO
                  DO 340 I = 1, 8
                     XNUM = XNUM*XM4 + P4(I)
                     XDEN = XDEN*XM4 + Q4(I)
  340             CONTINUE
                  RES = D4 + XM4*(XNUM/XDEN)
               ELSE 
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
                  RES = ZERO
                  IF (Y .LE. FRTBIG) THEN
                        RES = C(7)
                        YSQ = Y * Y
                        DO 450 I = 1, 6
                           RES = RES / YSQ + C(I)
  450                   CONTINUE
                  END IF
                  RES = RES/Y
                  CORR = LOG(Y)
                  RES = RES + SQRTPI - HALF*CORR
                  RES = RES + Y*(CORR-ONE)
            END IF
         ELSE
C----------------------------------------------------------------------
C  Return for bad arguments
C----------------------------------------------------------------------
            RES = XINF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
CS    ALGAMA = RES
      DLGAMA = RES
      RETURN
C ---------- Last line of DLGAMA ----------
      END
CS    REAL FUNCTION PSI(XX)
      DOUBLE PRECISION FUNCTION PSI(XX)
C----------------------------------------------------------------------
C
C This function program evaluates the logarithmic derivative of the
C   gamma function, 
C
C      psi(x) = d/dx (gamma(x)) / gamma(x) = d/dx (ln gamma(x))
C
C   for real x, where either
C
C          -xmax1 < x < -xmin (x not a negative integer), or
C            xmin < x.
C
C   The calling sequence for this function is 
C
C                  Y = PSI(X)
C
C   The main computation uses rational Chebyshev approximations
C   published in Math. Comp. 27, 123-127 (1973) by Cody, Strecok and
C   Thacher.  This transportable program is patterned after the
C   machine-dependent FUNPACK program PSI(X), but cannot match that
C   version for efficiency or accuracy.  This version uses rational
C   approximations that are theoretically accurate to 20 significant
C   decimal digits.  The accuracy achieved depends on the arithmetic
C   system, the compiler, the intrinsic functions, and proper selection
C   of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C The following machine-dependent constants must be declared in
C   DATA statements.  IEEE values are provided as a default.
C
C   XINF   = largest positive machine number
C   XMAX1  = beta ** (p-1), where beta is the radix for the
C            floating-point system, and p is the number of base-beta
C            digits in the floating-point significand.  This is an
C            upper bound on non-integral floating-point numbers, and
C            the negative of the lower bound on acceptable negative
C            arguments for PSI.  If rounding is necessary, round this
C            value down.
C   XMIN1  = the smallest in magnitude acceptable argument.  We
C            recommend XMIN1 = MAX(1/XINF,xmin) rounded up, where
C            xmin is the smallest positive floating-point number.
C   XSMALL = absolute argument below which  PI*COTAN(PI*X)  may be
C            represented by 1/X.  We recommend XSMALL < sqrt(3 eps)/pi,
C            where eps is the smallest positive number such that
C            1+eps > 1. 
C   XLARGE = argument beyond which PSI(X) may be represented by
C            LOG(X).  The solution to the equation
C               x*ln(x) = beta ** p
C            is a safe value.
C
C     Approximate values for some important machines are
C
C                        beta  p     eps     xmin       XINF  
C
C  CDC 7600      (S.P.)    2  48  7.11E-15  3.13E-294  1.26E+322
C  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
C  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
C  SUN 3/160     (D.P.)    2  53  1.11D-16  2.23D-308  1.79D+308
C  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
C                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
C   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307
C
C                         XMIN1      XMAX1     XSMALL    XLARGE
C
C  CDC 7600      (S.P.)  3.13E-294  1.40E+14  4.64E-08  9.42E+12
C  CRAY-1        (S.P.)  1.84E-2466 1.40E+14  4.64E-08  9.42E+12
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  1.18E-38   8.38E+06  1.90E-04  1.20E+06
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
C  IBM 3033      (D.P.)  1.39D-76   4.50D+15  5.80D-09  2.05D+15
C  SUN 3/160     (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
C  VAX 11/780    (S.P.)  5.89E-39   8.38E+06  1.35E-04  1.20E+06
C                (D.P.)  5.89D-39   3.60D+16  2.05D-09  2.05D+15
C   (G Format)   (D.P.)  1.12D-308  4.50D+15  5.80D-09  2.71D+14
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns XINF for  X < -XMAX1, for X zero or a negative
C    integer, or when X lies in (-XMIN1, 0), and returns -XINF
C    when X lies in (0, XMIN1).
C
C Intrinsic functions required are:
C
C     ABS, AINT, DBLE, INT, LOG, REAL, TAN
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division 
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 14, 1992
C
C----------------------------------------------------------------------
      INTEGER I,N,NQ
CS    REAL
      DOUBLE PRECISION
     1   AUG,CONV,DEN,FOUR,FOURTH,HALF,ONE,P1,P2,PIOV4,Q1,Q2,
     2   SGN,THREE,XLARGE,UPPER,W,X,XINF,XMAX1,XMIN1,XSMALL,X01,
     3   X01D,X02,XX,Z,ZERO
      DIMENSION P1(9),P2(7),Q1(8),Q2(6)
C----------------------------------------------------------------------
C  Mathematical constants.  PIOV4 = pi / 4
C----------------------------------------------------------------------
CS    DATA ZERO,FOURTH,HALF,ONE/0.0E0,0.25E0,0.5E0,1.0E0/
CS    DATA THREE,FOUR/3.0E0,4.0E0/,PIOV4/7.8539816339744830962E-01/
      DATA ZERO,FOURTH,HALF,ONE/0.0D0,0.25D0,0.5D0,1.0D0/
      DATA THREE,FOUR/3.0D0,4.0D0/,PIOV4/7.8539816339744830962D-01/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA XINF/3.40E+38/, XMIN1/1.18E-38/, XMAX1/8.38E+06/,
CS   1     XSMALL/1.90E-04/, XLARGE/1.20E+06/
      DATA XINF/1.79D+308/, XMIN1/2.23D-308/, XMAX1/4.50D+15/,
     1     XSMALL/5.80D-09/, XLARGE/2.71D+14/
C----------------------------------------------------------------------
C  Zero of psi(x)
C----------------------------------------------------------------------
CS    DATA X01/187.0E0/,X01D/128.0E0/,X02/6.9464496836234126266E-04/
      DATA X01/187.0D0/,X01D/128.0D0/,X02/6.9464496836234126266D-04/
C----------------------------------------------------------------------
C  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
C----------------------------------------------------------------------
CS    DATA P1/4.5104681245762934160E-03,5.4932855833000385356E+00,
CS   1        3.7646693175929276856E+02,7.9525490849151998065E+03,
CS   2        7.1451595818951933210E+04,3.0655976301987365674E+05,
CS   3        6.3606997788964458797E+05,5.8041312783537569993E+05,
CS   4        1.6585695029761022321E+05/
CS    DATA Q1/9.6141654774222358525E+01,2.6287715790581193330E+03,
CS   1        2.9862497022250277920E+04,1.6206566091533671639E+05,
CS   2        4.3487880712768329037E+05,5.4256384537269993733E+05,
CS   3        2.4242185002017985252E+05,6.4155223783576225996E-08/
      DATA P1/4.5104681245762934160D-03,5.4932855833000385356D+00,
     1        3.7646693175929276856D+02,7.9525490849151998065D+03,
     2        7.1451595818951933210D+04,3.0655976301987365674D+05,
     3        6.3606997788964458797D+05,5.8041312783537569993D+05,
     4        1.6585695029761022321D+05/
      DATA Q1/9.6141654774222358525D+01,2.6287715790581193330D+03,
     1        2.9862497022250277920D+04,1.6206566091533671639D+05,
     2        4.3487880712768329037D+05,5.4256384537269993733D+05,
     3        2.4242185002017985252D+05,6.4155223783576225996D-08/
C----------------------------------------------------------------------
C  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x) 
C     for  x > 3.0
C----------------------------------------------------------------------
CS    DATA P2/-2.7103228277757834192E+00,-1.5166271776896121383E+01,
CS   1        -1.9784554148719218667E+01,-8.8100958828312219821E+00,
CS   2        -1.4479614616899842986E+00,-7.3689600332394549911E-02,
CS   3        -6.5135387732718171306E-21/
CS    DATA Q2/ 4.4992760373789365846E+01, 2.0240955312679931159E+02,
CS   1         2.4736979003315290057E+02, 1.0742543875702278326E+02,
CS   2         1.7463965060678569906E+01, 8.8427520398873480342E-01/
      DATA P2/-2.7103228277757834192D+00,-1.5166271776896121383D+01,
     1        -1.9784554148719218667D+01,-8.8100958828312219821D+00,
     2        -1.4479614616899842986D+00,-7.3689600332394549911D-02,
     3        -6.5135387732718171306D-21/
      DATA Q2/ 4.4992760373789365846D+01, 2.0240955312679931159D+02,
     1         2.4736979003315290057D+02, 1.0742543875702278326D+02,
     2         1.7463965060678569906D+01, 8.8427520398873480342D-01/
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      X = XX
      W = ABS(X)
      AUG = ZERO
C----------------------------------------------------------------------
C  Check for valid arguments, then branch to appropriate algorithm
C----------------------------------------------------------------------
      IF ((-X .GE. XMAX1) .OR. (W .LT. XMIN1)) THEN
            GO TO 410
         ELSE IF (X .GE. HALF) THEN
            GO TO 200
C----------------------------------------------------------------------
C  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
C     Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.  
C----------------------------------------------------------------------
         ELSE IF (W .LE. XSMALL) THEN
            AUG = -ONE / X
            GO TO 150
      END IF
C----------------------------------------------------------------------
C  Argument reduction for cot
C----------------------------------------------------------------------
      IF (X .LT. ZERO) THEN
            SGN = PIOV4
         ELSE
            SGN = -PIOV4
      END IF
      W = W - AINT(W)
      NQ = INT(W * FOUR)
      W = FOUR * (W - CONV(NQ) * FOURTH)
C----------------------------------------------------------------------
C  W is now related to the fractional part of  4.0 * X.
C     Adjust argument to correspond to values in the first
C     quadrant and determine the sign.
C----------------------------------------------------------------------
      N = NQ / 2
      IF ((N+N) .NE. NQ) W = ONE - W
      Z = PIOV4 * W
      IF (MOD(N,2) .NE. 0) SGN = - SGN
C----------------------------------------------------------------------
C  determine the final value for  -pi * cotan(pi*x)
C----------------------------------------------------------------------
      N = (NQ + 1) / 2
      IF (MOD(N,2) .EQ. 0) THEN
C----------------------------------------------------------------------
C  Check for singularity
C----------------------------------------------------------------------
            IF (Z .EQ. ZERO) GO TO 410
            AUG = SGN * (FOUR / TAN(Z))
         ELSE
            AUG = SGN * (FOUR * TAN(Z))
      END IF
  150 X = ONE - X
  200 IF (X .GT. THREE) GO TO 300
C----------------------------------------------------------------------
C  0.5 <= X <= 3.0
C----------------------------------------------------------------------
      DEN = X
      UPPER = P1(1) * X
      DO 210 I = 1, 7
         DEN = (DEN + Q1(I)) * X
         UPPER = (UPPER + P1(I+1)) * X
  210 CONTINUE
      DEN = (UPPER + P1(9)) / (DEN + Q1(8))
      X = (X-X01/X01D) - X02
      PSI = DEN * X + AUG
      GO TO 500
C----------------------------------------------------------------------
C  3.0 < X 
C----------------------------------------------------------------------
  300 IF (X .LT. XLARGE) THEN
         W = ONE / (X * X)
         DEN = W
         UPPER = P2(1) * W
         DO 310 I = 1, 5
            DEN = (DEN + Q2(I)) * W
            UPPER = (UPPER + P2(I+1)) * W
  310    CONTINUE
         AUG = (UPPER + P2(7)) / (DEN + Q2(6)) - HALF / X + AUG
      END IF
      PSI = AUG + LOG(X)
      GO TO 500
C----------------------------------------------------------------------
C  Error return
C----------------------------------------------------------------------
  410 PSI = XINF
      IF (X .GT. ZERO) PSI = -XINF
  500 RETURN
C---------- Last card of PSI ----------
      END
      SUBROUTINE RIBESL(X,ALPHA,NB,IZE,B,NCALC)
C-------------------------------------------------------------------
C
C  This routine calculates Bessel functions I SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA,
C  with or without exponential scaling.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         I's or exponentially scaled I's (I*EXP(-X))
C         are to be calculated.  If I's are to be calculated,
C         X must be less than EXPARG (see below).
C ALPHA - Working precision fractional part of order for which
C         I's or exponentially scaled I's (I*EXP(-X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated,
C         and 2 if exponentially scaled I's are to be calculated.
C B     - Working precision output vector of length NB.  If the routine
C         terminates normally (NCALC=NB), the vector B contains the
C         functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
C         corresponding exponentially scaled functions.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector B, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   it     = Number of bits in the mantissa of a working precision
C            variable
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear
C            in mind that if ABS(X)=N, then at least N iterations
C            of the backward recursion will be executed.  The value
C            of 10.0 ** 4 is used on every machine.
C   EXPARG = Largest working precision argument that the library
C            EXP routine can handle and upper limit on the
C            magnitude of X when IZE=1; approximately
C            LOG(beta**maxexp)
C
C
C     Approximate values for some important machines are:
C
C                        beta       minexp      maxexp       it
C
C  CRAY-1        (S.P.)    2        -8193        8191        48
C  Cyber 180/855
C    under NOS   (S.P.)    2         -975        1070        48
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2         -126         128        24
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2        -1022        1024        53
C  IBM 3033      (D.P.)   16          -65          63        14
C  VAX           (S.P.)    2         -128         127        24
C  VAX D-Format  (D.P.)    2         -128         127        56
C  VAX G-Format  (D.P.)    2        -1024        1023        53
C
C
C                        NSIG       ENTEN       ENSIG      RTNSIG
C
C CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
C Cyber 180/855
C   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
C IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
C VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
C VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
C VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
C
C
C                         ENMTEN      XLARGE   EXPARG
C
C CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
C Cyber 180/855
C   under NOS   (S.P.)   1.25E-293    1.0E+4     741
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
C IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
C VAX           (S.P.)   1.17E-38     1.0E+4      88
C VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
C VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error,  NCALC .NE. NB, and not all I's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. 0:  An argument is out of range. For example,
C     NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. EXPARG.
C     In this case, the B-vector is not calculated, and NCALC is
C     set to MIN0(NB,0)-1 so that NCALC .NE. NB.
C
C  NB .GT. NCALC .GT. 0: Not all requested function values could
C     be calculated accurately.  This usually occurs because NB is
C     much larger than ABS(X).  In this case, B(N) is calculated
C     to the desired accuracy for N .LE. NCALC, but precision
C     is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C     for N .GT. NCALC (because it is too small to be represented),
C     and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C     significant figures of B(N) can be trusted.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, DGAMMA, GAMMA, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program is based on a program written by David J.
C  Sookne (2) that computes values of the Bessel functions J or
C  I of real argument and integer order.  Modifications include
C  the restriction of the computation to the I Bessel function
C  of non-negative real argument, the extension of the computation
C  to arbitrary positive order, the inclusion of optional
C  exponential scaling, and the elimination of most underflow.
C  An earlier version was published in (3).
C
C References: "A Note on Backward Recurrence Algorithms," Olver,
C              F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C              pp 941-947.
C
C             "Bessel Functions of Real Argument and Integer Order,"
C              Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C              125-132.
C
C             "ALGORITHM 597, Sequence of Modified Bessel Functions
C              of the First Kind," Cody, W. J., Trans. Math. Soft.,
C              1983, pp. 242-245.
C
C  Latest modification: March 14, 1992
C
C  Modified by: W. J. Cody and L. Stoltz
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C-------------------------------------------------------------------
      INTEGER IZE,K,L,MAGX,N,NB,NBMX,NCALC,NEND,NSIG,NSTART
CS    REAL              GAMMA,XX,
      DOUBLE PRECISION DGAMMA,XX,
     1 ALPHA,B,CONST,CONV,EM,EMPAL,EMP2AL,EN,ENMTEN,ENSIG,
     2 ENTEN,EXPARG,FUNC,HALF,HALFX,ONE,P,PLAST,POLD,PSAVE,PSAVEL,
     3 RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,X,XLARGE,ZERO
      DIMENSION B(NB)
C-------------------------------------------------------------------
C  Mathematical constants
C-------------------------------------------------------------------
CS    DATA ONE,TWO,ZERO,HALF,CONST/1.0E0,2.0E0,0.0E0,0.5E0,1.585E0/
      DATA ONE,TWO,ZERO,HALF,CONST/1.0D0,2.0D0,0.0D0,0.5D0,1.585D0/
C-------------------------------------------------------------------
C  Machine-dependent parameters
C-------------------------------------------------------------------
CS    DATA NSIG,XLARGE,EXPARG /8,1.0E4,88.0E0/
CS    DATA ENTEN,ENSIG,RTNSIG/1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN/4.7E-38/
      DATA NSIG,XLARGE,EXPARG /16,1.0D4,709.0D0/
      DATA ENTEN,ENSIG,RTNSIG/1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN/8.9D-308/
C-------------------------------------------------------------------
C  Statement functions for conversion
C-------------------------------------------------------------------
CS    CONV(N) = REAL(N)
CS    FUNC(XX) = GAMMA(XX)
      CONV(N) = DBLE(N)
      FUNC(XX) = DGAMMA(XX)
C-------------------------------------------------------------------
C Check for X, NB, OR IZE out of range.
C-------------------------------------------------------------------
      IF ((NB.GT.0) .AND. (X .GE. ZERO) .AND.
     1    (ALPHA .GE. ZERO) .AND. (ALPHA .LT. ONE) .AND.
     2    (((IZE .EQ. 1) .AND. (X .LE. EXPARG)) .OR.
     3     ((IZE .EQ. 2) .AND. (X .LE. XLARGE)))) THEN
C-------------------------------------------------------------------
C Use 2-term ascending series for small X
C-------------------------------------------------------------------
            NCALC = NB
            MAGX = INT(X)
            IF (X .GE. RTNSIG) THEN
C-------------------------------------------------------------------
C Initialize the forward sweep, the P-sequence of Olver
C-------------------------------------------------------------------
                  NBMX = NB-MAGX
                  N = MAGX+1
                  EN = CONV(N+N) + (ALPHA+ALPHA)
                  PLAST = ONE
                  P = EN / X
C-------------------------------------------------------------------
C Calculate general significance test
C-------------------------------------------------------------------
                  TEST = ENSIG + ENSIG
                  IF (2*MAGX .GT. 5*NSIG) THEN
                        TEST = SQRT(TEST*P)
                     ELSE
                        TEST = TEST / CONST**MAGX
                  END IF
                  IF (NBMX .GE. 3) THEN
C-------------------------------------------------------------------
C Calculate P-sequence until N = NB-1.  Check for possible overflow.
C-------------------------------------------------------------------
                     TOVER = ENTEN / ENSIG
                     NSTART = MAGX+2
                     NEND = NB - 1
                     DO 100 K = NSTART, NEND
                        N = K
                        EN = EN + TWO
                        POLD = PLAST
                        PLAST = P
                        P = EN * PLAST/X + POLD
                        IF (P .GT. TOVER) THEN
C-------------------------------------------------------------------
C To avoid overflow, divide P-sequence by TOVER.  Calculate
C P-sequence until ABS(P) .GT. 1.
C-------------------------------------------------------------------
                           TOVER = ENTEN
                           P = P / TOVER
                           PLAST = PLAST / TOVER
                           PSAVE = P
                           PSAVEL = PLAST
                           NSTART = N + 1
   60                      N = N + 1
                              EN = EN + TWO
                              POLD = PLAST
                              PLAST = P
                              P = EN * PLAST/X + POLD
                           IF (P .LE. ONE) GO TO 60
                           TEMPB = EN / X
C-------------------------------------------------------------------
C Calculate backward test, and find NCALC, the highest N
C such that the test is passed.
C-------------------------------------------------------------------
                           TEST = POLD*PLAST / ENSIG
                           TEST = TEST*(HALF-HALF/(TEMPB*TEMPB))
                           P = PLAST * TOVER
                           N = N - 1
                           EN = EN - TWO
                           NEND = MIN0(NB,N)
                           DO 80 L = NSTART, NEND
                              NCALC = L
                              POLD = PSAVEL
                              PSAVEL = PSAVE
                              PSAVE = EN * PSAVEL/X + POLD
                              IF (PSAVE*PSAVEL .GT. TEST) GO TO 90
   80                      CONTINUE
                           NCALC = NEND + 1
   90                      NCALC = NCALC - 1
                           GO TO 120
                        END IF
  100                CONTINUE
                     N = NEND
                     EN = CONV(N+N) + (ALPHA+ALPHA)
C-------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C-------------------------------------------------------------------
                     TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
                  END IF
C-------------------------------------------------------------------
C Calculate P-sequence until significance test passed.
C-------------------------------------------------------------------
  110             N = N + 1
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN * PLAST/X + POLD
                  IF (P .LT. TEST) GO TO 110
C-------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C-------------------------------------------------------------------
  120             N = N + 1
                  EN = EN + TWO
                  TEMPB = ZERO
                  TEMPA = ONE / P
                  EM = CONV(N) - ONE
                  EMPAL = EM + ALPHA
                  EMP2AL = (EM - ONE) + (ALPHA + ALPHA)
                  SUM = TEMPA * EMPAL * EMP2AL / EM
                  NEND = N - NB
                  IF (NEND .LT. 0) THEN
C-------------------------------------------------------------------
C N .LT. NB, so store B(N) and set higher orders to zero.
C-------------------------------------------------------------------
                        B(N) = TEMPA
                        NEND = -NEND
                        DO 130 L = 1, NEND
  130                      B(N+L) = ZERO
                     ELSE
                        IF (NEND .GT. 0) THEN
C-------------------------------------------------------------------
C Recur backward via difference equation, calculating (but
C not storing) B(N), until N = NB.
C-------------------------------------------------------------------
                           DO 140 L = 1, NEND
                              N = N - 1
                              EN = EN - TWO
                              TEMPC = TEMPB
                              TEMPB = TEMPA
                              TEMPA = (EN*TEMPB) / X + TEMPC
                              EM = EM - ONE
                              EMP2AL = EMP2AL - ONE
                              IF (N .EQ. 1) GO TO 150
                              IF (N .EQ. 2) EMP2AL = ONE
                              EMPAL = EMPAL - ONE
                              SUM = (SUM + TEMPA*EMPAL) * EMP2AL / EM
  140                      CONTINUE
                        END IF
C-------------------------------------------------------------------
C Store B(NB)
C-------------------------------------------------------------------
  150                   B(N) = TEMPA
                        IF (NB .LE. 1) THEN
                           SUM = (SUM + SUM) + TEMPA
                           GO TO 230
                        END IF
C-------------------------------------------------------------------
C Calculate and Store B(NB-1)
C-------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N)  = (EN*TEMPA) / X + TEMPB
                        IF (N .EQ. 1) GO TO 220
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
                  END IF
                  NEND = N - 2
                  IF (NEND .GT. 0) THEN
C-------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C-------------------------------------------------------------------
                     DO 200 L = 1, NEND
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*B(N+1)) / X +B(N+2)
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
  200                CONTINUE
                  END IF
C-------------------------------------------------------------------
C Calculate B(1)
C-------------------------------------------------------------------
                  B(1) = TWO*EMPAL*B(2) / X + B(3)
  220             SUM = (SUM + SUM) + B(1)
C-------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C-------------------------------------------------------------------
  230             IF (ALPHA .NE. ZERO)
     1               SUM = SUM * FUNC(ONE+ALPHA) * (X*HALF)**(-ALPHA)
                  IF (IZE .EQ. 1) SUM = SUM * EXP(-X)
                  TEMPA = ENMTEN
                  IF (SUM .GT. ONE) TEMPA = TEMPA * SUM
                  DO 260 N = 1, NB
                     IF (B(N) .LT. TEMPA) B(N) = ZERO
                     B(N) = B(N) / SUM
  260             CONTINUE
                  RETURN
C-------------------------------------------------------------------
C Two-term ascending series for small X.
C-------------------------------------------------------------------
               ELSE
                  TEMPA = ONE
                  EMPAL = ONE + ALPHA
                  HALFX = ZERO
                  IF (X .GT. ENMTEN) HALFX = HALF * X
                  IF (ALPHA .NE. ZERO) TEMPA = HALFX**ALPHA /FUNC(EMPAL)
                  IF (IZE .EQ. 2) TEMPA = TEMPA * EXP(-X)
                  TEMPB = ZERO
                  IF ((X+ONE) .GT. ONE) TEMPB = HALFX * HALFX
                  B(1) = TEMPA + TEMPA*TEMPB / EMPAL
                  IF ((X .NE. ZERO) .AND. (B(1) .EQ. ZERO)) NCALC = 0
                  IF (NB .GT. 1) THEN
                     IF (X .EQ. ZERO) THEN
                           DO 310 N = 2, NB
                              B(N) = ZERO
  310                      CONTINUE
                        ELSE
C-------------------------------------------------------------------
C Calculate higher-order functions.
C-------------------------------------------------------------------
                           TEMPC = HALFX
                           TOVER = (ENMTEN + ENMTEN) / X
                           IF (TEMPB .NE. ZERO) TOVER = ENMTEN / TEMPB
                           DO 340 N = 2, NB
                              TEMPA = TEMPA / EMPAL
                              EMPAL = EMPAL + ONE
                              TEMPA = TEMPA * TEMPC
                              IF (TEMPA .LE. TOVER*EMPAL) TEMPA = ZERO
                              B(N) = TEMPA + TEMPA*TEMPB / EMPAL
                              IF ((B(N) .EQ. ZERO) .AND. (NCALC .GT. N))
     1                             NCALC = N-1
  340                      CONTINUE
                     END IF
                  END IF
            END IF
         ELSE
            NCALC = MIN0(NB,0)-1
      END IF
      RETURN
C---------- Last line of RIBESL ----------
      END
      SUBROUTINE RJBESL(X, ALPHA, NB, B, NCALC)
C---------------------------------------------------------------------
C This routine calculates Bessel functions J sub(N+ALPHA) (X)
C   for non-negative argument X, and non-negative order N+ALPHA.
C
C
C  Explanation of variables in the calling sequence.
C
C   X     - working precision non-negative real argument for which
C           J's are to be calculated.
C   ALPHA - working precision fractional part of order for which
C           J's or exponentially scaled J'r (J*exp(X)) are
C           to be calculated.  0 <= ALPHA < 1.0.
C   NB  - integer number of functions to be calculated, NB > 0.
C           The first function calculated is of order ALPHA, and the
C           last is of order (NB - 1 + ALPHA).
C   B  - working precision output vector of length NB.  If RJBESL
C           terminates normally (NCALC=NB), the vector B contains the
C           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
C           corresponding exponentially scaled functions.
C   NCALC - integer output variable indicating possible errors.
C           Before using the vector B, the user should check that
C           NCALC=NB, i.e., all orders have been calculated to
C           the desired accuracy.  See Error Returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
C            then at least N iterations of the backward recursion
C            will be executed.  The value of 10.0 ** 4 is used on
C            every machine.
C
C
C     Approximate values for some important machines are:
C
C
C                            it    NSIG    ENTEN       ENSIG
C
C   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
C   Cyber 180/855
C     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
C   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
C   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
C   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
C   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16
C
C
C                           RTNSIG      ENMTEN      XLARGE
C
C   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
C   Cyber 180/855
C     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
C   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
C   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
C   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
C   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4
C
C*******************************************************************
C*******************************************************************
C
C  Error returns
C
C    In case of an error,  NCALC .NE. NB, and not all J's are
C    calculated to the desired accuracy.
C
C    NCALC .LT. 0:  An argument is out of range. For example,
C       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
C       In this case, B(1) is set to zero, the remainder of the
C       B-vector is not calculated, and NCALC is set to
C       MIN(NB,0)-1 so that NCALC .NE. NB.
C
C    NB .GT. NCALC .GT. 0: Not all requested function values could
C       be calculated accurately.  This usually occurs because NB is
C       much larger than ABS(X).  In this case, B(N) is calculated
C       to the desired accuracy for N .LE. NCALC, but precision
C       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C       for N .GT. NCALC (because it is too small to be represented),
C       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C       significant figures of B(N) can be trusted.
C
C
C  Intrinsic and other functions required are:
C
C     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,
C
C     REAL, SIN, SQRT
C
C
C  Acknowledgement
C
C   This program is based on a program written by David J. Sookne
C   (2) that computes values of the Bessel functions J or I of real
C   argument and integer order.  Modifications include the restriction
C   of the computation to the J Bessel function of non-negative real
C   argument, the extension of the computation to arbitrary positive
C   order, and the elimination of most underflow.
C
C  References: "A Note on Backward Recurrence Algorithms," Olver,
C               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C               pp 941-947.
C
C              "Bessel Functions of Real Argument and Integer Order,"
C               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C               125-132.
C  Latest modification: March 15, 1992
C
C  Author: W. J. Cody
C          Applied Mathematics Division
C          Argonne National Laboratory
C          Argonne, IL  60439
C
C---------------------------------------------------------------------
      INTEGER I,J,K,L,M,MAGX,N,NB,NBMX,NCALC,NEND,NSTART
CS    REAL               GAMMA,XX,
      DOUBLE PRECISION  DGAMMA,XX,
     1 ALPHA,ALPEM,ALP2EM,B,CAPP,CAPQ,CONV,EIGHTH,EM,EN,ENMTEN,ENSIG,
     2 ENTEN,FACT,FOUR,FUNC,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST,
     3 POLD,PSAVE,PSAVEL,RTNSIG,S,SUM,T,T1,TEMPA,TEMPB,TEMPC,TEST,
     4 THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1,TWOPI2,X,XC,XIN,XK,XLARGE,
     5 XM,VCOS,VSIN,Z,ZERO
      DIMENSION B(NB), FACT(25)
C---------------------------------------------------------------------
C  Mathematical constants
C
C   PI2    - 2 / PI
C   TWOPI1 - first few significant digits of 2 * PI
C   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
C            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
C---------------------------------------------------------------------
CS    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,
CS   1 1.935307179586476925286767E-3/
CS    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/
CS    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/
CS    DATA ONE30, THREE5 /130.0E0,35.0E0/
      DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535D0,6.28125D0,
     1 1.935307179586476925286767D-3/
      DATA ZERO, EIGHTH, HALF, ONE /0.0D0,0.125D0,0.5D0,1.0D0/
      DATA TWO, THREE, FOUR, TWOFIV /2.0D0,3.0D0,4.0D0,25.0D0/
      DATA ONE30, THREE5 /130.0D0,35.0D0/
C---------------------------------------------------------------------
C  Machine-dependent parameters
C---------------------------------------------------------------------
CS    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN, XLARGE /4.70E-38,1.0E4/
      DATA ENTEN, ENSIG, RTNSIG /1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN, XLARGE /8.90D-308,1.0D4/
C---------------------------------------------------------------------
C     Factorial(N)
C---------------------------------------------------------------------
CS    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,
CS   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,
CS   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,
CS   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,
CS   4 5.109094217170944E19,1.12400072777760768E21,
CS   5 2.585201673888497664E22,6.2044840173323943936E23/
      DATA FACT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,
     1 4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,
     2 8.71782912D10,1.307674368D12,2.0922789888D13,3.55687428096D14,
     3 6.402373705728D15,1.21645100408832D17,2.43290200817664D18,
     4 5.109094217170944D19,1.12400072777760768D21,
     5 2.585201673888497664D22,6.2044840173323943936D23/
C---------------------------------------------------------------------
C Statement functions for conversion and the gamma function.
C---------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CS    FUNC(XX) = GAMMA(XX)
      CONV(I) = DBLE(I)
      FUNC(XX) = DGAMMA(XX)
C---------------------------------------------------------------------
C Check for out of range arguments.
C---------------------------------------------------------------------
      MAGX = INT(X)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE)
     1       .AND. (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE))
     2   THEN
C---------------------------------------------------------------------
C Initialize result array to zero.
C---------------------------------------------------------------------
            NCALC = NB
            DO 20 I=1,NB
              B(I) = ZERO
   20       CONTINUE
C---------------------------------------------------------------------
C Branch to use 2-term ascending series for small X and asymptotic
C form for large X when NB is not too large.
C---------------------------------------------------------------------
            IF (X.LT.RTNSIG) THEN
C---------------------------------------------------------------------
C Two-term ascending series for small X.
C---------------------------------------------------------------------
               TEMPA = ONE
               ALPEM = ONE + ALPHA
               HALFX = ZERO
               IF (X.GT.ENMTEN) HALFX = HALF*X
               IF (ALPHA.NE.ZERO)
     1            TEMPA = HALFX**ALPHA/(ALPHA*FUNC(ALPHA))
               TEMPB = ZERO
               IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX
               B(1) = TEMPA + TEMPA*TEMPB/ALPEM
               IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
               IF (NB .NE. 1) THEN
                  IF (X .LE. ZERO) THEN
                        DO 30 N=2,NB
                          B(N) = ZERO
   30                   CONTINUE
                     ELSE
C---------------------------------------------------------------------
C Calculate higher order functions.
C---------------------------------------------------------------------
                        TEMPC = HALFX
                        TOVER = (ENMTEN+ENMTEN)/X
                        IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
                        DO 50 N=2,NB
                          TEMPA = TEMPA/ALPEM
                          ALPEM = ALPEM + ONE
                          TEMPA = TEMPA*TEMPC
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM
                          IF ((B(N).EQ.ZERO) .AND. (NCALC.GT.N))
     1                       NCALC = N-1
   50                   CONTINUE
                  END IF
               END IF
            ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN
C---------------------------------------------------------------------
C Asymptotic series for X .GT. 21.0.
C---------------------------------------------------------------------
               XC = SQRT(PI2/X)
               XIN = (EIGHTH/X)**2
               M = 11
               IF (X.GE.THREE5) M = 8
               IF (X.GE.ONE30) M = 4
               XM = FOUR*CONV(M)
C---------------------------------------------------------------------
C Argument reduction for SIN and COS routines.
C---------------------------------------------------------------------
               T = AINT(X/(TWOPI1+TWOPI2)+HALF)
               Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2
               VSIN = SIN(Z)
               VCOS = COS(Z)
               GNU = ALPHA + ALPHA
               DO 80 I=1,2
                 S = ((XM-ONE)-GNU)*((XM-ONE)+GNU)*XIN*HALF
                 T = (GNU-(XM-THREE))*(GNU+(XM-THREE))
                 CAPP = S*T/FACT(2*M+1)
                 T1 = (GNU-(XM+ONE))*(GNU+(XM+ONE))
                 CAPQ = S*T1/FACT(2*M+2)
                 XK = XM
                 K = M + M
                 T1 = T
                 DO 70 J=2,M
                   XK = XK - FOUR
                   S = ((XK-ONE)-GNU)*((XK-ONE)+GNU)
                   T = (GNU-(XK-THREE))*(GNU+(XK-THREE))
                   CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN
                   CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN
                   K = K - 2
                   T1 = T
   70            CONTINUE
                 CAPP = CAPP + ONE
                 CAPQ = (CAPQ+ONE)*(GNU*GNU-ONE)*(EIGHTH/X)
                 B(I) = XC*(CAPP*VCOS-CAPQ*VSIN)
                 IF (NB.EQ.1) GO TO 300
                 T = VSIN
                 VSIN = -VCOS
                 VCOS = T
                 GNU = GNU + TWO
   80         CONTINUE
C---------------------------------------------------------------------
C If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1
C---------------------------------------------------------------------
               IF (NB .GT. 2) THEN
                  GNU = ALPHA + ALPHA + TWO
                  DO 90 J=3,NB
                    B(J) = GNU*B(J-1)/X - B(J-2)
                    GNU = GNU + TWO
   90             CONTINUE
               END IF
C---------------------------------------------------------------------
C Use recurrence to generate results.  First initialize the
C calculation of P*S.
C---------------------------------------------------------------------
            ELSE
               NBMX = NB - MAGX
               N = MAGX + 1
               EN = CONV(N+N) + (ALPHA+ALPHA)
               PLAST = ONE
               P = EN/X
C---------------------------------------------------------------------
C Calculate general significance test.
C---------------------------------------------------------------------
               TEST = ENSIG + ENSIG
               IF (NBMX .GE. 3) THEN
C---------------------------------------------------------------------
C Calculate P*S until N = NB-1.  Check for possible overflow.
C---------------------------------------------------------------------
                  TOVER = ENTEN/ENSIG
                  NSTART = MAGX + 2
                  NEND = NB - 1
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA)
                  DO 130 K=NSTART,NEND
                     N = K
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN*PLAST/X - POLD
                     IF (P.GT.TOVER) THEN
C---------------------------------------------------------------------
C To avoid overflow, divide P*S by TOVER.  Calculate P*S until
C ABS(P) .GT. 1.
C---------------------------------------------------------------------
                        TOVER = ENTEN
                        P = P/TOVER
                        PLAST = PLAST/TOVER
                        PSAVE = P
                        PSAVEL = PLAST
                        NSTART = N + 1
  100                   N = N + 1
                           EN = EN + TWO
                           POLD = PLAST
                           PLAST = P
                           P = EN*PLAST/X - POLD
                        IF (P.LE.ONE) GO TO 100
                        TEMPB = EN/X
C---------------------------------------------------------------------
C Calculate backward test and find NCALC, the highest N such that
C the test is passed.
C---------------------------------------------------------------------
                        TEST = POLD*PLAST*(HALF-HALF/(TEMPB*TEMPB))
                        TEST = TEST/ENSIG
                        P = PLAST*TOVER
                        N = N - 1
                        EN = EN - TWO
                        NEND = MIN(NB,N)
                        DO 110 L=NSTART,NEND
                           POLD = PSAVEL
                           PSAVEL = PSAVE
                           PSAVE = EN*PSAVEL/X - POLD
                           IF (PSAVE*PSAVEL.GT.TEST) THEN
                              NCALC = L - 1
                              GO TO 190
                           END IF
  110                   CONTINUE
                        NCALC = NEND
                        GO TO 190
                     END IF
  130             CONTINUE
                  N = NEND
                  EN = CONV(N+N) + (ALPHA+ALPHA)
C---------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C---------------------------------------------------------------------
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
               END IF
C---------------------------------------------------------------------
C Calculate P*S until significance test passes.
C---------------------------------------------------------------------
  140          N = N + 1
                  EN = EN + TWO
                  POLD = PLAST
                  PLAST = P
                  P = EN*PLAST/X - POLD
               IF (P.LT.TEST) GO TO 140
C---------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C---------------------------------------------------------------------
  190          N = N + 1
               EN = EN + TWO
               TEMPB = ZERO
               TEMPA = ONE/P
               M = 2*N - 4*(N/2)
               SUM = ZERO
               EM = CONV(N/2)
               ALPEM = (EM-ONE) + ALPHA
               ALP2EM = (EM+EM) + ALPHA
               IF (M .NE. 0) SUM = TEMPA*ALPEM*ALP2EM/EM
               NEND = N - NB
               IF (NEND .GT. 0) THEN
C---------------------------------------------------------------------
C Recur backward via difference equation, calculating (but not
C storing) B(N), until N = NB.
C---------------------------------------------------------------------
                  DO 200 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     TEMPC = TEMPB
                     TEMPB = TEMPA
                     TEMPA = (EN*TEMPB)/X - TEMPC
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        IF (N.EQ.1) GO TO 210
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+TEMPA*ALP2EM)*ALPEM/EM
                     END IF
  200             CONTINUE
               END IF
C---------------------------------------------------------------------
C Store B(NB).
C---------------------------------------------------------------------
  210          B(N) = TEMPA
               IF (NEND .GE. 0) THEN
                  IF (NB .LE. 1) THEN
                        ALP2EM = ALPHA
                        IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE
                        SUM = SUM + B(1)*ALP2EM
                        GO TO 250
                     ELSE
C---------------------------------------------------------------------
C Calculate and store B(NB-1).
C---------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*TEMPA)/X - TEMPB
                        IF (N.EQ.1) GO TO 240
                        M = 2 - M
                        IF (M .NE. 0) THEN
                           EM = EM - ONE
                           ALP2EM = (EM+EM) + ALPHA
                           ALPEM = (EM-ONE) + ALPHA
                           IF (ALPEM.EQ.ZERO) ALPEM = ONE
                           SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                        END IF
                  END IF
               END IF
               NEND = N - 2
               IF (NEND .NE. 0) THEN
C---------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C---------------------------------------------------------------------
                  DO 230 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     B(N) = (EN*B(N+1))/X - B(N+2)
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                     END IF
  230             CONTINUE
               END IF
C---------------------------------------------------------------------
C Calculate B(1).
C---------------------------------------------------------------------
               B(1) = TWO*(ALPHA+ONE)*B(2)/X - B(3)
  240          EM = EM - ONE
               ALP2EM = (EM+EM) + ALPHA
               IF (ALP2EM.EQ.ZERO) ALP2EM = ONE
               SUM = SUM + B(1)*ALP2EM
C---------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C---------------------------------------------------------------------
  250          IF ((ALPHA+ONE).NE.ONE)
     1              SUM = SUM*FUNC(ALPHA)*(X*HALF)**(-ALPHA)
               TEMPA = ENMTEN
               IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
               DO 260 N=1,NB
                 IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO
                 B(N) = B(N)/SUM
  260          CONTINUE
            END IF
C---------------------------------------------------------------------
C Error return -- X, NB, or ALPHA is out of range.
C---------------------------------------------------------------------
         ELSE
            B(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
C---------------------------------------------------------------------
C Exit
C---------------------------------------------------------------------
  300 RETURN
C ---------- Last line of RJBESL ----------
      END
      SUBROUTINE RKBESL(X,ALPHA,NB,IZE,BK,NCALC)
C-------------------------------------------------------------------
C
C  This FORTRAN 77 routine calculates modified Bessel functions
C  of the second kind, K SUB(N+ALPHA) (X), for non-negative
C  argument X, and non-negative order N+ALPHA, with or without
C  exponential scaling.
C
C  Explanation of variables in the calling sequence
C
C  Description of output values ..
C
C X     - Working precision non-negative real argument for which
C         K's or exponentially scaled K's (K*EXP(X))
C         are to be calculated.  If K's are to be calculated,
C         X must not be greater than XMAX (see below).
C ALPHA - Working precision fractional part of order for which
C         K's or exponentially scaled K's (K*EXP(X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
C         and 2 if exponentially scaled K's are to be calculated.
C BK    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BK
C         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
C         or the corresponding exponentially scaled functions.
C         If (0 .LT. NCALC .LT. NB), BK(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BK, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   EPS    = The smallest positive floating-point number such that
C            1.0+EPS .GT. 1.0
C   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution
C            to equation:
C               W(X) * (1-1/8X+9/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C   SQXMIN = Square root of beta**minexp
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMIN   = Smallest positive machine number; approximately
C            beta**minexp
C
C
C     Approximate values for some important machines are:
C
C                          beta       minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)      2        -8193        8191    7.11E-15
C  Cyber 180/185
C    under NOS   (S.P.)      2         -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)      2         -126         128    1.19E-7
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)      2        -1022        1024    2.22D-16
C  IBM 3033      (D.P.)     16          -65          63    2.22D-16
C  VAX           (S.P.)      2         -128         127    5.96E-8
C  VAX D-Format  (D.P.)      2         -128         127    1.39D-17
C  VAX G-Format  (D.P.)      2        -1024        1023    1.11D-16
C
C
C                         SQXMIN       XINF        XMIN      XMAX
C
C CRAY-1        (S.P.)  6.77E-1234  5.45E+2465  4.59E-2467 5674.858
C Cyber 180/855
C   under NOS   (S.P.)  1.77E-147   1.26E+322   3.14E-294   672.788
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.08E-19    3.40E+38    1.18E-38     85.337
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.49D-154   1.79D+308   2.23D-308   705.342
C IBM 3033      (D.P.)  7.35D-40    7.23D+75    5.40D-79    177.852
C VAX           (S.P.)  5.42E-20    1.70E+38    2.94E-39     86.715
C VAX D-Format  (D.P.)  5.42D-20    1.70D+38    2.94D-39     86.715
C VAX G-Format  (D.P.)  7.46D-155   8.98D+307   5.57D-309   706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all K's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, the B-vector is not calculated,
C       and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
C       K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
C       the B-vector is not calculated.  Note that again
C       NCALC .NE. NB.
C
C  0 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BK(I) contains correct function
C       values for I .LE. NCALC, and contains the ratios
C       K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT
C
C
C Acknowledgement
C
C  This program is based on a program written by J. B. Campbell
C  (2) that computes values of the Bessel functions K of real
C  argument and real order.  Modifications include the addition
C  of non-scaled functions, parameterization of machine
C  dependencies, and the use of more accurate approximations
C  for SINH and SIN.
C
C References: "On Temme's Algorithm for the Modified Bessel
C              Functions of the Third Kind," Campbell, J. B.,
C              TOMS 6(4), Dec. 1980, pp. 581-586.
C
C             "A FORTRAN IV Subroutine for the Modified Bessel
C              Functions of the Third Kind of Real Order and Real
C              Argument," Campbell, J. B., Report NRC/ERB-925,
C              National Research Council, Canada.
C
C  Latest modification: March 15, 1992
C
C  Modified by: W. J. Cody and L. Stoltz
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C-------------------------------------------------------------------
      INTEGER I,IEND,ITEMP,IZE,J,K,M,MPLUS1,NB,NCALC
CS    REAL
      DOUBLE PRECISION
     1    A,ALPHA,BLPHA,BK,BK1,BK2,C,D,DM,D1,D2,D3,ENU,EPS,ESTF,ESTM,
     2    EX,FOUR,F0,F1,F2,HALF,ONE,P,P0,Q,Q0,R,RATIO,S,SQXMIN,T,TINYX,
     3    TWO,TWONU,TWOX,T1,T2,WMINF,X,XINF,XMAX,XMIN,X2BY4,ZERO
      DIMENSION BK(NB),P(8),Q(7),R(5),S(4),T(6),ESTM(6),ESTF(7)
C---------------------------------------------------------------------
C  Mathematical constants
C    A = LOG(2.D0) - Euler's constant
C    D = SQRT(2.D0/PI)
C---------------------------------------------------------------------
CS    DATA HALF,ONE,TWO,ZERO/0.5E0,1.0E0,2.0E0,0.0E0/
CS    DATA FOUR,TINYX/4.0E0,1.0E-10/
CS    DATA A/ 0.11593151565841244881E0/,D/0.797884560802865364E0/
      DATA HALF,ONE,TWO,ZERO/0.5D0,1.0D0,2.0D0,0.0D0/
      DATA FOUR,TINYX/4.0D0,1.0D-10/
      DATA A/ 0.11593151565841244881D0/,D/0.797884560802865364D0/
C---------------------------------------------------------------------
C  Machine dependent parameters
C---------------------------------------------------------------------
CS    DATA EPS/1.19E-7/,SQXMIN/1.08E-19/,XINF/3.40E+38/
CS    DATA XMIN/1.18E-38/,XMAX/85.337E0/
      DATA EPS/2.22D-16/,SQXMIN/1.49D-154/,XINF/1.79D+308/
      DATA XMIN/2.23D-308/,XMAX/705.342D0/
C---------------------------------------------------------------------
C  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
C                                         + Euler's constant
C         Coefficients converted from hex to decimal and modified
C         by W. J. Cody, 2/26/82
C  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
C  T    - Approximation for SINH(Y)/Y
C---------------------------------------------------------------------
CS    DATA P/ 0.805629875690432845E00,    0.204045500205365151E02,
CS   1        0.157705605106676174E03,    0.536671116469207504E03,
CS   2        0.900382759291288778E03,    0.730923886650660393E03,
CS   3        0.229299301509425145E03,    0.822467033424113231E00/
CS    DATA Q/ 0.294601986247850434E02,    0.277577868510221208E03,
CS   1        0.120670325591027438E04,    0.276291444159791519E04,
CS   2        0.344374050506564618E04,    0.221063190113378647E04,
CS   3        0.572267338359892221E03/
CS    DATA R/-0.48672575865218401848E+0,  0.13079485869097804016E+2,
CS   1       -0.10196490580880537526E+3,  0.34765409106507813131E+3,
CS   2        0.34958981245219347820E-3/
CS    DATA S/-0.25579105509976461286E+2,  0.21257260432226544008E+3,
CS   1       -0.61069018684944109624E+3,  0.42269668805777760407E+3/
CS    DATA T/ 0.16125990452916363814E-9, 0.25051878502858255354E-7,
CS   1        0.27557319615147964774E-5, 0.19841269840928373686E-3,
CS   2        0.83333333333334751799E-2, 0.16666666666666666446E+0/
CS    DATA ESTM/5.20583E1, 5.7607E0, 2.7782E0, 1.44303E1, 1.853004E2,
CS   1          9.3715E0/
CS    DATA ESTF/4.18341E1, 7.1075E0, 6.4306E0, 4.25110E1, 1.35633E0,
CS   1          8.45096E1, 2.0E1/
      DATA P/ 0.805629875690432845D00,    0.204045500205365151D02,
     1        0.157705605106676174D03,    0.536671116469207504D03,
     2        0.900382759291288778D03,    0.730923886650660393D03,
     3        0.229299301509425145D03,    0.822467033424113231D00/
      DATA Q/ 0.294601986247850434D02,    0.277577868510221208D03,
     1        0.120670325591027438D04,    0.276291444159791519D04,
     2        0.344374050506564618D04,    0.221063190113378647D04,
     3        0.572267338359892221D03/
      DATA R/-0.48672575865218401848D+0,  0.13079485869097804016D+2,
     1       -0.10196490580880537526D+3,  0.34765409106507813131D+3,
     2        0.34958981245219347820D-3/
      DATA S/-0.25579105509976461286D+2,  0.21257260432226544008D+3,
     1       -0.61069018684944109624D+3,  0.42269668805777760407D+3/
      DATA T/ 0.16125990452916363814D-9, 0.25051878502858255354D-7,
     1        0.27557319615147964774D-5, 0.19841269840928373686D-3,
     2        0.83333333333334751799D-2, 0.16666666666666666446D+0/
      DATA ESTM/5.20583D1, 5.7607D0, 2.7782D0, 1.44303D1, 1.853004D2,
     1          9.3715D0/
      DATA ESTF/4.18341D1, 7.1075D0, 6.4306D0, 4.25110D1, 1.35633D0,
     1          8.45096D1, 2.0D1/
C---------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      NCALC = MIN(NB,0)-2
      IF ((NB .GT. 0) .AND. ((ENU .GE. ZERO) .AND. (ENU .LT. ONE))
     1     .AND. ((IZE .GE. 1) .AND. (IZE .LE. 2)) .AND.
     2     ((IZE .NE. 1) .OR. (EX .LE. XMAX)) .AND.
     3     (EX .GT. ZERO))  THEN
            K = 0
            IF (ENU .LT. SQXMIN) ENU = ZERO
            IF (ENU .GT. HALF) THEN
                  K = 1
                  ENU = ENU - ONE
            END IF
            TWONU = ENU+ENU
            IEND = NB+K-1
            C = ENU*ENU
            D3 = -C
            IF (EX .LE. ONE) THEN
C---------------------------------------------------------------------
C  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
C                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
C---------------------------------------------------------------------
                  D1 = ZERO
                  D2 = P(1)
                  T1 = ONE
                  T2 = Q(1)
                  DO 10 I = 2,7,2
                     D1 = C*D1+P(I)
                     D2 = C*D2+P(I+1)
                     T1 = C*T1+Q(I)
                     T2 = C*T2+Q(I+1)
   10             CONTINUE
                  D1 = ENU*D1
                  T1 = ENU*T1
                  F1 = LOG(EX)
                  F0 = A+ENU*(P(8)-ENU*(D1+D2)/(T1+T2))-F1
                  Q0 = EXP(-ENU*(A-ENU*(P(8)+ENU*(D1-D2)/(T1-T2))-F1))
                  F1 = ENU*F0
                  P0 = EXP(F1)
C---------------------------------------------------------------------
C  Calculation of F0 =
C---------------------------------------------------------------------
                  D1 = R(5)
                  T1 = ONE
                  DO 20 I = 1,4
                     D1 = C*D1+R(I)
                     T1 = C*T1+S(I)
   20             CONTINUE
                  IF (ABS(F1) .LE. HALF) THEN
                        F1 = F1*F1
                        D2 = ZERO
                        DO 30 I = 1,6
                           D2 = F1*D2+T(I)
   30                   CONTINUE
                        D2 = F0+F0*F1*D2
                     ELSE
                        D2 = SINH(F1)/ENU
                  END IF
                  F0 = D2-ENU*D1/(T1*P0)
                  IF (EX .LE. TINYX) THEN
C--------------------------------------------------------------------
C  X.LE.1.0E-10
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                        BK(1) = F0+EX*F0
                        IF (IZE .EQ. 1) BK(1) = BK(1)-EX*BK(1)
                        RATIO = P0/F0
                        C = EX*XINF
                        IF (K .NE. 0) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
C  ALPHA .GE. 1/2
C--------------------------------------------------------------------
                              NCALC = -1
                              IF (BK(1) .GE. C/RATIO) GO TO 500
                              BK(1) = RATIO*BK(1)/EX
                              TWONU = TWONU+TWO
                              RATIO = TWONU
                        END IF
                        NCALC = 1
                        IF (NB .EQ. 1) GO TO 500
C--------------------------------------------------------------------
C  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
C--------------------------------------------------------------------
                        NCALC = -1
                        DO 80 I = 2,NB
                           IF (RATIO .GE. C) GO TO 500
                           BK(I) = RATIO/EX
                           TWONU = TWONU+TWO
                           RATIO = TWONU
   80                   CONTINUE
                        NCALC = 1
                        GO TO 420
                     ELSE
C--------------------------------------------------------------------
C  1.0E-10 .LT. X .LE. 1.0
C--------------------------------------------------------------------
                        C = ONE
                        X2BY4 = EX*EX/FOUR
                        P0 = HALF*P0
                        Q0 = HALF*Q0
                        D1 = -ONE
                        D2 = ZERO
                        BK1 = ZERO
                        BK2 = ZERO
                        F1 = F0
                        F2 = P0
  100                   D1 = D1+TWO
                        D2 = D2+ONE
                        D3 = D1+D3
                        C = X2BY4*C/D2
                        F0 = (D2*F0+P0+Q0)/D3
                        P0 = P0/(D2-ENU)
                        Q0 = Q0/(D2+ENU)
                        T1 = C*F0
                        T2 = C*(P0-D2*F0)
                        BK1 = BK1+T1
                        BK2 = BK2+T2
                        IF ((ABS(T1/(F1+BK1)) .GT. EPS) .OR.
     1                     (ABS(T2/(F2+BK2)) .GT. EPS))  GO TO 100
                        BK1 = F1+BK1
                        BK2 = TWO*(F2+BK2)/EX
                        IF (IZE .EQ. 2) THEN
                              D1 = EXP(EX)
                              BK1 = BK1*D1
                              BK2 = BK2*D1
                        END IF
                        WMINF = ESTF(1)*EX+ESTF(2)
                  END IF
               ELSE IF (EPS*EX .GT. ONE) THEN
C--------------------------------------------------------------------
C  X .GT. ONE/EPS
C--------------------------------------------------------------------
                  NCALC = NB
                  BK1 = ONE / (D*SQRT(EX))
                  DO 110 I = 1, NB
                     BK(I) = BK1
  110             CONTINUE
                  GO TO 500
               ELSE
C--------------------------------------------------------------------
C  X .GT. 1.0
C--------------------------------------------------------------------
                  TWOX = EX+EX
                  BLPHA = ZERO
                  RATIO = ZERO
                  IF (EX .LE. FOUR) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 .LE. X .LE. 4.0
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(1)/EX+ESTM(2))
                        M = INT(D2)
                        D1 = D2+D2
                        D2 = D2-HALF
                        D2 = D2*D2
                        DO 120 I = 2,M
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
  120                   CONTINUE
C--------------------------------------------------------------------
C  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
C    recurrence and K(ALPHA,X) from the wronskian
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(3)*EX+ESTM(4))
                        M = INT(D2)
                        C = ABS(ENU)
                        D3 = C+C
                        D1 = D3-ONE
                        F1 = XMIN
                        F0 = (TWO*(C+D2)/EX+HALF*EX/(C+D2+ONE))*XMIN
                        DO 130 I = 3,M
                           D2 = D2-ONE
                           F2 = (D3+D2+D2)*F0
                           BLPHA = (ONE+D1/D2)*(F2+BLPHA)
                           F2 = F2/EX+F1
                           F1 = F0
                           F0 = F2
  130                   CONTINUE
                        F1 = (D3+TWO)*F0/EX+F1
                        D1 = ZERO
                        T1 = ONE
                        DO 140 I = 1,7
                           D1 = C*D1+P(I)
                           T1 = C*T1+Q(I)
  140                   CONTINUE
                        P0 = EXP(C*(A+C*(P(8)-C*D1/T1)-LOG(EX)))/EX
                        F2 = (C+HALF-RATIO)*F1/EX
                        BK1 = P0+(D3*F0-F2+F0+BLPHA)/(F2+F1+F0)*P0
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(3)*EX+ESTF(4)
                     ELSE
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
C  recurrence, for  X .GT. 4.0
C--------------------------------------------------------------------
                        DM = AINT(ESTM(5)/EX+ESTM(6))
                        M = INT(DM)
                        D2 = DM-HALF
                        D2 = D2*D2
                        D1 = DM+DM
                        DO 160 I = 2,M
                           DM = DM-ONE
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
                           BLPHA = (RATIO+RATIO*BLPHA)/DM
  160                   CONTINUE
                        BK1 = ONE/((D+D*BLPHA)*SQRT(EX))
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(5)*(EX-ABS(EX-ESTF(7)))+ESTF(6)
                  END IF
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
C    K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                  BK2 = BK1+BK1*(ENU+HALF-RATIO)/EX
            END IF
C--------------------------------------------------------------------
C  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
C  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
C--------------------------------------------------------------------
            NCALC = NB
            BK(1) = BK1
            IF (IEND .EQ. 0) GO TO 500
            J = 2-K
            IF (J .GT. 0) BK(J) = BK2
            IF (IEND .EQ. 1) GO TO 500
            M = MIN(INT(WMINF-ENU),IEND)
            DO 190 I = 2,M
               T1 = BK1
               BK1 = BK2
               TWONU = TWONU+TWO
               IF (EX .LT. ONE) THEN
                     IF (BK1 .GE. (XINF/TWONU)*EX) GO TO 195
                     GO TO 187
                  ELSE
                     IF (BK1/EX .GE. XINF/TWONU) GO TO 195
               END IF
  187          CONTINUE
               BK2 = TWONU/EX*BK1+T1
               ITEMP = I
               J = J+1
               IF (J .GT. 0) BK(J) = BK2
  190       CONTINUE
  195       M = ITEMP
            IF (M .EQ. IEND) GO TO 500
            RATIO = BK2/BK1
            MPLUS1 = M+1
            NCALC = -1
            DO 410 I = MPLUS1,IEND
               TWONU = TWONU+TWO
               RATIO = TWONU/EX+ONE/RATIO
               J = J+1
               IF (J .GT. 1) THEN
                     BK(J) = RATIO
                  ELSE
                     IF (BK2 .GE. XINF/RATIO) GO TO 500
                     BK2 = RATIO*BK2
               END IF
  410       CONTINUE
            NCALC = MAX(MPLUS1-K,1)
            IF (NCALC .EQ. 1) BK(1) = BK2
            IF (NB .EQ. 1) GO TO 500
  420       J = NCALC+1
            DO 430 I = J,NB
               IF (BK(NCALC) .GE. XINF/BK(I)) GO TO 500
               BK(I) = BK(NCALC)*BK(I)
               NCALC = I
  430       CONTINUE
      END IF
  500 RETURN
C---------- Last line of RKBESL ----------
      END
      SUBROUTINE RYBESL(X,ALPHA,NB,BY,NCALC)
C----------------------------------------------------------------------
C
C  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision positive real argument for which
C         Y's are to be calculated.
C ALPHA - Working precision fractional part of order for which
C         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the
C         last is of order (NB - 1 + ALPHA).
C BY    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BY
C         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
C         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BY, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.  Let
C
C   beta   = Radix for the floating-point system
C   p      = Number of significant base-beta digits in the
C            significand of a floating-point number
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C
C Then the following machine-dependent constants must be declared
C   in DATA statements.  IEEE values are provided as a default.
C
C   EPS    = beta ** (-p)
C   DEL    = Machine number below which sin(x)/x = 1; approximately
C            SQRT(EPS).
C   XMIN   = Smallest acceptable argument for RBESY; approximately
C            max(2*beta**minexp,2/XINF), rounded up
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   THRESH = Lower bound for use of the asymptotic form; approximately
C            AINT(-LOG10(EPS/2.0))+1.0
C   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
C            and cosine functions have lost about half of their
C            precision at that point.
C
C
C     Approximate values for some important machines are:
C
C                        beta    p     minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15
C  Cyber 180/185
C    under NOS   (S.P.)    2    48      -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16
C  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17
C  VAX           (S.P.)    2    24      -128         127    5.96E-8
C  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17
C  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16
C
C
C                         DEL      XMIN      XINF     THRESH  XLARGE
C
C CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7
C Cyber 180/855
C   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8
C IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8
C VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4
C VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9
C VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all Y's are
C  calculated to the desired accuracy.
C
C  NCALC .LE. -1:  An argument is out of range. For example,
C       NB .LE. 0, or ABS(X) .GE. XLARGE.  In this case,
C       BY(1) = 0.0, the remainder of the BY-vector is not
C       calculated, and NCALC is set to MIN0(NB,0)-1  so that
C       NCALC .NE. NB.
C  1 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BY(I) contains correct function
C       values for I .LE. NCALC, and and the remaining NB-NCALC
C       array elements contain 0.0.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, COS, INT, LOG, SIN, SQRT
C
C
C Acknowledgement
C
C  This program draws heavily on Temme's Algol program for Y(a,x)
C  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
C  scheme is used for  x < THRESH, and Campbell's scheme is used
C  in the asymptotic region.  Segments of code from both sources
C  have been translated into Fortran 77, merged, and heavily modified.
C  Modifications include parameterization of machine dependencies,
C  use of a new approximation for ln(gamma(x)), and built-in
C  protection against overflow and destructive underflow.
C
C References: "Bessel functions J_nu(x) and Y_nu(x) of real
C              order and real argument," Campbell, J. B.,
C              Comp. Phy. Comm. 18, 1979, pp. 133-142.
C
C             "On the numerical evaluation of the ordinary
C              Bessel function of the second kind," Temme,
C              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.
C
C  Latest modification: March 15, 1992
C
C  Modified by: W. J. Cody
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C----------------------------------------------------------------------
      INTEGER I,K,NA,NB,NCALC
CS    REAL
      DOUBLE PRECISION
     1  ALFA,ALPHA,AYE,B,BY,C,CH,COSMU,D,DEL,DEN,DDIV,DIV,DMU,D1,D2,
     2  E,EIGHT,EN,ENU,EN1,EPS,EVEN,EX,F,FIVPI,G,GAMMA,H,HALF,ODD,
     3  ONBPI,ONE,ONE5,P,PA,PA1,PI,PIBY2,PIM5,Q,QA,QA1,Q0,R,S,SINMU,
     4  SQ2BPI,TEN9,TERM,THREE,THRESH,TWO,TWOBYX,X,XINF,XLARGE,XMIN,
     5  XNA,X2,YA,YA1,ZERO
      DIMENSION BY(NB),CH(21)
C----------------------------------------------------------------------
C  Mathematical constants
C    FIVPI = 5*PI
C    PIM5 = 5*PI - 15
C    ONBPI = 1/PI
C    PIBY2 = PI/2
C    SQ2BPI = SQUARE ROOT OF 2/PI
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/
CS    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/
CS    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/
CS    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/
CS    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA EIGHT,ONE5,TEN9/8.0D0,15.0D0,1.9D1/
      DATA FIVPI,PIBY2/1.5707963267948966192D1,1.5707963267948966192D0/
      DATA PI,SQ2BPI/3.1415926535897932385D0,7.9788456080286535588D-1/
      DATA PIM5,ONBPI/7.0796326794896619231D-1,3.1830988618379067154D-1/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/
CS    DATA THRESH,XLARGE/8.0E0,1.0E4/
      DATA DEL,XMIN,XINF,EPS/1.0D-8,4.46D-308,1.79D308,1.11D-16/
      DATA THRESH,XLARGE/16.0D0,1.0D8/
C----------------------------------------------------------------------
C  Coefficients for Chebyshev polynomial expansion of
C         1/gamma(1-x), abs(x) .le. .5
C----------------------------------------------------------------------
CS    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22,
CS   1         0.29017595056104745456E-20, 0.13639417919073099464E-18,
CS   2         0.23826220476859635824E-17,-0.90642907957550702534E-17,
CS   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13,
CS   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11,
CS   5         0.24230957900482704055E-09, 0.17451364971382984243E-08,
CS   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06,
CS   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04,
CS   8         0.12719271366545622927E-02, 0.17063050710955562222E-02,
CS   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00,
CS   A         0.92187029365045265648E+00/
      DATA CH/-0.67735241822398840964D-23,-0.61455180116049879894D-22,
     1         0.29017595056104745456D-20, 0.13639417919073099464D-18,
     2         0.23826220476859635824D-17,-0.90642907957550702534D-17,
     3        -0.14943667065169001769D-14,-0.33919078305362211264D-13,
     4        -0.17023776642512729175D-12, 0.91609750938768647911D-11,
     5         0.24230957900482704055D-09, 0.17451364971382984243D-08,
     6        -0.33126119768180852711D-07,-0.86592079961391259661D-06,
     7        -0.49717367041957398581D-05, 0.76309597585908126618D-04,
     8         0.12719271366545622927D-02, 0.17063050710955562222D-02,
     9        -0.76852840844786673690D-01,-0.28387654227602353814D+00,
     A         0.92187029365045265648D+00/
C----------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      IF ((NB .GT. 0) .AND. (X .GE. XMIN) .AND. (EX .LT. XLARGE)
     1       .AND. (ENU .GE. ZERO) .AND. (ENU .LT. ONE))  THEN
            XNA = AINT(ENU+HALF)
            NA = INT(XNA)
            IF (NA .EQ. 1) ENU = ENU - XNA
            IF (ENU .EQ. -HALF) THEN
                  P = SQ2BPI/SQRT(EX)
                  YA = P * SIN(EX)
                  YA1 = -P * COS(EX)
               ELSE IF (EX .LT. THREE) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for small X
C----------------------------------------------------------------------
                  B = EX * HALF
                  D = -LOG(B)
                  F = ENU * D
                  E = B**(-ENU)
                  IF (ABS(ENU) .LT. DEL) THEN
                        C = ONBPI
                     ELSE
                        C = ENU / SIN(ENU*PI)
                  END IF
C----------------------------------------------------------------------
C  Computation of sinh(f)/f
C----------------------------------------------------------------------
                  IF (ABS(F) .LT. ONE) THEN
                        X2 = F*F
                        EN = TEN9
                        S = ONE
                        DO 80 I = 1, 9
                           S = S*X2/EN/(EN-ONE)+ONE
                           EN = EN - TWO
   80                   CONTINUE
                     ELSE
                        S = (E - ONE/E) * HALF / F
                  END IF
C----------------------------------------------------------------------
C  Computation of 1/gamma(1-a) using Chebyshev polynomials
C----------------------------------------------------------------------
                  X2 = ENU*ENU*EIGHT
                  AYE = CH(1)
                  EVEN = ZERO
                  ALFA = CH(2)
                  ODD = ZERO
                  DO 40 I = 3, 19, 2
                     EVEN = -(AYE+AYE+EVEN)
                     AYE = -EVEN*X2 - AYE + CH(I)
                     ODD = -(ALFA+ALFA+ODD)
                     ALFA = -ODD*X2 - ALFA + CH(I+1)
   40             CONTINUE
                  EVEN = (EVEN*HALF+AYE)*X2 - AYE + CH(21)
                  ODD = (ODD+ALFA)*TWO
                  GAMMA = ODD*ENU + EVEN
C----------------------------------------------------------------------
C  End of computation of 1/gamma(1-a)
C----------------------------------------------------------------------
                  G = E * GAMMA
                  E = (E + ONE/E) * HALF
                  F = TWO*C*(ODD*E+EVEN*S*D)
                  E = ENU*ENU
                  P = G*C
                  Q = ONBPI / G
                  C = ENU*PIBY2
                  IF (ABS(C) .LT. DEL) THEN
                        R = ONE
                     ELSE
                        R = SIN(C)/C
                  END IF
                  R = PI*C*R*R
                  C = ONE
                  D = - B*B
                  H = ZERO
                  YA = F + R*Q
                  YA1 = P
                  EN = ZERO
  100             EN = EN + ONE
                  IF (ABS(G/(ONE+ABS(YA)))
     1                      + ABS(H/(ONE+ABS(YA1))) .GT. EPS) THEN
                        F = (F*EN+P+Q)/(EN*EN-E)
                        C = C * D/EN
                        P = P/(EN-ENU)
                        Q = Q/(EN+ENU)
                        G = C*(F+R*Q)
                        H = C*P - EN*G
                        YA = YA + G
                        YA1 = YA1+H
                        GO TO 100
                  END IF
                  YA = -YA
                  YA1 = -YA1/B
               ELSE IF (EX .LT. THRESH) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for moderate X
C----------------------------------------------------------------------
                  C = (HALF-ENU)*(HALF+ENU)
                  B = EX + EX
                  E = (EX*ONBPI*COS(ENU*PI)/EPS)
                  E = E*E
                  P = ONE
                  Q = -EX
                  R = ONE + EX*EX
                  S = R
                  EN = TWO
  200             IF (R*EN*EN .LT. E) THEN
                        EN1 = EN+ONE
                        D = (EN-ONE+C/EN)/S
                        P = (EN+EN-P*D)/EN1
                        Q = (-B+Q*D)/EN1
                        S = P*P + Q*Q
                        R = R*S
                        EN = EN1
                        GO TO 200
                  END IF
                  F = P/S
                  P = F
                  G = -Q/S
                  Q = G
  220             EN = EN - ONE
                  IF (EN .GT. ZERO) THEN
                        R = EN1*(TWO-P)-TWO
                        S = B + EN1*Q
                        D = (EN-ONE+C/EN)/(R*R+S*S)
                        P = D*R
                        Q = D*S
                        E = F + ONE
                        F = P*E - G*Q
                        G = Q*E + P*G
                        EN1 = EN
                        GO TO 220
                  END IF
                  F = ONE + F
                  D = F*F + G*G
                  PA = F/D
                  QA = -G/D
                  D = ENU + HALF -P
                  Q = Q + EX
                  PA1 = (PA*Q-QA*D)/EX
                  QA1 = (QA*Q+PA*D)/EX
                  B = EX - PIBY2*(ENU+HALF)
                  C = COS(B)
                  S = SIN(B)
                  D = SQ2BPI/SQRT(EX)
                  YA = D*(PA*S+QA*C)
                  YA1 = D*(QA1*S-PA1*C)
               ELSE
C----------------------------------------------------------------------
C  Use Campbell's asymptotic scheme.
C----------------------------------------------------------------------
                  NA = 0
                  D1 = AINT(EX/FIVPI)
                  I = INT(D1)
                  DMU = ((EX-ONE5*D1)-D1*PIM5)-(ALPHA+HALF)*PIBY2
                  IF (I-2*(I/2) .EQ. 0) THEN
                        COSMU = COS(DMU)
                        SINMU = SIN(DMU)
                     ELSE
                        COSMU = -COS(DMU)
                        SINMU = -SIN(DMU)
                  END IF
                  DDIV = EIGHT * EX
                  DMU = ALPHA
                  DEN = SQRT(EX)
                  DO 350 K = 1, 2
                     P = COSMU
                     COSMU = SINMU
                     SINMU = -P
                     D1 = (TWO*DMU-ONE)*(TWO*DMU+ONE)
                     D2 = ZERO
                     DIV = DDIV
                     P = ZERO
                     Q = ZERO
                     Q0 = D1/DIV
                     TERM = Q0
                     DO 310 I = 2, 20
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = -TERM*D1/DIV
                        P = P + TERM
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = TERM*D1/DIV
                        Q = Q + TERM
                        IF (ABS(TERM) .LE. EPS) GO TO 320
  310                CONTINUE
  320                P = P + ONE
                     Q = Q + Q0
                     IF (K .EQ. 1) THEN
                           YA = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                        ELSE
                           YA1 = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                     END IF
                     DMU = DMU + ONE
  350             CONTINUE
            END IF
            IF (NA .EQ. 1) THEN
               H = TWO*(ENU+ONE)/EX
               IF (H .GT. ONE) THEN
                  IF (ABS(YA1) .GT. XINF/H) THEN
                     H = ZERO
                     YA = ZERO
                  END IF
               END IF
               H = H*YA1 - YA
               YA = YA1
               YA1 = H
            END IF
C----------------------------------------------------------------------
C  Now have first one or two Y's
C----------------------------------------------------------------------
            BY(1) = YA
            BY(2) = YA1
            IF (YA1 .EQ. ZERO) THEN
                  NCALC = 1
               ELSE
                  AYE = ONE + ALPHA
                  TWOBYX = TWO/EX
                  NCALC = 2
                  DO 400 I = 3, NB
                     IF (TWOBYX .LT. ONE) THEN
                           IF (ABS(BY(I-1))*TWOBYX .GE. XINF/AYE)
     1                                                     GO TO 450
                        ELSE
                           IF (ABS(BY(I-1)) .GE. XINF/AYE/TWOBYX )
     1                                                     GO TO 450
                     END IF
                     BY(I) = TWOBYX*AYE*BY(I-1) - BY(I-2)
                     AYE = AYE + ONE
                     NCALC = NCALC + 1
  400             CONTINUE
            END IF
  450       DO 460 I = NCALC+1, NB
               BY(I) = ZERO
  460       CONTINUE
         ELSE
            BY(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
      RETURN
C---------- Last line of RYBESL ----------
      END
