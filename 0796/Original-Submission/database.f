C ***********************************************
C                                              **
C FUNCTION FZ CALLED BY THE INVERSION SOFTWARE **
C                                              **
C ***********************************************
      COMPLEX*16 FUNCTION FZ(Z)
C ***********************************************
C ************************************************
C     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
C     ..
C     .. Scalars in Common ..
      INTEGER NFUN
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX B,CI,CVAR,ARG,Z2,Z4,A1,A2,BB
      DOUBLE PRECISION R,C
C     ..
C     .. External Functions ..
      INTRINSIC CDEXP,CDLOG,CDSQRT
C     ..
C     .. Common blocks ..
      COMMON /NF/NFUN
C     ..
      GOTO(1,2,3,4,5,6,7,8,10,11,12,13,14,15,
     *  16,17,18,19,20,
     *  21,22,23,24,25,26,27,28,29,30,31,32,
     *  33,34,35,36,37,38) NFUN
C ************************************************
C FUNCTIONS ORDER AND NUMBERING REFLECT FUNCTIONS 
C ORDER AND NUMBERING AS IN THE PAPER:
C 
C   D'AMORE L., LACCETTI G., MURLI A.,  - 
C
C  "ALGORITHM XXX: A FORTRAN SOFTWARE
C   PACKAGE FOR THE NUMERICAL INVERSION OF THE 
C   LAPLACE TRANSFORM BASE ON FOURIER SERIES' METHOD"
C
C   ACM TRANS. MATH. SOFTWARE, VOL. ##,
C   NO. #, MONTH YEAR, PP. ##-##.
C *******************************************************
  1   FZ = (1.D0,0.D0)/Z
      GO TO 350
C
C ************************************************
  2   FZ = (2.D0,0.D0)* (CDSQRT(Z+ (1.D0,0.D0))-CDSQRT(Z))
      GO TO 350
C ************************************************
  3   FZ = (1.D0,0.D0)/CDSQRT(Z)
      GO TO 350
C ************************************************
C
  4   FZ = (Z*Z- (1.D0,0.D0))/ ((Z*Z+ (1.D0,0.D0))**2)
      GO TO 350
C ************************************************
C
  5   FZ = (1.D0,0.D0)/ (Z+ (1.D0,0.D0))**2
      GO TO 350
C ************************************************
C
  6   FZ = (1.D0,0.D0)/ (Z**2)
      GO TO 350
C ************************************************
C
  7   FZ = (1.D0,0.D0)/ (Z**2+ (1.D0,0.D0))
      GO TO 350
C ************************************************
C
  8   FZ = (1.D0,0.D0)/ (Z+ (0.5D0,0.D0))
      GO TO 350
C ************************************************
C
   9  FZ = 1.D0/CDSQRT(Z*Z+ (1.D0,0.D0))
      GO TO 350
C ************************************************
C
   10 FZ = CDEXP((-1.D0,0.D0)/Z)/CDSQRT(Z)
      GO TO 350
C ************************************************
C
  11  FZ = CDEXP((-4.D0,0.D0)*CDSQRT(Z))
      GO TO 350
C*************************************************
C
  12  CI = (0.D0,1.D0)
      CVAR = 1.D0/ (2* (0.D0,1.D0))
      FZ = CVAR*CDLOG((Z+CI)/ (Z-CI))
      GO TO 350
C *************************************************
C
  13  FZ = (1.D0,0.D0)/ ((Z+.2)**2+ (1.D0,0.D0))
      GO TO 350
C ************************************************
C
  14  FZ =  1.d0/Z**3
      GO TO 350
C ************************************************
C
   15 FZ = CDEXP(-2*Z)/Z
      GO TO 350
C *************************************************
C
  16  FZ = (1.D0,0.D0)/ (Z* ((1.D0,0.D0)+CDEXP(-Z)))
      GO TO 350
C ************************************************
   17 FZ = (1.D0,0.D0)/ (Z*Z+Z+ (1.D0,0.D0))
      GO TO 350
C *************************************************
C
  18  FZ = (3.D0,0.D0)/ (Z**2- (9.D0,0.D0))
      GO TO 350
C
C *************************************************
C
  19  FZ = (120.D0,0.D0)/Z**6
      GO TO 350
C
C ************************************************
C
   20 FZ = Z/ (Z**2+ (1.D0,0.D0))**2
      GO TO 350
C
C ************************************************
C
  21  FZ = (1.D0,0.D0)/ (Z+ (1.D0,0.D0)) -
     +     (1.D0,0.D0)/ (Z+ (1000.D0,0.D0))
      GO TO 350
C
C ************************************************
C
  22  FZ = Z/ (Z*Z+ (1.D0,0.D0))
      GO TO 350
C ************************************************
C
  23  FZ = (1.D0,0.D0)/ ((Z- (0.25D0,0.D0))**2)
      GO TO 350
C ************************************************
C
  24  FZ = (1.D0,0.D0)/ (Z*CDSQRT(Z))
      GO TO 350
C ************************************************
C
  25  FZ = (1.D0,0.D0)/CDSQRT(Z+ (1.D0,0.D0))
      GO TO 350
C ************************************************
C
  26  FZ = (Z+ (2.D0,0.D0))/ (Z*CDSQRT(Z))
      GO TO 350
C ************************************************
C
  27  FZ = (1.D0,0.D0)/ ((Z*Z+ (1.D0,0.D0))**2)
      GO TO 350
C ************************************************
C
  28  FZ = (1.D0,0.D0)/ (Z* (Z+ (1.D0,0.D0))**2)
      GO TO 350
C ************************************************
C
  29  FZ = (1.D0,0.D0)/ (Z**3- (8.D0,0.D0))
      GO TO 350
C ************************************************
C
  30  FZ = CDLOG((Z*Z+ (1.D0,0.D0))/ (Z*Z+ (4.D0,0.D0)))
      GO TO 350
C
C ************************************************
C
  31  FZ = CDLOG((Z+ (1.D0,0.D0))/Z)
      GO TO 350
C ************************************************
  32  FZ= CDLOG(Z)/Z
      GO TO 350
C ************************************************
C
  33  FZ = ((1.D0,0.D0)-CDEXP(-Z))/ (Z*Z)
      GO TO 350
C ************************************************
C
  34  FZ = (1.D0,0.D0)/ (Z* ((1.D0,0.D0)+CDEXP(Z)))
      GO TO 350
C
C *************************************************
C
  35  B = (1.D0,0.D0)/ (2.*Z) - (CDEXP(-2.*Z)/ (1.-CDEXP(-2.*Z)))
      FZ = (1./ (Z*Z+Z))*B
      GO TO 350
C
C **************************************************
  36   R=0.5D0
       C=0.4D0
       B=-R*CDSQRT((Z*(1+Z))/(1+C*Z))
       FZ=1./Z*CDEXP(B)
       GO TO 350
C **************************************************
  37    CONTINUE
C
C         FZ=CDEXP(-2.*PSI)/Z
C
C   COSH(PSI)= SQRT(1+Z**2+((Z**4)/16))
C
         Z2=Z*Z
         Z4=Z2*Z2
         ARG=1.d0+Z2+Z4/16
         A1=CDSQRT(ARG)
         A2=(Z/4.d0)*CDSQRT(16.d0+Z2)
         FZ=1.d0/(Z*(A1+A2)**2)
         GO TO 350
C **************************************************
  38  continue
      B=Z-CDSQRT(Z*Z-1.d0)
      BB=CDSQRT(Z)*CDSQRT(Z*Z-1.D0)*CDSQRT(Z-0.5D0*(CDSQRT(Z*Z-1.D0)))
      FZ=B/BB
      GOTO 350
C **************************************************

350   RETURN

      END

c *********************************************************************
c
C
C    FZ's COMPANION FUNCTION FEX TO COMPUTE THE EXACT VALUE OF THE   **
C    INVERSE TRANSFORM                                               **
C
c
c *********************************************************************
      DOUBLE PRECISION FUNCTION FEX(X)
C
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Scalars in Common ..
      INTEGER NFUN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,PI,PI2,SUM,T
      DOUBLE PRECISION U1,U2,PARTE1,DUFFY6,COMX,EULERO
      INTEGER K,N
      LOGICAL TROV
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DCOS,DEXP,DSIN,DSINH,DSQRT,DFLOAT
C     ..
C     ..
      INTEGER          LW, LIW
      PARAMETER        (LW=1500,LIW=LW/4)
      INTEGER          NOUT
      PARAMETER        (NOUT=6)
      INTEGER          KOUNT
      DOUBLE PRECISION ABSERR, EPSABS, EPSREL, RESULT
      INTEGER          IFAIL, INF
      DOUBLE PRECISION W(LW)
      INTEGER          IW(LIW)
      DOUBLE PRECISION FST3,FST4,FST52,FST51
      EXTERNAL         FST3,FST4,FST52,FST51
C     .. Common blocks ..
      COMMON /NF/NFUN
      COMMON           /TELNUM/COMX,KOUNT

      PI = 4.D0*DATAN(1.D0)
      COMX=X
      EULERO = .5772156D0

      GOTO(1,2,3,4,5,6,7,8,10,11,12,13,14,
     +     15,16,17,18,19,20,21,22,
     +     23,24,25,26,27,28,29,30,31,32,33,34,
c    +     35,36,37,38)NFUN
     +     35)NFUN
C ********************************************************************
C FUNCTIONS ORDER AND NUMBERING REFLECT FUNCTIONS 
C ORDER AND NUMBERING AS IN THE PAPER:
C 
C   D'AMORE L., LACCETTI G., MURLI A.,  - 
C
C  "ALGORITHM XXX: A FORTRAN SOFTWARE
C   PACKAGE FOR THE NUMERICAL INVERSION OF THE 
C   LAPLACE TRANSFORM BASE ON FOURIER SERIES' METHOD"
C
C   ACM TRANS. MATH. SOFTWARE, VOL. ##,
C   NO. #, MONTH YEAR, PP. ##-##.
C ********************************************************************
C
   1  FEX = 1.D0
      GO TO 380
C ********************************************************************
C
   2  FEX = (1.D0-DEXP(-X))/ (X*DSQRT(PI*X))
      GO TO 380
C ********************************************************************
C
  3   FEX = 1.D0/DSQRT(PI*X)
      GO TO 380
C ********************************************************************
C
  4   FEX = X*DCOS(X)
      GO TO 380
C ********************************************************************
C
  5   FEX = X*DEXP(-X)
      GO TO 380
C ********************************************************************
C
  6   FEX = X
      GO TO 380
C
C ********************************************************************
C
  7   FEX = DSIN(X)
      GO TO 380
C ****************************************************
C
  8   FEX = DEXP(-.5D0*X)
      GO TO 380
C ********************************************************************
C       Such inverse function is the Bessel function J_0.
C       Actually we compute it by using the Nag library. That's why
C       it does not appear here .
C
c 9   FEX = S17AEF(X,IFAIL)
c     GO TO 380
C ********************************************************************
C
   10 FEX = DCOS(2.D0*DSQRT(X))/DSQRT(PI*X)
      GO TO 380
C
C ********************************************************************
C
   11 FEX = 2.D0*DEXP(-4.D0/X)/ (X*DSQRT(PI*X))
      GO TO 380
C ****************************************************************
C
  12  FEX = DSIN(X)/X
      GO TO 380
C ****************************************************
C
  13  FEX = DEXP(-.2D0*X)*DSIN(X)
      GO TO 380
C ********************************************************************
C
  14  FEX = 0.5d0 * X**2
      GO TO 380
C ********************************************************************
C
   15 IF (X.GT.2.D0) THEN
          FEX = 1.D0

      ELSE IF (X.LT.2.D0) THEN
          FEX = 0.D0

      ELSE
          FEX = 0.5D0
      END IF

      GO TO 380
C ****************************************************************
C
  16  TROV = .FALSE.
      K = 0
  310 CONTINUE
      IF (X.GT.2*K .AND. X.LT.2*K+1) THEN
          FEX = 1.D0
          TROV = .TRUE.

      ELSE IF (X.GT.2*K+1 .AND. X.LT.2*K+2) THEN
          FEX = 0.D0
          TROV = .TRUE.
      END IF

      IF (X.EQ.DFLOAT(K)) THEN
          FEX = 0.5D0
          TROV = .TRUE.
      END IF

      K = K + 1
      IF (.NOT.TROV .AND. K.LE.49) GO TO 310
      GO TO 380
C
C ********************************************************************
C
   17 FEX = 2.D0/DSQRT(3.D0)*DEXP(-X/2.D0)*DSIN(X*DSQRT(3.D0)/2.D0)
      GO TO 380
C
C *****************************************************
C
  18  FEX = DSINH(3.D0*X)
      GO TO 380
C ****************************************************
C
  19  FEX = X**5
      GO TO 380
C ********************************************************************
C
   20 FEX = X/2.D0*DSIN(X)
      GO TO 380
C
C ********************************************************************
C
  21  FEX = DEXP(-X) - DEXP(-1000.D0*X)
      GO TO 380
C
C ********************************************************************
C
   22 FEX = DCOS(X)
      GO TO 380
C ********************************************************************
C
  23  FEX = X*DEXP(X/4.)
      GO TO 380
C ********************************************************************
C
  24  FEX = 2.D0*DSQRT(X/PI)
      GO TO 380
C ********************************************************************
C
  25  FEX = DEXP(-X)/DSQRT(PI*X)
      GO TO 380
C ********************************************************************
C
  26  FEX = (1.D0+4.D0*X)/DSQRT(PI*X)
      GO TO 380
C ********************************************************************
C
  27  FEX = (DSIN(X)-X*DCOS(X))/2.D0
      GO TO 380
C ********************************************************************
C
  28  FEX = 1 - DEXP(-X)* (1+X)
      GO TO 380
C ****************************************************************
C
  29  FEX = DEXP(-X)/12.D0* (DEXP(3.D0*X)-DCOS(DSQRT(3.D0)*X)-
     +      DSQRT(3.D0)*DSIN(DSQRT(3.D0)*X))
      GO TO 380
c ********************************************************************
C
   30 FEX = 2.D0* (DCOS(2.D0*X)-DCOS(X))/X
      GO TO 380
C ********************************************************************
C
  31  FEX = (1-DEXP(-X))/X
      GO TO 380

c *******************************************************************
  32  FEX=-EULERO - DLOG(X)
      GO TO 380

c *******************************************************************
C
  33  IF (X.GE.0.D0 .AND. X.LE.1.D0) THEN
          FEX = X
      ELSE
          FEX = 1.D0
      END IF
      GO TO 380
C ********************************************************************
C
  34  TROV = .FALSE.
      K = 0
  210 CONTINUE
      IF (X.GT.2*K .AND. X.LT.2*K+1) THEN
          FEX = 0.D0
          TROV = .TRUE.

      ELSE IF (X.GT.2*K+1 .AND. X.LT.2*K+2) THEN
          FEX = 1.D0
          TROV = .TRUE.
      END IF

      IF (X.EQ.DFLOAT(K)) THEN
          FEX = 0.5D0
          TROV = .TRUE.
      END IF

      K = K + 1
      IF (.NOT.TROV .AND. K.LE.49) GO TO 210
      GO TO 380
C *******************************************************************
C
  35  T = X
      PI2 = PI*PI
      A = (0.5D0- (DEXP(2.D0)/ (DEXP(2.D0)-1.))*DEXP(-T)) + 0.5D0
      N = 0.D0
      SUM = 0.D0
  370 CONTINUE
      N = N + 1.D0
      B = DSIN((N*PI*T)-DATAN(N*PI))
      B = B/ (N* (DSQRT((N*N*PI2)+1.D0)))
      SUM = SUM + B
      IF (N.LT.2.D+4) GO TO 370
      SUM = SUM/PI
      SUM = SUM*DEXP(-T)
      FEX = A - SUM
      GO TO 380
C ******************************************************
C36    EPSABS = 0.0D0
C       EPSREL = 1.0D-04
C       A = 0.0D0
C       INF = 1
C       KOUNT = 0
C       IFAIL = -1
C       CALL D01AMF(FST3,A,INF,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
C     +            IFAIL)
C       IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
C       FEX=(RESULT/PI)+0.5d0
C       goto 380
C99996 FORMAT (1X,A,I4)
C *****************************************************

C
C
C37    EPSABS = 0.0D0
C      EPSREL = 1.0D-04
C      U1=2.*DSQRT(2.D0-DSQRT(3.D0))
C      U2=2.*DSQRT(2.D0+DSQRT(3.D0))
C      A = .0D0
C      B =U1
C      KOUNT = 0
C      IFAIL = -1
C      CALL D01AKF(FST4,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
C     *            IFAIL)
C      IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
C      DUFFY6=RESULT
C      A=U2
C      B=4.D0
C      IFAIL=-1
C       KOUNT=0
C      CALL D01AKF(FST4,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
C     *            IFAIL)
C      IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
C      FEX=(-DUFFY6+RESULT)/PI +1.D0
C      GO TO 380
C ******************************************************
C
C38    EPSABS = 0.0D0
C      EPSREL = 1.0D-04
C      A = 0.0D0
C      B = (1.-0.5d0)/0.5d0
C      KOUNT = 0
C      IFAIL = -1
C      CALL D01AJF(FST51,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
C    * IFAIL)
C      IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
C      PARTE1=RESULT
C      A=0.d0
C      B=DSQRT((1.-0.5d0)/1.5d0)
C      KOUNT=0
C      IFAIL=-1
C      CALL D01AJF(FST52,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
C    * IFAIL)
C      FEX=(RESULT+PARTE1)*(2.D0/PI)


C ******************************************************
  380 RETURN
C ******************************************************

      END
C ******************************************************
C ******************************************************
      DOUBLE PRECISION FUNCTION FST51(U)
      DOUBLE PRECISION              U,N,C,B
      DOUBLE PRECISION              COMX, R
      INTEGER                       KOUNT
      INTRINSIC                     SIN, SQRT
      COMMON                        /TELNUM/COMX,KOUNT
      KOUNT = KOUNT + 1
      N=0.5D0
      C=(1.-N)/N
      B=DSQRT((1.-N)/(1.+N))
      R=DSQRT(U*U + N*N*(C*C-U*U))
      FST51 = DCOSH(COMX*U)*
     *      ((U*DSQRT((R+U)/2))+
     *      DSQRT(C**2-U**2)*DSQRT(
     *      (R-U)/2))/
     *      (R*DSQRT(C**2-U**2)
     *            *DSQRT(U))
      RETURN
      END
C ******************************************************
      DOUBLE PRECISION FUNCTION FST52(u)
      DOUBLE PRECISION              U,N,C,B
      DOUBLE PRECISION              COMX
      INTEGER                       KOUNT
      INTRINSIC                     SIN, SQRT
      COMMON                        /TELNUM/COMX,KOUNT
      KOUNT = KOUNT + 1
      N=0.5d0
      C=(1.-N)/N
      B=DSQRT((1.-N)/(1.+N))
      FST52=DCOS(COMX*U)*((U-DSQRT(C**2+U**2))/(DSQRT(U)*
     * DSQRT(C**2+U**2)*DSQRT(N*DSQRT(C**2+U**2)-U)))
      RETURN
      END
C ******************************************************

      DOUBLE PRECISION FUNCTION FST3(U)
      DOUBLE PRECISION  M,COMX,THETA,ARG1,ARG2,R,C,U
      INTEGER                       KOUNT
      INTRINSIC                     SQRT
      COMMON                        /TELNUM/COMX,KOUNT
      KOUNT = KOUNT + 1
      R=0.5d0
      c=0.4d0
      M=(1.D0+U*U)/(1.D0+C*C*U*U)
      M=M**0.25
      THETA=ATAN(U)-ATAN(U*C)
      THETA=THETA/2
      ARG1=-R*M*DSQRT(U/2)*(DCOS(THETA)-DSIN(THETA))
      ARG2=COMX*u-R*M*DSQRT(U/2)*(DCOS(THETA)+DSIN(THETA))
      FST3=DEXP(ARG1)*DSIN(ARG2)/u
      RETURN
      END
C ******************************************************
      DOUBLE PRECISION FUNCTION FST4(U)
      DOUBLE PRECISION              COMX,U,U1,U2,K
      INTEGER                       KOUNT
      INTRINSIC                     SIN, SQRT,ACOS
      COMMON                        /TELNUM/comx,KOUNT
      KOUNT = KOUNT + 1
      U1=4.*(2.D0-DSQRT(3.D0))
      U2=4.*(2.D0+DSQRT(3.D0))
      IF((U1 - U**2 )*( U2 - U**2) .GE. 0.D0) THEN
      K=ACOS(.25D0*DSQRT((U1-U**2)*(U2-U**2)))
      FST4=(DSIN(U*COMX+2.*K)-DSIN(U*COMX-2*K))/U
      endif
      RETURN
      END
C ******************************************************

