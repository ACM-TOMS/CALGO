C ===== 4. DOUBLE PRECISION TESTING AIDS FOR UNCONSTRAINED NONLINEAR
C =====     OPTIMIZATION.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
C     FUNCTIONS DEFINED BY SUBROUTINE OBJFCN. THE SUBROUTINE RETURNS
C     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
C     THE SEVENTH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
C     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
C     THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
C         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
C         THE STANDARD STARTING POINT. IF FACTOR IS UNITY, NO
C         MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
C     .. Scalar Arguments ..
      DOUBLE PRECISION FACTOR
      INTEGER N,NPROB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C1,C2,C3,C4,FIVE,H,HALF,ONE,TEN,THREE,TWENTY,
     +                 TWNTF,TWO,ZERO
      INTEGER IVAR,J
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION DFLOAT
C     ..
C     .. Data statements ..
      DATA ZERO,HALF,ONE,TWO,THREE,FIVE,TEN,TWENTY,TWNTF/0.0D0,0.5D0,
     +     1.0D0,2.0D0,3.0D0,5.0D0,1.0D1,2.0D1,2.5D1/
      DATA C1,C2,C3,C4/4.0D-1,2.5D0,1.5D-1,1.2D0/
C     ..
C     .. Statement Function definitions ..
      DFLOAT(IVAR) = IVAR
C     ..
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,20,30,40,50,60,80,100,120,
     +       140,150,160,170,190,210,230,240,
     +       250) NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 270
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      X(1) = ONE
      X(2) = TWO
      X(3) = ONE
      X(4) = ONE
      X(5) = ONE
      X(6) = ONE
c      do i = 1, 6
c        X(i) = 5.0d0
c      enddo
      GO TO 270
C
C     GAUSSIAN FUNCTION.
C
   30 CONTINUE
      X(1) = C1
      X(2) = ONE
      X(3) = ZERO
      GO TO 270
C
C     POWELL BADLY SCALED FUNCTION.
C
   40 CONTINUE
      X(1) = ZERO
      X(2) = ONE
      GO TO 270
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   50 CONTINUE
      X(1) = ZERO
      X(2) = TEN
      X(3) = TWENTY
      GO TO 270
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   60 CONTINUE
      H = ONE/DFLOAT(N)
      DO 70 J = 1,N
          X(J) = ONE - DFLOAT(J)*H
   70 CONTINUE
      GO TO 270
C
C     WATSON FUNCTION.
C
   80 CONTINUE
      DO 90 J = 1,N
          X(J) = ZERO
   90 CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION I.
C
  100 CONTINUE
      DO 110 J = 1,N
          X(J) = DFLOAT(J)
  110 CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION II.
C
  120 CONTINUE
      DO 130 J = 1,N
          X(J) = HALF
  130 CONTINUE
      GO TO 270
C
C     BROWN BADLY SCALED FUNCTION.
C
  140 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     BROWN AND DENNIS FUNCTION.
C
  150 CONTINUE
      X(1) = TWNTF
      X(2) = FIVE
      X(3) = -FIVE
      X(4) = -ONE
      GO TO 270
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  160 CONTINUE
      X(1) = FIVE
      X(2) = C2
      X(3) = C3
      GO TO 270
C
C     TRIGONOMETRIC FUNCTION.
C
  170 CONTINUE
      H = ONE/DFLOAT(N)
      DO 180 J = 1,N
          X(J) = H
  180 CONTINUE
      GO TO 270
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  190 CONTINUE
      DO 200 J = 1,N,2
          X(J) = -C4
          X(J+1) = ONE
  200 CONTINUE
      GO TO 270
C
C     EXTENDED POWELL SINGULAR FUNCTION.
C
  210 CONTINUE
      DO 220 J = 1,N,4
          X(J) = THREE
          X(J+1) = -ONE
          X(J+2) = ZERO
          X(J+3) = ONE
  220 CONTINUE
      GO TO 270
C
C     BEALE FUNCTION.
C
  230 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     WOOD FUNCTION.
C
  240 CONTINUE
      X(1) = -THREE
      X(2) = -ONE
      X(3) = -THREE
      X(4) = -ONE
      GO TO 270
C
C     CHEBYQUAD FUNCTION.
C
  250 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 260 J = 1,N
          X(J) = DFLOAT(J)*H
  260 CONTINUE
  270 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR.EQ.ONE) GO TO 320
      IF (NPROB.EQ.7) GO TO 290
      DO 280 J = 1,N
          X(J) = FACTOR*X(J)
  280 CONTINUE
      GO TO 310

  290 CONTINUE
      DO 300 J = 1,N
          X(J) = FACTOR
  300 CONTINUE
  310 CONTINUE
  320 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE OBJFCN(N,X,F,NPROB)
C     **********
C
C     SUBROUTINE OBJFCN
C
C     THIS SUBROUTINE DEFINES THE OBJECTIVE FUNCTIONS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE VALUES
C     OF N FOR FUNCTIONS 1,2,3,4,5,10,11,12,16 AND 17 ARE
C     3,6,3,2,3,2,4,3,2 AND 4, RESPECTIVELY.
C     FOR FUNCTION 7, N MAY BE 2 OR GREATER BUT IS USUALLY 6 OR 9.
C     FOR FUNCTIONS 6,8,9,13,14,15 AND 18 N MAY BE VARIABLE,
C     HOWEVER IT MUST BE EVEN FOR FUNCTION 14, A MULTIPLE OF 4 FOR
C     FUNCTION 15, AND NOT GREATER THAN 50 FOR FUNCTION 18.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE OBJFCN(N,X,F,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       F IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
C         THE NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,
C                            DSQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER N,NPROB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AP,ARG,BP,C100,C10000,C1P5,C1PD6,C25,C29,C2P25,
     +                 C2P625,C2PDM6,C3P5,C90,CP0001,CP1,CP2,CP25,CP5,
     +                 D1,D2,EIGHT,FIFTY,FIVE,FOUR,ONE,R,S1,S2,S3,T,T1,
     +                 T2,T3,TEN,TH,THREE,TPI,TWO,ZERO
      INTEGER I,IEV,IVAR,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FVEC(50),Y(15)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,DSQRT
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION DFLOAT
C     ..
C     .. Data statements ..
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,FIFTY/0.0D0,1.0D0,
     +     2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,5.0D1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,C25,
     +     C29,C90,C100,C10000,C1PD6/2.0D-6,1.0D-4,1.0D-1,2.0D-1,2.5D-1,
     +     5.0D-1,1.5D0,2.25D0,2.625D0,3.5D0,2.5D1,2.9D1,9.0D1,1.0D2,
     +     1.0D4,1.0D6/
      DATA AP,BP/1.0D-5,1.0D0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     +     Y(12),Y(13),Y(14),Y(15)/9.0D-4,4.4D-3,1.75D-2,5.4D-2,
     +     1.295D-1,2.42D-1,3.521D-1,3.989D-1,3.521D-1,2.42D-1,1.295D-1,
     +     5.4D-2,1.75D-2,4.4D-3,9.0D-4/
C     ..
C     .. Statement Function definitions ..
      DFLOAT(IVAR) = IVAR
C     ..
C
C     FUNCTION ROUTINE SELECTOR.
C
      GO TO (10,20,40,60,70,90,110,150,170,
     +       200,210,230,250,280,300,320,330,
     +       340) NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TH = DSIGN(CP25,X(2))
      IF (X(1).GT.ZERO) TH = DATAN(X(2)/X(1))/TPI
      IF (X(1).LT.ZERO) TH = DATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = DSQRT(ARG)
      T = X(3) - TEN*TH
      F = C100* (T**2+ (R-ONE)**2) + X(3)**2
      GO TO 390
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      F = ZERO
      DO 30 I = 1,13
c      DO 30 I = 1, 6
          D1 = DFLOAT(I)/TEN
          D2 = DEXP(-D1) - FIVE*DEXP(-TEN*D1) + THREE*DEXP(-FOUR*D1)
          S1 = DEXP(-D1*X(1))
          S2 = DEXP(-D1*X(2))
          S3 = DEXP(-D1*X(5))
          T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
          F = F + T**2
   30 CONTINUE
      GO TO 390
C
C     GAUSSIAN FUNCTION.
C
   40 CONTINUE
      F = ZERO
      DO 50 I = 1,15
          D1 = CP5*DFLOAT(I-1)
          D2 = C3P5 - D1 - X(3)
          ARG = -CP5*X(2)*D2**2
          R = DEXP(ARG)
          T = X(1)*R - Y(I)
          F = F + T**2
   50 CONTINUE
      GO TO 390
C
C     POWELL BADLY SCALED FUNCTION.
C
   60 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = DEXP(-X(1))
      S2 = DEXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      F = T1**2 + T2**2
      GO TO 390
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   70 CONTINUE
      F = ZERO
      DO 80 I = 1,10
          D1 = DFLOAT(I)
          D2 = D1/TEN
          S1 = DEXP(-D2*X(1))
          S2 = DEXP(-D2*X(2))
          S3 = DEXP(-D2) - DEXP(-D1)
          T = S1 - S2 - S3*X(3)
          F = F + T**2
   80 CONTINUE
      GO TO 390
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   90 CONTINUE
      T1 = ZERO
      T2 = ZERO
      DO 100 J = 1,N
          T1 = T1 + DFLOAT(J)* (X(J)-ONE)
          T2 = T2 + (X(J)-ONE)**2
  100 CONTINUE
      F = T2 + T1**2* (ONE+T1**2)
      GO TO 390
C
C     WATSON FUNCTION.
C
  110 CONTINUE
      F = ZERO
      DO 140 I = 1,29
          D1 = DFLOAT(I)/C29
          S1 = ZERO
          D2 = ONE
          DO 120 J = 2,N
              S1 = S1 + DFLOAT(J-1)*D2*X(J)
              D2 = D1*D2
  120     CONTINUE
          S2 = ZERO
          D2 = ONE
          DO 130 J = 1,N
              S2 = S2 + D2*X(J)
              D2 = D1*D2
  130     CONTINUE
          T = S1 - S2**2 - ONE
          F = F + T**2
  140 CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      F = F + X(1)**2 + T1**2
      GO TO 390
C
C     PENALTY FUNCTION I.
C
  150 CONTINUE
      T1 = -CP25
      T2 = ZERO
      DO 160 J = 1,N
          T1 = T1 + X(J)**2
          T2 = T2 + (X(J)-ONE)**2
  160 CONTINUE
      F = AP*T2 + BP*T1**2
      GO TO 390
C
C     PENALTY FUNCTION II.
C
  170 CONTINUE
      T1 = -ONE
      T2 = ZERO
      T3 = ZERO
      D1 = DEXP(CP1)
      D2 = ONE
      DO 190 J = 1,N
          T1 = T1 + DFLOAT(N-J+1)*X(J)**2
          S1 = DEXP(X(J)/TEN)
          IF (J.EQ.1) GO TO 180
          S3 = S1 + S2 - D2* (D1+ONE)
          T2 = T2 + S3**2
          T3 = T3 + (S1-ONE/D1)**2
  180     CONTINUE
          S2 = S1
          D2 = D1*D2
  190 CONTINUE
      F = AP* (T2+T3) + BP* (T1**2+ (X(1)-CP2)**2)
      GO TO 390
C
C     BROWN BADLY SCALED FUNCTION.
C
  200 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     BROWN AND DENNIS FUNCTION.
C
  210 CONTINUE
      F = ZERO
      DO 220 I = 1,20
          D1 = DFLOAT(I)/FIVE
          D2 = DSIN(D1)
          T1 = X(1) + D1*X(2) - DEXP(D1)
          T2 = X(3) + D2*X(4) - DCOS(D1)
          T = T1**2 + T2**2
          F = F + T**2
  220 CONTINUE
      GO TO 390
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  230 CONTINUE
      F = ZERO
      D1 = TWO/THREE
      DO 240 I = 1,99
          ARG = DFLOAT(I)/C100
          R = DABS((-FIFTY*DLOG(ARG))**D1+C25-X(2))
          T1 = R**X(3)/X(1)
          T2 = DEXP(-T1)
          T = T2 - ARG
          F = F + T**2
  240 CONTINUE
      GO TO 390
C
C     TRIGONOMETRIC FUNCTION.
C
  250 CONTINUE
      S1 = ZERO
      DO 260 J = 1,N
          S1 = S1 + DCOS(X(J))
  260 CONTINUE
      F = ZERO
      DO 270 J = 1,N
          T = DFLOAT(N+J) - DSIN(X(J)) - S1 - DFLOAT(J)*DCOS(X(J))
          F = F + T**2
  270 CONTINUE
      GO TO 390
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  280 CONTINUE
      F = ZERO
      DO 290 J = 1,N,2
          T1 = ONE - X(J)
          T2 = TEN* (X(J+1)-X(J)**2)
          F = F + T1**2 + T2**2
  290 CONTINUE
      GO TO 390
C
C     EXTENDED POWELL FUNCTION.
C
  300 CONTINUE
      F = ZERO
      DO 310 J = 1,N,4
          T = X(J) + TEN*X(J+1)
          T1 = X(J+2) - X(J+3)
          S1 = FIVE*T1
          T2 = X(J+1) - TWO*X(J+2)
          S2 = T2**3
          T3 = X(J) - X(J+3)
          S3 = TEN*T3**3
          F = F + T**2 + S1*T1 + S2*T2 + S3*T3
  310 CONTINUE
      GO TO 390
C
C     BEALE FUNCTION.
C
  320 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     WOOD FUNCTION.
C
  330 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      F = C100*S1**2 + S2**2 + C90*T1**2 + T2**2 + TEN* (S3+T3)**2 +
     +    (S3-T3)**2/TEN
      GO TO 390
C
C     CHEBYQUAD FUNCTION.
C
  340 CONTINUE
      DO 350 I = 1,N
          FVEC(I) = ZERO
  350 CONTINUE
      DO 370 J = 1,N
          T1 = ONE
          T2 = TWO*X(J) - ONE
          T = TWO*T2
          DO 360 I = 1,N
              FVEC(I) = FVEC(I) + T2
              TH = T*T2 - T1
              T1 = T2
              T2 = TH
  360     CONTINUE
  370 CONTINUE
      F = ZERO
      D1 = ONE/DFLOAT(N)
      IEV = -1
      DO 380 I = 1,N
          T = D1*FVEC(I)
          IF (IEV.GT.0) T = T + ONE/ (DFLOAT(I)**2-ONE)
          F = F + T**2
          IEV = -IEV
  380 CONTINUE
  390 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE OBJFCN.
C
      END
      SUBROUTINE GRDFCN(N,X,G,NPROB)
C     **********
C
C     SUBROUTINE GRDFCN
C
C     THIS SUBROUTINE DEFINES THE GRADIENT VECTORS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN THE PROLOGUE COMMENTS OF OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE GRDFCN(N,X,G,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       G IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE COMPONENTS
C         OF THE GRADIENT VECTOR OF THE NPROB OBJECTIVE FUNCTION
C         EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,
C                            DSQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
C     .. Scalar Arguments ..
      INTEGER N,NPROB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION G(N),X(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AP,ARG,BP,C100,C10000,C180,C19P8,C1P5,C1PD6,C200,
     +                 C20P2,C25,C29,C2P25,C2P625,C2PDM6,C3P5,CP0001,
     +                 CP1,CP2,CP25,CP5,D1,D2,EIGHT,FIFTY,FIVE,FOUR,ONE,
     +                 R,S1,S2,S3,T,T1,T2,T3,TEN,TH,THREE,TPI,TWENTY,
     +                 TWO,ZERO
      INTEGER I,IEV,IVAR,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FVEC(50),Y(15)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,DSQRT
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION DFLOAT
C     ..
C     .. Data statements ..
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,TWENTY,FIFTY/0.0D0,
     +     1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,2.0D1,5.0D1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,C19P8,
     +     C20P2,C25,C29,C100,C180,C200,C10000,C1PD6/2.0D-6,1.0D-4,
     +     1.0D-1,2.0D-1,2.5D-1,5.0D-1,1.5D0,2.25D0,2.625D0,3.5D0,
     +     1.98D1,2.02D1,2.5D1,2.9D1,1.0D2,1.8D2,2.0D2,1.0D4,1.0D6/
      DATA AP,BP/1.0D-5,1.0D0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     +     Y(12),Y(13),Y(14),Y(15)/9.0D-4,4.4D-3,1.75D-2,5.4D-2,
     +     1.295D-1,2.42D-1,3.521D-1,3.989D-1,3.521D-1,2.42D-1,1.295D-1,
     +     5.4D-2,1.75D-2,4.4D-3,9.0D-4/
C     ..
C     .. Statement Function definitions ..
      DFLOAT(IVAR) = IVAR
C     ..
C
C     GRADIENT ROUTINE SELECTOR.
C
      GO TO (10,20,50,70,80,100,130,190,220,
     +       260,270,290,310,350,370,390,400,
     +       410) NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TH = DSIGN(CP25,X(2))
      IF (X(1).GT.ZERO) TH = DATAN(X(2)/X(1))/TPI
      IF (X(1).LT.ZERO) TH = DATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = DSQRT(ARG)
      T = X(3) - TEN*TH
      S1 = TEN*T/ (TPI*ARG)
      G(1) = C200* (X(1)-X(1)/R+X(2)*S1)
      G(2) = C200* (X(2)-X(2)/R-X(1)*S1)
      G(3) = TWO* (C100*T+X(3))
      GO TO 490
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      DO 30 J = 1,N
          G(J) = ZERO
   30 CONTINUE
      DO 40 I = 1,13
c      DO 40 I = 1, 6
          D1 = DFLOAT(I)/TEN
          D2 = DEXP(-D1) - FIVE*DEXP(-TEN*D1) + THREE*DEXP(-FOUR*D1)
          S1 = DEXP(-D1*X(1))
          S2 = DEXP(-D1*X(2))
          S3 = DEXP(-D1*X(5))
          T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
          TH = D1*T
          G(1) = G(1) - S1*TH
          G(2) = G(2) + S2*TH
          G(3) = G(3) + S1*T
          G(4) = G(4) - S2*T
          G(5) = G(5) - S3*TH
          G(6) = G(6) + S3*T
   40 CONTINUE
      G(1) = TWO*X(3)*G(1)
      G(2) = TWO*X(4)*G(2)
      G(3) = TWO*G(3)
      G(4) = TWO*G(4)
      G(5) = TWO*X(6)*G(5)
      G(6) = TWO*G(6)
      GO TO 490
C
C     GAUSSIAN FUNCTION.
C
   50 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 60 I = 1,15
          D1 = CP5*DFLOAT(I-1)
          D2 = C3P5 - D1 - X(3)
          ARG = -CP5*X(2)*D2**2
          R = DEXP(ARG)
          T = X(1)*R - Y(I)
          S1 = R*T
          S2 = D2*S1
          G(1) = G(1) + S1
          G(2) = G(2) - D2*S2
          G(3) = G(3) + S2
   60 CONTINUE
      G(1) = TWO*G(1)
      G(2) = X(1)*G(2)
      G(3) = TWO*X(1)*X(2)*G(3)
      GO TO 490
C
C     POWELL BADLY SCALED FUNCTION.
C
   70 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = DEXP(-X(1))
      S2 = DEXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      G(1) = TWO* (C10000*X(2)*T1-S1*T2)
      G(2) = TWO* (C10000*X(1)*T1-S2*T2)
      GO TO 490
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   80 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 90 I = 1,10
          D1 = DFLOAT(I)
          D2 = D1/TEN
          S1 = DEXP(-D2*X(1))
          S2 = DEXP(-D2*X(2))
          S3 = DEXP(-D2) - DEXP(-D1)
          T = S1 - S2 - S3*X(3)
          TH = D2*T
          G(1) = G(1) - S1*TH
          G(2) = G(2) + S2*TH
          G(3) = G(3) - S3*T
   90 CONTINUE
      G(1) = TWO*G(1)
      G(2) = TWO*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  100 CONTINUE
      T1 = ZERO
      DO 110 J = 1,N
          T1 = T1 + DFLOAT(J)* (X(J)-ONE)
  110 CONTINUE
      T = T1* (ONE+TWO*T1**2)
      DO 120 J = 1,N
          G(J) = TWO* (X(J)-ONE+DFLOAT(J)*T)
  120 CONTINUE
      GO TO 490
C
C     WATSON FUNCTION.
C
  130 CONTINUE
      DO 140 J = 1,N
          G(J) = ZERO
  140 CONTINUE
      DO 180 I = 1,29
          D1 = DFLOAT(I)/C29
          S1 = ZERO
          D2 = ONE
          DO 150 J = 2,N
              S1 = S1 + DFLOAT(J-1)*D2*X(J)
              D2 = D1*D2
  150     CONTINUE
          S2 = ZERO
          D2 = ONE
          DO 160 J = 1,N
              S2 = S2 + D2*X(J)
              D2 = D1*D2
  160     CONTINUE
          T = S1 - S2**2 - ONE
          S3 = TWO*D1*S2
          D2 = TWO/D1
          DO 170 J = 1,N
              G(J) = G(J) + D2* (DFLOAT(J-1)-S3)*T
              D2 = D1*D2
  170     CONTINUE
  180 CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      G(1) = G(1) + X(1)* (TWO-FOUR*T1)
      G(2) = G(2) + TWO*T1
      GO TO 490
C
C     PENALTY FUNCTION I.
C
  190 CONTINUE
      T1 = -CP25
      DO 200 J = 1,N
          T1 = T1 + X(J)**2
  200 CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      DO 210 J = 1,N
          G(J) = D1* (X(J)-ONE) + X(J)*TH
  210 CONTINUE
      GO TO 490
C
C     PENALTY FUNCTION II.
C
  220 CONTINUE
      T1 = -ONE
      DO 230 J = 1,N
          T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  230 CONTINUE
      D1 = DEXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      DO 250 J = 1,N
          G(J) = DFLOAT(N-J+1)*X(J)*TH
          S1 = DEXP(X(J)/TEN)
          IF (J.EQ.1) GO TO 240
          S3 = S1 + S2 - D2* (D1+ONE)
          G(J) = G(J) + AP*S1* (S3+S1-ONE/D1)/FIVE
          G(J-1) = G(J-1) + AP*S2*S3/FIVE
  240     CONTINUE
          S2 = S1
          D2 = D1*D2
  250 CONTINUE
      G(1) = G(1) + TWO*BP* (X(1)-CP2)
      GO TO 490
C
C     BROWN BADLY SCALED FUNCTION.
C
  260 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      G(1) = TWO* (T1+X(2)*T3)
      G(2) = TWO* (T2+X(1)*T3)
      GO TO 490
C
C     BROWN AND DENNIS FUNCTION.
C
  270 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      G(4) = ZERO
      DO 280 I = 1,20
          D1 = DFLOAT(I)/FIVE
          D2 = DSIN(D1)
          T1 = X(1) + D1*X(2) - DEXP(D1)
          T2 = X(3) + D2*X(4) - DCOS(D1)
          T = T1**2 + T2**2
          S1 = T1*T
          S2 = T2*T
          G(1) = G(1) + S1
          G(2) = G(2) + D1*S1
          G(3) = G(3) + S2
          G(4) = G(4) + D2*S2
  280 CONTINUE
      G(1) = FOUR*G(1)
      G(2) = FOUR*G(2)
      G(3) = FOUR*G(3)
      G(4) = FOUR*G(4)
      GO TO 490
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  290 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      D1 = TWO/THREE
      DO 300 I = 1,99
          ARG = DFLOAT(I)/C100
          R = DABS((-FIFTY*DLOG(ARG))**D1+C25-X(2))
          T1 = R**X(3)/X(1)
          T2 = DEXP(-T1)
          T = T2 - ARG
          S1 = T1*T2*T
          G(1) = G(1) + S1
          G(2) = G(2) + S1/R
          G(3) = G(3) - S1*DLOG(R)
  300 CONTINUE
      G(1) = TWO*G(1)/X(1)
      G(2) = TWO*X(3)*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     TRIGONOMETRIC FUNCTION.
C
  310 CONTINUE
      S1 = ZERO
      DO 320 J = 1,N
          G(J) = DCOS(X(J))
          S1 = S1 + G(J)
  320 CONTINUE
      S2 = ZERO
      DO 330 J = 1,N
          TH = DSIN(X(J))
          T = DFLOAT(N+J) - TH - S1 - DFLOAT(J)*G(J)
          S2 = S2 + T
          G(J) = (DFLOAT(J)*TH-G(J))*T
  330 CONTINUE
      DO 340 J = 1,N
          G(J) = TWO* (G(J)+DSIN(X(J))*S2)
  340 CONTINUE
      GO TO 490
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  350 CONTINUE
      DO 360 J = 1,N,2
          T1 = ONE - X(J)
          G(J+1) = C200* (X(J+1)-X(J)**2)
          G(J) = -TWO* (X(J)*G(J+1)+T1)
  360 CONTINUE
      GO TO 490
C
C     EXTENDED POWELL FUNCTION.
C
  370 CONTINUE
      DO 380 J = 1,N,4
          T = X(J) + TEN*X(J+1)
          T1 = X(J+2) - X(J+3)
          S1 = FIVE*T1
          T2 = X(J+1) - TWO*X(J+2)
          S2 = FOUR*T2**3
          T3 = X(J) - X(J+3)
          S3 = TWENTY*T3**3
          G(J) = TWO* (T+S3)
          G(J+1) = TWENTY*T + S2
          G(J+2) = TWO* (S1-S2)
          G(J+3) = -TWO* (S1+S3)
  380 CONTINUE
      GO TO 490
C
C     BEALE FUNCTION.
C
  390 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      G(1) = -TWO* (S1*T1+S2*T2+S3*T3)
      G(2) = TWO*X(1)* (T1+X(2)* (TWO*T2+THREE*X(2)*T3))
      GO TO 490
C
C     WOOD FUNCTION.
C
  400 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      G(1) = -TWO* (C200*X(1)*S1+S2)
      G(2) = C200*S1 + C20P2*S3 + C19P8*T3
      G(3) = -TWO* (C180*X(3)*T1+T2)
      G(4) = C180*T1 + C20P2*T3 + C19P8*S3
      GO TO 490
C
C     CHEBYQUAD FUNCTION.
C
  410 CONTINUE
      DO 420 I = 1,N
          FVEC(I) = ZERO
  420 CONTINUE
      DO 440 J = 1,N
          T1 = ONE
          T2 = TWO*X(J) - ONE
          T = TWO*T2
          DO 430 I = 1,N
              FVEC(I) = FVEC(I) + T2
              TH = T*T2 - T1
              T1 = T2
              T2 = TH
  430     CONTINUE
  440 CONTINUE
      D1 = ONE/DFLOAT(N)
      IEV = -1
      DO 450 I = 1,N
          FVEC(I) = D1*FVEC(I)
          IF (IEV.GT.0) FVEC(I) = FVEC(I) + ONE/ (DFLOAT(I)**2-ONE)
          IEV = -IEV
  450 CONTINUE
      DO 470 J = 1,N
          G(J) = ZERO
          T1 = ONE
          T2 = TWO*X(J) - ONE
          T = TWO*T2
          S1 = ZERO
          S2 = TWO
          DO 460 I = 1,N
              G(J) = G(J) + FVEC(I)*S2
              TH = FOUR*T2 + T*S2 - S1
              S1 = S2
              S2 = TH
              TH = T*T2 - T1
              T1 = T2
              T2 = TH
  460     CONTINUE
  470 CONTINUE
      D2 = TWO*D1
      DO 480 J = 1,N
          G(J) = D2*G(J)
  480 CONTINUE
  490 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE GRDFCN.
C
      END
      SUBROUTINE HESFCN(N,X,HESD,HESL,NPROB)
C     **********
C
C     SUBROUTINE HESFCN
C
C     THIS SUBROUTINE DEFINES THE HESSIAN MATRICES OF 18
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.  THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE HESFCN (N, X, HESD, HESL, NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       HESD IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         DIAGONAL COMPONENTS OF THE HESSIAN MATRIX OF THE NPROB
C         OBJECTIVE FUNCTION EVALUATED AT X.
C
C       HESL IS AN OUTPUT ARRAY OF LENGTH N*(N-1)/2 WHICH CONTAINS
C         THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX OF THE
C         NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM.  NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS, ATAN, COS, EXP, LOG, SIGN, SIN,
C                            SQRT
C
C       INTEGER INLINE FUNCTION IX GIVES THE LOCATION OF A HESSIAN
C       ELEMENT (I,J), I>J, IN HESL
C
C     VICTORIA Z. AVERBUKH, SAMUEL A. FIGUEROA, AND
C     TAMAR SCHLICK, 1993.
C     **********
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,EIGHT,NINE,TEN,
     +                 FIFTY,CP0001,CP1,CP25,CP5,C1P5,C2P25,C2P625,C3P5,
     +                 C12,C19P8,C25,C29,C50,C100,C120,C200,C200P2,C202,
     +                 C220P2,C360,C400,C1000,C1080,C1200,C2000,C20000,
     +                 C2E8,C4E8,AP,BP,PI
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0,
     +          FIVE=5.0D0,SIX=6.0D0,EIGHT=8.0D0,NINE=9.0D0,TEN=1.0D1,
     +          FIFTY=5.0D1,CP0001=1.0D-4,CP1=1.0D-1,CP25=2.5D-1,
     +          CP5=5.0D-1,C1P5=1.5D0,C2P25=2.25D0,C2P625=2.625D0,
     +          C3P5=3.5D0,C12=1.2D1,C19P8=1.98D1,C25=2.5D1,C29=2.9D1,
     +          C50=5.0D1,C100=1.0D2,C120=1.2D2,C200=2.0D2,
     +          C200P2=2.002D2,C202=2.02D2,C220P2=2.202D2,C360=3.6D2,
     +          C400=4.0D2,C1000=1.0D3,C1080=1.08D3,C1200=1.2D3,
     +          C2000=2.0D3,C20000=2.0D4,C2E8=2.0D8,C4E8=4.0D8,
     +          AP=1.0D-5,BP=ONE,PI=3.141592653589793D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NPROB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION HESD(N),HESL(*),X(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ARG,D1,D2,D3,LOGR,P1,P2,PIARG,PIARG2,R,R3INV,S1,
     +                 S1S2,S1S3,S2,S2S3,S3,SS1,SS2,T,T1,T2,T3,TH,TT,
     +                 TT1,TT2,TTH
      INTEGER I,II,IVAR,J,JJ,K,M
      LOGICAL IEV
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FVEC(50),GVEC(50),Y(15)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,EXP,FLOAT,LOG,SIGN,SIN,SQRT
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION DFLOAT
      INTEGER IX
C     ..
C     .. Statement Function definitions ..

      IX(II,JJ) = (II-1)* (II-2)/2 + JJ
      DFLOAT(IVAR) = IVAR
C     ..
C     .. Data statements ..
      DATA Y/9.0D-4,4.4D-3,1.75D-2,5.4D-2,1.295D-1,2.42D-1,3.521D-1,
     +     3.989D-1,3.521D-1,2.42D-1,1.295D-1,5.4D-2,1.75D-2,4.4D-3,
     +     9.0D-4/
C     ..

C
C     HESSIAN ROUTINE SELECTOR.
C
      GO TO (10,20,80,100,110,130,170,260,300,
     +       340,350,390,430,480,510,540,550,
     +       560) NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
C
      IF (X(1).EQ.ZERO) THEN
          TH = SIGN(CP25,X(2))

      ELSE
          TH = ATAN(X(2)/X(1))/ (TWO*PI)
          IF (X(1).LT.ZERO) TH = TH + CP5
      END IF

      ARG = X(1)**2 + X(2)**2
      PIARG = PI*ARG
      PIARG2 = PIARG*ARG
      R3INV = ONE/SQRT(ARG)**3
      T = X(3) - TEN*TH
      S1 = FIVE*T/PIARG
      P1 = C2000*X(1)*X(2)*T/PIARG2
      P2 = (FIVE/PIARG)**2
      HESD(1) = C200 - C200* (R3INV-P2)*X(2)**2 - P1
      HESD(2) = C200 - C200* (R3INV-P2)*X(1)**2 + P1
      HESD(3) = C202
      HESL(1) = C200*X(1)*X(2)*R3INV +
     +          C1000/PIARG2* (T* (X(1)**2-X(2)**2)-FIVE*X(1)*X(2)/PI)
      HESL(2) = C1000*X(2)/PIARG
      HESL(3) = -C1000*X(1)/PIARG
      RETURN
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      DO 30 I = 1,6
          HESD(I) = ZERO
   30 CONTINUE
      DO 40 I = 1,15
          HESL(I) = ZERO
   40 CONTINUE
c      DO 230 I = 1, 6
      DO 50 I = 1,13
          D1 = DFLOAT(I)/TEN
          D2 = EXP(-D1) - FIVE*EXP(-TEN*D1) + THREE*EXP(-FOUR*D1)
          S1 = EXP(-D1*X(1))
          S2 = EXP(-D1*X(2))
          S3 = EXP(-D1*X(5))
          T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
          D2 = D1**2
          S1S2 = S1*S2
          S1S3 = S1*S3
          S2S3 = S2*S3
          HESD(1) = HESD(1) + D2*S1* (T+X(3)*S1)
          HESD(2) = HESD(2) - D2*S2* (T-X(4)*S2)
          HESD(3) = HESD(3) + S1**2
          HESD(4) = HESD(4) + S2**2
          HESD(5) = HESD(5) + D2*S3* (T+X(6)*S3)
          HESD(6) = HESD(6) + S3**2
          HESL(1) = HESL(1) - D2*S1S2
          HESL(2) = HESL(2) - D1*S1* (T+X(3)*S1)
          HESL(3) = HESL(3) + D1*S1S2
          HESL(4) = HESL(4) + D1*S1S2
          HESL(5) = HESL(5) + D1*S2* (T-X(4)*S2)
          HESL(6) = HESL(6) - S1S2
          HESL(7) = HESL(7) + D2*S1S3
          HESL(8) = HESL(8) - D2*S2S3
          HESL(9) = HESL(9) - D1*S1S3
          HESL(10) = HESL(10) + D1*S2S3
          HESL(11) = HESL(11) - D1*S1S3
          HESL(12) = HESL(12) + D1*S2S3
          HESL(13) = HESL(13) + S1S3
          HESL(14) = HESL(14) - S2S3
          HESL(15) = HESL(15) - D1*S3* (T+X(6)*S3)
   50 CONTINUE
      HESD(1) = X(3)*HESD(1)
      HESD(2) = X(4)*HESD(2)
      HESD(5) = X(6)*HESD(5)
      HESL(1) = X(3)*X(4)*HESL(1)
      HESL(3) = X(4)*HESL(3)
      HESL(4) = X(3)*HESL(4)
      HESL(7) = X(3)*X(6)*HESL(7)
      HESL(8) = X(4)*X(6)*HESL(8)
      HESL(9) = X(6)*HESL(9)
      HESL(10) = X(6)*HESL(10)
      HESL(11) = X(3)*HESL(11)
      HESL(12) = X(4)*HESL(12)
      DO 60 I = 1,6
          HESD(I) = TWO*HESD(I)
   60 CONTINUE
      DO 70 I = 1,15
          HESL(I) = TWO*HESL(I)
   70 CONTINUE
      RETURN
C
C     GAUSSIAN FUNCTION.
C
   80 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 90 I = 1,15
          D1 = CP5*DFLOAT(I-1)
          D2 = C3P5 - D1 - X(3)
          ARG = CP5*X(2)*D2**2
          R = EXP(-ARG)
          T = X(1)*R - Y(I)
          T1 = TWO*X(1)*R - Y(I)
          HESD(1) = HESD(1) + R**2
          HESD(2) = HESD(2) + R*T1*D2**4
          HESD(3) = HESD(3) + R* (X(2)*T1*D2**2-T)
          HESL(1) = HESL(1) - R*T1*D2**2
          HESL(2) = HESL(2) + D2*R*T1
          HESL(3) = HESL(3) + D2*R* (T-ARG*T1)
   90 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = CP5*X(1)*HESD(2)
      HESD(3) = TWO*X(1)*X(2)*HESD(3)
      HESL(2) = TWO*X(2)*HESL(2)
      HESL(3) = TWO*X(1)*HESL(3)
      RETURN
C
C     POWELL BADLY SCALED FUNCTION.
C
  100 CONTINUE
      S1 = EXP(-X(1))
      S2 = EXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      HESD(1) = C2E8*X(2)**2 + TWO*S1* (S1+T2)
      HESD(2) = C2E8*X(1)**2 + TWO*S2* (S2+T2)
      HESL(1) = C4E8*X(1)*X(2) + TWO*S1*S2 - C20000
      RETURN
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  110 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 120 I = 1,10
          D1 = DFLOAT(I)
          D2 = D1/TEN
          S1 = EXP(-D2*X(1))
          S2 = EXP(-D2*X(2))
          S3 = EXP(-D2) - EXP(-D1)
          T = S1 - S2 - S3*X(3)
          TH = T*D2**2
          HESD(1) = HESD(1) + TH*S1 + (D2*S1)**2
          HESD(2) = HESD(2) - TH*S2 + (D2*S2)**2
          HESD(3) = HESD(3) + S3**2
          HESL(1) = HESL(1) - S1*S2*D2**2
          HESL(2) = HESL(2) + D2*S1*S3
          HESL(3) = HESL(3) - D2*S2*S3
  120 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = TWO*HESD(2)
      HESD(3) = TWO*HESD(3)
      HESL(1) = TWO*HESL(1)
      HESL(2) = TWO*HESL(2)
      HESL(3) = TWO*HESL(3)
      RETURN
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  130 CONTINUE
      T1 = ZERO
      DO 140 J = 1,N
          T1 = T1 + DFLOAT(J)* (X(J)-ONE)
  140 CONTINUE
      T = ONE + SIX*T1**2
      M = 0
      DO 160 J = 1,N
          HESD(J) = TWO + TWO*T*DFLOAT(J)**2
          DO 150 K = 1,J - 1
              M = M + 1
              HESL(M) = TWO*T*DFLOAT(J*K)
  150     CONTINUE
  160 CONTINUE
      RETURN
C
C     WATSON FUNCTION.
C
  170 CONTINUE
      DO 180 J = 1,N
          HESD(J) = ZERO
  180 CONTINUE
      DO 190 J = 1,N* (N-1)/2
          HESL(J) = ZERO
  190 CONTINUE
      DO 230 I = 1,29
          D1 = DFLOAT(I)/C29
          D2 = ONE
          S1 = ZERO
          S2 = X(1)
          DO 200 J = 2,N
              S1 = S1 + DFLOAT(J-1)*D2*X(J)
              D2 = D1*D2
              S2 = S2 + D2*X(J)
  200     CONTINUE
          T = TWO* (S1-S2**2-ONE)*D1**2
          S3 = TWO*D1*S2
          D2 = ONE/D1
          M = 0
          DO 220 J = 1,N
              T1 = DFLOAT(J-1) - S3
              HESD(J) = HESD(J) + (T1**2-T)*D2**2
              D3 = ONE/D1
              DO 210 K = 1,J - 1
                  M = M + 1
                  HESL(M) = HESL(M) + (T1* (DFLOAT(K-1)-S3)-T)*D2*D3
                  D3 = D1*D3
  210         CONTINUE
              D2 = D1*D2
  220     CONTINUE
  230 CONTINUE
      T3 = X(2) - X(1)**2 - ONE
      HESD(1) = HESD(1) + ONE - TWO* (T3-TWO*X(1)**2)
      HESD(2) = HESD(2) + ONE
      HESL(1) = HESL(1) - TWO*X(1)
      DO 240 J = 1,N
          HESD(J) = TWO*HESD(J)
  240 CONTINUE
      DO 250 J = 1,N* (N-1)/2
          HESL(J) = TWO*HESL(J)
  250 CONTINUE
      RETURN
C
C     PENALTY FUNCTION I.
C
  260 CONTINUE
      T1 = -CP25
      DO 270 J = 1,N
          T1 = T1 + X(J)**2
  270 CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      M = 0
      DO 290 J = 1,N
          HESD(J) = D1 + TH + EIGHT*X(J)**2
          DO 280 K = 1,J - 1
              M = M + 1
              HESL(M) = EIGHT*X(J)*X(K)
  280     CONTINUE
  290 CONTINUE
      RETURN
C
C     PENALTY FUNCTION II.
C
  300 CONTINUE
      T1 = -ONE
      DO 310 J = 1,N
          T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  310 CONTINUE
      D1 = EXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      M = 0
      DO 330 J = 1,N
          HESD(J) = EIGHT*BP* (DFLOAT(N-J+1)*X(J))**2 + DFLOAT(N-J+1)*TH
          S1 = EXP(X(J)/TEN)
          IF (J.GT.1) THEN
              S3 = S1 + S2 - D2* (D1+ONE)
              HESD(J) = HESD(J) + AP*S1* (S3+S1-ONE/D1+TWO*S1)/C50
              HESD(J-1) = HESD(J-1) + AP*S2* (S2+S3)/C50
              DO 320 K = 1,J - 1
                  M = M + 1
                  T1 = EXP(DFLOAT(K)/TEN)
                  HESL(M) = EIGHT*DFLOAT(N-J+1)*DFLOAT(N-K+1)*X(J)*X(K)
  320         CONTINUE
              HESL(M) = HESL(M) + AP*S1*S2/C50
          END IF

          S2 = S1
          D2 = D1*D2
  330 CONTINUE
      HESD(1) = HESD(1) + TWO*BP
      RETURN
C
C     BROWN BADLY SCALED FUNCTION.
C
  340 CONTINUE
      HESD(1) = TWO + TWO*X(2)**2
      HESD(2) = TWO + TWO*X(1)**2
      HESL(1) = FOUR*X(1)*X(2) - FOUR
      RETURN
C
C     BROWN AND DENNIS FUNCTION.
C
  350 CONTINUE
      DO 360 I = 1,4
          HESD(I) = ZERO
  360 CONTINUE
      DO 370 I = 1,6
          HESL(I) = ZERO
  370 CONTINUE
      DO 380 I = 1,20
          D1 = DFLOAT(I)/FIVE
          D2 = SIN(D1)
          T1 = X(1) + D1*X(2) - EXP(D1)
          T2 = X(3) + D2*X(4) - COS(D1)
          T = EIGHT*T1*T2
          S1 = C12*T1**2 + FOUR*T2**2
          S2 = C12*T2**2 + FOUR*T1**2
          HESD(1) = HESD(1) + S1
          HESD(2) = HESD(2) + S1*D1**2
          HESD(3) = HESD(3) + S2
          HESD(4) = HESD(4) + S2*D2**2
          HESL(1) = HESL(1) + S1*D1
          HESL(2) = HESL(2) + T
          HESL(4) = HESL(4) + T*D2
          HESL(3) = HESL(3) + T*D1
          HESL(5) = HESL(5) + T*D1*D2
          HESL(6) = HESL(6) + S2*D2
  380 CONTINUE
      RETURN
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  390 CONTINUE
      DO 400 I = 1,3
          HESD(I) = ZERO
          HESL(I) = ZERO
  400 CONTINUE
      D1 = TWO/THREE
      DO 410 I = 1,99
          ARG = DFLOAT(I)/C100
          R = (-FIFTY*LOG(ARG))**D1 + C25 - X(2)
          T1 = ABS(R)**X(3)/X(1)
          T2 = EXP(-T1)
          T3 = T1*T2* (T1*T2+ (T1-ONE)* (T2-ARG))
          T = T1*T2* (T2-ARG)
          LOGR = LOG(ABS(R))
          HESD(1) = HESD(1) + T3 - T
          HESD(2) = HESD(2) + (T+X(3)*T3)/R**2
          HESD(3) = HESD(3) + T3*LOGR**2
          HESL(1) = HESL(1) + T3/R
          HESL(2) = HESL(2) - T3*LOGR
          HESL(3) = HESL(3) + (T-X(3)*T3*LOGR)/R
  410 CONTINUE
      HESD(1) = HESD(1)/X(1)**2
      HESD(2) = HESD(2)*X(3)
      HESL(1) = HESL(1)*X(3)/X(1)
      HESL(2) = HESL(2)/X(1)
      DO 420 I = 1,3
          HESD(I) = TWO*HESD(I)
          HESL(I) = TWO*HESL(I)
  420 CONTINUE
      RETURN
C
C     TRIGONOMETRIC FUNCTION.
C
  430 CONTINUE
      S1 = ZERO
      DO 440 J = 1,N
          HESD(J) = SIN(X(J))
          S1 = S1 + COS(X(J))
  440 CONTINUE
      S2 = ZERO
      M = 0
      DO 460 J = 1,N
          TH = COS(X(J))
          T = DFLOAT(N+J) - HESD(J) - S1 - DFLOAT(J)*TH
          S2 = S2 + T
          DO 450 K = 1,J - 1
              M = M + 1
              HESL(M) = SIN(X(K))* (DFLOAT(N+J+K)*HESD(J)-TH) -
     +                  HESD(J)*COS(X(K))
              HESL(M) = TWO*HESL(M)
  450     CONTINUE
          HESD(J) = DFLOAT(J* (J+2)+N)*HESD(J)**2 +
     +              TH* (TH-DFLOAT(2*J+2)*HESD(J)) +
     +              T* (DFLOAT(J)*TH+HESD(J))
  460 CONTINUE
      DO 470 J = 1,N
          HESD(J) = TWO* (HESD(J)+COS(X(J))*S2)
  470 CONTINUE
      RETURN
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  480 CONTINUE
      DO 490 J = 1,N* (N-1)/2
          HESL(J) = ZERO
  490 CONTINUE
      DO 500 J = 1,N,2
          HESD(J+1) = C200
          HESD(J) = C1200*X(J)**2 - C400*X(J+1) + TWO
          HESL(IX(J+1,J)) = -C400*X(J)
  500 CONTINUE
      RETURN
C
C     EXTENDED POWELL FUNCTION.
C
  510 CONTINUE
      DO 520 J = 1,N* (N-1)/2
          HESL(J) = ZERO
  520 CONTINUE
      DO 530 J = 1,N,4
          T2 = X(J+1) - TWO*X(J+2)
          T3 = X(J) - X(J+3)
          S1 = C12*T2**2
          S2 = C120*T3**2
          HESD(J) = TWO + S2
          HESD(J+1) = C200 + S1
          HESD(J+2) = TEN + FOUR*S1
          HESD(J+3) = TEN + S2
          HESL(IX(J+1,J)) = TWO*TEN
          HESL(IX(J+2,J)) = ZERO
          HESL(IX(J+2,J+1)) = -TWO*S1
          HESL(IX(J+3,J)) = -S2
          HESL(IX(J+3,J+1)) = ZERO
          HESL(IX(J+3,J+2)) = -TEN
  530 CONTINUE
      RETURN
C
C     BEALE FUNCTION.
C
  540 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      HESD(1) = TWO* (S1**2+S2**2+S3**2)
      HESD(2) = TWO*X(1)* (X(1)+TWO*T2+FOUR*X(1)*X(2)**2+SIX*X(2)*T3+
     +          NINE*X(1)*X(2)**4)
      HESL(1) = TWO* (T1-X(1)*S1) + FOUR*X(2)* (T2-X(1)*S2) +
     +          SIX* (T3-X(1)*S3)*X(2)**2
      RETURN
C
C     WOOD FUNCTION.
C
  550 CONTINUE
      HESD(1) = C1200*X(1)**2 - C400*X(2) + TWO
      HESD(2) = C220P2
      HESD(3) = C1080*X(3)**2 - C360*X(4) + TWO
      HESD(4) = C200P2
      HESL(1) = -C400*X(1)
      HESL(2) = ZERO
      HESL(3) = ZERO
      HESL(4) = ZERO
      HESL(5) = C19P8
      HESL(6) = -C360*X(3)
      RETURN
C
C     CHEBYQUAD FUNCTION.
C
  560 CONTINUE
      DO 570 I = 1,N
          FVEC(I) = ZERO
  570 CONTINUE
      DO 590 J = 1,N
          T1 = ONE
          T2 = TWO*X(J) - ONE
          T = TWO*T2
          DO 580 I = 1,N
              FVEC(I) = FVEC(I) + T2
              TH = T*T2 - T1
              T1 = T2
              T2 = TH
  580     CONTINUE
  590 CONTINUE
      D1 = ONE/FLOAT(N)
      IEV = .FALSE.
      DO 600 I = 1,N
          FVEC(I) = D1*FVEC(I)
          IF (IEV) FVEC(I) = FVEC(I) + ONE/ (DFLOAT(I)**2-ONE)
          IEV = .NOT. IEV
  600 CONTINUE
      D2 = TWO*D1
      M = 0
      DO 640 J = 1,N
          HESD(J) = FOUR*D1
          T1 = ONE
          T2 = TWO*X(J) - ONE
          T = TWO*T2
          S1 = ZERO
          S2 = TWO
          P1 = ZERO
          P2 = ZERO
          GVEC(1) = S2
          DO 610 I = 2,N
              TH = FOUR*T2 + T*S2 - S1
              S1 = S2
              S2 = TH
              TH = T*T2 - T1
              T1 = T2
              T2 = TH
              TH = EIGHT*S1 + T*P2 - P1
              P1 = P2
              P2 = TH
              GVEC(I) = S2
              HESD(J) = HESD(J) + FVEC(I)*TH + D1*S2**2
  610     CONTINUE
          HESD(J) = D2*HESD(J)
          DO 630 K = 1,J - 1
              M = M + 1
              HESL(M) = ZERO
              TT1 = ONE
              TT2 = TWO*X(K) - ONE
              TT = TWO*TT2
              SS1 = ZERO
              SS2 = TWO
              DO 620 I = 1,N
                  HESL(M) = HESL(M) + SS2*GVEC(I)
                  TTH = FOUR*TT2 + TT*SS2 - SS1
                  SS1 = SS2
                  SS2 = TTH
                  TTH = TT*TT2 - TT1
                  TT1 = TT2
                  TT2 = TTH
  620         CONTINUE
              HESL(M) = D2*D1*HESL(M)
  630     CONTINUE
  640 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE HESFCN.
C
      END
