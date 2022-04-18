! Algorithm 421 is coded by FORTRAN 77. Therefore, Algorithm 421 used here was
! modified for accommodation to Fortran 95

! In addition, Algorithm 421 was a little modified for increasing portability of
! the program.
!                                                        M. Kodama
!
! Added a few undeclared integer declarations and a real declaration for error
!                                                        Tim Hopkins
!*************************************************************************
      SUBROUTINE CDLGAM(CARG,CANS,ERROR,LF0)
! COMPLEX GAMMA AND LOGGAMMA FUNCTIONS WITH ERROR ESTIMATE
!
! CARG = A COMPLEX ARGUMENT, GIVEN AS A VECTOR OF 2 DOUBLE
!        PRECISION ELEMENTS CONSISTING OF THE REAL COMPONENT
!        FOLLOWED BY THE IMAGINARY COMPONENT
! CANS = THE COMPLEX ANSWER, OF THE SAME TYPE AS CARG
! ERROR = A REAL VARIABLE.  IT STANDS FOR AN ESTIMATE OF THE
!        ABSOLUTE ERROR OF THE ARGUMENT AS INPUT.  AS OUTPUT
!        IT GIVES AN ESTIMATE OF THE ABSOLUTE (FOR LOGGAMMA)
!        OR THE RELATIVE (FOR GAMMA) ERROR OF THE ANSWER
! LF0  = FLAG.  SET IT TO 0 FOR LOGGAMMA, AND 1 FOR GAMMA
      DOUBLE PRECISION CARG(2),CANS(2),COEF(7),F0,F1,G0,G1,&
         &PI,DPI,HL2P,AL2P,DELTA,DE0,DE1,Z1,Z2,ZZ1,W1,W2,Y1,&
         &A,B,U,U1,U2,UU1,UU2,UUU1,UUU2,V1,V2,VV1,VV2,T1,T2,&
         &H,H1,H2,AL1,AL2,DN,EPS,OMEGA
      REAL ERROR
      INTEGER LF1, LF2, LF3, LF0, K, N, J
      DATA COEF(1)/+0.641025641025641026D-2/
      DATA COEF(2)/-0.191752691752691753D-2/
      DATA COEF(3)/+0.841750841750841751D-3/
      DATA COEF(4)/-0.595238095238095238D-3/
      DATA COEF(5)/+0.793650793650793651D-3/
      DATA COEF(6)/-0.277777777777777778D-2/
      DATA COEF(7)/+0.833333333333333333D-1/
      DATA F0/840.07385296052619D0/,F1/20.001230821894200D0/
      DATA G0/1680.1477059210524D0/,G1/180.01477047052042D0/
      DATA PI/3.14159265358979324D0/
      DATA DPI/6.28318530717958648D0/
      DATA HL2P/0.918938533204672742D0/   !  HL2P=AL2P/2
      DATA AL2P/1.83787706640934548D0/    !  AL2P=LOG(2*PI)
! CONSTANTS EPS AND OMEGA ARE MACHINE DEPENDENT.
! EPS   IS THE BASIC ROUND-OFF UNIT.  FOR S/360 MACHINES,
! IT IS CHOSEN TO BE 16**-13.  FOR BINARY MACHINE OF N-
! BIT ACCURACY SET IT TO 2**(-N+1), AND INITIALIZE DE0
! AS 5.0 RATHER THAN AS 2.0
! OMEGA   IS THE LARGEST NUMBER REPRESENTABLE BY THE FLOAT
! POINT REPRESENTATION OF THE MACHINE.  FOR S/360
!     MACHINES, IT IS SLIGHTLY LESS THAN 16**63.
!      DATA EPS/2.20D-16/
!      DATA OMEGA/7.23700538D75/
! The above two lines were modified to the following two lines 
! for portability by M. Kodama.
      PARAMETER(EPS=EPSILON(1D0))
      PARAMETER(OMEGA=HUGE(1D0))
! The following two lines were added for portability by M. Kodama.
      REAL, PARAMETER:: OMEGA1=HUGE(1.)
      INTEGER, PARAMETER:: IRADIX=RADIX(1D0)
      Z1 = CARG(1)
      Z2 = CARG(2)
      DELTA = ABS(ERROR)
      DE0 = 2.0D0
      IF(IRADIX == 2) DE0 = 5.0D0
! The line above is added according to the comment sentence 19 lines above
! for portability by M. Kodama.
      DE1 = 0.0
! FORCE SIGN OF IMAGINARY PART OF ARG TO NON-NEGATIVE
      LF1 = 0
      IF (Z2 .GE. 0.0) GO TO 20
      LF1 = 1
      Z2 = -Z2
   20 LF2 = 0
      IF (Z1 .GE. 0.0) GO TO 100
! CASE WHEN REAL PART OF ARG IS NEGATIVE
      LF2 = 1
      LF1 = LF1-1
      T1 = AL2P - PI*Z2
      T2 = PI*(0.5D0 - Z1)
      U  = -DPI*Z2
      IF (U .GE. -0.1054D0) GO TO 40
      A  = 0.0D0
! IF E**U .LT. 10**(-17), IGNOR IT TO SAVE TIME AND TO AVOID
! IRRELEVANT UNDERFLOW
      IF (U .LE. -39.15D0) GO TO 30
      A = DEXP(U)
   30 H1 = 1.0D0 - A
      GO TO 50
   40 U2 = U*U
      A = -U*(F1*U2 + F0)
      H1 = (A + A)/((U2 + G1)*U2 + G0 + A)
      A = 1.0D0 - H1
! DINT IS THE DOUBLE PRECISION VERSION OF AINT, INTEGER EX-
!   TRACTION.  THIS FUNCTION IS NOT INCLUDED IN ANSI FORTRAN
!   .  WHEN THIS FUNCTION IS NOT PROVIDED BY THE SYSTEM,
!   EITHER SUPPLY IT AS AN EXTERNAL SUBROUTINE (AND TYPE THE
!   NAME DINT AS DOUBLE PRECISION), OR MODIFY THE NEXT
!   STATEMENT AS THE EXAMPLE FOR S/340 INDICATES.  FOR S/360
!   REPLACE IT WITH
!
!      DOUBLE PRECISION SCALE
!      DATA SCALE/Z4F00000000000000/
!  50 B = Z1 - ((Z1 - 0.5D0) + SCALE)
   50 B = Z1 - DINT(Z1 - 0.5D0)
      H2 = A*DSIN(DPI*B)
      B  = DSIN(PI*B)
      H1 = H1 + (B+B)*B*A
      H = DABS(H2) + H1 - DPI*A*DELTA
      IF (H .LE. 0.0) GO TO 500
      DE0 = DE0 + DABS(T1) + T2
      DE1 = PI + DPI*A/H
      Z1 = 1.0D0 - Z1
! CASE WHEN NEITHER REAL PART NOR IMAGINARY PART OF ARG IS
! NEGATIVE.  DEFINE THERSHOLD CURVE TO BE THE BROKEN LINES
! CONNECTING POINTS 10F0*I, 10F4.142*I, 0.1F14.042*I,AND
! 0.1FOMEGA*I
  100 LF3 = 0
      Y1 = Z1 - 0.5D0
      W1 = 0.0
      W2 = 0.0
      K  = 0
      B  = DMAX1(0.1D0, DMIN1(10.0D0, 14.142D0-Z2)) - Z1
      IF (B .LE. 0.0) GO TO 200
! CASE WHEN REAL PART OF ARG IS BETWEEN 0 AND THRESHOLD
      LF3 = 1
      ZZ1 = Z1
      N  = B + 1.0D0
      DN = N
      Z1 = Z1 + DN
      A  = Z1*Z1 + Z2*Z2
      V1 = Z1/A
      V2 = -Z2/A
! INITIALIZE U1+U2*I AS THE RIGHTMOST FACTOR 1-1/(Z+M)
      U1 = 1.0D0 - V1
      U2 = -V2
      K  = 6.0D0 - Z2*0.6D0 - ZZ1
      IF (K .LE. 0) GO TO 120
! FORWORD ASSEMBLY OF FACTORS (Z+J-1)/(Z+M)
      N  = N - K
      UU1 = (ZZ1*Z1 + Z2*Z2) / A
      UU2 = DN*Z2/A
      VV1 = 0.0
      VV2 = 0.0
      DO 110 J = 1,K
        B  = U1*(UU1+VV1) - U2*(UU2+VV2)
        U2 = U1*(UU2+VV2) + U2*(UU1+VV1)
        U1 = B
        VV1 = VV1 + V1
        VV2 = VV2 + V2
110   CONTINUE
  120 IF (N .LE. 1) GO TO 140
! BACKWARD ASSEMBLY OF FACTORS 1-J/(Z+N)
      VV1 = V1
      VV2 = V2
      DO 130  J = 2,N
        VV1 = VV1 + V1
        VV2 = VV2 + V2
        B  = U1*(1.0D0 - VV1) + U2*VV2
        U2 = -U1*VV2 + U2*(1.0D0 - VV1)
        U1 = B
130   CONTINUE
  140 U  = U1*U1 + U2*U2
      IF (U .EQ. 0.0) GO TO 500
      IF (LF0 .EQ. 0) GO TO 150
      IF (K .LE. 0) GO TO 200
  150 AL1 = DLOG(U)*0.5D0
      IF (LF0 .NE. 0) GO TO 160
      W1 = AL1
      W2 = DATAN2(U2,U1)
      IF (W2 .LT. 0.0) W2 = W2 + DPI
      IF (K .LE. 0) GO TO 200
  160 A = ZZ1 + Z2 - DELTA
      IF (A .LT. 0.0) GO TO 500
      DE0 = DE0 - AL1
      DE1 = DE1 + 2.0D0 + 1.0D0/A
! CASE WHEN REAL PART OF ARG IS GREATER THAN THRESHOLD
  200 A = Z1*Z1 + Z2*Z2
      AL1 = DLOG(A)*0.5D0
      AL2 = DATAN2(Z2,Z1)
      V1 = Y1*AL1 - Z2*AL2
      V2 = Y1*AL2 + Z2*AL1
! EVALUATE ASYMTOTIC TERMS.  IGNORE THIS TERM,IF ABS VAL(ARG) .GT.
! 10**9, TO SAVE TIME AND TO AVOID IRRELEVANT UNDERFLOW
      VV1 = 0.0
      VV2 = 0.0
      IF (A .GT. 1.0D18) GO TO 220
      UU1 = Z1/A
      UU2 = -Z2/A
      UUU1 = UU1*UU1 - UU2*UU2
      UUU2 = UU1*UU2*2.0D0
      VV1 = COEF(1)
      DO 210  J = 2,7
        B  = VV1*UUU1 - VV2*UUU2
        VV2 = VV1*UUU2 + VV2*UUU1
        VV1 = B + COEF(J)
210   CONTINUE
      B  = VV1*UU1 -VV2*UU2
      VV2 = VV1*UU2 + VV2*UU1
      VV1 = B
  220 W1 = (((VV1 + HL2P) - W1) - Z1) + V1
      W2 = ((VV2 - W2) -Z2) + V2
      DE0 = DE0 + DABS(V1) + DABS(V2)
      IF (K .LE. 0) DE1 = DE1 + AL1
! FINAL ASSEMBLY
      IF (LF2 .NE. 0) GO TO 310
      IF (LF0 .EQ. 0) GO TO 400
      A = DEXP(W1)
      W1 = A*DCOS(W2)
      W2 = A*DSIN(W2)
      IF (LF3 .EQ. 0) GO TO 400
      B  = (W1*U1 + W2*U2) / U
      W2 = (W2*U1 - W1*U2) / U
      W1 = B
      GO TO 400
  310 H = H1*H1 + H2*H2
      IF (H .EQ. 0.0) GO TO 500
      IF (LF0 .EQ. 0) GO TO 320
      IF (H .GT. 1.0D-2) GO TO 330
  320 A = DLOG(H)*0.5D0
      IF (H .LE. 1.0D-2) DE0 = DE0 - A
      IF (LF0 .NE. 0) GO TO 330
      W1 = (T1 - A) - W1
      W2 = (T2 - DATAN2(H2,H1)) - W2
      GO TO 400
  330 T1 = T1 - W1
      T2 = T2 - W2
      A  = DEXP(T1)
      T1 = A*DCOS(T2)
      T2 = A*DSIN(T2)
      W1 = (T1*H1 + T2*H2)/H
      W2 = (T2*H1 - T1*H2)/H
      IF (LF3 .EQ. 0) GO TO 400
      B  = W1*U1 - W2*U2
      W2 = W1*U2 + W2*U1
      W1 = B
  400 IF (LF1 .NE. 0) W2 = -W2
! TRUNCATION ERROR OF STIRLINGS FORMULA IS UP TO 3*10**-17.
      DE1 = DE0*EPS + 3.0D-17 + DE1*DELTA
      GO TO 600
! CASE WHEN ARGUMENT IS TOO CLOSE TO A SINGURARITY
!
  500 W1 = OMEGA
      W2 = OMEGA
!      DE1 = OMEGA
! The above line was replaced to the following line to avoid overflows by
! M. Kodama.
      DE1 = OMEGA1
!
  600 CANS(1) = W1
      CANS(2) = W2
      ERROR = DE1
      RETURN
      END
! The end of Algorithm 421

