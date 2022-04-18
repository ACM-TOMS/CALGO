C======================================================================
C                                                                     c
C     This file contains the FORTRAN program which demonstrates the   c
C use of the suite of routines for generating box constrained         c
C nonlinear optimization  problems.                                   c
C                                                                     c
C     Another purpose of this program is to test the integrity of the c
C routines in the suite after installation. It is by no means a       c
C thorough test.                                                      c
C                                                                     c
C     The other routines that should be linked to this main program   c
C are in files bcpgen.f and bcpex.f                                   c
C                                                                     c
C     The user is referred to the documentation                       c
C                                                                     c
C     Soares, J., Judice, J. and Facchinei, F., "Generating box       c
C     constrained optimization problems",                             c
C     ACM TOMS, (??)                                                  c
C and                                                                 c
C     Soares, J., Judice, J. and Facchinei, F., "FORTRAN subroutines  c
C     for generating box constrained optimization problems",          c
C     ACM TOMS, (??).                                                 c
C                                                                     c
C======================================================================


      PROGRAM VERIFY

C     Authors:          F. Facchinei, J. Soares and J. Judice
C     Date Written:     September 10, 1995
C     Date Modified:    September  3, 1996
C
C  When testing BCP01 (see comments in subroutine BCP01 for details)
C  we make use of the following variables:
C  ===================

C     Input parameters:
C     ( N, *, *, BND, LOWBND, ACTBND, PAIR, DEG, LOWM1, UPPM1,
C       LOWM2, UPPM2, WIDTH1, WIDTH2, HIF, ISEED, * )


C     Output parameters:
C     ( *, L, U, *, *, *, *, *, *, *, *, *, *, *, *, *, INFRM )


C  When testing BCP10 (see comments in subroutine BCP10 for details)
C  we make use of the following variables:
C  ====================

C     Input parameters:
C     ( N, X, *)


C     Output parameters:
C     ( *, *, F)



C  When testing BCP11 (see comments in subroutine BCP11 for details)
C  we make use of the following variables:
C  ====================

C     Input parameters:
C     ( N, X, *)

C     Output parameters:
C     ( *, *, GF)


C  When testing BCP12 (see comments in subroutine BCP12 for details)
C  we make use of the following variables:
C  ====================

C     Input parameters:
C     ( N, X, D, *)


C     Output parameters:
C     ( *, *, *, HFD)
C     .. See the comments in BCP01 for a description of the variables
C     in the common blocks.

C     .. Parameters ..
      INTEGER MAXN
      PARAMETER (MAXN=100)
      INTEGER NN
      PARAMETER (NN=10000)
C     ..
C     .. Scalars in Common ..
      INTEGER ACTL,ACTU,HI,INL,INU,NM1L,NM1U,NM2L,NM2U
C     ..
C     .. Arrays in Common ..
      REAL BETA(NN),XBAR(NN)
      INTEGER PTACT(NN)
C     ..
C     .. Local Scalars ..
      REAL ACUM,F,LOWM1,LOWM2,UPPM1,UPPM2,WIDTH1,WIDTH2
      INTEGER ACTBND,BND,DEG,HIF,I,INFRM,ISEED,LOWBND,N,PAIR
C     ..
C     .. Local Arrays ..
      REAL BETA2(MAXN),D(MAXN),GF(MAXN),HFD(MAXN),L(MAXN),L2(MAXN),
     +     U(MAXN),U2(MAXN),X(MAXN)
      INTEGER PTACT2(MAXN)
C     ..
C     .. External Subroutines ..
      EXTERNAL BCP01,BCP02,BCP10,BCP11,BCP12
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /BCPB1/INL,INU,ACTL,ACTU,NM1L,NM1U,NM2L,NM2U,HI
      COMMON /BCPB2/BETA,XBAR,PTACT
C     ..

C     N = 20 meaning that the problems have 20 variables
C     (N need to be less or equal to NN)

      N = 20

C     Initialize X and D for the calls to BCP10, BCP11 and BCP12

      DO 10 I = 1,N
          X(I) = 1.0 + 1.0/I
          D(I) = 1.0/I
   10 CONTINUE

C  *******************   Generate a problem with 10 lower bounds and
C  VERIFICATION TEST 1   10 upper bounds, distributed randomly.
C  *******************
C                        At the unconstrained solution, 8 bounds should
C                        be active. Among these active bounds, 4 should
C                        have a null multiplier. The remaining 4
C                        multipliers should be uniformly distributed in
C                        [1.0,2.0].
C
C                        Moreover, at the unconstrained solution the
C                        distance between each nonactive bound and the
C                        correspondent component should be
C                        uniformly distributed in [3.0,50.0].
C
C                        Use the simplest function h_i,
C                             h_i(x_i)= beta_i * ( x_i - xbar_i ),
C                        as described in the companion paper, to
C                        define the objective function.

      PRINT *,'**         VERIFICATION TEST 1               **'
      PRINT *

C     BND = 50, meaning that 50% of the all possible bounds (40)
C     are finite.

      BND = 50

C     LOWBND = 50, meaning that lower and upper bounds are
C     equally distributed among the finite bounds.

      LOWBND = 50

C     PAIR = 0, meaning that finite bounds need not exist in pairs.

      PAIR = 0

C     ACTBND = 40, meaning that 40% of the finite bounds are active
C     at the unconstrained solution.

      ACTBND = 40

C     Set the two ranges for the multipliers

      LOWM1 = 0.0
      UPPM1 = 0.0

      LOWM2 = 1.0
      UPPM2 = 2.0

C     DEG = 50, meaning that at the unconstrained solution 50% of the
C     multipliers associated with the active constraints are in the
C     first range while the others are in the second range.

      DEG = 50

C     Set the range for the nonactive but finite bounds.

      WIDTH1 = 3.0
      WIDTH2 = 50.0

C     Select the function h_i

      HIF = 1

C     Set ISEED

      ISEED = 1994

C     Generate the problem

      CALL BCP01(N,L,U,BND,LOWBND,ACTBND,PAIR,DEG,LOWM1,UPPM1,LOWM2,
     +           UPPM2,WIDTH1,WIDTH2,HIF,ISEED,INFRM)

C     Check INFRM

      IF (INFRM.NE.0) THEN
          CALL BCP02(INFRM)
          PRINT *,'**         VERIFICATION TEST 1 INTERRUPTED   **'
          STOP

      END IF

C     The output from BCP01 should have been the following in
C     single precision.

      L2(1) = -0.1000000E+21
      L2(2) = -23.95055
      L2(3) = 1.000000
      L2(4) = -21.96506
      L2(5) = 1.000000
      L2(6) = -0.1000000E+21
      L2(7) = -0.1000000E+21
      L2(8) = -34.71118
      L2(9) = -0.1000000E+21
      L2(10) = -0.1000000E+21
      L2(11) = -21.71301
      L2(12) = -0.1000000E+21
      L2(13) = -0.1000000E+21
      L2(14) = -14.06864
      L2(15) = -15.64741
      L2(16) = 1.000000
      L2(17) = -0.1000000E+21
      L2(18) = -0.1000000E+21
      L2(19) = -0.1000000E+21
      L2(20) = 1.000000

      U2(1) = 1.000000
      U2(2) = 1.000000
      U2(3) = 0.1000000E+21
      U2(4) = 0.1000000E+21
      U2(5) = 11.12918
      U2(6) = 1.000000
      U2(7) = 0.1000000E+21
      U2(8) = 4.876650
      U2(9) = 0.1000000E+21
      U2(10) = 26.85025
      U2(11) = 0.1000000E+21
      U2(12) = 0.1000000E+21
      U2(13) = 39.96765
      U2(14) = 1.000000
      U2(15) = 49.16487
      U2(16) = 0.1000000E+21
      U2(17) = 45.27851
      U2(18) = 0.1000000E+21
      U2(19) = 0.1000000E+21
      U2(20) = 0.1000000E+21

      PTACT2(1) = 3
      PTACT2(2) = 20
      PTACT2(3) = 16
      PTACT2(4) = 5
      PTACT2(5) = 2
      PTACT2(6) = 14
      PTACT2(7) = 1
      PTACT2(8) = 6
      PTACT2(9) = 10
      PTACT2(10) = 5
      PTACT2(11) = 4
      PTACT2(12) = 20
      PTACT2(13) = 12
      PTACT2(14) = 9
      PTACT2(15) = 7
      PTACT2(16) = 19
      PTACT2(17) = 11
      PTACT2(18) = 16
      PTACT2(19) = 3
      PTACT2(20) = 18

      BETA2(1) = 0.0
      BETA2(2) = 0.0
      BETA2(3) = 1.409947
      BETA2(4) = 1.979443
      BETA2(5) = 0.0
      BETA2(6) = 0.0
      BETA2(7) = 1.459439
      BETA2(8) = 1.787866
      BETA2(9) = 0.0
      BETA2(10) = 0.0
      BETA2(11) = 0.0
      BETA2(12) = 0.0
      BETA2(13) = 0.0
      BETA2(14) = 0.0
      BETA2(15) = 0.0
      BETA2(16) = 0.0
      BETA2(17) = 0.0
      BETA2(18) = 0.0
      BETA2(19) = 0.0
      BETA2(20) = 0.0

C     Determine discrepancy

      ACUM = 0.0
      DO 20 I = 1,N
          ACUM = ACUM + ABS(L(I)-L2(I))
          ACUM = ACUM + ABS(U(I)-U2(I))
          ACUM = ACUM + ABS(PTACT(I)-PTACT2(I))
          ACUM = ACUM + ABS(BETA(I)-BETA2(I))
   20 CONTINUE

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in generation = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 1 INTERRUPTED   **'
          STOP

      END IF

C     Evaluate f(X) and compute discrepancy.

      CALL BCP10(N,X,F)

      ACUM = ABS(F-670.363)
      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in    f(X)    = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 1 INTERRUPTED   **'
          STOP

      END IF

C     Evaluate grad f(X) and compute discrepancy.

      CALL BCP11(N,X,GF)

      ACUM = 0.0
      DO 30 I = 1,N
          ACUM = ACUM + GF(I)*GF(I)
   30 CONTINUE
      ACUM = ABS(ACUM-4389455.7521496)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in || grad f(X) ||^2 = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 1 INTERRUPTED   **'
          STOP

      END IF

C     Evaluation of the product Hessian of f at X by the vector D.
C     and compute discrepancy

      CALL BCP12(N,X,D,HFD)

      ACUM = 0.0
      DO 40 I = 1,N
          ACUM = ACUM + HFD(I)*HFD(I)
   40 CONTINUE
      ACUM = ABS(ACUM-15192289.986351)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in || H(X)d ||^2 = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 1 INTERRUPTED   **'
          STOP

      END IF

      PRINT *,'**  VERIFICATION TEST 1 ENDED SUCCESSFULLY   **'
      PRINT *

C  *******************   Generate a problem with 9 lower bounds and
C  VERIFICATION TEST 2   3 upper bounds, distributed ramdomly.
C  *******************
C                        At the unconstrained solution, 6 bounds should
C                        be active. All of them should have a null
C                        multiplier.
C
C                        Moreover, at the unconstrained solution the
C                        distance between each nonactive bound and the
C                        correspondent component should be
C                        uniformly distributed in [100.0,200.0].
C
C                        Use the following function h_i,
C         h_i(x_i)= ( x_i - xbar_i )^3 + beta_i ( x_i - xbar_i )
C                        as described in the companion paper, to
C                        define the objective function.

      PRINT *,'**         VERIFICATION TEST 2               **'
      PRINT *

C     BND = 30, meaning that 30% of the all possible bounds (40)
C     are finite.

      BND = 30

C     LOWBND = 75, meaning that 75% of the finite bounds are lower
C     bounds

      LOWBND = 75

C     PAIR = 0, meaning that finite bounds need not exist in pairs.

      PAIR = 0

C     ACTBND = 50, meaning that 50% of the finite bounds are active
C     at the unconstrained solution.

      ACTBND = 50

C     Set the two ranges for the multipliers

      LOWM1 = 0.0
      UPPM1 = 0.0

      LOWM2 = 0.0
      UPPM2 = 0.0

C     DEG = 100, meaning that at the unconstrained solution all the
C     multipliers associated with the active constraints are in the
C     first range.

      DEG = 100

C     Set the range for the nonactive but finite bounds.

      WIDTH1 = 100.0
      WIDTH2 = 200.0

C     Select the function h_i

      HIF = 2

C     Set ISEED

      ISEED = 1994

C     Generate the problem

      CALL BCP01(N,L,U,BND,LOWBND,ACTBND,PAIR,DEG,LOWM1,UPPM1,LOWM2,
     +           UPPM2,WIDTH1,WIDTH2,HIF,ISEED,INFRM)

C     Check INFRM

      IF (INFRM.NE.0) THEN
          CALL BCP02(INFRM)
          PRINT *,'**         VERIFICATION TEST 2 INTERRUPTED   **'
          STOP

      END IF

C     The output from BCP01 should have been the following in
C     single precision.

      L2(1) = -0.1000000E+21
      L2(2) = -145.7033
      L2(3) = 1.000000
      L2(4) = -141.4788
      L2(5) = 1.000000
      L2(6) = -0.1000000E+21
      L2(7) = -0.1000000E+21
      L2(8) = -0.1000000E+21
      L2(9) = -0.1000000E+21
      L2(10) = -0.1000000E+21
      L2(11) = -140.9426
      L2(12) = -0.1000000E+21
      L2(13) = -0.1000000E+21
      L2(14) = -124.6780
      L2(15) = -128.0370
      L2(16) = 1.000000
      L2(17) = -0.1000000E+21
      L2(18) = -0.1000000E+21
      L2(19) = -0.1000000E+21
      L2(20) = 1.000000

      U2(1) = 0.1000000E+21
      U2(2) = 0.1000000E+21
      U2(3) = 0.1000000E+21
      U2(4) = 0.1000000E+21
      U2(5) = 0.1000000E+21
      U2(6) = 177.5269
      U2(7) = 0.1000000E+21
      U2(8) = 0.1000000E+21
      U2(9) = 0.1000000E+21
      U2(10) = 0.1000000E+21
      U2(11) = 0.1000000E+21
      U2(12) = 0.1000000E+21
      U2(13) = 0.1000000E+21
      U2(14) = 1.000000
      U2(15) = 1.000000
      U2(16) = 0.1000000E+21
      U2(17) = 0.1000000E+21
      U2(18) = 0.1000000E+21
      U2(19) = 0.1000000E+21
      U2(20) = 0.1000000E+21

      PTACT2(1) = 5
      PTACT2(2) = 16
      PTACT2(3) = 3
      PTACT2(4) = 20
      PTACT2(5) = 14
      PTACT2(6) = 15
      PTACT2(7) = 9
      PTACT2(8) = 4
      PTACT2(9) = 16
      PTACT2(10) = 2
      PTACT2(11) = 10
      PTACT2(12) = 3
      PTACT2(13) = 7
      PTACT2(14) = 18
      PTACT2(15) = 13
      PTACT2(16) = 11
      PTACT2(17) = 5
      PTACT2(18) = 19
      PTACT2(19) = 8
      PTACT2(20) = 20

      BETA2(1) = 0.0
      BETA2(2) = 0.0
      BETA2(3) = 0.0
      BETA2(4) = 0.0
      BETA2(5) = 0.0
      BETA2(6) = 0.0
      BETA2(7) = 1.459439
      BETA2(8) = 1.787866
      BETA2(9) = 0.0
      BETA2(10) = 0.0
      BETA2(11) = 0.0
      BETA2(12) = 0.0
      BETA2(13) = 0.0
      BETA2(14) = 0.0
      BETA2(15) = 0.0
      BETA2(16) = 0.0
      BETA2(17) = 0.0
      BETA2(18) = 0.0
      BETA2(19) = 0.0
      BETA2(20) = 0.0

C     Determine discrepancy

      ACUM = 0.0
      DO 50 I = 1,N
          ACUM = ACUM + ABS(L(I)-L2(I))
          ACUM = ACUM + ABS(U(I)-U2(I))
          ACUM = ACUM + ABS(PTACT(I)-PTACT2(I))
          ACUM = ACUM + ABS(BETA(I)-BETA2(I))
   50 CONTINUE

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in generation = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 2 INTERRUPTED   **'
          STOP

      END IF

C     Evaluate f(X) and compute discrepancy.

      CALL BCP10(N,X,F)

      ACUM = ABS(F-669.924)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in    f(X)    = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 2 INTERRUPTED   **'
          STOP

      END IF

C     Evaluate grad f(X) and compute discrepancy.

      CALL BCP11(N,X,GF)

      ACUM = 0.0
      DO 60 I = 1,N
          ACUM = ACUM + GF(I)*GF(I)
   60 CONTINUE

      ACUM = ABS(ACUM-4389190.5930229)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in || grad f(X) ||^2 = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 2 INTERRUPTED   **'
          STOP

      END IF

C     Evaluation of the product Hessian of f at X by the vector D.
C     and compute discrepancy

      CALL BCP12(N,X,D,HFD)

      ACUM = 0.0
      DO 70 I = 1,N
          ACUM = ACUM + HFD(I)*HFD(I)
   70 CONTINUE

      ACUM = ABS(ACUM-15192921.294180)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in || H(X)d ||^2 = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 2 INTERRUPTED   **'
          STOP

      END IF

      PRINT *,'**  VERIFICATION TEST 2 ENDED SUCCESSFULLY   **'
      PRINT *

C  *******************   Generate a problem with 20 lower bounds and
C  VERIFICATION TEST 3   20 upper bounds.
C  *******************
C                        At the unconstrained solution, 20 bounds should
C                        be active. Among these active bounds, 10 should
C                        have a multiplier uniformly distributed in
C                        [0.001,0.1]. The other 10 are uniformly
C                        distributed in [10.0,20.0].
C
C                        Moreover, at the unconstrained solution, the
C                        distance between each nonactive bound and the
C                        correspondent component of such point should be
C                        uniformly distributed in [1.0,2.0].
C
C                        Use the following function h_i,
C      h_i(x_i)= (x_i-xbar_i)^(7/3)+ beta_i(x_i-xbar_i)
C                        as described in the companion paper.

      PRINT *,'**         VERIFICATION TEST 3               **'
      PRINT *

C     BND = 100, meaning that all possible bounds (40) are finite.

      BND = 100

C     LOWBND = 50, meaning that 50% of the finite bounds are lower
C     bounds. (Note: Any other value would result in an error)

      LOWBND = 50

C     PAIR = 0, meaning that finite bounds need not exist in pairs.
C     (Note: In this example its value is indifferent)

      PAIR = 0

C     ACTBND = 50, meaning that 50% of the finite bounds are active
C     at the unconstrained solution.

      ACTBND = 50

C     Set the two ranges for the multipliers

      LOWM1 = 0.001
      UPPM1 = 0.1

      LOWM2 = 10.0
      UPPM2 = 20.0

C     DEG = 50, meaning that at the unconstrained solution 50% of the
C     multipliers associated with the active constraints are in the
C     first range.

      DEG = 50

C     Set the range for the nonactive but finite bounds.

      WIDTH1 = 1.0
      WIDTH2 = 2.0

C     Select the function h_i

      HIF = 3

C     Set ISEED

      ISEED = 1994

C     Generate the problem

      CALL BCP01(N,L,U,BND,LOWBND,ACTBND,PAIR,DEG,LOWM1,UPPM1,LOWM2,
     +           UPPM2,WIDTH1,WIDTH2,HIF,ISEED,INFRM)

C     Check INFRM

      IF (INFRM.NE.0) THEN
          CALL BCP02(INFRM)
          PRINT *,'**         VERIFICATION TEST 3 INTERRUPTED   **'
          STOP

      END IF

C     The output from BCP01 should have been the following in
C     single precision.

      L2(1) = -0.1740453
      L2(2) = 1.000000
      L2(3) = 1.000000
      L2(4) = 1.000000
      L2(5) = 1.000000
      L2(6) = -0.5027176
      L2(7) = -0.4247884
      L2(8) = 1.000000
      L2(9) = -0.4194258
      L2(10) = -0.3779922
      L2(11) = 1.000000
      L2(12) = -0.2903705
      L2(13) = -0.9144108
      L2(14) = 1.000000
      L2(15) = 1.000000
      L2(16) = 1.000000
      L2(17) = -0.2567797
      L2(18) = -0.4670331
      L2(19) = -0.6959825
      L2(20) = 1.000000

      U2(1) = 1.000000
      U2(2) = 2.765269
      U2(3) = 2.115662
      U2(4) = 2.374369
      U2(5) = 2.164240
      U2(6) = 1.000000
      U2(7) = 1.000000
      U2(8) = 2.015358
      U2(9) = 1.000000
      U2(10) = 1.000000
      U2(11) = 2.597139
      U2(12) = 1.000000
      U2(13) = 1.000000
      U2(14) = 2.960955
      U2(15) = 2.390587
      U2(16) = 2.926635
      U2(17) = 1.000000
      U2(18) = 1.000000
      U2(19) = 1.000000
      U2(20) = 2.122221

      PTACT2(1) = 20
      PTACT2(2) = 14
      PTACT2(3) = 15
      PTACT2(4) = 11
      PTACT2(5) = 8
      PTACT2(6) = 4
      PTACT2(7) = 3
      PTACT2(8) = 5
      PTACT2(9) = 2
      PTACT2(10) = 16
      PTACT2(11) = 9
      PTACT2(12) = 7
      PTACT2(13) = 10
      PTACT2(14) = 1
      PTACT2(15) = 19
      PTACT2(16) = 6
      PTACT2(17) = 17
      PTACT2(18) = 12
      PTACT2(19) = 18
      PTACT2(20) = 13

      BETA2(1) = 0.4158475E-01
      BETA2(2) = 0.9796487E-01
      BETA2(3) = 0.5060987E-01
      BETA2(4) = 0.1604595E-01
      BETA2(5) = 0.3231400E-01
      BETA2(6) = 11.04975
      BETA2(7) = 13.14028
      BETA2(8) = 18.67464
      BETA2(9) = 14.59439
      BETA2(10) = 17.87866
      BETA2(11) = 0.9906817E-01
      BETA2(12) = 0.8076739E-01
      BETA2(13) = 0.9251605E-01
      BETA2(14) = 0.4730846E-01
      BETA2(15) = 0.6836787E-01
      BETA2(16) = 18.87231
      BETA2(17) = 16.89056
      BETA2(18) = 19.72215
      BETA2(19) = 10.18641
      BETA2(20) = 12.96863


C     Determine discrepancy

      ACUM = 0.0
      DO 80 I = 1,N
          ACUM = ACUM + ABS(L(I)-L2(I))
          ACUM = ACUM + ABS(U(I)-U2(I))
          ACUM = ACUM + ABS(PTACT(I)-PTACT2(I))
          ACUM = ACUM + ABS(BETA(I)-BETA2(I))
   80 CONTINUE

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in generation = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 3 INTERRUPTED   **'
          STOP

      END IF

C     Evaluate f(X) and compute discrepancy.

      CALL BCP10(N,X,F)

      ACUM = ABS(F-689.545)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in    f(X)    = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 3 INTERRUPTED   **'
          STOP

      END IF

C     Evaluate grad f(X) and compute discrepancy.

      CALL BCP11(N,X,GF)

      ACUM = 0.0
      DO 90 I = 1,N
          ACUM = ACUM + GF(I)*GF(I)
   90 CONTINUE

      ACUM = ABS(ACUM-4384443.0353111)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in || grad f(X) ||^2 = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 3 INTERRUPTED   **'
          STOP

      END IF

C     Evaluation of the product Hessian of f at X by the vector D.
C     and compute discrepancy

      CALL BCP12(N,X,D,HFD)

      ACUM = 0.0
      DO 100 I = 1,N
          ACUM = ACUM + HFD(I)*HFD(I)
  100 CONTINUE

      ACUM = ABS(ACUM-15191153.5)

      IF (ACUM.GT.1.0) THEN
          PRINT *,' Discrepancy in || H(X)d ||^2 = ',ACUM
          PRINT *
          PRINT *,'**         VERIFICATION TEST 3 INTERRUPTED   **'
          STOP

      END IF

      PRINT *,'**  VERIFICATION TEST 3 ENDED SUCCESSFULLY   **'
      PRINT *

      PRINT *,' END OF VERIFICATION '

      END

C     Extended Rosenbrock function. See More, J., Garbow, B. and
C     Hillstrom, K. "Testing Unconstrained Optimization Software",
C     ACM TOMS, Vol. 7, N. 1, March 1981, 17-41.

C******************************************************
C            EXTENDED ROSENBROCK
C******************************************************

C     -- UNCONSTRAINED SOLUTION

      SUBROUTINE SOLUTG(N,XBAR)

C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL XBAR(N)
C     ..
C     .. Local Scalars ..
      INTEGER J
C     ..
      DO 10 J = 1,N
          XBAR(J) = 1.0E+0
   10 CONTINUE

      RETURN

      END

C     -- OBJECTIVE FUNCTION

      SUBROUTINE FUNCTG(N,X,G)

C     .. Scalar Arguments ..
      REAL G
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL X(N)
C     ..
C     .. Local Scalars ..
      REAL T1,T2
      INTEGER J
C     ..
      G = 0.0
      DO 10 J = 1,N,2
          T1 = 1.0E+0 - X(J)
          T2 = 10.0E+0* (X(J+1)-X(J)*X(J))
          G = G + T1*T1 + T2*T2
   10 CONTINUE

      RETURN

      END

C     -- GRADIENT

      SUBROUTINE GRADG(N,X,GG)

C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL GG(N),X(N)
C     ..
C     .. Local Scalars ..
      REAL T1
      INTEGER J
C     ..
      DO 10 J = 1,N,2
          T1 = X(J+1) - X(J)*X(J)
          GG(J) = -400.0E+0*X(J)*T1 - 2.0E+0* (1.0E+0-X(J))
          GG(J+1) = 200.0E+0*T1
   10 CONTINUE

      RETURN

      END

C     -- HESSIAN BY VECTOR PRODUCT

      SUBROUTINE HESGD(N,X,D,HGD)

C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL D(N),HGD(N),X(N)
C     ..
C     .. Local Scalars ..
      REAL T1,T2
      INTEGER J
C     ..
      DO 10 J = 1,N,2
          T1 = -400.0E+0* (X(J+1)-3.0E+0*X(J)*X(J)) + 2.0E+0
          T2 = -400.0E+0*X(J)
          HGD(J) = T1*D(J) + T2*D(J+1)
          HGD(J+1) = T2*D(J) + 200.0E+0*D(J+1)
   10 CONTINUE

      RETURN

      END
