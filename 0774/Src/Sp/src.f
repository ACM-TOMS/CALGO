C======================================================================
C                                                                     c
C     This file contains the FORTRAN routines that the user should    c
C link to his/her code in order to employ a box constrained           c
C generated problem.                                                  c
C     Other routines that should be linked, besides the main program, c
C include subroutines to evaluate the unconstrained function          c
C as well as its derivatives.                                         c
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
C     To change the precision of this set of subroutines the user     c
C should do the following:                                            c
C                                                                     c
C     1. From double to real:                                         c
C         1.1. Replace all 'DOUBLE PRECISION' with 'REAL'             c
C         1.2.       >>      'D+'    >>   'E+'                        c
C                                                                     c
C     2. From real to double:                                         c
C         1.1. Replace all 'REAL' with 'DOUBLE PRECISION'             c
C         1.2.       >>     'E+'    >>   'D+'                         c
C                                                                     c
C     IMPORTANT EXCEPTION: Do not replace lower case 'real'. In case  c
C     of doubt make sure that RAND in its definition and calls is     c
C     always single precision.                                        c
C                                                                     c
C======================================================================

C    ..
C    .. Begin of subroutine BCP01 ..

      SUBROUTINE BCP01(N,L,U,BND,LOWBND,ACTBND,PAIR,DEG,LOWM1,UPPM1,
     +                 LOWM2,UPPM2,WIDTH1,WIDTH2,HIF,ISEED,INFRM)

C     _Authors:          F. Facchinei, J. Soares and J. Judice
C     _Date Written:     September 10, 1995
C     _Date Modified:    September  3, 1996
C     _Purpose:
C
C        BCP01 generates box constrained nonlinear optimization
C     problems of the form
C
C                      Minimize      f(x)
C                          s.t.   l <= x <= u
C
C     where x belongs to R^n and f is a continuously differentiable
C     function R^n -> R.
C
C     More precisely, BCP01 sets up the bounds, and the data required
C     to evaluate the objective function f and its derivatives.
C     Unsuccessful generation is detected by having a nonzero value
C     in INFRM as output.
C
C     _Arguments in parameter list:
C
C  N       (input) INTEGER
C          The number of variables in the problem. (Must be <= NN,
C          see common block below)
C
C  L, U    (output) DOUBLE PRECISION vectors of dimension N
C          The bounds of the new problem.
C
C  BND     (input) INTEGER
C          The number of bound constraints (or simply, bounds),
C          expressed as the percentage of the 2*n possible constraints,
C          constraints,that are really existing (Must be in [0,100]).
C          Example: if N is 500 and BND=40 then the problem has 400
C          bounds.
C
C  LOWBND  (input) INTEGER
C          The number of lower bounds, expressed as the percentage of
C          the number of constraints as determined by BND (Must be in
C          [0,100]). This quantity also determines the number of upper
C          bounds.
C          Example: if N is 500, BND is 40 and LOWBND=75 then there
C          are 300 (100) lower (upper) bounds.
C
C  ACTBND  (input) INTEGER
C          The number of active bounds at the unconstrained solution,
C          expressed as the percentage of the number of constraints as
C          determined by BND (Must be in [0,100]). These active
C          bounds are distributed proportionally to LOWBND.
C          Example: Continuing the previous example, if ACTBND
C          is 50 then 150 (50) of the lower (upper) bounds should
C          be active at the unconstrained solution.
C
C  PAIR    (input) INTEGER
C          If this parameter is set to 1 then the constraints as
C          determined by BND must come in pairs, i.e. each variable
C          has either no constraint or both lower and upper bounds (in
C          which case LOWBND must be set to 50).  If PAIR is
C          set to 0 then the indices of variables having upper bounds
C          are chosen randomly and independently of the lower bounds
C          (Must be 0 or 1).
C
C  DEG   (input) INTEGER
C          The number of Lagrange multipliers (at the unconstrained
C          solution) uniformly distributed in the interval
C          [LOWM1, UPPM1] (see below), expressed as the percentage of
C          the active constraints (Must be in [0,100]).
C          The remaining active constraints have Lagrange multipliers
C          uniformly distributed in [LOWM2, UPPM2] (see below).
C          Example:
C          Continuing the example, if  DEG=10 then 15 (5) active
C          lower (upper) bounds (at the unconstrained solution)  have
C          Lagrange multipliers uniformly distributed in [LOWM1,
C          UPPM1], while the remaining active bounds have Lagrange
C          multipliers uniformly distributed in [LOWM2, UPPM2].
C
C  LOWM1, UPPM1, LOWM2, UPPM2   (input) DOUBLE PRECISION
C          These parameters specify the ranges for the
C          Lagrange multipliers according to what explained above.
C          (Must be >= 0)
C
C  WIDTH1, WIDTH2   (input) DOUBLE PRECISION
C          These two parameters control the value of the bounds that
C          are not active at the unconstrained solution. More
C          precisely, if the i-th component has a lower bound not
C          active at the unconstrained solution xbar then this bound
C          is randomly chosen in the interval
C          [xbar_i - WIDTH2, xbar_i - WIDTH1]. Analogously,
C          if the i-th component has an upper bound not active at
C          the unconstrained solution xbar then this bound is
C          randomly chosen in the interval
C          [xbar_i + WIDTH1, xbar_i + WIDTH2] (Must be >= 0).
C
C  HIF    (input) INTEGER
C          This parameter determines which function h_i is to be used
C          by the generator. If HIF is equal to
C
C          1: h_i(x_i)= beta_i * ( x_i - xbar_i );
C
C          2: h_i(x_i)=( x_i - xbar_i )^3 + beta_i*( x_i - xbar_i );
C
C          3: h_i(x_i)=(x_i - xbar_i)^(7/3) + beta_i*(x_i - xbar_i).
C
C          (Must be  1, 2 or 3)
C
C  ISEED   (input) INTEGER
C          A seed for the random number generator. By running the box
C          constrained problems generator with the same choice of
C          parameters but two different seeds, two different box
C          constrained problems are generated with the same overall
C          characteristics (Must be >= 0).
C
C  INFRM   (output) INTEGER
C          Determines if a successful call to BCP01 has occurred
C          (INFRM=0) or not (INFRM<>0).
C          If INFRM is equal to
C            1: N is less than 1 or greater than NN
C            2, 3, 4, 5, 6, 11: Some integer parameter is defined
C               out of its range.
C            7: PAIR=1 but LOWBND is not 50.
C            8: Either LOWM1 or LOWM2 is less than zero.
C            9: Either LOWM1 is greater than UPPM1 or LOWM2 is greater
C               than UPPM2.
C            10: Either WIDTH1 is less than zero or WIDTH1 is greater
C               than WIDTH2.
C            12: It is impossible to make a selection because LOWBND
C               is too large.
C            13: It is impossible to make a selection because LOWBND
C               is too small.
C            14: It is impossible to make a selection because ACTBND
C               is too large.
C
C     _Explanation of the variables used in the common blocks:
C
C     NN     INTEGER
C            Maximum N allowed.
C
C     INL    INTEGER
C            Number of lower bounds.
C
C     INU    INTEGER
C            Number of upper bounds.
C
C     ACTL   INTEGER
C            Number of lower bounds active at the unconstrained solution
C
C     ACTU   INTEGER
C            Number of upper bounds active at the unconstrained solution
C
C     NM1L   INTEGER
C            Number of lower bounds active at the unconstrained solution
C            which have a Lagrange multiplier randomly chosen in
C            [LOWM1,UPPM1].
C
C     NM1U   INTEGER
C            Number of upper bounds active at the unconstrained solution
C            which have a Lagrange multiplier randomly chosen in
C            [LOWM1,UPPM1].
C
C     NM2L   INTEGER
C            Number of lower bounds active at the unconstrained solution
C            which have a Lagrange multiplier randomly chosen in
C            [LOWM2,UPPM2].
C
C     NM2U   INTEGER
C            Number of upper bounds active at the unconstrained solution
C            which have a Lagrange multiplier randomly chosen in
C            [LOWM2,UPPM2].
C
C     HI     INTEGER
C            Characterizes the function h_i used. HI=HIF.
C
C     PTACT  INTEGER vector of dimension NN
C            A vector of pointers that characterizes the bounds that are
C            active at the unconstrained solution, i.e.,
C            PTACT(II) (for II=1,ACTL) pinpoints an active lower bound
C            PTACT(II) (for II=ACTL+1,ACTL+ACTU) pinpoints an active
C            upper bound.
C
C     BETA   DOUBLE PRECISION vector of dimension NN
C            At the unconstrained solution, BETA characterizes the
C            multipliers of the bounds that are active, i.e.,
C            BETA(II) (for II=1,ACTL) is the multiplier associated
C            with the lower bound in the PTBETA(II) variable;
C            BETA(II) (for II=ACTL+1,ACTL+ACTU) is the multiplier
C            associated with the upper bound in the PTBETA(II) variable.
C
C     XBAR   DOUBLE PRECISION vector of dimension NN
C            The unconstrained solution.
C     .. Parameters ..
      INTEGER NN
      PARAMETER (NN=10000)
      REAL BIGINF
      PARAMETER (BIGINF=1.0E+20)
C     ..
C     .. Scalar Arguments ..
      REAL LOWM1,LOWM2,UPPM1,UPPM2,WIDTH1,WIDTH2
      INTEGER ACTBND,BND,DEG,HIF,INFRM,ISEED,LOWBND,N,PAIR
C     ..
C     .. Array Arguments ..
      REAL L(N),U(N)
C     ..
C     .. Scalars in Common ..
      INTEGER ACTL,ACTU,HI,INL,INU,NM1L,NM1U,NM2L,NM2U
C     ..
C     .. Arrays in Common ..
      REAL BETA(NN),XBAR(NN)
      INTEGER PTACT(NN)
C     ..
C     .. Local Scalars ..
      REAL AUX
      INTEGER I,II,INBND,J,J1,J2,JJ
C     ..
C     .. External Functions ..
      REAL RAND
      EXTERNAL RAND
C     ..
C     .. External Subroutines ..
      EXTERNAL BCP04,SOLUTG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /BCPB1/INL,INU,ACTL,ACTU,NM1L,NM1U,NM2L,NM2U,HI
      COMMON /BCPB2/BETA,XBAR,PTACT
C     ..

C     Test the input parameters.

      INFRM = 0

      IF ((N.GT.NN) .OR. (N.LE.0)) INFRM = 1
      IF ((BND.LT.0) .OR. (BND.GT.100)) INFRM = 2
      IF ((LOWBND.LT.0) .OR. (LOWBND.GT.100)) INFRM = 3
      IF ((ACTBND.LT.0) .OR. (ACTBND.GT.100)) INFRM = 4
      IF ((DEG.LT.0) .OR. (DEG.GT.100)) INFRM = 5
      IF ((PAIR.NE.0) .AND. (PAIR.NE.1)) INFRM = 6
      IF ((PAIR.EQ.1) .AND. (LOWBND.NE.50)) INFRM = 7
      IF ((LOWM1.LT.0) .OR. (LOWM2.LT.0)) INFRM = 8
      IF ((LOWM1.GT.UPPM1) .OR. (LOWM2.GT.UPPM2)) INFRM = 9
      IF ((WIDTH1.LE.0) .OR. (WIDTH1.GT.WIDTH2)) INFRM = 10
      IF ((HIF.LE.0) .OR. (HIF.GE.4)) INFRM = 11

      IF (INFRM.NE.0) RETURN

C     Determine the number of constraints

      INBND = (BND*2*N)/100

C     Determine the number of lower bounds

      INL = (LOWBND*INBND)/100

C     Safeguard the maximum number of lower bounds

      IF (INL.GT.N) THEN
          INFRM = 12
          RETURN

      END IF

C     Determine the number of upper bounds

      INU = INBND - INL

C     Safeguard the maximum number of upper bounds

      IF (INU.GT.N) THEN
          INFRM = 13
          RETURN

      END IF

C     Determine the number of active constraints at XBAR

      I = (ACTBND*INBND)/100
      ACTL = (ACTBND*INL)/100
      ACTU = I - ACTL

C     Safeguard some settings related to the parameter PAIR ..

      IF (PAIR.EQ.0) THEN

C         .. Since the bounds are supposed to be chosen randomly,
C        safeguard the total number of constraints
C        active at XBAR.

          IF (ACTL+ACTU.GT.N) THEN
              INFRM = 14
              RETURN

          END IF

      ELSE

C        .. Since the bounds are supposed to come in pairs,
C        safeguard that the number of lower bounds must
C        equal the number of upper bounds ..

          J = MIN(INL,INU)
          INL = J
          INU = J

C        .. in which case the number of bounds and active bounds at XBAR
C        should be corrected ..

          INBND = INL + INU
          ACTL = (ACTBND*INL)/100
          ACTU = ACTL

C        .. there should not exist one pair of bounds that are both
C        active at XBAR (Note: INL=INU).

          IF (ACTL+ACTU.GT.INL) THEN
              INFRM = 14
              RETURN

          END IF

      END IF

C     Read the optimal unconstrained point

      CALL SOLUTG(N,XBAR)

C     Set default values to nonexistent bounds and make some random
C     settings. For the moment, PTACT is used as an auxiliary vector.
C     After a call to BCP04 it will contain the set of integer values
C     1..N in a random order.

      DO 10 I = 1,N
          L(I) = -BIGINF
          U(I) = BIGINF
          PTACT(I) = I
   10 CONTINUE

      CALL BCP04(N,PTACT(1),ISEED)

C     Defining the lower bounds ..

C     .. actives at XBAR (Indices are in PTACT(1)..PTACT(ACTL)) ..

      J = 1
      DO 20 II = 1,ACTL
          I = PTACT(J)
          L(I) = XBAR(I)
          J = J + 1
   20 CONTINUE

C     .. and non actives at XBAR
C     (Indices are in PTACT(ACTL+1)..PTACT(INL)) ..

      DO 30 II = ACTL + 1,INL
          I = PTACT(J)
          AUX = WIDTH1 + RAND(ISEED)* (WIDTH2-WIDTH1)
          L(I) = XBAR(I) - AUX
          J = J + 1
   30 CONTINUE

C     .. (Note that PTACT(INL+1)..PTACT(N) contain indices of variables
C     that do not have lower bounds).

C     Defining the upper bounds ..

      IF (PAIR.EQ.0) THEN

C        .. since the bounds are supposed to be chosen randomly,
C        disorder PTACT(ACTL+1)..PTACT(N) ..

          J = N - ACTL
          CALL BCP04(J,PTACT(ACTL+1),ISEED)

C        .. and determine the upper bounds that are active at XBAR
C        (Indices are in PTACT(ACTL+1)..PTACT(ACTL+ACTU) ) ..

          J = 1
          JJ = ACTL + ACTU
          DO 40 II = 1,ACTU
              I = PTACT(JJ)
              U(I) = XBAR(I)
              PTACT(JJ) = PTACT(J)
              PTACT(J) = I
              J = J + 1
              JJ = JJ - 1
   40     CONTINUE

C        .. (PTACT(1)..PTACT(ACTU) contain the indices of upper bounds
C        that are active at XBAR)
C        Now, disorder PTACT(ACTU+1)..PTACT(N) ..

          J = N - ACTU
          CALL BCP04(J,PTACT(ACTU+1),ISEED)
C
C        .. to define the upper bounds that are not active at XBAR
C        (Indices are in PTACT(ACTU+1)..PTACT(INU)) ..

          J = ACTU + 1
          DO 50 II = ACTU + 1,INU
              I = PTACT(J)
              AUX = WIDTH1 + RAND(ISEED)* (WIDTH2-WIDTH1)
              U(I) = XBAR(I) + AUX
              J = J + 1
   50     CONTINUE

      ELSE

C        .. since the bounds are supposed to come in pairs,
C        disorder PTACT(ACTL+1)..PTACT(INL) ..

          J = INL - ACTL
          CALL BCP04(J,PTACT(ACTL+1),ISEED)

C        .. and determine the upper bounds that are actives at XBAR
C       (Indices are in PTACT(ACTL+1)..PTACT(ACTL+ACTU). Also note that
C       ACTL+ACTU.LE.INL=INU) ..

          J = 1
          JJ = ACTL + ACTU
          DO 60 II = 1,ACTU
              I = PTACT(JJ)
              U(I) = XBAR(I)
              PTACT(JJ) = PTACT(J)
              PTACT(J) = I
              J = J + 1
              JJ = JJ - 1
   60     CONTINUE

C        .. (PTACT(1)..PTACT(ACTU) contain the indices of upper bounds
C        that are active at XBAR).
C        Now, disorder PTACT(ACTU+1)..PTACT(INU) ..

          J = INU - ACTU
          CALL BCP04(J,PTACT(ACTU+1),ISEED)

C        .. to define the bounds that are not active at XBAR
C        (Indices are in PTACT(ACTU+1)..PTACT(INU)).

          J = ACTU + 1
          DO 70 II = ACTU + 1,INU
              I = PTACT(J)
              AUX = WIDTH1 + RAND(ISEED)* (WIDTH2-WIDTH1)
              U(I) = XBAR(I) + AUX
              J = J + 1
   70     CONTINUE
      END IF

C     Below PTACT is finally defined for which it was intended for.
C     That is ..

C     .. PTACT(1)..PTACT(ACTL) contains the indices of lower bounds
C     that are active at XBAR

C     .. PTACT(ACTL+1)..PTACT(ACTL+ACTU) contains the indices of upper
C     bounds that are active at XBAR

      J1 = 1
      J2 = ACTL + 1
      DO 80 I = 1,N
          IF (L(I).EQ.XBAR(I)) THEN
              PTACT(J1) = I
              J1 = J1 + 1

          ELSE IF (U(I).EQ.XBAR(I)) THEN
              PTACT(J2) = I
              J2 = J2 + 1
          END IF

   80 CONTINUE

C     Define the Lagrange Multipliers at XBAR ..

C     .. How many active constraints will have Lagrange
C     Multipliers in the interval [LOWM1,UPPM1] ..

      NM1L = ACTL*DEG/100
      NM1U = ACTU*DEG/100

C     .. How many active constraints will have Lagrange
C     Multipliers in the interval [LOWM2,UPPM2] ..

      NM2L = ACTL - NM1L
      NM2U = ACTU - NM1U

C     .. Now, disorder PTACT(1)..PTACT(ACTL) ..

      J = ACTL
      CALL BCP04(J,PTACT(1),ISEED)

C     .. to determine randomly the values of the multipliers at XBAR
C     for the active lower bounds ..

      DO 90 J = 1,NM1L
          BETA(J) = LOWM1 + RAND(ISEED)* (UPPM1-LOWM1)
   90 CONTINUE

      DO 100 J = NM1L + 1,NM1L + NM2L
          BETA(J) = LOWM2 + RAND(ISEED)* (UPPM2-LOWM2)
  100 CONTINUE

C     .. Now, disorder PTACT(ACTL+1)..PTACT(ACTL+ACTU) ..

      J = ACTU
      CALL BCP04(J,PTACT(ACTL+1),ISEED)

C     .. to determine randomly the values of the multipliers at XBAR
C     for the active upper bounds ..

      DO 110 J = ACTL + 1,ACTL + NM1U
          BETA(J) = LOWM1 + RAND(ISEED)* (UPPM1-LOWM1)
  110 CONTINUE

      DO 120 J = ACTL + NM1U + 1,ACTL + NM1U + NM2U
          BETA(J) = LOWM2 + RAND(ISEED)* (UPPM2-LOWM2)
  120 CONTINUE

C     Copy HIF into HI

      HI = HIF

      RETURN

      END

C     ..
C     .. End of subroutine BCP01 ..

C     ..
C     .. Begin of subroutine BCP02 ..

      SUBROUTINE BCP02(INFRM)

C     _Authors:         F. Facchinei, J. Soares and J. Judice
C     _Date Written:    September 10, 1995
C     _Date Modified:   September  3, 1996
C     _Purpose:
C
C       BCP02 prints an adequate error message on the screen according
C     to the integer value in INFRM.
C
C     .. Scalar Arguments ..

      INTEGER INFRM
C     ..

      IF (INFRM.EQ.0) THEN
          PRINT *,' **     INFRM = 0 **     EVERYTHING IS FINE      **'

      ELSE IF (INFRM.EQ.1) THEN
          PRINT *,' **     INFRM = 1 **       N<=0  OR N > NN       **'

      ELSE IF (INFRM.EQ.2) THEN
          PRINT *,' **     INFRM = 2 **     BND<0 OR BND>100        **'

      ELSE IF (INFRM.EQ.3) THEN
          PRINT *,' **     INFRM = 3 **   LOWBND<0 OR LOWBND>100    **'

      ELSE IF (INFRM.EQ.4) THEN
          PRINT *,' **     INFRM = 4 **   ACTBND<0 OR ACTBND>100    **'

      ELSE IF (INFRM.EQ.5) THEN
          PRINT *,' **     INFRM = 5 **      DEG<0 OR DEG>100       **'

      ELSE IF (INFRM.EQ.6) THEN
          PRINT *,' **     INFRM = 6 **       PAIR<> 0 AND 1        **'

      ELSE IF (INFRM.EQ.7) THEN
          PRINT *,' **     INFRM = 7 **  PAIR=1 AND LOWBND<>50      **'

      ELSE IF (INFRM.EQ.8) THEN
          PRINT *,' **     INFRM = 8 **       LOWM1<0 OR LOWM2<0    **'

      ELSE IF (INFRM.EQ.9) THEN
          PRINT *,' **     INFRM = 9 ** LOWM1>UPPM1 OR LOWM2>UPPM2  **'

      ELSE IF (INFRM.EQ.10) THEN
          PRINT *,' **     INFRM =10 **WIDTH1.LE.0 OR WIDTH1>WIDTH2 **'

      ELSE IF (INFRM.EQ.11) THEN
          PRINT *,' **     INFRM =11 **    HIF<1    OR HIF>3        **'

      ELSE IF (INFRM.EQ.12) THEN
          PRINT *,' **     INFRM =12 **     LOWBND TOO LARGE        **'

      ELSE IF (INFRM.EQ.13) THEN
          PRINT *,' **     INFRM =13 **     LOWBND TOO LITTLE       **'

      ELSE IF (INFRM.EQ.14) THEN
          PRINT *,' **     INFRM =14 **     ACTBND TOO LARGE        **'

      ELSE
          PRINT *,' **           INFRM OUT OF RANGE                 **'

      END IF

      END

C     ..
C     .. End of subroutine BCP02 ..

C     ..
C     .. Begin of subroutine BCP03 ..

      SUBROUTINE BCP03(N,L,U)

C     _Authors:         F. Facchinei, J. Soares and J. Judice
C     _Date Written:    September 10, 1995
C     _Date Modified:   September  3, 1996
C     _Purpose:
C
C       BCP03 prints on the screen all the data that is in the
C     common blocks. Such information was obtained by a
C     previous call to BCP01.
C
C     .. Scalar Arguments ..

      INTEGER N
C     ..
C     .. Array Arguments ..

      REAL L(N),U(N)
C     ..
C     .. Parameters ..

C     .. See the comments in BCP01 for a description of the variables
C     in the common blocks.
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
      INTEGER I
C     ..
C     .. Common blocks ..
      COMMON /BCPB1/INL,INU,ACTL,ACTU,NM1L,NM1U,NM2L,NM2U,HI
      COMMON /BCPB2/BETA,XBAR,PTACT
C     ..

      WRITE (*,FMT=9000)
      DO 10 I = 1,N
          WRITE (*,FMT=9010) I,L(I),XBAR(I),U(I),PTACT(I),BETA(I)
   10 CONTINUE

      PRINT *
      PRINT *,' #Lower Bounds = ',INL
      PRINT *,' #Upper Bounds = ',INU
      PRINT *,' #Active Lower Bounds = ',ACTL
      PRINT *,' #Active Upper Bounds = ',ACTU
      PRINT *,' #Active Lower Bounds with mult in the',
     +  ' first interval = ',NM1L
      PRINT *,' #Active Upper Bounds with mult in the',
     +  ' first interval = ',NM1U
      PRINT *,' #Active Lower Bounds with mult in the',
     +  ' second interval = ',NM2L
      PRINT *,' #Active Upper Bounds with mult in the',
     +  ' second interval = ',NM2U
      PRINT *,' Function h_i: ',HI


      RETURN

C     .. Formats ..

 9000 FORMAT (9X,'L',14X,'XBAR',12X,'U',8X,'PTACT',7X,'BETA')
 9010 FORMAT (I3,1X,G14.7,1X,G14.7,1X,G14.7,1X,I3,3X,G14.7)
      END

C     ..
C     .. End of subroutine BCP03 ..

C     ..
C     .. Begin of subroutine BCP04 ..

      SUBROUTINE BCP04(N,VEC,ISEED)

C     _Authors:         F. Facchinei, J. Soares and J. Judice
C     _Date Written:    September 10, 1995
C     _Date Modified:   September  3, 1996
C     _Purpose:
C
C       Given a vector of integer components, BCP04 alters the position
C     of its elements randomly.
C
C     _Arguments in parameter list:
C
C     N      (input) INTEGER
C        The dimension of the integer vector.
C
C     VEC    (input/output) INTEGER vector of dimension N
C        The vector whose components positions are to be randomly
C        changed.
C
C     ISEED  (input/output) INTEGER
C        A seed for the random number generator. It is modified.

C     .. Scalar Arguments ..

      INTEGER ISEED,N
C     ..
C     .. Array Arguments ..

      INTEGER VEC(N)
C     ..
C     .. External Functions ..


      REAL RAND
      EXTERNAL RAND
C     ..
C     .. Local Scalars ..
      INTEGER COPY,I,II
C     ..

      DO 10 I = N,1,-1

C        VEC(J), for J=I+1,N will remain unaltered

          II = 1 + I*RAND(ISEED)

C        Exchange VEC(II) with VEC(I)

          COPY = VEC(I)
          VEC(I) = VEC(II)
          VEC(II) = COPY

C        VEC(J), for J=I,N will remain unaltered

   10 CONTINUE

      RETURN

      END

C     ..
C     .. End of subroutine BCP04 ..

C     ..
C     .. Begin of subroutines HI0, HI1 and HI2 altogether ..

C     _Authors:         F. Facchinei, J. Soares and J. Judice
C     _Date Written:    September 10, 1995
C     _Date Modified:   September  3, 1996
C     _Purpose:
C
C     The routine:
C
C       HI0 evaluates the function h_i
C       HI1 evaluates the derivative h'_i
C       HI2 evaluates the second derivative h''_i
C
C     ..
C     .. begin of subroutine HI0 ..

      REAL FUNCTION HI0(XI,XIBAR,BETAI,HI)

C     .. Scalar Arguments ..
      REAL BETAI,XI,XIBAR
      INTEGER HI
C     ..
C     .. Local Scalars ..
      REAL AUX,AUX2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..

      IF (HI.EQ.1) THEN

          HI0 = BETAI* (XI-XIBAR)

      ELSE IF (HI.EQ.2) THEN

          AUX = XI - XIBAR
          HI0 = AUX*AUX*AUX + BETAI*AUX

      ELSE IF (HI.EQ.3) THEN

          AUX = XI - XIBAR
          IF (AUX.NE.0.0) THEN
              AUX2 = AUX* (ABS(AUX)** (4.0E+0/3.0E+0))

          ELSE
              AUX2 = 0.0E+0
          END IF

          HI0 = AUX2 + BETAI*AUX

      END IF

      END

C     ..
C     .. end of subroutine HI0 ..

C     ..
C     .. begin of subroutine HI1 ..

      REAL FUNCTION HI1(XI,XIBAR,BETAI,HI)

C     .. Scalar Arguments ..
      REAL BETAI,XI,XIBAR
      INTEGER HI
C     ..
C     .. Local Scalars ..
      REAL AUX,AUX2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..

      IF (HI.EQ.1) THEN

          HI1 = BETAI

      ELSE IF (HI.EQ.2) THEN

          AUX = XI - XIBAR
          HI1 = 3.0E+0*AUX*AUX + BETAI

      ELSE IF (HI.EQ.3) THEN

          AUX = XI - XIBAR
          IF (AUX.NE.0.0) THEN
              AUX2 = (ABS(AUX))** (4.0E+0/3.0E+0)

          ELSE
              AUX2 = 0.0E+0
          END IF

          HI1 = 7.0E+0*AUX2/3.0E+0 + BETAI

      END IF

      END

C     ..
C     .. end of subroutine HI1 ..

C     ..
C     .. begin of subroutine HI2 ..

      REAL FUNCTION HI2(XI,XIBAR,HI)

C     .. Scalar Arguments ..
      REAL XI,XIBAR
      INTEGER HI
C     ..
C     .. Local Scalars ..
      REAL AUX,AUX2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..

      IF (HI.EQ.1) THEN

          HI2 = 0.0E+0

      ELSE IF (HI.EQ.2) THEN

          HI2 = 6E+0* (XI-XIBAR)

      ELSE IF (HI.EQ.3) THEN

          AUX = XI - XIBAR
          IF (AUX.NE.0.0) THEN
              AUX2 = (ABS(AUX))** (4.0E+0/3.0E+0)/AUX

          ELSE
              AUX2 = 0.0E+0
          END IF

          HI2 = 28.0E+0*AUX2/9.0E+0

      END IF

      END

C     ..
C     .. end of subroutine HI2 ..

C     ..
C     .. No more routines HI.

C     ..
C     .. Begin of subroutine BCP10 ..

      SUBROUTINE BCP10(N,X,F)

C     _Authors:         F. Facchinei, J. Soares and J. Judice
C     _Date Written:    September 10, 1995
C     _Date Modified:   September  3, 1996
C     _Purpose:
C
C       BCP10 evaluates the objective function of the generated problem
C     at a given point.
C
C     _Arguments in parameter list:
C
C     N    (input) INTEGER
C          The number of variables in the problem.
C
C     X    (input) DOUBLE PRECISION vector of dimension N
C          The point where the objective function is to be evaluated.
C
C     F    (output) DOUBLE PRECISION
C          The objective function value.

C     .. Scalar Arguments ..

      REAL F
      INTEGER N
C     ..
C     .. Array Arguments ..

      REAL X(N)
C     ..
C     .. External Subroutines ..

      EXTERNAL FUNCTG
C     ..
C     .. External Functions ..

      REAL HI0
      EXTERNAL HI0
C     ..
C     .. Parameters ..

C     .. See the comments in BCP01 for a description of the variables
C     in the common blocks.
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
      INTEGER I,II
C     ..
C     .. Common blocks ..
      COMMON /BCPB1/INL,INU,ACTL,ACTU,NM1L,NM1U,NM2L,NM2U,HI
      COMMON /BCPB2/BETA,XBAR,PTACT
C     ..

C     Evaluate the unconstrained objective function

      CALL FUNCTG(N,X,F)

C     Use the function h_i according to the selection in HI

      DO 10 II = 1,ACTL
          I = PTACT(II)
          F = F + HI0(X(I),XBAR(I),BETA(II),HI)
   10 CONTINUE

      DO 20 II = ACTL + 1,ACTU
          I = PTACT(II)
          F = F - HI0(X(I),XBAR(I),BETA(II),HI)
   20 CONTINUE

      RETURN

      END

C     ..
C     .. End of subroutine BCP10.

C     ..
C     .. Begin of subroutine BCP11 ..

      SUBROUTINE BCP11(N,X,GF)

C     _Authors:          F. Facchinei, J. Soares and J. Judice
C     _Date Written:     September 10, 1995
C     _Date Modified:    September  3, 1996
C     _Purpose:
C
C       BCP11 evaluates the gradient of the objective function of the
C     generated problem at a given point.
C
C     _Arguments in parameter list:
C
C     N    (input) INTEGER
C          The number of variables in the problem.
C
C     X    (input) DOUBLE PRECISION vector of dimension N
C          The point where the gradient of the objective function is to
C          be evaluated.
C
C     GF   (output) DOUBLE PRECISION vector of dimension N
C          The gradient of objective function values.

C     .. Scalar Arguments ..

      INTEGER N
C     ..
C     .. Array Arguments ..

      REAL GF(N),X(N)
C     ..
C     .. External Subroutines ..

      EXTERNAL GRADG
C     ..
C     .. External Functions ..

      REAL HI1
      EXTERNAL HI1
C     ..
C     .. Parameters ..

C     .. See the comments in BCP01 for a description of the variables
C     in the common blocks.
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
      INTEGER I,II
C     ..
C     .. Common blocks ..
      COMMON /BCPB1/INL,INU,ACTL,ACTU,NM1L,NM1U,NM2L,NM2U,HI
      COMMON /BCPB2/BETA,XBAR,PTACT
C     ..

C     Evaluate the unconstrained objective function gradient

      CALL GRADG(N,X,GF)

      DO 10 II = 1,ACTL
          I = PTACT(II)
          GF(I) = GF(I) + HI1(X(I),XBAR(I),BETA(II),HI)
   10 CONTINUE

      DO 20 II = ACTL + 1,ACTU
          I = PTACT(II)
          GF(I) = GF(I) - HI1(X(I),XBAR(I),BETA(II),HI)
   20 CONTINUE

      RETURN

      END

C     ..
C     .. End of subroutine BCP11.

C     ..
C     .. Begin of subroutine BCP12 ..

      SUBROUTINE BCP12(N,X,D,HFD)

C     _Authors:         F. Facchinei, J. Soares and J. Judice
C     _Date Written:    September 10, 1995
C     _Date Modified:   September  3, 1996
C     _Purpose:
C
C       BCP12 evaluates the product between the Hessian matrix of
C     the objective function at a given point and a vector.
C
C     _Arguments in parameter list:
C
C     N    (input) INTEGER
C          The number of variables in the problem.
C
C     X    (input) DOUBLE PRECISION vector of dimension N
C          The point where the Hessian matrix of the objective function
C          is to be evaluated.
C
C     D    (input) DOUBLE PRECISION vector of dimension N
C          The vector that will be multiplied by the Hessian matrix.
C
C     HFD  (output) DOUBLE PRECISION vector of dimension N
C          The resulting vector.

C     .. Scalar Arguments ..

      INTEGER N
C     ..
C     .. Array Arguments ..

      REAL D(N),HFD(N),X(N)
C     ..
C     .. External Subroutines ..

      EXTERNAL HESGD
C     ..
C     .. External Functions ..

      REAL HI2
      EXTERNAL HI2
C     ..
C     .. Parameters ..
C     .. See the comments in BCP01 for a description of the variables in
C     the common blocks.
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
      INTEGER I,II
C     ..
C     .. Common blocks ..
      COMMON /BCPB1/INL,INU,ACTL,ACTU,NM1L,NM1U,NM2L,NM2U,HI
      COMMON /BCPB2/BETA,XBAR,PTACT
C     ..

C     Evaluate the unconstrained objective function Hessian X d

      CALL HESGD(N,X,D,HFD)

      DO 10 II = 1,ACTL
          I = PTACT(II)
          HFD(I) = HFD(I) + HI2(X(I),XBAR(I),HI)*D(I)
   10 CONTINUE

      DO 20 II = ACTL + 1,ACTU
          I = PTACT(II)
          HFD(I) = HFD(I) - HI2(X(I),XBAR(I),HI)*D(I)
   20 CONTINUE

      RETURN

      END

C     ..
C     .. End of subroutine BCP12


C     ..
C     .. Begin of subroutine RAND

      REAL FUNCTION RAND(ISEED)
C     **********
C
C     function rand
C
C     Rand is the portable random number generator of L. Schrage.
C
C     The generator is full cycle, that is, every integer from
C     1 to 2**31 - 2 is generated exactly once in the cycle.
C     It is completely described in TOMS 5(1979),132-138.
C
C     The function statement is
C
C       real function rand(iseed)
C
C     where
C
C       iseed is a positive integer variable which specifies
C         the seed to the random number generator. Given the
C         input seed, rand returns a random number in the
C         open interval (0,1). On output the seed is updated.
C
C     Argonne National Laboratory. MINPACK Project. March 1981.
C
C     **********
C     .. Scalar Arguments ..
      INTEGER ISEED
C     ..
C     .. Local Scalars ..
      REAL C
      INTEGER A,B15,B16,FHI,K,LEFTLO,P,XALO,XHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT
C     ..
C     .. Data statements ..
C
C     Set a = 7**5, b15 = 2**15, b16 = 2**16, p = 2**31-1, c = 1/p.
C
      DATA A/16807/,B15/32768/,B16/65536/,P/2147483647/
      DATA C/4.656612875E-10/
C     ..
C
C     There are 8 steps in rand.
C
C     1. Get 15 hi order bits of iseed.
C     2. Get 16 lo bits of iseed and form lo product.
C     3. Get 15 hi order bits of lo product.
C     4. Form the 31 highest bits of full product.
C     5. Get overflo past 31st bit of full product.
C     6. Assemble all the parts and pre-substract p.
C        The parentheses are essential.
C     7. Add p back if necessary.
C     8. Multiply by 1/(2**31-1).
C
      XHI = ISEED/B16
      XALO = (ISEED-XHI*B16)*A
      LEFTLO = XALO/B16
      FHI = XHI*A + LEFTLO
      K = FHI/B15
      ISEED = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
      IF (ISEED.LT.0) ISEED = ISEED + P
      RAND = C*FLOAT(ISEED)
      RETURN
C
C     Last card of function rand.
C
      END
