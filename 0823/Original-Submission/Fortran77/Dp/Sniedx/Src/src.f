
      SUBROUTINE GNCRML(MAXS,LSM,SHIFT)

C     GENERATING LOWER TRIAGULAR SCRMABLING MATRICES
C     AND SHIFT VECTORS.
C     INPUTS :
C       FROM INNEDX : MAXS
C       FROM BLOCK DATA "NIEDX" : S, MAXCOL
C
C     OUTPUTS :
C       TO INNEDX : LSM, SHIFT


C     .. Scalar Arguments ..
      INTEGER MAXS
C     ..
C     .. Array Arguments ..
      INTEGER LSM(40,31),SHIFT(40)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,L,LL,P,STEMP,TEMP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION UNI
      EXTERNAL UNI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INT,MOD
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RECIPD
      INTEGER COUNT,MAXCOL,S
C     ..
C     .. Arrays in Common ..
      INTEGER LASTQ(40),SV(40,31)
C     ..
C     .. Common blocks ..
      COMMON /NIEDX/S,MAXCOL,SV,COUNT,LASTQ,RECIPD
C     ..
C     .. Save statement ..
      SAVE /NIEDX/
C     ..
      DO 30 P = 1,S
          SHIFT(P) = 0
          L = 1
          DO 20 I = MAXS,1,-1
              LSM(P,I) = 0
              STEMP = MOD((INT(UNI()*1000.0)),2)
              SHIFT(P) = SHIFT(P) + STEMP*L
              L = 2*L
              LL = 1
              DO 10 J = 30,1,-1
                  IF (J.EQ.I) THEN
                      TEMP = 1

                  ELSE IF (J.LT.I) THEN
                      TEMP = MOD((INT(UNI()*1000.0)),2)

                  ELSE
                      TEMP = 0
                  END IF

                  LSM(P,I) = LSM(P,I) + TEMP*LL
                  LL = 2*LL
   10         CONTINUE
   20     CONTINUE
   30 CONTINUE

      RETURN

      END

      SUBROUTINE GNCRMU(USM,USHIFT)

C     GENERATING UPPER TRIANGULAR SCRMABLING MATRICES AND
C     SHIFT VECTORS.
C     INPUTS :
C       FROM BLOCK DATA "NIEDX" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSOBL : USM, USHIFT


C     .. Array Arguments ..
      INTEGER USHIFT(30),USM(30,30)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,STEMP,TEMP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION UNI
      EXTERNAL UNI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INT,MOD
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RECIPD
      INTEGER COUNT,MAXCOL,S
C     ..
C     .. Arrays in Common ..
      INTEGER LASTQ(40),SV(40,31)
C     ..
C     .. Common blocks ..
      COMMON /NIEDX/S,MAXCOL,SV,COUNT,LASTQ,RECIPD
C     ..
      SAVE /NIEDX/
C     ..
      DO 20 I = 1,MAXCOL
          STEMP = MOD((INT(UNI()*1000.0)),2)
          USHIFT(I) = STEMP
          DO 10 J = 1,MAXCOL
              IF (J.EQ.I) THEN
                  TEMP = 1

              ELSE IF (J.GT.I) THEN
                  TEMP = MOD((INT(UNI()*1000.0)),2)

              ELSE
                  TEMP = 0
              END IF

              USM(I,J) = TEMP
   10     CONTINUE
   20 CONTINUE
      RETURN

      END

      SUBROUTINE GNNEDX(DIMEN,ATMOST,IFLAG,MAXS,OUTS)
C
C      User Define:
C        DIMEN : dimension
C        ATMOST : sequence length
C        MAXS : Maximum Digits of Scrambling
C        IFLAG: User Choice of Sequences
C        IFLAG = 0 : No Scrambling
C        IFLAG = 1 : Owen type Scrambling
C        IFLAG = 2 : Faure-Tezuka type Scrambling
C        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
C

C     .. Local Scalars ..
      DOUBLE PRECISION F,SUM
      INTEGER ATMOST,DIMEN,I,IFLAG,J,MAXS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION QUASI(40),OUTS(40,10000)
      LOGICAL FLAG(2)

C     .. External Subroutines ..
      EXTERNAL GONEDX,INNEDX

C     ..

      CALL INNEDX(FLAG,DIMEN,ATMOST,QUASI,MAXS,IFLAG)

      DO 10 J = 1,DIMEN
          OUTS(J,1) = QUASI(J)
   10 CONTINUE
      DO 30 I = 2,ATMOST
          CALL GONEDX(QUASI)
          DO 20 J = 1,DIMEN
              OUTS(J,I) = QUASI(J)
   20     CONTINUE
   30 CONTINUE

      RETURN

      END

      SUBROUTINE GONEDX(QUASI)
C
C     THIS SUBROUTINE GENERATES A NEW
C     QUASIRANDOM VECTOR WITH EACH CALL
C
C     IT ADAPTS THE IDEAS OF ANTONOV AND SALEEV,
C     USSR COMPUT. MATHS. MATH. PHYS. 19 (1980),
C     252 - 256
C
C     THE USER MUST CALL "INNEDX" BEFORE CALLING
C     "GONEDX".  AFTER CALLING "INNEDX", TEST
C     FLAG(1) AND FLAG(2);  IF EITHER IS FALSE,
C     DO NOT CALL "GONEDX".  "GOGONEDXCHECKS
C     THAT THE USER DOES NOT MAKE MORE CALLS
C     THAN HE SAID HE WOULD : SEE THE COMMENTS
C     TO "INNEDX".
C
C     INPUTS:
C       FROM USER'S CALLING PROGRAM:
C         NONE
C
C       FROM LABELLED COMMON /NIEDX/:
C         SV       TABLE OF DIRECTION NUMBERS
C         S        DIMENSION
C         MAXCOL   LAST COLUMN OF V TO BE USED
C         COUNT    SEQUENCE NUMBER OF THIS CALL
C         LASTQ    NUMERATORS FOR LAST VECTOR GENERATED
C         RECIPD   (1/DENOMINATOR) FOR THESE NUMERATORS
C
C
C     FIND THE POSITION OF THE RIGHT-HAND ZERO IN COUNT
C
C     .. Array Arguments ..
      DOUBLE PRECISION QUASI(40)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RECIPD
      INTEGER COUNT,MAXCOL,S
C     ..
C     .. Arrays in Common ..
      INTEGER LASTQ(40),SV(40,31)
C     ..
C     .. Local Scalars ..
      INTEGER I,L
C     ..
C     .. External Functions ..
      INTEGER EXOR
      EXTERNAL EXOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common blocks ..
      COMMON /NIEDX/S,MAXCOL,SV,COUNT,LASTQ,RECIPD
C     ..
C     .. Save statement ..
      SAVE /NIEDX/
C     ..
      L = 0
      I = COUNT
   10 L = L + 1
      IF (MOD(I,2).EQ.1) THEN
          I = I/2
          GO TO 10

      END IF
C
C     CHECK THAT THE USER IS NOT CHEATING !
C

      IF (L.GT.MAXCOL) STOP ' TOO MANY CALLS ON GOSOBL'
C
C     CALCULATE THE NEW COMPONENTS OF QUASI,
C     FIRST THE NUMERATORS, THEN NORMALIZED
C
      DO 20 I = 1,S
          LASTQ(I) = EXOR(LASTQ(I),SV(I,L))
C
C     IF A FULL-WORD EXCLUSIVE-OR, SAY .XOR., IS AVAILABLE
C     THEN REPLACE THE PRECEDING STATEMENT BY
C
C         LASTQ(I) = LASTQ(I) .XOR. SV(I,L)
C
C     TO GET A FASTER, EXTENDED FORTRAN PROGRAM
C
          QUASI(I) = LASTQ(I)*RECIPD
   20 CONTINUE
C
      COUNT = COUNT + 1
C
      RETURN

      END

      SUBROUTINE INNEDX(FLAG,DIMEN,ATMOST,QUASI,MAXS,IFLAG)
C
C     THIS INITIALIZE NIED-XING SEQUENCE.
C     FIRST CHECK WHETHER THE USER-SUPPLIED
C     DIMENSION "DIMEN" OF THE QUASI-RANDOM
C     VECTORS IS STRICTLY BETWEEN 1 AND 16.
C     IF SO, FLAG(1) = .TRUE.
C
C     NEXT CHECK "ATMOST", AN UPPER BOUND ON THE NUMBER
C     OF CALLS THE USER INTENDS TO MAKE ON "GOSOBL".  IF
C     THIS IS POSITIVE AND LESS THAN 2**30, THEN FLAG(2) = .TRUE.
C     (WE ASSUME WE ARE WORKING ON A COMPUTER WITH
C     WORD LENGTH AT LEAST 31 BITS EXCLUDING SIGN.)
C     THE NUMBER OF COLUMNS OF THE ARRAY V WHICH
C     ARE INITIALIZED IS
C          MAXCOL = NUMBER OF BITS IN ATMOST.
C     IN "GONEDX" WE CHECK THAT THIS IS NOT EXCEEDED.
C
C     NX STORED NIEDERREiITER-XING SEQUENCE GENERATORS
C     WHICH IS BINARY FRACTIONS WITH RESPECT TO ROWS.
C     THE GENERATORS CAN BE FOUND URL OF WEBSITE
C     http://www.dismat.oeaw.ac.at/pirs/niedxing.html
C
C
C     THE NUMBERS IN V ARE ACTUALLY BINARY FRACTIONS WITH
C     RESPECT TO COLUMNS.
C
C     LSM ARE LOWER TRIANGULAR SCRAMBLING MATRICES.
C     USM ARE UPPER TRIANGULAR SCRAMBLING MATRICES.
C     SV ARE SCRAMBLING GENERATOR MATRICES AND THE NUMBERS
C     ARE BINARY FRACTIONS.
C     "RECIPD" HOLDS 1/(THE COMMON DENOMINATOR OF ALL
C     OF THEM).
C
C
C     "INNEDX" IMPLICITLY COMPUTES THE FIRST SHIFTED
C     VECTOR "LASTQ", AND RETURN IT TO THE CALLING
C     PROGRAM. SUBSEQUENT VECTORS COME FROM "GONEDX".
C     "LASTQ" HOLDS NUMERATORS OF THE LAST VECTOR GENERATED.
C
C
C     Non-Standard Intrinsic Funtion for f77
C     But Standard Intrinsic Fuction for f90 IBITS IS USED.
C
C     INPUTS :
C       FROM USER'S PROGRAM : DIMEN, ATMOST, MAXS,
C                             IFLAG, QUASI
C
C     OUTPUTS :
C       TO USER'S PROGRAM : FLAG, LASTQ
C       TO "GONEDX" VIA /NIEDX/ :
C         SV, S, MAXCOL, COUNT, LASTQ, RECIPD
C
C

C
C     CHECK PARAMETERS
C
C     .. Scalar Arguments ..
      INTEGER ATMOST,DIMEN,IFLAG,MAXS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION QUASI(40)
      LOGICAL FLAG(2)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RECIPD
      INTEGER COUNT,MAXCOL,S
C     ..
C     .. Arrays in Common ..
      INTEGER LASTQ(40),SV(40,31)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION LL
      INTEGER I,J,K,L,P,PP,TEMP1,TEMP2,TEMP3,TEMP4
      CHARACTER*20 NAME
C     ..
C     .. Local Arrays ..
      INTEGER LSM(40,31),NX(16,30),SHIFT(40),TV(40,31,31),USHIFT(30),
     +        USM(30,30),V(40,30)
C     ..
C     .. External Functions ..
      INTEGER EXOR
      EXTERNAL EXOR
C     ..
C     .. External Subroutines ..
      EXTERNAL GNCRML,GNCRMU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IBITS,MOD
C     ..
C     .. Common blocks ..
      COMMON /NIEDX/S,MAXCOL,SV,COUNT,LASTQ,RECIPD
C     ..
C     .. Save statement ..
      SAVE /NIEDX/
C     ..
      S = DIMEN
      FLAG(1) = (S.GE.1 .AND. S.LE.16)
      FLAG(2) = (ATMOST.GT.0 .AND. ATMOST.LT.2**30)
      IF (.NOT. (FLAG(1).AND.FLAG(2))) RETURN
      I = ATMOST
      MAXCOL = 0
   10 MAXCOL = MAXCOL + 1
      I = I/2
      IF (I.GT.0) GO TO 10
C
C   GET GENERATOR MATRICES
C
      IF (S.LE.4) THEN
          NAME = 'nxs4m30'

      ELSE IF (S.EQ.5) THEN
          NAME = 'nxs5m30'

      ELSE IF (S.EQ.6) THEN
          NAME = 'nxs6m30'

      ELSE IF (S.LE.7) THEN
          NAME = 'nxs7m30'

      ELSE IF (S.EQ.8) THEN
          NAME = 'nxs8m30'

      ELSE IF (S.EQ.9) THEN
          NAME = 'nxs9m30'

      ELSE IF (S.EQ.10) THEN
          NAME = 'nxs10m30'

      ELSE IF (S.LE.11) THEN
          NAME = 'nxs11m30'

      ELSE IF (S.EQ.12) THEN
          NAME = 'nxs12m30'

      ELSE IF (S.EQ.13) THEN
          NAME = 'nxs13m30'

      ELSE IF (S.LE.14) THEN
          NAME = 'nxs14m30'

      ELSE IF (S.EQ.15) THEN
          NAME = 'nxs15m30'

      ELSE IF (S.EQ.16) THEN
          NAME = 'nxs16m30'
      END IF

      OPEN (UNIT=1,FILE=NAME,STATUS='UNKNOWN')
      DO 20 I = 1,S
          READ (1,FMT=*) (NX(I,J),J=1,30)
   20 CONTINUE
      CLOSE (UNIT=1,STATUS='KEEP')

      DO 40 I = 1,S
          DO 30 J = 1,MAXCOL
              V(I,J) = 0
   30     CONTINUE
   40 CONTINUE

      L = 1
      DO 70 J = 30,1,-1
          DO 60 I = 1,S
              DO 50 K = 1,MAXCOL
                  V(I,K) = V(I,K) + IBITS(NX(I,J),30-K,1)*L
   50         CONTINUE
   60     CONTINUE
          L = 2*L
   70 CONTINUE

C
C   COMPUTING  GENERATOR MATRICES OF USER CHOICE
C

      IF (IFLAG.EQ.0) THEN
          DO 90 I = 1,S
              DO 80 J = 1,MAXCOL
                  SV(I,J) = V(I,J)
   80         CONTINUE
              SHIFT(I) = 0
   90     CONTINUE
          LL = 2.0** (30)

      ELSE
          IF ((IFLAG.EQ.1) .OR. (IFLAG.EQ.3)) THEN
              CALL GNCRML(MAXS,LSM,SHIFT)
              DO 130 I = 1,S
                  DO 120 J = 1,MAXCOL
                      L = 1
                      TEMP2 = 0
                      DO 110 P = MAXS,1,-1
                          TEMP1 = 0
                          DO 100 K = 1,30
                              TEMP1 = TEMP1 + (IBITS(LSM(I,P),K-1,1)*
     +                                IBITS(V(I,J),K-1,1))
  100                     CONTINUE
                          TEMP1 = MOD(TEMP1,2)
                          TEMP2 = TEMP2 + TEMP1*L
                          L = 2*L
  110                 CONTINUE
                      SV(I,J) = TEMP2
  120             CONTINUE
  130         CONTINUE
              LL = 2.0** (MAXS)
          END IF

          IF ((IFLAG.EQ.2) .OR. (IFLAG.EQ.3)) THEN
              CALL GNCRMU(USM,USHIFT)
              IF (IFLAG.EQ.2) THEN
                  MAXS = 30
              END IF

              DO 190 I = 1,S
                  DO 150 J = 1,MAXCOL
                      P = MAXS
                      DO 140 K = 1,MAXS
                          IF (IFLAG.EQ.2) THEN
                              TV(I,P,J) = IBITS(V(I,J),K-1,1)

                          ELSE
                              TV(I,P,J) = IBITS(SV(I,J),K-1,1)
                          END IF

                          P = P - 1
  140                 CONTINUE
  150             CONTINUE

                  DO 180 PP = 1,MAXCOL
                      TEMP2 = 0
                      TEMP4 = 0
                      L = 1
                      DO 170 J = MAXS,1,-1
                          TEMP1 = 0
                          TEMP3 = 0
                          DO 160 P = 1,MAXCOL
                              TEMP1 = TEMP1 + TV(I,J,P)*USM(P,PP)
                              IF (PP.EQ.1) THEN
                                  TEMP3 = TEMP3 + TV(I,J,P)*USHIFT(P)
                              END IF

  160                     CONTINUE
                          TEMP1 = MOD(TEMP1,2)
                          TEMP2 = TEMP2 + TEMP1*L
                          IF (PP.EQ.1) THEN
                              TEMP3 = MOD(TEMP3,2)
                              TEMP4 = TEMP4 + TEMP3*L
                          END IF

                          L = 2*L
  170                 CONTINUE
                      SV(I,PP) = TEMP2
                      IF (PP.EQ.1) THEN
                          IF (IFLAG.EQ.3) THEN
                              SHIFT(I) = EXOR(TEMP4,SHIFT(I))

                          ELSE
                              SHIFT(I) = TEMP4
                          END IF

                      END IF

  180             CONTINUE
  190         CONTINUE
              LL = 2.0** (MAXS)
          END IF

      END IF
C
C     RECIPD IS 1/(COMMON DENOMINATOR OF THE ELEMENTS IN V)
C
      RECIPD = 1.0/LL

C
C     SET UP FIRST VECTOR AND VALUES FOR "GONEDX"
C     AND RETURN THE FIRST SEQUENCE
C
      COUNT = 0
      DO 200 I = 1,S
          LASTQ(I) = SHIFT(I)
          QUASI(I) = LASTQ(I)*RECIPD
  200 CONTINUE
      RETURN

      END
