      SUBROUTINE GNCRML(LSM,SHIFT,MAXL)
C     .. Parameters ..
      INTEGER MAXLN
      PARAMETER (MAXLN=30)

C     .. Scalar Arguments ..
      INTEGER MAXL
C     ..
C     .. Array Arguments ..
      INTEGER LSM(1:500,0:MAXLN,0:MAXLN),SHIFT(1:500,0:MAXLN)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RQS
      INTEGER HISUM,NEXTN,QS,S,TESTN
C     ..
C     .. Arrays in Common ..
      INTEGER COEF(1:500,0:MAXLN,0:MAXLN),PGTEMP(0:MAXLN),
     +        SCOEF(1:500,0:MAXLN,0:MAXLN),ZTEMP(1:500,0:MAXLN)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,P,QSM1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION UNI
      EXTERNAL UNI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INT,MOD
C     ..
C     .. Common blocks ..
      COMMON /FAURE/S,QS,HISUM,SCOEF,NEXTN,TESTN,RQS,PGTEMP,ZTEMP,COEF
C     ..
C     .. Save statement ..
      SAVE /FAURE/
C     ..
      QSM1 = QS - 1
      DO 30 P = 1,S
          DO 20 I = 0,MAXL
              SHIFT(P,I) = MOD((INT(UNI()*1000.0)),QS)
              DO 10 J = 0,HISUM
                  IF (J.EQ.I) THEN
                      LSM(P,I,J) = MOD((INT(UNI()*1000.0)),QSM1) + 1

                  ELSE IF (J.LT.I) THEN
                      LSM(P,I,J) = MOD((INT(UNI()*1000.0)),QS)

                  ELSE
                      LSM(P,I,J) = 0
                  END IF

   10         CONTINUE
   20     CONTINUE
   30 CONTINUE

      RETURN

      END

      SUBROUTINE GNCRMU(USM,USHIFT,MAXL)

C     GENERATING LOWER TRIANGULAR SCRMABLING MATRICES AND
C     SHIFT VECTORS.
C     INPUTS :
C       FROM INSOBL : MAXL
C       FROM BLOCK DATA "SOBOL" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSOBL : USM, USHIFT

C     .. Parameters ..
      INTEGER MAXLN
      PARAMETER (MAXLN=30)

C     .. Scalar Arguments ..
      INTEGER MAXL
C     ..
C     .. Array Arguments ..
      INTEGER USHIFT(0:MAXLN),USM(0:MAXLN,0:MAXLN)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RQS
      INTEGER HISUM,NEXTN,QS,S,TESTN
C     ..
C     .. Arrays in Common ..
      INTEGER COEF(1:500,0:MAXLN,0:MAXLN),PGTEMP(0:MAXLN),
     +        SCOEF(1:500,0:MAXLN,0:MAXLN),ZTEMP(1:500,0:MAXLN)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,QSM1,STEMP,TEMP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION UNI
      EXTERNAL UNI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INT,MOD
C     ..
C     .. Common blocks ..
      COMMON /FAURE/S,QS,HISUM,SCOEF,NEXTN,TESTN,RQS,PGTEMP,ZTEMP,COEF
C     ..
C     .. Save statement ..
      SAVE /FAURE/
C     ..
      QSM1 = QS - 1
      DO 20 I = 0,HISUM
          STEMP = MOD((INT(UNI()*1000.0)),QS)
C               STEMP = 0
          USHIFT(I) = STEMP
          DO 10 J = 0,HISUM
              IF (J.EQ.I) THEN
                  TEMP = MOD((INT(UNI()*1000.0)),QSM1) + 1
C                  TEMP = 1
              ELSE IF (J.GT.I) THEN
C                 TEMP = MOD((int(UNI()*1000.0)),QS)
                  TEMP = 0

              ELSE
                  TEMP = 0
              END IF

              USM(I,J) = TEMP
   10     CONTINUE
   20 CONTINUE
      RETURN

      END

      SUBROUTINE GNFAUR(DIMEN,ATMOST,IFLAG,MAXS,OUTS)
C
C       User Define:
C        DIMEN : dimension
C        ATMOST : sequence length
C        MAXS : Maximum Digits of Scrambling Of Owen type Scrambling
C        IFLAG: User Choice of Sequences
C        IFLAG = 0 : No Scrambling
C        IFLAG = 1 : Owen type Scrambling
C        IFLAG = 2 : Faure-Tezuka type Scrambling
C        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
C

C     .. Scalar Arguments ..
      INTEGER ATMOST,DIMEN,IFLAG,MAXS
C     ..

C     .. Array Argument..
      DOUBLE PRECISION OUTS(500,10000)
C      ..

C     .. Local Scalars ..
      INTEGER I,J,MAXL,MAXX
C     ..

C     .. Local Arrays ..
      DOUBLE PRECISION QUASI(500)
      LOGICAL FLAG(2)
C     ..

C     .. External Subroutines ..
      EXTERNAL GSFAUR,INFAUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD
C     ..
      MAXL = MAXS - 1

      CALL INFAUR(FLAG,DIMEN,ATMOST,QUASI,MAXL,IFLAG,MAXX)
      IF (.NOT.FLAG(2)) THEN
          WRITE (*,FMT=*) 'ATMOST = ',ATMOST
          WRITE (*,FMT=*) 'ATMOST IS NOT OK'
          STOP

      END IF

      DO 10 J = 1,DIMEN
          OUTS(J,1) = QUASI(J)
   10 CONTINUE

      DO 30 I = 2,ATMOST
          CALL GSFAUR(QUASI,MAXX)
          DO 20 J = 1,DIMEN
              OUTS(J,I) = QUASI(J)
   20     CONTINUE
   30 CONTINUE

      RETURN

      END

      SUBROUTINE GSFAUR(QUASI,MAXL)
C
C       THIS SUBROUTINE GENERATES A NEW
C       QUASIRANDOM VECTOR WITH EACH CALL BY USING THE ALPHA-ARY
C       GRAY CODE METHOD OF LICHTNER "SIAM J. DISCRETE MATH(1998),
C       381-386". (SEE ESPECIALLY PAGE 381-382).
C
C       THE USER MUST CALL "INFAUR" BEFORE
C       CALLING "GSFAUR".
C       AFTER CALLING "INFAUR", TEST FLAG(1)
C       AND FLAG(2); IF EITHER IS FALSE, DO
C       NOT CALL GOFAUR. READ THE COMMENTS AT
C       THE BEGINNING OF INFAUR AND THEN
C       THOSE BELOW.
C
C       ALL INPUTS COME FROM "INFAUR" VIA
C       LABELLED COMMON "FAURE"; FOR THEIR
C       DEFINITIONS, SEE "INFAUR".
C
C       INPUTS:
C         S,QS,COEF,NEXTN,TESTN,HISUM,RQS,MAXL
C
C       OUTPUTS:
C         TO USER'S CALLING PROGRAM:
C         QUASI - A NEW SCRAMBLED QUASIRANDOM VECTOR
C

C
C
C
C
C       NEXTN HAS A REPRESENTATION IN BASE
C       QS OF THE FORM: SUM OVER J FROM ZERO
C       TO HISUM OF YTEMP(J)*(QS**J)
C
C       WE NOW COMPUTE THE YTEMP(J)'S SAME AS FAURE
C
C     .. Parameters ..
      INTEGER MAXLN
      PARAMETER (MAXLN=30)
C
C     .. Scalar Arguments ..
      INTEGER MAXL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION QUASI(500)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RQS
      INTEGER HISUM,NEXTN,QS,S,TESTN
C     ..
C     .. Arrays in Common ..
      INTEGER COEF(1:500,0:MAXLN,0:MAXLN),PGTEMP(0:MAXLN),
     +        SCOEF(1:500,0:MAXLN,0:MAXLN),ZTEMP(1:500,0:MAXLN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R
      INTEGER DIFF,I,K,KTEMP,LTEMP,MTEMP,POS,TEM
C     ..
C     .. Local Arrays ..
      INTEGER GTEMP(0:MAXLN),YTEMP(0:19)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common blocks ..
      COMMON /FAURE/S,QS,HISUM,SCOEF,NEXTN,TESTN,RQS,PGTEMP,ZTEMP,COEF
C     ..
C     .. Save statement ..
      SAVE /FAURE/
C     ..
      KTEMP = TESTN
      LTEMP = NEXTN
      DO 10 I = HISUM,0,-1
          KTEMP = KTEMP/QS
          MTEMP = MOD(LTEMP,KTEMP)
          YTEMP(I) = (LTEMP-MTEMP)/KTEMP
          LTEMP = MTEMP
   10 CONTINUE
C
C PERFORM THE CONVERT DIGIT FROM B-ARY TO B-ARY GRAY CODE
C
      GTEMP(HISUM) = YTEMP(HISUM)
      DO 20 I = HISUM - 1,0,-1
          TEM = YTEMP(I) - YTEMP(I+1)
          IF (TEM.LT.0) THEN
              GTEMP(I) = TEM + QS

          ELSE
              GTEMP(I) = MOD(TEM,QS)
          END IF

   20 CONTINUE
C
C FINDING THE POSITION OF COLUMN THAT IT'S VALUE
C HAS BEEN CHANGED.
C
      DO 30 I = 0,HISUM
          IF (GTEMP(I).NE.PGTEMP(I)) THEN
              POS = I
              DIFF = GTEMP(I) - PGTEMP(I)
          END IF

          PGTEMP(I) = GTEMP(I)
   30 CONTINUE
C
C UPDATE THE NEW VECTOR BY ADDING THE VALUE OF THE COLUMN THAT
C HAS BEEN CHANGED.
C
      DO 50 K = 1,S
          R = 0.0
          DO 40 I = MAXL,0,-1
              ZTEMP(K,I) = ZTEMP(K,I) + SCOEF(K,I,POS)*DIFF
              R = MOD(ZTEMP(K,I),QS) + RQS*R
   40     CONTINUE
          QUASI(K) = R*RQS
   50 CONTINUE

C
C       UPDATE NEXTN AND, IF NEEDED, TESTN AND
C       HISUM
C
      NEXTN = NEXTN + 1
      IF (NEXTN.EQ.TESTN) THEN
          TESTN = TESTN*QS
          HISUM = HISUM + 1
C
C       SINCE FLAG(2) IS TRUE,
C       HISUM STAYS UNDER 20
C
      END IF
C
      RETURN

      END


C
C      THIS MODIFIED ALGORITHM 659, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 1, P.88.

      SUBROUTINE INFAUR(FLAG,DIMEN,ATMOST,QUASI,MAXL,IFLAG,MAXX)
C
C       THIS SUBROUTINE FIRST CHECKS WHETHER
C       THE USER-SUPPLIED DIMENSION "DIMEN" OF THE
C       QUASIRANDOM VECTORS IS ACCEPTABLE
C       (STRICTLY BETWEEN 1 AND 41) : IF SO,
C       FLAG(1)=.TRUE.
C
C       THEN IT CALCULATES AN UPPER SUMMATION
C       LIMIT "HISUM" BASED ON "DIMEN" AND THE
C       USER-SUPPLIED NUMBER "ATMOST" OF QUASIRANDOM
C       VECTORS REQUIRED. FLAG(2)=.TRUE. IF
C       ATMOST IS OK.
C
C       IF FLAG(1) AND FLAG(2) ARE TRUE,
C       "INFAUR" NEXT PRODUCES THE OTHER
C       OUTPUTS LISTED BELOW PASSED TO
C       SUBROUTINE GOFAUR VIA LABELLED
C       COMMON "FAURE". THESE OUTPUTS ARE
C       IRRELEVANT TO THE USER.
C
C       FIRST CALL INFAUR. IF FLAG(1) AND
C       FLAG(2) ARE TRUE, EACH (SUBSEQUENT)
C       CALL TO GSFAUR GENERATES A NEW
C       QUASIRANDOM VECTOR.
C
C       INPUTS : DIMEN, ATMOST, MAXL
C
C       OUTPUTS
C          TO USERS CALLING PROGRAM:
C             QUASI :FIRST SEQUENCE
C             FLAG
C             QSS   : SAME AS QS - SEE BELOW
C
C
C          TO GSFAUR:
C             S      :DIMENSION
C             QS     :SMALLEST PRIME >=S
C             SCOEF  :SCRAMBLING GENERATING MATRICES
C             NEXTN  :THE NUMBER OF THE
C                     NEXT QUASIRANDOM
C                     VECTOR,INITIALIZED
C                     TO TESTN-1 HERE.
C             TESTN  :INITIALIZED TO 0
C             HISUM  :AFTER BEING USED TO
C                     PRODUCE COEF, INITIALIZED
C                     TO 3 FOR GOFAUR.
C             RQS    :1.0/QS.
C
C
C     .. Parameters ..
      INTEGER MAXLN
      PARAMETER (MAXLN=30)
C
C
C     .. Scalar Arguments ..
      INTEGER ATMOST,DIMEN,IFLAG,MAXL,MAXX
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION QUASI(500)
      LOGICAL FLAG(2)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION RQS
      INTEGER HISUM,NEXTN,QS,S,TESTN
C     ..
C     .. Arrays in Common ..
      INTEGER COEF(1:500,0:MAXLN,0:MAXLN),PGTEMP(0:MAXLN),
     +        SCOEF(1:500,0:MAXLN,0:MAXLN),ZTEMP(1:500,0:MAXLN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R
      INTEGER I,IN,J,K,K1,K2,L,TEMS,TEMS2
C     ..
C     .. Local Arrays ..
      INTEGER LSM(1:500,0:MAXLN,0:MAXLN),PRIMES(0:4,0:99),
     +        SHIFT(1:500,0:MAXLN),TSCOEF(1:500,0:MAXLN,0:MAXLN),
     +        USHIFT(0:MAXLN),USM(0:MAXLN,0:MAXLN)
C     ..
C     .. External Subroutines ..
      EXTERNAL GNCRML,GNCRMU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,MOD,NINT,REAL,INT
C     ..
C     .. Common blocks ..
      COMMON /FAURE/S,QS,HISUM,SCOEF,NEXTN,TESTN,RQS,PGTEMP,ZTEMP,COEF
C     ..
C     .. Save statement ..
      SAVE /FAURE/
C     ..
C     .. Data statements ..
C
      DATA (PRIMES(0,I),I=0,99)/0,1,2,3,5,5,7,7,11,11,11,11,13,13,17,17,
     +     17,17,19,19,23,23,23,23,29,29,29,29,29,29,31,31,37,37,37,37,
     +     37,37,41,41,41,41,43,43,47,47,47,47,53,53,53,53,53,53,59,59,
     +     59,59,59,59,61,61,67,67,67,67,67,67,71,71,71,71,73,73,79,79,
     +     79,79,79,79,83,83,83,83,89,89,89,89,89,89,97,97,97,97,97,97,
     +     97,97,101,101/
      DATA (PRIMES(1,I),I=0,99)/101,101,103,103,107,107,107,107,109,109,
     +     113,113,113,113,127,127,127,127,127,127,127,127,127,127,127,
     +     127,127,127,131,131,131,131,137,137,137,137,137,137,139,139,
     +     149,149,149,149,149,149,149,149,149,149,151,151,157,157,157,
     +     157,157,157,163,163,163,163,163,163,167,167,167,167,173,173,
     +     173,173,173,173,179,179,179,179,179,179,181,181,191,191,191,
     +     191,191,191,191,191,191,191,193,193,197,197,197,197,199,199/
      DATA (PRIMES(2,I),I=0,99)/211,211,211,211,211,211,211,211,211,211,
     +     211,211,223,223,223,223,223,223,223,223,223,223,223,223,227,
     +     227,227,227,229,229,233,233,233,233,239,239,239,239,239,239,
     +     241,241,251,251,251,251,251,251,251,251,251,251,257,257,257,
     +     257,257,257,263,263,263,263,263,263,269,269,269,269,269,269,
     +     271,271,277,277,277,277,277,277,281,281,281,281,283,283,293,
     +     293,293,293,293,293,293,293,293,293,307,307,307,307,307,307/
      DATA (PRIMES(3,I),I=0,99)/307,307,307,307,307,307,307,307,311,311,
     +     311,311,313,313,317,317,317,317,331,331,331,331,331,331,331,
     +     331,331,331,331,331,331,331,337,337,337,337,337,337,347,347,
     +     347,347,347,347,347,347,347,347,349,349,353,353,353,353,359,
     +     359,359,359,359,359,367,367,367,367,367,367,367,367,373,373,
     +     373,373,373,373,379,379,379,379,379,379,383,383,383,383,389,
     +     389,389,389,389,389,397,397,397,397,397,397,397,397,401,401/
      DATA (PRIMES(4,I),I=0,99)/401,401,409,409,409,409,409,409,409,409,
     +     419,419,419,419,419,419,419,419,419,419,421,421,431,431,431,
     +     431,431,431,431,431,431,431,433,433,439,439,439,439,439,439,
     +     443,443,443,443,449,449,449,449,449,449,457,457,457,457,457,
     +     457,457,457,461,461,461,461,463,463,467,467,467,467,479,479,
     +     479,479,479,479,479,479,479,479,479,479,487,487,487,487,487,
     +     487,487,487,491,491,491,491,499,499,499,499,499,499,499,499/
C
C       CHECK S
C
      S = DIMEN
      FLAG(1) = S .GT. 1 .AND. S .LT. 41
      IF (.NOT.FLAG(1)) RETURN
C
      TESTN = 0
      K1 = INT(S/100)
      K2 = (S-K1*100)
      QS = PRIMES(K1,K2)

C
C         COMPUTE LOG(ATMOST+TESTN) IN BASE QS
C         USING A RATIO OF NATURAL LOGS TO GET
C         AN UPPER BOUND ON (THE NUMBER OF
C         DIGITS IN THE BASE QS REPRESENTATION
C         OF ATMOST+TESTN) MINUS ONE.
C
      HISUM = NINT(LOG(REAL(ATMOST+TESTN))/LOG(REAL(QS)))
      FLAG(2) = HISUM .LT. 30
      IF (.NOT.FLAG(2)) RETURN
C
C        NOW FIND BINOMIAL COEFFICIENTS MOD QS
C        IN A UPPER-TRIANGULAR MATRIX "COEF"
C        USING RECURSION BINOM(I,J)=BINOM(I,J-1)
C        +BINOM(I-1,J-1) AND A=B+C IMPLIES MOD(A,D)=
C        MOD(MOD(B,D)+MOD(C,D),D)
C
C CONSTRUCTING UPPER-TIRANGULAR GENERATOR MATRICES
C UP TO DIMENSION S

      DO 20 I = 0,HISUM
          DO 10 J = 0,HISUM
              IF (I.EQ.J) THEN
                  COEF(1,I,J) = 1

              ELSE
                  COEF(1,I,J) = 0
              END IF

   10     CONTINUE
   20 CONTINUE

      COEF(2,0,0) = 1
      DO 30 J = 1,HISUM
          COEF(2,0,J) = 1
          COEF(2,J,J) = 1
   30 CONTINUE

      DO 50 I = 1,HISUM
          DO 40 J = I + 1,HISUM
              COEF(2,I,J) = MOD(COEF(2,I,J-1)+COEF(2,I-1,J-1),QS)
   40     CONTINUE
   50 CONTINUE


      DO 90 K = 3,S
          DO 80 J = 0,HISUM
              DO 70 I = 0,HISUM
                  COEF(K,I,J) = 0
                  IF (J.LT.I) THEN
                      GO TO 70

                  ELSE
                      TEMS = 0
                  END IF

                  DO 60 L = 0,J
                      TEMS = TEMS + COEF(K-1,I,L)*COEF(2,L,J)
   60             CONTINUE
                  COEF(K,I,J) = MOD(TEMS,QS)
   70         CONTINUE
   80     CONTINUE
   90 CONTINUE


C
C GENERATE USER CHOICE OF GENERATOR MATRICES
C
      IF (IFLAG.EQ.0) THEN
          MAXX = HISUM
          DO 120 K = 1,S
              DO 110 J = 0,HISUM
                  DO 100 I = 0,HISUM
                      SCOEF(K,I,J) = COEF(K,I,J)
  100             CONTINUE
                  SHIFT(K,J) = 0

  110         CONTINUE
  120     CONTINUE

      ELSE
          IF ((IFLAG.EQ.1) .OR. (IFLAG.EQ.3)) THEN
              MAXX = MAXL
              CALL GNCRML(LSM,SHIFT,MAXL)
              DO 160 K = 1,S
                  DO 150 J = 0,HISUM
                      DO 140 I = 0,MAXL
                          IF (K.EQ.1) THEN
                              SCOEF(K,I,J) = LSM(K,I,J)
                              IF (IFLAG.EQ.3) THEN
                                  TSCOEF(K,I,J) = SCOEF(K,I,J)
                              END IF

                          ELSE
                              TEMS = 0
                              IF (I.GE.HISUM) THEN
                                  IN = HISUM

                              ELSE
                                  IN = J
                              END IF

                              DO 130 L = 0,IN
                                  TEMS = TEMS + LSM(K,I,L)*COEF(K,L,J)
  130                         CONTINUE
                              SCOEF(K,I,J) = MOD(TEMS,QS)
                              IF (IFLAG.EQ.3) THEN
                                  TSCOEF(K,I,J) = SCOEF(K,I,J)
                              END IF

                          END IF

  140                 CONTINUE
  150             CONTINUE
  160         CONTINUE
          END IF

          IF ((IFLAG.EQ.2) .OR. (IFLAG.EQ.3)) THEN
              CALL GNCRMU(USM,USHIFT,HISUM)
              IF (IFLAG.EQ.2) THEN
                  MAXX = HISUM
              END IF

              DO 200 K = 1,S
                  DO 190 J = 0,HISUM
                      DO 180 I = 0,MAXX
                          TEMS = 0
                          TEMS2 = 0
                          IN = HISUM
                          DO 170 L = 0,IN
                              IF (IFLAG.EQ.2) THEN
                                  TEMS = TEMS + COEF(K,I,L)*USM(L,J)
                                  IF (J.EQ.0) THEN
                                      TEMS2 = TEMS2 +
     +                                        COEF(K,I,L)*USHIFT(L)
                                  END IF

                              END IF

                              IF (IFLAG.EQ.3) THEN
                                  TEMS = TEMS + TSCOEF(K,I,L)*USM(L,J)
                                  IF (J.EQ.0) THEN
                                      TEMS2 = TEMS2 +
     +                                        TSCOEF(K,I,L)*USHIFT(L)
                                  END IF

                              END IF

  170                     CONTINUE
                          SCOEF(K,I,J) = MOD(TEMS,QS)
                          IF (J.EQ.0) THEN
                              TEMS2 = MOD(TEMS2,QS)
                              IF (IFLAG.EQ.3) THEN
                                  SHIFT(K,I) = MOD((TEMS2+SHIFT(K,I)),
     +                                         QS)

                              ELSE
                                  SHIFT(K,I) = TEMS2
                              END IF

                          END IF

  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
          END IF

      END IF
C        CALCULATING THESE COEFFICIENTS
C        MOD QS AVOIDS POSSIBLE OVERFLOW
C        PROBLEMS WITH RAW BINOMIAL COEFFICIENTS
C

      RQS = 1.0/REAL(QS)

C     GENERATE FIRST SCRAMBLED QUASI RANDOM VECTOR

      DO 220 K = 1,S
          R = 0.0
          DO 210 I = MAXX,0,-1
              ZTEMP(K,I) = SHIFT(K,I)
              R = ZTEMP(K,I) + RQS*R
  210     CONTINUE
          QUASI(K) = R*RQS
  220 CONTINUE
C
C     SET UP FIRST VECTOR AND VALUES FOR GSFAUR
C
      DO 230 I = 0,HISUM
          PGTEMP(I) = 0
  230 CONTINUE

      TESTN = QS
      NEXTN = 1
      HISUM = 0

C     NOW COMPLETE INITIALIZATION
C     AS DESCRIBED IN SECTION 4.
C     NEXTN HAS 1 DIGITS IN BASE
C     QS, SO HISUM EQUALS 0.
C
      RETURN

      END
