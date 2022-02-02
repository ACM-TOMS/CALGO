
      DOUBLE PRECISION FUNCTION TESTF(N,DIMEN,QUASI)
C
C This version :  4 Mar 1992
C
C Provides a variety of test integrals for quasi-random
C sequences.  A call on TESTF computes an estimate of the
C integral ;  a call on EXACTF computes the exact value.
C
C     .. Scalar Arguments ..
      INTEGER DIMEN,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION QUASI(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,FLOAT,MOD,SIN
C     ..
C     .. Entry Points ..
      DOUBLE PRECISION EXACTF
C     ..
      GO TO (10,40,80,150) N
C
      ENTRY EXACTF(N,DIMEN)

      GO TO (30,60,140,170) N
C
C Test integral 1
C
   10 TESTF = 1.0
      DO 20 I = 1,DIMEN
          TESTF = TESTF*ABS(4*QUASI(I)-2)
   20 CONTINUE
      RETURN
C
   30 EXACTF = 1.0
      RETURN
C
C Test integral 2
C
   40 TESTF = 1.0
      DO 50 I = 1,DIMEN
          TESTF = TESTF*I*COS(I*QUASI(I))
   50 CONTINUE
      RETURN
C
   60 EXACTF = 1.0
      DO 70 I = 1,DIMEN
          EXACTF = EXACTF*SIN(FLOAT(I))
   70 CONTINUE
      RETURN
C
C Test integral 3
C
   80 TESTF = 1.0
      DO 130 I = 1,DIMEN
          X = 2*QUASI(I) - 1
          GO TO (90,100,110,120) MOD(I,4)

   90     TESTF = TESTF*X
          GO TO 130

  100     TESTF = TESTF* (2*X*X-1)
          GO TO 130

  110     TESTF = TESTF* (4*X*X-3)*X
          GO TO 130

  120     X = X*X
          TESTF = TESTF* (8*X*X-8*X+1)
  130 CONTINUE
      RETURN
C
  140 EXACTF = 0.0
      RETURN
C
C Test integral 4
C
  150 TESTF = 0
      X = 1
      DO 160 I = 1,DIMEN
          X = -X*QUASI(I)
          TESTF = TESTF + X
  160 CONTINUE
      RETURN
C
C
  170 X = 1.0/ (2** (DIMEN))
      IF (MOD(DIMEN,2).EQ.0) THEN
          EXACTF = (X-1)/3

      ELSE
          EXACTF = (X+1)/3
      END IF

      RETURN
C
      END
      PROGRAM TSNIED
C
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GNNIED"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
C      User Define:
C        DIMEN : dimension
C        ATMOST : sequence length
C        SAMS : Number of replications
C        IFLAG: User Choice of type Sequences
C        IFLAG = 0 : No Scrambling
C        IFLAG = 1 : Owen type Scrambling
C        IFLAG = 2 : Faure-Tezuka type Scrambling
C        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
C
C     .. Parameters ..
      INTEGER MAXDIM
      PARAMETER (MAXDIM=318)
C     .. Local Scalars ..
      DOUBLE PRECISION F,SUM
      INTEGER ATMOST,DIMEN,I,II,IFLAG,K,J,SAM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION OUTS(MAXDIM,10000)
      REAL SECOND, T1
C     ..

C     ..
C     .. External Subroutines ..
      EXTERNAL GNNIED,TIME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD
C     ..
      SAM = 1
      DIMEN = 2
      IFLAG = 1
      ATMOST = 2**12
      DO 30 II = 1,SAM
          WRITE (*,FMT=*) 'I = ITERATION NUMBER'
          WRITE (*,FMT=*) 'EI = ESTIMATED INTEGRAL'
          T1 = SECOND()
          CALL GNNIED(DIMEN,ATMOST,IFLAG,OUTS)
          SUM = 0.0
          K = 9
          DO 20 I = 1,ATMOST
              F = 1.0
              DO 10 J = 1,DIMEN
                  F = F*ABS(4.0*OUTS(J,I)-2.0)
   10         CONTINUE
              IF (MOD(I,2**K).EQ.0) THEN
                  K = K + 1
                  WRITE (*,FMT=*) 'I = ',I
                  WRITE (*,FMT=*) 'EI = ',SUM/I
              END IF

              SUM = SUM + F
   20     CONTINUE
          T1 = SECOND() - T1
          WRITE (*,FMT=*) 'Total time elapsed = ',T1
   30 CONTINUE
      STOP

      END
