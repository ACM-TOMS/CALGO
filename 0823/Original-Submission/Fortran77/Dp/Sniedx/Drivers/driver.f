      PROGRAM TNIEDX
C
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GNNIEDX"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
C      User Define:
C        DIMEN : dimension
C        ATMOST : sequence length
C        SAMS : Number of replications
C        MAXS : Maximum Digits of Scrambling Of Owen type Scrambling
C        IFLAG: User Choice of type Sequences
C        IFLAG = 0 : No Scrambling
C        IFLAG = 1 : Owen type Scrambling
C        IFLAG = 2 : Faure-Tezuka type Scrambling
C        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
C     .. Local Scalars ..
      DOUBLE PRECISION F,SUM
      INTEGER ATMOST,DIMEN,I,II,IFLAG,K,J,MAXS,SAM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION OUTS(40,10000)
      REAL SECOND, T1
C     ..

C     ..
C     .. External Subroutines ..
      EXTERNAL GNNEDX,TIME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD
C     ..
      SAM = 2
      MAXS = 31
      DIMEN = 2
      IFLAG = 1
      ATMOST = 2**12

      DO 30 II = 1,SAM
          WRITE (*,FMT=*) 'I = ITERATION NUMBER'
          WRITE (*,FMT=*) 'EI = ESTIMATED INTEGRAL'
          T1 = SECOND()
          CALL GNNEDX(DIMEN,ATMOST,IFLAG,MAXS,OUTS)
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
