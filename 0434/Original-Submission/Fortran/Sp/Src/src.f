      SUBROUTINE CONP ( MATRIX, NR, NC, PT, PS, PC )
C
C  INPUT ARGUMENTS.
C
C    MATRIX = SPECIFICATION OF THE CONTINGENCY TABLE.
C      THIS MATRIX IS PARTIONED AS FOLLOWS:
C
C      X(11).....X(IC)  R(1)
C        .  .....  .      .
C        .  .....  .      .
C      X(R1).....X(RC)  R(R)
C       C(1)..... C(C)    N
C
C      WHERE X(IJ) ARE THE OBSERVED CELL FREQUENCIES,
C      R(I) ARE THE ROW TOTALS, C(J) ARE THE COLUMN
C      TOTALS, AND N IS THE TOTAL SAMPLE SIZE.
C      NOTE THAT THE ORIGINAL CELL FREQUENCIES ARE
C      DESTROYED BY THIS SUBROUTINE.
C
C    NR = THE NUMBER OF ROWS IN MATRIX (R=NR-1).
C
C    NC = THE NUMBER OF COLUMNS IN MATRIX (C=NC-1).
C
C  OUTPUT ARGUMENTS.
C
C    PT = THE PROBABILITY OF OBTAINING THE GIVEN TABLE.
C
C    PS = THE PROBABILITY OF OBTAINING A TABLE AS PROBABLE
C      AS, OR LESS PROBABLE THAN, THE GIVEN TABLE.
C
C    PC = THE PROBABILITY OF OBTAINING SOME OF THE
C      TABLES POSSIBLE WITHIN THE CONSTRAINTS OF THE
C      MARGINAL TOTALS.  (THIS SHOULD BE 1.0.  DEVIATIONS
C      FROM 1.0 REFLECT THE ACCURACY OF THE COMPUTATION.)
C
C  EXTERNALS
C
C    FACLOG(N) = FUNCTION TO RETURN THE FLOATING POINT
C      VALUE OF LOG BASE 10 OF N FACTORIAL.
C
C     .. Scalar Arguments ..
      REAL PC,PS,PT
      INTEGER NC,NR
C     ..
C     .. Array Arguments ..
      INTEGER MATRIX(NR,NC)
C     ..
C     .. Local Scalars ..
      REAL PX,QXLOG,RXLOG
      INTEGER C,I,J,R,TEMP
C     ..
C     .. External Functions ..
      REAL FACLOG
      EXTERNAL FACLOG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C
      R = NR - 1
      C = NC - 1
C
C  COMPUTE LOG OF CONSTANT NUMERATOR.
C
      QXLOG = -FACLOG ( MATRIX(NR,NC) )
      DO 10 I = 1, R
        QXLOG = QXLOG + FACLOG ( MATRIX(I,NC) )
10    CONTINUE
      DO 20 J = 1, C
        QXLOG = QXLOG + FACLOG ( MATRIX(NR,J) )
20    CONTINUE
C
C  COMPUTE PROBABILITY OF GIVEN TABLE.
C
      RXLOG = 0.0
      DO 49 I  =  1, R
        DO 50 J = 1, C
          RXLOG = RXLOG + FACLOG ( MATRIX(I,J) )
50    CONTINUE
49    CONTINUE
      PT = 10.0**( QXLOG - RXLOG )
C
      PS = 0.0
      PC = 0.0
C
C  FILL LOWER RIGHT (R-1) X (C-1) CELLS WITH
C  MINIMIMUM OF ROW AND COLUMN TOTALS.
C
      DO 99 I  =  2, R
        DO 100 J = 2, C
          MATRIX(I,J) = MIN ( MATRIX(I,NC), MATRIX(NR,J) )
100   CONTINUE
99    CONTINUE
      GO TO 300
C
C  OBTAIN A NEW SET OF FREQUENCIES IN
C  LOWER RIGHT (R-1) X (C-1) CELLS.
C
200   DO 219 I  =  2, R
        DO 220 J = 2, C
          MATRIX(I,J) = MATRIX(I,J) - 1
          IF ( MATRIX(I,J) .GE. 0 ) GO TO 300
          MATRIX(I,J) = MIN ( MATRIX(I,NC), MATRIX(NR,J) )
220   CONTINUE
219   CONTINUE
      RETURN
C
C  FILL REMAINDER OF OBSERVED CELLS
C  ...COMPLETE COLUMN 1.
C
300   DO 320 I = 2, R
        TEMP = MATRIX(I,NC)
        DO 310 J = 2, C
          TEMP = TEMP - MATRIX(I,J)
310   CONTINUE
        IF ( TEMP .LT. 0 ) GO TO 200
        MATRIX(I,1) = TEMP
320   CONTINUE
C
C  ...COMPLETE ROW 1.
C
      DO 340 J = 1, C
        TEMP = MATRIX(NR,J)
        DO 330 I = 2, R
          TEMP = TEMP - MATRIX(I,J)
330   CONTINUE
        IF ( TEMP .LT. 0 ) GO TO 200
        MATRIX(1,J) = TEMP
340   CONTINUE
C
C  COMPUTE LOG OF THE DENOMINATOR.
C
      RXLOG = 0.0
      DO 349 I  =  1, R
        DO 350 J = 1, C
          RXLOG = RXLOG + FACLOG ( MATRIX(I,J) )
350   CONTINUE
349   CONTINUE
C
C  COMPUTE PX.  ADD TO PS IF PX .LE. PT
C  (ALLOW FOR ROUND-OFF ERROR).
C
      PX = 10.0** ( QXLOG - RXLOG )
      PC = PC + PX
      IF ( ( PT / PX ) .GT. 0.99999 ) PS = PS + PX
      GO TO 200
      END
      REAL FUNCTION FACLOG ( N )
C
C  INPUT ARGUMENT.
C
C    N = AN INTEGER GREATER THAN OR EQUAL TO ZERO.
C
C  FUNCTION RESULT.
C
C    FACLOG = THE LOG TO THE BASE 10 OF N FACTORIAL.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Local Scalars ..
      REAL ELOG,TPILOG,X
      INTEGER I,IFLAG
C     ..
C     .. Local Arrays ..
      REAL TABLE(101)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG10,FLOAT
C     ..
C     .. Save statement ..
      SAVE ELOG,IFLAG,TABLE,TPILOG
C     ..
      DATA TPILOG / 0.3990899342 /
      DATA ELOG / 0.4342944819 /
      DATA IFLAG / 0 /
C
C  USE STIRLINGS APPROXIMATION IF N GT 100.
C
      IF ( N .GT. 100 ) GO TO 50
C
C  LOOKUP UP ANSWER IF TABLE WAS GENERATED.
C
      IF ( IFLAG .EQ. 0 ) GO TO 100
10    FACLOG = TABLE(N+1)
      RETURN
C
C  HERE FOR STIRLINGS APPROXIMATION.
C
50    X = FLOAT ( N )
      FACLOG = ( X + 0.5 ) * ALOG10 ( X ) - X * ELOG + TPILOG
     &  + ELOG / ( 12.0 * X ) - ELOG / ( 360.0 * X * X * X )
      RETURN
C
C  HERE TO GENERATE LOG FACTORIAL TABLE.
C
100   TABLE(1) = 0.0
      DO 120 I = 2, 101
        X = FLOAT ( I - 1 )
        TABLE(I) = TABLE(I-1) + ALOG10 ( X )
120   CONTINUE
      IFLAG = 1
      GO TO 10
      END
