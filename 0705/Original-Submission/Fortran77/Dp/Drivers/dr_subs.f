C--**--CH2894--705--A:1--28:10:1999
C--**--CH2879--705--P:R--28:10:1999
      DOUBLE PRECISION FUNCTION DGAUS(MEAN, STDDEV)
C
      DOUBLE PRECISION MEAN, STDDEV
C
C     THIS ROUTINE GENERATES NORMALLY DISTRIBUTED PSEUDORANDOM NUMBERS
C     WITH GIVEN MEAN AND STANDARD DEVIATION.
C     YOU MUST INITIALIZE DXRAND() WITH "X = DXRAND(I)" WHERE I>0.
C
C     WRITTEN -
C       14NOV86 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     MODIFIED -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
      EXTERNAL DXRAND
      DOUBLE PRECISION RBUF, V1, V2, S, DXRAND
      LOGICAL RBGOOD
      SAVE RBUF, RBGOOD
      DATA RBGOOD/.FALSE./
C
C     CHECK FOR RANDOM NUMBER LEFT OVER IN BUFFER
      IF (RBGOOD) THEN
         DGAUS = STDDEV*RBUF + MEAN
         RBGOOD = .FALSE.
         RETURN
      ENDIF
C
C     GET RANDOM VECTOR IN UNIT CIRCLE
   20 CONTINUE
         V1 = 2.0D0*DXRAND(0) - 1.0D0
         V2 = 2.0D0*DXRAND(0) - 1.0D0
         S = V1*V1 + V2*V2
      IF (S .GE. 1.0D0) GOTO 20
C
      S = SQRT(-2.0D0*LOG(S)/S)
      RBUF = V1*S
      RBGOOD = .TRUE.
      DGAUS = STDDEV*V2*S + MEAN
      RETURN
C --- LAST LINE OF DGAUS ---
      END
      SUBROUTINE MADD (NA,NB,NC,M,N,A,B,C)
C
C     *****PARAMETERS:
      INTEGER NA,NB,NC,M,N
      DOUBLE PRECISION A(NA,N),B(NB,N),C(NC,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE MATRIX SUM A+B AND STORES THE
C     RESULT IN THE ARRAY C.  ALL MATRICES ARE M X N.  THE SUM
C     MAY BE OVERWRITTEN INTO A (B) BY DESIGNATING THE ARRAY C
C     TO BE A (B).
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA,NB,NC         ROW DIMENSIONS OF THE ARRAYS CONTAINING A,B,
C                         AND C,RESPECTIVELY, AS DECLARED IN THE CALLING
C                         PROGRAM DIMENSION STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRICES A, B, AND C;
C
C        N                NUMBER OF COLUMNS OF THE MATRICES A, B, AND C;
C
C        A                AN M X N MATRIX;
C
C        B                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        C                AN M X N ARRAY CONTAINING A+B.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      DO 20 J=1,N
         DO 10 I=1,M
            C(I,J)=A(I,J)+B(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF MADD
C
      END
      SUBROUTINE MSAVE (NA,NAS,M,N,A,ASAVE)
C
C     *****PARAMETERS:
      INTEGER NA,NAS,M,N
      DOUBLE PRECISION A(NA,N),ASAVE(NAS,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COPIES THE CONTENTS OF THE M X N ARRAY A INTO
C     THE M X N ARRAY ASAVE.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA,NAS           ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND
C                         AS, RESPECTIVELY, AS DECLARED IN THE CALLING
C                         PROGRAM DIMENSION STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRICES A AND ASAVE;
C
C        N                NUMBER OF COLUMNS OF THE MATRICES A AND ASAVE;
C
C        A                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        ASAVE            AN M X N ARRAY CONTAINING THE ARRAY A.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      DO 20 J=1,N
         DO 10 I=1,M
            ASAVE(I,J)=A(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF MSAVE
C
      END
      SUBROUTINE MSUB (NA,NB,NC,M,N,A,B,C)
C
C     *****PARAMETERS:
      INTEGER NA,NB,NC,M,N
      DOUBLE PRECISION A(NA,N),B(NB,N),C(NC,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE MATRIX DIFFERENCE A-B AND STORES
C     THE RESULT IN THE ARRAY C.  ALL MATRICES ARE M X N.  THE
C     DIFFERENCE MAY BE OVERWRITTEN INTO A BY DESIGNATING THE
C     ARRAY C TO BE A.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA,NB,NC         ROW DIMENSIONS OF THE ARRAYS CONTAINING A,B,
C                         AND C, RESPECTIVELY, AS DECLARED IN THE
C                         CALLING PROGRAM DIMENSION STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRICES A,B, AND C;
C
C        N                NUMBER OF COLUMNS OF THE MATRICES A,B, AND C;
C
C        A                AN M X N MATRIX;
C
C        B                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        C                AN M X N ARRAY CONTAINING A-B.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      DO 20 J=1,N
         DO 10 I=1,M
            C(I,J)=A(I,J)-B(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF MSUB
C
      END
      SUBROUTINE MULC (NA,NB,NC,M,N,L,A,B,C)
C
C     *****PARAMETERS:
      INTEGER NA,NB,NC,L,M,N
      DOUBLE PRECISION A(NA,N),B(NB,L),C(NC,L)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT A*B AND STORES THE
C        RESULT IN THE ARRAY C.  A IS M X N, B IS N X L, AND C IS
C        M X L.  THE ARRAY C MUST BE DISTINCT FROM BOTH A AND B.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA    ROW DIMENSION OF THE ARRAY CONTAINING A AS DECLARED
C              IN THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        NB    ROW DIMENSION OF THE ARRAY CONTAINING B AS DECLARED
C              IN THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        NC    ROW DIMENSION OF THE ARRAY CONTAINING C AS DECLARED
C              IN THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        L     NUMBER OF COLUMNS OF THE MATRICES B AND C;
C
C        M     NUMBER OF ROWS OF THE MATRICES A AND C;
C
C        N     NUMBER OF COLUMNS OF THE MATRIX A AND NUMBER OF ROWS
C              OF THE MATRIX B;
C
C        A     AN M X N MATRIX;
C
C        B     AN N X L MATRIX.
C
C     ON OUTPUT:
C
C        C     AN M X L ARRAY CONTAINING A*B.
C
C     *****HISTORY:
C     ORIGINAL BY ALAN J. LAUB
C     MODIFIED AND RENAMED BY JUDITH D. GARDINER
C
C     ------------------------------------------------------------------
C
      DO 40 J=1,L
         DO 10 I=1,M
            C(I,J)=0.0D0
10       CONTINUE
         DO 30 K=1,N
            DO 20 I=1,M
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
20          CONTINUE
30       CONTINUE
40    CONTINUE
      RETURN
C
C     LAST CARD OF MULC
C
      END
      SUBROUTINE RANDM(A, LDA, M, N, DPAR, JOB)
C
      INTEGER LDA,N,M,JOB
      DOUBLE PRECISION A(LDA,N),DPAR(3)
C
C     THIS SUBROUTINE IS A GENERAL FACILITY FOR PRODUCING RANDOM MATRI-
C     CES.  SEE JOB PARAMETER FOR OPTIONS.  TO USE THIS ROUTINE DXRAND
C     -MUST- BE INITIALIZED BY THE USER.
C
C     ON ENTRY -
C       LDA     INTEGER
C               LEADING DECLARED DIMENSION OF  A.
C
C       M,N     INTEGER
C               ROW AND COLUMN DIMENSION OF  A , RESPECTIVELY.
C
C       DPAR    DOUBLE PRECISION (3)
C               AN ARRAY OF PARAMETERS (SEE JOB BELOW).
C
C       JOB     INTEGER
C               INTEGER IN DECIMAL FORM  ABCDE  INDICATING THE FOLLOWING
C                  A,B,C      (NOT USED, SHOULD BE SET TO ZERO)
C                  D .EQ. 0   MAKE RANDOM DENSE MATRICES
C                  D .EQ. 1   MAKE RANDOM SPARSE MATRICES WITH EACH
C                             ELEMENT HAVING PROBABILITY DPAR(3) OF
C                             BEING  ZERO.
C                  E .EQ. 0   USE UNIFORMLY DISTRIBUTED NUMBERS ON THE
C                             INTERVAL  [ DPAR(1), DPAR(2) ) .
C                  E .EQ. 1   USE GAUSSIAN RANDOMS WITH MEAN  DPAR(1),
C                             AND VARIANCE  DPAR(2).
C                  E .EQ. 2   USE EXPONENTIALLY DISTRIBUTED VARIABLES
C                             IN THE RANGE [ DPAR(1), DPAR(2) ).
C
C     ON RETURN -
C       A       DOUBLE PRECISION(LDA,N)
C               M BY N MATRIX OF RANDOMS
C
C     SUBROUTINES AND FUNCTIONS USED -
C       (MATU) DGAUS; (FORTRAN) MOD IDINT DXRAND;
C
C     WRITTEN -
C       18FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-4691
C
      INTEGER I,J,ITYPE
      DOUBLE PRECISION DIFF,DM,DS,DA,DB,THRESH
      LOGICAL SPARSE
C
      EXTERNAL DXRAND
      DOUBLE PRECISION DXRAND,DGAUS
C
      ITYPE = MOD(JOB,10)
      THRESH = 0.0D0
      SPARSE = (MOD(JOB/10,10) .NE. 0)
C
      IF (ITYPE .NE. 0) GOTO 60
C        --- UNIFORM RANDOM ---
         DIFF = DPAR(2) - DPAR(1)
         DA = DPAR(1)
         IF (SPARSE) GOTO 30
            DO 20 J = 1,N
               DO 10 I = 1,M
                  A(I,J) = DIFF*DXRAND(0) + DA
   10          CONTINUE
   20       CONTINUE
         RETURN
   30    CONTINUE
            DO 50 J = 1,N
               DO 40 I = 1,M
                  A(I,J) = 0.0D0
                  IF (DXRAND(0) .GT. DPAR(3))
     *                A(I,J) = DIFF*DXRAND(0) + DA
   40          CONTINUE
   50       CONTINUE
         RETURN
   60 CONTINUE
      IF (ITYPE .NE. 1) GOTO 120
C        --- GAUSSIAN RANDOMS ---
         DM = DPAR(1)
         DS = DPAR(2)
         IF (SPARSE) GOTO 90
            DO 80 J = 1,N
               DO 70 I = 1,M
                     A(I,J) = DGAUS(DM,DS)
   70          CONTINUE
   80       CONTINUE
            RETURN
   90    CONTINUE
            DO 110 J = 1,N
               DO 100 I = 1,M
                  IF (DXRAND(0) .GT. DPAR(3))
     *               A(I,J) = DGAUS(DM,DS)
  100          CONTINUE
  110       CONTINUE
            RETURN
  120 CONTINUE
      IF (ITYPE .NE. 2) GOTO 180
C        --- EXPONENTIAL ---
         DA = DPAR(1)
         DB = LOG(DPAR(2)/DA)
         IF (SPARSE) GOTO 150
            DO 140 J = 1,N
               DO 130 I = 1,M
                  A(I,J) = DA*EXP(DB)
  130          CONTINUE
  140       CONTINUE
            RETURN
  150    CONTINUE
            DO 170 J = 1,N
               DO 160 I = 1,M
                  IF (DXRAND(0) .GT. THRESH)
     *               A(I,J) = DA*EXP(DB)
  160          CONTINUE
  170       CONTINUE
            RETURN
  180 CONTINUE
      RETURN
C
C --- LAST LINE OF RANDM ---
      END
      DOUBLE PRECISION FUNCTION DXRAND(I)
      INTEGER I
      REAL RAND
      IF (I .GT. 0) THEN
        CALL SEED(I)
        DXRAND = 0.0D0
      ELSE
        DXRAND = DBLE(RAND())
      ENDIF
      END
