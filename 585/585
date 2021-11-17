C     ALGORITHM 585, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.8, NO. 3,
C     SEP., 1982, P. 290.
C     PROGRAM TESTD (OUTPUT,TAPE6=OUTPUT)                               MAN   10
C     TEST DRIVER FOR THE EXTRAP SUBROUTINE IN THE                      MAN   20
C     MUHLBACH-NEVILLE-AITKEN ALGORITHM                                 MAN   30
C                                                                       MAN   40
      INTEGER I, INIT, INFO, K, L, M, MPI, N, NI                        MAN   50
      DOUBLE PRECISION A, GINIT, F, P, RESULT, SK, W, X, EPS, T         MAN   60
      DIMENSION A(10,10), GINIT(40), F(10), W(10)                       MAN   70
C                                                                       MAN   80
      DATA N /4/, NOUT /6/, EPS /1.D-30/                                MAN   90
C                                                                       MAN  100
C     N IS THE NUMBER OF TERMS IN THE LINEAR COMBINATION                MAN  110
C                                                                       MAN  120
C     W(I) MUST CONTAIN W                                               MAN  130
C                        I-1                                            MAN  140
C                                                                       MAN  150
      W(1) = 3.D0*DEXP(-1.D0) + 4.D0*DSIN(1.D0) - 1.D0                  MAN  160
      W(2) = -3.D0*DEXP(-2.D0) + 4.D0*DCOS(2.D0) - 3.5D0                MAN  170
      W(3) = 6.D0 - 1.D0/3.D0 - 3.D0*DEXP(-1.D0) - 4.D0*DCOS(1.D0)      MAN  180
      W(4) = 3.D0*DEXP(-1.D0) - 4.D0*DSIN(1.D0) - 3.D0                  MAN  190
C                                                                       MAN  200
C     A(I,J) MUST CONTAIN L   (F   )                                    MAN  210
C                          I-1  J-1                                     MAN  220
C                                                                       MAN  230
      A(1,1) = DEXP(-1.D0)                                              MAN  240
      A(2,1) = -DEXP(-2.D0)                                             MAN  250
      A(3,1) = 1.D0 - A(1,1)                                            MAN  260
      A(4,1) = A(1,1)                                                   MAN  270
      A(1,2) = DSIN(1.D0)                                               MAN  280
      A(2,2) = DCOS(2.D0)                                               MAN  290
      A(3,2) = 1.D0 - DCOS(1.D0)                                        MAN  300
      A(4,2) = -A(1,2)                                                  MAN  310
      A(1,3) = 1.D0                                                     MAN  320
      A(2,3) = 4.D0                                                     MAN  330
      A(3,3) = 1.D0/3.D0                                                MAN  340
      A(4,3) = 2.D0                                                     MAN  350
      A(1,4) = 0.D0                                                     MAN  360
      A(2,4) = 0.5D0                                                    MAN  370
      A(3,4) = -1.D0                                                    MAN  380
      A(4,4) = -1.D0                                                    MAN  390
C                                                                       MAN  400
      DO 50 L=1,10                                                      MAN  410
        T = L                                                           MAN  420
        X = T*0.2D0                                                     MAN  430
C                                                                       MAN  440
C     F(I) MUST CONTAIN F                                               MAN  450
C                        I-1                                            MAN  460
C                                                                       MAN  470
        F(1) = DEXP(-X)                                                 MAN  480
        F(2) = DSIN(X)                                                  MAN  490
        F(3) = X*X                                                      MAN  500
        F(4) = DLOG(X)                                                  MAN  510
        INIT = 0                                                        MAN  520
        DO 40 K=1,N                                                     MAN  530
          SK = W(K)*F(1)/A(K,1)                                         MAN  540
          M = MAX0(1,K-2)                                               MAN  550
C
C                     (K-1)
C     COMPUTATION OF F
C                     0,I
C
          DO 10 I=1,M
            GINIT(I) = A(K,I+1)*F(1)/A(K,1) - F(I+1)
   10     CONTINUE
          IF (K.LT.3) GO TO 30
C
C                     (I-1)
C     COMPUTATION OF F
C                     0,K-1
C
          DO 20 I=1,K
            MPI = M + I
            GINIT(MPI) = A(I,K)*F(1)/A(I,1) - F(K)
   20     CONTINUE
   30     CALL EXTRAP(INIT, EPS, N-1, SK, GINIT, RESULT, INFO)
   40   CONTINUE
        IF (INFO.LT.0) WRITE (NOUT,99998)
        IF (INFO.EQ.1) WRITE (NOUT,99997)
        P = 3.D0*F(1) + 4.D0*F(2) - F(3) + F(4)
        WRITE (NOUT,99999) X, P, RESULT
   50 CONTINUE
      STOP
99999 FORMAT (3D26.12)
99998 FORMAT (46H DIMENSIONS TOO SMALL IN THE SUBROUTINE EXTRAP)
99997 FORMAT (42H DIVISION BY ZERO IN THE SUBROUTINE EXTRAP)
      END
C     PROGRAM TESTD (OUTPUT,TAPE6=OUTPUT)                               MAN   10
C     TEST DRIVER FOR THE EXTRAP SUBROUTINE                             MAN   20
C     IN THE E-ALGORITHM                                                MAN   30
C                                                                       MAN   40
      INTEGER I, MAXCOL, INIT, INFO, K, M, MPI, NI                      MAN   50
      DOUBLE PRECISION KPI, SK, RESULT, S, GINIT, EPS, T                MAN   60
      DIMENSION GINIT(50)                                               MAN   70
C                                                                       MAN   80
      DATA NOUT /6/, INIT /0/, MAXCOL /5/, EPS /1.D-30/                 MAN   90
C                                                                       MAN  100
      DO 40 K=1,20                                                      MAN  110
C                                                                       MAN  120
C     COMPUTATION OF S                                                  MAN  130
C                     K-1                                               MAN  140
C                                                                       MAN  150
        T = K                                                           MAN  160
        SK = 1.D0/T                                                     MAN  170
C                                                                       MAN  180
C     COMPUTATION OF G (K-1)                                            MAN  190
C                     I                                                 MAN  200
C                                                                       MAN  210
        M = MAX0(1,K-2)                                                 MAN  220
        DO 10 I=1,M
          KPI = K + I
          GINIT(I) = -1.D0/(KPI-1.D0)/KPI
   10   CONTINUE
        IF (K.LT.3 .OR. K.GT.MAXCOL+1) GO TO 30
C
C     COMPUTATION OF G   (I-1)
C                     K-1
C
        DO 20 I=1,K
          KPI = K + I
          MPI = M + I
          GINIT(MPI) = -1.D0/(KPI-2.D0)/(KPI-1.D0)
   20   CONTINUE
   30   CALL EXTRAP(INIT, EPS, MAXCOL, SK, GINIT, RESULT, INFO)
        IF (INFO.LT.0) WRITE (NOUT,99998)
        IF (INFO.EQ.1) WRITE (NOUT,99997)
        T = K
        S = 1.D0/T**2
        T = (MAXCOL+1)*K
        IF (K.GT.MAXCOL+1) S = 1.D0/T
        WRITE (NOUT,99999) RESULT, S
   40 CONTINUE
      STOP
99999 FORMAT (3D26.12)
99998 FORMAT (46H DIMENSIONS TOO SMALL IN THE SUBROUTINE EXTRAP)
99997 FORMAT (42H DIVISION BY ZERO IN THE SUBROUTINE EXTRAP)
      END
      SUBROUTINE EXTRAP(INIT, EPS, MAXCOL, SK, GINIT, RESULT, INFO)     EXT   10
C
C   STATEMENT OF PURPOSE
C   ********************
C
C   THE SUBPROGRAM EXTRAP() IS A SUBROUTINE TO IMPLEMENT THE E-ALGORITHM
C   FOR SEQUENCE EXTRAPOLATION AND THE MUHLBACH-NEVILLE-AITKEN ALGORITHM
C   FOR THE GENERAL INTERPOLATION PROBLEM.
C
C
C   DESCRIPTION OF THE PARAMETERS OF THE SUBROUTINE EXTRAP()
C   ******************************************************
C
C
C INIT   INTEGER PARAMETER TO BE SET TO ZERO BEFORE THE FIRST
C        CALL OF THE SUBROUTINE.  DUE TO THE FACT THAT AFTER THE
C        COMPUTATION OF EACH ROW OF THE E-ARRAY ONE HAS TO EXIT THE
C        SUBROUTINE TO COMPUTE THE NEW VALUES OF THE DATA SK AND
C        GINIT(*), ALL THE VARIABLES INTERNAL TO THE SUBROUTINE MUST
C        REMAIN INTACT BETWEEN THE CALLS. THE SAME IS TRUE FOR THE
C        ARGUMENTS INIT , MAXCOL AND INFO.  THE VALUE OF EPS CAN BE
C        CHANGED.  DURING THE FIRST CALL OF THE SUBROUTINE INIT IS
C        CHANGED TO 1. TO USE THE SUBROUTINE FOR A NEW APPLICATION OF
C        THE ALGORITHM (AND NOT FOR THE SUBSEQUENT CALLS CORRESPONDING
C        TO THE SAME APPLICATION) INIT WILL HAVE TO BE SET AGAIN TO 0.
C
C EPS    THE AIM OF THIS PARAMETER IS TO AVOID A DIVISION BY ZERO.
C        WHEN THE ABSOLUTE VALUE OF A DENOMINATOR IS LESS THAN EPS THE
C        PARAMETER INFO IS SET TO 1.
C
C MAXCOL INTEGER PARAMETER GIVING THE INDEX OF THE LAST COLUMN
C        OF THE E-ARRAY THAT THE USER WANTS TO COMPUTE (SEE
C        BELOW IN THE EXPLANATIONS FOR THE PARAMETER RESULT.)
C
C SK     DOUBLE PRECISION ARGUMENT WHICH MUST CONTAIN THE VALUE
C        OF S    BEFORE THE K-TH CALL OF THE SUBROUTINE.
C            K-1
C
C GINIT  DOUBLE PRECISION ARRAY WHICH MUST CONTAIN THE VALUES OF
C        G (0) BEFORE THE FIRST CALL OF THE SUBROUTINE, G (1)
C         1                                              1
C        BEFORE THE SECOND CALL, G (K-1), G (K-1),...., G   (K-1),
C                                 1        2             K-2
C        G   (0), G   (1),...., G   (K-1) BEFORE THE K-TH CALL
C         K-1      K-1           K-1
C        IF MAXCOL+1 >= K >= 3 AND G (K-1),...., G      (K-1) IF
C                                   1             MAXCOL
C        K > MAXCOL + 1
C.
C RESULT  DOUBLE PRECISION ARGUMENT CONTAINING THE ANSWER GIVEN BY
C        THE SUBROUTINE AFTER THE CALL. RESULT SUCCESSIVELY CONTAINS
C         (0)    (0)         (0)       (1)       (2)
C        E    , E    ,...., E       , E       , E
C         0      1           MAXCOL    MAXCOL    MAXCOL
C
C INFO   INTEGER PARAMETER WHICH, ON NORMAL EXIT, HAS THE VALUE 0.
C        A NEGATIVE VALUE SIGNALS THAT DIMENSIONS ARE TOO SMALL IN
C        THE SUBROUTINE. THE VALUE 1 INDICATES A DIVISION BY A
C        QUANTITY SMALLER THAN EPS IN ABSOLUTE VALUE. IF INFO HAS A
C        VALUE DIFFERENT FROM ZERO IT IS IMPOSSIBLE TO ENTER AGAIN
C        INTO EXTRAP(). TO USE IT AGAIN THE PARAMETER INIT MUST BE
C        SET TO 0.
C
C
C THE USER MUST MODIFY THE VALUE OF NDIM ACCORDING TO THE PROBLEM SIZE.
C NDIM IS THE DIMENSION OF THE ARRAY E.
C THE DIMENSIONS OF G ARE FIXED AS ( NDIM , NDIM + 1 ).
C NDIM MUST BE GREATER OR EQUAL TO MAXCOL. THE DIMENSION OF GINIT(*) IN
C THE CALLING PROGRAM MUST BE GREATER OR EQUAL TO 2*(MAXCOL + 1).
C
C THE CODE IS IN DOUBLE PRECISION.
C TO CHANGE IT TO SINGLE PRECISION THE USER MUST MODIFY THE 'DOUBLE
C PRECISION' TYPE STATEMENTS TO 'REAL' AND CHANGE THE INTRINSIC FUNCTION
C 'DABS()' TO 'ABS()'.
C THE USER MUST BE CONSCIOUS THAT, AS USUAL IN SEQUENCES TRANSFORMATIONS
C DIFFICULTIES OFTEN OCCUR DUE TO CANCELLATION ERRORS AND THAT DOUBLE
C PRECISION IS NEEDED IN MANY OF THE CASES.
C
      INTEGER INIT, MAXCOL, INFO
      DOUBLE PRECISION GINIT, RESULT, SK
C
      INTEGER I, IM1, J, K, KN, KPI, M, NDI
      DOUBLE PRECISION EPS, S, R, E, G
      DIMENSION E(25), G(25,26), GINIT(1)
C
C
      IF (INIT.NE.0) GO TO 10
C
C     FIRST CALL OF THE SUBROUTINE
C
C        NDIM MUST BE EQUAL TO THE DIMENSION OF THE VECTOR E
C
      NDIM = 25
C
      INFO = 0
      IF (MAXCOL.GE.NDIM) GO TO 110
      K = 0
      INIT = 1
      E(1) = SK
      RESULT = SK
      G(1,1) = GINIT(1)
      RETURN
C
C     FOLLOWING CALLS OF THE SUBROUTINE
C
   10 IF (INFO.NE.0) RETURN
      K = K + 1
      KN = MIN0(K,MAXCOL)
      IF (K.GE.NDIM) GO TO 100
      G(1,K+1) = GINIT(1)
      IF (K.LE.2) GO TO 40
C
C     COMPUTATION OF THE NEXT ROW IN THE ARRAYS G  FOR I=2 TO M
C                                                I
C
      M = MIN0(K-1,MAXCOL)
      DO 30 I=2,M
        IM1 = I - 1
        DO 20 J=1,IM1
          R = G(J,K+1) - G(J,K)
          IF (DABS(R).LT.EPS) GO TO 120
          R = GINIT(I) + (G(I,J)-GINIT(I))/R*G(J,K+1)
          G(I,J) = GINIT(I)
          GINIT(I) = R
   20   CONTINUE
        G(I,K+1) = R
   30 CONTINUE
      G(I,K+1) = R
C
C     COMPUTATION OF THE WHOLE ARRAY G
C                                     K
C
   40 IF (K.EQ.1 .OR. K.GT.MAXCOL) GO TO 80
      G(K,1) = GINIT(K)
      DO 60 I=1,K
        KPI = K + I
        M = MIN0(I,K-1)
C
C         COMPUTATION OF THE NEXT ROW OF THE ARRAY
C
        DO 50 J=1,M
          R = G(J,I+1) - G(J,I)
          IF (DABS(R).LT.EPS) GO TO 120
          R = GINIT(KPI) + (G(K,J)-GINIT(KPI))/R*G(J,I+1)
          G(K,J) = GINIT(KPI)
          GINIT(KPI) = R
   50   CONTINUE
        IF (I.EQ.K) GO TO 70
        G(K,M+1) = R
   60 CONTINUE
   70 G(K,K+1) = R
C
C     COMPUTATION OF THE NEXT ROW IN THE E-ARRAY
C
   80 S = SK
      DO 90 I=1,KN
        R = G(I,K+1) - G(I,K)
        IF (DABS(R).LT.EPS) GO TO 120
        R = S + (E(I)-S)/R*G(I,K+1)
        E(I) = S
        S = R
   90 CONTINUE
      E(KN+1) = R
      RESULT = R
      RETURN
C
C     ERROR CONDITIONS
C
  100 INFO = -1
      RETURN
  110 INFO = -2
      RETURN
  120 INFO = 1
      RETURN
      END
