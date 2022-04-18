C     ALGORITHM 432 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 09,
C     P. 820.
      SUBROUTINE AXPXB(A,U,M,NA,NU,B,V,N,NB,NV,C,NC,EPSA,               AXPX  10
     1EPSB,FAIL)
C
C AXPXB IS A FORTRAN IV SUBROUTINE TO SOLVE THE REAL MATRIX
C EQUATION AX + XB - C.  THE MATRICES A AND B ARE TRANS-
C FORMED INTO REAL SCHUR FORM, AND THE TRANSFORMED SYSTEM IS
C SOLVED BY BACK SUBSTITUTION.  THE PROGRAM REQUIRES THE
C AUXILIARY SUBROUTINES HSHLDR, BCKMLT, SCHUR, AND SHRSLV.
C THE PARAMETERS IN THE ARGUMENT LIST ARE
C
C           A      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX A.  ON RETURN, THE LOWER TRIANGLE
C                  AND SUPERDIAGONAL OF THE ARRAY A CONTAIN
C                  A LOWER REAL SCHUR FORM OF A.  THE ARRAY
C                  A MUST BE DIMENSIONED AT LEAST M+1 BY
C                  M+1.
C           U      A DOUBLY SUBSCRIPTED ARRAY THAT, ON
C                  RETURN, CONTAINS THE ORTHOGONAL MATRIX
C                  THAT REDUCES A TO REAL SCHUR FORM.
C           M      THE ORDER OF THE MATRIX A.
C           NA     THE FIRST DIMENSION OF THE ARRAY A.
C           NU     THE FIRST DIMENSION OF THE ARRAY U.
C           B      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX B.  ON RETURN, THE UPPER TRIANGLE
C                  AND SUBDIAGONAL OF THE ARRAY B CONTAIN AN
C                  UPPER REAL SCHUR FORM OF B.  THE ARRAY B
C                  MUST BE DIMENSIONED AT LEAST M+1 BY M+1.
C           V      A DOUBLY SUBSCRIPTED ARRAY THAT, ON
C                  RETURN, CONTAINS THE ORTHOGONAL MATRIX
C                  THAT REDUCES B TO REAL SCHUR FORM.
C           N      THE ORDER OF THE MATRIX B.
C           NB     THE FIRST DIMENSION OF THE ARRAY B.
C           NV     THE FIRST DIMENSION OF THE ARRAY V.
C           C      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX C.  ON RETURN, C CONTAINS THE
C                  SOLUTION MATRIX X.
C           NC     THE FIRST DIMENSION OF THE ARRAY C.
C           EPSA   A CONVERGENCE CRITERION FOR THE REDUCTION
C                  OF A TO SCHUR FORM.  EPSA SHOULD BE SET
C                  SLIGHTLY SMALLER THAN 10.**(-N), WHERE N
C                  IS THE NUMBER OF SIGNIFICANT DIGITS IN
C                  THE ELEMENTS OF THE MATRIX A.
C           EPSB   A CONVERGENCE CRITERION FOR THE REDUCTION
C                  OF B TO REAL SCHUR FORM.
C           FAIL   AN INTEGER VARIABLE THAT, ON RETURN,
C                  CONTAINS AN ERROR SIGNAL.  IF FAIL IS
C                  POSITIVE (NEGATIVE) THEN THE PROGRAM WAS
C                  UNABLE TO REDUCE A (B) TO REAL SCHUR
C                  FORM.  IF FAIL IS ZERO, THE REDUCTIONS
C                  PROCEEDED WITHOUT MISHAP.
C
C WHEN EPSA IS NEGATIVE THE REDUCTION OF A TO REAL SCHUR
C FORM IS SKIPPED AND THE ARRAYS A AND U ARE ASSUMED TO
C CONTAIN THE SCHUR FORM AND ACCOMPANYING ORTHOGONAL MATRIX.
C THIS PERMITS THE EFFICIENT SOLUTION OF SEVERAL EQUATIONS
C OF THE FORM AX + BX = C WHEN A DOES NOT CHANGE.  LIKEWISE,
C IF EPSB IS NEGATIVE, THE REDUCTION OF B TO REAL SCHUR FORM
C IS SKIPPED.
C
      REAL
     1A(NA,1),U(NU,1),B(NB,1),V(NV,1),C(NC,1),EPSA,EPSB,TEMP
      INTEGER
     1M,NA,NU,N,NB,NV,NC,FAIL,M1,MM1,N1,NM1,I,J,K
      M1 = M+1
      MM1 = M-1
      N1 = N+1
      NM1 = N-1
C
C IF REQUIRED, REDUCE A TO LOWER REAL SCHUR FORM.
C
      IF(EPSA .LT. 0.) GO TO 35
      DO 10 I=1,M
        DO 10 J=I,M
          TEMP = A(I,J)
          A(I,J) = A(J,I)
          A(J,I) = TEMP
   10 CONTINUE
      CALL HSHLDR(A,M,NA)
      CALL BCKMLT(A,U,M,NA,NU)
      IF(MM1 .EQ. 0) GO TO 25
      DO 20 I=1,MM1
        A(I+1,I) = A(I,M1)
   20 CONTINUE
      CALL SCHUR(A,U,M,NA,NU,EPSA,FAIL)
      IF(FAIL .NE. 0) RETURN
   25 DO 30 I=1,M
        DO 30 J=I,M
          TEMP = A(I,J)
          A(I,J) = A(J,I)
          A(J,I) = TEMP
   30 CONTINUE
C
C IF REQUIRED, REDUCE B TO UPPER REAL SCHUR FORM.
C
   35 IF(EPSB .LT. 0.) GO TO 45
      CALL HSHLDR(B,N,NB)
      CALL BCKMLT(B,V,N,NB,NV)
      IF(NM1 .EQ. 0) GO TO 45
      DO 40 I=1,NM1
        B(I+1,I) = B(I,N1)
   40 CONTINUE
      CALL SCHUR(B,V,N,NB,NV,EPSB,FAIL)
      FAIL = -FAIL
      IF(FAIL .NE. 0) RETURN
C
C TRANSFORM C.
C
   45 DO 60 J=1,N
        DO 50 I=1,M
          A(I,M1) = 0.
          DO 50 K=1,M
            A(I,M1) = A(I,M1) + U(K,I)*C(K,J)
   50   CONTINUE
        DO 60 I=1,M
          C(I,J) = A(I,M1)
   60 CONTINUE
      DO 80 I=1,M
        DO 70 J=1,N
          B(N1,J) = 0.
          DO 70 K=1,N
            B(N1,J) = B(N1,J) + C(I,K)*V(K,J)
   70   CONTINUE
        DO 80 J=1,N
          C(I,J) = B(N1,J)
   80 CONTINUE
C
C SOLVE THE TRANSFORMED SYSTEM.
C
      CALL SHRSLV(A,B,C,M,N,NA,NB,NC)
C
C TRANSFORM C BACK TO THE SOLUTION.
C
      DO 100 J=1,N
        DO 90 I=1,M
          A(I,M1) = 0.
          DO 90 K=1,M
            A(I,M1) = A(I,M1) + U(I,K)*C(K,J)
   90   CONTINUE
        DO 100 I=1,M
          C(I,J) = A(I,M1)
  100 CONTINUE
      DO 120 I=1,M
        DO 110 J=1,N
          B(N1,J) = 0.
          DO 110 K=1,N
            B(N1,J) = B(N1,J) + C(I,K)*V(J,K)
  110   CONTINUE
        DO 120 J=1,N
          C(I,J) = B(N1,J)
  120 CONTINUE
      RETURN
      END
      SUBROUTINE SHRSLV(A,B,C,M,N,NA,NB,NC)                             SHRS1530
C SHRSLV IS A FORTRAN IV SUBROUTINE TO SOLVE THE REAL MATRIX
C EQATION AX + XB = C, WHERE A IS IN LOWER REAL SCHUR FORM
C AND B IS IN UPPER REAL SCHUR FORM.  SHRSLV USES THE AUX-
C ILIARY SUBROUTINE SYSSLV, WHICH IT COMMUNICATES WITH
C THROUGH THE COMMON BLOCK SLVBLK.  THE PARAMETERS IN THE
C ARGUMENT LIST ARE
C           A      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX A IN LOWER REAL SCHUR FORM.
C           B      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX B IN UPPER REAL SCHUR FORM.
C           C      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX C.
C           M      THE ORDER OF THE MATRIX A.
C           N      THE ORDER OF THE MATRIX B.
C           NA     THE FIRST DIMENSION OF THE ARRAY A.
C           NB     THE FIRST DIMENSION OF THE ARRAY B.
C           NC     THE FIRST DIMENSION OF THE ARRAY C.
C
      REAL
     1A(NA,1),B(NB,1),C(NC,1),T,P
      INTEGER
     1M,N,NA,NB,NC,K,KM1,DK,KK,L,LM1,DL,LL,I,IB,J,JA,NSYS
      COMMON/SLVBLK/T(5,5),P(5),NSYS
      L = 1
   10   LM1 = L-1
        DL = 1
        IF(L .EQ. N) GO TO 15
        IF(B(L+1,L) .NE. 0.) DL = 2
   15   LL = L+DL-1
        IF(L .EQ. 1) GO TO 30
        DO 20 J=L,LL
          DO 20 I=1,M
            DO 20 IB=1,LM1
              C(I,J) = C(I,J) - C(I,IB)*B(IB,J)
   20   CONTINUE
   30   K = 1
   40     KM1 = K-1
          DK = 1
          IF(K .EQ. M) GO TO 45
          IF(A(K,K+1) .NE. 0.) DK = 2
   45     KK = K+DK-1
          IF(K .EQ. 1) GO TO 60
          DO 50 I=K,KK
            DO 50 J=L,LL
              DO 50 JA=1,KM1
                C(I,J) = C(I,J) - A(I,JA)*C(JA,J)
   50     CONTINUE
   60     IF(DL .EQ. 2) GO TO 80
          IF(DK .EQ. 2) GO TO 70
          T(1,1) = A(K,K) + B(L,L)
          IF(T(1,1) .EQ. 0.) STOP
          C(K,L) = C(K,L)/T(1,1)
          GO TO 100
   70     T(1,1) = A(K,K) + B(L,L)
          T(1,2) = A(K,KK)
          T(2,1) = A(KK,K)
          T(2,2) = A(KK,KK) + B(L,L)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          GO TO 100
   80     IF(DK .EQ. 2) GO TO 90
          T(1,1) = A(K,K) + B(L,L)
          T(1,2) = B(LL,L)
          T(2,1) = B(L,LL)
          T(2,2) = A(K,K) + B(LL,LL)
          P(1) = C(K,L)
          P(2) = C(K,LL)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(K,LL) = P(2)
          GO TO 100
   90     T(1,1) = A(K,K) + B(L,L)
          T(1,2) = A(K,KK)
          T(1,3) = B(LL,L)
          T(1,4) = 0.
          T(2,1) = A(KK,K)
          T(2,2) = A(KK,KK) + B(L,L)
          T(2,3) = 0.
          T(2,4) = T(1,3)
          T(3,1) = B(L,LL)
          T(3,2) = 0.
          T(3,3) = A(K,K) + B(LL,LL)
          T(3,4) = T(1,2)
          T(4,1) = 0.
          T(4,2) = T(3,1)
          T(4,3) = T(2,1)
          T(4,4) = A(KK,KK) + B(LL,LL)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          P(3) = C(K,LL)
          P(4) = C(KK,LL)
          NSYS = 4
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          C(K,LL) = P(3)
          C(KK,LL) = P(4)
  100   K = K + DK
        IF(K .LE. M) GO TO 40
      L = L + DL
      IF(L .LE. N) GO TO 10
      RETURN
      END
      SUBROUTINE ATXPXA(A,U,C,N,NA,NU,NC,EPS,FAIL)                      ATXP2620
C
C ATXPXA IS A FORTRAN IV SUBROUTINE TO SOLVE THE REAL MATRIX
C EQUATION TRANS(A)*X + X*A = C, WHERE C IS SYMMETRIC AND
C TRANS(A) DENOTES THE TRANSPOSE OF A.  THE EQUATION IS
C TRANSFORMED SO THAT A IS IN UPPER REAL SCHUR FORM, AND THE
C TRANSFORMED EQUATION IS SOLVED BY A RECURSIVE PROCEDURE.
C THE PROGRAM REQUIRES THE AUXILIARY SUBROUTINES HSHLDR,
C BCKMLT, SCHUR, AND SYMSLV.  THE PARAMETERS IN THE ARGUMENT
C LIST ARE
C           A      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX A.  ON RETURN, THE UPPER TRIANGLE
C                  AND THE FIRST SUBDIAGONAL OF THE ARRAY A
C                  CONTAIN AN UPPER REAL SCHUR FORM OF A.
C                  THE ARRAY A MUST BE DIMENSIONED AT LEAST
C                  N+1 BY N+1.
C           U      A DOUBLY SUBSCRIPTED ARRAY THAT, ON
C                  RETURN, CONTAINS THE ORTHOGONAL MATRIX
C                  THAT REDUCES A TO UPPER REAL SCHUR FORM.
C           C      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX C.  ON RETURN, C CONTAINS THE
C                  SOLUTION MATRIX X.
C           N      THE ORDER OF THE MATRIX A.
C           NA     THE FIRST DIMENSION OF THE ARRAY A.
C           NU     THE FIRST DIMENSION OF THE ARRAY U.
C           NC     THE FIRST DIMENSION OF THE ARRAY C.
C           EPS    A CONVERGENCE CRITERION FOR THE REDUCTION
C                  OF A TO REAL SCHUR FORM.  EPS SHOULD BE
C                  SET SLIGHTLY SMALLER THAN 10.**(-N),
C                  WHERE N IS THE NUMBER OF SIGNIFICANT
C                  DIGITS IN THE ELEMENTS OF THE MATRIX A.
C           FAIL   AN INTEGER VARIABLE THAT, ON RETURN,
C                  CONTAINS AN ERROR SIGNAL.  IF FAIL IS
C                  NONZERO, THEN THE PROGRAM WAS UNABLE TO
C                  REDUCE A TO REAL SCHUR FORM.  IF FAIL IS
C                  ZERO, THE REDUCTION PROCEEDED WITHOUT
C                  MISHAP.
C
C WHEN EPS IS NEGATIVE, THE REDUCTION OF A TO REAL SCHUR
C FORM IS SKIPPED AND THE ARRAYS A AND U ARE ASSUMED TO
C CONTAIN THE SCHUR FORM AND ACCOMPANYING ORTHOGONAL MATRIX.
C THIS PERMITS THE EFFICIENT SOLUTION OF SEVERAL EQUATIONS
C WITH DIFFERENT RIGHT HAND SIDES.
C
      REAL
     1A(NA,1),U(NU,1),C(NC,1),EPS
      INTEGER
     1N,NA,NU,NC,FAIL,N1,NM1,I,J,K
      N1 = N+1
      NM1 = N-1
C
C IF REQUIRED, REDUCE A TO UPPER REAL SCHUR FORM.
C
      IF(EPS .LT. 0.) GO TO 15
      CALL HSHLDR(A,N,NA)
      CALL BCKMLT(A,U,N,NA,NU)
      DO 10 I=1,NM1
        A(I+1,I) = A(I,N1)
   10 CONTINUE
      CALL SCHUR(A,U,N,NA,NU,EPS,FAIL)
      IF(FAIL .NE. 0) RETURN
C
C TRANSFORM C.
C
   15 DO 20 I=1,N
        C(I,I) = C(I,I)/2.
   20 CONTINUE
      DO 40 I=1,N
        DO 30 J=1,N
          A(N1,J) = 0.
          DO 30 K=I,N
            A(N1,J) = A(N1,J) + C(I,K)*U(K,J)
   30   CONTINUE
        DO 40 J=1,N
          C(I,J) = A(N1,J)
   40 CONTINUE
      DO 60 J=1,N
        DO 50 I=1,N
          A(I,N1) = 0.
          DO 50 K=1,N
            A(I,N1) = A(I,N1) + U(K,I)*C(K,J)
   50   CONTINUE
        DO 60 I=1,N
          C(I,J) = A(I,N1)
   60 CONTINUE
      DO 70 I=1,N
        DO 70 J=I,N
          C(I,J) = C(I,J) + C(J,I)
          C(J,I) = C(I,J)
   70 CONTINUE
C
C SOLVE THE TRANSFORMED SYSTEM.
C
      CALL SYMSLV(A,C,N,NA,NC)
C
C TRANSFORM C BACK TO THE SOLUTION.
C
      DO 80 I=1,N
        C(I,I) = C(I,I)/2.
   80 CONTINUE
      DO 100 I=1,N
        DO 90 J=1,N
          A(N1,J) = 0.
          DO 90 K=I,N
            A(N1,J) = A(N1,J) + C(I,K)*U(J,K)
   90   CONTINUE
        DO 100 J=1,N
          C(I,J) = A(N1,J)
  100 CONTINUE
      DO 120 J=1,N
        DO 110 I=1,N
          A(I,N1) = 0.
          DO 110 K=1,N
            A(I,N1) = A(I,N1) + U(I,K)*C(K,J)
  110   CONTINUE
        DO 120 I=1,N
          C(I,J) = A(I,N1)
  120 CONTINUE
      DO 130 I=1,N
        DO 130 J=I,N
          C(I,J) = C(I,J) + C(J,I)
          C(J,I) = C(I,J)
  130 CONTINUE
      RETURN
      END
      SUBROUTINE SYMSLV(A,C,N,NA,NC)                                    SYMS3870
C
C SYMSLV IS A FORTRAN IV SUBROUTINE TO SOLVE THE REAL MATRIX
C EQUATION TRANS(A)*X + X*A = C, WHERE C IS SYMMETRIC, A IS
C IN UPPER REAL SCHUR FORM, AND TRANS(A) DENOTES THE TRANS-
C POSE OF A.  SYMSLV USES THE AUXILIARY SUBROUTINE SYSSLV,
C WHICH IT COMMUNICATES WITH THROUGH THE COMMON BLOCK
C SLVBLK.  THE PARAMETERS IN THE ARGUMENT LIST ARE
C           A      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX A IN UPPER REAL SCHUR FORM.
C           C      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX C.
C           N      THE ORDER OF THE MATRIX A.
C           NA     THE FIRST DIMENSION OF THE ARRAY A.
C           NC     THE FIRST DIMENSION OF THE ARRAY C.
C
      REAL
     1A(NA,1),C(NC,1),T,P
      INTEGER
     1N,NA,NC,K,KK,DK,KM1,L,LL,DL,LDL,I,IA,J,NSYS
      COMMON/SLVBLK/T(5,5),P(5),NSYS
      L = 1
   10   DL = 1
        IF(L .EQ. N) GO TO 20
        IF(A(L+1,L) .NE. 0.) DL = 2
   20   LL = L+DL-1
        K = L
   30     KM1 = K-1
          DK = 1
          IF(K .EQ. N) GO TO 35
          IF(A(K+1,K) .NE. 0.) DK = 2
   35     KK = K+DK-1
          IF(K .EQ. L) GO TO 45
          DO 40 I=K,KK
            DO 40 J=L,LL
              DO 40 IA=L,KM1
                C(I,J) = C(I,J) - A(IA,I)*C(IA,J)
   40     CONTINUE
   45     IF(DL .EQ. 2) GO TO 60
          IF(DK .EQ. 2 ) GO TO 50
          T(1,1) = A(K,K) + A(L,L)
          IF(T(1,1) .EQ. 0.) STOP
          C(K,L) = C(K,L)/T(1,1)
          GO TO 90
   50     T(1,1) = A(K,K) + A(L,L)
          T(1,2) = A(KK,K)
          T(2,1) = A(K,KK)
          T(2,2) = A(KK,KK) + A(L,L)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          GO TO 90
   60     IF(DK .EQ. 2) GO TO 70
          T(1,1) = A(K,K) + A(L,L)
          T(1,2) = A(LL,L)
          T(2,1) = A(L,LL)
          T(2,2) = A(K,K) + A(LL,LL)
          P(1) = C(K,L)
          P(2) = C(K,LL)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(K,LL) = P(2)
          GO TO 90
   70     IF(K .NE. L) GO TO 80
          T(1,1) = A(L,L)
          T(1,2) = A(LL,L)
          T(1,3) = 0.
          T(2,1) = A(L,LL)
          T(2,2) = A(L,L) + A(LL,LL)
          T(2,3) = T(1,2)
          T(3,1) = 0.
          T(3,2) = T(2,1)
          T(3,3) = A(LL,LL)
          P(1) = C(L,L)/2.
          P(2) = C(LL,L)
          P(3) = C(LL,LL)/2.
          NSYS = 3
          CALL SYSSLV
          C(L,L) = P(1)
          C(LL,L) = P(2)
          C(L,LL) = P(2)
          C(LL,LL) = P(3)
          GO TO 90
   80     T(1,1) = A(K,K) + A(L,L)
          T(1,2) = A(KK,K)
          T(1,3) = A(LL,L)
          T(1,4) = 0.
          T(2,1) = A(K,KK)
          T(2,2) = A(KK,KK) + A(L,L)
          T(2,3) = 0.
          T(2,4) = T(1,3)
          T(3,1) = A(L,LL)
          T(3,2) = 0.
          T(3,3) = A(K,K) + A(LL,LL)
          T(3,4) = T(1,2)
          T(4,1) = 0.
          T(4,2) = T(3,1)
          T(4,3) = T(2,1)
          T(4,4) = A(KK,KK) + A(LL,LL)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          P(3) = C(K,LL)
          P(4) = C(KK,LL)
          NSYS = 4
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          C(K,LL) = P(3)
          C(KK,LL) = P(4)
   90   K = K + DK
        IF(K .LE. N) GO TO 30
        LDL = L + DL
        IF(LDL .GT. N) RETURN
        DO 120 J=LDL,N
          DO 100 I=L,LL
            C(I,J) = C(J,I)
  100     CONTINUE
          DO 120 I=J,N
            DO 110 K=L,LL
              C(I,J) = C(I,J) - C(I,K)*A(K,J) - A(K,I)*C(K,J)
  110       CONTINUE
            C(J,I) = C(I,J)
  120   CONTINUE
      L = LDL
      GO TO 10
      END
      SUBROUTINE HSHLDR(A,N,NA)                                         HSHL5170
C
C HSHLDR IS A FORTRAN IV SUBROUTINE TO REDUCE A MATRIX TO
C UPPER HESSENBERG FORM BY ELEMENTARY HERMITIAN TRANSFORMA-
C TIONS (THE METHOD OF HOUSEHOLDER).  THE PARAMETERS IN THE
C ARGUMENT LIST ARE
C           A      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  MATRIX A.  ON RETURN, THE UPPER TRIANGLE
C                  OF THE ARRAY A MATRIX AND THE (N+1)-TH
C                  COLUMN CONTAIN THE SUBDIAGONAL ELEMENTS
C                  OF THE TRANSFORMED MATRIX.  ON RETURN,
C                  THE LOWER TRIANGLE AND THE (N+1)-TH ROW
C                  OF THE ARRAY A CONTAIN A HISTORY OF THE
C                  TRANSFORMATIONS.
C           N      THE ORDER OF THE MATRIX A.
C           NA     THE FIRST DIMENSION OF THE ARRAY A.
C
      REAL
     1A(NA,1),MAX,SUM,S,P
      INTEGER
     1N,NA,NM2,N1,L,L1,I,J
      NM2 = N-2
      N1 = N+1
      IF(N .EQ. 1) RETURN
      IF(N .GT. 2) GO TO 5
      A(1,N1) = A(2,1)
      RETURN
    5 DO 80 L=1,NM2
        L1 = L+1
        MAX = 0.
        DO 10 I=L1,N
          MAX = AMAX1(MAX,ABS(A(I,L)))
   10   CONTINUE
        IF(MAX .NE. 0.) GO TO 20
        A(L,N1) = 0.
        A(N1,L) = 0.
        GO TO 80
   20   SUM = 0.
        DO 30 I=L1,N
          A(I,L) = A(I,L)/MAX
          SUM = SUM + A(I,L)**2
   30   CONTINUE
        S = SIGN(SQRT(SUM),A(L1,L))
        A(L,N1) = -MAX*S
        A(L1,L) = S + A(L1,L)
        A(N1,L) = S*A(L1,L)
        DO 50 J=L1,N
          SUM = 0.
          DO 40 I=L1,N
            SUM = SUM + A(I,L)*A(I,J)
   40     CONTINUE
          P = SUM/A(N1,L)
          DO 50 I=L1,N
            A(I,J) = A(I,J) - A(I,L)*P
   50   CONTINUE
        DO 70 I=1,N
          SUM = 0.
          DO 60 J=L1,N
            SUM = SUM + A(I,J)*A(J,L)
   60     CONTINUE
          P = SUM/A(N1,L)
          DO 70 J=L1,N
            A(I,J) = A(I,J) - P*A(J,L)
   70   CONTINUE
   80 CONTINUE
      A(N-1,N1) = A(N,N-1)
      RETURN
      END
      SUBROUTINE BCKMLT(A,U,N,NA,NU)                                    BCKM5850
C
C BCKMLT IS A FORTRAN IV SUBROUTINE THAT, GIVEN THE OUTPUT
C OF THE SUBROUTINE HSHLDR, COMPUTES THE ORTHOGONAL MATRIX
C THAT REDUCES A TO UPPER HESSENBERG FORM.  THE PARAMETERS
C IN THE ARGUMENT LIST ARE
C           A      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  OUTPUT FROM HSHLDR.
C           U      A DOUBLY SUBSCRIPTED ARRAY THAT, ON
C                  RETURN, CONTAINS THE ORTHOGONAL MATRIX.
C           N      THE ORDER OF THE MATRIX A IN HSHLDR.
C           NA     THE FIRST DIMENSION OF THE ARRAY A.
C           NU     THE FIRST DIMENSION OF THE ARRAY U.
C
C THE ARRAYS A AND U MAY BE IDENTIFIED IN THE CALLING
C SEQUENCE.  IF THIS IS DONE, THE ELEMENTS OF THE ORTHOGONAL
C MATRIX WILL OVERWRITE THE OUTPUT OF HSHLDR.
C
      REAL
     1A(NA,1),U(NU,1),SUM,P
      INTEGER
     1N,NA,N1,NM1,NM2,LL,L,L1,I,J
      N1 = N+1
      NM1 = N-1
      NM2 = N-2
      U(N,N) = 1.
      IF(NM1 .EQ. 0) RETURN
      U(NM1,N) = 0.
      U(N,NM1) = 0.
      U(NM1,NM1) = 1.
      IF(NM2 .EQ. 0) RETURN
      DO 40 LL=1,NM2
        L = NM2-LL+1
        L1 = L+1
        IF(A(N1,L) .EQ. 0.) GO TO 25
        DO 20 J=L1,N
          SUM = 0.
          DO 10 I=L1,N
            SUM = SUM + A(I,L)*U(I,J)
   10     CONTINUE
          P = SUM/A(N1,L)
          DO 20 I=L1,N
            U(I,J) = U(I,J) - A(I,L)*P
   20   CONTINUE
   25   DO 30 I=L1,N
           U(I,L) = 0.
           U(L,I) = 0.
   30   CONTINUE
        U(L,L) = 1.
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SCHUR(H,U,NN,NH,NU,EPS,FAIL)                           SCHU6370
C
C SCHUR IS A FORTRAN IV SUBROUTINE TO REDUCE AN UPPER
C HESSENBERG MATRIX TO REAL SCHUR FORM BY THE QR METHOD WITH
C IMPLICIT ORIGIN SHIFTS.  THE PRODUCT OF THE TRANSFORMA-
C TIONS USED IN THE REDUCTION IS ACCUMULATED.  SCHUR IS AN
C ADAPTATION OF THE ALGOL PROGRAM HQR BY MARTIN, PETERS, AND
C WILKINSON (NUMER. MATH. 14 (1970) 219-231).  THE PARA-
C METERS IN THE ARGUMENT LIST ARE
C           H      A DOUBLY SUBSCRIPTED ARRAY CONTAINING THE
C                  UPPER HESSENBERG MATRIX H.  ON RETURN, H
C                  CONTAINS AN UPPER REAL SCHUR FORM OF H.
C                  THE ELEMENTS OF THE ARRAY H BELOW THE
C                  THIRD SUBDIAGONAL ARE UNDISTURBED.
C           U      A DOUBLY SUBSCRIPTED ARRAY CONTAINING ANY
C                  MATRIX.  ON RETURN, U CONTAINS THE MATRIX
C                  U*R(1)*R(2)..., WHERE R(I) ARE THE TRANS-
C                  FORMATIONS USED IN THE REDUCTION OF H.
C           NN     THE ORDER OF THE MATRICES H AND U.
C           NH     THE FIRST DIMENSION OF THE ARRAY H.
C           NU     THE FIRST DIMENSION OF THE ARRAY U.
C           EPS    A NUMBER USED IN DETERMINING WHEN AN
C                  ELEMENT OF H IS NEGLIGIBLE.  H(I,J) IS
C                  NEGLIGIBLE IF ABS(H(I,J)) IS LESS THAN OR
C                  EQUAL TO EPS TIMES THE INFINITY NORM OF
C                  H.
C           FAIL   AN INTEGER VARIABLE THAT, ON RETURN,
C                  CONTAINS AN ERROR SIGNAL. IF FAIL IS
C                  POSITIVE, THEN THE PROGRAM FAILED TO MAKE
C                  THE FAIL-1 OR FAIL-2 SUBDIAGONAL ELEMENT
C                  NEGLIGIBLE AFTER 30 ITERATIONS.
C
      REAL
     1H(NH,1),U(NU,1),EPS,HN,RSUM,TEST,P,Q,R,S,W,X,Y,Z
      INTEGER
     1NN,NA,NH,FAIL,I,ITS,J,JL,K,L,LL,M,MM,M2,M3,N,NA
      LOGICAL
     1LAST
      N = NN
      HN = 0.
      DO 20 I=1,N
        JL = MAX0(1,I-1)
        RSUM = 0.
        DO 10 J=JL,N
          RSUM = RSUM + ABS(H(I,J))
   10   CONTINUE
        HN = AMAX1(HN,RSUM)
   20 CONTINUE
      TEST = EPS*HN
      IF(HN .EQ. 0.) GO TO 230
   30 IF(N .LE. 1) GO TO 230
      ITS = 0
      NA = N-1
      NM2 = N-2
   40 DO 50 LL=2,N
      L = N-LL+2
        IF(ABS(H(L,L-1)) .LE. TEST) GO TO 60
   50 CONTINUE
      L = 1
      GO TO 70
   60 H(L,L-1) = 0.
   70 IF(L .LT. NA) GO TO 72
      N = L-1
      GO TO 30
   72 X = H(N,N)/HN
      Y = H(NA,NA)/HN
      R = (H(N,NA)/HN)*(H(NA,N)/HN)
      IF(ITS .LT. 30) GO TO 75
      FAIL = N
      RETURN
   75 IF(ITS.EQ.10 .OR. ITS.EQ.20) GO TO 80
      S = X + Y
      Y = X*Y - R
      GO TO 90
   80 Y = (ABS(H(N,NA)) + ABS(H(NA,NM2)))/HN
      S = 1.5*Y
      Y = Y**2
   90 ITS = ITS + 1
      DO 100 MM=L,NM2
        M = NM2-MM+L
        X = H(M,M)/HN
        R = H(M+1,M)/HN
        Z = H(M+1,M+1)/HN
        P = X*(X-S) + Y + R*(H(M,M+1)/HN)
        Q = R*(X+Z-S)
        R = R*(H(M+2,M+1)/HN)
        W = ABS(P) + ABS(Q) + ABS(R)
        P = P/W
        Q = Q/W
        R = R/W
        IF(M .EQ. L) GO TO 110
      IF(ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. ABS(P)*TEST)
     1GO TO 110
  100 CONTINUE
  110 M2 = M+2
      M3 = M+3
      DO 120 I=M2,N
        H(I,I-2) = 0.
  120 CONTINUE
      IF(M3 .GT. N) GO TO 140
      DO 130 I=M3,N
        H(I,I-3) = 0.
  130 CONTINUE
  140 DO 220 K=M,NA
        LAST = K.EQ.NA
        IF(K .EQ. M) GO TO 150
        P = H(K,K-1)
        Q = H(K+1,K-1)
        R = 0.
        IF(.NOT.LAST) R = H(K+2,K-1)
        X = ABS(P) + ABS(Q) + ABS(R)
        IF(X .EQ. 0.) GO TO 220
        P = P/X
        Q = Q/X
        R = R/X
  150   S = SQRT(P**2 + Q**2 + R**2)
        IF(P .LT. 0.) S = -S
        IF(K .NE. M) H(K,K-1) = -S*X
        IF(K.EQ.M .AND. L.NE.M) H(K,K-1) = -H(K,K-1)
        P = P + S
        X = P/S
        Y = Q/S
        Z = R/S
        Q = Q/P
        R = R/P
        DO 170 J=K,NN
          P = H(K,J) + Q*H(K+1,J)
          IF(LAST) GO TO 160
          P = P + R*H(K+2,J)
          H(K+2,J) = H(K+2,J) - P*Z
  160     H(K+1,J) = H(K+1,J) - P*Y
          H(K,J) = H(K,J) - P*X
  170   CONTINUE
        J = MIN0(K+3,N)
        DO 190 I=1,J
          P = X*H(I,K) + Y*H(I,K+1)
          IF(LAST) GO TO 180
          P = P + Z*H(I,K+2)
          H(I,K+2) = H(I,K+2) - P*R
  180     H(I,K+1) = H(I,K+1) - P*Q
          H(I,K) = H(I,K) - P
  190   CONTINUE
        DO 210 I=1,NN
          P = X*U(I,K) + Y*U(I,K+1)
          IF(LAST) GO TO 200
          P = P + Z*U(I,K+2)
          U(I,K+2) = U(I,K+2) - P*R
  200     U(I,K+1) = U(I,K+1) - P*Q
          U(I,K) = U(I,K) - P
  210   CONTINUE
  220 CONTINUE
      GO TO 40
  230 FAIL = 0
      RETURN
      END
      SUBROUTINE SYSSLV                                                 SYSS7920
C
C SYSSLV IS A FORTRAN IV SUBROUTINE THAT SOLVES THE LINEAR
C SYSTEM AX = B OF ORDER N LESS THAN 5 BY CROUT REDUCTION
C FOLLOWED BY BACK SUBSTITUTION.  THE MATRIX A, THE VECTOR
C B, AND THE ORDER N ARE CONTAINED IN THE ARRAYS A,B, AND
C THE VARIABLE N OF THE COMMON BLOCK SLVBLK.  THE SOLUTION
C IS RETURNED IN THE ARRAY B.
C
      COMMON/SLVBLK/A(5,5),B(5),N
      REAL MAX
    1 NM1 = N-1
      N1 = N+1
C
C COMPUTE THE LU FACTORIZATION OF A.
C
      DO 80 K=1,N
        KM1 = K-1
        IF(K.EQ.1) GO TO 20
        DO 10 I=K,N
          DO 10 J=1,KM1
            A(I,K) = A(I,K) - A(I,J)*A(J,K)
   10   CONTINUE
   20   IF(K.EQ.N) GO TO 100
        KP1 = K+1
      MAX = ABS(A(K,K))
        INTR = K
        DO 30 I=KP1,N
          AA = ABS(A(I,K))
          IF(AA .LE. MAX) GO TO 30
          MAX = AA
          INTR = I
   30   CONTINUE
        IF(MAX .EQ. 0.) STOP
        A(N1,K) = INTR
        IF(INTR .EQ. K) GO TO 50
        DO 40 J=1,N
          TEMP = A(K,J)
          A(K,J) = A(INTR,J)
          A(INTR,J) = TEMP
   40   CONTINUE
   50   DO 80 J=KP1,N
          IF(K.EQ.1) GO TO 70
          DO 60 I=1,KM1
            A(K,J) = A(K,J) - A(K,I)*A(I,J)
   60     CONTINUE
   70     A(K,J) = A(K,J)/A(K,K)
   80 CONTINUE
C
C INTERCHANGE THE COMPONENTS OF B.
C
  100 DO 110 J=1,NM1
        INTR = A(N1,J)
        IF(INTR .EQ. J) GO TO 110
        TEMP = B(J)
        B(J) = B(INTR)
        B(INTR) = TEMP
  110 CONTINUE
C
C SOLVE LX = B.
C
  200 B(1) = B(1)/A(1,1)
      DO 220 I=2,N
        IM1 = I-1
        DO 210 J=1,IM1
          B(I) = B(I) - A(I,J)*B(J)
  210   CONTINUE
        B(I) = B(I)/A(I,I)
  220 CONTINUE
C
C SOLVE UX = B.
C
  300 DO 310 II=1,NM1
        I = NM1-II+1
        I1 = I+1
        DO 310 J=I1,N
           B(I) = B(I) - A(I,J)*B(J)
  310 CONTINUE
      RETURN
      END
