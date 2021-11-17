      SUBROUTINE DONEST (N, V, X, ISGN, EST, KASE)
      INTEGER N, ISGN(N), KASE
      DOUBLE PRECISION V(N), X(N), EST

C
C     DONEST ESTIMATES THE 1-NORM OF A SQUARE, DOUBLE PRECISION MATRIX  A.
C     REVERSE COMMUNICATION IS USED FOR EVALUATING
C     MATRIX-VECTOR PRODUCTS.
C
C     ON ENTRY
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX.  N .GE. 1.
C
C        ISGN    INTEGER(N)
C                USED AS WORKSPACE.
C
C        KASE    INTEGER
C                = 0.
C
C     ON INTERMEDIATE RETURNS
C
C        KASE    = 1 OR 2.
C
C        X       DOUBLE PRECISION(N)
C                MUST BE OVERWRITTEN BY
C
C                     A*X,             IF KASE=1,
C                     TRANSPOSE(A)*X,  IF KASE=2,
C
C                AND DONEST MUST BE RE-CALLED, WITH ALL THE OTHER
C                PARAMETERS UNCHANGED.
C
C     ON FINAL RETURN
C
C        KASE    = 0.
C
C        EST     DOUBLE PRECISION
C                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C        V       DOUBLE PRECISION(N)
C                = A*W,   WHERE  EST = NORM(V)/NORM(W)
C                         (W  IS NOT RETURNED).
C
C     THIS VERSION DATED MARCH 16, 1988.
C     NICK HIGHAM, UNIVERSITY OF MANCHESTER.
C
C     MODIFIED FOR DOUBLE PRECISION ON JUNE 11, 1996.
C
C     REFERENCE
C     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C     SUBROUTINES AND FUNCTIONS
C     BLAS     IDAMAX, DASUM, DCOPY
C     GENERIC  ABS, NINT, DFLOAT, SIGN
C
        INTRINSIC DFLOAT
        DOUBLE PRECISION DFLOAT

        INTRINSIC ABS
        DOUBLE PRECISION ABS

        INTRINSIC SIGN
        DOUBLE PRECISION SIGN


      DOUBLE PRECISION DASUM
      INTEGER IDAMAX

      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO = 0.0D0)
      PARAMETER (ONE = 1.0D0)
      PARAMETER (TWO = 2.0D0)
C
C     INTERNAL VARIABLES
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN, ESTOLD, TEMP
C
      SAVE
C
      IF (KASE .EQ. 0) THEN
         DO 10,I = 1,N
            X(I) = ONE/DFLOAT(N)
 10             CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      ENDIF
C
      GOTO (100, 200, 300, 400, 500) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
C
 100   CONTINUE
      IF (N .EQ. 1) THEN
         V(1) = X(1)
         EST = ABS(V(1))
C        ... QUIT
         GOTO 510
      ENDIF
      EST = DASUM(N,X,1)
C
      DO 110,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I))
 110      CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
 200   CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
 220   CONTINUE
      DO 230,I = 1,N
         X(I) = ZERO
 230      CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
 300   CONTINUE
      CALL DCOPY(N,X,1,V,1)
      ESTOLD = EST
      EST = DASUM(N,V,1)
      DO 310,I = 1,N
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. ISGN(I) ) GOTO 320
 310      CONTINUE
C     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GOTO 410
C
 320   CONTINUE
C     TEST FOR CYCLING.
      IF (EST .LE. ESTOLD) GOTO 410
C
      DO 330,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I))
 330      CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
 400   CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF (   (  X(JLAST) .NE. ABS(X(J))  ) .AND.
     +       (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
C
C     ITERATION COMPLETE.  FINAL STAGE.
C
 410   CONTINUE
      ALTSGN = ONE
      DO 420,I = 1,N
         X(I) = ALTSGN * (ONE + DFLOAT(I-1)/DFLOAT(N-1))
         ALTSGN = -ALTSGN
 420      CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
C
C     ................ ENTRY   (JUMP = 5)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
 500   CONTINUE
      TEMP = TWO*DASUM(N,X,1)/DFLOAT(3*N)
      IF (TEMP. GT. EST) THEN
         CALL DCOPY(N,X,1,V,1)
         EST = TEMP
      ENDIF
C
 510   KASE = 0
      RETURN
C
      END
