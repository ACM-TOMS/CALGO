
      LOGICAL FUNCTION CKCOLR(MP, NP, DP, P, LDP1, LDP2, ZERCOL, Q, LDQ,
     *                        W, TOL, IERR)
C
C     PURPOSE
C
C     To check whether the polynomial matrix
C                                                    dp-1            dp
C        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s     + P(dp) * s  ,
C
C     is column reduced.
C
C     ARGUMENTS IN
C
C        MP - INTEGER.
C            The number of rows of the polynomial matrix P(s).
C            MP >= 1.
C        NP - INTEGER.
C            The number of columns of the polynomial matrix P(s).
C            NP >= 1.
C        DP - INTEGER.
C            The degree of the polynomial matrix P(s).
C            DP >= 0.
C        P - DOUBLE PRECISION array of DIMENSION (LDP1,LDP2,DP+1).
C            The leading MP by NP by (DP+1) part of this array must
C            contain the coefficients of the polynomial matrix P(s).
C            Specifically, P(i,j,k) must contain the coefficient of
C            s**(k-1) of the polynomial which is the (i,j)-th element
C            of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C            k = 1,2,...,DP+1.
C        LDP1 - INTEGER.
C            The leading dimension of array P as declared in the calling
C            program.
C            LDP1 >= MP.
C        LDP2 - INTEGER.
C            The second dimension of array P as declared in the calling
C            program.
C            LDP2 >= NP.
C
C        ARGUMENTS OUT
C
C        ZERCOL - LOGICAL array of DIMENSION at least (NP).
C            If ZERCOL(j) = .TRUE. then the j-th column of P(s) is zero;
C            otherwise the j-th column belongs to P1(s) (see METHOD).
C
C        WORKSPACE
C
C        Q - DOUBLE PRECISION array of DIMENSION (LDQ,NP).
C        LDQ - INTEGER.
C            The leading dimension of array Q as declared by the calling
C            program.
C            LDQ >= MP.
C        W - DOUBLE PRECISION array of DIMENSION (MP).
C
C     TOLERANCES
C
C        TOL - DOUBLE PRECISION.
C            A tolerance below which matrix elements are considered to
C            be zero. If the user sets TOL to be less than
C            EPS * (((DP+1)*MP)**2 * MAX(P(i,j,k))), then the tolerance
C            is taken as EPS * (((DP+1)*MP)**2 * MAX(P(i,j,k))),
C            i = 1,...,MP, j = 1,...,NP, k = 1,...,DP+1, where EPS is
C            the machine precision.
C
C     ERROR INDICATOR
C
C        IERR - INTEGER.
C           Unless the routine detects an error (see next section),
C           IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = 1 : Invalid input parameter(s).
C
C     METHOD
C
C     Let GAMC(P) be the contant matrix such that each of its columns
C     contains the coefficients of the highest power of s occurring
C     in the corresponding column of P(s), the so-called leading column
C     coefficient matrix. Then P(s) is called column reduced if there
C     exists a permutation matrix T such that P(s) = ( Z , P1(s) ) * T,
C     where Z is a zero matrix and GAMC(P1) has full column rank.
C
C     The algorithm used, which is in fact the QR decomposition of the
C     leading column coefficient matrix, is as follows:
C     The columns of the leading column coefficient matrix of P1(s) are
C     determined one by one, where a column is considered zero if its
C     Euclidean norm is less than TOL. To each new column the
C     Householder transformations are applied that have transformed the
C     submatrix of the former columns in upper triangular form. If the
C     new column is independent of its predecessors, then a new
C     Householder transformation is generated and applied such that the
C     augmented matrix is upper triangular.
C     The routine terminates after the last column of P(s) has been
C     treated or when the new column is dependent of its predecessors.
C
C     CONTRIBUTOR
C
C        A.J. Geurts (Eindhoven University of Technology).
C        C. Praagman (University of Groningen).
C
C     REVISIONS
C
C        1994, February 11.
C        1996, May 28.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
C     .. Scalar Arguments ..
      INTEGER MP, NP, DP, LDP1, LDP2, LDQ, IERR
      DOUBLE PRECISION TOL
C     .. Array Arguments ..
      DOUBLE PRECISION P(LDP1,LDP2,*), Q(LDQ,*), W(*)
      LOGICAL ZERCOL(*)
C     .. Local Scalars ..
      INTEGER H, J, J1, K
      DOUBLE PRECISION EPS, TOLER, NORM, TAU
      LOGICAL FULLRK, NOTCJ
C     .. External Subroutines/Functions ..
      EXTERNAL DCOPY, DLAMCH, DLANGE, DLARFG, DLARFV, DLASET, DNRM2
      DOUBLE PRECISION DLAMCH, DLANGE, DNRM2
C     .. Intrinsic Functions ..
      INTRINSIC DBLE, MAX, MIN
C
C     .. Executable Statements ..
C
C     Check input parameters
C
      IF (MP.LT.1 .OR. NP.LT.1 .OR. DP.LT.0 .OR. LDP1.LT.MP .OR.
     *    LDP2.LT.NP) THEN
          IERR = 1
          RETURN
      END IF
C
C     Computation of the tolerance. EPS is the machine precision.
C
      EPS = DLAMCH('E')
      TOLER = ZERO
      DO 10 K = 1, DP + 1
         TOLER = MAX(TOLER, DLANGE('M', MP, NP, P(1,1,K), LDP1, W))
   10 CONTINUE

      TOLER = DBLE(((DP+1) * MP))**2 * TOLER * EPS
      IF (TOLER .LT. TOL) TOLER = TOL
C
      CALL DLASET('G', MP, NP, ZERO, ZERO, Q, LDQ)
      J = 1
      J1 = 1
      FULLRK = .TRUE.
C     WHILE (FULLRK and J1 <= NP) DO
   20 IF (FULLRK .AND. J1.LE.NP) THEN
C
C        Find the j-th column of the leading column coefficient matrix
C        of P1(s) and put it in W.
C
         K = DP + 1
         NOTCJ = .TRUE.
C        WHILE (j-th column not found) DO
   30    IF (NOTCJ .AND. K.GE.1) THEN
            NORM = DNRM2(MP, P(1,J1,K), 1)
            IF (NORM .GE. TOLER) THEN
               CALL DCOPY(MP, P(1,J1,K), 1, W, 1)
               NOTCJ = .FALSE.
            END IF
            K = K - 1
            GO TO 30
         END IF
C        END WHILE 30
C
C        Check whether the j-th column is linearly independent of the
C        preceeding columns.
C
         IF (NOTCJ) THEN
            ZERCOL(J1) = .TRUE.
            J1 = J1 + 1
         ELSE
            ZERCOL(J1) = .FALSE.
C
C           Apply the Householder transformations Qh,
C           h = 1,...,min(mp,j) - 1, to W.
C
            DO 40 H = 1, MIN(MP,J) - 1
               CALL DLARFV(MP-H+1, Q(H+1,H), 1, Q(H,H), W(H), W(H+1), 1)
   40       CONTINUE
            NORM = DNRM2(MP-J+1, W(J), 1)
            IF (NORM .LT. TOLER) THEN
               FULLRK = .FALSE.
            ELSE
C
C              Generate the Householder transformation Qj.
C
               IF (J .LT. MP) THEN
                  CALL DLARFG(MP-J+1, W(J), W(J+1), 1, TAU)
                  CALL DCOPY(MP-J, W(J+1), 1, Q(J+1,J), 1)
                  Q(J,J) = TAU
               ELSE
                  Q(J,J) = ZERO
               END IF
            END IF
            J = J + 1
            J1 = J1 + 1
         END IF
         GO TO 20
      END IF
C     END WHILE 20
      CKCOLR = FULLRK
      RETURN
C *** Last line of CKCOLR ***
      END





C
      LOGICAL FUNCTION CKGAMC(MP, INV, GAMC, LDG, Q, LDQ, RWORK, TOL)
C
C     PURPOSE
C
C     To check whether the leading coefficient matrix has still full
C     column rank after a new column has been appended.
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C        MP - INTEGER.
C            The number of rows of the matrix GAMC.
C            MP >= 1.
C        INV - INTEGER.
C            The number of columns of the matrix GAMC.
C            INV >= 1.
C        GAMC - DOUBLE PRECISION array of DIMENSION (LDG,INV)
C            The leading MP by INV part of this array must contain the
C            leading coefficient matrix GAMC of which the first INV - 1
C            columns have been transformed in upper triangular form.
C            Note: this array is overwritten
C        LDG - INTEGER.
C            The leading dimension of array GAMC as declared in the
C            calling program.
C            LDG >= MP.
C        Q - DOUBLE PRECISION array of DIMENSION (LDQ,INV)
C            The leading MP by INV - 1 part of this array must contain
C            the constants and the vectors of the elementary Householder
C            transformations, as supplied by DLARFG, by which the first
C            INV - 1 columns of GAMC have been transformed into an upper
C            triangular matrix.
C        LDQ - INTEGER.
C            The leading dimension of array Q as declared in the calling
C            program.
C            LDQ >= MP.
C
C     ARGUMENTS OUT
C
C        GAMC - DOUBLE PRECISION array of DIMENSION (LDG,INV)
C            The leading MP by INV part of this array contains the
C            leading coefficient matrix GAMC transformed in upper
C            triangular form.
C        Q - DOUBLE PRECISION array of DIMENSION (LDQ,INV)
C            The leading MP by INV part of this array contains the
C            constants and the vectors of the elementary Householder
C            transformations by which GAMC has been transformed into an
C            upper triangular matrix.
C
C        WORKSPACE
C
C        RWORK - DOUBLE PRECISION array of DIMENSION (MP)
C
C     TOLERANCES
C
C        TOL - DOUBLE PRECISION.
C            A tolerance below which matrix elements are considered to
C            be zero.
C
C     METHOD
C
C     Let the first INV - 1 columns of GAMC be linearly independent,
C     which has been checked by former calls of CKGAMC. A new column is
C     appended. To this column the Householder transformations are
C     applied that have transformed the matrix of the former columns in
C     upper triangular form. If the new column is independent of its
C     predecessors, then a new Householder transformation is generated
C     and applied such that the augmented matrix is upper triangular.
C
C     CONTRIBUTOR
C
C        A.J. Geurts.
C
C     REVISIONS
C
C        1993, October 29.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
C     .. Scalar Arguments ..
      INTEGER MP, INV, LDG, LDQ
      DOUBLE PRECISION TOL
C     .. Array Arguments ..
      DOUBLE PRECISION GAMC(LDG,*), Q(LDQ,*), RWORK(*)
C     .. Local Scalars ..
      INTEGER J
      DOUBLE PRECISION NORM, TAU
C     .. External Subroutines/Functions ..
      EXTERNAL DCOPY, DLARFG, DLARFV, DLAVEC, DNRM2
      DOUBLE PRECISION DNRM2
C
C     .. Executable Statements ..
C
      IF (INV .LE. MP) THEN
         CALL DLAVEC(MP, ZERO, Q(1,INV), 1)
C
C        Check whether the INV-th column is linearly independent of the
C        preceeding columns by applying the Householder transformations
C        Q(h), h = 1,...,INV-1, to the INV-th column.
C
         CALL DCOPY(MP, GAMC(1,INV), 1, RWORK, 1)
         DO 10 J = 1, INV - 1
            CALL DLARFV(MP-J+1, Q(J+1,J), 1, Q(J,J), RWORK(J),
     *                  RWORK(J+1), 1)
   10    CONTINUE
         NORM = DNRM2(MP-INV+1, RWORK(INV), 1)
         IF (NORM .LT. TOL) THEN
            CKGAMC = .FALSE.
         ELSE
C
C           Generate the Householder transformation Q(INV) if INV < MP.
C
            CALL DCOPY(MP, RWORK, 1, GAMC(1,INV), 1)
            IF (INV .LT. MP) THEN
               CALL DLARFG(MP-INV+1, RWORK(INV), RWORK(INV+1), 1, TAU)
               GAMC(INV,INV) = RWORK(INV)
               CALL DLAVEC(MP-INV, ZERO, GAMC(INV+1,INV), 1)
               CALL DCOPY(MP-INV, RWORK(INV+1), 1, Q(INV+1,INV), 1)
               Q(INV,INV) = TAU
            ELSE
               Q(INV,INV) = ZERO
            END IF
            CKGAMC = .TRUE.
         END IF
      ELSE
         CKGAMC = .FALSE.
      END IF
      RETURN
C *** Last line of CKGAMC ***
      END






C
      SUBROUTINE COLRD1(MP, NP, DP, P, LDP1, LDP2, DR, DU, R, LDR1,
     *                  LDR2, U, LDU1, LDU2, ZERCOL, MU, S, SK, A, LDA,
     *                  AB, LDAB, Q, LDQ, Y, LDY, YI, GAMC, LDG, W,
     *                  TOL, IERR)
C
C     PURPOSE
C
C     To compute for a given polynomial matrix
C                                                    dp-1            dp
C        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s     + P(dp) * s  ,
C
C     which is not column reduced, a unimodular polynomial matrix U(s)
C     such that R(s) = P(s) * U(s) is column reduced.
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C        MP - INTEGER.
C            The number of rows of the polynomial matrix P(s).
C            MP >= 1.
C        NP - INTEGER.
C            The number of columns of the polynomial matrix P(s).
C            NP >= 1.
C        DP - INTEGER.
C            The degree of the polynomial matrix P(s).
C            DP >= 0.
C        P - DOUBLE PRECISION array of DIMENSION (LDP1,LDP2,DP+1).
C            The leading MP by NP by (DP+1) part of this array must
C            contain the coefficients of the polynomial matrix P(s).
C            Specifically, P(i,j,k) must contain the coefficient of
C            s**(k-1) of the polynomial which is the (i,j)-th element
C            of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C            k = 1,2,...,DP+1.
C        LDP1 - INTEGER.
C            The leading dimension of array P as declared in the calling
C            program.
C            LDP1 >= MP.
C        LDP2 - INTEGER.
C            The second dimension of array P as declared in the calling
C            program.
C            LDP2 >= NP.
C
C     ARGUMENTS OUT
C
C        DR - INTEGER.
C            The degree of the column reduced polynomial matrix R(s).
C        DU - INTEGER.
C            The degree of the unimodular polynomial matrix U(s).
C        R - DOUBLE PRECISION array of DIMENSION (LDR1,LDR2,DP+1).
C            The leading MP by NP by (DR+1) part of this array contains
C            the coefficients of the column reduced polynomial matrix
C            R(s). Specifically, R(i,j,k) contains the coefficient of
C            s**(k-1) of the polynomial which is the (i,j)-th element of
C            R(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C            k = 1,2,...,DR+1.
C        LDR1 - INTEGER.
C            The leading dimension of array R as declared in the calling
C            program.
C            LDR1 >= MP.
C        LDR2 - INTEGER.
C            The second dimension of array R as declared in the calling
C            program.
C            LDR2 >= NP.
C        U - DOUBLE PRECISION array of DIMENSION (LDU1,LDU2,NP*DP+1).
C            The leading NP by NP by (DU+1) part of this array contains
C            the coefficients of the unimodular polynomial matrix U(s).
C            Specifically, U(i,j,k) contains the coefficient of s**(k-1)
C            of the polynomial which is the (i,j)-th element of U(s),
C            where i = 1,2,...,NP, j = 1,2,...,NP and k = 1,2,...,DU+1.
C        LDU1 - INTEGER.
C            The leading dimension of array U as declared in the calling
C            program.
C            LDU1 >= NP.
C        LDU2 - INTEGER.
C            The second dimension of array U as declared in the calling
C            program.
C            LDU2 >= NP.
C        ZERCOL - LOGICAL array of DIMENSION at least (NP).
C            If ZERCOL(j) = .TRUE. then the j-th column of R(s) is zero;
C            otherwise the j-th column belongs to R1(s) (see METHOD).
C
C        WORKSPACE
C
C        MU - INTEGER array of DIMENSION at least (mamax+1),
C            where mamax = (NP * DP + 1) * MP.
C            On exit, this array contains the row indices of the left
C            upper elements of the right invertible diagonal submatrices
C            of A'.
C        S - INTEGER array of DIMENSION at least (mamax).
C            On exit, this array contains the column indices of the
C            pivots of A', that is A transformed in upper staircase
C            form.
C        SK - INTEGER array of DIMENSION at least (2*NP).
C        A - DOUBLE PRECISION array of DIMENSION (LDA,namax),
C            where namax = mamax + NP.
C            On exit, the upper block diagonal part of this array
C            contains the transformed matrix A', the lower part contains
C            the vectors of the Householder transformations.
C        LDA - INTEGER.
C            The leading dimension of array A as declared in the calling
C            program.
C            LDA >= mamax + 1.
C        AB - DOUBLE PRECISION array of DIMENSION (LDAB,namax).
C            Array in which the transformed matrix Ab is saved.
C        LDAB - INTEGER.
C            The leading dimension of array AB as declared in the
C            calling program.
C            LDAB >= mamax + 1.
C        Q - DOUBLE PRECISION array of DIMENSION (LDQ,NP).
C        LDQ - INTEGER.
C            The leading dimension of array Q as declared in the calling
C            program.
C            LDQ >= MP.
C        Y - DOUBLE PRECISION array of DIMENSION (LDY,NP*DP+2).
C            Array in which a null vector of sA - E or sA' - E is
C            stored.
C        LDY - INTEGER.
C            The leading dimension of array Y as declared in the calling
C            program.
C            LDA >= namax.
C        YI - DOUBLE PRECISION array of DIMENSION at least (NP).
C        GAMC - DOUBLE PRECISION array of DIMENSION (LDG,NP)
C            On exit, this array contains the leading column coefficient
C            matrix of R(s), which has full column rank.
C        LDG - INTEGER.
C            The leading dimension of array GAMC as declared in the
C            calling program.
C            LDG >= MP.
C        W - DOUBLE PRECISION array of DIMENSION (mamax+namax).
C
C     TOLERANCES
C
C        TOL - DOUBLE PRECISION.
C            A tolerance below which matrix elements are considered to
C            be zero.
C
C     ERROR INDICATOR
C
C        IERR - INTEGER.
C           Unless the routine detects an error (see next section),
C           IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = 2 : No column reduced R(s) has been found for the
C                   maximum b, (see METHOD).
C        IERR = 3 : The computation of a null vector has failed,
C                   because a diagonal block of A' is not right
C                   invertible.
C        IERR = 4 : The computation of R(s) has failed, because a
C                   computed null vector is not s**B times a polynomial
C                   vector.
C
C     METHOD
C
C     Let GAMC(P) be the constant matrix such that each of its columns
C     contains the coefficients of the highest power of s occurring in
C     the corresponding column of P(s), the so-called leading column
C     coefficient matrix. Then P(s) is called column reduced if there
C     exists a permutation matrix T such that P(s) = ( Z , P1(s) ) * T,
C     where Z is a zero matrix and GAMC(P1) has full column rank.
C
C     Let (U(s),Z(s))' be a minimal polynomial basis (MPB) for
C          b
C     Ker(s P(s), -I), for some b > 0. It has been proved, see [1], that
C     if b is greater than d'c(P), the sum of all but the smallest
C                                                                -b
C     column degrees of P(s), then U(s) is unimodular and R(s) = s  Z(s)
C     is column reduced and P(s) * U(s) = R(s).
C     The routine computes an MPB for b = 1,2,... and checks for each b
C     whether R(s) is column reduced, i.e. whether GAMC(R1) has full
C     column rank. The algorithm finishes with U(s) and R(s) as soon as
C     R(s) is column reduced.
C
C     REFERENCES
C
C     [1] Neven, W.H.L. and Praagman, C.
C         Column Reduction of Polynomial Matrices.
C         Linear Algebra and its Applications 188, 189, pp. 569-589,
C         1993.
C
C     CONTRIBUTORS
C
C        A.J. Geurts (Eindhoven University of Technology).
C        C. Praagman (University of Groningen).
C
C     REVISIONS
C
C        1993, November 8.
C        1996, May 21.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER MP, NP, DP, LDP1, LDP2, DR, DU, LDR1, LDR2, LDU1, LDU2,
     *        LDA, LDAB, LDQ, LDY, LDG, IERR
      DOUBLE PRECISION TOL
C     .. Array Arguments ..
      INTEGER MU(*), S(*), SK(*)
      DOUBLE PRECISION P(LDP1,LDP2,*), R(LDR1,LDR2,*), U(LDU1,LDU2,*),
     *                 A(LDA,*), AB(LDAB,*), Q(LDQ,*), Y(LDY,*), YI(*),
     *                 GAMC(LDG,*),  W(*)
      LOGICAL ZERCOL(*)
C     .. Local Scalars ..
      INTEGER BMAX, B, I, IC, IR, IND, NNV, NNVB, NCR1, J, K, KK,
     *        L, C1AI, DUB, DUR1, MA, NA, MAI, NAI, MUI
      DOUBLE PRECISION NORM, TAU
      LOGICAL COLRDC
C     .. External Subroutines/Functions ..
      EXTERNAL CKGAMC, COMPTV, COMPTY, COMPYI, DCOPY, DGEMV, DGER,
     *         DLACPY, DLAVEC, DLASET, DNRM2, HHTRAN, MKPENC
      DOUBLE PRECISION DNRM2
      LOGICAL CKGAMC
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C
C     .. Executable Statements ..
C
      BMAX = (NP-1) * DP + 1
      CALL MKPENC(MP, NP, DP, P, LDP1, LDP2, AB, LDAB, MA, NA)
      NNVB = 0
      DUB = -1
      DO 10 K = 1, NP * DP + 1
         CALL DLASET('G', NP, NP, ZERO, ZERO, U(1,1,K), LDU1)
   10 CONTINUE
      COLRDC = .FALSE.
      MU(1) = 1
      B = 0
C     WHILE (.NOT.COLRDC and B < BMAX) DO
   20 IF (.NOT.COLRDC .AND. B.LT.BMAX) THEN
         B = B + 1
C
C        Initialization of A.
C
         IF (B .EQ. 1) THEN
            CALL DLACPY('G', MA, NA, AB, LDAB, A, LDA)
            CALL DLAVEC(NA, ZERO, A(MA+1,1), LDA)
         ELSE
            CALL DLACPY('G', MA+1, NA, AB, LDAB, A, LDA)
            J = NP + MU(B-1)
            L = MA - MP - MU(B-1) + 1
            CALL DLASET('G', MP, J-1, ZERO, ZERO, A(MA+2,1), LDA)
            CALL DLASET('G', MP+1, L, ZERO, ZERO, A(MA+1,J), LDA)
            CALL DLASET('G', MP+1, MP, ZERO, ONE, A(MA+1,L+J), LDA)
            MA = MA + MP
            CALL DLASET('G', MA+1, MP, ZERO, ZERO, A(1,NA+1), LDA)
            NA = NA + MP
            DO 30 K = MU(B-1), MU(B) - 1
               TAU = A(K+1,S(K))
               IF (TAU .NE. ZERO) THEN
                  J = NP + K
                  L = MA - MP - K + 1
                  W(MA+1) = ONE
                  CALL DCOPY(L-1, A(K+2,S(K)), 1, W(MA+2), 1)
                  CALL DGEMV('N', MP, L, ONE, A(MA-MP+1,J), LDA,
     *                        W(MA+1), 1, ZERO, W, 1)
                  CALL DGER(MP, L, -TAU, W, 1, W(MA+1), 1,
     *                      A(MA-MP+1,J), LDA)
               END IF
   30       CONTINUE
         END IF
C
         I = B
         DU = DUB
         NNV = NNVB
         DR = -1
         DO 40 K = 1, DU + 1
            CALL DLASET('G', NP, NP-NNV, ZERO, ZERO, U(1,NNV+1,K),
     *                  LDU1)
   40    CONTINUE
         DO 50 K = 1, DP + 1
            CALL DLASET('G', MP, NP, ZERO, ZERO, R(1,1,K), LDR1)
   50    CONTINUE
         NCR1 = 0
         COLRDC = .TRUE.
C
C        NNV is the number of already found null vectors. With these
C        null vectors correspond the MP by NNV column reduced matrix
C        R(s). NCR1 is the number of columns in R1(s).
C
C        WHILE (NNV < NP, which means that not all null vectors of
C               sA - E  have been found, and R(s) is column reduced) DO
   60    IF ((NNV.LT.NP) .AND. COLRDC) THEN
C
C           Determine A(i,i), the i(th) right invertible diagonal block.
C           The left upper element of A(i,i) is A(MUI,C1AI) and the row
C           and column dimensions are MAI and NAI, respectively.
C
            MUI = MU(I)
            MAI = 0
            IF (I .EQ. 1) THEN
               C1AI = 1
               NAI = NP
            ELSE
               C1AI = MU(I-1) + NP
               NAI = MU(I) - MU(I-1)
            END IF
C
C           Compute the null vectors of A(i,i) and the corresponding
C           nullvectors of sA - E one by one, append the appropriate
C           parts of a computed null vector to U(s) and R(s) and check
C           whether R(s) remains column reduced.
C
            K = 1
C           WHILE (K <= NAI and R(s) still column reduced) DO
   70       IF (K.LE.NAI .AND. COLRDC) THEN
               IR = MUI + MAI
               IC = C1AI + K - 1
               L = MA - IR + 1
               IF (L .GT. 0) THEN
                  NORM = DNRM2(L, A(IR,IC), 1)
               ELSE
                  NORM = ZERO
               END IF
               IF (NORM .GE. TOL) THEN
C
C                 Generate and apply the Householder transformation for
C                 the IC-th column of A.
C
                  CALL HHTRAN(MA, NA, IC, IR, A, LDA, W, W(MA+1))
                  CALL DCOPY(L, W, 1, A(IR+1,IC), 1)
                  S(MUI+MAI) = IC
                  MAI = MAI + 1
                  SK(MAI) = IC - C1AI + 1
               ELSE
C
C                 Compute the null vector Y of sA - E corresponding
C                 to the IC-th column of A.
C
                  IF (L .GT. 0) THEN
                     CALL DLAVEC(L, ZERO, A(IR,IC), 1)
                  END IF
                  IF (MAI .GT. 0) THEN
                     CALL DCOPY(MAI, A(MUI,IC), 1, W, 1)
                     CALL COMPYI(MAI, K-1, A(MUI,C1AI), LDA, SK, W,
     *                           YI, IERR)
                     IF (IERR .NE. 0) RETURN
                  ELSE
                     CALL DLAVEC(K-1, ZERO, YI, 1)
                  END IF
                  YI(K) = -ONE
                  CALL COMPTY(NP, NA, I, K, YI, MU, A, LDA, S, Y, LDY,
     *                        SK(NP+1), W, IERR)
                  IF (IERR .NE. 0) THEN
                     IERR = 3
                     RETURN
                  END IF
                  CALL COMPTV(MA, NA, I, IR, A, LDA, S, Y, LDY, W)
C
C                 Append the first NP by I block of Y to the unimodular
C                 U, and the last MP by I block to the column reduced R
C                 and the last MP elements of the I-th column of Y to
C                 the leading column coefficient matrix GAMC.
C
                  NNV = NNV + 1
                  DUR1 = 0
                  DO 80 KK = 1, I
                     NORM = DNRM2(NP, Y(1,KK), 1)
                     IF (NORM .GE. TOL) THEN
                        DUR1 = KK - 1
                        CALL DCOPY(NP, Y(1,KK), 1, U(1,NNV,KK), 1)
                     ELSE
                        CALL DLAVEC(NP, ZERO, U(1,NNV,KK), 1)
                     END IF
   80             CONTINUE
                  DU = MAX(DU, DUR1)
                  IND = NA - MP + 1
                  DO 90 KK = 1, B
                     NORM = DNRM2(MP, Y(IND,KK), 1)
                     IF (NORM .NE. ZERO) THEN
                        IERR = 4
                        RETURN
                     END IF
   90             CONTINUE
                  DUR1 = -1
                  DO 100 KK = B + 1, I
                     NORM = DNRM2(MP, Y(IND,KK), 1)
                     IF (NORM .GE. TOL) THEN
                        DUR1 = KK - B - 1
                        CALL DCOPY(MP, Y(IND,KK), 1, R(1,NNV,KK-B), 1)
                     ELSE
                        CALL DLAVEC(MP, ZERO, R(1,NNV,KK-B), 1)
                     END IF
  100             CONTINUE
                  IF (DUR1 .NE. -1) THEN
                     NCR1 = NCR1 + 1
                     ZERCOL(NNV) = .FALSE.
                     DR = MAX(DR, DUR1)
                     CALL DCOPY(MP, Y(IND,DUR1+B+1), 1, GAMC(1,NCR1), 1)
C
C                    Check whether R is still column reduced.
C
                     COLRDC = CKGAMC(MP, NCR1, GAMC, LDG, Q, LDQ, W,
     *                               TOL)
                  ELSE
                     ZERCOL(NNV) = .TRUE.
                  END IF
               END IF
               K = K + 1
               GO TO 70
            END IF
C           END WHILE 70
C
C           Save the transformed A if I = B.
C
            IF (I .EQ. B) THEN
               CALL DLACPY('G', MA+1, NA, A, LDA, AB, LDAB)
               NNVB = NNV
               DUB = DU
            END IF
C
            IF (NNV.LE.NP .AND. COLRDC) THEN
C
C              If i > b + dp + 1, then the degree of the next computed
C              column in R(s) will be greater than dp. This will occur
C              if a null vector is not found, due to rounding errors.
C
               COLRDC = I .LE. (B + DP + 1)
               I = I + 1
               MU(I) = MU(I-1) + MAI
            END IF
            GO TO 60
         END IF
C        END WHILE 60
         GO TO 20
      END IF
C     END WHILE 20
C
      IF (.NOT. COLRDC) IERR = 2
      RETURN
C *** Last line of COLRD1 ***
      END






C
      SUBROUTINE COLRED(MP, NP, DP, P, LDP1, LDP2, DR, DU, R, LDR1,
     *                  LDR2, U, LDU1, LDU2, ZERCOL, IWORK, RWORK,
     *                  TOL, IERR)
C
C     PURPOSE
C
C     To compute for a given polynomial matrix
C                                                     dp-1            dp
C        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s     + P(dp) * s  ,
C
C     a unimodular polynomial matrix U(s) such that R(s) = P(s) * U(s)
C     is column reduced.
C
C     ARGUMENTS IN
C
C        MP - INTEGER.
C            The number of rows of the polynomial matrix P(s).
C            MP >= 1.
C        NP - INTEGER.
C            The number of columns of the polynomial matrix P(s).
C            NP >= 1.
C        DP - INTEGER.
C            The degree of the polynomial matrix P(s).
C            DP >= 0.
C        P - DOUBLE PRECISION array of DIMENSION (LDP1,LDP2,DP+1).
C            The leading MP by NP by (DP+1) part of this array must
C            contain the coefficients of the polynomial matrix P(s).
C            Specifically, P(i,j,k) must contain the coefficient of
C            s**(k-1) of the polynomial which is the (i,j)-th element
C            of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C            k = 1,2,...,DP+1.
C        LDP1 - INTEGER.
C            The leading dimension of array P as declared in the calling
C            program.
C            LDP1 >= MP.
C        LDP2 - INTEGER.
C            The second dimension of array P as declared in the calling
C            program.
C            LDP2 >= NP.
C
C     ARGUMENTS OUT
C
C        DR - INTEGER.
C            The degree of the column reduced polynomial matrix R(s).
C        DU - INTEGER.
C            The degree of the unimodular polynomial matrix U(s).
C        R - DOUBLE PRECISION array of DIMENSION (LDR1,LDR2,DP+1).
C            The leading MP by NP by (DR+1) part of this array contains
C            the coefficients of the column reduced polynomial matrix
C            R(s). Specifically, R(i,j,k) contains the coefficient of
C            s**(k-1) of the polynomial which is the (i,j)-th element of
C            R(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C            k = 1,2,...,DR+1.
C        LDR1 - INTEGER.
C            The leading dimension of array R as declared in the calling
C            program.
C            LDR1 >= MP.
C        LDR2 - INTEGER.
C            The second dimension of array R as declared in the calling
C            program.
C            LDR2 >= NP.
C        U - DOUBLE PRECISION array of DIMENSION (LDU1,LDU2,NP*DP+1).
C            The leading NP by NP by (DU+1) part of this array contains
C            the coefficients of the unimodular polynomial matrix U(s).
C            Specifically, U(i,j,k) contains the coefficient of s**(k-1)
C            of the polynomial which is the (i,j)-th element of U(s),
C            where i = 1,2,...,NP, j = 1,2,...,NP and k = 1,2,...,DU+1.
C        LDU1 - INTEGER.
C            The leading dimension of array U as declared in the calling
C            program.
C            LDU1 >= NP.
C        LDU2 - INTEGER.
C            The second dimension of array U as declared in the calling
C            program.
C            LDU2 >= NP.
C        ZERCOL - LOGICAL array of DIMENSION at least (NP).
C            If ZERCOL(j) = .TRUE. then the j-th column of R(s) is zero;
C            otherwise the j-th column belongs to R1(s) (see METHOD).
C
C     WORKSPACE
C
C        IWORK - INTEGER array of DIMENSION at least (liwork),
C            where liwork = 2*mamax + 2 * NP + 1,
C            and mamax = (NP*DP+1) * MP.
C        RWORK - DOUBLE PRECISION array of DIMENSION at least (lrwork),
C            where lrwork = (2*mamax+NP*DP+5) * namax + mamax + (2*MP+1) * NP,
C            and namax = mamax + NP.
C                           .
C     TOLERANCES
C
C        TOL - DOUBLE PRECISION.
C            A tolerance below which matrix elements are considered to
C            be zero. If the user sets TOL to be less than
C            EPS * (((DP+1)*MP)**2 * MAX(P(i,j,k))), then the tolerance
C            is taken as EPS * (((DP+1)*MP)**2 * MAX(P(i,j,k))),
C            i = 1,...,MP, j = 1,...,NP, k = 1,..., DP + 1, where EPS is
C            the machine precision.
C
C     ERROR INDICATOR
C
C        IERR - INTEGER.
C           Unless the routine detects an error (see next section),
C           IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = 1 : On entry, MP < 1 or NP < 1 or DP < 0 or
C                   LDP1 < MP or LDP2 < NP or
C                   LDR1 < MP or LDR2 < NP or
C                   LDU1 < NP or LDU2 < NP.
C        IERR = 2 : P(s) is not column reduced, and no column reduced 
C                   R(s) has been found.
C
C     METHOD
C
C     Let GAMC(P) be the constant matrix such that each of its columns
C     contains the coefficients of the highest power of s occurring in
C     the corresponding column of P(s), the so-called leading column
C     coefficient matrix. Then P(s) is called column reduced if there
C     exists a permutation matrix T such that P(s) = ( Z , P1(s) ) * T,
C     where Z is a zero matrix and GAMC(P1) has full column rank.
C
C     Let (U(s),Z(s))' be a minimal polynomial basis (MPB) for
C          b
C     Ker(s P(s), -I), for some b > 0. It has been proved, see [1], that
C     if b is greater than d'c(P), the sum of all but the smallest
C                                                                -b
C     column degrees of P(s), then U(s) is unimodular and R(s) = s  Z(s)
C     is column reduced and P(s) * U(s) = R(s).
C                                           b
C     The routine uses a linearization of (s P(s), -I) to compute an MPB
C     for b = 1,2,... and checks for each b whether R(s) is column
C     reduced, i.e. whether GAMC(R1) has full column rank. The algorithm
C     finishes with U(s) and R(s) as soon as R(s) is column reduced.
C
C     REFERENCES
C
C     [1] Neven, W.H.L. and Praagman, C.
C         Column Reduction of Polynomial Matrices.
C         Linear Algebra and its Applications 188, 189, pp. 569-589,
C         1993.
C     [2] Geurts, A.J. and Praagman, C.
C         A Fortran subroutine for column reduction of polynomial
C         matrices. EUT Report 94-WSK-01, Eindhoven, June 1994.
C
C     NUMERICAL ASPECTS
C
C     The algorithm used by the routine involves the construction of a
C                                                    b
C     special staircase form of a linearization of (s  P(s), -I) with
C     pivots considered to be non-zero when they are greater than or
C     equal to TOL. These pivots are then inverted in order to construct
C                         b
C     the columns of ker(s  P(s), -I).
C     The user is recommended to choose TOL of the order of the relative
C     error in the elements of P(s). If TOL is chosen to be too small,
C     then a very small element of insignificant value may be taken as
C     pivot. As a consequence, the correct null-vectors, and hence R(s),
C     may not be found. In the case that R(s) has not been found and in
C     the case that the elements of the computed U(s) and R(s) are large
C     relative to the elements of P(s) the user should consider trying
C     several values of TOL.
C
C     CONTRIBUTORS
C
C        A.J. Geurts (Eindhoven University of Technology).
C        C. Praagman (University of Groningen).
C
C     REVISIONS
C
C        1994, February 11.
C        1996, June 5.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER MP, NP, DP, LDP1, LDP2, DR, DU, LDR1, LDR2, LDU1, LDU2,
     *        IERR
      DOUBLE PRECISION TOL
C     .. Array Arguments ..
      INTEGER IWORK(*)
      DOUBLE PRECISION P(LDP1,LDP2,*), R(LDR1,LDR2,*), U(LDU1,LDU2,*),
     *                 RWORK(*)
      LOGICAL ZERCOL(*)
C     .. Local Scalars ..
      INTEGER LDA, LDAB, LDQ, LDY, LDG, MU, S, SK, A, AB, Q, Y, YI,
     *        GAMC, BMAX, DP1, K, MAMAX, NAMAX
      DOUBLE PRECISION EPS, TOLER
      LOGICAL COLRDC, PKZERO
C     .. External Subroutines/Functions ..
      EXTERNAL COLRD1, CKCOLR, DLACPY, DLAMCH, DLANGE, DLASET
      DOUBLE PRECISION DLAMCH, DLANGE
      LOGICAL CKCOLR
C     .. Intrinsic Functions ..
      INTRINSIC DBLE, MAX
C
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IF (MP.LT.1 .OR. NP.LT.1 .OR. DP.LT.0 .OR. LDP1.LT.MP .OR.
     *   LDP2.LT.NP .OR. LDR1.LT.MP .OR. LDR2.LT.NP .OR.
     *   LDU1.LT.NP .OR. LDU2.LT.NP) THEN
         IERR = 1
         RETURN
      END IF
C
      IERR = 0
C
C     Computation of the tolerance. EPS is the machine precision.
C
      EPS = DLAMCH('E')
      TOLER = ZERO
      DO 10 K = 1, DP + 1
         TOLER = MAX(TOLER, DLANGE('M', MP, NP, P(1,1,K), LDP1, RWORK))
   10 CONTINUE
C
      TOLER = DBLE(((DP+1) * MP)**2) * TOLER * EPS
      IF (TOLER .LT. TOL) TOLER = TOL
C
C     Computation of the true degree of P(s).
C
      K = DP + 2
      PKZERO = .TRUE.
C     WHILE (P(k) is a zero matrix) DO
   20 IF (PKZERO .AND. (K.GT.1)) THEN
         K = K - 1
         PKZERO = (DLANGE('M', MP, NP, P(1,1,K), LDP1, RWORK) .EQ. ZERO)
         GO TO 20
      END IF
C     END WHILE 20
      DP1 = K - 1
C
C     Check whether P(s) is already column reduced.
C
      Q = MP + 1
      LDQ = MP
      COLRDC = CKCOLR(MP, NP, DP1, P, LDP1, LDP2, ZERCOL, RWORK(Q), LDQ,
     *                RWORK, TOLER, IERR)
      IF (COLRDC) THEN
         DR = DP1
         DO 30 K = 1, DR + 1
            CALL DLACPY('G', MP, NP, P(1,1,K), LDP1, R(1,1,K), LDR1)
   30    CONTINUE
         DU = 0
         CALL DLASET('G', NP, NP, ZERO, ONE, U, LDU1)
         RETURN
      END IF
C
      BMAX = (NP - 1) * DP1 + 1
      MAMAX = (DP1 + BMAX) * MP
      NAMAX = MAMAX + NP
      LDA = MAMAX + 1
      LDAB = LDA
      LDY = NAMAX
      LDG = MP
      MU = 1
      S = MU + MAMAX + 1
      SK = S + MAMAX
      A = 1 + MAMAX + NAMAX 
      AB = A + (MAMAX + 1) * NAMAX
      Q = AB + (MAMAX + 1) * NAMAX
      Y = Q + MP * NP
      YI = Y + (BMAX + DP1 + 1) * NAMAX
      GAMC = YI + NP
C
      CALL COLRD1(MP, NP, DP1, P, LDP1, LDP2, DR, DU, R, LDR1, LDR2,
     *            U, LDU1, LDU2, ZERCOL, IWORK(MU), IWORK(S), IWORK(SK),
     *            RWORK(A), LDA, RWORK(AB), LDAB, RWORK(Q), LDQ,
     *            RWORK(Y), LDY, RWORK(YI), RWORK(GAMC), LDG, RWORK,
     *            TOLER, IERR)
C
C     Check whether the computed R(s) is column reduced.
C
      IF (IERR .EQ. 0) THEN
         Q = MP + 1
         COLRDC = CKCOLR(MP, NP, DR, R, LDR1, LDR2, ZERCOL, RWORK(Q),
     *                   LDQ, RWORK, TOL, IERR)
         IF (.NOT. COLRDC) IERR = 2
      ELSE
         IERR = 2
      END IF
      RETURN
C *** Last line of COLRED ***
      END






C
      SUBROUTINE COMPTV(MA, NA, I, J, A, LDA, S, Y, LDY, RWORK)
C
C     PURPOSE
C
C     To apply the Householder reflections, thus far applied to sA - E,
C     to the null vector Y(s) corresponding to a given column of the
C     transformed sA' - E, which transforms this null vector into a null
C     vector V(s) of the original pencil sA - E.
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C        MA - INTEGER.
C            The number of rows of matrix A.
C            MA >= 1.
C        NA - INTEGER.
C            The number of columns of matrix A.
C            NA >= MA.
C        I - INTEGER.
C            The actual number of columns in Y, which is also the index
C            of the diagonal block in A corresponding to the null vector
C            Y(s).
C            I >= 1.
C        J - INTEGER.
C            The row index in A which corresponds to the null vector
C            Y(s).
C            J >= 1.
C        A - DOUBLE PRECISION array of DIMENSION (LDA,NA).
C            The leading (MA + 1) by NA part of this array must contain
C            the transformed matrix A' in the upper part and the
C            Householder transformation vectors in the lower part.
C        LDA - INTEGER.
C            The leading dimension of array A as declared in the calling
C            program.
C            LDA >= MA + 1.
C        S - INTEGER ARRAY of DIMENSION at least (J-1).
C            The leading J - 1 elements of this array must contain the
C            indices of the pivots of the right invertible diagonal
C            submatrices, i.e., the pivot of A(m,m) is A(m,S(m)),
C            M=1,...,J-1. S(m) is also the index of the column in array
C            A in which the m-th non-trivial Householder transformation
C            vector is stored.
C        Y - DOUBLE PRECISION array of DIMENSION (LDY,I).
C            The leading NA by I part of this array must contain the
C            polynomial null vector Y(s) of sA' - E to be transformed,
C            where the t-th column (t = 1,...,I) must contain the
C            coefficient of s**(t-1).
C            Note: this array is overwritten.
C        LDY - INTEGER.
c            The leading dimension of array Y as declared in the calling
C            program.
C            LDY >= NA.
C
C     ARGUMENTS OUT
C
C        Y - DOUBLE PRECISION array of DIMENSION (LDY,I).
C            The leading NA by I part of this array contains the
C            transformed null vector V(s) = Q * Y(s), where Q  is the
C            product of the (J-1) Householder transformations Q(m).
C
C     WORKSPACE
C
C        RWORK - DOUBLE PRECISION array of DIMENSION (2*MA).
C
C     METHOD
C
C     Let
C         Q(m)  = ( I    0  )
C                 ( 0  P(m) )
C     be the elementary Householder transformation corresponding to the
C     pivot A(m,S(m)), augmented such that Q(m) is NA by NA, then
C     Q(1) Q(2) ... Q(J-1) Y ==> Y is computed.
C
C     CONTRIBUTOR
C
C        A.J. Geurts.
C
C     REVISIONS
C
C        1992, October 27.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER MA, NA, I, J, LDA, LDY
C     .. Array Arguments ..
      INTEGER S(*)
      DOUBLE PRECISION A(LDA,NA), Y(LDY,I), RWORK(2*NA)
C     .. Local Scalars ..
      INTEGER LEN, M, M1, MY
      DOUBLE PRECISION TAU
C     .. External Subroutines/Functions ..
      EXTERNAL DCOPY, DGEMV, DGER, DLAVEC
C
C     .. Executable Statements ..
C
C     Compute the matrix P(m) Y by w = Y'u and next Y = Y - tau * uw'.
C     for m = j-1, ... , 1.
C
      DO 20 M = J - 1, 1, -1
         IF (A(M+1,S(M)) .NE. ZERO) THEN
            LEN = MA - M + 1
            MY = NA - MA + M
            M1 = MA + 1
            CALL DCOPY(LEN, A(M+1,S(M)), 1, RWORK, 1)
            CALL DLAVEC(LEN, ZERO, RWORK(M1), 1)
            TAU = RWORK(1)
            RWORK(1) = ONE
            CALL DGEMV('T', LEN, I, ONE, Y(MY,1), LDY, RWORK, 1,
     *                 ZERO, RWORK(M1), 1)
            CALL DGER(LEN, I, -TAU, RWORK, 1, RWORK(M1), 1, Y(MY,1),
     *                LDY)
         END IF
   20 CONTINUE
      RETURN
C *** Last line of COMPTV ***
      END






C
      SUBROUTINE COMPTY(NP, NA, I, NYI, YI, MU, A, LDA, S, Y, LDY, SK,
     *                  RWORK, IERR)
C
C     PURPOSE
C
C     To compute a right null vector of the pencil sA' - E, where the
C     left part of A' is in staircase form. Actually, the computed
C     vector is the null vector of the corresponding left part of the
C     pencil.
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C        NP - INTEGER.
C            The number of columns of the polynomial matrix.
C            NP >= 1.
C        NA - INTEGER.
C            The number of columns of the matrix A.
C            NA >= 1.
C        I - INTEGER.
C            The index of the current diagonal block A(i,i) of the
C            matrix A being transformed into staircase form.
C            I >= 1.
C        NYI - INTEGER.
C            The length of the righthandside vector YI.
C            1 <= NYI <= NP.
C        YI - DOUBLE PRECISION array of DIMENSION at least (NP).
C            The right null vector of A(i,i).
C        MU - INTEGER array of DIMENSION at least (MA).
C            MU(k), k = 1, ..., i must contain the row index of the left
C            upper element of A(k,k).
C        A - DOUBLE PRECISION array of DIMENSION (LDA,NA).
C            The leading MU(i) - 1 by MU(i) + NP - 1 part of this array
C            must contain the part of the matrix A' which is in
C            staircase form.
C        LDA - INTEGER.
C            The leading dimension of array A as declared in the calling
C            program.
C            LDA >= MU(i) - 1.
C        S - INTEGER array of DIMENSION at least (MA).
C            The leading MU(i) - 1 elements of this array must contain
C            the column indices of the pivots of the right invertible
C            diagonal matrices A(k,k), k = 1, ..., i-1.
C
C     ARGUMENTS OUT
C
C        Y - DOUBLE PRECISION array of DIMENSION (LDY,NA).
C            The leading NA by i part of this array contains the
C            computed polynomial right null vector Y(s) of sA' - E,
C            where the j-th column contains the coefficient of s**(j-1).
C            The last NA - MU(i) - NP + 1 components of Y(s) are zero.
C        LDY - INTEGER.
C            The leading dimension of array Y as declared in the calling
C            program.
C            LDY >= NA.
C
C     WORK SPACE
C
C        SK - INTEGER array of DIMENSION at least (NP).
C        RWORK - DOUBLE PRECISION array of DIMENSION at least (NP).
C
C     ERROR INDICATOR
C
C        IERR - INTEGER.
C           Unless the routine detects an error (see next section),
C           IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C     IERR = k : A(k,k) is not right invertible.
C
C     METHOD
C
C     Let the pencil sA' - E, partially transformed up to the i-th
C     block, be
C
C        (  sA(1,1) sA(1,2)-E(1) sA(1,3)      . . .  sA(1,i)  . . .    )
C        (          sA(2,2)      sA(2,3)-E(2) . . .  sA(2,i)  . . .    )
C        (                         .                   .               )
C        (                                .            .               )
C        (                                       .     .               )
C        (                                           sA(i,i)  . . .    )
C        (                                                    . . .    )
C        (                                                    . . .    )
C
C     where A(k,k), k = 1,..., i is right invertible, A(k,k+1) is square
C     and E(k) = I of appropriate size.
C     Let Y(i,i), the (constant) right null vector of A(i,i), be given.
C     Then the routine computes a right null vector of sA' - E of the
C     form
C
C        ( Y(1,1) .     Y(1,k) . . . Y(1,j) . . Y(1,i)   )
C        (           .
C        (              .
C        (              Y(k,k) . . . Y(k,j) . . Y(k,i)   )
C        (                 .
C        (                         .
C        (                                .     Y(i-1,i) )
C        (                                      Y(i,i)   )
C        (   0                                     0     )
C
C     by comparing coefficients of equal power of s in the formula
C             i            j                 k
C     A(k,k) SUM Y(k,j) * s  = Y(k+1,k+1) * s  +
C            j=k
C     i-1                j                    j    i                   i
C     SUM (Y(k+1,j+1) - SUM A(k,l) Y(l,j)) * s  - SUM A(k,l) Y(l,i) * s
C    j=k+1             l=k+1                     l=k+1
C
C     Y(k,j), j = k,..., i, is a vector of length MU(k) - MU(k-1).
C
C     CONTRIBUTOR
C
C        A.J. Geurts.
C
C     REVISIONS
C
C        1992, October 27.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER NP, NA, I, NYI, LDA, LDY, IERR
C     .. Array Arguments ..
      INTEGER MU(*), S(*), SK(*)
      DOUBLE PRECISION YI(*), A(LDA,*), Y(LDY,*), RWORK(*)
C     .. Local Scalars ..
      INTEGER J, K, M, INDEXK, MUK, MUK1, NUK, MUK1NP
C     .. External Subroutines ..
      EXTERNAL COMPYI, DCOPY, DGEMV, DLAVEC, DLASET
C
C     .. Executable Statements ..
C
      IERR = 0
C
C     Initialization of the polynomial null vector Y(s).
C
      CALL DLASET('G', NA, I, ZERO, ZERO, Y, LDY)
      K = I
      IF (I .EQ. 1) THEN
         CALL DCOPY(NYI, YI, 1, Y(1,1), 1)
      ELSE
         INDEXK = NP + MU(I-1)
         CALL DCOPY(NYI, YI, 1, Y(INDEXK,I), 1)
         DO 30 K = I - 1, 1, -1
            MUK = MU(K)
            INDEXK = NP + MUK
            NUK = MU(K+1) - MUK
            DO 20 J = K, I
C
C              Compute the righthandside for Y(k,j) and store the
C              result in RWORK.
C
               IF (J .LT. I) THEN
                  CALL DCOPY(NUK, Y(INDEXK,J+1), 1, RWORK, 1)
               ELSE
                  CALL DLAVEC(NUK, ZERO, RWORK, 1)
               END IF
               IF (J .GT. K) THEN
                  CALL DGEMV('N', NUK, MU(J)-MUK, -ONE, A(MUK,INDEXK),
     *                       LDA, Y(INDEXK,J), 1, ONE, RWORK, 1)
               END IF
C
C              Solve A(k,k) * Y(k,j) = RWORK for Y(k,j).
C
               IF (K .EQ. 1) THEN
                  CALL COMPYI(NUK, NP, A(1,1), LDA, S, RWORK, Y(1,J),
     *                       IERR)
               ELSE
                  MUK1 = MU(K-1)
                  MUK1NP = MUK1 + NP
                  DO 10 M = 1, NUK
                     SK(M) = S(MUK-1+M) - MUK1NP + 1
   10             CONTINUE
                  CALL COMPYI(NUK, MUK-MUK1, A(MUK,MUK1NP), LDA, SK,
     *                       RWORK, Y(MUK1NP,J), IERR)
               END IF
               IF (IERR .NE. 0) THEN
                  IERR = K
                  RETURN
               END IF
   20       CONTINUE
   30    CONTINUE
      END IF
      RETURN
C *** Last line of COMPTY ***
      END






C
      SUBROUTINE COMPYI(M, N, A, LDA, S, V, Y, IERR)
C
C     PURPOSE
C
C     To compute a null vector yi of the right invertible diagonal
C     submatrix A(i,i), by solving an appropriate M by N system of
C     linear equations A y = v, where A is in staircase form.
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C     M - INTEGER.
C         The number of rows of matrix A.
C         M >= 1.
C     N - INTEGER.
C         The number of columns of matrix A.
C         N >= M.
C     A - DOUBLE PRECISION array of DIMENSION (LDA,N).
C         The leading M by N part of this array must contain the
C         matrix A.
C     LDA - INTEGER.
C         The leading dimension of array A as declared by the calling
C         program.
C         LDA >= M.
C     S - INTEGER array of DIMENSION at least (M).
C         S(i), i = 1,...,M, must contain the column index of the corner
C         in the i-th row of A.
C     V - DOUBLE PRECISION array of DIMENSION at least (M).
C         The righthand-side of the system of linear equations.
C
C     ARGUMENTS OUT
C
C     Y - DOUBLE PRECISION array of DIMENSION at least (N).
C         The computed solution of the system of linear equation.
C
C     ERROR INDICATOR
C
C     IERR - INTEGER.
C         Unless the routine detects an error (see next section),
C         IERR contains 0 on exit.
C
C     WARNINGS AMD ERRORS DETECTED BY THE ROUTINE
C
C     IERR = 3 : The matrix A is not right invertible.
C
C     METHOD
C
C     Let A * P = ( B | Z ) where P is a permutation matrix such that B
C     is nonsingular upper triangular. Z contains the remaining columns
C     of A. Then the system B x = v is solved and y = P ( x | 0 )'.
C
C     CONTRIBUTOR
C
C        A.J. Geurts.
C
C     REVISIONS
C
C        1992, October 27.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
C     .. Scalar Arguments ..
      INTEGER M, N, LDA, IERR
C     .. Array Arguments ..
      INTEGER S(M)
      DOUBLE PRECISION A(LDA,N), V(M), Y(N)
C     .. Local Scalars ..
      INTEGER J, K, SI
      DOUBLE PRECISION SUM
      LOGICAL FAIL
C     .. External Subroutines ..
      EXTERNAL DLARDV, DLAVEC
      DOUBLE PRECISION DLARDV
C
C     .. Executable Statements ..
C
C     Check input parameters.
C
      IF (N .LT. M) THEN
         IERR = 3
         RETURN
      END IF
C
      IERR = 0
      CALL DLAVEC(N, ZERO, Y, 1)
      SI = S(M)
      Y(SI) = DLARDV(V(M), A(M,SI), FAIL)
      IF (FAIL) THEN
         IERR = 3
         RETURN
      END IF
      DO 20 K = M - 1, 1, -1
         SUM = V(K)
         DO 10 J = K + 1, M
            SI = S(J)
            SUM = SUM - A(K,SI) * Y(SI)
   10    CONTINUE
         SI = S(K)
         Y(SI) = DLARDV(SUM, A(K,SI), FAIL)
         IF (FAIL) THEN
            IERR = 3
            RETURN
         END IF
   20 CONTINUE
      RETURN
C *** Last line of COMPYI ***
      END







      DOUBLE PRECISION FUNCTION DLARDV(A, B, FAIL)
C
C  -- LAPACK like auxiliary routine
C     December 8, 1994
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
      LOGICAL                           FAIL
C     ..
C
C  Purpose
C  =======
C
C  DLARDV returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
C
C  where  flmax  is a large value, via the function name. In addition if
C  a/b would overflow then  fail is returned as true, otherwise  fail is
C  returned as false.
C
C  Note that when  a and b  are both zero, fail is returned as true, but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that  abs( div ) = flmax.
C
C  When  b = 0  then  sign( a/b )  is taken as  sign( a ).
C
C  NB. This routine is derived from the Nag Fortran 77 O( 1 ) basic
C      linear algebra routine routine F06BLF.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE, ZERO
      PARAMETER             ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSB, DIV, FLMAX, FLMIN
      LOGICAL               FIRST
C     .. External Functions ..
      DOUBLE PRECISION      DLAMCH
      EXTERNAL              DLAMCH
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SIGN
C     .. Save statement ..
      SAVE                  FIRST, FLMIN, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
C
      IF( A.EQ.ZERO )THEN
         DIV = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            FLMIN =  DLAMCH('S')
            FLMAX =  1/FLMIN
         END IF
C
         IF( B.EQ.ZERO )THEN
            DIV  =  SIGN( FLMAX, A )
            FAIL = .TRUE.
         ELSE
            ABSB = ABS( B )
            IF( ABSB.GE.ONE )THEN
               FAIL = .FALSE.
               IF( ABS( A ).GE.ABSB*FLMIN )THEN
                  DIV = A/B
               ELSE
                  DIV = ZERO
               END IF
            ELSE
               IF( ABS( A ).LE.ABSB*FLMAX )THEN
                  FAIL = .FALSE.
                  DIV  =  A/B
               ELSE
                  FAIL = .TRUE.
                  DIV  = FLMAX
                  IF( ( ( A.LT.ZERO ).AND.( B.GT.ZERO ) ).OR.
     $                ( ( A.GT.ZERO ).AND.( B.LT.ZERO ) )     )
     $               DIV = -DIV
               END IF
            END IF
         END IF
      END IF
C
      DLARDV = DIV
      RETURN
C
C     End of DLARDV
C
      END






      SUBROUTINE DLARFV( N, V, INCV, TAU, ALPHA, X, INCX )
C
C  -- LAPACK like auxiliary routine
C     December 8, 1994
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TAU
      INTEGER            INCX, INCV, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), V( * )
C     ..
C
C  Purpose
C  =======
C
C  DLARFV performs a Householder reflection given by
C
C     ( alpha ) = H * ( alpha ) ,
C     (   x   )       (   x   )
C
C  where the orthogonal matrix H is given in the form
C
C     H = I - tau * ( 1 ) * ( 1 v') ,
C                   ( v )
C
C  where tau is a real scalar and v is a real (n-1)-element.
C  If tau is zero then H is assumed to be the unit matrix and the
C  transformation is skipped, otherwise tau must be in the range
C  ( 1.0, 2.0 ). tau and v will usually be supplied by routine DLARFG.
C
C  NB. This routine is derived from the Nag Fortran 77 O( n ) basic
C      linear agebra routine F06FUF.
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA
C     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
C     .. External Subroutines ..
      EXTERNAL           DAXPY
C     ..
C     .. Executable Statements ..
      IF( TAU.NE.ZERO ) THEN
         BETA  = ALPHA + DDOT( N, V, INCV, X, INCX )
         ALPHA = ALPHA - TAU * BETA
         CALL DAXPY( N-1, -TAU * BETA, V, INCV, X, INCX )
      END IF
C
      RETURN
C
C     End of DLARFV
C
      END






      SUBROUTINE DLAVEC( N, CONST, X, INCX )
C
C  -- LAPACK like auxiliary routine
C     August 31, 1994
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  Purpose
C  =======
C
C  DLAVEC performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C  NB. This routine is derived from the Nag Fortran 77 O( n ) basic
C      linear algebra routine F06FBF.
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C
C     .. Executable Statements ..
      IF( N.GT.0 ) THEN
         IF( CONST.NE.ZERO ) THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of DLAVEC
C
      END





C
      SUBROUTINE HHTRAN(MA, NA, L, K, A, LDA, Q, RWORK)
C
C     PURPOSE
C
C     To compute the Householder reflection Q which transforms a
C     subcolumn of A into the first unit vector and to apply Q left and
C     right to appropriate submatrices of A.
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C        MA - INTEGER.
C            The number of rows of matrix A.
C            MA >= 1.
C        NA - INTEGER.
C            The number of columns of matrix A.
C            NA >= MA.
C        L - INTEGER.
C            The index of the column of A to be transformed.
C            L >= 1.
C        K - INTEGER.
C            The row index from which the column is to be transformed.
C            K >= 1.
C        A - DOUBLE PRECISION array of DIMENSION (LDA,NA).
C            The leading MA by NA part of this array must contain the
C            matrix A. The left lower MA - K + 1 by L - 1 block of the
C            matrix A is understood to be zero made by former
C            transformations.
C            Note: this array is overwritten.
C        LDA - INTEGER.
C            The leading dimension of array A as declared in the calling
C            program.
C            LDA >= MA.
C
C     ARGUMENTS OUT
C
C        A - DOUBLE PRECISION array of DIMENSION (LDA,NA).
C            The leading MA by NA part of this array contains the
C            transformed matrix Q A Q , where Q  is the Householder
C                                1   2         i
C            transformation appropriately augmented with an identity
C            matrix.
C        Q - DOUBLE PRECISION array of DIMENSION at least (MA-K+1).
C            The leading MA - K + 1 elements of this array contain the
C            constant tau and the vector v which define the Householder
C            reflection (see DLARGF, i.e., Q(1) contains the value tau
C            and Q(i), i=2,...MA-K+1 contain the vector v).
C
C     WORKSPACE
C
C        RWORK - DOUBLE PRECISION array of DIMENSION at least (MA-K+1).
C
C     METHOD
C
C     Let w be the subcolumn A(i,L), i=K,...,MA. Then the Householder
C     reflection P = I - tau * u * u', where u = (1, v')', such that
C     Pw = e1 (the first unit vector) is computed.
C     Let
C         Q  = ( I  0 )       Q  = ( I  0 )
C          1   ( 0  P )        2   ( 0  P )
C     be such that Q  is MA by MA and Q  is NA by NA.
C                   1                2
C     Then Q  A Q   is computed and stored in A.
C           1    2
C
C     CONTRIBUTOR
C
C        A.J. Geurts.
C
C     REVISIONS
C
C        1992, October 27.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER MA, NA, K, L, LDA
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,NA), Q(MA-K+1), RWORK(MA-K+1)
C     .. Local Scalars ..
      INTEGER LEN, NMK, NL1
      DOUBLE PRECISION TAU
C     .. External Subroutines/Functions ..
      EXTERNAL DCOPY, DGEMV, DGER, DLARFG, DLAVEC
C
C     .. Executable Statements ..
C
      LEN = MA - K + 1
      IF (LEN .GT. 1) THEN
C
C        Generate the Householder transformation.
C
         CALL DCOPY(LEN, A(K,L), 1, Q, 1)
         CALL DLARFG(LEN, Q(1), Q(2), 1, TAU)
C
         IF (TAU .NE. ZERO) THEN
C
C           Compute matrix Q A by w = A'u and next A = A - tau * uw'.
C                           1
            Q(1) = ONE
            NL1 = NA - L + 1
            CALL DLAVEC(MA, ZERO, RWORK, 1)
            CALL DGEMV('T', LEN, NL1, ONE, A(K,L), LDA, Q, 1, ZERO,
     *                  RWORK, 1)
            CALL DGER(LEN, NL1, -TAU, Q, 1, RWORK, 1, A(K,L), LDA)
C
C           Compute matrix Q A Q  by w = Au and next A = A - tau * wu'.
C                           1   2
            NMK = NA - MA + K
            CALL DGEMV('N', MA, LEN, ONE, A(1,NMK), LDA, Q, 1, ZERO,
     *                  RWORK, 1)
            CALL DGER(MA, LEN, -TAU, RWORK, 1, Q, 1, A(1,NMK), LDA)
            Q(1) = TAU
         END IF
      ELSE
         Q(1) = ZERO
      END IF
      RETURN
C *** Last line of HHTRAN ***
      END






C
      SUBROUTINE MKPENC(MP, NP, DP, P, LDP1, LDP2, A, LDA, MA, NA)
C
C     PURPOSE
C
C     Given an MP x NP polynomial matrix of degree dp
C                                    dp-1            dp
C     P(s) = P(0) + ... + P(dp-1) * s     + P(dp) * s            (1)
C
C     the subroutine MKPENC constructs the first degree part A of the
C     linearization of the polynomial matrix P(s), where
C
C            | P(dp)    0  .        .   0 |
C            | P(dp-1)  I  0        .   0 |
C        A = |   .      0  I  .     .   0 |                      (2)
C            |   .            .  .  .   . |
C            |   .               I  0   0 |
C            | P(0)     0           I   0 |
C
C     REMARK: This auxiliary routine is intended to be called only from
C             the routine COLRED.
C
C     ARGUMENTS IN
C
C        MP - INTEGER.
C            The number of rows of the polynomial matrix P(s).
C            MP >= 1.
C        NP - INTEGER.
C            The number of columns of the polynomial matrix P(s).
C            NP >= 1.
C        DP - INTEGER.
C            The degree of the polynomial matrix P(s).
C            DP >= 1.
C        P - DOUBLE PRECISION array of DIMENSION (LDP1,LP2,DP+1).
C            The leading MP by NP by (DP+1) part of this array must
C            contain the coefficients of the polynomial matrix P(s).
C            Specifically, P(i,j,k) must contain the coefficient of
C            s**(k-1) of the polynomial which is the (i,j)-th element
C            of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C            k = 1,2,...,DP+1.
C        LDP1 - INTEGER.
C            The leading dimension of array P as declared in the calling
C            program.
C            LDP1 >= MP.
C        LDP2 - INTEGER.
C            The second dimension of array P as declared in the calling
C            program.
C            LDP2 >= NP.
C
C     ARGUMENTS OUT
C
C        A - DOUBLE PRECISION array of DIMENSION (LDA,(DP+1)*MP+NP).
C            The leading (DP+1)*MP by (DP+1)*MP+NP part of this array
C            contains the matrix A as described in (2).
C        LDA - INTEGER.
C            The leading dimension of array A as declared in the calling
C            program.
C            LDA >= (DP+1)*MP.
C        MA - INTEGER.
C            The number of rows of matrix A.
C        NA - INTEGER.
C            The number of columns of matrix A.
C
C     CONTRIBUTOR
C
C        A.J.Geurts.
C
C     REVISIONS
C
C        1992, October 27.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER MP, NP, DP, LDP1, LDP2, LDA, MA, NA
C     .. Array Arguments ..
      DOUBLE PRECISION P(LDP1,LDP2,*), A(LDA,*)
C     .. Local Scalars ..
      INTEGER J, M1, NJ
C     .. External Subroutines/Functions ..
      EXTERNAL DLACPY, DLASET
C
C     .. Executable Statements ..
C
C     Initialisation of the matrix A.
C
      M1 = DP * MP
      MA = M1 + MP
      NA = MA + NP
      CALL DLASET('G', MP, MA, ZERO, ZERO, A(1,NP+1), LDA)
      CALL DLASET('G', M1, MA, ZERO, ONE, A(MP+1,NP+1), LDA)
C
C     Insert the matrices P(0), P(1), ..., P(pd) in array A.
C
      NJ = M1 + 1
      DO 20 J = 1, DP + 1
         CALL DLACPY('G', MP, NP, P(1,1,J), LDP1, A(NJ,1), LDA)
         NJ = NJ - MP
   20 CONTINUE
C
      RETURN
C *** Last line of MKPENC ***
      END


