C
C     COLRED TEST DRIVER TEXT
C
C     Version dd February 24, 1995.
C     Revised May 21, 1996.
C
C     IMPLICIT NONE
C     .. Parameters ..
      INTEGER NIN, NOUT
      PARAMETER (NIN = 5, NOUT = 6)
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
C
      INTEGER DPMAX, MPMAX, NPMAX
      PARAMETER (DPMAX = 6, MPMAX = 8, NPMAX = 10)
      INTEGER MAMAX, NAMAX, LIWORK, LRWORK, DRMAX
      PARAMETER (MAMAX = (NPMAX * DPMAX + 1) * MPMAX, 
     *           NAMAX = MAMAX + NPMAX,
     *           LIWORK = 2 * MAMAX + 2 * NPMAX + 1,
     *           LRWORK = (2*MAMAX + NPMAX*DPMAX + 5) * NAMAX +
     *                    MAMAX + (2*MPMAX + 1) * NPMAX,
     *            DRMAX = (NPMAX+1)*DPMAX + 1)
      INTEGER LDP1, LDP2, LDR1, LDR2, LDU1, LDU2
      PARAMETER (LDP1 = MPMAX, LDP2 = NPMAX, LDR1 = MPMAX,
     *           LDR2 = NPMAX, LDU1 = NPMAX, LDU2 = NPMAX)
C     .. Local Scalars ..
      INTEGER MP, NP, DP, DR, DU, I, J, K, L, NUM, IERR
      DOUBLE PRECISION TOL
C     .. Local Arrays ..
      INTEGER IWORK(LIWORK)
      DOUBLE PRECISION P(LDP1,LDP2,DPMAX+1), R(LDR1,LDR2,DRMAX),
     *                 U(LDU1,LDU2,NPMAX*DPMAX+1), RWORK(LRWORK)
      LOGICAL ZERCOL(NPMAX)
C     .. External Subroutines/Functions ..
      EXTERNAL COLRED, MULTPM, PRMAPO
C
C     .. Executable Statements ..
C
C     Read NUM, that is the number of testcases to be run.
C
      READ(NIN, FMT = *) NUM
C
      DO 30 L = 1, NUM
C
         WRITE(NOUT, FMT = 99999) L
         READ(NIN, FMT = '()')
C
C        Read MP, NP, DP, TOL and next P(k), k = 0,....,DP row after row.
C
         READ(NIN, FMT = *) MP, NP, DP, TOL
C
         DO 20 K = 1, DP + 1
            READ(NIN, FMT = '()')
            DO 10 I = 1, MP
               READ(NIN, FMT = *) (P(I,J,K), J = 1, NP)
 10         CONTINUE
 20      CONTINUE
C
         WRITE(NOUT, FMT = 99998) DP, MP, NP, TOL
         CALL PRMAPO(MP, NP, DP, 5, NOUT, P, LDP1, LDP2, 'P', IERR)
C
         CALL COLRED(MP, NP, DP, P, LDP1, LDP2, DR, DU, R, LDR1, LDR2,
     *               U, LDU1, LDU2, ZERCOL, IWORK, RWORK, TOL, IERR)
C
         IF (IERR .EQ. 0) THEN
            WRITE (NOUT, FMT = 99997)
            CALL PRMAPO(NP, NP, DU, 5, NOUT, U, LDU1, LDU2, 'U', IERR)
C
            WRITE (NOUT, FMT = 99996)
            CALL PRMAPO(MP, NP, DR, 5, NOUT, R, LDR1, LDR2, 'R', IERR)
            WRITE (NOUT, FMT = 99995) (ZERCOL(J), J = 1, NP)
C
            CALL MULTPM(-ONE, MP, NP, NP, DP, DU, DR, P, LDP1, LDP2,
     *                  U, LDU1, LDU2, R, LDR1, LDR2, RWORK, IERR)
            IF (DR .GE. 0) THEN
               WRITE (NOUT, FMT = 99994)
               CALL PRMAPO(MP, NP, DR, 5, NOUT, R, LDR1, LDR2,
     *                     '(PU - R)', IERR)
            ELSE
               WRITE (NOUT, FMT = 99993)
            END IF
         ELSE
            WRITE (NOUT, FMT = 99992) IERR
         END IF
 30   CONTINUE
      STOP
C
99999 FORMAT (//' COLRED TEST DRIVER RESULTS, EXAMPLE(', I2, ')', /1X)
99998 FORMAT (' The input polynomial matrix:', //,
     *        ' P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s**(dp-1)',
     *        '  + P(dp) * s**dp', //, ' with degree DP =', I2,
     *        ' and size MP =', I2, ', NP =', I2, '.', //,
     *        ' The tolerance is:', D10.3, /1X)
99997 FORMAT (' The unimodular polynomial matrix U(s):')
99996 FORMAT (' The column reduced polynomial matrix R(s):')
99995 FORMAT (' ZERCOL(j), j = 1, NP:', 10(L2))
99994 FORMAT (' The residual matrix P(s) * U(s) - R(s):')
99993 FORMAT (' PU - R is the ZERO polynomial matrix.')
99992 FORMAT (' COLRED has failed: IERR =', I2)
      END






C
      SUBROUTINE MULTPM(ALPHA, RP1, CP1, CP2, DP1, DP2, DP3, P1, LDP11,
     *                  LDP12, P2, LDP21, LDP22, P3, LDP31, LDP32,
     *                  RWORK, IERR)
C
C     PURPOSE
C
C     To compute the coefficients of the real polynomial matrix
C
C        P(s) = P1(s) * P2(s) + alpha * P3(s),                      (1)
C
C     where P1(s), P2(s) and P3(s) are real polynomial matrices and
C     alpha is a real scalar.
C     Each of the polynomial matrices may by the zero matrix.
C
C     ARGUMENTS IN
C
C        ALPHA - DOUBLE PRECISION.
C           The value of the scalar factor alpha of the problem.
C        RP1 - INTEGER.
C           The number of rows of the matrices P1(s) and P3(s).
C           RP1 >= 1.
C        CP1 - INTEGER.
C           The number of columns of matrix P1(s) and the number of rows
C           of matrix P2(s).
C           CP1 >= 1.
C        CP2 - INTEGER.
C           The number of columns of the matrices P2(s) and P3(s).
C           CP2 >= 1.
C        DP1 - INTEGER.
C           The degree of the polynomial matrix P1(s).
C           DP1 >= -1.
C        DP2 - INTEGER.
C           The degree of the polynomial matrix P2(s).
C           DP2 >= -1.
C        DP3 - INTEGER.
C           The degree of the polynomial matrix P3(s).
C           DP3 >= -1.
C           Note: DP3 is overwritten.
C        P1 - DOUBLE PRECISION array of (LDP11,LDP12,*).
C           If DP1 >= 0, then the leading RP1 by CP1 by (DP1+1) part of
C           this array must contain the coefficients of the polynomial
C           matrix P1(s). Specifically, P1(i,j,k) must contain the
C           coefficient of s**(k-1) of the polynomial which is the
C           (i,j)-th element of P1(s), where i = 1,2, ...,RP1,
C           j = 1,2,...,CP1 and k = 1,2,...,DP1+1.
C           If DP1 = -1, then P1(s) is taken to be the zero polynomial
C           matrix and the array P1 is not referenced.
C        LDP11 - INTEGER.
C           The leading dimension of array P1 as declared in the calling
C           program.
C           LDP11 >= RP1 if DP1 >= 0.
C           LDP11 >= 1 if DP1 = -1.
C        LDP12 - INTEGER.
C           The second dimension of array P1 as declared in the calling
C           program.
C           LDP12 >= CP1 if DP1 >= 0,
C           LDP12 >= 1 if DP1 = -1.
C        P2 - DOUBLE PRECISION array of (LDP21,LDP22,*).
C           If DP2 >= 0, then the leading CP1 by CP2 by (DP2+1) part of
C           this array must contain the coefficients of the polynomial
C           matrix P2(s). Specifically, P2(i,j,k) must contain the
C           coefficient of s**(k-1) of the polynomial which is the
C           (i,j)-th element of P2(s), where i = 1,2, ...,CP1,
C           j = 1,2,...,CP2 and k = 1,2,...,DP2+1.
C           If DP2 = -1, then P2(s) is taken to be the zero polynomial
C           matrix and the array P2 is not referenced.
C        LDP21 - INTEGER.
C           The leading dimension of array P2 as declared in the calling
C           program.
C           LDP21 >= CP1 if DP2 >= 0.
C           LDP21 >= 1 if DP2 = -1.
C        LDP22 - INTEGER.
C           The second dimension of array P2 as declared in the calling
C           program.
C           LDP22 >= CP2 if DP2 >= 0,
C           LDP22 >= 1 if DP2 = -1.
C        P3 - DOUBLE PRECISION array of (LDP31,LDP32,lenp3),
C           where lenp3 = MAX(DP1+DP2,DP3,0) + 1.
C           If DP3 >= 0, then the leading RP1 by CP2 by (DP3+1) part of
C           this array must contain the coefficients of the polynomial
C           matrix P3(s). Specifically, P3(i,j,k) must contain the
C           coefficient of s**(k-1) of the polynomial which is the
C           (i,j)-th element of P3(s), where i = 1,2, ...,RP1,
C           j = 1,2,...,CP2 and k = 1,2,...,DP3+1.
C           If DP3 = -1, then P3(s) is taken to be the zero polynomial
C           matrix.
C           Note: this array is overwritten.
C        LDP31 - INTEGER.
C           The leading dimension of array P3 as declared in the calling
C           program.
C           LDP31 >= RP1.
C        LDP32 - INTEGER.
C           The second dimension of array P3 as declared in the calling
C           program.
C           LDP32 >= CP2.
C
C     ARGUMENTS OUT
C
C        DP3 - INTEGER.
C           The degree of the resulting polynomial matrix P(s).
C        P3 - DOUBLE PRECISION array of DIMENSION (LDP31,LDP32,lenp3).
C           If DP3 >= 0 on exit, then the leading RP1 by CP2 by (DP3+1)
C           part of this array contains the coefficients of P(s).
C           Specifically, P3(i,j,k) contains the coefficient of s**(k-1)
C           of the polynomial which is the (i,j)-th element of P(s),
C           where i = 1,2, ...,RP1, j = 1,2,...,CP2 and
C           k = 1,2,...,DP3+1.
C           If DP3 = -1 on exit, then P(s) is the zero polynomial
C           matrix and the contents of the array P3 are undefined.
C
C     WORK SPACE
C
C        RWORK - DOUBLE PRECISION array of DIMENSION at least (CP1).
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
C     Given the real polynomial matrices
C                DP1            i           DP2            i
C        P1(s) = SUM (a(i+1) * s ), P2(s) = SUM (b(i+1) * s ),
C                i=0                        i=0
C                DP3            i
C        P3(s) = SUM (c(i+1) * s ).
C                i=0
C     and a real scalar alpha, the routine computes the coefficients
C     d(1), d(2), ... of the polynomial matrix (1) from the formula
C                   s
C        d(i+1) := SUM (a(k+1) * b(i-k+1)) + alpha * c(i+1),
C                  k=r
C     where i = 0,1,...,DP1+DP2 and r and s depend on the value of i,
C     i.e. for r <= k <= s both a(k+1) and b(i-k+1) must exist.
C
C     CONTRIBUTOR
C
C        A.J. Geurts (Eindhoven University of Technology).
C
C     REVISIONS
C
C        1992, October 28.
C        1994, December 8.
C
C     IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
C     .. Scalar Arguments ..
      INTEGER RP1, CP1, CP2, DP1, DP2, DP3, LDP11, LDP12,
     *        LDP21, LDP22, LDP31, LDP32, IERR
      DOUBLE PRECISION ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION P1(LDP11, LDP12, *), P2(LDP21, LDP22, *),
     *                 P3(LDP31, LDP32, *), RWORK(*)
C     .. Local Scalars ..
      INTEGER H, I, J, K, DPOL3, E
      DOUBLE PRECISION W
      LOGICAL CFZERO
C     .. External Subroutines/Functions ..
      EXTERNAL DCOPY, DDOT, DLAVEC, DSCAL
      DOUBLE PRECISION DDOT
C
C     .. Executable Statements ..
C
C     Check input parameters.
C
      IF ((RP1.LT.1) .OR. (CP1.LT.1) .OR. (CP2.LT.1)
     *   .OR. (DP1.LT.-1) .OR. (DP2.LT.-1) .OR. (DP3.LT.-1)
     *   .OR. ((LDP11.LT.RP1) .AND. (DP1.GE.0))
     *   .OR. ((LDP11.LT.1) .AND. (DP1.EQ.-1))
     *   .OR. ((LDP12.LT.CP1) .AND. (DP1.GE.0))
     *   .OR. ((LDP12.LT.1) .AND. (DP1.EQ.-1))
     *   .OR. ((LDP21.LT.CP1) .AND. (DP2.GE.0))
     *   .OR. ((LDP21.LT.1) .AND. (DP2.EQ.-1))
     *   .OR. ((LDP22.LT.CP2) .AND. (DP2.GE.0))
     *   .OR. ((LDP22.LT.1) .AND. (DP2.EQ.-1))
     *   .OR. (LDP31.LT.RP1) .OR. (LDP32.LT.CP2)) THEN
         IERR = 1
         RETURN
      END IF
C
      IERR = 0
      IF (ALPHA .EQ. ZERO) THEN
         DP3 = -1
      END IF
C
      IF (DP3 .GE. 0) THEN
C
C        P3(s) := ALPHA * P3(s).
C
         DO 20 K = 1, DP3 + 1
            DO 10 J = 1, CP2
               CALL DSCAL(RP1, ALPHA, P3(1,J,K), 1)
   10       CONTINUE
   20    CONTINUE
      END IF
C
      IF ((DP1 .EQ. -1) .OR. (DP2 .EQ. -1)) RETURN
C
C     Neither of P1(s) and P2(s) is the zero polynomial.
C
      DPOL3 = DP1 + DP2
      IF (DPOL3 .GT. DP3) THEN
C
C        Initialize the additional part of P3(s) to zero.
C
         DO 40 K = DP3 + 2, DPOL3 + 1
            DO 30 J = 1, CP2
               CALL DLAVEC(RP1, ZERO, P3(1,J,K), 1)
   30       CONTINUE
   40    CONTINUE
         DP3 = DPOL3
      END IF
C                                                              k-1
C     The inner product of the j-th row of the coefficient of s    of
C                                                      i-1
C     P1(s) and the h-th column of the coefficient of s    of P2(s)
C                                                                k+i-2
C     contributes to the (j,h)-th element of the coefficient of s
C     of P3(s).
C
      DO 80 K = 1, DP1 + 1
         DO 70 J = 1, RP1
            CALL DCOPY(CP1, P1(J,1,K), LDP11, RWORK, 1)
            DO 60 I = 1, DP2 + 1
               E = K + I - 1
               DO 50 H = 1, CP2
                  W = DDOT(CP1, RWORK, 1, P2(1,H,I), 1)
                  P3(J,H,E) = W + P3(J,H,E)
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
C
C     Computation of the exact degree of P3(s).
C
      CFZERO = .TRUE.
C     WHILE (DP3 >= 0 and CFZERO) DO
   90 IF ((DP3 .GE. 0) .AND. CFZERO) THEN
         DPOL3 = DP3 + 1
         DO 110 I = 1, RP1
            DO 100 J = 1, CP2
               IF (P3(I,J,DPOL3) .NE. ZERO) CFZERO = .FALSE.
  100       CONTINUE
  110    CONTINUE
         IF (CFZERO) DP3 = DP3 - 1
         GO TO 90
      END IF
C     END WHILE 90
C
      RETURN
C *** Last line of MULTPM ***
      END






C
      SUBROUTINE PRMAPO(MP, NP, DP, L, NOUT, P, LDP1, LDP2, TEXT, IERR)
C
C     PURPOSE
C
C     To print the MP by NP coefficient matrices of a matrix polynomial
C                                                  dp-1           dp
C        P(s) = P(0) + P(1) * s + . . . P(dp-1) * s    + P(dp) * s  .
C
C     The elements of the matrices are output to 7 significant figures.
C
C     ARGUMENTS IN
C
C        MP - INTEGER.
C            The number of rows of the matrix polynomial P(s).
C            MP >= 1.
C        NP - INTEGER.
C            The number of columns of the matrix polynomial P(s).
C            NP >= 1.
C        DP - INTEGER.
C            The degree of the matrix polynomial P(s).
C            DP >= 0.
C        L - INTEGER.
C            The number of elements of the coefficient matrices to be
C            printed per line.
C            1 <= L <= 5.
C        NOUT - INTEGER.
C            The output channel to which the results are sent.
C            NOUT >= 0.
C        P - DOUBLE PRECISION array of DIMENSION (LDP1,LDP2,DP+1).
C            The leading MP by NP by (DP+1) part of this array must
C            contain the coefficients of the matrix polynomial P(s).
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
C        TEXT - CHARACTER*72.
C            Title caption of the coefficient matrices to be printed.
C            TEXT is followed by the degree of the coefficient matrix,
C            within brackets. If TEXT = ' ', then the coefficient
C            matrices are separated by an empty line.
C
C     ERROR INDICATOR
C
C        IERR - INTEGER.
C            Unless the routine detects an error (see next section),
C            IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = 1 : Invalid input parameter(s).
C
C     METHOD
C
C     For i = 1, 2, ..., DP + 1 the routine first prints the contents of
C     TEXT followed by (i-1) as a title, followed by the elements of the
C     MP by NP coefficient matrix P(i) such that
C     (i)  if NP < L, then the leading MP by NP part is printed;
C     (ii) if NP = k*L + p (where k, p > 0), then k MP by L blocks of
C          consecutive columns of P(i) are printed one after another
C          followed by one MP by p block containing the last p columns of
C          P(i).
C     Row numbers are printed on the left of each row and a column
C     number on top of each column.
C
C     CONTRIBUTOR
C
C        A.J. Geurts (Eindhoven University of Technology).
C
C     REVISIONS
C
C        1992, October 28.
C        1995, February 6.
C
C     IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER MP, NP, DP, L, NOUT, LDP1, LDP2, IERR
      CHARACTER*(*) TEXT
C     .. Array Arguments ..
      DOUBLE PRECISION P(LDP1,LDP2,*)
C     .. Local Scalars ..
      INTEGER I, J, J1, J2, JJ, K, LENTXT, LTEXT, N1
C     .. Intrinsic Functions ..
      INTRINSIC LEN, MIN
C
C     .. Executable Statements ..
C
C     Check input parameters.
C
      LENTXT = LEN(TEXT)
      LTEXT = MIN(72,LENTXT)
C     WHILE (TEXT(LTEXT:LTEXT) =  ' ') DO
   10 IF (TEXT(LTEXT:LTEXT) .EQ. ' ') THEN
         LTEXT = LTEXT - 1
         GO TO 10
      END IF
C     END WHILE 10
C
      IF (MP.LT.1 .OR. NP.LT.1 .OR. DP.LT.0 .OR. L.LT.1 .OR. L.GT.5
     *   .OR. NOUT.LT.0 .OR. LDP1.LT.MP .OR. LDP2.LT.NP) THEN
         IERR = 1
         RETURN
      END IF
C
      IERR = 0
C
      DO 50 K = 1, DP + 1
         IF (LTEXT .EQ. 0) THEN
            WRITE (NOUT, FMT = 99999)
         ELSE
            WRITE (NOUT, FMT = 99998) TEXT(1:LTEXT), K - 1
         END IF
         N1 = (NP - 1)/L
         J1 = 1
         J2 = L
         DO 30 J = 1, N1
            WRITE (NOUT, FMT = 99997) (JJ, JJ = J1, J2)
            DO 20 I = 1, MP
               WRITE (NOUT, FMT = 99996) I, (P(I,JJ,K), JJ = J1, J2)
   20       CONTINUE
            J1 = J1 + L
            J2 = J2 + L
   30    CONTINUE
         WRITE (NOUT, FMT = 99997) (J, J = J1, NP)
         DO 40 I = 1, MP
            WRITE (NOUT, FMT = 99996) I, (P(I,JJ,K), JJ = J1, NP)
   40    CONTINUE
   50 CONTINUE
      WRITE (NOUT, FMT = 99999)
C
      RETURN
99999 FORMAT (' ')
99998 FORMAT (/, 1X, A, '(', I2, ')')
99997 FORMAT (5X, 5(6X, I2, 7X))
99996 FORMAT (1X, I2, 2X, 5F15.7)
C *** Last line of PRMAPO ***
      END





