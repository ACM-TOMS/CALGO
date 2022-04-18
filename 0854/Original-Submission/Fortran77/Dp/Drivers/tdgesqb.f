*     DGESQB/DOSGSB EXAMPLE PROGRAM TEXT
*
*     .. Parameters ..
      IMPLICIT NONE
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
      INTEGER          NIN, NOUT
      PARAMETER        ( NIN = 5, NOUT = 6 )
      INTEGER          MMAX, NBMAX, NMAX
      PARAMETER        ( MMAX = 200, NBMAX = 64, NMAX = 200 )
      INTEGER          LDA, LDB, LDQ1, LDQ2, LDRES, LDWORK
      PARAMETER        ( LDA = MMAX, LDB = MMAX, LDQ1 = MMAX,
     $                   LDQ2 = MMAX, LDRES = MMAX,
     $                   LDWORK = 8*NMAX*NBMAX + 15*NBMAX*NBMAX )
*     .. Local Scalars .. 
      INTEGER          I, INFO, J, K, M, N
      DOUBLE PRECISION TEMP
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA, NMAX), B(LDB, NMAX),  CS(2*NMAX),
     $                 DWORK(LDWORK), Q1(LDQ1, MMAX), Q2(LDQ2,MMAX),
     $                 RES(LDRES,MMAX+2*NMAX),  TAU(NMAX)
*     .. External Subroutines ..
      EXTERNAL         DGEMM, DGESQB, DLACPY, DLASET, DOSGSB
*     .. External Functions .
      DOUBLE PRECISION DLANGE, DLAORS, DLAPY2
      EXTERNAL         DLANGE, DLAORS, DLAPY2
*     .. Executable Statements ..
      WRITE ( NOUT, FMT = 99999 )
*     Skip the heading in the data file and read the data.
      READ ( NIN, FMT = '()' )
      READ ( NIN, FMT = * )  M, N
      IF( M.LE.0 .OR. M.GT.MMAX ) THEN
         WRITE ( NOUT, FMT = 99993 ) M
      ELSE IF( N.LE.0 .OR. N.GT.NMAX ) THEN
         WRITE ( NOUT, FMT = 99992 ) N
      ELSE
         READ ( NIN, FMT = * ) ( ( A(I,J), J = 1,N ), I = 1,M )
         CALL DLACPY( 'All', M, N, A, LDA, RES, LDRES )
         READ ( NIN, FMT = * ) ( ( B(I,J), J = 1,N ), I = 1,M )
         CALL DLACPY( 'All', M, N, B, LDB, RES(1,N+1), LDRES )
         CALL DGESQB( M, N, A, LDA, B, LDB, CS, TAU, DWORK, LDWORK,
     $                INFO )
         INFO = 0
         IF ( INFO.NE.0 ) THEN
            WRITE ( NOUT, FMT = 99998 ) INFO
         ELSE
            K = MIN(M,N)
            CALL DLACPY( 'Lower', M-1, K, A(2,1), LDA, Q1(2,1), LDQ1 )
            CALL DLACPY( 'Lower', M, K, B, LDB, Q2, LDQ2 )
            CALL DOSGSB( 'No Transpose', 'No Transpose', M, M, K, Q1,
     $                   LDQ1, Q2, LDQ2, CS, TAU, DWORK, LDWORK, INFO )
            IF ( INFO.NE.0 ) THEN
               WRITE ( NOUT, FMT = 99997 ) INFO
            ELSE
               CALL DLASET( 'Lower', M-1, K, ZERO, ZERO, A(2,1), LDA )
               CALL DLASET( 'Lower', M, K, ZERO, ZERO, B, LDB )
C
               CALL DGEMM( 'No Transpose', 'No Transpose', M, N, M,
     $                     -ONE, Q1, LDQ1, A, LDA, ONE, RES, LDRES )
               CALL DGEMM( 'No Transpose', 'No Transpose', M, N, M, 
     $                     -ONE, Q2, LDQ2, B, LDB, ONE, RES, LDRES )
               CALL DGEMM( 'No Transpose', 'No Transpose', M, N, M,
     $                     ONE, Q2, LDQ2, A, LDA, ONE, RES(1,N+1),
     $                     LDRES )
               CALL DGEMM( 'No Transpose', 'No Transpose', M, N, M,
     $                     -ONE, Q1, LDQ1, B, LDB, ONE, RES(1,N+1),
     $                     LDRES )
               TEMP = DLAPY2( DLANGE( 'Frobenius', M, N, RES, LDRES,
     $                        DWORK ), DLANGE( 'Frobenius', M, N,
     $                        RES(1,N+1), LDRES, DWORK ) )
C
               WRITE ( NOUT, FMT = 99996 )
               DO 10  I = 1, M
                  WRITE (NOUT, FMT = 99994)
     $                  ( Q1(I,J), J = 1,M ), ( Q2(I,J), J = 1,M )
10             CONTINUE
               DO 20  I = 1, M
                  WRITE (NOUT, FMT = 99994)
     $                  ( -Q2(I,J), J = 1,M ), ( Q1(I,J), J = 1,M )
20             CONTINUE
               WRITE ( NOUT, FMT = 99991 ) DLAORS( .FALSE., .FALSE., M,
     $                  Q1, LDQ1, Q2, LDQ2, RES, LDRES ) 
               WRITE ( NOUT, FMT = 99995 )
               DO 30  I = 1, M
                  WRITE (NOUT, FMT = 99994) ( A(I,J), J = 1,N )
30             CONTINUE
               DO 40  I = 1, M
                  WRITE (NOUT, FMT = 99994) ( B(I,J), J = 1,N )
40             CONTINUE
               WRITE ( NOUT, FMT = 99990 ) TEMP
            END IF
         END IF
      END IF    
*
      STOP
*
99999 FORMAT (' TDGESQB EXAMPLE PROGRAM RESULTS',/1X)
99998 FORMAT (' INFO on exit from DGESQB = ',I2)
99997 FORMAT (' INFO on exit from DOSGSB = ',I2)
99996 FORMAT (' The orthogonal symplectic factor Q is ')
99995 FORMAT (/' The factor R is ')
99994 FORMAT (20(1X,F9.4))
99993 FORMAT (/' M is out of range.',/' M = ',I5)
99992 FORMAT (/' N is out of range.',/' N = ',I5)
99991 FORMAT (/' Orthogonality of Q: || Q''*Q - I ||_F = ',G7.2)
99990 FORMAT (/' Residual: || [A;B] - Q*R ||_F = ',G7.2)
      END
