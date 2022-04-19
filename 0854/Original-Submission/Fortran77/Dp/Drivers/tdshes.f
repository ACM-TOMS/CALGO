*     DSHES EXAMPLE PROGRAM TEXT
*
*     .. Parameters ..
      IMPLICIT NONE
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )      
      INTEGER          NIN, NOUT
      PARAMETER        ( NIN = 5, NOUT = 6 )
      INTEGER          NMAX
      PARAMETER        ( NMAX = 100 )
      INTEGER          LDA, LDQG, LDRES, LDU1, LDU2, LDWORK
      PARAMETER        ( LDA = NMAX, LDQG = NMAX, LDRES = NMAX,
     $                   LDU1 = NMAX, LDU2 = NMAX,
     $                   LDWORK = 2*NMAX*NMAX )
*     .. Local Scalars .. 
      CHARACTER*1      JOBU
      INTEGER          I, INFO, J, N
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA, NMAX), DWORK(LDWORK), QG(LDQG, NMAX+1),
     $                 RES(LDRES,3*NMAX+1), U1(LDU1,NMAX),
     $                 U2(LDU2, NMAX), WI(NMAX), WR(NMAX)
*     .. External Functions ..
      EXTERNAL         DLANHA, DLAORS, LSAME
      LOGICAL          LSAME
      DOUBLE PRECISION DLANHA, DLAORS
*     .. External Subroutines ..
      EXTERNAL         DGEMM, DLACPY, DLASET, DSHES, DSKMV, DSKR2K,
     $                 DTRMM
*     .. Executable Statements ..
      WRITE ( NOUT, FMT = 99999 )
*     Skip the heading in the data file and read the data.
      READ ( NIN, FMT = '()' )
      READ ( NIN, FMT = * )  N, JOBU
      IF ( N.LE.0 .OR. N.GT.NMAX ) THEN
         WRITE ( NOUT, FMT = 99993 ) N
      ELSE
         READ ( NIN, FMT = * ) ( ( A(I,J), J = 1,N ), I = 1,N )
         READ ( NIN, FMT = * ) ( ( QG(I,J), J = 1,N+1 ), I = 1,N )
         IF ( LSAME( JOBU, 'U' ) ) THEN
            CALL DLACPY( 'All', N, N, A, LDA, RES(1,N+1), LDRES )
            CALL DLACPY( 'All', N, N+1, QG, LDQG, RES(1,2*N+1), LDRES )
         END IF
         CALL DSHES( JOBU, N, A, LDA, QG, LDQG, U1, LDU1, U2, LDU2,
     $               WR, WI, DWORK, LDWORK, INFO )
         INFO = 0
         IF ( INFO.NE.0 ) THEN
            WRITE ( NOUT, FMT = 99998 ) INFO
         ELSE
            IF ( LSAME( JOBU, 'U' ) ) THEN
               WRITE ( NOUT, FMT = 99997 )
               DO 10  I = 1, N
                  WRITE ( NOUT, FMT = 99994 )
     $                  ( U1(I,J), J = 1,N ), ( U2(I,J), J = 1,N )
10             CONTINUE
               DO 20  I = 1, N
                  WRITE ( NOUT, FMT = 99994 )
     $                  ( -U2(I,J), J = 1,N ), ( U1(I,J), J = 1,N )
20             CONTINUE
               WRITE ( NOUT, FMT = 99992 ) DLAORS( .FALSE., .FALSE., N,
     $                 U1, LDU1, U2, LDU2, RES, LDRES )
            END IF
            WRITE ( NOUT, FMT = 99996 )
            DO 30  I = 1, N
               WRITE ( NOUT, FMT = 99994 ) ( A(I,J), J = 1,N )
30          CONTINUE
            WRITE ( NOUT, FMT = 99995 )
            DO 40  I = 1, N
               WRITE ( NOUT, FMT = 99994 ) ( QG(I,J), J = 1,N+1 )
40          CONTINUE
C
            IF ( LSAME( JOBU, 'U' ) ) THEN
C
               CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $                     U1, LDU1, A, LDA, ZERO, RES, LDRES )
               CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, -ONE,
     $                     RES, LDRES, U1, LDU1, ONE, RES(1,N+1),
     $                     LDRES )
               CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, ONE, 
     $                     U2, LDU2, A, LDA, ZERO, RES, LDRES )
               CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, -ONE,
     $                     RES, LDRES, U2, LDU2, ONE, RES(1,N+1),
     $                     LDRES )
               DO 50 I = 1, N
                  CALL DSKMV( 'Upper', N, ONE, QG(1,2), LDQG, U2(I,1),
     $                        LDU2, ZERO, RES(1,I), 1 )
50             CONTINUE
               CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N,
     $                     -ONE, U1, LDU1, RES, LDRES, ONE, RES(1,N+1),
     $                     LDRES )
               CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $                     U2, LDU2, A, LDA, ZERO, RES, LDRES )
               CALL DSKR2K( 'Lower', 'No Transpose', N, N, ONE, RES,
     $                      LDRES, U1, LDU1, ONE, RES(1,2*N+1), LDRES,
     $                      INFO )
               CALL DLACPY( 'Full', N, N-1, U2, LDU2, RES, LDRES )
               CALL DTRMM(  'Right', 'Upper' , 'No Transpose',
     $                      'Not unit', N, N-1, ONE, QG(1,3), LDQG,
     $                       RES, LDRES )
               CALL DSKR2K( 'Lower', 'No Transpose', N, N-1, ONE, RES,
     $                      LDRES, U2(1,2), LDU2, ONE, RES(1,2*N+1),
     $                      LDRES, INFO )
               CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $                     U1, LDU1, A, LDA, ZERO, RES, LDRES )
               CALL DSKR2K( 'Upper', 'No Transpose', N, N, ONE, RES,
     $                      LDRES, U2, LDU2, ONE, RES(1,2*N+2), LDRES,
     $                      INFO )

               CALL DLACPY( 'Full', N, N-1, U1, LDU1, RES, LDRES )
               CALL DTRMM(  'Right', 'Upper' , 'No Transpose',
     $                      'Not unit', N, N-1, ONE, QG(1,3), LDQG,
     $                       RES, LDRES )
               CALL DSKR2K( 'Upper', 'No Transpose', N, N-1, -ONE, RES,
     $                      LDRES, U1(1,2), LDU1, ONE, RES(1,2*N+2),
     $                      LDRES, INFO )
C
               WRITE ( NOUT, FMT = 99991 )  DLANHA( 'Skew-Hamiltonian',
     $                'Frobenius', N, RES(1,N+1), LDRES, RES(1,2*N+1),
     $                LDRES, DWORK )
            END IF
C
         END IF
      END IF    
*
99999 FORMAT (' DSHES EXAMPLE PROGRAM RESULTS',/1X)
99998 FORMAT (' INFO on exit from DSHES = ',I2)
99997 FORMAT (' The orthogonal symplectic factor U is ')
99996 FORMAT (/' The matrix A in Schur canonical form is ')
99995 FORMAT (/' The reduced matrix QG is ')
99994 FORMAT (20(1X,F9.4))
99993 FORMAT (/' N is out of range.',/' N = ',I5)
99992 FORMAT (/' Orthogonality of U: || U''*U - I ||_F = ',G7.2)
99991 FORMAT (/' Residual: || W - U*R*U'' ||_F = ',G7.2)
      END
