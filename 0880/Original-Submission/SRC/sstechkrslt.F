SUBROUTINE SSTECHKRSLT( N, D, E, M, W, Z, LDZ, HBANDR, WORK, RESULT, INFO )
!
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: HBANDR, INFO, LDZ, M, N
!
!.. Array Arguments ..
REAL( KIND=PREC ) :: D( * ), E( * ), RESULT( 2 ), Z( LDZ, * ), W( * ), WORK( * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTECHKRSLT checks the results of the LAPACK symmetric tridiagonal          !
!  eigensolvers. Two tests are performed:                                      !
!                                                                              !
!  RESULT( 1 ) = | T - Z*W*Z' | / ( |T|*N*ULP ) if M = N                       !
!              = | Z'*T*Z - W | / ( |T|*N*ULP ) otherwise                      !
!                                                                              !
!  RESULT( 2 ) = | I - Z*Z' | / ( N*ULP ) if M = N                             !
!              = | I - Z'*Z | / ( N*ULP ) otherwise                            !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          The dimension of the matrix.                                        !
!                                                                              !
!  D       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The N diagonal elements of the tridiagonal matrix.                  !
!                                                                              !
!  E       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The (N-1) off-diagonal elements of the tridiagonal matrix in        !
!          elements 1 to N-1, E(N) is not used.                                !
!                                                                              !
!  M       (input) INTEGER                                                     !
!          The total number of eigenvalues found, 0 <= M <= N.                 !
!                                                                              !
!  W       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The first M elements contain the selected eigenvalues in            !
!          ascending order.                                                    !
!                                                                              !
!  Z       (input) REAL( KIND=PREC ) array, dimension ( LDZ, M )               !
!          The first M columns of Z contain the orthonormal eigenvectors of    !
!          the tridiagonal matrix corresponding to the selected eigenvalues,   !
!          with the i-th column of Z holding the eigenvector associated        !
!          with W(i).                                                          !
!                                                                              !
!  LDZ     (input) INTEGER                                                     !
!          The leading dimension of the array Z.                               !
!                                                                              !
!  HBANDR  (input) INTEGER                                                     !
!          Sets the halfbandwidth of the matrices used in the tests            !
!          RESULT( 1 ) = | Z'*T*Z - W | / ( |T|*N*ULP )                        !
!          RESULT( 2 ) = | I - Z'*Z | / ( N*ULP )                              !
!          i.e. max(1,N*(HBANDR/100)) subdiagonals of Z'*T*Z and Z'*Z          !
!          are computed. This option should be used mainly for very            !
!          large problems in order to save time.                               !
!                                                                              !
!  WORK    (workspace) REAL( KIND=PREC ) array, dimension ( (N+1)*N )          !
!          Workspace.                                                          !
!                                                                              !
!  RESULT  (output) REAL( KIND=PREC ) array, dimension ( 2 )                   !
!          The values computed by the two tests described above. The           !
!          values are currently limited to 1/ulp, to avoid overflow.           !
!                                                                              !
!  INFO    (input) INTEGER                                                     !
!          Exit status from routine used to compute W and Z.                   !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, J, K, L, MINMN, NDIAG, NEIGV
REAL( KIND=PREC ) :: AUKJ, RNORM, SUM, TNORM, ULP, UNFL
!
!.. External Subroutines ..
REAL( KIND=PREC ), EXTERNAL :: SLAMCH, SLANSB, SLANST, SLANSY
EXTERNAL SGEMV, SLASET, SSYR, SSYRK
!
!.. Intrinsic Functions ..
INTRINSIC INT, MAX, MIN, REAL
!
!.. Executable Statements ......................................................
!
! Early return.           
!
IF ( HBANDR == 0 .OR. INFO /= 0 .OR. M == 0 ) THEN
   RESULT = ZERO
   RETURN
END IF
!
NEIGV = MIN( M, N )
NDIAG = MAX( 1, MIN( M, INT( N*(HBANDR/HNDRD) ) ) )
MINMN = N  !** MIN( M, N ) can be pessimistic **!
!
! Get machine parameters. 
!
ULP = SLAMCH( 'Precision' )
UNFL = SLAMCH( 'Safe minimum' )
RESULT = ONE / ULP
!
! Compute the 1-norm of T.
!
TNORM = SLANST( '1', N, D, E )
TNORM = MAX( TNORM, UNFL )
!
! Test the residuals.
!
IF      ( NDIAG /= N ) THEN
!
!       | Z'*T*Z - W | / ( |T|*N*ULP )         
!
        DO I = 1, NEIGV
           K = MIN( NDIAG, NEIGV-I+1 ) 
           L = N + NDIAG*(I-1) + 1
           WORK( 1 ) = D( 1 )*Z( 1,I ) + E( 1 )*Z( 2,I ) 
           DO J = 2, N-1
              WORK( J ) = E( J-1 )*Z( J-1,I ) + D( J )*Z( J,I ) + &
                          E( J )*Z( J+1,I ) 
           END DO
           WORK( N ) = E( N-1 )*Z( N-1,I ) + D( N )*Z( N,I )  
           CALL SGEMV( 'T', N, K, ONE, Z( 1,I ), LDZ, WORK, 1, &
                       ZERO, WORK( L ), 1 )
           WORK( L ) = WORK( L ) - W( I )
        END DO
!
        RNORM = SLANSB( 'M', 'L', NEIGV, NDIAG-1, WORK( N+1 ), NDIAG, WORK )
!
ELSE IF ( M == N ) THEN
!
!       | T - Z*W*Z' | / ( |T|*N*ULP )
!
        CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
!
        DO I = 1, N-1
           J = ( N+1 )*( I-1 )
           WORK( J+1 ) = D( I )
           WORK( J+2 ) = E( I )
        END DO
        WORK( N**2 ) = D( N )
        DO I = 1, N
           CALL SSYR( 'L', N, -W( I ), Z( 1,I ), 1, WORK, N )
        END DO
!
        RNORM = SLANSY( 'M', 'L', N, WORK, N, WORK( N**2+1 ) )
!
ELSE
!
!       | Z'*T*Z - W | / ( |T|*N*ULP )
!
        CALL SLASET( 'Full', M, M, ZERO, ZERO, WORK, M )
!
        DO I = 1, M
           DO J = 1, I
              L = M*( J-1 ) + I
              DO K = 1, N
                 AUKJ = D( K )*Z( K, J )
                 IF ( K /= 1 ) AUKJ = AUKJ + E( K-1 )*Z( K-1, J )
                 IF ( K /= N ) AUKJ = AUKJ + E( K )*Z( K+1, J )
                 WORK( L ) = WORK( L ) + Z( K, I )*AUKJ
              END DO
           END DO
           L = M*( I-1 ) + I
           WORK( L ) = WORK( L ) - W( I )
        END DO
!
        RNORM = SLANSY( 'M', 'L', M, WORK, M, WORK( M**2+1 ) )
!
END IF
!
IF ( TNORM > RNORM ) THEN
   RESULT( 1 ) = ( RNORM / TNORM ) / ( MINMN*ULP )
ELSE
   IF ( TNORM < ONE ) THEN
      RESULT( 1 ) = ( MIN( RNORM, MINMN*TNORM ) / TNORM )/ ( MINMN*ULP )
   ELSE
      RESULT( 1 ) = MIN( RNORM / TNORM, REAL( MINMN,PREC ) ) / ( MINMN*ULP )
   END IF
END IF
!
! Test the orthogonality of the computed eigenvectors.
!
IF      ( NDIAG /= N ) THEN
!
!       | I - Z'*Z | / ( N*ULP )        
!
        DO I = 1, NEIGV
           K = MIN( NDIAG, N-I+1 ) 
           L = N + NDIAG*(I-1) + 1
           CALL SGEMV( 'T', N, K, ONE, Z( 1,I ), LDZ, Z( 1,I ), 1, &
                       ZERO, WORK( L ), 1 )
           WORK( L ) = ONE - WORK( L )
        END DO
!
        RNORM = SLANSB( 'M', 'L', NEIGV, NDIAG-1, WORK( N+1 ), NDIAG, WORK )
!
ELSE IF ( M == N ) THEN
!
!       | I - Z*Z' | / ( N*ULP ) 
!
        CALL SSYRK( 'L', 'N', N, N, ONE, Z, LDZ, ZERO, WORK, N )
!
        DO I = 1, N
           J = ( N+1 )*( I-1 ) + 1
           WORK( J ) = WORK( J ) - ONE
        END DO
!
        RNORM = SLANSY( 'M', 'L', N, WORK, N, WORK( N**2+1 ) )
!
ELSE
!
!       | I - Z'*Z | / ( N*ULP )
!
        CALL SSYRK( 'L', 'T', M, N, ONE, Z, LDZ, ZERO, WORK, M )
!
        DO I = 1, M
           J = ( M+1 )*( I-1 ) + 1
           WORK( J ) = WORK( J ) - ONE
        END DO
!
        RNORM = SLANSY( 'M', 'L', M, WORK, M, WORK( M**2+1 ) )
!
END IF
!
RESULT( 2 ) = MIN( REAL( MINMN,PREC ), RNORM ) / ( MINMN*ULP )
!
END SUBROUTINE SSTECHKRSLT
