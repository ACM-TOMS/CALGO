SUBROUTINE SSTEBNDGAP( N, W, GAPMIN, ILO, IUO, VLO, VUO )
!
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: N
INTEGER, INTENT( IN ), OPTIONAL :: ILO, IUO
REAL( KIND=PREC ) :: GAPMIN
REAL( KIND=PREC ), INTENT( IN ), OPTIONAL :: VLO, VUO
!
!.. Array Argument ..
REAL( KIND=PREC ) :: W( * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEBNDGAP computes the minimum gap between                                 !
!     |W(ILO-1)-W(ILO)| and |W(IUO-W(IUO+1)|                                   !
!  or between                                                                  !
!     |W(right_of_VLO)-W(left_of_VLO)| and |W(right_of_VUO)-W(left_of_VUO)|.   !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          The dimension of the matrix.                                        !
!                                                                              !
!  W       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          Eigenvalues in ascending order.                                     !
!                                                                              !
!  ILO     (input,optional) INTEGER                                            !
!          Index of the smallest computed eigenvalue.                          !
!                                                                              !
!  IUO     (input,optional) INTEGER                                            !
!          Index of the largest computed eigenvalue.                           !
!                                                                              !
!  VLO     (input,optional) REAL( KIND=PREC )                                  !
!          Lower bound of the interval searched for eigenvalues.               !
!                                                                              !
!  VUO     (input,optional) REAL( KIND=PREC )                                  !
!          Upper bound of the interval searched for eigenvalues.               !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, J
REAL( KIND=PREC ) :: GAPL, GAPR, TNORM, ULP
!
!.. External Function ..
REAL( KIND=PREC ), EXTERNAL :: SLAMCH
!
!.. Intrinsic Functions ..
INTRINSIC ABS, MAX, MIN
!
!.. Executable Statements ......................................................
!
ULP = SLAMCH( 'Precision' )
TNORM = MAX( ABS( W( 1 ) ), ABS( W( N ) ), ULP )
!
IF      ( PRESENT( ILO ) .AND. PRESENT( IUO ) ) THEN
!
!       Gaps for case RANGE='I'.
!
        IF ( ILO == 1 ) THEN
           GAPL = ZERO
        ELSE
           GAPL = W( ILO ) - W( ILO-1 )
        END IF
        IF ( IUO == N ) THEN
           GAPR = ZERO
        ELSE
           GAPR = W( IUO+1 ) - W( IUO )
        END IF
!
ELSE IF ( PRESENT( VLO ) .AND. PRESENT( VUO ) ) THEN
!
!       Gaps for case RANGE='V'.
!
        J = 0
        DO I = 1, N
           IF ( W( I ) > VLO ) THEN
              J = I
              EXIT
           END IF
        END DO
        IF ( J > 1 ) THEN
           GAPL = W( J ) - W( J-1 )
        ELSE
           GAPL = ZERO
        END IF
        J = N
        DO I = 1, N
           IF ( W( N-I+1 ) < VUO ) THEN
              J = N-I+1
              EXIT
           END IF
        END DO
        IF ( J < N ) THEN
           GAPR = W( J+1 ) - W( J )
        ELSE
           GAPR = ZERO
        END IF
!
END IF
!
GAPMIN = MIN( MAX( ULP, ABS( GAPL ) ), MAX( ULP, ABS( GAPR ) ) ) / TNORM
!
END SUBROUTINE SSTEBNDGAP
