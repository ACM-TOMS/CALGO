SUBROUTINE DSTEDEFEIGV( ECOND, EDIST, ESIGN, ETYPE, ISEED, N, W )
!
USE DSTEDEFINITIONS
!
!.. Scalar Arguments ..
INTEGER :: ECOND, EDIST, ESIGN, ETYPE, N
!
!.. Array Arguments ..
INTEGER :: ISEED( 4 )
REAL( KIND=PREC ) :: W( N )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEDEFEIGV defines eigenvalue distributions to be used in the generation   !
!  of test matrices. See also subroutine DLATMS in LAPACK/TESTING/MATGEN.      !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  ECOND  (input) INTEGER                                                      !
!         Sets the condition number for ETYPE = 1:5,                           !
!         = 1, COND = 1 / SQRT(ULP), default                                   !
!         = 2, COND = 1 / (N*SQRT(ULP)),                                       !
!         = 3, COND = 1 / (10*N*SQRT(ULP)),                                    !
!         = 4, COND = 1 / ULP,                                                 !
!         = 5, COND = 1 / (N*ULP),                                             !
!         = 6, COND = 1 / (10*N*ULP),                                          !
!                                                                              !
!  EDIST  (input) INTEGER                                                      !
!         Specifies the type of the distribution to be used in the generation  !
!         of random eigenvalues                                                !
!         = 1, UNIFORM(  0, 1 )                                                !
!         = 2, UNIFORM( -1, 1 )                                                !
!         = 3, NORMAL( 0, 1 )                                                  !
!                                                                              !
!  ESIGN  (input) INTEGER                                                      !
!         Attributes signs to the entries in W                                 !
!         = 0, the entries of W are unchanged                                  !
!         = 1, multiplies each entry of W by +/- 1 with probability 0.5        !
!                                                                              !
!  ISEED  (input/output) INTEGER array, dimension ( 4 )                        !
!         Seed for the random number generator. Each entry of ISEED should     !
!         lie between 0 and 4095 inclusive and ISEED(4) should be odd.         !
!                                                                              !
!  ETYPE  (input) INTEGER                                                      !
!         Defines the eigenvalue distribution,                                 !
!         =  1, W( I ) = 1, W( 2:N ) = 1.0/COND;                               !
!         =  2, W( 1:N-1 ) = 1 and W( N ) = 1.0/COND;                          !
!         =  3, W( I ) = COND**(-(I-1)/(N-1));                                 !
!         =  4, W( I ) = 1-(I-1)/(N-1)*(1-1/COND);                             !
!         =  5, W is set to random numbers in the range (1/COND,1),            !
!               their logarithms are uniformly distributed;                    !
!         =  6, W is set to random numbers from the same distribution          !
!               as of the rest of the matrix;                                  !
!         =  7, W( I ) = ULP*I, I=1,2,...N-1, W( N ) = 1;                      !
!         =  8, W( 1 ) = ULP, W( I ) = 1+SQRT(ULP)*I, I=2,3,...N-1,            !
!               W( N ) = 2;                                                    !
!         =  9, W( 1 ) = 1, W( I ) = W( I-1 )+100*ULP, I=2,3,...N              !
!         <  0, has the same meaning as ABS(ETYPE), except that the            !
!               order of the elements of W is reversed.                        !
!                                                                              !
!         ULP = ( relative machine precision ) * ( base of the machine )       !
!                                                                              !
!  N      (input) INTEGER                                                      !
!         The order of the matrix, N > 0.                                      !
!                                                                              !
!  W      (output) REAL( KIND=PREC ) array, dimension ( N )                    !
!         The eigenvalue distribution.                                         !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I
REAL( KIND=PREC ) :: ALPHA, COND, TEMP, ULP
!
!.. External Functions ..
REAL( KIND=PREC ), EXTERNAL :: DLAMCH, DLARAN
!
!.. External Subroutine ..
EXTERNAL DLARNV
!
!.. Intrinsic Functions ..
INTRINSIC ABS, REAL, SQRT 
!
!.. Executable Statements ......................................................
!
ULP = DLAMCH( 'Precision' )
!
! Set the condition number for ETYPE = 1:5.        
!
SELECT CASE ( ABS( ETYPE ) )
   CASE ( 1:6 )
      SELECT CASE ( ECOND )
         CASE ( 1 ); COND = ONE / SQRT( ULP )
         CASE ( 2 ); COND = ONE / ( SQRT( ULP )*REAL( N,PREC ) )
         CASE ( 3 ); COND = ONE / ( SQRT( ULP )*REAL( N*10,PREC ) )
         CASE ( 4 ); COND = ONE / ULP
         CASE ( 5 ); COND = ONE / ( ULP*REAL( N,PREC ) )
         CASE ( 6 ); COND = ONE / ( ULP*REAL( N*10,PREC ) )
         CASE DEFAULT; COND = ONE 
      END SELECT
   CASE DEFAULT; COND = ONE 
END SELECT
!
! Generate the eigenvalue distribution.
!
SELECT CASE ( ABS( ETYPE ) )
!
CASE ( 1 )
!
!    W( 1 ) = 1, W( I ) = 1/COND, I=2,3,...N
!
     TEMP = ONE / COND
     W( 1 ) = ONE
     DO I = 2, N
        W( I ) = TEMP
     END DO
!
CASE ( 2 )
!
!    W( I ) = 1, I=1,2,...N-1, W( N ) = 1/COND
!
     DO I = 1, N-1
        W( I ) = ONE 
     END DO
     W( N ) = ONE / COND
!
CASE ( 3 )
!
!    W( I ) = COND**(-(I-1)/(N-1))
!
     W( 1 ) = ONE
     ALPHA = COND**( -ONE / REAL( N-1,PREC ) )
     DO I = 2, N
        W( I ) = ALPHA**( I-1 )
     END DO

CASE ( 4 )
!
!    W( I ) = 1-(I-1)/(N-1)*(1-1/COND)
!
     W( 1 ) = ONE
     TEMP = ONE / COND
     ALPHA = ( ONE-TEMP ) / REAL( N-1,PREC )
     DO I = 2, N
        W( I ) = REAL( N-I,PREC )*ALPHA + TEMP
     END DO
!
CASE ( 5 ) 
!
!    W = randomly distributed values on ( 1/COND , 1)
!
     ALPHA = LOG( ONE / COND )
     DO I = 1, N
        W( I ) = EXP( ALPHA*DLARAN( ISEED ) )
     END DO
!
CASE ( 6 ) 
!
!    W = randomly distributed values from EDIST
!
     CALL DLARNV( EDIST, ISEED, N, W )
!
CASE ( 7 ) 
!
!    W( I ) = ULP*I, I=1,2,...N-1, W( N ) = 1
!
     W( 1 ) = ULP
     DO I = 2, N-1
        W( I ) = W( I-1 ) + ULP
     END DO
     W( N ) = ONE
!
CASE ( 8 ) 
!
!    W( 1 ) = ULP, W( I ) = 1+SQRT(ULP)*I, I=2,3,...N-1, W( N ) = 2
!
     TEMP = SQRT( ULP )

     W( 2 ) = ONE + TEMP
     DO I = 3, N-1
        W( I ) = W( I-1 ) + TEMP
     END DO
     W( N ) = TWO
!
CASE ( 9 )
!
!    W( 1 ) = 1, W( I ) = W( I-1 )+100*ULP, I=2,3,...N
!
     TEMP = HNDRD*ULP
     W( 1 ) = ONE
     DO I = 2, N
        W( I ) = W( I-1 ) + TEMP
     END DO
!
CASE DEFAULT
!
     W( 1:N ) = ZERO
!
END SELECT
!
! Assign random signs to W.
!
IF ( ESIGN == 1 ) THEN 
   DO I = 1, N
      TEMP = DLARAN( ISEED )
      IF ( TEMP > HALF ) W( I ) = -W( I )
   END DO
END IF
!
! Reverse if ETYPE < 0. 
! 
IF ( ETYPE < 0 ) THEN
   DO I = 1, N / 2
      TEMP = W( I )
      W( I ) = W( N+1-I )
      W( N+1-I ) = TEMP
   END DO
END IF
! 
END SUBROUTINE DSTEDEFEIGV
