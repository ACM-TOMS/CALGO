SUBROUTINE SSTEDEFMTRX( MTYPE, N, D, E )
!
USE SSTEDEFINITIONS
!
!.. Scalar Arguments ..
INTEGER :: MTYPE, N
!
!.. Array Arguments ..
REAL( KIND=PREC ) :: D( N ), E( N )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEDEFMTRX generates a symmetric tridiagonal matrix as indicated by MTYPE. !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  MTYPE  (input) INTEGER                                                      !
!         Sets the matrix type:                                                !
!         = 0, zero matrix                                                     !
!         = 1, identity matrix                                                 !
!         = 2, (1,2,1) matrix                                                  !
!         = 3, Wilkinson matrix                                                !
!         = 4, Clement matrix                                                  !
!         = 5, Legendre orthogonal polynomials                                 !
!         = 6, Laguerre orthogonal polynomials                                 !
!         = 7, Hermit orthogonal polynomials                                   !
!                                                                              !
!  D      (output) REAL( KIND=PREC ) array, dimension ( N )                    !
!         The N diagonal elements of the tridiagonal matrix.                   !
!                                                                              !
!  E      (output) REAL( KIND=PREC ) array, dimension ( N )                    !
!         The (N-1) off-diagonal elements of the tridiagonal                   !
!         matrix in elements 1 to N-1, E(N) is set to 0.                       !
!                                                                              !
!  N      (input) INTEGER                                                      !
!         The order of the matrix. N > 0.                                      !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, J, K
REAL( KIND=PREC ) :: S
!
!.. Intrinsic Functions ..
INTRINSIC  CEILING, REAL
!
!.. Executable Statements ......................................................
!
SELECT CASE ( MTYPE )
!
CASE ( 0 )
!
!    Zero matrix
!
     D = ZERO
     E = ZERO
!
CASE ( 1 )
!
!    Identity matrix
!
     D = ONE
     E = ZERO
!
CASE ( 2 )
!
!    (1,2,1)-type matrix, for example 
!
!                   1   1   1
!    T = tridiag( 2   2   2   2 )
!                   1   1   1
!
     E( N ) = ZERO
     D( 1:N ) = TWO
     IF ( N > 1 ) E( 1:N-1 ) = ONE
!
CASE ( 3 )
!
!    Wilkinson-type matrix, for example
!
!                   1   1   1   1   1   1
!    T = tridiag( 3   2   1   0   1   2   3 )
!                   1   1   1   1   1   1
!
     IF ( N == 1 ) THEN
        D( 1 ) = ZERO
        E( 1 ) = ZERO
     ELSE
        J = N / 2
        K = CEILING( REAL( N,PREC ) / TWO ) + 1
        S = REAL( N+1,PREC )/TWO
        DO I = 1, J
           S = S - ONE
           D( I ) = S 
        END DO
        D( J + 1 ) = ZERO
        DO I = K, N 
           D( I ) = S
           S = S + ONE
        END DO
        DO I = 1, N - 1
           E( I ) = ONE
        END DO
        E( N ) = ZERO
     END IF
!
CASE ( 4 )
!
!    Clement-type matrix, for example
!
!                   2.2361   2.8284   3.0000   2.8284   2.2361
!    T = tridiag( 0        0        0        0        0        0 )
!                   2.2361   2.8284   3.0000   2.8284   2.2361
!
     D( 1:N ) = ZERO
     DO I = 1, N-1
        E( I ) = SQRT( FLOAT( I*(N-I) ) )
     END DO
!
CASE ( 5 )
!
!    Legendre orthogonal polynomials
!
     D = ZERO
     E( N ) = ZERO
     DO I = 2, N
        E( I-1 ) = FLOAT( I ) / SQRT( FLOAT( (2*I-1)*(2*I+1) ) )
     END DO
!
CASE ( 6 )
!
!    Laguerre orthogonal polynomials
!
     DO I = 1, N-1
        D( I ) = FLOAT( 2*I+1 )
        E( I ) = FLOAT( I+1 )
     END DO
     D( N ) = FLOAT( 2*N+1 )
     E( N ) = ZERO
!
CASE ( 7 )
!
!    Hermit orthogonal polynomials
!
     D = ZERO
     E( 1 ) = ONE
     E( N ) = ZERO
     DO I = 2, N-1
        E( I ) = SQRT( FLOAT( I ) )
     END DO
!
CASE DEFAULT
!
     D = ZERO
     E = ZERO
!
END SELECT
! 
END SUBROUTINE SSTEDEFMTRX
