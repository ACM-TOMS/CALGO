PROGRAM TESTF
!
!
! This program tests the Linear Shepard scattered data
! interpolation package.  The robust option and Ripple are
! tested also.
!

  USE REAL_PRECISION
  USE SHEPPACK, ONLY: LSHEP, LSHEPVAL, RIPPLE
  IMPLICIT NONE
  INTEGER :: I, IER, J, K, M, N, N1, NT, POWER, TEMP
  REAL(KIND=R8) :: E_bar, E_r, FA, FET, RND, X_MIN, X_MAX
  REAL(KIND=R8), ALLOCATABLE :: F(:), A(:,:), RW(:), X(:,:), XT(:)
  DATA X_MIN, X_MAX /0.0E+00_R8, 1.0E+00_R8/
  
  N1 = 8
  M = 5
!
! M is the dimension of the test data
! N1 is the number of samples per dimension for the test points
! 
  NT = N1**M
  DO K = 0, 4
     N = 100*(2**K)
     WRITE(*, *) '========================================================='
     WRITE(*, *) 'dimension                     :', M
     WRITE(*, *) 'number of sample points       :', N
     ALLOCATE ( A(M,N), X(M,N), F(N), RW(N), XT(M) )
!
! Generate N m-dimensional points for the initial sample.
!
     CALL RANDOM_SEED( )
     DO I = 1, N
        DO J = 1, M
           CALL RANDOM_NUMBER( RND )
           X(J, I) = X_MIN + (X_MAX - X_MIN) * RND
        END DO
     END DO
!
! Generate the function values for the initial sample.
!
        DO I = 1, N
           F(I) = F_EXACT( X(1:M, I), M ) 
        END DO

!
! Perform Linear Shepard interpolation on a uniform grid
! of N1**M points.  Calculate and report the error.
!
        CALL LSHEP ( M, N, X, F, A, RW, IER )
        E_r = 0.0_R8
        E_bar = 0.0_R8
        DO I = 0, NT-1
           TEMP = I
           DO J = 1, M
              POWER = N1**(M - J)
              XT(J) = 0.1_R8 + (TEMP / POWER) * (0.8_R8 / (N1 - 1))
              TEMP = TEMP - (TEMP / POWER) * POWER
           END DO
           FA =  LSHEPVAL( XT, M, N, X, F, A, RW, IER )
           FET = F_EXACT( XT, M )
           E_r = E_r + (FA - FET)**2
           E_bar = E_bar + ABS(FA - FET)
        END DO
        WRITE(*, *) '---------------LSHEP-------------------------------------'
        WRITE(*, *) 'The root mean squared error   :', SQRT( E_r / NT )
        WRITE(*, *) 'The mean absolute error       :', E_bar / NT
     
!
! Perform Robust Linear Shepard interpolation on a uniform grid
! of N1**M points.  Calculate and report the error.
!
        CALL LSHEP ( M, N, X, F, A, RW, IER, .TRUE. )
        E_r = 0.0_R8
        E_bar = 0.0_R8
        DO I = 0, NT-1
           TEMP = I
           DO J = 1, M
              POWER = N1**(M - J)
              XT(J) = 0.1_R8 + (TEMP / POWER) * (0.8_R8 / (N1 - 1))
              TEMP = TEMP - (TEMP / POWER) * POWER
           END DO
           FA =  LSHEPVAL( XT, M, N, X, F, A, RW, IER )
           FET = F_EXACT( XT, M )
           E_r = E_r + (FA - FET)**2
           E_bar = E_bar + ABS(FA - FET)
        END DO
        WRITE(*, *) '---------------Robust LSHEP------------------------------'
        WRITE(*, *) 'The root mean squared error   :', SQRT( E_r / NT )
        WRITE(*, *) 'The mean absolute error       :', E_bar / NT    

!
! Perform Ripple interpolation on a uniform grid
! of N1**M points.  Calculate and report the error.
!
        CALL RIPPLE(M, N, X, F, A, RW, IER)
        E_r = 0.0_R8
        E_bar = 0.0_R8
        DO I = 0, NT-1
           TEMP = I
           DO J = 1, M
              POWER = N1**(M - J)
              XT(J) = 0.1_R8 + (TEMP / POWER) * (0.8_R8 / (N1 - 1))
              TEMP = TEMP - (TEMP / POWER) * POWER
           END DO
           FA =  LSHEPVAL( XT, M, N, X, F, A, RW, IER )
           FET = F_EXACT( XT, M )
           E_r = E_r + (FA - FET)**2
           E_bar = E_bar + ABS(FA - FET)
        END DO
        WRITE(*, *) '---------------Ripple------------------------------------'
        WRITE(*, *) 'The root mean squared error   :', SQRT( E_r / NT )
        WRITE(*, *) 'The mean absolute error       :', E_bar / NT   
        WRITE(*, *) '---------------------------------------------------------'
     DEALLOCATE ( A, X, F, RW, XT )
  END DO
  
CONTAINS
  FUNCTION F_EXACT(X, M)
    USE REAL_PRECISION
    IMPLICIT NONE
    INTEGER :: M
    REAL(KIND=R8), DIMENSION(M) :: X        
    REAL(KIND=R8) :: F_EXACT
    REAL(KIND=R8), DIMENSION(M) :: TEMP
    
       TEMP(1:M) = 2.0_R8 * MIN( X(1:M), 1.0_R8 - X(1:M) )
       F_EXACT = PRODUCT( TEMP(1:M) )
    RETURN 
  END FUNCTION F_EXACT
     
END PROGRAM TESTF
