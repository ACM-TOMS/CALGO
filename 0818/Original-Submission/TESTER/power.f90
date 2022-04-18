      PROGRAM MAIN_POWER
!     Purpose: Build sparse matrix and 
!              perform power iteration via calls to Sparse BLAS

!     Header file of prototypes and named constants
      USE BLAS_SPARSE
      INTEGER, PARAMETER :: NMAX = 4, NNZ = 6
      INTEGER A,I,N,NITERS,ISTAT
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDX,JNDX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: VAL,Q,WORK
      REAL(KIND=DP) LAMBDA

      ALLOCATE (VAL(NNZ),INDX(NNZ),JNDX(NNZ))
      ALLOCATE (Q(NMAX),WORK(NMAX))

!     -----------------------------------
!     Define matrix, in coordinate format
!     -----------------------------------
      VAL  = (/ 1.1_dp, 2.2_dp, 2.4_dp, 3.3_dp, 4.1_dp, 4.4_dp/)
      INDX = (/  1,  2,  2,  3,  4,  4/)
      JNDX = (/  1,  2,  4,  3,  1,  4/)

      N = NMAX

!     ----------------------------------
!     Step 1:  Create Sparse BLAS handle
!     ----------------------------------
      CALL DUSCR_BEGIN( N, N, A, ISTAT)

!     -----------------------------------
!     Step 2:  Insert entries all at once
!     -----------------------------------
      CALL USCR_INSERT_ENTRIES(A, VAL, INDX, JNDX, ISTAT)

!     -----------------------------------------------
!     Step 3:  Complete construction of sparse matrix
!     -----------------------------------------------
      CALL USCR_END(A, ISTAT)


!     -----------------------------------------------
!     Step 4:  Call Power Method Routine
!     -----------------------------------------------
!     q      - eigenvector approximation. 
!     lambda - eigenvalue approximation.
      NITERS = 100
      CALL POWER_METHOD(A, Q, LAMBDA, N, NITERS, WORK, ISTAT)
      IF (ISTAT.NE.0) THEN
         WRITE(*,*) 'ERROR IN POWER_METHOD = ',ISTAT
      ELSE
         WRITE(*,*) 'NUMBER OF ITERATIONS = ',NITERS
         WRITE(*,*) 'APPROXIMATE DOMINANT EIGENVALUE = ',LAMBDA
      ENDIF

!     -----------------------------------------------
!     Step 5:  Release Matrix Handle
!     -----------------------------------------------
      CALL USDS(A,ISTAT)


    CONTAINS

      
        SUBROUTINE POWER_METHOD(A, Q, LAMBDA, N, NITERS, Z, ISTAT)
        USE BLAS_SPARSE
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(:),INTENT(INOUT) :: Q(N),Z(N)
        REAL(KIND=DP), INTENT(OUT) :: LAMBDA
        INTEGER, INTENT(IN) :: A,N,NITERS
        INTEGER, INTENT(OUT):: ISTAT
        INTEGER I,ITER,ISEED
        REAL(KIND=DP):: NORMZ
        REAL Y
        INTRINSIC RANDOM_NUMBER,DOT_PRODUCT
        
        !     Fill Z by random numbers
        DO I = 1, N
           CALL RANDOM_NUMBER(HARVEST=Y)
           Z(I)=DBLE(Y)
        END DO
        
        DO ITER = 1, NITERS
           !Compute 2-norm of Z
           NORMZ = SQRT(DOT_PRODUCT(Z(1:N),Z(1:N)))
           !Normalize Z
           if (NORMZ.NE.0) Z(1:N) = Z(1:N)/NORMZ
           !Copy Z to Q
           Q=Z
           !Set Z to 0
           Z=0.D0
           !Compute new Z
           CALL USMV(A, Q, Z, ISTAT)
           !Test error flag
           IF (ISTAT.NE.0) RETURN
           !New LAMBDA
           LAMBDA = DOT_PRODUCT(Q,Z)
        END DO

        RETURN
        END SUBROUTINE POWER_METHOD

      END PROGRAM MAIN_POWER
