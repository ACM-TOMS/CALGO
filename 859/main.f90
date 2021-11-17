PROGRAM MAIN

   USE BABDCR
!
!  This is a driver program to solve babd linear systems Ax = b whose 
!  coefficient matrix A has square blocks of the same size: 
!
!         ( Ba                                            Bb     )
!         ( S(0)   R(1)                                          )
!         (        S(1)   R(2)                                   )
!         (                                                      )
!     A = (               .       .                              )
!         (                      .      .                        )
!         (                            .      .                  )
!         (                                                      )
!         (                               S(NBLOKS-1)  R(NBLOKS) )
!
!  The matrix A must be given in three arrays: 
!
!   MATR_A  = [S(0) R(1) S(1) R(2),...., S(NBLOKS-1) R(NBLOKS)],
!
!   LFTBLK = [Ba],  RGTBLK = [Bb],
!
!  where each block has size NRWBLK by NRWBLK. It could 
!  be stored in a file containing MATR_A followed by LFTBLK and
!  RGTBLK. The r.h.s. is given in the array
!
!     VECT_B  = [f(0), f(1), ...., f(NBLOKS)],
!
!  where each element has length NRWBLK.
!
!  In order to use the subroutines BABDCR_FACT and BABDCR_SOLV, it is 
!  necessary to insert the USE BABDCR instruction. The other
!  considered subroutines are BABDCR_FACTSOLV to solve linear systems
!  without saving the matrices of the factorization, and BABDCR_SOLVT
!  that solves a system with the transpose of the factorized matrix.
! 
!  The factorization of the coefficient matrix with BABDCR_FACT
!  requires a fill-in array of size NRWBLK*NRWBLK*(NBLOKS-1) and an
!  integer array of size NRWBLK*NBLOKS to save the permutation
!  matrices. A call to BABDCR_FACT or BABDCR_FACTSOLV modifies the
!  original coefficient matrix. A call to BABDCR_SOLV or BABDCR_SOLVT
!  or BABDCR_FACTSOLV gives the solution of the considered linear
!  system in place of the r.h.s.
!
!  This program has the following features:
!    1) factorization and solution separately:
!       the factorization is saved and it is possible to compute the
!       condition number of the coefficient matrix and to perform
!       an iterative refinement
!    2) solution of a system previously factorized
!    3) factorization and solution with the same subroutine:
!       it is not possible to reuse the factorization and the 
!       coefficient matrix itself
!    4) solution of the transposed system with the matrix previously
!       factorized
!

   REAL(DP), ALLOCATABLE :: MATR_A( :, :, : ), VECT_B(:,:), &
                            LFTBLK( :, : ), RGTBLK( :, : ), &
                            FILL_IN( :, :, : ), &
                            VECT_R( :, : ), VECT_T( :, : ), &
                            MATR_O( :, :, : )
   INTEGER, ALLOCATABLE :: PERM( :, : )

   INTEGER :: NRWBLK, NBLOKS, FACTORED, INFO
   REAL(DP) :: NRM1, NRM1INV, BABDCR_NORM1, BABDCR_NORM1INV
   CHARACTER :: ANS, ANSCN, ANSIR

!   REAL(DP) :: TIME
!   INTEGER :: START_TIME, END_TIME, RATE_TIME, MAX_TIME


    FACTORED= 0
    ANS= INSERT_VALUE( FACTORED )
    DO WHILE ( ANS.NE.'5' )

      CALL BABDCR_INPUT( ANS, NRWBLK, NBLOKS, ANSCN, ANSIR )

      IF ( (ANSIR.EQ.'y').OR.(ANSIR.EQ.'Y') ) THEN
        ALLOCATE( VECT_T( NRWBLK, NBLOKS+1 ) )
        VECT_T = VECT_B
      END IF

      SELECT CASE ( ANS )
      CASE ('1')

        IF ( (ANSIR.EQ.'y').OR.(ANSIR.EQ.'Y') ) THEN
          IF ( ALLOCATED(MATR_O) ) DEALLOCATE( MATR_O )
          ALLOCATE( MATR_O( NRWBLK, NRWBLK, (NBLOKS+1)*2 ) )
          MATR_O(:,:,1:NBLOKS*2) = MATR_A
          MATR_O(:,:,NBLOKS*2+1) = LFTBLK
          MATR_O(:,:,NBLOKS*2+2) = RGTBLK 
        END IF

        IF ( (ANSCN.EQ.'y').OR.(ANSCN.EQ.'Y') ) THEN
          NRM1= BABDCR_NORM1( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK )
        END IF

! Start timing
!        CALL SYSTEM_CLOCK( COUNT=START_TIME, COUNT_RATE=RATE_TIME, &
!                           COUNT_MAX=MAX_TIME)

! factorization and solution
        CALL BABDCR_FACT( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                          PERM, FILL_IN, INFO )
        IF ( INFO.GT.0 ) THEN
          PRINT *, 'The coefficient matrix is singular'
          FACTORED= 0
        ELSE
          CALL BABDCR_SOLV( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                            PERM, FILL_IN, VECT_B )
          FACTORED= 1
        END IF

! Stop timing and compute the elapsed time in seconds 
!        CALL SYSTEM_CLOCK( END_TIME )
!        TIME= REAL(END_TIME-START_TIME)/REAL(RATE_TIME)
!        PRINT *, 'BABDCR_FACT + BABDCR_SOLV ', TIME
        
        IF ( (ANSCN.EQ.'y').OR.(ANSCN.EQ.'Y') ) THEN
          NRM1INV= BABDCR_NORM1INV( NRWBLK, NBLOKS, MATR_A, LFTBLK, &
                                    RGTBLK, PERM, FILL_IN )
          PRINT *, 'The condition number of the babd matrix is ', &
                   NRM1*NRM1INV
        END IF

      CASE ('2')
!
! Start timing
!        CALL SYSTEM_CLOCK( COUNT=START_TIME, COUNT_RATE=RATE_TIME, &
!                           COUNT_MAX=MAX_TIME)

! solution
        CALL BABDCR_SOLV( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                          PERM, FILL_IN, VECT_B )

! Stop timing and compute the elapsed time in seconds 
!        CALL SYSTEM_CLOCK( END_TIME )
!        TIME= REAL(END_TIME-START_TIME)/REAL(RATE_TIME)
!        PRINT *, 'BABDCR_FACT + BABDCR_SOLV ', TIME

      CASE ('3')
!
! Start timing
!        CALL SYSTEM_CLOCK( COUNT=START_TIME, COUNT_RATE=RATE_TIME, &
!                           COUNT_MAX=MAX_TIME)

! factorization and solution
        CALL BABDCR_FACTSOLV( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                             VECT_B, INFO )
        FACTORED= 0
        IF ( INFO.GT.0 ) THEN
          PRINT *, 'The coefficient matrix is singular'
          PRINT *, 'the solution can not be computed'
        ENDIF

! Stop timing and compute the elapsed time in seconds 
!        CALL SYSTEM_CLOCK( END_TIME )
!        TIME= REAL(END_TIME-START_TIME)/REAL(RATE_TIME)
!        PRINT *, 'BABDCR_FACT + BABDCR_SOLV ', TIME

      CASE ('4')
!
! Start timing
!        CALL SYSTEM_CLOCK( COUNT=START_TIME, COUNT_RATE=RATE_TIME, &
!                           COUNT_MAX=MAX_TIME)

! solution
        CALL BABDCR_SOLVT( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                           PERM, FILL_IN, VECT_B )

! Stop timing and compute the elapsed time in seconds 
!        CALL SYSTEM_CLOCK( END_TIME )
!        TIME= REAL(END_TIME-START_TIME)/REAL(RATE_TIME)
!        PRINT *, 'BABDCR_FACT + BABDCR_SOLV ', TIME

      END SELECT

      IF ( (ANSIR.EQ.'y').OR.(ANSIR.EQ.'Y') ) THEN
        ALLOCATE( VECT_R( NRWBLK, NBLOKS+1 ) )
        CALL BABD_MATMULT( 'N', NRWBLK, NBLOKS, &
                           MATR_O(:,:,1:NBLOKS*2), &
                           MATR_O(:,:,NBLOKS*2+1), &
                           MATR_O(:,:,NBLOKS*2+2), &
                           VECT_B, VECT_R )
        VECT_R= VECT_T-VECT_R
        PRINT *, 'norm of the residual ', MAXVAL(ABS(VECT_R))
        CALL BABDCR_SOLV( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                          PERM, FILL_IN, VECT_R )
        VECT_B= VECT_B+VECT_R
        CALL BABD_MATMULT( 'N', NRWBLK, NBLOKS, &
                           MATR_O(:,:,1:NBLOKS*2), &
                           MATR_O(:,:,NBLOKS*2+1), &
                           MATR_O(:,:,NBLOKS*2+2), &
                           VECT_B, VECT_R )
        VECT_R= VECT_T-VECT_R
        PRINT *, 'norm of the residual (after the recursion step) ', &
                 MAXVAL(ABS(VECT_R))
        DEALLOCATE( VECT_T, VECT_R )
      ENDIF

      IF ( INFO.EQ.0 ) THEN
        PRINT *, 'Solution vector (in VECT_B)'
        PRINT *, VECT_B
      ENDIF

      ANS= INSERT_VALUE( FACTORED )

    ENDDO

   IF ( ALLOCATED(VECT_B) ) DEALLOCATE( VECT_B )
   IF ( ALLOCATED(MATR_A) ) DEALLOCATE( MATR_A )
   IF ( ALLOCATED(LFTBLK) ) DEALLOCATE( LFTBLK )
   IF ( ALLOCATED(RGTBLK) ) DEALLOCATE( RGTBLK )
   IF ( ALLOCATED(PERM) ) DEALLOCATE( FILL_IN, PERM )
   IF ( ALLOCATED(MATR_O) ) DEALLOCATE( MATR_O )

   STOP

   CONTAINS

     FUNCTION INSERT_VALUE( FACTORED )
!
! Purpose
! =======
!
!   This function specifies the operation to perform on the linear
!   system among factorization of the coefficient matrix, solution,
!   or both.
!
! Parameters
! ==========
!
! Input variables:
!   FACTORED integer,           flag variable;
!                               if FACTORED=0, the factorization is not
!                               available
!                               if FACTORED=1, the factorization has been
!                               previously computed and is available
!
     INTEGER, INTENT(IN) :: FACTORED
     CHARACTER :: INSERT_VALUE, ANS

      ANS= '0'
      DO WHILE ( (ANS.NE.'1').AND.(ANS.NE.'2').AND.(ANS.NE.'3') &
                 .AND.(ANS.NE.'4').AND.(ANS.NE.'5') )
        WRITE(*,9999)
        READ  *, ANS
        print *, ANS
        IF ( (FACTORED.EQ.0).AND.((ANS.EQ.'2').OR.(ANS.EQ.'4')) ) THEN
          WRITE(*,9998)
          ANS= '0'
        END IF
      ENDDO

9999  FORMAT ('----Menu---------------- '/, &
             ' 1.      Factorization (BABDCR_FACT) and Solution ', &
             '(BABDCR_SOLV)'/, &
             ' 2.      Solution (BABDCR_SOLV)'/, &
             ' 3.      Factorization and Solution (BABDCR_FACTSOLV)'/, &
             ' 4.      Solution of the transposed system ', &
             '(BABDCR_SOLVT)'/, &
             ' 5.      Exit'  )

9998  FORMAT ( 'A factorized coefficient matrix does not exist'/, &
              'It is necessary to select a coefficient matrix '/, &
              ' that must be previously factorized by choosing '/, &
              ' 1 (Factorization) in the Menu' )

      INSERT_VALUE= ANS

     RETURN
     END FUNCTION INSERT_VALUE



     SUBROUTINE BABDCR_INPUT( ANS, NRWBLK, NBLOKS, ANSCN, ANSIR )
!
!  == subroutine BABDCR_INPUT =========
!  Giuseppe Romanazzi, Pierluigi Amodio
!  Dipartimento di Matematica
!  Universita' di Bari
!  October 20, 2005
!  ====================================
!
! Purpose
! ======
!
!  This subroutine defines the input data for the linear system.
!  It uses the global variables previously defined. If these are 
!  already defined, they are deallocated and then again allocated
!  with the new dimensions.
!  
!  The coefficient matrix may be defined by means of an intrinsic 
!  random fortran procedure or, read from an input file containing 
!  MATR_A followed by LFTBLK and RGTBLK.  The r.h.s. may be read
!  from an input file, defined as a random vector or such that the
!  solution has all the components equal to 1.
!  
!
! Parameters
! ==========
!
! Input variables:
!   ANS character,    the kind of operation that may be performed.
!                     If ANS='1' or ANS='3' MATR_A is read from file
!                     or using an intrinsic random fortran procedure;
!                     if ANS='2' or ANS='3' VECT_B is read from file
!                     or using an intrinsic random fortran procedure
!
! Input/Output variables:
!   NRWBLK  integer,  the dimension of each block of the babd matrix, 
!                     NRWBLK>0
!
!   NBLOKS  integer,  the number of blocks S(i) and R(i) in the babd
!                     matrix, NBLOKS>0
!
! Output variables:
!   ANSCN character,  if ANSCN='1', it is required to compute the
!                     condition number of the babd matrix
!
!   ANSIR character,  if ANSIR='1', it is required to compute one step
!                     of iterative refinement
!
     CHARACTER, INTENT(IN)   :: ANS
     INTEGER, INTENT(IN OUT) :: NRWBLK, NBLOKS
     CHARACTER, INTENT(OUT)  :: ANSCN, ANSIR
! local variables
     REAL(DP), ALLOCATABLE :: SOL_X(:,:)
     INTEGER               :: I, J, K
     CHARACTER             :: ANS1
     CHARACTER(LEN=20)     :: FILE_NAME

      IF ( (ANS.EQ.'1').OR.(ANS.EQ.'3') ) THEN
 10     CONTINUE
!
!  Define the size of each block
!

        NRWBLK=0
        DO WHILE ( NRWBLK.LE.0 )
          PRINT *, 'Insert the dimension of each block, NRWBLK>0 '
          READ  *, NRWBLK
        ENDDO
!
!  Define the number of blocks S(i) and R(i)
!
        NBLOKS=0
        DO WHILE ( NBLOKS.LE.0 )
          PRINT *, 'Insert the number of block rows, NBLOKS>0 '
          READ  *, NBLOKS
        ENDDO
!
!  Definition the new problem
!    
        IF ( ALLOCATED(MATR_A) ) DEALLOCATE( MATR_A )
        IF ( ALLOCATED(LFTBLK) ) DEALLOCATE( LFTBLK )
        IF ( ALLOCATED(RGTBLK) ) DEALLOCATE( RGTBLK )
        IF ( ALLOCATED(VECT_B) ) DEALLOCATE( VECT_B )
        ALLOCATE( MATR_A( NRWBLK, NRWBLK, NBLOKS*2 ), &
                  LFTBLK( NRWBLK, NRWBLK ), &
                  RGTBLK( NRWBLK, NRWBLK ), &
                  VECT_B( NRWBLK, NBLOKS+1 ) )
!
!  Input of the coefficient matrices
!
        ANS1= '0'
        DO WHILE ( (ANS1.NE.'r').AND.(ANS1.NE.'f').AND. &
                   (ANS1.NE.'R').AND.(ANS1.NE.'F') )
          PRINT *, 'Input matrix : random or from file? (r/f)'
          READ  *, ANS1
        ENDDO

        SELECT CASE ( ANS1 )
        CASE('r', 'R')
!
! random coefficient matrix
!
          CALL RANDOM_NUMBER( MATR_A )
          MATR_A= (MATR_A-0.5)*20
          CALL RANDOM_NUMBER( LFTBLK )
          LFTBLK= (LFTBLK-0.5)*20
          CALL RANDOM_NUMBER( RGTBLK )
          RGTBLK= (RGTBLK-0.5)*20
        CASE('f', 'F')
!
! the coefficient matrix is stored in a file containing MATR_A
! followed by LFTBLK and RGTBLK
!
 20       CONTINUE
          FILE_NAME= READ_FILE( 1 )
          OPEN( 10, FILE=FILE_NAME, STATUS='OLD' )
          READ( 10, *, END=5, ERR=15 ) (((MATR_A(I,J,K), I=1,NRWBLK), &
                                          J=1,NRWBLK), K=1,NBLOKS*2)
          READ( 10, *, END=5, ERR=15 ) ((LFTBLK(I,J), I=1,NRWBLK), &
                                         J=1,NRWBLK)
          READ( 10, *, END=5, ERR=15 ) ((RGTBLK(I,J), I=1,NRWBLK), &
                                         J=1,NRWBLK)
          CLOSE( 10 )
          GO TO 30
 5        PRINT *, 'The specified file has not the required length'
          PRINT *, 'Please, insert again size and number of blocks'
          CLOSE( 10 )
          GO TO 10
 15       PRINT *, 'The considered file is corrupted'
          CLOSE( 10 )
          GO TO 20
 30       CONTINUE
        END SELECT
!
!  Input of the r.h.s.
!
        ANS1= '0'
        DO WHILE ( (ANS1.NE.'r').AND.(ANS1.NE.'R').AND. &
                   (ANS1.NE.'k').AND.(ANS1.NE.'K').AND. &
                   (ANS1.NE.'f').AND.(ANS1.NE.'F') )
          PRINT *, 'Input rhs of size ', (NBLOKS+1)*NRWBLK, &
                   ' : random, with known solution (1,...,1) or ', &
                   'from file? ',  '(r/k/f)'
          READ *, ANS1
        ENDDO

        SELECT CASE ( ANS1 )
        CASE('r', 'R')
          CALL RANDOM_NUMBER( VECT_B )
          VECT_B= (VECT_B-0.5)*20
        CASE('k', 'K')
          ALLOCATE( SOL_X( NRWBLK, NBLOKS+1 ) )
          SOL_X = 1D0
          CALL BABD_MATMULT( 'N', NRWBLK, NBLOKS, MATR_A, LFTBLK, &
                             RGTBLK, SOL_X, VECT_B )
          DEALLOCATE( SOL_X )
        CASE('f', 'F')
 45       CONTINUE
          FILE_NAME= READ_FILE( 2 )
          OPEN ( 20, FILE=FILE_NAME, STATUS='OLD' )
          READ ( 20, *, END=65, ERR=55 ) (( VECT_B(I,J), &
                                           I=1,NRWBLK ), J=1,NBLOKS+1 )
          GO TO 40
 55       CONTINUE
          PRINT *, 'The considered file is corrupted'
          VECT_B = 0
          CLOSE( 20 )
          GO TO 45
 65       WRITE(*,899) 
 40       CONTINUE
          CLOSE( 20 )
        END SELECT

      ELSEIF ( (ANS.EQ.'2').OR.(ANS.EQ.'4') ) THEN

        IF ( ALLOCATED(VECT_B) ) DEALLOCATE(VECT_B)
        ALLOCATE( VECT_B( NRWBLK, NBLOKS+1 ) )

        ANS1= '0'
        DO WHILE ( (ANS1.NE.'r').AND.(ANS1.NE.'f').AND. &
                   (ANS1.NE.'R').AND.(ANS1.NE.'F') )
          PRINT *, 'Input rhs of size ', (NBLOKS+1)*NRWBLK, &
                   ' : random or from file? (r/f)'
          READ *, ANS1
        ENDDO

        SELECT CASE ( ANS1 )
        CASE('r', 'R')
          CALL RANDOM_NUMBER( VECT_B )
          VECT_B= (VECT_B-0.5)*20
        CASE('f', 'F')
 75       CONTINUE
          FILE_NAME= READ_FILE( 2 )
          OPEN( 20, FILE=FILE_NAME, STATUS='OLD' )
          READ( 20, *, END=95, ERR=85 ) (( VECT_B(I,J), &
                                          I=1,NRWBLK ), J=1,NBLOKS+1 )
          GO TO 50
 85       CONTINUE
          PRINT *, 'The considered file is corrupted'
          VECT_B = 0
          CLOSE( 20 )
          GO TO 75
 95       WRITE(*,899) 
 50       CONTINUE
          CLOSE( 20 )
        END SELECT

      ENDIF

      IF ( ANS.EQ.'1' ) THEN
        IF ( ALLOCATED(PERM) ) DEALLOCATE( FILL_IN, PERM )
        ALLOCATE( PERM( NRWBLK*2,NBLOKS ), &
                  FILL_IN( NRWBLK,NRWBLK,NBLOKS-1 ) )
!
!  Computation of the condition number of the matrix
!
        ANSCN= '0'
        DO WHILE ( (ANSCN.NE.'y').AND.(ANSCN.NE.'Y').AND. &
                   (ANSCN.NE.'n').AND.(ANSCN.NE.'N') )
          PRINT *, 'Do you want to estimate the condition number? (y/n)'
          READ *, ANSCN
        ENDDO
!
!  Computation of one step of iterative refinement
!
        ANSIR= '0'
        DO WHILE ( (ANSIR.NE.'y').AND.(ANSIR.NE.'Y').AND. &
                   (ANSIR.NE.'n').AND.(ANSIR.NE.'N') )
          PRINT *, 'Do you want to compute one step of iterative ', &
                 'refinement? (y/n)'
          READ *, ANSIR
        ENDDO
      ELSE
        ANSCN= 'n'
        IF ( ANS.NE.'2' ) THEN
          ANSIR= 'n'
        ENDIF
      ENDIF

 899  FORMAT( 'The specified file has not the required length'/, &
              'A number of zeros has been added to obtain a '/, &
              'vector of the desired length' )

     RETURN
     
     END SUBROUTINE BABDCR_INPUT
   
 
     FUNCTION READ_FILE( VAL )
!
! Purpose
! =======
!
!  This function specify the name of the function containing the 
!  input matrix or r.h.s.
!
! Parameters
! ==========
!
! Input variables:
!   VAL integer,            flag variable;
!                           if VAL=1 we need to specify the name
!                           of the input matrix;
!                           if VAL=2, we need to specify the name
!                           of the input r.h.s.
!
     INTEGER, INTENT(IN) :: VAL
     CHARACTER(LEN=20) :: READ_FILE, FILE_NAME
     LOGICAL :: EXT

      IF ( VAL.EQ.1 ) THEN
        WRITE(*,999) 'matrix'
      ELSE
        WRITE(*,999) 'r.h.s.'
      ENDIF
      READ  *, FILE_NAME
      INQUIRE ( FILE=FILE_NAME, EXIST=EXT )
      DO WHILE ( .NOT.EXT )
        PRINT *, 'A file with this name does not exist'
        IF ( VAL.EQ.1 ) THEN
          WRITE(*,999) 'matrix'
        ELSE
          WRITE(*,999) 'r.h.s.'
        ENDIF
        READ  *, FILE_NAME
        INQUIRE ( FILE=FILE_NAME, EXIST=EXT )
      END DO

 999  FORMAT ( 'Insert the name of the input file ', &
               'containing the ', A6 )

      READ_FILE = FILE_NAME

     RETURN
     
     END FUNCTION READ_FILE

   END PROGRAM MAIN
