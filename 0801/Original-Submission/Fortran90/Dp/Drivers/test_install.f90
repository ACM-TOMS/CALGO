! This file contains a main program to test the correctness of the
! compiled code; it is uncommented and has no further use beyond testing
! the installation.  Author: Layne T. Watson, 10/1999.

! Compile this file (free form Fortran 90) and link it to the object
! files from the compiles of polsys_plp.f90 (free form) and lapack_plp.f
! (fixed format).  Then run the executable with input file INPUT.DAT
! (upper case).  A message indicating apparent success or failure of the
! installation is written to standard out.

PROGRAM TEST_INSTALL

USE POLSYS

IMPLICIT NONE
INTEGER, PARAMETER:: NN=30, MMAXT=50
INTEGER:: BPLP, I, IFLAG1, J, K, M, MAXT, N, NUMRR=1
INTEGER, DIMENSION(NN):: NUM_TERMS, NUM_SETS
INTEGER, DIMENSION(NN,NN):: NUM_INDICES
INTEGER, DIMENSION(NN,NN,NN):: INDEX
INTEGER, DIMENSION(NN,MMAXT,NN):: DEG
INTEGER, DIMENSION(:), POINTER:: IFLAG2, NFE
REAL (KIND=R8):: TRACKTOL, FINALTOL, SINGTOL
REAL (KIND=R8), DIMENSION(8):: SSPAR
REAL (KIND=R8), DIMENSION(NN):: SCALE_FACTORS
REAL (KIND=R8), DIMENSION(:), POINTER:: ARCLEN, LAMBDA
COMPLEX (KIND=R8), DIMENSION(NN,MMAXT):: COEF
COMPLEX (KIND=R8), DIMENSION(:,:), POINTER:: ROOTS
COMPLEX (KIND=R8), DIMENSION(2,4):: EROOTS = RESHAPE(SOURCE=(/  &
  (  2.34233851959121E+03_R8,  0.0E00_R8),   &
  ( -7.88344824094120E-01_R8,  0.0E00_R8),   &
  (  9.08921229615388E-02_R8,  0.0E00_R8),   &
  ( -9.11497098197499E-02_R8,  0.0E00_R8),   &
  (  1.61478579234357E-02_R8,  1.68496955498881E+00_R8),     &
  (  2.67994739614461E-04_R8,  4.42802993973661E-03_R8),     &
  (  1.61478579234359E-02_R8, -1.68496955498881E+00_R8),     &
  (  2.67994739614461E-04_R8, -4.42802993973661E-03_R8)  /), &
  SHAPE=(/ 2,4 /) )
CHARACTER (LEN=80):: TITLE
CHARACTER (LEN=80), DIMENSION(NN):: P
LOGICAL:: NEW_PROBLEM, ROOT_COUNT_ONLY

NAMELIST /PROBLEM/ COEF, DEG, FINALTOL, NEW_PROBLEM, N, NUMRR, NUM_TERMS,  &
  SINGTOL, SSPAR, TITLE, TRACKTOL
NAMELIST /SYSPARTITION/ INDEX, NUM_INDICES, NUM_SETS, P, ROOT_COUNT_ONLY

NULLIFY(IFLAG2, NFE, ARCLEN, LAMBDA, ROOTS) ! Disassociate pointers.

OPEN (UNIT=3,FILE='INPUT.DAT',ACTION='READ',POSITION='REWIND',   &
     DELIM='APOSTROPHE',STATUS='OLD')

SSPAR(1:8) = 0.0_R8 ; DEG = 0 ; COEF = (0.0_R8,0.0_R8)

READ (3,NML=PROBLEM)

CALL CLEANUP_POL
ALLOCATE(POLYNOMIAL(N))
DO I=1,N
  POLYNOMIAL(I)%NUM_TERMS=NUM_TERMS(I)
  ALLOCATE(POLYNOMIAL(I)%TERM(NUM_TERMS(I)))
  DO J=1,NUM_TERMS(I)
    ALLOCATE(POLYNOMIAL(I)%TERM(J)%DEG(N+1))
    POLYNOMIAL(I)%TERM(J)%COEF=COEF(I,J) 
    POLYNOMIAL(I)%TERM(J)%DEG(1:N)=DEG(I,J,1:N)
  END DO
END DO

READ (3,NML=SYSPARTITION)

CALL CLEANUP_PAR
ALLOCATE(PARTITION_SIZES(N))
PARTITION_SIZES(1:N)=NUM_SETS(1:N)
ALLOCATE(PARTITION(N))
DO I=1,N
  ALLOCATE(PARTITION(I)%SET(PARTITION_SIZES(I)))
  DO J=1,PARTITION_SIZES(I)
    PARTITION(I)%SET(J)%NUM_INDICES=NUM_INDICES(I,J)
    ALLOCATE(PARTITION(I)%SET(J)%INDEX(NUM_INDICES(I,J)))
    PARTITION(I)%SET(J)%INDEX(1:NUM_INDICES(I,J)) = &
      INDEX(I,J,1:NUM_INDICES(I,J))
  END DO 
END DO

CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,IFLAG1,IFLAG2, &
    ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS)

SINGTOL = 0.0_R8
DO I=1,BPLP
  SINGTOL = MAX(SINGTOL, MINVAL(SUM(ABS(SPREAD(   &
    EROOTS(1:2,I),DIM=2,NCOPIES=BPLP) - ROOTS(1:2,1:BPLP)), DIM=1)))
END DO

IF (SINGTOL < 1.0E-6_R8) THEN
  WRITE (*,*) 'Test problem was solved correctly. The installation ', &
    'appears correct.'
ELSE
  WRITE (*,*) 'Warning!  Test problem was not solved correctly.'
END IF

CLOSE (UNIT=3)
CALL CLEANUP_POL
CALL CLEANUP_PAR
STOP

CONTAINS

SUBROUTINE CLEANUP_POL

! Deallocates structure POLYNOMIAL.

IF (.NOT. ALLOCATED(POLYNOMIAL)) RETURN
DO I=1,SIZE(POLYNOMIAL)
  DO J=1,NUMT(I)
    DEALLOCATE(POLYNOMIAL(I)%TERM(J)%DEG)
  END DO
  DEALLOCATE(POLYNOMIAL(I)%TERM)
END DO
DEALLOCATE(POLYNOMIAL)
RETURN
END SUBROUTINE CLEANUP_POL

SUBROUTINE CLEANUP_PAR

! Deallocates structure PARTITION.

IF (.NOT. ALLOCATED(PARTITION)) RETURN
DO I=1,SIZE(PARTITION)
  DO J=1,PARTITION_SIZES(I)
    DEALLOCATE(PARTITION(I)%SET(J)%INDEX)
  END DO
  DEALLOCATE(PARTITION(I)%SET)
END DO
DEALLOCATE(PARTITION)
DEALLOCATE(PARTITION_SIZES)	
RETURN
END SUBROUTINE CLEANUP_PAR

END PROGRAM TEST_INSTALL

                                                                      !!!
SUBROUTINE TARGET_SYSTEM_USER(N,PROJ_COEF,XC,F,DF)
! Template for user written subroutine to evaluate the (complex) target
! system F(XC) and its (complex) N x N Jacobian matrix DF(XC).  XC(1:N+1)
! is in complex projective coordinates, and the homogeneous coordinate
! XC(N+1) is explicitly eliminated from F(XC) and DF(XC) using the
! projective transformation (cf. the comments in START_POINTS_PLP).  The
! comments in the internal subroutine TARGET_SYSTEM should be read before
! attempting to write this subroutine; pay particular attention to the
! handling of the homogeneous coordinate XC(N+1).  DF(:,N+1) is not
! referenced by the calling program.

USE REAL_PRECISION
USE GLOBAL_PLP
IMPLICIT NONE
INTEGER, INTENT(IN):: N
COMPLEX (KIND=R8), INTENT(IN), DIMENSION(N+1):: PROJ_COEF,XC
COMPLEX (KIND=R8), INTENT(OUT):: F(N), DF(N,N+1)

! For greater efficiency, replace the following code (which is just the
! internal POLSYS_PLP subroutine TARGET_SYSTEM) with hand-crafted code.

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
INTEGER:: DEGREE, I, J, K, L
COMPLEX (KIND=R8):: T, TS
DO I=1,N
  TS = (0.0_R8, 0.0_R8)
  DO J=1,POLYNOMIAL(I)%NUM_TERMS
    T = POLYNOMIAL(I)%TERM(J)%COEF
    DO K=1,N+1
      DEGREE = POLYNOMIAL(I)%TERM(J)%DEG(K)
      IF (DEGREE == 0) CYCLE
      T = T * XC(K)**DEGREE
    END DO
    TS = TS + T
  END DO
  F(I)=TS
END DO 

DF=(0.0_R8,0.0_R8)

DO I=1,N
  DO J=1,N+1
    TS = (0.0_R8,0.0_R8)
    DO K=1,POLYNOMIAL(I)%NUM_TERMS
      DEGREE = POLYNOMIAL(I)%TERM(K)%DEG(J)
      IF (DEGREE == 0) CYCLE
      T = POLYNOMIAL(I)%TERM(K)%COEF * DEGREE * (XC(J)**(DEGREE - 1))
      DO L=1,N+1
        DEGREE = POLYNOMIAL(I)%TERM(K)%DEG(L)
        IF ((L == J) .OR. (DEGREE == 0)) CYCLE
        T = T * (XC(L)**DEGREE)
      END DO
      TS = TS + T
    END DO
    DF(I,J) = TS
  END DO
END DO

DO I=1,N
  DF(I,1:N) = DF(I,1:N) + PROJ_COEF(1:N) * DF(I,N+1)
END DO
! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

RETURN
END SUBROUTINE TARGET_SYSTEM_USER
