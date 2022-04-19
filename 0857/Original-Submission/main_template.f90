!*****************************************************************************
! Authors: HaiJun Su, J. Michael McCarthy, Masha Sosonkina, Layne T. Watson.
!
! Date:    August, 2004. 
!*****************************************************************************
! This file contains a sample main program and user written subroutine
! for the POLSYS_GLP package. 
!
! This is a modified version of MAIN_TEMPLATE for POLSYS_PLP adapted
! for parallel computation using MPI and general linear product set
! structures.

PROGRAM MAIN_TEMPLATE
!
! MAIN_TEMPLATE is a template for calling BEZOUT_GLP and POLSYS_GLP.
! There are two options provided by MAIN_TEMPLATE:  (1) MAIN_TEMPLATE
! returns only the generalized GLP Bezout number ("root count") of the
! target polynomial system based on a system partition provided by the
! user (calls BEZOUT_GLP) or (2) MAIN_TEMPLATE returns the root count,
! homotopy path tracking statistics, error flags, and the roots (calls
! POLSYS_GLP).  For the first option set the SYSGLPSET logical switch
! ROOT_COUNT_ONLY = .TRUE., and for the second option set ROOT_COUNT_ONLY
! = .FALSE..
! 
! The input file contains data that defines the target system, system
! covering (set structure), and set degrees of the start system.
!
! This main program illustrates how to find the root count and solve a 
! polynomial system.  The data is read in using NAMELISTs, which makes the 
! input data file self-explanatory.  The problem definition is given in
! the NAMELIST /PROBLEM/ and the GLP system set structure is defined in the
! NAMELIST /SYSGLPSET/.  A new polynomial system definition is
! signalled by setting the variable NEW_PROBLEM=.TRUE. in the /PROBLEM/
! namelist.  Thus, a data file describing a polynomial system to
! solve might look like:
! 
! &PROBLEM NEW_PROBLEM=.TRUE.  data  /
! &SYSGLPSET ROOT_COUNT_ONLY=.TRUE.  data  /  count roots
! 
! &PROBLEM NEW_PROBLEM=.TRUE.  data  /
! &SYSGLPSET ROOT_COUNT_ONLY=.FALSE.  data  /  find roots
! 
! Note that static arrays are used below only to support NAMELIST input;
! the actual storage of the polynomial system and covering information
! in the data structures in the module GLOBAL_GLP is very compact.
!
! INCLUDE 'mpif.h' ! Exists in POLSYS2.

USE POLSYS2
IMPLICIT NONE

! Local variables.
INTEGER, PARAMETER:: MMAXT = 500, NN = 20
INTEGER:: BGLP, I, IFLAG1, J, K, M, MAXT, N, NUMRR = 1
INTEGER:: NUM_REAL, NUM_FAILED, NUM_SOL
INTEGER, DIMENSION(NN):: NUM_SETS, NUM_TERMS
INTEGER, DIMENSION(NN,NN):: NUM_INDICES, SET_DEG
INTEGER, DIMENSION(NN,NN,NN):: INDEX
INTEGER, DIMENSION(NN,MMAXT,NN):: DEG
INTEGER, DIMENSION(:), POINTER:: IFLAG2, INDEX_PATH_TRACKED, NFE
REAL (KIND=R8):: FINALTOL, SINGTOL, TRACKTOL
REAL (KIND=R8), DIMENSION(8):: SSPAR
REAL (KIND=R8), DIMENSION(NN):: SCALE_FACTORS
REAL (KIND=R8), DIMENSION(:), POINTER:: ARCLEN, LAMBDA
COMPLEX (KIND=R8), DIMENSION(NN,MMAXT):: COEF
COMPLEX (KIND=R8), DIMENSION(:,:), POINTER:: ROOTS
CHARACTER (LEN=80):: TITLE
CHARACTER (LEN=80), DIMENSION(NN):: DG, P
CHARACTER (LEN=80):: INPUT_FILE_NAME, OUTPUT_FILE_NAME

LOGICAL:: NEW_PROBLEM, NO_SCALING, RECALL, ROOT_COUNT_ONLY, USER_F_DF

! MPI variables.
INTEGER, PARAMETER:: MASTER_PROC = 0 !Process 0 is the master process.
INTEGER:: IERR, RC
INTEGER:: NUM_PROC  ! The number of processes.
INTEGER:: RANK_PROC ! The process RANK_PROC.

INTEGER, DIMENSION(:), POINTER:: PATH_COUNT,PATH_COUNT_DISP
REAL (KIND=R8):: END_TIME, START_TIME
REAL (KIND=R8), DIMENSION(:), POINTER:: RUN_TIME


NAMELIST /PROBLEM/ COEF,DEG,FINALTOL,N,NEW_PROBLEM,NUMRR,NUM_TERMS,&
                   TITLE,TRACKTOL,SINGTOL,SSPAR 
NAMELIST /SYSGLPSET/ DG,INDEX,NUM_INDICES,NUM_SETS,P,ROOT_COUNT_ONLY,SET_DEG

!Disassociate pointers.
NULLIFY(IFLAG2, NFE, ARCLEN, LAMBDA, ROOTS, INDEX_PATH_TRACKED) 

! Initialize MPI.                               
CALL MPI_INIT(IERR)

IF (IERR .NE. 0) THEN
 WRITE (*,*) 'Error starting MPI program. Terminating.'
 CALL MPI_ABORT(MPI_COMM_WORLD, RC, IERR)
 STOP
END IF

! Get my process number, RANK_PROC.
CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK_PROC, IERR)
! Get total number of processes used.
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUM_PROC, IERR)

! MAIN_TEMPLATE reads the target polynomial system definition and the
! system covering (set structure) specification from the file INPUT.DAT.
INPUT_FILE_NAME = 'INPUT.DAT'
OUTPUT_FILE_NAME = 'OUTPUT.DAT'

IF (RANK_PROC .EQ. MASTER_PROC) THEN
  WRITE (*,*) 'Total of ', NUM_PROC, ' processes have been initialized.'
  WRITE (*,*) 'Processors are reading POLSYS_GLP data from "', &
    TRIM(INPUT_FILE_NAME), '".'
END IF

ALLOCATE(PATH_COUNT(NUM_PROC))
ALLOCATE(PATH_COUNT_DISP(NUM_PROC))
ALLOCATE(RUN_TIME(NUM_PROC))

SSPAR(1:8) = 0.0_R8 ; DEG = 0 ; COEF = (0.0_R8,0.0_R8)

IF (RANK_PROC .EQ. MASTER_PROC) THEN
  OPEN (UNIT=7,FILE=OUTPUT_FILE_NAME,ACTION='WRITE',STATUS='REPLACE',& 
      DELIM='NONE')
END IF
OPEN (UNIT=3,FILE=INPUT_FILE_NAME,ACTION='READ',POSITION='REWIND',   &
      DELIM='APOSTROPHE',STATUS='OLD')


MAIN_LOOP: DO

! Zero out  various counters.
NUM_SOL = 0      ! Number of finite solutions.
NUM_REAL = 0     ! Number of finite real solutions.
NUM_FAILED = 0   ! Number of homotopy path tracking failures.

READ (3,NML=PROBLEM,END=1000)

! Allocate storage for the target system in POLYNOMIAL.
IF (NEW_PROBLEM) THEN
  CALL CLEANUP_POL
  ALLOCATE(POLYNOMIAL(N))
  DO I=1,N
    POLYNOMIAL(I)%NUM_TERMS = NUM_TERMS(I)
    ALLOCATE(POLYNOMIAL(I)%TERM(NUM_TERMS(I)))
    DO J=1,NUM_TERMS(I)
      ALLOCATE(POLYNOMIAL(I)%TERM(J)%DEG(N+1))
      POLYNOMIAL(I)%TERM(J)%COEF = COEF(I,J) 
      POLYNOMIAL(I)%TERM(J)%DEG(1:N) = DEG(I,J,1:N)
    END DO
  END DO
END IF

READ (3,NML=SYSGLPSET)

! Allocate storage for the system set structure in COVER. 
CALL CLEANUP_PAR
ALLOCATE(COVER_SIZES(N))
COVER_SIZES(1:N) = NUM_SETS(1:N)
ALLOCATE(COVER(N))
DO I=1,N
  ALLOCATE(COVER(I)%SET(COVER_SIZES(I)))
  DO J=1,COVER_SIZES(I)
    COVER(I)%SET(J)%NUM_INDICES = NUM_INDICES(I,J)
    COVER(I)%SET(J)%SET_DEG = SET_DEG(I,J)
    ALLOCATE(COVER(I)%SET(J)%INDEX(NUM_INDICES(I,J)))
    COVER(I)%SET(J)%INDEX(1:NUM_INDICES(I,J)) = &
                    INDEX(I,J,1:NUM_INDICES(I,J))
  END DO 
END DO

IF (ROOT_COUNT_ONLY) THEN
  ! Have the master compute GLP Bezout number.
  IF (RANK_PROC .EQ. MASTER_PROC) THEN
    MAXT = MAXVAL(NUM_TERMS(1:N))
    CALL BEZOUT_GLP(N,MAXT,SINGTOL,BGLP)
  END IF
ELSE
  ! Compute all BGLP roots of the target polynomial system.
  IF (RANK_PROC .EQ. MASTER_PROC) THEN
    WRITE (*,*) 'Path tracking started.'
  END IF
  START_TIME = MPI_WTIME()

  ! Compute roots of the target polynomial system.
  CALL POLSYS_GLP(INDEX_PATH_TRACKED, PATH_COUNT(RANK_PROC+1),N,&
     TRACKTOL,FINALTOL,SINGTOL,SSPAR,BGLP,IFLAG1,IFLAG2,&
     ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS,NUMRR=NUMRR)

  END_TIME = MPI_WTIME()
  RUN_TIME(RANK_PROC+1) = END_TIME - START_TIME
  IF (RANK_PROC .EQ. MASTER_PROC) THEN
    WRITE (*,*) 'Path tracking finished.'
  END IF 
  IF (IFLAG1 .NE. 0) THEN
    IF (RANK_PROC .EQ. MASTER_PROC) THEN
      WRITE (*,*) 'Program aborted with error flag IFLAG1 =', IFLAG1,'.'
    END IF
    GOTO 1000
  END IF

  ! Gather PATH_COUNT and timings from slave processes to the master. 
  CALL MPI_GATHER(PATH_COUNT(RANK_PROC+1),1,MPI_INTEGER,PATH_COUNT, &
       1,MPI_INTEGER,MASTER_PROC,MPI_COMM_WORLD,IERR)

  CALL MPI_GATHER(RUN_TIME(RANK_PROC+1),1,MPI_DOUBLE_PRECISION,RUN_TIME, &
       1,MPI_DOUBLE_PRECISION,MASTER_PROC,MPI_COMM_WORLD,IERR)

  IF (RANK_PROC .EQ. MASTER_PROC) THEN
    PATH_COUNT_DISP(1) = 0
    DO I=2,NUM_PROC
      PATH_COUNT_DISP(I) = PATH_COUNT_DISP(I-1) + PATH_COUNT(I-1)
    END DO
  END IF

! Gather solution information from slave processes to the master.
  CALL MPI_GATHERV(INDEX_PATH_TRACKED, PATH_COUNT(RANK_PROC+1), MPI_INTEGER,&
     INDEX_PATH_TRACKED, PATH_COUNT, PATH_COUNT_DISP, MPI_INTEGER, &
     MASTER_PROC, MPI_COMM_WORLD,IERR)

  CALL MPI_GATHERV(ARCLEN, PATH_COUNT(RANK_PROC+1), MPI_DOUBLE_PRECISION, &
     ARCLEN, PATH_COUNT, PATH_COUNT_DISP, MPI_DOUBLE_PRECISION, &
     MASTER_PROC, MPI_COMM_WORLD,IERR)

  CALL MPI_GATHERV(NFE,PATH_COUNT(RANK_PROC+1),MPI_INTEGER,NFE,PATH_COUNT,&
     PATH_COUNT_DISP,MPI_INTEGER,MASTER_PROC,MPI_COMM_WORLD,IERR)

  CALL MPI_GATHERV(IFLAG2,PATH_COUNT(RANK_PROC+1),MPI_INTEGER,IFLAG2,&
     PATH_COUNT,PATH_COUNT_DISP,MPI_INTEGER,MASTER_PROC,MPI_COMM_WORLD,IERR)

  CALL MPI_GATHERV(LAMBDA, PATH_COUNT(RANK_PROC+1),MPI_DOUBLE_PRECISION, &
     LAMBDA,PATH_COUNT,PATH_COUNT_DISP,MPI_DOUBLE_PRECISION, MASTER_PROC, &
     MPI_COMM_WORLD, IERR)

! Gather roots, each of size N padded with one homogeneous
! variable, from slave processes.
  CALL MPI_GATHERV(ROOTS(1:N+1, 1:PATH_COUNT(RANK_PROC+1)), &
     (N+1)*PATH_COUNT(RANK_PROC+1),MPI_DOUBLE_COMPLEX,&
     ROOTS,(N+1)*PATH_COUNT,(N+1)*PATH_COUNT_DISP,MPI_DOUBLE_COMPLEX, &
     MASTER_PROC, MPI_COMM_WORLD,IERR)
END IF

! Master process controls solutions' output.
IF (RANK_PROC .EQ. MASTER_PROC) THEN
  IF (NEW_PROBLEM) THEN
    WRITE (7,190) TITLE,TRACKTOL,FINALTOL,SINGTOL,SSPAR(5),N
    190 FORMAT(///A80//'TRACKTOL, FINALTOL =',2ES22.14, &
      /,'SINGTOL (0 SETS DEFAULT) =',ES22.14, &
      /,'SSPAR(5) (0 SETS DEFAULT) =',ES22.14, &
      /,'NUMBER OF EQUATIONS =',I3)
  END IF
  IF (.NOT. ROOT_COUNT_ONLY) THEN
    WRITE (*,*) 'Master process is writing solutions to "',  &
      TRIM(OUTPUT_FILE_NAME), '".'
    M = 1
    DO I=1,NUM_PROC
      WRITE (7,500) I,PATH_COUNT(I), RUN_TIME(I)
      500 FORMAT(/'===== PROCESSOR ', I4,' TRACKED ', I6, &
                  ' PATHS IN ', ES12.3, ' secs =====')
      DO K=1,PATH_COUNT(I)
        WRITE (7,600) INDEX_PATH_TRACKED(M),ARCLEN(M),NFE(M),IFLAG2(M)
        600 FORMAT(/'PATH NUMBER =',I10//'ARCLEN =',ES22.14/'NFE =',I5/ &
                   'IFLAG2 =',I3)
        ! Designate solutions as "FAILED" or "NORMAL."
        IF (MOD(IFLAG2(M),10) == 1) THEN
          ! Normal error return.
          WRITE (7,610) 1.0_R8,LAMBDA(M)
          610 FORMAT('LAMBDA =',ES22.14,', ESTIMATED ERROR =',ES22.14)
          ! Designate solutions as "FINITE" or "INFINITE."
          IF (ABS(ROOTS(N+1,M)) < 1.0E-6_R8) THEN
            WRITE (7,620) 
            620 FORMAT('SOLUTION AT INFINITY'/)
          ELSE
            NUM_SOL = NUM_SOL + 1
            ! Designate solutions as "REAL" or "COMPLEX."
            IF (ANY(ABS(AIMAG(ROOTS(1:N,M))) >= 1.0E-4_R8)) THEN
              WRITE (7,630)
              630 FORMAT('FINITE COMPLEX SOLUTION'/)
            ELSE
              NUM_REAL = NUM_REAL + 1
              WRITE (7,640)
              640 FORMAT('FINITE REAL SOLUTION'/)
            END IF
          END IF
        ELSE
          NUM_FAILED = NUM_FAILED + 1 ! Homotopy path tracking failed.
          WRITE (7,650) LAMBDA(M)
          650 FORMAT('LAMBDA =',ES22.14/)
        END IF

        WRITE (7,660) (J,ROOTS(J,M),J=1,N)
        660 FORMAT(('X(',I2,') = (',ES22.14,',',ES22.14,')'))
        WRITE (7,670) N + 1, ROOTS(N+1,M)
        670 FORMAT(/,'X(',I2,') = (',ES22.14,',',ES22.14,')')

        M = M + 1
      END DO
    END DO

    WRITE (7,910) NUM_PROC, BGLP, NUM_SOL, NUM_REAL, NUM_SOL-NUM_REAL, &
          BGLP-NUM_SOL-NUM_FAILED, NUM_FAILED, MAXVAL(RUN_TIME)
    910 FORMAT(  &
          /'=========Number of processors used: ',I6,'  ========', &
          /'Bezout GLP number (BGLP)          : ',I6 &
          /'Number of finite solutions        : ',I6 &
          /'Number of finite real solutions   : ',I6 &
          /'Number of finite complex solutions: ',I6 &
          /'Number of solutions at infinity   : ',I6 &
          /'Number of homotopy path failures  : ',I6 &
          /'Maximum running time              : ',ES11.3,' secs',/52('='))
   ELSE  ! ROOT_COUNT_ONLY=.TRUE.
     WRITE (7,1010) BGLP,(J,TRIM(P(J)),J,TRIM(DG(J)),J=1,N)
     1010 FORMAT(/60('='),/ &
            'Bezout GLP number (BGLP) =',I8,' for the system covering:', &
            /('P(',I2,') = ',A,',   DG(',I2,') = ',A))
     WRITE (7,FMT="(60('='))")
   END IF
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

END DO MAIN_LOOP

DEALLOCATE(PATH_COUNT)
DEALLOCATE(PATH_COUNT_DISP)
DEALLOCATE(RUN_TIME)

1000 IF (RANK_PROC .EQ. MASTER_PROC) CLOSE (UNIT=7)  
CLOSE (UNIT=3)  
CALL CLEANUP_POL
CALL CLEANUP_PAR

! Shut down MPI.
CALL MPI_FINALIZE(IERR)

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

! Deallocates structure COVER.

IF (.NOT. ALLOCATED(COVER)) RETURN
DO I=1,SIZE(COVER)
  DO J=1,COVER_SIZES(I)
    DEALLOCATE(COVER(I)%SET(J)%INDEX)
  END DO
  DEALLOCATE(COVER(I)%SET)
END DO
DEALLOCATE(COVER)
DEALLOCATE(COVER_SIZES)  
RETURN
END SUBROUTINE CLEANUP_PAR

END PROGRAM MAIN_TEMPLATE

SUBROUTINE TARGET_SYSTEM_USER(N,PROJ_COEF,XC,F,DF)
! Template for user written subroutine to evaluate the (complex) target
! system F(XC) and its (complex) N x N Jacobian matrix DF(XC).  XC(1:N+1)
! is in complex projective coordinates, and the homogeneous coordinate
! XC(N+1) is explicitly eliminated from F(XC) and DF(XC) using the
! projective transformation (cf. the comments in START_POINTS_GLP).  The
! comments in the internal subroutine TARGET_SYSTEM should be read before
! attempting to write this subroutine; pay particular attention to the
! handling of the homogeneous coordinate XC(N+1).  DF(:,N+1) is not
! referenced by the calling program.

USE REAL_PRECISION
USE GLOBAL_GLP
IMPLICIT NONE
INTEGER, INTENT(IN):: N
COMPLEX (KIND=R8), INTENT(IN), DIMENSION(N+1):: PROJ_COEF,XC
COMPLEX (KIND=R8), INTENT(OUT):: F(N), DF(N,N+1)

! For greater efficiency, replace the following code (which is just the
! internal POLSYS_GLP subroutine TARGET_SYSTEM) with hand-crafted code.

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
  F(I) = TS
END DO 

DF = (0.0_R8,0.0_R8)

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
