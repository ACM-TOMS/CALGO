! This file contains a sample main program and user written subroutine
! for the POLSYS_PLP package.  Layne T. Watson, Steven M. Wise, Andrew
! J. Sommese, August, 1998.  Cosmetic changes, 10/1999.

PROGRAM MAIN_TEMPLATE
!
! MAIN_TEMPLATE is a template for calling BEZOUT_PLP and POLSYS_PLP.
! There are two options provided by MAIN_TEMPLATE:  (1) MAIN_TEMPLATE
! returns only the generalized PLP Bezout number ("root count") of the
! target polynomial system based on a system partition provided by the
! user (calls BEZOUT_PLP) or (2) MAIN_TEMPLATE returns the root count,
! homotopy path tracking statistics, error flags, and the roots (calls
! POLSYS_PLP).  For the first option set the logical switch
! ROOT_COUNT_ONLY = .TRUE., and for the second option set ROOT_COUNT_ONLY
! = .FALSE..
! 
! The file INPUT.DAT contains data for several sample target systems
! and system partitions.  This main program illustrates how to find the
! root count for several different partitions for the same polynomial
! system, and also how to solve more than one polynomial system in the
! same run.  The data is read in using NAMELISTs, which makes the data
! file INPUT.DAT self-explanatory.  The problem definition is given in
! the NAMELIST /PROBLEM/ and the PLP system partition is defined in the
! NAMELIST /SYSPARTITION/.  A new polynomial system definition is
! signalled by setting the variable NEW_PROBLEM=.TRUE. in the /PROBLEM/
! namelist.  Thus a data file describing several different polynomial
! systems to solve, and exploring different system partitions for the
! same polynomial system, might look like
! 
! &PROBLEM NEW_PROBLEM=.TRUE.  data  /
! &SYSPARTITION ROOT_COUNT_ONLY=.FALSE.  data  /  finds roots
! 
! &PROBLEM NEW_PROBLEM=.TRUE.  data  /
! &SYSPARTITION ROOT_COUNT_ONLY=.TRUE.  data  /  finds root count only
! &PROBLEM NEW_PROBLEM=.FALSE.  /
! &SYSPARTITION ROOT_COUNT_ONLY=.TRUE.  data  /  a different root count
! &PROBLEM NEW_PROBLEM=.FALSE.  /
! &SYSPARTITION ROOT_COUNT_ONLY=.TRUE.  data  /  another root count
! 
! Note that static arrays are used below only to support NAMELIST input;
! the actual storage of the polynomial system and partition information
! in the data structures in the module GLOBAL_PLP is very compact.


USE POLSYS

! Local variables.
IMPLICIT NONE
INTEGER, PARAMETER:: NN = 30, MMAXT = 50
INTEGER:: BPLP, I, IFLAG1, J, K, M, MAXT, N, NUMRR = 1
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
CHARACTER (LEN=80):: TITLE
CHARACTER (LEN=80), DIMENSION(NN):: P
LOGICAL:: NEW_PROBLEM, NO_SCALING, RECALL, ROOT_COUNT_ONLY, USER_F_DF

NAMELIST /PROBLEM/ COEF, DEG, FINALTOL, NEW_PROBLEM, N, NUMRR, NUM_TERMS,  &
  SINGTOL, SSPAR, TITLE, TRACKTOL
NAMELIST /SYSPARTITION/ INDEX, NUM_INDICES, NUM_SETS, P, ROOT_COUNT_ONLY

NULLIFY(IFLAG2, NFE, ARCLEN, LAMBDA, ROOTS) ! Disassociate pointers.

! MAIN_TEMPLATE reads the target polynomial system definition and the
! system partition specification from the file INPUT.DAT.

OPEN (UNIT=3,FILE='INPUT.DAT',ACTION='READ',POSITION='REWIND',   &
     DELIM='APOSTROPHE',STATUS='OLD')

! All output is to the file OUTPUT.DAT, which is overwritten.

OPEN (UNIT=7,FILE='OUTPUT.DAT',ACTION='WRITE',STATUS='REPLACE',DELIM='NONE')

SSPAR(1:8) = 0.0_R8 ; DEG = 0 ; COEF = (0.0_R8,0.0_R8)

MAIN_LOOP: &
DO

  READ (3,NML=PROBLEM,END=500)

  IF (NEW_PROBLEM) THEN
    WRITE (7,190) TITLE,TRACKTOL,FINALTOL,SINGTOL,SSPAR(5),N
    190 FORMAT(///A80//'TRACKTOL, FINALTOL =',2ES22.14, &
      /,'SINGTOL (0 SETS DEFAULT) =',ES22.14, &
      /,'SSPAR(5) (0 SETS DEFAULT) =',ES22.14, &
      /,'NUMBER OF EQUATIONS =',I3)
    WRITE (7,200)
    200 FORMAT(/'****** COEFFICIENT TABLEAU ******')
    DO I=1,N
      WRITE (7,210) I,NUM_TERMS(I)
      210 FORMAT(/,'POLYNOMIAL(',I2,')%NUM_TERMS =',I3)
      DO J=1,NUM_TERMS(I)
        WRITE (7,220) (I,J,K,DEG(I,J,K), K=1,N)
        220 FORMAT('POLYNOMIAL(',I2,')%TERM(',I2,')%DEG(',I2,') =',I2)
        WRITE (7,230) I,J,COEF(I,J)
        230 FORMAT('POLYNOMIAL(',I2,')%TERM(',I2,')%COEF = (',ES22.14, &
            ',',ES22.14,')')
      END DO
    END DO

    ! Allocate storage for the target system in POLYNOMIAL.

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

  READ (3,NML=SYSPARTITION)

  ! Allocate storage for the system partition in PARTITION. 

  CALL CLEANUP_PAR
  ALLOCATE(PARTITION_SIZES(N))
  PARTITION_SIZES(1:N) = NUM_SETS(1:N)
  ALLOCATE(PARTITION(N))
  DO I=1,N
    ALLOCATE(PARTITION(I)%SET(PARTITION_SIZES(I)))
    DO J=1,PARTITION_SIZES(I)
      PARTITION(I)%SET(J)%NUM_INDICES = NUM_INDICES(I,J)
      ALLOCATE(PARTITION(I)%SET(J)%INDEX(NUM_INDICES(I,J)))
      PARTITION(I)%SET(J)%INDEX(1:NUM_INDICES(I,J)) = &
        INDEX(I,J,1:NUM_INDICES(I,J))
    END DO 
  END DO

  IF (ROOT_COUNT_ONLY) THEN

    ! Compute only the PLP Bezout number BPLP for this partition.

    MAXT = MAXVAL(NUM_TERMS(1:N))
    CALL BEZOUT_PLP(N,MAXT,SINGTOL,BPLP)
  ELSE

    ! Compute all BPLP roots of the target polynomial system.

    CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,IFLAG1,IFLAG2, &
      ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS)
  END IF

  WRITE (7,240) BPLP, (K,TRIM(P(K)),K=1,N)
  240 FORMAT(//,'GENERALIZED PLP BEZOUT NUMBER (BPLP) =',I10, &
    /'BASED ON THE FOLLOWING SYSTEM PARTITION:',/('P(',I2,') = ',A))

  IF (.NOT. ROOT_COUNT_ONLY) THEN
    DO M=1,BPLP
      WRITE (7,260) M,ARCLEN(M),NFE(M),IFLAG2(M)
      260 FORMAT(/'PATH NUMBER =',I10//'ARCLEN =',ES22.14/'NFE =',I5/ &
        'IFLAG2 =',I3)

      ! Designate solutions as "REAL" or "COMPLEX."

      IF (ANY(ABS(AIMAG(ROOTS(1:N,M))) >= 1.0E-4_R8)) THEN
        WRITE (7,270,ADVANCE='NO')
        270 FORMAT('COMPLEX, ')
      ELSE
        WRITE (7,280,ADVANCE='NO')
        280 FORMAT('REAL, ')
      END IF

      ! Designate solutions as "FINITE" or "INFINITE."

      IF (ABS(ROOTS(N+1,M)) < 1.0E-6_R8) THEN
        WRITE (7,290)
        290 FORMAT('INFINITE SOLUTION')
      ELSE
        WRITE (7,300)
        300 FORMAT('FINITE SOLUTION')
      END IF
      IF (MOD(IFLAG2(M),10) == 1) THEN
        WRITE (7,310) 1.0_R8,LAMBDA(M)
        310 FORMAT('LAMBDA =',ES22.14,', ESTIMATED ERROR =',ES22.14/)
      ELSE
        WRITE (7,315) LAMBDA(M)
        315 FORMAT('LAMBDA =',ES22.14/)
      END IF
      WRITE (7,320) (J,ROOTS(J,M),J=1,N)
      320 FORMAT(('X(',I2,') = (',ES22.14,',',ES22.14,')'))
      WRITE (7,330) N + 1, ROOTS(N+1,M)
      330 FORMAT(/,'X(',I2,') = (',ES22.14,',',ES22.14,')')
    END DO
  END IF

END DO MAIN_LOOP

500 CALL TEST_OPTIONS ! This tests various options, and is not part of a
                      ! typical main program.
CLOSE (UNIT=3) ; CLOSE (UNIT=7)
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

SUBROUTINE TEST_OPTIONS
IMPLICIT NONE

! Illustrate use of optional arguments NUMRR, NO_SCALING, USER_F_DF:

TRACKTOL = 1.0E-6_R8;  FINALTOL = 1.0E-8_R8
CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,IFLAG1,IFLAG2, &
  ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS, NUMRR=1, NO_SCALING=.TRUE.,  &
  USER_F_DF=.TRUE.)

M = 13
WRITE (7,FMT="(//'Testing optional arguments.')")
WRITE (7,260) M,ARCLEN(M),NFE(M),IFLAG2(M)
IF (MOD(IFLAG2(M),10) == 1) THEN
  WRITE (7,310) 1.0_R8,LAMBDA(M)
ELSE
  WRITE (7,315) LAMBDA(M)
END IF
WRITE (7,320) (J,ROOTS(J,M),J=1,N)
WRITE (7,330) N + 1, ROOTS(N+1,M)

! Now retrack one of these paths (#13) using the RECALL option:

IFLAG2(13) = -2
TRACKTOL = 1.0E-10_R8;  FINALTOL = 1.0E-14_R8
CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,IFLAG1,IFLAG2, &
  ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS, NUMRR=3, NO_SCALING=.TRUE.,  &
  USER_F_DF=.TRUE., RECALL=.TRUE.)

M = 13
WRITE (7,FMT="(//'Statistics for retracked path.')")
WRITE (7,260) M,ARCLEN(M),NFE(M),IFLAG2(M)
IF (MOD(IFLAG2(M),10) == 1) THEN
  WRITE (7,310) 1.0_R8,LAMBDA(M)
ELSE
  WRITE (7,315) LAMBDA(M)
END IF
WRITE (7,320) (J,ROOTS(J,M),J=1,N)
WRITE (7,330) N + 1, ROOTS(N+1,M)
RETURN

260 FORMAT(/'PATH NUMBER =',I10//'ARCLEN =',ES22.14/'NFE =',I5/ &
    'IFLAG2 =',I3)
310 FORMAT('LAMBDA =',ES22.14,', ESTIMATED ERROR =',ES22.14/)
315 FORMAT('LAMBDA =',ES22.14/)
320 FORMAT(('X(',I2,') = (',ES22.14,',',ES22.14,')'))
330 FORMAT(/,'X(',I2,') = (',ES22.14,',',ES22.14,')')
END SUBROUTINE TEST_OPTIONS

END PROGRAM MAIN_TEMPLATE

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
