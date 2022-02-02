! This file (sample_pmain.f95) contains a test program that calls pVTdirect,
! the parallel implementation of DIRECT, to optimize five test 
! objective functions defined in objfunc.f95.

PROGRAM SAMPLE_PMAIN
USE PVTdirect_MOD  ! The module for pVTdirect.
IMPLICIT NONE

! Local variables.
CHARACTER(len=13) :: nmlfile ! Name list file name.  
CHARACTER(len=2), DIMENSION(5) :: objnm ! Objective function short names.
INTEGER :: bin ! Number of evaluations per task.
INTEGER :: chkpt_start ! Checkpointing method option given 
  ! in the namelist file.
INTEGER :: c_switch ! Convex hull processing option. 
INTEGER :: eval_lim ! Limit on number of function evaluations.
INTEGER :: i ! Loop counter. 
INTEGER :: ierr ! Error status.
INTEGER :: iter_lim ! Limit on number of iterations.
INTEGER, PARAMETER :: MAX_N=154 ! The maximum problem dimension due to the 
  ! limited array sizes for 'LB' and 'UB' (lower and upper bounds) given in 
  ! the namelist file. Users can change it to accommodate a particular 
  ! application.
INTEGER :: N ! Problem dimension given in the namelist file.
INTEGER :: n_masters ! Number of masters.
INTEGER :: n_optbox ! Number of best boxes given in the namelist file.
INTEGER :: n_subdomains ! Number of subdomains.
INTEGER, DIMENSION(5) :: Nv ! Problem dimensions for verification.
INTEGER:: proc_id ! Processor ID.
INTEGER, DIMENSION(:), ALLOCATABLE:: status ! Return statuses from pVTdirect().
INTEGER :: test_opt ! A user-specified test option for selecting an objective 
  ! function. 
REAL(KIND = R8) :: diam_lim ! Minimum box diameter given in the namelist file.
REAL(KIND = R8) :: eps_fmin ! Tolerance (in the namelist file) defining the 
  ! minimum acceptable potential improvement in a potentially optimal box.
REAL(KIND = R8) :: FMIN ! Objective function value at X.
REAL(KIND = R8), DIMENSION(5) :: Fresults ! FMIN results for verification.
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: L, U ! Lower and upper bounds 
  ! defining n-dimensional search space.
REAL(KIND = R8), DIMENSION(MAX_N) :: LB, UB ! Lower and upper bounds given in 
  ! the namelist file.
REAL(KIND = R8) :: min_sep ! Minimum separation (in the namelist file) between 
  ! center points of the best boxes.
REAL(KIND = R8) :: objf_conv ! The smallest acceptable relative improvement 
  ! in FMIN between iterations, specified in the namelist file.
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: Wa ! Weight array.
REAL(KIND = R8), DIMENSION(MAX_N) :: weight ! A positive real array for 
  ! computing the distance between two points X and Y as: 
  ! SQRT(SUM( (X-Y)*diag(weight)*(X-Y) )).
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: X ! Minimizing vector.
REAL(KIND = R8), DIMENSION(5,5) :: Xresults ! X results for verification.
TYPE(HyperBox), DIMENSION(:), ALLOCATABLE :: box_set ! An empty array of TYPE 
  ! (Hyperbox) to be filled with as many as SIZE(box_set) best boxes if the 
  ! optional argument BOX_SET is set.

! Procedure interfaces for objective functions.
INTERFACE
  FUNCTION  Obj_GR(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END
  FUNCTION  Obj_QU(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END
  FUNCTION  Obj_RO(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END
  FUNCTION  Obj_SC(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END
  FUNCTION  Obj_MI(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END
END INTERFACE

! Name lists for parallel scheme 'PSCHEME', problem configuration 'PROBLEM', 
! optimization parameters 'OPTPARM', and checkpointing options 'CHKPTOP'.
NAMELIST /PSCHEME/ n_subdomains, n_masters, bin
NAMELIST /PROBLEM/ N, LB, UB
NAMELIST /OPTPARM/ iter_lim, eval_lim, diam_lim, objf_conv, eps_fmin, &
  c_switch, min_sep, weight, n_optbox
NAMELIST /CHKPTOP/ chkpt_start

! Manually assign 'test_opt'. 
! 1: Obj_GR
! 2: Obj_QU
! 3: Obj_RO
! 4: Obj_SC
! 5: Obj_MI
test_opt = 1

! Prepare the objective function names, problem dimensions, test results 
! for verification.
objnm(:)= (/'GR', 'QU', 'RO', 'SC', 'MI'/)
Nv(:)=(/2, 3, 4, 2, 5/)
Fresults(:) = (/0.0, -8.7558E+01, 5.4171E-08, -8.3796E+02, -4.6876/)
Xresults(:,:) = reshape((/0.0, 0.0, 0.0,0.0, 0.0,&
                3.0, 3.0, 3.0, 0.0, 0.0, &
                1.0, 1.0, 1.0, 1.0, 0.0, &
                420.9687, 420.9687, 0.0, 0.0, 0.0, &
                2.20284785, 1.57079499, 1.28493473, 1.92273352, 1.72018931/), &
                (/5,5/))
! Call pVTdirect subroutine to optimize the specified test function.
! Read the namelists for the current 'test_opt'. This serves as an
! example of varying input parameters at run-time for calling pVTdirect.
nmlfile = 'pdirect'//objnm(test_opt)//'.nml'
OPEN(UNIT=111, FILE=nmlfile, IOSTAT=ierr)
IF (ierr /= 0) THEN
  CLOSE(UNIT=111)
  WRITE (*,FMT='(3A,I3)') "Error opening file ",nmlfile,":",ierr
  STOP
END IF
READ(UNIT=111, NML=PSCHEME, IOSTAT=ierr)
IF (ierr /= 0) THEN
  CLOSE(UNIT=111)
  WRITE (*,FMT='(3A,I3)') "Error reading name list PSCHEME in ",nmlfile, &
    ":",ierr
  STOP
END IF
READ(UNIT=111, NML=PROBLEM, IOSTAT=ierr)
IF (ierr /= 0) THEN
  CLOSE(UNIT=111)
  WRITE (*,FMT='(3A,I3)') "Error reading name list PROBLEM in ",nmlfile, &
    ":",ierr
  STOP
END IF
READ(UNIT=111, NML=OPTPARM, IOSTAT=ierr)
IF (ierr /= 0) THEN
  CLOSE(UNIT=111)
  WRITE (*,FMT='(3A,I3)') "Error reading name list OPTPARM in ",nmlfile, &
    ":", ierr
  STOP
END IF
READ(UNIT=111, NML=CHKPTOP, IOSTAT=ierr)
IF (ierr /= 0) THEN
  CLOSE(UNIT=111)
  WRITE (*,FMT='(3A,I3)') "Error reading name list CHKPTOP in ",nmlfile, &
    ":", ierr
  STOP
END IF
CLOSE(UNIT=111)
! Check the sanity of the input parameters from the namelist.
IF (N <= MAX_N) THEN
  ALLOCATE(L(N))
  ALLOCATE(U(N))
  ALLOCATE(X(N))
  ALLOCATE(Wa(N))
ELSE
  STOP "Please make MAX_N larger in this main program."
END IF 
L(1:N) = LB(1:N)
U(1:N) = UB(1:N)
Wa(1:N) = weight(1:N)
IF ( (n_subdomains < 0) .OR. (n_masters < 0) ) THEN
  STOP "ERROR: n_subdomains or n_masters cannot be negative."
END IF
! Allocate 'status', an array to hold return statuses for all subdomains.
ALLOCATE(status(n_subdomains))
IF ( (objf_conv < 0.0_R8) .OR. (objf_conv > 1.0_R8) &
  .OR. (eps_fmin < 0.0_R8) .OR. (eps_fmin > 1.0_R8) ) THEN
  STOP "ERROR: objf_conv or eps_fmin is out of range."
END IF
IF ( (c_switch < 0) .OR. (c_switch > 1) ) THEN
  STOP "ERROR: c_switch is out of range."
END IF
IF ( (n_optbox < 1) .OR. ( (n_optbox > 1) .AND. &
                           (min_sep < EPSILON(0.0_R8)) ) ) THEN
  STOP "ERROR: min_sep or n_optbox is out of range."
END IF
IF ( (n_optbox > 1) .AND. (n_masters*n_subdomains > 1) ) THEN
  STOP "ERROR: multiple masters cannot output multiple best boxes."
END IF
IF ( (chkpt_start < 0) .OR. (chkpt_start > 2) ) THEN
  STOP "ERROR: chkpt_start is out of range."
END IF
SELECT CASE (test_opt)
  CASE (1) ! Griewank function (GR).
    IF (n_optbox == 1) THEN ! A single best box output is specified.
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN 
        CALL pVTdirect(N, L, U, Obj_GR, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          RESTART=chkpt_start) 
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    ELSE ! The output of multiple best boxes is specified.
      ALLOCATE(box_set(n_optbox)) ! Allocate 'box_set'.
      DO i = 1, SIZE(box_set)
        ALLOCATE(box_set(i)%c(N))
        ALLOCATE(box_set(i)%side(N))
      END DO
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_GR, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          BOX_SET=box_set, NUM_BOX=n_optbox, MIN_SEP=min_sep, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    END IF
  CASE (2) ! Quartic function (QU).
    IF (n_optbox == 1) THEN ! A single best box output is specified.
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_QU, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    ELSE ! The output of multiple best boxes is specified.
      ALLOCATE(box_set(n_optbox)) ! Allocate 'box_set'.
      DO i = 1, SIZE(box_set)
        ALLOCATE(box_set(i)%c(N))
        ALLOCATE(box_set(i)%side(N))
      END DO
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_QU, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          BOX_SET=box_set, NUM_BOX=n_optbox, MIN_SEP=min_sep, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    END IF
  CASE (3) ! Rosenbrock's valley (RO).
    IF (n_optbox == 1) THEN ! A single best box output is specified.
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_RO, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    ELSE ! The output of multiple best boxes is specified.
      ALLOCATE(box_set(n_optbox)) ! Allocate 'box_set'.
      DO i = 1, SIZE(box_set)
        ALLOCATE(box_set(i)%c(N))
        ALLOCATE(box_set(i)%side(N))
      END DO
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_RO, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          BOX_SET=box_set, NUM_BOX=n_optbox, MIN_SEP=min_sep, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    END IF
  CASE (4) ! Schwefel's function (SC).
    IF (n_optbox == 1) THEN  ! A single best box output is specified.
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_SC, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    ELSE ! The output of multiple best boxes is specified.    
      ALLOCATE(box_set(n_optbox)) ! Allocate 'box_set'.
      DO i = 1, SIZE(box_set)
        ALLOCATE(box_set(i)%c(N))
        ALLOCATE(box_set(i)%side(N))
      END DO
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_SC, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          BOX_SET=box_set, NUM_BOX=n_optbox, MIN_SEP=min_sep, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    END IF
  CASE (5) ! Michalewicz's function (MI). 
    IF (n_optbox == 1) THEN ! A single best box output is specified.
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_MI, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    ELSE ! The output of multiple best boxes is specified.
      ALLOCATE(box_set(n_optbox)) ! Allocate 'box_set'.
      DO i = 1, SIZE(box_set)
        ALLOCATE(box_set(i)%c(N))
        ALLOCATE(box_set(i)%side(N))
      END DO
      ! Initialization.
      CALL pVTdirect_init(ierr)
      IF (ierr == 0) THEN
        CALL pVTdirect(N, L, U, Obj_MI, X, FMIN, proc_id, status, &
          SWITCH=c_switch, MAX_ITER=iter_lim, MAX_EVL=eval_lim, &
          MIN_DIA=diam_lim, OBJ_CONV=objf_conv, EPS=eps_fmin, W=Wa, &
          N_SUB=n_subdomains, N_MASTER=n_masters, BINSIZE=bin, &
          BOX_SET=box_set, NUM_BOX=n_optbox, MIN_SEP=min_sep, &
          RESTART=chkpt_start)
        ! Finalization.
        CALL pVTdirect_finalize()
      ELSE
        STOP "Error in the MPI initialization."
      END IF
    END IF
END SELECT
! Check the returned 'status' at the root processor.
IF (proc_id == 0) THEN
  ! Print out the header of test results.
  WRITE (*,*) "pVTdirect test results:"  
  IF(ALL(status < 10)) THEN ! Normal return.
    WRITE (*,113) objnm(test_opt),n_masters,FMIN,diam_lim,&
                  iter_lim,eval_lim
    113 FORMAT(' Test problem',A3,                          &
      ' completed successfully with ',I3,' master(s) per subdomain.'/ &
      ' Minimum objective function value =',ES15.7,'.',/    &
      ' The minimum box diameter is',ES15.7,',',/           &
      ' the number of iterations is',I8,',',/               &
      ' and the number of objective function evaluations is',I10,'.'/)
    DO i = 1, n_subdomains
      WRITE(*, 114, ADVANCE='NO') i, status(i)
      114 FORMAT(/, " Subdomain", I3," completed with stopping rule", &
      I2, '.'/)
    END DO  
    IF ( N == Nv(test_opt) ) THEN
      ! Verify 'FMIN'. 
      IF ( ABS(FMIN - Fresults(test_opt))/&
        (1.0_R8 + ABS(Fresults(test_opt))) <= 1.0E-03_R8 ) THEN
        WRITE (*,*) "The global minimum value was found successfully."
      ELSE
        WRITE (*,*) "ERROR: Failed to locate the global minimum point."
      END IF
    ELSE
      WRITE (*,*) "N is different from the one expected for verification."
    END IF
    IF (n_optbox == 1) THEN ! Print the result of the single best box.
      WRITE (*,123,ADVANCE='NO') FMIN,X
      123 FORMAT(/," The best box with value",ES15.7," is at X ="/ &
        " (",(T3,5ES15.7))
      WRITE (*,124)
      124 FORMAT(")")
    ELSE ! Print the results of multiple best boxes.
      DO i = 1, n_optbox
        WRITE (*,133,ADVANCE='NO') i, box_set(i)%val, box_set(i)%c
        133 FORMAT(/," Best box",I3," has value",ES15.7," at X ="/ &
        " (",(T3,5ES15.7))
        WRITE (*,124)
      END DO
      DO i = 1, SIZE(box_set)
        DEALLOCATE(box_set(i)%c)
        DEALLOCATE(box_set(i)%side)
      END DO
      DEALLOCATE(box_set) ! Deallocate 'box_set'.
    END IF
    IF (N == Nv(test_opt)) THEN
      ! Verify 'X'. 
      IF ( ANY( ABS((X(1:N) - Xresults(1:N,test_opt))/(1.0_R8 + &
        ABS(Xresults(1:N,test_opt)))) > 1.0E-03_R8 ) ) THEN
        WRITE (*,FMT='(/A)') &
          " WARNING: Failed to locate the global minimum point."
      ELSE
        WRITE (*,FMT='(/A)') &
          " The global minimum point was located successfully."
      END IF
    END IF
  ELSE ! Error return.
    DO i = 1, n_subdomains
      IF (status(i) >= 10) THEN
        WRITE (*,143) objnm(test_opt), status(i), i
        143 FORMAT(' Test problem ',A2,' failed with error flag',I3, &
        ' for subdomain', I3, '.'/ &
        ' See comments in pVTdirect to interpret the flag.'//)
      END IF
    END DO
  END IF
  ! Put a line at the end.
  WRITE (*,153) test_opt
  153 FORMAT(" ----- End of test results for function",I2,". -----"//)
END IF
! Deallocate arrays.
DEALLOCATE(L)
DEALLOCATE(U)
DEALLOCATE(X)
DEALLOCATE(Wa)
DEALLOCATE(status)

END PROGRAM SAMPLE_PMAIN 

! Include the objective function definitions.
INCLUDE "objfunc.f95"
