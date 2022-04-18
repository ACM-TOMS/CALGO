! This file (simple_main.f95) contains a simple main program that calls
! VTdirect, the sequential implementation of DIRECT, with hardwired inputs
! to optimize the Griewank test function defined in the file objfunc.f95.
! The output produced by this program should be close to:
!
! VTdirect simple main program results:
! Test problem GR completed successfully with stopping rule 1.
! Minimum objective function value =  0.0000000E+00.
! The minimum box diameter is  6.9115407E-12,
! the number of iterations is      49,
! and the number of objective function evaluations is      3835.
!
!
! The best box with value  0.0000000E+00 is at X =
! ( -1.0514799E-08 -1.4761463E-08)

PROGRAM SIMPLE_MAIN
USE VTdirect_MOD  ! The module for VTdirect.
IMPLICIT NONE

! Local variables.
INTEGER :: c_switch ! Convex hull processing option.
INTEGER :: eval_lim ! Limit on number of function evaluations.
INTEGER :: i ! Loop counter. 
INTEGER :: iter_lim ! Limit on number of iterations.
INTEGER :: N ! Problem dimension.
INTEGER :: n_optbox ! Number of best boxes.
INTEGER :: status ! Return status from VTdirect().
REAL(KIND = R8) :: diam_lim ! Minimum box diameter.
REAL(KIND = R8) :: eps_fmin ! Tolerance defining the minimum 
  ! acceptable potential improvement in a potentially optimal box.
REAL(KIND = R8) :: FMIN ! Objective function value at X.
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: L, U ! Lower and upper bounds
  ! defining n-dimensional search space.
REAL(KIND = R8) :: min_sep ! Minimum separation between center points 
  ! of the best boxes.
REAL(KIND = R8) :: obj_dec ! Termination condition based on relative
  ! decrease in objective function.
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: Wa ! Weight array for
  ! computing the distance between two points X and Y as:
  ! SQRT(SUM( (X-Y)*diag(Wa)*(X-Y) )).
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: X ! Minimizing vector.
TYPE(HyperBox), DIMENSION(:), ALLOCATABLE :: box_set ! An empty array of TYPE
  ! (Hyperbox) to be filled with as many as SIZE(box_set) best boxes if the
  ! optional argument BOX_SET is set.

! Procedure interface for the Griewank function.
INTERFACE
  FUNCTION  Obj_GR(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END
END INTERFACE

! Prepare the hardwired inputs.
N = 2
ALLOCATE(L(N))
ALLOCATE(U(N))
ALLOCATE(X(N))
ALLOCATE(Wa(N))
L(1:N) = (/-20.0_R8, -20.0_R8/) 
U(1:N) = (/30.0_R8, 30.0_R8/)
iter_lim = 49
eval_lim = 0
diam_lim = 0.0_R8

! Initialization of useful optional arguments.
c_switch = 0; obj_dec = 0.0_R8; eps_fmin = 0.0_R8

! Optional arguments needed to return multiple best boxes.
min_sep = 1.0_R8
Wa(1:N) = (/1, 1/)
n_optbox = 1

! Print out the header of simple main program results.
WRITE (*,*) "VTdirect simple main program results:"
! Call VTdirect subroutine to optimize the Griewank function.
IF (n_optbox == 1) THEN ! A single best box output is specified.
  CALL VTdirect(N, L, U, Obj_GR, X, FMIN, status, &
    MAX_ITER=iter_lim, MAX_EVL=eval_lim, MIN_DIA=diam_lim, &
    ! The remaining optional arguments can be omitted.
    SWITCH=c_switch, OBJ_CONV=obj_dec, EPS=eps_fmin, RESTART=0)
ELSE ! The output of multiple best boxes is specified.
  ALLOCATE(box_set(n_optbox)) ! Allocate 'box_set'.
  DO i = 1, SIZE(box_set)
    ALLOCATE(box_set(i)%c(N))
    ALLOCATE(box_set(i)%side(N))
  END DO
  CALL VTdirect(N, L, U, Obj_GR, X, FMIN, status, &
    MAX_ITER=iter_lim, MAX_EVL=eval_lim, MIN_DIA=diam_lim, &
    BOX_SET=box_set, NUM_BOX=n_optbox, MIN_SEP=min_sep, W=Wa, &
    ! The remaining optional arguments can be omitted.
    SWITCH=c_switch, OBJ_CONV=obj_dec, EPS=eps_fmin, RESTART=0)
END IF

! Check the returned 'status'.
IF(status < 10) THEN ! Normal return.
  WRITE (*,113) status,FMIN,diam_lim,iter_lim,eval_lim
  113 FORMAT(' Test problem GR completed successfully with stopping rule', &
        I2,'.', / ' Minimum objective function value =',ES15.7,'.',/ &
      ' The minimum box diameter is',ES15.7,',',/           &
      ' the number of iterations is',I8,',',/               &
      ' and the number of objective function evaluations is',I10,'.'/)
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
ELSE ! Error return.
  WRITE (*,143) status
  143 FORMAT(' Test problem GR failed with error flag',I3,'.'/  &
      ' See comments in VTdirect to interpret flag.'//)
END IF

! Deallocate arrays.
DEALLOCATE(L)
DEALLOCATE(U)
DEALLOCATE(X)
DEALLOCATE(Wa)

END PROGRAM SIMPLE_MAIN

! Include the objective function definition.
INCLUDE "objfunc.f95"
