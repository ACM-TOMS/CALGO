PROGRAM SAMPLES_MAIN
! Test code that runs the VTMOP driver on one of the test problems from
! VTMOP_FUNC_MOD in serial execution mode.
!
! The outputs are printed to the file samples_out.txt

! Import the library of VTMOP interfaces.
USE VTMOP_LIB
! Import the library of sample multiobjective functions.
USE VTMOP_FUNC_MOD
IMPLICIT NONE

! Output file.
CHARACTER(LEN=20), PARAMETER :: OUTFILE = "samples_out.txt"

! Adjust the problem dimensions below. D is the number of design variables,
! and P is the number of objectives. Any values P >= 2 and D >= P are allowed.
INTEGER, PARAMETER :: D=5, P=3

! Set the following to FALSE if changing the objective function is changed.
LOGICAL, PARAMETER :: CHECK_MAE = .TRUE.

! Pointer to the objective function, defined in VTMOP_FUNC_MOD.
! Uncomment one the following.
!PROCEDURE(VTMOP_MOD_OBJ_INT), POINTER :: OBJ_F => CONVEX ! Convex Pareto front.
!PROCEDURE(VTMOP_MOD_OBJ_INT), POINTER :: OBJ_F => DTLZ1 ! Planar Pareto front.
PROCEDURE(VTMOP_MOD_OBJ_INT), POINTER :: OBJ_F => DTLZ2 ! Sphere Pareto front.
!PROCEDURE(VTMOP_MOD_OBJ_INT), POINTER :: OBJ_F => DTLZ3 ! Sphere Pareto front.
!PROCEDURE(VTMOP_MOD_OBJ_INT), POINTER :: OBJ_F => DTLZ5 ! Low-dim Pareto curve.
!PROCEDURE(VTMOP_MOD_OBJ_INT), POINTER :: OBJ_F => DTLZ7 ! Discont. Par. front.

! Local variables used for defining the problem.
INTEGER :: I ! Loop index.
INTEGER :: IERR ! Error flag.
REAL(KIND=R8) :: MU ! Working precision.
REAL(KIND=R8) :: LB(D), UB(D) ! Lower/upper bound constraints.
REAL(KIND=R8) :: OBJ_BOUNDS(P,2) ! Lower/upper objective bounds.
REAL(KIND=R8), ALLOCATABLE :: DB_X(:,:), DB_F(:,:) ! Full database
REAL(KIND=R8), ALLOCATABLE :: PARETO_X(:,:), PARETO_F(:,:) ! Solution set.

! Get the unit roundoff.
MU = EPSILON(0.0_R8)

! Set the bounds constraints.
LB(:) = 0.0_R8
UB(:) = 1.0_R8

! Don't use any objective bounds.
OBJ_BOUNDS(:,1) = -HUGE(0.0_R8)
OBJ_BOUNDS(:,2) = HUGE(0.0_R8)

! Call VTMOP driver to solve the problem.

! Do 1000 function evaluations with checkpointing.
CALL VTMOP_SOLVE( D, P, LB, UB, OBJ_F, PARETO_X, PARETO_F, IERR,   &
                  ADAPTIVE_SEARCH=.TRUE., BB_BUDGET=1000,          &
                  MAXITERS=1000, SEARCH_BUDGET=D,                  &
                  INITIAL_SBUDGET=2*D, LOPT_BUDGET=2500,           &
                  DECAY=0.5_R8, DES_TOL=SQRT(MU), EPS=SQRT(MU),    &
                  EPSW=MU**0.25_R8, OBJ_TOL=SQRT(MU),              &
                  MIN_RADF=0.02_R8, TRUST_RADF=0.2_R8,             &
                  DES_PTS=DB_X, OBJ_PTS=DB_F, LOCAL_OPT=GPS,       &
                  OBJ_BOUNDS=OBJ_BOUNDS, FIT_SURROGATES=LSHEP_FIT, &
                  EVAL_SURROGATES=LSHEP_EVAL, PFLAG=0, ICHKPT=1    )
IF (IERR .GE. 10) THEN
   WRITE(*,11) "An error occurred. IERR=", IERR; STOP; END IF

! Do another 1000 iterations in recovery mode to test the checkpoint recovery.
CALL VTMOP_SOLVE( D, P, LB, UB, OBJ_F, PARETO_X, PARETO_F, IERR,   &
                  ADAPTIVE_SEARCH=.TRUE., BB_BUDGET=2000,          &
                  MAXITERS=2000, SEARCH_BUDGET=D,                  &
                  INITIAL_SBUDGET=2*D, LOPT_BUDGET=2500,           &
                  DECAY=0.5_R8, DES_TOL=SQRT(MU), EPS=SQRT(MU),    &
                  EPSW=MU**0.25_R8, OBJ_TOL=SQRT(MU),              &
                  MIN_RADF=0.02_R8, TRUST_RADF=0.2_R8,             &
                  DES_PTS=DB_X, OBJ_PTS=DB_F, LOCAL_OPT=GPS,       &
                  OBJ_BOUNDS=OBJ_BOUNDS, FIT_SURROGATES=LSHEP_FIT, &
                  EVAL_SURROGATES=LSHEP_EVAL, PFLAG=0, ICHKPT=-1   )
IF (IERR .GE. 10) THEN
   WRITE(*,11) "An error occurred. IERR=", IERR; STOP; END IF

! Check the solution.
CALL CHECK_SOLUTION()

! Write the outputs to OUTFILE.
OPEN(99, FILE=OUTFILE, STATUS="REPLACE")

! Print the summary statistics.
WRITE(99,*) "Summary:"
WRITE(99,12) "Number of nondominated/efficient points: ", SIZE(PARETO_X,2)
WRITE(99,12) "Total number of function evaluations: ", SIZE(DB_X,2)

! Print the nondominated objective set.
WRITE(99,*)
WRITE(99,*) "Nondominated point set:"
DO I = 1, SIZE(PARETO_F, 2)
   WRITE(99,10) PARETO_F(:,I)
END DO

! Print the efficient set.
WRITE(99,*)
WRITE(99,*) "Efficient point set:"
DO I = 1, SIZE(PARETO_X, 2)
   WRITE(99,10) PARETO_X(:,I)
END DO

! Print the full database of observed objective points.
WRITE(99,*)
WRITE(99,*) "Full database of objective values observed:"
DO I = 1, SIZE(DB_F,2)
   WRITE(99,10) DB_F(:,I)
END DO

! Print the full database of design points evaluated.
WRITE(99,*)
WRITE(99,*) "Full database of evaluated design points:"
DO I = 1, SIZE(DB_X,2)
   WRITE(99,10) DB_X(:,I)
END DO

! Close the output file.
CLOSE(99)

! Free the output arrays.
DEALLOCATE(DB_F, DB_X, PARETO_F, PARETO_X)

! Output formats.
10 FORMAT(1X,5ES15.7)
11 FORMAT(1X,A,I3)
12 FORMAT(1X,A,I5)
13 FORMAT(2X,A)

CONTAINS

SUBROUTINE CHECK_SOLUTION
! Auxiliary subroutine for checking residuals against a tolerance,
! in the special case of DTLZ2.

! Local variables.
INTEGER :: N
REAL(KIND=R8) :: MAE, TOLERANCE
! N is the number of points on the Pareto front approximation.
N = SIZE(PARETO_F,2)
! The tolerance is a small number proportional to (P * D).
TOLERANCE = MAX(0.01_R8, SQRT(MU)) * REAL(P*D, KIND=R8)
! Compute the mean absolute error (MAE).
MAE = 0.0_R8
DO I = 1, N
   MAE = MAE + ABS(NORM2(PARETO_F(:,I)) - 1.0_R8)
END DO
MAE = MAE / REAL(N, KIND=R8)
! Perform sanity checks for correctness.
IF (N < P + 1) THEN
   WRITE(*,*) "The number of solution points is unreasonably low."
   WRITE(*,*) "An installation error may have occurred."
ELSE IF (CHECK_MAE .AND. MAE > TOLERANCE) THEN
   WRITE(*,*) "The MAE did not meet the tolerance."
   WRITE(*,*) "An installation error may have occurred."
ELSE
   WRITE(*,*) "The serial installation appears to be correct."
   WRITE(*,*) "The full output is contained in the file ", OUTFILE
END IF
RETURN
END SUBROUTINE CHECK_SOLUTION

END PROGRAM SAMPLES_MAIN
