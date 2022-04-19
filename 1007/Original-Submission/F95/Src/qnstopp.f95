! This file (qnstopp.f95) contains the module QNSTOPP_MOD that declares
! the subroutine (QNSTOPP) for the parallel QNSTOP implementation, which
! assumes OpenMP 3.1.
!
! The following LAPACK routines are directly called:
!   + DGELS, DGESV, DGETRI, DGETRF,
!   + DPOTRF,
!   + DSYEV, DSYSV,
!   + ILAENV.
!
! The following BLAS routines are directly called:
!   + DGEMV, DGEMM, DNRM2,
!   + DSYMM, DSYMV, DSYR, DSYR2.
!
MODULE QNSTOPP_MOD
USE REAL_PRECISION, ONLY: R8
USE ZIGARRAY, ONLY: RANDOM_NORMAL
!$ USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
IMPLICIT NONE

PRIVATE
PUBLIC:: LATINDESIGN, QNSTOPP, R8

CONTAINS
SUBROUTINE QNSTOPP(P, LB, UB, OBJ_FUNC, XI, MODE, FMIN, STATUS, SWITCH, &
   NSTART, START_PTS, N, MAX_ITER, MAX_EVAL, MIN_TAU, TAU, GAIN,        &
   GAMMAV, GAMMAW, ETA, TRACE, ASP_GEN, OMP)
! This is a parallel implementation of the bound constrained quasi-Newton
! stochastic optimization algorithm QNSTOP described in:
!
! B. S. Castle, Quasi-Newton Methods for Stochastic Optimization
! and Proximity-Based Methods for Disparate Information Fusion,
! Ph.D.thesis, Indiana University, Bloomington, IN, 2012.
!
! The algorithm minimizes the stochastic function f(x) inside the box
! LB <= x <= UB, returning an ellipsoid center \xi estimating the minimum
! point.
!
!
! On input:
!
! P is the dimension of LB, UB, and XI.
!
! LB(1:P) is a real array giving lower bounds on XI.
!
! UB(1:P) is a real array giving upper bounds on XI.
!
! OBJ_FUNC is the name of the real function procedure defining the
!    objective function f(x) to be minimized. OBJ_FUNC(C,IFLAG) returns
!    the value f(C) with IFLAG=0, or IFLAG/=0 if f(C) is not defined.
!    OBJ_FUNC is precisely defined in the INTERFACE block below.
!
! XI(1:P) is a real array containing the initial ellipsoid center.
!
! MODE is a CHARACTER*1 variable indicating how QNSTOPP is to be used:
!    'G' or 'g' for global optimization.
!    'S' or 's' for stochastic optimization.
!
! Optional arguments:
!
! SWITCH =
!    1   Use a single starting point for QNSTOPP, given as input XI.
!        This is the default if SWITCH is omitted.
!    2   Use multiple starting points given in the optional array argument
!        START_PTS(1:P,:), which must then be present.  XI(1:P) is not used
!        on input.
!    3   Use the number of starting points given in the optional integer
!        argument NSTART, which must then be present; the points are
!        derived by Latin hypercube sampling in the feasible box.  The input
!        XI(1:P) is the last of these points.
!    4   Use the number of starting points given in the optional integer
!        argument NSTART, which must then be present; the points are
!        computed adaptively by the user supplied real vector-valued
!        function  ASP_GEN(FIRST,P,XISTART,XIFINAL,FMIN)  described below,
!        whose name must be given as an argument.  The input XI(1:P) is the
!        first of these points.
!
! NSTART is the integer number of start points to be automatically
!    generated, as indicated by the value of SWITCH.
!
! START_PTS(1:P,:) is the assumed shape real array of start points defined
!    by the user, corresponding to SWITCH = 2.  The number of such
!    start points is inferred from the column dimension of START_PTS.
!
! N is the size of the experimental design used in each quasi-Newton
!    iteration (number of evaluations of f(x) per iteration).  The default
!    value is N = (3P+4)/2, and N must satisfy N >= P+2.
!
! MAX_ITER is the maximum number of iterations (quasi-Newton steps) per
!    start point allowed; defines stopping rule 1. If MAX_ITER is present
!    but <= 0 on input, there is no iteration limit and the total number
!    of iterations executed for all start points is returned in MAX_ITER.
!
! MAX_EVAL is the maximum number of function evaluations per start point
!    allowed; defines stopping rule 2. If MAX_EVAL is present but <= 0 on
!    input, there is no limit on the number of function evaluations per
!    start point, and the total number of function evaluations from all
!    start points is returned in MAX_EVAL.
!
! MIN_TAU is the minimum acceptable design ellipsoid radius TAU allowed;
!    defines stopping rule 3.  If MIN_TAU is present but <= 0 on input,
!    no minimum radius TAU is enforced, and the minimum final TAU from
!    all runs (one run per start point) is returned in MIN_TAU.
!
! TAU is the initial real radius for the ellipsoidal design region
!    { X \in R^p | (X - X_k)^t W_k (X - X_k) <= \tau_k^2 },
!    where X_k is the current point and W_k is a symmetric positive
!    definite scaling matrix.  If omitted, TAU defaults to 1/10 the
!    diameter of the feasible box LB <= x <= UB.
!
! GAIN is a positive real number, relevant only for MODE = 'G', that
!    defines the decay factor such that the design ellipsoid radius at
!    iteration k is \tau_k = GAIN/(GAIN + k - 1)*TAU.  If GAIN is
!    omitted, TAU is not decayed in mode 'G'.
!
! GAMMAV is a real eccentricity bound for the matrix V, an intermediate
!    matrix used for updating the scaling matrix W.  The default value is
!    20.
!
! GAMMAW is a real eccentricity bound for the scaling matrix W, precisely
!    (1/GAMMAW)*I_p <= W <= GAMMAW*I_p.  The default value is 20.
!
! ETA is a nonnegative real number, relevant only for MODE = 'S', used to
!    bound the Hessian matrix update,
!    -ETA*I <= H_k - H_{k-1} <= ETA*I.  The default value is 1.
!
! TRACE is a positive integer specifying the logical I/O unit for
!    archiving the optimization history.  If TRACE is present, for each
!    iteration for each run (different start point), the ellipsoid center
!    XI, the objective function value f(XI), the best experimental design
!    point XS, and its function value f(XS) are written to the Fortran I/O
!    unit TRACE.
!
! ASP_GEN is the name of the real vector-valued function returning the
!    next starting point Xnew = ASP_GEN(FIRST,P,XISTART,XIFINAL,FMIN).
!    FIRST = .TRUE. on the first call, and .FALSE. on successive calls.
!    P is the vector dimension.  XISTART, XIFINAL are the initial, final
!    points, respectively, from the previous completed optimization
!    with minimum FMIN.  ASP_GEN is precisely defined in the INTERFACE
!    block below.
!
! OMP =
!   0 or absent corresponds to no parallelization,
!   1 corresponds to parallelizing just the loop over the start points,
!   2 corresponds to parallelizing just the loop over the sampling point
!     objective function evaluations,
!   3 corresponds to the (nested) parallelization of both loops (maximum
!     of OMP_NUM_THREADS**2 threads, where OMP_NUM_THREADS is the OpenMP
!     variable specifying the maximum number of threads launched per
!     parallel region).
!
! On output:
!
! XI(1:P) is a real vector containing the final ellipsoid center
!    corresponding to the best predicted objective function value (for
!    mode 'S') or the actual best objective function value (for mode 'G')
!    from all runs (different start points).
!
! FMIN is the real value of the objective function at the above best
!    center XI.
!
! STATUS is a return status flag. The units decimal digit specifies
!    which stopping rule was satisfied on a successful return. The tens
!    decimal digit indicates a successful return, or an error condition
!    with the cause of the error condition reported in the units digit.
!
! Tens digit =
!  0 Normal return.
!    Units digit =
!     1   Stopping rule 1 (iteration limit) satisfied.
!     2   Stopping rule 2 (function evaluation limit) satisfied.
!     3   Stopping rule 3 (minimum TAU value reached) satisfied.
!  1 Input data error.
!    Units digit =
!     0   P < 1.
!     1   Assumed shape array LB, UB, XI, or START_PTS does not have
!         dimension P.
!     2   Some lower bound is >= the corresponding upper bound.
!     3   The value of NSTART, N, TAU, GAIN, GAMMAV, GAMMAW, or ETA is
!         invalid (unreasonable).
!     4   None of MAX_EVAL, MAX_ITER, and MIN_TAU (with MODE='S') are
!         specified; there is no stopping rule.
!     5   Invalid MODE, OMP, or SWITCH value.
!     6   The I/O unit TRACE must be positive.
!     7   START_PTS must be present when SWITCH=2, NSTART must be present
!         when SWITCH=3 or 4, and START_PTS and NSTART together are
!         incompatible.  When SWITCH=4, ASP_GEN must be present and
!         MOD(OMP,2) .NE. 1.
!  2 Memory allocation error or failure.
!    Units digit =
!     0   Work array allocation.
!     1   Latin hypercube experimental design array.
!  3 IFLAG returned from a call to OBJ_FUNC(C,IFLAG) was not zero,
!      indicating the objective function was not defined at the point C,
!      which is returned in XI.  The units digit = 0.
!  4 Some start point is not in the box [LB,UB].  The invalid start point
!      is returned in XI.  The units digit = 0.
!  5 A LAPACK routine returned an error.
!    Units digit =
!     1   DSYEV failed to converge.
!     2   DGELS was passed a matrix that does not have full rank.
!     3   DGESV performed a factorization resulting in a singular U.
!     4   DPOTRF could not complete the factorization.
!     5   DSYSV produced a factorization with a singular block diagonal.
!     6   DGETRF performed a factorization resulting in a singular U.
!     7   DGETRI was passed a singular matrix.
!
!  For example,
!     03 indicates a normal return (tens digit = 0) with "stopping rule 3
!        satisfied" (units digit = 3), and
!     12 indicates an input error (tens digit = 1) when "some lower bound
!        is >= the corresponding upper bound" (units digit = 2).
!
! Optional arguments:
!
! MAX_ITER (if present and <= 0 on input) contains the total number of
!    iterations for all start points.
!
! MAX_EVAL (if present and <= 0 on input) contains the total number of
!    function evaluations for all start points.
!
! MIN_TAU (if present and <= 0 on input) contains the minimum final TAU from
!    all runs (one run per start point).
!
!
INTEGER, INTENT(IN):: P
REAL(KIND=R8), DIMENSION(P), INTENT(IN):: LB, UB
INTERFACE
  FUNCTION OBJ_FUNC(C, IFLAG) RESULT(F)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND=R8), DIMENSION(:), INTENT(IN):: C
    INTEGER, INTENT(OUT):: IFLAG
    REAL(KIND=R8):: F
  END FUNCTION OBJ_FUNC
END INTERFACE
REAL(KIND=R8), DIMENSION(:), INTENT(IN OUT):: XI
CHARACTER(LEN=1), INTENT(IN):: MODE
REAL(KIND=R8), INTENT(OUT):: FMIN
INTEGER, INTENT(OUT):: STATUS
INTEGER, INTENT(IN), OPTIONAL:: SWITCH, NSTART
REAL(KIND=R8), DIMENSION(:,:), INTENT(IN), OPTIONAL:: START_PTS
INTEGER, INTENT(IN), OPTIONAL:: N
INTEGER, INTENT(IN OUT), OPTIONAL:: MAX_ITER, MAX_EVAL
REAL(KIND=R8), INTENT(IN OUT), OPTIONAL:: MIN_TAU
REAL(KIND=R8), INTENT(IN), OPTIONAL:: TAU, GAIN, GAMMAV, GAMMAW, ETA
INTEGER, INTENT(IN), OPTIONAL:: TRACE
INTERFACE
  FUNCTION ASP_GEN(FIRST, P, XISTART, XIFINAL, FMIN) RESULT(XNEW)
    USE REAL_PRECISION, ONLY : R8
    LOGICAL, INTENT(IN):: FIRST
    INTEGER, INTENT(IN):: P
    REAL(KIND=R8), DIMENSION(P), INTENT(IN):: XISTART, XIFINAL
    REAL(KIND=R8), INTENT(IN):: FMIN
    REAL(KIND=R8), DIMENSION(P):: XNEW
  END FUNCTION ASP_GEN
END INTERFACE
OPTIONAL:: ASP_GEN
INTEGER, INTENT(IN), OPTIONAL:: OMP

! Auxiliary functions.
INTEGER, EXTERNAL:: ILAENV ! LAPACK.
REAL(KIND=R8), EXTERNAL:: DNRM2 ! BLAS.

! Local versions of optional arguments.
INTEGER:: LMAX_EVAL, LMAX_ITER, LN, LNSTART, LOMP, LSWITCH, LTRACE
REAL(KIND=R8):: LETA, LGAMMAV, LGAMMAW, LMIN_TAU, LTAU

! Local variables.
INTEGER:: CurrentIter ! The current iteration.
INTEGER:: I, J, K, M ! Loop iteration variables.
INTEGER:: N5choose2 ! Convenience variable for (5n choose 2).
INTEGER:: SharedStatus ! Status variable shared between all threads.

REAL(KIND=R8):: Mu ! Solution to the trust region subproblem.
REAL(KIND=R8):: Mu_c, Mu_d ! Constants for the \mu update formula.
REAL(KIND=R8):: Tau_a, TauSave ! Constant and variable for updating \tau.
REAL(KIND=R8):: Tol ! Relative tolerance.

REAL(KIND=R8), DIMENSION(P):: BoxCenter ! Center of feasible box.
REAL(KIND=R8), DIMENSION(P):: BoxDiag ! Diagonal of feasible box.
REAL(KIND=R8), DIMENSION(P):: CRay ! Ray from box center to farthest vertex.
REAL(KIND=R8), DIMENSION(P):: EigVal ! Work array for eigenvalues.
REAL(KIND=R8), DIMENSION(P):: InitialXI ! The initial start point.
REAL(KIND=R8), DIMENSION(P):: InvalidXI ! Storage for failed point.
REAL(KIND=R8), DIMENSION(P+1), TARGET:: t_k ! Regression experiment.
REAL(KIND=R8), DIMENSION(P):: XiUpdate ! The update to \xi_k.
REAL(KIND=R8), DIMENSION(:), ALLOCATABLE:: YValues ! Observed function values.

REAL(KIND=R8), DIMENSION(P,P):: Hessian
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE:: LHSdesign ! Latin hypercube.
REAL(KIND=R8), DIMENSION(P,P):: W ! Scaling matrix.

! Parameters.
REAL(KIND=R8), PARAMETER:: Kappa = 0.01_R8, Plus = 1.1_R8 ! Constants for
! trust region subproblem solution.
REAL(KIND=R8), PARAMETER:: Tau_b = 7.0_R8/16.0_R8 ! Parameter for updating
! experimental design region and trust region radius.

! Temporary variables and work matrices.
INTEGER:: IWork, ThreadNum
INTEGER, DIMENSION(:,:), ALLOCATABLE:: IWork_2D
REAL(KIND=R8):: RWork_1, RWork_4, RWork_5
REAL(KIND=R8), DIMENSION(P), TARGET:: RWork1D_1, RWork1D_2
REAL(KIND=R8), DIMENSION(:), ALLOCATABLE:: RWork1D_3, RWork1D_4, ThreadFMIN
REAL(KIND=R8), DIMENSION(P,P+1):: RWork2D_1
REAL(KIND=R8), DIMENSION(P+1,P+1):: RWork2D_2
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE:: RWork2D_3
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE, TARGET:: RWork2D_4, RWork2D_5, &
  RWork2D_6, ThreadMinXi

! Pointers.
REAL(KIND=R8), DIMENSION(:), POINTER:: Gradient ! Pointer to t_k(2:P+1).
REAL(KIND=R8), DIMENSION(:), POINTER:: PrevGradient
REAL(KIND=R8), DIMENSION(:), POINTER:: SclGradient

REAL(KIND=R8), DIMENSION(:,:), POINTER:: D_kT
REAL(KIND=R8), DIMENSION(:,:), POINTER:: SclHessian
REAL(KIND=R8), DIMENSION(:,:), POINTER:: XValues

! Perform sanity check of input arguments and set local variables derived
! from input arguments.
IF (P < 1) THEN
  STATUS = 10; RETURN
END IF
IF ( (SIZE(LB) .NE. P) .OR. (SIZE(UB) .NE. P) .OR. (SIZE(XI) .NE. P) ) &
  THEN; STATUS = 11; RETURN
ELSE IF (PRESENT(START_PTS)) THEN
  IF (SIZE(START_PTS,DIM=1) .NE. P) THEN; STATUS = 11; RETURN; END IF
END IF
IF (ANY(LB .GE. UB)) THEN
  STATUS = 12; RETURN
END IF
! For all optional variables <var> except START_PTS, save a local copy in
! L<var>, or initialize L<var> with the default value for <var>.
BoxCenter(1:P) = 0.5_R8 ! Scaled box center (1/2)(LB + UB).
BoxDiag = UB - LB ! For mapping feasible box to unit hypercube.
CRay(1:P) = 0.5_R8 ! Scaled radial ray (1/2)(UB - LB).
RWork_1 = DNRM2(P,BoxDiag(1:P),1) ! Box diameter.
Tol = EPSILON(1.0_R8)*RWork_1
IF (PRESENT(NSTART)) THEN
  IF (NSTART < 2) THEN
    STATUS = 13; RETURN
  ELSE; LNSTART = NSTART
  END IF
END IF
IF (PRESENT(N)) THEN
  IF (N < P+2) THEN
    STATUS = 13; RETURN
  ELSE; LN = N
  END IF
ELSE; LN = (3*P + 4)/2
END IF
IF (PRESENT(TAU)) THEN
  IF (TAU < Tol) THEN
    STATUS = 13; RETURN
  ELSE; TauSave = TAU*SQRT(REAL(P,R8))/RWork_1
  END IF
ELSE; TauSave = 0.1_R8*RWork_1 ! 1/10 box diameter.
END IF
IF (PRESENT(GAIN)) THEN
  IF (GAIN < Tol) THEN; STATUS = 13; RETURN; END IF
END IF
IF (PRESENT(GAMMAV)) THEN
  IF (GAMMAV < 1.0_R8) THEN
    STATUS = 13; RETURN
  ELSE; LGAMMAV = GAMMAV
  END IF
ELSE; LGAMMAV = 20.0_R8
END IF
IF (PRESENT(GAMMAW)) THEN
  IF (GAMMAW < 1.0_R8) THEN
    STATUS = 13; RETURN
  ELSE; LGAMMAW = GAMMAW
  END IF
ELSE; LGAMMAW = 20.0_R8
END IF
IF (PRESENT(ETA)) THEN
  IF (ETA < 0.0_R8) THEN; STATUS = 13; RETURN; END IF
  LETA = ETA
ELSE
  LETA = 1.0_R8
END IF
IF ( ALL(MODE .NE. (/'G','g','S','s'/) ) ) THEN
  STATUS = 15; RETURN
END IF
IF (PRESENT(MAX_ITER)) THEN
  IF (MAX_ITER .LE. 0) THEN
    LMAX_ITER = 0
  ELSE; LMAX_ITER = MAX_ITER
  END IF
ELSE; LMAX_ITER = 0
END IF
IF (PRESENT(MAX_EVAL)) THEN
  IF (MAX_EVAL .LE. 0) THEN
    LMAX_EVAL = 0
  ELSE; LMAX_EVAL = MAX_EVAL
  END IF
ELSE; LMAX_EVAL = 0
END IF
IF (PRESENT(MIN_TAU)) THEN
  IF (MIN_TAU .LE. 0.0_R8) THEN
    LMIN_TAU = 0.0_R8
  ELSE; LMIN_TAU = MIN_TAU*SQRT(REAL(P,R8))/RWork_1
  END IF
ELSE; LMIN_TAU = 0.0_R8
END IF
IF ((LMAX_ITER .EQ. 0) .AND. (LMAX_EVAL .EQ. 0) .AND. &
    ((LMIN_TAU .EQ. 0.0_R8) .OR. &
    ((LMIN_TAU > 0.0_R8) .AND. (.NOT. PRESENT(GAIN)) .AND. &
    (ANY(MODE .EQ. (/'G','g'/)))))) THEN
  STATUS = 14; RETURN ! Infinite iteration loop.
END IF
IF (PRESENT(SWITCH)) THEN
  IF (SWITCH < 1 .OR. SWITCH > 4) THEN
    STATUS = 15; RETURN
  ELSE; LSWITCH = SWITCH
  END IF
ELSE; LSWITCH = 1
END IF
IF (PRESENT(TRACE)) THEN
  IF (TRACE .LE. 0) THEN
    STATUS = 16; RETURN
  ELSE; LTRACE = TRACE
  END IF
ELSE; LTRACE = 0 ! Corresponds to no historical output.
END IF
IF (PRESENT(OMP)) THEN
  IF (OMP < 0 .OR. OMP > 3) THEN; STATUS = 15; RETURN; END IF
  IF ((LSWITCH .EQ. 4) .AND. (MOD(OMP,2) .EQ. 1)) THEN
    STATUS = 17; RETURN
  END IF
!$ IF (OMP .EQ. 3) CALL OMP_SET_NESTED(1)
  LOMP = OMP
ELSE; LOMP = 0 ! No parallelization.
END IF
IF ( (LSWITCH .EQ. 2 .AND. .NOT. PRESENT(START_PTS)) .OR. &
  (ANY(LSWITCH .EQ. (/3,4/)) .AND. .NOT. PRESENT(NSTART)) .OR. &
  (PRESENT(NSTART) .AND. PRESENT(START_PTS)) .OR. &
  ((LSWITCH .EQ. 4) .AND. .NOT. PRESENT(ASP_GEN)) ) THEN
  STATUS = 17; RETURN
END IF ! Ensures compatibility between SWITCH, NSTART, and START_PTS.
IF (ANY(LSWITCH .EQ. (/1,3,4/)) .AND. &
   (ANY(XI < LB) .OR. ANY(XI > UB))) THEN
  STATUS = 40; RETURN
END IF
IF (LSWITCH .EQ. 1) LNSTART = 1
IF (LSWITCH .EQ. 2) LNSTART = SIZE(START_PTS, DIM=2)

! End of argument sanity checking and initialization of local variables
! for the optional arguments.

CALL RANDOM_SEED(SIZE=IWork)
ALLOCATE(IWork_2D(IWork,1))
IWork_2D(1:IWork,1) = 32749 ! Prime seed from POLSYS_PLP homotopy code.
CALL RANDOM_SEED(PUT=IWork_2D(1:IWork,1))
DEALLOCATE(IWork_2D)

! Allocate work arrays.
N5choose2 = (5*LN*(5*LN-1))/2
IWork = MAX(N5choose2, (ILAENV(1, 'DSYEV', 'L', LN, -1, -1, -1) + 2)*LN)
ThreadNum = 1
!$OMP PARALLEL SHARED(LNSTART,LOMP,ThreadNum) DEFAULT(none)
!$ IF ((OMP_GET_THREAD_NUM() .EQ. 0) .AND. (MOD(LOMP,2) .EQ. 1)) THEN
!$   ThreadNum = MIN(OMP_GET_NUM_THREADS(),LNSTART)
!$ END IF
!$OMP END PARALLEL
ALLOCATE(ThreadFMIN(ThreadNum), ThreadMinXi(P,ThreadNum), &
  YValues(LN), IWork_2D(N5choose2,4), RWork1D_3(5*IWork), &
  RWork1D_4(5*LN), RWork2D_3(LN,LN), RWork2D_4(P,5*LN), &
  RWork2D_5(P+1,LN), RWork2D_6(5*LN,P), STAT=STATUS)
IF (STATUS .NE. 0) THEN; STATUS = 20; RETURN; END IF

IF (LSWITCH .EQ. 3) THEN
! Initialize and generate the Latin hypercube.
  ALLOCATE(LHSdesign(1:P,1:LNSTART), STAT=STATUS)
  IF (STATUS .NE. 0) THEN; CALL CLEANUP; STATUS = 21; RETURN; END IF
  CALL LATINDESIGN(P, LNSTART, LB, UB, XI, LHSdesign)
END IF

IF (ANY(MODE .EQ. (/'S','s'/))) THEN
  Tau_a = TauSave*2.0_R8**Tau_b ! For \tau update formula.
  Mu_d = LETA*LGAMMAW + SQRT(EPSILON(1.0_R8)) ! For \mu update formula.
END IF
LTAU = TauSave

! Initialization for threads.
ThreadFMIN(:) = Huge(0.0_R8)
SharedStatus = 0

!$OMP PARALLEL DO &
!
!The SHARED list contains shared variables across all threads.
!$OMP& SHARED(BoxCenter,BoxDiag,CRay,GAIN,InitialXI,InvalidXI,LB,LHSdesign, &
!$OMP&   LETA,LGAMMAV,LGAMMAW,LMAX_EVAL,LMAX_ITER,LMIN_TAU,LN,LNSTART,      &
!$OMP&   LOMP,LSWITCH,LTRACE,MODE,Mu_d,N5choose2,P,SharedStatus,            &
!$OMP&   START_PTS,Tau_a,TauSave,ThreadFMIN,ThreadMinXi,Tol,UB)             &
!
!The FIRSTPRIVATE list contains initialized variables every thread should
!have a private copy of.
!$OMP& FIRSTPRIVATE(LTAU,Hessian,RWork_5,W,XI) &
!
!The PRIVATE list contains uninitialized variables every thread should
!have a private copy of.
!$OMP& PRIVATE(D_kT,EigVal,Gradient,IWork,IWork_2D,K,Mu,Mu_c,PrevGradient, &
!$OMP&   RWork_1,RWork_4,RWork1D_1,RWork1D_2,RWork1D_3,RWork1D_4,          &
!$OMP&   RWork2D_1,RWork2D_2,RWork2D_3,RWork2D_4,RWork2D_5,RWork2D_6,      &
!$OMP&   SclGradient,SclHessian,t_k,ThreadNum,XiUpdate,XValues,YValues)    &
!
!The LASTPRIVATE list contains private variables needed after the
!parallel region finishes.
!$OMP& LASTPRIVATE(CurrentIter,STATUS,LTAU) &
!
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& IF (MOD(LOMP,2) .EQ. 1) &
!Explicitly define the variables in the parallel region.
!$OMP& DEFAULT (NONE)
!
StartPtLoop: DO M=1,LNSTART
!$OMP ATOMIC READ
  STATUS = SharedStatus
  IF (STATUS .NE. 0) CYCLE StartPtLoop
  ThreadNum = 1
!$  ThreadNum = MOD(OMP_GET_THREAD_NUM(),LNSTART)+1

  SELECT CASE (LSWITCH)
  CASE (2)
    XI = (START_PTS(1:P,M) - LB)/BoxDiag
  CASE (3)
    XI = (LHSdesign(1:P,M) - LB)/BoxDiag
  CASE (4) ! RWork_5 = f(XI).
    SELECT CASE (M)
    CASE (2)
      XI = (ASP_GEN(.TRUE., P, LB + BoxDiag*InitialXi, LB + BoxDiag*XI, &
        RWork_5) - LB)/BoxDiag ! ASP_GEN generates points in [LB,UB].
    CASE (3:)
      XI = (ASP_GEN(.FALSE., P, LB + BoxDiag*InitialXi, LB + BoxDiag*XI, &
        RWork_5) - LB)/BoxDiag ! ASP_GEN generates points in [LB,UB].
    END SELECT
  CASE (1)
    XI = (XI - LB)/BoxDiag
  END SELECT
! Feasibility test ANY(XI<LB) .OR. ANY(XI>UB) mapped to unit hypercube.
  IF (ANY(XI < 0.0_R8) .OR. ANY(XI > 1.0_R8)) THEN
!$OMP CRITICAL
    SharedStatus = 40; InvalidXI(1:P) = XI(1:P)
!$OMP END CRITICAL
    CYCLE StartPtLoop
  END IF
! Initialize the minimum value from the thread.
  RWork_5 = OBJ_FUNC(LB + BoxDiag*XI, IWork)
  IF (IWork .NE. 0) THEN
!$OMP CRITICAL
    SharedStatus = 30; InvalidXI(1:P) = XI(1:P)
!$OMP END CRITICAL
    CYCLE StartPtLoop
  END IF
  IF (RWork_5 < ThreadFMIN(ThreadNum)) THEN
    ThreadFMIN(ThreadNum) = RWork_5
    ThreadMinXi(1:P,ThreadNum) = XI
  END IF
! Initialization for QNSTOP iteration.
  W = 0.0_R8
  FORALL (I=1:P) W(I,I) = 1.0_R8
  Hessian = W; InitialXI = XI
  CurrentIter = 0

IterationLoop: DO
! Check the stopping conditions.
  IF ((LMAX_ITER > 0 .AND. CurrentIter + 1 > LMAX_ITER)) THEN
    STATUS = 1; EXIT IterationLoop
  END IF
  IF (LMAX_EVAL > 0 .AND. (CurrentIter + 1)*(LN + 1) + 1 > LMAX_EVAL) THEN
    STATUS = 2; EXIT IterationLoop
  END IF
  IF (LTAU < LMIN_TAU) THEN; STATUS = 3; EXIT IterationLoop; END IF

! Step 1. Update iteration counts and \tau.
  CurrentIter = CurrentIter + 1
  IF (ANY(MODE .EQ. (/'S','s'/))) THEN
! For the stochastic optimization case, the radius is
! \tau_k = a(k+1)^{-b}, where a>0 and b\in(0,0.5).
    LTAU = Tau_a/(1.0_R8 + REAL(CurrentIter, R8))**Tau_b
  ELSE IF (PRESENT(GAIN)) THEN
! For the global optimization case with decay, the radius is
! \tau_k = g/(g + k - 1) \tau_0, where g>0 is the gain.
    LTAU = GAIN/(GAIN + REAL(CurrentIter - 1, R8))*TauSave
  END IF

! Step 2. Regression experiment.
!
! Step 2a.
! Construct a \Pi-poised design region
!
!   {X_{k1}, ... X_{kN}} \subset E_k(\tau_k) \cap \Theta
!
! where
!
!   E_k(\tau_k) = { X \in R^p : (X-\xi_k)^T W_k (X-\xi_k) <= \tau_k^2 }
!
! is the experimental design region and
!
!   W_k = {W_k \in R^{p \times p} : W = W_k^T, \det(W_k) = 1,
!           \gamma^{-1} I_p \preceq W_k \preceq \gamma I_p}
!
! is a scaling matrix for some \gamma \geq 1.
!
! Mathematical details are given in:
!
! B. S. Castle, Quasi-Newton Methods for Stochastic Optimization
! and Proximity-Based Methods for Disparate Information Fusion,
! Ph.D.thesis, Indiana University, Bloomington, IN, 2012.

! Generate a set X_{5N} of 5N design points uniformly sampled from the
! design region and trim 4N points near each other to give the set X_N.
! This is a heuristic to avoid solving a computationally intensive
! problem to generate N good sample points.
  RWork2D_1(1:P,1:P) = W
  CALL DSYEV('V', 'U', P, RWork2D_1(1:P,1:P), P, EigVal, RWork1D_3, &
    SIZE(RWork1D_3), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 51
    CYCLE StartPtLoop
  END IF
  FORALL (I=1:P) RWork1D_2(I) = SQRT(LTAU**2/EigVal(I)) ! = \sigma.

  DO I=1,P
    CALL RANDOM_NORMAL(RWork2D_6(1:5*LN,I))
    RWork2D_6(1:5*LN,I) = RWork2D_6(1:5*LN,I) * RWork1D_2(I)
  END DO

  RWork1D_2 = 1.0_R8/(RWork1D_2*RWork1D_2)
  CALL DGEMV('N', 5*LN, P, 1.0_R8, RWork2D_6**2, 5*LN, RWork1D_2, 1, 0.0_R8, &
    RWork1D_4, 1)
  RWork1D_4 = SQRT(RWork1D_4)

  CALL RANDOM_NUMBER(RWork1D_3(1:5*LN))
  RWork1D_3(1:5*LN) = RWork1D_3(1:5*LN)**(1.0_R8/REAL(P, KIND=R8))

  FORALL (J=1:5*LN,I=1:P)
    RWork2D_4(I,J)=(RWork2D_6(J,I)/RWork1D_4(J))*RWork1D_3(J)
  END FORALL
! RWork2D_4 = X_{5N}.

! Shift the sampled points to the design region.
  DO I=1,5*LN
    RWork1D_3(1:P) = XI
    CALL DGEMV('N', P, P, 1.0_R8, RWork2D_1(1:P,1:P), P, RWork2D_4(1:P,I), &
      1, 1.0_R8, RWork1D_3(1:P), 1)
    RWork2D_4(:,I) = RWork1D_3(1:P)
  END DO

! Calculate the Euclidean pairwise distance between elements I and J
! in X_{5N}, which is treated as m column vectors x1, x2, ..., xm.
! RWork1D_3 is arranged in the order
!   ((1,2), (1,3), ..., (1,m), (2,3), ..., (2,m), ..., (m-1,m))
  K = 1
  DO I=1,5*LN
    DO J=I+1,5*LN
      RWork1D_2(1:P) = RWork2D_4(1:P,I) - RWork2D_4(1:P,J)
      RWork1D_3(K) = DNRM2(P,RWork1D_2(1:P), 1)
      K = K + 1
    END DO
  END DO

  CALL QSORTC(RWork1D_3(1:N5choose2), IWork_2D(1:N5choose2,1))

  K = 1
  DO I=2,5*LN
    IWork_2D(K:K+5*LN-I,3) = I-1
    DO J=I,5*LN
      IWork_2D(K,2) = J
      K = K+1
    END DO
  END DO

  IWork_2D(1:5*LN,4) = 1

  I = 5*LN; J = 1
  DO WHILE (I > LN)
    IF (IWork_2D(IWork_2D(IWork_2D(J,1),2),4) .NE. 0 .AND. &
        IWork_2D(IWork_2D(IWork_2D(J,1),3),4) .NE. 0) THEN
      CALL RANDOM_NUMBER(RWork1D_3(1))
      IF (RWork1D_3(1) > 0.5) THEN
        IWork_2D(IWork_2D(IWork_2D(J,1),2),4) = 0
      ELSE
        IWork_2D(IWork_2D(IWork_2D(J,1),3),4) = 0
      END IF
      I = I-1
    END IF
    J = J+1
  END DO

  K = 1
  DO I=1,5*LN
    IF (IWork_2D(I,4) .EQ. 1) THEN
      RWork2D_4(:,K) = RWork2D_4(:,I)
      K = K+1
    END IF
  END DO

! Project to the feasible set \Theta = [0,1]^P, to which [LB,UB] has been
! mapped, thus (UB-LB)_i = 1 below.
! If X \not \in \Theta, the projection is computed by
! C + \bar t (X-C), where C is the box center and
!
!   \bar t = min{ min_{(X-C)_i>0} { (UB-LB)_i / (2 (X-C)_i) },
!                 min_{(X-C)_i<0} {-(UB-LB)_i / (2 (X-C)_i) } }.
  DO I=1,LN
    IF (ANY(RWork2D_4(:,I) > 1.0_R8) .OR. ANY(RWork2D_4(:,I) < 0.0_R8)) THEN
      RWork_4 = 1.0_R8
      DO K=1,P
        IF (ABS(RWork2D_4(K,I) - BoxCenter(K)) .GE. Tol) THEN
          RWork_4 = MIN(RWork_4,CRay(K)/ABS(RWork2D_4(K,I) - BoxCenter(K)))
        END IF
      END DO
      RWork2D_4(:,I) = BoxCenter(:) + RWork_4*(RWork2D_4(:,I) - BoxCenter(:))
    END IF
  END DO
  XValues => RWork2D_4(1:P,1:LN) ! Array of LN sampled points.

! Step 2b.
! Observe the response vector Y_k = (y_{k1}, ..., y_{kN})^T,
! where y_{ki} = f(X_{ki}) + noise.

!$OMP PARALLEL DO SHARED(BoxDiag,InvalidXI,LB,LN,P,    &
!$OMP&                   SharedStatus,XValues,YValues) &
!$OMP& PRIVATE(IWork) IF (LOMP/2 .EQ. 1) SCHEDULE(DYNAMIC) DEFAULT(NONE)
EvalLoop: DO I=1,LN
!$OMP ATOMIC READ
    IWork = SharedStatus
    IF (IWork .NE. 0) THEN
      CYCLE EvalLoop
    END IF
    YValues(I) = OBJ_FUNC(LB + BoxDiag*XValues(1:P,I), IWork)
    IF (IWork .NE. 0) THEN
!$OMP CRITICAL
      SharedStatus = 30; InvalidXI(1:P) = XValues(1:P,I)
!$OMP END CRITICAL
    END IF
  END DO EvalLoop
!$OMP END PARALLEL DO
!$OMP ATOMIC READ
  IWork = SharedStatus
  IF (IWork .NE. 0) CYCLE StartPtLoop

  IF (LTRACE > 0) THEN
    IWork = MINLOC(YValues, 1)
!$OMP CRITICAL
    WRITE (LTRACE,117) M,CurrentIter,RWork_5,YValues(IWork),LB+BoxDiag*XI(1:P)
117 FORMAT(//'Start point:',I7,', iteration:',I8,','/'f(XI) =',ES15.7, &
          ', f(XS) =',ES15.7/'Ellipsoid center XI ='/(1X,6ES12.4))
    WRITE (LTRACE,119) LB+BoxDiag*XValues(1:P,IWork)
119 FORMAT(/'Best experimental design point XS ='/(1X,6ES12.4))
!$OMP END CRITICAL
  END IF

! Step 2c.
! The response surface is modeled by the linear model
! y_{ki} = \hat f_k + (X_{ki} - \xi_k)^T \hat g_k + \epsilon_{ki},
! where \epsilon_{ki} accounts for lack of fit.
!
! The gradient \hat g_k = ((t_k)_2, ..., (t_k)_{p+1})^T is obtained by
! observing the responses and solving
!
!   D_k^T D_k t_k = D_k^T Y_k,
!
! where
!
!   D_k = [ 1, (X_{k1}-\xi_k)^T
!           ...
!           1, (X_{kN}-\xi_k)^T ].
!
  RWork2D_5(1,1:LN) = 1.0_R8
  FORALL (I=1:LN) RWork2D_5(2:P+1,I) = XValues(:,I) - XI
  NULLIFY(XValues)
  D_kT => RWork2D_5 ! Gramian matrix D_k transposed.
  RWork2D_3(1:P+1,1:LN) = D_kT
  RWork2D_6(1:LN,1) = YValues
  CALL DGELS('T', P+1, LN, 1, RWork2D_3(1:P+1,1:LN), P+1, RWork2D_6(1:LN,1), &
    LN, RWork1D_3, SIZE(RWork1D_3), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 52
    CYCLE StartPtLoop
  END IF
  IF (CurrentIter > 1) THEN
    RWork1D_2 = Gradient
    PrevGradient => RWork1D_2 ! Gradient from previous iteration.
  END IF
  t_k = RWork2D_6(1:P+1,1)
  Gradient => t_k(2:P+1) ! Current gradient estimate.

! Step 3. Secant update.
!    a. (MODE 'G') Update {\hat H_k} with the BFGS update.
!    b. (MODE 'S') Update {\hat H_k} with the constrained SR1 update.
!
! For the following, let
!   \nu_k = {\hat g_k} - {\hat g_{k-1}} and
!   s_k = \xi_k - \xi_{k-1}.
  Iteration1IF: IF (CurrentIter > 1) THEN
  IF (ANY(MODE .EQ. (/'G','g'/))) THEN
! In the global context, QNSTOP methods update the Hessian matrix
! with the BFGS update
!
!   {\hat H_k} = {\hat H_{k-1}} + U1 - U2,
!
! where
!
!   U1 = {\nu_k \nu_k^T} / {\nu_k^T s_k} and
!   U2 = {{\hat H_{k-1}} s_k s_k^T {\hat H_{k-1}}}/{s_k^T {\hat H_{k-1}} s_k}.
    RWork1D_2 = Gradient - PrevGradient
    NULLIFY(PrevGradient)
    RWork_1 = DOT_PRODUCT(RWork1D_2,XiUpdate) ! \nu_k^T s_k
! Skip Hessian update if s_k^T s_k or \nu_k^T s_k is small.
Hupdate: IF ((DOT_PRODUCT(XiUpdate,XiUpdate) > EPSILON(1.0_R8)) .AND. &
      (ABS(RWork_1) > EPSILON(1.0_R8))) THEN
    CALL DSYMV('U', P, 1.0_R8, Hessian, P, XiUpdate, 1, 0.0_R8, &
      RWork1D_3(1:P), 1)
! RWork1D_3 = \hat H_{k-1} s_k.
    CALL DSYR('U', P, -1.0_R8/DOT_PRODUCT(RWork1D_3(1:P), XiUpdate), &
      RWork1D_3(1:P), 1, Hessian, P)
! Hessian = Hessian - U2.
    CALL DSYR('U', P, 1.0_R8/RWork_1, RWork1D_2, 1, &
      Hessian, P)
! Hessian = Hessian + U1. End BFGS update.
    END IF Hupdate
  ELSE
! In the stochastic context, QNSTOP methods constrain the Hessian
! matrix update to satisfy
!
!   -\eta I_p \preceq {\hat H_k} - {\hat H_{k-1}} \preceq \eta I_p
!
! for some \eta \geq 0.
!
! Conceptually, this prevents the quadratic model from changing
! drastically from one iteration to the next.
    RWork1D_2 = Gradient - PrevGradient
    NULLIFY(PrevGradient)
    RWork2D_3(1:P,1:P) = 0.0_R8 ! RWork2D_3 will contain the Hessian update.

    IF (DOT_PRODUCT(XiUpdate, XiUpdate) &
      > EPSILON(1.0_R8)*MAX(1.0_R8, DOT_PRODUCT(XI, XI))) THEN
! Apply the SR1 update
!
!   H_k = H_{k-1} + z_k z_k^T / (z_k^T s_k),
!
! where z_k = \nu_k - H_{k-1} s_k.
      RWork1D_1 = RWork1D_2
      CALL DSYMV('U', P, -1.0_R8, Hessian, P, XiUpdate, 1, 1.0_R8, &
        RWork1D_1, 1)
! z_k = RWork1D_1(1:P).
      RWork_1 = DOT_PRODUCT(XiUpdate, RWork1D_1) ! = s_k^T z_k.
      IF (ABS(RWork_1) < EPSILON(1.0_R8) * &
        DOT_PRODUCT(RWork1D_1, RWork1D_1)) THEN
        RWork_1 = 1.0_R8/EPSILON(1.0_R8) ! SR1 update norm too large.
      ELSE
        CALL DSYR('U', P, 1.0_R8/RWork_1, RWork1D_1, 1, RWork2D_3(1:P,1:P), P)
! RWork2D_3(1:P,1:P) = z_k z_k^T /(z_k^T s_k).
        RWork_1 = DOT_PRODUCT(RWork1D_1, RWork1D_1)/ABS(RWork_1)
      END IF
      IF (RWork_1 > LETA) THEN
! If the LMI bounds are not satisfied, compute the best constrained
! solution.  The mathematical details are in:
!
! B. S. Castle, Quasi-Newton Methods for Stochastic Optimization
! and Proximity-Based Methods for Disparate Information Fusion,
! Ph.D.thesis, Indiana University, Bloomington, IN, 2012.
        RWork2D_6(1:P,1:P) = 0.0_R8
        CALL DSYR('U', P, 1.0_R8, LETA*XiUpdate, 1, RWork2D_6(1:P,1:P), P)
        CALL DSYR2('U', P, -1.0_R8, LETA*XiUpdate, 1, RWork1D_1, 1, &
          RWork2D_6(1:P,1:P), P)
! Qp = RWork2D_6(1:P,1:P) = r r^T - r z_k^T - z_k r^T,
! where r = \eta s_k.
        CALL DSYEV('V', 'U', P, RWork2D_6(1:P,1:P), P, EigVal, RWork1D_3, &
          SIZE(RWork1D_3), STATUS)
        IF (STATUS .NE. 0) THEN
!$OMP ATOMIC WRITE
          SharedStatus = 51
          CYCLE StartPtLoop
        END IF
        IWork = MINLOC(EigVal, 1)
        RWork2D_4(1:P,1:P) = 0.0_R8
        CALL DSYR('U',P,LETA,RWork2D_6(1:P,IWork),1,RWork2D_4(1:P,1:P),P)
! U_1 = RWork2D_4(1:P,1:P) = \eta v v^T, where v is the eigenvector
! corresponding to the smallest eigenvalue of Qp.

        RWork2D_6(1:P,1:P) = 0.0_R8
        CALL DSYR('U', P, 1.0_R8, LETA*XiUpdate, 1, RWork2D_6(1:P,1:P), P)
        CALL DSYR2('U', P, 1.0_R8, LETA*XiUpdate, 1, RWork1D_1, 1, &
          RWork2D_6(1:P,1:P), P)
! Qn = RWork2D_6(1:P,1:P) = r r^T + r z_k^T + z_k r^T.
        CALL DSYEV('V', 'U', P, RWork2D_6(1:P,1:P), P, EigVal, RWork1D_3, &
          SIZE(RWork1D_3), STATUS)
        IF (STATUS .NE. 0) THEN
!$OMP ATOMIC WRITE
          SharedStatus = 51
          CYCLE StartPtLoop
        END IF
        IWork = MINLOC(EigVal, 1)
        RWork2D_1 = 0.0_R8
        CALL DSYR('U', P, -LETA, RWork2D_6(1:P,IWork), 1, RWork2D_1, P)
! U_2 = RWork2D_1(1:P,1:P)  = -\eta v v^T, where v is the eigenvector
! corresponding to the smallest eigenvalue of Qn.

        RWork2D_2(1:P,1:2) = SPREAD(-RWork1D_2, DIM=2, NCOPIES=2)
        CALL DSYMV('U', P, 1.0_R8, Hessian+RWork2D_4(1:P,1:P), P, XiUpdate, &
          1, 1.0_R8, RWork2D_2(1:P,1), 1)
        CALL DSYMV('U', P, 1.0_R8, Hessian+RWork2D_1(1:P,1:P), P, XiUpdate, &
          1, 1.0_R8, RWork2D_2(1:P,2), 1)
! Choose the update U_i such that ||(\hat H_{k-1}+U_i)s_k - \nu_k|| is minimal.
        IF (DOT_PRODUCT(RWork2D_2(1:P,1), RWork2D_2(1:P,1)) .LE. &
            DOT_PRODUCT(RWork2D_2(1:P,2), RWork2D_2(1:P,2))) THEN
          RWork2D_3(1:P,1:P) = RWork2D_4(1:P,1:P)
        ELSE
          RWork2D_3(1:P,1:P) = RWork2D_1(1:P,1:P)
        END IF
      END IF
    END IF
    Hessian = Hessian + RWork2D_3(1:P,1:P)
  END IF
  END IF Iteration1IF

! Step 4. Update iterate.
!    a. (MODE 'G') Compute \mu_k by using
!       [\hat H_k + \mu_k W_k] s_k = -{\hat g_k}.
!    b. (MODE 'S') Compute \mu_k by the formula
!       \mu_k = d(c + k - 1), where c\ge 0 and d > \eta\gamma.
!
! QNSTOP methods utilize an ellipsoidal trust region concentric
! with the design region for controlling step length.
  IF (ANY(MODE .EQ. (/'G','g'/)) .OR. (CurrentIter .EQ. 1) ) THEN
! Solve the trust region subproblem with Algorithm 7.3.6 in:
!
! Conn, A. R., Gould, N. I., & Toint, P. L.,
! Trust Region Methods, SIAM, Philadelphia, 2000.

! Compute sqrt(W^{-1}) = v sqrt((diag(\lambda))^{-1}) v^T,
! where v contains the eigenvectors of W and
! \lambda contains corresponding eigenvalues.
    RWork2D_1(1:P,1:P) = W
    CALL DSYEV('V', 'U', P, RWork2D_1(1:P,1:P), P, EigVal, RWork1D_3, &
      SIZE(RWork1D_3), IWork)
    IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
      SharedStatus = 51
      CYCLE StartPtLoop
    END IF
    WHERE (EigVal(:) > EPSILON(1.0_R8))
      EigVal(:) = SQRT(1.0_R8/EigVal(:))
    END WHERE
    RWork2D_3(1:P,1:P) = RWork2D_1(1:P,1:P)
    FORALL (I=1:P) RWork2D_4(1:P,I) = RWork2D_1(1:P,I) * EigVal(I)
    CALL DGEMM('N', 'T', P, P, P, 1.0_R8, RWork2D_4(1:P,1:P), P, &
      RWork2D_3(1:P,1:P), P, 0.0_R8, &
      RWork2D_1(1:P,1:P), P) ! RWork2D_1=sqrt(W^{-1}).

! Compute scaled Hessian \bar H_k = sqrt(W^{-1}) \hat H_k sqrt(W^{-1}).
    CALL DSYMM('R', 'U', P, P, 1.0_R8, Hessian, P, RWork2D_1, P, 0.0_R8, &
      RWork2D_3(1:P, 1:P), P)
    CALL DSYMM('R', 'U', P, P, 1.0_R8, RWork2D_1, P, RWork2D_3(1:P,1:P), P, &
      0.0_R8, RWork2D_6(1:P,1:P), P)
    SclHessian => RWork2D_6(1:P,1:P) ! Scaled Hessian \bar H_k.

! Compute scaled gradient \bar g_k = sqrt(W^{-1}) \hat g_k.
    CALL DSYMV('U', P, 1.0_R8, RWork2D_1, P, Gradient, 1, 0.0_R8, &
      RWork1D_1, 1)
    SclGradient => RWork1D_1(1:P) ! Scaled gradient \bar g_k.

! Let \lambda_1 be the smallest eigenvalue of {\bar H_k}.
! If \lambda_1 < 0, set \mu_k = -Plus \lambda_1.
! Otherwise, set \mu_k = 0.
    RWork2D_4(1:P,1:P) = SclHessian
    CALL DSYEV('V', 'U', P, RWork2D_4(1:P,1:P), P, EigVal, RWork1D_3, &
      SIZE(RWork1D_3), IWork)
    IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
      SharedStatus = 51
      CYCLE StartPtLoop
    END IF

    IWork = MINLOC(EigVal, 1)
    RWork_1 = EigVal(IWork)
    IF (RWork_1 < 0.0_R8) THEN
      Mu = -Plus * RWork_1
    ELSE
      Mu = 0.0_R8
    END IF

! Let L = v sqrt(E + \mu_k).
! s is the solution to the linear system L^T X = b, where
! b is the solution to the linear system L X = - {\bar g_k}.
    EigVal = SQRT(EigVal + Mu)
    RWork2D_3(1:P,1:P) = 0.0_R8
    FORALL (I=1:P) RWork2D_3(I,I) = EigVal(I)
    CALL DGEMM('N', 'N', P, P, P, 1.0_R8, RWork2D_4(1:P,1:P), P, &
      RWork2D_3(1:P,1:P), P, 0.0_R8, RWork2D_2(1:P,1:P), P) ! RWork2D_2 = L
    RWork1D_4(1:P) = -SclGradient
    RWork2D_3(1:P,1:P) = RWork2D_2(1:P,1:P)
    CALL DGESV(P, 1, RWork2D_3(1:P,1:P), P, IWork_2D(1:P,2), RWork1D_4(1:P), &
      P, STATUS)
! RWork1D_4 = b.
    IF (STATUS .NE. 0) THEN
!$OMP ATOMIC WRITE
      SharedStatus = 53
      CYCLE StartPtLoop
    END IF
    RWork2D_3(1:P,1:P) = TRANSPOSE(RWork2D_2(1:P,1:P))
    CALL DGESV(P, 1, RWork2D_3(1:P,1:P), P, IWork_2D(1:P,2), RWork1D_4(1:P), &
      P, STATUS)
    IF (STATUS .NE. 0) THEN
!$OMP ATOMIC WRITE
      SharedStatus = 53
      CYCLE StartPtLoop
    END IF
! RWork1D_4 = s.

    RWork_1 = DNRM2(P, RWork1D_4(1:P), 1)
    IF (RWork_1 > LTAU) THEN
! Update the approximated solution to satisfy ||s|| - \tau <= Kappa \tau.
      DO WHILE (ABS(RWork_1 - LTAU) > Kappa*LTAU)
! Update \mu such that
!
!   \mu = \mu + ((||s||-\tau)/\tau)*(||s||^2)/(||w||^2),
!
! where w is the solution to the linear system LX = s.
        RWork2D_3(1:P,1:P) = RWork2D_2(1:P,1:P)
        RWork1D_2 = RWork1D_4(1:P)
        CALL DGESV(P, 1, RWork2D_3(1:P,1:P), P, IWork_2D(1:P,1), RWork1D_2, &
          P, IWork)
! w = RWork1D_2.
        IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
          SharedStatus = 53
          CYCLE StartPtLoop
        END IF
        RWork_4 = DNRM2(P, RWork1D_2(1:P), 1)
        Mu = Mu + ((RWork_1 - LTAU)/LTAU)*(RWork_1**2)/(RWork_4**2)

! Compute the Cholesky factorization L L^T = \bar H_k + \mu I.
        RWork2D_2(1:P,1:P) = SclHessian
        FORALL (I=1:P) RWork2D_2(I,I) = RWork2D_2(I,I) + Mu
        CALL DPOTRF('L', P, RWork2D_2(1:P,1:P), P, IWork)
        IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
          SharedStatus = 54
          CYCLE StartPtLoop
        END IF
        FORALL (J=2:P)
          FORALL (I=1:J-1) RWork2D_2(I,J) = 0.0_R8
        END FORALL
! L = RWork2D_2.

! Update s as the solution to the linear system L^T x = c, where
! c is the solution to the linear system L x = -{\bar g_k}.
        RWork1D_4(1:P) = -SclGradient
        RWork2D_3(1:P,1:P) = RWork2D_2(1:P,1:P)
        CALL DGESV(P, 1, RWork2D_3(1:P,1:P), P, IWork_2D(1:P,2), &
          RWork1D_4(1:P), P, IWork)
        IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
          SharedStatus = 53
          CYCLE StartPtLoop
        END IF
! c = RWork1D_4.
        RWork2D_3(1:P,1:P) = TRANSPOSE(RWork2D_2(1:P,1:P))
        CALL DGESV(P, 1, RWork2D_3(1:P,1:P), P, IWork_2D(1:P,2), &
          RWork1D_4(1:P), P, IWork)
        IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
          SharedStatus = 53
          CYCLE StartPtLoop
        END IF
! s = RWork1D_4.

        RWork_1 = DNRM2(P, RWork1D_4(1:P), 1)
      END DO
    END IF
    NULLIFY(SclHessian, SclGradient)
! Use MODE 'G' trust region value of \mu for first mode 'S' step.
    IF (ANY(MODE .EQ. (/'S','s'/))) Mu_c = Mu/Mu_d
  ELSE
    Mu = Mu_d*(Mu_c + REAL(CurrentIter - 1, KIND=R8))
  END IF

! Step 4c.
! Compute the step and update \xi_k by solving the linear system
!
!   ({\hat H_k} + \mu_k W_k) (\xi_{k+1} - \xi_k) = -{\hat g_k}.
!
  XiUpdate = -Gradient
  CALL DSYSV('U', P, 1, Hessian + Mu*W, P, IWork_2D(1:P,2), XiUpdate, &
    P, RWork1D_3, SIZE(RWork1D_3), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 55
    CYCLE StartPtLoop
  END IF

  XI = XI + XiUpdate
! Project to the feasible set \Theta = [0,1]^P, to which [LB,UB] has been
! mapped, thus (UB-LB)_i = 1 below.
! If X \not \in \Theta, the projection is computed by
! C + \bar t (X-C), where C is the box center and
!
!   \bar t = min{ min_{(X-C)_i>0} { (UB-LB)_i / (2 (X-C)_i) },
!                 min_{(X-C)_i<0} {-(UB-LB)_i / (2 (X-C)_i) } }.
  IF (ANY(XI(1:P) > 1.0_R8) .OR. ANY(XI(1:P) < 0.0_R8)) THEN
    RWork_4 = 1.0_R8
    DO K=1,P
      IF (ABS(XI(K) - BoxCenter(K)) .GE. Tol) THEN
        RWork_4 = MIN(RWork_4, CRay(K)/ABS(XI(K) - BoxCenter(K)))
      END IF
    END DO
    XI(1:P) = BoxCenter(1:P) + RWork_4*(XI(1:P) - BoxCenter(1:P))
  END IF

  RWork_5 = OBJ_FUNC(LB + BoxDiag*XI, IWork)
  IF (IWork .NE. 0) THEN
!$OMP CRITICAL
    SharedStatus = 30; InvalidXI(1:P) = XI(1:P)
!$OMP END CRITICAL
    CYCLE StartPtLoop
  END IF
  IF (RWork_5 < ThreadFMIN(ThreadNum)) THEN
    ThreadFMIN(ThreadNum) = RWork_5
    ThreadMinXi(1:P,ThreadNum) = XI
  END IF

! Step 5. Update the subsequent design ellipsoid.
!
! With the index set J={2, ..., p+1}, compute an approximation
!
!   V_k = 4 \sigma^2 (\tilde V_k)_{J,J}
!                              ...
! for the covariance matrix of \nabla \hat m_k(\xi_{k+1}-\xi_k),
! where \tilde V_k = (D_k^T D_k)^{-1}, and
!
!   \sigma^2 = {(\hat y_k - y_k)^T (\hat y_k - y_k) \over (N-p+1)}
!
! is the ordinary least squares estimate of the variance
! and \hat y_k = D_k t_k.
!
  CALL DGEMV('T', P+1, LN, 1.0_R8, D_kT, P+1, t_k, 1, 0.0_R8, &
    RWork2D_6(1:LN,1), 1)
  RWork1D_3(1:LN) = RWork2D_6(1:LN,1) - YValues
  RWork_1 = DOT_PRODUCT(RWork1D_3(1:LN), RWork1D_3(1:LN))/REAL(LN-(P+1),R8)
! RWork_1 = \sigma^2.

! Skip W update if \sigma^2 is too small.
Wupdate: IF (RWork_1 > EPSILON(1.0_R8)) THEN
  CALL DGEMM('N', 'T', P+1, P+1, LN, 1.0_R8, D_kT, P+1, D_kT, P+1, 0.0_R8, &
    RWork2D_2, P+1)
! RWork2D_2 = D_k^T D_k.
  NULLIFY(D_kT)
  CALL DGETRF(P+1, P+1, RWork2D_2, P+1, IWork_2D(1:P+1,2), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 56
    CYCLE StartPtLoop
  END IF
  CALL DGETRI(P+1, RWork2D_2, P+1, IWork_2D(1:P+1,2), RWork1D_3, &
    SIZE(RWork1D_3), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 57
    CYCLE StartPtLoop
  END IF
  RWork2D_4(1:P,1:P) = 4.0_R8*RWork_1*RWork2D_2(2:P+1,2:P+1) ! = V_k.

! Modify the eigenvalues of V_k to lie within a factor \gamma_v of the
! spectral radius \rho(V_k) by resetting small eigenvalues to
! \rho(V_k)/\gamma_v.
  RWork2D_2(1:P,1:P) = RWork2D_4(1:P,1:P)
  CALL DSYEV('V', 'U', P, RWork2D_2(1:P,1:P), P, EigVal, RWork1D_3, &
    SIZE(RWork1D_3), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 51
    CYCLE StartPtLoop
  END IF
  RWork_1 = MAXVAL(EigVal)/LGAMMAV
  WHERE (EigVal(:) < RWork_1) EigVal(:) = RWork_1

  EigVal(1:P) = 1.0_R8/EigVal(1:P)
  FORALL (I=1:P) RWork2D_3(1:P,I) = RWork2D_2(1:P,I)*EigVal(I)
  CALL DGEMM('N', 'T', P, P, P, 1.0_R8, RWork2D_3(1:P,1:P), P, &
    RWork2D_2(1:P,1:P), P, 0.0_R8, &
    RWork2D_6(1:P,1:P), P) ! RWork2D_6(1:P,1:P) = V_k^{-1}.

! Let
!
!   W_{k+1} = (\hat H_k + \mu_k W_k)^T V_k^{-1} (\hat H_k + \mu_k W_k).
!
! The ellipsoidal approximation of the 1 - \alpha percentile
! confidence set for the minimizer is given by
!
!   E_{k+1}(\chi_{p,1-\alpha}) = {X \in \realp : (X-\xi_{k+1})^T
!     W_{k+1} (X-\xi_{k+1}) \leq \chi_{p, 1-\alpha}^2 }.
!
  RWork2D_3(1:P,1:P) = Hessian + Mu*W
  CALL DSYMM('L', 'U', P, P, 1.0_R8, RWork2D_3(1:P,1:P), P, &
    RWork2D_6(1:P,1:P), P, 0.0_R8, RWork2D_4(1:P,1:P), P)
  CALL DSYMM('R', 'U', P, P, 1.0_R8, RWork2D_3(1:P,1:P), P, &
    RWork2D_4(1:P,1:P), P, 0.0_R8, W, P)
! W = W_{k+1}.

! Strictly using the updates for W_{k+1} can lead to degenerate
! ellipsoids. To obtain useful design ellipsoids, the constraints
! \gamma_w^{-1}I_p \preceq W_{k+1} \preceq \gamma_w I_p and
! \det(W_{k+1}) = 1 are enforced by modifying the eigenvalues.
!
  RWork2D_4(1:P,1:P) = W
  CALL DSYEV('V', 'U', P, RWork2D_4(1:P,1:P), P, EigVal, RWork1D_3, &
    SIZE(RWork1D_3), IWork)
  IF (IWork .NE. 0) THEN
!$OMP ATOMIC WRITE
    SharedStatus = 51
    CYCLE StartPtLoop
  END IF
  RWork_1 = MAXVAL(EigVal)/LGAMMAW
  WHERE (EigVal(:) < RWork_1) EigVal(:) = RWork_1
! Scale eigenvalues so their product is one.
  RWork_1 = 1.0_R8/PRODUCT(EigVal(:)**(1.0_R8/REAL(P,KIND=R8)))
  EigVal(:) = EigVal(:)*RWork_1

  FORALL (I=1:P) RWork2D_3(1:P,I) = RWork2D_4(1:P,I)*EigVal(I)
  CALL DGEMM('N', 'T', P, P, P, 1.0_R8, RWork2D_3(1:P,1:P), P, &
    RWork2D_4(1:P,1:P), P, 0.0_R8, W, P)
  END IF Wupdate
END DO IterationLoop

END DO StartPtLoop
!$OMP END PARALLEL DO

IF (SharedStatus .EQ. 0) THEN
  IWork = MINLOC(ThreadFMIN,1)
  FMIN = ThreadFMIN(IWork)
! Convert from [0,1]^P back to [LB,UB].
  XI = LB + BoxDiag*ThreadMinXi(1:P,IWork)
  IF (PRESENT(MAX_ITER)) THEN
    IF (MAX_ITER .LE. 0) MAX_ITER = LNSTART*CurrentIter
  END IF
  IF (PRESENT(MAX_EVAL)) THEN
    IF (MAX_EVAL .LE. 0) MAX_EVAL = LNSTART*(CurrentIter*(LN + 1) + 1)
  END IF
  IF (PRESENT(MIN_TAU)) THEN
    IF (MIN_TAU .LE. 0.0_R8) THEN
      MIN_TAU = LTAU*DNRM2(P,BoxDiag(1:P),1)/SQRT(REAL(P,R8))
    END IF
  END IF
ELSE
  STATUS = SharedStatus
! Convert from [0,1]^P back to [LB,UB].
  IF (ANY(STATUS .EQ. [30, 40])) XI(1:P) = LB + BoxDiag*InvalidXI(1:P)
END IF

CALL CLEANUP
RETURN

CONTAINS
SUBROUTINE CLEANUP ! Clean up allocated arrays before exiting.
IF (ALLOCATED(YValues)) THEN
  DEALLOCATE(YValues, IWork_2D, RWork1D_3, RWork1D_4, &
    RWork2D_3, RWork2D_4, RWork2D_5, RWork2D_6, ThreadFMIN, ThreadMinXi)
END IF
IF (ALLOCATED(LHSdesign)) DEALLOCATE(LHSdesign)
RETURN
END SUBROUTINE CLEANUP

END SUBROUTINE QNSTOPP


SUBROUTINE QSORTC(A, IDX)
! This is a QuickSort routine adapted from Orderpack 2.0.
!
! Also, this implementation incorporates ideas from "A Practical Introduction
! to Data Structures and Algorithm Analysis", by Clifford Shaffer.
!
! It sorts real numbers into ascending numerical order
! and keeps an index of the value's original array position.
!
! Author: Will Thacker, Winthrop University, July 2013.
!
! QSORTC sorts the real array A and keeps the index of the value's
! original array position along with the value (integer array IDX).
!
! On input:
!
! A(:) is the array to be sorted.
!
! On output:
!
! A(:) is sorted.
!
! IDX(1:SIZEOF(A)) contains the original positions of the sorted values.
!    I.e., sorted(i) = orginal_unsorted(IDX(i)).
!
!
REAL(KIND=R8), DIMENSION(:), INTENT(IN OUT):: A
INTEGER, DIMENSION(SIZE(A)), INTENT(OUT):: IDX

! Local variables

INTEGER:: I   ! Loop iteration variable.

! Initialize the array of original positions.
FORALL (I=1:SIZE(A)) IDX(I)=I

CALL QSORTC_HELPER(A, IDX, 1, SIZE(A))
RETURN

CONTAINS
RECURSIVE SUBROUTINE QSORTC_HELPER(A, IDX, ISTART, ISTOP)
! This internal recursive subroutine performs the recursive quicksort
! algorithm.  It is needed because the boundaries of the part of the
! array being sorted change and the initial call to the sort routine
! does not need to specify boundaries since, generally, the user will
! want to sort the entire array passed.
!
! On input:
!
! A(:) contains a subpart to be sorted.
!
! IDX(i) contains the initial position of the value A(i) before sorting.
!
! ISTART is the starting position of the subarray to be sorted.
!
! ISTOP is the ending position of the subarray to be sorted.
!
! On output:
!
! A(ISTART:ISTOP) will be sorted.
!
! IDX(i) contains the original position for the value at A(i).
!
!
REAL(KIND=R8), DIMENSION(:), INTENT (IN OUT):: A
INTEGER, DIMENSION(SIZE(A)), INTENT(IN OUT):: IDX
INTEGER, INTENT(IN):: ISTART, ISTOP

!  Local variables
INTEGER:: ILEFT ! A position on the left to be swapped with value at IRIGHT.
INTEGER:: IMID ! The middle position used to select the pivot.
INTEGER:: IRIGHT ! A position on the right to be swapped with value at ILEFT.
INTEGER:: ITEMP  ! Used for swapping within IDX.
REAL(KIND=R8):: ATEMP ! Used for swapping.
REAL(KIND=R8):: PIVOT ! Holds the temporary pivot.

! INSMAX is used to stop recursively dividing the array and to instead
! use a sort that is more efficient for small arrays than quicksort.
!
! The best cutoff point is system dependent.

INTEGER, PARAMETER:: INSMAX=24

! Check to see if we have enough values to make quicksort useful.
! Otherwise let the insertion sort handle it.

IF ((ISTOP - ISTART) < INSMAX) THEN
  CALL INSERTION(A, IDX, ISTART, ISTOP)
ELSE

! Use the median of the first, middle and last items for the pivot
! and place the median (pivot) at the end of the list.
! Putting it at the end of the list allows for a guard value to keep
! the loop from falling off the right end of the array (no need to
! check for at the end of the subarray EACH time through the loop).

  IMID = (ISTART + ISTOP)/2

  IF (A(ISTOP) < A(ISTART)) THEN
    ATEMP = A(ISTART)
    A(ISTART) = A(ISTOP)
    A(ISTOP) = ATEMP

    ITEMP = IDX(ISTART)
    IDX(ISTART) = IDX(ISTOP)
    IDX(ISTOP) = ITEMP
  END IF

  IF (A(IMID) < A(ISTOP)) THEN
    ATEMP = A(ISTOP)
    A(ISTOP) = A(IMID)
    A(IMID) = ATEMP

    ITEMP = IDX(ISTOP)
    IDX(ISTOP) = IDX(IMID)
    IDX(IMID) = ITEMP

    IF (A(ISTOP) < A(ISTART)) THEN
      ATEMP = A(ISTOP)
      A(ISTOP) = A(ISTART)
      A(ISTART) = ATEMP

      ITEMP = IDX(ISTOP)
      IDX(ISTOP) = IDX(ISTART)
      IDX(ISTART) = ITEMP
    END IF
  END IF

! Now, the first position has a value that is less or equal to the
! partition. So, we know it belongs in the left side of the partition
! and we can skip it. Also, the pivot is at the end.  So, there is
! no need to compare the pivot with itself.

  PIVOT = A(ISTOP)
  ILEFT = ISTART + 1
  IRIGHT = ISTOP - 1

  DO WHILE (ILEFT < IRIGHT)
! Find a value in the left side that is bigger than the pivot value.
! Pivot is at the right end so ILEFT will not fall off the end
! of the subarray.

    DO WHILE (A(ILEFT) < PIVOT)
      ILEFT = ILEFT + 1
    END DO

    DO WHILE (IRIGHT .NE. ILEFT)
        IF (A(IRIGHT) .LT. PIVOT) EXIT
            IRIGHT = IRIGHT - 1
    END DO

! Now we have a value bigger than pivot value on the left side that can be
! swapped with a value smaller than the pivot on the right side.
!
! This gives us all values less than pivot on the left side of the
! array and all values greater than the pivot on the right side.

    ATEMP = A(IRIGHT)
    A(IRIGHT) = A(ILEFT)
    A(ILEFT) = ATEMP

    ITEMP = IDX(IRIGHT)
    IDX(IRIGHT) = IDX(ILEFT)
    IDX(ILEFT) = ITEMP

  END DO
!
! The last swap was in error (since the while condition is not checked
! until after the swap is done) so we swap again to fix it.

! This is done (once) rather than having an if (done many times) in the
! loop to prevent the swapping.


  ATEMP = A(IRIGHT)
  A(IRIGHT) = A(ILEFT)
  A(ILEFT) = ATEMP

  ITEMP = IDX(IRIGHT)
  IDX(IRIGHT) = IDX(ILEFT)
  IDX(ILEFT) = ITEMP

! Put the pivot value in its correct spot (between the 2 partitions)
! When the WHILE condition finishes, ILEFT is greater than IRIGHT.
! So, ILEFT has the position of the first value in the right side.
! This is where we can put the pivot (and where it will finally rest,
! so no need to look at it again).  Also, place the first value of
! the right side (being displaced by the pivot) at the end of the
! subarray (since it is bigger than the pivot).

  ATEMP = A(ISTOP)
  A(ISTOP) = A(ILEFT)
  A(ILEFT) = ATEMP

  ITEMP = IDX(ISTOP)
  IDX(ISTOP) = IDX(ILEFT)
  IDX(ILEFT) = ITEMP

  CALL QSORTC_HELPER(A, IDX, ISTART, ILEFT-1)
  CALL QSORTC_HELPER(A, IDX, ILEFT+1, ISTOP)
END IF

END SUBROUTINE QSORTC_HELPER


SUBROUTINE INSERTION(A, IDX, ISTART, ISTOP)
! This subroutine performs an insertion sort used for sorting
! small subarrays efficiently.
!
! This subroutine sorts a subarray of A (between positions ISTART
! and ISTOP) keeping a record of the original position (array IDX).

! On input:
!
! A(:) contains a subpart to be sorted.
!
! IDX(i) contains the initial position of the value A(i) before sorting.
!
! ISTART is the starting position of the subarray to be sorted.
!
! ISTOP is the ending position of the subarray to be sorted.
!
! On output:
!
! A(ISTART:ISTOP) will be sorted.
!
! IDX(i) contains the original position for the value at A(i).
!

REAL(KIND=R8), DIMENSION(:), INTENT (IN OUT):: A
INTEGER, DIMENSION(SIZE(A)), INTENT(IN OUT)::IDX
INTEGER, INTENT(IN):: ISTART, ISTOP

! Local variables.

REAL(KIND=R8):: AMIN  ! Temporary minimum.
REAL(KIND=R8):: ATEMP  ! The value to be inserted.
INTEGER:: I    ! Index variable.
INTEGER:: IABOVE ! Index to find insertion point.
INTEGER:: IMIN ! Temporary minimum position.
INTEGER:: ITEMP ! Temporary for swapping.

IF (ISTOP .EQ. ISTART) THEN
  RETURN
END IF

! Find the smallest and put it at the top as a "guard" so there is
! no need for the DO WHILE to check if it is going past the top.

AMIN = A(ISTART)
IMIN = ISTART

DO I=ISTOP,ISTART+1,-1
  IF (A(I) < AMIN) THEN
    AMIN = A(I)
    IMIN = I
  END IF
END DO

A(IMIN) = A(ISTART)
A(ISTART) = AMIN

ITEMP = IDX(ISTART)
IDX(ISTART) = IDX(IMIN)
IDX(IMIN) = ITEMP

! Insertion sort the rest of the array.
DO I=ISTART+2,ISTOP
  ATEMP = A(I)
  ITEMP = IDX(I)
  IABOVE = I - 1
  IF (ATEMP < A(IABOVE)) THEN
    A(I) = A(IABOVE)
    IDX(I) = IDX(IABOVE)
    IABOVE = IABOVE - 1

! Stop moving items down when the position for "insertion" is found.
!
! Do not have to check for "falling off" the beginning of the
! array since the smallest value is a guard value in the first position.

    DO WHILE (ATEMP < A(IABOVE))
      A(IABOVE+1) = A(IABOVE)
      IDX(IABOVE+1) = IDX(IABOVE)
      IABOVE = IABOVE - 1
    END DO
  END IF
  A(IABOVE+1) = ATEMP
  IDX(IABOVE+1) = ITEMP
END DO
END SUBROUTINE INSERTION

END SUBROUTINE QSORTC

SUBROUTINE LATINDESIGN(P, LNSTART, LB, UB, XI, LHSDESIGN)
! LATINDESIGN generates a Latin hypercube experimental design with
! LNSTART points, one of which is XI, in the P-dimensional box
! B = { x | LB <= x <= UB }.  Partition each interval [LB(J),UB(J)]
! into LNSTART equal length subintervals, which partitions the box B
! into cells. A Latin hypercube experimental design has the property
! that if the LNSTART points are projected onto any coordinate
! direction K, for 1 <= K <= P, the projection has exactly one point
! in each of the subintervals of [LB(K),UB(K)].  The cells centered
! along the diameter [LB,UB] of the box B are called the diagonal cells.
!
! The algorithm is to first generate LNSTART points in the diagonal
! cells, stored as the columns of the array LHSDESIGN(1:P,1:LNSTART).
! Then swap appropriate coordinates with the last column so that the
! given point XI can be inserted as the last column.  Finally, permute
! the coordinates within rows of LHSDESIGN(1:P,1:LNSTART-1) by doing
! LNSTART-2 random interchanges in each row.
!
! On input:
!
! P is the dimension.
!
! LNSTART is the integer number of points in the experimental design.
!
! LB(1:P) is a real array giving the lower bounds.
!
! UB(1:P) is a real array giving the upper bounds, defining the box
!   { x | LB <= x <= UB } in which design is generated.
!
! XI(1:P) is a real array containing the start point to be inserted in the
!   Latin hypercube design.
!
! On output:
!
! LHSDESIGN(1:P,1:LNSTART) is a real array whose columns are the Latin
!   hypercube experimental design, with XI being the last column.

INTEGER, INTENT(IN):: P, LNSTART
REAL(KIND=R8), DIMENSION(P), INTENT(IN):: LB, UB, XI
REAL(KIND=R8), DIMENSION(P,LNSTART), INTENT(OUT):: LHSDESIGN

! Local variables.
INTEGER:: I, J, L   ! Temporary loop variable.
REAL(KIND=R8):: TEMP  ! Used for swapping coordinates.
REAL(KIND=R8), DIMENSION(P):: WORK  ! Holds the bin sizes.
! SWAP is an array of random values for picking swap rows.
REAL(KIND=R8), DIMENSION(P,LNSTART-2):: SWAP


CALL RANDOM_NUMBER(HARVEST=LHSDESIGN(1:P,1:LNSTART))
WORK = (UB - LB)/REAL(LNSTART, KIND=R8)    ! Bin sizes.

! Start with all the points in diagonal cells.
DO J=1,LNSTART
  LHSDESIGN(1:P,J) = LB(1:P) + WORK(1:P)*(REAL(J-1, KIND=R8) + &
    LHSDESIGN(1:P,J))
END DO

! Now place the given start point XI into the last column by swapping
! coordinates.
DO J=1,P
  L = MIN(INT((XI(J) - LB(J))/WORK(J)) + 1, LNSTART )
  LHSDESIGN(J,L) = LHSDESIGN(J,LNSTART)
  LHSDESIGN(J,LNSTART) = XI(J)
END DO

CALL RANDOM_NUMBER(HARVEST=SWAP(1:P, 1:LNSTART-2))
! Permute coordinates within each row of LHSDESIGN(1:P,1:LNSTART-1) by
! doing LNSTART - 2 interchanges.
DO I=1,P
  DO J=LNSTART-1,2,-1
    L = INT(SWAP(I,J-1)*REAL(J,KIND=R8)) + 1
    TEMP = LHSDESIGN(I,J)
    LHSDESIGN(I,J) = LHSDESIGN(I,L)
    LHSDESIGN(I,L) = TEMP
  END DO
END DO
RETURN
END SUBROUTINE LATINDESIGN

END MODULE QNSTOPP_MOD
