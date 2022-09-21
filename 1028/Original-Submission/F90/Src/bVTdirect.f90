! This file (bVTdirect.f90) contains the module bVTdirect_MOD that
! declares subroutines and functions used by the batched (parallel)
! implementation of VTDIRECT, which is used by VTMOP. Several modifications
! were made to the interface for the serial VTdirect code, allowing
! VTMOP to pass a P-dimensional scalarizing weight vector in a way
! that is threadsafe.
!
MODULE bVTdirect_MOD
USE VTDIRECT_COMMSUB ! Module (shared_modules.f95) for subroutines
  ! commonly used by VTdirect and pVTdirect.
CONTAINS

SUBROUTINE bVTdirect(D, P, L, U, WEIGHTS, OBJ_FUNC, X, FMIN, STATUS, SWITCH, &
                    MAX_ITER, MAX_EVL, MIN_DIA, OBJ_CONV, EPS, &
                    MIN_SEP, W, BOX_SET, NUM_BOX )
IMPLICIT NONE
! This is an OpenMP parallel implementation of the DIRECT global unconstrained
! optimization algorithm described in:
!
!    D.R. Jones, C.D. Perttunen, and B.E. Stuckman, Lipschitzian
!    optimization without the Lipschitz constant, Journal of Optimization
!    Theory and Applications, Vol. 79, No. 1, 1993, pp. 157-181.
!
! The algorithm to minimize f(x) inside the box L <= x <= U is as follows:
!
!    1. Normalize the search space to be the unit hypercube. Let c_1 be
!       the center point of this hypercube and evaluate f(c_1).
!    2. Identify the set S of potentially optimal rectangles.
!    3. For all rectangles j in S:
!       3a. Identify the set I of dimensions with the maximum side length.
!           Let delta equal one-third of this maximum side length.
!       3b. Sample the function at the points c +- delta*e_i for all i
!           in I, where c is the center of the rectangle and e_i is the
!           ith unit vector.
!       3c. Divide the rectangle containing c into thirds along the
!           dimensions in I, starting with the dimension with the lowest
!           value of f(c +- delta*e_i) and continuing to the dimension
!           with the highest f(c +- delta*e_i).
!    4. Repeat 2.-3. until stopping criterion is met.
!
! On input:
!
! D is the dimension of L, U, and X.
!
! P is the dimension of WEIGHTS.
!
! L(1:D) is a real array giving lower bounds on X.
!
! U(1:D) is a real array giving upper bounds on X.
!
! WEIGHTS(1:P) is a scalarizing weight vector that is passed to OBJ_FUNC.
!
! OBJ_FUNC is the name of the real function procedure defining the
!    objective function w^T f(x) to be minimized. OBJ_FUNC(C, SCAL_W, IFLAG)
!    returns the value SCAL_W^T f(C) with IFLAG=0, or IFLAG/=0 if f(C) is
!    not defined. OBJ_FUNC is precisely defined in the INTERFACE block below.
!
! Optional arguments:
!
! SWITCH =
!    1   select potentially optimal boxes on the convex hull of the
!        (box diameter, function value) points (default).
!    0   select as potentially optimal the box with the smallest function
!        value for each diameter that is above the roundoff level.
!        This is an aggressive selection procedure that generates many more
!        boxes to subdivide.
!
! MAX_ITER is the maximum number of iterations (repetitions of Steps 2-3)
!    allowed; defines stopping rule 1. If MAX_ITER is present but <= 0
!    on input, there is no iteration limit and the number of iterations
!    executed is returned in MAX_ITER.
!
! MAX_EVL is the maximum number of function evaluations allowed; defines
!    stopping rule 2. If MAX_EVL is present but <= 0 on input, there is no
!    limit on the number of function evaluations, which is returned in
!    MAX_EVL.
!
! MIN_DIA is the minimum box diameter allowed; defines stopping rule 3.
!    If MIN_DIA is present but <= 0 on input, a minimum diameter below
!    the roundoff level is not permitted, and the box diameter of the
!    box containing the smallest function value FMIN is returned in
!    MIN_DIA.
!
! OBJ_CONV is the smallest acceptable relative improvement in the minimum
!    objective function value 'FMIN' between iterations; defines
!    stopping rule 4. OBJ_CONV must be positive and greater than the round
!    off level.  If absent, it is taken as zero.
!
! EPS is the tolerance defining the minimum acceptable potential
!    improvement in a potentially optimal box.  Larger EPS values
!    eliminate more boxes from consideration as potentially optimal,
!    and bias the search toward exploration.  EPS must be positive and
!    greater than the roundoff level.  If absent, it is taken as
!    zero.  EPS > 0 is incompatible with SWITCH = 0.
!
! MIN_SEP is the specified minimal (weighted) distance between the
!    center points of the boxes returned in the optional array BOX_SET.
!    If absent or invalid, MIN_SEP is taken as 1/2 the (weighted) diameter
!    of the box [L, U].
!
! W(1:D) is a positive real array.  The distance between two points X
!    and Y is defined as SQRT(SUM( (X-Y)*W*(X-Y) )).  If absent, W is
!    taken as all ones.
!
! BOX_SET is an empty array (TYPE HyperBox) allocated to hold the desired
!    number of boxes.
!
! On output:
!
! X(1:D) is a real vector containing the sampled box center with the
!    minimum objective function value FMIN.
!
! FMIN is the minimum function value.
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
!     3   Stopping rule 3 (minimum diameter reached) satisfied. The
!         minimum diameter corresponds to the box for which X and
!         FMIN are returned.
!     4   Stopping rule 4 (relative change in 'FMIN') satisfied.
!  1 Input data error.
!    Units digit =
!     0   D < 2.
!     1   Assumed shape array L, U, W, or X does not have size D.
!     2   Some lower bound is >= the corresponding upper bound.
!     3   MIN_DIA, OBJ_CONV, or EPS is invalid or below the roundoff level.
!     4   None of MAX_EVL, MAX_ITER, MIN_DIA, and OBJ_CONV are specified;
!         there is no stopping rule.
!     5   Invalid SWITCH value.
!     6   SWITCH = 0 and EPS > 0 are incompatible.
!  2 Memory allocation error or failure.
!    Units digit =
!     0   BoxMatrix type allocation.
!     1   BoxLink or BoxLine type allocation.
!     2   int_vector or real_vector type allocation.
!     3   HyperBox type allocation.
!     4   BOX_SET is allocated with a wrong problem dimension.
!
!  For example, 
!     03 indicates a normal return (tens digit = 0) with "stopping rule 3
!        satisfied" (units digit = 3), and
!     12 indicates an input error (tens digit = 1) when "some lower bound
!        is >= the corresponding upper bound" (units digit = 2).
!
! Optional arguments:
!
! MAX_ITER (if present) contains the number of iterations.
!
! MAX_EVL (if present) contains the number of function evaluations.
!
! MIN_DIA (if present) contains the diameter of the box associated with
!    X and FMIN.
!
! MIN_SEP (if present) is unchanged if it was a reasonable value on
!    input. Otherwise, it is reset to the default value.
!
! W (if present) is unchanged if it was positive on input. Any
!    non-positive component is reset to one.
!
! BOX_SET (if present) is an array of TYPE (HyperBox) containing the
!    best boxes with centers separated by at least MIN_SEP.
!    The number of returned boxes NUM_BOX <= SIZE(BOX_SET) is as
!    large as possible given the requested separation.
!
! NUM_BOX (if present) is the number of boxes returned in the array
!    BOX_SET(1:).
!
INTEGER, INTENT(IN):: D, P
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: L
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: U
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: WEIGHTS
INTERFACE
  FUNCTION OBJ_FUNC(C, SCAL_W, IFLAG) RESULT(F)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: C
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: SCAL_W ! Scalarizing weights.
    INTEGER, INTENT(OUT):: IFLAG
    REAL(KIND = R8):: F
  END FUNCTION OBJ_FUNC
END INTERFACE
REAL(KIND = R8), DIMENSION(:), INTENT(OUT):: X
REAL(KIND = R8), INTENT(OUT):: FMIN
INTEGER, INTENT(OUT):: STATUS
INTEGER, INTENT(IN), OPTIONAL:: SWITCH
INTEGER, INTENT(INOUT), OPTIONAL:: MAX_ITER
INTEGER, INTENT(INOUT), OPTIONAL:: MAX_EVL
REAL(KIND = R8), INTENT(INOUT), OPTIONAL:: MIN_DIA
REAL(KIND = R8), INTENT(IN), OPTIONAL:: OBJ_CONV
REAL(KIND = R8), INTENT(IN), OPTIONAL:: EPS
REAL(KIND = R8), INTENT(INOUT), OPTIONAL:: MIN_SEP
REAL(KIND = R8), DIMENSION(:), INTENT(INOUT), OPTIONAL:: W
TYPE(HyperBox), DIMENSION(:), INTENT(INOUT), OPTIONAL:: BOX_SET
INTEGER, INTENT(OUT), OPTIONAL:: NUM_BOX

! Local variables.
INTEGER:: alloc_err ! Allocation error status.
INTEGER:: b_id ! Box matrix identifier.
INTEGER:: boxset_ind ! BOX_SET array index counter.
INTEGER:: col ! Local column index.
INTEGER:: eval_c ! Function evaluation counter.
INTEGER:: i, j, k ! Loop counters.
INTEGER:: ierr ! Error status for file I/O.
INTEGER:: iflag ! Error flag for subroutine calls.
INTEGER:: i_start ! Records the start index for searching in a node of
  ! 'setInd'.
INTEGER:: lbc ! If 1, LBC (limiting box columns) is enabled (default).
  ! If 0, LBC is disabled when the size of 'BOX_SET' is greater than 1 or
  ! MAX_ITER is not specified.
INTEGER:: stop_rule ! Bits 0, 1, 2, 3 being set correspond to stopping
  ! rules 1 (iteration limit), 2 (function evaluation limit), 3 (minimum
  ! box diameter),4 (relative change in 'FMIN') respectively.
INTEGER:: SWITCH_I ! Local copy of argument SWITCH.
INTEGER:: t ! Iteration counter.
LOGICAL:: do_it ! Sign to process first box in each column of BoxMatrix.
REAL(KIND = R8), DIMENSION(D):: current_center ! Center coordinates of
  ! the best unmarked box in the current box data structure.
REAL(KIND = R8):: dia ! Diameter squared associated with 'FMIN'.
REAL(KIND = R8):: dia_limit ! Minimum diameter permitted.
REAL(KIND = R8):: EPS_I ! Local copy of argument EPS.
REAL(KIND = R8):: fmin_old ! FMIN backup.
REAL(KIND = R8):: MIN_SEP_I ! Local copy of argument MIN_SEP.
REAL(KIND = R8), DIMENSION(D):: tmp_x ! Temporary variable for 'unit_x'.
REAL(KIND = R8), DIMENSION(D):: UmL ! An array containing U(:)-L(:).
REAL(KIND = R8), DIMENSION(D):: unit_x ! X normalized to unit hypercube.
REAL(KIND = R8), DIMENSION(D):: W_I ! Local copy of weights W.
TYPE(BoxMatrix), POINTER:: m_head ! The first box matrix.
TYPE(BoxLink), POINTER:: p_l ! Pointer to the current box link.
TYPE(BoxMatrix), POINTER:: p_b ! Pointer to box matrix.
TYPE(Hyperbox), POINTER:: p_box ! Box for the removed parent box.
TYPE(HyperBox), POINTER:: p_save ! Pointer to the saved best box.
TYPE(int_vector), POINTER:: p_start ! Pointer to the start node for
  ! searching the column with CONVEX_BIT set in 'setInd'.
TYPE(BoxLine):: setB ! Set of newly sampled boxes.
TYPE(int_vector), POINTER:: p_setInd ! Pointer to a node of 'setInd'.
TYPE(int_vector), POINTER:: setFcol ! A linked list. Each node holds free
  ! column indices in BoxMatrices.
TYPE(int_vector):: setI ! Set I of dimensions with the maximum side length.
TYPE(int_vector), POINTER:: setInd ! A linked list. Each node holds column
  ! indices corresponding to different squared diameters in 'setDia'.
TYPE(real_vector), POINTER:: setDia ! A linked list. Each node holds
  ! current different squared diameters from largest to smallest.
TYPE(ValList):: setW ! Function values for newly sampled center points.

! Perform sanity check of input arguments and set local variables derived
! from input arguments.
STATUS = 0
CALL sanitycheck(STATUS)
IF (STATUS /= 0) RETURN
! Assign 'row_w' and 'col_w' in terms of 'D'.
IF (D <= 10) THEN
  row_w = MAX(10, 2*D)
ELSE
  row_w = 17 + CEILING(LOG(REAL(D))/LOG(2.0))
END IF
IF (lbc == 1) THEN
  ! If LBC is used, limit 'row_w' to be not greater than 'MAX_ITER'.
  ! In case of MAX_ITER==1, force row_w = 2.
  IF (row_w > MAX(MAX_ITER, 2)) row_w = MAX(MAX_ITER, 2)
END IF
col_w = 35*D
! Tolerance for REAL number equality tests.
EPS4N = REAL(4*D, KIND=R8)*EPSILON(1.0_R8)

! Allocate 'setI', 'setB' and 'setW'.
ALLOCATE(setI%elements(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
NULLIFY(setI%flags)
setI%dim = 0
ALLOCATE(setB%Line(2*D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+1;RETURN;END IF
ALLOCATE(setB%dir(2*D),STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+1;RETURN;END IF
DO i = 1, 2*D
  ALLOCATE(setB%Line(i)%c(D), STAT = alloc_err)
  IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+1;RETURN;END IF
  ALLOCATE(setB%Line(i)%side(D), STAT = alloc_err)
  IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+1;RETURN;END IF
END DO
setB%ind = 0
ALLOCATE(setW%val(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
ALLOCATE(setW%dir(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
setW%dim = 0

! Allocate 'setDia', 'setInd', and 'setFcol' for the first box matrix.
ALLOCATE(setDia)
ALLOCATE(setDia%elements(col_w), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
NULLIFY(setDia%next)
setDia%id = 1
setDia%dim = 0
ALLOCATE(setInd)
ALLOCATE(setInd%elements(col_w), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
ALLOCATE(setInd%flags(col_w), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
setInd%flags(:) = 0
NULLIFY(setInd%next)
NULLIFY(setInd%prev)
setInd%id = 1
setInd%dim = 0
ALLOCATE(setFcol)
ALLOCATE(setFcol%elements(col_w), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS=ALLOC_ERROR+2;RETURN;END IF
NULLIFY(setFcol%next)
NULLIFY(setFcol%prev)
NULLIFY(setFcol%flags)
setFcol%id = 1
setFcol%dim = 0

! Allocate p_box.
ALLOCATE(p_box)
ALLOCATE(p_box%c(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS = ALLOC_ERROR+3;RETURN;END IF
ALLOCATE(p_box%side(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS = ALLOC_ERROR+3;RETURN;END IF

! Allocate tempbox.
ALLOCATE(tempbox)
ALLOCATE(tempbox%c(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS = ALLOC_ERROR+3;RETURN;END IF
ALLOCATE(tempbox%side(D), STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS = ALLOC_ERROR+3;RETURN;END IF

! Step 1: Normalization of the search space, and initialization
!         of first hyperbox.
ALLOCATE(m_head, STAT = alloc_err)
IF (alloc_err /= 0) THEN;STATUS = ALLOC_ERROR;RETURN;END IF

! Initialize 't'.
t = 0
! Initialize the data structures with the initial unit box.
iflag = 0
FMIN = HUGE(1.0_R8) ! Just so come compilers can verify INTENT(OUT)
                    ! attribute of FMIN, initialized in 'init'.
CALL init(m_head, iflag)
! Check the returned 'iflag'.
IF (iflag /= 0) THEN;STATUS = iflag;RETURN;END IF

! Initialization for main loop.
t = 1
eval_c = 1
MAIN_LOOP: DO

  ! Step 2: Identify the set S of potentially optimal boxes (rectangles).
  ! Preprocess for identifying potentially optimal hyperboxes. Find and
  ! process the hyperboxes that are on the convex hull if 'SWITCH_I'== 1;
  ! otherwise, process the first box of each column until reaching
  ! the one with 'FMIN'. Set the CONVEX_BIT in 'flags', starting
  ! from the first one in 'setInd', until reaching the column with 'FMIN'.
  p_setInd => setInd
  MARKLOOP: DO
    DO i = 1, p_setInd%dim
      p_setInd%flags(i) = IBSET(p_setInd%flags(i), CONVEX_BIT)
      ! Check if the column has reached the one with 'FMIN' that has box
      ! diameter recorded in 'dia'.
      b_id = (p_setInd%elements(i) - 1)/col_w + 1
      col = MOD(p_setInd%elements(i) - 1, col_w) + 1
      p_b => m_head
      DO j = 1, b_id - 1
        p_b => p_b%child
      END DO
      IF (dia == p_b%M(1,col)%diam) EXIT MARKLOOP
    END DO
    ! Move to the next linked node of 'setInd' if present.
    IF (ASSOCIATED(p_setInd%next)) THEN
      p_setInd => p_setInd%next
    ELSE
      ! Exit when all linked nodes of 'setInd' have been checked.
      EXIT MARKLOOP
    END IF
  END DO MARKLOOP
  IF (SWITCH_I == 1) THEN ! Convex hull processing is on.
    ! Remove the columns that are not on the convex hull of potentially
    ! optimal curve. 'p_setInd' and 'i' point to the box with 'FMIN'.
    ! Pass them to findconvex for terminating the loop of identifying
    ! boxes on convex hull when EPS_I /= 0.
    CALL findconvex(m_head, p_setInd, i, setInd)
  END IF
  ! Process the set of potentially optimal boxes. They are the first
  ! boxes of all columns with CONVEX_BIT set in 'flags' in 'setInd'.
  ! Initialize 'i_start' and 'p_start' in order to search such
  ! columns in 'setInd' starting from the smaller diameters to larger
  ! diameters so sifting up newly generated boxes is not going to override
  ! any convex hull boxes that reside on the first positions of heaps
  ! of box columns (the first boxes on columns with CONVEX_BIT set
  ! in 'flags' are potentially optimal, also called "convex hull boxes").

  ! Locate the last element in the last link of 'setInd'.
  p_start => setInd
  DO WHILE (ASSOCIATED(p_start%next))
    IF (p_start%next%dim > 0) THEN
      p_start => p_start%next
    ELSE ! Exit when the next link is empty.
      EXIT
    END IF
  END DO
  i_start = p_start%dim
  ! If OBJ_CONV is present and not zero, save 'FMIN' in 'fmin_old' to be
  ! compared with the updated 'FMIN' in subroutine sampleF.
  IF (PRESENT(OBJ_CONV)) THEN
    IF (OBJ_CONV /= 0.0_R8) fmin_old = FMIN
  END IF
  ! Loop processing any boxes in the columns with CONVEX_BIT set in
  ! 'flags' in 'setInd'.
  INNER:  DO
    do_it = .FALSE.
    ! Find such a box column in linked list 'setInd' starting from
    ! position 'i_start' in the node 'p_start'. If found, 'do_it' will
    ! be set TRUE and index 'i' and node 'p_setInd' will be returned.
    p_setInd => findcol(i_start, p_start, i, do_it)
    IF (do_it) THEN
      ! Step 3:
      ! Step 3a: Obtain the 'setI' of dimensions with the maximum
      !          side length for the first box on column
      !          'p_setInd%elements(i)', where 'i' is the index
      !          in 'setInd' for the column holding the hyperbox to
      !          subdivide.
      CALL findsetI(m_head, p_setInd%elements(i), setI)

      ! Step 3b: Sample new center points at c+delta*e_i and
      !          c-delta*e_i for all dimensions in 'setI', where
      !          c is the center of the parent hyperbox being processed,
      !          and e_i is the ith unit vector. Evaluate the objective
      !          function at new center points and keep track of current
      !          global minimum 'FMIN' and its associated 'unit_x'.
      CALL sampleP(p_setInd%elements(i), setI, m_head, setB)

      ! Obtain function evaluations for all new center points.
      CALL sampleF(setB, eval_c, ierr)
      IF (ierr /= 0) THEN;STATUS = ierr;RETURN;END IF

      ! Step 3c: Divide the hyperbox containing c into thirds along the
      !          dimensions in 'setI', starting with the dimension with
      !          the lowest function value of f(c +- delta*e_i) and
      !          continuing to the dimension with the highest function
      !          value f(c +- delta*e_i).
      CALL divide(i, p_setInd%id, m_head, setB, setDia, setInd, setFcol, &
                  p_box, setW, setI, iflag)
      IF (iflag /= 0) THEN;STATUS = ALLOC_ERROR+iflag-1;RETURN;END IF
    ELSE
      ! There are no more boxes to divide for this iteration. Release
      ! empty box columns and apply 'LBC' (if enabled) to squeeze out
      ! unnecessary memory storage.
      CALL squeeze()
      EXIT
    END IF
  END DO INNER
  ! Check stop rules:
  ! Stop rule 1: maximum iterations.
  IF (BTEST(stop_rule, 0)) THEN
    IF (t >= MAX_ITER) THEN;STATUS = 1;EXIT MAIN_LOOP;END IF
  END IF
  ! Stop rule 2: maximum evaluations.
  IF (BTEST(stop_rule, 1)) THEN
    IF (eval_c >= MAX_EVL ) THEN;STATUS = 2;EXIT MAIN_LOOP;END IF
  END IF
  ! Stop rule 3: minimum diameter.
  ! Check if minimum diameter has been reached regardless of whether
  ! MIN_DIA was specified.
  IF (sqrt(dia) <= dia_limit) THEN;STATUS = 3;EXIT MAIN_LOOP;END IF
  ! Stop rule 4: objective function convergence
  IF (PRESENT(OBJ_CONV)) THEN
    IF ((OBJ_CONV /= 0.0_R8) .AND. (fmin_old /= FMIN)) THEN
      ! 'FMIN' has been updated.
      IF (fmin_old-FMIN < (1.0_R8 + ABS(fmin_old))*OBJ_CONV) THEN
        STATUS = 4
        EXIT MAIN_LOOP
      END IF
    END IF
  END IF

  ! Update iteration counter.
  t = t + 1
END DO MAIN_LOOP

! Preparation for return to the caller.
! Scale 'unit_x' back to 'X' in original coordinates.
X = L + unit_x*UmL
! Return current diameter of the box with 'FMIN' in the original frame.
IF (PRESENT(MIN_DIA)) MIN_DIA = SQRT(dia)
! Return the total iterations and evaluations.
IF (PRESENT(MAX_ITER)) MAX_ITER = t
IF (PRESENT(MAX_EVL)) MAX_EVL = eval_c

! Find as many as SIZE(BOX_SET) best boxes.
BOXSET: IF (PRESENT(BOX_SET)) THEN
  ! Put the function value and the center coordinates of the best box
  ! into 'BOX_SET'. Initialize the 'current_center' as the best box.
  current_center = unit_x
  boxset_ind = boxset_ind + 1
  BOX_SET(boxset_ind)%val= FMIN
  ! Scale the center coordinates to the original coordinate system.
  BOX_SET(boxset_ind)%c = L + unit_x*UmL
  ! Loop to find SIZE(BOX_SET)-1 best boxes.
  OUTER1: DO k = 1, SIZE(BOX_SET)-1
    ! Set an initial value to be compared with box function values.
    ! Reuse the real variable 'fmin_old' (FMIN backup).
    fmin_old = HUGE(0.0_R8)
    ! Loop over the entire data structure starting from the box matrix
    ! 'm_head'.
    p_b => m_head
    ! Initialize the number of marked boxes (with zero diameter) 't' to be
    ! 0. Reuse the integer variable 't' (loop counter for main loop).
    t = 0
    INNER1: DO WHILE(ASSOCIATED(p_b))
      ! Check all the columns in 'p_b'.
      INNER2: DO i = 1, col_w
        DO j = 1, MOD(p_b%ind(i)-1, row_w) + 1
          ! Only in the first pass, locate and mark the box with 'FMIN'.
          IF (k == 1) THEN
            IF (ALL(unit_x == p_b%M(j,i)%c)) THEN
              ! Fill in the 'side' and 'diam' of the first best box in
              ! BOX_SET. Scale them back to the original coordinate
              ! system. Mark this box by setting its diameter zero and
              ! update 't'.
              BOX_SET(1)%side = p_b%M(j,i)%side*UmL
              BOX_SET(1)%diam = SUM(BOX_SET(1)%side**2)
              p_b%M(j,i)%diam = 0.0_R8
              t = t + 1
              CYCLE
            END IF
          END IF
          ! Process unmarked boxes (with non-zero diameter) and
          ! update 't'.
          IF (p_b%M(j,i)%diam /= 0) THEN
            ! Compute the weighted separation between 'current_center'
            ! and the center of this box 'p_b%M(j,i)'.
            tmp_x = (current_center-p_b%M(j,i)%c)*UmL
            IF (DOT_PRODUCT(tmp_x*W_I, tmp_x) < MIN_SEP_I) THEN
              ! If the separation is less than 'MIN_SEP_I', mark this
              ! box by setting its diameter zero.
              p_b%M(j,i)%diam = 0.0_R8
              t = t + 1
            ELSE
              ! If the separation >= MIN_SEP_I, compare the function value
              ! of the box with the current best value ('fmin_old').
              IF (p_b%M(j,i)%val < fmin_old) THEN
                fmin_old = p_b%M(j,i)%val
                ! Backup the pointer to the current best box.
                p_save => p_b%M(j,i)
              END IF
            END IF
          ELSE
            t = t + 1
          END IF
        END DO
        ! Repeat the above steps if any box link exists.
        IF (p_b%ind(i) > row_w) THEN
          ! There must be box link(s).
          p_l => p_b%sibling(i)%p
          DO WHILE(ASSOCIATED(p_l))
            DO j = 1, p_l%ind
              IF (p_l%Line(j)%diam /= 0) THEN
                tmp_x = (current_center - p_l%Line(j)%c)*UmL
                IF (DOT_PRODUCT(tmp_x*W_I, tmp_x)< MIN_SEP_I) THEN
                  p_l%Line(j)%diam = 0.0_R8
                  t = t + 1
                ELSE
                  IF (p_l%Line(j)%val < fmin_old) THEN
                    fmin_old = p_l%Line(j)%val
                    p_save => p_l%Line(j)
                  END IF
                END IF
              ELSE
                t = t + 1
              END IF
            END DO
            ! Move to the next box link, if any.
            p_l => p_l%next
          END DO
        END IF
      END DO INNER2
      ! Move to the next box matrix, if any.
      p_b => p_b%child
    END DO INNER1
    IF (ASSOCIATED(p_save)) THEN
      IF(p_save%diam /= 0) THEN
        ! Found the next best box. Put it into BOX_SET and mark it.
        boxset_ind = boxset_ind + 1
        BOX_SET(boxset_ind) = p_save
        ! Scale the coordinates back to the original coordinate system.
        BOX_SET(boxset_ind)%c = L + BOX_SET(boxset_ind)%c*UmL
        BOX_SET(boxset_ind)%side = BOX_SET(boxset_ind)%side*UmL
        BOX_SET(boxset_ind)%diam = SUM(BOX_SET(boxset_ind)%side**2)
        p_save%diam = 0
        t = t + 1
        ! Update 'current_center'.
        current_center = p_save%c
      END IF
    ELSE ! Exit when the next best box is not available.
      EXIT OUTER1
    END IF
    ! If all boxes in the data structure are marked, exit.
    IF (eval_c <= t) EXIT OUTER1
  END DO OUTER1
  ! If NUM_BOX or MIN_SEP is specified as an optional input
  ! argument, assign a value to it.
  IF(PRESENT(NUM_BOX)) NUM_BOX = boxset_ind
  IF(PRESENT(MIN_SEP)) MIN_SEP = SQRT(MIN_SEP_I)
END IF BOXSET

! Deallocate all the data structures explicitly allocated, including
! box matrices, box links and setI, setB, setW, setInd, setFcol, setDia,
! and p_box.
CALL cleanup()

RETURN
CONTAINS

SUBROUTINE cleanup()
IMPLICIT NONE
! Cleans up all data structures allocated to prevent a memory leak.
!
! On input: None.
!
! On output: None.
!
! Local variables.
INTEGER:: i, j  ! Loop counters.
TYPE(BoxMatrix), POINTER:: p_b, p_bm
TYPE(BoxLink), POINTER:: p_l
TYPE(int_vector), POINTER:: p_seti
TYPE(real_vector), POINTER:: p_setr

! Deallocate box links and box matrices starting from the first box
! matrix. First deallocate all box links associated with each box
! matrix, and finally deallocate the box matrix.
p_b => m_head
! Check all columns with box links that will be deallocated one by one
! starting from the last box link.
DO WHILE(ASSOCIATED(p_b))
  ! Check all the columns in 'p_b'.
  DO i = 1, col_w
    IF (p_b%ind(i) > row_w) THEN
      ! There must be box link(s). Chase to the last one and start
      ! deallocating them one by one.
      p_l => p_b%sibling(i)%p
      DO WHILE(ASSOCIATED(p_l%next))
       p_l => p_l%next
      END DO
      ! Found the last box link 'p_l'. Trace back and deallocate all
      ! links.
      DO WHILE(ASSOCIATED(p_l))
        IF (ASSOCIATED(p_l%prev)) THEN
          ! 'p_l's previous link is still a box link.
          p_l => p_l%prev
        ELSE
          ! There is no box link before 'p_l'. This is the first box link
          ! of this column.
          DO j = 1, row_w
            DEALLOCATE(p_l%Line(j)%c)
            DEALLOCATE(p_l%Line(j)%side)
          END DO
          DEALLOCATE(p_l%Line)
          DEALLOCATE(p_l)
          EXIT
        END IF
        DO j = 1, row_w
          DEALLOCATE(p_l%next%Line(j)%c)
          DEALLOCATE(p_l%next%Line(j)%side)
        END DO
        DEALLOCATE(p_l%next%Line)
        DEALLOCATE(p_l%next)
      END DO
    END IF
  END DO
  ! Save the pointer of this box matrix for deallocation.
  p_bm => p_b
  ! Before it's deallocated, move to the next box matrix.
  p_b => p_b%child
  ! Deallocate this box matrix with all box links cleaned up.
  DEALLOCATE(p_bm%ind)
  DEALLOCATE(p_bm%sibling)
  DO i = 1, row_w
    DO j = 1, col_w
      DEALLOCATE(p_bm%M(i,j)%c)
      DEALLOCATE(p_bm%M(i,j)%side)
    END DO
  END DO
  DEALLOCATE(p_bm%M)
  DEALLOCATE(p_bm)
END DO

! Deallocate 'setI', 'setB' and 'setW'.
DEALLOCATE(setI%elements)
DO i = 1, 2*D
  DEALLOCATE(setB%Line(i)%c)
  DEALLOCATE(setB%Line(i)%side)
END DO
DEALLOCATE(setB%Line)
DEALLOCATE(setB%dir)
DEALLOCATE(setW%val)
DEALLOCATE(setW%dir)

! Deallocate nodes of 'setDia', 'setInd' and 'setFcol' starting from
! the last node.
p_setr => setDia
DO WHILE(ASSOCIATED(p_setr%next))
  p_setr => p_setr%next
END DO
! Found the last link pointed to by 'p_setr' of 'setDia', so deallocate
! links one by one until reaching the head node that has a null 'prev'.
DO
  DEALLOCATE(p_setr%elements)
  IF (p_setr%id /= 1) THEN
    p_setr => p_setr%prev
    DEALLOCATE(p_setr%next)
  ELSE
    DEALLOCATE(setDia)
    EXIT
  END IF
END DO
p_seti => setInd
DO WHILE(ASSOCIATED(p_seti%next))
  p_seti => p_seti%next
END DO
! Found the last link pointed to by 'p_seti' of 'setInd', so deallocate
! links one by one until reaching the head node that has a null 'prev'.
DO
  DEALLOCATE(p_seti%elements)
  DEALLOCATE(p_seti%flags)
  IF (p_seti%id /= 1) THEN
    p_seti => p_seti%prev
    DEALLOCATE(p_seti%next)
  ELSE
    DEALLOCATE(setInd)
    EXIT
  END IF
END DO
p_seti => setFcol
DO WHILE(ASSOCIATED(p_seti%next))
  p_seti => p_seti%next
END DO
! Found the last link pointed to by 'p_seti' of 'setFcol', so deallocate
! links one by one until reaching the head node that has a null 'prev'.
DO
  DEALLOCATE(p_seti%elements)
  IF (p_seti%id /= 1) THEN
    p_seti => p_seti%prev
    DEALLOCATE(p_seti%next)
  ELSE
    DEALLOCATE(setFcol)
    EXIT
  END IF
END DO

! Deallocate p_box.
DEALLOCATE(p_box%c)
DEALLOCATE(p_box%side)
DEALLOCATE(p_box)

! Deallocate tempbox
DEALLOCATE(tempbox%c)
DEALLOCATE(tempbox%side)
DEALLOCATE(tempbox)
RETURN
END SUBROUTINE cleanup

SUBROUTINE divide(parent_i, id, b, setB, setDia, setInd, setFcol, &
                  p_box, setW, setI, iflag)
IMPLICIT NONE
! Divide the first box on a column of one of box matrices 'b', starting
! from the dimension with minimum w to the one with maximum w, where w is
! min{f(c+delta), f(c-delta)}.
!
! On input:
! parent_i - The index in one of nodes of 'setInd' and 'setDia' for the
!            column holding the parent box to divide. Each element in
!            'setInd' has the same index as the one in 'setDia'.
! id       - The identifier of the node of type 'setInd'.
! b        - The head link of box matrices.
! setB     - A set of 'HyperBox' type structures, each with newly sampled
!            center point coordinates and the corresponding function
!            value. After dividing, it contains complete boxes with
!            associated side lengths and the squared diameters.
! setDia   - A linked list of current different squared diameters of box
!            matrices. It's sorted from the biggest to the smallest.
! setInd   - A linked list of column indices corresponding to the
!            different squared diameters in 'setDia'.
! setFcol  - A linked list of free columns in box matrices.
! p_box    - A 'HyperBox' type structure to hold removed parent box to
!            subdivide.
! setW     - A set of type 'ValList' used to sort wi values, where wi is
!            defined as min{f(c+delta*e_i), f(c-delta*e_i)}, the minimum of
!            function values at the two newly sampled points.
!
! On output:
! b        - 'b' has the parent box removed and contains the newly formed
!            boxes after dividing the parent box.
! setB     - Cleared set of type 'BoxLine'. All newly formed boxes have
!            been inserted to 'b'.
! setDia   - Updated linked list 'setDia' with new squared diameters of
!            boxes, if any.
! setInd   - Updated linked list 'setInd' with new column indices
!            corresponding to newly added squared diameters in 'setDia'.
! setFcol  - Updated linked list 'setFcol' with current free columns in
!           'b'.
! p_box    - A 'HyperBox' structure holding removed parent box to
!            subdivide.
! setW     - 'setW' becomes empty after dividing.
! setI     - A set of dimensions with the order of dimensions for
!            dividing. It is cleared after dividing.
! iflag    - Status to return.
!            0    Normal return.
!            1    Allocation failures.
!
INTEGER, INTENT(IN):: parent_i
INTEGER, INTENT(IN):: id
TYPE(BoxMatrix), INTENT(INOUT),TARGET:: b
TYPE(BoxLine), INTENT(INOUT):: setB
TYPE(real_vector),INTENT(INOUT):: setDia
TYPE(int_vector), INTENT(INOUT), TARGET:: setInd
TYPE(int_vector), INTENT(INOUT):: setFcol
TYPE(HyperBox), INTENT(INOUT):: p_box
TYPE(ValList), INTENT(INOUT):: setW
TYPE(int_vector), INTENT(INOUT):: setI
INTEGER, INTENT(OUT):: iflag

! Local variables.
INTEGER:: b_id, b_j, i, j, k, status
INTEGER, DIMENSION(2*D):: sortInd
TYPE(BoxLink), POINTER:: p_l
TYPE(BoxMatrix), POINTER:: p_b
TYPE(int_vector), POINTER:: p_i, p_setInd
REAL(KIND = R8):: temp

! Initialize 'iflag' for a normal return.
iflag = 0
! Find the desired node of 'setInd'.
p_setInd => setInd
DO i = 1, id - 1
  p_setInd => p_setInd%next
END DO

! Clear the CONVEX_BIT of 'flags' as being processed.
p_setInd%flags(parent_i) = IBCLR(p_setInd%flags(parent_i), CONVEX_BIT)
IF (p_setInd%elements(parent_i) <= col_w) THEN
  ! This column is in the head link of box matrices.
  p_b => b
  b_j = p_setInd%elements(parent_i)
ELSE
  ! Find the box matrix that contains this column.
  b_id = (p_setInd%elements(parent_i)-1)/col_w + 1
  b_j = MOD(p_setInd%elements(parent_i)-1, col_w) + 1
  p_b => b
  DO i = 1, b_id-1
    p_b => p_b%child
  END DO
END IF

! Fill out 'setW'.
DO i = 1, setB%ind, 2
  ! Add minimum 'val' of a pair of newly sampled center points
  ! into 'setW'.
  setW%val((i+1)/2) = MIN(setB%Line(i)%val, setB%Line(i+1)%val)
  setW%dir((i+1)/2) = setB%dir(i)
END  DO
setW%dim = setB%ind/2

! Find the order of dimensions for further dividing by insertion
! sorting wi's in 'setW'.
DO i = 2, setW%dim
  DO j = i, 2, -1
    IF (setW%val(j) < setW%val(j-1)) THEN
      ! Element j is smaller than element j-1, so swap them. Also,
      ! the associated directions are swapped.
      temp = setW%val(j)
      k = setW%dir(j)
      setW%val(j) = setW%val(j-1)
      setW%dir(j) = setW%dir(j-1)
      setW%val(j-1) = temp
      setW%dir(j-1) = k
    ELSE
      EXIT
    END IF
  END DO
END DO

! Sort the indices of boxes in setB according to the dividing order in
! 'setW%dir'. Record the sorted indices in 'sortInd'.
DO i = 1, setW%dim
  DO j = 1, setB%ind, 2
    IF (setB%dir(j) == setW%dir(i)) THEN
      sortInd(2*i-1) = j
      sortInd(2*i) = j + 1
    END IF
  END DO
END DO
! 'setW%dir' contains the order of dimensions to divide the parent box.
! Loop dividing on all dimensions in 'setW%dir' by setting up the new
! side lengths as 1/3 of parent box side lengths for each newly
! sampled box center.
DO i = 1, setW%dim
  temp = p_b%M(1,b_j)%side(setW%dir(i))/3.0_R8
  DO j = i, setW%dim
   setB%Line(sortInd(2*j-1))%side(setW%dir(i)) = temp
   setB%Line(sortInd(2*j))%side(setW%dir(i)) = temp
  END DO
  ! Modify the parent's side lengths.
  p_b%M(1,b_j)%side(setW%dir(i))= temp
END DO
! Clear 'setW' for next time.
setW%dim = 0

! Remove the parent box from box matrix 'p_b'.
p_box = p_b%M(1,b_j)
! Move the last box to the first position.
IF (p_b%ind(b_j) <= row_w) THEN
  ! There are no box links.
  p_b%M(1,b_j) = p_b%M(p_b%ind(b_j),b_j)
ELSE
  ! There are box links. Chase to the last box link.
  p_l => p_b%sibling(b_j)%p
  DO i = 1, (p_b%ind(b_j)-1)/row_w-1
    p_l => p_l%next
  END DO
  p_b%M(1, b_j) = p_l%Line(p_l%ind)
  p_l%ind = p_l%ind-1
END IF
p_b%ind(b_j) = p_b%ind(b_j)-1

! Adjust this box column to a heap by calling siftdown.
CALL siftdown(p_b, b_j, 1)

! Update 'setDia', 'setInd' and 'setFcol' if this column is empty.
! Find which node 'setInd' is associated with by checking 'setInd%id'.
IF (p_b%ind(b_j) == 0)THEN
  ! This column is empty. Remove this diameter squared from a
  ! corresponding node of 'setDia'.
  CALL rmNode(col_w, p_setInd%id-1, parent_i, setDia)
  ! Push the released column back to top of 'setFcol'.
  IF (setFcol%dim < col_w) THEN
    ! The head node of 'setFcol' is not full.
    CALL insNode(col_w, p_setInd%elements(parent_i), &
                    setFcol%dim + 1, setFcol)
  ELSE
    ! The head node is full. There must be at least one more node
    ! for 'setFcol'. Find the last non-full node of 'setFcol' to
    ! insert the released column.
    p_i => setFcol%next
    DO
      IF (p_i%dim < col_w) THEN
        ! Found it.
        CALL insNode(col_w, p_setInd%elements(parent_i), &
                        p_i%dim + 1, p_i)
        EXIT
      END IF
      ! Go to the next node.
      p_i=> p_i%next
    END DO
  END IF
  ! Remove the column index from a corresponding node of 'setInd'.
  CALL rmNode(col_w, 0, parent_i, p_setInd)
END IF

! Modify the diameter squared for the parent box temporarily saved in
! 'p_box'. Compute the diameter in the original frame.
p_box%diam = DOT_PRODUCT(p_box%side*UmL, p_box%side*UmL)

! Update 'dia' associated with 'FMIN', which has coordinates in 'unit_x'.
IF (ALL(unit_x == p_box%c)) dia = p_box%diam

! Compute squared diameters for all new boxes in 'setB'.
DO i = 1, setB%ind
  setB%Line(i)%diam = DOT_PRODUCT(setB%Line(i)%side*UmL, &
                                  setB%Line(i)%side*UmL)
  ! Update 'dia' if needed.
  IF (ALL(unit_x == setB%Line(i)%c)) dia = setB%Line(i)%diam
END DO

! Add all new boxes in 'setB' and 'p_box' to 'b' according to different
! squared diameters and different function values.
DO i = 1, setB%ind
  CALL insMat(setB%Line(i), b, setDia, setInd, setFcol, status, 0)
  IF (status /=0) THEN;iflag = status;RETURN;END IF
END DO
CALL insMat(p_box, b, setDia, setInd, setFcol, status, 0)
IF (status /=0) THEN;iflag = status;RETURN;END IF
! Clear 'setB' and 'setI' for calling divide() next time.
setB%ind = 0
setI%dim = 0
RETURN
END SUBROUTINE divide

FUNCTION findcol(i_start, p_start, index, do_it) RESULT(p_iset)
IMPLICIT NONE
! Find the leftmost column (setInd%elements(index)) on the lower convex
! hull, in the plot of (box diameter, function value) points, with
! CONVEX_BIT of 'flags' set in linked list 'setInd', which indicates a
! potentially optimal box to be subdivided.
!
! On input:
! i_start - The index to start searching in node 'p_start'.
! p_start - The pointer to the node at which to start searching.
!
! On output:
! index    - The found index in node 'p_iset' of the linked list 'setInd'.
! do_it    - The returned sign to continue processing or not.
! i_start  - The index at which to resume searching in node 'p_start'
!            next time.
! p_start  - The pointer to the node at which to resume searching next time.
! p_iset   - The returned node, which contains the next box to subdivide.
!
INTEGER, INTENT(INOUT) :: i_start
TYPE(int_vector), POINTER :: p_start
INTEGER, INTENT(OUT) :: index
LOGICAL, INTENT(OUT) :: do_it
TYPE(int_vector), POINTER :: p_iset

! Local variables.
INTEGER :: i, start
TYPE(int_vector), POINTER :: p_set

do_it = .FALSE.
start = i_start
p_set => p_start
! Search from the last to the first element in the link 'p_set'.
DO WHILE(ASSOCIATED(p_set))
  DO i = start, 1, -1
    ! Find the last column with CONVEX_BIT set in 'flags' of 'p_set'.
    IF (BTEST(p_set%flags(i), CONVEX_BIT)) THEN
      ! Clear the CONVEX_BIT of 'flags' as being processed.
      p_set%flags(i) = IBCLR(p_set%flags(i), CONVEX_BIT)
      do_it = .TRUE.
      index = i
      p_iset => p_set
      ! Save them to i_start and p_start for resuming searching
      ! next time.
      i_start = i
      p_start => p_set
      EXIT
    END IF
  END DO
  ! There are no more box column with CONVEX_BIT set in 'flags' in
  ! this node. Go to the previous one (if any).
  IF (.NOT. do_it) THEN
    IF(ASSOCIATED(p_set%prev)) THEN
      p_set => p_set%prev
      ! Reset 'start' to be the last one for all the following iterations
      ! except the first one which resumed from 'i_start'.
      start = p_set%dim
    ELSE
      EXIT
    END IF
  ELSE ! Exit when the convex hull box is found.
    EXIT
  END IF
END DO

RETURN
END FUNCTION findcol

SUBROUTINE findconvex(b, p_fmin, i_fmin, setInd)
IMPLICIT NONE
! In 'setInd', clear CONVEX_BIT of columns if the first boxes on these
! columns are not on convex hull. Bit CONVEX_BIT with value 0 indicates
! the first box on the column is not one of potentially optimal boxes.
! This is determined by comparing slopes. If EPS_I is 0, starting from the
! first column, find the maximum slope from the first box on that column
! to the first boxes on all other columns until reaching the box with
! 'FMIN'. Then, starting from the next column with the first box on
! convex hull, repeat the procedure until no more columns before the
! column with 'FMIN' to check. If EPS_I is greater than 0, the outer loop
! breaks out when the maximum slope is less than the value:
! (val-(FMIN-(|FMIN|+1)EPS_I))/diam.
!
! On input:
! b      - The head link of box matrices.
! p_fmin - Pointer to the node holding the column index of the box with
!          'FMIN'.
! i_fmin - Index of the column in the node 'p_fmin'.
! setInd - A linked list holding column indices of box matrices.
!
! On output:
! setInd - 'setInd' has the modified column indices.
!
TYPE(BoxMatrix), INTENT(IN), TARGET:: b
TYPE(int_vector), POINTER:: p_fmin
INTEGER, INTENT(IN):: i_fmin
TYPE(int_vector), INTENT(INOUT), TARGET:: setInd

! Local variables.
INTEGER:: b_id1, b_id2, col1, col2, i, j, k, target_i
LOGICAL:: stop_fmin
REAL(KIND = R8):: slope, slope_max
TYPE(BoxMatrix), POINTER:: p_b1, p_b2
TYPE(int_vector), POINTER:: p_setInd1, p_setInd2, target_set

! Initialize the first node pointer.
p_setInd1 => setInd

! Initialization for outer loop that processes all columns before
! the column containing 'FMIN' in order to find a convex hull curve.
stop_fmin = .FALSE.
i = 1
k = 1
OUTLOOP: DO WHILE((.NOT. stop_fmin) .AND. ASSOCIATED(p_setInd1))
  ! Initialization for inner loop, which computes the slope from the
  ! first box on the fixed column 'i' to the first boxes on all the
  ! other columns, before reaching the column containing a box with
  ! 'FMIN', to locate the target column with maximum slope. Mark off
  ! any columns in between the fixed first column and the target column.
  NULLIFY(target_set)
  slope_max = -HUGE(slope)
  p_setInd2 => p_setInd1
  ! Fix the first convex hull column as column 'i' in 'p_setInd1'.
  ! The second column used to calculate the slope has index 'k' in
  ! 'p_setInd2'.  'k' is incremented up to the column index corresponding
  ! to 'FMIN'.  Find the box matrix 'p_b1' and the local column index
  ! 'col1'.
  b_id1 =(p_setInd1%elements(i)-1)/col_w + 1
  col1 = MOD(p_setInd1%elements(i)-1, col_w) + 1
  p_b1 => b
  DO j = 1, b_id1 - 1
    p_b1 => p_b1%child
  END DO
  ! Check if the first column has reached the column with 'FMIN'. If
  ! so, break out of the outer loop.
  IF (dia == p_b1%M(1,col1)%diam) EXIT OUTLOOP
  k = i + 1
  INLOOP: DO
    IF (k > p_setInd2%dim) THEN
      ! Move to the next node as k increments beyond the maximum
      ! length for each node of 'setInd'.
      p_setInd2 => p_setInd2%next
      IF (.NOT. ASSOCIATED(p_setInd2)) EXIT INLOOP
      IF (p_setInd2%dim == 0) EXIT INLOOP
      k = 1
    END IF
    ! To compute the slope from the first box on column 'i' of 'p_setInd1'
    ! to the first box on column 'k' of 'p_setInd2', find the local
    ! column index 'col2' and the corresponding box matrix 'p_b2'.
    b_id2 = (p_setInd2%elements(k) - 1)/col_w + 1
    col2 = MOD(p_setInd2%elements(k) - 1, col_w) + 1
    p_b2 => b
    DO j = 1, b_id2 - 1
      p_b2 => p_b2%child
    END DO
    ! Use the slope formula (f1-f2)/(d1-d2), where f1 and f2 are the
    ! function values at the centers of the two boxes with diameters
    ! d1 and d2.
    slope = (p_b1%M(1,col1)%val - p_b2%M(1,col2)%val) / &
            (SQRT(p_b1%M(1,col1)%diam) - SQRT(p_b2%M(1,col2)%diam))
    ! Compare the new slope with the current maximum slope. Keep track
    ! of the target column index and the target node.
    IF (slope > slope_max) THEN
      slope_max = slope
      target_i = k
      target_set => p_setInd2
      ! Check if this target column contains 'FMIN'.
      IF (dia == p_b2%M(1,col2)%diam) stop_fmin = .TRUE.
    END IF
    ! If the target column contains 'FMIN', break out of the inner loop.
    IF (dia == p_b2%M(1,col2)%diam) EXIT INLOOP
    ! Move on to the next column.
    k = k + 1
  END DO INLOOP
  ! Mark off boxes in between.
  IF (ASSOCIATED(target_set)) CALL markoff(i, target_i, p_setInd1, target_set)
  ! Check if EPS_I /= 0. If so, stop if the found 'slope_max' from
  ! the first box is less than the desired accuracy of the solution.
  IF ((EPS_I /= 0) .AND. ASSOCIATED(target_set)) THEN
    IF ((p_b1%M(1,col1)%val - (FMIN - (ABS(FMIN) + 1.0_R8)*EPS_I))/ &
       SQRT(p_b1%M(1,col1)%diam) > slope_max ) THEN
      ! Mark off the first boxes on the columns from the column target_i
      ! to the one with 'FMIN'.
      target_set%flags(target_i) = IBCLR(target_set%flags(target_i), &
                                         CONVEX_BIT)
      CALL markoff(target_i, i_fmin, target_set, p_fmin)
      ! Mark off the first box on the colume with 'FMIN' if it is not marked
      ! off yet as the target_set.
      IF (BTEST(p_fmin%flags(i_fmin), CONVEX_BIT)) p_fmin%flags(i_fmin) = &
          IBCLR(p_fmin%flags(i_fmin), CONVEX_BIT)
      EXIT OUTLOOP
    END IF
  END IF
  ! To start the next pass, the first fixed column jumps to the target column
  ! just found which is the next column on convex hull.
  i = target_i
  p_setInd1 => target_set
END DO OUTLOOP
RETURN
END SUBROUTINE findconvex

SUBROUTINE findsetI(b, col, setI)
IMPLICIT NONE
! Fill out 'setI', holding dimensions with the maximum side length
! of the first box on 'col' in box matrix links 'b'.
!
! On input:
! b   - The head link of box matrices.
! col - The global column index of box matrix links.
!
! On output:
! setI - The set of dimensions with the maximum side length.
!
TYPE(BoxMatrix), INTENT(IN), TARGET:: b
INTEGER, INTENT(IN):: col
TYPE(int_vector), INTENT(INOUT):: setI

! Local variables.
INTEGER:: b_id, i, j
REAL(KIND = R8) :: temp
TYPE(BoxMatrix), POINTER:: p_b

! Find the box matrix link that 'col' is associated with.
IF (col <= col_w) THEN
  p_b => b
  j = col
ELSE
  b_id = (col-1)/col_w + 1
  j = MOD(col-1, col_w) + 1
  p_b => b
  DO i = 1, b_id-1
    p_b => p_b%child
  END DO
END IF
! Search for the maximum side length.
temp = MAXVAL(p_b%M(1,j)%side(:))
! Find all the dimensions with the maximum side length.
DO i = 1, D
  IF ((ABS(p_b%M(1,j)%side(i)-temp)/temp) <= EPS4N) THEN
    ! Add dimension index to 'setI'.
    CALL insNode(D, i, setI%dim + 1, setI)
  END IF
END DO
RETURN
END SUBROUTINE findsetI

SUBROUTINE init(b, status)
IMPLICIT NONE
! Allocate the arrays and initialize the first center point.
! Evaluates the function value at the center point and initializes
! 'FMIN' and 'unit_x'.
!
! On input: None.
!
! On output:
! b      - The first box matrix to initialize.
! status - Status of return.
!          =0  Successful.
!          >0  Allocation or restart related error.
!
TYPE(BoxMatrix), INTENT(OUT), TARGET:: b
INTEGER, INTENT(OUT):: status

! Local variables.
INTEGER:: i, iflag, j

! Normal status.
status = 0

! Allocate arrays.
ALLOCATE(b%M(row_w, col_w), STAT = iflag)
IF (iflag /= 0) THEN;status = 1;RETURN;END IF
ALLOCATE(b%ind(col_w), STAT = iflag)
IF (iflag /= 0) THEN;status = 1;RETURN;END IF
! Clear the box counter for each column.
b%ind(:) = 0
! Nullify the child link to the next box matrix.
NULLIFY(b%child)
ALLOCATE(b%sibling(col_w), STAT = iflag)
IF (iflag /= 0) THEN;status = 1;RETURN;END IF
DO i = 1, col_w
  NULLIFY(b%sibling(i)%p)
END DO
DO i = 1, row_w
  DO j = 1, col_w
    ALLOCATE(b%M(i,j)%c(D), STAT = iflag)
    IF (iflag /= 0) THEN;status = 1;RETURN;END IF
    ALLOCATE(b%M(i,j)%side(D), STAT = iflag)
    IF (iflag /= 0) THEN;status = 1;RETURN;END IF
  END DO
END DO

! Starting from the last column, push free columns to 'setFcol'.
DO i = 1, col_w-1
  setFcol%elements(i) = col_w-(i-1)
END DO
! Use the first column.
setFcol%dim = col_w-1
! Initialize the center of the first unit hypercube in box matrix 'b'
! and 'unit_x' in the normalized coordinate system.
b%M(1,1)%c(:) = 0.5_R8
b%M(1,1)%side(:) = 1.0_R8
unit_x(:) = 0.5_R8
! Evaluate objective function at 'c'.
! Store the function value and initialize 'FMIN'.
iflag = 0
! For usage with VTMOP, this evaluation is performed before calling
! bVTdirect, and so OBJ_FUNC only querries VTMOP's internal database.
FMIN = OBJ_FUNC(L + b%M(1,1)%c(:)*UmL, WEIGHTS, iflag)
! Check 'iflag' to deal with undefined function values.
IF (iflag /= 0) THEN
  ! Assign a huge value to 'FMIN' for an infeasible region.
  FMIN = HUGE(1.0_R8)
END IF
b%M(1,1)%val = FMIN

! Initialize the diameter squared for this box and 'dia', the diameter
! squared associated with 'FMIN'. Convert 'side' back to the orginal frame
! for computing the real 'dia' that will be compared with 'dia_limit' or
! 'MIN_DIA' given by the user.
dia = DOT_PRODUCT(b%M(1,1)%side*UmL, b%M(1,1)%side*UmL)
b%M(1,1)%diam = dia
! Initialize the 'ind' for the first column and 'id' for this box matrix.
b%ind(1) = 1
b%id = 1
! Initialize 'setDia', 'setInd', and 'setFcol'.
setDia%dim = 1
setDia%elements(1) = b%M(1,1)%diam
! Set the first box as being on convex hull by setting the CONVEX_BIT
! of the 'flags'.
setInd%dim = 1
setInd%elements(1) = 1
setInd%flags(1) = IBSET(setInd%flags(1), CONVEX_BIT)
RETURN
END SUBROUTINE init

SUBROUTINE markoff(i, target_i, p_iset1, target_set)
IMPLICIT NONE
! Mark off columns in between the column 'i' of 'p_iset1' and column
! 'target_i' of 'target_set' by clearing the CONVEX_BIT in 'flags'.
!
! On input:
! i          - Column index of the first box for computing slope in findconvex.
! target_i   - Column index of the second box for computing slope in
!              findconvex.
! p_iset1  - The node of 'setInd' holding the column 'i'.
! target_set - The node of 'setInd' holding the column 'target_i'.
!
! On output:
! p_iset1  - 'p_iset1' has changed column indices.
! target_set - 'target_set' has changed column indices.
!
INTEGER, INTENT(IN) :: i
INTEGER, INTENT(IN) :: target_i
TYPE(int_vector), INTENT(INOUT), TARGET :: p_iset1
TYPE(int_vector), POINTER :: target_set

! Local variables.
INTEGER :: j
TYPE(int_vector), POINTER :: p_set

! Check if any columns in between.
IF (ASSOCIATED(target_set, p_iset1)) THEN
  ! If 'target_i' is next to column 'i' or no any columns in between, return.
  IF ((i == target_i) .OR. (i+1 == target_i)) RETURN
END IF

! Clear all CONVEX_BITs in 'flags' in between.
j = i
p_set => p_iset1
DO
  j = j + 1
  IF (j > p_set%dim) THEN
    p_set => p_set%next
    IF (.NOT. ASSOCIATED(p_set)) EXIT
    IF (p_set%dim == 0) EXIT
    j = 1
  END IF
  ! Check if at the target node.
  IF (ASSOCIATED(target_set, p_set)) THEN
    ! If 'j' has reached 'target_i', exit.
    IF (j == target_i) EXIT
  END IF
  ! Clear the CONVEX_BIT of 'flags' for column 'j' of 'p_set'.
  p_set%flags(j) = IBCLR(p_set%flags(j), CONVEX_BIT)
END DO
END SUBROUTINE markoff

SUBROUTINE sampleF(setB, eval_c, ierr)
IMPLICIT NONE
! Evaluate the objective function at each newly sampled center point.
! Keep updating 'FMIN' and 'unit_x'.
!
! On input:
! setB   - The set of newly sampled boxes with their center points'
!          coordinates.
! eval_c - The counter of evaluations.
!
! On output:
! setB   - The set of newly sampled boxes with added function values
!          at center points.
! eval_c - The updated counter of evaluations.
! ierr   - Return status.
!          0  Normal return.
!          >0 Error return.
!
TYPE(BoxLine), INTENT(INOUT):: setB
INTEGER, INTENT(INOUT):: eval_c
INTEGER, INTENT(OUT):: ierr

! Local variables.
INTEGER:: i, iflag

! Initialize 'ierr'.
ierr = 0

! Normal evaluation in the original coordinate system.
! Loop evaluating all the new center points of boxes in 'setB'.

! Execute function evaluations in an out-of-order parallel loop.
!$OMP PARALLEL DO SCHEDULE(DYNAMIC), DEFAULT(SHARED), PRIVATE(i, iflag)
DO i = 1, setB%ind
  iflag = 0
  setB%Line(i)%val = OBJ_FUNC(L + setB%Line(i)%c(:)*UmL, WEIGHTS, iflag)
  ! Check 'iflag'.
  IF (iflag /= 0) THEN
    ! Assign a huge value for an infeasible region.
    setB%Line(i)%val = HUGE(1.0_R8)
  END IF
END DO
!$OMP END PARALLEL DO

! Commit the function evaluations to the database in order.
DO i = 1, setB%ind
  ! Update evaluation counter.
  eval_c = eval_c + 1
  IF (FMIN > setB%Line(i)%val) THEN
    ! Update 'FMIN' and 'unit_x'.
    FMIN = setB%Line(i)%val
    unit_x(:) = setB%Line(i)%c(:)
  END IF
  ! Update FMIN and 'unit_x' in Lexicographical order.
  IF (FMIN == setB%Line(i)%val) THEN
    IF (setB%Line(i)%c .lexLT. unit_x) THEN
      ! The new point is smaller lexicographically than 'unit_x'.
      unit_x(:) = setB%Line(i)%c(:)
    END IF
  END IF
END DO
RETURN
END SUBROUTINE sampleF

SUBROUTINE sampleP(col, setI, b, setB)
IMPLICIT NONE
! In each dimension 'i' in 'setI', samples at the two points c+delta*e_i
! and c-delta*e_i. e_i is the 'i'th unit vector. In 'setB', records all
! the new points as the centers of boxes that will be formed completely
! through subroutines 'sampleF' and 'divide'.
!
! On input:
! col  - The global column index of the box to subdivide.
! setI - The set of dimensions with maximum side length.
! b    - The head link of box matrices.
! setB - The empty set of type 'HyperBox' that will hold newly sampled
!        points as the centers of new boxes.
!
! On output:
! setI - The set of dimensions with maximum side length.
! b    - The head link of box matrices.
! setB - The set of boxes that contains the newly sampled center points.
!
INTEGER, INTENT(IN):: col
TYPE(int_vector), INTENT(INOUT):: setI
TYPE(BoxMatrix), INTENT(INOUT), TARGET:: b
TYPE(BoxLine), INTENT(INOUT):: setB

! Local variables.
INTEGER:: b_id  ! Identifier of the associated box matrix.
INTEGER:: i     ! Loop counters.
INTEGER:: j     ! Local column index converted from the global one 'col'.
INTEGER:: new_i ! Index of new points in setB.
REAL (KIND=R8):: delta ! 1/3 of the maximum side length.
TYPE (BoxMatrix), POINTER:: p_b ! Pointer to the associated box matrix.

! Find the box matrix that 'col' is associated with. Store the pointer
! to box matrix in 'p_b'. The local column index 'j' will be converted from
! 'col'.
IF (col <= col_w) THEN
  p_b => b
  j = col
ELSE
  b_id = (col-1)/col_w + 1
  j = MOD(col-1, col_w) + 1
  p_b => b
  DO i = 1, b_id-1
    p_b => p_b%child
  END DO
END IF

! Find the maximum side length by obtaining the dimension in 'setI'.
! Then, extract the maximum side length from the first box on column 'j'
! of box matrix 'p_b'. Calculate 'delta', 1/3 of the maximum side length.
delta = p_b%M(1,j)%side(setI%elements(setI%dim)) / 3.0_R8

! Loop sampling two new points in all dimensions 'i' in 'setI'.
! c+delta*e_i => newpt_1; c-delta*e_i => newpt_2, where e_i is the 'i'th
! unit vector.
DO i = 1, setI%dim
  new_i = setB%ind + 1
  ! Copy the coordinates of parent box to the two new boxes
  setB%Line(new_i)%c(:) = p_b%M(1, j)%c(:)
  setB%Line(new_i+1)%c(:) = p_b%M(1, j)%c(:)
  ! Assign changed coordinates to the two new points in 'setB'.
  setB%Line(new_i)%c(setI%elements(i)) = &
    p_b%M(1, j)%c(setI%elements(i)) + delta
  setB%Line(new_i+1)%c(setI%elements(i)) = &
  p_b%M(1, j)%c(setI%elements(i))-delta

  ! Record the directions with changes in 'setB%dir' for further
  ! processing to find the dividing order of dimensions.
  setB%dir(new_i) = setI%elements(i)
  setB%dir(new_i+1) = setI%elements(i)
  ! Update 'ind' of 'setB'.
  setB%ind = setB%ind + 2
  ! Initialize side lengths of new points for further dividing
  ! by copying the sides from the parent box.
  setB%Line(new_i)%side(:) = p_b%M(1,j)%side(:)
  setB%Line(new_i+1)%side(:) = p_b%M(1,j)%side(:)
END DO
! Clear 'setI'.
setI%dim = 0
RETURN
END SUBROUTINE sampleP

SUBROUTINE sanitycheck(iflag)
IMPLICIT NONE
! Check the sanity of the input arguments, and set all local variables
! derived from input arguments.
!
! On input: None.
!
! On output:
! iflag  - The sanity check result.
!
INTEGER, INTENT(INOUT):: iflag

! Initialize 'iflag'.
iflag = 0
!
! Check required arguments.
!
IF (D < 2) THEN
  iflag = INPUT_ERROR
  RETURN
ELSE
  N_I = D
END IF

IF ((SIZE(X) /= D) .OR. (SIZE(L) /= D) .OR. (SIZE(U) /= D)) THEN
  iflag = INPUT_ERROR + 1
  RETURN
END IF
IF ((SIZE(WEIGHTS) /= P)) THEN ! Check WEIGHTS(:) against P.
  iflag = INPUT_ERROR + 1
  RETURN
END IF
IF (ANY(L >= U)) THEN
  iflag = INPUT_ERROR + 2
  RETURN
END IF
IF (ANY(WEIGHTS < 0.0_R8)) THEN ! Check that WEIGHTS(:) is nonnegative.
  iflag = INPUT_ERROR + 2
  RETURN
END IF
!
! Check optional arguments.
!
IF (PRESENT(W)) THEN
  IF (SIZE(W) /= D) THEN
    iflag = INPUT_ERROR + 1
    RETURN
  END IF
END IF
! Default: processing boxes only on convex hull.
SWITCH_I = 1
IF (PRESENT(SWITCH)) THEN
  IF ((SWITCH < 0) .OR. (SWITCH > 1)) THEN
    iflag = INPUT_ERROR + 5
    RETURN
  END IF
  IF (SWITCH == 0) THEN
    IF (PRESENT(EPS)) THEN
      IF (EPS > 0.0_R8) THEN
        iflag = INPUT_ERROR + 6
        RETURN
      END IF
    END IF
  END IF
  ! Assign the local copy 'SWITCH_I'.
  SWITCH_I = SWITCH
END IF
! Enable LBC (limiting box columns) by default.
lbc = 1
stop_rule = 0
! When MAX_ITER <= 0, the number of iterations will be returned on exit.
IF (PRESENT(MAX_ITER)) THEN
  IF (MAX_ITER > 0) THEN
    ! Set bit 0 of stop_rule.
    stop_rule = IBSET(stop_rule, STOP_RULE1)
  ELSE ! Disable LBC (limiting box columns).
    lbc = 0
  END IF
ELSE ! Disable LBC (limiting box columns).
  lbc = 0
END IF
! When MAX_EVL <=0, the number of evaluations will be returned on exit.
IF (PRESENT(MAX_EVL)) THEN
  IF (MAX_EVL > 0) THEN
    ! Set bit 1 of stop_rule.
    stop_rule = IBSET(stop_rule, STOP_RULE2)
    ! When 'MAX_ITER' is positive and 'MAX_EVL' is
    ! sufficiently small, disable LBC (limiting box columns).
    IF (PRESENT(MAX_ITER)) THEN
      IF ((MAX_ITER > 0) .AND. (MAX_EVL*(2*D+2) < 2D+6)) lbc = 0
    END IF
  END IF
END IF

! Even if user doesn't specify 'MIN_DIA', a diameter smaller than
! SQRT(SUM(UmL*UmL))*EPSILON(1.0_R8) is not permitted to occur.
! When MIN_DIA <=0, the diameter associated with X and FMIN will be
! returned on exit. Assign 'UmL'.
UmL(:) = U(:) - L(:)
dia_limit =  SQRT(SUM(UmL*UmL))*EPSILON(1.0_R8)
IF (PRESENT(MIN_DIA)) THEN
  IF (MIN_DIA > 0.0_R8) THEN
    IF (MIN_DIA < dia_limit) THEN
      iflag = INPUT_ERROR + 3
      RETURN
    ELSE
      dia_limit = MIN_DIA
      ! Set bit 2 of stop_rule.
      stop_rule = IBSET(stop_rule, STOP_RULE3)
    END IF
  END IF
END IF

! When OBJ_CONV is present a minimum relative change in the minimum
! objective function value will be enforced.
IF (PRESENT(OBJ_CONV)) THEN
  IF (OBJ_CONV /= 0.0_R8) THEN
    IF ((OBJ_CONV < EPSILON(1.0_R8)*REAL(D,KIND=R8))&
       .OR.(OBJ_CONV >= 1.0_R8)) THEN
      iflag = INPUT_ERROR + 3
      RETURN
    ELSE
      ! Set bit 3 of stop_rule.
      stop_rule = IBSET(stop_rule, STOP_RULE4)
    END IF
  END IF
END IF

! When EPS is present a test involving EPS is used to define potentially
! optimal boxes.  The absence of this test is equivalent to EPS=0.
EPS_I = 0.0_R8
IF (PRESENT(EPS)) THEN
  IF (EPS /= 0.0_R8) THEN
    IF ((EPS < EPSILON(1.0_R8)) .OR. (EPS > 1.0_R8)) THEN
      iflag = INPUT_ERROR + 3
      RETURN
    ELSE
      EPS_I = EPS
    END IF
  END IF
END IF

! Check if stop_rule has at least at 1 bit set. Otherwise no stopping rule
! has been given.
IF (stop_rule == 0) THEN;iflag = INPUT_ERROR+4;RETURN;END IF

! Check if MIN_SEP, W, and BOX_SET are correctly set.
IF (PRESENT(BOX_SET)) THEN
  ! Set default weights for distance definition.
  W_I(1:D) = 1.0_R8
  ! Verify and set weights for distance definition.
  IF (PRESENT(W)) THEN
    WHERE (W > 0)
      W_I = W
    ELSEWHERE
      W = 1.0_R8
    END WHERE
  END IF
  ! Compute the weighted diameter of the original design space. Reuse
  ! the variable 'dia' (the diameter squared associated with 'FMIN').
  dia = SQRT(SUM(UmL*W_I*UmL))
  ! Check the optional argument MIN_SEP. Set a default value if MIN_SEP
  !  is not present or not correctly assigned.
  MIN_SEP_I = (0.5_R8*dia)**2
  IF (PRESENT(MIN_SEP)) THEN
    IF ((MIN_SEP < dia*EPSILON(1.0_R8)) .OR. (MIN_SEP > dia)) &
      MIN_SEP = 0.5_R8*dia
    MIN_SEP_I = MIN_SEP**2
  END IF
  ! Check the size of its component arrays c(:) and side(:).
  DO i = 1, SIZE(BOX_SET)
    IF (ASSOCIATED(BOX_SET(i)%c)) THEN
      IF (SIZE(BOX_SET(i)%c(:)) /= D) THEN
        iflag = ALLOC_ERROR + 4
        RETURN
      END IF
    ELSE ! Allocate component 'c'.
      ALLOCATE(BOX_SET(i)%c(D))
    END IF
    IF (ASSOCIATED(BOX_SET(i)%side)) THEN
      IF (SIZE(BOX_SET(i)%side) /= D) THEN
        iflag = ALLOC_ERROR + 4
        RETURN
      END IF
    ELSE ! Allocate component 'side'.
      ALLOCATE(BOX_SET(i)%side(D))
    END IF
  END DO
  ! Initialize the index counter 'boxset_ind'.
  boxset_ind = 0
  ! Disable LBC (limiting box columns) that removes boxes needed for
  ! finding 'BOX_SET'.
  lbc = 0
END IF

RETURN
END SUBROUTINE sanitycheck

SUBROUTINE squeeze()
IMPLICIT NONE
! Scan through box columns to remove empty box columns and squeeze box
! column lengths to MAX_ITER-t+1 if LBC is enabled.
!
! On input: None.
!
! On output: None.
!
! Local variables.
INTEGER:: b_id ! Box matrix index.
INTEGER:: i, i1, j, j1 ! Loop counters.
INTEGER:: leaf ! Index to a leaf node in a heap(a box column).
INTEGER:: maxid ! Global index to a box with the largest function value.
INTEGER:: maxpos ! Local index of the box (ranked 'maxid' globally in the
  ! heap) in the box matrix or box link.
INTEGER:: mycol ! Local box column index.
TYPE(BoxLink), POINTER:: p_l, p_l1 ! Pointer to a box link.
TYPE(BoxMatrix), POINTER:: p_b ! Pointer to a box matrix.
TYPE(int_vector), POINTER:: p_seti ! Pointer to a link of 'int_vector'.
REAL(KIND = R8):: maxf ! The largest function value in a heap.
REAL(KIND = R8), DIMENSION(D):: maxc ! Coordinates of the box with 'maxf'.

! The scan starts from the first box column indexed in 'setInd'.
p_seti => setInd
DO WHILE(ASSOCIATED(p_seti))
  i = 1
  DO
    ! Exit when reaching the last box column index in the setInd link
    ! pointed to by 'p_seti'.
    IF (i > p_seti%dim) EXIT
    IF (p_seti%elements(i) <= col_w) THEN
      ! This box column is in the head link of the linked list of box
      ! matrices. Let 'p_b' point to the box matrix and find the local
      ! column index.
      p_b => m_head
      mycol = p_seti%elements(i)
    ELSE
      ! Find the box matrix link that contains this box column.
      b_id = (p_seti%elements(i)-1)/col_w + 1
      mycol = MOD(p_seti%elements(i)-1, col_w) + 1
      p_b => m_head
      DO j = 1, b_id-1
        p_b => p_b%child
      END DO
    END IF
    IF (p_b%ind(mycol) == 0)  THEN
      ! This column is empty. Remove this diameter squared from the
      ! corresponding node of 'setDia'.
      CALL rmNode(col_w, p_seti%id-1, i, setDia)
      ! Push the released box column back to top of 'setFcol'.
      IF (setFcol%dim < col_w) THEN
        ! The head node of 'setFcol' is not full.
        CALL insNode(col_w, p_seti%elements(i), setFcol%dim + 1, setFcol)
      ELSE
        ! The head node is full. There must be at least one more node
        ! for 'setFcol'. Find the last non-full node of 'setFcol' to
        ! insert the released column.
        p_seti => setFcol%next
        DO
          IF (p_seti%dim < col_w) THEN
            ! Found it.
            CALL insNode(col_w, p_seti%elements(i), p_seti%dim + 1, p_seti)
            EXIT
          END IF
          ! Go to the next node.
          p_seti=> p_seti%next
        END DO
      END IF
      ! Remove the column index from a corresponding node of 'setInd'.
      CALL rmNode(col_w, 0, i, p_seti)
      ! Decrease 'i' by 1 to compensate for this removal.
      i = i - 1
    ELSE ! This box column is not empty.
      IF (lbc == 1) THEN ! Adjust the box column length to
                         ! MAX(MAX_ITER - t + 1, row_w).
        IF (p_b%ind(mycol) > MAX(MAX_ITER - t + 1, row_w)) THEN
          ! Loop removing the largest element until the goal is reached.
          DO
            IF (p_b%ind(mycol) == MAX(MAX_ITER - t + 1,row_w)) EXIT
            ! Look for the largest element starting from the first leaf.
            ! Initialize and keep update variables ('maxid', 'maxpos',
            ! 'maxf', and 'maxc').
            leaf = p_b%ind(mycol)/2 + 1
            IF (leaf <= row_w) THEN
              ! The first leaf node starts inside 'M', part of the box
              ! matrix.
              maxid = 0
              maxpos = leaf
              maxf = p_b%M(leaf,mycol)%val
              maxc(:) = p_b%M(leaf,mycol)%c(:)
            ELSE ! The first leaf node starts from a box link.
              maxid = (leaf - 1)/row_w
              maxpos = MOD(leaf - 1, row_w) + 1
              p_l => p_b%sibling(mycol)%p
              DO k = 1, maxid - 1
                p_l => p_l%next
              END DO
              maxf = p_l%Line(maxpos)%val
              maxc(:) = p_l%Line(maxpos)%c(:)
            END IF
            IF (maxid == 0) THEN ! Search starts from 'M'.
              DO k = maxpos + 1, row_w ! Compare with the boxes positioned
                ! later than the current 'maxpos' in 'M'.
                ! Exit when reaching the end.
                IF (k > p_b%ind(mycol)) EXIT
                ! Update the variables with the newly found box that has
                ! larger function value or larger lexicographical order.
                IF((p_b%M(k,mycol)%val > maxf) .OR.  &
                  ((p_b%M(k,mycol)%val == maxf) .AND.  &
                  (maxc .lexLT. p_b%M(k,mycol)%c))) THEN
                  maxf = p_b%M(k,mycol)%val
                  maxc(:) = p_b%M(k,mycol)%c(:)
                  maxpos = k
                END IF
              END DO
              ! Search moves on to the box links (if any) and reset
              ! counters.
              p_l => p_b%sibling(mycol)%p
              k = 0
              j = 0
              i1 = 0
            ELSE
              ! Start search from a box link and initialize the counters.
              k = maxid - 1
              j = maxpos
              i1 = maxid - 1
            END IF
            ! Search through all the box links after the current box with
            ! the largest value, which is stored in the box link indexed
            ! by 'i1' (if 0, the box is in 'M').
            DO j1 = 1, (p_b%ind(mycol) - 1)/row_w - i1
              ! Exit if there are no more box links.
              IF (.NOT. ASSOCIATED(p_l)) EXIT
              IF (p_l%ind == 0) THEN
                ! Deallocate an empty box link and exit.
                DEALLOCATE(p_l)
                EXIT
              END IF
              k = k + 1 ! Increase the global index of a box link.
              DO
                j = j + 1 ! Increase the local counter inside a box link.
                IF (j > p_l%ind) THEN
                  ! Reset 'j' at the end of the box link.
                  j = 0
                  EXIT
                END IF
                ! Update the varaibles with the newly found box that has
                ! larger function valu or larger lexicographical order.
                IF ((p_l%Line(j)%val > maxf) .OR.  &
                   ((p_l%Line(j)%val == maxf) .AND.  &
                   (maxc .lexLT. p_l%Line(j)%c))) THEN
                  maxf = p_l%Line(j)%val
                  maxc = p_l%Line(j)%c
                  maxid = k
                  maxpos = j
                END IF
              END DO
              p_l => p_l%next ! Move to the next link.
            END DO
            ! Found the largest fval box at link maxid and position maxpos.
            ! If maxpos is not the last element position in the heap,
            ! replace position maxpos with the last element of the heap
            ! and sift it up.
            IF (maxpos + maxid*row_w < p_b%ind(mycol)) THEN
              ! Locate the last heap element.
              IF (p_b%ind(mycol) <= row_w) THEN
                ! The last element is inside M.
                p_b%M(maxpos,mycol) = p_b%M(p_b%ind(mycol),mycol)
                ! Update the box counter for the box column.
                p_b%ind(mycol) = p_b%ind(mycol) - 1
                CALL siftup (p_b, mycol, maxpos)
              ELSE ! The last element is inside a box link.
                IF (maxid == 0) THEN
                  ! The largest element is inside M.
                  p_l1 => p_b%sibling(mycol)%p
                  DO k = 1, (p_b%ind(mycol) - 1)/row_w - 1
                    p_l1 => p_l1%next
                  END DO
                  p_b%M(maxpos,mycol) = &
                    p_l1%Line(MOD(p_b%ind(mycol) - 1, row_w) + 1)
                ELSE ! The largest element is inside a box link.
                  p_l => p_b%sibling(mycol)%p
                  DO k = 1, maxid - 1
                    p_l => p_l%next
                  END DO
                  p_l1 => p_l
                  DO k = 1, (p_b%ind(mycol) - 1)/row_w - maxid
                    p_l1 => p_l1%next
                  END DO
                  p_l%Line(maxpos) = &
                    p_l1%Line(MOD(p_b%ind(mycol) - 1,row_w) + 1)
                END IF
                ! Update counters for the box link and the box column.
                p_l1%ind = p_l1%ind - 1
                p_b%ind(mycol) = p_b%ind(mycol) - 1
                ! Siftup the box that replaced the largest fval box.
                CALL siftup(p_b, mycol, maxid*row_w + maxpos)
              END IF
            ELSE ! The largest element is the last element.
              IF (maxid > 0) THEN
                ! Locate the box link that holds the largest element.
                p_l => p_b%sibling(mycol)%p
                DO k = 1, (p_b%ind(mycol) - 1)/row_w - 1
                  p_l => p_l%next
                END DO
                p_l%ind = p_l%ind - 1
              END IF
              ! Update the counter.
              p_b%ind(mycol) = p_b%ind(mycol) - 1
            END IF
          END DO
        END IF
      END IF
    END IF
    ! Go to the next box column.
    i = i + 1
  END DO
  ! Go to the next 'setInd' link.
  p_seti => p_seti%next
END DO
RETURN
END SUBROUTINE squeeze

END SUBROUTINE bVTdirect

END MODULE bVTdirect_MOD
