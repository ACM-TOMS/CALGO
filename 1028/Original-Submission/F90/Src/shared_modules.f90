! This file (shared_modules.f95) contains modules REAL_PRECISION,
! VTDIRECT_GLOBAL that defines data types, parameters, and module
! procedures used by both VTdirect and pVTdirect, VTDIRECT_COMMSUB
! that declares the subroutines and functions shared by VTdirect
! and pVTdirect, and VTDIRECT_CHKPT that defines data types and module
! procedules for the checkpointing feature used in both VTdirect and
! pVTdirect.

MODULE REAL_PRECISION  ! From HOMPACK90.
  ! This is for 64-bit arithmetic.
  INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

MODULE VTDIRECT_GLOBAL
USE REAL_PRECISION
IMPLICIT NONE
!
!HyperBox: Defines an n-dimensional box.
!  val  - Function value at the box center.
!  c    - The center point coordinates.
!  side - Box side lengths for all dimensions.
!  diam - Box diameter squared.
!
TYPE HyperBox
  REAL(KIND = R8) :: val
  REAL(KIND = R8), DIMENSION(:), POINTER :: c
  REAL(KIND = R8), DIMENSION(:), POINTER :: side
  REAL(KIND = R8) :: diam
END TYPE HyperBox
!
!BoxLink: Contains 1-D array of hyperboxes, linked to each column of
!         BoxMatrix when needed.
!  Line - 1-D array of boxes.
!  ind  - Index of last box in array 'Line'.
!  next - The pointer to the next BoxLink.
!  prev - The pointer to the previous BoxLink.
!
TYPE BoxLink
  TYPE(HyperBox), DIMENSION(:), POINTER :: Line
  INTEGER :: ind
  TYPE(BoxLink), POINTER :: next
  TYPE(BoxLink), POINTER :: prev
END TYPE BoxLink
!
!P_BoxLink: Contains a pointer to a BoxLink for a column in BoxMatrix.
!  p - Pointer to a BoxLink.
!
TYPE P_BoxLink
  TYPE(BoxLink), POINTER :: p
END TYPE P_BoxLink
!
!BoxLine: Contains 1-D array of newly sampled hyperboxes.
!  Line - 1-D array of boxes.
!  ind  - Index of last box in array 'Line'.
!  dir  - Directions in which box centers are sampled.
!
TYPE BoxLine
  TYPE(HyperBox), DIMENSION(:), POINTER :: Line
  INTEGER :: ind
  INTEGER, DIMENSION(:), POINTER :: dir
END TYPE BoxLine
!
!BoxMatrix: Contains 2-D array of hyperboxes.
!  M       - 2-D array of boxes.
!  ind     - An array holding the number of boxes in all the columns in 'M'.
!  child   - The pointer to the next BoxMatrix with the smaller diameters.
!  sibling - The pointer array for all columns, pointing to the next
!            BoxLinks with the same diameters.
!  id      - Identifier of this box matrix among all box matrices.
!
TYPE BoxMatrix
  TYPE(HyperBox), DIMENSION(:,:), POINTER :: M
  INTEGER, DIMENSION(:), POINTER :: ind
  TYPE(BoxMatrix), POINTER :: child
  TYPE(P_BoxLink), DIMENSION(:), POINTER :: sibling
  INTEGER :: id
END TYPE BoxMatrix
!
! P_Box: Contains a pointer to a hyperbox.
!
TYPE P_Box
  TYPE (HyperBox), POINTER :: p_box
END TYPE P_Box
!
!int_vector : A list holding integer values.
!  dim      - The index of the last element in the list.
!  elements - The integer array.
!  flags    - The integer array holding the status for 'elements'.
!             Bit 0: status of convex hull processing.
!  next     - The pointer to the next integer vector list.
!  prev     - The pointer to the previous integer vector list.
!  id       - Identifier of this integer vector list among all lists.
!
TYPE int_vector
  INTEGER :: dim
  INTEGER, DIMENSION(:), POINTER :: elements
  INTEGER, DIMENSION(:), POINTER :: flags
  TYPE(int_vector), POINTER :: next
  TYPE(int_vector), POINTER :: prev
  INTEGER :: id
END TYPE int_vector
!
!real_vector: A list holding real values.
!  dim      - The index of the last element in the list.
!  elements - The real array.
!  next     - The pointer to the next real vector list.
!  prev     - The pointer to the previous real vector list.
!  id       - Identifier of this real vector list among all lists.
!
TYPE real_vector
  INTEGER :: dim
  REAL(KIND = R8), DIMENSION(:), POINTER :: elements
  TYPE(real_vector), POINTER :: next
  TYPE(real_vector), POINTER :: prev
  INTEGER :: id
END TYPE real_vector
!
!ValList: a list for sorting the wi for all dimensions i corresponding
!  to the maximum side length. wi = min {f(c+delta*ei), f(c-delta*ei)},
!  the minimum of objective function values at the center point c +
!  delta*ei and the center point c - delta*ei, where delta is one-third of
!  this maximum side length and ei is the ith standard basis vector.
!  dim - The index of the last element in the list.
!  val - An array holding the minimum function values.
!  dir - An array holding the sampling directions corresponding to the
!        function values in array 'val'.
!
TYPE ValList
  INTEGER :: dim
  REAL(KIND = R8), DIMENSION(:), POINTER :: val
  INTEGER, DIMENSION(:), POINTER :: dir
END TYPE

! Parameters.

! Argument input error.
INTEGER, PARAMETER :: INPUT_ERROR = 10
! Allocation failure error.
INTEGER, PARAMETER :: ALLOC_ERROR = 20
! File I/O error.
INTEGER, PARAMETER :: FILE_ERROR = 30
! Stop rule 1: maximum iterations.
INTEGER, PARAMETER :: STOP_RULE1 = 0
! Stop rule 2: maximum evaluations.
INTEGER, PARAMETER :: STOP_RULE2 = 1
! Stop rule 3: minimum diameter.
INTEGER, PARAMETER :: STOP_RULE3 = 2
! Stop rule 4: minimum relative change in objective function.
INTEGER, PARAMETER :: STOP_RULE4 = 3

! 'flags' bits for 'setInd' of type 'int_vector'.
INTEGER, PARAMETER :: CONVEX_BIT = 0
! Bit 0: if set, the first box on the corresponding column is in the convex
!        hull box set.

! Global variables.

INTEGER :: col_w, row_w ! Factors defining reasonable memory space for
  ! each box matrix allocation.
INTEGER :: N_I          ! Local copy of 'N'.
REAL(KIND = R8) :: EPS4N
TYPE(Hyperbox), POINTER :: tempbox ! Box for swapping heap elements.

!$OMP THREADPRIVATE(col_w, row_w, N_I, EPS4N, tempbox)

! Interfaces.
INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE AssgBox
END INTERFACE

INTERFACE insNode
  MODULE PROCEDURE insNodeI
  MODULE PROCEDURE insNodeR
END INTERFACE insNode

INTERFACE rmNode
  MODULE PROCEDURE rmNodeI
  MODULE PROCEDURE rmNodeR
END INTERFACE rmNode

CONTAINS

SUBROUTINE AssgBox(x, y)
IMPLICIT NONE
! Copy the contents of box 'y' to box 'x'.
!
! On input:
! y - Box with type 'HyperBox'.
!
! On output:
! x - Box with type 'HyperBox' having contents of box 'y'.
!
TYPE(HyperBox), INTENT(IN) :: y
TYPE(HyperBox), INTENT(INOUT) :: x
x%val = y%val
x%diam = y%diam
x%c = y%c
x%side = y%side
RETURN
END SUBROUTINE AssgBox

SUBROUTINE insNodeR(n, pt, index, Set)
IMPLICIT NONE
! Insert a real number 'pt' at the indexed position 'index' of 'Set'.
!
! On input:
! n     - The maximum length of the real array in each node of 'Set'.
! pt    - The real number to be inserted to 'Set'.
! index - The position at which to insert 'pt' in a node of 'Set'.
! Set   - A linked list of type(real_vector) nodes.
!
! On output:
! Set   - 'Set' has an added real number and modified 'dim' component.
!
INTEGER, INTENT(IN) :: n
REAL(KIND = R8), INTENT(IN) :: pt
INTEGER, INTENT(IN) :: index
TYPE(real_vector), INTENT(INOUT), TARGET :: Set

! Local variables.
TYPE(real_vector), POINTER :: p_set ! Pointer to a node of 'Set'.

! Insert 'pt' into 'Set'.
IF (Set%dim < n ) THEN
  ! The head node is not full. There are no other nodes.
  ! Update 'dim'.
  Set%dim = Set%dim + 1
  IF (index == Set%dim) THEN
    ! The desired position is at end, so insert 'pt' at end.
    Set%elements(Set%dim) = pt
  ELSE
    ! The desired position is not at end, so shift elements before
    ! insertion.
    Set%elements(index+1:Set%dim) = Set%elements(index:Set%dim-1)
    Set%elements(index) = pt
  END IF
ELSE
  ! The head node is full. Check other nodes.
  p_set => Set%next
  ! To shift elements, find the last node which is not full.
  DO WHILE(p_set%dim == n)
    p_set => p_set%next
  END DO
  ! Found the last node 'p_set' which is not full.
  ! Update 'dim' of 'p_set'. Shift element(s) inside this node, if any.
  p_set%dim = p_set%dim + 1
  ! Loop shifting until reaching the node 'Set' at which to insert 'pt'.
  DO WHILE(.NOT. ASSOCIATED(p_set, Set))
    ! Shift element(s) inside this node,  if any.
    p_set%elements(2:p_set%dim) = p_set%elements(1:p_set%dim-1)
    ! Shift the last element from the previous node to this one.
    p_set%elements(1) = p_set%prev%elements(n)
    ! Finished shifting this node. Go to the previous node.
    p_set => p_set%prev
  END DO
  ! Reached the original node 'Set'.
  IF (index == Set%dim) THEN
    ! The desired position is at end, so insert 'pt' at end.
    Set%elements(Set%dim) = pt
  ELSE
    ! The desired position is not at end, so shift elements before
    ! insertion.
    Set%elements(index+1:Set%dim) = Set%elements(index:Set%dim-1)
    Set%elements(index) = pt
  END IF
END IF
RETURN
END SUBROUTINE insNodeR

SUBROUTINE insNodeI(n, pt, index, Set)
IMPLICIT NONE
! Insert an integer number 'pt' at the indexed position 'index' of 'Set'.
!
! On input:
! n        - The maximum length of the integer array in each node of 'Set'.
! pt       - The integer number to be inserted to 'Set'.
! index    - The position at which to insert 'pt'.
! Set      - A linked list of type(int_vector) nodes.
!
! On output:
! Set   - 'Set' has an added integer number and modified 'dim' component.

INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: pt
INTEGER, INTENT(IN) :: index
TYPE(int_vector), INTENT(INOUT), TARGET :: Set

! Local variables.
TYPE(int_vector), POINTER :: p_set ! Pointer to a node of 'Set'.
LOGICAL :: setflags ! The flag to operate on component 'flags'.

! Assign 'setflags'.
setflags = ASSOCIATED(Set%flags)
! Insert 'pt' into 'Set'.
IF (Set%dim < n ) THEN
  ! The head node is not full. There are no other nodes.
  Set%dim = Set%dim + 1
  IF (index == Set%dim) THEN
    ! The desired position is at end, so insert 'pt' at end.
    Set%elements(Set%dim) = pt
    ! Clear the 'flags'.
    IF (setflags) Set%flags(Set%dim) = 0
ELSE
    ! The desired position is not at end, so shift elements before
    ! insertion.
    Set%elements(index+1:Set%dim) = Set%elements(index:Set%dim-1)
    IF (setflags) THEN
      ! Shift the 'flags'.
      Set%flags(index+1:Set%dim) = Set%flags(index:Set%dim-1)
      ! Clear the 'flags'.
      Set%flags(index) = 0
    END IF
    ! Insert 'pt'.
    Set%elements(index) = pt
  END IF
ELSE
  ! The head node is full. There must be other nodes.
  p_set => Set%next
  ! To shift elements, find the last node which is not full.
  DO WHILE(p_set%dim == n)
    p_set => p_set%next
  END DO
  ! Found the last node 'p_set' which is not full.
  ! Update 'dim' of 'p_set'. Shift element(s), if any.
  p_set%dim = p_set%dim + 1
  ! Loop shifting until reaching the original node 'Set'.
  DO WHILE(.NOT. ASSOCIATED(p_set, Set))
    ! Shift element(s) inside this node, if any.
    p_set%elements(2:p_set%dim) = p_set%elements(1:p_set%dim-1)
    ! Shift the last element from the previous node to this one.
    p_set%elements(1) = p_set%prev%elements(n)
    IF (setflags) THEN
      ! Shift the 'flags'.
      p_set%flags(2:p_set%dim) = p_set%flags(1:p_set%dim-1)
      p_set%flags(1) = p_set%prev%flags(n)
    END IF
    ! Finished shifting this node. Go to the previous node.
    p_set => p_set%prev
  END DO
  ! Reached the original node 'Set'.
  IF (index == Set%dim) THEN
    ! The desired position is at end, so insert 'pt' at end.
    Set%elements(Set%dim) = pt
    ! Clear the 'flags'.
    IF (setflags) Set%flags(Set%dim) = 0
  ELSE
    ! The desired position is not at end, so shift elements before
    ! insertion.
    Set%elements(index+1:Set%dim) = Set%elements(index:Set%dim-1)
    IF (setflags) THEN
      ! Shift the 'flags'.
      Set%flags(index+1:Set%dim) = Set%flags(index:Set%dim-1)
      ! Clear the 'flags'.
      Set%flags(index) = 0
    END IF
    ! Insert 'pt'.
    Set%elements(index) = pt
  END IF
END IF
RETURN
END SUBROUTINE insNodeI

SUBROUTINE rmNodeI(n, offset, index, Set)
IMPLICIT NONE
! Remove an integer entry at position 'index' from the integer array
! in the node at 'offset' links from the beginning of the linked list
! 'Set'.
!
! On input:
! n        - The maximum length of the integer array in each node of 'Set'.
! offset   - The offset of the desired node from the first node of 'Set'.
! index    - The position at which to delete an integer from the integer
!            array in the node.
! Set      - A linked list of type(int_vector) nodes.
!
! On output:
! Set    - The desired node of 'Set' has the indexed integer entry removed
!          and the 'dim' component modified.
!
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: offset
INTEGER, INTENT(IN) :: index
TYPE(int_vector), INTENT(INOUT), TARGET :: Set

! Local variables.
INTEGER :: i  ! Loop counter.
TYPE(int_vector), POINTER :: p_set  ! Pointer to a node of 'Set'.
LOGICAL :: setflags ! The flag to operate on component 'flags'.

! Assign 'setflags'.
setflags = ASSOCIATED(Set%flags)
! Find the desired node.
p_set => Set
DO i = 1, offset
  p_set => p_set%next
END DO
IF (index < p_set%dim) THEN
  ! It's not the last entry in 'p_set', so shift elements.
  p_set%elements(index:p_set%dim-1) = p_set%elements(index+1:p_set%dim)
  ! Shift the 'flags'.
  IF (setflags) p_set%flags(index:p_set%dim-1) = p_set%flags(index+1:p_set%dim)
END IF
IF (p_set%dim < n) THEN
  ! There are not other elements in next node, so remove the indexed
  ! entry directly from 'Set' by updating 'dim'.
  p_set%dim = p_set%dim - 1
ELSE
  ! There might be nodes in which to shift elements.
  ! Check if any element(s) in next node to shift.
  IF (ASSOCIATED(p_set%next)) THEN
    p_set => p_set%next
    DO
      IF (p_set%dim > 0) THEN
        ! There are elements to shift.
        ! Shift one element from p_next into its previous node.
        p_set%prev%elements(n) = p_set%elements(1)
        ! Shift elements inside p_next.
        p_set%elements(1:p_set%dim-1) = p_set%elements(2:p_set%dim)
        IF (setflags) THEN
          ! Shift the 'flags'.
          p_set%prev%flags(n) = p_set%flags(1)
          p_set%flags(1:p_set%dim-1) = p_set%flags(2:p_set%dim)
        END IF
      ELSE
        ! There are no elements to shift. Update 'dim' of previous node.
        p_set%prev%dim = p_set%prev%dim - 1
        EXIT
      END IF
      ! Move on to the next node, if any. If there are no more nodes, update
      ! 'dim'.
      IF (ASSOCIATED(p_set%next)) THEN
        p_set => p_set%next
      ELSE
        p_set%dim = p_set%dim - 1
        EXIT
      END IF
    END DO
  ELSE
    ! There are no more nodes. Update 'dim' of 'p_set'.
    p_set%dim = p_set%dim - 1
  END IF
END IF
RETURN
END SUBROUTINE rmNodeI

SUBROUTINE rmNodeR(n, offset, index, Set)
IMPLICIT NONE
! Remove a real entry at position 'index' from the real array
! in the node at 'offset' links from the beginning of the linked list
! 'Set'.
!
! On input:
! n      - The maximum length of the real array in each node of 'Set'.
! offset - The offset of the desired node from the first node of 'Set'.
! index  - The position at which to delete a real entry from the real
!          array in the node.
! Set    - A linked list of type(real_vector) nodes.
!
! On output:
! Set    - The desired node of 'Set' has the indexed real entry removed
!          and the 'dim' component modified.
!
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: offset
INTEGER, INTENT(IN) :: index
TYPE(real_vector), INTENT(INOUT), TARGET :: Set

! Local variables.
INTEGER :: i  ! Loop counter.
TYPE(real_vector), POINTER :: p_set  ! Pointer to a node of 'Set'.

! Find the desired node.
p_set => Set
DO i = 1, offset
  p_set => p_set%next
END DO
IF (index < p_set%dim)THEN
  ! It's not the last entry in 'p_set', so shift elements.
  p_set%elements(index:p_set%dim-1) = p_set%elements(index+1:p_set%dim)
END IF
IF (p_set%dim < n) THEN
  ! There are not other elements in next node, so remove the indexed
  ! entry directly from 'Set' by updating 'dim'.
  p_set%dim = p_set%dim - 1
ELSE
  ! There might be nodes in which to shift elements.
  ! Check if any element(s) in next node to shift.
  IF (ASSOCIATED(p_set%next)) THEN
    p_set => p_set%next
    DO
      IF (p_set%dim > 0) THEN
        ! There are elements to shift.
        ! Shift one element from p_next into its previous node.
        p_set%prev%elements(n) = p_set%elements(1)
        ! Shift elements inside p_next.
        p_set%elements(1:p_set%dim-1) = p_set%elements(2:p_set%dim)
      ELSE
        ! There are no elements to shift. Update 'dim' of previous node.
        p_set%prev%dim = p_set%prev%dim - 1
        EXIT
      END IF
      ! Move on to the next node if any. If there are no more nodes, update
      ! 'dim'.
      IF (ASSOCIATED(p_set%next)) THEN
        p_set => p_set%next
      ELSE
        p_set%dim = p_set%dim - 1
        EXIT
      END IF
    END DO
  ELSE
    ! There are no  more nodes. Update 'dim' of 'p_set'.
    p_set%dim = p_set%dim - 1
  END IF
END IF
RETURN
END SUBROUTINE rmNodeR

END MODULE VTDIRECT_GLOBAL

MODULE VTDIRECT_COMMSUB
USE VTDIRECT_GLOBAL
INTERFACE OPERATOR (.lexLT.) ! Lexicographic comparison operator.
  MODULE PROCEDURE lexOrder
END INTERFACE
CONTAINS

SUBROUTINE binaryS(p_rset, diam, pos, found)
IMPLICIT NONE
! Using a binary search, match the squared diameter 'diam' with an
! element in the node 'p_rset' of 'setDia', and return the position
! 'pos' if a match is found. If there is no match, return the right
! 'pos' at which to insert 'diam'. When 'pos' is returned as 0, 'diam'
! should be inserted after all others. If 'pos' is not 0, 'diam' could
! be inserted at the position 'pos' of 'p_rset' depending on 'found'.
! See insMat().

! On input:
! p_rset - A pointer to one of the nodes in 'setDia'.
! diam   - The diameter squared to match against.
!
! On output:
! pos   - The position in 'p_rset' for a match or insertion.
! found - Status indicating whether 'diam' is found in 'p_rest' or not.
!
TYPE(real_vector), INTENT(IN) :: p_rset
REAL(KIND = R8), INTENT(IN) :: diam
INTEGER, INTENT(OUT) :: pos
LOGICAL, INTENT(OUT) :: found

! Local variables.
INTEGER :: low, mid, up

! Initialization for searching.
found = .FALSE.

! Initialize limits outside bounds of array for binary search.
low = 0
up = p_rset%dim + 1
IF (p_rset%dim > 0) THEN
  ! Check with the first and the last.
  IF (p_rset%elements(1) <= diam) THEN
    ! 'diam' is the biggest.
    IF (ABS(p_rset%elements(1) - diam)/MAX(diam, p_rset%elements(1)) <= EPS4N) THEN
      ! 'diam' is the same as the first one.
      found = .TRUE.
    END IF
    pos = 1
    RETURN
  ELSE
    IF (p_rset%elements(up-1) >= diam) THEN
      IF (ABS(p_rset%elements(up-1) - diam)/MAX(p_rset%elements(up-1), diam) &
         <= EPS4N) THEN
        ! 'diam' is the smallest one. Same as the last one in 'p_rset'.
        found = .TRUE.
        pos = up-1
      ELSE
        ! 'diam' is smaller than all in 'p_rset'. Set 'pos' 0 to insert
        ! 'diam' after all others.
        found = .FALSE.
        pos = 0
      END IF
      RETURN
    ELSE
      ! 'diam' falls in between the biggest and the smallest. Apply binary
      ! search.
      DO WHILE((low + 1) < up)
        mid = (low + up) / 2
        IF (ABS(diam - p_rset%elements(mid))/            &
          MAX(diam, p_rset%elements(mid)) <= EPS4N) THEN
          ! 'diam' found.
          up = mid
          EXIT
        END IF
        IF (diam < p_rset%elements(mid))THEN
          low = mid
        ELSE
          up = mid
        END IF
      END DO
      ! Check if it's found.
      mid = up
      IF (ABS(diam - p_rset%elements(mid))/MAX(diam, p_rset%elements(mid)) &
         <= EPS4N) THEN
        ! Found it, so assign 'mid' to 'pos' in order to insert the
        ! associated box in the same column as 'mid'.
        found = .TRUE.
        pos = mid
      ELSE
        found = .FALSE.
        IF (diam > p_rset%elements(mid) )THEN
          ! 'diam' is bigger than the one at 'mid'. Set 'pos' to be 'mid'
          ! in order to insert 'diam' before 'mid'.
          pos = mid
        ELSE
          ! 'diam' is smaller than the one at 'mid'. Set 'pos' to be one
          ! after 'mid' in order to insert 'diam' after 'mid'.
          pos = mid + 1
        END IF
      END IF
    END IF
  END IF
ELSE
  ! 'p_rset' is empty. Set 'pos'=0 to insert 'diam' at the end of 'p_rset'.
  found = .FALSE.
  pos = 0
END IF
RETURN
END SUBROUTINE binaryS

SUBROUTINE checkblinks(col, b, status)
IMPLICIT NONE
! Check if this column needs a new box link. Create one if needed.
!
! On input:
! col - The local column index.
! b   - The current link of box matrices.
!
! On output:
! b      - 'b' has the newly added box link for 'col' if needed.
! status - Return status.
!           0   Successful.
!           1   Allocation error.
!
INTEGER, INTENT(IN) :: col
TYPE(BoxMatrix), INTENT(INOUT) :: b
INTEGER, INTENT(OUT) :: status

! Local variables.
TYPE(BoxLink), POINTER :: newBoxLink, p_link
INTEGER :: i

! Set normal status.
status = 0

IF ((b%ind(col) == row_w) .AND. (.NOT. ASSOCIATED(b%sibling(col)%p))) THEN
  ! When the M part of the box column is full and there are no box links:
  ! Make a new box link linked to it. Allocate a new one.
  ALLOCATE(newBoxLink, STAT = status)
  IF (status /= 0) RETURN
  CALL initLink(newBoxLink, status)
  IF (status /= 0) RETURN
  ! This is the first box link that does not have previous link.
  NULLIFY(newBoxLink%prev)
  ! Link it as 'sibling' of this column.
  b%sibling(col)%p => newBoxLink
ELSE ! There is at least one box link.
  IF (b%ind(col) > row_w) THEN
    ! Find the last box link 'p_blink'.
    p_link => b%sibling(col)%p
    DO i=1,(b%ind(col)-1)/row_w-1
      p_link => p_link%next
    END DO
    IF ((p_link%ind == row_w) .AND. (.NOT. ASSOCIATED(p_link%next))) THEN
      ! It's full.
      ! Need a new box link. Allocate a new one.
      ALLOCATE(newBoxLink, STAT = status)
      IF (status /= 0) RETURN
      CALL initLink(newBoxLink, status)
      IF (status /= 0) RETURN
      ! Link the new box link to the last one as the next link.
      p_link%next => newBoxLink
      ! Link the last box link to the new one as the previous link.
      newBoxLink%prev => p_link
    END IF
  END IF
END IF
RETURN
END SUBROUTINE checkblinks

FUNCTION findpt(b, col, i_last, index, p_last) RESULT(p_index)
IMPLICIT NONE
! Find the pointer for the box 'index' in the column 'col' of box
! matrix 'b'. If this box 'index' is in one of the box links, record
! the pointer to the box link holding this box 'index' in 'p_last'
! and compute the box position offset 'i_last'. These two records will
! be used to resume chasing the pointers for the heap elements closer
! to the bottom.
!
! On input:
! b      - Box matrix holding the box column 'col' with the box 'index'.
! col    - Column index.
! i_last - Box position offset used in finding the starting box position
!          from the box link 'p_last'.
! index  - Box index.
! p_last - Pointer to the last box link that has been chased so far.
!
! On output:
! i_last - Updated 'i_last'.
! p_last - Updated 'p_last'.
!
TYPE(BoxMatrix), INTENT(IN), TARGET :: b
INTEGER, INTENT(IN) :: col
INTEGER, INTENT(INOUT) :: i_last
INTEGER, INTENT(IN) :: index
TYPE(BoxLink), POINTER :: p_last
TYPE(HyperBox), POINTER :: p_index

! Local variables.
TYPE(BoxLink), POINTER :: p_l   ! Pointer to a box link.
INTEGER :: i

IF (.NOT. ASSOCIATED(p_last)) THEN
  ! 'p_last' has not been set, so start from the first box matrix 'b'.
  IF (index <= row_w) THEN
    ! The box 'index' is in 'M' array of 'b'.
    p_index => b%M(index,col)
  ELSE
    ! Chase to the box link that this box belongs to.
    p_l => b%sibling(col)%p
    DO i = 1, (index-1)/row_w -1
      p_l => p_l%next
    END DO
    ! Found the box link that holds the box 'index'.
    p_index => p_l%Line(MOD(index-1, row_w)+1)
    ! Set 'p_last' and 'i_last'.
    p_last => p_l
    i_last = ((index-1)/row_w)*row_w
  END IF
ELSE
  ! Start from 'p_last', because it is the last box link that has been
  ! processed.
  p_l => p_last
  DO i = 1, (index-i_last-1)/row_w
    p_l => p_l%next
  END DO
  ! Found the box link that holds the box 'index'.
  p_index => p_l%Line(MOD(index-1,row_w)+1)
  ! Set 'p_last' and 'i_last'.
  p_last => p_l
  i_last = ((index-1)/row_w)*row_w
END IF

RETURN
END FUNCTION findpt

SUBROUTINE initLink(newBoxLink, iflag)
IMPLICIT NONE
! Initialize a new box link.
!
! On input:
!   newBoxLink - A new box link.
!
! On output:
!   newBoxLink - 'newBoxLink' with initialized structures.
!   iflag      - Return status.
!               0   Normal.
!               1   Allocation failure.
!
TYPE(BoxLink), INTENT(INOUT) :: newBoxLink
INTEGER, INTENT(OUT) :: iflag

! Local variables.
INTEGER :: i

! Initialize 'iflag'.
iflag = 0
! Allocate 'Line' of the new BoxLink.
ALLOCATE(newBoxLink%Line(row_w),STAT = iflag)
IF (iflag /= 0) RETURN
DO i = 1, row_w
  ALLOCATE(newBoxLink%Line(i)%c(N_I),STAT = iflag)
  IF (iflag /= 0) RETURN
  ALLOCATE(newBoxLink%Line(i)%side(N_I),STAT = iflag)
  IF (iflag /= 0) RETURN
END DO
! Nullify pointers 'next' and 'prev'.
NULLIFY(newBoxLink%next)
NULLIFY(newBoxLink%prev)
! Initialize the counter for boxes.
newBoxLink%ind = 0
RETURN
END SUBROUTINE initLink

SUBROUTINE  insBox(box, col, b, iflag)
IMPLICIT NONE
! Insert 'box' in column 'col' of box matrices 'b'. If all positions
! in 'col' are full, make a new box link linked to the end of
! this column.
!
! On input:
! box   - The box to be inserted.
! col   - The global column index at which to insert 'box'.
! b     - The head link of box matrices.
!
! On output:
! b     - 'b' has a newly added 'box'.
! iflag - Return status.
!          0   Normal.
!          1   Allocation failure.
!
TYPE(HyperBox), INTENT(IN) :: box
INTEGER, INTENT(IN) :: col
TYPE(BoxMatrix), INTENT(INOUT), TARGET :: b
INTEGER, INTENT(OUT) :: iflag

! Local variables.
INTEGER :: b_id, i, mycol, status
TYPE(BoxLink), POINTER :: p_blink
TYPE(BoxMatrix), POINTER :: p_b

! Initialize 'iflag' as a normal return.
iflag = 0

! Locate the box matrix in which to insert 'box'.
mycol = col
IF (mycol <= col_w) THEN
  p_b => b
ELSE
  b_id = (mycol-1)/col_w + 1
  mycol = MOD(mycol-1, col_w) + 1
  p_b => b
  DO i=1, b_id-1
    p_b => p_b%child
  END DO
END IF

! Insert the box at the end of column 'mycol' of box matrix 'p_b'.
IF (p_b%ind(mycol) < row_w) THEN
  ! There is no box links.
  p_b%M(p_b%ind(mycol)+1, mycol) = box
ELSE
  ! There are box links. Chase to the last box link.
  p_blink => p_b%sibling(mycol)%p
  DO i = 1, p_b%ind(mycol)/row_w - 1
    p_blink => p_blink%next
  END DO
  p_blink%ind = p_blink%ind + 1
  p_blink%Line(p_blink%ind) = box
END IF
! Update 'ind' of the column ('ind' of 'p_b' counts all the boxes in
! this column including the ones in its box links.).
p_b%ind(mycol) = p_b%ind(mycol) + 1

! Siftup the box in the heap.
CALL siftup(p_b, mycol, p_b%ind(mycol))

! Add a new box link if needed.
CALL checkblinks(mycol, p_b, status)
IF (status /=0) THEN ; iflag = 1 ; END IF

RETURN
END SUBROUTINE insBox

SUBROUTINE insMat(box, b, setDia, setInd, setFcol, status, isConvex)
IMPLICIT NONE
! Retrieve all box matrices to find the place at which to insert 'box'.
! Insert it in the column with the same squared diameter, or a new
! column if the squared diameter is new. In the same column, the smaller
! 'val', the earlier the position.
!
! On input:
! box     - The box to be inserted.
! b       - The head link of box matrices.
! setDia  - A linked list holding different squared diameters.
! setInd  - A linked list holding the column indices corresponding to
!           'setDia'.
! setFcol - A linked list holding free columns of box matrices.
! isConvex- A flag to indicate if 'box' is on the convex hull.
!
! On output:
! b       - 'b' has the newly added 'box' and updated counters.
! setDia  - 'setDia' has a newly added squared diameter if any and
!           updated 'dim' if modified.
! setInd  - 'setInd' has a newly added column index if needed and
!           updated 'dim' if modified.
! setFcol - 'setFcol' has a column index removed if needed and updated
!           'dim' if modified.
! status  - Status of processing.
!            0    Normal.
!            1    Allocation failures of type 'BoxMatrix'.
!            2    Allocation failures of type 'BoxLink'.
!            3    Allocation failures of type 'real_vector' or 'int_vector'.
!
!
TYPE(HyperBox), INTENT(IN) :: box
TYPE(BoxMatrix), INTENT(INOUT), TARGET :: b
TYPE(real_vector), INTENT(INOUT), TARGET :: setDia
TYPE(int_vector), INTENT(INOUT), TARGET :: setInd
TYPE(int_vector), INTENT(INOUT), TARGET :: setFcol
INTEGER, INTENT(OUT) :: status
INTEGER, INTENT(IN) :: isConvex

! Local variables.
INTEGER :: b_id, b_pos, i, iflag, j, pos
LOGICAL :: found
TYPE(BoxMatrix), POINTER :: p_b
TYPE(int_vector), POINTER :: p_iset
TYPE(real_vector), POINTER :: p_rset

! Initialization for msearchSet().
pos = 0
found = .FALSE.
NULLIFY(p_rset)
iflag = 0
status = 0

! Locate a node of 'setDia' into which 'diam' can be inserted.
p_rset => msearchSet(setDia, box%diam)
CALL binaryS(p_rset, box%diam, pos, found)
IF (found) THEN
  ! A match is found in 'p_rset' of 'setDia'.
  ! Find the corresponding node in 'setInd' to match 'p_rset'.
  p_iset => setInd
  DO i = 1, p_rset%id-1
    p_iset => p_iset%next
  END DO
  ! Insert 'box' in the column indexed by 'pos' in a node of 'setInd'.
  CALL insBox(box, p_iset%elements(pos), b, iflag)
  IF (iflag == 1 ) THEN ; status = 2 ; RETURN ; END IF
  IF (isConvex == 1) THEN ! Set 'CONVEX_BIT'.
    p_iset%flags(pos) = IBSET(p_iset%flags(pos), CONVEX_BIT)
  END IF
ELSE
  ! No match is found. It's a new squared diameter.
  IF (pos == 0)THEN
    ! Obtain a free column from 'setFcol' to insert 'box' after all
    ! other columns.
    IF (setFcol%dim > 0)THEN
      ! 'setFcol' is not empty, so pop a column from the top of 'setFcol'
      ! nodes.
      IF (setFcol%dim < col_w) THEN
        ! The head node is not full, therefore it must be the top node.
        i = setFcol%elements(setFcol%dim)
        ! Update 'dim'.
        setFcol%dim = setFcol%dim - 1
      ELSE
        ! There might be other nodes with element(s).
        p_iset => setFcol%next
        IF (ASSOCIATED(p_iset)) THEN
          ! Chase to the top node of 'setFcol' with element(s).
          DO WHILE(p_iset%dim == col_w)
            p_iset => p_iset%next
          END DO
          ! The top node could be 'p_iset' or its 'prev'.
          IF (p_iset%dim /= 0) THEN
            i = p_iset%elements(p_iset%dim)
            p_iset%dim = p_iset%dim - 1
          ELSE
            i = p_iset%prev%elements(p_iset%prev%dim)
            p_iset%prev%dim = p_iset%prev%dim - 1
          END IF
        ELSE
          ! There are no more nodes with elements. Pop a column from the
          ! head node of 'setFcol'.
          i = setFcol%elements(setFcol%dim)
          ! Update 'dim'.
          setFcol%dim = setFcol%dim -1
        END IF
      END IF
    ELSE
      ! There are no free columns, so make a new box matrix.
      CALL newMat(b, setDia, setInd, setFcol, iflag)
      IF (iflag /= 0) THEN ; status = iflag ; RETURN ; END IF
      ! Pop a column from the top of 'setFcol' for use.
      i = setFcol%elements(setFcol%dim)
      CALL rmNode(col_w, 0, setFcol%dim, setFcol)
    END IF
    ! Found the global column index 'i' at which to insert 'box'. Convert
    ! it to a local column index 'b_pos' and locate the box
    ! matrix 'p_b' at which to insert 'box'.
    IF (i <= col_w) THEN
      p_b => b
      b_pos = i
    ELSE
      b_id = (i-1)/col_w + 1
      b_pos = MOD(i-1, col_w) + 1
      p_b => b
      DO j = 1, b_id -1
        p_b => p_b%child
      END DO
    END IF
    ! Insert 'box' at the beginning of the new column 'b_pos'.
    p_b%M(1,b_pos) = box
    ! Locate the nodes in both 'setDia' and 'setInd' at which to insert
    ! the new squared diameter and the column index 'i' at the end of
    ! both linked lists ('pos' is 0).
    IF (setDia%dim < col_w)THEN
      ! There are no more nodes with elements. Assign the head node
      ! to 'p_rset'.
      p_rset => setDia
    ELSE
      ! There are other nodes to check.
      p_rset => setDia%next
      IF (ASSOCIATED(p_rset)) THEN
        ! Chase to the end of the linked list.
        DO WHILE(p_rset%dim == col_w)
          p_rset => p_rset%next
        END DO
      ELSE
        ! There are no more nodes. Assign the head node to 'p_rset'
        p_rset => setDia
      END IF
    END IF
    ! Found the node 'p_rset' of 'setDia' to insert.
    CALL insNode(col_w, box%diam, p_rset%dim+1, p_rset)
    ! Find the corresponding node in 'setInd' at which to insert the
    ! column index 'i'.
    p_iset => setInd
    DO j =1, p_rset%id -1
      p_iset => p_iset%next
    END DO
    CALL insNode(col_w, i, p_iset%dim+1, p_iset)
    ! Set 'CONVEX_BIT'.
    IF (isConvex == 1) p_iset%flags(p_iset%dim) = &
      IBSET(p_iset%flags(p_iset%dim), CONVEX_BIT)
    ! Update 'ind' of col 'b_pos' in 'p_b'.
    p_b%ind(b_pos) = 1
  ELSE
    ! 'pos' is not 0. 'p_rset' points to the right node of 'setDia' to
    ! insert the new squared diameter.
    ! Obtain a free column from 'setFcol' to insert a new column before the
    ! column indexed by the returned 'pos'.
    IF (setFcol%dim > 0)THEN
      ! 'setFcol' is not empty, so pop a column from the top of 'setFcol'
      ! nodes.
      IF (setFcol%dim < col_w) THEN
      ! The head node is not full, so it must be the top.
        i = setFcol%elements(setFcol%dim)
        setFcol%dim = setFcol%dim - 1
      ELSE
        ! There might be nodes with free columns.
        p_iset => setFcol%next
        IF (ASSOCIATED(p_iset)) THEN
          ! Chase to the top of 'setFcol' links.
          DO WHILE(p_iset%dim == col_w)
            p_iset => p_iset%next
          END DO
          ! The top node could be 'p_iset' or its 'prev'.
          IF (p_iset%dim /= 0) THEN
            i = p_iset%elements(p_iset%dim)
            p_iset%dim = p_iset%dim - 1
          ELSE
            i = p_iset%prev%elements(p_iset%prev%dim)
            p_iset%prev%dim = p_iset%prev%dim - 1
          END IF
        ELSE
          ! There are no more nodes with elements. Pop a column from the
          ! head node of 'setFcol'.
          i = setFcol%elements(setFcol%dim)
          setFcol%dim = setFcol%dim -1
        END IF
      END IF
    ELSE
      ! There are no free columns, so make a new box matrix.
      CALL newMat(b, setDia, setInd, setFcol, iflag)
      IF (iflag /= 0) THEN ; status = iflag ; RETURN ; END IF
      ! Pop a column for use.
      i = setFcol%elements(setFcol%dim)
      CALL rmNode(col_w, 0, setFcol%dim, setFcol)
    END IF
    ! Found the global column index 'i' at which to insert 'box'.
    ! Convert it to a local column index 'b_pos' and locate the
    ! box matrix 'p_b' in which to insert 'box'.
    IF (i <= col_w) THEN
      p_b => b
      b_pos = i
    ELSE
      b_id = (i-1)/col_w + 1
      b_pos = MOD(i-1, col_w) + 1
      p_b => b
      DO j = 1, b_id -1
        p_b => p_b%child
      END DO
    END IF
    ! Add 'box' to be the first on column 'b_pos' of 'p_b'.
    p_b%M(1,b_pos) = box
    ! Insert the new squared diameter at the position 'pos' in 'p_rset'
    ! of 'setDia'.
    CALL insNode(col_w, box%diam, pos, p_rset)
    ! Insert the corresponding column index 'i' at the same position 'pos'
    ! in a node of 'setInd'.
    p_iset => setInd
    DO j = 1, p_rset%id -1
      p_iset => p_iset%next
    END DO
    CALL insNode(col_w, i, pos, p_iset)
    IF (isConvex == 1) p_iset%flags(pos) = &
      IBSET(p_iset%flags(pos), CONVEX_BIT)
    ! Update 'ind' of box matrix 'p_b'.
    p_b%ind(b_pos) = 1
  END IF
END IF

RETURN
END SUBROUTINE insMat

FUNCTION lexOrder(a, b) RESULT (smaller)
IMPLICIT NONE
! Compare box centers 'a' and 'b' lexicographically. If 'a' is
! lexicographically smaller than 'b', return .TRUE. .  Define operator
! .lexLT. .
!
! On input:
! a - center of a box.
! b - center of a box.
!
! On output:
! smaller - .TRUE.   'a' is smaller than 'b' lexicographically.
!           .FALSE.  Otherwise.
!
REAL(KIND = R8), DIMENSION(:), INTENT(IN) :: a
REAL(KIND = R8), DIMENSION(:), INTENT(IN) :: b
LOGICAL :: smaller

! Local variables.
INTEGER :: i

! Initialize 'smaller'.
smaller = .FALSE.

! Loop for comparing coordinates of 'a' and 'b'.
DO i = 1,SIZE(a)
  IF (a(i) < b(i)) THEN
    smaller = .TRUE.
    EXIT
  ELSE IF (a(i) == b(i)) THEN
    CONTINUE
  ELSE ! a(i) > b(i)
    EXIT
  END IF
END DO
RETURN
END FUNCTION lexOrder

FUNCTION msearchSet(setDia, diam) RESULT(p_rset)
IMPLICIT NONE
! Find the right node in 'setDia' in which to insert 'diam'.
!
! On input:
! setDia - A linked list holding different squared diameters.
! diam   - A diameter squared to be inserted in a node in 'setDia'.
!
! On output:
! p_rset - Pointer to the right node of 'setDia'.
!
TYPE(real_vector), INTENT(IN), TARGET :: setDia
REAL(KIND = R8), INTENT(IN) :: diam
TYPE(real_vector), POINTER :: p_rset

! Local variables.
TYPE(real_vector), POINTER :: p_setDia

! Initialize 'p_setDia'.
p_setDia => setDia
DO
  IF (p_setDia%dim > 0) THEN
    ! There are elements to be compared with 'diam'.
    IF (diam >= p_setDia%elements(1)) THEN
      ! 'diam' is the biggest. Return this node as 'p_rset'.
      p_rset => p_setDia
      EXIT
    ELSE
      IF ((diam >= p_setDia%elements(p_setDia%dim)) .OR.     &
         ((ABS(diam - p_setDia%elements(p_setDia%dim))        &
         /MAX(diam, p_setDia%elements(p_setDia%dim))) <= EPS4N) ) THEN
        ! 'diam' is within the range of elements in this node.
        p_rset => p_setDia
        EXIT
      ELSE
        ! 'diam' is smaller than the last element. Go on to the next
        ! node, if any.
        IF (ASSOCIATED(p_setDia%next)) THEN
          p_setDia => p_setDia%next
        ELSE
          ! There are no more nodes. Return this pointer.
          p_rset => p_setDia
          EXIT
        END IF
      END IF
    END IF
  ELSE
    ! It's empty. Return this pointer.
    p_rset => p_setDia
    EXIT
  END IF
END DO
RETURN
END FUNCTION msearchSet

SUBROUTINE newMat(b, setDia, setInd, setFcol, iflag)
IMPLICIT NONE
! Make a new box matrix and its associated linked lists for holding
! different squared diameters, column indices, and free columns.
! Link them to existing structures.
!
! On input:
! b       - The head link of box matrices.
! setDia  - Linked list holding different squared diameters of
!           current boxes.
! setInd  - Linked list holding column indices of box matrices.
! setFcol - Linked list holding free columns of box matrices.
!
! On output:
! b       - 'b' has a box matrix link added at the end.
! setDia  - 'setDia' has a node added at the end.
! setInd  - 'setInd' has a node added at the end.
! setFcol - 'setFcol' has a node added at the end.
! iflag   - Return status.
!           0   Normal.
!           1   Allocation failures of type 'BoxMatrix'.
!           3   Allocation failures of type 'real_vector' or 'int_vector'.
!
TYPE(BoxMatrix), INTENT(INOUT), TARGET :: b
TYPE(real_vector), INTENT(INOUT), TARGET :: setDia
TYPE(int_vector), INTENT(INOUT), TARGET :: setInd
TYPE(int_vector), INTENT(INOUT), TARGET :: setFcol
INTEGER, INTENT(OUT) :: iflag

! Local variables.
INTEGER :: alloc_err, i, j
TYPE(BoxMatrix), POINTER :: new_b, p_b
TYPE(int_vector), POINTER :: n_setFcol, n_setInd, p_setFcol, p_setInd
TYPE(real_vector), POINTER :: n_setDia, p_setDia

! Initialize iflag for normal return.
iflag = 0
! Allocate a new box matrix.
ALLOCATE(new_b, STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 1; RETURN; END IF
! Allocate its associated arrays.
ALLOCATE(new_b%M(row_w, col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 1; RETURN; END IF
ALLOCATE(new_b%ind(col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 1; RETURN; END IF

! Clear the box counter for each column.
new_b%ind(:) = 0
! Nullify pointers for a child link and box links.
NULLIFY(new_b%child)
ALLOCATE(new_b%sibling(col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 1; RETURN; END IF
DO i = 1, col_w
  NULLIFY(new_b%sibling(i)%p)
END DO
DO i = 1, row_w
  DO j = 1, col_w
    ALLOCATE(new_b%M(i,j)%c(N_I), STAT=alloc_err)
    IF (alloc_err /= 0) THEN; iflag = 1; RETURN; END IF
    ALLOCATE(new_b%M(i,j)%side(N_I), STAT=alloc_err)
    IF (alloc_err /= 0) THEN; iflag = 1; RETURN; END IF
  END DO
END DO

! Find the last box matrix to link with the new one.
p_b => b
DO WHILE(ASSOCIATED(p_b%child))
  p_b => p_b%child
END DO
! Found the last box matrix p_b. Link it to new box matrix 'b'.
p_b%child => new_b
! Set up 'id' for new_b.
new_b%id = p_b%id + 1

! Allocate new 'setDia', 'setInd' and 'setFcol' for 'new_b'.
! Find the corresponding nodes of 'setDia', 'setInd' and 'setFcol'.
p_setDia => setDia
p_setInd => setInd
p_setFcol => setFcol
DO i=1, p_b%id-1
  p_setDia => p_setDia%next
  p_setInd => p_setInd%next
  p_setFcol => p_setFcol%next
END DO
! Allocate a new node for 'setDia'. Initialize its structure.
ALLOCATE(n_setDia, STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
ALLOCATE(n_setDia%elements(col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
NULLIFY(n_setDia%next)
NULLIFY(n_setDia%prev)
n_setDia%id = p_setDia%id +1
n_setDia%dim = 0
! Allocate a new node for 'setInd'. Initialize its structure.
ALLOCATE(n_setInd, STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
ALLOCATE(n_setInd%elements(col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
ALLOCATE(n_setInd%flags(col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
n_setInd%flags(:)= 0
NULLIFY(n_setInd%next)
NULLIFY(n_setInd%prev)
n_setInd%id = p_setInd%id + 1
n_setInd%dim = 0
! Allocate a new node for 'setFcol'. Initialize its structure.
ALLOCATE(n_setFcol, STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
ALLOCATE(n_setFcol%elements(col_w), STAT=alloc_err)
IF (alloc_err /= 0) THEN; iflag = 3; RETURN; END IF
NULLIFY(n_setFcol%next)
NULLIFY(n_setFcol%prev)
n_setFcol%id = p_setFcol%id + 1
n_setFcol%dim = 0
! Link them to the end of existing sets.
p_setDia%next => n_setDia
n_setDia%prev => p_setDia
p_setInd%next => n_setInd
n_setInd%prev => p_setInd
p_setFcol%next => n_setFcol
n_setFcol%prev => p_setFcol
! Fill up 'setFcol' with new columns from the new box matrix.
! Starting from the last column, push free columns to the top of 'setFcol'.
DO i = 1, col_w
  setFcol%elements(i) = new_b%id*col_w - (i-1)
END DO
setFcol%dim = col_w
RETURN
END SUBROUTINE newMat

SUBROUTINE siftdown(b, col, index)
IMPLICIT NONE
! Sift down the heap element at 'index' through the heap column 'col' in
! the box matrix 'b'.
!
! On input:
! b     - Box matrix holding the box column 'col' for siftdown.
! col   - Column index.
! index - Index of the box to be sifted down.
!
! On output:
! b     - Box matrix with the elements rearranged by siftdown.
!
TYPE(BoxMatrix), INTENT(INOUT), TARGET :: b
INTEGER, INTENT(IN) :: col, index

! Local variables.
INTEGER :: i           ! Loop counters.
INTEGER :: i_last      ! Index of the last box that has been processed
                       ! in the previous iteration.
INTEGER :: i_last_backup  ! i_last's backup used to go back to i_last,
                          ! which may be updated in the current iteration.
INTEGER :: left, right    ! Indices for the left and right children.
TYPE(BoxLink), POINTER :: p_last        ! Pointer to the last box link
                                        ! that has been processed.
TYPE(BoxLink), POINTER :: p_last_backup ! Pointer backup for 'p_last'.
TYPE(HyperBox), POINTER :: p_i          ! Pointer to the heap parent box.
TYPE(HyperBox), POINTER :: p_left, p_right  ! Pointers to the left and right
                                            ! child boxes.
! Exit if it is an empty column.
IF (b%ind(col) == 0) RETURN

! Initialization.
NULLIFY(p_last)
i_last=0

! Starting siftdown operation from the box 'index'
i = index
DO
! Find the indices for the left and right children.
left = 2*i
right = 2*i + 1
! If 'i' is a leaf, exit.
IF (left > b%ind(col)) EXIT
! Find the pointers to the the ith box and its left child. Record the
! pointer 'p_last' and index 'i_last' for the currently processed box
! link.
  p_i => findpt(b, col, i_last, i, p_last)
  p_left => findpt(b, col, i_last, left, p_last)
  IF (left < b%ind(col)) THEN
    ! Backup the pointer 'p_last' and the index 'i_last', because they will
    ! be updated when finding the pointer to the right child. If the right
    ! child is in the correct place in the heap, restore the pointer
    ! 'p_last' and the index 'i_last'.
    p_last_backup => p_last
    i_last_backup = i_last
    p_right => findpt(b, col, i_last, right, p_last)
    IF (p_left%val > p_right%val) THEN
      p_left => p_right
      left = left + 1
    ELSE
      IF (p_left%val == p_right%val) THEN
        ! When values are equal, smaller lex order wins.
        IF (p_right%c .lexLT. p_left%c) THEN
          p_left => p_right
          left = left + 1
        ELSE
          ! Restore 'p_last' and 'i_last' to the values before finding
          ! the pointer for the right child.
          p_last => p_last_backup
          i_last = i_last_backup
        END IF
      ELSE
        ! Restore 'p_last' and 'i_last' to the values before finding
        ! the pointer for the right child.
        p_last => p_last_backup
        i_last = i_last_backup
      END IF
    END IF
  END IF
  ! If the boxes are in the correct order, exit.
  IF (p_i%val < p_left%val) EXIT
  IF ((p_i%val == p_left%val) .AND. (p_i%c .lexLT. p_left%c)) EXIT
  ! Swap the boxes pointed to by 'pi' and 'p_left'.
  tempbox = p_i
  p_i = p_left
  p_left = tempbox
  ! Continue siftdown operation from the left child of 'i'.
  i = left
END DO

RETURN
END SUBROUTINE siftdown

SUBROUTINE siftup(b, col, index)
IMPLICIT NONE
! Sift up the heap element at 'index' through the heap column 'col' in
! the box matrix 'b'.
!
! On input:
! b     - Box matrix holding the box column 'col' for siftup.
! col   - Column index.
! index - Index of the box to be sifted up.
!
! On output:
! b     - Box matrix with the elements rearranged by siftup.
!
TYPE(BoxMatrix), INTENT(INOUT), TARGET :: b
INTEGER, INTENT(IN) :: col, index

! Local variables.
INTEGER :: i, j      ! Loop counters.
INTEGER :: parent    ! Index for the parent.
INTEGER :: e_idi, e_idp, l_idi, l_idp
TYPE(BoxLink), POINTER :: p_blinki, p_blinkp ! Pointers for the box link.
TYPE(HyperBox), POINTER :: p_i  ! Pointer to the box indexed as 'i'.
TYPE(HyperBox), POINTER :: p_p  ! Pointer to the parent box.

! Exit if it is an empty column.
IF (b%ind(col) == 0) RETURN

! Starting siftup operation from the box 'index'.
i = index
IF (i == 1) RETURN ! Exit if there is only one box in the column.
l_idi = (i-1)/row_w  ! If 0, it's inside the box matrix M.
e_idi = MOD(i-1,row_w)+1
! Locate the box 'i'.
IF (l_idi == 0) THEN
  p_i => b%M(i, col)
ELSE
  ! Chase to the box link that this box belongs to.
  p_blinki => b%sibling(col)%p
  DO j = 1, l_idi-1
    p_blinki => p_blinki%next
  END DO
  ! Found the box link that holds the box 'i'.
  p_i => p_blinki%Line(e_idi)
END IF

DO WHILE (i /= 1) ! If the root has been reached, exit.
  ! Find the index for the parent.
  parent = i/2
  ! Compute link id and element id for the parent.
  l_idp = (parent-1)/row_w
  e_idp = MOD(parent-1,row_w)+1

  ! Locate the box 'parent'.
  IF (l_idp == 0) THEN
    p_p => b%M(parent, col)
  ELSE
    p_blinkp => b%sibling(col)%p
    DO j = 1, l_idp-1
      p_blinkp => p_blinkp%next
    END DO
    p_p => p_blinkp%Line(e_idp)
  END IF

  ! Found the boxes 'parent' and 'i'. Compare their 'val's. If box 'i' has
  ! a smaller 'val', swap box 'i' with box 'parent'. Otherwise, box 'i'
  ! moves up to point to the box 'parent' and repeats the comparison and
  ! swaps (if needed), until box 'i' becomes the root.
  IF ((p_i%val < p_p%val) .OR. &
     ((p_i%val == p_p%val) .AND. (p_i%c .lexLT. p_p%c))) THEN
    tempbox = p_i
    p_i = p_p
    p_p = tempbox
    l_idi = l_idp
    e_idi = e_idp
    i = parent
    p_i => p_p
    IF (l_idi /= 0) p_blinki => p_blinkp
  ELSE
    i = 1
  END IF
END DO

RETURN
END SUBROUTINE siftup
END MODULE VTDIRECT_COMMSUB

