SUBROUTINE SROT_RMD(n, x, incx, y, incy, c, s, xa, ya, ca, sa, sel)
  implicit none
  integer, intent(in) :: n, incx, incy
  real, intent(in) :: x(*), y(*), c, s
  real, intent(inout) :: xa(*), ya(*), ca, sa
  character(2), intent(in) :: sel

  ! PURPOSE
  !    Calculates the reverse mode derivative of SROT from BLAS.
  !
  ! ARGUMENTS
  !    If SROT was called with the arguments
  !
  !       n, x, incx, y, incy, c, s
  !
  !    then the corresponding call to SROTG_RMD should begin with the arguments
  !
  !       n, x, incx, y, incy, c, s
  !
  !    with the same values as they had on exit from the SROT-call. All these
  !    arguments will remain unchanged on exit. In addition the following
  !    arguments should be provided:
  !
  !       xa
  !          (input, output, real vector of the same dimension and increment as x)
  !          On entry: The adjoint of the x produced by SROT
  !          On exit: The adjoint of the x supplied to SROT
  !
  !       ya
  !          (input, output, real vector of the same dimension and increment as y)
  !          On entry: The adjoint of the y produced by SROT
  !          On exit: The adjoint of the y supplied to SROT
  !
  !       ca
  !          (input, output, real scalar)
  !          ca += the adjoint of c due to the SROT-call
  !
  !       sa
  !          (input, output, real scalar)
  !          ca += the adjoint of c due to the SROT-call
  !
  !       sel
  !          (input, character*3)
  !          Used to select which adjoints to update/compute:
  !            sel(1:1) = '1' if (xa,ya) should be computed, else sel(1:1) = '0'
  !            sel(2:2) = '1' if (ca,sa) should be updated, else sel(2:2) = '0'
  !          For example, to update only (xa,ya), set sel = '10'.
  !
  ! OPERATIONS
  !    BLAS: [x'; y'] := G*[x'; y']
  !    RMD:  [xa'; ya'] := G'*[xa' ya']
  !          [ca; sa] += G'*[c1; s1]
  !          where:
  !             c1 = dot(xa,x) + dot(ya,y)  (xa, ya are values on entry)
  !             s1 = dot(xa,y) - dot(ya,x)  (xa, ya are values on entry)
  !             and G = [c s; -s c]

  ! Notes:
  !   Ga = [xa';ya']*[x y] where x, y are values supplied to SROT
  ! or:
  !   Ga = [xa';ya']*[x,y]G where x, y are values produced by SROT
  ! Also:
  !   ca = Ga(1,1) + Ga(2,2) and
  !   sa = Ga(1,2) - Ga(2,1)
  ! Combined, we obtain the RMD-operations provided above.
  
  real, external :: sdot
  real :: c1, s1;

  if (sel(2:2) == '1') then
    ! Here the updating of ca and sa has been rewritten in terms of the values
    ! that x and y have on exit from srot. 
    c1 = sdot(n, xa, incx, x, incx) + sdot(n, ya, incy, y, incy)
    s1 = sdot(n, xa, incx, y, incy) - sdot(n, ya, incy, x, incx)
    ca = ca + c1*c - s1*s
    sa = sa + s1*c + c1*s
  endif

  if (sel(1:1) == '1') then
    call srot(n, xa, incx, ya, incy, c, -s)
  endif
  
END SUBROUTINE SROT_RMD
