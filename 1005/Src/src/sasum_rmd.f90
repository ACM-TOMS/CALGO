SUBROUTINE SASUM_RMD(n, x, incx, xa, ba)
  implicit none
  integer :: n, incx
  real :: ba
  real :: x(*), xa(*)

  ! PURPOSE
  !    Calculate the reverse mode derivative of SASUM from BLAS.
  !
  ! ARGUMENTS
  !    If SASUM was called with the statement
  !
  !       b = SASUM(n, x, incx)
  !
  !    then SASUM_RMD should be called with the arguments:
  !    
  !       n, x, incx
  !    
  !    with the values which they had on the SASUM call. These arguments will
  !    remain unchanged on exit. Note that b is omitted. In addition the
  !    following arguments should be provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SASUM call.
  !
  !    ba
  !       (input, real scalar)
  !       The adjoint of b
  !
  ! NOTES
  !    SASUM is not differentiable for x-elements which are 0. The adjoints
  !    of such elements is returned as 0.
  !
  ! OPERATIONS
  !    BLAS:  b := sum |x(i)|  (1-norm)
  !    RMD:   xa += |ba|*sign(x)    
  !           ba unchanged
  !           where sign(x) = 1 where x > 0, -1 where x < 0 and 0 where x = 0

  integer :: i, k
  k = 1
  do i = 1,n
    if (x(k) /= 0) xa(k) = xa(k) + sign(ba, x(k))
    k = k + incx
  enddo

END SUBROUTINE SASUM_RMD
