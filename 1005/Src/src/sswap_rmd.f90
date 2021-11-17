SUBROUTINE SSWAP_RMD(n, incx, incy, xa, ya)
  implicit none
  integer :: n, incx, incy
  real :: xa(*), ya(*)

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSWAP from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSWAP was called with the arguments
  !
  !       n, x, incx, y, incy
  !
  !    then the corresponding call to SSWAP_RMD should begin with the arguments
  !
  !       n, incx, incy
  !
  !    with the same values. These arguments will remain unchanged on exit. In
  !    addition the following arguments should be provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa := the value of ya on entry
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       ya := the value of xa on entry
  !  
  ! OPERATIONS
  !    BLAS: x <--> y
  !    RMD:  xa <--> ya

  call sswap(n, ya, incy, xa, incx)
END SUBROUTINE SSWAP_RMD
