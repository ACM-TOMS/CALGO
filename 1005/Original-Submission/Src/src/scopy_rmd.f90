SUBROUTINE SCOPY_RMD(n, incx, incy, xa, ya)
  implicit none
  integer :: n, incx, incy
  real :: xa(*), ya(*)

  ! PURPOSE
  !    Calculate the reverse mode derivative of SCOPY from BLAS.
  !
  ! ARGUMENTS
  !    If SCOPY was called with the arguments
  !
  !       n, x, incx, y, incy,
  !
  !    then the call to SCOPY_RMD should begin with the arguments
  !
  !       n, incx, incy
  !
  !    with the same values. These arguments will remain unchanged on exit. In
  !    addition the following arguments should be provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += adjoint of x
  !
  !    ya
  !       (input, real vector of the same dimension and increment as y)
  !       ya += adjoint of y
  !
  ! NOTE
  !    As there is only one output there is no need for a sel parameter
  !
  ! OPERATIONS
  !    BLAS: y = x
  !    RMD:  xa += ya
  !          ya unchanged
  
  call saxpy(n, 1.0, ya, incy, xa, incx)
END SUBROUTINE SCOPY_RMD
