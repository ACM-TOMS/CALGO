SUBROUTINE SAXPY_RMD(n, alpha, incx, incy, xa, ya)
  implicit none
  integer :: n, incx, incy
  real    :: alpha, xa(*), ya(*)

  ! PURPOSE
  !    Calculates the reverse mode derivative of SAXPY from BLAS.
  !
  ! ARGUMENTS
  !    If SAXPY was called with the arguments
  !
  !       n, alpha, x, incx, y, incy
  !
  !    then the corresponding call to SAXPY_RMD should begin with the arguments
  !
  !       n, alpha, incx, incy
  !
  !    with the same values. These arguments will remain unchanged on exit. In
  !    addition the following arguments should be provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SAXPY call.
  !
  !    ya
  !       (input, real vector of the same dimension and increment as y)
  !       The adjoint of y.
  !
  ! NOTE
  !    SAXPY computes y := alpha*x + y, so that the adjoint of y is unchanged,
  !    and needs no update. Therefore a sel parameter is not needed (it is
  !    implicitly assumed to be '1X', to update xa, X may be 0 or 1 because ya
  !    is unchanged)
  !
  ! OPERATIONS
  !    BLAS: y := alpha*x + y
  !    RMD:  xa += alpha*ya
  !          ya unchanged

  call saxpy(n, alpha, ya, incy, xa, incx)
END SUBROUTINE SAXPY_RMD
