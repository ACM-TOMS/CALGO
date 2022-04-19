SUBROUTINE SSCAL_RMD(n, alpha, incx, xa)
  integer :: n, incx
  real    :: alpha, xa(*)

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSCAL from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSCAL was called with the arguments
  !
  !       n, alpha, x, incx
  !
  !    then SSCAL_RMD should be called with the arguments
  !
  !       n, alpha, incx
  !
  !    with the same values. These arguments will remain unchanged on exit. Note
  !    that x is omitted. In addition the following parameter should be provided:
  !
  !    xa
  !       (output, real vector of the same dimension and increment as x)
  !       xa := adjoint of x
  !
  ! OPERATIONS
  !    BLAS: x := alpha*x
  !    RMD:  xa := alpha*xa

  call sscal(n, alpha, xa, incx)
END SUBROUTINE SSCAL_RMD
