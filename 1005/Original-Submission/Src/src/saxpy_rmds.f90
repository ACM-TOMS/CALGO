SUBROUTINE SAXPY_RMDS(n, x, incx, incy, alphaa, ya)
  implicit none
  integer :: n, incx, incy
  real    :: alphaa, x(*), ya(*)

  ! PURPOSE
  !    Calculates the adjoint of alpha for SAXPY from BLAS.
  !
  ! ARGUMENTS
  !    If SAXPY was called with the arguments
  !
  !       n, alpha, x, incx, y, incy
  !
  !    then the corresponding call to SAXPY_RMDS should begin with the arguments
  !
  !       n, x, incx, incy
  !
  !    with the same values. These arguments will remain unchanged on exit. Note
  !    that alpha and y are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SAXPY call.
  !
  !    ya
  !       (input, real vector of the same dimension and increment as y)
  !       The adjoint of the y produced by SAXPY.
  !
  ! OPERATIONS
  !    BLAS: y := alpha*x + y
  !    RMD:  alphaa += ya'*x

  real, external :: sdot
  
  alphaa = alphaa + sdot(n, ya, incy, x, incx)

END SUBROUTINE SAXPY_RMDS
