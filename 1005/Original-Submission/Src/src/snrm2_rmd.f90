SUBROUTINE SNRM2_RMD(n, x, incx, b, xa, ba)
  implicit none
  integer :: n, incx
  real :: b, ba
  real :: x(*), xa(*)

  ! PURPOSE
  !    Calculate the reverse mode derivative of SNRM2 from BLAS.
  !
  ! ARGUMENTS
  !    If SNRM2 was called with the statement
  !
  !       b = SNRM2(n, x, incx)
  !
  !    then SNRM2_RMD should be called with the same arguments:
  !    
  !       n, x, incx
  !    
  !    with the same values. These arguments will remain unchanged on exit. In
  !    addition the following arguments should be provided:
  !
  !    b
  !       (input, real scalar)
  !       The result of the SNRM2 call.
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SDOT call.
  !
  !    ba
  !       (input, real scalar)
  !       The adjoint of b
  !
  ! NOTE
  !    As there is only one vector input a sel parameter is not needed.
  !
  ! OPERATIONS
  !    BLAS:  b := norm(x)  (2-norm)
  !    RMD:   xa += ba*x/b
  !           ba unchanged

  call saxpy(n, ba/b, x, incx, xa, incx)

END SUBROUTINE SNRM2_RMD
