SUBROUTINE SSCAL_RMDS(n, x0, incx, alphaa, xa)
  integer :: n, incx
  real    :: x0(*), xa(*), alphaa

  ! PURPOSE
  !    Calculate the adjoint of alpha for SSCAL from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSCAL was called with the arguments
  !
  !       n, alpha, x, incx
  !
  !    then SSCAL_RMDS should be called with the arguments
  !
  !       n, x0, incx
  !
  !    which all except x0 should have the same values as they had on the SSCAL-
  !    call, and x0 should have the value that x had on entry to the SSCAL-call
  !    (SSCAL only changes x). All these arguments will remain unchanged on
  !    exit. Note that alpha is omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += adjoint of alpha due to the SSCAL-call
  !
  !    xa
  !       (input, real vector if the same dimension and increment as x)
  !       The adjoint of the x produced by SSCAL
  !
  ! OPERATIONS
  !    BLAS: x := alpha*x0
  !    RMD:  alphaa += xa'*x0

  real, external :: sdot

  alphaa = alphaa + sdot(n, xa, incx, x0, incx)

END SUBROUTINE SSCAL_RMDS
