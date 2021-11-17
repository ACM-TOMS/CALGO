SUBROUTINE SSPR2_RMDS(uplo, n, x, incx, y, incy, alphaa, APa)
  integer :: n, incx, incy
  real :: alphaa, x(*), y(*), APa(*)
  character :: uplo

  ! PURPOSE
  !    Calculate the adjoint of alpha for SSPR2 from BLAS
  !
  !
  ! ARGUMENTS
  !    If SSPR2 was called with the arguments:
  !
  !       uplo, n, alpha, x, incx, y, incy, AP
  !
  !    then the corresponding call to SSPR2_RMD should begin with the arguments
  !
  !       uplo, n, x, incx, y, incy
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that alpha and AP are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SSPR2 call.
  !
  !    APa
  !       (input, real packed triangular matrix of the same dimensions as AP,
  !       and stored in the same half according to uplo)
  !       The adjoint of AP.
  !
  ! OPERATIONS
  !    BLAS: (with uplo = 'L')
  !          AP += tril(alpha*x*y' + alpha*y*x'), where AP is lower triangluar packed
  !          i.e. sym(AP) := alpha*(x*y' + y*x') + sym(AP)
  !          (with uplo = 'U')
  !          AP += triu(alpha*x*y' + alpha*y*x'), where AP is upper triangluar packed
  !          i.e. sym(AP') := alpha*(x*y' + y*x') + sym(AP')
  !    RMD:  (with uplo = 'L' or 'U':)
  !          alphaa += x'*APa*y + y'*APa*x

  integer :: i, kx, ky, kA
  real, external :: sdot
    
  kx = 1
  ky = 1
  kA = 1
  if (uplo == 'l' .or. uplo == 'L') then
    do i = 1, n
      alphaa = alphaa + sdot(n+1-i, x(kx), incx, APa(kA), 1)*y(ky)
      alphaa = alphaa + sdot(n+1-i, y(ky), incy, APa(kA), 1)*x(kx)
      kx = kx + incx
      ky = ky + incy
      kA = kA + n+1-i
    enddo
  else
    do i = 1, n
      alphaa = alphaa + sdot(i, x, incx, APa(kA), 1)*y(ky)
      alphaa = alphaa + sdot(i, y, incy, APa(kA), 1)*x(kx)
      kx = kx + incx
      ky = ky + incy
      kA = kA + i
    enddo
  endif

END SUBROUTINE SSPR2_RMDS
