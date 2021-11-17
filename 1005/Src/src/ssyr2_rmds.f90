SUBROUTINE SSYR2_RMDS(uplo, n, x, incx, y, incy, lda, alphaa, Aa)
  integer :: n, incx, incy, lda
  real :: alphaa, x(*), y(*), Aa(lda,*)
  character :: uplo

  ! PURPOSE
  !    Calculate the adjoint of alpha for SSYR2 from BLAS
  !
  !
  ! ARGUMENTS
  !    If SSYR2 was called with the arguments:
  !
  !       uplo, n, alpha, x, incx, y, incy, A, lda
  !
  !    then the corresponding call to SSYR2_RMD should begin with the arguments
  !
  !       uplo, n, x, incx, y, incy, lda
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that alpha and A are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SSYR2 call.
  !
  !    Aa
  !       (input, real triangular matrix of the same dimensions as A, and stored
  !       in the same half according to uplo)
  !       The adjoint of A.
  !
  ! OPERATIONS
  !    BLAS: (with uplo = 'L')
  !          A += alpha*tril(x*y' + y*x'), where A is lower triangluar
  !          i.e. sym(A) := alpha*(x*y' + y*x') + sym(A)
  !          (with uplo = 'U')
  !          A += alpha*triu(x*y' + y*x'), where A is upper triangluar
  !          i.e. sym(A') := alpha*(x*y' + y*x') + sym(A')
  !    RMD:  (with uplo = 'L' or 'U':)
  !          alphaa += x'*Aa*y + y'*Aa*x

  integer :: i, kx, ky
  real, external :: sdot
  
  kx = 1
  ky = 1
  if (uplo == 'l' .or. uplo == 'L') then
    do i = 1, n
      alphaa = alphaa + sdot(n+1-i, x(kx), incx, Aa(i,i), 1)*y(ky)
      alphaa = alphaa + sdot(n+1-i, y(ky), incy, Aa(i,i), 1)*x(kx)
      kx = kx + incx
      ky = ky + incy
    enddo
  else
    do i = 1, n
      alphaa = alphaa + sdot(i, x, incx, Aa(1,i), 1)*y(ky)
      alphaa = alphaa + sdot(i, y, incy, Aa(1,i), 1)*x(kx)
      kx = kx + incx
      ky = ky + incy
    enddo
  endif

END SUBROUTINE SSYR2_RMDS
