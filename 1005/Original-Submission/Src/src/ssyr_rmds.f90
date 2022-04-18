SUBROUTINE SSYR_RMDS(uplo, n, x, incx, lda, alphaa, Aa)
  integer :: n, incx, lda
  real :: alphaa, x(*), Aa(lda,n)
  character :: uplo

  ! PURPOSE
  !    Calculate the adjoint of alpha for SSYR from BLAS
  !
  !
  ! ARGUMENTS
  !    If SSYR was called with the arguments
  !
  !       uplo, n, alpha, x, incx, A, lda
  !    
  !    then the corresponding call to SSYR_RMD should begin with the arguments
  !
  !       uplo, n, x, incx, lda
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that alpha and A are omitted. In addition the following arguments
  !    should be provided:
  !
  !    alphaa
  !       (input, output, real scalar)
  !       alphaa += the adjoint of alpha due to the SSYR call.
  !
  !    Aa
  !       (input, real triangular matrix of the same dimensions as A, and stored
  !       in the same half according to uplo)
  !       The adjoint of A.
  !
  ! OPERATIONS
  !    (with uplo = 'L')
  !    BLAS: A += alpha*tril(x*x') i.e. sym(A) := alpha*x*x' + sym(A)
  !    RMD:  alphaa += x'*Aa*x
  !    (with uplo = 'U')
  !    BLAS: A += alpha*triu(x*x') i.e. sym(A') := alpha*x*x' + sym(A')
  !    RMD:  alphaa += x'*Aa*x

  integer :: i, k
  real, external :: sdot

  k = 1
  if (uplo == 'l' .or. uplo == 'L') then
    do i = 1, n
      alphaa = alphaa + sdot(n+1-i, x(k), incx, Aa(i,i), 1)*x(k)
      k = k + incx
    enddo
  else
    do i = 1, n
      alphaa = alphaa + sdot(i, x, incx, Aa(1,i), 1)*x(k)
      k = k + incx
    enddo
  endif

END SUBROUTINE SSYR_RMDS
