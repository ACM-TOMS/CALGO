SUBROUTINE SSYR_RMD(uplo, n, alpha, x, incx, lda, xa, Aa)
  integer :: n, incx, lda
  real :: alpha, x(*), xa(*), Aa(lda,n)
  character :: uplo

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSYR from BLAS. 
  !
  !
  ! ARGUMENTS
  !    If SSYR was called with the arguments
  !
  !       uplo, n, alpha, x, incx, A, lda
  !    
  !    then the corresponding call to SSYR_RMD should begin with the arguments
  !
  !       uplo, n, alpha, x, incx, lda
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that A is omitted. In addition the following arguments should be
  !    provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SSYR call.
  !
  !    Aa
  !       (input, real triangular matrix of the same dimensions as A, and stored
  !       in the same half according to uplo)
  !       The adjoint of A.
  !
  ! NOTE
  !    Aa is always unchanged so that a sel parameter is not needed.
  !
  ! OPERATIONS
  !    BLAS: A += alpha*tril(x*x') i.e. sym(A) := alpha*x*x' + sym(A)
  !    RMD:  xa += alpha*(Aa + Aa')*x, equiv.to: xa += alpha*(diag(Aa)+sym(Aa))*x
  !          Aa unchanged

  call sgbmv('n', n, n, 0, 0, alpha, Aa, lda+1, x, incx, 1.0, xa, incx)
  call ssymv(uplo, n, alpha, Aa, lda, x, incx, 1.0, xa, incx)
END SUBROUTINE SSYR_RMD
