SUBROUTINE SSYR2_RMD(uplo, n, alpha, x, incx, y, incy, lda, xa, ya, Aa, sel)
  implicit none
  character :: uplo
  integer :: n, incx, incy, lda
  real :: alpha, x(*), y(*), xa(*), ya(*), Aa(lda,*)
  character(3) :: sel

  ! PURPOSE
  !    SSYR2_RMD calculates the reverse mode derivative of the SSYR2 routine from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSYR2 was called with the arguments:
  !
  !       uplo, n, alpha, x, incx, y, incy, A, lda
  !
  !    then the corresponding call to SSYR2_RMD should begin with the arguments
  !
  !       uplo, n, alpha, x, incx, y, incy, lda
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that A is omitted. In addition the following arguments should be
  !    provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SSYR2 call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       ya += the adjoint of y due to the SSYR2 call.
  !
  !    Aa
  !       (input, real triangular matrix of the same dimensions as A, and stored
  !       in the same half according to uplo)
  !       The adjoint of A.
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !          sel(1:1) = '1' if xa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if ya should be updated, else sel(2:2) = '0'
  !       For example, to update only xa, set sel = '10'.
  !
  ! OPERATIONS
  !    BLAS: (with uplo = 'L')
  !          A += tril(alpha*x*y' + alpha*y*x'), where A is lower triangluar
  !          i.e. sym(A) := alpha*(x*y' + y*x') + sym(A)
  !          (with uplo = 'U')
  !          A += triu(alpha*x*y' + alpha*y*x'), where A is upper triangluar
  !          i.e. sym(A') := alpha*(x*y' + y*x') + sym(A')
  !    RMD:  (with uplo = 'L' or 'U':)
  !          xa += alpha*(Aa + Aa')*y (eqiv.to: xa += alpha*(diag(Aa)+sym(Aa))*y)
  !          ya += alpha*(Aa + Aa')*x (eqiv.to: ya += alpha*(diag(Aa)+sym(Aa))*x)
  !          Aa unchanged

  ! Local variables
  logical selx, sely

  selx = sel(1:1) == '1'
  sely = sel(2:2) == '1'

  if(selx) then
    call sgbmv('n', n, n, 0, 0, alpha, Aa, lda+1, y, incy, 1.0, xa, incx)
    call ssymv(uplo, n, alpha, Aa, lda, y, incy, 1.0, xa, incx)
  end if
  if(sely) then
    call sgbmv('n', n, n, 0, 0, alpha, Aa, lda+1, x, incx, 1.0, ya, incy)
    call ssymv(uplo, n, alpha, Aa, lda, x, incx, 1.0, ya, incy)
  end if
END SUBROUTINE SSYR2_RMD
