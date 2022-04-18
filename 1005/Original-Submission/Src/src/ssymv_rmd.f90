SUBROUTINE SSYMV_RMD(uplo, n, alpha, A, lda, x, incx, beta, incy, Aa, xa, ya, sel)
  implicit none
  character :: uplo
  integer :: n, lda, incx, incy
  real :: alpha, beta
  real :: A(lda, *), x(*), Aa(lda, *), xa(*), ya(*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSYMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSYMV was called with the arguments
  !
  !       uplo, n, alpha, A, lda, x, incx, beta, y, incy
  !
  !    then the corresponding call to SSYMV_RMD should begin with the same
  !    arguments
  !
  !       uplo, n, alpha, A, lda, x, incx, beta, incy
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that y is omitted. In addition the following arguments should be
  !    provided:
  !
  !    Aa
  !       (input, output, real triangular matrix of the same dimensions as A,
  !       and stored in the same half according to uplo)
  !       Aa += the adjoint of A due to the SSYMV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SSYMV call
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       On entry: the adjoint of the y produced by SSYMV
  !       On exit: the adjoint of the y supplied to SSYMV
  !
  !    sel
  !       (input, character*3)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be updated, else sel(2:2) = '0'
  !          sel(3:3) = '1' if ya should be computed, else sel(3:3) = '0'
  !       For example, to update only Aa, set sel = '100'.
  !       
  ! OPERATIONS
  !    (with uplo = 'L')
  !    BLAS: y = alpha*sym(A)*x + beta*y                  A lower triangular
  !    RMD:  Aa += alpha*(tril(x*ya'+ya*x') - diag(x*ya'))  where ya is value on entry
  !          xa += alpha*sym(A)*ya                        do. 
  !          ya := beta*ya
  !    (with uplo = 'U')
  !    BLAS: y = alpha*sym(A')*x + beta*y                 A upper triangular
  !    RMD:  Aa += alpha*(triu(x*ya'+ya*x') - diag(x*ya'))  where ya is value on entry
  !          xa += alpha*sym(A')*ya                       do. 
  !          ya := beta*ya
  
  ! Local variables
  logical :: sela, selx, sely
  integer :: i
  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  sely = sel(3:3) == '1'

  if(sela) then
    ! Update Aa
    call ssyr2(uplo, n, alpha, x, incx, ya, incy, Aa, lda)
    do i=1,n
      Aa(i,i) = Aa(i,i) - alpha*x(1 + (i-1)*incx)*ya(1 + (i-1)*incy)
    enddo
  end if

  if(selx) then
    ! Update xa
    call ssymv(uplo, n, alpha, A, lda, ya, incy, 1.0, xa, incx)
  end if

  if(sely) then 
    ! Update ya
    call sscal(n, beta, ya, incy)
  end if
END SUBROUTINE SSYMV_RMD

