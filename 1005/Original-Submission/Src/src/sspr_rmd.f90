SUBROUTINE SSPR_RMD(uplo, n, alpha, x, incx, xa, APa)
  implicit none
  integer :: n, incx
  real :: alpha, x(*), xa(*), APa(*)
  character :: uplo

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSPR from BLAS. 
  !
  !
  ! ARGUMENTS
  !    If SSPR was called with the arguments
  !
  !       uplo, n, alpha, x, incx, AP
  !    
  !    then the corresponding call to SSPR_RMD should begin with the arguments
  !
  !       uplo, n, alpha, x, incx
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that AP is omitted. In addition the following arguments should be
  !    provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SSPR operation
  !
  !    APa
  !       (input, real packed triangular matrix stored in a vector with
  !       n*(n+1)/2 elements in the same way as AP)
  !       The adjoint of AP
  !
  ! NOTE
  !    A sel parameter is not needed because APa is unchanged.
  !
  ! OPERATIONS
  !    The same as for SSYR_RMD except that packed storage is used.

  ! Local variables
  integer :: i, ix, k

  ! xa += alpha*diag(APa)*x
  k  = 1
  ix = 1
  if(uplo == 'u' .or. uplo == 'U') then
    do i = 1,n
      xa(ix) = xa(ix) + alpha*APa(k)*x(ix)
      k = k+i+1
      ix = ix + incx
    end do
  else
    do i = 1,n
      xa(ix) = xa(ix) + alpha*APa(k)*x(ix)
      k = k+n-i+1
      ix = ix + incx
    end do
  end if
      
  ! xa += alpha*sym(APa)*x
  call sspmv(uplo, n, alpha, APa, x, incx, 1.0, xa, incx)
END SUBROUTINE SSPR_RMD
