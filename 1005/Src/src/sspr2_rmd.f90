SUBROUTINE SSPR2_RMD(uplo, n, alpha, x, incx, y, incy, xa, ya, APa, sel)
  implicit none
  integer :: n, incx, incy
  real :: alpha, x(*), y(*), xa(*), ya(*), APa(*)
  character :: uplo
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the reverse mode derivative of SSPR2 from BLAS. 
  !
  !
  ! ARGUMENTS
  !    If SSPR2 was called with the arguments
  !
  !       uplo, n, alpha, x, incx, y, incy, AP
  !    
  !    then the corresponding call to SSPR2_RMD should begin with the arguments
  !
  !       uplo, n, alpha, x, incx, y, incy
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that AP is omitted. In addition the following arguments should be
  !    provided:
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += the adjoint of x due to the SSPR2 call.
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       ya += the adjoint of y due to the SSPR2 call.
  !
  !    APa
  !       (input, real packed triangular matrix of the same dimensions as AP,
  !       and stored in a vector with n*(n+1)/2 elements in the same way as AP)
  !       The adjoint of AP.
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update:
  !          sel(1:1) = '1' if xa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if ya should be updated, else sel(2:2) = '0'
  !       For example, to update only xa, set sel = '10'.
  !
  ! OPERATIONS
  !    The same as for SSYR2_RMD except that packed storage is used
  
  ! Local variables
  integer :: i, ix, iy, k
  logical selx, sely

  selx = sel(1:1) == '1'
  sely = sel(2:2) == '1'

  if(selx) then
    k  = 1
    ix = 1
    iy = 1
    ! xa += alpha*diag(APa)*y
    if(uplo == 'u' .or. uplo == 'U') then
      do i = 1,n
        xa(ix) = xa(ix) + alpha*APa(k)*y(iy)
        k = k+i+1
        ix = ix + incx
        iy = iy + incy
      end do
    else
      do i = 1,n
        xa(ix) = xa(ix) + alpha*APa(k)*y(iy)
        k = k+n-i+1
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    ! xa += alpha*sym(APa)*y
    call sspmv(uplo, n, alpha, APa, y, incy, 1.0, xa, incx)
  end if
  
  if(sely) then
    k  = 1
    ix = 1
    iy = 1
    ! ya += alpha*diag(APa)*x
    if(uplo == 'u' .or. uplo == 'U') then
      do i = 1,n
        ya(iy) = ya(iy) + alpha*APa(k)*x(ix)
        k = k+i+1
        ix = ix + incx
        iy = iy + incy
      end do
    else
      do i = 1,n
        ya(iy) = ya(iy) + alpha*APa(k)*x(ix)
        k = k+n-i+1
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    ! ya += alpha*sym(APa)*x
    call sspmv(uplo, n, alpha, APa, x, incx, 1.0, ya, incy)
  end if

END SUBROUTINE SSPR2_RMD
