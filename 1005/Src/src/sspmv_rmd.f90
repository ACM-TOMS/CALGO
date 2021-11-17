SUBROUTINE SSPMV_RMD(uplo, n, alpha, AP, x, incx, beta, incy, APa, xa, ya, sel)
  implicit none
  character :: uplo
  integer :: n, incx, incy
  real :: alpha, beta
  real :: AP(*), x(*), APa(*), xa(*), ya(*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SSPMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If SSPMV was called with the arguments
  !
  !       uplo, n, alpha, AP, x, incx, beta, y, incy
  !
  !    then the corresponding SSPMV_RMD call should begin with the arguments
  !
  !       uplo, n, alpha, AP, x, incx, beta, incy
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    Note that y is omitted. In addition the following arguments should be
  !    provided:
  !
  !    APa
  !       (input, output, real triangular packed matrix stored in a vector with
  !       (n*(n+1))/2 elements in the same way as AP)
  !       APa += adjoint of AP due to the SSPMV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       xa += adjoint of x due to the SSPMV call
  !
  !    ya
  !       (input, output, real vector of the same dimension and increment as y)
  !       On entry: the adjoint of the y produced by SSPMV
  !       On exit: the adjoint of the y supplied to SSPMV
  !
  !    sel
  !       (input, character*3)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if APa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be updated, else sel(2:2) = '0'
  !          sel(3:3) = '1' if ya should be computed, else sel(3:3) = '0'
  !       For example, to update only APa, set sel = '100'.
  !       
  ! OPERATIONS
  !    Same as for SSYMV_RMD except that A and Aa use packed storage

  ! Local variables
  integer :: i, k
  logical sela, selx, sely
  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  sely = sel(3:3) == '1'

  if(sela) then
    ! Update APa
    call sspr2(uplo, n, alpha, x, incx, ya, incy, APa)
    if(uplo == 'u' .or. uplo == 'U') then
      k = 0
      do i = 1,n
        k = k+i
        APa(k) = APa(k) - alpha*x(1 + (i-1)*incx)*ya(1 + (i-1)*incy)
      end do
    else
      k = 1
      do i = 1,n
        APa(k) = APa(k) - alpha*x(1 + (i-1)*incx)*ya(1 + (i-1)*incy)
        k = k+n-i+1
      end do
    end if
  end if
    
  if(selx) then
    ! Update xa
    call sspmv(uplo, n, alpha, AP, ya, incy, 1.0, xa, incx)
  end if

  if(sely) then
    ! Update ya
    call sscal(n, beta, ya, incy)
  end if
END SUBROUTINE SSPMV_RMD

