SUBROUTINE STPMV_RMD(uplo, trans, diag, n, AP, x0, incx, APa, xa, sel)
  implicit none
  character :: uplo, trans, diag
  integer :: n, incx
  real :: AP(*), x0(*), APa(*), xa(*)
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the reverse mode derivative of STPMV from BLAS.
  !
  !
  ! ARGUMENTS
  !    If STPMV was called with the following arguments:
  !
  !       uplo, trans, diag, n, AP, x, incx
  !
  !    then the corresponding call to STPMV_RMD should begin with the arguments
  !
  !       uplo, trans, diag, n, AP, x0, incx
  !
  !    which all except x0 should have the same values as they had on the STRMV
  !    call, and x0 should have the value that x had on entry to the STPMV-call
  !    (STPMV only changes the x-argument). All these arguments will remain
  !    unchanged on exit. In addition the following arguments should be
  !    provided:
  !
  !    APa
  !       (input, output, real triangular matrix of the same dimensions as AP,
  !       stored in a vector with (n*(n+1))/2 elements in the same way as AP)
  !       APa += the adjoint of AP due to the STPMV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       On entry: the adjoint of the x produced by STPMV
  !       On exit: the adjoint of the x supplied to STPMV
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if APa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be computed, else sel(2:2) = '0'
  !       For example, if sel = '01', then only xa will be computed.
  !       
  ! OPERATIONS
  !    Same as for for STRMV_RMD except that packed storage is used

  ! Local variables
  integer :: i, k, d
  logical sela, selx
  character :: ntrans, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  ntrans = trno(1 + (index('NnTt', trans)-1)/2) ! N-->T, T-->N
  
  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  if(sela) then
    ! Calculate APa
    if(trans == 'n' .or. trans == 'N') then
      ! BLAS operation
      !   x = AP*x0  where AP is triangular packed
      ! RMD operation
      !   APa += xa*x0'
      if(uplo == 'u' .or. uplo == 'U') then
        ! AP is upper triangular
        k = 1 + d
        do i = 1+d,n
          call saxpy(i-d, x0(1+(i-1)*incx), xa, incx, APa(k), 1)
          k = k + i
        end do
      else
        ! A is lower triangular
        k = 1
        do i = 1,n-d
          call saxpy(n-i+1-d, x0(1+(i-1)*incx), xa(1+(i+d-1)*incx), incx, &
              APa(k+d*i), 1)
          k = k + n-i+1-d
        end do
      end if
    else
      ! BLAS operation
      !   x = A'*x  where A is triangular packed
      ! RMD operation
      !   Aa += x*xa'
      if(uplo == 'u' .or. uplo == 'U') then
        ! A is upper triangular
        k = 1 + d
        do i = 1+d,n
          call saxpy(i-d, xa(1+(i-1)*incx), x0, incx, APa(k), 1)
          k = k + i
        end do
      else
        ! A is lower triangular
        k = 1
        do i = 1,n-d
          call saxpy(n-i+1-d, xa(1+(i-1)*incx), x0(1+(i+d-1)*incx), incx, &
              APa(k+d*i), 1)
          k = k + n-i+1-d
        end do
      end if
    end if
  end if

  if(selx) then
    call stpmv(uplo, ntrans, diag, n, AP, xa, incx)
  end if
  
END SUBROUTINE STPMV_RMD
