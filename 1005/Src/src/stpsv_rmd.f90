SUBROUTINE STPSV_RMD(uplo, trans, diag, n, AP, x, incx, APa, xa, wrk, sel)
  implicit none
  character :: uplo, trans, diag
  integer :: n, incx
  real :: AP(*), x(*), APa(*), xa(*), wrk(*)
  character(2) :: sel
  
  ! PURPOSE
  !    Calculate the reverse mode derivative of STPSV from BLAS.
  !
  ! ARGUMENTS
  !    If STPSV was called with the following arguments:
  !
  !       uplo, trans, diag, n, AP, x, incx
  !
  !    then the corresponding call to STPSV_RMD should begin with the same arguments
  !
  !       uplo, trans, diag, n, AP, x, incx
  !
  !    with the same values. All these arguments will remain unchanged on exit.
  !    In addition the following arguments should be provided:
  !
  !    APa
  !       (input, output, real triangular matrix of the same dimensions as AP,
  !       stored in a vector with (n*(n+1))/2 elements in the same way as AP)
  !       APa += adjoint of AP due to the STPSV call
  !
  !    xa
  !       (input, output, real vector of the same dimension and increment as x)
  !       On entry: The adjoint of the x produced by STPSV
  !       On exit: The adjoint of the x supplied to STPSV
  !
  !    wrk
  !       (output, real vector of dimension at least n)
  !       When sel(2:2) = '0' so that a new xa should not be computed it is
  !       necessary to supply STPSV_RMD with a workspace vector. When sel(2:2) =
  !       '1', wrk is not referenced, because xa serves its purpose. In this case
  !       a dummy value may be given instead
  !
  !    sel
  !       (input, character*2)
  !       Used to select which adjoints to update/compute:
  !          sel(1:1) = '1' if APa should be updated, else sel(1:1) = '0'
  !          sel(2:2) = '1' if xa should be computed, else sel(2:2) = '0'
  !       For example, if sel = '01', then only xa will be computed
  !
  ! OPERATIONS
  !    Same as for STRSV_RMD except that packed storage is used

  ! Local variables
  integer :: i, k, d
  logical sela, selx, ww, wx
  character :: ntrans, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selx = sel(2:2) == '1'
  wx = selx
  ww = .not.selx
  
  ntrans = trno(1 + (index('NnTt', trans)-1)/2) ! N-->T, T-->N
  
  if(diag == 'n' .or. diag == 'N') then
    d = 0
  else
    d = 1
  end if

  if (selx) then
    call stpsv(uplo, ntrans, diag, n, AP, xa, incx)
  else
    call scopy(n, xa, incx, wrk, 1)
    call stpsv(uplo, ntrans, diag, n, AP, wrk, 1)
  endif
  if(sela) then
    ! Calculate APa
    if(trans == 'n' .or. trans == 'N') then
      if(uplo == 'u' .or. uplo == 'U') then
        ! AP is upper triangular
        k = 1 + d
        do i = 1+d,n
          if (wx) call saxpy(i-d, -x(1+(i-1)*incx), xa, incx, APa(k), 1)
          if (ww) call saxpy(i-d, -x(1+(i-1)*incx), wrk, 1, APa(k), 1)
          k = k + i
        end do
      else
        ! A is lower triangular
        k = 1
        do i = 1,n-d
          if (wx) call saxpy(n-i+1-d, -x(1+(i-1)*incx), xa(1+(i+d-1)*incx), incx, APa(k +d*i), 1)
          if (ww) call saxpy(n-i+1-d, -x(1+(i-1)*incx), wrk(i+d)          , 1   , APa(k +d*i), 1)
          k = k + n-i+1-d
        end do
      end if
    else
      ! BLAS operation
      !   x := AP^(-T)*x
      ! RMD operation
      !   APa -= tril(x*xa')
      if(uplo == 'u' .or. uplo == 'U') then
        ! AP is upper triangular
        k = 1 + d
        do i = 1+d,n
          if (wx) call saxpy(i-d, -xa(1+(i-1)*incx), x, incx, APa(k), 1)
          if (ww) call saxpy(i-d, -wrk(i)          , x, incx, APa(k), 1)
          k = k + i
        end do
      else
        ! AP is lower triangular
        k = 1
        do i = 1,n-d
          if (wx) call saxpy(n-i+1-d, -xa(1+(i-1)*incx), x(1+(i+d-1)*incx), incx, APa(k+d*i), 1)
          if (ww) call saxpy(n-i+1-d, -wrk(i)          , x(1+(i+d-1)*incx), incx, APa(k+d*i), 1)
          k = k + n-i+1-d
        end do
      end if
    end if
  end if
  
END SUBROUTINE STPSV_RMD
