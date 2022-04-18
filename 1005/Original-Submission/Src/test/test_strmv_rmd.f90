PROGRAM TEST_STRMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer     :: n, sizeIn, sizeOut, u, t, d
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1, esa, incx
  !
  character :: uplo(2), trans(2), diag(2)
  !
  uplo  = ['L', 'U']
  trans = ['N', 'T']
  diag  = ['N', 'U']
  !
  do n = 1, 5
    do u = 1, 2
    do t = 1, 2
    do d = 1, 2
      do esa = 0, 3, 3
      do incx = 1, 3, 2
        sizeIn  = n + (n*(n+1))/2
        sizeOut = n
        testcase = fdata([real::], [n, esa, incx], [uplo(u), trans(t), diag(d)])
        call rmd_stestrandom(F, F_rmd, sizeIn, sizeOut, tol, rep, testcase)
      end do
      end do
    end do
    end do
    end do
  end do
END PROGRAM TEST_STRMV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: n, sizeA, esa, incx
  character uplo, trans, diag
  !  
  real, allocatable :: A(:, :), x(:)
  !
  uplo  = testcase % characters(1)
  trans = testcase % characters(2)
  diag  = testcase % characters(3)
  n     = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  !
  sizeA = (n*(n+1))/2
  !  
  allocate(A(n+esa, n), x(n*incx))
  !
  call rmd_strishape_v2m(uplo, n, u, A, n+esa) ! A
  call scopy(n, u(sizeA+1), 1, x, incx)              ! x
  !
  ! BLAS operation
  call strmv(uplo, trans, diag, n, A, n+esa, x, incx)
  !  
  ! Copy back from x to v
  call scopy(n, x, incx, v, 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)     :: u(*), va(*), v(*)
  real, intent(inout)  :: ua(*)
  !  
  ! Local variables
  real, allocatable :: A(:, :), x0(:), Aa(:, :), xa(:)
  integer :: n, sizeA, esa, incx
  character uplo, trans, diag
  !
  uplo   = testcase % characters(1)
  trans  = testcase % characters(2)
  diag   = testcase % characters(3)
  n      = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  !  
  allocate(A(n+esa, n), x0(n*incx), Aa(n+esa, n), xa(n*incx))
  !
  sizeA = (n*(n+1))/2
  !
  ! Copy inital values to local variables
  call rmd_strishape_v2m(uplo, n, u, A, n+esa) ! A
  call rmd_strishape_v2m(uplo, n, ua, Aa, n+esa) ! Aa
  call scopy(n, u(sizeA+1), 1, x0, incx)          ! x0
  call scopy(n, va, 1, xa, incx)           ! xa
  !
  ! RMD operation
  call strmv_rmd(uplo, trans, diag, n, A, n+esa, x0, incx, Aa, xa, '11')
  !  
  ! Copy result to output variables
  call rmd_strishape_m2v(uplo, n, ua, Aa, n+esa) ! Aa to ua
  call scopy(n, xa, incx, ua(sizeA+1), 1)           ! xa to ua(sizeA+1)
  !  
END SUBROUTINE F_RMD
