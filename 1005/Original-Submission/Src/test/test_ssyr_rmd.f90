PROGRAM TEST_SSYR_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer :: n, i, esa, incx, dimu
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  do n = 1, 5
    do esa = 0, 3, 3
    do incx = 1, 3, 2
      i = (n*(n+1))/2
      dimu = i+n+1
      testcase = fdata([real::], [n, esa, incx, dimu], ['L'])
      call rmd_stestrandom(F, F_rmd, dimu, i, tol, rep, testcase)
      testcase % characters(1) = 'U'
      call rmd_stestrandom(F, F_rmd, dimu, i, tol, rep, testcase)
    end do 
    end do 
  end do
END PROGRAM TEST_SSYR_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  !  
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: n, esa, incx, dimu
  real :: alpha
  character uplo
  !  
  real, allocatable :: A(:, :), x(:)
  !
  uplo  = testcase % characters(1)
  n     = testcase % integers(1)
  esa   = testcase % integers(2)
  incx  = testcase % integers(3)
  dimu  = testcase % integers(4)
  !
  allocate(A(n+esa, n), x(n*incx))
  !
  ! Copy from u to upper or lower triangle of A
  ! depending on the value of uplo
    call rmd_strishape_v2m(uplo, n, u, A, n+esa)
  ! Copy from u to x
    call scopy(n, u(n*(n+1)/2+1), 1, x, incx)
    alpha = u(dimu)
  !
  ! Actual blas operation
  call ssyr(uplo, n, alpha, x, incx, A, n+esa) ! A = alpha*x*x' + A
  !
  ! Copy back from A to v
    call rmd_strishape_m2v(uplo, n, v, A, n+esa)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in)     :: u(*), va(*), v(*)
  real, intent(inout)  :: ua(*)
  !
  real, allocatable :: x(:), Aa(:, :), xa(:)
  integer :: n, k, esa, incx, dimu
  real :: alpha, alphaa
  character uplo
  !
  uplo  = testcase % characters(1)
  n     = testcase % integers(1)
  esa   = testcase % integers(2)
  incx  = testcase % integers(3)
  dimu  = testcase % integers(4)
  !
  allocate(x(n*incx), Aa(n+esa, n), xa(n*incx))
  !
  call rmd_strishape_v2m(uplo, n, va, Aa, n+esa)
  !
  k = (n*(n+1))/2
  !  
  call scopy(n, u(k+1), 1, x, incx)
  call scopy(n, ua(k+1), 1, xa, incx)
  alpha = u(dimu)
  alphaa = ua(dimu)
  !
  call ssyr_rmds(uplo, n, x, incx, n+esa, alphaa, Aa)
  call ssyr_rmd(uplo, n, alpha, x, incx, n+esa, xa, Aa)
  !  
  call scopy(k, va, 1, ua, 1)
  call scopy(n, xa, incx, ua(k+1), 1)
  ua(dimu) = alphaa
  !
END SUBROUTINE F_RMD
