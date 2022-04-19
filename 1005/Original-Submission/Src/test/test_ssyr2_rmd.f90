PROGRAM TEST_SSYR2_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer :: n, esa, incx, incy, dimu, dimv
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  n = 2
  dimv = n*(n+1)/2
  dimu = n*(n+1)/2 + 2*n + 1
  testcase = fdata([real::], [n, 0, 1, 1, dimu], ['U'])
  call rmd_stest123(F, F_rmd, dimu, dimv, tol, testcase)
  !
  do n = 1,5
    do esa = 0, 3, 3
      do incx = 1, 3
        do incy = 1, 3
          dimv = n*(n+1)/2
          dimu = n*(n+1)/2 + 2*n + 1
          testcase = fdata([real::], [n, esa, incx, incy, dimu], ['L'])
          call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
          testcase % characters(1) = 'U'
          call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
        end do
      end do
    end do
  end do
END PROGRAM TEST_SSYR2_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: n, k, esa, incx, incy, dimu
  real :: alpha
  character uplo
  !  
  real, allocatable :: A(:, :), x(:), y(:)
  !
  uplo  = testcase % characters(1)
  n     = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  incy = testcase % integers(4)
  dimu = testcase % integers(5)
  !
  allocate(A(n+esa, n), x(n*incx), y(n*incy))
  !
  ! Copy from u to upper or lower triangle of A
  ! depending on the value of uplo
  call rmd_strishape_v2m(uplo, n, u, A, n+esa)
  !
  k = (n*(n+1))/2
  call scopy(n, u(k+1), 1, x, incx)
  call scopy(n, u(k+n+1), 1, y, incy)
  alpha = u(dimu)
  !
  ! Actual blas operation
  ! A = alpha*(x*y' + y*x') + A
  call ssyr2(uplo, n, alpha, x, incx, y, incy, A, n+esa)   
  !
  ! Copy back from A to v
  call rmd_strishape_m2v(uplo, n, v, A, n+esa)
  !
  ! call disp('A = ', A)
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! 
  type(fdata), intent(inout) :: testcase
  real, intent(in)     :: u(*), va(*), v(*)
  real, intent(inout)  :: ua(*)
  !
  real, allocatable :: x(:), y(:), Aa(:, :), xa(:), ya(:)
  integer :: n, k, esa, incx, incy, dimu
  real :: alpha, alphaa
  character uplo
  !
  uplo  = testcase % characters(1)
  n     = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  incy = testcase % integers(4)
  dimu = testcase % integers(5)
  !
  allocate(x(n*incx), y(n*incy), Aa(n+esa, n), xa(n*incx), ya(n*incy))
  !
  call rmd_strishape_v2m(uplo, n, va, Aa, n+esa)
  !
  k = (n*(n+1))/2
  !  
  call scopy(n, u(k+1), 1, x, incx)
  call scopy(n, u(k+n+1), 1, y, incy)
  call scopy(n, ua(k+1), 1, xa, incx)
  call scopy(n, ua(k+n+1), 1, ya, incy)
  alpha = u(dimu)
  alphaa = ua(dimu)
  !
  call ssyr2_rmds(uplo, n, x, incx, y, incy, n+esa, alphaa, Aa)
  call ssyr2_rmd(uplo, n, alpha, x, incx, y, incy, n+esa, xa, ya, Aa, '111')
  !  
  call rmd_strishape_m2v(uplo, n, ua, Aa, n+esa)
  call scopy(n, xa, incx, ua(k+1), 1)
  call scopy(n, ya, incy, ua(k+n+1), 1)
  ua(dimu) = alphaa
  !
END SUBROUTINE F_RMD
