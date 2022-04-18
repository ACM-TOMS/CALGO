PROGRAM TEST_SSPR2_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer :: n, dimu, dimv, incx, incy
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  do n = 1, 5
    do incx = 1, 3, 2
      do incy = 1, 3, 2
        dimu  = n*(n+1)/2 + 2*n + 1
        dimv = n*(n+1)/2
        testcase = fdata([real::], [n, incx, incy, dimu, dimv], ['L'])
        call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
        testcase % characters(1) = 'U'
        call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
      end do
    end do
  end do
END PROGRAM TEST_SSPR2_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: n, incx, incy, dimu, dimv
  real :: alpha
  character uplo
  real, allocatable :: x(:), y(:), AP(:)
  !
  uplo   = testcase % characters(1)
  n      = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  dimu = testcase % integers(4)
  dimv = testcase % integers(5)
  alpha = u(dimu)
  !
  allocate(x(n*incx), y(n*incy), AP(dimv))
  !
  call scopy(dimv, u, 1, AP, 1)
  call scopy(n, u(dimv + 1), 1, x, incx)
  call scopy(n, u(dimv + n + 1), 1, y, incy)
  !
  ! Actual blas operation
  call sspr2(uplo, n, alpha, x, incx, y, incy, AP) ! AP = alpha*(x*y' + y*x') + AP
  call scopy(dimv, AP, 1, v, 1)
  !  
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)        :: u(*), va(*), v(*)
  real, intent(inout)     :: ua(*)
  !
  ! Local variables
  real, allocatable :: x(:), y(:), xa(:), ya(:), AP(:), APa(:)
  integer :: n, incx, incy, dimu, dimv
  real :: alpha, alphaa
  character uplo
  !
  uplo   = testcase % characters(1)
  n      = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  dimu = testcase % integers(4)
  dimv = testcase % integers(5)
  !
  allocate(x(n*incx), y(n*incy), xa(n*incx), ya(n*incy), AP(dimv), APa(dimv))
  !
  call scopy(n, u(dimv+1), 1, x, incx)
  call scopy(n, u(dimv+n+1), 1, y, incy)
  call scopy(n, ua(dimv+1), 1, xa, incx)
  call scopy(n, ua(dimv+n+1), 1, ya, incy)
  call scopy(dimv, va, 1, APa, 1)
  alpha = u(dimu)
  alphaa = ua(dimu)
  !
  call sspr2_rmds(uplo, n, x, incx, y, incy, alphaa, APa)
  call sspr2_rmd(uplo, n, alpha, x, incx, y, incy, xa, ya, APa, '111')
  !  
  call scopy(dimv, APa, 1, ua, 1)
  call scopy(n, xa, incx, ua(dimv+1), 1)
  call scopy(n, ya, incy, ua(dimv+n+1), 1)
  ua(dimu) = alphaa
  !
END SUBROUTINE F_RMD
