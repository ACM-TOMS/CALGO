PROGRAM TEST_SSPMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  type(fdata) :: testcase
  integer :: n, dimu, incx, incy
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  n = 1
  dimu = n*(n+1)/2 + 2*n + 2
  testcase = fdata([real::], [n, 1, 1, dimu], ['L'])
  call rmd_stest123(F, F_rmd, dimu, n, tol, testcase)
  do n = 1, 5
    do incx = 1, 3, 2
      do incy = 1, 3, 2
        dimu = n*(n+1)/2 + 2*n + 2
        testcase = fdata([real::], [n, incx, incy, dimu], ['L'])
        call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
        testcase = fdata([real::], [n, incx, incy, dimu], ['U'])
        call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
      end do
    end do
  end do
END PROGRAM TEST_SSPMV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  ! 
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: n, sizeAP, incx, incy, dimu
  real :: alpha, beta
  character uplo
  !  
  real, allocatable :: AP(:), x(:), y(:)
  !
  uplo   = testcase % characters(1)
  n      = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  dimu = testcase % integers(4)
  sizeAP = n*(n+1)/2
  !  
  allocate(AP(sizeAP), x(n*incx), y(n*incy))
  !  
  call scopy(sizeAP, u, 1, AP, 1)
  call scopy(n, u(sizeAP+1), 1, x, incx)
  call scopy(n, u(sizeAP+n+1), 1, y, incy)
  alpha = u(dimu-1)
  beta = u(dimu)
  !
  ! Actual BLAS operation
  ! y = alpha*AP*x + beta*y
  call sspmv(uplo, n, alpha, AP, x, incx, beta, y, incy) 
  !  
  ! Copy back from y to v
  call scopy(n, y, incy, v, 1)
  !  
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in)        :: u(*), va(*), v(*)
  real, intent(inout)     :: ua(*)
  !
  real, allocatable :: AP(:), x(:), APa(:), xa(:), ya(:), y0(:)
  integer :: n, sizeAP, incx, incy, dimu
  real :: alpha, beta, alphaa, betaa
  character uplo
  !
  uplo   = testcase % characters(1)
  n      = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  sizeAP = (n*(n+1))/2
  !
  allocate(AP(sizeAP), x(n*incx), APa(sizeAP), xa(n*incx), ya(n*incy), y0(n*incy))
  !
  call scopy(sizeAP, u, 1, AP, 1)
  call scopy(sizeAP, ua, 1, APa, 1)    
  call scopy(n, u(sizeAP+1), 1, x, incx)
  call scopy(n, ua(sizeAP+1), 1, xa, incx)
  call scopy(n, va, 1, ya, incy)
  call scopy(n, u(sizeAP+n+1), 1, y0, incy)
  dimu = testcase % integers(4)
  alpha = u(dimu-1)
  beta = u(dimu)
  alphaa = ua(dimu-1)
  betaa = ua(dimu)
  !
  call sspmv_rmds(uplo, n, AP, x, incx, y0, incy, alphaa, betaa, ya, '11')
  call sspmv_rmd(uplo, n, alpha, AP, x, incx, beta, incy, APa, xa, ya, '111')
  !
  call scopy(sizeAP, APa, 1, ua, 1)
  call scopy(n, xa, incx, ua(sizeAP+1), 1)
  call scopy(n, ya, incy, ua(sizeAP+n+1), 1)
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !  
END SUBROUTINE F_RMD
