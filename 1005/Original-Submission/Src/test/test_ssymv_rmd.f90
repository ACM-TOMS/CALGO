PROGRAM TEST_SSYMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer :: n, esa, incx, incy, dimu
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  n = 1
  dimu = n*(n+1)/2 + 2*n + 2
  testcase = fdata([real::], [n, 0, 1, 1, dimu], ['L'])
  call rmd_stest123(F, F_rmd, dimu, n, tol, testcase)
  do n = 1, 5
    do esa = 0, 3, 3
      do incx = 1, 3, 2
        do incy = 1, 3, 2
          dimu = n*(n+1)/2 + 2*n + 2
          testcase = fdata([real::], [n, esa, incx, incy, dimu], ['L'])
          call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
          testcase % characters(1) = 'U'
          call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
        end do
      end do
    end do
  end do
END PROGRAM TEST_SSYMV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  ! 
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: n, k, esa, incx, incy, dimu
  real :: alpha, beta
  character uplo
  !  
  real, allocatable :: A(:, :), x(:), y(:)
  !
  uplo  = testcase % characters(1)
  n     = testcase % integers(1)
  esa  = testcase % integers(2)
  incx = testcase % integers(3)
  incy = testcase % integers(4)
  dimu   = testcase % integers(5)
  !
  allocate(A(n+esa, n), x(n*incx), y(n*incy))
  !  
  ! Copy from u to upper or lower triangle of A
  ! depending on the value of uplo
  call rmd_strishape_v2m(uplo, n, u, A, n+esa)
  !
  k = n*(n+1)/2
  call scopy(n, u(k+1), 1, x, incx)
  call scopy(n, u(k+n+1), 1, y, incy)
  alpha = u(dimu - 1)
  beta = u(dimu)

  !
  ! Actual blas operation
  ! y = alpha*A*x + beta*y
  call ssymv(uplo, n, alpha, A, n+esa, x, incx, beta, y, incy)
  !  
  ! Copy back from y to v
  call scopy(n, y, incy, v, 1)
  !  
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in)    :: u(*), va(*), v(*)
  real, intent(inout) :: ua(*)
  !
  real, allocatable :: A(:, :), x(:), y0(:), Aa(:, :), xa(:), ya(:)
  integer :: n, k, esa, incx, incy, lda, dimu
  real :: alpha, beta, alphaa, betaa
  character uplo
  !
  uplo  = testcase % characters(1)
  n     = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  incy = testcase % integers(4)
  !
  allocate(A(n+esa, n), x(n*incx), Aa(n+esa, n), xa(n*incx), ya(n*incy))
  allocate(y0(n*incy))
  !
  call rmd_strishape_v2m(uplo, n, u, A, n+esa)
  call rmd_strishape_v2m(uplo, n, ua, Aa, n+esa)
  !
  k = n*(n+1)/2
  !  
  call scopy(n, u(k+1), 1, x, incx)
  call scopy(n, ua(k+1), 1, xa, incx)
  call scopy(n, u(k+n+1), 1, y0, incy)
  call scopy(n, va, 1, ya, incy)
  dimu   = testcase % integers(5)
  alpha = u(dimu-1)
  beta = u(dimu)
  alphaa = ua(dimu-1)
  betaa = ua(dimu)
  !
  lda = n + esa
  call ssymv_rmds(uplo, n, A, lda, x, incx, y0, incy, alphaa, betaa, ya, '11')
  call ssymv_rmd(uplo, n, alpha, A, lda, x, incx, beta, incy, Aa, xa, ya, '111')
  !
  call rmd_strishape_m2v(uplo, n, ua, Aa, n+esa)
  call scopy(n, xa, incx, ua(k+1), 1)
  call scopy(n, ya, incy, ua(k+n+1), 1)
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !  
END SUBROUTINE F_RMD
