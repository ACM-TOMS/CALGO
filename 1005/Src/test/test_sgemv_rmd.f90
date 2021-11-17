PROGRAM TEST_SGEMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  integer :: m, n, nn, rep, esa, incx, incy
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  rep = 1
  m = 1
  n = 1
  nn = m*n + m + n + 2
  testcase = fdata([real::], [m, n, 0, 1, 1, nn], ['n'])
  call rmd_stest123(F, F_rmd, nn, m, tol, testcase)
  do m = 1, 3
    do n = 1, 3
      do esa = 0, 3, 3
        do incx = 1, 3, 2
          do incy = 1, 3, 2
            nn = m*n + m + n + 2
            testcase = fdata([real::], [m, n, esa, incx, incy, nn], ['n'])
            call rmd_stestrandom(F, F_rmd, nn, m, tol, rep, testcase)
            testcase = fdata([real::], [m, n, esa, incx, incy, nn], ['t'])
            call rmd_stestrandom(F, F_rmd, nn, n, tol, rep, testcase)
          end do
        end do
      end do
    end do
  end do
END PROGRAM

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), x(:), y(:)
  integer :: m, n, nx, ny, esa, incx, incy, nn
  real :: alpha, beta
  character trans
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  trans = testcase % characters(1)
  esa = testcase % integers(3)
  incx = testcase % integers(4)
  incy = testcase % integers(5)
  nn = testcase % integers(6)
  !
  if(trans == 'n' .or. trans == 'N') then
    nx = n
    ny = m
  else
    nx = m
    ny = n
  end if
  !
  allocate(A(m+esa, n), x(nx*incx), y(ny*incy))
  !
  A(1:m, 1:n) = reshape(u(1:m*n), [m, n])
  call scopy(nx, u(m*n+1), 1, x, incx)
  call scopy(ny, u(m*n+nx+1), 1, y, incy)
  alpha = u(nn - 1)
  beta = u(nn)
  !
  call sgemv(trans, m, n, alpha, A, m+esa, x, incx, beta, y, incy)
  !
  call scopy(ny, y, incy, v, 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  !use dispmodule
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), x(:), Aa(:, :), xa(:), ya(:), y0(:)
  integer :: m, n, nx, ny, esa, incx, incy, nn
  real :: alpha, beta, alphaa, betaa
  character trans
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  trans = testcase % characters(1)
  esa = testcase % integers(3)
  incx = testcase % integers(4)
  incy = testcase % integers(5)
  nn = testcase % integers(6)
  !
  if(trans == 'n' .or. trans == 'N') then
    nx = n
    ny = m
  else
    nx = m
    ny = n
  end if
  !
  allocate(A(m+esa, n), Aa(m+esa, n), x(nx*incx), xa(nx*incx), ya(ny*incy), y0(ny*incy))
  !
  A(1:m, 1:n) = reshape(u(1:m*n), [m, n])
  Aa(1:m, 1:n) = reshape(ua(1:m*n), [m, n])
  call scopy(nx, u(m*n+1), 1, x, incx)
  call scopy(nx, ua(m*n+1), 1, xa, incx)
  call scopy(ny, va(1), 1, ya, incy)
  call scopy(ny, u(m*n+nx+1), 1, y0, incy)
  alpha = u(nn-1)
  beta = u(nn)
  alphaa = ua(nn-1)
  betaa = ua(nn)
  !
  call sgemv_rmds(trans, m, n, A, m+esa, x, incx, y0, incy, alphaa, betaa, ya, '11')
  call sgemv_rmd(trans, m, n, alpha, A, m+esa, x, incx, beta, incy, Aa, xa, ya, '111')
  !
  ua(1:m*n) = reshape(Aa(1:m, 1:n), [m*n])
  call scopy(nx, xa, incx, ua(m*n+1), 1)
  call scopy(ny, ya, incy, ua(m*n+nx+1), 1)
  ua(nn-1) = alphaa
  ua(nn) = betaa
  !
END SUBROUTINE F_RMD
