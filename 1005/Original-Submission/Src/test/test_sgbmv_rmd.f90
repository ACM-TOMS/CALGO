PROGRAM TEST_SGBMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  integer :: m, n, kl, ku, rep, dimu, esa, incx, incy
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  m = 1
  n = 2
  kl = 0
  ku = 1
  
  dimu = rmd_dimband(m, n, kl, ku) + m + n + 2
  testcase = fdata([real::], [m, n, kl, ku, 0, 1, 1, dimu], ['n']);
  call rmd_stest123(F, F_rmd, dimu, m, tol, testcase)
  !
  rep = 1
  ! Note: Must have 
  do m = 1, 5
    do n = 1, 5
      do kl = 0, m-1
        do ku = 0, n-1
          do esa = 0, 3, 3
            do incx = 1, 3, 2
              do incy = 1, 3, 2
                dimu = rmd_dimband(m, n, kl, ku) + m + n + 2
                testcase = fdata([real::], [m, n, kl, ku, esa, incx, incy, dimu], ['n'])
                call rmd_stestrandom(F, F_rmd, dimu, m, tol, rep, testcase)
                testcase % characters(1) = 't'
                call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  !
END PROGRAM TEST_SGBMV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), x(:), y(:)
  integer :: m, n, nx, ny, kl, ku, heightA, sizeA, esa, incx, incy, lda, dimu
  real :: alpha, beta
  character trans
  !  
  m = testcase % integers(1)
  n = testcase % integers(2)
  kl = testcase % integers(3)
  ku = testcase % integers(4)
  trans = testcase % characters(1)
  esa = testcase % integers(5)
  incx = testcase % integers(6)
  incy = testcase % integers(7)
  dimu = testcase % integers(8)
  !
  if(trans == 'n' .or. trans == 'N') then
    nx = n
    ny = m
  else
    nx = m
    ny = n
  end if
  heightA = kl + ku + 1
  sizeA = dimu - n - m - 2
  lda = heightA + esa
  !
  allocate(A(lda, n), x(nx*incx), y(ny*incy))
  call rmd_sbandshape_v2m(m, n, kl, ku, u, A, lda)
  call scopy(nx, u(sizeA+1), 1, x, incx)
  call scopy(ny, u(sizeA+nx+1), 1, y, incy)
  alpha = u(dimu-1)
  beta = u(dimu)
  !
  call sgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
  !
  call scopy(ny, y, incy, v, 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  use printutil
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), x(:), y(:), Aa(:, :), xa(:), ya(:), y0(:)
  integer :: m, n, kl, ku, nx, ny, heightA, sizeA, esa, incx, incy, lda, dimu
  real :: alpha, beta, alphaa, betaa
  character trans
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  kl = testcase % integers(3)
  ku = testcase % integers(4)
  trans = testcase % characters(1)
  esa = testcase % integers(5)
  incx = testcase % integers(6)
  incy = testcase % integers(7)
  dimu = testcase % integers(8)
  !
  if(trans == 'n' .or. trans == 'N') then
    nx = n
    ny = m
  else
    nx = m
    ny = n
  end if
  heightA = kl + ku + 1
  sizeA = dimu - n - m - 2
  lda = heightA + esa
  !
  allocate(A(heightA+esa, n), x(nx*incx), y(ny*incy))
  allocate(Aa(heightA+esa, n), xa(nx*incx), ya(ny*incy), y0(ny*incy))
  !
  call rmd_sbandshape_v2m(m, n, kl, ku, u, A, lda)
  call rmd_sbandshape_v2m(m, n, kl, ku, ua, Aa, lda)
  call scopy(nx, u(sizeA + 1), 1, x, incx)
  call scopy(ny, u(sizeA + nx + 1), 1, y0, incy)
  call scopy(nx, ua(sizeA + 1), 1, xa, incx)
  call scopy(ny, v, 1, y, incy)
  call scopy(ny, va, 1, ya, incy)
  alpha = u(dimu-1)
  beta  = u(dimu)
  alphaa = ua(dimu-1)
  betaa  = ua(dimu)
  !
  call sgbmv_rmds(trans, m, n, kl, ku, A, lda, x, incx, y0, incy, alphaa, betaa, ya, '11')
  call sgbmv_rmd(trans, m, n, kl, ku, alpha, A, lda, x, incx, &
    beta, incy, Aa, xa, ya, '111')
  !
  call rmd_sbandshape_m2v(m, n, kl, ku, ua, Aa, lda)
  call scopy(nx, xa, incx, ua(sizeA+1), 1)
  call scopy(ny, ya, incy, ua(sizeA+nx+1), 1)
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !  
END SUBROUTINE F_RMD
