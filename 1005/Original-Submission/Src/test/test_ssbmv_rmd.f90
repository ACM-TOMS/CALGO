PROGRAM TEST_SSBMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  integer :: n, k, ul, rep, esa, incx, incy, dimu
  real :: tol
  type(fdata) :: testcase
  character :: uplo(2)
  !
  tol = rmd_stolerance
  rep = 1
  uplo = ['u', 'l']
  !
  do n = 1, 5
    do k = 0, n-1
      do ul = 1, 2
        do esa = 0, 3, 3
          do incx = 1, 3, 2
            do incy = 1, 3, 2
              dimu = rmd_dimband(n, n, k, 0) + n + n + 2
              !dimu = (k+1)*n + 2*n + 2
              testcase = fdata([real::], [n, k, esa, incx, incy, dimu], [uplo(ul)])
              call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
            end do
          end do
        end do
      end do
    end do
  end do
  !  
END PROGRAM TEST_SSBMV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  !  
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), x(:), y(:)
  integer :: n, k, esa, incx, incy, dimu, lda, sizeA
  real :: alpha, beta
  character uplo
  !
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  esa = testcase % integers(3)
  incx = testcase % integers(4)
  incy = testcase % integers(5)
  dimu = testcase % integers(6)
  alpha = u(dimu - 1)
  beta = u(dimu)
  lda = k + 1 + esa
  !
  allocate(A(lda, n), x(n*incx), y(n*incy))
  !
  sizeA = dimu - 2 - 2*n
  if(uplo == 'u' .or. uplo == 'U') then
    call rmd_sbandshape_v2m(n, n, 0, k, u, A, lda)
  else
    call rmd_sbandshape_v2m(n, n, k, 0, u, A, lda)
  end if
  call scopy(n, u(sizeA + 1), 1, x, incx)
  call scopy(n, u(sizeA + n + 1), 1, y, incy)
  !
  call ssbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
  !
  call scopy(n, y, incy, v, 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  !  
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), x(:), y(:), Aa(:, :), xa(:), ya(:), y0(:)
  integer :: n, k, esa, incx, incy, dimu, lda, sizeA
  real :: alpha, beta, alphaa, betaa
  character uplo
  !
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  esa = testcase % integers(3)
  incx = testcase % integers(4)
  incy = testcase % integers(5)
  dimu = testcase % integers(6)
  alpha = u(dimu - 1)
  beta = u(dimu)
  alphaa = ua(dimu - 1)
  betaa = ua(dimu)
  lda = k + 1 + esa
  !
  allocate(A(lda, n), x(n*incx), y(n*incy), y0(n*incy))
  allocate(Aa(lda, n), xa(n*incx), ya(n*incy))
  !
  if(uplo == 'u' .or. uplo == 'U') then
    call rmd_sbandshape_v2m(n, n, 0, k, u, A, lda)
    call rmd_sbandshape_v2m(n, n, 0, k, ua, Aa, lda)
  else
    call rmd_sbandshape_v2m(n, n, k, 0, u, A, lda)
    call rmd_sbandshape_v2m(n, n, k, 0, ua, Aa, lda)
  end if
  sizeA = dimu - 2 - 2*n
  call scopy(n, u(sizeA + 1), 1, x, incx)
  call scopy(n, ua(sizeA + 1), 1, xa, incx)
  call scopy(n, v, 1, y, incy)
  call scopy(n, va, 1, ya, incy)
  call scopy(n, u(sizeA + n + 1), 1, y0, incy)
  !
  call ssbmv_rmds(uplo, n, k, A, lda, x, incx, y0, incy, alphaa, betaa, ya, '11')
  call ssbmv_rmd(uplo, n, k, alpha, A, lda, x, incx, beta, incy, Aa, xa, ya, '111')
  !
  if(uplo == 'u' .or. uplo == 'U') then
    call rmd_sbandshape_m2v(n, n, 0, k, ua, Aa, lda)
  else
    call rmd_sbandshape_m2v(n, n, k, 0, ua, Aa, lda)
  end if
  call scopy(n, xa, incx, ua(sizeA + 1), 1)
  call scopy(n, ya, incy, ua(sizeA + n + 1), 1)
  ua(dimu - 1) = alphaa
  ua(dimu) = betaa
  !
END SUBROUTINE F_RMD
