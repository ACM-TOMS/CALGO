PROGRAM TEST_STBMV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  integer :: n, k, ul, tr, di, rep, esa, incx, dimu
  type(fdata) :: testcase
  character :: uplo(2), trans(2), diag(2)
  !
  tol = rmd_stolerance
  rep = 1
  uplo = ['u', 'l']
  trans = ['n', 't']
  diag = ['n', 'u']
  !
  do n = 1, 5
    do k = 0, n-1
      do ul = 1, 2
        do tr = 1, 2
          do di = 1, 2
            do esa = 0, 3, 3
              do incx = 1, 3, 2
                dimu = rmd_dimband(n, n, k, 0) + n
                testcase = fdata([real::], [n, k, esa, incx], [uplo(ul), trans(tr), diag(di)])
                call rmd_stestrandom(F, F_rmd, dimu, n, tol, rep, testcase)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
END PROGRAM

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), x(:)
  integer :: n, k, esa, incx, dimu
  character uplo, trans, diag
  !
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  trans = testcase % characters(2)
  diag = testcase % characters(3)
  esa = testcase % integers(3)
  incx = testcase % integers(4)
  dimu = rmd_dimband(n, n, k, 0) + n
  !
  allocate(A(k+1+esa, n), x(n*incx))
  A = 777
  !
  if(uplo == 'u' .or. uplo == 'U') then
    call rmd_sbandshape_v2m(n, n, 0, k, u, A, k+1+esa)
  else
    call rmd_sbandshape_v2m(n, n, k, 0, u, A, k+1+esa)
  end if
  call scopy(n, u(dimu - n + 1), 1, x, incx)
  !
  call stbmv(uplo, trans, diag, n, k, A, k+1+esa, x, incx)
  !
  call scopy(n, x, incx, v, 1)
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
  real, allocatable :: A(:, :), x(:), Aa(:, :), xa(:)
  integer :: n, k, esa, incx, dimu
  character uplo, trans, diag
  !
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  trans = testcase % characters(2)
  diag = testcase % characters(3)
  esa = testcase % integers(3)
  incx = testcase % integers(4)
  dimu = rmd_dimband(n, n, k, 0) + n
  !
  allocate(A(k+1+esa, n), x(n*incx))
  allocate(Aa(k+1+esa, n), xa(n*incx))
  !
  if(uplo == 'u' .or. uplo == 'U') then
    call rmd_sbandshape_v2m(n, n, 0, k, u, A, k+1+esa)
    call rmd_sbandshape_v2m(n, n, 0, k, ua, Aa, k+1+esa)
  else
    call rmd_sbandshape_v2m(n, n, k, 0, u, A, k+1+esa)
    call rmd_sbandshape_v2m(n, n, k, 0, ua, Aa, k+1+esa)
  end if
  call scopy(n, u(dimu - n + 1), 1, x, incx) ! x0
  call scopy(n, va, 1, xa, incx)
  !
  call stbmv_rmd(uplo, trans, diag, n, k, A, k+1+esa, x, incx, Aa, xa, '11')
  !
  if(uplo == 'u' .or. uplo == 'U') then
    call rmd_sbandshape_m2v(n, n, 0, k, ua, Aa, k+1+esa)
  else
    call rmd_sbandshape_m2v(n, n, k, 0, ua, Aa, k+1+esa)
  end if
  call scopy(n, xa, incx, ua(dimu - n + 1), 1)
  !
END SUBROUTINE F_RMD
