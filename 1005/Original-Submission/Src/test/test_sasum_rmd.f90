PROGRAM TEST_SASUM_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol, u(3)
  integer :: n, incx
  type(fdata) :: testcase
  tol = rmd_stolerance
  u = 0
  do n = 1, 3
    do incx = 1, 3, 2
      testcase = fdata([real::], [n, incx], ['0'])
      call random_number(u(1:n))
      where(abs(u) < 0.1) u = sign(0.1, u)
      call rmd_stestf(F, F_rmd, u, n, 1, tol, testcase)
    end do
  end do
END PROGRAM TEST_SASUM_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  real, external :: sasum
  real, allocatable :: x(:)
  integer :: n, incx
  n = testcase % integers(1)
  incx = testcase % integers(2)
  allocate(x(n*incx))
  call scopy(n, u, 1, x, incx)
  v(1) = sasum(n, x, incx)
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  real, allocatable :: x(:), xa(:)
  real :: asuma
  integer :: n, incx
  n = testcase % integers(1)
  incx = testcase % integers(2)
  allocate(x(n*incx), xa(n*incx))
  call scopy(n, u, 1, x, incx)
  call scopy(n, ua, 1, xa, incx)
  asuma = va(1)
  call sasum_rmd(n, x, incx, xa, asuma)
  call scopy(n, xa, incx, ua, 1)
END SUBROUTINE F_RMD
