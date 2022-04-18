PROGRAM TEST_SNRM2_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  integer :: n, rep, incx
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  rep = 1
  !
  do n = 3, 3
    do incx = 1, 3, 2
      testcase = fdata([real::], [n, incx], ['0'])
      call rmd_stestrandom(F, F_rmd, n, 1, tol, rep, testcase)
    end do
  end do
END PROGRAM TEST_SNRM2_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, external :: snrm2
  !
  real, allocatable :: x(:)
  integer :: n, incx
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  !
  allocate(x(n*incx))
  call scopy(n, u(1:n), 1, x, incx)
  !
  v(1) = snrm2(n, x, incx)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  !
  real, allocatable :: x(:), xa(:)
  real :: nrm2, nrm2a
  integer :: n, incx
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  !
  allocate(x(n*incx), xa(n*incx))
  call scopy(n, u(1:n), 1, x, incx)
  call scopy(n, ua(1:n), 1, xa, incx)
  nrm2 = v(1)
  nrm2a = va(1)
  !
  call snrm2_rmd(n, x, incx, nrm2, xa, nrm2a)
  !
  call scopy(n, xa, incx, ua(1:n), 1)
  !
END SUBROUTINE F_RMD
