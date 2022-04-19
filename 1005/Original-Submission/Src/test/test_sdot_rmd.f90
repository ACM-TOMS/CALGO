PROGRAM TEST_SDOT_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  integer :: n, rep, incx, incy
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  rep = 1
  !   
  do n = 1, 5
    do incx = 1, 3, 2
    do incy = 1, 3, 2
      testcase = fdata([real::], [n, incx, incy], ['0'])
      call rmd_stestrandom(F, F_rmd, 2*n, 1, tol, rep, testcase)
    end do
    end do
  end do
  !
END PROGRAM TEST_SDOT_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, external :: sdot
  !
  real, allocatable :: x(:), y(:)
  integer :: n, incx, incy
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !
  allocate(x(n*incx), y(n*incy))
  call scopy(n, u(1:n), 1, x, incx)
  call scopy(n, u(n+1:2*n), 1, y, incy)
  !
  v(1) = sdot(n, x, incx, y, incy)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: x(:), y(:), xa(:), ya(:)
  real :: dota !, dot
  integer :: n, incx, incy
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !
  allocate(x(n*incx), y(n*incy))
  allocate(xa(n*incx), ya(n*incy))
  !
  call scopy(n, u(1:n), 1, x, incx)
  call scopy(n, u(n+1:2*n), 1, y, incy)
  call scopy(n, ua(1:n), 1, xa, incx)
  call scopy(n, ua(1:n), 1, ya, incy)
  ! dot = v(1)
  dota = va(1)
  !
  call sdot_rmd(n, x, incx, y, incy, dota, xa, ya, '11')
  !
  call scopy(n, xa, incx, ua(1:n), 1)
  call scopy(n, ya, incy, ua(n+1:n*2), 1)
  !
END SUBROUTINE F_RMD
