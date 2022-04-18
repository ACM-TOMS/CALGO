PROGRAM TEST_SDSDOT_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol, u(3)
  integer :: n, rep, incx, incy
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  rep = 1
  !
  do n = 1, 4
    do incx = 1, 3, 2
      do incy = 1, 3, 2
        testcase = fdata([real::], [n, incx, incy], ['1', '1', '1'])
        call rmd_stestrandom(F, F_rmd, 2*n + 1, 1, tol, rep, testcase)
      end do
    end do
  end do
  !
  u = [0.333, 0.444, 0.555]
  !
  ! Check that sdsdot accumulates adjoints correctly
  testcase = fdata([real::], [1, 1, 1], ['1', '1', '1'])
  call rmd_stestf(F, F_rmd, u, 3, 1, tol, testcase,  accum = [1, 2, 3])
  !
  ! Check that sel selects
  testcase % characters = ['1', '0', '0']
  !call rmd_stestf(F, F_rmd, u, 3, 1, tol, testcase,  mask = [.true., .false., .false.])
  testcase % characters = ['0', '1', '0']
  !call rmd_stestf(F, F_rmd, u, 3, 1, tol, testcase,  mask = [.false., .true., .false.])
  testcase % characters = ['0', '0', '1']
  !call rmd_stestf(F, F_rmd, u, 3, 1, tol, testcase,  &
  !  &             mask = [.false., .false., .true.])
END PROGRAM TEST_SDSDOT_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, external :: sdsdot
  !
  real, allocatable :: x(:), y(:)
  integer :: n, incx, incy
  real :: b
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !
  allocate(x(n*incx), y(n*incy))
  b = u(1)
  call scopy(n, u(2 : n+1), 1, x, incx)
  call scopy(n, u(n+2 : 2*n+1), 1, y, incy)
  !
  v(1) = sdsdot(n, b, x, incx, y, incy)
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
  real :: dota, ba !, dot
  integer :: n, incx, incy
  character(3) :: sel
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  sel = testcase % characters(1) // testcase % characters(2) &
    &                            // testcase % characters(3)
  !
  allocate(x(n*incx), y(n*incy))
  allocate(xa(n*incx), ya(n*incy))
  !
  call scopy(n, u(2:n+1), 1, x, incx)
  call scopy(n, u(n+2:2*n+1), 1, y, incy)
  ba = ua(1)
  call scopy(n, ua(2:n+1), 1, xa, incx)
  call scopy(n, ua(n+2:2*n+1), 1, ya, incy)
  ! dot = v(1)
  dota = va(1)
  !
  call sdsdot_rmd(n, x, incx, y, incy, dota, ba, xa, ya, sel)
  !
  ua(1) = ba
  call scopy(n, xa, incx, ua(2 : n+1), 1)
  call scopy(n, ya, incy, ua(n+2 : 2*n+1), 1)
  !
END SUBROUTINE F_RMD
