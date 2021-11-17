PROGRAM TEST_SROT_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd

  integer, parameter :: nmax = 4
  integer :: n, dimu, dimv, incx, incy
  type(fdata) :: testcase
  real :: tol = rmd_stolerance, c, s
  real :: u(2*nmax + 2)
  logical :: mask(2*nmax + 2)
  !
  do n = 1,nmax
    do incx = 1,2,3
      do incy = 1,2,3
        dimu = 2*n + 2
        dimv = 2*n
        call random_number(u(1:2*n))
        call random_number(c)
        s = sqrt(1.0 - c**2)
        u(2*n+1:2*n+2) = [c, s]
        testcase = fdata([real::], [n, incx, incy, 1, 1], [character::])
        ! All selected:
        call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase)
        ! accum gives indices in u for which to check that F_rmd updates
        ! correctly
        call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase, &
          & accum = [dimu-1, dimu])
        ! Only ca, sa selected:
        mask = .false.
        mask(dimu-1:dimu) = .true.
        testcase % integers(4:5) = [0,1]
        call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase, mask = mask(1:dimu))
        ! Only xa, ya selected:
        mask = .true.
        mask(dimu-1:dimu) = .false.
        testcase % integers(4:5) = [1,0]
        call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase, mask = mask(1:dimu))
      enddo
    enddo
  enddo
END PROGRAM TEST_SROT_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  real :: c, s
  integer :: n, incx, incy
  real, allocatable :: x(:), y(:)
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  allocate(x(n*incx), y(n*incy))
  call scopy(n, u, 1, x, incx)
  call scopy(n, u(n+1), 1, y, incy)
  c = u(2*n + 1)
  s = u(2*n + 2)
  call srot(n, x, incx, y, incy, c, s)
  call scopy(n, x, incx, v, 1)
  call scopy(n, y, incy, v(n+1), 1)
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  real :: c, s, ca, sa
  integer :: n, incx, incy, selxy, selcs
  real, allocatable :: x(:), y(:), xa(:), ya(:)
  character(2) :: sel
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  selxy = testcase % integers(4)
  selcs = testcase % integers(5)
  sel = '00'
  if (selxy == 1) sel(1:1) = '1'
  if (selcs == 1) sel(2:2) = '1'
  allocate(x(n*incx), y(n*incy), xa(n*incx), ya(n*incy))
  call scopy(n, v, 1, x, incx)
  call scopy(n, v(n+1), 1, y, incy)
  call scopy(n, va, 1, xa, incx)
  call scopy(n, va(n+1), 1, ya, incy)
  c = u(2*n + 1)
  s = u(2*n + 2)
  ca = ua(2*n + 1)
  sa = ua(2*n + 2)
  call srot_rmd(n, x, incx, y, incy, c, s, xa, ya, ca, sa, sel)
  call scopy(n, xa, incx, ua, 1)
  call scopy(n, ya, incy, ua(n+1), 1)
  ua(2*n + 1) = ca
  ua(2*n + 2) = sa
END SUBROUTINE F_RMD
