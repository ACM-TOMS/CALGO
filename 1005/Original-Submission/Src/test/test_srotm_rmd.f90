PROGRAM TEST_SROTM_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd

  integer, parameter :: nmax = 4
  integer :: i, n, dimu, dimv, incx, incy
  type(fdata) :: testcase
  real :: tol = rmd_stolerance
  real :: u(2*nmax + 4)
  logical :: mask(2*nmax + 4)
  integer :: flag
  !
  do n = 1,nmax
    do incx = 1,2,3
      do incy = 1,2,3
        do flag = 0,0
          select case(flag)
          case(-2)
            dimu = 2*n
          case(-1)
            dimu = 2*n + 4
          case(0:1)
            dimu = 2*n + 2
          endselect
          dimv = 2*n
          call random_number(u(1:dimu))
          testcase = fdata([real::], [n, incx, incy, flag, dimu, 1, 1], [character::])
          ! All selected:
          call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase)
          ! accum gives indices in u for which to check that F_rmd updates
          ! correctly
          call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase, &
            &             accum = [(i, i = 2*n+1, dimu)])
          ! Only param selected:
          mask = .false.
          mask(2*n+1:dimu) = .true.
          testcase % integers(6:7) = [0,1]
          call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase, mask = mask(1:dimu))
          ! Only xa, ya selected:
          mask = .true.
          mask(2*n+1:dimu) = .false.
          testcase % integers(6:7) = [1,0]
          call rmd_stestf(F, F_rmd, u(1:dimu), dimu, dimv, tol, testcase, mask = mask(1:dimu))
        enddo
      enddo
    enddo
  enddo
END PROGRAM TEST_SROTM_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  real :: param(5)
  integer :: n, incx, incy, flag, dimu
  real, allocatable :: x(:), y(:)
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  flag = testcase % integers(4)
  dimu = testcase % integers(5)
  allocate(x(n*incx), y(n*incy))
  call scopy(n, u, 1, x, incx)
  call scopy(n, u(n+1), 1, y, incy)
  param(1) = flag
  select case(flag)
  case(-1)
    param(2:5) = u(2*n+1:dimu)
  case(0)
    param(3:4) = u(2*n+1:dimu)
  case(1)
    param(2) = u(2*n+1)
    param(5) = u(dimu)
  end select
  call srotm(n, x, incx, y, incy, param)
  call scopy(n, x, incx, v, 1)
  call scopy(n, y, incy, v(n+1), 1)
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  real :: param(5), parama(5)
  integer :: n, incx, incy, selxy, selparam, flag, dimu
  real, allocatable :: x(:), y(:), xa(:), ya(:)
  character(2) :: sel
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  flag = testcase % integers(4)
  dimu = testcase % integers(5)
  selxy = testcase % integers(6)
  selparam = testcase % integers(7)
  sel = '00'
  if (selxy == 1) sel(1:1) = '1'
  if (selparam == 1) sel(2:2) = '1'
  allocate(x(n*incx), y(n*incy), xa(n*incx), ya(n*incy))
  call scopy(n, v, 1, x, incx)
  call scopy(n, v(n+1), 1, y, incy)
  call scopy(n, va, 1, xa, incx)
  call scopy(n, va(n+1), 1, ya, incy)
  param(1) = flag
  select case(flag)
  case(-1)
    param(2:5) = u(2*n+1:dimu)
    parama(2:5) = ua(2*n+1:dimu)
  case(0)
    param(3:4) = u(2*n+1:dimu)
    parama(3:4) = ua(2*n+1:dimu)
  case(1)
    param([2,5]) = u(2*n+1:dimu)
    parama([2,5]) = ua(2*n+1:dimu)
  end select
  call srotm_rmd(n, x, incx, y, incy, param, xa, ya, parama, sel)
  call scopy(n, xa, incx, ua, 1)
  call scopy(n, ya, incy, ua(n+1), 1)
  select case(flag)
  case(-1)
    ua(2*n+1:dimu) = parama(2:5)
  case(0)
    ua(2*n+1:dimu) = parama(3:4) 
  case(1)
    ua(2*n+1:dimu) = parama([2,5])
  end select
END SUBROUTINE F_RMD
