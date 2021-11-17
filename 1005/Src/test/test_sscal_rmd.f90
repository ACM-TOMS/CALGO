PROGRAM TEST_SSCAL_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  integer :: m, rep, incx
  real :: tol
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  rep = 1
  !
  do m = 1, 5
    do incx = 1, 4
      testcase = fdata([real::], [m, incx], [character::])
      call rmd_stestrandom(F, F_rmd, m+1, m, tol, rep, testcase)
    end do
  end do
END PROGRAM TEST_SSCAL_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: x(:)
  integer :: m, incx
  real :: alpha
  !
  m = testcase % integers(1)
  incx = testcase % integers(2)
  !  
  allocate(x(m*incx))
  call scopy(m, u, 1, x, incx)
  alpha = u(m+1)
  !  
  call sscal(m, alpha, x, incx)
  !  
  call scopy(m, x, incx, v, 1)
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
  real, allocatable :: xa(:), x0(:)
  integer :: m, incx
  real :: alpha, alphaa
  !
  m = testcase % integers(1)
  incx = testcase % integers(2)
  !
  allocate(xa(m*incx), x0(m*incx))
  call scopy(m, u, 1, x0, incx)
  call scopy(m, va, 1, xa, incx)
  alpha = u(m+1)
  alphaa = ua(m+1)
  !  
  call sscal_rmds(m, x0, incx, alphaa, xa)
  call sscal_rmd(m, alpha, incx, xa)
  !
  call scopy(m, xa, incx, ua, 1)
  ua(m+1) = alphaa
  !
END SUBROUTINE F_RMD
