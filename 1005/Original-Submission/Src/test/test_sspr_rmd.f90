PROGRAM TEST_SSPR_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer :: n, dimu, dimv, incx
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  do n = 1, 5
    do incx = 1, 3, 2
      dimu  = (n*(n+1))/2 + n + 1
      dimv = (n*(n+1))/2
      testcase = fdata([real::], [n, incx, dimu, dimv], ['L'])
      call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
      testcase % characters(1) = 'U'
      call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
    end do
  end do
END PROGRAM TEST_SSPR_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: n, incx, dimu, dimv
  real :: alpha
  character uplo
  real, allocatable ::x(:), AP(:)
  !
  uplo = testcase % characters(1)
  n    = testcase % integers(1)
  incx = testcase % integers(2)
  dimu = testcase % integers(3)
  dimv = testcase % integers(4)
  !
  allocate(x(n*incx), AP(dimv))
  !
  AP = u(1:dimv)
  call scopy(n, u(dimv+1), 1, x, incx)
  alpha = u(dimu)
  !
  call sspr(uplo, n, alpha, x, incx, AP)
  v(1:dimv) = AP
  !  
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)     :: u(*), va(*), v(*)
  real, intent(inout)  :: ua(*)
  !
  ! Local variables
  real, allocatable :: x(:), xa(:), APa(:)
  integer :: n, incx, dimu, dimv
  real :: alpha, alphaa
  character uplo
  !
  uplo = testcase % characters(1)
  n    = testcase % integers(1)
  incx = testcase % integers(2)
  dimu = testcase % integers(3)
  dimv = testcase % integers(4)
  !
  allocate(x(n*incx), xa(n*incx), APa(dimv))
  !
  call scopy(n, u(dimv+1), 1, x, incx)
  call scopy(n, ua(dimv+1), 1, xa, incx)
  APa = va(1:dimv)
  alpha = u(dimu)
  alphaa = ua(dimu)
  !
  call sspr_rmds(uplo, n, x, incx, alphaa, APa) 
  call sspr_rmd(uplo, n, alpha, x, incx, xa, APa)
  !
  ua(1:dimv) = APa
  call scopy(n, xa, incx, ua(dimv+1), 1)
  ua(dimu) = alphaa
  !
END SUBROUTINE F_RMD
