PROGRAM TEST_SROTG_RMD
  use rmd_stesttools
  implicit none
  integer :: i
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd  
  type(fdata) :: testcase
  real :: tol = rmd_stolerance, u(2)
  testcase = fdata([real::], [integer::], [character::])
  !
  do i = 1,10
    call random_number(u)
    u = u*0.8 + 0.1
    call rmd_stestf(F, F_rmd, u, 2, 3, tol, testcase)
  enddo
END PROGRAM TEST_SROTG_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  real :: a, b, d, c, s
  !
  a = u(1)
  b = u(2)
  call srotg(a, b, c, s)
  d = a
  v(1) = c
  v(2) = s
  v(3) = d
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  real :: c, s, d, aa, ba, ca, sa, da
  !
  c = v(1)
  s = v(2)
  d = v(3)
  aa = ua(1)
  ba = ua(2)
  ca = va(1)
  sa = va(2)
  da = va(3)
  call srotg_rmd(c, s, d, aa, ba, ca, sa, da)
  ua(1) = aa
  ua(2) = ba
END SUBROUTINE F_RMD
