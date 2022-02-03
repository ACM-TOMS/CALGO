
program test

use alg133

integer :: i
real :: xval

integer, parameter :: nvals = 5000
integer :: seed(2) = (/0, 32767/)


call setseed(seed)

do i = 1, nvals
  xval = random(0.0e0, 1.0e0)
end do

! Check 5000th number is as expected
if (x(1) == 23089 .and. x(2) == 808991) then
  write(6,*) "5000th number determined correctly"
else
  write(6,*) "IMPLEMENTATION ERROR"
endif

end
