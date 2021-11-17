!     ALGORITHM 266, COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN CACM, VOL 8, NO. 10,
!     OCT., 1965, P. 605--606.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!  UNDER NO CIRCUMSTANCES SHOULD THIS RANDOM NUMBER
!!!!  GENERATOR BE USED FOR SERIOUS COMPUTING
!!!!
!!!!  Code provided solely for archival purposes
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module alg266
  integer :: x

contains
  subroutine setseed(y)
  integer, intent(in) :: y
  
! set the initial seed of the generator. 
!
  x = y

  end subroutine setseed

  real function random(a, b)
  real, intent(in) :: a, b

!
! Implementation of acm algorithm 266. This is provided
! purely for historical reasons and should never be used
! in practice.
!
! The implementation assumes that the default integer is
! at least 31 bits (excluding the sign).
!
! The returned value is in the range (a,b).
!
! The generator used is
!   x(n+1) = 3125*x(n) mod 2^26
! which has a period of 2^24.
!
  integer, parameter :: multiplier(3) = (/25,25,5/)
  integer, parameter :: modval = 67108864
  real, parameter :: mrecip = 2.0**(-26)
  integer :: i

  do i = 1, 3
    x = mod (multiplier(i) * x, modval)
  end do


! Generate a real value in the required range
  random = real(x)*mrecip*(b-a) + a

  end function random

end module

