!     ALGORITHM 133, COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN CACM, VOL 5, NO. 11,
!     NOV., 1962, P. 553.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!  UNDER NO CIRCUMSTANCES SHOULD THIS RANDOM NUMBER
!!!!  GENERATOR BE USED FOR SERIOUS COMPUTING
!!!!
!!!!  Code provided solely for archival purposes
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module alg133
  integer :: x(2)

contains
  subroutine setseed(y)
  integer, intent(in) :: y(2)
  
! set the initial seed of the generator. This value is
! stored as y(1)*2**20 + y(2).
!
! Thus the value 28 395 423 107 is stored as
! y(1) = 27079; y(2) = 1033603
!
  x = y

  end subroutine setseed

  real function random(a, b)
  real, intent(in) :: a, b

!
! Implementation of acm algorithm 133. This is provided
! purely for historical reasons and should never be used
! in practice.
!
! The implementation assumes that the default integer is
! at least 23 bits (excluding the sign).
!
! The returned value is in the range (a,b).
!
! The generator used is
!   x(n+1) = 5*x(n) mod 2^35
! which has a period of 2^33.
!
! mscale = 2^20, topscale = 2^(35-20) = 2^15;
! mrecip = 1.0/2^35; toprecip = 1.0/2^15
!
  integer, parameter :: mult = 5, mscale = 1048576, topscale = 32768
  real, parameter :: mrecip = 2.0**(-35), toprecip = 2.0**(-15)

  x = mult * x

! Add overflow from x(2) and remove it from x(2)
  if (x(2) >= mscale) then
    x(1) = x(1) + x(2)/mscale
    x(2) = mod(x(2), mscale)
  endif

! Possibly need to tidy up x(1)
  if (x(1) >= topscale) then
    x(1) = mod(x(1), topscale)
  endif

! Generate a real value in the required range
  random = (x(1)*toprecip + x(2)*mrecip)*(b-a) + a

  end function random

end module
