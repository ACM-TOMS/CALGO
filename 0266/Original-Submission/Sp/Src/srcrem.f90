!     REMARK ON ALGORITHM 266, COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN CACM, VOL 9, NO. 9,
!     SEP., 1966, P. 687--687.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!  UNDER NO CIRCUMSTANCES SHOULD THIS RANDOM NUMBER
!!!!  GENERATOR BE USED FOR SERIOUS COMPUTING
!!!!
!!!!  Code provided solely for archival purposes
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module alg266r
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
! Implementation of acm remark on algorithm 266. This is provided
! purely for historical reasons and should never be used
! in practice.
!
! The implementation assumes that the default integer is
! at least 29 bits (excluding the sign).
!
! The returned value is in the range (a,b).
!
! The generator used is
!   x(n+1) = 125*x(n) mod 2796203
! which has a period of ?????
!
  integer, parameter :: modval = 2796203
  real, parameter :: mrecip = 1.0/2796203

  x = mod (125 * x, modval)


! Generate a real value in the required range
  random = real(x)*mrecip*(b-a) + a

  end function random

end module
