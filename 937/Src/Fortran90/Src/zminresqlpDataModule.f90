!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zminresqlpDataModule.f90
!
! Defines precision complex(kind=dp) and a few complex constants for use
! in other modules.
!
! 20 Aug 2012: Created from minresqlpDataModule.f90 for complex constants
!              zzero and zone.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zminresqlpDataModule

  use   minresqlpDataModule, only  : dp, sp, ip, prcsn, zero, realmin, eps, one, debug, realmin, realmax, testMtx, testSymortho

  implicit none

  intrinsic :: cmplx

  complex(dp), parameter, public :: zzero = cmplx(0.0_dp,0.0_dp,dp), zone = cmplx(1.0_dp,0.0_dp,dp)

end module zminresqlpDataModule
