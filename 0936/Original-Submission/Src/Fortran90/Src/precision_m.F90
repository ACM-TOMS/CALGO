! Module define_precision provides the kind value needed
! to define the precision of a complete package along
! with definitions of all commonly used precisions
module precision_m
! Some Fortran compilers don't yet support int32, and int64.
!  use, intrinsic :: iso_fortran_env, only: int32, int64
  intrinsic kind, selected_real_kind
  ! .. Parameters .. to define the standard precisions
  integer, parameter :: kinds = kind(0.0e0)
  integer, parameter :: kindd = kind(0.0d0)
  integer, parameter :: kindq = selected_real_kind(30) ! Quad precision
  integer, parameter :: kind80 = selected_real_kind(18) ! 80 bit reals
! The last 2 above may return a kind of -1 if not supported.
!  integer, parameter :: ikind = int32 ! Integer kind
!  integer, parameter :: ikind64 = int64 ! To use 64 bit integers
end module precision_m
