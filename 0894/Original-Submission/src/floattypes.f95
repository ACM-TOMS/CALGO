! This module defines the suffices of the intrinsic floating-point data types.
!
! The suffix ``sp'' is the single    precision, (unit round-off = 2^{ -24}),
!            ``wp'' is the double    precision, (unit round-off = 2^{ -53}),
!            ``tp'' is the triple    precision, (unit round-off = 2^{ -64}), and
!            ``qp'' is the quadruple precision. (unit round-off = 2^{-113}).
!
! The program works correctly even if the actual unit round-off is smaller,
! though we cannot fully exploit the corresponding higher precision.
!
! The two extended precisions, ``tp'' and ``qp'', are disabled in some of
! the makefiles provided herewith, because they are not always supported.
! When one of them is available, please add either ``-D__USE_TPREC'' or
! ``-D__USE_QPREC'' to the command line arguments of your FORTRAN compiler.
!
! This module is intended for internal use only.

module floattypes

  implicit none

  integer, parameter :: sp = selected_real_kind( 6,  37)     ! u := 2^{ -24}
  integer, parameter :: wp = selected_real_kind(15, 307)     ! u := 2^{ -53}
#ifdef __USE_TPREC
  integer, parameter :: tp = selected_real_kind(18,4900)     ! u := 2^{ -64}
#endif
#ifdef __USE_QPREC
  integer, parameter :: qp = selected_real_kind(33,4900)     ! u := 2^{-113}
#endif

  public

end module floattypes

