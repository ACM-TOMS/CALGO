module LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  May 13, 2016
!
   interface LA_XVALS

   module procedure SXVALS
   module procedure DXVALS

   end interface

contains
   
   function SXVALS( k, xx )
   use LA_CONSTANTS32, only: wp, zero
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#endif
   real(wp) :: SXVALS
!
!  .. Scalar Arguments ..
   integer :: k
   real(wp) :: xx
!  ..
!
!  Purpose
!  =======
!
!  SXVALS returns one of the IEEE exceptional values:
!    1:  -Infinity
!    2:  +Infinity
!    3:  NaN
!  The intrinsic IEEE_VALUE is used if available, otherwise
!  the values are computed.
!
!  Arguments
!  =========
!
!  K       (input) INTEGER
!          If 1 <= K <= 3, SXVALS returns one of the IEEE exceptional
!          values as shown above; otherwise, SXVALS returns 0.
!
!  XX      (input) REAL
!          A variable of the desired KIND, used only for interface
!          matching.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: x
#ifndef USE_IEEE_INTRINSIC
   real(wp) :: y, z
#endif
!  ..
!
   x = zero
#ifdef USE_IEEE_INTRINSIC
   if( k == 1 ) then
     x = IEEE_VALUE( xx, ieee_negative_inf )
   else if( k == 2 ) then
     x = IEEE_VALUE( xx, ieee_positive_inf )
   else if( k == 3 ) then
     x = IEEE_VALUE( xx, ieee_quiet_nan )
   end if
#else
   y = HUGE(xx)
   z = y*y
   if( k == 1 ) then
     x = -z
   else if( k == 2 ) then
     x = z
   else if( k == 3 ) then
     x = z/z
   end if
#endif
   SXVALS = x
end function SXVALS

   function DXVALS( k, xx )
   use LA_CONSTANTS, only: wp, zero
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#endif
   real(wp) :: DXVALS
!
!  .. Scalar Arguments ..
   integer :: k
   real(wp) :: xx
!  ..
!
!  Purpose
!  =======
!
!  DXVALS returns one of the IEEE exceptional values:
!    1:  -Infinity
!    2:  +Infinity
!    3:  NaN
!  The intrinsic IEEE_VALUE is used if available, otherwise
!  the values are computed.
!
!  Arguments
!  =========
!
!  K       (input) INTEGER
!          If 1 <= K <= 3, DXVALS returns one of the IEEE exceptional
!          values as shown above; otherwise, DXVALS returns 0.
!
!  XX      (input) REAL
!          A variable of the desired KIND, used only for interface
!          matching.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: x
#ifndef USE_IEEE_INTRINSIC
   real(wp) :: y, z
#endif
!  ..
!
   x = zero
#ifdef USE_IEEE_INTRINSIC
   if( k == 1 ) then
     x = IEEE_VALUE( xx, ieee_negative_inf )
   else if( k == 2 ) then
     x = IEEE_VALUE( xx, ieee_positive_inf )
   else if( k == 3 ) then
     x = IEEE_VALUE( xx, ieee_quiet_nan )
   end if
#else
   y = HUGE( xx )
   z = y*y
   if( k == 1 ) then
     x = -z
   else if( k == 2 ) then
     x = z
   else if( k == 3 ) then
     x = z/z
   end if
#endif
   DXVALS = x
end function DXVALS

end module LA_XXVALS
