function SCASUM( n, x, incx )
   use LA_CONSTANTS32, only: wp, zero
   real(wp) :: SCASUM
!
!  Updated Level 1 BLAS
!  E. Anderson
!  March 17, 2015
!
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!
!  Purpose
!  =======
!
!  SCASUM computes the 1-norm of an n-element complex vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vector x.
!
!  X       (input) COMPLEX array, dimension (1+(N-1)*abs(INCX))
!          The n-element vector x.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector x.
!          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!          If INCX = 0, x isn't a vector so there is no need to call
!          this subroutine.  If you call it anyway, it will count x(1)
!          in the vector norm N times.
!
! =====================================================================
!
!  .. Local Scalars ..
   integer :: i, ix
   real(wp) :: stmp
   complex(wp) :: zdum
!  ..
!  .. Statement Functions ..
   real(wp) :: CABS1
   CABS1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
!  ..
!
!  Quick return if possible
!
   SCASUM = zero
   if( n <= 0 ) return
!
   stmp = zero
   if( incx == 1 ) then
!
!     Special case for incx = 1
!
      do i = 1, n
         stmp = stmp + CABS1(x(i))
      end do
   else
!
!     General increment case
!
      ix = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      do i = 1, n
         stmp = stmp + CABS1(x(ix))
         ix = ix + incx
      end do
   end if
   SCASUM = stmp
   return
end function
