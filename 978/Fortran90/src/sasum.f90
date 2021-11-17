function SASUM( n, x, incx )
   use LA_CONSTANTS32, only: wp, zero
   real(wp) :: SASUM
!
!  Updated Level 1 BLAS
!  E. Anderson
!  March 5, 2015
!
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!
!  Purpose
!  =======
!
!  SASUM computes the 1-norm of an n-element real vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vector x.
!
!  X       (input) REAL array, dimension (1+(N-1)*abs(INCX))
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
!  ..
!
!  Quick return if possible
!
   SASUM = zero
   if( n <= 0 ) return
!
   stmp = zero
   if( incx == 1 ) then
!
!     Special case for incx = 1
!
      do i = 1, n
         stmp = stmp + abs(x(i))
      end do
   else
!
!     General increment case
!
      ix = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      do i = 1, n
         stmp = stmp + abs(x(ix))
         ix = ix + incx
      end do
   end if
   SASUM = stmp
   return
end function
