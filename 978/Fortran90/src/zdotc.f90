function ZDOTC( n, x, incx, y, incy )
   use LA_CONSTANTS, only: wp, czero
   complex(wp) :: ZDOTC
!
!  Updated Level 1 BLAS
!  E. Anderson
!  March 17, 2015
!
!  .. Scalar Arguments ..
   integer :: incx, incy, n
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*), y(*)
!  ..
!
!  Purpose
!  =======
!
!  ZDOTC computes the dot product of two vectors x and y,
!  conjugating the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vectors x and y.
!
!  X       (input) COMPLEX array, dimension (1+(N-1)*abs(INCX))
!          The n-element vector x.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector x.
!          If INCX >= 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!
!  Y       (input) COMPLEX array, dimension (1+(N-1)*abs(INCY))
!          The n-element vector y.
!
!  INCY    (input) INTEGER
!          The increment between successive values of the vector y.
!          If INCY >= 0, Y(1+(i-1)*INCY) = y(i) for 1 <= i <= n
!          If INCY < 0, Y(1-(n-i)*INCY) = y(i) for 1 <= i <= n
!
! =====================================================================
!
!  .. Local Scalars ..
   integer :: i, ix, iy
   complex(wp) :: stmp
!  ..
!
!  Quick return if possible
!
   ZDOTC = czero
   if( n <= 0 ) return
!
   stmp = czero
   if( incx == 1 .and. incy == 1 ) then
!
!     Special case for incx = incy = 1
!
      do i = 1, n
         stmp = stmp + conjg(x(i))*y(i)
      end do
   else
!
!     General increment case
!
      ix = 1
      iy = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      if( incy < 0 ) iy = 1 - (n-1)*incy
      do i = 1, n
         stmp = stmp + conjg(x(ix))*y(iy)
         ix = ix + incx
         iy = iy + incy
      end do
   end if
   ZDOTC = stmp
   return
end function
