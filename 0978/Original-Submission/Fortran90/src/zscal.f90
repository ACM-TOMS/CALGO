subroutine ZSCAL( n, alpha, x, incx )
   use LA_CONSTANTS, only: wp, cone
!
!  Updated Level 1 BLAS
!  E. Anderson
!  August 11, 2016
!
!  .. Scalar Arguments ..
   integer :: incx, n
   complex(wp) :: alpha
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!
!  Purpose
!  =======
!
!  ZSCAL multiplies an n-element vector x by a complex scalar alpha.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vector x.
!
!  ALPHA   (input) COMPLEX
!          The scalar alpha.
!
!  X       (input/output) COMPLEX array, dimension (1+(N-1)*abs(INCX))
!          On entry, the n-element vector x.
!          On exit, the scaled vector alpha*x.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector x.
!          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!          If INCX = 0, x isn't a vector so there is no need to call
!          this subroutine.  If you call it anyway, it will scale x(1)
!          by alpha N times.
!
! =====================================================================
!
!  .. Local Scalars ..
   integer :: i, ix
!  ..
!
!  Quick return if possible
!
   if( n <= 0 .or. alpha == cone ) return
!
   if( incx == 1 ) then
!
!     Special case for incx = 1
!
      do i = 1, n
         x(i) = alpha*x(i)
      end do
   else
!
!     General increment case
!
      ix = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      do i = 1, n
         x(ix) = alpha*x(ix)
         ix = ix + incx
      end do
   end if
   return
end subroutine
