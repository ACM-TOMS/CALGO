subroutine ZDSCAL( n, alpha, x, incx )
   use LA_CONSTANTS, only: wp, one
!
!  Updated Level 1 BLAS
!  E. Anderson
!  August 11, 2016
!
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: alpha
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!
!  Purpose
!  =======
!
!  ZDSCAL multiplies an n-element vector x by a real scalar alpha.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vector x.
!
!  ALPHA   (input) REAL
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
   if( n <= 0 .or. alpha == one ) return
!
   if( incx == 1 ) then
!
!     Special case for incx = 1
!
      do i = 1, n
         x(i) = cmplx(alpha*real(x(i)),alpha*aimag(x(i)),wp)
      end do
   else
!
!     General increment case
!
      ix = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      do i = 1, n
         x(ix) = cmplx(alpha*real(x(ix)),alpha*aimag(x(ix)),wp)
         ix = ix + incx
      end do
   end if
   return
end subroutine
