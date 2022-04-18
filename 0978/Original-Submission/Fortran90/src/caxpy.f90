subroutine CAXPY( n, alpha, x, incx, y, incy )
   use LA_CONSTANTS32, only: wp, czero
!
!  Updated Level 1 BLAS
!  E. Anderson
!  March 17, 2015
!
!  .. Scalar Arguments ..
   integer :: incx, incy, n
   complex(wp) :: alpha
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*), y(*)
!  ..
!
!  Purpose
!  =======
!
!  CAXPY computes the vector sum y <- alpha*x + y
!  (constant times a vector plus a vector).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vectors x and y.
!
!  ALPHA   (input) COMPLEX
!          The scalar alpha.
!
!  X       (input) COMPLEX array, dimension (1+(N-1)*abs(INCX))
!          The n-element vector x.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector x.
!          If INCX >= 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!
!  Y       (input/output) COMPLEX array, dimension (1+(N-1)*abs(INCY))
!          On entry, the n-element vector y.
!          On exit, the vector sum alpha*x + y.
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
!  ..
!
!  Quick return if possible
!
   if( n <= 0 .or. alpha == czero ) return
!
   if( incx == 1 .and. incy == 1 ) then
!
!     Special case for incx = incy = 1
!
      do i = 1, n
         y(i) = y(i) + alpha*x(i)
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
         y(iy) = y(iy) + alpha*x(ix)
         ix = ix + incx
         iy = iy + incy
      end do
   end if
   return
end subroutine
