integer function ISAMAX( n, x, incx )
   use LA_CONSTANTS32, only: wp, one
!
!  Updated Level 1 BLAS
!  E. Anderson
!  February 9, 2017
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
!  ISAMAX finds the index of the element of maximum absolute value
!  in the n-element real vector x. If the vector contains NaNs,
!  ISAMAX returns the index of the largest element excluding the
!  NaN entries. If all entries are NaN, ISAMAX returns 1.
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
!          If INCX = 0, all elements of x are the same so ISAMAX will
!          return 1.
!
! =====================================================================
!
!  .. Local Scalars ..
   integer :: i, imax, ix
   real(wp) :: smax
!  ..
!
!  Quick return if possible
!
   ISAMAX = 0
   if( n <= 0 ) return
   ISAMAX = 1
   if( n == 1 .or. incx == 0 ) return
!
   imax = 1
   smax = -one
   if( incx == 1 ) then
!
!     Special case for incx = 1
!
      do i = 1, n
         if( abs(x(i)) > smax ) then
            imax = i
            smax = abs(x(i))
         end if
      end do
   else
!
!     General increment case
!
      ix = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      do i = 1, n
         if( abs(x(ix)) > smax ) then
            imax = i
            smax = abs(x(ix))
         end if
         ix = ix + incx
      end do
   end if
   ISAMAX = imax
   return
end function
