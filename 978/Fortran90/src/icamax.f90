integer function ICAMAX( n, x, incx )
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
   complex(wp) :: x(*)
!  ..
!
!  Purpose
!  =======
!
!  ICAMAX finds the index of the element of maximum absolute value
!  in the n-element complex vector x. If the vector contains NaNs,
!  ICAMAX returns the index of the largest element excluding the
!  NaN entries. If all entries are NaN, ICAMAX returns 1.
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
!          If INCX = 0, all elements of x are the same so ICAMAX will
!          return 1.
!
! =====================================================================
!
!  .. Local Scalars ..
   integer :: i, imax, ix
   real(wp) :: smax
   complex(wp) :: zdum
!  ..
!  .. Statement Functions ..
   real(wp) :: CABS1
   CABS1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
!  ..
!
!  Quick return if possible
!
   ICAMAX = 0
   if( n <= 0 ) return
   ICAMAX = 1
   if( n == 1 .or. incx == 0 ) return
!
   imax = 1
   smax = -one
   if( incx == 1 ) then
!
!     Special case for incx = 1
!
      do i = 1, n
         if( CABS1(x(i)) > smax ) then
            imax = i
            smax = CABS1(x(i))
         end if
      end do
   else
!
!     General increment case
!
      ix = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      do i = 1, n
         if( CABS1(x(ix)) > smax ) then
            imax = i
            smax = CABS1(x(ix))
         end if
         ix = ix + incx
      end do
   end if
   ICAMAX = imax
   return
end function
