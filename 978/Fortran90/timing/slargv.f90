subroutine SLARGV( n, x, incx, y, incy, c, incc )
   use LA_CONSTANTS32, only: wp, zero, one, safmin, safmax, rtmin, rtmax
!
!  LAPACK auxiliary routine
!  E. Anderson
!  July 30, 2016
!
!  .. Scalar Arguments ..
   integer            incc, incx, incy, n
!     ..
!  .. Array Arguments ..
   real(wp)           c( * ), x( * ), y( * )
!  ..
!
!  Purpose
!  =======
!
!  SLARGV generates a vector of real plane rotations, determined by
!  elements of the real vectors x and y. For i = 1,2,...,n
!
!     (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
!     ( -s(i)  c(i) ) ( y(i) ) = (   0  )
!
!  where c**2 + s**2 = 1.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of plane rotations to be generated.
!
!  X       (input/output) REAL array,
!                         dimension (1+(N-1)*INCX)
!          On entry, the vector x.
!          On exit, x(i) is overwritten by a(i), for i = 1,...,n.
!
!  INCX    (input) INTEGER
!          The increment between elements of X.
!
!  Y       (input/output) REAL array,
!                         dimension (1+(N-1)*INCY)
!          On entry, the vector y.
!          On exit, the sines of the plane rotations.
!
!  INCY    (input) INTEGER
!          The increment between elements of Y.
!
!  C       (output) REAL array, dimension (1+(N-1)*INCC)
!          The cosines of the plane rotations.
!
!  INCC    (input) INTEGER
!          The increment between elements of C.
!
!  =====================================================================
!
!  .. Local Scalars ..
   integer            i, ic, ix, iy
   real(wp)           d, f, f1, fs, g, g1, gs, p, u, uu
!  ..
!  .. Intrinsic Functions ..
   intrinsic          abs, sign, sqrt
!  ..
!  .. Executable Statements ..
!
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   iy = 1
   if( incy < 0 ) iy = 1 - (n-1)*incy
   ic = 1
   if( incc < 0 ) ic = 1 - (n-1)*incc
   do i = 1, n
!
!     Use identical algorithm to SLARTG
!
      f = x( ix )
      g = y( iy )
      f1 = abs( f )
      g1 = abs( g )
      if( g == zero ) then
         c( ic ) = one
         x( ix ) = f
      else if( f == zero ) then
         c( ic ) = zero
         y( iy ) = sign( one, g )
         x( ix ) = g1
      else if( f1 > rtmin .and. f1 < rtmax .and. &
               g1 > rtmin .and. g1 < rtmax ) then
         d = sqrt( f*f + g*g )
         p = one / d
         c( ic ) = f1*p
         y( iy ) = g*sign( p, f )
         x( ix ) = sign( d, f )
      else
         u = min( safmax, max( safmin, f1, g1 ) )
         uu = one / u
         fs = f*uu
         gs = g*uu
         d = sqrt( fs*fs + gs*gs )
         p = one / d
         c( ic ) = abs( fs )*p
         y( iy ) = gs*sign( p, f )
         x( ix ) = sign( d, f )*u
      end if
      ic = ic + incc
      iy = iy + incy
      ix = ix + incx
   end do
   return
end subroutine
