subroutine CLARGV( n, x, incx, y, incy, c, incc )
   use LA_CONSTANTS32, only: wp, zero, one, safmin, safmax, rtmin, rtmax, czero
!
!  LAPACK auxiliary routine
!  E. Anderson
!  August 4, 2016
!
!  .. Scalar Arguments ..
   integer            incc, incx, incy, n
!     ..
!  .. Array Arguments ..
   real(wp)           c( * )
   complex(wp)        x( * ), y( * )
!  ..
!
!  Purpose
!  =======
!
!  CLARGV generates a vector of complex plane rotations with real
!  cosines, determined by elements of the complex vectors x and y.
!  For i = 1,2,...,n
!
!     (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
!     ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
!
!  where c**2 + conjg(s)*s = 1.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of plane rotations to be generated.
!
!  X       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)
!          On entry, the vector x.
!          On exit, x(i) is overwritten by r(i), for i = 1,...,n.
!
!  INCX    (input) INTEGER
!          The increment between elements of X.
!
!  Y       (input/output) COMPLEX array, dimension (1+(N-1)*INCY)
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
   real(wp)           d, f1, g1, f2, g2, h2, p, u, uu, v, vv, w
   complex(wp)        f, g, fs, gs, t
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, conjg, max, min, real, sqrt
!  ..
!  .. Statement Functions ..
   real(wp) :: ABSSQ
!  ..
!  .. Statement Function definitions ..
   ABSSQ( t ) = real( t )**2 + aimag( t )**2
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
!     Use identical algorithm to CLARTG
!
      f = x( ix )
      g = y( iy )
      if( g == czero ) then
         c( ic ) = one
         x( ix ) = f
      else if( f == czero ) then
         c( ic ) = zero
         g1 = max( abs(real(g)), abs(aimag(g)) )
         if( g1 > rtmin .and. g1 < rtmax ) then
!
!           Use unscaled algorithm
!
            g2 = ABSSQ( g )
            d = sqrt( g2 )
            y( iy ) = conjg( g ) / d
            x( ix ) = d
         else
!
!        Use scaled algorithm
!
            u = min( safmax, max( safmin, g1 ) )
            uu = one / u
            gs = g*uu
            g2 = ABSSQ( gs )
            d = sqrt( g2 )
            y( iy ) = conjg( gs ) / d
            x( ix ) = d*u
         end if
      else
         f1 = max( abs(real(f)), abs(aimag(f)) )
         g1 = max( abs(real(g)), abs(aimag(g)) )
         if( f1 > rtmin .and. f1 < rtmax .and. &
             g1 > rtmin .and. g1 < rtmax ) then
!
!           Use unscaled algorithm
!
            f2 = ABSSQ( f )
            g2 = ABSSQ( g )
            h2 = f2 + g2
            if( f2 > rtmin .and. h2 < rtmax ) then
               d = sqrt( f2*h2 )
            else
               d = sqrt( f2 )*sqrt( h2 )
            end if
            p = 1 / d
            c( ic ) = f2*p
            y( iy ) = conjg( g )*( f*p )
            x( ix ) = f*( h2*p )
         else
!
!           Use scaled algorithm
!
            u = min( safmax, max( safmin, f1, g1 ) )
            uu = one / u
            gs = g*uu
            g2 = ABSSQ( gs )
            if( f1*uu < rtmin ) then
!
!              f is not well-scaled when scaled by g1.
!              Use a different scaling for f.
!
               v = min( safmax, max( safmin, f1 ) )
               vv = one / v
               w = v * uu
               fs = f*vv
               f2 = ABSSQ( fs )
               h2 = f2*w**2 + g2
            else
!
!              Otherwise use the same scaling for f and g.
!
               w = one
               fs = f*uu
               f2 = ABSSQ( fs )
               h2 = f2 + g2
            end if
            if( f2 > rtmin .and. h2 < rtmax ) then
               d = sqrt( f2*h2 )
            else
               d = sqrt( f2 )*sqrt( h2 )
            end if
            p = 1 / d
            c( ic ) = ( f2*p )*w
            y( iy ) = conjg( gs )*( fs*p )
            x( ix ) = ( fs*( h2*p ) )*u
         end if
      end if
      ix = ix + incx
      iy = iy + incy
      ic = ic + incc
   end do
   return
end subroutine
