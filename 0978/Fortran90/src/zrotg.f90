subroutine ZROTG( a, b, c, s )
   use LA_CONSTANTS, only: wp, zero, one, czero, rtmin, rtmax, safmin, safmax
!
!  Updated Level 1 BLAS
!  E. Anderson
!  August 4, 2016
!
!  .. Scalar Arguments ..
   real(wp) :: c
   complex(wp) :: a, b, s
!  ..
!
!  Purpose
!  =======
!
!  ZROTG constructs a plane rotation
!     [  c         s ] [ a ] = [ r ]
!     [ -conjg(s)  c ] [ b ]   [ 0 ]
!  where c is real, s ic complex, and c**2 + conjg(s)*s = 1.
!
!  Further Details
!  ===============
!
!  The computation uses the formulas
!     |x| = sqrt( Re(x)**2 + Im(x)**2 )
!     sgn(x) = x / |x|  if x /= 0
!            = 1        if x  = 0
!     c = |a| / sqrt(|a|**2 + |b|**2)
!     s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
!  When a and b are real and r /= 0, the formulas simplify to
!     r = sgn(a)*sqrt(|a|**2 + |b|**2)
!     c = a / r
!     s = b / r
!  the same as in SROTG when |a| > |b|.  When |b| >= |a|, the
!  sign of c and s will be different from those computed by SROTG
!  if the signs of a and b are not the same.
!
!  Arguments
!  =========
!
!  A       (input/output) COMPLEX
!          On entry, the scalar a.
!          On exit, the scalar r.
!
!  B       (input) COMPLEX
!          The scalar b.
!
!  C       (output) REAL
!          The scalar c.
!
!  S       (output) REAL
!          The scalar s.
!
! =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
   complex(wp) :: f, fs, g, gs, r, t
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
   f = a
   g = b
   if( g == czero ) then
      c = one
      s = czero
      r = f
   else if( f == czero ) then
      c = zero
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         g2 = ABSSQ( g )
         d = sqrt( g2 )
         s = conjg( g ) / d
         r = d
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         d = sqrt( g2 )
         s = conjg( gs ) / d
         r = d*u
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
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
         c = f2*p
         s = conjg( g )*( f*p )
         r = f*( h2*p )
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, f1, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         if( f1*uu < rtmin ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            vv = one / v
            w = v * uu
            fs = f*vv
            f2 = ABSSQ( fs )
            h2 = f2*w**2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
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
         c = ( f2*p )*w
         s = conjg( gs )*( fs*p )
         r = ( fs*( h2*p ) )*u
      end if
   end if
   a = r
   return
end subroutine
