subroutine DROTG( a, b, c, s )
   use LA_CONSTANTS, only: wp, zero, one, safmin, safmax
!
!  Updated Level 1 BLAS
!  E. Anderson
!  July 26, 2016
!
!  .. Scalar Arguments ..
   real(wp) :: a, b, c, s
!  ..
!
!  Purpose
!  =======
!
!  DROTG constructs a plane rotation
!     [  c  s ] [ a ] = [ r ]
!     [ -s  c ] [ b ]   [ 0 ]
!  satisfying c**2 + s**2 = 1.
!
!  Further Details
!  ===============
!
!  The computation uses the formulas
!     sigma = sgn(a)    if |a| >  |b|
!           = sgn(b)    if |b| >= |a|
!     r = sigma*sqrt( a**2 + b**2 )
!     c = 1; s = 0      if r = 0
!     c = a/r; s = b/r  if r != 0
!  The subroutine also computes
!     z = s    if |a| > |b|,
!       = 1/c  if |b| >= |a| and c != 0
!       = 1    if c = 0
!  This allows c and s to be reconstructed from z as follows:
!     If z = 1, set c = 0, s = 1.
!     If |z| < 1, set c = sqrt(1 - z**2) and s = z.
!     If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).
!
!  Arguments
!  =========
!
!  A       (input) REAL
!          On entry, the scalar a.
!          On exit, the scalar r.
!
!  B       (input/output) REAL
!          On entry, the scalar b.
!          On exit, the scalar z.
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
   real(wp) :: anorm, bnorm, scl, sigma, r, z
!  ..
   anorm = abs(a)
   bnorm = abs(b)
   if( anorm == zero ) then
      c = zero
      s = one
      a = b
      b = one
   else if( bnorm == zero ) then
      c = one
      s = zero
      b = zero
   else
      scl = min( safmax, max( safmin, anorm, bnorm ) )
      if( anorm > bnorm ) then
         sigma = sign(one,a)
      else
         sigma = sign(one,b)
      end if
      r = sigma*( scl*sqrt((a/scl)**2 + (b/scl)**2) )
      c = a/r
      s = b/r
      if( anorm > bnorm ) then
         z = s
      else if( c /= zero ) then
         z = one/c
      else
         z = one
      end if
      a = r
      b = z
   end if
   return
end subroutine
