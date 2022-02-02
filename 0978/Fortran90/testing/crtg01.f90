subroutine CRTG01( f, g, cs, sn, r, result )
   use LA_CONSTANTS32, only: wp, zero, one, ulp, safmin, safmax, czero, cone
   use LA_XISNAN
   use LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  May 26, 2016
!
!  .. Scalar Arguments ..
   real(wp) :: cs
   complex(wp) :: sn, f, g, r
!  ..
!  .. Array Arguments ..
   real(wp) :: result(*)
!  ..
!
!  Purpose
!  =======
!
!  CRTG01 is an auxiliary test routine for CB1TRTG.  The input values
!  are components of the rotation
!     [     CS      SN ] [ F ] = [ R ]
!     [ -conjg(SN)  CS ] [ G ] = [ 0 ]
!  where CS is real, SN is complex, and CS**2 + conjg(SN)*SN = 1.
!
!  If f and g are finite, the identities checked are
!     1:  CS*F + SN*G = R
!     2:  -conjg(SN)*F + CS*G = 0
!     3:  CS**2 + conjg(SN)*SN = 1
!     4:  CS >= 0
!
!  If either f or g is NaN, the only check is that r is NaN.
!  If either f or g is infinite, the only check is that r is infinite
!  or NaN.
!
!  Arguments
!  =========
!
!  F       (input) COMPLEX
!          The first component of the vector to be rotated.
!
!  G       (input) COMPLEX
!          The second component of the vector to be rotated.
!
!  CS      (input) REAL
!          The cosine of the rotation.
!
!  SN      (input) COMPLEX
!          The sine of the rotation.
!
!  R       (input) COMPLEX
!          The first component of the rotated vector.
!
!  RESULT  (output) REAL array, dimension(4)
!          The computed test ratios.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: d, fi, fmax, fr, gi, gmax, gr, p, rr, ri, rnorm
   complex(wp) :: fs, gs, rs
!  ..

   d = LA_XVALS( 2, d )
   fr = real(f)
   fi = aimag(f)
   gr = real(g)
   gi = aimag(g)
   rr = real(r)
   ri = aimag(r)
   if( LA_ISNAN( fr ) .or. LA_ISNAN( fi ) .or. &
       LA_ISNAN( gr ) .or. LA_ISNAN( gi ) ) then
!
!     If either f or g is NaN, the only check is that r is NaN.
!
      if( LA_ISNAN( rr ) .or. LA_ISNAN( ri ) ) then
         result(1) = zero
      else
         result(1) = one / ulp
      end if
      result(2) = zero
      result(3) = zero
      result(4) = zero
   else if( abs(fr) == d .or. abs(fi) == d .or. &
            abs(gr) == d .or. abs(gi) == d ) then
!
!     If either f or g is infinite, the only check is that r is infinite
!     or NaN.
!
      if( LA_ISNAN( rr ) .or. LA_ISNAN( ri ) .or. abs(rr) == d .or. abs(ri) == d ) then
         result(1) = zero
      else
         result(1) = one / ulp
      end if
      result(2) = zero
      result(3) = zero
      result(4) = zero
   else
!
!     Otherwise, check that the rotation solves the equation above.
!
      fmax = max( abs(fr), abs(fi) )
      gmax = max( abs(gr), abs(gi) )
      d = min( safmax, max( safmin, fmax, gmax ) )
      p = one / d
!
!     Test 1:  CS*F + SN*G = R
!
      rs = r*p
      fs = f*p
      gs = g*p
      rnorm = abs(rs)
      if( rnorm == zero ) then
         rnorm = max( abs(fs), abs(gs) )
         if( rnorm == zero ) then
           rnorm = one
         end if
      end if
      result(1) = ( abs( rs - ( cs*fs + sn*gs ) ) / rnorm ) / ulp
!
!     Test 2:  -conjg(SN)*F + CS*G = 0
!
      rnorm = max( abs(fs), abs(gs) )
      if( rnorm == zero ) then
         rnorm = one
      end if
      result(2) = ( abs( czero - ( -conjg(sn)*fs + cs*gs ) ) / rnorm ) / ulp
!
!     Test 3:  CS**2 + conjg(SN)*SN = 1
!
      result(3) = abs( cone - ( cs*cs + conjg(sn)*sn ) ) / ulp
!
!     Test 4:  CS >= 0
!
      if( cs >= zero ) then
         result(4) = zero
      else
         result(4) = one / ulp
      end if
   end if
end subroutine
