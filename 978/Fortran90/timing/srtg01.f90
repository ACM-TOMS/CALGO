subroutine SRTG01( likbla, f, g, cs, sn, r, result )
   use LA_CONSTANTS32, only: wp, zero, one, ulp, safmax, safmin
   use LA_XISNAN
   use LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  June 10, 2016
!
!  .. Scalar Arguments ..
   integer :: likbla
   real(wp) :: cs, sn, f, g, r
!  ..
!  .. Array Arguments ..
   real(wp) :: result(*)
!  ..
!
!  Purpose
!  =======
!
!  SRTG01 is an auxiliary test routine for SB1TRTG.  The input values
!  are components of the rotation
!     [  CS  SN ] [ F ] = [ R ]
!     [ -SN  CS ] [ G ] = [ 0 ]
!  where CS**2 + SN**2 = 1.
!
!  If f and g are finite, the identities checked are
!     1:  CS*F + SN*G = R
!     2:  -SN*F + CS*G = 0
!     3:  CS**2 + SN**2 = 1
!     4:  CS >= 0
!
!  If either f or g is NaN, the only check is that r is NaN.
!  If either f or g is infinite, the only check is that r is infinite
!  or NaN.
!
!  Arguments
!  =========
!
!  LIKBLA  (input) INTEGER
!          = 0: LAPACK-style rotation, in which CS >= 0.
!          = 1: BLAS-style rotation, in which CS >= 0 if |F| > |G|.
!
!  F       (input) REAL
!          The first component of the vector to be rotated.
!
!  G       (input) REAL
!          The second component of the vector to be rotated.
!
!  CS      (input) REAL
!          The cosine of the rotation.
!
!  SN      (input) REAL
!          The sine of the rotation.
!
!  R       (input) REAL
!          The first component of the rotated vector.
!
!  RESULT  (output) REAL array, dimension(4)
!          The computed test ratios.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: d, fs, gs, p, rnorm, rs
!  ..
!
   d = LA_XVALS( 2, d )
   if( LA_ISNAN( f ) .or. LA_ISNAN( g ) ) then
!
!     If either f or g is NaN, the only check is that r is NaN.
!
      if( LA_ISNAN( r ) ) then
         result(1) = zero
      else
         result(1) = one / ulp
      end if
      result(2) = zero
      result(3) = zero
      result(4) = zero
   else if( abs(f) == d .or. abs(g) == d ) then
!
!     If either f or g is infinite, the only check is that r is infinite
!     or NaN.
!
      if( LA_ISNAN( r ) .or. abs(r) == d ) then 
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
      d = min( safmax, max( safmin, abs(f), abs(g) ) )
      p = one / d
!
!     Test 1:  CS*F + SN*G = R
!
      fs = f*p
      gs = g*p
      rs = r*p
      rnorm = abs(rs)
      if( rnorm == zero ) then
         rnorm = max( abs(fs), abs(gs) )
         if( rnorm == zero ) then
           rnorm = one
         end if
      end if
      result(1) = ( abs( rs - ( cs*fs + sn*gs ) ) / rnorm ) / ulp
!
!     Test 2:  -SN*F + CS*G = 0
!
      result(2) = abs( zero - ( -sn*fs + cs*gs ) ) / ulp
!
!     Test 3:  CS**2 + SN**2 = 1
!
      result(3) = abs( one - ( cs*cs + sn*sn ) ) / ulp
!
!     Test 4:  CS >= 0
!
      if( likbla == 0 .or. ( abs(f) > abs(g) ) ) then
         if( cs >= zero ) then
            result(4) = zero
         else
            result(4) = one / ulp
         end if
      else
         if( cs*sign( one, f )*sign( one, g ) >= zero ) then
            result(4) = zero
         else
            result(4) = one / ulp
         end if
      end if
   end if
end subroutine
