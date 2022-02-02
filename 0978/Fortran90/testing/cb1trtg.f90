subroutine CB1TRTG( rtname, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, one, three, ulp, safmin, safmax
   use LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  May 26, 2016
!
!  .. Scalar Arguments ..
   character*(*) :: rtname
   integer :: nout
   real(wp) :: thresh
!  ..
!
!  Purpose
!  =======
!
!  CB1TRTG tests the generation of plane rotations in LINPACK CROTG
!  and LAPACK CLARTG.  Tests are done with +/- combinations of the
!  following values:
!
!  0, very small, ulp, 1/3, 1, 1/ulp, very big, infinity, NaN
!
!  Arguments
!  =========
!
!  RTNAME  (input) CHARACTER*(*)
!          The subroutine to test, either CROTG or CLARTG.
!
!  THRESH  (input) REAL
!          The threshold value for the test ratios.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
!  .. Parameters ..
   integer, parameter :: nvals = 16, ntests = 4
!  ..
!  .. Local Scalars ..
   integer :: i1, i2, j1, j2, k, likbla, nfail, nrun
   real(wp) :: cs
   complex(wp) :: sn, f, g, r, z
!  ..
!  .. Local Arrays ..
   real(wp) :: tstval(nvals), result(ntests)
!  ..
!  .. External Functions ..
   logical :: LSAMEN
   external :: LSAMEN
!  ..
!
!  Test for a valid routine name
!
   if( LSAMEN( 5, rtname, 'CROTG' ) ) then
      likbla = 1
   else if( LSAMEN( 6, rtname, 'CLARTG' ) ) then
      likbla = 0
   else
      return
   end if
!
!  Values to use for F and G
!
   tstval(1) = -safmax
   tstval(2) = -( one / ulp )
   tstval(3) = -one
   tstval(4) = -( one / three )
   tstval(5) = -ulp
   tstval(6) = -safmin
   tstval(7) = zero
   tstval(8) = safmin
   tstval(9) = ulp
   tstval(10) = ( one / three )
   tstval(11) = one
   tstval(12) = ( one / ulp )
   tstval(13) = safmax
   tstval(14) = LA_XVALS( 1, cs )
   tstval(15) = LA_XVALS( 2, cs )
   tstval(16) = LA_XVALS( 3, cs )
!
!  Do for all pairs of these values
!
   nfail = 0
   nrun = 0
   do i1 = 1, nvals
   do i2 = 1, nvals
      f = cmplx( tstval(i1), tstval(i2), wp )
      do j1 = 1, nvals
      do j2 = 1, nvals
         g = cmplx( tstval(j1), tstval(j2), wp )
         if( LSAMEN( 5, rtname, 'CROTG' ) ) then
            r = f
            z = g
            call CROTG( r, z, cs, sn )
         else
            call CLARTG( f, g, cs, sn, r )
         end if 

!        Test the resulting values of CS, SN, and R.
!
         call CRTG01( f, g, cs, sn, r, result )
!
!        Compare tests to the threshold
!
         do k = 1, ntests
            if( result( k ).lt.thresh ) then
!
!              Do the test this way to detect NaNs
!
            else
               if( nfail == 0 ) then
                  write( nout, 99 ) thresh
               END IF
               write( nout, 98 ) rtname, real(f), aimag(f), &
                  real(g), aimag(g), k, result(k)
               nfail = nfail + 1
            end if
         end do
         nrun = nrun + 4
      end do
      end do
   end do
   end do
   if( nfail == 0 ) then
      WRITE(nout,97) rtname, nrun
   else
      write(nout,96) rtname, nfail, nrun
   end if
99 format( /, 'Tests for CROTG/CLARTG:', / &
           '   1:  C*F + S*G   = R', / &
           '   2:  -conjg(S)*F + C*G  = 0', / &
           '   3:  C**2 + conjg(S)*S = 1', / &
           '   4:  C >= 0',/ &
           'Tests pass the threshold if <', G13.5, / )
98 FORMAT( 1X, A6, ': F = (', G13.5, ',', G13.5, '), G = (', G13.5, &
           ',', G13.5, '), test(', I1, ') =', G13.5 )
97 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6, ' tests failed ', &
           'to pass the threshold', / )
end subroutine
