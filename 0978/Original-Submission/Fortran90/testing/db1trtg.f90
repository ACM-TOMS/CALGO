subroutine DB1TRTG( rtname, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, one, three, ulp, safmin, safmax
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
!  DB1TRTG tests the generation of plane rotations in the BLAS DROTG
!  and LAPACK DLARTG.  Tests are done with +/- combinations of the
!  following values:
!
!  0, very small, ulp, 1/3, 1, 1/ulp, very big, infinity, NaN
!
!  Arguments
!  =========
!
!  RTNAME  (input) CHARACTER*(*)
!          The subroutine to test, either DROTG or DLARTG.
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
   integer :: i, j, k, likbla, nfail, nrun
   real(wp) :: cs, sn, f, g, r, z
!  ..
!  .. Local Arrays ..
   real(wp) :: tstval(nvals), result(ntests)
!  ..
!  .. External Functions ..
   logical :: LSAMEN
   external :: LSAMEN
!  ..
!
!  Check for a valid routine name
!
   if( LSAMEN( 5, rtname, 'DROTG' ) ) then
      likbla = 1
   else if( LSAMEN( 6, rtname, 'DLARTG' ) ) then
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
   tstval(14) = LA_XVALS( 1, f )
   tstval(15) = LA_XVALS( 2, f )
   tstval(16) = LA_XVALS( 3, f )
!
!  Do for all pairs of these values
!
   nfail = 0
   nrun = 0
   do i = 1, nvals
      f = tstval(i)
      do j = 1, nvals
         g = tstval(j)
         if( LSAMEN( 5, rtname, 'DROTG' ) ) then
            r = f
            z = g
            call DROTG( r, z, cs, sn )
         else
            call DLARTG( f, g, cs, sn, r )
         end if 

!        Test the resulting values of CS, SN, and R.
!
         call DRTG01( likbla, f, g, cs, sn, r, result )
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
               write( nout, 98 ) rtname, f, g, cs, sn, k, result(k)
               nfail = nfail + 1
            end if
         end do
         nrun = nrun + 4
      end do
   end do
   if( nfail == 0 ) then
      WRITE(nout,97) rtname, nrun
   else
      write(nout,96) rtname, nfail, nrun
   end if
99 format( /, 'Tests for DROTG/DLARTG:', / &
           '   1:  C*F + S*G   = R', / &
           '   2:  -S*F + C*G  = 0', / &
           '   3:  C**2 + S**2 = 1', / &
           '   4:  C >= 0',/ &
           'Tests pass the threshold if <', G13.5, / )
98 format( 1X, A6, ': F = ', G13.5, ', G = ', G13.5, ', C = ', G13.5, ', S = ', &
           G13.5, ', test(', I1, ') =', G13.5 )
97 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6, ' tests failed ', &
           'to pass the threshold', / )
end subroutine
