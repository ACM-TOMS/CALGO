subroutine ZB1TSSQ( nn, nvals, nix, ixvals, lenx, x, xs, z, &
      work, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, ulp, half, one, two, &
                             safmin, safmax, smlnum, bignum
   use LA_XISNAN
   use LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  February 17, 2017
!
!  .. Scalar Arguments ..
   integer :: nn, nix, lenx, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   integer :: ixvals(*), nvals(*)
   real(wp) :: work(*)
   complex(wp) :: x(*), xs(*), z(*)
!  ..
!
!  Purpose
!  =======
!
!  ZB1TSSQ tests the sum of squares routine ZLASSQ.  Tests are done
!  with combinations of the following values:
!
!  0, very small, small, ulp, 1, 1/ulp, big, very big, infinity, NaN
!
!  One of these values is used to initialize SCL and SUMSQ using
!  one of the following scaling schemes:
!     1) SCL = 1.0, SUMSQ = value
!     2) SCL = value, SUMSQ = 1.0
!     3) SCL = s, SUMSQ = value / s for some s
!  Then the real and complex parts of a random vector of length N are
!  filled with random values from [-1,1] and scaled by another of these
!  values.
!
!  Arguments
!  =========
!
!  NN      (input) INTEGER
!          The number of values of N contained in the vector NVALS.
!
!  NVALS   (input) INTEGER array, dimension (NN)
!          The values of the vector length N.
!
!  NIX     (input) INTEGER
!          The number of values of INCX contained in the vector IXVALS.
!
!  IXVALS  (input) INTEGER array, dimension (NIX)
!          The values of the vector increment INCX.
!
!  LENX    (input) INTEGER
!          The length of the workspace vectors. LENX >= (1+(NMAX-1)*IXMAX),
!          where NMAX is the largest value of N and IXMAX is the largest
!          value of INCX in absolute value.
!
!  X       (workspace) COMPLEX array, dimension (LENX)
!
!  XS      (workspace) COMPLEX array, dimension (LENX)
!
!  Z       (workspace) COMPLEX array, dimension (NMAX)
!
!  WORK    (workspace) REAL array, dimension (2*NMAX)
!
!  THRESH  (input) REAL
!          The threshold value for the test ratios.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!  =====================================================================
!
!  .. Parameters ..
   integer, parameter :: nv = 10
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, in, incx, incxs, is, iv, ix, iw, ixmax, j, ks, &
      n, nerrs, nfail, nmax, ns, ntests
   real(wp) :: s, scl, snrm, ssq, sumsq, trat, v0, v1, workssq, &
      y1, y2, ymax, ymin, ynrm, znrm
   complex(wp) :: rogue
!  ..
!  .. Local Arrays ..
   real(wp) :: values(nv)
!  ..
!  .. External Functions ..
   real(wp) :: DLAPY2
!  external :: DLAPY2
!  ..
!
   values(1) = zero
   values(2) = two*safmin
   values(3) = smlnum
   values(4) = ulp
   values(5) = one
   values(6) = one / ulp
   values(7) = bignum
   values(8) = safmax
   values(9) = LA_XVALS( 2, v0 )
   values(10) = LA_XVALS( 3, v0 )
   nfail = 0
   ntests = 0
   rogue = cmplx( 1234.5678_wp, -1234.5678_wp, wp )
   subnam = 'ZLASSQ'
!
!  Check that the arrays are large enough
!
   nmax = 0
   do i = 1, nn
      nmax = max( nmax, nvals(i) )
   end do
   ixmax = 0
   do i = 1, nix
      ixmax = max( ixmax, ixvals(i) )
   end do
   if( lenx < nmax*ixmax ) then
      write(nout,99) nmax, ixmax, lenx, nmax*ixmax
      return
   end if
!
!  Test for each value of n
!
   do in = 1, nn
!
!    Generate 2*n values in (-1,1).
!
     ns = max( 0, nvals(in) )
     ks = 2*ns
     if( ks > 0 ) then
       call random_number(work(1:ks))
     end if
     do i = 1, ks
       work(i) = one - two*work(i)
     end do
!
!    Compute the sum of squares of the random values
!    by an unscaled algorithm.
!
     workssq = zero
     do i = 1, ks
       workssq = workssq + work(i)*work(i)
     end do
!
!    Construct the test vector with one known value
!    and the rest from the random work array multiplied
!    by a scaling factor
!
     do iv = 1, nv
       v0 = values(iv)
       if( v0 > one ) then
         v0 = v0*half
       end if
       do is = 1, 3
         z(1) = cmplx( 0.6_wp*sqrt(v0), -0.8_wp*sqrt(v0), wp )
         if( is == 1 ) then
           s = one
           ssq = v0
         else if( is == 2 ) then
           s = sqrt(v0)
           ssq = one
         else if( v0 < one ) then
           s = one / 3.0_wp
           ssq = v0*9.0_wp
         else
           s = 3.0_wp
           ssq = v0/9.0_wp
         end if
         do iw = 1, nv
           v1 = values(iw)
           if( v1 > one ) then
             v1 = ( v1*half ) / sqrt( real(ks+1,wp) )
           end if
           do i = 1, ns
             z(i+1) = cmplx( v1*work(2*i-1), v1*work(2*i), wp )
           end do
!
!          Compute the expected value of the sum of
!          squares using the z vector of length n+1.
!
           if( ns >= 0 ) then
             y1 = sqrt(abs(v0))
           else
             y1 = zero
           end if
           if( ns >= 1 ) then
             y2 = abs(v1)*sqrt(workssq)
           else
             y2 = zero
           end if
           ymin = min( y1, y2 )
           ymax = max( y1, y2 )
!
!          Expected value is NaN if either is NaN. The test
!          for ymin == ymax avoids further computation if both
!          are infinity.
!
           if( LA_ISNAN( y1 ) .or. LA_ISNAN( y2 ) ) then
             ynrm = LA_XVALS( 3, v0 )
           else if( ymin == ymax ) then
             ynrm = sqrt(two)*ymax
           else if( ymax == zero ) then
             ynrm = zero
           else
             ynrm = ymax * sqrt( one + ( ymin / ymax )**2 )
           end if
!
!          Test for each value of incx
!
           do j = 1, nix
             xs(1:lenx) = rogue
             incxs = ixvals(j)
             ix = 1
             if( incxs < 0 ) ix = 1 - (ns-1)*incxs
             do i = 1, ns
               xs(ix) = z(i+1)
               ix = ix + incxs
             end do
!
!            Call ZLASSQ to compute the sum of squares
!
             scl = s
             sumsq = ssq
             n = ns
             x(1:lenx) = xs(1:lenx)
             incx = incxs
             call ZLASSQ( n, x, incx, scl, sumsq )
             snrm = scl * sqrt(sumsq)
!
!            See what values changed inside the subroutine
!
             nerrs = 0
             call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
             call ZARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
             call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!            Compare snrm and zrnm
!
             if( ns > 0 .and. incxs == 0 ) then
               y1 = sqrt(abs(v0))
               y2 = sqrt( real(ns,wp) )*DLAPY2( real(x(1)), aimag(x(1)) )
               ymin = min( y1, y2 )
               ymax = max( y1, y2 )
               if( LA_ISNAN( y1 ) .or. LA_ISNAN( y2 ) ) then
                 znrm = LA_XVALS( 3, v0 )
               else if( ymin == ymax ) then
                 znrm = sqrt(two)*ymax
               else if( ymax == zero ) then
                 znrm = zero
               else
                 znrm = ymax * sqrt( one + ( ymin / ymax )**2 )
               end if
             else
               znrm = ynrm
             end if
             if( LA_ISNAN( snrm ) .or. LA_ISNAN( znrm ) ) then
               if( LA_ISNAN( snrm ) .neqv. LA_ISNAN( znrm ) ) then
                 trat = one / ulp
               else
                 trat = zero
               end if
             else if( snrm == znrm ) then
               trat = zero
             else if( znrm == zero ) then
               trat = snrm / ulp
             else
               trat = ( abs( snrm - znrm ) / znrm ) / &
                      ( two*real(max(1,ns),wp)*ulp )
             end if
             if( LA_ISNAN( trat ) .or. trat >= thresh ) then
               write(nout,97) subnam, ns, incxs, iv, is, iw, trat
               nfail = nfail + 1
             end if
             ntests = ntests + 1
           end do
         end do
       end do
     end do
   end do
   if( nfail == 0 ) then
     write(nout,95) subnam, ntests
   else
     write(nout,94) subnam, nfail, ntests
   end if
   return
99 format( ' Not enough space to test ', A6, ': NMAX = ',I6, &
           ', IXMAX = ',I6,/,'   LENX = ',I6,', must be at least ',I6 )
97 format( 1X, A6, ': N=', I6,', INCX=', I4, ', IV=', I1, ', IS=', I1, &
           ', IW=', I1, ', test=', E15.8 )
95 format( 1X, A6, ' passed the computational tests (', I6, &
           ' tests run)', / )
94 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
