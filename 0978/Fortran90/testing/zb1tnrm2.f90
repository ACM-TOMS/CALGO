subroutine ZB1TNRM2( nn, nvals, nix, ixvals, lenx, x, xs, z, &
      work, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, ulp, half, one, two, three, four, &
                             safmin, safmax, smlnum, bignum
   use LA_XISNAN
   use LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  February 16, 2017
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
!  ZB1TNRM2 tests the 2-norm routine DZNRM2.  Tests are done with
!  combinations of the following values:
!
!  0, very small, small, ulp, 1, 1/ulp, big, very big, infinity, NaN
!
!  One of these values is used to initialize the real and imaginary
!  parts of x(1) and the real and imaginary parts of x(2:N) are
!  filled with random values from [-1,1] scaled by another of
!  these values.
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
!  WORK    (workspace) REAL array, dimension (NMAX)
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
   integer :: i, in, incx, incxs, iv, ix, iw, ixmax, j, ks, &
      n, nerrs, nfail, nm1, nmax, ns, ntests
   real(wp) :: snrm, trat, workssq, v0, v1, y1, y2, ymax, ymin, ynrm, znrm
   complex(wp) :: rogue
!  ..
!  .. Local Arrays ..
   real(wp) :: values(nv)
!  ..
!  .. External Functions ..
   real(wp) :: DZNRM2
   external :: DZNRM2
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
   subnam = 'DZNRM2'
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
!    Generate 2*(n-1) values in (-1,1).
!
     ns = nvals(in)
     nm1 = max( 0, ns - 1 )
     ks = 2*nm1
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
       if( abs(v0) > one ) then
          v0 = v0 / four
       end if
       z(1) = cmplx( v0, -three*v0, wp )
       do iw = 1, nv
         v1 = values(iw)
         if( abs(v1) > one ) then
            v1 = ( v1*half ) / sqrt( real(ks+1,wp) )
         end if
         do i = 1, nm1
           z(i+1) = cmplx( v1*work(2*i-1), v1*work(2*i), wp )
         end do
!
!        Compute the expected value of the 2-norm
!
         if( ns > 0 ) then
           y1 = abs(v0)*sqrt(10.0_wp)
         else
           y1 = zero
         end if
         if( ns > 1 ) then
           y2 = abs(v1)*sqrt(workssq)
         else
           y2 = zero
         end if
         ymin = min( y1, y2 )
         ymax = max( y1, y2 )
!
!        Expected value is NaN if either is NaN. The test
!        for ymin == ymax avoids further computation if both
!        are infinity.
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
!        Test for each value of incx
!
         do j = 1, nix
           xs(1:lenx) = rogue
           incxs = ixvals(j)
           ix = 1
           if( incxs < 0 ) ix = 1 - (ns-1)*incxs
           do i = 1, ns
             xs(ix) = z(i)
             ix = ix + incxs
           end do
!
!          Call DZNRM2 to compute the 2-norm
!
           n = ns
           x(1:lenx) = xs(1:lenx)
           incx = incxs
           snrm = DZNRM2( n, x, incx )
!
!          See what values changed inside the subroutine
!
           nerrs = 0
           call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
           call ZARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
           call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!          Compare snrm and zrnm.  Roundoff error grows like O(n)
!          in this implementation so we scale the test ratio accordingly.
!
           if( incxs == 0 ) then
              y1 = abs( real(x(1)) )
              y2 = abs( aimag(x(1)) )
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
              znrm = sqrt( real(ns,wp) )*znrm
           else
              znrm = ynrm
           end if
           if( LA_ISNAN( snrm ) .or. LA_ISNAN( znrm ) ) then
             if( LA_ISNAN( snrm ) .neqv. LA_ISNAN( znrm ) ) then
                trat = one / ulp
             else
                trat = zero
             end if
           else if( znrm == zero ) then
             trat = snrm / ulp
           else
             trat = ( abs( snrm - znrm ) / znrm ) / ( two*real(max(1,ns),wp)*ulp )
           end if
           if (trat >= thresh) then
             write(nout,98) subnam, ns, incxs, iv, iw, trat
           end if
           if( nerrs > 0 ) nfail = nfail + 1
           ntests = ntests + 1
         end do
       end do
     end do
   end do
   if( nfail == 0 ) then
     write(nout,97) subnam, ntests
   else
     write(nout,96) subnam, nfail, ntests
   end if
   return
99 format( ' Not enough space to test ', A6, ': NMAX = ',I6, &
           ', IXMAX = ',I6,/,'   LENX = ',I6,', must be at least ',I6 )
98 format( 1X, A6, ': N=', I6,', INCX=', I4, ', IV=', I2, ', IW=', I2, &
           ', test=', E15.8 )
97 format( 1X, A6, ' passed the computational tests (', I6, &
           ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
