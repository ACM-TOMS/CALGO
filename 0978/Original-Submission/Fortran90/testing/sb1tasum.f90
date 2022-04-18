subroutine SB1TASUM( nn, nvals, nix, ixvals, lenx, x, xs, z, &
      work, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, one, two, ulp, safmin, safmax, &
                             smlnum, bignum
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
   real(wp) :: x(*), xs(*), z(*), work(*)
!  ..
!
!  Purpose
!  =======
!
!  SB1TASUM tests the 1-norm routine SASUM. Tests are done with
!  combinations of the following values:
!
!  0, very small, small, ulp, 1, 1/ulp, big, very big, infinity, NaN
!
!  One of these values is used to initialize x(1) and x(2:N) is
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
!  X       (workspace) REAL array, dimension (LENX)
!
!  XS      (workspace) REAL array, dimension (LENX)
!
!  Z       (workspace) REAL array, dimension (NMAX)
!
!  WORK    (workspace) REAL array, dimension (NMAX)
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
   integer, parameter :: nv = 10
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, in, incx, incxs, iv, ix, iw, ixmax, j, n, &
      nerrs, nfail, nmax, ns, ntests
   real(wp) :: snrm, trat, v0, v1, wnrm, y1, y2, ynrm, znrm
!  ..
!  .. Local Arrays ..
   real(wp) :: values(nv)
!  ..
!  .. External Functions ..
   real(wp) :: SASUM
   external :: SASUM
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
   subnam = 'SASUM '
!
!  Check that the arrays are large enough
!
   nmax = 1
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
     ns = nvals(in)
     if( ns > 1 ) then
       call random_number(work(1:ns-1))
     end if
     do i = 1, ns-1
       work(i) = one - two*work(i)
     end do
     do iv = 1, nv
       v0 = values(iv)
       if( abs(v0) > bignum ) v0 = v0 / real(ns+1,wp)
       z(1) = v0
       do iw = 1, nv
         v1 = values(iw)
         if( abs(v1) > bignum ) v1 = v1 / real(ns+1,wp)
         wnrm = zero
         do i = 1, ns-1
           z(i+1) = v1*work(i)
           wnrm = wnrm + abs(z(i+1))
         end do
!
!        Compute the expected value of the 1-norm
!
         if( ns > 0 ) then
           y1 = abs(v0)
         else
           y1 = zero
         end if
         if( ns > 1 ) then
           y2 = wnrm
         else
           y2 = zero
         end if
         ynrm = y1 + y2
!
!        Test for each value of incx
!
         do j = 1, nix
           xs(1:lenx) = bignum
           incxs = ixvals(j)
           ix = 1
           if( incxs < 0 ) ix = 1 - (ns-1)*incxs
           do i = 1, ns
             xs(ix) = z(i)
             ix = ix + incxs
           end do
!
!          Call SASUM to compute the 1-norm
!
           n = ns
           x(1:lenx) = xs(1:lenx)
           incx = incxs
           snrm = SASUM( n, x, incx )
!
!          See what values changed inside the subroutine
!
           nerrs = 0
           call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
           call SARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
           call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!          Compare snrm and znrm
!
           if( incxs == 0 ) then
              znrm = real(ns,wp)*abs( xs(1) )
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
             trat = ( abs( snrm - znrm ) / znrm ) / ( real(max(1,ns),wp)*ulp )
           end if
           if( LA_ISNAN( trat ) .or. trat >= thresh ) then
             write(nout,98) subnam, ns, incxs, iv, iw, trat
             nerrs = nerrs + 1
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
97 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
