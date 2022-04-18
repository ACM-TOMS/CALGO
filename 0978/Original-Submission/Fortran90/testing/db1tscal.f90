subroutine DB1TSCAL( nn, nvals, nix, ixvals, nal, alvals, lenx, x, xs, &
      work, iwork, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, one, two, ulp, bignum
!
!  Level 1 BLAS test program
!  E. Anderson
!  February 15, 2017
!
!  .. Scalar Arguments ..
   integer :: nn, nix, nal, lenx, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   integer :: ixvals(*), nvals(*), iwork(*)
   real(wp) :: alvals(*), x(*), xs(*), work(*)
!  ..
!
!  Purpose
!  =======
!
!  DB1TSCAL tests the scaling routine DSCAL.
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
!  NAL     (input) INTEGER
!          The number of values of ALPHA contained in the vector ALVALS.
!
!  ALVALS  (input) REAL array, dimension (NAL)
!          The values of the scaling constant ALPHA.
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
!  WORK    (workspace) REAL array, dimension (NMAX)
!
!  IWORK   (workspace) INTEGER array, dimension (LENX)
!
!  THRESH  (input) REAL
!          The threshold value for the test ratios.
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, in, incx, incxs, ix, ixmax, j, k, n, &
      nerrs, nfail, nmax, ns, ntests
   real(wp) :: alpha, alphas, rogue, tmax, trat, xexp, xmax, &
      xmex, xnrm
!  ..
!  .. External Subroutines ..
   external :: DSCAL
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'DSCAL '
   rogue = bignum
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
     ns = nvals(in)
     call random_number(work(1:ns))
     do i = 1, ns
       work(i) = one - two*work(i)
     end do
!
!    Test for each value of incx
!
     do j = 1, nix
       incxs = ixvals(j)
       xs(1:lenx) = rogue
       iwork(1:lenx) = 0
       ix = 1
       if( incxs < 0 ) ix = 1 - (ns-1)*incxs
       do i = 1, ns
         iwork(ix) = 1
         xs(ix) = work(i)
         ix = ix + incxs
       end do
!
!      Test for each value of alpha
!
       do k = 1, nal
         alphas = alvals(k)
!
!        Call DSCAL to scale the vector x
!
         n = ns
         alpha = alphas
         x(1:lenx) = xs(1:lenx)
         incx = incxs
         call DSCAL( n, alpha, x, incx )
!
!        See what values changed inside the subroutine
!
         nerrs = 0
         call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
         call DARGCKS( nerrs, subnam, 'ALPHA', alpha, alphas, nout )
         call DARGCKX( nerrs, subnam, 'X', lenx, iwork, x, xs, nout )
         call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!        Check that the values of x were scaled correctly
!
         ix = 1
         if( incxs < 0 ) ix = 1 - (ns-1)*incxs
         if( ns > 0 .and. incxs == 0 ) then
           xexp = (alphas**ns)*work(ns)
           xnrm = abs( xexp )
           if( xnrm == zero ) then
             trat = abs( x(ix) ) / ulp
           else
             trat = ( abs( x(ix) - xexp ) / xnrm ) / ulp
           end if
           imax = 1
           xmax = x(ix)
           xmex = xexp
           tmax = trat
         else
           imax = 1
           xmax = zero
           tmax = -one
           do i = 1, ns
             xexp = alphas*work(i)
             xnrm = abs( xexp )
             if( xnrm == zero ) then
               trat = abs( x(ix) ) / ulp
             else
               trat = ( abs( x(ix) - xexp ) / xnrm ) / ulp
             end if
             if( trat > tmax ) then
               imax = i
               xmax = x(ix)
               xmex = xexp
               tmax = trat
             end if
             ix = ix + incxs
           end do
         end if
         if( tmax >= thresh ) then
           write(nout,98) subnam, ns, alphas, incxs, trat
           write(nout,97) imax, xmax, xmex
           nerrs = nerrs + 1
         end if
         if( nerrs > 0 ) nfail = nfail + 1
         ntests = ntests + 1
       end do
     end do
   end do
   if( nfail == 0 ) then
     write(nout,96) subnam, ntests
   else
     write(nout,95) subnam, nfail, ntests
   end if
   return
99 format( ' Not enough space to test ', A6, ': NMAX = ',I6, &
           ', IXMAX = ',I6,/,'   LENX = ',I6,', must be at least ',I6 )
98 format( 1X, A6, ': N=', I6,', ALPHA=', E15.8, ', INCX=', I4, &
           ', test=', E15.8 )
97 format( '   occurs at X(', I6, ')=', E15.8, ', expected=', E15.8 )
96 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
