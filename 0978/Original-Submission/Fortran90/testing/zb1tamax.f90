subroutine ZB1TAMAX( nn, nvals, nix, ixvals, lenx, x, xs, &
      work, nout )
   use LA_CONSTANTS, only: wp, zero, one, two, bignum
   use LA_XXVALS
!
!  Level 1 BLAS test program
!  E. Anderson
!  February 16, 2017
!
!  .. Scalar Arguments ..
   integer :: nn, nix, lenx, nout
!  ..
!  .. Array Arguments ..
   integer :: ixvals(*), nvals(*)
   real(wp) :: work(*)
   complex(wp) :: x(*), xs(*)
!  ..
!
!  Purpose
!  =======
!
!  ZB1TAMAX tests the infinity-norm routine IZAMAX. Tests are done
!  with the largest element in absolute value appearing in positions
!
!  1, 2, n/3, 2n/3, n-1, n
!
!  and the largest element being both positive and negative.
!  The real and imaginary parts of the other values of the vector are
!  set to random values in [-1,1] or to alternating NaN and random
!  values.
!
!  Arguments
!  =========
!
!  NN      (input) INTGER
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
!          The length of the workspace vectors. LENX >= (1+(NMAX-1)*IMAX),
!          where NMAX is the largest value of N and IXMAX is the largest
!          value of INCX in absolute value.
!
!  X       (workspace) COMPLEX array, dimension (LENX)
!
!  XS      (workspace) COMPLEX array, dimension (LENX)
!
!  WORK    (workspace) REAL array, dimension (2*NMAX)
!
!  NOUT    (input) INTEGER
!          The unit number for output.
!
!  =====================================================================
!
!  .. Parameters ..
   integer, parameter :: nvmax = 8
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, ibla, iexp, imax, imax2, in, incx, incxs, iu, iv, &
      iw, ix, ixmax, j, k, n, nerrs, nfail, nmax, ns, ntests, nv, nw
   real(wp) :: pnan, rogue
   complex(wp) :: xmax
!  ..
!  .. Local Arrays ..
   integer :: ivalues(nvmax)
!  ..
!  .. External Functions ..
   integer :: IZAMAX
   external :: IZAMAX
!  ..
!
   nfail = 0
   ntests = 0
   rogue = cmplx( 666.666_wp, -666.666_wp, wp )
   subnam = 'IZAMAX'
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
     call random_number(work(1:2*ns))
     do i = 1, 2*ns
       work(i) = one - two*work(i)
     end do
     nv = 1
     ivalues(nv) = 1
     if( ns >= 2 ) then
       nv = nv + 1
       ivalues(nv) = 2
     end if
     if( ns >= 9 ) then
        nv = nv + 1
        ivalues(nv) = ns/3
     end if
     if( ns >= 5 ) then
        nv = nv + 1
        ivalues(nv) = (2*ns)/3
     end if
     if( ns >= 4 ) then
       nv = nv + 1
       ivalues(nv) = ns-1
     end if
     if( ns >= 3 ) then
       nv = nv + 1
       ivalues(nv) = ns
     end if
!
!    Test with filler values set from the work array or
!    alternating between NaN and work values.
!
     do iu = 1, 2
!
!      Test each position in ivalues
!
       do iv = 1, nv
         imax = ivalues(iv)
!
!        Test for xmax positive or negative
!
         do j = 1, 2
           if( j == 1 ) then
             xmax = cmplx( bignum, bignum, wp )
           else
             xmax = cmplx( -bignum, -bignum, wp )
           end if
!
!          Do up to 3 additional tests with a second maximum to see
!          if the routine finds the first one
!
           nw = 1
           if( ns >= imax+1 ) nw = nw + 1
           if( ns >= imax+2 ) nw = nw + 1
           if( ns >= imax+3 ) nw = nw + 1
           do iw = 1, nw
             if( iw == 1 ) then
               imax2 = imax
             else if( iw == 2 ) then
               imax2 = ns
             else if( iw == 3 ) then
               imax2 = imax + 1
             else if( iw == 3 ) then
               imax2 = imax + 2
             end if
!  
!            Test for each value of incx
!
             do k = 1, nix
               incxs = ixvals(k)
!  
!              Fill x from work or set to NaN except for imax and imax2
!
               xs(1:lenx) = rogue
               ix = 1
               if( incxs < 0 ) ix = 1 - (ns-1)*incxs
               do i = 1, ns
                 if( i == imax .or. i == imax2 ) then
                   xs(ix) = xmax
                 else if( iu == 2 .and. mod(i,2) == 1 ) then
                   pnan = LA_XVALS( 3, pnan )
                   xs(ix) = cmplx( pnan, pnan, wp )
                 else
                   xs(ix) = cmplx( work(2*i-1), work(2*i), wp )
                 end if
                 ix = ix + incxs
               end do
!
!              Call IZAMAX to find the maximum element
!
               n = ns
               x(1:lenx) = xs(1:lenx)
               incx = incxs
               ibla = IZAMAX( n, x, incx )
!
!              See what values changed inside the subroutine
!
               nerrs = 0
               call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
               call ZARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
               call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!              Check that we found the right maximum
!
               if( ns == 0 ) then
                 iexp = 0
               else
                 if( incxs == 0 ) then
                   iexp = 1
                 else
                   iexp = imax
                 end if
               end if
               if( ibla /= iexp ) then
                 write(nout,98) subnam, ns, incxs, iu, iv, iw, ibla, iexp
               end if
               if( nerrs > 0 ) nfail = nfail + 1
               ntests = ntests + 1
             end do
           end do
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
98 format( 1X, A6, ': N=', I6,', INCX=', I4, ', IU=', I1, ', IV=', I1, &
           ', IW=', I1, ', IMAX=', I6, ', expected=', I6 )
97 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
