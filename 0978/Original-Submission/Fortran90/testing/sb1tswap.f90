subroutine SB1TSWAP( nn, nvals, nix, ixvals, lenx, x, y, xs, ys, &
                      work, iwork1, iwork2, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, one, ulp
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
   integer :: ixvals(*), nvals(*), iwork1(*), iwork2(*)
   real(wp) :: x(*), xs(*), y(*), ys(*), work(*)
!  ..
!
!  Purpose
!  =======
!
!  SB1TSWAP tests the vector swap routine SSWAP.
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
!  Y       (workspace) REAL array, dimension (LENX)
!
!  XS      (workspace) REAL array, dimension (LENX)
!
!  YS      (workspace) REAL array, dimension (LENX)
!
!  WORK    (workspace) REAL array, dimension (NMAX)
!
!  IWORK1  (workspace) INTEGER array, dimension (LENX)
!
!  IWORK2  (workspace) INTEGER array, dimension (LENX)
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
   integer :: i, imax, in, incx, incxs, incy, incys, ix, ixmax, iy, &
      j, jmax, k, n, nerrs, nfail, nmax, ns, ntests
   real(wp) :: rogue1, rogue2, tmax, trat, umax, urat, xexp, xmax, &
      xmex, xnrm, yexp, ymax, ymex, ynrm
!  ..
!
   nfail = 0
   ntests = 0
   rogue1 = -666.666_wp
   rogue2 = -999.999_wp
   subnam = 'SSWAP '
!
!  Check that the x and y arrays are large enough
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
!
!    Noncomputational routine, set x(i) = real(i) and
!    y(i) = -real(i)
!
     do i = 1, ns
       work(i) = real(i)
     end do
!
!    Test for each value of incx
!
     do j = 1, nix
       incxs = ixvals(j)
       xs(1:lenx) = rogue1
       iwork1(1:lenx) = 0
       ix = 1
       if( incxs < 0 ) ix = 1 - (ns-1)*incxs
       do i = 1, ns
         xs(ix) = work(i)
         iwork1(ix) = 1
         ix = ix + incxs
       end do
!
!      Test for each value of incy
!
       do k = 1, nix
         incys = ixvals(k)
         ys(1:lenx) = rogue2
         iwork2(1:lenx) = 0
         iy = 1
         if( incys < 0 ) iy = 1 - (ns-1)*incys
         do i = 1, ns
           ys(iy) = -work(i)
           iwork2(iy) = 1
           iy = iy + incys
         end do
         n = ns
         x(1:lenx) = xs(1:lenx)
         incx = incxs
         y(1:lenx) = ys(1:lenx)
         incy = incys
         call SSWAP( n, x, incx, y, incy )
!
!        See what values changed inside the subroutine
!
         nerrs = 0
         call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
         call SARGCKX( nerrs, subnam, 'X', lenx, iwork1, x, xs, nout )
         call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
         call SARGCKX( nerrs, subnam, 'Y', lenx, iwork2, y, ys, nout )
         call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!        Check that the expected values are in x and y
!
         imax = 1
         jmax = 1
         xmax = zero
         ymax = zero
         tmax = -one
         umax = -one
         ix = 1
         if( incxs < 0 ) ix = 1 - (ns-1)*incxs
         iy = 1
         if( incys < 0 ) iy = 1 - (ns-1)*incys
         do i = 1, ns
!
!          Determine what should be in x(i) and y(i).
!          The possibility of zero increments complicates this.
!
           if( incxs /= 0 ) then
             if( incys /= 0 ) then
               xexp = -work(i)
               yexp = work(i)
             else
               if( i == 1 ) then
                 xexp = -work(ns)
               else
                 xexp = work(i-1)
               end if
               yexp = work(ns)
             end if
           else
             if( incys /= 0 ) then
               xexp = -work(ns)
               if( i == 1 ) then
                 yexp = work(ns)
               else
                 yexp = -work(i-1)
               end if
             else
               if( mod( n, 2 ) == 0 ) then
                 xexp = work(ns)
                 yexp = -work(ns)
               else
                 xexp = -work(ns)
                 yexp = work(ns)
               end if
             end if
           end if
!
!          Check that the expected values are in x
!
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
!
!          Check that the expected values are in y
!
           ynrm = abs( yexp )
           if( ynrm == zero ) then
             urat = abs( y(iy) ) / ulp
           else
             urat = ( abs( y(iy) - yexp ) / ynrm ) / ulp
           end if
           if( urat > umax ) then
             jmax = i
             ymax = y(iy)
             ymex = yexp
             umax = urat
           end if
           iy = iy + incys
         end do
         if( max(tmax,umax) >= thresh ) then
           write(nout,98) subnam, ns, incxs, incys, max(tmax,umax)
           write(nout,97) imax, xmax, xmex, jmax, ymax, ymex
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
!
99 format( ' Not enough space to test ', A6, ': NMAX = ',I6, &
           ', IXMAX = ',I6,/,'   LENX = ',I6,', must be at least ',I6 )
98 format( 1X, A6, ': N=', I6, ', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at X(', I6, ')=', E15.8, ', expected ', E15.8, /, &
           '   and/or at Y(', I6, ')=', E15.8, ', expected ', E15.8 )
96 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
