subroutine CB1TSROT( nn, nvals, nix, ixvals, lenx, x, y, xs, ys, &
                     work1, work2, iwork1, iwork2, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, half, one, two, three, ulp
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
   integer :: ixvals(*), nvals(*), iwork1(*), iwork2(*)
   real(wp) :: work1(*), work2(*)
   complex(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  Purpose
!  =======
!
!  CB1TSROT tests the Level 1 BLAS subroutine CSROT to apply
!  a real plane rotation to a pair of complex vectors x and y.
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
!  LENX    (input) INTEGER
!          The length of the workspace vectors. LENX >= (1+(NMAX-1)*IXMAX),
!          where NMAX is the largest value of N and IXMAX is the largest
!          value of INCX in absolute value.
!
!  X       (workspace) COMPLEX array, dimension (LENX)
!
!  Y       (workspace) COMPLEX array, dimension (LENX)
!
!  XS      (workspace) COMPLEX array, dimension (LENX)
!
!  YS      (workspace) COMPLEX array, dimension (LENX)
!
!  WORK1   (workspace) REAL array, dimension (2*NMAX)
!
!  WORK2   (workspace) REAL array, dimension (2*NMAX)
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
!  .. Parameters ..
   integer, parameter :: nr = 8
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, in, incx, incxs, incy, incys, ir, ix, &
      ixmax, iy, j, jmax, k, n, nerrs, nfail, nmax, ns, ntests
   real(wp) :: c, cs, s, ss, tmax, trat, umax, urat, xnrm, ynrm
   complex(wp) :: rogue1, rogue2, xexp, xmax, xmex, xtmp, yexp, &
      ymax, ymex, ytmp
!  ..
!  .. Local Arrays ..
   real(wp) :: cvals(nr), svals(nr)
!  ..
!
   nfail = 0
   ntests = 0
   rogue1 = cmplx( 666.666_wp, -666.666_wp, wp )
   rogue2 = cmplx( 999.999_wp, -999.999_wp, wp )
   subnam = 'CSROT '
!
!  Test for each of the following rotations:
!     1) c = 0, s = 1
!     2) c = 0, s = -1
!     3) c = 1/2, s = sqrt(3)/2
!     4) c = 1/2, s = -sqrt(3)/2
!     5) c = -sqrt(3)/2, s = 1/2
!     6) c = -sqrt(3)/2, s = -1/2
!     7) c = 1, s = 0
!     8) c = -1, s = 0
!
   cvals(1) = zero; svals(1) = one
   cvals(2) = zero; svals(2) = -one
   cvals(3) = half; svals(3) = sqrt(three)*half
   cvals(4) = half; svals(4) = -sqrt(three)*half
   cvals(5) = -sqrt(three)*half; svals(5) = half
   cvals(6) = -sqrt(three)*half; svals(6) = -half
   cvals(7) = one; svals(7) = zero
   cvals(8) = -one; svals(8) = zero
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
     call random_number(work1(1:2*ns))
     do i = 1, 2*ns
       work1(i) = one - two*work1(i)
     end do
     call random_number(work2(1:2*ns))
     do i = 1, 2*ns
       work2(i) = one - two*work2(i)
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
         xs(ix) = cmplx( work1(2*i-1), work1(2*i), wp )
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
           ys(iy) = cmplx( work2(2*i-1), work2(2*i), wp )
           iwork2(iy) = 1
           iy = iy + incys
         end do
!
!        Test for each of the predetermined rotations
!
         do ir = 1, nr
           cs = cvals(ir)
           ss = svals(ir)
!
!          Call CSROT to apply the rotation described by c and s
!
           n = ns
           x(1:lenx) = xs(1:lenx)
           incx = incxs
           y(1:lenx) = ys(1:lenx)
           incy = incys
           c = cs
           s = ss
           call CSROT( n, x, incx, y, incy, c, s )
!
!          See what values changed inside the subroutine
!
           nerrs = 0
           call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
           call CARGCKX( nerrs, subnam, 'X', lenx, iwork1, x, xs, nout )
           call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
           call CARGCKX( nerrs, subnam, 'Y', lenx, iwork2, y, ys, nout )
           call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
           call SARGCKS( nerrs, subnam, 'C', c, cs, nout )
           call SARGCKS( nerrs, subnam, 'S', s, ss, nout )
!
!          Check that the expected values are in x and y.
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
           if( ns > 0 ) then
             if( incxs == 0 ) xexp = cmplx( work1(2*ns-1), work1(2*ns), wp )
             if( incys == 0 ) yexp = cmplx( work2(2*ns-1), work2(2*ns), wp )
           end if
           do i = 1, ns
             if( incxs /= 0 ) then
               xexp = cmplx( work1(2*i-1), work1(2*i), wp )
             end if
             if( incys /= 0 ) then
               yexp = cmplx( work2(2*i-1), work2(2*i), wp )
             end if
             xtmp = cs*xexp + ss*yexp
             xnrm = abs( xtmp )
             ytmp = -ss*xexp + cs*yexp
             ynrm = abs( ytmp )
             xexp = xtmp
             yexp = ytmp
!
!            Check value in x(i) against xexp.
!
             if( incxs /= 0 .or. i == ns ) then
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
             end if
             ix = ix + incxs
!
!            Check value in y(i) against yexp.
!
             if( incys /= 0 .or. i == ns ) then
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
             end if
           end do
           if( max(tmax,umax) >= thresh ) then
             write(nout,98) subnam, ns, incxs, incys, ir, max(tmax,umax)
             write(nout,97) imax, xmax, xmex, jmax, ymax, ymex
             nerrs = nerrs + 1
           end if
           if( nerrs > 0 ) nfail = nfail + 1
           ntests = ntests + 1
         end do
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
98 format( 1X, A6, ': N=', I6, ', INCX=', I4, ', INCY=', I4, ', IR=', I2, &
           ', test=', E15.8 )
97 format( '   occurs at X(', I6, ')=(', E15.8, ',', E15.8, &
           '), expected (', E15.8, ',', E15.8, ')', /, &
           '   and/or at Y(', I6, ')=(', E15.8, ',', E15.8, &
           '), expected (', E15.8, ',', E15.8, ')' )  
96 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
