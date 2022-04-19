subroutine DB1TAXPY( nn, nvals, nix, ixvals, nal, alvals, lenx, &
                     x, y, xs, ys, work1, work2, iwork, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, ulp, one, bignum
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
   real(wp) :: alvals(*), x(*), xs(*), y(*), ys(*), work1(*), work2(*)
!  ..
!
!  Purpose
!  =======
!
!  DB1TAXPY tests the vector sum routine DAXPY.
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
!  Y       (workspace) REAL array, dimension (LENX)
!
!  XS      (workspace) REAL array, dimension (LENX)
!
!  YS      (workspace) REAL array, dimension (LENX)
!
!  WORK1   (workspace) REAL array, dimension (NMAX)
!
!  WORK2   (workspace) REAL array, dimension (NMAX)
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
   integer :: i, imax, in, incx, incxs, incy, incys, ix, ixmax, iy, &
      j, k, l, n, nerrs, nfail, nmax, ns, ntests
   real(wp) :: alpha, alphas, rogue, tmax, trat, xexp, yexp, ymax, &
      ymex, ynrm
!  ..
!
   nfail = 0
   ntests = 0
   rogue = bignum
   subnam = 'DAXPY '
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
!    Set x(i) = real(i) and y(i) = real(i)/1000
!
     do i = 1, ns
       work1(i) = real(i)
       work2(i) = real(i) / 1000._wp
     end do
!
!    Test for each value of incx
!
     do j = 1, nix
       incxs = ixvals(j)
       ix = 1
       if( incxs < 0 ) ix = 1 - (ns-1)*incxs
       xs(1:lenx) = rogue
       do i = 1, ns
         xs(ix) = work1(i)
         ix = ix + incxs
       end do
!
!      Test for each value of incy
!
       do k = 1, nix
         incys = ixvals(k)
         ys(1:lenx) = rogue
         iwork(1:lenx) = 0
         iy = 1
         if( incys < 0 ) iy = 1 - (ns-1)*incys
         do i = 1, ns
           ys(iy) = work2(i)
           iwork(iy) = 1
           iy = iy + incys
         end do
!
!        Test for each value of alpha
!
         do l = 1, nal
           alphas = alvals(l)
           n = ns
           alpha = alphas
           x(1:lenx) = xs(1:lenx)
           incx = incxs
           y(1:lenx) = ys(1:lenx)
           incy = incys
           call DAXPY( n, alpha, x, incx, y, incy )
!
!          See what values changed inside the subroutine
!
           nerrs = 0
           call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
           call DARGCKS( nerrs, subnam, 'ALPHA', alpha, alphas, nout )
           call DARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
           call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
           call DARGCKX( nerrs, subnam, 'Y', lenx, iwork, y, ys, nout )
           call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!          Check that the expected values are in y.
!          In case of failure, only report the largest one.
!
           imax = 1
           ymax = zero
           tmax = -one
           iy = 1
           if( incys < 0 ) iy = 1 - (ns-1)*incys
           xexp = work1(ns)
           yexp = work2(ns)
           do i = 1, ns
             if( incxs /= 0 ) then
               xexp = work1(i)
             end if
             if( incys == 0 ) then
               yexp = alphas*xexp+yexp
             else
               yexp = alphas*xexp+work2(i)
             end if
             if( incys /= 0 .or. i == ns ) then
               ynrm = abs( yexp )
               if( ynrm == zero ) then
                 trat = abs( y(iy) ) / ulp
               else
                 trat = ( abs( y(iy) - yexp ) / ynrm ) / ulp
               end if
               if( trat > tmax ) then
                 imax = i
                 ymax = y(iy)
                 ymex = yexp
                 tmax = trat
               end if
             end if
             iy = iy + incys
           end do
           if( tmax >= thresh ) then
             write(nout,98) subnam, ns, alphas, incxs, incys, tmax
             write(nout,97) imax, ymax, ymex
             nerrs = nerrs + 1
           end if
           if( nerrs > 0 ) nfail = nfail + 1
           ntests = ntests + 1
         end do
       end do
     end do
   end do
   if( nfail.eq.0 ) then
      write(*,96) subnam, ntests
   else
     write(*,95) subnam, nfail, ntests
   end if
   return
!
99 format( ' Not enough space to test ', A6, ': NMAX = ',I6, &
           ' , IXMAX = ',I6,/,'   LENX = ',I6,', must be at least ',I6 )
98 format( 1X, A6, ': N=', I6, ', ALPHA=', E15.8, ', INCX=', I4, &
           ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at Y(', I6, ')=', E15.8, ', expected=', E15.8 )  
96 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
