subroutine ZB1TDOTC( nn, nvals, nix, ixvals, lenx, x, y, xs, ys, &
      work1, work2, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, one, two, ulp, czero
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
   real(wp) :: work1(*), work2(*)
   complex(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  Purpose
!  =======
!
!  ZB1TDOTC tests the dot product routine ZDOTC.
!
!  No special scaling is applied to the vectors x and y except for
!  Re(x(1)) and Re(x(m)) for some value of m near the middle of the range
!  1:n, which are scaled by 1.0, -Infinity, Infinity, or NaN.
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
!  .. Parameters ..
   integer, parameter :: nv = 4
!  ..
!  .. Local Scalars ..
   logical :: setexp
   character*6 :: subnam
   integer :: i, in, incx, incxs, incy, incys, iv, iw, ix, ixmax, iy, &
      j, k, m, n, nerrs, nfail, nmax, ns, ntests
   real(wp) :: pinf, pnan, snrm, trat, v0, v1
   complex(wp) :: rogue1, rogue2, scom, sexp, sexp1, sexp2, sexp3, &
      sexp4, xn, xtmp, yn, ytmp
!  ..
!  .. Local Arrays ..
   real(wp) :: values(nv)
!  ..
!  .. External Functions ..
   complex(wp) :: ZDOTC
   external :: ZDOTC
!  ..
!
   values(1) = one
   values(2) = LA_XVALS( 1, v0 )
   values(3) = LA_XVALS( 2, v0 )
   values(4) = LA_XVALS( 3, v0 )
   nfail = 0
   ntests = 0
   rogue1 = cmplx( 666.666_wp, -666.666_wp, wp )
   rogue2 = cmplx( 999.888_wp, -999.999_wp, wp )
   subnam = 'ZDOTC '
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
     call random_number(work1(1:2*ns))
     do i = 1, 2*ns
       work1(i) = one - two*work1(i)
     end do
     call random_number(work2(1:2*ns))
     do i = 1, 2*ns
       work2(i) = one - two*work2(i)
     end do
     m = int(ns/2) + 1
!
!    Make signs the same in positions 1 and m
!
     if( ns > 0 ) then
       work1(2*m-1) = sign( work1(2*m-1), work1(1) )
       work1(2*m) = sign( work1(2*m), work1(2) )
       work2(2*m-1) = sign( work2(2*m-1), work2(1) )
       work2(2*m) = sign( work2(2*m), work2(2) )
     end if
!
!    Choose a scaling for Re(x(1))
!
     do iv = 1, nv
       v0 = values(iv)
!
!      Choose a scaling for Re(x(m)) for m near the
!      middle of the range 1:n.
!
       do iw = 1, nv
         v1 = values(iw)
!
!        Calculate the expected value of the dot product
!
         setexp = .true.
         pinf = LA_XVALS( 2, pinf )
         pnan = LA_XVALS( 3, pnan )
         if( ns <= 0 ) then
           sexp = zero
         else if( ns == 1 ) then
           xtmp = cmplx( v0*work1(1), work1(2), wp )
           ytmp = cmplx( work2(1), work2(2), wp )
           sexp = conjg(xtmp)*ytmp
         else if( LA_ISNAN( v0 ) ) then
           sexp = cmplx( pnan, pnan )
         else
!
!          Calculate the dot product four ways:
!          sexp1: incx /= 0, incy /= 0
!          sexp2: incx == 0, incy /= 0
!          sexp3: incx /= 0, incy == 0
!          sexp4: incx == 0, incy == 0
!
           sexp1 = zero
           sexp2 = zero
           sexp3 = zero
           sexp4 = zero
           xn = cmplx( v0*work1(2*ns-1), work1(2*ns), wp )
           yn = cmplx( work2(2*ns-1), work2(2*ns), wp )
           do i = 1, ns
             ytmp = cmplx( work2(2*i-1), work2(2*i), wp )
             if( i == 1 ) then
               xtmp = cmplx( v0*work1(2*i-1), work1(2*i), wp )
               sexp1 = sexp1 + conjg( xtmp )*ytmp
               sexp2 = sexp2 + conjg( xn )*ytmp
               sexp3 = sexp3 + conjg( xtmp )*yn
               sexp4 = sexp4 + conjg( xn )*yn
             else if( i == m ) then
               xtmp = cmplx( v1*work1(2*i-1), work1(2*i), wp )
               sexp1 = sexp1 + conjg( xtmp )*ytmp
               sexp2 = sexp2 + conjg( xn )*ytmp
               sexp3 = sexp3 + conjg( xtmp )*yn
               sexp4 = sexp4 + conjg( xn )*yn
             else
               xtmp = cmplx( work1(2*i-1), work1(2*i), wp )
               sexp1 = sexp1 + conjg( xtmp )*ytmp
               sexp2 = sexp2 + conjg( xn )*ytmp
               sexp3 = sexp3 + conjg( xtmp )*yn
               sexp4 = sexp4 + conjg( xn )*yn
             end if
           end do
           setexp = .false.
         end if
         if( setexp ) then
           sexp1 = sexp
           sexp2 = sexp
           sexp3 = sexp
           sexp4 = sexp
         end if
!
!        Test for each value of incx
!
         do j = 1, nix
           incxs = ixvals(j)
           xs(1:lenx) = rogue1
           ix = 1
           if( incxs < 0 ) ix = 1 - (ns-1)*incxs
           do i = 1, ns
             xs(ix) = cmplx( work1(2*i-1), work1(2*i), wp )
             ix = ix + incxs
           end do
!
!          Scale Re(x(1)) and Re(x(m))
!
           if( ns > 0 ) then
             ix = 1
             if( incxs < 0 ) ix = 1 - (ns-1)*incxs
             xs(ix) = cmplx( v0*real(xs(ix)), aimag(xs(ix)), wp )
           end if
           if( ns > 1 .and. incxs /= 0 ) then
             if( incxs > 0 ) then
               ix = 1 + (m-1)*incxs
             else
               ix = 1 - (ns-m)*incxs
             end if
             xs(ix) = cmplx( v1*real(xs(ix)), aimag(xs(ix)), wp )
           end if
!
!          Test for each value of incy
!
           do k = 1, nix
             incys = ixvals(k)
             ys(1:lenx) = rogue2
             iy = 1
             if( incys < 0 ) iy = 1 - (ns-1)*incys
             do i = 1, ns
               ys(iy) = cmplx( work2(2*i-1), work2(2*i), wp )
               iy = iy + incys
             end do
!
!            Call ZDOTC to compute the dot product of x and y
!
             n = ns
             x(1:lenx) = xs(1:lenx)
             incx = incxs
             y(1:lenx) = ys(1:lenx)
             incy = incys
             scom = ZDOTC( n, x, incx, y, incy )
!
!            See what values changed inside the subroutine
!
             nerrs = 0
             call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
             call ZARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
             call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
             call ZARGCKV( nerrs, subnam, 'Y', lenx, y, ys, nout )
             call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!            Compare scom and sexp
!
             if( incxs /= 0 ) then
               if( incys /= 0 ) then
                 sexp = sexp1
               else
                 sexp = sexp3
               end if
             else
               if( incys /= 0 ) then
                 sexp = sexp2
               else
                 sexp = sexp4
               end if
             end if
             snrm = abs( sexp )
             if( snrm == zero ) then
               trat = abs( scom ) / ulp
             else
               trat = ( abs( scom - sexp ) / snrm ) / ( real(max(1,ns),wp)*ulp )
             end if
             if (trat >= thresh) then
               write(nout,98) subnam, ns, incxs, incys, iv, iw, trat
               nerrs = nerrs + 1
             end if
             if( nerrs > 0 ) nfail = nfail + 1
             ntests = ntests + 1
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
98 format( 1X, A6, ': N=', I6,', INCX=', I4, ', INCY=', I4, ', IV=', I1, &
           ', IW=', I1, ', test=', E15.8 )
97 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
