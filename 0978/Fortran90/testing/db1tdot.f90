subroutine DB1TDOT( nn, nvals, nix, ixvals, lenx, x, y, xs, ys, &
      work1, work2, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, one, two, ulp
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
   real(wp) :: x(*), xs(*), y(*), ys(*), work1(*), work2(*)
!  ..
!
!  Purpose
!  =======
!
!  DB1TDOT tests the dot product routine DDOT.
!
!  No special scaling is applied to the vectors x and y except for
!  x(1) and x(m) for some value of m near the middle of the range 1:n,
!  which are scaled by 1.0, -Infinity, Infinity, or NaN.
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
!  .. Parameters ..
   integer, parameter :: nv = 4
!  ..
!  .. Local Scalars ..
   logical :: setexp
   character*6 :: subnam
   integer :: i, in, incx, incxs, incy, incys, iv, iw, ix, ixmax, &
      iy, j, k, m, n, nerrs, nfail, nmax, ns, ntests
   real(wp) :: pinf, rogue1, rogue2, scom, sexp, sexp1, sexp2, sexp3, &
      sexp4, snrm, trat, v0, v1, xn, yn
!  ..
!  .. Local Arrays ..
   real(wp) :: values(nv)
!  ..
!  .. External Functions ..
   real(wp) :: DDOT
   external :: DDOT
!  ..
!
   values(1) = one
   values(2) = LA_XVALS( 1, v0 )
   values(3) = LA_XVALS( 2, v0 )
   values(4) = LA_XVALS( 3, v0 )
   nfail = 0
   ntests = 0
   rogue1 = -666.666_wp
   rogue2 = -999.999_wp
   subnam = 'DDOT  '
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
     call random_number(work1(1:ns))
     do i = 1, ns
       work1(i) = one - two*work1(i)
     end do
     call random_number(work2(1:ns))
     do i = 1, ns
       work2(i) = one - two*work2(i)
     end do
     m = int(ns/2) + 1
!
!    Make signs the same in positions 1 and m of both vectors
!
     if( ns > 0 ) then
       work1(m) = sign( work1(m), work1(1) )
       work2(m) = sign( work2(m), work2(1) )
     end if
!
!    Choose a scaling for x(1)
!
     do iv = 1, nv
       v0 = values(iv)
!
!      Choose a scaling for x(m) for m near the
!      middle of the range 1:n.
!
       do iw = 1, nv
         v1 = values(iw)
!
!        Calculate the expected value of the dot product.
!
         setexp = .true.
         pinf = LA_XVALS( 2, pinf )
         if( ns <= 0 ) then
           sexp = zero
         else if( ns == 1 ) then
           sexp = v0*work1(1)*work2(1)
         else if( LA_ISNAN( v0 ) ) then
           sexp = LA_XVALS( 3, v0 )
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
           xn = v0*work1(ns)
           yn = work2(ns)
           do i = 1, ns
             if( i == 1 ) then
               sexp1 = sexp1 + v0*work1(i)*work2(i)
               sexp2 = sexp2 + xn*work2(i)
               sexp3 = sexp3 + v0*work1(i)*yn
               sexp4 = sexp4 + xn*yn
             else if( i == m ) then
               sexp1 = sexp1 + v1*work1(i)*work2(i)
               sexp2 = sexp2 + xn*work2(i)
               sexp3 = sexp3 + v1*work1(i)*yn
               sexp4 = sexp4 + xn*yn
             else
               sexp1 = sexp1 + work1(i)*work2(i)
               sexp2 = sexp2 + xn*work2(i)
               sexp3 = sexp3 + work1(i)*yn
               sexp4 = sexp4 + xn*yn
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
             xs(ix) = work1(i)
             ix = ix + incxs
           end do
!
!          Scale x(1) and x(m)
!
           if( ns > 0 ) then
             ix = 1
             if( incxs < 0 ) ix = 1 - (ns-1)*incxs
             xs(ix) = v0*xs(ix)
           end if
           if( ns > 1 .and. incxs /= 0 ) then
             if( incxs > 0 ) then
               ix = 1 + (m-1)*incxs
             else
               ix = 1 - (ns-m)*incxs
             end if
             xs(ix) = v1*xs(ix)
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
               ys(iy) = work2(i)
               iy = iy + incys
             end do
!
!            Call DDOT to compute the dot product of x and y
!
             n = ns
             x(1:lenx) = xs(1:lenx)
             incx = incxs
             y(1:lenx) = ys(1:lenx)
             incy = incys
             scom = DDOT( n, x, incx, y, incy )
!
!            See what values changed inside the subroutine
!
             nerrs = 0
             call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
             call DARGCKV( nerrs, subnam, 'X', lenx, x, xs, nout )
             call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
             call DARGCKV( nerrs, subnam, 'Y', lenx, y, ys, nout )
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
             if( LA_ISNAN( scom ) .or. LA_ISNAN( sexp ) ) then
               if( LA_ISNAN( scom ) .neqv. LA_ISNAN( sexp ) ) then
                  trat = one / ulp
               else
                  trat = zero
               end if
             else if( scom == sexp ) then
               trat = zero
             else if( snrm == zero ) then
               trat = abs( scom ) / ulp
             else
               trat = ( abs( scom - sexp ) / snrm ) / ( real(max(1,ns),wp)*ulp )
             end if
             if( LA_ISNAN( trat ) .or. trat >= thresh ) then
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
