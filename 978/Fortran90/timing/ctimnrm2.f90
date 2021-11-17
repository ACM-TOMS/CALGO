subroutine CTIMNRM2( n, nrep, x, xs, work, thresh, rtim, nout )
   use LA_CONSTANTS32, only: wp, zero, half, one, two, ulp, smlnum, bignum, &
                           sbig, ssml, tbig, tsml, czero
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: rtim(*), work(*)
   complex(wp) :: x(*), xs(*)
!  ..
!
!  CTIMNRM2 times the 2-norm routine SCNRM2 from the Level 1 BLAS.
!  This is the in-cache test with the same vector used for each trial.
!
!  E. Anderson
!  July 13, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, itype, j, n4, nerrs, nfail, ns, ntests
   real(wp) :: ssq1, ssq2, t1, t2, tottim, trat, v1, v2, y1, y2, &
      ymin, ymax, znrm
!  ..
!  .. Local Arrays ..
   real(wp), allocatable :: sbla(:)
!  ..
!  .. External Functions ..
   real(wp) :: SCNRM2, SECOND
   external :: SCNRM2, SECOND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'SCNRM2'
   ns = n
   incxs = 1
   allocate(sbla(nrep))
!
!  Do for each of 8 vector types:
!   1. Values from the safe range that do not require scaling
!   2. Values from the range that would cause overflow
!   3. Values from the range that would cause underflow
!   4. All zero
!   5. Half safe, half large
!   6. Half safe, half small
!   7. Half safe, half zero
!   8. One fourth each of safe, large, small, and zero
!
   do itype = 1, 8 
!
!    Initialize x with random values in the range (0,1)
!
     if( itype /= 4 ) then
       call random_number(work(1:2*n))
     end if
     v1 = one
     v2 = one
     ssq1 = zero
     ssq2 = zero
     if( itype == 1 ) then
!
!      Set x to values in (-1,1)
!
       do i = 1, n
         y1 = one - two*work(2*i-1)
         y2 = one - two*work(2*i)
         xs(i) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
       end do
     else if( itype == 2 ) then
!
!      Scale x by a large number
!
       do i = 1, n
         y1 = work(2*i-1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i)
         if( y2 < half ) y2 = y2 - one
         ssq1 = ssq1 + y1*y1 + y2*y2
         xs(i) = cmplx( bignum*y1, bignum*y2, wp )
       end do
       v1 = bignum
     else if( itype == 3 ) then
!
!      Scale x by a small number
!
       do i = 1, n
         y1 = work(2*i-1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i)
         if( y2 < half ) y2 = y2 - one
         ssq1 = ssq1 + y1*y1 + y2*y2
         xs(i) = cmplx( smlnum*y1, smlnum*y2, wp )
       end do
       v1 = smlnum
     else if( itype == 4 ) then
!
!      Set x to zero
!
       do i = 1, n
         xs(i) = czero
       end do
     else if( itype == 5 ) then
!
!      Set x to alternating safe and large values
!
       do i = 1, n-1, 2
         y1 = one - two*work(2*i-1)
         y2 = one - two*work(2*i)
         xs(i) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
!
         y1 = work(2*i+1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i+2)
         if( y2 < half ) y2 = y2 - one
         xs(i+1) = cmplx( bignum*y1, bignum*y2, wp )
         ssq2 = ssq2 + y1*y1 + y2*y2
       end do
       if( mod( n, 2 ) == 1 ) then
         y1 = one - two*work(2*n-1)
         y2 = one - two*work(2*n)
         xs(n) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
       end if
       v2 = bignum
     else if( itype == 6 ) then
!
!      Set x to alternating safe and small values
!
       do i = 1, n-1, 2
         y1 = one - two*work(2*i-1)
         y2 = one - two*work(2*i)
         xs(i) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
!
         y1 = work(2*i+1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i+2)
         if( y2 < half ) y2 = y2 - one
         xs(i+1) = cmplx( smlnum*y1, smlnum*y2, wp )
         ssq2 = ssq2 + y1*y1 + y2*y2
       end do
       if( mod( n, 2 ) == 1 ) then
         y1 = one - two*work(2*n-1)
         y2 = one - two*work(2*n)
         xs(n) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
       end if
       v2 = smlnum
     else if( itype == 7 ) then
!
!      Set x to alternating safe and zero values
!
       do i = 1, n-1, 2
         y1 = one - two*work(2*i-1)
         y2 = one - two*work(2*i)
         xs(i) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
         xs(i+1) = czero
       end do
       if( mod( n, 2 ) == 1 ) then
         y1 = one - two*work(2*n-1)
         y2 = one - two*work(2*n)
         xs(n) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
       end if
     else if( itype == 8 ) then
!
!      Set x to alternating safe, large, small, and zero values
!
       do i = 1, n-3, 4
         y1 = one - two*work(2*i-1)
         y2 = one - two*work(2*i)
         xs(i) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
!
         y1 = work(2*i+1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i+2)
         if( y2 < half ) y2 = y2 - one
         xs(i+1) = cmplx( bignum*y1, bignum*y2, wp )
         ssq2 = ssq2 + y1*y1 + y2*y2
!
         y1 = work(2*i+3)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i+4)
         if( y2 < half ) y2 = y2 - one
         xs(i+2) = cmplx( smlnum*y1, smlnum*y2, wp )
!
         xs(i+3) = czero
       end do
!
!      Clean up if n is not a multiple of 4
!
       n4 = n - mod(n,4)
       i = n4+1
       if( i <= n ) then
         y1 = one - two*work(2*i-1)
         y2 = one - two*work(2*i)
         xs(i) = cmplx( y1, y2, wp )
         ssq1 = ssq1 + y1*y1 + y2*y2
       end if
       i = n4+2
       if( i <= n ) then
         y1 = work(2*i-1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i)
         if( y2 < half ) y2 = y2 - one
         ssq2 = ssq2 + y1*y1 + y2*y2
         xs(i) = cmplx( bignum*y1, bignum*y2, wp )
       end if
       i = n4+3
       if( i <= n ) then
         y1 = work(2*i-1)
         if( y1 < half ) y1 = y1 - one
         y2 = work(2*i)
         if( y2 < half ) y2 = y2 - one
         xs(i) = cmplx( smlnum*y1, smlnum*y2, wp )
       end if
       v2 = bignum
     end if
!
!    Compute the expected value of the sum of squares
!
     y1 = v1*sqrt(ssq1)
     y2 = v2*sqrt(ssq2)
     ymin = min( y1, y2 )
     ymax = max( y1, y2 )
     if( ymax == zero ) then
       znrm = zero
     else
       znrm = ymax * sqrt( one + ( ymin / ymax )**2 )
     end if
!
!    Do the test nrep times but only test the final one
!
     n = ns
     x(1:ns) = xs(1:ns)
     incx = incxs
     t1 = SECOND()
     do j = 1, nrep
       sbla(j) = SCNRM2( n, x, incx )
     end do
     t2 = SECOND()
     tottim = t2 - t1
     rtim(itype) = tottim / nrep
!    write(nout,99) subnam, ns, itype, rtim(itype)
!
!    Compare sbla and zrnm
!
     nerrs = 0
     if( znrm == zero ) then
       trat = sbla(nrep) / ulp
     else
       trat = ( abs( sbla(nrep) - znrm ) / znrm ) / ( real(max(1,n),wp)*ulp )
     end if
     if (trat >= thresh) then
       write(nout,98) subnam, ns, itype, trat
       nerrs = 1
     end if
!
!    See what values changed inside the subroutine
!
     call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
     call CARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
     call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
     if( nerrs > 0 ) nfail = nfail + 1
     ntests = ntests + 1
   end do
!
   if( nfail /= 0 ) then
     write(nout,97) subnam, nfail, ntests
   end if
   deallocate(sbla)
   return
99 format( A6, ': N=', I8, ', type=', I4, ', time=', E15.8 )
98 format( A6, ': N=', I8,', type=', I4, ', test=', E15.8 )
97 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
