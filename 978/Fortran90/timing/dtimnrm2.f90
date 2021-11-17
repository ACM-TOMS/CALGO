subroutine DTIMNRM2( n, nrep, x, xs, thresh, rtim, nout )
   use LA_CONSTANTS, only: wp, zero, half, one, two, ulp, smlnum, bignum, &
                           sbig, ssml, tbig, tsml
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: rtim(*), x(*), xs(*)
!  ..
!
!  DTIMNRM2 times the 2-norm routine DNRM2 from the Level 1 BLAS.
!  This is the in-cache test with the same vector used for each trial.
!
!  E. Anderson
!  July 8, 2016
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
   real(wp) :: DNRM2, DSECND
   external :: DNRM2, DSECND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'DNRM2 '
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
       call random_number(xs(1:n))
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
         xs(i) = one - two*xs(i)
         ssq1 = ssq1 + xs(i)*xs(i)
       end do
     else if( itype == 2 ) then
!
!      Scale x by a large number
!
       do i = 1, n
         if( xs(i) < half ) xs(i) = xs(i) - one
         ssq1 = ssq1 + xs(i)*xs(i)
         xs(i) = bignum*xs(i)
       end do
       v1 = bignum
     else if( itype == 3 ) then
!
!      Scale x by a small number
!
       do i = 1, n
         if( xs(i) < half ) xs(i) = xs(i) - one
         ssq1 = ssq1 + xs(i)*xs(i)
         xs(i) = smlnum*xs(i)
       end do
       v1 = smlnum
     else if( itype == 4 ) then
!
!      Set x to zero
!
       do i = 1, n
         xs(i) = zero
       end do
     else if( itype == 5 ) then
!
!      Set x to alternating safe and large values
!
       do i = 1, n-1, 2
         xs(i) = one - two*xs(i)
         ssq1 = ssq1 + xs(i)*xs(i)
         if( xs(i+1) < half ) xs(i+1) = xs(i+1) - one
         ssq2 = ssq2 + xs(i+1)*xs(i+1)
         xs(i+1) = bignum*xs(i+1)
       end do
       if( mod( n, 2 ) == 1 ) then
         xs(n) = one - two*xs(n)
         ssq1 = ssq1 + xs(n)*xs(n)
       end if
       v2 = bignum
     else if( itype == 6 ) then
!
!      Set x to alternating safe and small values
!
       do i = 1, n-1, 2
         xs(i) = one - two*xs(i)
         ssq1 = ssq1 + xs(i)*xs(i)
         if( xs(i+1) < half ) xs(i+1) = xs(i+1) - one
         ssq2 = ssq2 + xs(i+1)*xs(i+1)
         xs(i+1) = smlnum*xs(i+1)
       end do
       if( mod( n, 2 ) == 1 ) then
         xs(n) = one - two*xs(n)
         ssq1 = ssq1 + xs(n)*xs(n)
       end if
       v2 = smlnum
     else if( itype == 7 ) then
!
!      Set x to alternating safe and zero values
!
       do i = 1, n-1, 2
         xs(i) = one - two*xs(i)
         ssq1 = ssq1 + xs(i)*xs(i)
         xs(i+1) = zero
       end do
       if( mod( n, 2 ) == 1 ) then
         xs(n) = one - two*xs(n)
         ssq1 = ssq1 + xs(n)*xs(n)
       end if
     else if( itype == 8 ) then
!
!      Set x to alternating safe, large, small, and zero values
!
       do i = 1, n-3, 4
         xs(i) = one - two*xs(i)
         ssq1 = ssq1 + xs(i)*xs(i)
         if( xs(i+1) < half ) xs(i+1) = xs(i+1) - one
         ssq2 = ssq2 + xs(i+1)*xs(i+1)
         xs(i+1) = bignum*xs(i+1)
         if( xs(i+2) < half ) xs(i+2) = xs(i+2) - one
         xs(i+2) = smlnum*xs(i+2)
         xs(i+3) = zero
       end do
!
!      Clean up if n is not a multiple of 4
!
       n4 = n - mod(n,4)
       i = n4+1
       if( i <= n ) then
         xs(i) = one - two*xs(i)
         ssq1 = ssq1 + xs(i)*xs(i)
       end if
       i = n4+2
       if( i <= n ) then
         if( xs(i) < half ) xs(i) = xs(i) - one
         ssq2 = ssq2 + xs(i)*xs(i)
         xs(i) = bignum*xs(i)
       end if
       i = n4+3
       if( i <= n ) then
         if( xs(i) < half ) xs(i) = xs(i) - one
         xs(i) = smlnum*xs(i)
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
     x(1:n) = xs(1:n)
     incx = incxs
     t1 = DSECND()
     do j = 1, nrep
       sbla(j) = DNRM2( n, x, incx )
     end do
     t2 = DSECND()
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
     call DARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
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
