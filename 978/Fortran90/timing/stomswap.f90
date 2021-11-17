subroutine STOMSWAP( n, nrep, x, ldx, y, ldy, xs, ys, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, one, ulp
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: x(ldx,*), xs(*), y(ldy,*), ys(*)
!  ..
!
!  STOMSWAP times SSWAP from the Level 1 BLAS, a subroutine to swap
!  two vectors.  This is the out-of-cache test with different vectors
!  used for each trial.
!
!  E. Anderson
!  July 12, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, ix, iy, &
      j, jmax, nerrs, nfail, ns, ntests
   real(wp) :: t1, t2, tottim, tmax, trat, umax, urat, &
      xexp, xmax, xnrm, yexp, ymax, ynrm
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
   external :: SECOND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'SSWAP '
   ns = n
   incxs = 1
   incys = 1
!
!  Set x(i) = real(i) and y(i) = real(i)/1000
!
   do i = 1, n
     xs(i) = real(i)
     ys(i) = real(i) / 1000._wp
   end do
!
!  Do the test nrep times but only test the final one
!
   n = ns
   do j = 1, nrep
     x(1:n,j) = xs(1:n)
     y(1:n,j) = ys(1:n)
   end do
   incx = incxs
   incy = incys
   t1 = SECOND()
   do j = 1, nrep
     call SSWAP( n, x(1,j), incx, y(1,j), incy )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   write(nout,99) 'SSWAP ', ns, tottim / nrep
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!  Check that the expected values are in x
!
   imax = 1
   xmax = zero
   tmax = -one
   ix = 1
   if( incxs < 0 ) ix = 1 - (ns-1)*incxs
   do i = 1, ns
     xexp = ys(i)
     xnrm = abs( xexp )
     if( xnrm == zero ) then
       trat = abs( x(ix,nrep) ) / ulp
     else
       trat = ( abs( x(ix,nrep) - xexp ) / xnrm ) / ulp
     end if
     if( trat > tmax ) then
       imax = i
       xmax = x(ix,nrep)
       tmax = trat
     end if
     ix = ix + incxs
   end do
!
!  Check that the expected values are in y
!
   jmax = 1
   ymax = zero
   umax = -one
   iy = 1
   if( incys < 0 ) iy = 1 - (ns-1)*incys
   do i = 1, ns
     yexp = xs(i)
     ynrm = abs( yexp )
     if( ynrm == zero ) then
       urat = abs( y(iy,nrep) ) / ulp
     else
       urat = ( abs( y(iy,nrep) - yexp ) / ynrm ) / ulp
     end if
     if( urat > umax ) then
       jmax = i
       ymax = y(iy,nrep)
       umax = urat
     end if
     iy = iy + incys
   end do
   if( max(tmax,umax) >= thresh ) then
     write(nout,98) subnam, ns, incxs, incys, max(tmax,umax)
     xexp = ys(imax)
     yexp = xs(jmax)
     write(nout,97) imax, xmax, xexp, jmax, ymax, yexp
     nerrs = nerrs + 1
   end if
   if( nerrs > 0 ) nfail = nfail + 1
   ntests = ntests + 1
!
   if( nfail /= 0 ) then
     write(nout,95) subnam, nfail, ntests
   end if
   return
!
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( A6, ': N=', I6, ', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at X(', I6, ')=', E15.8, ', expected ', E15.8, /, &
           '   and/or at Y(', I6, ')=', E15.8, ', expected ', E15.8 )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
