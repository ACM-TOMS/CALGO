subroutine CTOMCOPY( n, nrep, x, ldx, y, ldy, xs, ys, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, one, ulp, czero
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(ldx,*), xs(*), y(ldy,*), ys(*)
!  ..
!
!  CTOMCOPY times CCOPY from the Level 1 BLAS, a subroutine to copy
!  a vector x to a vector y.  This is the out-of-cache test with
!  different vectors used for each trial.
!
!  E. Anderson
!  July 13, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, iy, &
      j, nerrs, nfail, ns, ntests
   real(wp) :: rogue, t1, t2, tottim, tmax, trat, ynrm
   complex(wp) :: yexp, ymax
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
!  ..
!
   nfail = 0
   ntests = 0
   rogue = -999.999_wp
   subnam = 'CCOPY '
   ns = n
   incxs = 1
   incys = 1
!
!  Set x(i) = cmplx( real(i), -real(i) )
!
   do i = 1, ns
     xs(i) = cmplx( real(i), -real(i), wp )
   end do
   ys(1:ns) = cmplx( rogue, -rogue, wp )
!
!  Do the test nrep times but only test the final one
!
   tottim = zero
   n = ns
   do j = 1, nrep
     x(1:n,j) = xs(1:n)
     y(1:n,j) = ys(1:n)
   end do
   t1 = SECOND()
   do j = 1, nrep
     n = ns
     incx = incxs
     incy = incys
     call CCOPY( n, x(1,j), incx, y(1,j), incy )
   end do
   t2 = SECOND()
   tottim = tottim + t2 - t1
   write(nout,99) 'CCOPY ', ns, tottim / nrep
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call CARGCKV( nerrs, subnam, 'X', ns, x(1,nrep), xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!  Check that the expected values are in y
!
   iy = 1
   if( incys < 0 ) iy = 1 - (ns-1)*incys
   imax = 1
   ymax = zero
   tmax = -one
   do i = 1, ns
     yexp = xs(i)
     ynrm = abs( yexp )
     if( ynrm == zero ) then
       trat = abs( y(iy,nrep) ) / ulp
     else
       trat = ( abs( y(iy,nrep) - yexp ) / ynrm ) / ulp
     end if
     if( trat > tmax ) then
       imax = i
       ymax = y(iy,nrep)
       tmax = trat
     end if
     iy = iy + incys
   end do
   if( tmax >= thresh ) then
     write(nout,98) subnam, ns, incxs, incys, tmax
     yexp = xs(imax)
     write(nout,97) imax, ymax, yexp
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
98 format( A6, ': N=', I8, ', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at Y(', I8, ')=', E15.8, ', expected=', E15.8 )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
