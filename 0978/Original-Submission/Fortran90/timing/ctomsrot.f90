subroutine CTOMSROT( n, nrep, x, ldx, y, ldy, xs, ys, work, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, half, one, two, three, ulp
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp ) :: work(*)
   complex(wp) :: x(ldx,*), xs(*), y(ldy,*), ys(*)
!  ..
!
!  CTOMSROT times CSROT from the Level 1 BLAS, a subroutine to apply
!  a plane rotation to a pair of vectors x and y.  This is the
!  out-of-cache test with different vectors used for each trial.
!
!  E. Anderson
!  July 12, 2016
!
!  .. Parameters ..
   integer, parameter :: nr = 8
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, ir, ix, &
      iy, j, jmax, nerrs, nfail, ns, ntests
   real(wp) :: c, cs, rogue1, rogue2, s, ss, t1, t2, tottim, tmax, &
      trat, umax, urat, xnrm, ynrm
   complex(wp) :: xexp, xmax, yexp, ymax
!  ..
!  .. Local Arrays ..
   real(wp) :: cvals(nr), svals(nr)
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
   external :: SECOND
!  ..
!
   nfail = 0
   ntests = 0
   rogue1 = -666.666_wp
   rogue2 = -999.999_wp
   subnam = 'CSROT '
   ns = n
   incxs = 1
   incys = 1
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
!  Set x and y to random values in (-1,1).
!
   call random_number(work(1:2*n))
   do i = 1, n
     xs(i) = cmplx( one - two*work(2*i-1), one - two*work(2*i), wp )
   end do
   call random_number(work(1:2*n))
   do i = 1, n
     ys(i) = cmplx( one - two*work(2*i-1), one - two*work(2*i), wp )
   end do
!
!  Just test one of the rotations
!
   ir = 3
   cs = cvals(ir)
   ss = svals(ir)
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
   c = cs
   s = ss
   t1 = SECOND()
   do j = 1, nrep
     call CSROT( n, x(1,j), incx, y(1,j), incy, c, s )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   write(nout,99) 'CSROT ', ns, tottim / nrep
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
   call SARGCKS( nerrs, subnam, 'C', c, cs, nout )
   call SARGCKS( nerrs, subnam, 'S', s, ss, nout )
!
!  Check that the expected values are in x
!
   imax = 1
   xmax = zero
   tmax = -one
   ix = 1
   if( incxs < 0 ) ix = 1 - (ns-1)*incxs
   do i = 1, ns
     xexp = cs*xs(i) + ss*ys(i)
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
     yexp = -ss*xs(i) + cs*ys(i)
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
     write(nout,98) subnam, ns, incxs, incys, ir, max(tmax,umax)
     xexp = cs*xs(imax) + ss*ys(imax)
     yexp = -ss*xs(jmax) + cs*ys(jmax)
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
98 format( A6, ': N=', I8, ', INCX=', I4, ', INCY=', I4, ', IR=', I2, &
           ', test=', E15.8 )
97 format( '   occurs at X(', I8, ')=(', E15.8, ',', E15.8, &
           '), expected=(', E15.8, ',', E15.8, ')', /, &
           '   and/or at Y(', I8, ')=(', E15.8, ',', E15.8, &
           '), expected=(', E15.8, ',', E15.8, ')' )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
