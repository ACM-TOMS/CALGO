subroutine DTOMAXPY( n, nrep, nal, alvals, x, ldx, y, ldy, xs, ys, &
                     thresh, nout )
   use LA_CONSTANTS, only: wp, zero, ulp, one
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, n, nrep, nal, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: alvals(*), x(ldx,*), xs(*), y(ldy,*), ys(*)
!  ..
!
!  DTOMAXPY times DAXPY from the Level 1 BLAS, a subroutine to compute
!  the vector sum y <- alpha*x + y.  This is the out-of-cache test with
!  different vectors used for each trial.
!
!  E. Anderson
!  July 8, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, iy, &
      j, l, nerrs, nfail, ns, ntests
   real(wp) :: alpha, alphas, rogue, tmax, trat, yexp, ymax, ynrm, &
      t1, t2, tottim
!  ..
!  .. External Functions ..
   real(wp) :: DSECND
!  ..
!
   nfail = 0
   ntests = 0
   rogue = -ONE
   subnam = 'DAXPY '
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
!  Test for each value of alpha
!
   do l = 1, nal
     alphas = alvals(l)
!
!    Do the test nrep times but only test the final one
!
     n = ns
     alpha = alphas
     do j = 1, nrep
       x(1:n,j) = xs(1:n)
       y(1:n,j) = ys(1:n)
     end do
     incx = incxs
     incy = incys
     t1 = DSECND()
     do j = 1, nrep
       call DAXPY( n, alpha, x(1,j), incx, y(1,j), incy )
     end do
     t2 = DSECND()
     tottim = t2 - t1
     write(nout,99) 'DAXPY ', ns, alphas, tottim / nrep
!
!    See what values changed inside the subroutine
!
     nerrs = 0
     call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
     call DARGCKS( nerrs, subnam, 'ALPHA', alpha, alphas, nout )
     call DARGCKV( nerrs, subnam, 'X', n, x(1,nrep), xs, nout )
     call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
     call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!    Check that the expected values are in y.
!    In case of failure, only report the largest one.
!
     imax = 1
     ymax = zero
     tmax = -one
     iy = 1
     if( incys < 0 ) iy = 1 - (ns-1)*incys
     do i = 1, ns
       yexp = alphas*xs(i)+ys(i)
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
       write(nout,98) subnam, ns, alphas, incxs, incys, tmax
       yexp = alphas*xs(imax)+ys(imax)
       write(nout,97) imax, ymax, yexp
       nerrs = nerrs + 1
     end if
     if( nerrs > 0 ) nfail = nfail + 1
     ntests = ntests + 1
   end do
   if( nfail /= 0 ) then
     write(*,95) subnam, nfail, ntests
   end if
   return
!
99 format( A6, ': N=', I8, ', ALPHA=', E15.8, ', time=', E15.8 )
98 format( A6, ': N=', I8, ', ALPHA=', E15.8, ', INCX=', I4, &
           ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at Y(', I8, ')=', E15.8, ', expected=', E15.8 )  
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
