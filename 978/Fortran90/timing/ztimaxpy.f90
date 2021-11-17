subroutine ZTIMAXPY( n, nrep, nal, alvals, x, y, xs, ys, &
                     work, thresh, rtim, nout )
   use LA_CONSTANTS, only: wp, zero, ulp, one, two, czero
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nal, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: rtim(*), work(*)
   complex(wp) :: alvals(*), x(*), xs(*), y(*), ys(*)
!  ..
!
!  ZTIMAXPY times ZAXPY from the Level 1 BLAS, a subroutine to compute
!  the vector sum y <- alpha*x + y.  This is the in-cache test with the
!  same vectors used for each trial.
!
!  E. Anderson
!  July 18, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, iy, &
      j, jmax, l, nerrs, nfail, ns, ntests
   real(wp) :: tmax, trat, ynrm, t1, t2, tottim
   complex(wp) :: alpha, alphas, yexp, ymax
!  ..
!  .. External Functions ..
   real(wp) :: DSECND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'ZAXPY '
   ns = n
   incxs = 1
   incys = 1
!
!  Fill the real and imginary parts of the vector x with random
!  values from (-1,1).
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
!  Test for each value of alpha
!
   do l = 1, nal
     alphas = alvals(l)
!
!    Do the test nrep times but only test the final one
!
     n = ns
     alpha = alphas
     x(1:n) = xs(1:n)
     incx = incxs
     y(1:n) = ys(1:n)
     incy = incys
     t1 = DSECND()
     if( nrep == 1 ) then
        jmax = 2
     else
        jmax = nrep - mod(nrep,2)
     end if
     incxs = 1
     incys = 1
     do j = 1, jmax-1, 2
       call ZAXPY( n, alpha, x, incx, y, incy )
       call ZAXPY( n, -alpha, x, incx, y, incy )
     end do
     t2 = DSECND()
     tottim = t2 - t1
     rtim(l) = tottim / jmax
!
     n = ns
     alpha = alphas
     x(1:n) = xs(1:n)
     incx = incxs
     y(1:n) = ys(1:n)
     incy = incys
     call ZAXPY( n, alpha, x, incx, y, incy )
!    write(nout,99) 'ZAXPY ', ns, alphas, rtim
!
!    See what values changed inside the subroutine
!
     nerrs = 0
     call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
     call ZARGCKS( nerrs, subnam, 'ALPHA', alpha, alphas, nout )
     call ZARGCKV( nerrs, subnam, 'X', n, x, xs, nout )
     call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
     call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!    Check that the expected values are in y.
!    In case of failure, only report the largest one.
!
     imax = 1
     ymax = czero
     tmax = -one
     iy = 1
     if( incys < 0 ) iy = 1 - (ns-1)*incys
     do i = 1, ns
       yexp = alphas*xs(i)+ys(i)
       ynrm = abs( yexp )
       if( ynrm == zero ) then
         trat = abs( y(iy) ) / ulp
       else
         trat = ( abs( y(iy) - yexp ) / ynrm ) / ulp
       end if
       if( trat > tmax ) then
         imax = i
         ymax = y(iy)
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
99 format( A6, ': N=', I8, ', ALPHA=(', E15.8, ',', E15.8, &
           '), time=', E15.8 )
98 format( A6, ': N=', I8, ', ALPHA=(', E15.8, ',', E15.8, &
           '), INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at Y(', I8, ')=(', E15.8, ',', E15.8, &
           '), expected=(', E15.8, ',', E15.8, ')' )  
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
