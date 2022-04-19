subroutine CTOMAXPY( n, nrep, nal, alvals, x, ldx, y, ldy, xs, ys, &
                     work, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, ulp, one, two
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, n, nrep, nal, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: work(*)
   complex(wp) :: alvals(*), x(ldx,*), xs(*), y(ldy,*), ys(*)
!  ..
!
!  CTOMAXPY times CAXPY from the Level 1 BLAS, a subroutine to compute
!  the vector sum y <- alpha*x + y.  This is the out-of-cache test with
!  different vectors used for each trial.
!
!  E. Anderson
!  July 13, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, iy, &
      j, l, nerrs, nfail, ns, ntests
   real(wp) :: tmax, trat, ynrm, t1, t2, tottim
   complex(wp) :: alpha, alphas, yexp, ymax
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'CAXPY '
   ns = n
   incxs = 1
   incys = 1
   call random_number(work(1:2*ns))
   do i = 1, ns
     xs(i) = cmplx( one - two*work(2*i-1), one - two*work(2*i), wp )
   end do
   call random_number(work(1:2*ns))
   do i = 1, ns
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
     tottim = zero
     do j = 1, nrep
       x(1:ns,j) = xs(1:ns)
       y(1:ns,j) = ys(1:ns)
     end do
     t1 = SECOND()
     do j = 1, nrep
       n = ns
       alpha = alphas
       incx = incxs
       incy = incys
       call CAXPY( n, alpha, x(1,j), incx, y(1,j), incy )
     end do
     t2 = SECOND()
     tottim = tottim + t2 - t1
     write(nout,99) 'CAXPY ', ns, alphas, tottim / nrep
!
!    See what values changed inside the subroutine
!
     nerrs = 0
     call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
     call CARGCKS( nerrs, subnam, 'ALPHA', alpha, alphas, nout )
     call CARGCKV( nerrs, subnam, 'X', n, x(1,nrep), xs, nout )
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
99 format( A6, ': N=', I8, ', ALPHA=(', E15.8, ',', E15.8, &
           '), time=', E15.8 )
98 format( A6, ': N=', I8, ', ALPHA=(', E15.8, ',', E15.8, &
           '), INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at Y(', I8, ')=(', E15.8, ',', E15.8, &
           '), expected=(', E15.8, ',', E15.8, ')' )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
