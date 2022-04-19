subroutine CTOMDOTU( n, nrep, x, ldx, y, ldy, xs, ys, work, thresh, nout )
   use LA_CONSTANTS32, only: wp, zero, one, two, ulp, czero
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: work(*)
   complex(wp) :: x(ldx,*), xs(*), y(ldy,*), ys(*)
!  ..
!
!  CTOMDOTU times CDOTU from the Level 1 BLAS, a function to compute the
!  dot product of a vector x and a vector y.  This is the out-of-cache
!  test with different vectors used for each trial.
!
!  E. Anderson
!  July 13, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, incy, incys, j, nerrs, nfail, ns, ntests
   real(wp) :: rogue1, rogue2, snrm, trat, t1, t2, tottim
   complex(wp) :: cexp
!  ..
!  .. Local Arrays ..
   complex(wp), allocatable :: cbla(:)
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
   complex(wp) :: CDOTU
   external :: CDOTU, SECOND
!  ..
!
   nfail = 0
   ntests = 0
   rogue1 = -666.666_wp
   rogue2 = -999.999_wp
   subnam = 'CDOTU '
   ns = n
   incxs = 1
   incys = 1
   allocate(cbla(nrep))
!
!  Set xs and ys to random values from (-1,1).
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
!  Calculate the expected value of the dot product
!
   cexp = czero
   do i = 1, n
     cexp = cexp + xs(i)*ys(i)
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
     cbla(j) = CDOTU( n, x(1,j), incx, y(1,j), incy )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   write(nout,99) 'CDOTU ', ns, tottim / nrep
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call CARGCKV( nerrs, subnam, 'X', ns, x(1,nrep), xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call CARGCKV( nerrs, subnam, 'Y', ns, y(1,nrep), ys, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!  Compare cbla and cexp
!
   snrm = abs( cexp )
   if( snrm == zero ) then
     trat = abs( cbla(nrep) ) / ulp
   else
     trat = ( abs( cbla(nrep) - cexp ) / snrm ) / ( real(max(1,n),wp)*ulp )
   end if
   if (trat >= thresh) then
     write(nout,98) subnam, ns, incxs, incys, trat
     nerrs = nerrs + 1
   end if
   if( nerrs > 0 ) nfail = nfail + 1
   ntests = ntests + 1
!
   if( nfail /= 0 ) then
     write(nout,96) subnam, nfail, ntests
   end if
   deallocate(cbla)
   return
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( A6, ': N=', I8,', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
