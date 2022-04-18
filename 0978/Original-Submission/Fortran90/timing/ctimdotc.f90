subroutine CTIMDOTC( n, nrep, x, y, xs, ys, work, thresh, rtim, nout )
   use LA_CONSTANTS32, only: wp, zero, one, two, ulp, czero
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim, thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: work(*)
   complex(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  CTIMDOTC times CDOTC from the Level 1 BLAS, a function to compute the
!  dot product of a vector x and a vector y.  This is the in-cache test
!  with the same vectors used for each trial.
!
!  E. Anderson
!  July 13, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, incy, incys, j, nerrs, nfail, ns, ntests
   real(wp) :: snrm, trat, t1, t2, tottim
   complex(wp) :: cexp
!  ..
!  .. Local Arrays ..
   complex(wp), allocatable :: cbla(:)
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
   complex(wp) :: CDOTC
   external :: CDOTC, SECOND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'CDOTC '
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
     cexp = cexp + conjg(xs(i))*ys(i)
   end do
!
!  Do the test nrep times but only test the final one
!
   n = ns
   x(1:n) = xs(1:n)
   incx = incxs
   y(1:n) = ys(1:n)
   incy = incys
   t1 = SECOND()
   do j = 1, nrep
     cbla(j) = CDOTC( n, x, incx, y, incy )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   rtim = tottim / nrep
!  write(nout,99) 'CDOTC ', ns, rtim
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call CARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call CARGCKV( nerrs, subnam, 'Y', ns, y, ys, nout )
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
     write(nout,*) ' trat = ',trat
     write(nout,*) ' thresh = ',thresh
     write(nout,*) ' cexp = ',cexp
     write(nout,*) ' ccom = ',cbla(nrep)
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
