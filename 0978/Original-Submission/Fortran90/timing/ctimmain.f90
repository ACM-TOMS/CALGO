program CTIMMAIN
   use LA_CONSTANTS32, only: wp
!
!  Time the Level 1 BLAS in cache and out of cache.
!
!  The purpose of this timing is to compare algorithms, not to
!  draw pretty performance curves.  Post processing would involve
!  comparing the in-cache and out-of-cache performance of two or
!  more algorithm variants.
!
!  A sample input file is as follows:
!  10000          Value of N for in-cache (IC) test
!  128            Offset between X and Y in the IC work array
!  100            Number of repetitions for IC test
!  100000         Value of N for out-of-cache (OC) test
!  128            Offset between X and Y in the OC work array
!  100            Number of repetitions for OC test
!  3              Number of values of alpha
!  1.0 -1.0 0.3   Values of alpha
!  all            Or list the algorithms to time one per line
!
!  The BLAS routines ICAMAX, SCASUM, CAXPY, CCOPY, CDOTC/CDOTU,
!  SCNRM2, CROT/CSROT, CSCAL/CSSCAL, and CSWAP are timed in this
!  procedure, as well as the LAPACK routine CLARGV (the vector version
!  of the BLAS CROTG or LAPACK CLARTG).  The algorithm types that
!  can be specified in the input file are amax, asum, axpy, copy,
!  dot, nrm2, rot, rotg, scal, and swap.
!
!  Output takes the form
!  In-cache test with N=10000, NOFF=100, NREP=100000
!  -------------------------------------------------
!  CAXPY : N=   10000, ALPHA= 0.10000000E+01, time= 0.99999993E-05
!  CAXPY : N=   10000, ALPHA=-0.10000000E+01, time= 0.99999993E-05
!  CAXPY : N=   10000, ALPHA= 0.30000000E+00, time= 0.99999993E-05
!  ...
!  Out-of-cache test with N=100000, NOFF=128, NREP=100
!  ---------------------------------------------------
!  ...
!
!  E. Anderson
!  July 26, 2016
!  =================================================================
!
!  Read the sizes of n, noff (offset), and nrep (no. of repetitions)
!  for the in-cache and out-of-cache tests.
!
!  .. Parameters ..
   integer, parameter :: nin = 5
   integer, parameter :: nout = 6
   integer, parameter :: ixamax = 1
   integer, parameter :: ixasum = 2
   integer, parameter :: ixaxpy = 3
   integer, parameter :: ixcopy = 4
   integer, parameter :: ixdotc = 5
   integer, parameter :: ixdotu = 6
   integer, parameter :: ixnrm2 = 7
   integer, parameter :: ixrot  = 8
   integer, parameter :: ixsrot = 9
   integer, parameter :: ixrotg = 10
   integer, parameter :: ixscal = 11
   integer, parameter :: ixsscal = 12
   integer, parameter :: ixswap = 13
   integer, parameter :: nsub   = 13
!  ..
!  .. Local Scalars ..
   logical :: fatal
   character*6 :: path
   character*80 :: aline
   integer :: i, ix, ixs, iy, iys, j, ldx, ldy, ldz, n, nal, ncal, &
      noff, nrep, noc, nocoff, nocrep
   real(wp) :: thresh
!  ..
!  .. Local Arrays ..
   logical :: dosub(nsub)
   character*6 :: subnam(nsub)
   real(wp), allocatable, dimension(:) :: alvals, rtim, work
   complex(wp), allocatable, dimension(:) :: calvals
   complex(wp), allocatable, dimension(:,:) :: wic, woc, voc
!  ..
!  .. External Functions ..
   logical :: lsamen
!  ..
!
   subnam(ixamax) = 'ICAMAX'
   subnam(ixasum) = 'SCASUM'
   subnam(ixaxpy) = 'CAXPY '
   subnam(ixcopy) = 'CCOPY '
   subnam(ixdotc) = 'CDOTC '
   subnam(ixdotu) = 'CDOTU '
   subnam(ixnrm2) = 'SCNRM2'
   subnam(ixrot)  = 'CROT  '
   subnam(ixsrot) = 'CSROT '
   subnam(ixrotg) = 'CLARGV'
   subnam(ixscal) = 'CSCAL '
   subnam(ixsscal) = 'CSSCAL'
   subnam(ixswap) = 'CSWAP '
!
!  Read the parameter values
!
   read( nin, * ) n
   read( nin, * ) noff
   read( nin, * ) nrep
   read( nin, * ) noc
   read( nin, * ) nocoff
   read( nin, * ) nocrep
!
!  Read the values of ALPHA
!
   read( nin, * ) nal
   if( nal < 1 ) then
      write( nout, 99 )'nal ', nal, 1
      fatal = .true.
      read( nin, * )
   else
      allocate( alvals( nal ) )
      read( nin, * )( alvals( i ), i = 1, nal )
   end if
!
!  Read the values of CALPHA
!
   read( nin, * ) ncal
   if( ncal < 1 ) then
      write( nout, 99 )'ncal', ncal, 1
      fatal = .true.
      read( nin, * )
   else
      allocate( calvals( ncal ) )
      read( nin, * )( calvals( i ), i = 1, ncal )
   end if
   allocate( rtim( max(nal,ncal,9) ) )
   read( nin, * ) thresh
!
!  Allocate storage for the in-cache and out-of-cache tests
!  We need storage for one or more of x, xsave, y, and ysave,
!  where the requirements are
!    amax: x
!    asum: x
!    axpy: x, y, ysave
!    copy: x, y, ysave
!    dot : x, y
!    nrm2: x
!    rot : x, y, xsave, ysave
!    rotg : x, y, xsave, ysave, z
!    scal: x, xsave
!    swap: x, y
!  Copies in xsave and ysave are only required for the in-cache
!  tests, the out-of-cache tests do not reuse vectors.
!  The test of rotg actually tests the LAPACK routine xlargv.
!
   allocate(work(1:2*max(n,noc,1)))
   allocate(wic(n+noff,5))
   allocate(woc(nocrep*(noc+nocoff),3))
   allocate(voc(noc+nocoff,2))
   ldx = noc+nocoff
   ldy = noc+nocoff
   ldz = noc+nocoff
   do i = 1, nsub
      dosub(i) = .false.
   end do
!
!  Read the test paths for later reference
!
10 continue
   read( nin, '(a)', end = 20 ) aline
   path = aline(1:6)
   if( lsamen( 4, path, 'all ' ) ) then 
      do i = 1, nsub
         dosub(i) = .true.
      end do
   else if( lsamen( 4, path, 'amax' ) .or. lsamen( 6, path, 'ICAMAX' ) ) then
      dosub(ixamax) = .true.
   else if( lsamen( 4, path, 'asum' ) .or. lsamen( 6, path, 'SCASUM' ) ) then
      dosub(ixasum) = .true.
   else if( lsamen( 4, path, 'axpy' ) .or. lsamen( 5, path, 'CAXPY' ) ) then
      dosub(ixaxpy) = .true.
   else if( lsamen( 4, path, 'copy' ) .or. lsamen( 5, path, 'CCOPY' ) ) then
      dosub(ixcopy) = .true.
   else if( lsamen( 4, path, 'dotc' ) .or. lsamen( 5, path, 'CDOTC' ) ) then
      dosub(ixdotc) = .true.
   else if( lsamen( 4, path, 'dotu' ) .or. lsamen( 5, path, 'CDOTU' ) ) then
      dosub(ixdotu) = .true.
   else if( lsamen( 4, path, 'nrm2' ) .or. lsamen( 6, path, 'SCNRM2' ) ) then
      dosub(ixnrm2) = .true.
   else if( lsamen( 4, path, 'rot ' ) .or. lsamen( 5, path, 'CROT ' ) ) then
      dosub(ixrot ) = .true.
   else if( lsamen( 4, path, 'srot' ) .or. lsamen( 5, path, 'CSROT' ) ) then
      dosub(ixsrot) = .true.
   else if( lsamen( 4, path, 'rotg' ) .or. lsamen( 5, path, 'CROTG' ) .or. &
        lsamen( 6, path, 'CLARTG' ) .or. lsamen( 6, path, 'CLARGV' ) ) then
      dosub(ixrotg) = .true.
   else if( lsamen( 4, path, 'scal' ) .or. lsamen( 5, path, 'CSCAL' ) ) then
      dosub(ixscal) = .true.
   else if( lsamen( 5, path, 'sscal' ) .or. lsamen( 6, path, 'CSSCAL' ) ) then
      dosub(ixsscal) = .true.
   else if( lsamen( 4, path, 'swap' ) .or. lsamen( 5, path, 'CSWAP' ) ) then
      dosub(ixswap) = .true.
   else
      write( nout, 98 ) path
   end if
   go to 10
20 continue
!
!  Do the in-cache tests
!
   ix = 1
   iy = 1 + n + noff
   ixs = iy + n + noff
   iys = ixs + n + noff
   write(nout,97) 'In-cache', n, noff, nrep
   do i = 1, nsub
      select case(i)
         case(ixamax)
            if( dosub(ixamax) ) then
            call CTIMAMAX( n, nrep, wic(1,1), wic(1,2), work, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixasum)
            if( dosub(ixasum) ) then
            call CTIMASUM( n, nrep, wic(1,1), wic(1,2), work, thresh, &
                           rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixaxpy)
            if( dosub(ixaxpy) ) then
            call CTIMAXPY( n, nrep, ncal, calvals, wic(1,1), wic(1,2), &
                           wic(1,3), wic(1,4), work, thresh, rtim, nout )
            do j = 1, ncal
              write(nout,95) subnam(i), n, calvals(j), rtim(j)
            end do
            end if
         case(ixcopy)
            if( dosub(ixcopy) ) then
            call CTIMCOPY( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixdotc)
            if( dosub(ixdotc) ) then
            call CTIMDOTC( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), work, thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixdotu)
            if( dosub(ixdotu) ) then
            call CTIMDOTU( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), work, thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixnrm2)
            if( dosub(ixnrm2) ) then
            call CTIMNRM2( n, nrep, wic(1,1), wic(1,2), work, thresh, &
                           rtim, nout )
            do j = 1, 8
               write(nout,93) subnam(i), n, j, rtim(j)
            end do
            end if
         case(ixrot )
            if( dosub(ixrot ) ) then
            call CTIMROT( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                          wic(1,4), work, thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixsrot)
            if( dosub(ixsrot) ) then
            call CTIMSROT( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), work, thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixrotg)
            if( dosub(ixrotg) ) then
            call CTIMRTG( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                          wic(1,4), wic(1,5), work, thresh, rtim, nout )
            do j = 1, 9
               write(nout,93) subnam(i), n, j, rtim(j)
            end do
            end if
         case(ixscal)
            if( dosub(ixscal) ) then
            call CTIMSCAL( n, nrep, ncal, calvals, wic(1,1), wic(1,2), &
                           work, thresh, rtim, nout )
            do j = 1, ncal
              write(nout,95) subnam(i), n, calvals(j), rtim(j)
            end do
            end if
         case(ixsscal)
            if( dosub(ixsscal) ) then
            call CTIMSSCAL( n, nrep, nal, alvals, wic(1,1), wic(1,2), &
                            work, thresh, rtim, nout )
            do j = 1, nal
              write(nout,94) subnam(i), n, alvals(j), rtim(j)
            end do
            end if
         case(ixswap)
            if( dosub(ixswap) ) then
            call CTIMSWAP( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
      end select
   end do
!
!  Do the out-of-cache tests
!
   write(nout,*)
   write(nout,97) 'Out-of-cache', noc, nocoff, nocrep
   do i = 1, nsub
      select case(i)
         case(ixamax)
            if( dosub(ixamax) ) then
            call CTOMAMAX( noc, nocrep, woc(1,1), ldx, voc(1,1), work, nout )
            end if
         case(ixasum)
            if( dosub(ixasum) ) then
            call CTOMASUM( noc, nocrep, woc(1,1), ldx, voc(1,1), work, &
                           thresh, nout )
            end if
         case(ixaxpy)
            if( dosub(ixaxpy) ) then
            call CTOMAXPY( noc, nocrep, ncal, calvals, woc(1,1), ldx, woc(1,2), &
                           ldy, voc(1,1), voc(1,2), work, thresh, nout )
            end if
         case(ixcopy)
            if( dosub(ixcopy) ) then
            call CTOMCOPY( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixdotc)
            if( dosub(ixdotc) ) then
            call CTOMDOTC( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), work, thresh, nout )
            end if
         case(ixdotu)
            if( dosub(ixdotu) ) then
            call CTOMDOTU( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), work, thresh, nout )
            end if
         case(ixnrm2)
            if( dosub(ixnrm2) ) then
            call CTOMNRM2( noc, nocrep, woc(1,1), ldx, voc(1,1), work, &
                           thresh, rtim, nout )
            end if
         case(ixrot )
            if( dosub(ixrot ) ) then
            call CTOMROT( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                          voc(1,1), voc(1,2), work, thresh, nout )
            end if
         case(ixsrot)
            if( dosub(ixsrot) ) then
            call CTOMSROT( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), work, thresh, nout )
            end if
         case(ixrotg)
            if( dosub(ixrotg) ) then
            call CTOMRTG( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                          woc(1,3), ldz, voc(1,1), voc(1,2), work, &
                          thresh, nout )
            end if
         case(ixscal)
            if( dosub(ixscal) ) then
            call CTOMSCAL( noc, nocrep, ncal, calvals, woc(1,1), ldx, &
                           voc(1,1), work, thresh, nout )
            end if
         case(ixsscal)
            if( dosub(ixsscal) ) then
            call CTOMSSCAL( noc, nocrep, nal, alvals, woc(1,1), ldx, &
                            voc(1,1), work, thresh, nout )
            end if
         case(ixswap)
            if( dosub(ixswap) ) then
            call CTOMSWAP( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), thresh, nout )
            end if
      end select
   end do
99 format( 'Invalid input value: ', A4, '=', I6, '; must be >=', I6 )
98 format( / A4, ':  Unrecognized path name' )
97 format( A, ' tests for N=', I8, ', NOFF=', I5, ', NREP=', I6, &
           /'----------------------------------------------------------' )
96 format( A6, ': N=', I8, ', time=', E15.8 )
95 format( A6, ': N=', I8, ', ALPHA=(', E15.8, ',', E15.8, '), time=', E15.8 )
94 format( A6, ': N=', I8, ', ALPHA=', E15.8, ', time=', E15.8 )
93 format( A6, ': N=', I8, ', type=', I4, ', time=', E15.8 )
end program
