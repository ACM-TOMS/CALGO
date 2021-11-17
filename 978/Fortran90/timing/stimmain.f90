program STIMMAIN
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
!  The BLAS routines ISAMAX, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
!  SROT, SSCAL, and SSWAP are timed in this procedure, as well
!  as the LAPACK routines SLAPY2 and SLARGV (the vector version
!  of the BLAS SROTG or LAPACK SLARTG).  The algorithm types that
!  can be specified in the input file are amax, asum, axpy, copy,
!  dot, nrm2, py2, rot, rotg, scal, and swap.
!
!  Output takes the form
!  In-cache test with N=10000, NOFF=100, NREP=100000
!  -------------------------------------------------
!  SAXPY : N=   10000, ALPHA= 0.10000000E+01, time= 0.99999993E-05
!  SAXPY : N=   10000, ALPHA=-0.10000000E+01, time= 0.99999993E-05
!  SAXPY : N=   10000, ALPHA= 0.30000000E+00, time= 0.99999993E-05
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
   integer, parameter :: ixdot  = 5
   integer, parameter :: ixnrm2 = 6
   integer, parameter :: ixrot  = 7
   integer, parameter :: ixrotg = 8
   integer, parameter :: ixscal = 9
   integer, parameter :: ixswap = 10
   integer, parameter :: ixpy2  = 11
   integer, parameter :: nsub   = 11
!  ..
!  .. Local Scalars ..
   logical :: fatal
   character*6 :: path
   character*80 :: aline
   integer :: i, ix, ixs, iy, iys, j, ldx, ldy, ldz, n, nal, noff, nrep, &
      noc, nocoff, nocrep
   real(wp) :: thresh
!  ..
!  .. Local Arrays ..
   logical :: dosub(nsub)
   character*6 :: subnam(nsub)
   real(wp), allocatable :: alvals(:), rtim(:), wic(:,:), woc(:,:), voc(:,:)
!  ..
!  .. External Functions ..
   logical :: lsamen
!  ..
!
   subnam(ixamax) = 'ISAMAX'
   subnam(ixasum) = 'SASUM '
   subnam(ixaxpy) = 'SAXPY '
   subnam(ixcopy) = 'SCOPY '
   subnam(ixdot)  = 'SDOT  '
   subnam(ixnrm2) = 'SNRM2 '
   subnam(ixrot)  = 'SROT  '
   subnam(ixrotg) = 'SLARGV'
   subnam(ixscal) = 'SSCAL '
   subnam(ixswap) = 'SSWAP '
   subnam(ixpy2)  = 'SLAPY2'
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
!  Read the values of alpha
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
   allocate( rtim( max(nal,9) ) )
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
   else if( lsamen( 4, path, 'amax' ) .or. lsamen( 6, path, 'ISAMAX' ) ) then
      dosub(ixamax) = .true.
   else if( lsamen( 4, path, 'asum' ) .or. lsamen( 5, path, 'SASUM' ) ) then
      dosub(ixasum) = .true.
   else if( lsamen( 4, path, 'axpy' ) .or. lsamen( 5, path, 'SAXPY' ) ) then
      dosub(ixaxpy) = .true.
   else if( lsamen( 4, path, 'copy' ) .or. lsamen( 5, path, 'SCOPY' ) ) then
      dosub(ixcopy) = .true.
   else if( lsamen( 4, path, 'dot ' ) .or. lsamen( 5, path, 'SDOT ' ) ) then
      dosub(ixdot ) = .true.
   else if( lsamen( 4, path, 'nrm2' ) .or. lsamen( 5, path, 'SNRM2' ) ) then
      dosub(ixnrm2) = .true.
   else if( lsamen( 4, path, 'rot ' ) .or. lsamen( 5, path, 'SROT ' ) ) then
      dosub(ixrot ) = .true.
   else if( lsamen( 4, path, 'rotg' ) .or. lsamen( 5, path, 'SROTG' ) .or. &
        lsamen( 6, path, 'SLARTG' ) .or. lsamen( 6, path, 'SLARGV' ) ) then
      dosub(ixrotg) = .true.
   else if( lsamen( 4, path, 'scal' ) .or. lsamen( 5, path, 'SSCAL' ) ) then
      dosub(ixscal) = .true.
   else if( lsamen( 4, path, 'swap' ) .or. lsamen( 5, path, 'SSWAP' ) ) then
      dosub(ixswap) = .true.
   else if( lsamen( 4, path, 'py2 ' ) .or. lsamen( 6, path, 'SLAPY2' ) ) then
      dosub(ixpy2) = .true.
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
            call STIMAMAX( n, nrep, wic(1,1), wic(1,2), rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixasum)
            if( dosub(ixasum) ) then
            call STIMASUM( n, nrep, wic(1,1), wic(1,2), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixaxpy)
            if( dosub(ixaxpy) ) then
            call STIMAXPY( n, nrep, nal, alvals, wic(1,1), wic(1,2), &
                           wic(1,3), wic(1,4), thresh, rtim, nout )
            do j = 1, nal
              write(nout,95) subnam(i), n, alvals(j), rtim(j)
            end do
            end if
         case(ixcopy)
            if( dosub(ixcopy) ) then
            call STIMCOPY( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixdot )
            if( dosub(ixdot ) ) then
            call STIMDOT( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                          wic(1,4), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixnrm2)
            if( dosub(ixnrm2) ) then
            call STIMNRM2( n, nrep, wic(1,1), wic(1,2), thresh, rtim, nout )
            do j = 1, 8
               write(nout,94) subnam(i), n, j, rtim(j)
            end do
            end if
         case(ixrot )
            if( dosub(ixrot ) ) then
            call STIMROT( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                          wic(1,4), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixrotg)
            if( dosub(ixrotg) ) then
            call STIMRTG( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                          wic(1,4), wic(1,5), thresh, rtim, nout )
            do j = 1, 9
               write(nout,94) subnam(i), n, j, rtim(j)
            end do
            end if
         case(ixscal)
            if( dosub(ixscal) ) then
            call STIMSCAL( n, nrep, nal, alvals, wic(1,1), wic(1,2), &
                           thresh, rtim, nout )
            do j = 1, nal
              write(nout,95) subnam(i), n, alvals(j), rtim(j)
            end do
            end if
         case(ixswap)
            if( dosub(ixswap) ) then
            call STIMSWAP( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                           wic(1,4), thresh, rtim, nout )
            write(nout,96) subnam(i), n, rtim(1)
            end if
         case(ixpy2)
            if( dosub(ixpy2) ) then
            call STIMPY2( n, nrep, wic(1,1), wic(1,2), wic(1,3), &
                          wic(1,4), wic(1,5), thresh, rtim, nout )
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
            call STOMAMAX( noc, nocrep, woc(1,1), ldx, voc(1,1), nout )
            end if
         case(ixasum)
            if( dosub(ixasum) ) then
            call STOMASUM( noc, nocrep, woc(1,1), ldx, voc(1,1), thresh, nout )
            end if
         case(ixaxpy)
            if( dosub(ixaxpy) ) then
            call STOMAXPY( noc, nocrep, nal, alvals, woc(1,1), ldx, woc(1,2), &
                           ldy, voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixcopy)
            if( dosub(ixcopy) ) then
            call STOMCOPY( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixdot )
            if( dosub(ixdot ) ) then
            call STOMDOT( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                          voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixnrm2)
            if( dosub(ixnrm2) ) then
            call STOMNRM2( noc, nocrep, woc(1,1), ldx, voc(1,1), thresh, &
                           rtim, nout )
            end if
         case(ixrot )
            if( dosub(ixrot ) ) then
            call STOMROT( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                          voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixrotg)
            if( dosub(ixrotg) ) then
            call STOMRTG( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                          woc(1,3), ldz, voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixscal)
            if( dosub(ixscal) ) then
            call STOMSCAL( noc, nocrep, nal, alvals, woc(1,1), ldx, &
                           voc(1,1), thresh, nout )
            end if
         case(ixswap)
            if( dosub(ixswap) ) then
            call STOMSWAP( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                           voc(1,1), voc(1,2), thresh, nout )
            end if
         case(ixpy2)
            if( dosub(ixpy2) ) then
            call STOMPY2( noc, nocrep, woc(1,1), ldx, woc(1,2), ldy, &
                          woc(1,3), ldz, voc(1,1), voc(1,2), thresh, &
                          nout )
            end if
      end select
   end do
99 format( 'Invalid input value: ', A4, '=', I6, '; must be >=', I6 )
98 format( / A4, ':  Unrecognized path name' )
97 format( A, ' tests for N=', I8, ', NOFF=', I5, ', NREP=', I6, &
           /'----------------------------------------------------------' )
96 format( A6, ': N=', I8, ', time=', E15.8 )
95 format( A6, ': N=', I8, ', ALPHA=', E15.8, ', time=', E15.8 )
94 format( A6, ': N=', I8, ', type=', I4, ', time=', E15.8 )
end program
