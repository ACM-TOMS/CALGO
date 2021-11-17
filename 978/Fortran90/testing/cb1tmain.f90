program CB1TMAIN
!
!  Test program for the complex Level 1 BLAS.
!
!  The program must be driven by a short data file.  The first 8
!  records of the file are read using list-directed input, while the
!  remaining records are read as a single character string.  An
!  annotated example of a data file can be obtained by deleting the
!  first 3 characters from the following 19 lines (the case of the
!  character string in the last 11 lines is unimportant):
!  Level 1 BLAS test program
!  4                      Number of values of N
!  1 2 3 4                Values of N
!  4                      Number of values of INCX and INCY
!  1 -1 5 -5              Values of INCX and INCY
!  5                      Number of values of real ALPHA
!  0.0 1.0 -1.0 0.3 -0.3  Values of ALPHA
!  3                      Number of values of complex ALPHA
!  (1.0,0.0) (0.0,1.0) (0.3,-0.7)  Values of CALPHA
!  10.0                   Threshold value of test ratio
!  amax                   Test IxAMAX (index of largest element)
!  asum                   Test xASUM (1-norm of a vector)
!  axpy                   Test xAXPY (y <- alpha*x + y)
!  copy                   Test xCOPY (y <- x)
!  dot                    Test xDOTC and xDOTU (dot products)
!  nrm2                   Test xNRM2 (2-norm of a vector)
!  rot                    Test xROT (application of Givens rotation)
!  rotg                   Test xROTG (computation of Givens rotation)
!  scal                   Test xSCAL (scale a vector by a scalar)
!  ssq                    Test xLASSQ (scaled sum of squares)
!  swap                   Test xSWAP (swap two vectors)
!
!  The test of the computation of Givens rotations separately tests the
!  BLAS routine xROTG and the LAPACK routine xLARTG, producing two
!  lines of output.  The test of xLASSQ is also a test of an LAPACK
!  auxiliary routine not in the original BLAS.  The Level 1 BLAS
!  routines to compute and apply modified Givens rotations, xROTMG and
!  xROTM, are not included in this test program.
!
!  E. Anderson
!  May 3, 2016
!
!  .. Parameters ..
   use LA_CONSTANTS32, only: wp, one
   integer, parameter :: nin = 5, nout = 6
!  ..
!  .. Local Scalars ..
   character*4 :: path
   character*80 :: aline
   logical :: fatal
   integer :: i, in, ixmax, lenw, lenx, nal, ncal, nmax, nn, nix
   real(wp) :: thresh
!  ..
!  .. Local Arrays ..
   integer, allocatable, dimension(:) :: nvals, ixvals, iwork, iwork2
   real(wp), allocatable, dimension(:) :: alvals, work, work2
   complex(wp), allocatable, dimension(:) :: calvals, x1, x2, y1, y2
!  ..
!  .. External Functions ..
   logical :: LSAMEN
   external :: LSAMEN
!  ..
!
!  Read and skip the first line of the input file
!
   fatal = .false.
   read( nin, * )
   write( nout, 99 )
!
!  Read the values of N
!
   read( nin, * ) nn
   if( nn < 1 ) then
      write( nout, 98 )' nn ', nn, 1
      fatal = .true.
      read( nin, * )
   else
      allocate( nvals( nn ) )
      read( nin, * )( nvals( i ), i = 1, nn )
   end if
!
!  Skip any negative values of N
!
   nmax = 0
   in = 0
   do i = 1, nn
      if( nvals(i) < 0 ) then
         write( nout, 98 )' n  ', nvals(i), 0
      else
         in = in + 1
         if( in < i ) nvals(in) = nvals(i)
         nmax = max( nmax, nvals(i) )
      end if
   end do
   nmax = nmax + 1
   nn = in
   if( .not.fatal .and. nn <= 0 ) then
      write( nout, 98 )' nn ', nn, 1
      nn = 0
      fatal = .true.
   end if
   if( nn > 0 ) then
      write( nout, 97 )'N   ', ( nvals(i), i = 1, nn )
   end if
!
!  Read the values of INCX
!
   read( nin, * ) nix
   if( nix < 1 ) then
      write( nout, 98 )'nix ', nix, 1
      fatal = .true.
      read( nin, * )
   else
      allocate( ixvals( nix ) )
      read( nin, * )( ixvals( i ), i = 1, nix )
   end if
   ixmax = 0
   do i = 1, nix
      ixmax = max( ixmax, abs(ixvals(i)) )
   end do
   if( nix > 0 ) then
      write( nout, 97 )'INCX', ( ixvals(i), i = 1, nix )
   end if
!
!  Read the values of ALPHA
!
   read( nin, * ) nal
   if( nal < 1 ) then
      write( nout, 98 )'nal ', nal, 1
      fatal = .true.
      read( nin, * )
   else
      allocate( alvals( nal ) )
      read( nin, * )( alvals( i ), i = 1, nal )
   end if
   if( nal > 0 ) then
      write( nout, 93 )'ALPHA', ( alvals(i), i = 1, min(5,nal) )
   end if
   if( nal > 5 ) then
     write( nout, 92 )( alvals(i), i = 6, nal )
   end if
!
!  Read the values of CALPHA
!
   read( nin, * ) ncal
   if( ncal < 1 ) then
      write( nout, 98 )'ncal', ncal, 1
      fatal = .true.
      read( nin, * )
   else
      allocate( calvals( ncal ) )
      read( nin, * )( calvals( i ), i = 1, ncal )
   end if
   if( ncal > 0 ) then
      write( nout, 91 )'CALPHA', ( calvals(i), i = 1, min(3,ncal) )
   end if
   if( ncal > 5 ) then
     write( nout, 90 )( calvals(i), i = 4, ncal )
   end if
!
!  Allocate space for the work arrays
!
   lenx = max( nmax*ixmax, 1 )
   lenw = max( nmax, 1 )
   allocate( x1(1:lenx) )
   allocate( x2(1:lenx) )
   allocate( y1(1:lenx) )
   allocate( y2(1:lenx) )
   allocate( work(1:2*lenw) )
   allocate( work2(1:2*lenw) )
   allocate( iwork(1:lenx) )
   allocate( iwork2(1:lenx) )
!
!  Read the threshold value for the test ratios.
!
   read( nin, * )thresh
   write( nout, 96 )thresh
!
!  Quit if fatal errors were encountered
!
   if( fatal ) then
      write( nout, 94 )
      go to 30
   end if
!
!  Read a test path and call the appropriate test routine
!
10 continue
   read( nin, '(a)', end = 20 ) aline
   path = aline(1:4)
   if( lsamen( 4, path, 'amax' ) ) then
      call CB1TAMAX( nn, nvals, nix, ixvals, lenx, x1, x2, &
                     work, nout )
   else if( lsamen( 4, path, 'asum' ) ) then
      call CB1TASUM( nn, nvals, nix, ixvals, lenx, x1, y1, x2, &
                     work, thresh, nout )
   else if( lsamen( 4, path, 'axpy' ) ) then
      call CB1TAXPY( nn, nvals, nix, ixvals, ncal, calvals, lenx, &
                     x1, y1, x2, y2, work, work2, iwork, thresh, nout )
   else if( lsamen( 4, path, 'copy' ) ) then
      call CB1TCOPY( nn, nvals, nix, ixvals, lenx, x1, y1, x2, y2, &
                     work, iwork, thresh, nout )
   else if( lsamen( 4, path, 'dot ' ) ) then
      call CB1TDOTC( nn, nvals, nix, ixvals, lenx, x1, y1, x2, y2, &
                     work, work2, thresh, nout )
      call CB1TDOTU( nn, nvals, nix, ixvals, lenx, x1, y1, x2, y2, &
                     work, work2, thresh, nout )
   else if( lsamen( 4, path, 'nrm2' ) ) then
      call CB1TNRM2( nn, nvals, nix, ixvals, lenx, x1, x2, y1, &
                     work, thresh, nout )
   else if( lsamen( 4, path, 'rot ' ) ) then
      call CB1TROT( nn, nvals, nix, ixvals, lenx, x1, y1, x2, y2, &
                    work, work2, iwork, iwork2, thresh, nout )
      call CB1TSROT( nn, nvals, nix, ixvals, lenx, x1, y1, x2, y2, &
                     work, work2, iwork, iwork2, thresh, nout )
   else if( lsamen( 4, path, 'rotg' ) ) then
      call CB1TRTG( 'CROTG ', thresh, nout )
      call CB1TRTG( 'CLARTG', thresh, nout )
   else if( lsamen( 4, path, 'scal' ) ) then
      call CB1TSCAL( nn, nvals, nix, ixvals, ncal, calvals, lenx, &
                     x1, x2, work, iwork, thresh, nout )
      call CB1TSSCAL( nn, nvals, nix, ixvals, nal, alvals, lenx, &
                      x1, x2, work, iwork, thresh, nout )
   else if( lsamen( 4, path, 'ssq ' ) ) then
      call CB1TSSQ( nn, nvals, nix, ixvals, lenx, x1, y1, x2, &
                    work, thresh, nout )
   else if( lsamen( 4, path, 'swap' ) ) then
      call CB1TSWAP( nn, nvals, nix, ixvals, lenx, x1, y1, x2, y2, &
                     work, iwork, iwork2, thresh, nout )
   else
      write( nout, 95 ) path
   end if
   go to 10
20 continue
   write( nout, 89 )
30 continue
!
!  Format statements
!
99 format( ' Tests of the COMPLEX Level 1 BLAS', &
           // ' The following parameter values will be used:' )
98 format( ' Invalid input value: ', A4, '=', I6, '; must be >=', I6 )
97 format( 4X, A4, ':  ', 10I7, / 11X, 10I7 )
96 format( /, ' Routines pass computational tests if test ratio is ', &
           'less than', F8.2, // )
95 format( /1X, A4, ':  Unrecognized path name' )
94 format( ' Execution not attempted due to input errors' )
93 format( 2X, A6, ': ', 5( 1X, G12.5 ) )
92 format( 10X, 5( 1X, G12.5 ) )
91 format( 2X, A6, ': ', 3( 1X, '(', G12.5, ',', G12.5, ')' ) )
90 format( 9X, 3( 1X, '(', G12.5, ',', G12.5, ')' ) )
89 format( ' End of tests' )
end program
