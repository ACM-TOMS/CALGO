program pdttrandutv
  !
  implicit none
  !
  ! Include MPI headers.
  include 'mpif.h'
  !
  ! .. Local Scalars ..
  integer            ictxt, iam, nprocs, mycol, myrow, npcol, nprow, &
                     m, n, nb, q_factor, print_matrices, &
                     ierr, len_processor_name
  double precision   wtime, resid
  character*( MPI_MAX_PROCESSOR_NAME ) processor_name
  ! ..
  ! .. Local Arrays ..
  integer            iseed( 4 )
  !
  ! .. External Subroutines ..
  external           blacs_barrier, blacs_exit, blacs_get, blacs_gridexit, & 
                     blacs_gridinfo, blacs_gridinit, blacs_pinfo, &
                     mpi_get_processor_name, pd_test_and_time_randutv
  ! ..
  ! .. Executable Statements ..

  !
  ! Get starting information.
  !
  call blacs_pinfo( iam, nprocs )

  !
  ! Set initial values for execution.
  ! =================================
  !
  ! Grid dimensions.
  nprow = 2
  npcol = 2

  ! Matrix and block size dimensions (m >= n).
  m   =  6
  n   =  6
  nb  =  2

  ! Q factor for power iterations.
  q_factor = 2

  ! Print matrices (1=yes;0=no).
  print_matrices = 1

  ! Seeds for random number generation.
  iseed( 1 ) = 2 
  iseed( 2 ) = 5
  iseed( 3 ) = 8
  iseed( 4 ) = 1

  ! 
  ! Check matrix dimensions.
  ! 
  if( m < n ) then
    if( iam == 0 ) then
      write( *, * )
      write( *, * ) '*** ERROR in pdttrandutv:  m should be >= n.'
      write( *, * ) 'm:  ', m
      write( *, * ) 'n:  ', n
    end if
    stop
  end if

  ! 
  ! Check grid size.
  ! 
  if( ( nprow * npcol ) > nprocs ) then
    if( iam == 0 ) then
      write( *, * )
      write( *, * ) '*** ERROR in pdttrandutv: ', & 
                    'Grid size requires too many processes.'
      write( *, * ) 'nprow:  ', nprow
      write( *, * ) 'npcol:  ', npcol
      write( *, * ) 'nprocs: ', nprocs
    end if
    stop
  end if

  !
  ! Define process grid.
  !
  call blacs_get( -1, 0, ictxt )
  call blacs_gridinit( ictxt, 'row-major', nprow, npcol )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

  !
  ! Processes must skip the processing if they are not included in the grid.
  !
  if( ( myrow == -1 ).or.( mycol == -1 ) ) then
    write( *, '(1x, a, 3(i4), a)' ) 'Process ', iam, myrow, mycol, &
        ' is excluded from the grid.'
    go to 20
  end if

  !
  ! Print initial messages.
  !
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,78("="))' )
    write( *, * ) 'Testing and timing randUTV'
    write( *, '(1x,78("="))' )
    write( *, '(/,1x,a,i3,a,i3)' ) 'Grid size:   ', nprow, ' x ', npcol
    write( *, * ) 'm:  ', m
    write( *, * ) 'n:  ', n
    write( *, * ) 'nb: ', nb
    write( *, * )
  end if

  !
  ! Every process prints its processor name.
  !
  call blacs_barrier( ictxt, 'all' )
  call mpi_get_processor_name( processor_name, len_processor_name, ierr )
  write( *, '(1x, a, 3(i4), 2(a) )' ) 'Proc: ', iam, myrow, mycol, &
      ' on: ', trim( processor_name )


  call blacs_barrier( ictxt, 'all' )
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, * )
  end if
  call blacs_barrier( ictxt, 'all' )

  !
  ! Test and time the randUTV factorization.
  !
  call pd_test_and_time_randutv( ictxt, m, n, nb, q_factor, iseed, &
           print_matrices, wtime, resid )

  call blacs_gridexit( ictxt )

 20 continue
  !
  ! Print final message.
  !
  if( iam == 0 ) then
     write( *, '(/,1x,a)' ) 'End of program'
  end if

  call blacs_exit( 0 )
  stop

  !
  ! End of pdttrandutv
  !
end

! =============================================================================
subroutine pd_test_and_time_randutv( ictxt, m, n, nb, q_factor, iseed, &
               print_matrices, wtime, resid )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            ictxt, q_factor, m, n, nb, print_matrices
  double precision   wtime, resid
  ! ..
  ! .. Array Arguments ..
  integer            iseed( * )
  ! ..
  !
  ! Purpose
  ! =======
  !
  ! It tests and times the randUTV factorization.
  !
  ! =====================================================================
  !
  ! .. Parameters ..
  integer            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                     LLD_, MB_, M_, NB_, N_, RSRC_
  parameter          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                       CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                       RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
  double precision   ONE, ZERO, TWO, BILLION
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0, TWO = 2.0D+0, & 
                       BILLION = 1.0D+9 )
  ! ..
  ! .. Local Arrays ..
  integer            desca( DLEN_ ), descac( DLEN_ ), &
                     descu( DLEN_ ), descv( DLEN_ ), &
                     niseed( 4 )
  double precision, dimension ( : ), allocatable :: a, ac, u, v
  ! ..
  ! .. Local Scalars ..
  integer            nprow, npcol, myrow, mycol, &
                     npa, nqa, npac, nqac, npu, nqu, npv, nqv, &
                     j, info, allocstat
  ! ..
  ! .. External Functions ..
  integer            numroc
  double precision   check_utv_resids
  external           numroc, check_utv_resids
  ! ..
  ! .. External Subroutines ..
  external           blacs_barrier, blacs_gridinfo, descinit, &
                     generate_random_matrix, pdgeutv, pdlacpy, &
                     print_distributed_matrix, slboot, slcombine, sltimer
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          int, max, dble
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

  !
  ! Some initializations.
  !
  npa   = numroc( m, nb, myrow, 0, nprow )
  nqa   = numroc( n, nb, mycol, 0, npcol )
  npac  = numroc( m, nb, myrow, 0, nprow )
  nqac  = numroc( n, nb, mycol, 0, npcol )
  npu   = numroc( m, nb, myrow, 0, nprow )
  nqu   = numroc( m, nb, mycol, 0, npcol )
  npv   = numroc( n, nb, myrow, 0, nprow )
  nqv   = numroc( n, nb, mycol, 0, npcol )

  !
  ! Initialize array descriptors for the matrices.
  !
  call descinit( desca,  m, n, nb, nb, 0, 0, ictxt, max( 1, npa  ), info )
  call descinit( descac, m, n, nb, nb, 0, 0, ictxt, max( 1, npac ), info )
  call descinit( descu,  m, m, nb, nb, 0, 0, ictxt, max( 1, npu  ), info )
  call descinit( descv,  n, n, nb, nb, 0, 0, ictxt, max( 1, npv  ), info )

  !
  ! Allocate local arrays.
  !
  allocate( a( npa * nqa ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for a ***'

  allocate( ac( npac * nqac ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for ac ***'

  allocate( u( npu * nqu ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for u ***'

  allocate( v( npv * nqv ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for v ***'

  !
  ! Generate matrices A and AC.
  !
  do j = 1, 4
    niseed( j ) = iseed( j )
  end do
  call generate_random_matrix( 1, m, n, niseed, a, desca )
  call pdlacpy( 'all', m, n, a, 1, 1, desca, &
                             ac, 1, 1, descac )

  !
  ! Print initial matrix A.
  !
  if( print_matrices == 1 ) then
    call print_distributed_matrix( m, n, a, 1, 1, desca, 'ai' )
  end if

  !
  ! Factorize matrix and generate orthonormal matrices U and V. Get times.
  !
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, * ) 'Starting pdgeutv...'
  end if

  call slboot()
  call blacs_barrier( ictxt, 'all' )
  call sltimer( 1 )

  call pdgeutv( m, n, a, desca, 1, u, descu, 1, v, descv, & 
                q_factor, info )

  call sltimer( 1 )
  call slcombine( ictxt, 'all', '>', 'w', 1, 1, wtime )

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, * ) 'End of pdgeutv.'
  end if

  !
  ! Check info argument returned.
  !
  if( info /= 0 ) then
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write( *, * )
      write( *, * ) 'Info code returned by pdgeutv = ', info
      write( *, * )
    end if
  end if

  !
  ! Print final matrices.
  !
  if( print_matrices == 1 ) then
    call print_distributed_matrix( m, n, a, 1, 1, desca, 'af' )
    call print_distributed_matrix( m, m, u, 1, 1, descu, 'uf' )
    call print_distributed_matrix( n, n, v, 1, 1, descv, 'vf' )
  end if

  !
  ! Check residuals.
  !
  if( info == 0 ) then
    resid = check_utv_resids( m, n, ac, descac, u, descu, a, desca, v, descv )
  else 
    resid = - ONE
  end if 

  ! 
  ! Print time.
  ! 
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, * )
    write( *, '(1x,a,f12.5)' ) 'Time (in s.) spent in the factorization: ', &
                               wtime
  end if

  !
  ! Deallocate local arrays and vectors.
  !
  deallocate( a )
  deallocate( ac )
  deallocate( u )
  deallocate( v )

  return
  !
  ! End of pd_test_and_time_randutv
  !
end

! =============================================================================
double precision function check_utv_resids( m, n, ac, descac, u, descu, & 
                                            t, desct, v, descv )
  !
  implicit none
  !
  ! .. Scalar arguments ..
  integer            m, n
  ! ..
  ! .. Array arguments ..
  integer            descac( * ), descu( * ), desct( * ), descv( * )
  double precision   ac( * ), u( * ), t( * ), v( * )
  ! ..
  !
  ! Purpose
  ! =======
  !
  ! It checks the residuals of the randUTV factorization.
  !
  ! =====================================================================
  !
  ! .. Parameters ..
  integer            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, &
                     LLD_, MB_, M_, NB_, N_, RSRC_
  parameter          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1, &
                       CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                       RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
  double precision   ONE, ZERO
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  ! ..
  ! .. Local arrays ..
  integer            descut( DLEN_ ), desctt( DLEN_ ), &
                     descmm( DLEN_ ), descmn( DLEN_ ), descnn( DLEN_ )
  double precision, dimension ( : ), allocatable :: ut, tt, mm, mn, nn
  ! ..
  ! .. Local scalars ..
  double precision   nrm, resida, residb, residc, residd, reside, residf, &
                     maxres
  integer            ictxt, nb, nprow, npcol, myrow, mycol, &
                     info, allocstat, &
                     nput, nqut, nptt, nqtt, & 
                     npmm, nqmm, npmn, nqmn, npnn, nqnn, dummy
  ! ..
  ! .. External Functions ..
  integer            numroc
  double precision   pdlange, compar_svd
  external           numroc, pdlange, compar_svd
  ! ..
  ! .. External Subroutines ..
  external           blacs_gridinfo, descinit, pdgemm, pdlacpy, pdlaset
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          int, max
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = descac( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
  nb    = descac( NB_ )

  !
  ! Some initializations.
  !
  nput  = numroc( m, nb, myrow, 0, nprow )
  nqut  = numroc( n, nb, mycol, 0, npcol )
  nptt  = numroc( m, nb, myrow, 0, nprow )
  nqtt  = numroc( n, nb, mycol, 0, npcol )
  npmm  = numroc( m, nb, myrow, 0, nprow )
  nqmm  = numroc( m, nb, mycol, 0, npcol )
  npmn  = numroc( m, nb, myrow, 0, nprow )
  nqmn  = numroc( n, nb, mycol, 0, npcol )
  npnn  = numroc( n, nb, myrow, 0, nprow )
  nqnn  = numroc( n, nb, mycol, 0, npcol )

  !
  ! Initialize array descriptors for the matrices MM, and MN.
  !
  call descinit( descut, m, n, nb, nb, 0, 0, ictxt, max( 1, nput ), info )
  call descinit( desctt, m, n, nb, nb, 0, 0, ictxt, max( 1, nptt ), info )
  call descinit( descmm, m, m, nb, nb, 0, 0, ictxt, max( 1, npmm ), info )
  call descinit( descmn, m, n, nb, nb, 0, 0, ictxt, max( 1, npmn ), info )
  call descinit( descnn, n, n, nb, nb, 0, 0, ictxt, max( 1, npnn ), info )

  !
  ! Allocate local arrays.
  !
  allocate( ut( nput * nqut ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for ut ***'

  allocate( tt( nptt * nqtt ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for tt ***'

  allocate( mm( npmm * nqmm ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for mm ***'

  allocate( mn( npmn * nqmn ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for mn ***'

  allocate( nn( npnn * nqnn ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for nn ***'

  !
  ! Print initial message.
  !
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, * )
    write( *, '(/,1x,a)' ) 'Residuals of randUTV factorization'
  end if

  !
  ! Compute residual: || Ac - U * triu( T ) * V' || / || Ac ||.
  !
  call pdlacpy( 'All', m, n, t, 1, 1, desct, &
                             tt, 1, 1, desctt )
  call pdlaset( 'Lower', m - 1, n, ZERO, ZERO, tt, 2, 1, desctt )

  call pdgemm( 'No transpose', 'No transpose', m, n, m, &
               ONE, u, 1, 1, descu, tt, 1, 1, desctt, &
               ZERO,  ut, 1, 1, descut )

  call pdlacpy( 'All', m, n, ac, 1, 1, descac, &
                             mn, 1, 1, descmn )
  call pdgemm( 'No transpose', 'Transpose', m, n, n, &
               -ONE, ut, 1, 1, descut, v, 1, 1, descv, &
               ONE,  mn, 1, 1, descmn )
  !
  resida = pdlange( 'Frobenius', m, n, mn, 1, 1, descmn, dummy )
  nrm = pdlange( 'Frobenius', m, n, ac, 1, 1, descac, dummy )
  if( nrm /= ZERO ) then
    resida = resida / nrm
  end if 

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,a,e12.4)' ) &
        '|| Ac - U * triu( T ) * Vt || / || Ac || =  ', resida
  end if

  !
  ! Compute residual: || I - U' * U || / || U ||.
  !
  call pdlaset( 'All', m, m, ZERO, ONE, mm, 1, 1, descmm )
  call pdgemm( 'Transpose', 'No transpose', m, m, m, &
               -ONE, u, 1, 1, descu, u, 1, 1, descu, &
               ONE,  mm, 1, 1, descmm )
  !
  residb = pdlange( 'Frobenius', m, m, mm, 1, 1, descmm, dummy )
  nrm = pdlange( 'Frobenius', m, m, u, 1, 1, descu, dummy )
  if( nrm /= ZERO ) then
    residb = residb / nrm
  end if 

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,a,e12.4)' ) &
        '|| I - Ut * U || / || U || =                ', residb
  end if

  !
  ! Compute residual: || I - U * U' || / || U ||.
  !
  call pdlaset( 'All', m, m, ZERO, ONE, mm, 1, 1, descmm )
  call pdgemm( 'No transpose', 'Transpose', m, m, m, &
               -ONE, u, 1, 1, descu, u, 1, 1, descu, &
               ONE,  mm, 1, 1, descmm )
  !
  residc = pdlange( 'Frobenius', m, m, mm, 1, 1, descmm, dummy )
  nrm = pdlange( 'Frobenius', m, m, u, 1, 1, descu, dummy )
  if( nrm /= ZERO ) then
    residc = residc / nrm
  end if 

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,a,e12.4)' ) &
        '|| I - U * Ut || / || U || =                ', residc
  end if

  !
  ! Compute residual: || I - V' * V || / || V ||.
  !
  call pdlaset( 'All', n, n, ZERO, ONE, nn, 1, 1, descnn )
  call pdgemm( 'Transpose', 'No transpose', n, n, n, &
               -ONE, v, 1, 1, descv, v, 1, 1, descv, &
               ONE,  nn, 1, 1, descnn )
  !
  residd = pdlange( 'Frobenius', n, n, nn, 1, 1, descnn, dummy )
  nrm = pdlange( 'Frobenius', n, n, v, 1, 1, descv, dummy )
  if( nrm /= ZERO ) then
    residd = residd / nrm
  end if 

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,a,e12.4)' ) &
        '|| I - Vt * V || / || V || =                ', residd
  end if

  !
  ! Compute residual: || I - V * V' || / || V ||.
  !
  call pdlaset( 'All', n, n, ZERO, ONE, nn, 1, 1, descnn )
  call pdgemm( 'No transpose', 'Transpose', n, n, n, &
               -ONE, v, 1, 1, descv, v, 1, 1, descv, &
               ONE,  nn, 1, 1, descnn )
  !
  reside = pdlange( 'Frobenius', n, n, nn, 1, 1, descnn, dummy )
  nrm = pdlange( 'Frobenius', n, n, v, 1, 1, descv, dummy )
  if( nrm /= ZERO ) then
    reside = reside / nrm
  end if 

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,a,e12.4)' ) &
        '|| I - V * Vt || / || V || =                ', reside
  end if

  !
  ! Compute residual: || singular_values( ac ) - singular_values( triu( t ) ||
  !                   / || singular_values( ac ) ||.
  !
  residf = compar_svd( m, n, ac, descac, t, desct )

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, '(1x,a,e12.4)' ) &
        '|| svd( a ) - svd( t )|| / || svd( a ) || = ', residf
  end if

  !
  ! Deallocate local arrays and vectors.
  !
  deallocate( ut )
  deallocate( tt )
  deallocate( mm )
  deallocate( mn )
  deallocate( nn )

  !
  ! Compute and print maximum residual.
  !
  maxres = max( resida, residb, residc, residd, reside, residf )

  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, * )
    write( *, '(1x,a,e12.4)' ) 'Maximum residual: ', maxres
  end if

  check_utv_resids = maxres
  return
  !
  ! End of check_utv_resids
  !
end

! =============================================================================
double precision function compar_svd( m, n, a1, desca1, a2, desca2 )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            m, n
  ! ..
  ! .. Array Arguments ..
  integer            desca1( * ), desca2( * )
  double precision   a1( * ), a2( * )
  ! ..
  !
  ! Purpose
  ! =======
  !
  ! It computes the following:
  !   || singular_values( a1 ) - singular_values( a2 ) || / 
  !   || singular_values( a1 ) ||.
  !
  ! =====================================================================
  !
  ! .. Parameters ..
  integer            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                     LLD_, MB_, M_, NB_, N_, RSRC_
  parameter          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                       CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                       RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
  double precision   ONE, ZERO
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  ! ..
  ! .. Local Arrays ..
  double precision, dimension ( : ), allocatable :: ac, vs1, vs2, work
  integer            descac( dlen_ )
  ! ..
  ! .. Local Scalars ..
  integer            ictxt, nprow, npcol, myrow, mycol, &
                     nb, npac, nqac, allocstat, len_work, info
  double precision   scalar_work, nrm, res
  ! ..
  ! .. External Functions ..
  integer            numroc
  double precision   dnrm2
  external           numroc, dnrm2
  ! ..
  ! .. External Subroutines ..
  external           blacs_gridinfo, daxpy, descinit, pdgesvd, pdlacpy
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          int, max
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = desca1( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
  nb    = desca1( NB_ )

  !
  ! Some initializations.
  !
  npac  = numroc( m, nb, myrow, 0, nprow )
  nqac  = numroc( n, nb, mycol, 0, npcol )
  call descinit( descac, m, n, nb, nb, 0, 0, ictxt, max( 1, npac ), info )

  !
  ! Compute optimal real workspace length.
  !
  call pdgesvd( 'None', 'None', m, n, a1, 1, 1, desca1, &
                vs1, a1, 1, 1, desca1, a1, 1, 1, desca1, &
                scalar_work, -1, info )
  if( info /= 0 ) then
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write ( *, * ) '*** Error in compar_svd:', &
                     'Info after query-request  pdgesvd: ', info
    end if
  end if
  len_work = int( scalar_work )

  !
  ! Allocate local arrays and vectors.
  !
  allocate( ac( npac * nqac ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for ac ***'

  allocate( vs1( min( m, n ) ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for vs1 ***'

  allocate( vs2( min( m, n ) ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for vs2 ***'

  allocate( work( len_work ), stat = allocstat )
  if ( allocstat /= 0 ) stop '*** Not enough memory for work ***'

  !
  ! Compute singular values of a1.
  !
  call pdlacpy( 'All', m, n, a1, 1, 1, desca1, &
                             ac, 1, 1, descac )
  call pdgesvd( 'None', 'None', m, n, ac, 1, 1, descac, &
                vs1, ac, 1, 1, descac, ac, 1, 1, descac, &
                work, len_work, info )
  if( info /= 0 ) THEN
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write ( *, * ) '*** Error in compar_svd:', &
                     'Info after first pdgesvd: ', info
    end if
  end if

  !
  ! Compute singular values of a2.
  !
  call pdlacpy( 'All', m, n, a2, 1, 1, desca2, &
                             ac, 1, 1, descac )
  call pdgesvd( 'None', 'None', m, n, ac, 1, 1, descac, &
                vs2, ac, 1, 1, descac, ac, 1, 1, descac, &
                work, len_work, info )
  if( info /= 0 ) then
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write ( *, * ) '*** Error in compar_svd:', &
                     'Info after second pdgesvd: ', info
    end if
  end if

  !
  ! Compute the norm of the singular values.
  !
  nrm = dnrm2( min( m, n ), vs1, 1 )

  !
  ! Compute the norm of the difference between the singular values.
  !
  call daxpy( min( m, n ), - ONE, vs2, 1, vs1, 1 )
  res = dnrm2( min( m, n ), vs1, 1 )

  !
  ! Compute the relative residual.
  !
  if( nrm /= ZERO ) then
    compar_svd = res / nrm
  else 
    compar_svd = res
  end if 
  !
  ! Deallocate local arrays and vectors.
  !
  deallocate( ac )
  deallocate( vs1 )
  deallocate( vs2 )
  deallocate( work )

  return
  !
  ! End of compar_svd
  !
end

! =============================================================================
subroutine generate_random_matrix( matrix_type, m, n, iseed, a, desca )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            matrix_type, m, n
  ! ..
  ! .. Array Arguments ..
  integer            desca( * ), iseed( * )
  double precision   a( * )
  ! ..
  ! 
  ! Purpose
  ! =======
  ! 
  ! It generates random double-precision m-by-n matrix a.
  ! This code should not be used for large matrices since it is not
  ! efficient.
  !
  ! =====================================================================
  !
  ! .. Parameters ..
  integer            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                     LLD_, MB_, M_, NB_, N_, RSRC_
  parameter          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                       CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                       RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
  double precision   ONE, ZERO
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  ! ..
  ! .. Local Scalars ..
  integer            ictxt, mycol, myrow, npcol, nprow, i, j
  double precision   scale_factor, elem
  ! ..
  ! .. External Subroutines ..
  external           blacs_gridinfo, dlarnv, pdelset
  ! ..
  ! .. Executable Statements ..

  ! 
  ! Get grid parameters.
  ! 
  ictxt = desca( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

  !
  ! Check matrix_type.
  !
  if( matrix_type == 1 ) then
    ! 
    ! Fill A with random numbers.
    ! 
    do j = 1, n
      do i = 1, m
        call dlarnv( 1, iseed, 1, elem )
        call pdelset( a, i, j, desca, elem )
      end do
    end do

  else
    ! 
    ! Fill A with integer entries converted to double precision, and 
    ! scaled down.
    ! 
    if( ( m == 0 ).or.( n == 0 ) ) then
      scale_factor = ONE
    else
      scale_factor = ONE / ( dble( m ) * dble( n ) )
    end if
    elem = ZERO
    do j = 1, n
      do i = j, m
        elem = elem + ONE
        call pdelset( a, i, j, desca, elem * scale_factor )
      end do
      do i = 1, j - 1
        elem = elem + ONE
        call pdelset( a, i, j, desca, elem * scale_factor )
      end do
    end do
    !
    if( ( m > 0 ).and.( n > 0 ) ) then
      elem = ONE + 0.1
      call pdelset( a, 1, 1, desca, elem * scale_factor )
    end if

  end if
  ! 
  return
  ! 
  ! End of generate_random_matrix
  ! 
end

! =============================================================================
subroutine print_distributed_matrix( m, n, a, ia, ja, desca, cname )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            m, n, ia, ja
  character*( * )    cname
  ! ..
  ! .. Array Arguments ..
  integer            desca( * )
  double precision   a( * )
  ! ..
  !
  ! Purpose
  ! =======
  !
  ! It prints the distributed matrix.
  !
  ! =====================================================================
  !
  ! .. Parameters ..
  integer            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                     LLD_, MB_, M_, NB_, N_, RSRC_
  parameter          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                       CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                       RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
  double precision   ONE, ZERO
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  ! ..
  ! .. Local Scalars ..
  integer            ictxt, myrow, mycol, nprow, npcol, &
                     i, j
  double precision   alpha
  ! ..
  ! .. External Subroutines ..
  external           blacs_barrier, blacs_gridinfo, pdelget
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = desca( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

  ! Initial barrier.
  call blacs_barrier( ictxt, 'All' )

  ! Print matrix heading.
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, fmt = '(a,a)' ) cname, ' = [ '
  end if

  ! Main loops.
  do i = ia, ia + m - 1
    do j = ja, ja + n - 1

      ! Get and print element (i,j).
      call pdelget( 'All', ' ', alpha, a, i, j, desca )
      if( ( myrow == 0 ).and.( mycol == 0 ) ) then
        write( *, fmt = '(1x,e14.7)', advance = 'no' ) alpha
      end if

    end do
    ! Print end of line for every row.
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write ( *, * )
    end if
  end do

  ! Print end of matrix.
  if( ( myrow == 0 ).and.( mycol == 0 ) ) then
    write( *, fmt = '(a)' ) '];'
  end if

  ! Final barrier.
  call blacs_barrier( ictxt, 'All' )

  return
  !
  ! End of print_distributed_matrix
  !
end

! =============================================================================
subroutine print_local_matrix( m, n, a, lda, cname )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            m, n, lda
  character*(*)      cname
  ! ..
  ! .. Array Arguments ..
  double precision   a( lda, n )
  ! ..
  !
  ! Purpose
  ! =======
  !
  ! It prints the local matrix in Matlab/Octave format.
  !
  ! =====================================================================
  !
  ! .. Local Scalars ..
  integer            i, j
  ! ..
  ! .. Executable Statements ..

  ! Print matrix heading.
  write( *, fmt = '(a,a)' ) cname, ' = [ '

  ! Main loops.
  do i = 1, m
    do j = 1, n
        write( *, fmt = '(1x,e12.5)', advance = 'no' ) a( i, j )
    end do
    write ( *, * )
  end do

  ! Print end of matrix.
  write( *, fmt = '(a)' ) '];'

  return
  !
  ! End of print_local_matrix
  !
end

