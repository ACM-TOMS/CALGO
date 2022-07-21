! =============================================================================
subroutine pdgeutv( m, n, a, desca, build_u, u, descu, build_v, v, descv, &
                    q, info )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            m, n, q, info
  logical            build_u, build_v
  ! ..
  ! .. Array Arguments ..
  integer            desca( * ), descu( * ), descv( * )
  double precision   a( * ), u( * ), v( * )
  ! ..
  !
  ! PURPOSE
  !
  ! To compute a randUTV factorization
  !     A = U * T * V'.
  !
  ! ARGUMENTS
  !
  ! Input/Output Parameters
  !
  ! m       (global input) integer
  !         The number of rows to be operated on, i.e. the number of
  !         rows of the distributed matrix A. m >= 0. m >= n.
  !
  ! n       (global input) integer
  !         The number of columns to be operated on, i.e. the number
  !         of columns of the distributed matrix A. n >= 0. m >= n.
  !         NOTE: In case m >> n, an initial triangularization of matrix A 
  !         with the QR factorization and then the randUTV factorization of 
  !         the triangular factor could improve performances. This technique 
  !         is used in some linear factorizations such as the SVD when m >> n.
  !
  ! a       (local input/local output) double precision pointer into
  !         the local memory to an array, dimension (LLD_a,LOCc(n))
  !         On entry, the local pieces of the m-by-n distributed
  !         matrix A which is to be factored.
  !         On exit, the elements on and above the diagonal of A
  !         contain the min(m,n) by n upper trapezoidal matrix T (T
  !         is upper triangular if m >= n); the elements below the
  !         diagonal are zeros.
  !
  ! desca   (global and local input) integer array, dimension DLEN_
  !         The array descriptor for the distributed matrix A.
  !
  ! build_u (global input) logical
  !         If it is true, matrix U must be build.
  !
  ! u       (local output) double precision pointer into
  !         the local memory to an array, dimension (LLD_u,LOCc(m))
  !         On exit, matrix u contains the orthonormal m x m matrix U.
  !
  ! descu   (global and local input) integer array, dimension DLEN_
  !         The array descriptor for the distributed matrix U.
  !
  ! build_v (global input) logical
  !         If it is true, matrix V must be build.
  !
  ! v       (local output) double precision pointer into
  !         the local memory to an array, dimension (LLD_v,LOCc(n))
  !         On exit, matrix v contains the orthonormal n x n matrix V.
  !
  ! descv   (global and local input) integer array, dimension DLEN_
  !         The array descriptor for the distributed matrix V.
  !
  ! q       (global input) integer
  !         The number of iterations for the iterative power process.
  !
  ! info    (global input) integer
  !
  ! Error Indicator
  !
  ! info    (global output) integer
  !         = 0: Successful exit.
  !         < 0: If the i-th argument is an array and the j-entry had
  !              an illegal value, then info = -(i*100+j), if the
  !              i-th argument is a scalar and had an illegal value,
  !              then info = -i.
  !
  ! PARALLEL EXECUTION RECOMMENDATIONS
  !
  ! Restrictions:
  ! o The distribution blocks must be square (MB=NB).
  !
  ! ******************************************************************
  !
  ! .. Parameters ..
  integer            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                     LLD_, MB_, M_, NB_, N_, RSRC_
  parameter          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                       CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                       RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
  double precision   ONE, ZERO
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  logical            PROFILE, PRINT_INTERMEDIATE_MATRICES
  parameter          ( PROFILE = .FALSE., & 
                       PRINT_INTERMEDIATE_MATRICES = .FALSE. )
  ! ..
  ! .. Local Scalars ..
  integer            ictxt, myrow, mycol, nprow, npcol, & 
                     nb, lda, j, k, mn, nbrow, nbcol, &
                     mpa, nqa, mpg, mpy, nqg, nqy, &
                     idx_row_of_j, idx_col_of_j, & 
                     proc_row_of_j, proc_col_of_j, &
                     ii, offjj, allocstat, seed_size
  double precision   t_normal, t_power, t_rgeqr, t_rormq, t_lgeqr, t_lormq, &
                     t_svd, t_mma12, t_mma01
  ! ..
  ! .. Local Arrays ..
  integer            descg( DLEN_ ), descy( DLEN_ ), &
                     descsvdu( DLEN_ ), descsvdvt( DLEN_ ), idum( 1 )
  double precision, dimension( : ), allocatable :: g, y, vtau, & 
                     svda, svdu, svds, svdvt, s, work_larft
  integer, dimension( : ), allocatable :: seeds
  ! ..
  ! .. External Functions ..
  integer            numroc
  external           numroc
  ! ..
  ! .. External subroutines ..
  external           blacs_barrier, blacs_gridinfo, chk1mat, & 
                     compute_local_svd, descset, dlacpy, dlaset, & 
                     generate_normal_matrix, infog2l, Multiply_BAB, &
                     Multiply_BBA, pchk1mat, pdgemm, pdlarft, pdlaset, & 
                     pdmygeqr2, pdmylarfb_fc, pxerbla, slcombine, sltimer
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          max, min, random_seed
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = desca( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
  lda   = desca( LLD_ )
  nb    = desca( NB_ )

  ! 
  ! Check matrix dimensions.
  ! 
  if( m < n ) then 
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write( *, '(/,1x,a)' ) '*** ERROR in pdgeutv: m should be >= n.'
      write( *, * ) 'm:  ', m
      write( *, * ) 'n:  ', n
    end if
    info = -1
    return
  end if

  !
  ! Test the input parameters.
  !
  info = 0
  if( nprow == -1 ) then
    info = -(400+CTXT_)
  else
    call chk1mat( M, 1, N, 2, 1, 1, desca, 4, info )
    if( ( info == 0 ).and. build_u ) then
      call chk1mat( M, 1, M, 1, 1, 1, descu, 7, info )
    end if 
    if( ( info == 0 ).and. build_v ) then
      call chk1mat( N, 2, N, 2, 1, 1, descv, 10, info )
    end if 
    !
    ! Check block sizes and alignments.
    !
    if( info == 0 ) then
      if( desca( MB_ ) /= desca( NB_ ) ) then
        info = -(400+NB_)
      end if
      !
      if( ( info == 0 ) .and. build_u ) then
        if( ictxt /= descu( CTXT_ ) ) then
          info = -(700+CTXT_)
        else if( descu( MB_ ) /= descu( NB_ ) ) then
          info = -(700+NB_)
        else if( desca( MB_ ) /= descu( MB_ ) ) then
          info = -(700+MB_)
        else if( desca( rsrc_ ) /= descu( rsrc_ ) ) then
          info = -(700+rsrc_)
        else if( desca( csrc_ ) /= descu( csrc_ ) ) then
          info = -(700+csrc_)
        end if
      end if
      !
      if( ( info == 0 ) .and. build_v ) then
        if( ictxt /= descv( CTXT_ ) ) then
          info = -(1000+CTXT_)
        else if( descv( MB_ ) /= descv( NB_ ) ) then
          info = -(1000+NB_)
        else if( desca( MB_ ) /= descv( MB_ ) ) then
          info = -(1000+MB_)
        else if( desca( rsrc_ ) /= descv( rsrc_ ) ) then
          info = -(1000+rsrc_)
        else if( desca( csrc_ ) /= descv( csrc_ ) ) then
          info = -(1000+csrc_)
        end if
      end if
    end if
    !
    if( info == 0 ) then
      call pchk1mat( M, 1, N, 2, 1, 1, desca,  4, 0, idum, idum, info )
      if ( build_u ) then
        call pchk1mat( M, 1, M, 1, 1, 1, descu,  7, 0, idum, idum, info )
      end if
      if ( build_v ) then
        call pchk1mat( N, 2, N, 2, 1, 1, descv, 10, 0, idum, idum, info )
      end if
    end if
  end if
  !
  if( info /= 0 ) then
    call pxerbla( ictxt, 'pdgeutv', -info )
    return
  end if

  !
  ! Quick return if possible.
  !
  mn = min( m, n )
  if( mn == 0 ) return

  !
  ! Allocate local arrays and vectors.
  !
  allocate( vtau( n ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for vtau ***'

  allocate( svda( nb * nb ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for sa ***'

  allocate( svdu( nb * nb ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for su ***'

  allocate( svds( nb ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for vecs ***'

  allocate( svdvt( nb * nb ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for svt ***'

  allocate( s( nb * nb ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for s ***'

  allocate( work_larft( nb * nb ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for work_larft ***'

  !
  ! Some initializations.
  !
  mpa = numroc( m, nb, myrow, desca( RSRC_ ), nprow )
  nqa = numroc( n, nb, mycol, desca( CSRC_ ), npcol )

  !
  ! Set different random seeds in every process.
  !
  call random_seed( size = seed_size )
  allocate( seeds( seed_size ) )
  do k = 1, seed_size
    !!!! seeds( k ) = k
    seeds( k ) = 101 + k + myrow + mycol * nprow
  end do
  !!!! print *, myrow, mycol, ( seeds( k ), k = 1, seed_size )
  call random_seed( put = seeds )
  deallocate( seeds )

  !
  ! Set matrices U and V to the identity, if they must be built.
  !
  if( build_u ) then
    call pdlaset( 'All', m, m, ZERO, ONE, u, 1, 1, descu )
  end if
  if( build_v ) then
    call pdlaset( 'All', n, n, ZERO, ONE, v, 1, 1, descv )
  end if

  !
  ! *************************
  ! * Compute factorization *
  ! *************************
  !
  do j = 1, n, nb
    nbrow = min( nb, m - j + 1 )
    nbcol = min( nb, n - j + 1 )
    call infog2l( j, j, desca, nprow, npcol, myrow, mycol, &
                  idx_row_of_j, idx_col_of_j, proc_row_of_j, proc_col_of_j )

    if( PRINT_INTERMEDIATE_MATRICES ) then
      if( ( myrow == 0 ).AND.( mycol == 0 ) ) then
        write ( *, '(/, 40("-"))' )
        write ( *, * ) 'Iter: ', j, nb, nbrow, nbcol
        write ( *, '(40("-"), /)' )
      end if
    end if
  
    !
    ! =================================================================
    ! Rotate maximal mass of A( :, j:n ) into the current column block.
    ! =================================================================
    !
    ! Perform this processing only if there are more columns to the right of
    ! the current column block.
    !
    if( ( n - j - nbcol + 1 ) > 0 ) then
      !
      ! Create matrices G and Y.
      ! ========================
      !
      mpg = numroc( m - j + 1,  nb, myrow, proc_row_of_j, nprow )
      nqg = numroc( nb, nb, mycol, proc_col_of_j, npcol )
      mpy = numroc( n - j + 1,  nb, myrow, proc_row_of_j, nprow )
      nqy = numroc( nb, nb, mycol, proc_col_of_j, npcol )
  
      call descset( descg, m - j + 1, nb, nb, nb, &
                    proc_row_of_j, proc_col_of_j, ictxt, max( 1, mpg ) )
      call descset( descy, n - j + 1, nb, nb, nb, &
                    proc_row_of_j, proc_col_of_j, ictxt, max( 1, mpy ) )
  
      allocate( g( max( 1, mpg * nqg ) ), stat = allocstat )
      if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for g ***'
  
      allocate( y( max( 1, mpy * nqy ) ), stat = allocstat )
      if ( allocstat /= 0 ) stop 'pdgeutv: *** Not enough memory for y ***'
  
      ! 
      ! Generate normal random matrix G.
      ! ================================
      !
      if( PROFILE ) then
        call blacs_barrier( ictxt, 'All' )
        call sltimer( 2 )
      end if

      !!!! call pdlaset( 'All', m - j + 1, nbcol, ZERO, ONE, g, 1, 1, descg )
      call generate_normal_matrix( m - j + 1, nbcol, g, 1, 1, descg )

      if( PROFILE ) then
        call sltimer( 2 )
      end if

      if( PRINT_INTERMEDIATE_MATRICES ) then
        call print_distributed_matrix( m - j + 1, nbcol, g, 1, 1, descg, 'g0' )
      end if

      !         
      ! Compute the sampling matrix Y.
      ! ==============================
      !         
      ! Y = Aloc' * G.
      !
      if( PROFILE ) then
        call blacs_barrier( ictxt, 'All' )
        call sltimer( 3 )
      end if

      call pdgemm( 'Transpose', 'No transpose', n - j + 1, nbcol, m - j + 1, &
                   ONE, a, j, j, desca, g, 1, 1, descg, &
                   ZERO, y, 1, 1, descy )

      if( PRINT_INTERMEDIATE_MATRICES ) then
        call print_distributed_matrix( n - j + 1, nbcol, y, 1, 1, descy, 'y0' )
      end if

      !
      ! Perform the 'power iteration' process.
      ! ======================================
      !
      ! for i_iter = 1:n_iter
      !   Y = Aloc' * ( Aloc * Y );
      ! end
      !
      do k = 1, q
        ! Reuse matrix G.
        call pdgemm( 'No transpose', 'No transpose', &
                     m - j + 1, nbcol, n - j + 1, &
                     ONE, a, j, j, desca, y, 1, 1, descy, &
                     ZERO, g, 1, 1, descg )
        call pdgemm( 'Transpose', 'No transpose', &
                     n - j + 1, nbcol, m - j + 1, &
                     ONE, a, j, j, desca, g, 1, 1, descg, &
                     ZERO, y, 1, 1, descy )
      end do

      if( PROFILE ) then
        call sltimer( 3 )
      end if

      if( PRINT_INTERMEDIATE_MATRICES ) then
        call print_distributed_matrix( n - j + 1, nbcol, y, 1, 1, descy, 'y1' )
      end if

      !         
      ! Update A from the right side. Update V if asked.
      ! ================================================
      !
      ! Construct the Householder transformations to be applied from the right.
      !
      if( PROFILE ) then
        call blacs_barrier( ictxt, 'All' )
        call sltimer( 4 )
      end if

      call pdmygeqr2( n - j + 1, nbcol, y, 1, 1, descy, vtau )

      if( PROFILE ) then
        call sltimer( 4 )
      end if

      if( PRINT_INTERMEDIATE_MATRICES ) then
        call print_distributed_matrix( n - j + 1, nbcol, y, 1, 1, descy, 'y2' )
      end if

      ! 
      ! Apply the Householder transformations to rotate maximal mass into the 
      ! current column block.
      ! T(:,[J2,J3])  = T(:,[J2,J3])*Vloc;
      !
      if( PROFILE ) then
        call blacs_barrier( ictxt, 'All' )
        call sltimer( 5 )
      end if

      !
      ! Form the triangular factor of the block reflector.
      !
      call pdlarft( 'Forward', 'Columnwise', n - j + 1, nbcol, &
                    y, 1, 1, descy, vtau, s, work_larft )
      !
      ! Apply H to matrix A from the right.
      !
      call pdmylarfb_fc( 'Right', 'No transpose', m, n - j + 1, nbcol, & 
                         y, 1, 1, descy, s, a, 1, j, desca )

      !
      ! Apply H to matrix V from the right.
      !
      if( build_v ) then
        call pdmylarfb_fc( 'Right', 'No transpose', n, n - j + 1, nbcol, & 
                           y, 1, 1, descy, s, v, 1, j, descv )
      end if

      ! Deallocate matrices G and Y.
      deallocate( g )
      deallocate( y )

      if( PROFILE ) then
        call sltimer( 5 )
      end if
    end if

    if( PRINT_INTERMEDIATE_MATRICES ) then
      call print_distributed_matrix( m, n, a, 1, 1, desca, 'a1' )
    end if

    !
    ! ========================================================
    ! Annihilate elements below the diagonal in current block.
    ! ========================================================
    ! 
    ! Perform this processing only if there are more rows below the diagonal 
    ! block.
    !
    if( ( m - j - nbrow + 1 ) > 0 ) then
      !
      ! Update A from the left side. Update U if asked.
      ! ===============================================
      !
      ! Determine the rotations to be applied from the left.
      ! [Uloc,Dloc] = LOCAL_nonpiv_QR(T([J2,I3],J2));
      !
      if( PROFILE ) then
        call blacs_barrier( ictxt, 'All' )
        call sltimer( 6 )
      end if

      call pdmygeqr2( m - j + 1, nbcol, a, j, j, desca, vtau )

      if( PROFILE ) then
        call sltimer( 6 )
      end if

      if( PRINT_INTERMEDIATE_MATRICES ) then
        call print_distributed_matrix( m, n, a, 1, 1, desca, 'a2' )
      end if

      ! 
      ! Update the rest of matrix A with transformations from the second QR.
      !
      if( PROFILE ) then
        call blacs_barrier( ictxt, 'All' )
        call sltimer( 7 )
      end if

      !
      ! Form the triangular factor of the block reflector.
      !
      call pdlarft( 'Forward', 'Columnwise', m - j + 1, nbcol, & 
                    a, j, j, desca, vtau, s, work_larft )

      !
      ! Apply H' to matrix A from the left.
      !
      call pdmylarfb_fc( 'Left', 'Transpose', &
                         m - j + 1, n - j - nbcol + 1, nbcol, &
                         a, j, j, desca, s, a, j, j + nbcol, desca )

      !
      ! Apply H to matrix U from the right.
      !
      if( build_u ) then
        call pdmylarfb_fc( 'Right', 'No transpose', m, m - j + 1, nbcol, & 
                           a, j, j, desca, s, u, 1, j, descu )
      end if

      if( PROFILE ) then
        call sltimer( 7 )
      end if

      if( PRINT_INTERMEDIATE_MATRICES ) then
        call print_distributed_matrix( m, n, a, 1, 1, desca, 'a3' )
      end if

      !
      ! Set to zero elements below the diagonal in current block.
      !
      call pdlaset( 'Lower', m - j, nbcol, ZERO, ZERO, a, j + 1, j, desca )

    end if

    if( PRINT_INTERMEDIATE_MATRICES ) then
      call print_distributed_matrix( m, n, a, 1, 1, desca, 'a4' )
    end if

    ! 
    ! ===============================================================
    ! Compute SVD of diagonal block, and update matrices A, U, and V.
    ! ===============================================================
    ! 
    if( PROFILE ) then
      call blacs_barrier( ictxt, 'All' )
      call sltimer( 8 )
    end if

    !
    ! The SVD of the diagonal block is computed locally by the process owner 
    ! of the current diagonal block of A.
    !
    if( ( myrow == proc_row_of_j ).AND.( mycol == proc_col_of_j ) ) then

      ! Copy current diagonal block of A into local matrix SA.
      offjj = idx_row_of_j + ( idx_col_of_j - 1 ) * lda
      call dlacpy( 'All', nbrow, nbcol, a( offjj ), lda, svda, nbrow )

      ! Set current diagonal block to zero.
      call dlaset( 'All', nbrow, nbcol, ZERO, ZERO, a( offjj ), lda )

      ! Compute svd of local matrix SA.
      !!!! call print_local_matrix( nbrow, nbcol, sa, nbrow, 'sai' )
      call compute_local_svd( nbrow, nbcol, svda, nbrow, svdu, nbrow, svds, &
                              svdvt, nbcol, info )
      !!!! call print_local_matrix( nbrow, nbcol, sa, nbrow, 'saf' )
      !!!! call print_local_matrix( nbrow, nbrow, su, nbrow, 'suf' )
      !!!! call print_local_matrix( nbcol, nbcol, svt, nbcol, 'svtf' )
      !!!! call print_local_matrix( min( nbrow, nbcol ), 1, vecs, nbrow, 'sf' )

      ! Copy back the singular values into the diagonal of the diagonal block.
      do ii = 1, min( nbrow, nbcol )
        a( offjj - 1 + ii + ( ii - 1 ) * lda ) = svds( ii )
      end do
    end if

    if( PROFILE ) then
      call sltimer( 8 )
    end if

    if( PROFILE ) then
      call blacs_barrier( ictxt, 'All' )
      call sltimer( 9 )
    end if

    call descset( descsvdu, nbrow, nbrow, nb, nb, &
                  proc_row_of_j, proc_col_of_j, ictxt, max( 1, nbrow ) )

    if( PRINT_INTERMEDIATE_MATRICES ) then
      call print_distributed_matrix( m, n, a, 1, 1, desca, 'a5' )
    end if

    !
    ! Apply U of miniSVD to A.
    !
    if( ( j + nbcol - 1 ) < n ) then
      call Multiply_BAB( 'Transpose', 'No transpose', & 
                         nbrow, n - j - nbcol + 1, &
                         svdu, descsvdu, a, j, j + nbcol, desca )
    end if

    if( PRINT_INTERMEDIATE_MATRICES ) then
      call print_distributed_matrix( m, n, a, 1, 1, desca, 'a6' )
    end if

    !
    ! Apply U of miniSVD to U.
    !
    if( build_u ) then
      call Multiply_BBA( 'No transpose', 'No transpose', m, nbrow, &
                         svdu, descsvdu, u, 1, j, descu )
    end if

    if( PROFILE ) then
      call sltimer( 9 )
    end if

    if( PROFILE ) then
      call blacs_barrier( ictxt, 'All' )
      call sltimer( 10 )
    end if

    call descset( descsvdvt, nbcol, nbcol, nb, nb, &
                  proc_row_of_j, proc_col_of_j, ictxt, max( 1, nbcol ) )

    if( PRINT_INTERMEDIATE_MATRICES ) then
      call print_distributed_matrix( m, n, a, 1, 1, desca, 'a7' )
    end if

    !
    ! Apply V of miniSVD to A.
    !
    if( j > 1 ) then
      call Multiply_BBA( 'No transpose', 'Transpose', j - 1, nbcol, &
                         svdvt, descsvdvt, a, 1, j, desca )
    end if

    if( PRINT_INTERMEDIATE_MATRICES ) then
       call print_distributed_matrix( m, n, a, 1, 1, desca, 'a8' )
    end if

    !
    ! Apply V of miniSVD to V.
    !
    if( build_v ) then
      call Multiply_BBA( 'No transpose', 'Transpose', n, nbcol, &
                         svdvt, descsvdvt, v, 1, j, descv )
    end if

    if( PROFILE ) then
      call sltimer( 10 )
    end if

    if( PRINT_INTERMEDIATE_MATRICES ) then
      call print_distributed_matrix( m, n, a, 1, 1, desca, 'a9' )
    end if
  end do

  ! Deallocate local arrays and vectors.
  deallocate( vtau )
  deallocate( svda )
  deallocate( svdu )
  deallocate( svds )
  deallocate( svdvt )
  deallocate( s )
  deallocate( work_larft )
 
  ! Compute and print final timings, if needed.
  if( PROFILE ) then
    ! Gather maximum of all CPU and WALL clock timings.
    call slcombine( ictxt, 'All', '>', 'w', 1, 2, t_normal )
    call slcombine( ictxt, 'All', '>', 'w', 1, 3, t_power )
    call slcombine( ictxt, 'All', '>', 'w', 1, 4, t_rgeqr )
    call slcombine( ictxt, 'All', '>', 'w', 1, 5, t_rormq )
    call slcombine( ictxt, 'All', '>', 'w', 1, 6, t_lgeqr )
    call slcombine( ictxt, 'All', '>', 'w', 1, 7, t_lormq )
    call slcombine( ictxt, 'All', '>', 'w', 1, 8, t_svd )
    call slcombine( ictxt, 'All', '>', 'w', 1, 9, t_mma12 )
    call slcombine( ictxt, 'All', '>', 'w', 1, 10, t_mma01 )
  
    if( ( myrow == 0 ).and.( mycol == 0 ) ) then
      write ( *, '(/,1x,a)' ) 'Profiling of pdgeutv (time in s.): '
      write ( *, '(1x,a,f12.4)' ) '  t_normal: ', t_normal
      write ( *, '(1x,a,f12.4)' ) '  t_power:  ', t_power
      write ( *, '(1x,a,f12.4)' ) '  t_rgeqr:  ', t_rgeqr
      write ( *, '(1x,a,f12.4)' ) '  t_rormq:  ', t_rormq
      write ( *, '(1x,a,f12.4)' ) '  t_lgeqr:  ', t_lgeqr
      write ( *, '(1x,a,f12.4)' ) '  t_lormq:  ', t_lormq
      write ( *, '(1x,a,f12.4)' ) '  t_svd:    ', t_svd
      write ( *, '(1x,a,f12.4)' ) '  t_mma12:  ', t_mma12
      write ( *, '(1x,a,f12.4)' ) '  t_mma01:  ', t_mma01
      write ( *, '(1x,a,f12.4)' ) '  t_total:  ', & 
          t_normal + t_power + t_rgeqr + t_rormq + t_lgeqr + t_lormq + & 
          t_svd + t_mma12 + t_mma01
      write ( *, '(1x,a)' ) 'End of profiling'
    end if
  end if

  return
  ! *** Last line of pdgeutv ***
end

! =============================================================================
subroutine pdmygeqr2( m, n, a, ia, ja, desca, tau )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            ia, ja, m, n
  ! ..
  ! .. Array Arguments ..
  integer            desca( * )
  double precision   a( * ), tau( * )
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It computes the QR factorization of a real distributed m-by-n matrix
  !  sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = Q * R.
  !
  !  =====================================================================
  !
  ! .. Local Scalars ..
  double precision   scalar_work
  integer            info, len_work, allocstat
  ! ..
  ! .. Local Arrays ..
  double precision, dimension ( : ), allocatable :: work
  ! ..
  ! .. External Subroutines ..
  external           pdgeqr2
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          int
  ! ..
  ! .. Executable Statements ..

  !
  ! Compute workspace length.
  !
  call pdgeqr2( m, n, a, ia, ja, desca, tau, scalar_work, -1, info )
  if ( info /= 0 ) then
    write ( *, * ) &
        '*** ERROR in pdmygeqr2: Info of call to compute wk length: ', info
    stop
  end if
  len_work = int( scalar_work )

  !
  ! Allocate workspace.
  !
  allocate( work( len_work ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdmygeqr2: *** Not enough memory for work ***'

  !
  ! Call to pdgeqr2.
  !
  call pdgeqr2( m, n, a, ia, ja, desca, tau, work, len_work, info )
  if ( info /= 0 ) then
    write ( *, * ) '*** ERROR in pdmygeqr2: Info of pdgeqr2: ', info
    stop
  end if

  !
  ! Deallocate workspace.
  !
  deallocate( work )

  return
  !
  ! End of pdmygeqr2
  !
end

! =============================================================================
subroutine pdmylarfb_fc( side, trans, m, n, k, v, iv, jv, descv, t, & 
                         c, ic, jc, descc )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  character          side, trans
  integer            ic, iv, jc, jv, k, m, n
  ! ..
  ! .. Array Arguments ..
  integer            descc( * ), descv( * )
  double precision   c( * ), t( * ), v( * )
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It applies a real block reflector Q or its transpose Q**T to a real 
  !  distributed M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) from the 
  !  left or from the right.
  !  Matrix Q is intrisically stored in matrix v (Householder vectors) and 
  !  vector t (tau factors).
  !
  !  Workspace required inside this function
  !  =======================================
  !
  !  WORK    (local workspace) double precision array, dimension (LWORK)
  !          If STOREV = 'C',
  !            if SIDE = 'L',
  !              LWORK >= ( NqC0 + MpC0 ) * K
  !            else if SIDE = 'R',
  !              LWORK >= ( NqC0 + MAX( NpV0 + NUMROC( NUMROC( N+ICOFFC,
  !                         NB_V, 0, 0, NPCOL ), NB_V, 0, 0, LCMQ ),
  !                         MpC0 ) ) * K
  !            end if
  !          else if STOREV = 'R',
  !            if SIDE = 'L',
  !              LWORK >= ( MpC0 + MAX( MqV0 + NUMROC( NUMROC( M+IROFFC,
  !                         MB_V, 0, 0, NPROW ), MB_V, 0, 0, LCMP ),
  !                         NqC0 ) ) * K
  !            else if SIDE = 'R',
  !              LWORK >= ( MpC0 + NqC0 ) * K
  !            end if
  !          end if
  !
  !          where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ),
  !
  !          IROFFV = MOD( IV-1, MB_V ), ICOFFV = MOD( JV-1, NB_V ),
  !          IVROW = INDXG2P( IV, MB_V, MYROW, RSRC_V, NPROW ),
  !          IVCOL = INDXG2P( JV, NB_V, MYCOL, CSRC_V, NPCOL ),
  !          MqV0 = NUMROC( M+ICOFFV, NB_V, MYCOL, IVCOL, NPCOL ),
  !          NpV0 = NUMROC( N+IROFFV, MB_V, MYROW, IVROW, NPROW ),
  !
  !          IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
  !          ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
  !          ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
  !          MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
  !          NpC0 = NUMROC( N+ICOFFC, MB_C, MYROW, ICROW, NPROW ),
  !          NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
  !
  !          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
  !          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
  !          the subroutine BLACS_GRIDINFO.
  !
  !  =====================================================================
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
  double precision   scalar_work
  integer            ictxt, myrow, mycol, nprow, npcol, &
                     INFO, LWORK, LCM, LCMQ, & 
                     IROFFV, ICOFFV, IVROW, IVCOL, MqV0, NpV0, & 
                     IROFFC, ICOFFC, ICROW, ICCOL, MpC0, NpC0, NqC0, &
                     RSRC_C, CSRC_C, RSRC_V, CSRC_V, &
                     MB_C, NB_C, MB_V, NB_V, &
                     len_work, allocstat
  ! ..
  ! .. Local Arrays ..
  double precision, dimension ( : ), allocatable :: work
  ! ..
  ! .. External Functions ..
  logical            lsame
  integer            ilcm, indxg2p, numroc
  external           ilcm, indxg2p, numroc, lsame
  ! ..
  ! .. External subroutines ..
  external           blacs_gridinfo, pdlarfb, pdormqr
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          max, min, mod
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = descc( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

  !
  ! Compute workspace for pdlarfb with arguments 'Forward' and 'Columnwise'.
  ! 
  RSRC_C = DESCC( RSRC_ )
  CSRC_C = DESCC( CSRC_ )
  RSRC_V = DESCV( RSRC_ )
  CSRC_V = DESCV( CSRC_ )
  MB_C   = DESCC( MB_ )
  NB_C   = DESCC( NB_ )
  MB_V   = DESCV( MB_ )
  NB_V   = DESCV( NB_ )
  ! 
  LCM = ILCM( NPROW, NPCOL )
  LCMQ = LCM / NPCOL
  !
  IROFFV = MOD( IV-1, MB_V )
  ICOFFV = MOD( JV-1, NB_V )
  IVROW = INDXG2P( IV, MB_V, MYROW, RSRC_V, NPROW )
  IVCOL = INDXG2P( JV, NB_V, MYCOL, CSRC_V, NPCOL )
  MqV0 = NUMROC( M+ICOFFV, NB_V, MYCOL, IVCOL, NPCOL )
  NpV0 = NUMROC( N+IROFFV, MB_V, MYROW, IVROW, NPROW )
  !
  IROFFC = MOD( IC-1, MB_C ) 
  ICOFFC = MOD( JC-1, NB_C )
  ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW )
  ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL )
  MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW )
  NpC0 = NUMROC( N+ICOFFC, MB_C, MYROW, ICROW, NPROW )
  NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL )
  !
  if( lsame( SIDE, 'L' ) ) then
    LWORK = ( NqC0 + MpC0 ) * K
  elseif( lsame( SIDE, 'R' ) ) then
    LWORK = ( NqC0 + MAX( NpV0 + NUMROC( NUMROC( N+ICOFFC, &
              NB_V, 0, 0, NPCOL ), NB_V, 0, 0, LCMQ ), &
              MpC0 ) ) * K
  else
    LWORK = 0
  end if

  !
  ! Compute workspace required by pdormqr, and compare with the one of pdlarfb.
  ! The workspace length of pdormqr should be larger than that of pdlarfb.
  !
  call pdormqr( SIDE, TRANS, M, N, K, &
                V, IV, JV, DESCV, T, C, IC, JC, DESCC, &
                scalar_work, -1, INFO )
  len_work = int( scalar_work )
  if ( len_work < LWORK ) stop '*** ERROR in pdmylarfb_fc: len_work < LWORK.'

  !
  ! Allocate workspace.
  !
  len_work = LWORK
  allocate( work( len_work ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'pdmylarfb_fc: *** Not enough memory for work ***'

  !
  ! Call to pdlarfb.
  !
  call pdlarfb( SIDE, TRANS, 'Forward', 'Columnwise', M, N, K, V, IV, &
                JV, DESCV, T, C, IC, JC, DESCC, work )

  !
  ! Deallocate workspace.
  !
  deallocate( work )

  return
  !
  ! End of pdmylarfb_fc
  !
end

! =============================================================================
subroutine compute_local_svd( m, n, a, lda, u, ldu, s, vt, ldvt, info )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            m, n, lda, ldu, ldvt, info
  ! ..
  ! .. Array Arguments ..
  double precision   a( lda, * ), u( ldu, * ), s( * ), vt( ldvt, * )
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It computes the SVD factorization of local matrix A, and returns 
  !  orthonormal matrices U and VT, and the singular values in s.
  !  This operation is performed on local matrices and vectors.
  !
  !  =====================================================================
  !
  ! .. Local Scalars ..
  integer            len_work, allocstat
  double precision   scalar_work
  ! ..
  ! .. Local Arrays ..
  double precision, dimension ( : ), allocatable :: work
  ! ..
  ! .. External subroutines ..
  external           dgesvd
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic          int
  ! ..
  ! .. Executable Statements ..

  !
  ! Obtain the optimal real workspace length.
  !
  len_work = -1
  call dgesvd( 'All', 'All', m, n, a, lda, s, u, ldu, vt, ldvt, &
               scalar_work, len_work, info )
  if( info /= 0 ) then
    write ( *, * ) '*** ERROR in compute_local_svd: ', & 
        'Info of call to compute wk length: ', info
  end if
  len_work = int( scalar_work )

  !
  ! Allocate workspace.
  !
  allocate( work( len_work ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'compute_local_svd: *** Not enough memory for work ***'

  !
  ! Compute singular values and vectors.
  !
  call dgesvd( 'All', 'All', m, n, a, lda, s, u, ldu, vt, ldvt, &
               work, len_work, info )
  if( info /= 0 ) then
    write ( *, * ) '*** ERROR in compute_local_svd: Info of dgesvd: ', info
  end if

  !
  ! Deallocate workspace.
  !
  deallocate( work )

  !
  ! End of compute_local_svd
  !
  return
end

! =============================================================================
subroutine Multiply_BBA( transa, transb, m, n, a, desca, b, ib, jb, descb )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  character*( * )    transa, transb
  integer            m, n, ib, jb
  ! ..
  ! .. Array Arguments ..
  integer            desca( * ), descb( * )
  double precision   a( * ), b( * )
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It computes the following operation:  B := B * A, where A and B are 
  !  transposed according to arguments "transa" and "transb".
  !  Matrix B is m x n.
  !
  !  =====================================================================
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
  integer            ictxt, nprow, npcol, myrow, mycol, nb, & 
                     proc_row_of_ib_b, proc_col_of_jb_b, &
                     mpbc, nqbc, len_bc, allocstat
  ! ..
  ! .. Local Arrays ..
  integer            descbc( dlen_ )
  double precision, dimension ( : ), allocatable :: bc
  ! ..
  ! .. External Functions ..
  integer            indxg2p, numroc
  external           indxg2p, numroc
  ! ..
  ! .. External Subroutines ..
  external           blacs_gridinfo, descset, pdgemm, pdlacpy
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = descb( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
  nb    = descb( NB_ )

  ! Prepare descriptor for copy of matrix b.
  proc_row_of_ib_b = indxg2p( ib, descb( MB_ ), myrow, descb( RSRC_ ), nprow )
  proc_col_of_jb_b = indxg2p( jb, descb( NB_ ), mycol, descb( CSRC_ ), npcol )
  mpbc = numroc( m, nb, myrow, proc_row_of_ib_b, nprow )
  nqbc = numroc( n, nb, mycol, proc_col_of_jb_b, npcol )
  call descset( descbc, m, n, nb, nb, &
                proc_row_of_ib_b, proc_col_of_jb_b, ictxt, max( 1, mpbc ) )

  ! Allocate copy of matrix b.
  len_bc = max( 1, mpbc * nqbc )
  allocate( bc( len_bc ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'Multiply_BBA: *** Not enough memory for bc ***'

  ! Copy b into bc.
  call pdlacpy( 'All', m, n, b, ib, jb, descb, &
                             bc, 1, 1, descbc )

  ! Compute:  b := bc * a.
  call pdgemm( transa, transb, m, n, n, &
               ONE, bc, 1, 1, descbc, a, 1, 1, desca, &
               ZERO, b, ib, jb, descb )

  ! Deallocate copy of matrix B.
  deallocate( bc )

  return
  !
  ! End of Multiply_BBA
  !
end

! =============================================================================
subroutine Multiply_BAB( transa, transb, m, n, a, desca, b, ib, jb, descb )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  character*( * )    transa, transb
  integer            m, n, ib, jb
  ! ..
  ! .. Array Arguments ..
  integer            desca( * ), descb( * )
  double precision   a( * ), b( * )
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It computes the following operation:  B := A * B, where A and B are 
  !  transposed according to arguments "transa" and "transb".
  !  Matrix B is m x n.
  !
  !  =====================================================================
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
  integer            ictxt, nprow, npcol, myrow, mycol, nb, & 
                     proc_row_of_ib_b, proc_col_of_jb_b, &
                     mpbc, nqbc, len_bc, allocstat
  ! ..
  ! .. Local Arrays ..
  integer            descbc( DLEN_ )
  double precision, dimension ( : ), allocatable :: bc
  ! ..
  ! .. External Functions ..
  integer            indxg2p, numroc
  external           indxg2p, numroc
  ! ..
  ! .. External Subroutines ..
  external           blacs_gridinfo, descset, pdgemm, pdlacpy
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = descb( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
  nb    = descb( NB_ )

  !
  ! Prepare descriptor for copy of matrix b.
  !
  proc_row_of_ib_b = indxg2p( ib, descb( MB_ ), myrow, descb( RSRC_ ), nprow )
  proc_col_of_jb_b = indxg2p( jb, descb( NB_ ), mycol, descb( CSRC_ ), npcol )
  mpbc = numroc( m, nb, myrow, proc_row_of_ib_b, nprow )
  nqbc = numroc( n, nb, mycol, proc_col_of_jb_b, npcol )
  call descset( descbc, m, n, nb, nb, &
                proc_row_of_ib_b, proc_col_of_jb_b, ictxt, max( 1, mpbc ) )

  !
  ! Allocate copy of matrix b.
  !
  len_bc = max( 1, mpbc * nqbc )
  allocate( bc( len_bc ), stat = allocstat )
  if ( allocstat /= 0 ) stop 'Multiply_BAB: *** Not enough memory for bc ***'

  !
  ! Copy b into bc.
  !
  call pdlacpy( 'All', m, n, b, ib, jb, descb, &
                             bc, 1, 1, descbc )

  !
  ! Compute:  b := a * bc.
  !
  call pdgemm( transa, transb, m, n, m, &
               ONE, a, 1, 1, desca, bc, 1, 1, descbc, &
               ZERO, b, ib, jb, descb )

  !
  ! Deallocate copy of matrix B.
  !
  deallocate( bc )

  return
  !
  ! End of Multiply_BAB
  !
end

! =============================================================================
subroutine generate_normal_matrix( m, n, a, ia, ja, desca )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  integer            m, n, ia, ja
  ! ..
  ! .. Array Arguments ..
  integer            desca( * )
  double precision   a( * )
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It generates a distributed random matrix with elements in a normal
  !  distribuition.
  !
  !  =====================================================================
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
  integer            ictxt, nprow, npcol, myrow, mycol, lda, nb, &
                     i, j, mp, nq, ip, jq, iarow, jacol
  ! ..
  ! .. External Functions ..
  integer            numroc
  double precision   generate_normal_random_number
  external           numroc, generate_normal_random_number
  ! ..
  ! .. External subroutines ..
  external           blacs_gridinfo, infog2l
  ! ..
  ! .. Executable Statements ..

  !
  ! Get grid parameters.
  !
  ictxt = desca( CTXT_ )
  call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
  lda   = desca( LLD_ )
  nb    = desca( NB_ )

  !
  ! Compute local indices.
  !
  call infog2l( ia, ja, desca, nprow, npcol, myrow, mycol, &
                ip, jq, iarow, jacol )
  mp = numroc( ia + m - 1, nb, myrow, desca( RSRC_ ), nprow )
  nq = numroc( ja + n - 1, nb, mycol, desca( CSRC_ ), npcol )

  !
  ! Process elements in local matrix.
  !

  do j = jq, nq
    do i = ip, mp
      !!!! num = generate_normal_random_number( ZERO, ONE )
      !!!! print *, myrow, mycol, i, j, num
      !!!! a( i + ( j - 1 ) * lda ) = num
      a( i + ( j - 1 ) * lda ) = generate_normal_random_number( ZERO, ONE )
    end do
  end do

  !
  ! End of generate_normal_matrix
  !
  return
end

! =============================================================================
double precision function generate_normal_random_number( mu, sigma )
  !
  implicit none
  !
  ! .. Scalar Arguments ..
  double precision   mu, sigma
  ! ..
  !
  !  Purpose
  !  =======
  !
  !  It returns a double-precision random number with a normal distribution.
  !
  !  =====================================================================
  !
  ! .. Parameters ..
  double precision   ONE, ZERO
  parameter          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  ! ..
  ! .. Local Scalars ..
  logical, save ::          alternate = .FALSE.
  double precision, save :: b1, b2
  double precision ::       u, c1, c2, a, factor
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic                 log, sqrt, random_number
  ! ..
  ! .. Executable Statements ..

  !
  ! Quick return.
  !
  if( alternate ) then
    alternate = .not. alternate
    generate_normal_random_number = mu + sigma * b2
    return
  end if

  !
  ! Main loop.
  !
  do
    call random_number( u )
    c1 = -1.0 + 2.0 * u
    call random_number( u )
    c2 = -1.0 + 2.0 * u
    a  = c1 * c1 + c2 * c2;
    if( .not.( ( a == ZERO ).or.( a >= ONE ) ) ) exit
  end do
  factor    = sqrt( ( -2.0 * log( a ) ) / a );
  b1        = c1 * factor;
  b2        = c2 * factor;
  alternate = .not. alternate
  generate_normal_random_number = mu + sigma * b1

  !
  ! End of generate_normal_random_number
  !
  return
end

