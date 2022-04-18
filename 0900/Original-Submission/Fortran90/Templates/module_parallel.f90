!*********************************************************************************************************************************
!*
!* MODULE: PARALLEL
!*
!* PURPOSE: provides procedures for parallelization
!*
!* DEPENDENCIES: 
!*               - PRECISION
!*               - MPI
!*               - MODEL
!*               - OBSERVATIONS
!* 
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
!*          
!*          - VARIABLES DEFINED BY MPIF.H FILE (from MPI)
!*          
!*          - MODESMODEL (integer): number of model modes (from MODEL)
!*          - MODELERROR (real array): covariance matrix of model errors (from MODEL)
!*          - SQRTMODELERROR (real array): square root of covariance matrix of model errors (from MODEL)
!*
!*          - MODESOBS (integer): number of observation modes (from OBSERVATIONS)
!*          - NUMBEROBS (integer): number of observations (from OBSERVATIONS)
!*          - IFOBS (logical): flag to determine if there are observations (from OBSERVATIONS)
!*          - OBSVALUE (real array): observation values (from OBSERVATIONS)
!*          - COVOBS (real array): covariance matrix of observation errors (from OBSERVATIONS)
!*          - SQRTCOVOBS (real array): square root of the covariance matrix of observation errors (from OBSERVATIONS)
!*
!*          - IERROR (integer): variable to detect errors
!*          - RANK (integer): which machine
!*          - NPROC (integer): number of processors
!*          - BLOCK_SIZE (integer): block size
!*          - NPROW (integer): number of rows in the processor grid
!*          - NPCOL (integer): number of columns in the processor grid
!*          - CONTEXT (integer): it is a universe where messages exist and do not interact with other context's messages
!*          - MYROW (integer): calling processor's row number in the processor grid 
!*          - MYCOL (integer): calling processor's column number in the processor grid 
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module parallel
  !*
  use precision
  use model
  use observations
  !*
  implicit none
  !*
  include "mpif.h"
  !*
  integer( kind = pin ) , public :: ierror !* variable to detect errors
  integer( kind = pin ) , public :: rank !* which machine
  integer( kind = pin ) , public :: nproc !* number of processors
  integer( kind = pin ) , public :: block_size !* block size
  integer( kind = pin ) , public :: nprow !* number of rows in the processor grid
  integer( kind = pin ) , public :: npcol !* number of columns in the processor grid
  integer( kind = pin ) , public :: context !* it is a universe where messages exist and do not interact with other context's mes&
                                            !* &sages
  integer( kind = pin ) , public :: myrow !* calling processor's row number in the processor grid 
  integer( kind = pin ) , public :: mycol !* calling processor's column number in the processor grid 
  !*
contains
  !*
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_INIT
  !*
  !* PURPOSE: initialization of MPI
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - MPI_INIT (from MPI)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the MPI_INIT subroutine
  !*
  !*******************************************************************************************************************************
  subroutine parallel_init()
    implicit none
    call mpi_init( ierror )
  end subroutine parallel_init
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_RANKSIZE
  !*
  !* PURPOSE: get rank and size 
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS: none
  !*
  !* CALLS: 
  !*        - MPI_COMM_RANK (from MPI)
  !*        - MPI_COMM_SIZE (from MPI)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the MPI_COMM_RANK and MPI_COMM_SIZE subroutines
  !*
  !*******************************************************************************************************************************
  subroutine parallel_ranksize()
    implicit none
    call mpi_comm_rank( mpi_comm_world , rank , ierror )
    call mpi_comm_size( mpi_comm_world , nproc , ierror )
  end subroutine parallel_ranksize
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_FINALIZE
  !*
  !* PURPOSE: finalization of MPI
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - MPI_FINALIZE (from MPI)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the MPI_FINALIZE subroutine
  !*
  !*******************************************************************************************************************************
  subroutine parallel_finalize()
    implicit none
    call mpi_finalize( ierror )
  end subroutine parallel_finalize
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_PREGRID
  !*
  !* PURPOSE: preparation of the grid
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - SL_INIT (from SCALAPACK)
  !*        - BLACS_GRIDINFO (from BLACS)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the SL_INIT and BLACS_GRIDINFO subroutines
  !*
  !*******************************************************************************************************************************
  subroutine parallel_pregrid()
    implicit none
    block_size = 2_pin
    nprow = int( sqrt( nproc * 1.0_pre ) )
    npcol = nproc / nprow
    call sl_init( context , nprow , npcol )
    call blacs_gridinfo( context , nprow , npcol , myrow , mycol )
  end subroutine parallel_pregrid
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_POSTGRID
  !*
  !* PURPOSE: finalization of the grid
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - BLACS_GRIDEXIT (from BLACS)  
  !*
  !* COMMENTS: this subroutine is a wrapper of the BLACS_GRIDEXIT subroutine
  !*
  !*******************************************************************************************************************************
  subroutine parallel_postgrid()
    implicit none
    call blacs_gridexit( context )
  end subroutine parallel_postgrid
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_PREDISTRIBUTE
  !*
  !* PURPOSE: distribution a global matrix and getting descriptors
  !*
  !* INPUTS: 
  !*         - M_GLOBAL (integer): number of rows of the global matrix
  !*         - N_GLOBAL (integer): number of columns of the global matrix
  !*         - A_GLOBAL (real array of size M_GLOBAL x N_GLOBAL): global matrix
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*         - A_DESC (integer array of size 9): BLACS context handle identifying the created process grid
  !*         - M_LOCAL (integer): number of rows of the local matrix
  !*         - N_LOCAL (integer): number of columns of the local matrix
  !*         - A_LOCAL (pointer to a real array of dimension M_LOCAL x N_LOCAL): local matrix
  !*
  !* CALLS:
  !*        - DESCINIT (from SCALAPACK)
  !*        - BLACS_BARRIER (from BLACS)
  !*  
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine parallel_predistribute( m_global , n_global , a_global , a_desc , m_local , n_local , a_local )
    implicit none
    integer( kind = pin ) , intent( in ) :: m_global
    integer( kind = pin ) , intent( in ) :: n_global
    real( kind = pre ) , dimension( m_global , n_global ) , intent( in ) :: a_global
    integer( kind = pin ) , dimension( 9_pin ) , intent( out ) :: a_desc
    integer( kind = pin ) , intent( out ) :: m_local
    integer( kind = pin ) , intent( out ) :: n_local
    real( kind = pre ) , dimension( : , : ) , pointer :: a_local
    integer( kind = pin ) :: info
    integer( kind = pin ) :: i_local
    integer( kind = pin ) :: j_local
    integer( kind = pin ) :: i_global
    integer( kind = pin ) :: j_global
    integer( kind = pin ) :: numroc
    integer( kind = pin ) :: indxl2g
    m_local = numroc( m_global , block_size , myrow , 0_pin , nprow )
    n_local = numroc( n_global , block_size , mycol , 0_pin , npcol )
    call descinit( a_desc , m_global , n_global , block_size , block_size , 0_pin , 0_pin , context , m_local , info )
    if ( info /= 0_pin ) then
       print*,'error: info /= 0'
       print*,'info = ' , info
       stop
    end if
    allocate( a_local( m_local , n_local ) )
    do i_local = 1_pin , m_local
       i_global = indxl2g( i_local , block_size , myrow , 0_pin , nprow )
       do j_local = 1_pin , n_local
          j_global = indxl2g( j_local , block_size , mycol , 0_pin , npcol )
          a_local( i_local , j_local ) = a_global( i_global , j_global )
       end do
    end do
    call blacs_barrier( context , 'all' )
  end subroutine parallel_predistribute
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_POSTDISTRIBUTE
  !*
  !* PURPOSE: gathering local pieces into a global matrix
  !*
  !* INPUTS:
  !*         - M_GLOBAL (integer): number of rows of the global matrix
  !*         - N_GLOBAL (integer): number of columns of the global matrix
  !*
  !* INPUTS/OUTPUTS: 
  !*                 - A_LOCAL (pointer to a real array): pointer to a local matrix
  !*
  !* OUTPUTS:
  !*         - A_GLOBAL (real array of size M_GLOBAL x N_GLOBAL): global matrix
  !*
  !* CALLS:
  !*        - INDXL2G (from PBLAS)
  !*        - SGSUM2D (from PBLAS)
  !*        - DGSUM2D (from PBLAS)
  !*  
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine parallel_postdistribute( m_global , n_global , a_global , a_local )
    implicit none
    integer( kind = pin ) , intent( in ) :: m_global
    integer( kind = pin ) , intent( in ) :: n_global
    real( kind = pre ) , dimension( m_global , n_global ) , intent( out ) :: a_global
    real( kind = pre ) , dimension( : , : ) , pointer :: a_local
    integer :: m_local
    integer :: n_local
    integer :: i_local
    integer :: i_global
    integer :: j_local
    integer :: j_global
    integer :: indxl2g
    m_local = size( a_local , 1_pin )
    n_local = size( a_local , 2_pin )
    a_global = 0.0_pre
    do i_local = 1_pin , m_local
       i_global = indxl2g( i_local , block_size , myrow , 0_pin , nprow )
       do j_local = 1_pin , n_local
          j_global = indxl2g( j_local , block_size , mycol , 0_pin , npcol )
          a_global( i_global , j_global ) = a_local( i_local , j_local )
       end do
    end do
    if ( pre == low ) then
       call sgsum2d( context , 'A' , 'T' , m_global , n_global , a_global , m_global , - 1_pin , - 1_pin )
    else if ( pre == high ) then
       call dgsum2d( context , 'A' , 'T' , m_global , n_global , a_global , m_global , - 1_pin , - 1_pin )
    end if
    if ( associated( a_local ) ) deallocate( a_local )
  end subroutine parallel_postdistribute
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_GETRF
  !* 
  !* PURPOSE: performs a LU factorization with partial pivoting
  !*
  !* INPUTS: 
  !*         - MA_GLOBAL (integer): number or rows of the global matrix
  !*         - NA_GLOBAL (integer): number of columns of the global matrix
  !*
  !* INPUTS/OUTPUTS:
  !*                 - A_GLOBAL (real array of size MA_GLOBAL x NA_GLOBAL): global matrix
  !*
  !* OUTPUTS:
  !*          - IPIV_GLOBAL (integer array of size MA_GLOBAL): permutations
  !*
  !* CALLS:
  !*        - PARALLEL_PREGRID (from PARALLEL)
  !*        - PARALLEL_PREDISTRIBUTE (from PARALLEL)
  !*        - PARALLEL_POSTDISTRIBUTE (from PARALLEL)
  !*        - PARALLEL_POSTGRID (from PARALLEL)
  !*        - PSGETRF (from SCALAPACK)
  !*        - PDGETRF (from SCALAPACK)
  !*        - DESCSET (from SCALAPACK)
  !*        - PSLAPIV (from SCALAPACK)
  !*        - PDLAPIV (from SCALAPACK)
  !*        - DGETRF (from LAPACK)
  !*        - SGETRF (from LAPACK)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the PDGETRF and PSGETRF subroutines
  !*
  !*******************************************************************************************************************************
  subroutine parallel_getrf( ma_global , na_global , a_global , ipiv_global )
    implicit none
    integer( kind = pin ) , intent( in ) :: ma_global
    integer( kind = pin ) , intent( in ) :: na_global
    real( kind = pre ) , dimension( ma_global , na_global ) , intent( inout ) :: a_global
    integer( kind = pin ) , dimension( ma_global ) , intent( out ) :: ipiv_global
    integer( kind = pin ) :: ma_local
    integer( kind = pin ) :: na_local
    integer( kind = pin ) :: mp_local
    integer( kind = pin ) :: np_local
    integer( kind = pin ) :: i
    integer( kind = pin ) :: info
    integer( kind = pin ) :: numroc
    integer( kind = pin ) , dimension( 9_pin ) :: a_desc
    integer( kind = pin ) , dimension( 9_pin ) :: piv_desc
    integer( kind = pin ) , dimension( 9_pin ) :: ipiv_desc
    integer( kind = pin ) , allocatable , dimension( : ) :: iwork
    integer( kind = pin ) , allocatable , dimension( : ) :: ipiv_local
    real( kind = pre ) , allocatable , dimension( : , : ) :: piv_global
    real( kind = pre ) , dimension( : , : ) , pointer :: a_local
    real( kind = pre ) , dimension( : , : ) , pointer :: piv_local 
    if ( nproc > 1_pin ) then
       allocate( piv_global( ma_global , 1_pin ) )
       do i = 1_pin , ma_global
          piv_global( i , 1_pin ) = i * 1.0_pre
       end do
       call parallel_pregrid()
       call parallel_predistribute( ma_global , na_global , a_global , a_desc , ma_local , na_local , a_local )
       call parallel_predistribute( ma_global , 1_pin , piv_global , piv_desc , mp_local , np_local , piv_local )
       allocate( ipiv_local( ma_local + block_size ) )
       if ( pre == low ) then
          call psgetrf( ma_global , na_global , a_local , 1_pin , 1_pin , a_desc , ipiv_local , info )
       else if ( pre == high ) then
          call pdgetrf( ma_global , na_global , a_local , 1_pin , 1_pin , a_desc , ipiv_local , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       allocate( iwork( nprow * npcol ) )
       call descset( ipiv_desc , a_desc( 3_pin ) + a_desc( 5_pin ) * nprow , 1_pin , a_desc( 5_pin ) , 1_pin , a_desc( 7_pin ) , &
            &mycol , context , a_desc( 5_pin ) + numroc( a_desc( 3_pin ) , a_desc( 5_pin ) , myrow , a_desc( 7_pin ) , nprow ) )
       if ( pre == low ) then
          call pslapiv( 'F' , 'R' , 'C' , ma_global , 1_pin , piv_local , 1_pin , 1_pin , piv_desc , ipiv_local , 1_pin , 1_pin ,&
               & ipiv_desc , iwork )
       else if ( pre == high ) then
          call pdlapiv( 'F' , 'R' , 'C' , ma_global , 1_pin , piv_local , 1_pin , 1_pin , piv_desc , ipiv_local , 1_pin , 1_pin ,&
               & ipiv_desc , iwork )
       end if
       deallocate( iwork )
       deallocate( ipiv_local )
       call parallel_postdistribute( ma_global , na_global , a_global , a_local )
       call parallel_postdistribute( ma_global , 1_pin , piv_global , piv_local )
       call parallel_postgrid()
       ipiv_global = int( piv_global( : , 1_pin ) )
       deallocate( piv_global )
    else
       if ( pre == high ) then
          call dgetrf( ma_global , na_global , a_global , ma_global , ipiv_global , info )
       else if ( pre == low ) then
          call sgetrf( ma_global , na_global , a_global , ma_global , ipiv_global , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
    end if
  end subroutine parallel_getrf
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_GESV
  !*
  !* PURPOSE: get a solution of a linear system
  !*
  !* INPUTS:
  !*         - NA_GLOBAL (integer): number of rows of the matrix
  !*         - NRHS_GLOBAL (integer): number of right hand sides
  !*
  !* INPUTS/OUTPUTS:
  !*                 - A_GLOBAL (real array of size NA_GLOBAL x NA_GLOBAL): on entry the coefficient matrix, on exit the factors&
  !*                                                                        & of the LU decomposition 
  !*                 - B_GLOBAL (real array of size NA_GLOBAL x NRHS_GLOBAL): on entry the right hand side matrix, on exit the s&
  !*                                                                          &olution of the system
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - PARALLEL_PREGRID (from PARALLEL)
  !*        - PARALLEL_PREDISTRIBUTE (from PARALLEL)
  !*        - PARALLEL_POSTDISTRIBUTE (from PARALLEL)
  !*        - PARALLEL_POSTGRID (from PARALLEL)
  !*        - PSGESV (from SCALAPACK)
  !*        - PDGESV (from SCALAPACK)
  !*        - DESCSET (from SCALAPACK)
  !*        - PSLAPIV (from SCALAPACK)
  !*        - PDLAPIV (from SCALAPACK)
  !*        - DGESV (from LAPACK)
  !*        - SGESV (from LAPACK)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the PDGESV and PSGESV subroutines
  !*
  !*******************************************************************************************************************************
  subroutine parallel_gesv( na_global , a_global , nrhs_global , b_global )
    implicit none
    integer( kind = pin ) , intent( in ) :: na_global
    real( kind = pre ) , dimension( na_global , na_global ) , intent( inout ) :: a_global
    integer( kind = pin ) , intent( in ) :: nrhs_global
    real( kind = pre ) , dimension( na_global , nrhs_global ) , intent( inout ) :: b_global
    integer( kind = pin ) , allocatable , dimension( : ) :: ipiv_global
    integer( kind = pin ) :: info
    integer( kind = pin ) :: ma_local
    integer( kind = pin ) :: na_local
    integer( kind = pin ) :: mb_local
    integer( kind = pin ) :: nb_local
    integer( kind = pin ) :: mp_local
    integer( kind = pin ) :: np_local
    integer( kind = pin ) :: i
    integer( kind = pin ) :: numroc
    integer( kind = pin ) , allocatable , dimension( : ) :: iwork
    integer( kind = pin ) , allocatable , dimension( : ) :: ipiv_local
    integer( kind = pin ) , dimension( 9_pin ) :: a_desc
    integer( kind = pin ) , dimension( 9_pin ) :: b_desc
    integer( kind = pin ) , dimension( 9_pin ) :: piv_desc
    integer( kind = pin ) , dimension( 9_pin ) :: ipiv_desc
    real( kind = pre ) , dimension( : , : ) , pointer :: a_local
    real( kind = pre ) , dimension( : , : ) , pointer :: b_local
    real( kind = pre ) , dimension( : , : ) , pointer :: piv_local
    real( kind = pre ) , allocatable , dimension( : , : ) :: piv_global
    if ( nproc > 1_pin ) then
       allocate( piv_global( na_global , 1_pin ) )
       do i = 1_pin , na_global
          piv_global( i , 1_pin ) = i * 1.0_pre
       end do
       call parallel_pregrid()
       call parallel_predistribute( na_global , na_global , a_global , a_desc , ma_local , na_local , a_local )
       call parallel_predistribute( na_global , nrhs_global , b_global , b_desc , mb_local , nb_local , b_local )
       call parallel_predistribute( na_global , 1_pin , piv_global , piv_desc , mp_local , np_local , piv_local )
       allocate( ipiv_local( ma_local + block_size ) )
       if ( pre == low ) then
          call psgesv( na_global , nrhs_global , a_local , 1_pin , 1_pin , a_desc , ipiv_local , b_local , 1_pin , 1_pin , b_desc&
               & , info )
       else if ( pre == high ) then
          call pdgesv( na_global , nrhs_global , a_local , 1_pin , 1_pin , a_desc , ipiv_local , b_local , 1_pin , 1_pin , b_desc&
               & , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       allocate( iwork( nprow * npcol ) )
       call descset( ipiv_desc , a_desc( 3_pin ) + a_desc( 5_pin ) * nprow , 1_pin , a_desc( 5_pin ) , 1_pin , a_desc( 7_pin ) , &
            &mycol , context , a_desc( 5_pin ) + numroc( a_desc( 3_pin ) , a_desc( 5_pin ) , myrow , a_desc( 7_pin ) , nprow ) )
       if ( pre == low ) then
          call pslapiv( 'F' , 'R' , 'C' , na_global , 1_pin , piv_local , 1_pin , 1_pin , piv_desc , ipiv_local , 1_pin , 1_pin ,&
               & ipiv_desc , iwork )
       else if ( pre == high ) then
          call pdlapiv( 'F' , 'R' , 'C' , na_global , 1_pin , piv_local , 1_pin , 1_pin , piv_desc , ipiv_local , 1_pin , 1_pin ,&
               & ipiv_desc , iwork )
       end if
       deallocate( iwork )
       deallocate( ipiv_local )
       deallocate( a_local )
       call parallel_postdistribute( na_global , nrhs_global , b_global , b_local )
       call parallel_postdistribute( na_global , 1_pin , piv_global , piv_local )
       call parallel_postgrid()
       deallocate( piv_global )
    else
       allocate( ipiv_global( na_global ) )
       if ( pre == high ) then
          call dgesv( na_global , nrhs_global , a_global , na_global , ipiv_global , b_global , na_global , info )
       else if ( pre == low ) then
          call sgesv( na_global , nrhs_global , a_global , na_global , ipiv_global , b_global , na_global , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       deallocate( ipiv_global )
    end if
  end subroutine parallel_gesv
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_GEMM
  !*
  !* PURPOSE: matrix multiplication C = alpha * op( A ) * op( B ) + beta * C
  !*
  !* INPUTS:
  !*         - TRANSA (character*1): indicates if A_GLOBAL has to be transposed or not
  !*         - TRANSB (character*1): indicates if B_GLOBAL has to be transposed or not
  !*         - MA_GLOBAL (integer): number of rows of the matrix A_GLOBAL
  !*         - NA_GLOBAL (integer): number of columns of the matrix B_GLOBAL
  !*         - A_GLOBAL (real array of size MA_GLOBAL x NA_GLOBAL): matrix A
  !*         - MB_GLOBAL (integer): number of rows of the matrix B_GLOBAL
  !*         - NB_GLOBAL (integer): number of columns of the matrix B_GLOBAL
  !*         - B_GLOBAL (real array of size MB_GLOBAL x NB_GLOBAL): matrix B
  !*         - MC_GLOBAL (integer): number of rows of the matrix C_GLOBAL
  !*         - NC_GLOBAL (integer): number of columns of the matrix C_GLOBAL
  !*         - ALPHA (real): alpha
  !*         - BETA (real): beta
  !*
  !* INPUTS/OUTPUTS:
  !*                 - C_GLOBAL (real array of size MC_GLOBAL x NC_GLOBAL): matrix C
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - PARALLEL_PREGRID (from PARALLEL)
  !*        - PARALLEL_PREDISTRIBUTE (from PARALLEL)
  !*        - PARALLEL_POSTDISTRIBUTE (from PARALLEL)
  !*        - PARALLEL_POSTGRID (from PARALLEL)
  !*        - PSGEMM (from PBLAS)
  !*        - PDGEMM (from PBLAS)
  !*        - DGEMM (from BLAS)
  !*        - SGEMM (from BLAS)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the PDGEMM and PSGEMM subroutines 
  !*
  !*******************************************************************************************************************************
  subroutine parallel_gemm( transa , transb , ma_global , na_global , a_global , mb_global , nb_global , b_global , mc_global , n&
       &c_global , c_global , alpha , beta )
    implicit none
    character( kind = pch , len = 1_pin ) , intent( in ) :: transa
    character( kind = pch , len = 1_pin ) , intent( in ) :: transb
    integer( kind = pin ) , intent( in ) :: ma_global
    integer( kind = pin ) , intent( in ) :: na_global
    real( kind = pre ) , dimension( ma_global , na_global ) , intent( in ) :: a_global
    integer( kind = pin ) , intent( in ) :: mb_global
    integer( kind = pin ) , intent( in ) :: nb_global
    real( kind = pre ) , dimension( mb_global , nb_global ) , intent( in ) :: b_global
    integer( kind = pin ) , intent( in ) :: mc_global
    integer( kind = pin ) , intent( in ) :: nc_global
    real( kind = pre ) , dimension( mc_global , nc_global ) , intent( inout ) :: c_global
    real( kind = pre ) , intent( in ) :: alpha
    real( kind = pre ) , intent( in ) :: beta
    integer( kind = pin ) :: ma_local
    integer( kind = pin ) :: na_local
    integer( kind = pin ) :: mb_local
    integer( kind = pin ) :: nb_local
    integer( kind = pin ) :: mc_local
    integer( kind = pin ) :: nc_local
    integer( kind = pin ) , dimension( 9_pin ) :: a_desc
    integer( kind = pin ) , dimension( 9_pin ) :: b_desc
    integer( kind = pin ) , dimension( 9_pin ) :: c_desc
    real( kind = pre ) , dimension( : , : ) , pointer :: a_local
    real( kind = pre ) , dimension( : , : ) , pointer :: b_local
    real( kind = pre ) , dimension( : , : ) , pointer :: c_local
    if ( nproc > 1_pin ) then
       call parallel_pregrid()
       call parallel_predistribute( ma_global , na_global , a_global , a_desc , ma_local , na_local , a_local )
       call parallel_predistribute( mb_global , nb_global , b_global , b_desc , mb_local , nb_local , b_local )
       call parallel_predistribute( mc_global , nc_global , c_global , c_desc , mc_local , nc_local , c_local )
       if ( transa == 'T' .and. transb == 'N' ) then
          if ( pre == low ) then
             call psgemm( transa , transb , na_global , nc_global , mb_global , alpha , a_local , 1_pin , 1_pin , a_desc , b_loca&
                  &l , 1_pin , 1_pin , b_desc , beta , c_local , 1_pin , 1_pin , c_desc )
          else if ( pre == high ) then
             call pdgemm( transa , transb , na_global , nc_global , mb_global , alpha , a_local , 1_pin , 1_pin , a_desc , b_loca&
                  &l , 1_pin , 1_pin , b_desc , beta , c_local , 1_pin , 1_pin , c_desc )
          end if
       else if ( transa == 'N' .and. transb == 'T' ) then
          if ( pre == low ) then
             call psgemm( transa , transb , ma_global , nc_global , nb_global , alpha , a_local , 1_pin , 1_pin , a_desc , b_loca&
                  &l , 1_pin , 1_pin , b_desc , beta , c_local , 1_pin , 1_pin , c_desc )
          else if ( pre == high ) then
             call pdgemm( transa , transb , ma_global , nc_global , nb_global , alpha , a_local , 1_pin , 1_pin , a_desc , b_loca&
                  &l , 1_pin , 1_pin , b_desc , beta , c_local , 1_pin , 1_pin , c_desc )
          end if
       else if ( transa == 'N' .and. transb == 'N' ) then
          if ( pre == low ) then
             call psgemm( transa , transb , ma_global , nc_global , mb_global , alpha , a_local , 1_pin , 1_pin , a_desc , b_loca&
                  &l , 1_pin , 1_pin , b_desc , beta , c_local , 1_pin , 1_pin , c_desc )
          else if ( pre == high ) then
             call pdgemm( transa , transb , ma_global , nc_global , mb_global , alpha , a_local , 1_pin , 1_pin , a_desc , b_loca&
                  &l , 1_pin , 1_pin , b_desc , beta , c_local , 1_pin , 1_pin , c_desc )
          end if
       end if
       deallocate( a_local )
       deallocate( b_local )
       call parallel_postdistribute( mc_global , nc_global , c_global , c_local )
       call parallel_postgrid()
    else
       if ( transa == 'T' .and. transb == 'N' ) then
          if ( pre == high ) then
             call dgemm( transa , transb , na_global , nc_global , mb_global , alpha , a_global , mb_global , b_global , mb_globa&
                  &l , beta , c_global , mc_global )
          else if ( pre == low ) then
             call sgemm( transa , transb , na_global , nc_global , mb_global , alpha , a_global , mb_global , b_global , mb_globa&
                  &l , beta , c_global , mc_global )
          end if
       else if ( transa == 'N' .and. transb == 'T' ) then
          if ( pre == high ) then
             call dgemm( transa , transb , ma_global , nc_global , nb_global , alpha , a_global , ma_global , b_global , nc_globa&
                  &l , beta , c_global , ma_global )
          else if ( pre == low ) then
             call sgemm( transa , transb , ma_global , nc_global , nb_global , alpha , a_global , ma_global , b_global , nc_globa&
                  &l , beta , c_global , ma_global )
          end if
       else if ( transa == 'N' .and. transb == 'N' ) then
          if ( pre == high ) then
             call dgemm( transa , transb , ma_global , nc_global , mb_global , alpha , a_global , ma_global , b_global , mb_globa&
                  &l , beta , c_global , ma_global )
          else if ( pre == low ) then
             call sgemm( transa , transb , ma_global , nc_global , mb_global , alpha , a_global , ma_global , b_global , mb_globa&
                  &l , beta , c_global , ma_global ) 
          end if
       end if
    end if
  end subroutine parallel_gemm
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_OPERATOR
  !*
  !* PURPOSE: parallelization of certain loops including operators using the master-slave strategy
  !*
  !* INPUTS:
  !*         - L (integer): step time index
  !*         - MATRIXIN (real array): matrix before operator is applied
  !*         - OPTION (character*): option to choose the operator
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - MATRIXOUT (real array): matrix after operator is applied
  !*
  !* CALLS:
  !*        - MPI_RECV (from MPI)
  !*        - MPI_SEND (from MPI)
  !*        - MPI_BCAST (from MPI)
  !*        - MODEL_TANGMODEL (from MODEL)
  !*        - MODEL_MODEL (from MODEL)
  !*        - OBSERVATIONS_TANGOBSOP (from OBSERVATIONS)
  !*        - OBSERVATIONS_OBSOP (from OBSERVATIONS)
  !*  
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine parallel_operator( l , matrixin , matrixout , option )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( : , : ) , intent( in ) :: matrixin
    real( kind = pre ) , dimension( : , : ) , intent( out ) :: matrixout
    character( kind = pch , len = * ) , intent( in ) :: option
    integer( kind = pin ) :: rowsin
    integer( kind = pin ) :: colsin
    integer( kind = pin ) :: rowsout
    integer( kind = pin ) :: colsout
    integer( kind = pin ) :: j
    integer( kind = pin ) :: down
    integer( kind = pin ) :: source
    integer( kind = pin ) , dimension( mpi_status_size ) :: status
    real( kind = pre ) , allocatable , dimension( : ) :: vectoraux
    rowsin = size( matrixin , 1_pin ) 
    colsin = size( matrixin , 2_pin ) 
    rowsout = size( matrixout , 1_pin )
    colsout = size( matrixout , 2_pin )
    allocate( vectoraux( rowsout + 1_pin ) ) 
    if ( nproc > 1_pin ) then 
       if ( rank == 0_pin ) then
          down = 0_pin
          j = 1_pin
          do 
             if ( pre == low ) then
                call mpi_recv( vectoraux , rowsout + 1_pin , mpi_real , mpi_any_source , mpi_any_tag , mpi_comm_world , status , &
                     &ierror )
             else if ( pre == high ) then
                call mpi_recv( vectoraux , rowsout + 1_pin , mpi_double_precision , mpi_any_source , mpi_any_tag , mpi_comm_world&
                     & , status , ierror )
             end if
             if ( int( vectoraux( rowsout + 1_pin ) ) /= 0_pin ) then
                matrixout( 1_pin : rowsout , int( vectoraux( rowsout + 1_pin ) ) ) = vectoraux( 1_pin : rowsout ) 
             end if
             source = status( mpi_source )
             if ( j <= colsout ) then
                call mpi_send( j , 1_pin , mpi_integer , source , 0_pin , mpi_comm_world , ierror )
                j = j + 1_pin
             else
                call mpi_send( colsout + 1_pin , 1_pin , mpi_integer , source , 0_pin , mpi_comm_world , ierror ) 
                down = down + 1_pin
             end if
             if ( down == nproc - 1_pin ) exit
          end do
       else
          vectoraux = 0.0_pre
          if ( pre == low ) then
             call mpi_send( vectoraux , rowsout + 1_pin , mpi_real , 0_pin , 0_pin , mpi_comm_world , ierror ) 
          else if ( pre == high ) then
             call mpi_send( vectoraux , rowsout + 1_pin , mpi_double_precision , 0_pin , 0_pin , mpi_comm_world , ierror ) 
          end if
          do
             call mpi_recv( j , 1_pin , mpi_integer , mpi_any_source , mpi_any_tag , mpi_comm_world , status , ierror )
             if ( j > colsout ) then 
                exit
             else
                if ( option == 'tangmodelT' ) then
                   call model_tangmodel( l , matrixin( j , 1_pin : colsin ) , vectoraux( 1_pin : rowsout ) )
                else if ( option == 'tangmodelN' ) then
                   call model_tangmodel( l , matrixin( 1_pin : rowsin , j ) , vectoraux( 1_pin : rowsout ) )
                else if ( option == 'tangobsopT' ) then
                   call observations_tangobsop( l , matrixin( j , 1_pin : colsin ) , vectoraux( 1_pin : rowsout ) )
                else if ( option == 'tangobsopN' ) then
                   call observations_tangobsop( l , matrixin( 1_pin : rowsin , j ) , vectoraux( 1_pin : rowsout ) )
                else if ( option == 'modelN' ) then
                   call model_model( l , matrixin( 1_pin : rowsin , j ) , vectoraux( 1_pin : rowsout ) )
                else if ( option == 'obsopN' ) then
                   call observations_obsop( l , matrixin( 1_pin : rowsin , j ) , vectoraux( 1_pin : rowsout ) )
                end if
                vectoraux( rowsout + 1_pin ) = j * 1.0_pre
             end if
             if ( pre == low ) then
                call mpi_send( vectoraux , rowsout + 1_pin , mpi_real , 0_pin , 0_pin , mpi_comm_world , ierror )
             else if ( pre == high ) then
                call mpi_send( vectoraux , rowsout + 1_pin , mpi_double_precision , 0_pin , 0_pin , mpi_comm_world , ierror )
             end if
          end do
       end if
    else
       do j = 1_pin , colsout
          if ( option == 'tangmodelT' ) then
             call model_tangmodel( l , matrixin( j , 1_pin : colsin ) , matrixout( 1_pin : rowsout , j ) )
          else if ( option == 'tangmodelN' ) then
             call model_tangmodel( l , matrixin( 1_pin : rowsin , j ) , matrixout( 1_pin : rowsout , j ) )
          else if ( option == 'tangobsopT' ) then
             call observations_tangobsop( l , matrixin( j , 1_pin : colsin ) , matrixout( 1_pin : rowsout , j ) )
          else if ( option == 'tangobsopN' ) then
             call observations_tangobsop( l , matrixin( 1_pin : rowsin , j ) , matrixout( 1_pin : rowsout , j ) )
          else if ( option == 'modelN' ) then
             call model_model( l , matrixin( 1_pin : rowsin , j ) , matrixout( 1_pin : rowsout , j ) )
          else if ( option == 'obsopN' ) then
             call observations_obsop( l , matrixin( 1_pin : rowsin , j ) , matrixout( 1_pin : rowsout , j ) )
          end if
       end do
    end if
    deallocate( vectoraux )
    if ( pre == low ) then
       call mpi_bcast( matrixout , rowsout * colsout , mpi_real , 0_pin , mpi_comm_world , ierror )
    else if ( pre == high ) then
       call mpi_bcast( matrixout , rowsout * colsout , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
    end if
  end subroutine parallel_operator
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_GESVD
  !*
  !* PURPOSE: performs a singular value decomposition A = U * SIGMA * VT
  !* 
  !* INPUTS:
  !*         - MA_GLOBAL (integer): number or rows of A
  !*         - NA_GLOBAL (integer): number or columns of A
  !*         - MU_GLOBAL (integer): number of rows of U
  !*         - NU_GLOBAL (integer): number of columns of U
  !*         - MVT_GLOBAL (integer): number of rows of VT
  !*         - NVT_GLOBAL (integer): number of columns of VT
  !*         - MSV_GLOBAL (integer): number of singular values
  !*
  !* INPUTS/OUTPUTS:
  !*                 - JOBU (character*1): flag to determine if parts of the decomposition are skipped
  !*                 - JOBVT (character*1): flag to determine if parts of the decomposition are skipped
  !*                 - A_GLOBAL (real array of size MA_GLOBAL x NA_GLOBAL): matrix A
  !*
  !* OUTPUTS:
  !*          - U_GLOBAL (real array of size MU_GLOBAL x NU_GLOBAL): matrix U
  !*          - VT_GLOBAL (real array of size MVT_GLOBAL x NVT_GLOBAL): matrix VT
  !*          - SV_GLOBAL (real array of size MSV_GLOBAL): vector of singular values
  !*
  !* CALLS:
  !*                   - PARALLEL_PREGRID (from PARALLEL)
  !*                   - PARALLEL_PREDISTRIBUTE (from PARALLEL)
  !*                   - PARALLEL_POSTDISTRIBUTE (from PARALLEL)
  !*                   - PARALLEL_POSTGRID (from PARALLEL)
  !*                   - PDGESVD (from SCALAPACK)
  !*                   - PSGESVD (from SCALAPACK)
  !*                   - DGESVD (from LAPACK)
  !*                   - SGESVD (from LAPACK)
  !*  
  !* COMMENTS: this subroutine is a wrapper of the PSGESVD and PDGESVD subroutines
  !*
  !*******************************************************************************************************************************
  subroutine parallel_gesvd( jobu , jobvt , ma_global , na_global , a_global , mu_global , nu_global , u_global , mvt_global , n&
       &vt_global , vt_global , msv_global , sv_global )
    implicit none
    character( kind = pch , len = 1_pin ) , intent( inout ) :: jobu
    character( kind = pch , len = 1_pin ) , intent( inout ) :: jobvt
    integer( kind = pin ) , intent( in ) :: ma_global
    integer( kind = pin ) , intent( in ) :: na_global
    real( kind = pre ) , dimension( ma_global , na_global ) , intent( inout ) :: a_global
    integer( kind = pin ) , intent( in ) :: mu_global
    integer( kind = pin ) , intent( in ) :: nu_global
    real( kind = pre ) , dimension( mu_global , nu_global ) , intent( out ) :: u_global
    integer( kind = pin ) , intent( in ) :: mvt_global
    integer( kind = pin ) , intent( in ) :: nvt_global
    real( kind = pre ) , dimension( mvt_global , nvt_global ) , intent( out ) :: vt_global
    integer( kind = pin ) , intent( in ) :: msv_global
    real( kind = pre ) , dimension( msv_global ) , intent( out ) :: sv_global
    integer( kind = pin ) :: lwork , info
    real( kind = pre ) , allocatable , dimension( : ) :: work
    integer( kind = pin ) :: ma_local
    integer( kind = pin ) :: na_local
    integer( kind = pin ) :: mu_local
    integer( kind = pin ) :: nu_local
    integer( kind = pin ) :: mvt_local
    integer( kind = pin ) :: nvt_local
    integer( kind = pin ) , dimension( 9_pin ) :: a_desc
    integer( kind = pin ) , dimension( 9_pin ) :: u_desc
    integer( kind = pin ) , dimension( 9_pin ) :: vt_desc
    real( kind = pre ) , dimension( : , : ) , pointer :: a_local
    real( kind = pre ) , dimension( : , : ) , pointer :: u_local
    real( kind = pre ) , dimension( : , : ) , pointer :: vt_local
    if ( nproc > 1_pin ) then
       call parallel_pregrid()
       if ( jobu == 'A' ) jobu = 'V'
       if ( jobvt == 'A' ) jobvt = 'V'
       call parallel_predistribute( ma_global , na_global , a_global , a_desc , ma_local , na_local , a_local )
       call parallel_predistribute( mu_global , nu_global , u_global , u_desc , mu_local , nu_local , u_local )
       call parallel_predistribute( mvt_global , nvt_global , vt_global , vt_desc , mvt_local , nvt_local , vt_local )
       allocate( work( 1_pin ) )
       !* jobu   = jobu
       !* jobvt  = jobvt
       !* m      = ma_global
       !* n      = na_global
       !* a      = a_local
       !* ia     = 1_pin
       !* ja     = 1_pin
       !* desca  = a_desc
       !* s      = sv_global
       !* u      = u_local
       !* iu     = 1_pin
       !* ju     = 1_pin
       !* descu  = u_desc
       !* vt     = vt_local
       !* ivt    = 1_pin
       !* jvt    = 1_pin
       !* descvt = vt_desc
       !* work   = work
       !* lwork  = - 1_pin
       !* info   = info   
       if ( pre == high ) then
          call pdgesvd( jobu , jobvt , ma_global , na_global , a_local , 1_pin , 1_pin , a_desc , sv_global , u_local , 1_pin , 1&
               &_pin , u_desc , vt_local , 1_pin , 1_pin , vt_desc , work , - 1_pin , info )
       else if ( pre == low ) then
          call psgesvd( jobu , jobvt , ma_global , na_global , a_local , 1_pin , 1_pin , a_desc , sv_global , u_local , 1_pin , 1&
               &_pin , u_desc , vt_local , 1_pin , 1_pin , vt_desc , work , - 1_pin , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       lwork = ceiling( work( 1_pin ) )
       deallocate( work )
       allocate( work( lwork ) )
       !* jobu   = jobu
       !* jobvt  = jobvt
       !* m      = ma_global
       !* n      = na_global
       !* a      = a_local
       !* ia     = 1_pin
       !* ja     = 1_pin
       !* desca  = a_desc
       !* s      = sv_global
       !* u      = u_local
       !* iu     = 1_pin
       !* ju     = 1_pin
       !* descu  = u_desc
       !* vt     = vt_local
       !* ivt    = 1_pin
       !* jvt    = 1_pin
       !* descvt = vt_desc
       !* work   = work
       !* lwork  = lwork
       !* info   = info
       if ( pre == high ) then
          call pdgesvd( jobu , jobvt , ma_global , na_global , a_local , 1_pin , 1_pin , a_desc , sv_global , u_local , 1_pin , 1&
               &_pin , u_desc , vt_local , 1_pin , 1_pin , vt_desc , work , lwork , info )
       else if ( pre == low ) then
          call psgesvd( jobu , jobvt , ma_global , na_global , a_local , 1_pin , 1_pin , a_desc , sv_global , u_local , 1_pin , 1&
               &_pin , u_desc , vt_local , 1_pin , 1_pin , vt_desc , work , lwork , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       deallocate( work )
       call parallel_postdistribute( ma_global , na_global , a_global , a_local )
       call parallel_postdistribute( mu_global , nu_global , u_global , u_local )
       call parallel_postdistribute( mvt_global , nvt_global , vt_global , vt_local )
       call parallel_postgrid()
    else
       if ( jobu == 'V' ) jobu = 'A'
       if ( jobvt == 'V' ) jobvt = 'A'
       allocate( work( 1_pin ) )
       !* jobu  = jobu
       !* jobvt = jobvt
       !* m     = ma_global
       !* n     = na_global
       !* a     = a_global
       !* lda   = ma_global
       !* s     = sv_global
       !* u     = u_global
       !* ldu   = mu_global
       !* vt    = vt_global
       !* ldvt  = mvt_global
       !* work  = work
       !* lwork = - 1_pin
       !* info  = info
       if ( pre == high ) then
          call dgesvd( jobu , jobvt , ma_global , na_global , a_global , ma_global , sv_global , u_global , mu_global , vt_global&
               & , mvt_global , work , - 1_pin , info )
       else if ( pre == low ) then
          call sgesvd( jobu , jobvt , ma_global , na_global , a_global , ma_global , sv_global , u_global , mu_global , vt_global&
               & , mvt_global , work , - 1_pin , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       lwork = ceiling( work( 1_pin ) )
       deallocate( work )
       allocate( work( lwork ) )
       !* jobu  = jobu
       !* jobvt = jobvt
       !* m     = ma_global
       !* n     = na_global
       !* a     = a_global
       !* lda   = ma_global
       !* s     = sv_global
       !* u     = u_global
       !* ldu   = mu_global
       !* vt    = vt_global
       !* ldvt  = mvt_global
       !* work  = work
       !* lwork = lwork
       !* info  = info
       if ( pre == high ) then
          call dgesvd( jobu , jobvt , ma_global , na_global , a_global , ma_global , sv_global , u_global , mu_global , vt_global&
               & , mvt_global , work , lwork , info )
       else if ( pre == low ) then
          call sgesvd( jobu , jobvt , ma_global , na_global , a_global , ma_global , sv_global , u_global , mu_global , vt_global&
               & , mvt_global , work , lwork , info )
       end if
       if ( info /= 0_pin ) then
          print*,'error: info /= 0'
          print*,'info = ' , info
          stop
       end if
       deallocate( work )
    end if
  end subroutine parallel_gesvd
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: PARALLEL_SVDCUT
  !*
  !* PURPOSE: cut a matrix of size ROWS x COLSIN to a matrix of size ROWS x COLSOUT ( COLSIN >= COLSOUT ) taking the directions a&
  !*          &ssociated with the leading singular values
  !*
  !* INPUTS:
  !*         - ROWS (integer): number of rows of the matrix MATRIXIN
  !*         - COLSIN (integer): number of columns of the matrix MATRIXIN
  !*         - COLSOUT (integer): number of columns of the reduced matrix MATRIXOUT
  !*         - MATRIXIN (real array of size ROWS x COLSIN): matrix to be cut
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - MATRIXOUT (real array of size ROWS x COLSOUT): matrix MATRIXIN reduced to COLSOUT columns
  !*
  !* CALLS:
  !*        - PARALLEL_GEMM
  !*        - PARALLEL_GESVD
  !*  
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine parallel_svdcut( rows , colsin , colsout , matrixin , matrixout )
    implicit none
    integer( kind = pin ) , intent( in ) :: rows 
    integer( kind = pin ) , intent( in ) :: colsin
    integer( kind = pin ) , intent( in ) :: colsout
    real( kind = pre ) , dimension( rows , colsin ) , intent( in ) :: matrixin
    real( kind = pre ) , dimension( rows , colsout ) , intent( out ) :: matrixout
    real( kind = pre ) , allocatable , dimension( : ) :: singvalues
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1u
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1vt
    character( kind = pch , len = 1_pin ) :: jobu
    character( kind = pch , len = 1_pin ) :: jobvt
    if ( colsin < colsout ) then
       print*,'error: colsin < colsout'
       print*,'colsin = ' , colsin
       print*,'colsout = ' , colsout
       stop
    end if
    allocate( matrixaux1( colsin , colsin ) )
    call parallel_gemm( 'T' , 'N' , rows , colsin , matrixin , rows , colsin , matrixin , colsin , colsin , matrixaux1 , 1.0_pre &
         &, 0.0_pre )
    allocate( singvalues( colsin ) )
    allocate( matrixaux1u( colsin , colsin ) )
    allocate( matrixaux1vt( colsin , colsin ) )
    jobu = 'A'
    jobvt = 'N'
    call parallel_gesvd( jobu , jobvt , colsin , colsin , matrixaux1 , colsin , colsin , matrixaux1u , colsin , colsin , matrixa&
         &ux1vt , colsin , singvalues )
    deallocate( singvalues )
    deallocate( matrixaux1vt )
    deallocate( matrixaux1 )
    call parallel_gemm( 'N' , 'N' , rows , colsin , matrixin , colsin , colsout , matrixaux1u( 1_pin : colsin , 1_pin : colsout )&
         & , rows , colsout , matrixout , 1.0_pre , 0.0_pre )
    deallocate( matrixaux1u )
  end subroutine parallel_svdcut
  !*
end module parallel
