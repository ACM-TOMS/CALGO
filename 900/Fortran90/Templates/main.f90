!*********************************************************************************************************************************
!*
!* PROGRAM: MAIN
!*
!* PURPOSE: main program 
!*
!* DEPENDENCIES:
!*                - PRECISION
!*                - RANDOM
!*                - MODEL
!*                - OBSERVATIONS
!*                - PARALLEL
!*                - INITIALIZE
!*                - EKF [ENKF] [RRSQRTKF] [RRSQRTENKF] 
!* 
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
!*
!*          - IDUM (integer): seed for random generators, it must be set to a negative value (from RANDOM)
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
!*          - IERROR (integer): variable to detect errors (from PARALLEL)
!*          - RANK (integer): which machine (from PARALLEL)
!*          - NPROC (integer): number of processors (from PARALLEL)
!*          - BLOCK_SIZE (integer): block size (from PARALLEL)
!*          - NPROW (integer): number of rows in the processor grid (from PARALLEL)
!*          - NPCOL (integer): number of columns in the processor grid (from PARALLEL)
!*          - CONTEXT (integer): it is a universe where messages exist and do not interact with other context's messages (from PA&
!*                               &RALLEL)
!*          - MYROW (integer): calling processor's row number in the processor grid (from PARALLEL)
!*          - MYCOL (integer): calling processor's column number in the processor grid (from PARALLEL)
!*
!*          - DIMSPACESTATE (integer): dimension of the space state (from INITIALIZE)
!*          - NUMBERSAMPLES (integer): number of samples (from INITIALIZE)
!*          - MODESANALYSIS (integer): number of analysis modes (from INITIALIZE)
!*          - FIRST (logical): flag to determine the first step (from INITIALIZE)
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
program main

   !* INCLUDING MODULES
   use precision
   use precision
   use random
   use model
   use observations
   use parallel
   use initialize
   use ekf
   ...
	
   implicit none

   !* DECLARATION OF VARIABLES
   ...

   !* INITIALIZATION OF PARALLELIZATION
   call parallel_init()

   !* GETTING PARALLELIZATION PARAMETERS
   call parallel_ranksize()

   !* SETTING PARAMETERS AND BROADCASTING
   if ( rank == 0_pin ) then
     call initialize_parametersup()
   end if
   call mpi_bcast(...)
   ...
 
   !* ALLOCATIONS
   allocate(...)
   ...   

   !* INITIALIZATION
   if ( rank == 0_pin ) then
      call initialize_initialize1( state , covariance )
   end if
   call mpi_bcast(...)

   !* LOOP

   do l = 1 , ...

      !* PREDICTION
      call ekf_predictor( l , state , covariance ) 

      !* SET IFOBS
      if ( rank == 0_pin ) then
         call observations_ifobservations( l + 1_pin )         
      end if
      call mpi_bcast(...)
   
      if ( ifobs ) then !* IN CASE THERE ARE OBSERVATIONS     
   
        !* SET NUMBEROBS    
        if ( rank == 0_pin ) then
           call observations_numberobs( l + 1_pin )
        end if
        call mpi_bcast(...)

        !* ALLOCATE OBSVALUE
        allocate( obsvalue( numberobs ) )

        !* SET OBSVALUE
        if ( rank == 0_pin ) then
           call observations_obsvalue( l + 1_pin )
        end if
        call mpi_bcast(...)
        
        !* CORRECTION
        call ekf_corrector( l , state , covariance )

        !* DEALLOCATE OBSVALUE
        deallocate( obsvalue )

      end if

   end do

   !* DEALLOCATIONS
   ...

   !* END PARALLELIZATION
   call parallel_finalize()

end program
