FUNCTION DSTEFLOPS( ) RESULT( FLOPS )
!
USE DSTEDEFINITIONS
INCLUDE 'f90papi.h'
!
!.. Function Result ..
REAL( KIND=PREC ) :: FLOPS
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEFLOPS is the interface for a flops counter.                             !
!                                                                              !
!==============================================================================!
!  
!.. Scalars ..
INTEGER*4 :: PAPI_INFO
INTEGER*8 :: PAPI_FLPINS
REAL :: PAPI_MFLOPS, PAPI_PTIME, PAPI_RTIME
!
!.. External Subroutine ..
EXTERNAL PAPIF_FLOPS
!.. Executable Statements ......................................................
!
!
!  Interface to PAPI FLOPS counter 
!
!  PAPI_RTIME  : the total real (wallclock) time since the first rate call
!  PAPI_PTIME  : the total process time since the first rate call         
!  PAPI_FLPINS : the total floating point instructions since the first rate call
!  PAPI_MFLOPS : millions of floating point operations per second achieved since 
!                the latest rate call
!  PAPI_INFO   : return code
!
CALL PAPIF_FLOPS( PAPI_RTIME, PAPI_PTIME, PAPI_FLPINS, &
                  PAPI_MFLOPS, PAPI_INFO )
FLOPS = REAL( PAPI_MFLOPS,PREC )
!
END FUNCTION DSTEFLOPS
SUBROUTINE DSTEPRNFLOP( N, FLOPS, ITEST, ICASE )
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: ICASE, ITEST, N
REAL( KIND=PREC ) :: FLOPS
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEPRNFLOP prints flops counting.                                          !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N      (input) INTEGER                                                      !
!         The dimension of the matrix.                                         !
!                                                                              !
!  FLOPS  (input) REAL( KIND=PREC )                                            !
!         Flops counting.                                                      !
!                                                                              !
!  ITEST  (input) INTEGER                                                      !
!         Identifies the routine associated to FLOPS.                          !
!                                                                              !
!  ICASE  (input) INTEGER                                                      !
!         Identifies the case associated to FLOPS.                             !
!                                                                              !
!==============================================================================!
!
!.. Executable Statements ......................................................
!
IF      ( ICASE == 0 ) THEN
        OPEN ( UNIT=40, FILE="dstetester.out.flops" )
ELSE IF ( ITEST == 0 ) THEN
        CLOSE( UNIT=40 )
ELSE
        WRITE( UNIT=40, FMT='(A,I4,3A,I6,A,1P,E12.4)' ) &
               'Case:', ICASE, ', ', IDTEST( ITEST ), ', N = ', N, &
               ', FLOPS =', FLOPS
END IF
!
END SUBROUTINE DSTEPRNFLOP
