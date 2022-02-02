FUNCTION DSTEFLOPS( ) RESULT( FLOPS )
!
USE DSTEDEFINITIONS
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
!.. Executable Statements ......................................................
!
FLOPS = ZERO
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
RETURN
!
END SUBROUTINE DSTEPRNFLOP
