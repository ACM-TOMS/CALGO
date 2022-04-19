! This is a modified version of
! Marsaglia & Tsang's generator for random normals.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! A reference to the method can be found at,
! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for 
! generating random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08&k=LzsPWH25gRR3YYS3VmwyEA%3D%3D%0A&r=
! xJs6G7HNRwaeoHbudTmgPMrHSQzD6DDzJbkF1xvfwL0%3D%0A&m=
! Q9Ud7OCB4Hp12uzLJM01fJsXqpWqhleLOtKfC1eZk6w%3D%0A&s=
! 937d576cf1f4b636ff342bc454c87e8733112f4135a17b83d01458876c28c688
! Latest version - 1 January 2001

! This version was modified on July 2013 by Will Thacker, 
! Winthrop University.

! The changes made are:
!
! removed the exponential random number generator,
!
! removed the need for integers to be exactly 32-bits,
!
! removed the internal integer and real uniformly distributed random 
! number generators,
!
! changed the name of the function rnorm to a subroutine RANDOM_NORMAL,
!
! uses the FORTRAN intrinsic RANDOM_NUMBER rather than an internal
! generator (that was machine dependent),
!
! now returns an array of normal variates rather than one at a time.
!
! The initializing subroutine ZIGSET no longer sets a random seed 
! (instead, a user can call the intrinsic, RANDOM_SEED, before calling 
! RANDOM_NORMAL if you want something other than the default seed).
!
! The system is now set up as a module to more easily allow 
! initializing of tables only once, rather than reinitializing for 
! each call to RANDOM_NORMAL, and to facilitate sharing the tables.

MODULE ZIGARRAY
USE REAL_PRECISION
IMPLICIT NONE

PRIVATE
REAL(KIND = R8), PARAMETER:: HALF=0.5_R8
REAL(KIND = R8), PARAMETER:: M1=2147483648.0_R8
REAL(KIND = R8):: DN=3.442619855899_R8 
REAL(KIND = R8):: Q
REAL(KIND = R8):: TN=3.442619855899_R8
REAL(KIND = R8):: VN=0.00991256303526217_R8 
INTEGER, SAVE:: KN(0:127)
REAL(KIND = R8), SAVE:: FN(0:127), WN(0:127) 
LOGICAL, SAVE:: INITIALIZED=.FALSE.

PUBLIC:: ZIGSET, RANDOM_NORMAL

CONTAINS
SUBROUTINE RANDOM_NORMAL(RV) 
!  Generate an array of real random normals with a mean of 0 and a 
!  standard deviation of 1, using the intrinsic function RANDOM_NUMBER.  

!  On input:
!     No input values are required.
!
!  On output:
!  
!  RV(:) is a real array with the computed random normal variates.
!
REAL(KIND=R8), INTENT(OUT), DIMENSION(:):: RV 

!  Local variables

REAL(KIND = R8), PARAMETER:: R = 3.442620_R8
REAL(KIND = R8):: TEMP, X, Y ! Temporaries. 
REAL(KIND = R8):: TRV(2) ! Temporary for getting 2 random numbers.
INTEGER:: HZ, I, IZ, N ! Temporary integers.
LOGICAL:: DONE ! Tells when a value finally meets criteria for a normal variate.

!  Initialize the tables needed for the algorithm if they have not been 
!    set up yet.
IF( .NOT. INITIALIZED ) CALL ZIGSET()
 
! Get many uniform random real values.

CALL RANDOM_NUMBER(RV)
N = SIZE(RV)

! Go through the random uniforms and change them to random normals.

DO I=1,N
!  Compute a random (positive or negative) random 32-bit integer

  HZ = INT(((RV(I)*2.0_R8) - 1.0_R8)*M1)  
  IZ = IAND( HZ, 127 )   ! Get a random index into the tables
  IF(ABS(HZ) < KN(IZ)) THEN
     RV(I) = HZ * WN(IZ)
  ELSE
     DONE = .FALSE.
     DO WHILE (.NOT. DONE)
        IF ( IZ == 0 ) THEN
           DO
              CALL RANDOM_NUMBER(TRV)
              X = -0.2904764_R8*LOG(TRV(1))
              Y = -LOG(TRV(2))
              IF (Y+Y .GE. X*X) EXIT
           END DO
           RV(I) = R + X
           IF (HZ .LE.  0) RV(I) = -RV(I)
           DONE = .TRUE.
        ELSE
           X = HZ*WN(IZ)
           CALL RANDOM_NUMBER(TEMP)
           IF (FN(IZ) + TEMP*(FN(IZ-1) - FN(IZ)) < EXP(-HALF*X*X) ) THEN
              RV(I) = X
              DONE = .TRUE.
           ELSE
              CALL RANDOM_NUMBER(TEMP)
              HZ = INT(((TEMP*2.0_R8) - 1.0_R8)*M1)
              IZ = IAND(HZ, 127)
              IF (ABS(HZ) < KN(IZ)) THEN
                 RV(I) = HZ*WN(IZ)
                 DONE = .TRUE.
              END IF
           END IF
        END IF
     END DO
  END IF
END DO
RETURN
END SUBROUTINE RANDOM_NORMAL

SUBROUTINE ZIGSET()

! Local variables.
INTEGER  :: I ! Loop index.

!  Set up the tables used by RANDOM_NORMAL

!  This modified method still uses the original implementation scheme 
!  where the KN table is set up as an integer array for efficiency 
!  when random integers were generated.
!
!  A possible improvement might be to change KN to real and not have to 
!  generate random integers.  However, random integers are still needed 
!  for indexing into the arrays.

   Q = VN*EXP(HALF*DN*DN)
   KN(0) = INT((DN/Q)*M1)
   KN(1) = 0
   WN(0) = Q/M1
   WN(127) = DN/M1
   FN(0) = 1.0_R8
   FN(127) = EXP(-HALF*DN*DN )
   DO  I = 126,1,-1
      DN = SQRT(-2.0_R8*LOG(VN/DN + EXP(-HALF*DN*DN)))
      KN(I+1) = INT((DN/TN)*M1)
      TN = DN
      FN(I) = EXP(-HALF*DN*DN)
      WN(I) = DN/M1
   END DO

   INITIALIZED = .TRUE.
   RETURN
END SUBROUTINE ZIGSET


END MODULE ZIGARRAY
