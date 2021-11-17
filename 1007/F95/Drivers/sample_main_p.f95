! This file (sample_main_p.f95) contains a sample main program that calls
! QNSTOPP to minimize the Griewank function, and verify the installation.

PROGRAM SAMPLE_MAIN_P
USE QNSTOPP_MOD
IMPLICIT NONE

INTERFACE
  FUNCTION Obj_GR(c, iflag) RESULT(f)
    USE REAL_PRECISION, ONLY : R8
    REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
    INTEGER, INTENT(OUT):: iflag
    REAL(KIND = R8):: f
  END FUNCTION Obj_GR
END INTERFACE

! Local variables.
INTEGER:: EVALlim, ITERlim, STATUS

REAL(KIND=R8), DIMENSION(50):: LB = -20.0_R8 ! Lower bounds.
REAL(KIND=R8), DIMENSION(50):: UB = 30.0_R8  ! Upper bounds.
REAL(KIND=R8), DIMENSION(50):: XI            ! Initial start point.
REAL(KIND=R8):: FMIN  ! Final objective function value.

WRITE (*,'(/A/)') "This test may take several minutes to finish."

! Use of QNSTOPP for deterministic global optimization, MODE = 'G'.
EVALlim = 0
ITERlim = 100
XI = 1.0_R8
CALL QNSTOPP(50, LB, UB, Obj_GR, XI, 'G', FMIN, STATUS, SWITCH=3, &
  NSTART=10, N=75, MAX_ITER=ITERlim, MAX_EVAL=EVALlim, TAU=1.0_R8,&
  GAIN=10.0_R8, OMP=1)
IF ((STATUS .GE. 10) .OR. (FMIN > 1.0E-2_R8)) THEN
  WRITE (*,'(/A/A/)') "The deterministic global optimization failed;", &
    "the installation appears incorrect."
  STOP
ELSE
  WRITE (*,113) STATUS, FMIN, ITERlim, EVALlim
  113 FORMAT('Deterministic global optimization completed successfully &
    &with status',I3.2,'.',/ &
    'Minimum objective function value =',ES15.7,',',/ &
    'the number of iterations is',I8,', and for all start points ',/ &
    'the total number of objective function evaluations is',I10,'.'/)
END IF

! Use of QNSTOPP for (local) stochastic optimization, MODE = 'S'.
EVALlim = 0
ITERlim = 500
LB(1:2) = -1.0_R8 ! Lower bounds.
UB(1:2) = 7.0_R8  ! Upper bounds.
XI(1:2) = (/ 1.5_R8, 2.0_R8 /)
CALL QNSTOPP(2, LB(1:2), UB(1:2), Obj_GR, XI(1:2), 'S', FMIN, STATUS, &
  SWITCH=1, N=75, MAX_ITER=ITERlim, MAX_EVAL=EVALlim, TAU=10.0_R8, &
  GAMMAW=2.0_R8, ETA=0.1_R8, OMP=2)
IF ((STATUS .GE. 10) .OR. (FMIN > 1.0E-2_R8)) THEN
  WRITE (*,'(/A/)') "The stochastic optimization failed; &
    &the installation appears incorrect."
  STOP
ELSE
  WRITE (*,114) STATUS, FMIN, ITERlim, EVALlim
  114 FORMAT('Stochastic optimization completed successfully &
    &with status',I3.2,'.',/ &
    'Minimum objective function value =',ES15.7,',',/ &
    'the number of iterations is',I8,', and for all start points ',/ &
    'the total number of objective function evaluations is',I10,'.'/)
END IF

WRITE (*,'(/A/)') "The installation appears correct."

END PROGRAM SAMPLE_MAIN_P

FUNCTION  Obj_GR(c, iflag) RESULT(f)
! On input:
! c     - Point coordinates.
!
! On output:
! f     - Function value at 'c'.
! iflag - A flag that is used to indicate the status of the
!         function evaluation. It is 0 for normal status.
!
! Obj_GR: Griewank function.
! The function formula is
! f(c) = 1+sum(c(i)^2/d)-product(cos(c(i)/sqrt(i))),
! where d = 500.0 and i is the summation index ranging from 1 to N (the number
! of dimensions).
! The global minimum is f(c)=0 at c(:)=(0,...,0), when c is in [-20, 30]^N.
!
! This implementation is provided in:
!
! J. He, L. T. Watson, and M. Sosonkina, Algorithm 897: VTDIRECT95:
! serial and parallel codes for the global optimization algorithm DIRECT,
! ACM Trans. Math. Software, 36, Article 17 (2009), pp. 1--24.
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Dummy variables.
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
INTEGER, INTENT(OUT):: iflag
REAL(KIND = R8):: f

! Local variables.
INTEGER:: i
REAL(KIND = R8):: d
REAL(KIND = R8)::prod

! Assign a value to 'd'.
d = 500.0_R8
! Compute the product.
prod = 1.0_R8
DO i = 1, SIZE(c)
  prod = prod * COS( c(i)/SQRT(REAL(i,KIND=R8)) )
END DO
f = 1.0_R8 + DOT_PRODUCT(c,c) / d - prod
iflag = 0

RETURN
END FUNCTION Obj_GR
