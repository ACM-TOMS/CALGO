! This file (objfunc.f95) contains five test objective functions for
! both VTdirect and pVTdirect.
!

! For all the objective functions:
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
! where d = 500.0 (a bigger d value gives more local minima) and
! i is the summation index ranging from 1 to N (the number of
! dimensions). The global minimum is f(c)=0 at c(:)=(0,...,0), when c
! is in [-20, 30]^N.
!
FUNCTION  Obj_GR(c, iflag) RESULT(f)
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

! Obj_QU: Quartic function.
! The function formula is f(c) = sum(2.2*(c(i)+e_i)^2-(c(i)-e_i)^4),
! where e_i is a uniformly distributed random number in [0.2,0.4]
! and i is the summation index ranging from 1 to N (the
! number of dimensions). Here, take all e_i = 0.3 for a
! deterministic result. The global minimum is f(c)=-29.186*N at
! c(:)=(3,...,3) when c is in [-2,3]^N.
!
FUNCTION Obj_QU(c, iflag) RESULT(f)
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Dummy variables.
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
INTEGER, INTENT(OUT):: iflag
REAL(KIND = R8):: f

! Local variables.
INTEGER:: i
REAL(KIND = R8):: e_i

! Assign a value to 'e_i'.
e_i = 0.3_R8
! Compute the sum.
f = 0.0_R8
DO i = 1, SIZE(c)
  f = f + 2.2_R8 * (c(i) + e_i)**2 - (c(i) - e_i)**4
END DO
iflag = 0

RETURN
END FUNCTION Obj_QU

! Obj_RO: Rosenbrock's valley function.
! The function is  f(c)=sum(100*(c(i+1)-c(i)^2)^2+(1-c(i))^2),
! where i is the summation index ranging from 1 to N-1 (the
! number of dimensions minus one). The global minimum is f(c)=0 at
! c(:)=(1,...,1), when c is in [-2.048, 2.048]^N.
!
FUNCTION Obj_RO(c, iflag) RESULT(f)
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Dummy variables.
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
INTEGER, INTENT(OUT):: iflag
REAL(KIND = R8):: f

! Local variables.
INTEGER:: i

! Compute the sum.
f = 0.0_R8
DO i = 1, SIZE(c)-1
  f = f + 100.0_R8 * (c(i+1) - c(i)**2)**2 + (1.0_R8 - c(i))**2
END DO
iflag = 0
RETURN
END FUNCTION Obj_RO

! Obj_SC: Schwefel's function.
! The function formula is  f(c)=sum(-c(i)*sin(sqrt(abs(c(i))))),
! where i is the summation index ranging from 1 to N (the number of
! dimensions). The global minimum is f(c)=-N*418.9829 at
! c(:)=(420.9687,...,420.9687), when c is in [-500,500]^N.
!
FUNCTION Obj_SC(c, iflag) RESULT(f)
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Dummy variables.
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
INTEGER, INTENT(OUT):: iflag
REAL(KIND = R8):: f

! Compute the sum.
f = -SUM( c(:)*SIN( SQRT(ABS(c(:))) ) )
iflag = 0

RETURN
END FUNCTION Obj_SC

! Obj_MI: Michalewicz's function.
! The function formula is
! f(c)=-sum(sin(c(i))*(sin(i*c(i)^2/pi))^(2*m)), where i is the
! summation index ranging form 1 to N (the number of dimensions),
! and m=10 (a bigger m value presents a narrower valley).
! The global minimum is f(c)=-4.687 (N=5) at
! c=(2.203, 1.571, 1.285, 1.923, 1.720),
! when c is in [0, pi]^N.
!
FUNCTION Obj_MI(c, iflag) RESULT(f)
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Dummy variables.
REAL(KIND = R8), DIMENSION(:), INTENT(IN):: c
INTEGER, INTENT(OUT):: iflag
REAL(KIND = R8):: f

! Local variables.
INTEGER:: i,m
REAL(KIND = R8), PARAMETER:: PI = 3.141592653589793_R8

! Assign a value to 'm'.
m = 10

! Compute the sum.
f = 0.0_R8
DO i = 1, SIZE(c)
  f = f + SIN(c(i)) * (SIN( i * c(i)**2 / PI ))**(2*m)
END DO
f = -f
iflag = 0

RETURN
END FUNCTION Obj_MI
