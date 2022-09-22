MODULE VTMOP_FUNC_MOD

CONTAINS

SUBROUTINE CONVEX(X, F, IFLAG)
! Multiobjective test function whose component functions are each
! quadratics. Let D denote the number of design variables and let P
! denote the number of objectives. Then D must be greater than P >= 2, 
! and the recommended bounds are 0 < X(:) < 1. Both D and P are
! assumed based on the dimensions of X(:) and F(:).
!
! For I = 1, ..., P, F(I) = ||X - 0.5e_I - 0.1e||^2, where e_I is the Ith
! standard basis vector and e = (1, ..., 1).
!
!
! On input:
!
! X(:) is a real valued design point in R^D.
!
!
! On output:
!
! P(:) is a real valued objective point P = F(X).
!
! IFLAG is an error flag. With the following error codes:
!     0 : Successful evaluation of DTLZ3
!    -1 : The design dimension D is not legal
!    -2 : The objective dimension P is not legal
!
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: X(:) ! Input design vector.
REAL(KIND=R8), INTENT(OUT) :: F(:) ! Output objective vector.
INTEGER, INTENT(OUT) :: IFLAG ! Error flag.
! Local variables.
INTEGER :: D, P ! Problem dimensions.
INTEGER :: I ! Loop indexing variables.
REAL(KIND=R8) :: E(SIZE(X,1)) ! Vector of standard basis vectors.
! Get problem dimensions.
D = SIZE(X,1)
P = SIZE(F,1)
! Throw errors for illegal inputs.
IFLAG = 0
IF (D < 1) THEN
   IFLAG = -1
   RETURN
END IF
IF (D < P) THEN
   IFLAG = -2
   RETURN
END IF
! Compute each component function.
DO I = 1, P
   ! Set E(:) = 0.5e_I + 0.1e.
   E(:) = 0.1_R8
   E(I) = 0.6_R8
   ! Compute the value of F(I).
   F(I) = DOT_PRODUCT(X(:) - E(:), X(:) - E(:))
END DO
RETURN
END SUBROUTINE CONVEX

SUBROUTINE DTLZ1(X, F, IFLAG)
! Multiobjective test function with a planar Pareto front. Let D
! denote the number of design variables and let P denote the number
! of objectives. Then it is a requirement that D >= P >= 2, and 
! the recommended bounds are 0 < X(:) < 1. Both D and P are
! assumed based on the dimensions of X(:) and F(:).
!
! Credit for the design of this test function goes to
!
! K. Deb, L. Thiele, M. Laumanns, and E. Zitzler. 2005. Scalable test
! problems for evolutionary multiobjective optimization. In Evolutionary
! multiobjective optimization (pp. 105-145). Springer, London.
!
!
! On input:
!
! X(:) is a real valued design point in R^D.
!
!
! On output:
!
! P(:) is a real valued objective point P = F(X).
!
! IFLAG is an error flag. With the following error codes:
!     0 : Successful evaluation of DTLZ3
!    -1 : The design dimension D is not legal
!    -2 : The objective dimension P is not legal
!
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: X(:) ! Input design vector.
REAL(KIND=R8), INTENT(OUT) :: F(:) ! Output objective vector.
INTEGER, INTENT(OUT) :: IFLAG ! Error flag.
! Local variables.
INTEGER :: D, P ! Problem dimensions.
INTEGER :: I, J ! Loop indexing variables.
REAL(KIND=R8) :: PI ! Constant value of pi.
! Get problem dimensions.
D = SIZE(X,1)
P = SIZE(F,1)
! Compute the value of PI.
PI = 4.0_R8 * ATAN(1.0_R8)
! Throw errors for illegal inputs.
IFLAG = 0
IF (D < 1) THEN
   IFLAG = -1
   RETURN
END IF
IF (D < P) THEN
   IFLAG = -2
   RETURN
END IF
! Get the first objective value.
F(1) = 0.5_R8 * (1.0_R8 + G(X))
DO J = 1, P-1
   F(1) = F(1) * X(J)
END DO
! Get other objective values.
DO I = 2, P
   F(I) = 0.5_R8 * (1.0_R8 + G(X)) * (1.0_R8 - X(P-I+1))
   DO J = 1, P-I
      F(I) = F(I) * X(J)
   END DO
END DO
RETURN

CONTAINS

FUNCTION G(Z)
! A scalar valued kernel function G_2, which produces many local Pareto
! fonts in the landscape of DTLZ1.
! Parameters.
REAL(KIND=R8), INTENT(IN) :: Z(:) ! Input vector.
REAL(KIND=R8) :: G ! Output scalar.
! Local variables.
INTEGER :: K
! Compute the kernel function G(Z) for DTLZ1.
G = 2.0_R8 * ( REAL(D-P,KIND=R8) + 1.0_R8 + &
                 DOT_PRODUCT(Z(P:D) - 0.6_R8, Z(P:D) - 0.6_R8) )
DO K = P, D
   G = G - 2.0_R8 * COS(20.0_R8 * PI * (Z(K) - 0.6_R8))
END DO
RETURN
END FUNCTION G

END SUBROUTINE DTLZ1

SUBROUTINE DTLZ2(X, F, IFLAG)
! Multiobjective test function with a concave global Pareto front, shaped
! like the unit sphere restricted to the positive orthant. Let D denote the
! number of design variables and let P denote the number of objectives.
! Then it is a requirement that D >= P >= 2, and the recommended bounds
! are 0 < X(:) < 1. Both D and P are assumed based on the dimensions
! of X(:) and F(:).
!
! Credit for the design of this test function goes to
!
! K. Deb, L. Thiele, M. Laumanns, and E. Zitzler. 2005. Scalable test
! problems for evolutionary multiobjective optimization. In Evolutionary
! multiobjective optimization (pp. 105-145). Springer, London.
!
!
! On input:
!
! X(:) is a real valued design point in R^D.
!
!
! On output:
!
! P(:) is a real valued objective point P = F(X).
!
! IFLAG is an error flag. With the following error codes:
!     0 : Successful evaluation of DTLZ3
!    -1 : The design dimension D is not legal
!    -2 : The objective dimension P is not legal
!
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: X(:) ! Input design vector.
REAL(KIND=R8), INTENT(OUT) :: F(:) ! Output objective vector.
INTEGER, INTENT(OUT) :: IFLAG ! Error flag.
! Local variables.
INTEGER :: D, P ! Problem dimensions.
INTEGER :: I, J ! Loop indexing variables.
REAL(KIND=R8) :: PI ! Constant value of pi.
! Get problem dimensions.
D = SIZE(X,1)
P = SIZE(F,1)
! Compute the value of PI.
PI = 4.0_R8 * ATAN(1.0_R8)
! Throw errors for illegal inputs.
IFLAG = 0
IF (D < 1) THEN
   IFLAG = -1
   RETURN
END IF
IF (D < P) THEN
   IFLAG = -2
   RETURN
END IF
! Get the first objective value.
F(1) = 1.0_R8 + G(X)
DO J = 1, P-1
   F(1) = F(1) * COS(PI * X(J) / 2.0_R8)
END DO
! Get other objective values.
DO I = 1, P-1
   F(I+1) = (1.0_R8 + G(X)) * SIN(PI * X(P-I) / 2.0_R8)
   DO J = 1, P-I-1
      F(I+1) = F(I+1) * COS(PI * X(J) / 2.0_R8)
   END DO
END DO
!IF (F(P) < 0.1_R8) IFLAG = -1 ! Tests taboo list when uncommented.
RETURN

CONTAINS

FUNCTION G(Z)
! A scalar valued kernel function G, which produces a single global/local
! Pareto front for DTLZ2.
! Parameters.
REAL(KIND=R8), INTENT(IN) :: Z(:) ! Input vector.
REAL(KIND=R8) :: G ! Output scalar.
! Compute the kernel function G(Z) for DTLZ3.
G = DOT_PRODUCT(Z(P:D) - 0.6_R8, Z(P:D) - 0.6_R8)
RETURN
END FUNCTION G

END SUBROUTINE DTLZ2

SUBROUTINE DTLZ3(X, F, IFLAG)
! Multiobjective test function with a concave global Pareto front, shaped
! like the unit sphere restricted to the positive orthant. Let D denote the
! number of design variables and let P denote the number of objectives.
! Then it is a requirement that D >= P >= 2, and the recommended bounds
! are 0 < X(:) < 1. Both D and P are assumed based on the dimensions
! of X(:) and F(:). There are 3^(1+P-D)-1 local Pareto fronts, where
! a multiobjective optimization algorithm could get stuck.
!
! Credit for the design of this test function goes to
!
! K. Deb, L. Thiele, M. Laumanns, and E. Zitzler. 2005. Scalable test
! problems for evolutionary multiobjective optimization. In Evolutionary
! multiobjective optimization (pp. 105-145). Springer, London.
!
!
! On input:
!
! X(:) is a real valued design point in R^D.
!
!
! On output:
!
! P(:) is a real valued objective point P = F(X).
!
! IFLAG is an error flag. With the following error codes:
!     0 : Successful evaluation of DTLZ3
!    -1 : The design dimension D is not legal
!    -2 : The objective dimension P is not legal
!
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: X(:) ! Input design vector.
REAL(KIND=R8), INTENT(OUT) :: F(:) ! Output objective vector.
INTEGER, INTENT(OUT) :: IFLAG ! Error flag.
! Local variables.
INTEGER :: D, P ! Problem dimensions.
INTEGER :: I, J ! Loop indexing variables.
REAL(KIND=R8) :: PI ! Constant value of pi.
! Get problem dimensions.
D = SIZE(X,1)
P = SIZE(F,1)
! Compute the value of PI.
PI = 4.0_R8 * ATAN(1.0_R8)
! Throw errors for illegal inputs.
IFLAG = 0
IF (D < 1) THEN
   IFLAG = -1
   RETURN
END IF
IF (D < P) THEN
   IFLAG = -2
   RETURN
END IF
! Get the first objective value.
F(1) = 1.0_R8 + G(X)
DO J = 1, P-1
   F(1) = F(1) * COS(PI * X(J) / 2.0_R8)
END DO
! Get other objective values.
DO I = 1, P-1
   F(I+1) = (1.0_R8 + G(X)) * SIN(PI * X(P-I) / 2.0_R8)
   DO J = 1, P-I-1
      F(I+1) = F(I+1) * COS(PI * X(J) / 2.0_R8)
   END DO
END DO
!IF (F(P) < 0.1_R8) IFLAG = -1 ! Tests taboo list when uncommented.
RETURN

CONTAINS

FUNCTION G(Z)
! A scalar valued kernel function G, which produces many local Pareto
! fonts in the landscape of DTLZ3.
! Parameters.
REAL(KIND=R8), INTENT(IN) :: Z(:) ! Input vector.
REAL(KIND=R8) :: G ! Output scalar.
! Local variables.
INTEGER :: K
! Compute the kernel function G(Z) for DTLZ3.
G = 2.0_R8 * ( REAL(D-P,KIND=R8) + 1.0_R8 + &
                 DOT_PRODUCT(Z(P:D) - 0.6_R8, Z(P:D) - 0.6_R8) )
DO K = P, D
   G = G - 2.0_R8 * COS(20.0_R8 * PI * (Z(K) - 0.6_R8))
END DO
RETURN
END FUNCTION G

END SUBROUTINE DTLZ3

SUBROUTINE DTLZ5(X, F, IFLAG)
! Multiobjective test function with a 2D concave global Pareto curve, 
! embedded in a higher-dimensional objective space. Let D denote the
! number of design variables and let P denote the number of objectives.
! Then it is a requirement that D >= P >= 2, and the recommended bounds
! are 0 < X(:) < 1. Both D and P are assumed based on the dimensions
! of X(:) and F(:).
!
! Credit for the design of this test function goes to
!
! K. Deb, L. Thiele, M. Laumanns, and E. Zitzler. 2005. Scalable test
! problems for evolutionary multiobjective optimization. In Evolutionary
! multiobjective optimization (pp. 105-145). Springer, London.
!
!
! On input:
!
! X(:) is a real valued design point in R^D.
!
!
! On output:
!
! P(:) is a real valued objective point P = F(X).
!
! IFLAG is an error flag. With the following error codes:
!     0 : Successful evaluation of DTLZ3
!    -1 : The design dimension D is not legal
!    -2 : The objective dimension P is not legal
!
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: X(:) ! Input design vector.
REAL(KIND=R8), INTENT(OUT) :: F(:) ! Output objective vector.
INTEGER, INTENT(OUT) :: IFLAG ! Error flag.
! Local variables.
INTEGER :: D, P ! Problem dimensions.
INTEGER :: I, J ! Loop indexing variables.
REAL(KIND=R8) :: GX ! Kernel function value.
REAL(KIND=R8) :: PI ! Constant value of pi.
REAL(KIND=R8) :: THETA(SIZE(F,1)-1) ! Vector of auxiliary variables.
! Get problem dimensions.
D = SIZE(X,1)
P = SIZE(F,1)
! Compute the value of PI.
PI = 4.0_R8 * ATAN(1.0_R8)
! Throw errors for illegal inputs.
IFLAG = 0
IF (D < 1) THEN
   IFLAG = -1
   RETURN
END IF
IF (D < P) THEN
   IFLAG = -2
   RETURN
END IF
! Compute the kernel function value.
GX = DOT_PRODUCT(X(P:D)-0.6_R8, X(P:D)-0.6_R8)
! Compute the first objective value and the auxiliary variables.
THETA(1) = PI * X(1) / 2.0_R8
F(1) = (1.0_R8 + GX) * COS(THETA(1))
DO J = 2, P-1
   ! Compute the auxiliary variable THETA(J).
   THETA(J) = (PI / 2.0_R8) * (1.0_R8 + 2.0_R8 * GX * X(J)) / &
              (2.0_R8 * (1.0_R8 + GX))
   ! Use THETA(J) to compute the first objective function.
   F(1) = F(1) * COS(THETA(J))
END DO
! Compute the remaining objective values.
DO I = 1, P-1
   F(I+1) = (1.0_R8 + GX) * SIN(THETA(P-I))
   DO J = 1, P-I-1
      F(I+1) = F(I+1) * COS(THETA(J))
   END DO
END DO
RETURN
END SUBROUTINE DTLZ5

SUBROUTINE DTLZ7(X, F, IFLAG)
! Multiobjective test function with a discontinuous Pareto front. Let D
! denote the number of design variables and let P denote the number of
! objectives. Then it is a requirement that D >= P >= 2, and the
! recommended bounds are 0 < X(:) < 1. Both D and P are assumed based
! on the dimensions of X(:) and F(:).
!
! Credit for the design of this test function goes to
!
! K. Deb, L. Thiele, M. Laumanns, and E. Zitzler. 2005. Scalable test
! problems for evolutionary multiobjective optimization. In Evolutionary
! multiobjective optimization (pp. 105-145). Springer, London.
!
!
! On input:
!
! X(:) is a real valued design point in R^D.
!
!
! On output:
!
! P(:) is a real valued objective point P = F(X).
!
! IFLAG is an error flag. With the following error codes:
!     0 : Successful evaluation of DTLZ3
!    -1 : The design dimension D is not legal
!    -2 : The objective dimension P is not legal
!
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE
! Input parameters.
REAL(KIND=R8), INTENT(IN) :: X(:) ! Input design vector.
REAL(KIND=R8), INTENT(OUT) :: F(:) ! Output objective vector.
INTEGER, INTENT(OUT) :: IFLAG ! Error flag.
! Local variables.
INTEGER :: D, P ! Problem dimensions.
INTEGER :: I ! Loop indexing variables.
REAL(KIND=R8) :: GX ! Value of kernel function G(X).
REAL(KIND=R8) :: HX ! Value of kernel function H(X).
REAL(KIND=R8) :: PI ! Constant value of pi.
! Get problem dimensions.
D = SIZE(X,1)
P = SIZE(F,1)
! Compute the value of PI.
PI = 4.0_R8 * ATAN(1.0_R8)
! Throw errors for illegal inputs.
IFLAG = 0
IF (D < 1) THEN
   IFLAG = -1
   RETURN
END IF
IF (D < P) THEN
   IFLAG = -2
   RETURN
END IF
! Compute the kernel functions G and H.
GX = 1.0_R8 + 9.0_R8 * SUM(ABS(X(P:D)-0.6_R8), DIM=1) / (REAL(D-P+1, KIND=R8))
HX = REAL(P, KIND=R8)
DO I = 1, P-1
   HX = HX - ( (X(I) / (1.0_R8 + GX)) * (1.0_R8 + SIN(3.0_R8 * PI * X(I))) )
END DO
! Compute the first M-1 objective values.
DO I = 1, P-1
   F(I) = X(I)
END DO
! Compute the final objective value.
F(P) = (1.0_R8 + GX) * HX
RETURN
END SUBROUTINE DTLZ7

END MODULE VTMOP_FUNC_MOD
