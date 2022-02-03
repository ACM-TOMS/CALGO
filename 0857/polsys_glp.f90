! This file contains all the modules and external subroutines for the
! parallel version of the package POLSYS_GLP, except for the LAPACK
! routines used, which are distributed in a separate file.  Layne T.
! Watson, Steven M. Wise, Andrew J. Sommese, August, 1998.  Cosmetic
! changes, 10/1999.  Extension from POLSYS_PLP to POLSYS_GLP with
! MPI-based parallelization by Hai-Jun Su, J. Michael McCarthy, Masha
! Sosonkina, Layne T. Watson, August, 2004.

      MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
      INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
      END MODULE REAL_PRECISION

                                      !!!
MODULE GLOBAL_GLP

! The module GLOBAL_GLP contains derived data types, arrays, and
! functions used in POLSYS_GLP and related subroutines.  GLOBAL_GLP uses
! the HOMPACK90 module REAL_PRECISION for 64-bit arithmetic.
USE REAL_PRECISION, ONLY: R8
IMPLICIT NONE
INTEGER, PARAMETER:: LARGE=SELECTED_INT_KIND(15)
REAL (KIND=R8), PARAMETER:: PI=3.1415926535897932384626433_R8

! TARGET SYSTEM: Let X be a complex N-dimensional vector.  POLSYS_GLP
! is used to solve the polynomial system, called the target system,
! F(X)=0, where F is represented by the following derived data types:

TYPE TERM_TYPE
   COMPLEX (KIND=R8):: COEF
   INTEGER, DIMENSION(:), POINTER:: DEG
END TYPE TERM_TYPE
TYPE POLYNOMIAL_TYPE
   TYPE(TERM_TYPE), DIMENSION(:), POINTER:: TERM
   INTEGER:: NUM_TERMS
END TYPE POLYNOMIAL_TYPE
TYPE(POLYNOMIAL_TYPE), DIMENSION(:), ALLOCATABLE:: POLYNOMIAL

! The mathematical representation of the target system F is, for I=1,...,N,
!
! F_I(X) = SUM_{J=1}^{POLYNOMIAL(I)%NUM_TERMS}
!          POLYNOMIAL(I)%TERM(J)%COEF * 
!          PRODUCT_{K=1}^N  X(K)**POLYNOMIAL(I)%TERM(J)%DEG(K).
!
! Any program calling POLSYS_GLP (such as the sample main program
! MAIN_TEMPLATE) must aquire data and allocate storage for the target
! system as illustrated below:  
!
! ALLOCATE(POLYNOMIAL(N))
! DO I=1,N
!   READ (*,*) POLYNOMIAL(I)%NUM_TERMS
!   ALLOCATE(POLYNOMIAL(I)%TERM(POLYNOMIAL(I)%NUM_TERMS))
!   DO J=1,POLYNOMIAL(I)%NUM_TERMS
!     ALLOCATE(POLYNOMIAL(I)%TERM(J)%DEG(N+1))
!     READ (*,*) POLYNOMIAL(I)%TERM(J)%COEF,POLYNOMIAL(I)%TERM(J)%DEG(1:N)
!   END DO
! END DO
!
! START SYSTEM/COVER:  In a generalized linear product (GLP)
! formulation the start system G(X)=0 must have the same set structure 
! as the target system P(X)=0.  G and P are represented by the derived
! data types:

INTEGER, DIMENSION(:), ALLOCATABLE:: COVER_SIZES
TYPE SET_TYPE 
   INTEGER, DIMENSION(:), POINTER:: INDEX
   INTEGER:: NUM_INDICES
   INTEGER:: SET_DEG
   COMPLEX (KIND=R8), DIMENSION(:), POINTER:: START_COEF
END TYPE SET_TYPE
TYPE COVER_TYPE
   TYPE(SET_TYPE), DIMENSION(:), POINTER:: SET 
END TYPE COVER_TYPE
TYPE(COVER_TYPE), DIMENSION(:), ALLOCATABLE:: COVER

! The mathematical representation of the start system G is, for I=1,...,N,
!
! G_I(X) = PRODUCT_{J=1}^{COVER_SIZES(I)}
!          ( L(I,J)**COVER(I)%SET(J)%SET_DEG - 1.0 ),
!
! where the linear factors L(I,J) are
!
! L(I,J) = SUM_{K=1}^{COVER(I)%SET(J)%NUM_INDICES}
!   COVER(I)%SET(J)%START_COEF(K) * X(COVER(I)%SET(J)%INDEX(K)).
!
! The system covering P=(P(1),...,P(N)) is comprised of the
! coverings P(I) = {S(I,1),...S(I, COVER_SIZES(I))}, where the sets
! of variables S(I,J) are defined by
!
! S(I,J) = UNION_{K=1}^{COVER(I)%SET(J)%NUM_INDICES}
!          { X(COVER(I)%SET(J)%INDEX(K)) }.
!
! The calling program must acquire data and allocate storage as
! illustrated below:
!
! ALLOCATE(COVER_SIZES(N))
! READ (*,*) COVER_SIZES(1:N)
! ALLOCATE(COVER(N))
! DO I=1,N
!   ALLOCATE(COVER(I)%SET(COVER_SIZES(I))
!   DO J=1, COVER_SIZES(I)
!     READ (*,*) COVER(I)%SET(J)%NUM_INDICES
!     ALLOCATE(COVER(I)%SET(J)%INDEX(COVER(I)%SET(J)%NUM_INDICES))
!     READ (*,*) COVER(I)%SET(J)%INDEX
!   END DO
! END DO
!
! START_COEF is generated randomly by POLSYS_GLP.
! SET_DEG is provided in the input data file.
! Subroutine CHECK_GLP is used to check the validity of a given generalized 
! linear decomposition. 


CONTAINS

! INDEXING FUNCTIONS FOR THE TARGET SYSTEM:
!
! C(I,J) retrieves the coefficient of the Jth term of the Ith polynomial
! component of the target system.

COMPLEX (KIND=R8) FUNCTION C(I,J) 
  IMPLICIT NONE
  INTEGER:: I,J
  C = POLYNOMIAL(I)%TERM(J)%COEF
END FUNCTION C

! D(I,J,K) retrieves the degree of the Kth variable in the Jth term of
! the Ith polynomial component of the target system.

INTEGER FUNCTION D(I,J,K)
  IMPLICIT NONE
  INTEGER:: I,J,K
  D = POLYNOMIAL(I)%TERM(J)%DEG(K)
END FUNCTION D

! NUMT(I) retrieves the number of terms in the Ith polynomial component of
! the target system F(X).

INTEGER FUNCTION NUMT(I)
  IMPLICIT NONE
  INTEGER:: I
  NUMT = POLYNOMIAL(I)%NUM_TERMS 
END FUNCTION NUMT

! The target system is succinctly specified with the retrieval functions:
!
! F_I(X) = SUM_{J=1}^{NUMT(I)} C(I,J) * PRODUCT_{K=1}^N  X(K)**D(I,J,K).
!
! INDEXING FUNCTIONS FOR THE START SYSTEM/COVER:
!
! PAR(I,J,K) retrieves the index of the Kth variable in the Jth set
! S(I,J) of the Ith covering P(I).

INTEGER FUNCTION PAR(I,J,K)
  IMPLICIT NONE
  INTEGER:: I,J,K
  PAR = COVER(I)%SET(J)%INDEX(K)
END FUNCTION PAR

! SC(I,J,K) retrieves the coefficient of the variable with index
! PAR(I,J,K) in the Jth factor of the Ith component of the start system
! G(X).

COMPLEX (KIND=R8) FUNCTION SC(I,J,K)
  IMPLICIT NONE
  INTEGER:: I,J,K
  SC = COVER(I)%SET(J)%START_COEF(K)
END FUNCTION SC

! SD(I,J) retrieves the set degree of the Jth set S(I,J) in the Ith
! covering P(I).

INTEGER FUNCTION SD(I,J)
  IMPLICIT NONE
  INTEGER:: I,J
  SD = COVER(I)%SET(J)%SET_DEG
END FUNCTION SD

! NUMV(I,J) retrieves the number of variables in the Jth set S(I,J) of
! the Ith covering P(I).

INTEGER FUNCTION NUMV(I,J)
  IMPLICIT NONE
  INTEGER:: I,J
  NUMV = COVER(I)%SET(J)%NUM_INDICES
END FUNCTION NUMV

! Both the start system and the set structure are succinctly specified with
! retrieval functions: 
!
! G_I(X) = PRODUCT_{J=1}^{COVER_SIZES(I)}
! ( [ SUM_{K=1}^{NUMV(I,J)} SC(I,J,K)*X(PAR(I,J,K)) ]**SD(I,J) - 1.0 ),
!
! and P(I) = { S(I,1),...,S(I,COVER_SIZES(I)) }, where
!
! S(I,J) = UNION_{K=1}^{NUMV(I,J)} { X(PAR(I,J,K)) }.

END MODULE GLOBAL_GLP

                                                                     !!!
MODULE POLSYS2

! This module contains the subroutines POLSYS_GLP (finds all or some of
! the roots of a polynomial system defined in the module GLOBAL_GLP),
! BEZOUT_GLP (computes the generalized Bezout number), and SINGSYS_GLP
! (checks the nonsingularity of a generic start point).  Typically a
! user would only call POLSYS_GLP, and thus include in their main
! program the statements:
!   USE GLOBAL_GLP
!   USE POLSYS2, ONLY: POLSYS_GLP
! An expert user might want to call BEZOUT_GLP or SINGSYS_GLP
! separately, and thus these routines are also provided as module
! procedures.

USE GLOBAL_GLP
INCLUDE 'mpif.h'

CONTAINS
                                                                     !!!
SUBROUTINE POLSYS_GLP(INDEX_PATH_TRACKED,PATH_COUNT,N,TRACKTOL,FINALTOL,&
           SINGTOL,SSPAR,BGLP,IFLAG1,IFLAG2,ARCLEN,LAMBDA,ROOTS,NFE,&
           SCALE_FACTORS,NUMRR,RECALL,NO_SCALING,USER_F_DF)

! Using a probability-one globally convergent homotopy method,
! POLSYS_GLP finds all finite isolated complex solutions to a system
! F(X) = 0 of N polynomial equations in N unknowns with complex
! coefficients.  A generalized linear product (GLP) formulation is used
! for the start system of the homotopy map.
!
! POLSYS_GLP uses the module GLOBAL_GLP, which contains the definition
! of the polynomial system to be solved, and also defines the notation
! used below.  The user may also find it beneficial at some point to
! refer to the documentation for STEPNX in the HOMPACK90 package.
!
! The representation of F(X) is stored in the module GLOBAL_GLP.  Using
! the same notation as GLOBAL_GLP, F(X) is defined mathematically by
!
! F_I(X)=SUM_{J=1}^{NUMT(I)} C(I,J) * PRODUCT_{K=1}^N X(K)**D(I,J,K),
!
! for I=1,...,N.
! 
! POLSYS_GLP features target system scaling, a projective
! transformation so that the homotopy zero curves are tracked in complex
! projective space, and a generalized linear product (GLP) formulation of
! the start system.  Scaling may be disabled by the optional argument
! NO_SCALING.  Whatever the case, the roots of F(X) are always returned
! unscaled and untransformed.  The GLP set structure, possibly different 
! for each component F_I(X), is defined in the module GLOBAL_GLP.
!
! Scaling is carried out in the internal subroutine SCALE_GLP, and is
! an independent preprocessing step.  SCALE_GLP modifies the polynomial
! coefficients and creates and stores unscaling factors SCALE_FACTORS
! for the variables X(I).  The problem is solved with the scaled
! coefficients and scaled variables.  The coefficients of the target
! polynomial system, which are contained in the global structure
! POLYNOMIAL, remain in modified form on return from POLSYS_GLP.
!
! With the projective transformation, the system is essentially recast in
! homogeneous coordinates, Z(1),...,Z(N+1), and solved in complex
! projective space.  The resulting solutions are untransformed via
! X(I) = Z(I)/Z(N+1), I=1,...N, unless this division would cause
! overflow, in which case Re(X(I)) = Im(X(I)) = HUGE(1.0_R8).
! On return, for the Jth path, ROOTS(I,J) = X(I) for I=1,...,N, and
! ROOTS(N+1,J) = Z(N+1), the homogeneous variable.
!
! In the GLP scheme the number of paths that must be tracked can be
! less, and commonly far less, than the "total degree" because of the
! specialized start system G(X) = 0.  The structure of the start system
! is determined by the system set structure P.  The representations of both
! are stored in the module GLOBAL_GLP, and following the comments there,
! are defined mathematically as follows:
!
! The system set structure P=(P(1),...,P(N)) is comprised of the 
! coverings P(I) = {S(I,1),...S(I, COVER_SIZES(I))}, where the sets
! of variables S(I,J) are defined by
!
! S(I,J) = UNION_{K=1}^{NUMV(I,J)} {X(PAR(I,J,K))}.
!
! The degree of each set is provided in the input data file.
! The only restriction on the system set structure P is that each monomial
! of the target system must be in the span of the set structure of the start
! system.  CHECK_GLP returns an error if the start system does not contain a
! a monomial that is in the target system.
! 
! The start system is defined mathematically, for I=1,...,N, by
! 
! G_I(X) = PRODUCT_{J=1}^{COVER_SIZES(I)} ( L(I,J)**SD(I,J)-1.0 ),
!
! where the linear factors L(I,J) are
!
! L(I,J) = SUM{K=1}^{NUMV(I,J)} SC(I,J,K)*X(PAR(I,J,K)).
!
! Contained in this module (POLSYS2) is the routine BEZOUT_GLP.  This
! routine calculates the generalized GLP Bezout number, based on the
! system set structure P and SET_DEG provided by the user, by counting the 
! number of solutions to the start system.  The user is encouraged to explore 
! several system set structures with BEZOUT_GLP before calling POLSYS_GLP.  
! See the sample calling program MAIN_TEMPLATE and the comments in 
! BEZOUT_GLP.
!
! Internal routines:  INIT_GLP, INTERP, OUTPUT_GLP, RHO, ROOT_OF_UNITY,
!   ROOT_GLP, SCALE_GLP, START_POINTS_GLP, START_SYSTEM, TANGENT_GLP,
!   TARGET_SYSTEM. 
!
! External routines called:  CHECK_GLP, BEZOUT_GLP, SINGSYS_GLP, STEPNX.
!
!
! On input:
!
! N is the dimension of the target polynomial system.
!
! TRACKTOL is the local error tolerance allowed the path tracker along
!   the path.  ABSERR and RELERR (of STEPNX) are set to TRACKTOL.
!
! FINALTOL is the accuracy desired for the final solution.  It is used
!   for both the absolute and relative errors in a mixed error criterion.
!
! SINGTOL is the singularity test threshold used by SINGSYS_GLP.  If
!   SINGTOL <= 0.0 on input, then SINGTOL is reset to a default value.
!
! SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P) is a
!   vector of parameters used for the optimal step size estimation. If
!   SSPAR(I) <= 0.0 on input, it is reset to a default value by STEPNX. 
!   See the comments in STEPNX for more information.
!
! Optional arguments:
!
! NUMRR is the number of multiples of 1000 steps that will be tried before
!   abandoning a path.  If absent, NUMRR is taken as 1.
!
! RECALL is used to retrack certain homotopy paths.  It's use assumes
!   BGLP contains the Bezout number (which is not recalculated),
!   SCALE_FACTORS contains the variable unscaling factors, and that
!   IFLAG2(1:BGLP) exists.  The Ith homotopy path is retracked if
!   IFLAG2(I) = -2, and skipped otherwise.
!
! NO_SCALING indicates that the target polynomial is not to be scaled.
!   Scaling is done by default when NO_SCALING is absent.
!
! USER_F_DF indicates (when present) that the user is providing a subroutine
!   TARGET_SYSTEM_USER to evaluate the (complex) target system F(XC) and
!   its (complex) N x N Jacobian matrix DF(XC).  XC(1:N+1) is in
!   complex projective coordinates, and the homogeneous coordinate XC(N+1)
!   is explicitly eliminated from F(XC) and DF(XC) using the projective
!   transformation (cf. the comments in START_POINTS_GLP).
!
!
! The following objects must be allocated and defined as described in
! GLOBAL_GLP:
! 
! POLYNOMIAL(I)%NUM_TERMS is the number of terms in the Ith component
!   F_I(X) of the target polynomial system, for I=1,...,N.
!
! POLYNOMIAL(I)%TERM(J)%COEF is the coefficient of the Jth term in the Ith
!   component of the target polynomial system, for J=1,...,NUMT(I), and
!   I=1,...,N.
!
! POLYNOMIAL(I)%TERM(J)%DEG(K) is the degree of the Kth variable in the
!   Jth term of the Ith component of the target polynomial system, for
!   K=1,...,N, J=1,...NUMT(I), and I=1,...,N.
!
! COVER_SIZES(I) is the number of sets in the Ith component
!   covering P(I), for I=1,...,N.
!
! COVER(I)%SET(J)%NUM_INDICES is the number of indices stored in the
!   Jth set S(I,J) of the Ith component covering P(I), for
!   J=1,...,COVER_SIZES(I), and I=1,...,N.
!
! COVER(I)SET(J)%INDEX(K) is the index of the Kth variable stored
!   in the Jth set S(I,J) of the Ith component covering P(I).
!
!
! On output:
!
! 
! PATH_COUNT is the number of paths tracked by a processor.
!
! BGLP is the generalized Bezout number corresponding to the
!   generalized linear product (GLP) formulation defined by the system
!   set structure P.
!
! IFLAG1
!   = 0  for a normal return.
!
!   = -1 if either POLYNOMIAL or COVER was improperly allocated.
!
!   = -2 if any POLYNOMIAL(I)%TERM(J)%DEG(K) is less than zero.
!
!   = -3 if F_I(X) = CONSTANT for some I.
!
!   = -4 if SUM_{J=1}^{COVER_SIZES(I)} 
!          COVER(I)SET(J)%NUM_INDICES < N, for some I.
! 
!   = -5 if UNION_{J=1}^{COVER_SIZES}
!          S(I,J) /= {1,2,...,N-1,N}, for some I.
!
!   = -6 if the optional argument RECALL was present but any of BGLP
!        or the arrays ARCLEN, IFLAG2, LAMBDA, NFE, ROOTS are
!        inconsistent with the previous call to POLSYS_GLP.
!
!   = -7 if the array SCALE_FACTORS is too small.
!
!   = -8 if the input GLP set structure is invalid (inconsistent with
!        the set structure of the given polynomial).
!
! IFLAG2(1:PATH_COUNT) is an integer array which returns information about
!   each path tracked by a processor.  Precisely, for each path I that was 
!   tracked, IFLAG2(I):
!   = 1 + 10*C, where C is the cycle number of the path, for a normal return.
!
!   = 2 if the specified error tolerance could not be met.  Increase
!       TRACKTOL and rerun.
!
!   = 3 if the maximum number of steps allowed was exceeded.  To track
!       the path further, increase NUMRR and rerun the path.
!
!   = 4 if the Jacobian matrix does not have full rank.  The algorithm has
!       failed (the zero curve of the homotopy map cannot be followed any
!       further).
!
!   = 5 if the tracking algorithm has lost the zero curve of the homotopy
!       map and is not making progress.  The error tolerances TRACKTOL and
!       FINALTOL were too lenient.  The problem should be restarted with
!       smaller error tolerances.
!
!   = 6 if the normal flow Newton iteration in STEPNX or ROOT_GLP failed
!       to converge.  The error error tolerances TRACKTOL or FINALTOL may
!       be too stringent.
!
!   = 7 if ROOT_GLP failed to find a root in 10*NUMRR iterations.
!
! ARCLEN( I) is the approximate arc length of the Ith path, for I=1,...,BGLP.
!
! LAMBDA(I), if MOD(IFLAG2(I),10) = 1, contains an error estimate of
!   the normalized residual of the scaled, transformed polynomial
!   system of equations at the scaled, transformed root for the Ith path
!   (LAMBDA for this path is assumed to be 1). Otherwise LAMBDA(I) is the
!   final value of the homotopy parameter lambda on the Ith path, for
!   I=1,...,PATH_COUNT.
!
! ROOTS(1:N,I) are the complex roots (untransformed and unscaled) of
!   the target polynomial corresonding to the Ith path, for I=1,...,BGLP.
!
! ROOTS(N+1,I) is the homogeneous variable of the target polynomial
!   system in complex projective space corresponding to ROOTS(1:N,I).
!
! NFE(I) is the number of Jacobian matrix evaluations required to track
!   the Ith path, for I=1,...,PATH_COUNT.
!
! SCALE_FACTORS(1:N) contains the unscaling factors for the variables X(I).
!   These are needed only on a recall when scaling was done on the original
!   call to POLSYS_GLP (NO_SCALING was absent).


USE GLOBAL_GLP
IMPLICIT NONE

! Processor path variables.
INTEGER, INTENT(OUT):: PATH_COUNT

! POLSYS_GLP variables.
INTEGER, INTENT(IN):: N
REAL (KIND=R8), INTENT(IN):: TRACKTOL, FINALTOL
REAL (KIND=R8), INTENT(IN OUT):: SINGTOL
REAL (KIND=R8), DIMENSION(8), INTENT(IN OUT):: SSPAR
INTEGER, INTENT(IN OUT):: BGLP, IFLAG1
INTEGER, DIMENSION(:), POINTER:: IFLAG2
REAL (KIND=R8), DIMENSION(:), POINTER:: ARCLEN, LAMBDA
COMPLEX (KIND=R8), DIMENSION(:,:), POINTER:: ROOTS
INTEGER, DIMENSION(:), POINTER:: NFE
REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: SCALE_FACTORS
INTEGER, OPTIONAL, INTENT(IN):: NUMRR
LOGICAL, OPTIONAL, INTENT(IN):: RECALL, NO_SCALING, USER_F_DF

INTEGER, DIMENSION(:), POINTER:: INDEX_PATH_TRACKED

INTERFACE
  SUBROUTINE STEPNX(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,   &
      ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)
    USE REAL_PRECISION
    INTEGER, INTENT(IN):: N
    INTEGER, INTENT(IN OUT):: NFE,IFLAG
    LOGICAL, INTENT(IN OUT):: START,CRASH
    REAL (KIND=R8), INTENT(IN OUT):: HOLD,H,RELERR,ABSERR,S,RHOLEN, &
      SSPAR(8)
    REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
    REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD, &
      TZ,W,WP
    REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: Z0,Z1
  END SUBROUTINE STEPNX
END INTERFACE

! Local variables.
INTEGER:: BTEMP, I, IFLAG, II, ITER, J, JJ, K, KK, LIMIT, MAXPS, &
  MAXT, NNFE, NUM_RERUNS
INTEGER, SAVE:: BGLP_SAVE
INTEGER, DIMENSION(N):: CHECK_PAR, DLEX_NUM, DLEX_SAVE, FLEX_NUM, FLEX_SAVE
INTEGER, DIMENSION(2*N+1):: PIVOT
REAL (KIND=R8):: ABSERR, H, HOLD, RELERR, RHOLEN, S
REAL (KIND=R8), DIMENSION(2*N):: A, DRHOL, RHOV, Z
REAL (KIND=R8), DIMENSION(2*N+1):: Y, YP, YOLD, YOLDS, YPOLD, TZ, W, WP
REAL (KIND=R8), DIMENSION(3*(2*N+1)):: ALPHA
REAL (KIND=R8), DIMENSION(2*N+1,12):: YS
REAL (KIND=R8), DIMENSION(N,N):: RAND_MAT
REAL, DIMENSION(N,N):: RANDNUMS
REAL (KIND=R8), DIMENSION(N+1,N):: MAT
REAL (KIND=R8), DIMENSION(2*N,2*N):: DRHOX
REAL (KIND=R8), DIMENSION(2*N,2*N+2):: QR
COMPLEX (KIND=R8), DIMENSION(N-1):: TAU
COMPLEX (KIND=R8), DIMENSION(N):: B, F, G, V
COMPLEX (KIND=R8), DIMENSION(N+1):: PROJ_COEF, XC
COMPLEX (KIND=R8), DIMENSION(N,N):: AA
COMPLEX (KIND=R8), DIMENSION(N,N+1):: DF, DG
COMPLEX (KIND=R8), DIMENSION(:,:), ALLOCATABLE:: TEMP1G, TEMP2G
LOGICAL:: CRASH, NONSING, START


! Parallel computing variables:
!
! MASTER_PROC is the processor with rank 0.
! NUM_PROC is the number of processors requested for the run.
! PATH_CYCLE is the current number of paths that have been tracked by
! a processor.
! PATH_TRACKING is a path index to be distributed by the master to a
! slave processor.
! RANK_PROC is an index that denotes a processor.

INTEGER, PARAMETER:: MASTER_PROC = 0 
INTEGER:: NUM_PROC, PATH_CYCLE, PATH_TRACKING, RANK_PROC
REAL (KIND=R8):: END_TIME, RUN_TIME, START_TIME

! MPI system variables:
!
! IERR is a return error code, not used.
! MPI_SOURCE is a variable that denotes a processor generating a message.
! SENDER is a variable indentifying the rank of a processor.
! STATUS is used to query processors via SENDER = STATUS(MPI_SOURCE).

INTEGER:: IERR, SENDER
INTEGER, DIMENSION(MPI_STATUS_SIZE):: STATUS

! Set MPI variables to initial values and get current processor rank.
CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK_PROC, IERR)

! Determine number of processors requested.
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUM_PROC, IERR)

! Begin input data check.

IFLAG1 = 0  ! Normal return.

! Check that dimensions are valid.
IF ((N <= 0) .OR. (SIZE(POLYNOMIAL) /= N)          & 
             .OR. ANY((/(NUMT(I),I=1,N)/) <= 0)    &
             .OR. (SIZE(COVER) /= N)           &
             .OR. ANY(COVER_SIZES <= 0)) THEN
  IFLAG1 = -1
  RETURN
END IF
DO I=1,N
  IF ((SIZE(POLYNOMIAL(I)%TERM) /= NUMT(I))                   &
    .OR. (SIZE(COVER(I)%SET) /= COVER_SIZES(I))       &
    .OR. ANY((/(NUMV(I,J),J=1,COVER_SIZES(I))/) <= 0)) THEN  
    IFLAG1 = -1
    RETURN
  END IF
END DO
DO I=1,N
  DO J=1,NUMT(I)
    IF (SIZE(POLYNOMIAL(I)%TERM(J)%DEG) /= N + 1) THEN
      IFLAG1 = -1
      RETURN
    END IF
  END DO
  DO J=1,COVER_SIZES(I)
    IF (SIZE(COVER(I)%SET(J)%INDEX) /= NUMV(I,J)) THEN
      IFLAG1 = -1
      RETURN
    END IF
  END DO
END DO

! Check that the target system has no negative powers.
DO I=1,N
  DO J=1,NUMT(I)
    IF (ANY(POLYNOMIAL(I)%TERM(J)%DEG(1:N) < 0)) THEN
      IFLAG1 = -2
      RETURN
    END IF
  END DO
END DO

! Check that the target system has no constant-valued components.
DO I=1,N
  IF (ALL( (/( SUM(POLYNOMIAL(I)%TERM(J)%DEG(1:N)),   & 
                                 J=1,NUMT(I) )/) == 0)) THEN
    IFLAG1 = -3
    RETURN
  END IF
END DO

DO I=1,N
  IF (SUM( (/(NUMV(I,J),J=1,COVER_SIZES(I))/) ) < N) THEN
    IFLAG1 = -4
    RETURN
  END IF
  CHECK_PAR(1:N) = 0
  DO J=1,COVER_SIZES(I)
    DO K=1,NUMV(I,J)
      CHECK_PAR(PAR(I,J,K)) = CHECK_PAR(PAR(I,J,K)) + 1
    END DO
  END DO
  IF (ANY(CHECK_PAR < 1)) THEN
    IFLAG1 = -5
    RETURN
  END IF
END DO


! Check consistency on a recall.
IF (PRESENT(RECALL)) THEN
  IF ( (BGLP /= BGLP_SAVE) .OR. (SIZE(ARCLEN) < BGLP)              &
                           .OR. (SIZE(IFLAG2) < BGLP)              &
                           .OR. (SIZE(LAMBDA) < BGLP)              &
                           .OR. (SIZE(NFE) < BGLP)                 &
                           .OR. (SIZE(ROOTS,DIM=2) < BGLP) ) THEN
    IFLAG1 = -6
    RETURN
  END IF
END IF

! Check SCALE_FACTORS array size.
IF (SIZE(SCALE_FACTORS) < N) THEN
  IFLAG1 = -7
  RETURN
END IF

! End input data check.

! Initialize the POINTER aguments of POLSYS_GLP.
MAXT = MAXVAL((/(NUMT(I),I=1,N)/))
IF ( .NOT. PRESENT(RECALL)) THEN
  
  ! Return error if GLP set structure not valid.
  CALL BEZOUT_GLP(N,MAXT,SINGTOL,BGLP)
  IF (BGLP < 0 ) THEN
    IFLAG1 = -8 ! Invalid GLP covering.
    RETURN
  END IF

  BGLP_SAVE = BGLP  ! Save Bezout number for recall check.
  
  IF (ASSOCIATED(ARCLEN)) THEN
    IF (SIZE(ARCLEN) < BGLP) THEN
      DEALLOCATE(ARCLEN) ; ALLOCATE(ARCLEN(BGLP))
    END IF
  ELSE
    ALLOCATE(ARCLEN(BGLP))
  END IF
  IF (ASSOCIATED(IFLAG2)) THEN
    IF (SIZE(IFLAG2) < BGLP) THEN
      DEALLOCATE(IFLAG2) ; ALLOCATE(IFLAG2(BGLP))
    END IF
  ELSE
    ALLOCATE(IFLAG2(BGLP))
  END IF
  
  IF(RANK_PROC==MASTER_PROC) THEN
    IFLAG2 = -2
  ELSE
    IFLAG2 = 0
  END IF
  
  IF (ASSOCIATED(NFE)) THEN
    IF (SIZE(NFE) < BGLP) THEN
      DEALLOCATE(NFE) ; ALLOCATE(NFE(BGLP))
    END IF
  ELSE
    ALLOCATE(NFE(BGLP))
  END IF
  IF (ASSOCIATED(LAMBDA)) THEN
    IF (SIZE(LAMBDA) < BGLP) THEN
      DEALLOCATE(LAMBDA) ; ALLOCATE(LAMBDA(BGLP))
    END IF
  ELSE
    ALLOCATE(LAMBDA(BGLP))
  END IF
  IF (ASSOCIATED(ROOTS)) THEN
      DEALLOCATE(ROOTS) ; ALLOCATE(ROOTS(N+1,BGLP))
  ELSE
    ALLOCATE(ROOTS(N+1,BGLP))
  END IF
  IF (ASSOCIATED(INDEX_PATH_TRACKED)) THEN
    IF (SIZE(INDEX_PATH_TRACKED) < BGLP) THEN
      DEALLOCATE(INDEX_PATH_TRACKED); ALLOCATE(INDEX_PATH_TRACKED(BGLP))
    ENDIF
  ELSE
   ALLOCATE(INDEX_PATH_TRACKED(BGLP))
  END IF
END IF

! Allocate storage for the start system.
DO I=1,N
  DO J=1,COVER_SIZES(I)
    ALLOCATE(COVER(I)%SET(J)%START_COEF(NUMV(I,J))) 
  END DO
END DO

! Allocate working space for homotopy map derivative calculation.
MAXPS = MAXVAL(COVER_SIZES)                  
ALLOCATE(TEMP1G(N,MAXPS), TEMP2G(N,MAXPS))        

! Get real random numbers uniformly distributed in [-1,-1/2] union
! [1/2,1] for RAND_MAT, which is used in SINGSYS_GLP. 
CALL RANDOM_NUMBER(HARVEST=RANDNUMS)
RANDNUMS = RANDNUMS - 0.5 + SIGN(0.5,RANDNUMS - 0.5)
RAND_MAT = REAL(RANDNUMS,KIND=R8)

! Set default value for singularity threshold SINGTOL in SINGSYS_GLP.
IF (SINGTOL <= REAL(N,KIND=R8)*EPSILON(1.0_R8))  &
    SINGTOL = SQRT(EPSILON(1.0_R8))

! Scale the target polynomial system as requested.
IF (PRESENT(NO_SCALING)) THEN
  SCALE_FACTORS = 0.0_R8
ELSE IF (.NOT. PRESENT(RECALL)) THEN
  CALL SCALE_GLP
END IF

! Initialize the start system for the homotopy map.
CALL INIT_GLP

! Set main loop initial values.
FLEX_NUM(1:N-1) = 1
FLEX_NUM(N) = 0
FLEX_SAVE = 0
IF (PRESENT(NUMRR)) THEN
  NUM_RERUNS = MAX(NUMRR,1)
ELSE
  NUM_RERUNS = 1
END IF

PATH_COUNT = 0
PATH_CYCLE = 0
PATH_TRACKING = 0
INDEX_PATH_TRACKED = 0

IF (NUM_PROC > 1 .AND. RANK_PROC==MASTER_PROC) THEN
  ! Executed by master processor.
  START_TIME = MPI_WTIME()
  J = MAX(1, (BGLP)/100)
  K = 1
  JJ = 1
  DO I=1,BGLP
    IF(IFLAG2(I) /= -2) CYCLE
    ! Receive message from slave processor.
    CALL MPI_RECV(PATH_TRACKING,1,MPI_INTEGER,MPI_ANY_SOURCE,&
       MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
    ! Get the rank of the sender.
    SENDER = STATUS(MPI_SOURCE)

    PATH_TRACKING = I
    ! Send next path index.
    CALL MPI_SEND(PATH_TRACKING,1,MPI_INTEGER,SENDER,SENDER,&
       MPI_COMM_WORLD,IERR)

    ! Print to standard out for each 1% of total paths completed by
    ! uncommenting the WRITE statements below.
    IF (K < J ) THEN 
      K = K + 1
    ELSE 
      K = 1
      END_TIME = MPI_WTIME()
      RUN_TIME = END_TIME - START_TIME

    ! WRITE (*,100) NINT(100.0*J*JJ/BGLP), J*JJ, BGLP, &
    !    RUN_TIME*(BGLP-I)/I
    ! 100 FORMAT(I3, "% or", I6, " out of ", I6, " paths completed.  ", &
    !            ES12.4, " seconds left.")
      JJ = JJ + 1
    END IF
  END DO
  
  ! Send end of job message to all slave processors.
  DO I=1,NUM_PROC-1
    CALL MPI_RECV(PATH_TRACKING,1,MPI_INTEGER,MPI_ANY_SOURCE,&
       MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
    SENDER = STATUS(MPI_SOURCE)
    CALL MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER,SENDER,0,&
       MPI_COMM_WORLD,IERR)
  END DO

ELSE ! Executed by slave processors.

! If only one processor is specified, then parallel code is not executed.
IF (NUM_PROC > 1) THEN
  ! Request a path from master processor.
  CALL MPI_SEND(PATH_TRACKING,1,MPI_INTEGER,MASTER_PROC,MASTER_PROC,&
                MPI_COMM_WORLD,IERR)
  ! Receive a path from master processor.
  CALL MPI_RECV(PATH_TRACKING,1,MPI_INTEGER,MASTER_PROC,&
                MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
  IF (STATUS(MPI_TAG) /= 0) THEN
    IFLAG2(PATH_TRACKING) = -2
  ELSE
    GO TO 1000
  END IF
END IF

! Main loop over all possible lexicographic vectors FLEX_NUM(1:N)
! corresponding to linear factors.

MAIN_LOOP: &
DO 

  DO J=N,1,-1
    IF (FLEX_NUM(J) < COVER_SIZES(J)) THEN
      K = J
      EXIT      
    END IF
  END DO 
  FLEX_NUM(K) = FLEX_NUM(K) + 1
  IF (K + 1 <= N) THEN
    FLEX_NUM(K+1:N) = 1
  END IF
  ! Check if the subsystem of the start system defined by the
  ! lexicographic vector FLEX_NUM is singular. 
  CALL SINGSYS_GLP(N,FLEX_NUM,FLEX_SAVE,SINGTOL,RAND_MAT,MAT,NONSING)

  ! If the subsystem is nonsingular, track a path.
  NONSING_START_POINT: IF (NONSING) THEN
    BTEMP = PRODUCT( (/(SD(I,FLEX_NUM(I)),I=1,N)/) )
    DLEX_NUM(1:N-1) = 1
    DLEX_NUM(N) = 0
    DLEX_SAVE = 0
  
    ! Cycle through all lexicographic vectors DLEX_NUM(1:N) corresponding
    ! to roots of unity, defined by the set degrees specified in
    ! (/(SD(I,FLEX_NUM(I)),I=1,N)/).
    SD_LEX_LOOP: DO II=1,BTEMP
      DO JJ=N,1,-1
        IF (DLEX_NUM(JJ) < SD(JJ,FLEX_NUM(JJ))) THEN
          KK = JJ
          EXIT
        END IF
      END DO
      DLEX_NUM(KK) = DLEX_NUM(KK) + 1
      IF (KK + 1 <= N) THEN
         DLEX_NUM(KK+1:N) = 1
      END IF
       
      PATH_CYCLE = PATH_CYCLE + 1
      
      IF (IFLAG2(PATH_CYCLE) /= -2) CYCLE SD_LEX_LOOP
      
      PATH_COUNT = PATH_COUNT + 1
      
      INDEX_PATH_TRACKED(PATH_COUNT) = PATH_CYCLE
      
      ! Get the start point for the homotopy path defined by FLEX_NUM and
      ! DLEX_NUM.
      CALL START_POINTS_GLP
  
      NNFE = 0
      IFLAG = -2
      Y(1) = 0.0_R8 ; Y(2:2*N+1) = Z(1:2*N)
      YP(1) = 1.0_R8 ;  YP(2:2*N+1) = 0.0_R8
      YOLD = Y ;  YPOLD = YP
      HOLD = 1.0_R8 ; H = 0.1_R8
      S = 0.0_R8
      LIMIT = 1000*NUM_RERUNS
      START = .TRUE.
      CRASH = .FALSE.
  
      ! Track the homotopy path.
  
      TRACKER: DO ITER=1,LIMIT
        IF (Y(1) < 0.0_R8) THEN
          IFLAG = 5
          EXIT TRACKER
        END IF
  
        ! Set different error tolerance if the trajectory Y(S) has any high
        ! curvature components.
        RELERR = TRACKTOL
        ABSERR = TRACKTOL
        IF (ANY(ABS(YP - YPOLD) > 10.0_R8*HOLD)) THEN
          RELERR = FINALTOL
          ABSERR = FINALTOL
        END IF
  
        ! Take a step along the homotopy zero curve.
        CALL STEP_GLP
        IF (IFLAG > 0) THEN
          EXIT TRACKER
        END IF
        IF (Y(1) >= .97_R8) THEN
          RELERR = FINALTOL
          ABSERR = FINALTOL
  
          ! Enter end game.
          CALL ROOT_GLP
          EXIT TRACKER
        END IF
  
        ! D LAMBDA/DS >= 0 necessarily.  This condition is forced here.
        IF (YP(1) < 0.0_R8) THEN
  
          ! Reverse the tangent direction so D LAMBDA/DS = YP(1) > 0.
          YP = -YP
          YPOLD = YP
  
          ! Force STEPNX to use the linear predictor for the next step only.
          START = .TRUE.
        END IF
      END DO TRACKER
  
      ! Set error flag if limit on number of steps exceeded.
      IF (ITER >= LIMIT) IFLAG = 3
  
      ARCLEN(PATH_COUNT) = S
      NFE(PATH_COUNT) = NNFE
      IFLAG2(PATH_COUNT) = IFLAG
      LAMBDA(PATH_COUNT) = Y(1)
  
      ! Convert from real to complex arithmetic.
      XC(1:N) = CMPLX(Y(2:2*N:2),Y(3:2*N+1:2),KIND=R8)
  
      ! Untransform and unscale solutions.
      CALL OUTPUT_GLP
      ROOTS(1:N,PATH_COUNT) = XC(1:N)
      ROOTS(N+1,PATH_COUNT) = XC(N+1)
      
      ! One path tracking finished, request another one.
      ! If only one processor is specified, then parallel code is not executed.
      IF (NUM_PROC > 1 ) THEN 
        ! Request a path from master processor.
        CALL MPI_SEND(PATH_TRACKING,1,MPI_INTEGER,MASTER_PROC,MASTER_PROC,&
                      MPI_COMM_WORLD,IERR)
        ! Receive a path from master processor.
        CALL MPI_RECV(PATH_TRACKING,1,MPI_INTEGER,MASTER_PROC,&
                      MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
        IF (STATUS(MPI_TAG) /= 0)  THEN 
          IFLAG2(PATH_TRACKING) = -2
        ELSE ! No more jobs.
          EXIT MAIN_LOOP
        END IF
      END IF
      
    END DO SD_LEX_LOOP
  END IF NONSING_START_POINT
  
  IF (ALL(FLEX_NUM == COVER_SIZES)) THEN
     EXIT MAIN_LOOP
  END IF      

END DO MAIN_LOOP

END IF 


! Clean up working storage in STEPNX.
1000 IFLAG = -42
CALL STEPNX (2*N,NNFE,IFLAG,START,CRASH,HOLD,H,RELERR, &
  ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)

! Deallocate the storage for the start system and working storage.
DO I=1,N
  DO J=1,COVER_SIZES(I)
    DEALLOCATE(COVER(I)%SET(J)%START_COEF)
  END DO
END DO
DEALLOCATE(TEMP1G,TEMP2G)
RETURN

CONTAINS

                                                                     !!!
SUBROUTINE SCALE_GLP
! SCALE_GLP scales the complex coefficients of a polynomial system of N
! equations in N unknowns, F(X)=0, where the Jth term of the Ith equation
! looks like 
!
!    C(I,J) * X(1)**D(I,J,1) ... X(N)**D(I,J,N).
!
! The Ith equation is scaled by 10**FACE(I).  The Kth variable is scaled
! by 10**FACV(K).  In other words, X(K)=10**FACV(K)*Y(K), where Y solves
! the scaled equation.  The scaled equation has the same form as the
! original, except that CSCL(I,J) replaces POLYNOMIAL(I)%TERM(J)%COEF,
! where 
!
! CSCL(I,J)=C(I,J)*10**(FACE(I)+FACV(1)*D(I,J,1)+...+FACV(N)*D(I,J,N)).
!
! The criterion for generating FACE and FACV is that of minimizing the
! sum of squares of the exponents of the scaled coefficients.  It turns
! out that this criterion reduces to solving a single linear system,
! ALPHA*X=BETA, as defined in the code below.  See Meintjas and Morgan,
! "A methodology for solving chemical equilibrium problems," General
! Motors Research Laboratories Technical Report GMR-4971.
!
! Calls the LAPACK routines DGEQRF, DORMQR, and the BLAS routines
! DTRSV and IDAMAX.
!
! On exit:
!
! SCALE_FACTORS(K) = FACV(K) is the scale factor for X(K), K=1,...,N.
!  Precisely, the unscaled solution
!  X(K) = 10**FACV(K) * (computed scaled solution).
!
! POLYNOMIAL(I)%TERM(J)%COEF = CSCL(I,J) is the scaled complex
!  coefficient, for J=1,...,NUMT(I), and I=1,...,N.

! Local variables.
IMPLICIT NONE
INTEGER:: COUNT, I, ICMAX, IRMAX, J, K, L, LENR
INTEGER, DIMENSION(N):: NNUMT
INTEGER, DIMENSION(N,MAXT,N):: DDEG 
REAL (KIND=R8):: DUM, RTOL, TUM
REAL (KIND=R8), DIMENSION(:), POINTER:: FACE, FACV
REAL (KIND=R8), DIMENSION(2*N), TARGET:: BETA, RWORK, XWORK
REAL (KIND=R8), DIMENSION(2*N,2*N):: ALPHA
REAL (KIND=R8), DIMENSION(N,MAXT):: CMAG

INTERFACE
  INTEGER FUNCTION IDAMAX(N,X,STRIDE)
    USE REAL_PRECISION
    INTEGER:: N,STRIDE
    REAL (KIND=R8), DIMENSION(N):: X
  END FUNCTION IDAMAX
END INTERFACE

LENR = N*(N+1)/2
SCALE_FACTORS(1:N) = 0.0_R8 ! This corresponds to no scaling.

! Delete exact zero coefficients, just for scaling.
NNUMT = 0
DO I=1,N
  COUNT = 0
  DO J=1,NUMT(I)
    IF (ABS(C(I,J)) > 0.0_R8) THEN
      COUNT = COUNT + 1
      NNUMT(I) = NNUMT(I) + 1
      CMAG(I,COUNT) = LOG10(ABS(C(I,J)))
      DDEG(I,COUNT,1:N) = (/(D(I,J,K),K=1,N)/)
    END IF
  END DO
END DO

! Generate the matrix ALPHA.
ALPHA(1:N,1:N) = 0.0_R8
DO I=1,N
  ALPHA(I,I) = REAL(NNUMT(I),KIND=R8)
END DO
DO I=1,N
  ALPHA(N+1:2*N,I) = REAL(SUM(DDEG(I,1:NNUMT(I),1:N),DIM=1),KIND=R8)
END DO
DO L=1,N
  DO K=1,L
    ICMAX = 0
    DO I=1,N
      ICMAX = ICMAX + DOT_PRODUCT(DDEG(I,1:NNUMT(I),L),DDEG(I,1:NNUMT(I),K))
    END DO
    ALPHA(N+L,N+K) = REAL(ICMAX,KIND=R8)
    ALPHA(N+K,N+L) = ALPHA(N+L,N+K)
  END DO
END DO
ALPHA(1:N,N+1:2*N) = TRANSPOSE(ALPHA(N+1:2*N,1:N))

! Compute the QR-factorization of the matrix ALPHA.
CALL DGEQRF(2*N,2*N,ALPHA,2*N,XWORK,BETA,2*N,I)

! Check for ill-conditioned scaling matrix.
IRMAX = 1
ICMAX = 1
DO J=2,N
  I = IDAMAX(J,ALPHA(1,J),1)
  IF (ABS(ALPHA(I,J)) > ABS(ALPHA(IRMAX,ICMAX))) THEN 
    IRMAX = I
    ICMAX = J
  END IF
END DO
RTOL = ABS(ALPHA(IRMAX,ICMAX))*EPSILON(1.0_R8)*REAL(N,KIND=R8)
DO I=1,N
  IF (ABS(ALPHA(I,I)) < RTOL) THEN  ! ALPHA is ill conditioned.
    RETURN  ! Default to no scaling at all.
  END IF
END DO

! Generate the column BETA.
DO K=1,N
  BETA(K) = -SUM(CMAG(K,1:NNUMT(K)))
  TUM = 0.0_R8
  DO I=1,N
    TUM = TUM + SUM(CMAG(I,1:NNUMT(I)) * REAL(DDEG(I,1:NNUMT(I),K),KIND=R8))
  END DO
  BETA(N+K) = -TUM
END DO

! Solve the linear system ALPHA*X=BETA.
CALL DORMQR('L','T',2*N,1,2*N-1,ALPHA,2*N,XWORK,BETA,2*N,RWORK,2*N,I)  
CALL DTRSV('U','N','N',2*N,ALPHA,2*N,BETA,1) 

! Generate FACE, FACV, and the scaled coefficients CSCL(I,J).
FACE => BETA(1:N)
FACV => BETA(N+1:2*N)
DO I=1,N
  DO J=1,NUMT(I)
    DUM = ABS(C(I,J))
    IF (DUM /= 0.0) THEN
      TUM = FACE(I) + LOG10(DUM) + DOT_PRODUCT(FACV(1:N), &
                                   POLYNOMIAL(I)%TERM(J)%DEG(1:N))
      POLYNOMIAL(I)%TERM(J)%COEF = (10.0_R8**TUM) * (C(I,J)/DUM)
    ENDIF
  END DO
END DO

SCALE_FACTORS(1:N) = FACV(1:N)
RETURN
END SUBROUTINE SCALE_GLP

                                                                     !!!
SUBROUTINE INIT_GLP
! INIT_GLP homogenizes the homotopy map, and harvests random complex
! numbers which define the start system and the projective transformation.
! 
!
! On exit:
!
! POLYNOMIAL(I)%TERM(J)%DEG(N+1) is the degree of the homogeneous variable
!  in the Jth term of the Ith component of the target system.
!
! COVER(I)%SET(J)%START_COEF(K) is the coefficient of X(PAR(I,J,K)) in
!  the linear factor L(I,J).  (L(I,J) is defined in GLOBAL_GLP.)  
!
! PROJ_COEF(I) is the coefficient of X(I) in the projective transformation,
!  when I=1,...,N, and the constant term in the projective transformation,
!  when I=N+1. 
!

! Local variables.
IMPLICIT NONE
INTEGER:: COUNT, I, J, K, SEED_SIZE
INTEGER, DIMENSION(:), ALLOCATABLE:: SEED
REAL, DIMENSION(N*N+N+1,2):: RANDS
REAL (KIND=R8), DIMENSION(N*N+N+1,2):: RANDSR8

! Construct the homogenization of the homotopy map.  Note: 
! Homogenization of the start system is implicit.
DO I=1,N
  DO J=1,NUMT(I)
    POLYNOMIAL(I)%TERM(J)%DEG(N+1) = SUM((/(SD(I,K),K=1, &
      COVER_SIZES(I))/)) - SUM(POLYNOMIAL(I)%TERM(J)%DEG(1:N))
  END DO
END DO

! Get the random coefficients START_COEF which define the start system
! and the random coefficients PROJ_COEF which define the projective
! transformation.

!CALL RANDOM_SEED(SIZE=SEED_SIZE)  ! Explicit seed used by POLSYS_PLP;
!ALLOCATE(SEED(SEED_SIZE))         ! default system seed will be used here.
!SEED(1:SEED_SIZE) = 32749
!CALL RANDOM_SEED(PUT=SEED(1:SEED_SIZE))
CALL RANDOM_SEED
CALL RANDOM_NUMBER(HARVEST=RANDS)
!DEALLOCATE(SEED)

RANDS = 2.0 * RANDS - 1.0
RANDSR8 = REAL(RANDS,KIND=R8)
COUNT = 1
DO I=1,N
  DO J=1,COVER_SIZES(I)
    DO K=1,NUMV(I,J)
      COVER(I)%SET(J)%START_COEF(K) = CMPLX(RANDSR8(COUNT,1),  &
        RANDSR8(COUNT,2),KIND=R8)
      COUNT = COUNT + 1
    END DO
  END DO
END DO
PROJ_COEF(1:N+1) = CMPLX(RANDSR8(COUNT:COUNT+N,1),  &
  RANDSR8(COUNT:COUNT+N,2),KIND=R8)

RETURN
END SUBROUTINE INIT_GLP

                                                                     !!!
SUBROUTINE START_POINTS_GLP
! START_POINTS_GLP finds a starting point for the homotopy map
! corresponding to the lexicographic vector FLEX_NUM (defining the
! variable sets) and the lexicographic vector DLEX_NUM (defining the
! particular start point among all those defined by FLEX_NUM).  The
! (complex) start point z is the solution to a nonsingular linear system
! AA z = B, defined by (cf. the notation in the module GLOBAL_GLP)
!
! L(1,FLEX_NUM(1)) - R(DLEX_NUM(1)-1,SD(1,FLEX_NUM(1))) * X(N+1) = 0,
! .
! .
! .
! L(N,FLEX_NUM(N)) - R(DLEX_NUM(N)-1,SD(N,FLEX_NUM(N))) * X(N+1) = 0,
! X(N+1) = SUM_{J=1}^N PROJ_COEF(J)*X(J) + PROJ_COEF(N+1),
!
! where the last equation is the projective transformation, X(N+1) is
! the homogeneous coordinate, and R(K,M)=e**(i*2*PI*K/M) is an Mth root
! of unity.  The homogeneous variable X(N+1) is explicitly eliminated,
! resulting in an N x N complex linear system for z=(X(1),...,X(N)).
!
! START_POINTS_GLP calculates a start point in an efficient way:  For each
! fixed lexicographic number LEX_NUM, the routine reuses, if possible,
! previous Householder reflections in the LQ decomposition of AA.
!
! Calls the LAPACK routines ZLARFG, ZLARFX, the BLAS routine ZTRSV, and the
! internal function ROOT_OF_UNITY. 
!
! On exit:
!
! Z(1:2N) is a real vector representing the (complex) start point z.

! Local variables.
IMPLICIT NONE
INTEGER:: I, J, K
COMPLEX (KIND=R8):: ROOT, WORK(1)

! (Re)set the coefficient matrix AA, and set B.
DO I=1,N
  IF (DLEX_SAVE(I) /= DLEX_NUM(I)) THEN
    DLEX_SAVE(I+1:N) = 0
    DO J=1,N
      ROOT = ROOT_OF_UNITY(DLEX_NUM(J)-1,SD(J,FLEX_NUM(J)))
      B(J) = ROOT * PROJ_COEF(N+1)
      IF (J >= I) THEN
        AA(J,1:N) = (0.0_R8,0.0_R8)
        K = NUMV(J,FLEX_NUM(J))
        AA(J,COVER(J)%SET(FLEX_NUM(J))%INDEX(1:K)) = &
          COVER(J)%SET(FLEX_NUM(J))%START_COEF(1:K)
        AA(J,1:N) = AA(J,1:N) - PROJ_COEF(1:N) * ROOT
      END IF
    END DO
    EXIT
  END IF
END DO

! Special code for the case N=1.
IF (N == 1) THEN
  WORK(1) = B(1)/AA(1,1)
  Z(1) = REAL(WORK(1))
  Z(2) = AIMAG(WORK(1))
  DLEX_SAVE = DLEX_NUM
  RETURN
END IF

! Update the LQ factorization of AA.
IF (DLEX_SAVE(1) /= DLEX_NUM(1)) THEN
  AA(1,1:N) = CONJG(AA(1,1:N))
  CALL ZLARFG(N,AA(1,1),AA(1,2:N),1,TAU(1))
END IF
DO I=2,N
  IF (DLEX_SAVE(I) /= DLEX_NUM(I)) THEN
    DO J=1,I-1
      V(J) = (1.0_R8,0.0_R8)
      V(J+1:N) = AA(J,J+1:N)
      CALL ZLARFX('R',1,N-J+1,V(J:N),TAU(J),AA(I,J:N),1,WORK)
    END DO
    IF (I < N) THEN
      AA(I,I:N) = CONJG(AA(I,I:N))
      CALL ZLARFG(N-I+1,AA(I,I),AA(I,I+1:N),1,TAU(I))
    END IF
  END IF
END DO
DLEX_SAVE = DLEX_NUM

! Solve the linear system AA Z = B, by solving L Q Z = B.

! L W = B.
CALL ZTRSV('L','N','N',N,AA(1:N,1:N),N,B(1:N),1)
! Z = CONJG(Q') W.
DO I=N-1,1,-1
  V(I) = (1.0_R8,0.0_R8)
  V(I+1:N) = AA(I,I+1:N)
  CALL ZLARFX('L',N-I+1,1,V(I:N),TAU(I),B(I:N),N,WORK)
END DO

! Convert the complex start point to a real vector.
Z(1:2*N:2) = REAL(B)
Z(2:2*N:2) = AIMAG(B)
RETURN
END SUBROUTINE START_POINTS_GLP

                                                                     !!!
COMPLEX (KIND=R8) FUNCTION ROOT_OF_UNITY(K,N) RESULT(RU)
! RU = e**(i*2*PI*K/N).
  IMPLICIT NONE
  INTEGER:: K, N
  REAL (KIND=R8):: ANGLE
  ANGLE = 2.0_R8*PI*(REAL(K,KIND=R8)/REAL(N,KIND=R8))
  RU = CMPLX(COS(ANGLE),SIN(ANGLE),KIND=R8)
  RETURN
END FUNCTION ROOT_OF_UNITY

                                                                     !!!
SUBROUTINE STEP_GLP

! Driver for reverse call external subroutine STEPNX from HOMPACK90.

IMPLICIT NONE
INTEGER:: FAIL,IFLAGS
FAIL=0
STEP: DO
  CALL STEPNX(2*N,NNFE,IFLAG,START,CRASH,HOLD,H,RELERR, &
    ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)
  IF (CRASH) THEN
    IFLAG = 2
    EXIT
  END IF
  IFLAGS = IFLAG
  SELECT CASE (IFLAGS)
    CASE (-2)      ! Successful step.
      EXIT
    CASE (-12)     ! Compute tangent vector.
      RHOLEN = 0.0_R8
      CALL TANGENT_GLP
      IF (IFLAG == 4) THEN
        IFLAG = IFLAGS - 100
        FAIL = FAIL + 1
      ENDIF
    CASE (-32,-22) ! Compute tangent vector and Newton step.
      RHOLEN = -1.0_R8
      CALL TANGENT_GLP(NEWTON_STEP=.TRUE.)
      IF (IFLAG == 4) THEN
        IFLAG = IFLAGS - 100
        FAIL = FAIL + 1
      ENDIF
    CASE (4,6) ! STEPNX failed.
      EXIT
  END SELECT
  IF (FAIL == 2) THEN
    IFLAG = 4 ; RETURN
  ENDIF
END DO STEP
RETURN
END SUBROUTINE STEP_GLP

                                                                     !!!
SUBROUTINE TANGENT_GLP(NEWTON_STEP)
! This subroutine builds the Jacobian matrix of the homotopy map,
! computes a QR decomposition of that matrix, and then calculates the
! (unit) tangent vector and (if NEWTON_STEP is present) the Newton
! step.
!
! On input:
!
! NEWTON_STEP is a logical optional argument which, if present,
!    indicates that the Newton step is also to be calculated.
!
! RHOLEN < 0 if the norm of the homotopy map evaluated at
!    (LAMBDA, X) is to be computed.  If  RHOLEN >= 0  the norm is not
!    computed and RHOLEN is not changed.
!
! W(1:2*N+1) = current point (LAMBDA(S), X(S)).
!
! YPOLD(1:2*N+1) = unit tangent vector at previous point on the zero
!    curve of the homotopy map.
!
! On output:
!
! RHOLEN = ||RHO(LAMBDA(S), X(S))|| if  RHOLEN < 0  on input.
!    Otherwise RHOLEN is unchanged.
!
! WP(1:2*N+1) = dW/dS = unit tangent vector to integral curve of
!    d(homotopy map)/dS = 0  at  W(S) = (LAMBDA(S), X(S)) .
!
! TZ = the Newton step = -(pseudo inverse of  (d RHO(W(S))/d LAMBDA ,
!    d RHO(W(S))/dX)) * RHO(W(S)) .
!
! IFLAG  is unchanged, unless the QR factorization detects a rank < N,
!    in which case the tangent and Newton step vectors are not computed
!    and TANGENT_GLP returns with  IFLAG = 4 .
!
!
! Calls  DGEQPF, DNRM2, DORMQR, RHO.

IMPLICIT NONE
LOGICAL, INTENT(IN), OPTIONAL:: NEWTON_STEP
REAL (KIND=R8):: LAMBDA, SIGMA, WPNORM
INTEGER:: I, J, K

INTERFACE
  FUNCTION DNRM2(N,X,STRIDE)
    USE REAL_PRECISION
    INTEGER:: N,STRIDE
    REAL (KIND=R8):: DNRM2,X(N)
  END FUNCTION DNRM2
END INTERFACE

! Compute the Jacobian matrix, store it and homotopy map in QR.
!
!  QR = ( D RHO(LAMBDA,X)/D LAMBDA , D RHO(LAMBDA,X)/DX ,
!                                              RHO(LAMBDA,X) )  .
!
! Force LAMBDA >= 0 for tangent calculation.
IF (W(1) < 0.0_R8) THEN
  LAMBDA = 0.0_R8
ELSE
  LAMBDA = W(1)
END IF

! RHO(W) evaluates the homotopy map and its Jacobian matrix at W,
! leaving the results in the arrays RHOV, DRHOL, and DRHOX.
CALL RHO(LAMBDA,W(2:2*N+1))
QR(1:2*N,1) = DRHOL(1:2*N)
QR(1:2*N,2:2*N+1) = DRHOX(1:2*N,1:2*N)
QR(1:2*N,2*N+2) = RHOV(1:2*N)

! Compute the norm of the homotopy map if it was requested.
IF (RHOLEN < 0.0_R8) RHOLEN = DNRM2(2*N,QR(:,2*N+2),1)

! Reduce the Jacobian matrix to upper triangular form.
PIVOT = 0
CALL DGEQPF(2*N,2*N+1,QR,2*N,PIVOT,WP,ALPHA,K)

IF (ABS(QR(2*N,2*N)) <= ABS(QR(1,1))*EPSILON(1.0_R8)) THEN 
  IFLAG = 4
  RETURN
ENDIF

! Apply Householder reflections to last column of QR (which contains
! RHO(A,W)).
CALL DORMQR('L','T',2*N,1,2*N-1,QR,2*N,WP,QR(:,2*N+2),2*N,  &
  ALPHA, 3*(2*N+1),K)

! Compute kernel of Jacobian matrix, yielding WP=dW/dS.
TZ(2*N+1) = 1.0_R8
DO I=2*N,1,-1
  J = I + 1
  TZ(I) = -DOT_PRODUCT(QR(I,J:2*N+1),TZ(J:2*N+1))/QR(I,I)
END DO
WPNORM = DNRM2(2*N+1,TZ,1)
WP(PIVOT) = TZ/WPNORM
IF (DOT_PRODUCT(WP,YPOLD) < 0.0_R8) WP = -WP

! WP is the unit tangent vector in the correct direction.
IF (.NOT. PRESENT(NEWTON_STEP)) RETURN

! Compute the minimum norm solution of [d RHO(W(S))] V = -RHO(W(S)).
! V is given by  P - (P,Q)Q  , where P is any solution of
! [d RHO] V = -RHO  and Q is a unit vector in the kernel of [d RHO].

ALPHA(2*N+1) = 1.0_R8
DO I=2*N,1,-1
  J = I + 1
  ALPHA(I) = -(DOT_PRODUCT(QR(I,J:2*N+1),ALPHA(J:2*N+1)) + QR(I,2*N+2)) &
            /QR(I,I)
END DO
TZ(PIVOT) = ALPHA(1:2*N+1)

! TZ now contains a particular solution P, and WP contains a vector Q
! in the kernel (the unit tangent).
SIGMA = DOT_PRODUCT(TZ,WP)
TZ = TZ - SIGMA*WP

! TZ is the Newton step from the current point W(S) = (LAMBDA(S), X(S)).
RETURN
END SUBROUTINE TANGENT_GLP

                                                                     !!!
SUBROUTINE ROOT_GLP
! In a deleted neighborhood of a solution (1,X(SBAR)), the homotopy zero
! curve (LAMBDA(S),X(S)) is assumed to safisfy X = X(LAMBDA), a consequence
! of the Implicit Function Theorem and the fact that the Jacobian matrix
! D RHO(A,LAMBDA(S),X(S))/DX is nonsingular in a sufficiently small
! deleted neighborhood of an isolated solution.  Let
!         TAU = 1 - LAMBDA = SIGMA**C,
! where the positive integer C is the cycle number of the root.  Then
!   X(LAMBDA) = X(1 - TAU) = X(1 - SIGMA**C) = Z(SIGMA)
! is an analytic function of SIGMA in a neighborhood of SIGMA=0.  This fact
! is exploited by guessing C and interpolating Z(SIGMA) within its
! Maclaurin series' radius of convergence, but far enough away from 0 to
! avoid numerical instability.  This annulus is called the "operating
! range" of the algorithm.  The interpolant to analytic Z(SIGMA) is then
! evaluated at SIGMA=0 to estimate the root X(1)=Z(0).

! Local variables.
IMPLICIT NONE
INTEGER, PARAMETER:: CHAT_MAX=8, LITFH = 7
INTEGER:: C, CHAT(1), CHAT_BEST, CHAT_OLD, GOING_BAD, I, & 
  J, ML_ITER, N2P1, RETRY
REAL (KIND=R8):: ACCURACY, FV(12), GM, H_SAVE, HC, HQ, HQ_BEST,    & 
  HQMHC(CHAT_MAX), L(-3:2), S_SAVE, SIGMA(-3:2), SHRINK, T, TOL_1, &   
  TOL_2, V(12) 
LOGICAL:: EVEN, FIRST_JUMP, REUSE

INTERFACE
  FUNCTION DNRM2(N,X,STRIDE)
    USE REAL_PRECISION
    INTEGER:: N, STRIDE
    REAL (KIND=R8):: DNRM2, X(N)
  END FUNCTION DNRM2
END INTERFACE

N2P1 = 2*N + 1
ACCURACY = MAX(FINALTOL,SQRT(EPSILON(1.0_R8))*10.0_R8**2)
HQ_BEST = 10.0_R8*ACCURACY
CHAT_BEST = 0 ; CHAT_OLD = 0 ; GOING_BAD = 0
FIRST_JUMP = .TRUE. ; REUSE = .FALSE.
YOLDS = 0.0_R8

! Save the first point.
H_SAVE = HOLD  
S_SAVE = S - HOLD
YS(:,1) = YOLD ; YS(:,2) = YPOLD

! If Y(1) >= 1 or if YP(1) <= 0 back up to YOLD and generate another point.

REFINE_Y: DO

  IF ((Y(1) >= 1.0_R8) .OR. (YP(1) <= 0.0_R8)) THEN
    SHRINK = 1.0_R8

    ! Try 3 times to get a point.
    DO I=1,3
      SHRINK = SHRINK * .75_R8
      S = S_SAVE
      H = MIN(H_SAVE, SHRINK*(1.0_R8 - YS(1,1))/YS(1,2))

      ! If Y(1)>=1 increase RELERR and ABSERR to prevent STEPNX from making
      ! the stepsize too small.
      IF (Y(1) >= 1.0_R8) THEN
        RELERR = TRACKTOL ; ABSERR = TRACKTOL
      END IF
      Y = YS(:,1) ; YP = YS(:,2)
      START = .TRUE.
      CALL STEP_GLP
      RELERR = FINALTOL ; ABSERR = FINALTOL
      IF (IFLAG > 0) THEN
        IFLAG = 4 ; RETURN
      ELSE IF ((Y(1) < 1.0_R8) .AND. (YP(1) > 0.0_R8) .AND.  &
              (Y(1) > YS(1,1))) THEN
        ITER = ITER + 1
        EXIT REFINE_Y
      ELSE IF (I == 3) THEN 
        IFLAG = 7 ; RETURN
      END IF
    END DO
  ELSE

    ! Refine the second point Y to FINALTOL accuracy.  If the refinement
    ! fails, back up and get another point.
    W = Y
    RHOLEN = 0.0_R8
    DO J=1,LITFH
      CALL TANGENT_GLP(NEWTON_STEP=.TRUE.)
      NNFE = NNFE + 1
      IF (IFLAG > 0) THEN
        IFLAG = -2
        YP(1) = -1.0_R8 ; CYCLE REFINE_Y
      END IF
      W = W + TZ

      ! Test for erratic LAMBDA.
      IF (W(1) >= 1.0_R8 .OR. WP(1) <= 0.0_R8 .OR. W(1) <= YS(1,1)) THEN
        YP(1) = -1.0_R8 ; CYCLE REFINE_Y 
      END IF
      IF (DNRM2(N2P1,TZ,1) <= FINALTOL * (DNRM2(N2P1,W,1) + 1.0_R8)) EXIT

      ! Test for lack of convergence.
      IF (J == LITFH) THEN
        YP(1) = -1.0_R8 ; CYCLE REFINE_Y
      END IF
    END DO
    Y = W ; YP = WP
    S = S - HOLD
    W = Y - YOLD
    HOLD = DNRM2(N2P1,W,1)
    S = S + HOLD
    EXIT REFINE_Y
  END IF

END DO REFINE_Y

! Save the second point.
YS(:,3) = Y ; YS(:,4) = YP
H_SAVE = H ; S_SAVE = S

! Try entire end game interpolation process RETRY=10*NUMRR times.
RETRY = 10*NUM_RERUNS

MAIN_LOOP: &
DO ML_ITER=1,RETRY

  ! Get close enough to SIGMA=0 (LAMBDA=1) so that a Hermite cubic
  ! interpolant is accurate to within TOL_1 (defined by CHAT).
  OPERATING_RANGE: DO

    ! Enforce LIMIT on the number of steps.
    IF (ITER >= LIMIT) THEN
      IFLAG = 3 ; EXIT MAIN_LOOP
    END IF

    SHRINK = 1.0_R8
    DO J=1,3
      SHRINK = .75_R8*SHRINK

      ! Get a third point Y with Y(1) < 1.
      H = MIN(H_SAVE, SHRINK*(1.0_R8 - Y(1))/YP(1))
      CALL STEP_GLP
      IF (IFLAG > 0) THEN
        IFLAG = 4 ; EXIT MAIN_LOOP
      ELSE IF ((Y(1) >= 1.0_R8) .OR. (YP(1) <= 0.0_R8) .OR.  &
        (Y(1) <= YS(1,3))) THEN
        ! Back up and try again with a smaller step.
        Y = YS(:,3) ; YP = YS(:,4) ; YOLD = YS(:,1) ; YPOLD = YS(:,2)
        S = S_SAVE
      ELSE
        ITER = ITER + 1
        EXIT
      END IF
      IF (J == 3) THEN 
        IFLAG = 7 ; EXIT MAIN_LOOP
      END IF
    END DO

    ! Save the third point.
    YS(:,5) = Y ; YS(:,6) = YP
    H_SAVE = H ; S_SAVE = S

    ! L(2) < L(1) < L(0) < 1.

    L(2) = YS(1,1) ; L(1) = YS(1,3) ; L(0) = YS(1,5)

    ! Test approximation quality for each cycle number C = 1,...,CHAT_MAX.

    SHRINK = 1.0_R8/(1.0_R8 + MAXVAL(ABS(YS(2:N2P1,5))))
    DO C=1,CHAT_MAX
      SIGMA(0:2) = (1.0_R8 - L(0:2))**(1.0_R8/REAL(C,KIND=R8))

      ! 0 < SIGMA(0) < SIGMA(1) < SIGMA(2).
      ! Compute difference between Hermite quintic HQ(SIGMA) interpolating at
      ! SIGMA(0:2) and Hermite cubic HC(SIGMA) interpolating at SIGMA(0:1).
      ! The interpolation points for the Newton form are (SIGMA(0), SIGMA(0),
      ! SIGMA(1), SIGMA(1), SIGMA(2), SIGMA(2)).  The function values are in
      ! YS(:,5:1:-2) and the derivatives YS(:,6:2:-2) = dX/dS have to be
      ! converted to dX/dSIGMA.

      T = 0.0_R8
      V(1:6) = (/ (SIGMA(J),SIGMA(J),J=0,2) /)
      DO J=2,N2P1
        FV(1:5:2) = YS(J,5:1:-2)
        FV(2:6:2) = (YS(J,6:2:-2)/YS(1,6:2:-2)) * (-REAL(C,KIND=R8)) * &
          SIGMA(0:2)**(C-1)
        CALL INTERP(V(1:6),FV(1:6))
        T = MAX(T,ABS(FV(5) - SIGMA(2)*FV(6)))
      END DO

      ! T*(SIGMA(1)*SIGMA(0))**2 = ||HQ(0) - HC(0)||_infty.

      HQMHC(C) = T*((SIGMA(1)*SIGMA(0))**2)*SHRINK
    END DO

    ! Find best estimate CHAT of cycle number.
    CHAT = MINLOC(HQMHC)

    ! If there has been one successful jump across the origin (with
    ! CHAT_BEST) and the cycle number prediction changes, then the process
    ! may be leaving the operating range.

    IF (( .NOT. FIRST_JUMP) .AND. (CHAT(1) /= CHAT_BEST)) THEN
      GOING_BAD = GOING_BAD + 1
      IF (GOING_BAD == 2) EXIT MAIN_LOOP
    END IF
    TOL_1 = ACCURACY*10.0_R8**(REAL(CHAT(1),KIND=R8)/2.0_R8)
    IF (HQMHC(CHAT(1)) <= TOL_1) THEN
      EXIT OPERATING_RANGE
    ELSE IF ( .NOT. FIRST_JUMP) THEN
      GOING_BAD = GOING_BAD + 1
      IF (GOING_BAD == 2) EXIT MAIN_LOOP
    END IF

    ! Shift point history, and try to get closer to SIGMA=0.
    YS(:,1:2) = YS(:,3:4) ; YS(:,3:4) = YS(:,5:6) ; REUSE = .FALSE.
  END DO OPERATING_RANGE

  ! Add 3 new points past SIGMA=0 such that
  ! SIGMA(2) > SIGMA(1) > SIGMA(0) > 0 > SIGMA(-1) > SIGMA(-2) > SIGMA(-3).
  ! If CHAT is odd then the corresponding LAMBDA are such that
  !   L(2)   <   L(1)   <   L(0)   < 1 <   L(-1)   <   L(-2)   <   L(-3),
  ! and if CHAT is even then
  !   L(2)   <   L(1)   <   L(0)   < 1
  !                                  1 >   L(-1)   >   L(-2)   >   L(-3).

  SIGMA(0:2) = (1.0_R8 - L(0:2))**(1.0_R8/REAL(CHAT(1),KIND=R8))
  DO I=1,3
    V(1:4+2*I) = (/ (SIGMA(J),SIGMA(J),J=2,1-I,-1) /)
    DO J=2,N2P1
      FV(1:3+2*I:2) = YS(J,1:3+2*I:2)
      FV(2:4+2*I:2) = (YS(J,2:4+2*I:2)/YS(1,2:4+2*I:2)) * &
        (-REAL(CHAT(1),KIND=R8)) * SIGMA(2:1-I:-1)**(CHAT(1)-1)
      CALL INTERP(V(1:4+2*I),FV(1:4+2*I))
      CALL INTERP(V(1:4+2*I),FV(1:4+2*I),-SIGMA(I-1),W(J))
    END DO
    IF (MOD(CHAT(1),2) == 0) THEN
      EVEN = .TRUE.
      W(1) = L(I-1)
    ELSE
      EVEN = .FALSE.
      W(1) = 2.0_R8 - L(I-1)
    END IF

    ! W now contains the (predicted) point symmetric to SIGMA(I-1) with
    ! respect to SIGMA=0.
    RHOLEN = 0.0_R8

    ! Correct the prediction.  If there has been one successful jump across
    ! the origin, correction failures may indicate that the process is
    ! leaving the operating range.
    DO J=1,LITFH
      CALL TANGENT_GLP(NEWTON_STEP=.TRUE.)
      NNFE = NNFE + 1

      ! Test for singular Jacobian matrix.
      IF (IFLAG > 0) EXIT MAIN_LOOP
      W = W + TZ

      ! Test for erratic LAMBDA.
      IF ((( .NOT. EVEN) .AND. (W(1) <= 1.0_R8)) .OR.  &
          (EVEN .AND. (W(1) >= 1.0_R8))) THEN
        IF ( .NOT. FIRST_JUMP) THEN
          GOING_BAD = GOING_BAD + 1
          IF (GOING_BAD == 2) EXIT MAIN_LOOP
        END IF  
        YS(:,1:2) = YS (:,3:4) ; YS(:,3:4) = YS(:,5:6) 
        REUSE = .FALSE. ; CYCLE MAIN_LOOP
      END IF
      IF (DNRM2(N2P1,TZ,1) <= FINALTOL * (DNRM2(N2P1,W,1) + 1.0_R8)) EXIT

      ! Test for lack of convergence.
      IF (J == LITFH) THEN
        IF ( .NOT. FIRST_JUMP) THEN
          GOING_BAD = GOING_BAD + 1
          IF (GOING_BAD == 2) EXIT MAIN_LOOP
        END IF
        YS(:,1:2) = YS (:,3:4) ; YS(:,3:4) = YS(:,5:6)
        REUSE = .FALSE. ; CYCLE MAIN_LOOP
      END IF
    END DO

    ! Ensure that the tangent vector has the correct direction.
    IF (EVEN) THEN
      IF (WP(1) > 0.0_R8) WP = -WP
    ELSE
      IF (WP(1) < 0.0_R8) WP = -WP
    END IF

    ! Update the lambda (L), sigma (SIGMA), and history (YS) arrays.
    L(-I) = W(1)
    SIGMA(-I) =  -(ABS(L(-I) - 1.0_R8))**(1.0_R8/REAL(CHAT(1),KIND=R8))
    YS(:,5+2*I) = W ; YS(:,6+2*I) = WP

    ! Reuse old points if the cycle number estimation has not changed
    ! from the last iteration, and the origin was successfully jumped in
    ! the last iteration.
    IF (REUSE .AND. (CHAT(1) == CHAT_OLD)) EXIT
  END DO

  ! Construct 12th order interpolant and estimate the root at SIGMA=0. 
  HC = 0.0_R8 ; HQ = 0.0_R8 ; T = 0.0_R8
  V(1:12) = (/ (SIGMA(J),SIGMA(J),J=-3,2) /)
  DO J=2,N2P1
    FV(1:11:2) = YS(J,11:1:-2)
    FV(2:12:2) = (YS(J,12:2:-2)/YS(1,12:2:-2)) * &
      (-REAL(CHAT(1),KIND=R8)) * SIGMA(-3:2)**(CHAT(1)-1)
    CALL INTERP(V(1:12),FV(1:12))
    CALL INTERP(V(1:12),FV(1:12),0.0_R8,W(J))

    ! Difference between 8th and 6th order Hermite interpolants.
    T  = MAX(T ,ABS(FV( 7) - SIGMA(0)*FV( 8)))

    ! Difference between 10th and 8th order Hermite interpolants.
    HC = MAX(HC,ABS(FV( 9) - SIGMA(1)*FV(10)))

    ! Difference between 12th and 10th order Hermite interpolants.
    HQ = MAX(HQ,ABS(FV(11) - SIGMA(2)*FV(12)))
  END DO
  SHRINK = 1.0_R8/(1.0_R8 + MAXVAL(ABS(W(2:N2P1))))
  T  =  T*((PRODUCT(SIGMA(-3:-1)))**2)*SHRINK  ! ||H_7  - H_5||/(1+||W||)
  HC = HC*((PRODUCT(SIGMA(-3: 0)))**2)*SHRINK  ! ||H_9  - H_7||/(1+||W||)
  HQ = HQ*((PRODUCT(SIGMA(-3: 1)))**2)*SHRINK  ! ||H_11 - H_9||/(1+||W||)

  ! Check both accuracy and consistency of Hermite interpolants.
  TOL_2 = FINALTOL * (10**(CHAT(1) - 1))
  GM = SQRT(TOL_1 * TOL_2)
  IF ((T <= TOL_1) .AND. (HC <= GM) .AND. (HQ <= TOL_2)) THEN

    ! Full convergence.
    IF (FIRST_JUMP) FIRST_JUMP = .FALSE.
    YOLDS(2:N2P1) = W(2:N2P1) ; HQ_BEST = HQ
    CHAT_BEST = CHAT(1)
    EXIT MAIN_LOOP
  ELSE IF (HQ > 1.01_R8*HQ_BEST) THEN
    IF ( .NOT. FIRST_JUMP) THEN
      GOING_BAD = GOING_BAD + 1
      IF (GOING_BAD == 2) EXIT MAIN_LOOP
    END IF
  ELSE

    ! Progress has been made.
    IF (FIRST_JUMP) FIRST_JUMP = .FALSE.
    GOING_BAD = 0
    YOLDS(2:N2P1) = W(2:N2P1) ; HQ_BEST = HQ
    CHAT_BEST = CHAT(1)
  END IF

  ! Shift point history.
  YS(:,1:2) = YS(:,3:4) ; YS(:,3:4) = YS(:,5:6)

  ! If the cycle number estimate does not change in the next iteration, the
  ! points found across the origin can be reused.
  REUSE = .TRUE. ; CHAT_OLD = CHAT(1)
  SIGMA(-3)   = SIGMA(-2)  ; SIGMA(-2)  = SIGMA(-1)
  YS(:,11:12) = YS(:,9:10) ; YS(:,9:10) = YS(:,7:8)

END DO MAIN_LOOP

IF (ML_ITER >= RETRY) IFLAG = 7

! Return final solution in Y.
IF ( .NOT. FIRST_JUMP) THEN
  Y(1) = HQ_BEST ; Y(2:N2P1) = YOLDS(2:N2P1)
  IFLAG = 1 + 10*CHAT_BEST
END IF
RETURN
END SUBROUTINE ROOT_GLP

                                                                     !!!
SUBROUTINE INTERP(T,FT,X,FX)
! Given data points T(:) and function values FT(:)=f(T(:)), INTERP
! computes the Newton form of the interpolating polynomial to f at T(:).
! T is assumed to be sorted, and if
! T(I-1) < T(I) = T(I+1) = ... = T(I+K) < T(I+K+1) then
! FT(I)=f(T(I)), FT(I+1)=f'(T(I)), ..., FT(I+K)=f^{(K)}(T(I)).
! On return FT(K) contains the divided difference f[T(1),...,T(K)], and
! FX contains the interpolating polynomial evaluated at X.  If X and FX
! are present, the divided differences are not calculated.

IMPLICIT NONE
REAL (KIND=R8), DIMENSION(:):: T, FT
REAL (KIND=R8), OPTIONAL:: X, FX

! Local variables.
REAL (KIND=R8):: FOLD,SAVE
INTEGER:: I,K,N

N = SIZE(T)
IF (.NOT. PRESENT(X)) THEN ! Calculate divided differences.
  DO K=1,N-1
    FOLD = FT(K)
    DO I=K+1,N
      IF (T(I) == T(I-K)) THEN
        FT(I) = FT(I)/REAL(K,KIND=R8)
      ELSE
        SAVE = FT(I)
        FT(I) = (FT(I) - FOLD)/(T(I) - T(I-K))
        FOLD = SAVE
      END IF
    END DO
  END DO
  RETURN
END IF
FX = FT(N) ! Evaluate Newton polynomial.
DO K=N-1,1,-1
  FX = FX*(X - T(K)) + FT(K)
END DO
RETURN
END SUBROUTINE INTERP

                                                                     !!!
SUBROUTINE RHO(LAMBDA,X)
! RHO evaluates the (complex) homotopy map
!
!   RHO(A,LAMBDA,X) = LAMBDA*F(X) + (1 - LAMBDA)*GAMMA*G(X),
!
! where GAMMA is a random complex constant, and the Jacobian
! matrix [ D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX ] at
! (A,LAMBDA,X), and updates the global arrays RHOV (the homotopy map),
! DRHOX (the derivative of the homotopy with repect to X) , and DRHOL
! (the derivative with respect to LAMBDA).  The vector A corresponds
! mathematically to all the random coefficients in the start system, and
! is not explicitly referenced by RHO.  X, on entry, is real, but since
! arithmetic in RHO is complex, X is converted to complex form.  Before
! return RHO converts the homotopy map and the two derivatives back to
! real.  Precisely, suppose XC is the complexification of X, i.e.,
!
! XC(1:N)=CMPLX(X(1:2*N-1:2),X(2:2*N:2)).
!
! Let CRHOV(A,LAMBDA,XC) be the (complex) homotopy map.  Then RHOV
! is just
!
! RHOV(1:2*N-1:2) = REAL( CRHOV(1:N)),
! RHOV(2:2*N  :2) = AIMAG(CRHOV(1:N)).
!
! Let CDRHOXC = D CRHOV(A,LAMBDA,XC)/D XC denote the (complex) derivative
! of the homotopy map with respect to XC, evaluated at (A,LAMBDA,XC).
! DRHOX is obtained by
!
! DRHOX(2*I-1,2*J-1) = REAL(CDRHOXC(I,J)),
! DRHOX(2*I  ,2*J  ) = DRHOX(2*I-1,2*J-1),
! DRHOX(2*I  ,2*J-1) = AIMAG(CDRHOXC(I,J)),
! DRHOX(2*I-1,2*J  ) = -DRHOX(2*I  ,2*J-1),
!
! for I, J = 1,...,N.  Let CDRHOL = D CRHOV(A,LAMBDA,XC)/D LAMBDA denote
! the (complex) derivative of the homotopy map with respect to LAMBDA,
! evaluated at (A,LAMBDA,XC).  Then DRHOL is obtained by
!
! DRHOL(1:2*N-1:2) = REAL( CDRHOL(1:N)),
! DRHOL(2:2*N  :2) = AIMAG(CDRHOL(1:N)).
!
! (None of CRHOV, CDRHOXC, or CDRHOL are in the code.)
!
! Internal subroutines:  START_SYSTEM, TARGET_SYSTEM.
! External (optional, user written) subroutine:  TARGET_SYSTEM_USER.
!
! On input:
!
! LAMBDA is the continuation parameter.
!
! X(1:2*N) is the real 2*N-dimensional evaluation point.
!
! On exit:
!
! LAMBDA and X are unchanged.
!
! RHOV(1:2*N) is the real (2*N)-dimensional representation of the
!   homotopy map RHO(A,LAMBDA,X).
!
! DRHOX(1:2*N,1:2*N) is the real (2*N)-by-(2*N)-dimensional
!   representation of D RHO(A,LAMBDA,X)/DX evaluated at (A,LAMBDA,X).
!
! DRHOL(1:2*N) is the real (2*N)-dimensional representation of
!   D RHO(A,LAMBDA,X)/D LAMBDA evaluated at (A,LAMBDA,X).

IMPLICIT NONE
REAL (KIND=R8), INTENT(IN):: LAMBDA
REAL (KIND=R8), DIMENSION(2*N), INTENT(IN):: X

INTERFACE
SUBROUTINE TARGET_SYSTEM_USER(N,PROJ_COEF,XC,F,DF)
  USE REAL_PRECISION
  INTEGER, INTENT(IN):: N
  COMPLEX (KIND=R8), INTENT(IN), DIMENSION(N+1):: PROJ_COEF,XC
  COMPLEX (KIND=R8), INTENT(OUT):: F(N), DF(N,N+1)
END SUBROUTINE TARGET_SYSTEM_USER
END INTERFACE

! Local variables.
INTEGER:: I, J
REAL (KIND=R8):: ONEML
COMPLEX (KIND=R8):: GAMMA

ONEML = 1.0_R8 - LAMBDA
GAMMA = (.0053292102547824_R8,.9793238462643383_R8)

! Convert the real-valued evaluation point X to a complex vector.
XC(1:N) = CMPLX(X(1:2*N-1:2),X(2:2*N:2),KIND=R8)

! Calculate the homogeneous variable.
XC(N+1) = SUM(PROJ_COEF(1:N)*XC(1:N)) + PROJ_COEF(N+1)

CALL START_SYSTEM              ! Returns G and DG.
IF (PRESENT(USER_F_DF)) THEN   ! Returns F and DF.
  CALL TARGET_SYSTEM_USER(N,PROJ_COEF,XC,F,DF) ! User written subroutine.
ELSE
  CALL TARGET_SYSTEM ! Internal subroutine.
END IF

! Convert complex derivatives to real derivatives via the Cauchy-Riemann
! equations.
DO I=1,N
  DO J=1,N  
    DRHOX(2*I-1,2*J-1) = LAMBDA*REAL(DF(I,J)) + ONEML*REAL(DG(I,J)*GAMMA)
    DRHOX(2*I  ,2*J  ) = DRHOX(2*I-1,2*J-1)
    DRHOX(2*I  ,2*J-1) = LAMBDA*AIMAG(DF(I,J)) + ONEML*AIMAG(DG(I,J)*GAMMA)
    DRHOX(2*I-1,2*J  ) = -DRHOX(2*I,2*J-1)
  END DO
END DO
  DRHOL(1:2*N-1:2) = REAL(F) - REAL(G*GAMMA)
  DRHOL(2:2*N:2  ) = AIMAG(F) - AIMAG(G*GAMMA) 
  RHOV(1:2*N-1:2)  = LAMBDA*REAL(F) + ONEML*REAL(G*GAMMA)
  RHOV(2:2*N:2  )  = LAMBDA*AIMAG(F) + ONEML*AIMAG(G*GAMMA)
RETURN
END SUBROUTINE RHO

                                                                     !!!
SUBROUTINE START_SYSTEM
! START_SYSTEM evaluates the start system G(XC) and the Jacobian matrix
!   DG(XC).  Arithmetic is complex.
!
! On exit:
!
! G(:) contains the complex N-dimensional start system evaluated at XC(:).
!
! DG(:,:) contains the complex N-by-N-dimensional Jacobian matrix of
!   the start system evaluted at XC(:). 

! Local variables.
IMPLICIT NONE
INTEGER:: I, J, K, L
COMPLEX (KIND=R8):: TEMP
COMPLEX (KIND=R8), DIMENSION(N):: TEMP2

! TEMP1G AND TEMP2G are employed to reduce recalculation in G and DG.
! Note:  If SD(I,J)=0, then the corresponding factor is 1, not 0.
TEMP1G = (0.0_R8,0.0_R8)
TEMP2G = (0.0_R8,0.0_R8)
DO I=1,N
  DO J=1,COVER_SIZES(I)
    IF (COVER(I)%SET(J)%SET_DEG == 0) THEN
      TEMP2G(I,J) = (1.0_R8,0.0_R8)
    ELSE
      K = COVER(I)%SET(J)%NUM_INDICES
      TEMP1G(I,J) = SUM( COVER(I)%SET(J)%START_COEF(1:K)*  &
                         XC(COVER(I)%SET(J)%INDEX(1:K)) )
      TEMP2G(I,J) = TEMP1G(I,J)**COVER(I)%SET(J)%SET_DEG - &
                        XC(N+1)**COVER(I)%SET(J)%SET_DEG
    END IF
  END DO
  G(I) = PRODUCT(TEMP2G(I,1:COVER_SIZES(I)))
END DO

! Calculate the derivative of G with respect to XC(1),...,XC(N)
!   in 3 steps.
! STEP 1:  First treat XC(N+1) as an independent variable. 

DG = (0.0_R8,0.0_R8)  
DO I=1,N
  DO J=1,COVER_SIZES(I)
    IF (COVER(I)%SET(J)%SET_DEG == 0) CYCLE
    K = COVER(I)%SET(J)%NUM_INDICES
    TEMP2(1:K) = COVER(I)%SET(J)%SET_DEG * &
      COVER(I)%SET(J)%START_COEF(1:K) * &
      (TEMP1G(I,J)**(COVER(I)%SET(J)%SET_DEG - 1))
    TEMP = (1.0_R8,0.0_R8)
    DO L=1,COVER_SIZES(I)
      IF (L == J) CYCLE
      TEMP = TEMP * TEMP2G(I,L)
    END DO

    DG(I,COVER(I)%SET(J)%INDEX(1:K)) = &
      DG(I,COVER(I)%SET(J)%INDEX(1:K)) + TEMP2(1:K) * TEMP
  END DO
END DO

! STEP 2:  Now calculate the N-by-1 Jacobian matrix of G with
!   respect to XC(N+1) using the product rule.

DO I=1,N
  DO J=1,COVER_SIZES(I)
    IF (COVER(I)%SET(J)%SET_DEG == 0) CYCLE            
    TEMP = -COVER(I)%SET(J)%SET_DEG * &
      (XC(N+1)**(COVER(I)%SET(J)%SET_DEG - 1))
    DO K=1,COVER_SIZES(I)
      IF (K == J) CYCLE
      TEMP = TEMP*TEMP2G(I,K)
    END DO
    DG(I,N+1) = DG(I,N+1) + TEMP
  END DO
END DO

! STEP 3:  Use the chain rule with XC(N+1) considered as a function
!   of XC(1),...,XC(N).

DO I=1,N
  DG(I,1:N) = DG(I,1:N) + DG(I,N+1) * PROJ_COEF(1:N)
END DO
RETURN
END SUBROUTINE START_SYSTEM

                                                                     !!!
SUBROUTINE TARGET_SYSTEM
! TARGET_SYSTEM calculates the target system F(XC) and the Jacobian matrix
!   DF(XC).  Arithmetic is complex.
!
! On exit:
!
! F(:) contains the complex N-dimensional target system evaluated
!   at XC(:).
!
! DF(:,:) is the complex N-by-N-dimensional Jacobian matrix of the
!   target system evaluated at XC(:).

! Local variables.
IMPLICIT NONE
INTEGER:: I, J, K, L
COMPLEX (KIND=R8):: T, TS

! Evaluate F(XC).  For efficiency, indexing functions and array sections
! are avoided.
DO I=1,N
  TS = (0.0_R8, 0.0_R8)
  DO J=1,POLYNOMIAL(I)%NUM_TERMS
    T = POLYNOMIAL(I)%TERM(J)%COEF
    DO K=1,N+1
      IF (POLYNOMIAL(I)%TERM(J)%DEG(K) == 0) CYCLE
      T = T * XC(K)**POLYNOMIAL(I)%TERM(J)%DEG(K)
    END DO
    TS = TS + T
  END DO
  F(I) = TS
END DO 

! Calulate the Jacobian matrix DF(XC).
DF = (0.0_R8,0.0_R8)
DO I=1,N
  DO J=1,N+1
    TS = (0.0_R8,0.0_R8)
    DO K=1,POLYNOMIAL(I)%NUM_TERMS
      IF (POLYNOMIAL(I)%TERM(K)%DEG(J) == 0) CYCLE
      T = POLYNOMIAL(I)%TERM(K)%COEF * POLYNOMIAL(I)%TERM(K)%DEG(J) * &
          (XC(J)**(POLYNOMIAL(I)%TERM(K)%DEG(J) - 1))
      DO L=1,N+1
        IF ((L == J) .OR. (POLYNOMIAL(I)%TERM(K)%DEG(L) == 0)) CYCLE
        T = T * (XC(L)**POLYNOMIAL(I)%TERM(K)%DEG(L))
      END DO
      TS = TS + T
    END DO
    DF(I,J) = TS
  END DO
END DO

! Convert DF to partials with respect to XC(1),...,XC(N) by
!   applying the chain rule with XC(N+1) considered as a function
!   of XC(1),...,XC(N).
DO I=1,N
  DF(I,1:N) = DF(I,1:N) + PROJ_COEF(1:N) * DF(I,N+1)
END DO
RETURN
END SUBROUTINE TARGET_SYSTEM

                                                                     !!!
SUBROUTINE OUTPUT_GLP
! OUTPUT_GLP first untransforms (converts from projective to affine
! coordinates) and then unscales a root.
!
! On entry:
!
! XC(1:N) contains a root in projective coordinates, with the (N+1)st
!   projective coordinate XC(N+1) implicitly defined by the
!   projective transformation.
!
! On exit:
!
! XC(1:N) contains the untransformed (affine), unscaled root.
!
! XC(N+1) is the homogeneous coordinate of the root of the scaled
!   target system, if scaling was performed. 

IMPLICIT NONE
INTEGER:: I
REAL (KIND=R8), PARAMETER:: BIG=HUGE(1.0_R8)

! Calculate the homogeneous coordinate XC(N+1) using the vector XC(1:N)
! with the projective transformation, then untransform XC(1:N) (convert
! to affine coordinates).
XC(N+1) = SUM(PROJ_COEF(1:N)*XC(1:N)) + PROJ_COEF(N+1)

! Deal carefully with solutions at infinity.
IF (ABS(XC(N+1)) < 1.0_R8) THEN
  DO I=1,N
    IF (ABS(XC(I)) >= BIG*ABS(XC(N+1))) THEN
      XC(I) = CMPLX(BIG,BIG,KIND=R8)  ! Solution at infinity.
    ELSE
      XC(I) = XC(I)/XC(N+1)
    END IF
  END DO
ELSE
  XC(1:N) = XC(1:N)/XC(N+1)
END IF

! Unscale the variables.
IF (.NOT. PRESENT(NO_SCALING)) THEN
  DO I=1,N
    IF (REAL(XC(I)) /= BIG) XC(I) = XC(I)*(10.0_R8**SCALE_FACTORS(I))
  END DO
END IF

RETURN
END SUBROUTINE OUTPUT_GLP
END SUBROUTINE POLSYS_GLP

                                                                     !!!
SUBROUTINE CHECK_GLP(N, VALID_FLAG, POL_INDEX, TERM_INDEX)
! CHECK_GLP checks that the monomials of the target system are in the span
!   of the set structure specified for the start system. 
! For each monomial of each target polynomial: 
! 1. Form the array FLAG_MTX with the degree of the monomial as NUM_ROWS,
!    and the degree of the set structure as NUM_COLS. 
! 2. Initialize FLAG_MTX so that an element is 1 if the variable is present 
!      in the set and 0 otherwise.
! 3. Scan FLAG_MTX to determine whether a valid assignment exists.  
!
! On entry:
! N: dimension of polynomial system.
!
! On exit:
! VALID_FLAG
!    = 0   for a normal return.
!    = -1  if a monomial is not in the specified set structure.
!    = -2  if the total degree of a term is greater than the total degree 
!          of the set structure.
!    = -3  if a variable is not in the set structure.
!
! POL_INDEX provides, when VALID_FLAG<0, the polynomial index of the monomial 
!   that is not in the specified GLP set structure. 
!
! TERM_INDEX provides, when VALID_FLAG<0, the term index of the monomial
!   that is not in the specified GLP set structure.


USE GLOBAL_GLP
IMPLICIT NONE
INTEGER, INTENT(IN):: N
INTEGER, INTENT(OUT):: VALID_FLAG
INTEGER, INTENT(OUT):: POL_INDEX,TERM_INDEX

! Local variables.
INTEGER:: COLS, FLAG, I, J, K, L, MAX_TOT_DEG, NUM_COLS, NUM_ROWS, ROWS
INTEGER, DIMENSION(N):: TOT_DEG

! FLAG_MTX(I,J) = 1 if the variable corresponding to row I is present in
! the set corresponding to column J, and 0 otherwise.
INTEGER, DIMENSION(:,:), ALLOCATABLE:: FLAG_MTX
! COL_TAKEN(I) = 1 if the Ith column is occupied and 0 otherwise.
INTEGER, DIMENSION(:), ALLOCATABLE:: COL_TAKEN
! LEX_NUM(I) is the column index that the Ith row is scanning.
INTEGER, DIMENSION(:), ALLOCATABLE:: LEX_NUM

DO I=1,N
  TOT_DEG(I) = SUM( (/ (SD(I,J),J=1,COVER_SIZES(I)) /) )
END DO
MAX_TOT_DEG = MAXVAL(TOT_DEG(1:N))

ALLOCATE(FLAG_MTX(MAX_TOT_DEG, MAX_TOT_DEG))
ALLOCATE(LEX_NUM(MAX_TOT_DEG))
ALLOCATE(COL_TAKEN(MAX_TOT_DEG))

FLAG_MTX = 0

POL_LOOP: DO POL_INDEX=1,N
  NUM_COLS = TOT_DEG(POL_INDEX)
  TERM_LOOP: DO TERM_INDEX=1, NUMT(POL_INDEX)    
    VALID_FLAG = 0 ! Valid.
    I = POL_INDEX ! I is the short name of POL_INDEX.
    L = TERM_INDEX ! L is the short name of TERM_INDEX.
    NUM_ROWS = SUM( (/(D(I,L,K),K=1,N)/) )
    IF (NUM_ROWS ==0 ) CYCLE TERM_LOOP
    IF (NUM_ROWS > NUM_COLS ) THEN
      VALID_FLAG = -2 ! The total degree is too high.
      EXIT POL_LOOP 
    END IF

    ! Set flag matrix.
    ROWS = 1
    DO K=1,N
      IF (D(I, L, K)==0) CYCLE
      COLS = 1
      DO J=1, COVER_SIZES(I)
        IF (SD(I,J) ==0) CYCLE
        FLAG = 0
        IF (ANY(COVER(I)%SET(J)%INDEX(1:NUMV(I,J)) == K )) FLAG = 1
        FLAG_MTX( ROWS:ROWS+D(I,L,K)-1, COLS:COLS+SD(I,J)-1) = FLAG
        COLS = COLS + SD(I,J)
      END DO
      ROWS = ROWS + D(I,L,K)
    END DO

    ! If all the elements in a row are all zero, a variable is present in
    ! the term that is not in any set in the cover.
    DO K=1,NUM_ROWS
      IF (ALL(FLAG_MTX(K,1:NUM_COLS) == 0 )) THEN
        VALID_FLAG = -3
        EXIT POL_LOOP
      END IF
    END DO
 
    CALL CHECK_ONE_TERM
    IF (VALID_FLAG < 0) EXIT POL_LOOP ! Invalid set covering.
  END DO TERM_LOOP
END DO POL_LOOP

DEALLOCATE(LEX_NUM,COL_TAKEN,FLAG_MTX)
RETURN

CONTAINS

SUBROUTINE CHECK_ONE_TERM()
! CHECK_ONE_TERM scans FLAG_MTX and attempts to construct the array LEX_NUM 
! of column indices of FLAG_MTX such that
!   FLAG_MTX(I,LEX_NUM(I)) = 1  for each I=1,...,NUM_ROWS,
! and all the elements of LEX_NUM are distinct.  Success means that the
! given term is in the span of the polynomials defined by the component
! set covering.
!
! On entry:
! FLAG_MTX(1:NUM_ROWS,1:NUM_COLS) is a matrix with entries that are 1 or 0.
!
! On exit:
! VALID_FLAG 
!    =  0  if a valid array LEX_NUM exists.
!    = -1  otherwise.

IMPLICIT NONE
INTEGER:: II, JJ, KK, LL 

LEX_NUM = 0
COL_TAKEN = 0
KK = 1

CHECK_LOOP: DO
  IF (KK==0) THEN
    VALID_FLAG = -1 ! Set covering is invalid.
    RETURN 
  END IF
  DO II=KK,NUM_ROWS
    JJ = 0
    DO LL=LEX_NUM(II)+1,NUM_COLS
      IF ((FLAG_MTX(II,LL)==1) .AND. (COL_TAKEN(LL)==0)) THEN
        JJ = LL
        EXIT
      END IF
    END DO
    IF ( JJ==0 )  THEN
      LEX_NUM(II:NUM_ROWS) = 0
        DO LL=II-1, NUM_ROWS
          IF (LL > 0 .AND. LEX_NUM(LL) > 0) COL_TAKEN(LEX_NUM(LL)) = 0
        END DO
      KK = II-1
      CYCLE CHECK_LOOP
    END IF
    COL_TAKEN(JJ) = 1
    LEX_NUM(II) = JJ
  END DO
  RETURN ! Set covering is valid (for this term).
END DO CHECK_LOOP
END SUBROUTINE CHECK_ONE_TERM
END SUBROUTINE CHECK_GLP

                                                                     !!!
SUBROUTINE BEZOUT_GLP(N,MAXT,TOL,BGLP,  POL_INDEX,TERM_INDEX)
!
! BEZOUT_GLP calculates and returns only the generalized Bezout number
! BGLP of the target polynomial system, based on the variable covering
! P defined in the module GLOBAL_GLP.  BEZOUT_GLP finds BGLP very
! quickly, which is useful for exploring alternative coverings.
!
! Calls SINGSYS_GLP.
! 
! On input:
!
! N is the dimension of the target system.
!
! MAXT is the maximum number of terms in any component of the target
!   system.  MAXT = MAX((/(NUMT(I),I=1,N)/)).
!
! TOL is the singularity test threshold used by SINGSYS_GLP.  If
!   TOL <= 0.0 on input, TOL is reset to the default value
!   SQRT(EPSILON(1.0_R8)).
!
! GLOBAL_GLP allocatable objects POLYNOMIAL, COVER_SIZES, and
!   COVER (see GLOBAL_GLP documentation) must be allocated and
!   defined in the calling program.
!
! On output:
!
! N and MAXT are unchanged, and TOL may have been changed as described
!   above.
!
! BGLP is the generalized Bezout number for the target system based on
!   the variable covering P defined in the module GLOBAL_GLP.  BGLP = -1
!   indicates that the provided system covering P is inconsistent
!   with the given target polynomial system.
!
! POL_INDEX, TERM_INDEX are optional variables that give the polynomial
!   index and term index, respectively, of the target system monomial that
!   is inconsistent with the given system set covering.  These are useful
!   for debugging an incorrect input system covering.

USE GLOBAL_GLP
IMPLICIT NONE
INTEGER, INTENT(IN):: N, MAXT
REAL (KIND=R8), INTENT(IN OUT):: TOL
INTEGER, INTENT(OUT):: BGLP
INTEGER, OPTIONAL:: POL_INDEX, TERM_INDEX

!INTERFACE
!  SUBROUTINE CHECK_GLP(N, VALID_FLAG, POL_INDEX, TERM_INDEX)
!  USE GLOBAL_GLP
!  INTEGER, INTENT(IN):: N
!  INTEGER, INTENT(OUT):: VALID_FLAG
!  INTEGER, INTENT(OUT):: POL_INDEX,TERM_INDEX
!  END SUBROUTINE CHECK_GLP
!  SUBROUTINE SINGSYS_GLP(N,LEX_NUM,LEX_SAVE,TOL,RAND_MAT,MAT,NONSING)
!  USE GLOBAL_GLP
!  INTEGER, INTENT(IN):: N
!  INTEGER, DIMENSION(N), INTENT(IN OUT):: LEX_NUM,LEX_SAVE
!  REAL (KIND=R8), INTENT(IN):: TOL
!  REAL (KIND=R8), DIMENSION(N,N), INTENT(IN):: RAND_MAT
!  REAL (KIND=R8), DIMENSION(N+1,N), INTENT(IN OUT):: MAT
!  LOGICAL, INTENT(OUT):: NONSING
!  END SUBROUTINE SINGSYS_GLP
!END INTERFACE

! Local variables.
INTEGER:: J, K, L, POL_IND, TERM_IND, VALID_GLP
INTEGER, DIMENSION(MAXT):: DHOLD
INTEGER, DIMENSION(N):: LEX_NUM, LEX_SAVE
REAL (KIND=R8), DIMENSION(N+1,N):: MAT
REAL (KIND=R8), DIMENSION(N,N):: RAND_MAT
REAL, DIMENSION(N,N):: RANDNUMS
LOGICAL:: NONSING

! Set default value for singularity threshold TOL.
IF (TOL <= REAL(N,KIND=R8)*EPSILON(1.0_R8)) TOL = SQRT(EPSILON(1.0_R8))

! Initialize RAND_MAT with random numbers uniformly distributed in
! [-1,-1/2] union [1/2,1].
CALL RANDOM_SEED
CALL RANDOM_NUMBER(HARVEST=RANDNUMS)
RANDNUMS = RANDNUMS - 0.5 + SIGN(0.5, RANDNUMS - 0.5)
RAND_MAT = REAL(RANDNUMS,KIND=R8)

! GLP does not calculate set degrees, these must be provided as input.

! Check that the system set covering is valid.
CALL CHECK_GLP(N, VALID_GLP, POL_IND, TERM_IND)
IF (VALID_GLP < 0)  THEN
  BGLP = -1
  IF (PRESENT(POL_INDEX)) POL_INDEX = POL_IND
  IF (PRESENT(TERM_INDEX)) TERM_INDEX = TERM_IND
  RETURN
END IF

! Compute Bezout number using lexicographic ordering.
BGLP = 0
LEX_NUM(1:N-1) = 1
LEX_NUM(N) = 0
LEX_SAVE = 0

MAIN_LOOP: DO 
  DO J=N,1,-1
    IF (LEX_NUM(J) < COVER_SIZES(J)) THEN
      L = J
      EXIT
    END IF
  END DO
  LEX_NUM(L) = LEX_NUM(L) + 1
  IF (L + 1 <= N) LEX_NUM(L+1:N) = 1

  ! Test singularity of start subsystem corresponding to lexicographic
  ! vector LEX_NUM.
  CALL SINGSYS_GLP(N,LEX_NUM,LEX_SAVE,TOL,RAND_MAT,MAT,NONSING)
  IF (NONSING) THEN
    BGLP = BGLP + PRODUCT((/(SD(K,LEX_NUM(K)),K=1,N)/))
  END IF
  IF (ALL(LEX_NUM == COVER_SIZES)) EXIT
END DO MAIN_LOOP
RETURN
END SUBROUTINE BEZOUT_GLP

                                                                     !!!
SUBROUTINE SINGSYS_GLP(N,LEX_NUM,LEX_SAVE,TOL,RAND_MAT,MAT,NONSING)
!
! SINGSYS_GLP determines if the subsystem of the start system
! corresponding to the lexicographic vector LEX_NUM is nonsingular,
! or if a family of subsystems of the start system defined by
! LEX_NUM and LEX_SAVE is singular, by using Householder reflections and
! tree pruning.  Using the notation defined in the module GLOBAL_GLP,
! the vector LEX_NUM defines a linear system of equations
!     L(1,LEX_NUM(1)) = constant_1
!                     .
!                     .
!                     .
!     L(N,LEX_NUM(N)) = constant_N
! which, if nonsingular for generic coefficients, defines
! PRODUCT((/ (SD(K,LEX_NUM(K)), K=1,N) /)) nonsingular starting points
! for homotopy paths.  Nonsingularity of a generic coefficient matrix is
! checked by computing a QR decomposition of the transpose of the
! coefficient matrix.  Observe that if the first J rows are rank
! deficient, then all lexicographic vectors (LEX_NUM(1:J), *) also
! correspond to singular systems, and thus the tree of all possible
! lexicographic orderings can be pruned.
!
! The QR factorization is maintained as a product of Householder
! reflections, and updated based on the difference between LEX_SAVE
! (the value of LEX_NUM returned from the previous call to SINGSYS_GLP)
! and the current input LEX_NUM.  LEX_SAVE and LEX_NUM together
! implicitly define a family of subsystems, namely, all those
! corresponding to lexicographic orderings with head LEX_NUM(1:J),
! where J is the smallest index such that LEX_SAVE(J) /= LEX_NUM(J).
!
! Calls LAPACK subroutines DLARFX and DLARFG.
!
! On input:
!
! N is the dimension of the start and target systems.
!
! LEX_NUM(1:N) is a lexicographic vector which specifies a particular
!   subsystem (and with LEX_SAVE a family of subsystems) of the start
!   system.
!
! LEX_SAVE(1:N) holds the value of LEX_NUM returned from the previous
!   call, and should not be changed between calls to SINGSYS_GLP.  Set
!   LEX_SAVE=0 on the first call to SINGSYS_GLP.
!
! TOL is the singularity test threshold.  The family of subsystems
!   corresponding to lexicographic vectors (LEX_NUM(1:J), *) is declared
!   singular if ABS(R(J,J)) < TOL for the QR factorization of a generic
!   start system coefficient matrix.
!
! RAND_MAT(N,N) is a random matrix with entries uniformly distributed
!   in [-1,-1/2] union [1/2,1], used to seed the random generic
!   coefficient matrix MAT.  RAND_MAT should not change between calls to
!   SINGSYS_GLP.
!
! On output:
!
! LEX_NUM is unchanged if NONSING=.TRUE.  If NONSING=.FALSE.,
!   LEX_NUM(1:J) is unchanged, and
!   LEX_NUM(J+1:N) = COVER_SIZES(J+1:N), where J is the smallest
!   index such that ABS(R(J,J)) < TOL for the QR factorization of the
!   generic start system coefficient matrix corresponding to LEX_NUM
!   (on input).
!
! LEX_SAVE = LEX_NUM.
!
! NONSING = .TRUE. if the subsystem of the start system defined by
!   LEX_NUM is nonsingular.  NONSING = .FALSE. otherwise, which means that
!   the entire family of subsystems corresponding to lexicographic vectors
!   (LEX_NUM(1:J), *) is singular, where J is the smallest index such that
!   ABS(R(J,J)) < TOL for the QR factorization of the generic start system
!   coefficient matrix corresponding to LEX_NUM (on input).
!
! Working storage:
!
! MAT(N+1,N) is updated on successive calls to SINGSYS_GLP, and should
! not be changed by the calling program.  MAT can be undefined on the
! first call to SINGSYS_GLP (when LEX_SAVE = 0).  Define J as the
! smallest index where LEX_SAVE(J) /= LEX_NUM(J).  Upon exit after a
! subsequent call, for some M >= J, MAT contains, in the first M columns,
! a partial QR factorization stored as a product of Householder
! reflections, and, in the last N-M columns, random numbers that define
! the subsystem of the start system corresponding to the lexicographic
! vector LEX_NUM.  For 1<=K<=M, V(2:N+1-K)=MAT(K+1:N,K), V(1)=1, together
! with TAU=MAT(N+1,K), define a Householder reflection of dimension
! N+1-K.

USE GLOBAL_GLP

IMPLICIT NONE
INTEGER, INTENT(IN):: N
INTEGER, DIMENSION(N), INTENT(IN OUT):: LEX_NUM, LEX_SAVE
REAL (KIND=R8), INTENT(IN):: TOL
REAL (KIND=R8), DIMENSION(N,N), INTENT(IN):: RAND_MAT
REAL (KIND=R8), DIMENSION(N+1,N), INTENT(IN OUT):: MAT
LOGICAL, INTENT(OUT):: NONSING

! Local variables.
INTEGER:: I, J, K
REAL (KIND=R8), DIMENSION(N):: V
REAL (KIND=R8):: WORK(1)

IF (N == 1) THEN
  LEX_SAVE = LEX_NUM
  NONSING = .TRUE.
  RETURN
END IF

! (Re)set MAT (in column form) from LEX_NUM.
DO I=1,N
  IF (LEX_SAVE(I) /= LEX_NUM(I)) THEN
    LEX_SAVE(I+1:N) = 0
    DO K=I,N
      MAT(1:N+1,K) = 0.0_R8
      DO J=1,NUMV(K,LEX_NUM(K))
        MAT(PAR(K,LEX_NUM(K),J),K) = RAND_MAT(PAR(K,LEX_NUM(K),J),K)
      END DO
    END DO
    EXIT
  END IF
END DO

! Recompute QR factorization of MAT starting where first change in
! LEX_NUM occurred.
NONSING = .FALSE.
IF (LEX_SAVE(1) /= LEX_NUM(1)) THEN
  ! Skip QR factorization and prune tree if this set degree = 0.
  IF (SD(1,LEX_NUM(1)) == 0) THEN
    LEX_NUM(2:N) = COVER_SIZES(2:N)
    LEX_SAVE = LEX_NUM
    RETURN
  ELSE
    CALL DLARFG(N,MAT(1,1),MAT(2:N,1),1,MAT(N+1,1))
  END IF
END IF
DO J=2,N
  IF (LEX_SAVE(J) /= LEX_NUM(J)) THEN

    ! Skip rest of QR factorization and prune tree if this set degree = 0.
    IF (SD(J,LEX_NUM(J)) == 0) THEN
      IF (J < N) LEX_NUM(J+1:N) = COVER_SIZES(J+1:N)
      EXIT
    END IF
    DO K=1,J-1
      V(K) = 1.0_R8
      V(K+1:N) = MAT(K+1:N,K)
      CALL DLARFX('L',N-K+1,1,V(K:N),MAT(N+1,K),MAT(K:N,J),N-K+1,WORK)
    END DO
    IF (J < N) CALL DLARFG(N-J+1,MAT(J,J),MAT(J+1:N,J),1,MAT(N+1,J))

    ! Check singularity of subsystem corresponding to lexicographic
    ! vector (LEX_NUM(1:J), *).
    IF (ABS(MAT(J,J)) < TOL) THEN
      IF (J < N) LEX_NUM(J+1:N) = COVER_SIZES(J+1:N)
      EXIT
    END IF
  END IF

  ! Subsystem corresponding to LEX_NUM is nonsingular when J==N here.
  IF (J == N) NONSING = .TRUE.
END DO

! Save updated LEX_NUM for next call.
LEX_SAVE = LEX_NUM
RETURN
END SUBROUTINE SINGSYS_GLP

END MODULE POLSYS2
                                                                     !!!


! ----------------------------------------------------------------------
!
! The following modules and external subroutines are from HOMPACK90.


!  This module provides global allocatable arrays used for the sparse
!  matrix data structures, and by the polynomial system solver.  The
!  MODULE HOMOTOPY uses this module.
!
      MODULE HOMPACK90_GLOBAL
      USE REAL_PRECISION
      INTEGER, DIMENSION(:), ALLOCATABLE:: COLPOS, IPAR, ROWPOS
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE:: PAR, PP, QRSPARSE
      END MODULE HOMPACK90_GLOBAL


      MODULE HOMOTOPY       ! Interfaces for user written subroutines.
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL
!
! Interface for subroutine that evaluates F(X) and returns it in the vector V.
      INTERFACE
         SUBROUTINE F(X,V)
         USE REAL_PRECISION
         REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
         REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
         END SUBROUTINE F
      END INTERFACE
!
! Interface for subroutine that returns in V the K-th column of the Jacobian 
! matrix of F(X) evaluated at X. 
      INTERFACE
         SUBROUTINE FJAC(X,V,K)
         USE REAL_PRECISION
         REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
         REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
         INTEGER, INTENT(IN):: K
         END SUBROUTINE FJAC
      END INTERFACE
!
! Interface for subroutine that evaluates RHO(A,LAMBDA,X) and returns it 
! in the vector V.
      INTERFACE
         SUBROUTINE RHO(A,LAMBDA,X,V)
         USE REAL_PRECISION
         REAL (KIND=R8), INTENT(IN):: A(:),X(:)
         REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
         REAL (KIND=R8), INTENT(OUT):: V(:)
         END SUBROUTINE RHO
      END INTERFACE
! The following code is specifically for the polynomial system driver
! POLSYS1H, and should be used verbatim with POLSYS1H in the external 
! subroutine RHO.  
!     USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR  ! FOR POLSYS1H ONLY.
!     INTERFACE
!       SUBROUTINE HFUNP(N,A,LAMBDA,X)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN):: N
!       REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
!       END SUBROUTINE HFUNP
!     END INTERFACE
!     INTEGER:: J,NPOL
! FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
!     IF (LAMBDA .LT. 0.0) LAMBDA=0.0
!     NPOL=IPAR(1)
!     CALL HFUNP(NPOL,A,LAMBDA,X)
!     DO J=1,2*NPOL
!       V(J)=PAR(IPAR(3 + (4-1)) + (J-1))
!     END DO
!     RETURN
! If calling FIXP?? or STEP?? directly, supply appropriate replacement
! code in the external subroutine RHO.
!
! Interface for subroutine that calculates and returns in A the vector
! Z such that RHO(Z,LAMBDA,X) = 0 .
      INTERFACE
         SUBROUTINE RHOA(A,LAMBDA,X)
         USE REAL_PRECISION
         REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: A
         REAL (KIND=R8), INTENT(IN):: LAMBDA,X(:)
         END SUBROUTINE RHOA
      END INTERFACE
!
! Interface for subroutine that returns in the vector V the Kth column
! of the Jacobian matrix [D RHO/D LAMBDA, D RHO/DX] evaluated at the
! point (A, LAMBDA, X).
      INTERFACE
         SUBROUTINE RHOJAC(A,LAMBDA,X,V,K)
         USE REAL_PRECISION
         REAL (KIND=R8), INTENT(IN):: A(:),X(:)
         REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
         REAL (KIND=R8), INTENT(OUT):: V(:)
         INTEGER, INTENT(IN):: K
         END SUBROUTINE RHOJAC
      END INTERFACE
! The following code is specifically for the polynomial system driver
! POLSYS1H, and should be used verbatim with POLSYS1H in the external 
! subroutine RHOJAC.  
!     USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR  ! FOR POLSYS1H ONLY.
!     INTERFACE
!       SUBROUTINE HFUNP(N,A,LAMBDA,X)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN):: N
!       REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
!       END SUBROUTINE HFUNP
!     END INTERFACE
!     INTEGER:: J,NPOL,N2
!     NPOL=IPAR(1)
!     N2=2*NPOL
!     IF (K .EQ. 1) THEN
! FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
!       IF (LAMBDA .LT. 0.0) LAMBDA=0.0
!       CALL HFUNP(NPOL,A,LAMBDA,X)
!       DO J=1,N2
!         V(J)=PAR(IPAR(3 + (6-1)) + (J-1))
!       END DO
!       RETURN
!     ELSE
!       DO J=1,N2
!         V(J)=PAR(IPAR(3 + (5-1)) + (J-1) + N2*(K-2))
!       END DO
!     ENDIF
!
!     RETURN
! If calling FIXP?? or STEP?? directly, supply appropriate replacement
! code in the external subroutine RHOJAC.
!
!
! Interface for subroutine that evaluates a sparse Jacobian matrix of
! F(X) at X, and operates as follows:
!
! If MODE = 1,
! evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
! the result in packed skyline storage format in QRSPARSE.  LENQR is the
! length of QRSPARSE, and ROWPOS contains the indices of the diagonal
! elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
! ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPOS
! is not used by this storage format.
!
! If MODE = 2,
! evaluate the N x N Jacobian matrix of F(X) at X, and return the result
! in sparse row storage format in QRSPARSE.  LENQR is the length of
! QRSPARSE, ROWPOS contains the indices of where each row begins within
! QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
! the corresponding elements in QRSPARSE.  Even if zero, the diagonal
! elements of the Jacobian matrix must be stored in QRSPARSE.
      INTERFACE
         SUBROUTINE FJACS(X)
         USE REAL_PRECISION
         USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
         REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
         END SUBROUTINE FJACS
      END INTERFACE
!
!
! Interface for subroutine that evaluates a sparse Jacobian matrix of
! RHO(A,X,LAMBDA) at (A,X,LAMBDA), and operates as follows:
!
! If MODE = 1,
! evaluate the N X N symmetric Jacobian matrix [D RHO/DX] at
! (A,X,LAMBDA), and return the result in packed skyline storage format in
! QRSPARSE.  LENQR is the length of QRSPARSE, and ROWPOS contains the
! indices of the diagonal elements of [D RHO/DX] within QRSPARSE.  PP
! contains -[D RHO/D LAMBDA] evaluated at (A,X,LAMBDA).  Note the minus
! sign in the definition of PP.  The allocatable array COLPOS is not used
! in this storage format.
!
! If MODE = 2,
! evaluate the N X (N+1) Jacobian matrix [D RHO/DX, D RHO/DLAMBDA] at
! (A,X,LAMBDA), and return the result in sparse row storage format in
! QRSPARSE.  LENQR is the length of QRSPARSE, ROWPOS contains the indices
! of where each row begins within QRSPARSE, and COLPOS (of length LENQR)
! contains the column indices of the corresponding elements in QRSPARSE.
! Even if zero, the diagonal elements of the Jacobian matrix must be
! stored in QRSPARSE.  The allocatable array PP is not used in this
! storage format.
!
      INTERFACE
         SUBROUTINE RHOJS(A,LAMBDA,X)
         USE REAL_PRECISION
         USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
         REAL (KIND=R8), INTENT(IN):: A(:),LAMBDA,X(:)
         END SUBROUTINE RHOJS
      END INTERFACE
      END MODULE HOMOTOPY

      SUBROUTINE STEPNX(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,   &
         ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)
!
!  STEPNX  takes one step along the zero curve of the homotopy map
! using a predictor-corrector algorithm.  The predictor uses a Hermite
! cubic interpolant, and the corrector returns to the zero curve along
! the flow normal to the Davidenko flow.  STEPNX  also estimates a
! step size H for the next step along the zero curve.  STEPNX  is an
! expert user version of STEPN(F|S), written using the reverse call
! protocol.  All matrix data structures and numerical linear algebra
! are the responsibility of the calling program.  STEPNX  indicates to
! the calling program, via flags, at which points  RHO(A,LAMBDA,X)  and
! [ D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX ]  must be
! evaluated, and what linear algebra must be done with these functions.
! Out of range arguments can also be signaled to  STEPNX , which will
! attempt to modify its steplength algorithm to reflect this
! information.
!
! The following interface block should be inserted in the calling
! program:
!
!     INTERFACE
!       SUBROUTINE STEPNX(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
!    &    ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)
!       USE HOMOTOPY
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN):: N
!       INTEGER, INTENT(IN OUT):: NFE,IFLAG
!       LOGICAL, INTENT(IN OUT):: START,CRASH
!       REAL (KIND=R8), INTENT(IN OUT):: HOLD,H,RELERR,ABSERR,S,RHOLEN,
!    &    SSPAR(8)
!       REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
!       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD,
!    &    TZ,W,WP
!       REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: Z0,Z1
!       END SUBROUTINE STEPNX
!     END INTERFACE
!
! ON INPUT:
!
! N = dimension of X and the homotopy map.
!
! NFE = number of Jacobian matrix evaluations.
!
! IFLAG = -2, -1, or 0, indicating the problem type, on the first
!         call to  STEPNX .  STEPNX  does not distinguish between
!         these values, but they are permitted for consistency with
!         the rest of HOMPACK.
!
!       = 0-10*R, -1-10*R, or -2-10*R, R = 1,2,3, indicate to  STEPNX
!         where to resume after a reverse call.  The calling program
!         must not modify  IFLAG  after a reverse call, except as
!         noted next.
!
!       = -40, -41, or -42, used for a final call to deallocate working
!         storage, after all path tracking is finished.  START  and
!         IFLAG  are reset on return.
!
!       = -100-10*R, -101-10*R, -102-10*R, R = 1,2,3, indicate to
!         STEPNX  where to resume after a reverse call, and that the
!         requested evaluation point was out of range.  STEPNX  will
!         reduce  H  and try again.
!
! START = .TRUE. on first call to  STEPNX , .FALSE. otherwise.
!
! HOLD = ||Y - YOLD||; should not be modified by the user.
!
! H = upper limit on length of step that will be attempted.  H  must be
!    set to a positive number on the first call to  STEPNX .
!    Thereafter  STEPNX  calculates an optimal value for  H , and  H
!    should not be modified by the user.
!
! RELERR, ABSERR = relative and absolute error values.  The iteration is
!    considered to have converged when a point W=(LAMBDA,X) is found 
!    such that
!
!    ||Z|| <= RELERR*||W|| + ABSERR  ,          where
!
!    Z is the Newton step to W=(LAMBDA,X).
!
! S = (approximate) arc length along the homotopy zero curve up to
!    Y(S) = (LAMBDA(S), X(S)).
!
! Y(1:N+1) = previous point (LAMBDA(S), X(S)) found on the zero curve of 
!    the homotopy map.
!
! YP(1:N+1) = unit tangent vector to the zero curve of the homotopy map
!    at  Y .
!
! YOLD(1:N+1) = a point before  Y  on the zero curve of the homotopy map.
!
! YPOLD(1:N+1) = unit tangent vector to the zero curve of the homotopy
!    map at  YOLD .
!
! A(:) = parameter vector in the homotopy map.
!
! TZ(1:N+1), W(1:N+1), and WP(1:N+1)  are work arrays used for the
!    Newton step calculation and the interpolation.  On reentry after
!    a reverse call,  WP  and  TZ  contain the tangent vector and
!    Newton step, respectively, at the point  W .  Precisely,
!    D RHO(A,W)/DW WP = 0,  WP^T YP > 0,  ||WP|| = 1,
!    and  TZ  is the minimum norm solution of
!    D RHO(A,W)/DW TZ = - RHO(A,W).
!
! RHOLEN = ||RHO(A,W)||_2 is required by some reverse calls.
!
! SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
!    a vector of parameters used for the optimal step size estimation.
!    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
!    by  STEPNX .  Otherwise the input value of  SSPAR(J)  is used.
!    See the comments below in  STEPNX  for more information about
!    these constants.
!
!
! ON OUTPUT:
!
! N  and  A  are unchanged.
!
! NFE  has been updated.
!
! IFLAG  
!    = -22, -21, -20, -32, -31, or -30 requests the calling program to
!      return the unit tangent vector in  WP , the normal flow Newton
!      step in  TZ , and the 2-norm of the homotopy map in  RHOLEN ,
!      all evaluated at the point  W .
!
!    = -12, -11, or -10 requests the calling program to return in  WP
!      the unit tangent vector at  W .
!
!    = -2, -1, or 0 (unchanged) on a normal return after a successful
!      step.
!
!    = 4 if a Jacobian matrix with rank < N has occurred.  The
!        iteration was not completed.
!
!    = 6 if the iteration failed to converge.  W  contains the last
!        Newton iterate.
!
!    = 7 if input arguments or array sizes are invalid, or  IFLAG  was
!        changed during a reverse call.
!
! START = .FALSE. on a normal return.
!
! CRASH 
!    = .FALSE. on a normal return.
!
!    = .TRUE. if the step size  H  was too small.  H  has been
!      increased to an acceptable value, with which  STEPNX  may be
!      called again.
!
!    = .TRUE. if  RELERR  and/or  ABSERR  were too small.  They have
!      been increased to acceptable values, with which  STEPNX  may
!      be called again.
!
! HOLD = ||Y - YOLD||.
!
! H = optimal value for next step to be attempted.  Normally  H  should
!    not be modified by the user.
!
! RELERR, ABSERR  are unchanged on a normal return.
!
! S = (approximate) arc length along the zero curve of the homotopy map 
!    up to the latest point found, which is returned in  Y .
!
! Y, YP, YOLD, YPOLD  contain the two most recent points and tangent
!    vectors found on the zero curve of the homotopy map.
!
! SSPAR  may have been changed to default values.
!
!
! Z0(1:N+1), Z1(1:N+1)  are allocatable work arrays used for the
!    estimation of the next step size  H .
!
! Calls  DNRM2 .
!
      USE HOMOTOPY
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: N
      INTEGER, INTENT(IN OUT):: NFE,IFLAG
      LOGICAL, INTENT(IN OUT):: START,CRASH
      REAL (KIND=R8), INTENT(IN OUT):: HOLD,H,RELERR,ABSERR,S,RHOLEN,   &
        SSPAR(8)
      REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD,   &
        TZ,W,WP
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: Z0,Z1
!
! ***** LOCAL VARIABLES. *****
!
      REAL (KIND=R8), SAVE:: DCALC,DELS,F0,F1,FOURU,FP0,FP1,   &
        HFAIL,HT,LCALC,RCALC,TEMP,TWOU
      INTEGER, SAVE:: IFLAGC,ITNUM,J,JUDY,NP1
      LOGICAL, SAVE:: FAIL
!
! ***** END OF SPECIFICATION INFORMATION. *****
!
! THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
! THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
! STATEMENT:
      INTEGER, PARAMETER:: LITFH=4
!
! DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
!
      REAL (KIND=R8):: DD001,DD0011,DD01,DD011,DNRM2,QOFS
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) -    &
                                  DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +   &
         DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
!
!
      NP1=N+1
      IF (IFLAG > 0) RETURN
      IF ((START .AND. IFLAG < -2) .OR. SIZE(Y) /= NP1 .OR.   &
        SIZE(YP) /= NP1 .OR. SIZE(YOLD) /= NP1 .OR.   &
        SIZE(YPOLD) /= NP1 .OR. SIZE(TZ) /= NP1 .OR.   &
        SIZE(W) /= NP1 .OR. SIZE(WP) /= NP1 .OR.   &
        (.NOT. START .AND. -MOD(-IFLAG,100) /= IFLAGC .AND.   &
        ABS(IFLAG)/10 /= 4)) THEN
        IFLAG=7
        RETURN
      ENDIF
      IFLAGC=-MOD(-IFLAG,10)
!
! PICK UP EXECUTION WEHRE IT LEFT OFF AFTER A REVERSE CALL.
!
      IF (IFLAG < -2) THEN
        GO TO (50,100,400,700), MOD(ABS(IFLAG),100)/10
      ENDIF
      TWOU=2.0*EPSILON(1.0_R8)
      FOURU=TWOU+TWOU
      CRASH=.TRUE.
! THE ARCLENGTH  S  MUST BE NONNEGATIVE.
      IF (S .LT. 0.0) RETURN
! IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
      IF (H .LT. FOURU*(1.0+S)) THEN
        H=FOURU*(1.0+S)
        RETURN
      ENDIF
! IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE VALUES.
      TEMP=DNRM2(NP1,Y,1)+1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. 0.0) THEN
          RELERR=FOURU*(1.0+FOURU)
          ABSERR=MAX(ABSERR,0.0_R8)
        ELSE
          ABSERR=FOURU*TEMP
        ENDIF
        RETURN
      ENDIF
      CRASH=.FALSE.
      IF (.NOT. START) GO TO 300
!
! *****  STARTUP SECTION (FIRST STEP ALONG ZERO CURVE).  *****
!
      FAIL=.FALSE.
      START=.FALSE.
      IF (ALLOCATED(Z0)) DEALLOCATE(Z0)
      IF (ALLOCATED(Z1)) DEALLOCATE(Z1)
      ALLOCATE(Z0(NP1),Z1(NP1))
!
! SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
! LET Z[K] DENOTE THE NEWTON ITERATES ALONG THE FLOW NORMAL TO THE
! DAVIDENKO FLOW AND Y THEIR LIMIT.
! IDEAL CONTRACTION FACTOR:  ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)= .5
! IDEAL RESIDUAL FACTOR:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= .01
! IDEAL DISTANCE FACTOR:  ||Z[1] - Y|| / ||Z[0] - Y||
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= .5
! MINIMUM STEP SIZE  HMIN .
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_R8)
! MAXIMUM STEP SIZE  HMAX .
      IF (SSPAR(5) .LE. 0.0) SSPAR(5)= 1.0
! MINIMUM STEP SIZE REDUCTION FACTOR  BMIN .
      IF (SSPAR(6) .LE. 0.0) SSPAR(6)= .1_R8
! MAXIMUM STEP SIZE EXPANSION FACTOR  BMAX .
      IF (SSPAR(7) .LE. 0.0) SSPAR(7)= 3.0
! ASSUMED OPERATING ORDER  P .
      IF (SSPAR(8) .LE. 0.0) SSPAR(8)= 2.0
!
! DETERMINE SUITABLE INITIAL STEP SIZE.
      H=MIN(H, .10_R8, SQRT(SQRT(RELERR*TEMP+ABSERR)))
! USE LINEAR PREDICTOR ALONG TANGENT DIRECTION TO START NEWTON ITERATION.
      YPOLD(1)=1.0
      YPOLD(2:NP1)=0.0
! REQUEST TANGENT VECTOR AT Y VIA REVERSE CALL.
      W=Y
      YP=YPOLD
      IFLAG=IFLAGC-10
      IFLAGC=IFLAG
      NFE=NFE+1
      RETURN
 50   YP=WP
! IF THE STARTING POINT IS OUT OF RANGE, GIVE UP.
      IF (IFLAG .LE. -100) THEN
        IFLAG=6
        RETURN
      ENDIF
 70   W=Y + H*YP
      Z0=W
      JUDY=1                                    ! DO JUDY=1,LITFH
 80   IF (JUDY > LITFH) GO TO 200
! REQUEST THE CALCULATION OF THE NEWTON STEP  TZ  AT THE CURRENT
! POINT  W  VIA REVERSE CALL.
        IFLAG=IFLAGC-20
        IFLAGC=IFLAG
        NFE=NFE+1
        RETURN
100     IF (IFLAG .LE. -100) GO TO 200
!
! TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
! COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
! GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)   &
                                                       GO TO 600
!
        JUDY=JUDY+1
      GO TO 80                                   ! END DO
!
! NO CONVERGENCE IN  LITFH  ITERATIONS.  REDUCE  H  AND TRY AGAIN.
200   IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      GO TO 70
!
! ***** END OF STARTUP SECTION. *****
!
! ***** PREDICTOR SECTION. *****
!
300   FAIL=.FALSE.
! COMPUTE POINT PREDICTED BY HERMITE INTERPOLANT.  USE STEP SIZE  H
! COMPUTED ON LAST CALL TO  STEPNX .
320   DO J=1,NP1
        W(J)=QOFS(YOLD(J),YPOLD(J),Y(J),YP(J),HOLD,HOLD+H)
      END DO
      Z0=W 
!
! ***** END OF PREDICTOR SECTION. *****
!
! ***** CORRECTOR SECTION. *****
!
      JUDY=1                          ! CORRECTOR: DO JUDY=1,LITFH
350   IF (JUDY > LITFH) GO TO 500
! REQUEST THE CALCULATION OF THE NEWTON STEP  TZ  AT THE CURRENT
! POINT  W  VIA REVERSE CALL.
        IFLAG=IFLAGC-30
        IFLAGC=IFLAG
        NFE=NFE+1
        RETURN
400     IF (IFLAG .LE. -100) GO TO 500
!
! TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
! COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
! GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)   &
                                                       GO TO 600
!
        JUDY=JUDY+1
      GO TO 350                              ! END DO CORRECTOR
!
! NO CONVERGENCE IN  LITFH  ITERATIONS.  RECORD FAILURE AT CALCULATED  H , 
! SAVE THIS STEP SIZE, REDUCE  H  AND TRY AGAIN.
500   FAIL=.TRUE.
      HFAIL=H
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      GO TO 320
!
! ***** END OF CORRECTOR SECTION. *****
!
! ***** MOP-UP SECTION. *****
!
! YOLD  AND  Y  ALWAYS CONTAIN THE LAST TWO POINTS FOUND ON THE ZERO
! CURVE OF THE HOMOTOPY MAP.  YPOLD  AND  YP  CONTAIN THE TANGENT
! VECTORS TO THE ZERO CURVE AT  YOLD  AND  Y , RESPECTIVELY.
!
600   YPOLD=YP
      YOLD=Y
      Y=W
      YP=WP
      W=Y - YOLD
! UPDATE ARC LENGTH.
      HOLD=DNRM2(NP1,W,1)
      S=S+HOLD
!
! ***** END OF MOP-UP SECTION. *****
!
! ***** OPTIMAL STEP SIZE ESTIMATION SECTION. *****
!
! CALCULATE THE DISTANCE FACTOR  DCALC .
      TZ=Z0 - Y
      W=Z1 - Y
      DCALC=DNRM2(NP1,TZ,1)
      IF (DCALC .NE. 0.0) DCALC=DNRM2(NP1,W,1)/DCALC
!
! THE OPTIMAL STEP SIZE HBAR IS DEFINED BY
!
!   HT=HOLD * [MIN(LIDEAL/LCALC, RIDEAL/RCALC, DIDEAL/DCALC)]**(1/P)
!
!     HBAR = MIN [ MAX(HT, BMIN*HOLD, HMIN), BMAX*HOLD, HMAX ]
!
! IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, SET THE CONTRACTION
! FACTOR  LCALC  TO ZERO.
      IF (ITNUM .EQ. 1) LCALC = 0.0
! FORMULA FOR OPTIMAL STEP SIZE.
      IF (LCALC+RCALC+DCALC .EQ. 0.0) THEN
        HT = SSPAR(7) * HOLD
      ELSE 
        HT = (1.0/MAX(LCALC/SSPAR(1), RCALC/SSPAR(2), DCALC/SSPAR(3)))   &
             **(1.0/SSPAR(8)) * HOLD
      ENDIF
!  HT  CONTAINS THE ESTIMATED OPTIMAL STEP SIZE.  NOW PUT IT WITHIN
! REASONABLE BOUNDS.
      H=MIN(MAX(HT,SSPAR(6)*HOLD,SSPAR(4)), SSPAR(7)*HOLD, SSPAR(5))
      IF (ITNUM .EQ. 1) THEN
! IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, DON'T DECREASE  H .
        H=MAX(H,HOLD)
      ELSE IF (ITNUM .EQ. LITFH) THEN
! IF CONVERGENCE REQUIRED THE MAXIMUM  LITFH  ITERATIONS, DON'T
! INCREASE  H .
        H=MIN(H,HOLD)
      ENDIF
! IF CONVERGENCE DID NOT OCCUR IN  LITFH  ITERATIONS FOR A PARTICULAR
! H = HFAIL , DON'T CHOOSE THE NEW STEP SIZE LARGER THAN  HFAIL .
      IF (FAIL) H=MIN(H,HFAIL)
!
!
      IFLAG=IFLAGC
      RETURN
! CLEAN UP ALLOCATED WORKING STORAGE.
 700  START=.TRUE.
      IFLAG=IFLAGC
      IF (ALLOCATED(Z0)) DEALLOCATE(Z0)
      IF (ALLOCATED(Z1)) DEALLOCATE(Z1)
      RETURN
      END SUBROUTINE STEPNX
