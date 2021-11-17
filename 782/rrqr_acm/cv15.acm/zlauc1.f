      LOGICAL FUNCTION ZLAUC1( K, X, SMIN, W, GAMMA, THRESH )
*
*     This code is part of a release of the package for computing
*     rank-revealing QR Factorizations written by:
*     ==================================================================
*     Christian H. Bischof        and   Gregorio Quintana-Orti
*     Math. and Comp. Sci. Div.         Departamento de Informatica
*     Argonne National Lab.             Universidad Jaime I
*     Argonne, IL 60439                 Campus P. Roja, 12071 Castellon
*     USA                               Spain
*     bischof@mcs.anl.gov               gquintan@inf.uji.es
*     ==================================================================
*     $Revision: 1.42 $
*     $Date: 96/12/30 16:59:41 $
*
*     .. Scalar Arguments ..
      INTEGER            K
      DOUBLE PRECISION   SMIN, THRESH
      COMPLEX*16         GAMMA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         W( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PREC_LAUC1 applies incremental condition estimation to determine whether
*  the K-th column of A, stored in vector W, would be acceptable as a pivot
*  column with respect to the threshold THRESH.
*
*  Arguments
*  =========
*
*  K       (input) INTEGER
*          Length of vector X.
*
*  X       (input/output) COMPLEX*16 array, dimension ( K )
*          On entry, X(1:K-1) contains an approximate smallest left singular
*          vector of the upper triangle of A(1:k-1,1:k-1).
*          On exit, if ZLAUC1 == .TRUE., X contains an approximate
*          smallest left singular vector of the upper triangle of A(1:k,1:k);
*          if ZLAUC1 == .FALSE., X is unchanged.
*
*  SMIN    (input/output) DOUBLE PRECISION
*          On entry, an estimate for the smallest singular value of the
*          upper triangle of A(1:k-1,1:k-1).
*          On exit, if ZLAUC1 == .TRUE., SMIN is an estimate of the
*          smallest singular value of the upper triangle of  A(1:k,1:k);
*          if ZLAUC1 == .FALSE., SMIN is unchanged.
*
*  W       (input) FLOATING_DECLARE array, dimension ( K-1 )
*          The K-th column of matrix A excluding the diagonal element.
*
*  GAMMA   (input) COMPLEX*16
*          Diagonal entry in k-th column of A if column k were to
*          be  accepted.
*
*  THRESH  (input) DOUBLE PRECISION
*          If the approximate smallest singular value for A(1:K,1:K)
*          is smaller than THRESH, the kth column is rejected.
*
*  (ZLAUC1) (output) LOGICAL
*          If the k-th column of A is found acceptable, ZLAUC1
*          returns .TRUE., otherwise it returns .FALSE.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   SMINPR
      COMPLEX*16         SINE, COSINE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLAIC1, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*
*     Try to use diagonal element as condition estimator
*
      IF( THRESH.GT.ABS( GAMMA ) ) THEN
         ZLAUC1 = .FALSE.
         RETURN
      END IF
*
*     Use incremental condition estimation to determine an estimate
*     SMINPR and an approximate singular vector [SINE*X,COSINE]'
*     for A(K,K).
*
      CALL ZLAIC1( 2, K-1, X, SMIN, W, GAMMA, SMINPR,
     $             SINE, COSINE )
      IF( THRESH.GT.SMINPR ) THEN
         ZLAUC1 =  .FALSE.
      ELSE
         CALL ZSCAL( K-1, SINE, X, 1 )
         X( K ) = COSINE
         SMIN = SMINPR
         ZLAUC1 = .TRUE.
      END IF
      RETURN
*
*     End of ZLAUC1
*
      END
