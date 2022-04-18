      SUBROUTINE DHABAK( JOB, SGN, N, ILO, SCALE, M, V1, LDV1, V2, LDV2,
     $                   INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     DHABAK applies the inverse of a balancing transformation, computed
C     by the routines DHABAL or DSHBAL, to a 2*N-by-M matrix
C
C               [   V1   ]
C               [        ],
C               [ op(V2) ]
C
C     where op(V2) = +/-V2.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the type of inverse transformation required:
C             = 'N': do nothing, return immediately;
C             = 'P': do inverse transformation for permutation only;
C             = 'S': do inverse transformation for scaling only;
C             = 'B': do inverse transformations for both permutation and
C                    scaling.
C             JOB must be the same as the argument JOB supplied to
C             DHABAL or DSHBAL.
C
C     SGN     CHARACTER*1
C             Specifies the sign of op(V2):
C             = 'P': op(V2) =  V2;
C             = 'N': op(V2) = -V2.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of rows of the matrices V1 and V2. N >= 0.
C
C     ILO     (input) INTEGER
C             The integer ILO determined by DHABAL or DSHBAL.
C             1 <= ILO <= N+1.
C
C     SCALE   (input) DOUBLE PRECISION array, dimension (N)
C             Details of the permutation and scaling factors, as
C             returned by DHABAL or DSHBAL.
C
C     M       (input) INTEGER
C             The number of columns of the matrices V1 and V2.  M >= 0.
C
C     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,M)
C             On entry, the leading N-by-M part of this array must
C             contain the matrix V1.
C             On exit, the leading N-by-N part of this array is
C             overwritten by the updated matrix V1 of the transformed
C             matrix.
C
C     LDV1    (input) INTEGER
C             The leading dimension of the array V1. LDV1 >= max(1,N).
C
C     V2      (input/output) DOUBLE PRECISION array, dimension (LDV2,M)
C             On entry, the leading N-by-M part of this array must
C             contain the matrix V2.
C             On exit, the leading N-by-N part of this array is
C             overwritten by the updated matrix V2 of the transformed
C             matrix.
C
C     LDV2    (input) INTEGER
C             The leading dimension of the array V2. LDV2 >= max(1,N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] P. BENNER:
C     Symplectic balancing of Hamiltonian matrices.
C     SIAM J. Sci. Comput., 22(5):1885--1904, 2000.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     KEYWORDS
C
C     Balancing, Hamiltonian matrix, skew-Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER*1       JOB, SGN
      INTEGER           ILO, INFO, LDV1, LDV2, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  SCALE(*), V1(LDV1,*), V2(LDV2,*)
C     .. Local Scalars ..
      LOGICAL           LPERM, LSCAL, LSGN, SYSW
      INTEGER           I, K
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DRSCL, DSCAL, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      LPERM = LSAME( JOB, 'P' ).OR.LSAME( JOB, 'B' )
      LSCAL = LSAME( JOB, 'S' ).OR.LSAME( JOB, 'B' )
      LSGN  = LSAME( SGN, 'N' )
      IF ( .NOT.LPERM .AND. .NOT.LSCAL
     $     .AND. .NOT.LSAME( JOB, 'N' ) ) THEN
         INFO = -1
      ELSE IF ( .NOT. LSGN .AND. .NOT.LSAME( SGN, 'P' ) ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.N+1 ) THEN
         INFO = -4
      ELSE IF ( M.LT.0 ) THEN
         INFO = -6
      ELSE IF ( LDV1.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF ( LDV2.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
C
C     Return if there were illegal values.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DHABAK', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF ( N.EQ.0 .OR. M.EQ.0 .OR. LSAME( JOB, 'N' ) )
     $   RETURN
C
C     Inverse scaling.
C
      IF ( LSCAL ) THEN
         DO 20 I = ILO, N
            CALL DRSCL( M, SCALE(I), V1(I,1), LDV1 )
   20    CONTINUE
         DO 30 I = ILO, N
            CALL DRSCL( M, SCALE(I), V2(I,1), LDV2 )
   30    CONTINUE
      END IF
C
C     Inverse permutation.
C
      IF ( LPERM ) THEN
         DO 40 I = ILO-1, 1, -1
            K = INT( SCALE( I ) )
            SYSW = ( K.GT.N )
            IF ( SYSW )
     $            K = K - N
C
            IF ( K.NE.I ) THEN
C
C              Exchange rows k <-> i.
C
               CALL DSWAP( M, V1(I,1), LDV1, V1(K,1), LDV1 )
               CALL DSWAP( M, V2(I,1), LDV2, V2(K,1), LDV2 )
            END IF
C
            IF ( SYSW ) THEN
C
C              Exchange V1(k,:) <-> V2(k,:).
C
               CALL DSWAP( M, V1(K,1), LDV1, V2(K,1), LDV2 )
               IF ( LSGN ) THEN
                  CALL DSCAL( M, -ONE, V2(K,1), LDV2 )
               ELSE
                  CALL DSCAL( M, -ONE, V1(K,1), LDV1 )
               END IF
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C *** Last line of DHABAK ***
      END
