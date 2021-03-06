      SUBROUTINE DLAQR6( JOB, WANTT, WANTZ, KACC22, N, KTOP, KBOT, 
     $                   NSHFTS, SR, SI, H, LDH, ILOZ, IHIZ, Z, 
     $                   LDZ, V, LDV, IHV, LDIHV, HV, LDHV, U, LDU, 
     $                   NV, WV, LDWV, NH, WH, LDWH )
*
*  -- ScaLAPACK auxiliary routine (version 1.7.x) --
*     University of Umea and HPC2N, Umea, Sweden
*     February 2008
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
     $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV, LDHV,
     $                   LDIHV
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      INTEGER            IHV( LDIHV, * )
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ),
     $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ),
     $                   Z( LDZ, * ), HV( LDHV, * )
*     ..
*
*     This auxiliary subroutine called by PDLAQR5 performs a
*     single small-bulge multi-shift QR sweep, moving the chain
*     of bulges from top to bottom in the submatrix 
*     H(KTOP:KBOT,KTOP:KBOT), collecting the transformations in the 
*     matrix HV *or* accumulating the transformations in the matrix
*     Z (see below).
*
*     This is a modified version of DLAQR5 from LAPACK 3.1.
*
* ======================================================================
*
*      JOB    (input) character scalar
*             Set the kind of job to do in DLAQR6, as follows:
*             JOB = 'I': Introduce and chase bulges in submatrix
*             JOB = 'C': Chase bulges from top to bottom of submatrix 
*             JOB = 'O': Chase bulges off submatrix
*
*      WANTT  (input) logical scalar
*             WANTT = .true. if the quasi-triangular Schur factor
*             is being computed.  WANTT is set to .false. otherwise.
*
*      WANTZ  (input) logical scalar
*             WANTZ = .true. if the orthogonal Schur factor is being
*             computed.  WANTZ is set to .false. otherwise. 
*
*      KACC22 (input) integer with value 0, 1, or 2.
*             Specifies the computation mode of far-from-diagonal
*             orthogonal updates.
*        = 0: DLAQR6 does not accumulate reflections and does not
*             use matrix-matrix multiply to update far-from-diagonal
*             matrix entries.
*        = 1: DLAQR6 accumulates reflections and uses matrix-matrix
*             multiply to update the far-from-diagonal matrix entries.
*        = 2: DLAQR6 accumulates reflections, uses matrix-matrix
*             multiply to update the far-from-diagonal matrix entries,
*             and takes advantage of 2-by-2 block structure during
*             matrix multiplies.
*
*      N      (input) integer scalar
*             N is the order of the Hessenberg matrix H upon which this
*             subroutine operates.
*
*      KTOP   (input) integer scalar
*      KBOT   (input) integer scalar
*             These are the first and last rows and columns of an
*             isolated diagonal block upon which the QR sweep is to be
*             applied. It is assumed without a check that
*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
*             and
*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
*
*      NSHFTS (input) integer scalar
*             NSHFTS gives the number of simultaneous shifts.  NSHFTS
*             must be positive and even.
*
*      SR     (input) DOUBLE PRECISION array of size (NSHFTS)
*      SI     (input) DOUBLE PRECISION array of size (NSHFTS)
*             SR contains the real parts and SI contains the imaginary
*             parts of the NSHFTS shifts of origin that define the
*             multi-shift QR sweep.
*
*      H      (input/output) DOUBLE PRECISION array of size (LDH,N)
*             On input H contains a Hessenberg matrix.  On output a
*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
*             to the isolated diagonal block in rows and columns KTOP
*             through KBOT.
*
*      LDH    (input) integer scalar
*             LDH is the leading dimension of H just as declared in the
*             calling procedure.  LDH.GE.MAX(1,N).
*
*      ILOZ   (input) INTEGER
*      IHIZ   (input) INTEGER
*             Specify the rows of Z to which transformations must be
*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
*
*      Z      (input/output) DOUBLE PRECISION array of size (LDZ,IHI)
*             If WANTZ = .TRUE., then the QR Sweep orthogonal
*             similarity transformation is accumulated into
*             Z(ILOZ:IHIZ,ILO:IHI) from the right.
*             If WANTZ = .FALSE., then Z is unreferenced.
*
*      LDZ    (input) integer scalar
*             LDA is the leading dimension of Z just as declared in
*             the calling procedure. LDZ.GE.N.
*
*      V      (workspace) DOUBLE PRECISION array of size (LDV,NSHFTS/2)
*
*      LDV    (input) integer scalar
*             LDV is the leading dimension of V as declared in the
*             calling procedure.  LDV.GE.3.
*
*      IHV    (output) INTEGER array of size (LDIHV,KBOT-KTOP+1)
*             If WANTZ = .true., this array is not referenced.
*             Otherwise on output, this array contains sizes and
*             the application indices of the reflectors in the array
*             HV: IHV(1,i) contains the size of the i:th reflector, and
*             k=IHV(2,i) contains the corresponding application position 
*             in H(k,k).
*
*      LDIHV  (input) integer scalar
*             LDIHV is the leading dimension of IHV as declared in the
*             calling procedure.  LDIVH.GE.2.
*
*      HV     (output) DOUBLE PRECISION array of size (LDHV,KBOT-KTOP+1)
*             If WANTZ = .true., this array is not referenced.
*             Otherwise on output, this array contains the 2x2 and 3x3 
*             reflectors defining the bulge-chase from top to bottom of 
*             the considered matrix H(KTOP:KBOT,KTOP:KBOT), where column
*             j in HV corresponds to the j:th applied reflector.
*
*      LDVH   (input) integer scalar
*             LDV is the leading dimension of V as declared in the
*             calling procedure.  LDV.GE.3.
*
*      U      (workspace) DOUBLE PRECISION array of size
*             (LDU,3*NSHFTS-3)
*
*      LDU    (input) integer scalar
*             LDU is the leading dimension of U just as declared in the
*             in the calling subroutine.  LDU.GE.3*NSHFTS-3.
*
*      NH     (input) integer scalar
*             NH is the number of columns in array WH available for
*             workspace. NH.GE.1 is required for usage of this 
*             workspace, otherwise the updates of the far-from-diagonal
*             elements will be updated without level 3 BLAS.
*
*      WH     (workspace) DOUBLE PRECISION array of size (LDWH,NH)
*
*      LDWH   (input) integer scalar
*             Leading dimension of WH just as declared in the
*             calling procedure.  LDWH.GE.3*NSHFTS-3.
*
*      NV     (input) integer scalar
*             NV is the number of rows in WV agailable for workspace.
*             NV.GE.1 is required for usage of this 
*             workspace, otherwise the updates of the far-from-diagonal
*             elements will be updated without level 3 BLAS.
*
*      WV     (workspace) DOUBLE PRECISION array of size
*             (LDWV,3*NSHFTS-3)
*
*      LDWV   (input) integer scalar
*             LDWV is the leading dimension of WV as declared in the
*             in the calling subroutine.  LDWV.GE.NV.
*
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*        Robert Granat, Department of Computing Science and HPC2N,
*        University of Umea, Sweden.
*
*     ============================================================
*     Reference:
*
*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*     Algorithm Part I: Maintaining Well Focused Shifts, and
*     Level 3 Performance, SIAM Journal of Matrix Analysis,
*     volume 23, pages 929--947, 2002.
*
*     ============================================================
*     .. Parameters ..
      LOGICAL            DEBUG, PRINT
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0, DEBUG = .FALSE.,
     $                     PRINT = .FALSE.)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM,
     $                   SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, TST1, TST2,
     $                   ULP
      INTEGER            I, I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN,
     $                   JROW, JTOP, K, K1, KDU, KMS, KNZ, KRCOL, KZS,
     $                   M, M22, MBOT, MEND, MSTART, MTOP, NBMPS, NDCOL,
     $                   NS, NU, SINCOL, EINCOL, UINCOL, IPHV, CHUNK,
     $                   THREADS, JLEN2, JCOL2, GCHUNK, JROW2, MAXCHUNK
      LOGICAL            ACCUM, BLK22, BMP22, INTRO, CHASE, OFF, ALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            OMP_GET_NUM_THREADS, PILAENVX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, OMP_GET_NUM_THREADS,
     $                   PILAENVX
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          ABS, DBLE, MAX, MIN, MOD
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   VT( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLABAD, DLACPY, DLAQR1, DLARFG, DLASET,
     $                   DTRMM
*     ..
*     .. Executable Statements ..
*
      IF( DEBUG ) WRITE(*,*) 'Entering dlaqr6:',JOB,N,KTOP,KBOT,NH,NV
*
*     ==== If there are no shifts, then there is nothing to do. ====
*
      IF( NSHFTS.LT.2 )
     $   RETURN
*
*     ==== If the active block is empty or 1-by-1, then there
*     .    is nothing to do. ====
*
      IF( KTOP.GE.KBOT )
     $   RETURN
*
*     Get chunk size in case of multithreading
*
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
      THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
#else
      THREADS = 1
#endif  
      IF( DEBUG ) WRITE(*,*) '#threads=',THREADS
*
*     ==== Shuffle shifts into pairs of real shifts and pairs
*     .    of complex conjugate shifts assuming complex
*     .    conjugate shifts are already adjacent to one
*     .    another. ====
*
      DO 10 I = 1, NSHFTS - 2, 2
         IF( SI( I ).NE.-SI( I+1 ) ) THEN
*
            SWAP = SR( I )
            SR( I ) = SR( I+1 )
            SR( I+1 ) = SR( I+2 )
            SR( I+2 ) = SWAP
*
            SWAP = SI( I )
            SI( I ) = SI( I+1 )
            SI( I+1 ) = SI( I+2 )
            SI( I+2 ) = SWAP
         END IF
   10 CONTINUE
*
*     ==== NSHFTS is supposed to be even, but if is odd,
*     .    then simply reduce it by one.  The shuffle above
*     .    ensures that the dropped shift is real and that
*     .    the remaining shifts are paired. ====
*
      NS = NSHFTS - MOD( NSHFTS, 2 )
*
*     ==== Machine constants for deflation ====
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries ? This is only performed if both NH and NV is 
*          greater than 1. ====
*
      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
      ACCUM = ACCUM .AND. NH.GE.1 .AND. NV.GE.1
*      ACCUM = .FALSE.
      IF( DEBUG ) WRITE(*,*) 'ACCUM=',ACCUM
*
*     ==== If so, exploit the 2-by-2 block structure? ====
*
      BLK22 = ( NS.GT.2 ) .AND. ( KACC22.EQ.2 )
*
*     ==== Decode JOB ====
*
      ALL = LSAME( JOB, 'A' )
      IF( .NOT. ALL )
     $     INTRO = LSAME( JOB, 'I' )
      IF( .NOT. ALL .AND. .NOT. INTRO )
     $     CHASE = LSAME( JOB, 'C' )
      IF( .NOT. ALL .AND. .NOT. INTRO .AND. .NOT. CHASE ) THEN
         OFF = LSAME( JOB, 'O' )
         IF( .NOT. OFF ) 
     $        RETURN
      END IF
*
*     ==== clear trash ====
*
      IF( INTRO.OR.ALL .AND. KTOP+2.LE.KBOT )
     $   H( KTOP+2, KTOP ) = ZERO
*
*     ==== NBMPS = number of 2-shift bulges in the chain ====
*
      NBMPS = NS / 2
      IF( DEBUG ) WRITE(*,*) 'NBMPS =',NBMPS
*
*     ==== KDU = width of slab ====
*
      KDU = 6*NBMPS - 3
      IF( DEBUG ) WRITE(*,*) 'KDU =',KDU
*
*     Set loop limits for bulge-chasing depending on working mode
*
      IF( ALL ) THEN
         SINCOL = 3*( 1-NBMPS ) + KTOP - 1
         EINCOL = KBOT - 2
         UINCOL = 3*NBMPS - 2
      ELSEIF( INTRO ) THEN
         SINCOL = 3*( 1-NBMPS ) + KTOP - 1
         EINCOL = KBOT - 3*NBMPS - 1
         UINCOL = 3*NBMPS - 2
      ELSEIF( CHASE ) THEN
         SINCOL = KTOP
         EINCOL = KBOT - 3*NBMPS - 1
         UINCOL = 3*NBMPS - 2
      ELSEIF( OFF ) THEN
         SINCOL = KTOP
         EINCOL = KBOT - 2
         UINCOL = 3*NBMPS - 2
      END IF
      IPHV = 0
      IF( DEBUG ) WRITE(*,*) 'SINCOL,EINCOL,UINCOL=',
     $     SINCOL,EINCOL,UINCOL  
*
*     ==== Create and/or chase chains of NBMPS bulges ====
*
      DO 220 INCOL = SINCOL, EINCOL, UINCOL
         IF( DEBUG ) WRITE(*,*) 'INCOL=',INCOL
         NDCOL = MIN( INCOL + KDU, EINCOL )
         IF( DEBUG ) WRITE(*,*) 'NDCOL=',NDCOL
         IF( ACCUM )
     $      CALL DLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
*
*        ==== Near-the-diagonal bulge chase.  The following loop
*        .    performs the near-the-diagonal part of a small bulge
*        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
*        .    chunk extends from column INCOL to column NDCOL
*        .    (including both column INCOL and column NDCOL). The
*        .    following loop chases a 3*NBMPS column long chain of
*        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
*        .    may be less than KTOP and and NDCOL may be greater than
*        .    KBOT indicating phantom columns from which to chase
*        .    bulges before they are actually introduced or to which
*        .    to chase bulges beyond column KBOT.)  ====
*
         DO 150 KRCOL = INCOL, MIN( EINCOL, INCOL+3*NBMPS-3, KBOT-2 )
            IF( DEBUG ) WRITE(*,*) 'KRCOL=',KRCOL
*
*           ==== Bulges number MTOP to MBOT are active double implicit
*           .    shift bulges.  There may or may not also be small
*           .    2-by-2 bulge, if there is room.  The inactive bulges
*           .    (if any) must wait until the active bulges have moved
*           .    down the diagonal to make room.  The phantom matrix
*           .    paradigm described above helps keep track.  ====
*
            MTOP = MAX( 1, ( ( KTOP-1 )-KRCOL+2 ) / 3+1 )
            MBOT = MIN( NBMPS, ( KBOT-KRCOL ) / 3 )
            M22 = MBOT + 1
            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+3*( M22-1 ) ).EQ.
     $              ( KBOT-2 )
*
*           ==== Generate reflections to chase the chain right
*           .    one column.  (The minimum value of K is KTOP-1.) ====
*
            DO 20 M = MTOP, MBOT
               K = KRCOL + 3*( M-1 )
               IF( K.EQ.KTOP-1 ) THEN
                  CALL DLAQR1( 3, H( KTOP, KTOP ), LDH, SR( 2*M-1 ),
     $                         SI( 2*M-1 ), SR( 2*M ), SI( 2*M ),
     $                         V( 1, M ) )
                  ALPHA = V( 1, M )
                  CALL DLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
                  IF( .NOT. WANTZ ) THEN
                     IPHV = IPHV + 1
                     CALL DCOPY( 3, V( 1, M ), 1, HV( 1, IPHV ), 1 )
                     IHV( 1, IPHV ) = 3
                     IHV( 2, IPHV ) = KTOP
                  END IF
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M ) = H( K+2, K )
                  V( 3, M ) = H( K+3, K )
                  CALL DLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
*
*                 ==== A Bulge may collapse because of vigilant
*                 .    deflation or destructive underflow.  (The
*                 .    initial bulge is always collapsed.) Use
*                 .    the two-small-subdiagonals trick to try
*                 .    to get it started again. If V(2,M).NE.0 and
*                 .    V(3,M) = H(K+3,K+1) = H(K+3,K+2) = 0, then
*                 .    this bulge is collapsing into a zero
*                 .    subdiagonal.  It will be restarted next
*                 .    trip through the loop.)
*
                  IF( V( 1, M ).NE.ZERO .AND.
     $                ( V( 3, M ).NE.ZERO .OR. ( H( K+3,
     $                K+1 ).EQ.ZERO .AND. H( K+3, K+2 ).EQ.ZERO ) ) )
     $                 THEN
*
*                    ==== Typical case: not collapsed (yet). ====
*
                     H( K+1, K ) = BETA
                     H( K+2, K ) = ZERO
                     H( K+3, K ) = ZERO
                  ELSE
*
*                    ==== Atypical case: collapsed.  Attempt to
*                    .    reintroduce ignoring H(K+1,K).  If the
*                    .    fill resulting from the new reflector
*                    .    is too large, then abandon it.
*                    .    Otherwise, use the new one. ====
*
                     CALL DLAQR1( 3, H( K+1, K+1 ), LDH, SR( 2*M-1 ),
     $                            SI( 2*M-1 ), SR( 2*M ), SI( 2*M ),
     $                            VT )
                     SCL = ABS( VT( 1 ) ) + ABS( VT( 2 ) ) +
     $                     ABS( VT( 3 ) )
                     IF( SCL.NE.ZERO ) THEN
                        VT( 1 ) = VT( 1 ) / SCL
                        VT( 2 ) = VT( 2 ) / SCL
                        VT( 3 ) = VT( 3 ) / SCL
                     END IF
*
*                    ==== The following is the traditional and
*                    .    conservative two-small-subdiagonals
*                    .    test.  ====
*                    .
                     IF( ABS( H( K+1, K ) )*( ABS( VT( 2 ) )+
     $                   ABS( VT( 3 ) ) ).GT.ULP*ABS( VT( 1 ) )*
     $                   ( ABS( H( K, K ) )+ABS( H( K+1,
     $                   K+1 ) )+ABS( H( K+2, K+2 ) ) ) ) THEN
*
*                       ==== Starting a new bulge here would
*                       .    create non-negligible fill.   If
*                       .    the old reflector is diagonal (only
*                       .    possible with underflows), then
*                       .    change it to I.  Otherwise, use
*                       .    it with trepidation. ====
*
                        IF( V( 2, M ).EQ.ZERO .AND. V( 3, M ).EQ.ZERO )
     $                       THEN
                           V( 1, M ) = ZERO
                        ELSE
                           H( K+1, K ) = BETA
                           H( K+2, K ) = ZERO
                           H( K+3, K ) = ZERO
                        END IF
                     ELSE
*
*                       ==== Stating a new bulge here would
*                       .    create only negligible fill.
*                       .    Replace the old reflector with
*                       .    the new one. ====
*
                        ALPHA = VT( 1 )
                        CALL DLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
                        REFSUM = H( K+1, K ) + H( K+2, K )*VT( 2 ) +
     $                           H( K+3, K )*VT( 3 )
                        H( K+1, K ) = H( K+1, K ) - VT( 1 )*REFSUM
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                        V( 1, M ) = VT( 1 )
                        V( 2, M ) = VT( 2 )
                        V( 3, M ) = VT( 3 )
                     END IF
                  END IF
                  IF( .NOT. WANTZ ) THEN
                     IPHV = IPHV + 1
                     CALL DCOPY( 3, V( 1, M ), 1, HV( 1, IPHV ), 1 )
                     IHV( 1, IPHV ) = 3
                     IHV( 2, IPHV ) = K+1
                  END IF
               END IF
   20       CONTINUE
*
*           ==== Generate a 2-by-2 reflection, if needed. ====
*
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 ) THEN
               IF( K.EQ.KTOP-1 ) THEN
                  CALL DLAQR1( 2, H( K+1, K+1 ), LDH, SR( 2*M22-1 ),
     $                         SI( 2*M22-1 ), SR( 2*M22 ), SI( 2*M22 ),
     $                         V( 1, M22 ) )
                  BETA = V( 1, M22 )
                  CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
                  IF( .NOT. WANTZ ) THEN
                     IPHV = IPHV + 1
                     CALL DCOPY( 2, V( 1, M22 ), 1, HV( 1, IPHV ), 1 )
                     IHV( 1, IPHV ) = 2
                     IHV( 2, IPHV ) = K+1
                  END IF
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M22 ) = H( K+2, K )
                  CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
                  IF( .NOT. WANTZ ) THEN
                     IPHV = IPHV + 1
                     CALL DCOPY( 2, V( 1, M22 ), 1, HV( 1, IPHV ), 1 )
                     IHV( 1, IPHV ) = 2
                     IHV( 2, IPHV ) = K+1
                  END IF
                  H( K+1, K ) = BETA
                  H( K+2, K ) = ZERO
               END IF
            ELSE
*
*              ==== Initialize V(1,M22) here to avoid possible undefined
*              .    variable problems later. ====
*
               V( 1, M22 ) = ZERO
               IF( .NOT. WANTZ ) THEN
                  IPHV = IPHV + 1
                  HV( 1, IPHV ) = ZERO
                  IHV( 1, IPHV ) = 2
                  IHV( 2, IPHV ) = K+1
               END IF
            END IF
*
*           ==== Multiply H by reflections from the left ====
*            
            IF( ACCUM ) THEN
               JBOT = MIN( MAX(INCOL+KDU,NDCOL), KBOT )
               IF( DEBUG ) WRITE(*,*) 'NDCOL,KBOT,JBOT=',
     $              NDCOL,KBOT,JBOT
            ELSE IF( WANTT ) THEN
               JBOT = N
            ELSE
               JBOT = KBOT
            END IF
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, MEND, M, K, REFSUM )
*$OMP DO
#endif
            DO 40 J = MAX( KTOP, KRCOL ), JBOT
               IF( DEBUG ) WRITE(*,*) 'Left mult on H: J=',J
               MEND = MIN( MBOT, ( J-KRCOL+2 ) / 3 )
               DO 30 M = MTOP, MEND
                  K = KRCOL + 3*( M-1 )
                  REFSUM = V( 1, M )*( H( K+1, J )+V( 2, M )*
     $                     H( K+2, J )+V( 3, M )*H( K+3, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M )
                  H( K+3, J ) = H( K+3, J ) - REFSUM*V( 3, M )
 30            CONTINUE
 40         CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
            IF( BMP22 ) THEN
               K = KRCOL + 3*( M22-1 )
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO
#endif
               DO 50 J = MAX( K+1, KTOP ), JBOT
                  REFSUM = V( 1, M22 )*( H( K+1, J )+V( 2, M22 )*
     $                     H( K+2, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M22 )
   50          CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
            END IF
             IF( PRINT )
     $              CALL DLAPRNT(N,N,H,1,1,LDH,'H5',6)
*
*           ==== Multiply H by reflections from the right.
*           .    Delay filling in the last row until the
*           .    vigilant deflation check is complete. ====
*
            IF( ACCUM ) THEN
               JTOP = MAX( KTOP, INCOL )
            ELSE IF( WANTT ) THEN
               JTOP = 1
            ELSE
               JTOP = KTOP
            END IF
            DO 90 M = MTOP, MBOT
               IF( V( 1, M ).NE.ZERO ) THEN
                  K = KRCOL + 3*( M-1 )
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO
#endif
                  DO 60 J = JTOP, MIN( KBOT, K+3 )
                     IF( DEBUG ) WRITE(*,*) 'Right mult on H: J=',J
                     REFSUM = V( 1, M )*( H( J, K+1 )+V( 2, M )*
     $                        H( J, K+2 )+V( 3, M )*H( J, K+3 ) )
                     H( J, K+1 ) = H( J, K+1 ) - REFSUM
                     H( J, K+2 ) = H( J, K+2 ) - REFSUM*V( 2, M )
                     H( J, K+3 ) = H( J, K+3 ) - REFSUM*V( 3, M )
   60             CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
*
                  IF( ACCUM ) THEN
*
*                    ==== Accumulate U. (If necessary, update Z later
*                    .    with with an efficient matrix-matrix
*                    .    multiply.) ====
*
                     KMS = K - INCOL
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO
#endif
                     DO 70 J = MAX( 1, KTOP-INCOL ), KDU
                        REFSUM = V( 1, M )*( U( J, KMS+1 )+V( 2, M )*
     $                           U( J, KMS+2 )+V( 3, M )*U( J, KMS+3 ) )
                        U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                        U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*V( 2, M )
                        U( J, KMS+3 ) = U( J, KMS+3 ) - REFSUM*V( 3, M )
   70                CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif                
                  ELSE IF( WANTZ ) THEN
*
*                    ==== U is not accumulated, so update Z
*                    .    now by multiplying by reflections
*                    .    from the right. ====
*
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO
#endif
                     DO 80 J = ILOZ, IHIZ
                        REFSUM = V( 1, M )*( Z( J, K+1 )+V( 2, M )*
     $                           Z( J, K+2 )+V( 3, M )*Z( J, K+3 ) )
                        Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                        Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*V( 2, M )
                        Z( J, K+3 ) = Z( J, K+3 ) - REFSUM*V( 3, M )
   80                CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
                  END IF
               END IF
   90       CONTINUE
             IF( PRINT )
     $           CALL DLAPRNT(N,N,H,1,1,LDH,'H6',6)
*
*           ==== Special case: 2-by-2 reflection (if needed) ====
*
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 .AND. ( V( 1, M22 ).NE.ZERO ) ) THEN
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO
#endif
               DO 100 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = V( 1, M22 )*( H( J, K+1 )+V( 2, M22 )*
     $                     H( J, K+2 ) )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*V( 2, M22 )
  100          CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
*
               IF( ACCUM ) THEN
                  KMS = K - INCOL
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO 
#endif
                  DO 110 J = MAX( 1, KTOP-INCOL ), KDU
                     REFSUM = V( 1, M22 )*( U( J, KMS+1 )+V( 2, M22 )*
     $                        U( J, KMS+2 ) )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*V( 2, M22 )
  110             CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
               ELSE IF( WANTZ ) THEN
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( J, REFSUM )
*$OMP DO
#endif
                  DO 120 J = ILOZ, IHIZ
                     REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )*
     $                        Z( J, K+2 ) )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*V( 2, M22 )
  120             CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
               END IF
            END IF
*
*           ==== Vigilant deflation check ====
*
            MSTART = MTOP
            IF( KRCOL+3*( MSTART-1 ).LT.KTOP )
     $         MSTART = MSTART + 1
            MEND = MBOT
            IF( BMP22 )
     $         MEND = MEND + 1
            IF( KRCOL.EQ.KBOT-2 )
     $         MEND = MEND + 1
            DO 130 M = MSTART, MEND
               K = MIN( KBOT-1, KRCOL+3*( M-1 ) )
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( H( K+1, K ).NE.ZERO ) THEN
                  TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) )
                  IF( TST1.EQ.ZERO ) THEN
                     IF( K.GE.KTOP+1 )
     $                  TST1 = TST1 + ABS( H( K, K-1 ) )
                     IF( K.GE.KTOP+2 )
     $                  TST1 = TST1 + ABS( H( K, K-2 ) )
                     IF( K.GE.KTOP+3 )
     $                  TST1 = TST1 + ABS( H( K, K-3 ) )
                     IF( K.LE.KBOT-2 )
     $                  TST1 = TST1 + ABS( H( K+2, K+1 ) )
                     IF( K.LE.KBOT-3 )
     $                  TST1 = TST1 + ABS( H( K+3, K+1 ) )
                     IF( K.LE.KBOT-4 )
     $                  TST1 = TST1 + ABS( H( K+4, K+1 ) )
                  END IF
                  IF( ABS( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) )
     $                 THEN
                     H12 = MAX( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
                     H21 = MIN( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
                     H11 = MAX( ABS( H( K+1, K+1 ) ),
     $                     ABS( H( K, K )-H( K+1, K+1 ) ) )
                     H22 = MIN( ABS( H( K+1, K+1 ) ),
     $                     ABS( H( K, K )-H( K+1, K+1 ) ) )
                     SCL = H11 + H12
                     TST2 = H22*( H11 / SCL )
*
                     IF( TST2.EQ.ZERO .OR. H21*( H12 / SCL ).LE.
     $                   MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                  END IF
               END IF
  130       CONTINUE
*
*           ==== Fill in the last row of each bulge. ====
*
            MEND = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 3 )
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( M, K, REFSUM )
*$OMP DO
#endif
            DO 140 M = MTOP, MEND
               K = KRCOL + 3*( M-1 )
               REFSUM = V( 1, M )*V( 3, M )*H( K+4, K+3 )
               H( K+4, K+1 ) = -REFSUM
               H( K+4, K+2 ) = -REFSUM*V( 2, M )
               H( K+4, K+3 ) = H( K+4, K+3 ) - REFSUM*V( 3, M )
  140       CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
*
*           ==== End of near-the-diagonal bulge chase. ====
*
  150    CONTINUE
*
*        ==== Use U (if accumulated) to update far-from-diagonal
*        .    entries in H.  If required, use U to update Z as
*        .    well. ====
*
         IF( ACCUM ) THEN
            IF( WANTT ) THEN
               JTOP = 1
               JBOT = N
            ELSE
               JTOP = KTOP
               JBOT = KBOT
            END IF
            IF( DEBUG ) WRITE(*,*) 'JTOP,JBOT=',JTOP,JBOT
            K1 = MAX( 1, KTOP-INCOL )
            NU = ( KDU-MAX( 0, MAX(INCOL+KDU,NDCOL)-KBOT ) ) - K1 + 1
            IF( DEBUG ) WRITE(*,*) 'K1,NU=',K1,NU
            IF( ( .NOT.BLK22 ) .OR. ( INCOL.LT.KTOP ) .OR.
     $          ( NDCOL.GT.KBOT ) .OR. ( NS.LE.2 ) .OR.
     $           NU.LT.KDU ) THEN
*
*              ==== Updates not exploiting the 2-by-2 block
*              .    structure of U.  K1 and NU keep track of
*              .    the location and size of U in the special
*              .    cases of introducing bulges and chasing
*              .    bulges off the bottom.  In these special
*              .    cases and in case the number of shifts
*              .    is NS = 2, there is no 2-by-2 block
*              .    structure to exploit.  ====           
*
*              ==== Horizontal Multiply ====
*
               IF( PRINT )
     $              CALL DLAPRNT(N,N,H,1,1,LDH,'H7',6)
               DO 160 JCOL = MIN(MAX(INCOL+KDU,NDCOL),KBOT)+ 1, JBOT, NH
                  IF( DEBUG ) WRITE(*,*) 'JCOL=',JCOL
                  JLEN = MIN( NH, JBOT-JCOL+1 )
                  IF( DEBUG ) WRITE(*,*) 'DGEMM 1:',NU,JLEN,NU
#ifdef USE_OMP
                  IF( THREADS.EQ.1 ) THEN
                     GCHUNK = JLEN
                  ELSE
                     GCHUNK = MAX( 1, MIN( MAXCHUNK, JLEN / THREADS  ) )
                  END IF
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( JCOL2, JLEN2 )
*$OMP DO
                  DO 165 JCOL2 = JCOL, JCOL+JLEN-1, GCHUNK
                     JLEN2 = MIN( GCHUNK, JCOL+JLEN-1-JCOL2+1 )
                     CALL DGEMM( 'C', 'N', NU, JLEN2, NU, ONE, 
     $                    U( K1, K1 ), LDU, H( INCOL+K1, JCOL2 ), LDH, 
     $                    ZERO, WH(1, JCOL2-JCOL+1 ), LDWH )
                     CALL DLACPY( 'ALL', NU, JLEN2, WH(1, JCOL2-JCOL+1), 
     $                    LDWH, H( INCOL+K1, JCOL2 ), LDH )
 165              CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                  CALL DGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ),
     $                        LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH,
     $                        LDWH )
                  CALL DLACPY( 'ALL', NU, JLEN, WH, LDWH,
     $                         H( INCOL+K1, JCOL ), LDH )
#endif
  160          CONTINUE
               IF( PRINT )
     $              CALL DLAPRNT(N,N,H,1,1,LDH,'H8',6)
*
*              ==== Vertical multiply ====
*
               DO 170 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
                  IF( DEBUG ) WRITE(*,*) 'JROW=',JROW
                  JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
                  IF( DEBUG ) WRITE(*,*) 'DGEMM 2:',JLEN,NU,NU
#ifdef USE_OMP
                   IF( THREADS.EQ.1 ) THEN
                     GCHUNK = JLEN
                  ELSE
                     GCHUNK = MAX( 1, MIN( MAXCHUNK, JLEN / THREADS  ) )
                  END IF
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( JROW2, JLEN2 )
*$OMP DO
                  DO 175 JROW2 = JROW, JROW+JLEN-1, GCHUNK
                     JLEN2 = MIN( GCHUNK, JROW+JLEN-1-JROW2+1 )
                     CALL DGEMM( 'N', 'N', JLEN2, NU, NU, ONE,
     $                           H( JROW2, INCOL+K1 ), LDH, U( K1, K1 ),
     $                           LDU, ZERO, WV( JROW2-JROW+1, 1), LDWV )
                     CALL DLACPY( 'ALL', JLEN2, NU, WV( JROW2-JROW+1,1), 
     $                            LDWV, H( JROW2, INCOL+K1 ), LDH )
 175              CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                  CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE,
     $                        H( JROW, INCOL+K1 ), LDH, U( K1, K1 ),
     $                        LDU, ZERO, WV, LDWV )
                  CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV,
     $                         H( JROW, INCOL+K1 ), LDH )
#endif
  170          CONTINUE
               IF( PRINT )
     $              CALL DLAPRNT(N,N,H,1,1,LDH,'H9',6)
*
*              ==== Z multiply (also vertical) ====
*
               IF( WANTZ ) THEN
                  DO 180 JROW = ILOZ, IHIZ, NV
                     IF( DEBUG ) WRITE(*,*) 'JROW=',JROW
                     JLEN = MIN( NV, IHIZ-JROW+1 )
                     IF( DEBUG ) WRITE(*,*) 'DGEMM 3:',JLEN,NU,NU
#ifdef USE_OMP
                     IF( THREADS.EQ.1 ) THEN
                        GCHUNK = JLEN
                     ELSE
                        GCHUNK = MAX( 1, MIN( MAXCHUNK, JLEN / THREADS))
                     END IF
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( JROW2, JLEN2 )
*$OMP DO
                     DO 185 JROW2 = JROW, JROW+JLEN-1, GCHUNK
                        JLEN2 = MIN( GCHUNK, JROW+JLEN-1-JROW2+1 )
                        CALL DGEMM( 'N', 'N', JLEN2, NU, NU, ONE,
     $                              Z( JROW2, INCOL+K1 ), LDZ, 
     $                              U( K1, K1 ), LDU, ZERO, 
     $                              WV( JROW2-JROW+1, 1), LDWV )
                        CALL DLACPY( 'ALL', JLEN2, NU, 
     $                               WV( JROW2-JROW+1, 1), LDWV,
     $                               Z( JROW2, INCOL+K1 ), LDZ )
 185                 CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                     CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE,
     $                           Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ),
     $                           LDU, ZERO, WV, LDWV )
                     CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV,
     $                            Z( JROW, INCOL+K1 ), LDZ )
#endif
  180             CONTINUE
               END IF
            ELSE
*
*              ==== Updates exploiting U's 2-by-2 block structure.
*              .    (I2, I4, J2, J4 are the last rows and columns
*              .    of the blocks.) ====
*
               I2 = ( KDU+1 ) / 2
               I4 = KDU
               J2 = I4 - I2
               J4 = KDU
*
*              ==== KZS and KNZ deal with the band of zeros
*              .    along the diagonal of one of the triangular
*              .    blocks. ====
*
               KZS = ( J4-J2 ) - ( NS+1 )
               KNZ = NS + 1
*
*              ==== Horizontal multiply ====
*
               DO 190 JCOL = MIN(MAX(INCOL+KDU,NDCOL),KBOT)+ 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
                  IF( DEBUG ) WRITE(*,*) 'JCOL,JLEN=',JCOL,JLEN
*
*                 ==== Copy bottom of H to top+KZS of scratch ====
*                  (The first KZS rows get multiplied by zero.) ====
*
                  IF( DEBUG ) WRITE(*,*) 'INCOL,J2=',INCOL,J2
                  CALL DLACPY( 'ALL', KNZ, JLEN, H( INCOL+1+J2, JCOL ),
     $                 LDH, WH( KZS+1, 1 ), LDWH )
#ifdef USE_OMP
                  IF( THREADS.EQ.1 ) THEN
                     GCHUNK = JLEN
                  ELSE
                     GCHUNK = MAX( 1, MIN( MAXCHUNK, JLEN / THREADS))
                  END IF
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( JCOL2, JLEN2 )
*$OMP DO
                  DO 195 JCOL2 = JCOL, JCOL+JLEN-1, GCHUNK
                     JLEN2 = MIN( GCHUNK, JCOL+JLEN-1-JCOL2+1 )
                     CALL DLASET( 'ALL', KZS, JLEN2, ZERO, ZERO, 
     $                    WH( 1, JCOL2-JCOL+1 ), LDWH )
                     CALL DTRMM( 'L', 'U', 'C', 'N', KNZ, JLEN2, ONE,
     $                    U( J2+1, 1+KZS ), LDU, 
     $                    WH( KZS+1, JCOL2-JCOL+1 ), LDWH )
                     CALL DGEMM( 'C', 'N', I2, JLEN2, J2, ONE, U, 
     $                    LDU, H( INCOL+1, JCOL2 ), LDH, ONE, 
     $                    WH(1, JCOL2-JCOL+1 ), LDWH )
                     CALL DLACPY( 'ALL', J2, JLEN2, 
     $                    H( INCOL+1, JCOL2 ), LDH, 
     $                    WH( I2+1, JCOL2-JCOL+1  ), LDWH )
                     CALL DTRMM( 'L', 'L', 'C', 'N', J2, JLEN2, ONE,
     $                    U( 1, I2+1 ), LDU, 
     $                    WH( I2+1, JCOL2-JCOL+1 ), LDWH )
                     CALL DGEMM( 'C', 'N', I4-I2, JLEN2, J4-J2, ONE,
     $                    U( J2+1, I2+1 ), LDU,
     $                    H( INCOL+1+J2, JCOL2 ), LDH, ONE,
     $                    WH( I2+1, JCOL2-JCOL+1 ), LDWH )
                     CALL DLACPY( 'ALL', KDU, JLEN2, 
     $                    WH( 1, JCOL2-JCOL+1 ), LDWH,
     $                    H( INCOL+1, JCOL2 ), LDH )
 195              CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                  CALL DLASET( 'ALL', KZS, JLEN, ZERO, ZERO, WH, LDWH )
*
*                 ==== Multiply by U21' ====
*
                  CALL DTRMM( 'L', 'U', 'C', 'N', KNZ, JLEN, ONE,
     $                        U( J2+1, 1+KZS ), LDU, WH( KZS+1, 1 ),
     $                        LDWH )
*
*                 ==== Multiply top of H by U11' ====
*
                  IF( DEBUG ) WRITE(*,*) 'DGEMM 4:',I2,JLEN,J2
                  CALL DGEMM( 'C', 'N', I2, JLEN, J2, ONE, U, LDU,
     $                        H( INCOL+1, JCOL ), LDH, ONE, WH, LDWH )
*
*                 ==== Copy top of H bottom of WH ====
*
                  CALL DLACPY( 'ALL', J2, JLEN, H( INCOL+1, JCOL ), LDH,
     $                         WH( I2+1, 1 ), LDWH )
*
*                 ==== Multiply by U21' ====
*
                  IF( DEBUG ) WRITE(*,*) 'DTRMM 2:',J2,JLEN
                  CALL DTRMM( 'L', 'L', 'C', 'N', J2, JLEN, ONE,
     $                        U( 1, I2+1 ), LDU, WH( I2+1, 1 ), LDWH )
*
*                 ==== Multiply by U22 ====
*
                  IF( DEBUG ) WRITE(*,*) 'DGEMM 5:',I4-I2,JLEN,J4-J2
                  CALL DGEMM( 'C', 'N', I4-I2, JLEN, J4-J2, ONE,
     $                        U( J2+1, I2+1 ), LDU,
     $                        H( INCOL+1+J2, JCOL ), LDH, ONE,
     $                        WH( I2+1, 1 ), LDWH )
*
*                 ==== Copy it back ====
*
                  CALL DLACPY( 'ALL', KDU, JLEN, WH, LDWH,
     $                         H( INCOL+1, JCOL ), LDH )
#endif
  190          CONTINUE
*
*              ==== Vertical multiply ====
*
               DO 200 JROW = JTOP, MAX( INCOL, KTOP ) - 1, NV
                  JLEN = MIN( NV, MAX( INCOL, KTOP )-JROW )
*
*                 ==== Copy right of H to scratch (the first KZS
*                 .    columns get multiplied by zero) ====
*
                  CALL DLACPY( 'ALL', JLEN, KNZ, H( JROW, INCOL+1+J2 ),
     $                         LDH, WV( 1, 1+KZS ), LDWV )
#ifdef USE_OMP
                   IF( THREADS.EQ.1 ) THEN
                     GCHUNK = JLEN
                  ELSE
                     GCHUNK = MAX( 1, MIN( MAXCHUNK, JLEN / THREADS  ) )
                  END IF
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( JROW2, JLEN2 )
*$OMP DO
                  DO 205 JROW2 = JROW, JROW+JLEN-1, GCHUNK
                     JLEN2 = MIN( GCHUNK, JROW+JLEN-1-JROW2+1 )
                     CALL DLASET( 'ALL', JLEN2, KZS, ZERO, ZERO, 
     $                    WV(JROW2-JROW+1, 1 ), LDWV )
                     CALL DTRMM( 'R', 'U', 'N', 'N', JLEN2, KNZ, ONE,
     $                    U( J2+1, 1+KZS ), LDU, 
     $                    WV( JROW2-JROW+1, 1+KZS ), LDWV )
                     CALL DGEMM( 'N', 'N', JLEN2, I2, J2, ONE,
     $                    H( JROW2, INCOL+1 ), LDH, U, LDU, ONE, 
     $                    WV( JROW2-JROW+1,1), LDWV )
                     CALL DLACPY( 'ALL', JLEN2, J2, H( JROW2, INCOL+1 ), 
     $                    LDH, WV( JROW2-JROW+1, 1+I2 ), LDWV )
                     CALL DTRMM( 'R', 'L', 'N', 'N', JLEN2, I4-I2, ONE,
     $                    U( 1, I2+1 ), LDU, WV( JROW2-JROW+1, 1+I2 ), 
     $                    LDWV )
                     CALL DGEMM( 'N', 'N', JLEN2, I4-I2, J4-J2, ONE,
     $                    H( JROW2, INCOL+1+J2 ), LDH,
     $                    U( J2+1, I2+1 ), LDU, ONE, 
     $                    WV( JROW2-JROW+1, 1+I2 ), LDWV )
                     CALL DLACPY( 'ALL', JLEN2, KDU, 
     $                    WV( JROW2-JROW+1, 1 ), LDWV,
     $                    H( JROW2, INCOL+1 ), LDH )
 205              CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                  CALL DLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV, LDWV )
*
*                 ==== Multiply by U21 ====
*
                  CALL DTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,
     $                        U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),
     $                        LDWV )
*
*                 ==== Multiply by U11 ====
*
                  IF( DEBUG ) WRITE(*,*) 'DGEMM 6:',JLEN,I2,J2
                  CALL DGEMM( 'N', 'N', JLEN, I2, J2, ONE,
     $                        H( JROW, INCOL+1 ), LDH, U, LDU, ONE, WV,
     $                        LDWV )
*
*                 ==== Copy left of H to right of scratch ====
*
                  CALL DLACPY( 'ALL', JLEN, J2, H( JROW, INCOL+1 ), LDH,
     $                         WV( 1, 1+I2 ), LDWV )
*
*                 ==== Multiply by U21 ====
*
                  IF( DEBUG ) WRITE(*,*) 'DTRMM 4:',JLEN,I4-I2
                  CALL DTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,
     $                        U( 1, I2+1 ), LDU, WV( 1, 1+I2 ), LDWV )
*
*                 ==== Multiply by U22 ====
*
                  IF( DEBUG ) WRITE(*,*) 'DGEMM 7:',JLEN,I4-I2,I4-J2
                  CALL DGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,
     $                        H( JROW, INCOL+1+J2 ), LDH,
     $                        U( J2+1, I2+1 ), LDU, ONE, WV( 1, 1+I2 ),
     $                        LDWV )
*
*                 ==== Copy it back ====
*
                  CALL DLACPY( 'ALL', JLEN, KDU, WV, LDWV,
     $                         H( JROW, INCOL+1 ), LDH )
#endif
  200          CONTINUE
*
*              ==== Multiply Z (also vertical) ====
*
               IF( WANTZ ) THEN
                  DO 210 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
*
*                    ==== Copy right of Z to left of scratch (first
*                    .     KZS columns get multiplied by zero) ====
*
                     CALL DLACPY( 'ALL', JLEN, KNZ,
     $                            Z( JROW, INCOL+1+J2 ), LDZ,
     $                            WV( 1, 1+KZS ), LDWV )
#ifdef USE_OMP
                     IF( THREADS.EQ.1 ) THEN
                        GCHUNK = JLEN
                     ELSE
                        GCHUNK = MAX( 1, MIN( MAXCHUNK, JLEN / THREADS))
                     END IF
*$OMP PARALLEL DEFAULT( SHARED ), PRIVATE( JROW2, JLEN2 )
*$OMP DO
                     DO 215 JROW2 = JROW, JROW+JLEN-1, GCHUNK
                        JLEN2 = MIN( GCHUNK, JROW+JLEN-1-JROW2+1 )
                        CALL DLASET( 'ALL', JLEN2, KZS, ZERO, ZERO, 
     $                       WV( JROW2-JROW+1, 1 ), LDWV )
                        CALL DTRMM( 'R', 'U', 'N', 'N', JLEN2, KNZ, ONE,
     $                       U( J2+1, 1+KZS ), LDU, 
     $                       WV( JROW2-JROW+1, 1+KZS ), LDWV )
                        CALL DGEMM( 'N', 'N', JLEN2, I2, J2, ONE,
     $                       Z( JROW2, INCOL+1 ), LDZ, U, LDU, ONE,
     $                       WV( JROW2-JROW+1, 1 ), LDWV )
                        CALL DLACPY( 'ALL', JLEN2, J2, 
     $                       Z( JROW2, INCOL+1 ), LDZ, 
     $                       WV( JROW2-JROW+1, 1+I2 ), LDWV )
                        CALL DTRMM( 'R', 'L', 'N', 'N', JLEN2, I4-I2, 
     $                       ONE, U( 1, I2+1 ), LDU, 
     $                       WV( JROW2-JROW+1, 1+I2 ), LDWV )
                        CALL DGEMM( 'N', 'N', JLEN2, I4-I2, J4-J2, ONE,
     $                       Z( JROW2, INCOL+1+J2 ), LDZ,
     $                       U( J2+1, I2+1 ), LDU, ONE,
     $                       WV( JROW2-JROW+1, 1+I2 ), LDWV )
                        CALL DLACPY( 'ALL', JLEN2, KDU, 
     $                       WV( JROW2-JROW+1, 1), LDWV, 
     $                       Z( JROW2, INCOL+1 ), LDZ )
 215                 CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
*
*                    ==== Multiply by U12 ====
*
                     CALL DLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV,
     $                            LDWV )
                     IF( DEBUG ) WRITE(*,*) 'DTRMM 5:',JLEN,KNZ
                     CALL DTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,
     $                           U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),
     $                           LDWV )
*
*                    ==== Multiply by U11 ====
*
                     IF( DEBUG ) WRITE(*,*) 'DGEMM 8:',JLEN,I2,J2
                     CALL DGEMM( 'N', 'N', JLEN, I2, J2, ONE,
     $                           Z( JROW, INCOL+1 ), LDZ, U, LDU, ONE,
     $                           WV, LDWV )
*
*                    ==== Copy left of Z to right of scratch ====
*
                     CALL DLACPY( 'ALL', JLEN, J2, Z( JROW, INCOL+1 ),
     $                            LDZ, WV( 1, 1+I2 ), LDWV )
*
*                    ==== Multiply by U21 ====
*
                     IF( DEBUG ) WRITE(*,*) 'DTRMM 6:',JLEN,I4-I2
                     CALL DTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,
     $                           U( 1, I2+1 ), LDU, WV( 1, 1+I2 ),
     $                           LDWV )
*
*                    ==== Multiply by U22 ====
*
                     IF( DEBUG ) WRITE(*,*) 'DGEMM 8:',JLEN,I4-I2,J4-J2
                     CALL DGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,
     $                           Z( JROW, INCOL+1+J2 ), LDZ,
     $                           U( J2+1, I2+1 ), LDU, ONE,
     $                           WV( 1, 1+I2 ), LDWV )
*
*                    ==== Copy the result back to Z ====
*
                     CALL DLACPY( 'ALL', JLEN, KDU, WV, LDWV,
     $                            Z( JROW, INCOL+1 ), LDZ )
#endif
  210             CONTINUE
               END IF
            END IF
         END IF
         IF( PRINT ) THEN
            CALL DLAPRNT(N,N,H,1,1,LDH,'H',6)
            CALL DLAPRNT(N,N,Z,1,1,LDZ,'Z',6)
         END IF
  220 CONTINUE
*
*     ==== Clear out workspaces and return
*
      IF( N.GE.5 ) 
     $     CALL DLASET( 'Lower', N-4, N-4, ZERO, ZERO, H(5,1), LDH )
*
      IF( DEBUG ) WRITE(*,*) 'Leaving dlaqr6:',JOB,N,KTOP,KBOT
*
*     ==== End of DLAQR6 ====
*
      END
