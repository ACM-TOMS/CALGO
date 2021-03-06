      SUBROUTINE PDROT( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, 
     $                  INCY, CS, SN, WORK, LWORK, INFO )
C   
C  -- ScaLAPACK-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Ume?, Sweden
C     written by Robert Granat, May 15, 2007.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER            N, IX, JX, INCX, IY, JY, INCY, LWORK, INFO
      DOUBLE PRECISION   CS, SN
C     ..
C     .. Array Arguments ..
      INTEGER            DESCX( * ), DESCY( * )
      DOUBLE PRECISION   X( * ), Y( * ), WORK( * )
C     ..    
C
C  Purpose
C  =======
C  PDROT applies a planar rotation defined by CS and SN to the
C  two distributed vectors sub(X) and sub(Y).   
C
C  Notes
C  =====
C
C  Each global data object is described by an associated description
C  vector.  This vector stores the information required to establish
C  the mapping between an object element and its corresponding process
C  and memory location.
C
C  Let A be a generic term for any 2D block cyclicly distributed array.
C  Such a global array has an associated description vector DESCA.
C  In the following comments, the character _ should be read as
C  "of the global array".
C
C  NOTATION        STORED IN      EXPLANATION
C  --------------- -------------- --------------------------------------
C  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
C                                 DTYPE_A = 1.
C  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
C                                 the BLACS process grid A is distribu-
C                                 ted over. The context itself is glo-
C                                 bal, but the handle (the integer
C                                 value) may vary.
C  M_A    (global) DESCA( M_ )    The number of rows in the global
C                                 array A.
C  N_A    (global) DESCA( N_ )    The number of columns in the global
C                                 array A.
C  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
C                                 the rows of the array.
C  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
C                                 the columns of the array.
C  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
C                                 row of the array A is distributed.
C  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
C                                 first column of the array A is
C                                 distributed.
C  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
C                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
C
C  Let K be the number of rows or columns of a distributed matrix,
C  and assume that its process grid has dimension p x q.
C  LOCr( K ) denotes the number of elements of K that a process
C  would receive if K were distributed over the p processes of its
C  process column.
C  Similarly, LOCc( K ) denotes the number of elements of K that a
C  process would receive if K were distributed over the q processes of
C  its process row.
C  The values of LOCr() and LOCc() may be determined via a call to the
C  ScaLAPACK tool function, NUMROC:
C          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
C  An upper bound for these quantities may be computed by:
C          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C     
C  Arguments
C  =========
C     
C  N       (global input) INTEGER
C          The number of elements to operate on when applying the planar
C          rotation to X and Y. N>=0.
C 
C  X       (local input/local output) DOUBLE PRECSION array of dimension 
C          ( (JX-1)*M_X + IX + ( N - 1 )*abs( INCX ) )
C          This array contains the entries of the distributed vector 
C          sub( X ).
C 
C  IX      (global input) INTEGER 
C          The global row index of the submatrix of the distributed 
C          matrix X to operate on. If INCX = 1, then it is required
C          that IX = IY. 1 <= IX <= M_X.
C
C  JX      (global input) INTEGER 
C          The global column index of the submatrix of the distributed 
C          matrix X to operate on. If INCX = M_X, then it is required
C          that JX = JY. 1 <= IX <= N_X.
C 
C  DESCX   (global and local input) INTEGER array of dimension 9
C          The array descriptor of the distributed matrix X. 
C  
C  INCX    (global input) INTEGER 
C          The global increment for the elements of X. Only two values 
C          of INCX are supported in this version, namely 1 and M_X. 
C          Moreover, it must hold that INCX = M_X if INCY = M_Y and
C          that INCX = 1 if INCY = 1.
C
C  Y       (local input/local output) DOUBLE PRECSION array of dimension 
C          ( (JY-1)*M_Y + IY + ( N - 1 )*abs( INCY ) )
C          This array contains the entries of the distributed vector 
C          sub( Y ).
C 
C  IY      (global input) INTEGER 
C          The global row index of the submatrix of the distributed 
C          matrix Y to operate on. If INCY = 1, then it is required
C          that IY = IX. 1 <= IY <= M_Y.
C
C  JY      (global input) INTEGER 
C          The global column index of the submatrix of the distributed 
C          matrix Y to operate on. If INCY = M_X, then it is required
C          that JY = JX. 1 <= JY <= N_Y.
C
C  DESCY   (global and local input) INTEGER array of dimension 9
C          The array descriptor of the distributed matrix Y. 
C
C  INCY    (global input) INTEGER 
C          The global increment for the elements of Y. Only two values 
C          of INCY are supported in this version, namely 1 and M_Y. 
C          Moreover, it must hold that INCY = M_Y if INCX = M_X and
C          that INCY = 1 if INCX = 1.
C
C  CS,SN (global input) DOUBLE PRECISION
C          The parameters defining the properties of the planar 
C          rotation. It must hold that 0 <= CS,SN <= 1 and that
C          SN**2 + CS**2 = 1. The latter is hardly checked in
C          finite precision arithmetics.   
C
C  WORK    (local input) DOUBLE PRECISION array of dimension LWORK           
C          Local workspace area.
C
C  LWORK   (local input) INTEGER
C          The length of the workspace array WORK. 
C          If INCX = 1 and INCY = 1, then LWORK = 2*MB_X
C
C          If LWORK = -1, then a workspace query is assumed; the
C          routine only calculates the optimal size of the WORK array,
C          returns this value as the first entry of the IWORK array, and
C          no error message related to LIWORK is issued by PXERBLA.
C  
C  INFO    (global output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          If the i-th argument is an array and the j-entry had
C          an illegal value, then INFO = -(i*100+j), if the i-th
C          argument is a scalar and had an illegal value, then INFO = -i.
C 
C  Additional requirements
C  =======================
C
C  The following alignment requirements must hold:
C  (a) DESCX( MB_ ) = DESCY( MB_ ) and DESCX( NB_ ) = DESCY( NB_ )
C  (b) DESCX( RSRC_ ) = DESCY( RSRC_ )
C  (c) DESCX( CSRC_ ) = DESCY( CSRC_ )
C
C  =====================================================================
C     ..
C     ..
C     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 ) 
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, LEFT, RIGHT
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, NPROCS, 
     $                   MB, NB, XYROWS, XYCOLS, RSRC1, RSRC2, CSRC1,
     $                   CSRC2, ICOFFXY, IROFFXY, MNWRK, LLDX, LLDY,
     $                   INDX, JXX, XLOC1, XLOC2, RSRC, CSRC, YLOC1,
     $                   YLOC2, JYY, IXX, IYY
C     ..
C     .. External Functions ..
      INTEGER            NUMROC, INDXG2P, INDXG2L
      EXTERNAL           NUMROC, INDXG2P, INDXG2L
C     ..
C     .. External Subroutines ..
      EXTERNAL           DROT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Local Functions .. 
      INTEGER            ICEIL
C     .. 
C     .. Executable Statements ..
C
C     Get grid parameters 
C     
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
C
C     Test and decode inparameters
C
      LQUERY = LWORK.EQ.-1
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSEIF( IX.LT.1 .OR. IX.GT.DESCX(M_) ) THEN
         INFO = -3
      ELSEIF( JX.LT.1 .OR. JX.GT.DESCX(N_) ) THEN
         INFO = -4
      ELSEIF( INCX.NE.1 .AND. INCX.NE.DESCX(M_) ) THEN
         INFO = -6
       ELSEIF( IY.LT.1 .OR. IY.GT.DESCY(M_) ) THEN
         INFO = -8
      ELSEIF( JY.LT.1 .OR. JY.GT.DESCY(N_) ) THEN
         INFO = -9
      ELSEIF( INCY.NE.1 .AND. INCY.NE.DESCY(M_) ) THEN
         INFO = -11
      ELSEIF( (INCX.EQ.DESCX(M_) .AND. INCY.NE.DESCY(M_)) .OR.
     $        (INCX.EQ.1 .AND. INCY.NE.1 ) ) THEN
         INFO = -11
      ELSEIF( (INCX.EQ.1 .AND. INCY.EQ.1) .AND.
     $        IX.NE.IY ) THEN
         INFO = -8
      ELSEIF( (INCX.EQ.DESCX(M_) .AND. INCY.EQ.DESCY(M_)) .AND.
     $        JX.NE.JY ) THEN
         INFO = -9
      END IF
C
C     Compute the direction of the planar rotation
C
      LEFT  = INCX.EQ.DESCX(M_) .AND. INCY.EQ.DESCY(M_)
      RIGHT = INCX.EQ.1 .AND. INCY.EQ.1
C
C     Check blocking factors and root processor
C
      IF( INFO.EQ.0 ) THEN
         IF( LEFT .AND. DESCX(NB_).NE.DESCY(NB_) ) THEN
            INFO = -(100*5 + NB_)
         END IF
         IF( RIGHT .AND. DESCX(MB_).NE.DESCY(NB_) ) THEN
            INFO = -(100*10 + MB_)
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LEFT .AND. DESCX(CSRC_).NE.DESCY(CSRC_) ) THEN
            INFO = -(100*5 + CSRC_)
         END IF
         IF( RIGHT .AND. DESCX(RSRC_).NE.DESCY(RSRC_) ) THEN
            INFO = -(100*10 + RSRC_)
         END IF
      END IF
C
C     Compute workspace
C
      MB = DESCX( MB_ )
      NB = DESCX( NB_ )
      IF( LEFT ) THEN
         RSRC1 = INDXG2P( IX, MB, MYROW, DESCX(RSRC_), NPROW )
         RSRC2 = INDXG2P( IY, MB, MYROW, DESCY(RSRC_), NPROW )
         CSRC  = INDXG2P( JX, NB, MYCOL, DESCX(CSRC_), NPCOL ) 
         ICOFFXY = MOD( JX - 1, NB )
         XYCOLS = NUMROC( N+ICOFFXY, NB, MYCOL, CSRC, NPCOL )
         IF( ( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC2 ) .AND.
     $         MYCOL.EQ.CSRC ) XYCOLS = XYCOLS - ICOFFXY
         IF( RSRC1.NE.RSRC2 ) THEN
            MNWRK = XYCOLS 
         ELSE
            MNWRK = 0
         END IF
      ELSEIF( RIGHT ) THEN
         CSRC1 = INDXG2P( JX, NB, MYCOL, DESCX(CSRC_), NPCOL )
         CSRC2 = INDXG2P( JY, NB, MYCOL, DESCY(CSRC_), NPCOL )
         RSRC  = INDXG2P( IX, MB, MYROW, DESCX(RSRC_), NPROW ) 
         IROFFXY = MOD( IX - 1, MB )
         XYROWS = NUMROC( N+IROFFXY, MB, MYROW, RSRC, NPROW )
         IF( ( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC2  ) .AND.
     $         MYROW.EQ.RSRC ) XYROWS = XYROWS - IROFFXY
         IF( CSRC1.NE.CSRC2 ) THEN
            MNWRK = XYROWS
         ELSE
            MNWRK = 0
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LWORK.LT.MNWRK ) INFO = -15
      END IF
C     
C     Return if some argument is incorrect
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDROT', -INFO )
         RETURN
      ELSEIF( LQUERY ) THEN
         WORK( 1 ) = DBLE(MNWRK)
         RETURN
      END IF
C    
C     Quick return if possible.
C    
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
C
C     Extract local leading dimensions
C
      LLDX = DESCX( LLD_ )
      LLDY = DESCY( LLD_ )
C
C     If we have only one process, use the corresponding LAPACK
C     routine and return
C
      IF( NPROCS.EQ.1 ) THEN
         IF( LEFT ) THEN
            CALL DROT( N, X((JX-1)*LLDX+IX), LLDX, Y((JY-1)*LLDY+IY), 
     $           LLDY, CS, SN )
         ELSEIF( RIGHT ) THEN
            CALL DROT( N, X((JX-1)*LLDX+IX), 1, Y((JY-1)*LLDY+IY), 
     $           1, CS, SN )
         END IF
         RETURN
      END IF
C
C     Exchange data between processors if necessary and perform planar
C     rotation
C
      IF( LEFT ) THEN
         DO 10 INDX = 1, NPCOL
            IF( MYROW.EQ.RSRC1 ) THEN
               IF( INDX.EQ.1 ) THEN
                  JXX = JX
               ELSE
                  JXX = JX-ICOFFXY + (INDX-1)*NB
               END IF
               CALL INFOG2L( IX, JXX, DESCX, NPROW, NPCOL, MYROW, 
     $                       MYCOL, XLOC1, XLOC2, RSRC, CSRC )
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IF( RSRC1.NE.RSRC2 ) THEN
                     CALL DGESD2D( ICTXT, 1, XYCOLS, 
     $                             X((XLOC2-1)*LLDX+XLOC1), LLDX, 
     $                             MOD(RSRC+1,NPROW), CSRC )
                     CALL DGERV2D( ICTXT, 1, XYCOLS, WORK, 1,
     $                             MOD(RSRC+1,NPROW), CSRC )
                     CALL DROT( XYCOLS, X((XLOC2-1)*LLDX+XLOC1),
     $                          LLDX, WORK, 1, CS, SN )
                  ELSE
                     CALL INFOG2L( IY, JXX, DESCY, NPROW, NPCOL, 
     $                             MYROW, MYCOL, YLOC1, YLOC2, RSRC, 
     $                             CSRC )
                     CALL DROT( XYCOLS, X((XLOC2-1)*LLDX+XLOC1),
     $                          LLDX, Y((YLOC2-1)*LLDY+YLOC1), LLDY, CS, 
     $                          SN )
                  END IF
               END IF
            END IF
            IF( MYROW.EQ.RSRC2 .AND. RSRC1.NE.RSRC2 ) THEN
               IF( INDX.EQ.1 ) THEN
                  JYY = JY
               ELSE
                  JYY = JY-ICOFFXY + (INDX-1)*NB
               END IF
               CALL INFOG2L( IY, JYY, DESCY, NPROW, NPCOL, MYROW, 
     $                       MYCOL, YLOC1, YLOC2, RSRC, CSRC )
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  CALL DGESD2D( ICTXT, 1, XYCOLS, 
     $                          Y((YLOC2-1)*LLDY+YLOC1), LLDY, 
     $                          MOD(RSRC-1+NPROW,NPROW), CSRC )
                  CALL DGERV2D( ICTXT, 1, XYCOLS, WORK, 1,
     $                          MOD(RSRC-1+NPROW,NPROW), CSRC )
                  CALL DROT( XYCOLS, WORK, 1, Y((YLOC2-1)*LLDY+YLOC1),
     $                       LLDY, CS, SN )
               END IF
            END IF
 10      CONTINUE
      ELSEIF( RIGHT ) THEN
         DO 20 INDX = 1, NPROW
            IF( MYCOL.EQ.CSRC1 ) THEN
               IF( INDX.EQ.1 ) THEN
                  IXX = IX
               ELSE
                  IXX = IX-IROFFXY + (INDX-1)*MB
               END IF
               CALL INFOG2L( IXX, JX, DESCX, NPROW, NPCOL, MYROW, 
     $                       MYCOL, XLOC1, XLOC2, RSRC, CSRC )
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IF( CSRC1.NE.CSRC2 ) THEN
                     CALL DGESD2D( ICTXT, XYROWS, 1, 
     $                             X((XLOC2-1)*LLDX+XLOC1), LLDX, 
     $                             RSRC, MOD(CSRC+1,NPCOL) )
                     CALL DGERV2D( ICTXT, XYROWS, 1, WORK, XYROWS,
     $                             RSRC, MOD(CSRC+1,NPCOL) )
                     CALL DROT( XYROWS, X((XLOC2-1)*LLDX+XLOC1),
     $                          1, WORK, 1, CS, SN )
                  ELSE
                     CALL INFOG2L( IXX, JY, DESCY, NPROW, NPCOL, 
     $                             MYROW, MYCOL, YLOC1, YLOC2, RSRC, 
     $                             CSRC )
                     CALL DROT( XYROWS, X((XLOC2-1)*LLDX+XLOC1),
     $                          1, Y((YLOC2-1)*LLDY+YLOC1), 1, CS, 
     $                          SN )
                  END IF
               END IF
            END IF
            IF( MYCOL.EQ.CSRC2 .AND. CSRC1.NE.CSRC2 ) THEN
               IF( INDX.EQ.1 ) THEN
                  IYY = IY
               ELSE
                  IYY = IY-IROFFXY + (INDX-1)*MB
               END IF
               CALL INFOG2L( IYY, JY, DESCY, NPROW, NPCOL, MYROW, 
     $                       MYCOL, YLOC1, YLOC2, RSRC, CSRC )
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  CALL DGESD2D( ICTXT, XYROWS, 1, 
     $                          Y((YLOC2-1)*LLDY+YLOC1), LLDY, 
     $                          RSRC, MOD(CSRC-1+NPCOL,NPCOL) )
                  CALL DGERV2D( ICTXT, XYROWS, 1, WORK, XYROWS,
     $                          RSRC, MOD(CSRC-1+NPCOL,NPCOL) )
                  CALL DROT( XYROWS, WORK, 1, Y((YLOC2-1)*LLDY+YLOC1),
     $                       1, CS, SN )
               END IF
            END IF
 20      CONTINUE
      END IF      
C     
C     Store minimum workspace requirements in WORK-array and return 
C
      WORK( 1 ) = DBLE(MNWRK)
      RETURN
C
C     End of PDROT
C
      END
