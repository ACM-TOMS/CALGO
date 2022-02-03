      SUBROUTINE PDLAHQR( WANTT, WANTZ, N, ILO, IHI, A, DESCA, WR, WI,
     $                    ILOZ, IHIZ, Z, DESCZ, WORK, LWORK, IWORK,
     $                    ILWORK, INFO )
*
*  -- ScaLAPACK driver routine (version 1.7.x) --
*     University of Umea and HPC2N, Umea, Sweden
*     February, 2008
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, ILWORK, INFO, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), WI( * ), WORK( * ), WR( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  This is a wrapper for codes calling PDLAHQR which redirects the call
*  to PDHSEQR and then returns. (By setting WRAPPER below to .FALSE.,
*  this subroutine instead calls PDLAQR1, our modified version of 
*  PDLAHQR.  
*
*  Restrictions: it is assumed that ILOZ.EQ.ILO .AND. IHIZ.EQ.IHI .
*
*  =====================================================================
*
*     .. Parameters ..
      LOGICAL            WRAPPER, DEBUG
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, 
     $                     ZERO = 0.0D+00, DEBUG = .FALSE. )
#ifdef USE_NEWPQR
      PARAMETER          ( WRAPPER = .TRUE. )
#else
      PARAMETER          ( WRAPPER = .FALSE. )
#endif
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I
      CHARACTER          JOB, COMPZ
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   TIMINGS( 10 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDHSEQR, PXERBLA
*     ..
*     .. Executable Statements ..
*
      LQUERY = LWORK.EQ.-1 .OR. ILWORK.EQ.-1
      IF( .NOT. LQUERY ) THEN
         IF( ILOZ.NE.ILO ) THEN
            INFO = -10
         END IF
         IF( IHIZ.NE.IHI ) THEN
            INFO = -11
         END IF
         IF( INFO.LT.0 ) THEN
            CALL PXERBLA( DESCA(CTXT_), 'PDLAHQR', -INFO )
            RETURN
         END IF
      END IF
*     
      IF( WANTT ) THEN
         JOB = 'S'
      ELSE
         JOB = 'E'
      END IF
      IF( WANTZ ) THEN
         COMPZ = 'V'
      ELSE
         COMPZ = 'N'
      END IF
*     
      IF( WRAPPER ) THEN
         DO 10 I = 1, 10
            TIMINGS( I ) = ZERO
 10      CONTINUE
         CALL PDHSEQR( JOB, COMPZ, N, ILO, IHI, A, DESCA, WR, WI, Z,
     $                 DESCZ, WORK, LWORK, IWORK, ILWORK, INFO,
     $                 TIMINGS )
         IF( .NOT. LQUERY .AND. INFO.LT.0 ) THEN
            CALL PXERBLA( DESCA(CTXT_), 'PDHSEQR', -INFO )
            RETURN
         ELSEIF( .NOT. LQUERY .AND. INFO.GT.0 ) THEN
            CALL PXERBLA( DESCA(CTXT_), 'PDHSEQR', -INFO )
            RETURN
         END IF
      ELSE
         CALL PDLAQR1( WANTT, WANTZ, N, ILO, IHI, A, DESCA, WR, WI,
     $                 ILOZ, IHIZ, Z, DESCZ, WORK, LWORK, 
     $                 IWORK, ILWORK, INFO )
         IF( .NOT. LQUERY .AND. INFO.LT.0 ) THEN
            CALL PXERBLA( DESCA(CTXT_), 'PDLAQR1', -INFO )
            RETURN
         ELSEIF( .NOT. LQUERY .AND. INFO.GT.0 ) THEN
            CALL PXERBLA( DESCA(CTXT_), 'PDLAQR1', -INFO )
            RETURN
         END IF
      END IF
*
      RETURN
*
*     END OF PDLAHQR
*
      END
