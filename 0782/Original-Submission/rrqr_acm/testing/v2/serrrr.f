      SUBROUTINE SERRRR( PATH, NUNIT )
*
*  -- LAPACK test routine (version 1.0b) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     Rewritten for RRQR subroutines.
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  SERRRR tests the error exits for SGEQPX and SGEQPY.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name for the routines to be tested.
*
*  NUNIT   (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 2 )
      REAL               ZERO, ONE

      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )



*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            INFO, RANK, LW
      REAL               IRCOND, ORCOND
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ),
     $                   SVLUES( 4 ), W( 2*NMAX+3*NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SGEQPX, SGEQPY
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      IRCOND = ZERO
      RANK = 1
      LW = 2*NMAX+3*NMAX
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      A( 1, 1 ) = 1.0E+0
      A( 1, 2 ) = 2.0E+0
      A( 2, 2 ) = 3.0E+0
      A( 2, 1 ) = 4.0E+0






      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'RR' ) ) THEN
*
*        Test error exits for Rank-Revealing QR factorization
*
*        SGEQPX
*
         SRNAMT = 'SGEQPX'
         INFOT = 1
         CALL SGEQPX( 4, 0, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGEQPX( 1, -1, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGEQPX( 1, 0, -1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SGEQPX( 1, 0, 0, -1, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SGEQPX( 1, 2, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SGEQPX( 2, 2, 0, 0, A, 2, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SGEQPX( 1, 0, 0, 0, A, 1, C, 1, IP, -ONE, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL SGEQPX( 1, 0, 1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, 2, INFO )
         CALL CHKXER( 'SGEQPX', INFOT, NOUT, LERR, OK )
*
*        SGEQPY
*
         SRNAMT = 'SGEQPY'
         INFOT = 1
         CALL SGEQPY( 4, 0, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGEQPY( 1, -1, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SGEQPY( 1, 0, -1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SGEQPY( 1, 0, 0, -1, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SGEQPY( 1, 2, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SGEQPY( 2, 2, 0, 0, A, 2, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SGEQPY( 1, 0, 0, 0, A, 1, C, 1, IP, -ONE, ORCOND,
     $                RANK, SVLUES, W, LW, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL SGEQPY( 1, 0, 1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, 2, INFO )
         CALL CHKXER( 'SGEQPY', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRRR
*
      END
