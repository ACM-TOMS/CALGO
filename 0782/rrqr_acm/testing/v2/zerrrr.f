      SUBROUTINE ZERRRR( PATH, NUNIT )
*
*  -- LAPACK test routine (version 1.0b) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     Rewritten for new least-squares solvers.
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  ZERRRR tests the error exits for ZGEQPX and ZGEQPY.
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
      DOUBLE PRECISION   ZERO, ONE



      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            INFO, RANK, LW
      DOUBLE PRECISION   IRCOND, ORCOND
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      DOUBLE PRECISION   SVLUES( 4 ), RW( 2*NMAX )
      COMPLEX*16         A( NMAX, NMAX ), C( NMAX, NMAX ),
     $                   W( 2*NMAX+3*NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, ZGEQPX
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






      A( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
      A( 1, 2 ) = ( 2.0D+0, 0.0D+0 )
      A( 2, 2 ) = ( 3.0D+0, 0.0D+0 )
      A( 2, 1 ) = ( 4.0D+0, 0.0D+0 )

      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'RR' ) ) THEN
*
*        Test error exits for Rank-Revealing QR factorization
*
*        ZGEQPX
*
         SRNAMT = 'ZGEQPX'
         INFOT = 1
         CALL ZGEQPX( 4, 0, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEQPX( 1, -1, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEQPX( 1, 0, -1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGEQPX( 1, 0, 0, -1, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGEQPX( 1, 2, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGEQPX( 2, 2, 0, 0, A, 2, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGEQPX( 1, 0, 0, 0, A, 1, C, 1, IP, -ONE, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL ZGEQPX( 1, 0, 3, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, 2, RW, INFO )
         CALL CHKXER( 'ZGEQPX', INFOT, NOUT, LERR, OK )
*
*        ZGEQPY
*
         SRNAMT = 'ZGEQPY'
         INFOT = 1
         CALL ZGEQPY( 4, 0, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEQPY( 1, -1, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEQPY( 1, 0, -1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGEQPY( 1, 0, 0, -1, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGEQPY( 1, 2, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGEQPY( 2, 2, 0, 0, A, 2, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGEQPY( 1, 0, 0, 0, A, 1, C, 1, IP, -ONE, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL ZGEQPY( 1, 0, 3, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, 2, RW, INFO )
         CALL CHKXER( 'ZGEQPY', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of ZERRRR
*
      END
