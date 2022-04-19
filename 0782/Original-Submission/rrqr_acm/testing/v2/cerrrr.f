      SUBROUTINE CERRRR( PATH, NUNIT )
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
*  CERRRR tests the error exits for CGEQPX and CGEQPY.
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
      REAL               SVLUES( 4 ), RW( 2*NMAX )
      COMPLEX            A( NMAX, NMAX ), C( NMAX, NMAX ),
     $                   W( 2*NMAX+3*NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, CGEQPX
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

      A( 1, 1 ) = ( 1.0E+0, 0.0E+0 )
      A( 1, 2 ) = ( 2.0E+0, 0.0E+0 )
      A( 2, 2 ) = ( 3.0E+0, 0.0E+0 )
      A( 2, 1 ) = ( 4.0E+0, 0.0E+0 )






      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'RR' ) ) THEN
*
*        Test error exits for Rank-Revealing QR factorization
*
*        CGEQPX
*
         SRNAMT = 'CGEQPX'
         INFOT = 1
         CALL CGEQPX( 4, 0, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEQPX( 1, -1, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGEQPX( 1, 0, -1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEQPX( 1, 0, 0, -1, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGEQPX( 1, 2, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGEQPX( 2, 2, 0, 0, A, 2, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGEQPX( 1, 0, 0, 0, A, 1, C, 1, IP, -ONE, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGEQPX( 1, 0, 3, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, 2, RW, INFO )
         CALL CHKXER( 'CGEQPX', INFOT, NOUT, LERR, OK )
*
*        CGEQPY
*
         SRNAMT = 'CGEQPY'
         INFOT = 1
         CALL CGEQPY( 4, 0, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEQPY( 1, -1, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGEQPY( 1, 0, -1, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEQPY( 1, 0, 0, -1, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGEQPY( 1, 2, 0, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGEQPY( 2, 2, 0, 0, A, 2, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGEQPY( 1, 0, 0, 0, A, 1, C, 1, IP, -ONE, ORCOND,
     $                RANK, SVLUES, W, LW, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGEQPY( 1, 0, 3, 0, A, 1, C, 1, IP, IRCOND, ORCOND,
     $                RANK, SVLUES, W, 2, RW, INFO )
         CALL CHKXER( 'CGEQPY', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRRR
*
      END
