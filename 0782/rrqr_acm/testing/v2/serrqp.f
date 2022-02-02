      SUBROUTINE SERRQP( PATH, NUNIT )
*
*  -- LAPACK test routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  SERRQP tests the error exits for SGEQPF.
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
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            INFO
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      REAL               A( NMAX, NMAX ), TAU( NMAX ), W( 3*NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SGEQPF
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
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = 1.
      A( 1, 2 ) = 2.
      A( 2, 2 ) = 3.
      A( 2, 1 ) = 4.
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'QP' ) ) THEN
*
*        Test error exits for QR factorization with pivoting
*
*        SGEQPF
*
         SRNAMT = 'SGEQPF'
         INFOT = 1
         CALL SGEQPF( -1, 0, A, 1, IP, TAU, W, INFO )
         CALL CHKXER( 'SGEQPF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SGEQPF( 0, -1, A, 1, IP, TAU, W, INFO )
         CALL CHKXER( 'SGEQPF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SGEQPF( 2, 0, A, 1, IP, TAU, W, INFO )
         CALL CHKXER( 'SGEQPF', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRQP
*
      END
