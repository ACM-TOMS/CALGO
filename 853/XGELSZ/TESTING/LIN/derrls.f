      SUBROUTINE DERRLS( PATH, NUNIT )
*
*     This code is part of a package for solving rank deficient least
*     squares problems, written by:
*     ==================================================================
*     L. Foster                   and   R. Kommu
*     Department of Mathematics         Department of Physics
*     San Jose State University         San Jose State University
*     San Jose, CA 95192                San Jose, CA 95192
*     foster@math.sjsu.edu              rkommu@email.sjsu.edu
*     ==================================================================
*     03/05/2004
*
*     This code is a modification of the corresponding
*     LAPACK (Version 3.0) routine
*
*  -- LAPACK test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  DERRLS tests the error exits for the DOUBLE PRECISION least squares
*  driver routines (DGELS, DGELSS, DGELSX, DGELSY, DGELSZ, DGELSD).
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
      INTEGER            INFO, IRNK
      DOUBLE PRECISION   RCOND
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ),
     $                   W( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DGELS, DGELSD, DGELSS, DGELSX,
     $                   DGELSY, DGELSZ
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
      A( 1, 1 ) = 1.0D+0
      A( 1, 2 ) = 2.0D+0
      A( 2, 2 ) = 3.0D+0
      A( 2, 1 ) = 4.0D+0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        Test error exits for the least squares driver routines.
*
*        DGELS
*
         SRNAMT = 'DGELS '
         INFOT = 1
         CALL DGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, LERR, OK )
*
*        DGELSS
*
         SRNAMT = 'DGELSS'
         INFOT = 1
         CALL DGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO )
         CALL CHKXER( 'DGELSS', INFOT, NOUT, LERR, OK )
*
*        DGELSX
*
         SRNAMT = 'DGELSX'
         INFOT = 1
         CALL DGELSX( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'DGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSX( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'DGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSX( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'DGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSX( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'DGELSX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSX( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, INFO )
         CALL CHKXER( 'DGELSX', INFOT, NOUT, LERR, OK )
*
*        DGELSY
*
         SRNAMT = 'DGELSY'
         INFOT = 1
         CALL DGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGELSY( 2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSY', INFOT, NOUT, LERR, OK )
*
*        DGELSZ
*
         SRNAMT = 'DGELSZ'
         INFOT = 1
         CALL DGELSZ( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSZ( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSZ( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSZ( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSZ( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10,
     $                INFO )
         CALL CHKXER( 'DGELSZ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGELSZ( 2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO )
         CALL CHKXER( 'DGELSZ', INFOT, NOUT, LERR, OK )
*
*        DGELSD
*
         SRNAMT = 'DGELSD'
         INFOT = 1
         CALL DGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP,
     $                INFO )
         CALL CHKXER( 'DGELSD', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRLS
*
      END
