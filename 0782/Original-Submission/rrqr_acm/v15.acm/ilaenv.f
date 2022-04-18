      INTEGER          FUNCTION ILAENV( ISPEC, SUBNAM, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     August 17, 1990
*
*     ** TEST VERSION **
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    OPTS, SUBNAM
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV returns machine and problem-dependent parameters.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies what quantity is to be returned (as the function's
*          value):
*          = 1: The optimum blocksize.  ILAENV(1,...)=1 is a flag that
*               no blocking should be done.
*          = 2: The minimum blocksize.
*          = 3: The "crossover point".  When two versions of a solution
*               method are implemented, one of which is faster for
*               problems and the other faster for smaller problems,
*               this is the largest problem for which the small-problem
*               method is preferred.
*          = 4: the number of shifts to use.  At present, this is only
*               appropriate for the nonsymmetric eigenvalue and
*               generalized eigenvalue routines, and -1 will be returned
*               if SUBNAM(2:6) is not 'HSEQR' or 'HGEQZ'.
*          = 5: The minimum second dimension for blocked updates.  This
*               is primarily intended for methods which update
*               rectangular blocks using Householder transformations
*               with short Householder vectors (e.g., xLAEBC and
*               xLAGBC).  If a  k x k  Householder transformation is
*               used to update a  k x m  block, then blocking (i.e.,
*               use of xGEMM) will only be done if k is at least
*               ILAENV(2,...) and m is at least ILAENV(5,...)
*          = 6: (Used only by the SVD drivers.)  When reducing an m x n
*               matrix to bidiagonal form, if the larger dimension is
*               less than ILAENV(6,...,m,n,,), then the usual procedure
*               is preferred.  If the larger dimension is larger than
*               ILAENV(6,...,m,n,,), then it is preferred to first use
*               a QR (or LQ) factorization to make it triangular, and
*               then reduce it to bidiagonal form.
*
*          = 7: Number of processors.
*
*  SUBNAM  (input) CHARACTER*(*)
*          The name of the calling routine, or the routine which is
*          expected to use the value returned.
*
*  OPTS    (input) CHARACTER*(*)
*          The values of the CHARACTER*1 options passed to the routine
*          whose name is in SUBNAM, all run together.  For example,
*          a subroutine called with "CALL SGEZZZ('T','Y',N4,'C',A,LDA)"
*          would pass 'SGEZZZ' as SUBNAM and 'TYC' as OPTS.
*
*  N1,N2,N3,N4 (input) INTEGER
*          The problem dimensions, in the order that they appear in
*          the calling sequence.  If there is only one dimension
*          (customarily called N), then that value should be passed
*          as N1 and N2, N3, and N4 will be ignored, etc.
*
* (ILAENV) (output) INTEGER
*          The function value returned will be the value specified
*          by ISPEC.  If no reasonable value is available, a negative
*          value will be returned, otherwise the returned value will
*          be non-negative.  Note that the value returned may be
*          unreasonably large for the problem, e.g., a blocksize of
*          32 for an 8 x 8 problem, and thus should be restricted
*          to whatever range is appropriate.
*
*  Further Details
*  ======= =======
*
*  The calling sequence is intended to match up to the arguments of
*  the routine needing the value in a simple and mindless way.
*  For example, if SGEZZZ were defined as
*
*       SUBROUTINE SGEZZZ( OPT1, OPT2, N4, OPT3, N1, A, LDA )
*       CHARACTER*1  OPT1, OPT2, OPT3
*       INTEGER      N1, LDA, N4
*       REAL         A(LDA,*)
*
*  then in SGEZZ, the blocksize would be found by a call like:
*
*       NBLOCK = ILAENV( 1, 'SGEZZZ', OPT1//OPT2//OPT3, N4, N1, 0, 0 )
*
*  It would be further checked and restricted by code like the
*  following:
*
*       NBLOCK = MAX( 1, MIN( N4, NBLOCK ) )
*       IF( NBLOCK.EQ.1 ) THEN
*  c
*  c    unblocked method
*  c
*    . . .
*
*=======================================================================
*
*     .. Local Scalars ..
      CHARACTER          NAME1
      CHARACTER*2        NAME23
      CHARACTER*3        NAME46
      CHARACTER*6        NAMWRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, REAL
*     ..
*     .. Arrays in Common ..
      INTEGER            NENVIR( 10 )
*     ..
*     .. Common blocks ..
      COMMON             / CENVIR / NENVIR
*     ..
*     .. Save statement ..
      SAVE               / CENVIR /
*     ..
*     .. Executable Statements ..
*
*     ISPEC > 7 or ISPEC < 1: Error (return -1)
*
      IF( ISPEC.LT.1 .OR. ISPEC.GT.7 ) THEN
         ILAENV = -1
         RETURN
      END IF
*
*     ISPEC=7: Number of processors
*
      IF( ISPEC.EQ.7 ) THEN
         ILAENV = MAX( 1, NENVIR( 7 ) )
         RETURN
      END IF
*
*     ISPEC=6: Jim's crossover.
*
      IF( ISPEC.EQ.6 ) THEN
         ILAENV = INT( REAL( MAX( N1, N2 ) )*1.6E0 )
         RETURN
      END IF
*
*     ISPEC=1 through ISPEC=6: split up name into components
*
      NAMWRK = SUBNAM
      NAME1 = NAMWRK( 1: 1 )
      NAME23 = NAMWRK( 2: 3 )
      NAME46 = NAMWRK( 4: 6 )
*
*     Test version: just use number from common block
*
      ILAENV = MAX( 0, NENVIR( ISPEC ) )
      RETURN
*
*     End of ILAENV
*
      END

