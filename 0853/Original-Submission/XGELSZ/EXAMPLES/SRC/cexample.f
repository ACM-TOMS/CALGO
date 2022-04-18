      PROGRAM EXAMPLE
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
*     Example -- a simple example that demonstrates the use of CGELSZ
*
      INTEGER            M, N, NRHS, LWORK
      PARAMETER          ( M = 5, N = 3, NRHS = 1, LWORK = 32 )
*
      INTEGER            RANK, INFO
      REAL               RCOND
      INTEGER            JPVT( N )
      COMPLEX            A( M, N ), B( M, N ), WORK( LWORK )
      REAL               RWORK( 2*N )
*
      INTEGER            I, J
*     external subroutines
      EXTERNAL           CGELSZ
*
*     zero out the JPVT vector
      DO 10 I = 1, N
         JPVT( I ) = 0
   10 CONTINUE
*
*     setup the A matrix
      A( 1, 1 ) = 3
      A( 1, 2 ) = 4
      A( 1, 3 ) = 5
      A( 2, 1 ) = 5
      A( 2, 2 ) = 6
      A( 2, 3 ) = 7
      A( 3, 1 ) = 6
      A( 3, 2 ) = 8
      A( 3, 3 ) = 10
      A( 4, 1 ) = 10
      A( 4, 2 ) = 12
      A( 4, 3 ) = 14
      A( 5, 1 ) = 11
      A( 5, 2 ) = 14
      A( 5, 3 ) = 17
*
*     setup the B (rhs) matrix
      B( 1, 1 ) = 300
      B( 2, 1 ) = 600
      B( 3, 1 ) = 900
      B( 4, 1 ) = 1200
      B( 5, 1 ) = 1500
*
      WRITE( *, * )
      WRITE( *, * )'Solve A x = b for a singular matrix A where '
      WRITE( *, * )
      WRITE( *, * )'A is '
*
      DO 20 I = 1, M
         WRITE( *, '(3(A,G8.2,A,G8.2,A))' )( '(', REAL( A( I, J ) ),
     $      ',', AIMAG( A( I, J ) ), ')  ', J = 1, N )
   20 CONTINUE
      WRITE( *, * )
      WRITE( *, * )'and b is '
*
      DO 30 I = 1, M
         WRITE( *, '(A,G13.5,A,G13.5,A)' )'(', REAL( B( I, 1 ) ), ',',
     $      AIMAG( B( I, 1 ) ), ')  '
   30 CONTINUE
*
*     call cgelsz
      INFO = 0
      RCOND = 1.0E-4
      CALL CGELSZ( M, N, NRHS, A, M, B, M, JPVT, RCOND, RANK, WORK,
     $             LWORK, RWORK, INFO )
*
*     print the results
      WRITE( *, * )
      WRITE( *, * )'The calculated rank is', RANK
      WRITE( *, * )'The true rank is      ', 2
*
      WRITE( *, * )
      WRITE( *, * )'The calculated solution is'
      DO 40 I = 1, N
         WRITE( *, * )( B( I, J ), J = 1, NRHS )
   40 CONTINUE
*
*      The true solution (the generalized inverse of A times B)
*
      B( 1, 1 ) = -13
      B( 2, 1 ) = 29
      B( 3, 1 ) = 71
      WRITE( *, * )
      WRITE( *, * )'The true solution ( pinv(A)*b ) is'
      DO 50 I = 1, N
         WRITE( *, * )( B( I, J ), J = 1, NRHS )
   50 CONTINUE
      WRITE( *, * )
*
      END
