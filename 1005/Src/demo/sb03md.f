      SUBROUTINE SB03MD( DICO, JOB, FACT, TRANA, N, C, LDC, A, LDA, U, 
     $                   LDU, SCALE, SEP, FERR, WR, WI, IWORK, DWORK,
     $                   LDWORK, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve for X either the real continuous-time Lyapunov equation
C
C        op(A)'*X + X*op(A) = scale*C                             (1)
C
C     or the real discrete-time Lyapunov equation
C
C        op(A)'*X*op(A) - X = scale*C                             (2)
C
C     and/or estimate an associated condition number, called separation,
C     where op(A) = A or A' (A**T) and C is symmetric (C = C').
C     (A' denotes the transpose of the matrix A.) A is N-by-N, the right
C     hand side C and the solution X are N-by-N, and scale is an output
C     scale factor, set less than or equal to 1 to avoid overflow in X.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the equation from which X is to be determined
C             as follows:
C             = 'C':  Equation (1), continuous-time case;
C             = 'D':  Equation (2), discrete-time case.
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'X':  Compute the solution only;
C             = 'S':  Compute the separation only;
C             = 'B':  Compute both the solution and the separation.
C
C     FACT    CHARACTER*1
C             Specifies whether or not the real Schur factorization
C             of the matrix A is supplied on entry, as follows:
C             = 'F':  On entry, A and U contain the factors from the
C                     real Schur factorization of the matrix A;
C             = 'N':  The Schur factorization of A will be computed
C                     and the factors will be stored in A and U.
C
C     TRANA   CHARACTER*1
C             Specifies the form of op(A) to be used, as follows:
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, X, and C.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A. If FACT = 'F', then A contains
C             an upper quasi-triangular matrix in Schur canonical form;
C             the elements below the upper Hessenberg part of the
C             array A are not referenced.
C             On exit, if INFO = 0 or INFO = N+1, the leading N-by-N
C             upper Hessenberg part of this array contains the upper
C             quasi-triangular matrix in Schur canonical form from the
C             Schur factorization of A. The contents of array A is not
C             modified if FACT = 'F'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     U       (input or output) DOUBLE PRECISION array, dimension
C             (LDU,N)
C             If FACT = 'F', then U is an input argument and on entry
C             the leading N-by-N part of this array must contain the
C             orthogonal matrix U of the real Schur factorization of A.
C             If FACT = 'N', then U is an output argument and on exit,
C             if INFO = 0 or INFO = N+1, it contains the orthogonal
C             N-by-N matrix from the real Schur factorization of A.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry with JOB = 'X' or 'B', the leading N-by-N part of
C             this array must contain the symmetric matrix C.
C             On exit with JOB = 'X' or 'B', if INFO = 0 or INFO = N+1,
C             the leading N-by-N part of C has been overwritten by the
C             symmetric solution matrix X.
C             If JOB = 'S', C is not referenced.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= 1,        if JOB = 'S';
C             LDC >= MAX(1,N), otherwise.
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C
C     SEP     (output) DOUBLE PRECISION
C             If JOB = 'S' or JOB = 'B', and INFO = 0 or INFO = N+1, SEP
C             contains the estimated separation of the matrices op(A)
C             and -op(A)', if DICO = 'C' or of op(A) and op(A)', if
C             DICO = 'D'.
C             If JOB = 'X' or N = 0, SEP is not referenced.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'B', and INFO = 0 or INFO = N+1, FERR contains an
C             estimated forward error bound for the solution X.
C             If XTRUE is the true solution, FERR bounds the relative
C             error in the computed solution, measured in the Frobenius
C             norm:  norm(X - XTRUE)/norm(XTRUE).
C             If JOB = 'X' or JOB = 'S', FERR is not referenced.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             If FACT = 'N', and INFO = 0 or INFO = N+1, WR and WI
C             contain the real and imaginary parts, respectively, of
C             the eigenvalues of A.
C             If FACT = 'F', WR and WI are not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N*N)
C             This array is not referenced if JOB = 'X'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the
C             optimal value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 1, and
C             If JOB = 'X' then
C                If FACT = 'F', LDWORK >= N*N,           for DICO = 'C';
C                               LDWORK >= MAX(N*N, 2*N), for DICO = 'D';
C                If FACT = 'N', LDWORK >= MAX(N*N, 3*N).
C             If JOB = 'S' or JOB = 'B' then
C                If FACT = 'F', LDWORK >= 2*N*N,       for DICO = 'C';
C                               LDWORK >= 2*N*N + 2*N, for DICO = 'D'.
C                If FACT = 'N', LDWORK >= MAX(2*N*N, 3*N), DICO = 'C';
C                               LDWORK >= 2*N*N + 2*N, for DICO = 'D'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, the QR algorithm failed to compute all
C                   the eigenvalues (see LAPACK Library routine DGEES);
C                   elements i+1:n of WR and WI contain eigenvalues
C                   which have converged, and A contains the partially
C                   converged Schur form;
C             = N+1:  if DICO = 'C', and the matrices A and -A' have
C                   common or very close eigenvalues, or
C                   if DICO = 'D', and matrix A has almost reciprocal
C                   eigenvalues (that is, lambda(i) = 1/lambda(j) for
C                   some i and j, where lambda(i) and lambda(j) are
C                   eigenvalues of A and i <> j); perturbed values were
C                   used to solve the equation (but the matrix A is
C                   unchanged).
C
C     METHOD
C
C     The Schur factorization of a square matrix  A  is given by
C
C        A = U*S*U'
C
C     where U is orthogonal and S is block upper triangular with 1-by-1
C     and 2-by-2 blocks on its diagonal, these blocks corresponding to
C     the eigenvalues of A, the 2-by-2 blocks being complex conjugate
C     pairs. This factorization is obtained by numerically stable
C     methods: first A is reduced to upper Hessenberg form (if FACT =
C     'N') by means of Householder transformations and then the
C     QR Algorithm is applied to reduce the Hessenberg form to S, the
C     transformation matrices being accumulated at each step to give U.
C     If A has already been factorized prior to calling the routine
C     however, then the factors U and S may be supplied and the initial
C     factorization omitted.
C                   _            _
C     If we now put C = U'CU and X = UXU' equations (1) and (2) (see
C     PURPOSE) become (for TRANS = 'N')
C          _   _    _
C        S'X + XS = C,                                               (3)
C     and
C          _    _   _
C        S'XS - X = C,                                               (4)
C
C     respectively. Partition S, C and X as
C                            _   _         _   _
C            (s    s')      (c   c')      (x   x')
C            ( 11    )  _   ( 11   )  _   ( 11   )
C        S = (       ), C = (      ), X = (      )
C            (       )      ( _    )      ( _    )
C            ( 0   S )      ( c  C )      ( x  X )
C                   1             1             1
C                _      _
C     where s  , c  and x  are either scalars or 2-by-2 matrices and s,
C            11   11     11
C     _     _
C     c and x are either (N-1) element vectors or matrices with two
C     columns. Equations (3) and (4) can then be re-written as
C           _     _        _
C        s' x   + x  s   = c                                       (3.1)
C         11 11    11 11    11
C
C          _   _           _    _
C        S'x + xs        = c - sx                                  (3.2)
C         1      11              11
C
C                                _    _
C        S'X + X S       = C - (sx' + xs')                         (3.3)
C         1 1   1 1         1
C     and
C           _       _       _
C        s' x  s  - x     = c                                      (4.1)
C         11 11 11   11      11
C
C          _     _          _    _
C        S'xs  - x        = c - sx  s                              (4.2)
C         1  11                   11 11
C
C                                _            _        _
C        S'X S - X        = C - sx  s' - [s(S'x)' + (S'x)s']       (4.3)
C         1 1 1   1          1    11         1        1
C                                                  _
C     respectively. If DICO = 'C' ['D'], then once x   has been
C                                                   11
C     found from equation (3.1) [(4.1)], equation (3.2) [(4.2)] can be
C                                        _
C     solved by forward substitution for x and then equation (3.3)
C     [(4.3)] is of the same form as (3) [(4)] but of the order (N-1) or
C     (N-2) depending upon whether s   is 1-by-1 or 2-by-2.
C                                   11
C                             _      _
C     When s   is 2-by-2 then x  and c   will be 1-by-2 matrices and s,
C           11                 11     11
C     _     _
C     x and c are matrices with two columns. In this case, equation
C     (3.1) [(4.1)] defines the three equations in the unknown elements
C        _
C     of x   and equation (3.2) [(4.2)] can then be solved by forward
C         11                 _
C     substitution, a row of x being found at each step.
C
C     REFERENCES
C
C     [1] Barraud, A.Y.                   T
C         A numerical algorithm to solve A XA - X = Q.
C         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977.
C
C     [2] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [3] Hammarling, S.J.
C         Numerical solution of the stable, non-negative definite
C         Lyapunov equation.
C         IMA J. Num. Anal., 2, pp. 303-325, 1982.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     FURTHER COMMENTS
C
C     If DICO = 'C', SEP is defined as the separation of op(A) and
C     -op(A)':
C
C            sep( op(A), -op(A)' ) = sigma_min( T )
C
C     and if DICO = 'D', SEP is defined as
C
C            sep( op(A), op(A)' ) = sigma_min( T )
C
C     where sigma_min(T) is the smallest singular value of the
C     N*N-by-N*N matrix
C
C       T = kprod( I(N), op(A)' ) + kprod( op(A)', I(N) )  (DICO = 'C'),
C
C       T = kprod( op(A)', op(A)' ) - I(N**2)              (DICO = 'D').
C
C     I(x) is an x-by-x identity matrix, and kprod denotes the Kronecker
C     product. The program estimates sigma_min(T) by the reciprocal of
C     an estimate of the 1-norm of inverse(T). The true reciprocal
C     1-norm of inverse(T) cannot differ from sigma_min(T) by more
C     than a factor of N.
C
C     When SEP is small, small changes in A, C can cause large changes
C     in the solution of the equation. An approximate bound on the
C     maximum relative error in the computed solution is
C
C                      EPS * norm(A) / SEP      (DICO = 'C'),
C
C                      EPS * norm(A)**2 / SEP   (DICO = 'D'),
C
C     where EPS is the machine precision.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, July 1997.
C     Supersedes Release 2.0 routine SB03AD by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C
C     KEYWORDS
C
C     Lyapunov equation, orthogonal transformation, real Schur form,
C     Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOB, TRANA
      INTEGER           INFO, LDA, LDC, LDU, LDWORK, N
      DOUBLE PRECISION  FERR, SCALE, SEP
C     .. Array Arguments ..
      INTEGER           IWORK( * )
      DOUBLE PRECISION  A( LDA, * ), C( LDC, * ), DWORK( * ),
     $                  U( LDU, * ), WI( * ), WR( * )
C     .. Local Scalars ..
      LOGICAL           CONT, NOFACT, NOTA, WANTBH, WANTSP, WANTX
      CHARACTER         NOTRA, NTRNST, TRANST, UPLO
      INTEGER           I, IERR, KASE, LWA, MINWRK, NN, NN2, SDIM
      DOUBLE PRECISION  EPS, EST, SCALEF
C     .. Local Arrays ..
      LOGICAL           BWORK( 1 )
C     .. External Functions ..
      LOGICAL           LSAME, SELECT
      DOUBLE PRECISION  DLAMCH, DLANHS
      EXTERNAL          DLAMCH, DLANHS, LSAME, SELECT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEES, DLACON, MB01RD, SB03MX, SB03MY,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      CONT   = LSAME( DICO,  'C' )
      WANTX  = LSAME( JOB,   'X' )
      WANTSP = LSAME( JOB,   'S' )
      WANTBH = LSAME( JOB,   'B' )
      NOFACT = LSAME( FACT,  'N' )
      NOTA   = LSAME( TRANA, 'N' )
      NN  = N*N
      NN2 = 2*NN
C
      INFO = 0
      IF( .NOT.CONT .AND. .NOT.LSAME( DICO, 'D' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTBH .AND. .NOT.WANTSP .AND. .NOT.WANTX ) THEN
         INFO = -2
      ELSE IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.NOTA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                         .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDU.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( WANTSP .AND. LDC.LT.1 .OR.
     $    .NOT.WANTSP .AND. LDC.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE
         IF ( WANTX ) THEN
            IF ( NOFACT ) THEN
               MINWRK = MAX( NN, 3*N )
            ELSE IF ( CONT ) THEN
               MINWRK = NN
            ELSE
               MINWRK = MAX( NN, 2*N )
            END IF
         ELSE
            IF ( CONT ) THEN
               IF ( NOFACT ) THEN
                  MINWRK = MAX( NN2, 3*N )
               ELSE
                  MINWRK = NN2
               END IF
            ELSE
               MINWRK = NN2 + 2*N
            END IF
         END IF
         IF( LDWORK.LT.MAX( 1, MINWRK ) )
     $      INFO = -19
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB03MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         SCALE = ONE
         IF( WANTBH )
     $      FERR  = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
      LWA = 0
C
      IF( NOFACT ) THEN
C
C        Compute the Schur factorization of A.
C        Workspace:  need   3*N;
C                    prefer larger.
C        (Note: Comments in the code beginning "Workspace:" describe the
C        minimal amount of real workspace needed at that point in the
C        code, as well as the preferred amount for good performance.
C        NB refers to the optimal block size for the immediately
C        following subroutine, as returned by ILAENV.)
C
         CALL DGEES( 'Vectors', 'Not ordered', SELECT, N, A, LDA, SDIM,
     $               WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         IF( INFO.GT.0 )
     $      RETURN
         LWA = INT( DWORK( 1 ) )
      END IF
C
      IF( .NOT.WANTSP ) THEN
C
C        Transform the right-hand side.
C        Workspace:  N*N.
C
         NTRNST = 'N'
         TRANST = 'T'
         UPLO   = 'U'
         CALL MB01RD( UPLO, TRANST, N, N, ZERO, ONE, C, LDC, U, LDU, C,
     $                LDC, DWORK, LDWORK, INFO )
C
         DO 10 I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
   10    CONTINUE
C
         LWA = MAX( LWA, NN )
C
C        Solve the transformed equation.
C        Workspace for DICO = 'D':  2*N.
C
         IF ( CONT ) THEN
            CALL SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO )
         ELSE
            CALL SB03MX( TRANA, N, A, LDA, C, LDC, SCALE, DWORK, INFO )
         END IF
         IF( INFO.GT.0 )
     $      INFO = N + 1
C
C        Transform back the solution.
C        Workspace:  N*N.
C
         CALL MB01RD( UPLO, NTRNST, N, N, ZERO, ONE, C, LDC, U, LDU, C,
     $                LDC, DWORK, LDWORK, IERR )
C
         DO 20 I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
   20    CONTINUE
C
      END IF
C
      IF( .NOT.WANTX ) THEN
C
C        Estimate the separation.
C        Workspace:  2*N*N       for DICO = 'C';
C                    2*N*N + 2*N for DICO = 'D'.
C
         IF( NOTA ) THEN
            NOTRA = 'T'
         ELSE
            NOTRA = 'N'
         END IF
C
         EST = ZERO
         KASE = 0
C        REPEAT
   30    CONTINUE
         CALL DLACON( NN, DWORK(NN+1), DWORK, IWORK, EST, KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
               IF( CONT ) THEN
                  CALL SB03MY( TRANA, N, A, LDA, DWORK, N, SCALEF,
     $                         IERR )
               ELSE
                  CALL SB03MX( TRANA, N, A, LDA, DWORK, N, SCALEF,
     $                         DWORK(NN2+1), IERR )
               END IF
            ELSE
               IF( CONT ) THEN
                  CALL SB03MY( NOTRA, N, A, LDA, DWORK, N, SCALEF,
     $                         IERR )
               ELSE
                  CALL SB03MX( NOTRA, N, A, LDA, DWORK, N, SCALEF,
     $                         DWORK(NN2+1), IERR )
               END IF
            END IF
            GO TO 30
         END IF
C        UNTIL KASE = 0
C
         SEP = SCALEF / EST
C
         IF( WANTBH ) THEN
C
C           Get the machine precision.
C
            EPS = DLAMCH( 'P' )
C
C           Compute the estimate of the relative error.
C
            IF ( CONT ) THEN
               FERR = EPS*DLANHS( 'Frobenius', N, A, LDA, DWORK )/SEP
            ELSE
               FERR = EPS*DLANHS( 'Frobenius', N, A, LDA, DWORK )**2/SEP
            END IF
         END IF
      END IF
C
      DWORK( 1 ) = DBLE( MAX( LWA, MINWRK ) )
      RETURN
C *** Last line of SB03MD ***
      END
      SUBROUTINE MB01RD( UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, A, LDA,
     $                   X, LDX, DWORK, LDWORK, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To compute the matrix formula
C        _
C        R = alpha*R + beta*op( A )*X*op( A )',
C                                                 _
C     where alpha and beta are scalars, R, X, and R are symmetric
C     matrices, A is a general matrix, and op( A ) is one of
C
C        op( A ) = A   or   op( A ) = A'.
C
C     The result is overwritten on R.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1                                         _
C             Specifies which triangles of the symmetric matrices R, R,
C             and X are given as follows:
C             = 'U':  the upper triangular part is given;
C             = 'L':  the lower triangular part is given.
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( A ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER           _
C             The order of the matrices R and R and the number of rows
C             of the matrix op( A ).  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrix X and the number of columns of the
C             the matrix op( A ).  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then R need not be
C             set before entry, except when R is identified with X in
C             the call (which is possible only in this case).
C
C     BETA    (input) DOUBLE PRECISION
C             The scalar beta. When beta is zero then A and X are not
C             referenced.
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M)
C             On entry with UPLO = 'U', the leading M-by-M upper
C             triangular part of this array must contain the upper
C             triangular part of the symmetric matrix R; the strictly
C             lower triangular part of the array is used as workspace.
C             On entry with UPLO = 'L', the leading M-by-M lower
C             triangular part of this array must contain the lower
C             triangular part of the symmetric matrix R; the strictly
C             upper triangular part of the array is used as workspace.
C             On exit, the leading M-by-M upper triangular part (if
C             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of
C             this array contains the corresponding triangular part of
C                                 _
C             the computed matrix R. If beta <> 0, the remaining
C             strictly triangular part of this array contains the
C             corresponding part of the matrix expression
C             beta*op( A )*T*op( A )', where T is the triangular matrix
C             defined in the Method section.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,M).
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,k)
C             where k is N when TRANS = 'N' and is M when TRANS = 'T' or
C             TRANS = 'C'.
C             On entry with TRANS = 'N', the leading M-by-N part of this
C             array must contain the matrix A.
C             On entry with TRANS = 'T' or TRANS = 'C', the leading
C             N-by-M part of this array must contain the matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,l),
C             where l is M when TRANS = 'N' and is N when TRANS = 'T' or
C             TRANS = 'C'.
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, if UPLO = 'U', the leading N-by-N upper
C             triangular part of this array must contain the upper
C             triangular part of the symmetric matrix X and the strictly
C             lower triangular part of the array is not referenced.
C             On entry, if UPLO = 'L', the leading N-by-N lower
C             triangular part of this array must contain the lower
C             triangular part of the symmetric matrix X and the strictly
C             upper triangular part of the array is not referenced.
C             On exit, each diagonal element of this array has half its
C             input value, but the other elements are not modified.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= MAX(1,N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, the leading M-by-N part of this
C             array (with the leading dimension MAX(1,M)) returns the
C             matrix product beta*op( A )*T, where T is the triangular
C             matrix defined in the Method section.
C             This array is not referenced when beta = 0.
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,M*N), if  beta <> 0;
C             LDWORK >= 1,          if  beta =  0.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -k, the k-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The matrix expression is efficiently evaluated taking the symmetry
C     into account. Specifically, let X = T + T', with T an upper or
C     lower triangular matrix, defined by
C
C        T = triu( X ) - (1/2)*diag( X ),  if UPLO = 'U',
C        T = tril( X ) - (1/2)*diag( X ),  if UPLO = 'L',
C
C     where triu, tril, and diag denote the upper triangular part, lower
C     triangular part, and diagonal part of X, respectively. Then,
C
C        op( A )*X*op( A )' = B + B',
C
C     where B := op( A )*T*op( A )'. Matrix B is not symmetric, but it
C     can be written as tri( B ) + stri( B ), where tri denotes the
C     triangular part specified by UPLO, and stri denotes the remaining
C     strictly triangular part. Let R = V + V', with V defined as T
C     above. Then, the required triangular part of the result can be
C     written as
C
C        alpha*V + beta*tri( B )  + beta*(stri( B ))' +
C                 alpha*diag( V ) + beta*diag( tri( B ) ).
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately
C
C                   2         2
C        3/2 x M x N + 1/2 x M
C
C     operations.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004,
C     Apr. 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANS, UPLO
      INTEGER           INFO, LDA, LDR, LDWORK, LDX, M, N
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), R(LDR,*), X(LDX,*)
C     .. Local Scalars ..
      CHARACTER*12      NTRAN
      LOGICAL           LTRANS, LUPLO
      INTEGER           J, JWORK, LDW, NROWA
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DLACPY, DLASCL, DLASET,
     $                  DSCAL, DTRMM, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      LUPLO  = LSAME( UPLO,  'U' )
      LTRANS = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
C
      IF ( LTRANS ) THEN
         NROWA = N
         NTRAN = 'No transpose'
      ELSE
         NROWA = M
         NTRAN = 'Transpose'
      END IF
C
      LDW = MAX( 1, M )
C
      IF(      ( .NOT.LUPLO  ).AND.( .NOT.LSAME( UPLO,  'L' ) ) )THEN
         INFO = -1
      ELSE IF( ( .NOT.LTRANS ).AND.( .NOT.LSAME( TRANS, 'N' ) ) )THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDR.LT.LDW ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( ( BETA.NE.ZERO .AND. LDWORK.LT.MAX( 1, M*N ) )
     $     .OR.( BETA.EQ.ZERO .AND. LDWORK.LT.1 ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01RD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      CALL DSCAL( N, HALF, X, LDX+1 )
      IF ( M.EQ.0 )
     $   RETURN
C
      IF ( BETA.EQ.ZERO .OR. N.EQ.0 ) THEN
         IF ( ALPHA.EQ.ZERO ) THEN
C
C           Special case alpha = 0.
C
            CALL DLASET( UPLO, M, M, ZERO, ZERO, R, LDR )
         ELSE
C
C           Special case beta = 0 or N = 0.
C
            IF ( ALPHA.NE.ONE )
     $         CALL DLASCL( UPLO, 0, 0, ONE, ALPHA, M, M, R, LDR, INFO )
         END IF
         RETURN
      END IF
C
C     General case: beta <> 0. Efficiently compute
C        _
C        R = alpha*R + beta*op( A )*X*op( A )',
C
C     as described in the Method section.
C
C     Compute W = beta*op( A )*T in DWORK.
C     Workspace: need M*N.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code.)
C
      IF( LTRANS ) THEN
         JWORK = 1
C
         DO 10 J = 1, N
            CALL DCOPY( M, A(J,1), LDA, DWORK(JWORK), 1 )
            JWORK = JWORK + LDW
 10      CONTINUE
C
      ELSE
         CALL DLACPY( 'Full', M, N, A, LDA, DWORK, LDW )
      END IF
C
      CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', M, N, BETA,
     $            X, LDX, DWORK, LDW )
C
C     Compute Y = alpha*V + W*op( A )' in R. First, set to zero the
C     strictly triangular part of R not specified by UPLO. That part
C     will then contain beta*stri( B ).
C
      IF ( ALPHA.NE.ZERO ) THEN
         IF ( M.GT.1 ) THEN
            IF ( LUPLO ) THEN
               CALL DLASET( 'Lower', M-1, M-1, ZERO, ZERO, R(2,1), LDR )
            ELSE
               CALL DLASET( 'Upper', M-1, M-1, ZERO, ZERO, R(1,2), LDR )
            END IF
         END IF
         CALL DSCAL( M, HALF, R, LDR+1 )
      END IF
C
      CALL DGEMM( 'No transpose', NTRAN, M, M, N, ONE, DWORK, LDW, A,
     $            LDA, ALPHA, R, LDR )
C
C     Add the term corresponding to B', with B = op( A )*T*op( A )'.
C
      IF( LUPLO ) THEN
C
         DO 20 J = 1, M
            CALL DAXPY( J, ONE, R(J,1), LDR, R(1,J), 1 )
   20    CONTINUE
C
      ELSE
C
         DO 30 J = 1, M
            CALL DAXPY( J, ONE, R(1,J), 1, R(J,1), LDR )
 30      CONTINUE
C
      END IF
C
      RETURN
C *** Last line of MB01RD ***
      END
      SUBROUTINE SB03MX( TRANA, N, A, LDA, C, LDC, SCALE, DWORK, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve the real discrete Lyapunov matrix equation
C
C            op(A)'*X*op(A) - X = scale*C
C
C     where op(A) = A or A' (A**T), A is upper quasi-triangular and C is
C     symmetric (C = C'). (A' denotes the transpose of the matrix A.)
C     A is N-by-N, the right hand side C and the solution X are N-by-N,
C     and scale is an output scale factor, set less than or equal to 1
C     to avoid overflow in X. The solution matrix X is overwritten
C     onto C.
C
C     A must be in Schur canonical form (as returned by LAPACK routines
C     DGEES or DHSEQR), that is, block upper triangular with 1-by-1 and
C     2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
C     diagonal elements equal and its off-diagonal elements of opposite
C     sign.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANA   CHARACTER*1
C             Specifies the form of op(A) to be used, as follows:
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, X, and C.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             upper quasi-triangular matrix A, in Schur canonical form.
C             The part of A below the first sub-diagonal is not
C             referenced.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading N-by-N part of this array must
C             contain the symmetric matrix C.
C             On exit, if INFO >= 0, the leading N-by-N part of this
C             array contains the symmetric solution matrix X.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,N).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if A has almost reciprocal eigenvalues; perturbed
C                   values were used to solve the equation (but the
C                   matrix A is unchanged).
C
C     METHOD
C
C     A discrete-time version of the Bartels-Stewart algorithm is used.
C     A set of equivalent linear algebraic systems of equations of order
C     at most four are formed and solved using Gaussian elimination with
C     complete pivoting.
C
C     REFERENCES
C
C     [1] Barraud, A.Y.                   T
C         A numerical algorithm to solve A XA - X = Q.
C         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977.
C
C     [2] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
C     Supersedes Release 2.0 routine SB03AZ by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom, October 1982.
C     Based on DTRLPD by P. Petkov, Tech. University of Sofia, September
C     1993.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000.
C     A. Varga, DLR Oberpfaffenhofen, March 2002.
C
C     KEYWORDS
C
C     Discrete-time system, Lyapunov equation, matrix algebra, real
C     Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOTRNA, LUPPER
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT,
     $                   MINK1N, MINK2N, MINL1N, MINL2N, NP1
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, P11, P12, P21, P22,
     $                   SCALOC, SMIN, SMLNUM, XNORM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION   VEC( 2, 2 ), X( 2, 2 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANHS
      EXTERNAL           DDOT, DLAMCH, DLANHS, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLABAD, DLALN2, DSCAL, DSYMV, SB03MV, SB04PX,
     $                   XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      NOTRNA = LSAME( TRANA, 'N' )
      LUPPER = .TRUE.
C
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                      .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03MX', -INFO )
         RETURN
      END IF
C
      SCALE = ONE
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Set constants to control overflow.
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( N*N ) / EPS
      BIGNUM = ONE / SMLNUM
C
      SMIN = MAX( SMLNUM, EPS*DLANHS( 'Max', N, A, LDA, DWORK ) )
      NP1  = N + 1
C
      IF( NOTRNA ) THEN
C
C        Solve    A'*X*A - X = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        upper-left corner column by column by
C
C          A(K,K)'*X(K,L)*A(L,L) - X(K,L) = C(K,L) - R(K,L),
C
C        where
C                    K           L-1
C          R(K,L) = SUM {A(I,K)'*SUM [X(I,J)*A(J,L)]} +
C                   I=1          J=1
C
C                    K-1
C                   {SUM [A(I,K)'*X(I,L)]}*A(L,L).
C                    I=1
C
C        Start column loop (index = L).
C        L1 (L2): column index of the first (last) row of X(K,L).
C
         LNEXT = 1
C
         DO 60 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 60
            L1 = L
            L2 = L
            IF( L.LT.N ) THEN
               IF( A( L+1, L ).NE.ZERO )
     $            L2 = L2 + 1
               LNEXT = L2 + 1
            END IF
C
C           Start row loop (index = K).
C           K1 (K2): row index of the first (last) row of X(K,L).
C
            DWORK( L1 )   = ZERO
            DWORK( N+L1 ) = ZERO
            CALL DSYMV( 'Lower', L1-1, ONE, C, LDC, A( 1, L1 ), 1, ZERO,
     $                  DWORK, 1 )
            CALL DSYMV( 'Lower', L1-1, ONE, C, LDC, A( 1, L2 ), 1, ZERO,
     $                  DWORK( NP1 ), 1 )
C
            KNEXT = L
C
            DO 50 K = L, N
               IF( K.LT.KNEXT )
     $            GO TO 50
               K1 = K
               K2 = K
               IF( K.LT.N ) THEN
                  IF( A( K+1, K ).NE.ZERO )
     $               K2 = K2 + 1
                  KNEXT = K2 + 1
               END IF
C
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1, A( 1, K1 ), 1, DWORK, 1 ) + A( L1, L1 )
     $                *DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) )
                  SCALOC = ONE
C
                  A11 = A( K1, K1 )*A( L1, L1 ) - ONE
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 10 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
C
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
C
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( K2 ) = DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ),
     $                                1 )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K2, A( 1, K1 ), 1, DWORK, 1 ) + A( L1, L1 )
     $                *DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K2, A( 1, K2 ), 1, DWORK, 1 ) + A( L1, L1 )
     $                *DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 ) )
C
                  CALL DLALN2( .TRUE., 2, 1, SMIN, A( L1, L1 ),
     $                         A( K1, K1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 20 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
C
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( N+K1 ) = DDOT( L1-1, C( K1, 1 ), LDC,
     $                                  A( 1, L2 ), 1 )
                  P11 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  P12 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1, A( 1, K1 ), 1, DWORK, 1 ) +
     $                 P11*A( L1, L1 ) + P12*A( L2, L1 ) )
C
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( K1, A( 1, K1 ), 1, DWORK( NP1 ), 1 ) +
     $                 P11*A( L1, L2 ) + P12*A( L2, L2 ) )
C
                  CALL DLALN2( .TRUE., 2, 1, SMIN, A( K1, K1 ),
     $                         A( L1, L1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 30 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   30                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
C
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( K2 ) = DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( N+K1 ) = DDOT( L1-1, C( K1, 1 ), LDC,
     $                                  A( 1, L2 ), 1 )
                  DWORK( N+K2 ) = DDOT( L1-1, C( K2, 1 ), LDC,
     $                                  A( 1, L2 ), 1 )
                  P11 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  P12 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  P21 = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  P22 = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K2, A( 1, K1 ), 1, DWORK, 1 ) +
     $                 P11*A( L1, L1 ) + P12*A( L2, L1 ) )
C
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( K2, A( 1, K1 ), 1, DWORK( NP1 ), 1 ) +
     $                 P11*A( L1, L2 ) + P12*A( L2, L2 ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K2, A( 1, K2 ), 1, DWORK, 1 ) +
     $                 P21*A( L1, L1 ) + P22*A( L2, L1 ) )
C
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( K2, A( 1, K2 ), 1, DWORK( NP1 ), 1 ) +
     $                 P21*A( L1, L2 ) + P22*A( L2, L2 ) )
C
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MV( .FALSE., LUPPER, A( K1, K1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL SB04PX( .TRUE., .FALSE., -1, 2, 2,
     $                            A( K1, K1 ), LDA, A( L1, L1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 40 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
C
   50       CONTINUE
C
   60    CONTINUE
C
      ELSE
C
C        Solve    A*X*A' - X = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        bottom-right corner column by column by
C
C            A(K,K)*X(K,L)*A(L,L)' - X(K,L) = C(K,L) - R(K,L),
C
C        where
C
C                    N            N
C          R(K,L) = SUM {A(K,I)* SUM [X(I,J)*A(L,J)']} +
C                   I=K         J=L+1
C
C                      N
C                   { SUM [A(K,J)*X(J,L)]}*A(L,L)'
C                    J=K+1
C
C        Start column loop (index = L)
C        L1 (L2): column index of the first (last) row of X(K,L)
C
         LNEXT = N
C
         DO 120 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 120
            L1 = L
            L2 = L
            IF( L.GT.1 ) THEN
               IF( A( L, L-1 ).NE.ZERO ) THEN
                  L1 = L1 - 1
                  DWORK( L1 ) = ZERO
                  DWORK( N+L1 ) = ZERO
               END IF
               LNEXT = L1 - 1
            END IF
            MINL1N = MIN( L1+1, N )
            MINL2N = MIN( L2+1, N )
C
C           Start row loop (index = K)
C           K1 (K2): row index of the first (last) row of X(K,L)
C
            IF( L2.LT.N ) THEN
               CALL DSYMV( 'Upper', N-L2, ONE, C( L2+1, L2+1 ), LDC,
     $                     A( L1, L2+1 ), LDA, ZERO, DWORK( L2+1 ), 1 )
               CALL DSYMV( 'Upper', N-L2, ONE, C( L2+1, L2+1 ), LDC,
     $                     A( L2, L2+1 ), LDA, ZERO, DWORK( NP1+L2 ), 1)
            END IF
C
            KNEXT = L
C
            DO 110 K = L, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 110
               K1 = K
               K2 = K
               IF( K.GT.1 ) THEN
                  IF( A( K, K-1 ).NE.ZERO )
     $               K1 = K1 - 1
                  KNEXT = K1 - 1
               END IF
               MINK1N = MIN( K1+1, N )
               MINK2N = MIN( K2+1, N )
C
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  DWORK( K1 ) = DDOT( N-L1, C( K1, MINL1N ), LDC,
     $                                A( L1, MINL1N ), LDA )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K1+1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L1 ), 1 )*A( L1, L1 ) )
                  SCALOC = ONE
C
                  A11 = A( K1, K1 )*A( L1, L1 ) - ONE
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 70 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   70                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
C
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
C
                  DWORK( K1 ) = DDOT( N-L1, C( K1, MINL1N ), LDC,
     $                                A( L1, MINL1N ), LDA )
                  DWORK( K2 ) = DDOT( N-L1, C( K2, MINL1N ), LDC,
     $                                A( L1, MINL1N ), LDA )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 )*A( L1, L1 ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( NP1-K1, A( K2, K1 ), LDA, DWORK( K1 ), 1 )
     $               + DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 )*A( L1, L1 ) )
C
                  CALL DLALN2( .FALSE., 2, 1, SMIN, A( L1, L1 ),
     $                         A( K1, K1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 80 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
C
                  DWORK( K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                A( L1, MINL2N ), LDA )
                  DWORK( N+K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                  A( L2, MINL2N ), LDA )
                  P11 = DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                        C( MINK1N, L1 ), 1 )
                  P12 = DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                        C( MINK1N, L2 ), 1 )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + P11*A( L1, L1 ) + P12*A( L1, L2 ) )
C
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( N+K1 ), 1)
     $               + P11*A( L2, L1 ) + P12*A( L2, L2 ) )
C
                  CALL DLALN2( .FALSE., 2, 1, SMIN, A( K1, K1 ),
     $                         A( L1, L1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 90 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
C
                  DWORK( K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                A( L1, MINL2N ), LDA )
                  DWORK( K2 ) = DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                                A( L1, MINL2N ), LDA )
                  DWORK( N+K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                  A( L2, MINL2N ), LDA )
                  DWORK( N+K2 ) = DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                                  A( L2, MINL2N ), LDA )
                  P11 = DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                        C( MINK2N, L1 ), 1 )
                  P12 = DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                        C( MINK2N, L2 ), 1 )
                  P21 = DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                        C( MINK2N, L1 ), 1 )
                  P22 = DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                        C( MINK2N, L2 ), 1 )
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + P11*A( L1, L1 ) + P12*A( L1, L2 ) )
C
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( N+K1 ),
     $                       1) + P11*A( L2, L1 ) + P12*A( L2, L2 ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( NP1-K1, A( K2, K1 ), LDA, DWORK( K1 ),
     $                       1) + P21*A( L1, L1 ) + P22*A( L1, L2 ) )
C
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( NP1-K1, A( K2, K1 ), LDA, DWORK( N+K1 ), 1)
     $               + P21*A( L2, L1 ) + P22*A( L2, L2 ) )
C
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MV( .TRUE., LUPPER, A( K1, K1 ), LDA, VEC,
     $                            2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL SB04PX( .FALSE., .TRUE., -1, 2, 2,
     $                            A( K1, K1 ), LDA, A( L1, L1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 100 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
C
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
C
  110       CONTINUE
C
  120    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of SB03MX ***
      END
      SUBROUTINE SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve the real Lyapunov matrix equation
C
C            op(A)'*X + X*op(A) = scale*C
C
C     where op(A) = A or A' (A**T), A is upper quasi-triangular and C is
C     symmetric (C = C'). (A' denotes the transpose of the matrix A.)
C     A is N-by-N, the right hand side C and the solution X are N-by-N,
C     and scale is an output scale factor, set less than or equal to 1
C     to avoid overflow in X. The solution matrix X is overwritten
C     onto C.
C
C     A must be in Schur canonical form (as returned by LAPACK routines
C     DGEES or DHSEQR), that is, block upper triangular with 1-by-1 and
C     2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
C     diagonal elements equal and its off-diagonal elements of opposite
C     sign.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANA   CHARACTER*1
C             Specifies the form of op(A) to be used, as follows:
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, X, and C.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             upper quasi-triangular matrix A, in Schur canonical form.
C             The part of A below the first sub-diagonal is not
C             referenced.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading N-by-N part of this array must
C             contain the symmetric matrix C.
C             On exit, if INFO >= 0, the leading N-by-N part of this
C             array contains the symmetric solution matrix X.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,N).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if A and -A have common or very close eigenvalues;
C                   perturbed values were used to solve the equation
C                   (but the matrix A is unchanged).
C
C     METHOD
C
C     Bartels-Stewart algorithm is used. A set of equivalent linear
C     algebraic systems of equations of order at most four are formed
C     and solved using Gaussian elimination with complete pivoting.
C
C     REFERENCES
C
C     [1] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
C     Supersedes Release 2.0 routine SB03AY by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom, October 1982.
C     Based on DTRLYP by P. Petkov, Tech. University of Sofia, September
C     1993.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C
C     KEYWORDS
C
C     Continuous-time system, Lyapunov equation, matrix algebra, real
C     Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOTRNA, LUPPER
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT,
     $                   MINK1N, MINK2N, MINL1N, MINL2N
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, SCALOC, SMIN,
     $                   SMLNUM, XNORM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANHS
      EXTERNAL           DDOT, DLAMCH, DLANHS, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLABAD, DLALN2, DLASY2, DSCAL, SB03MW, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      NOTRNA = LSAME( TRANA, 'N' )
      LUPPER = .TRUE.
C
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                      .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03MY', -INFO )
         RETURN
      END IF
C
      SCALE = ONE
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Set constants to control overflow.
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( N*N ) / EPS
      BIGNUM = ONE / SMLNUM
C
      SMIN = MAX( SMLNUM, EPS*DLANHS( 'Max', N, A, LDA, DUM ) )
C
      IF( NOTRNA ) THEN
C
C        Solve    A'*X + X*A = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        upper-left corner column by column by
C
C          A(K,K)'*X(K,L) + X(K,L)*A(L,L) = C(K,L) - R(K,L),
C
C        where
C                   K-1                    L-1
C          R(K,L) = SUM [A(I,K)'*X(I,L)] + SUM [X(K,J)*A(J,L)].
C                   I=1                    J=1
C
C        Start column loop (index = L).
C        L1 (L2): column index of the first (last) row of X(K,L).
C
         LNEXT = 1
C
         DO 60 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 60
            L1 = L
            L2 = L
            IF( L.LT.N ) THEN
               IF( A( L+1, L ).NE.ZERO )
     $            L2 = L2 + 1
               LNEXT = L2 + 1
            END IF
C
C           Start row loop (index = K).
C           K1 (K2): row index of the first (last) row of X(K,L).
C
            KNEXT = L
C
            DO 50 K = L, N
               IF( K.LT.KNEXT )
     $            GO TO 50
               K1 = K
               K2 = K
               IF( K.LT.N ) THEN
                  IF( A( K+1, K ).NE.ZERO )
     $               K2 = K2 + 1
                  KNEXT = K2 + 1
               END IF
C
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
                  SCALOC = ONE
C
                  A11 = A( K1, K1 ) + A( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 10 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
C
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ), 1 ) )
C
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 20 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
C
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L2 ), 1 ) )
C
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( L1, L1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 30 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   30                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
C
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L2 ), 1 ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ), 1 ) )
C
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 ) +
     $                 DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L2 ), 1 ) )
C
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MW( .FALSE., LUPPER, A( K1, K1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL DLASY2( .TRUE., .FALSE., 1, 2, 2, A( K1, K1 ),
     $                            LDA, A( L1, L1 ), LDA, VEC, 2, SCALOC,
     $                            X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 40 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
C
   50       CONTINUE
C
   60    CONTINUE
C
      ELSE
C
C        Solve    A*X + X*A' = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        bottom-right corner column by column by
C
C            A(K,K)*X(K,L) + X(K,L)*A(L,L)' = C(K,L) - R(K,L),
C
C        where
C                      N                     N
C            R(K,L) = SUM [A(K,I)*X(I,L)] + SUM [X(K,J)*A(L,J)'].
C                    I=K+1                 J=L+1
C
C        Start column loop (index = L).
C        L1 (L2): column index of the first (last) row of X(K,L).
C
         LNEXT = N
C
         DO 120 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 120
            L1 = L
            L2 = L
            IF( L.GT.1 ) THEN
               IF( A( L, L-1 ).NE.ZERO )
     $            L1 = L1 - 1
               LNEXT = L1 - 1
            END IF
            MINL1N = MIN( L1+1, N )
            MINL2N = MIN( L2+1, N )
C
C           Start row loop (index = K).
C           K1 (K2): row index of the first (last) row of X(K,L).
C
            KNEXT = L
C
            DO 110 K = L, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 110
               K1 = K
               K2 = K
               IF( K.GT.1 ) THEN
                  IF( A( K, K-1 ).NE.ZERO )
     $               K1 = K1 - 1
                  KNEXT = K1 - 1
               END IF
               MINK1N = MIN( K1+1, N )
               MINK2N = MIN( K2+1, N )
C
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L1 ), 1 ) +
     $                 DDOT( N-L1, C( K1, MINL1N ), LDC,
     $                       A( L1, MINL1N ), LDA ) )
                  SCALOC = ONE
C
                  A11 = A( K1, K1 ) + A( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 70 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   70                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
C
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                     C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                     A( L1, MINL2N ), LDA ) )
C
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 80 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
C
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L2 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L2, MINL2N ), LDA ) )
C
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( L1, L1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 90 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
C
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
C
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
C
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L2 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L2, MINL2N ), LDA ) )
C
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
C
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                       C( MINK2N, L2 ), 1 ) +
     $                 DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                       A( L2, MINL2N ), LDA ) )
C
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MW( .TRUE., LUPPER, A( K1, K1 ), LDA, VEC,
     $                            2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL DLASY2( .FALSE., .TRUE., 1, 2, 2, A( K1, K1 ),
     $                            LDA, A( L1, L1 ), LDA, VEC, 2, SCALOC,
     $                            X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
C
                  IF( SCALOC.NE.ONE ) THEN
C
                     DO 100 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
C
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
C
  110       CONTINUE
C
  120    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of SB03MY ***
      END
      SUBROUTINE SB03MV( LTRAN, LUPPER, T, LDT, B, LDB, SCALE, X, LDX,
     $                   XNORM, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve for the 2-by-2 symmetric matrix X in
C
C            op(T)'*X*op(T) - X = SCALE*B,
C
C     where T is 2-by-2, B is symmetric 2-by-2, and op(T) = T or T',
C     where T' denotes the transpose of T.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LTRAN   LOGICAL
C             Specifies the form of op(T) to be used, as follows:
C             = .FALSE.:  op(T) = T,
C             = .TRUE. :  op(T) = T'.
C
C     LUPPER  LOGICAL
C             Specifies which triangle of the matrix B is used, and
C             which triangle of the matrix X is computed, as follows:
C             = .TRUE. :  The upper triangular part;
C             = .FALSE.:  The lower triangular part.
C
C     Input/Output Parameters
C
C     T       (input) DOUBLE PRECISION array, dimension (LDT,2)
C             The leading 2-by-2 part of this array must contain the
C             matrix T.
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= 2.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,2)
C             On entry with LUPPER = .TRUE., the leading 2-by-2 upper
C             triangular part of this array must contain the upper
C             triangular part of the symmetric matrix B and the strictly
C             lower triangular part of B is not referenced.
C             On entry with LUPPER = .FALSE., the leading 2-by-2 lower
C             triangular part of this array must contain the lower
C             triangular part of the symmetric matrix B and the strictly
C             upper triangular part of B is not referenced.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= 2.
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor. SCALE is chosen less than or equal to 1
C             to prevent the solution overflowing.
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,2)
C             On exit with LUPPER = .TRUE., the leading 2-by-2 upper
C             triangular part of this array contains the upper
C             triangular part of the symmetric solution matrix X and the
C             strictly lower triangular part of X is not referenced.
C             On exit with LUPPER = .FALSE., the leading 2-by-2 lower
C             triangular part of this array contains the lower
C             triangular part of the symmetric solution matrix X and the
C             strictly upper triangular part of X is not referenced.
C             Note that X may be identified with B in the calling
C             statement.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= 2.
C
C     XNORM   (output) DOUBLE PRECISION
C             The infinity-norm of the solution.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if T has almost reciprocal eigenvalues, so T
C                   is perturbed to get a nonsingular equation.
C
C             NOTE: In the interests of speed, this routine does not
C                   check the inputs for errors.
C
C     METHOD
C
C     The equivalent linear algebraic system of equations is formed and
C     solved using Gaussian elimination with complete pivoting.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is stable and reliable, since Gaussian elimination
C     with complete pivoting is used.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
C     Based on DLALD2 by P. Petkov, Tech. University of Sofia, September
C     1993.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Discrete-time system, Lyapunov equation, matrix algebra.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     FOUR = 4.0D+0 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL            LTRAN, LUPPER
      INTEGER            INFO, LDB, LDT, LDX
      DOUBLE PRECISION   SCALE, XNORM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), T( LDT, * ), X( LDX, * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IP, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   EPS, SMIN, SMLNUM, TEMP, XMAX
C     ..
C     .. Local Arrays ..
      INTEGER            JPIV( 3 )
      DOUBLE PRECISION   BTMP( 3 ), T9( 3, 3 ), TMP( 3 )
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     ..
C     .. Executable Statements ..
C
C     Do not check the input parameters for errors.
C
      INFO = 0
C
C     Set constants to control overflow.
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
C
C     Solve equivalent 3-by-3 system using complete pivoting.
C     Set pivots less than SMIN to SMIN.
C
      SMIN = MAX( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ),
     $            ABS( T( 2, 1 ) ), ABS( T( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      T9( 1, 1 ) = T( 1, 1 )*T( 1, 1 ) - ONE
      T9( 2, 2 ) = T( 1, 1 )*T( 2, 2 ) + T( 1, 2 )*T( 2, 1 ) - ONE
      T9( 3, 3 ) = T( 2, 2 )*T( 2, 2 ) - ONE
      IF( LTRAN ) THEN
         T9( 1, 2 ) = T( 1, 1 )*T( 1, 2 ) + T( 1, 1 )*T( 1, 2 )
         T9( 1, 3 ) = T( 1, 2 )*T( 1, 2 )
         T9( 2, 1 ) = T( 1, 1 )*T( 2, 1 )
         T9( 2, 3 ) = T( 1, 2 )*T( 2, 2 )
         T9( 3, 1 ) = T( 2, 1 )*T( 2, 1 )
         T9( 3, 2 ) = T( 2, 1 )*T( 2, 2 ) + T( 2, 1 )*T( 2, 2 )
      ELSE
         T9( 1, 2 ) = T( 1, 1 )*T( 2, 1 ) + T( 1, 1 )*T( 2, 1 )
         T9( 1, 3 ) = T( 2, 1 )*T( 2, 1 )
         T9( 2, 1 ) = T( 1, 1 )*T( 1, 2 )
         T9( 2, 3 ) = T( 2, 1 )*T( 2, 2 )
         T9( 3, 1 ) = T( 1, 2 )*T( 1, 2 )
         T9( 3, 2 ) = T( 1, 2 )*T( 2, 2 ) + T( 1, 2 )*T( 2, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      IF ( LUPPER ) THEN
         BTMP( 2 ) = B( 1, 2 )
      ELSE
         BTMP( 2 ) = B( 2, 1 )
      END IF
      BTMP( 3 ) = B( 2, 2 )
C
C     Perform elimination.
C
      DO 50 I = 1, 2
         XMAX = ZERO
C
         DO 20 IP = I, 3
C
            DO 10 JP = I, 3
               IF( ABS( T9( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T9( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   10       CONTINUE
C
   20    CONTINUE
C
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 3, T9( IPSV, 1 ), 3, T9( I, 1 ), 3 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 3, T9( 1, JPSV ), 1, T9( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T9( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T9( I, I ) = SMIN
         END IF
C
         DO 40 J = I + 1, 3
            T9( J, I ) = T9( J, I ) / T9( I, I )
            BTMP( J ) = BTMP( J ) - T9( J, I )*BTMP( I )
C
            DO 30 K = I + 1, 3
               T9( J, K ) = T9( J, K ) - T9( J, I )*T9( I, K )
   30       CONTINUE
C
   40    CONTINUE
C
   50 CONTINUE
C
      IF( ABS( T9( 3, 3 ) ).LT.SMIN )
     $   T9( 3, 3 ) = SMIN
      SCALE = ONE
      IF( ( FOUR*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T9( 1, 1 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T9( 2, 2 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T9( 3, 3 ) ) ) THEN
         SCALE = ( ONE / FOUR ) / MAX( ABS( BTMP( 1 ) ),
     $               ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
      END IF
C
      DO 70 I = 1, 3
         K = 4 - I
         TEMP = ONE / T9( K, K )
         TMP( K ) = BTMP( K )*TEMP
C
         DO 60 J = K + 1, 3
            TMP( K ) = TMP( K ) - ( TEMP*T9( K, J ) )*TMP( J )
  60     CONTINUE
C
  70  CONTINUE
C
      DO 80 I = 1, 2
         IF( JPIV( 3-I ).NE.3-I ) THEN
            TEMP = TMP( 3-I )
            TMP( 3-I ) = TMP( JPIV( 3-I ) )
            TMP( JPIV( 3-I ) ) = TEMP
         END IF
  80  CONTINUE
C
      X( 1, 1 ) = TMP( 1 )
      IF ( LUPPER ) THEN
         X( 1, 2 ) = TMP( 2 )
      ELSE
         X( 2, 1 ) = TMP( 2 )
      END IF
      X( 2, 2 ) = TMP( 3 )
      XNORM = MAX( ABS( TMP( 1 ) ) + ABS( TMP( 2 ) ),
     $             ABS( TMP( 2 ) ) + ABS( TMP( 3 ) ) )
C
      RETURN
C *** Last line of SB03MV ***
      END
      SUBROUTINE SB03MW( LTRAN, LUPPER, T, LDT, B, LDB, SCALE, X, LDX,
     $                   XNORM, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve for the 2-by-2 symmetric matrix X in
C
C            op(T)'*X + X*op(T) = SCALE*B,
C
C     where T is 2-by-2, B is symmetric 2-by-2, and op(T) = T or T',
C     where T' denotes the transpose of T.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LTRAN   LOGICAL
C             Specifies the form of op(T) to be used, as follows:
C             = .FALSE.:  op(T) = T,
C             = .TRUE. :  op(T) = T'.
C
C     LUPPER  LOGICAL
C             Specifies which triangle of the matrix B is used, and
C             which triangle of the matrix X is computed, as follows:
C             = .TRUE. :  The upper triangular part;
C             = .FALSE.:  The lower triangular part.
C
C     Input/Output Parameters
C
C     T       (input) DOUBLE PRECISION array, dimension (LDT,2)
C             The leading 2-by-2 part of this array must contain the
C             matrix T.
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= 2.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,2)
C             On entry with LUPPER = .TRUE., the leading 2-by-2 upper
C             triangular part of this array must contain the upper
C             triangular part of the symmetric matrix B and the strictly
C             lower triangular part of B is not referenced.
C             On entry with LUPPER = .FALSE., the leading 2-by-2 lower
C             triangular part of this array must contain the lower
C             triangular part of the symmetric matrix B and the strictly
C             upper triangular part of B is not referenced.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= 2.
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor. SCALE is chosen less than or equal to 1
C             to prevent the solution overflowing.
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,2)
C             On exit with LUPPER = .TRUE., the leading 2-by-2 upper
C             triangular part of this array contains the upper
C             triangular part of the symmetric solution matrix X and the
C             strictly lower triangular part of X is not referenced.
C             On exit with LUPPER = .FALSE., the leading 2-by-2 lower
C             triangular part of this array contains the lower
C             triangular part of the symmetric solution matrix X and the
C             strictly upper triangular part of X is not referenced.
C             Note that X may be identified with B in the calling
C             statement.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= 2.
C
C     XNORM   (output) DOUBLE PRECISION
C             The infinity-norm of the solution.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if T and -T have too close eigenvalues, so T
C                   is perturbed to get a nonsingular equation.
C
C             NOTE: In the interests of speed, this routine does not
C                   check the inputs for errors.
C
C     METHOD
C
C     The equivalent linear algebraic system of equations is formed and
C     solved using Gaussian elimination with complete pivoting.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is stable and reliable, since Gaussian elimination
C     with complete pivoting is used.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
C     Based on DLALY2 by P. Petkov, Tech. University of Sofia, September
C     1993.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Continuous-time system, Lyapunov equation, matrix algebra.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     FOUR = 4.0D+0 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL            LTRAN, LUPPER
      INTEGER            INFO, LDB, LDT, LDX
      DOUBLE PRECISION   SCALE, XNORM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), T( LDT, * ), X( LDX, * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IP, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   EPS, SMIN, SMLNUM, TEMP, XMAX
C     ..
C     .. Local Arrays ..
      INTEGER            JPIV( 3 )
      DOUBLE PRECISION   BTMP( 3 ), T9( 3, 3 ), TMP( 3 )
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     ..
C     .. Executable Statements ..
C
C     Do not check the input parameters for errors
C
      INFO = 0
C
C     Set constants to control overflow
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
C
C     Solve equivalent 3-by-3 system using complete pivoting.
C     Set pivots less than SMIN to SMIN.
C
      SMIN = MAX( MAX( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ),
     $                 ABS( T( 2, 1 ) ), ABS( T( 2, 2 ) ) )*EPS,
     $            SMLNUM )
      T9( 1, 3 ) = ZERO
      T9( 3, 1 ) = ZERO
      T9( 1, 1 ) = T( 1, 1 )
      T9( 2, 2 ) = T( 1, 1 ) + T( 2, 2 )
      T9( 3, 3 ) = T( 2, 2 )
      IF( LTRAN ) THEN
         T9( 1, 2 ) = T( 1, 2 )
         T9( 2, 1 ) = T( 2, 1 )
         T9( 2, 3 ) = T( 1, 2 )
         T9( 3, 2 ) = T( 2, 1 )
      ELSE
         T9( 1, 2 ) = T( 2, 1 )
         T9( 2, 1 ) = T( 1, 2 )
         T9( 2, 3 ) = T( 2, 1 )
         T9( 3, 2 ) = T( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )/TWO
      IF ( LUPPER ) THEN
         BTMP( 2 ) = B( 1, 2 )
      ELSE
         BTMP( 2 ) = B( 2, 1 )
      END IF
      BTMP( 3 ) = B( 2, 2 )/TWO
C
C     Perform elimination
C
      DO 50 I = 1, 2
         XMAX = ZERO
C
         DO 20 IP = I, 3
C
            DO 10 JP = I, 3
               IF( ABS( T9( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T9( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   10       CONTINUE
C
   20    CONTINUE
C
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 3, T9( IPSV, 1 ), 3, T9( I, 1 ), 3 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 3, T9( 1, JPSV ), 1, T9( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T9( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T9( I, I ) = SMIN
         END IF
C
         DO 40 J = I + 1, 3
            T9( J, I ) = T9( J, I ) / T9( I, I )
            BTMP( J ) = BTMP( J ) - T9( J, I )*BTMP( I )
C
            DO 30 K = I + 1, 3
               T9( J, K ) = T9( J, K ) - T9( J, I )*T9( I, K )
   30       CONTINUE
C
   40    CONTINUE
C
   50 CONTINUE
C
      IF( ABS( T9( 3, 3 ) ).LT.SMIN )
     $   T9( 3, 3 ) = SMIN
      SCALE = ONE
      IF( ( FOUR*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T9( 1, 1 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T9( 2, 2 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T9( 3, 3 ) ) ) THEN
         SCALE = ( ONE / FOUR ) / MAX( ABS( BTMP( 1 ) ),
     $               ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
      END IF
C
      DO 70 I = 1, 3
         K = 4 - I
         TEMP = ONE / T9( K, K )
         TMP( K ) = BTMP( K )*TEMP
C
         DO 60 J = K + 1, 3
            TMP( K ) = TMP( K ) - ( TEMP*T9( K, J ) )*TMP( J )
  60     CONTINUE
C
  70  CONTINUE
C
      DO 80 I = 1, 2
         IF( JPIV( 3-I ).NE.3-I ) THEN
            TEMP = TMP( 3-I )
            TMP( 3-I ) = TMP( JPIV( 3-I ) )
            TMP( JPIV( 3-I ) ) = TEMP
         END IF
  80  CONTINUE
C
      X( 1, 1 ) = TMP( 1 )
      IF ( LUPPER ) THEN
         X( 1, 2 ) = TMP( 2 )
      ELSE
         X( 2, 1 ) = TMP( 2 )
      END IF
      X( 2, 2 ) = TMP( 3 )
      XNORM = MAX( ABS( TMP( 1 ) ) + ABS( TMP( 2 ) ),
     $             ABS( TMP( 2 ) ) + ABS( TMP( 3 ) ) )
C
      RETURN
C *** Last line of SB03MW ***
      END
      SUBROUTINE SB04PX( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
     $                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve for the N1-by-N2 matrix X, 1 <= N1,N2 <= 2, in
C
C            op(TL)*X*op(TR) + ISGN*X = SCALE*B,
C
C     where TL is N1-by-N1, TR is N2-by-N2, B is N1-by-N2, and ISGN = 1
C     or -1.  op(T) = T or T', where T' denotes the transpose of T.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LTRANL  LOGICAL
C             Specifies the form of op(TL) to be used, as follows:
C             = .FALSE.:  op(TL) = TL,
C             = .TRUE. :  op(TL) = TL'.
C
C     LTRANR  LOGICAL
C             Specifies the form of op(TR) to be used, as follows:
C             = .FALSE.:  op(TR) = TR,
C             = .TRUE. :  op(TR) = TR'.
C
C     ISGN    INTEGER
C             Specifies the sign of the equation as described before.
C             ISGN may only be 1 or -1.
C
C     Input/Output Parameters
C
C     N1      (input) INTEGER
C             The order of matrix TL.  N1 may only be 0, 1 or 2.
C
C     N2      (input) INTEGER
C             The order of matrix TR.  N2 may only be 0, 1 or 2.
C
C     TL      (input) DOUBLE PRECISION array, dimension (LDTL,N1)
C             The leading N1-by-N1 part of this array must contain the
C             matrix TL.
C
C     LDTL    INTEGER
C             The leading dimension of array TL.  LDTL >= MAX(1,N1).
C
C     TR      (input) DOUBLE PRECISION array, dimension (LDTR,N2)
C             The leading N2-by-N2 part of this array must contain the
C             matrix TR.
C
C     LDTR    INTEGER
C             The leading dimension of array TR.  LDTR >= MAX(1,N2).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,N2)
C             The leading N1-by-N2 part of this array must contain the
C             right-hand side of the equation.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N1).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor. SCALE is chosen less than or equal to 1
C             to prevent the solution overflowing.
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,N2)
C             The leading N1-by-N2 part of this array contains the
C             solution of the equation.
C             Note that X may be identified with B in the calling
C             statement.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= MAX(1,N1).
C
C     XNORM   (output) DOUBLE PRECISION
C             The infinity-norm of the solution.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if TL and -ISGN*TR have almost reciprocal
C                   eigenvalues, so TL or TR is perturbed to get a
C                   nonsingular equation.
C
C             NOTE: In the interests of speed, this routine does not
C                   check the inputs for errors.
C
C     METHOD
C
C     The equivalent linear algebraic system of equations is formed and
C     solved using Gaussian elimination with complete pivoting.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is stable and reliable, since Gaussian elimination
C     with complete pivoting is used.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 2000.
C     This is a modification and slightly more efficient version of
C     SLICOT Library routine SB03MU.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Discrete-time system, Sylvester equation, matrix algebra.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                     TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
     $                   X( LDX, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1,
     $                   TEMP, U11, U12, U22, XMAX
C     ..
C     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ),
     $                   LOCU22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
C     ..
C     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     ..
C     .. Data statements ..
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / ,
     $                   LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
C     ..
C     .. Executable Statements ..
C
C     Do not check the input parameters for errors.
C
      INFO  = 0
      SCALE = ONE
C
C     Quick return if possible.
C
      IF( N1.EQ.0 .OR. N2.EQ.0 ) THEN
         XNORM = ZERO
         RETURN
      END IF
C
C     Set constants to control overflow.
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SGN = ISGN
C
      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K
C
C     1-by-1: TL11*X*TR11 + ISGN*X = B11.
C
   10 CONTINUE
      TAU1 = TL( 1, 1 )*TR( 1, 1 ) + SGN
      BET  = ABS( TAU1 )
      IF( BET.LE.SMLNUM ) THEN
         TAU1 = SMLNUM
         BET  = SMLNUM
         INFO = 1
      END IF
C
      GAM = ABS( B( 1, 1 ) )
      IF( SMLNUM*GAM.GT.BET )
     $   SCALE = ONE / GAM
C
      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN
C
C     1-by-2:
C     TL11*[X11 X12]*op[TR11 TR12] + ISGN*[X11 X12] = [B11 B12].
C                      [TR21 TR22]
C
   20 CONTINUE
C
      SMIN = MAX( MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
     $                 ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
     $                *ABS( TL( 1, 1 ) )*EPS,
     $            SMLNUM )
      TMP( 1 ) = TL( 1, 1 )*TR( 1, 1 ) + SGN
      TMP( 4 ) = TL( 1, 1 )*TR( 2, 2 ) + SGN
      IF( LTRANR ) THEN
         TMP( 2 ) = TL( 1, 1 )*TR( 2, 1 )
         TMP( 3 ) = TL( 1, 1 )*TR( 1, 2 )
      ELSE
         TMP( 2 ) = TL( 1, 1 )*TR( 1, 2 )
         TMP( 3 ) = TL( 1, 1 )*TR( 2, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40
C
C     2-by-1:
C     op[TL11 TL12]*[X11]*TR11 + ISGN*[X11] = [B11].
C       [TL21 TL22] [X21]             [X21]   [B21]
C
   30 CONTINUE
      SMIN = MAX( MAX( ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
     $                 ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
     $                *ABS( TR( 1, 1 ) )*EPS,
     $            SMLNUM )
      TMP( 1 ) = TL( 1, 1 )*TR( 1, 1 ) + SGN
      TMP( 4 ) = TL( 2, 2 )*TR( 1, 1 ) + SGN
      IF( LTRANL ) THEN
         TMP( 2 ) = TL( 1, 2 )*TR( 1, 1 )
         TMP( 3 ) = TL( 2, 1 )*TR( 1, 1 )
      ELSE
         TMP( 2 ) = TL( 2, 1 )*TR( 1, 1 )
         TMP( 3 ) = TL( 1, 2 )*TR( 1, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
   40 CONTINUE
C
C     Solve 2-by-2 system using complete pivoting.
C     Set pivots less than SMIN to SMIN.
C
      IPIV = IDAMAX( 4, TMP, 1 )
      U11  = TMP( IPIV )
      IF( ABS( U11 ).LE.SMIN ) THEN
         INFO = 1
         U11  = SMIN
      END IF
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      IF( ABS( U22 ).LE.SMIN ) THEN
         INFO = 1
         U22  = SMIN
      END IF
      IF( BSWAP ) THEN
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      ELSE
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      END IF
      IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR.
     $    ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      END IF
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      IF( XSWAP ) THEN
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      END IF
      X( 1, 1 ) = X2( 1 )
      IF( N1.EQ.1 ) THEN
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X2( 1 ) ) + ABS( X2( 2 ) )
      ELSE
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X2( 1 ) ), ABS( X2( 2 ) ) )
      END IF
      RETURN
C
C     2-by-2:
C     op[TL11 TL12]*[X11 X12]*op[TR11 TR12] + ISGN*[X11 X12] = [B11 B12]
C       [TL21 TL22] [X21 X22]   [TR21 TR22]        [X21 X22]   [B21 B22]
C
C     Solve equivalent 4-by-4 system using complete pivoting.
C     Set pivots less than SMIN to SMIN.
C
   50 CONTINUE
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
     $            ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
      SMIN = MAX( ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
     $            ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )*SMIN
      SMIN = MAX( EPS*SMIN, SMLNUM )
      T16( 1, 1 ) = TL( 1, 1 )*TR( 1, 1 ) + SGN
      T16( 2, 2 ) = TL( 2, 2 )*TR( 1, 1 ) + SGN
      T16( 3, 3 ) = TL( 1, 1 )*TR( 2, 2 ) + SGN
      T16( 4, 4 ) = TL( 2, 2 )*TR( 2, 2 ) + SGN
      IF( LTRANL ) THEN
         T16( 1, 2 ) = TL( 2, 1 )*TR( 1, 1 )
         T16( 2, 1 ) = TL( 1, 2 )*TR( 1, 1 )
         T16( 3, 4 ) = TL( 2, 1 )*TR( 2, 2 )
         T16( 4, 3 ) = TL( 1, 2 )*TR( 2, 2 )
      ELSE
         T16( 1, 2 ) = TL( 1, 2 )*TR( 1, 1 )
         T16( 2, 1 ) = TL( 2, 1 )*TR( 1, 1 )
         T16( 3, 4 ) = TL( 1, 2 )*TR( 2, 2 )
         T16( 4, 3 ) = TL( 2, 1 )*TR( 2, 2 )
      END IF
      IF( LTRANR ) THEN
         T16( 1, 3 ) = TL( 1, 1 )*TR( 1, 2 )
         T16( 2, 4 ) = TL( 2, 2 )*TR( 1, 2 )
         T16( 3, 1 ) = TL( 1, 1 )*TR( 2, 1 )
         T16( 4, 2 ) = TL( 2, 2 )*TR( 2, 1 )
      ELSE
         T16( 1, 3 ) = TL( 1, 1 )*TR( 2, 1 )
         T16( 2, 4 ) = TL( 2, 2 )*TR( 2, 1 )
         T16( 3, 1 ) = TL( 1, 1 )*TR( 1, 2 )
         T16( 4, 2 ) = TL( 2, 2 )*TR( 1, 2 )
      END IF
      IF( LTRANL .AND. LTRANR ) THEN
         T16( 1, 4 ) = TL( 2, 1 )*TR( 1, 2 )
         T16( 2, 3 ) = TL( 1, 2 )*TR( 1, 2 )
         T16( 3, 2 ) = TL( 2, 1 )*TR( 2, 1 )
         T16( 4, 1 ) = TL( 1, 2 )*TR( 2, 1 )
      ELSE IF( LTRANL .AND. .NOT.LTRANR ) THEN
         T16( 1, 4 ) = TL( 2, 1 )*TR( 2, 1 )
         T16( 2, 3 ) = TL( 1, 2 )*TR( 2, 1 )
         T16( 3, 2 ) = TL( 2, 1 )*TR( 1, 2 )
         T16( 4, 1 ) = TL( 1, 2 )*TR( 1, 2 )
      ELSE IF( .NOT.LTRANL .AND. LTRANR ) THEN
          T16( 1, 4 ) = TL( 1, 2 )*TR( 1, 2 )
          T16( 2, 3 ) = TL( 2, 1 )*TR( 1, 2 )
          T16( 3, 2 ) = TL( 1, 2 )*TR( 2, 1 )
          T16( 4, 1 ) = TL( 2, 1 )*TR( 2, 1 )
      ELSE
          T16( 1, 4 ) = TL( 1, 2 )*TR( 2, 1 )
          T16( 2, 3 ) = TL( 2, 1 )*TR( 2, 1 )
          T16( 3, 2 ) = TL( 1, 2 )*TR( 1, 2 )
          T16( 4, 1 ) = TL( 2, 1 )*TR( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )
C
C     Perform elimination.
C
      DO 100 I = 1, 3
         XMAX = ZERO
C
         DO 70 IP = I, 4
C
            DO 60 JP = I, 4
               IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   60       CONTINUE
C
   70    CONTINUE
C
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T16( I, I ) = SMIN
         END IF
C
         DO 90 J = I + 1, 4
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
C
            DO 80 K = I + 1, 4
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
   80       CONTINUE
C
   90    CONTINUE
C
  100 CONTINUE
C
      IF( ABS( T16( 4, 4 ) ).LT.SMIN )
     $   T16( 4, 4 ) = SMIN
      IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
         SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ),
     $                ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ),
     $                ABS( BTMP( 4 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
         BTMP( 4 ) = BTMP( 4 )*SCALE
      END IF
C
      DO 120 I = 1, 4
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
C
         DO 110 J = K + 1, 4
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
  110    CONTINUE
C
  120 CONTINUE
C
      DO 130 I = 1, 3
         IF( JPIV( 4-I ).NE.4-I ) THEN
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         END IF
  130 CONTINUE
C
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) ) + ABS( TMP( 3 ) ),
     $             ABS( TMP( 2 ) ) + ABS( TMP( 4 ) ) )
C
      RETURN
C *** Last line of SB04PX ***
      END
      LOGICAL FUNCTION  SELECT( PAR1, PAR2 )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     Void logical function for DGEES.
C
      DOUBLE PRECISION  PAR1, PAR2
C
      SELECT = .TRUE.
      RETURN
      END
