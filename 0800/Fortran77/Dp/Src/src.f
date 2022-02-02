      SUBROUTINE DCROOT(XR,XI,YR,YI)
C
C     PURPOSE
C
C     To compute the complex square root YR + i*YI of a complex number
C     XR + i*XI  in real arithmetic.  The branch is chosen so that
C     YR .GE. 0.0  and  SIGN(YI) .EQ. SIGN(XI).
C
C     ARGUMENTS
C
C     XR      (input) DOUBLE PRECISION
C     XI      (input) DOUBLE PRECISION
C             These scalars define the real and imaginary part of the
C             complex number of which the root is sought.
C
C     YR      (input) DOUBLE PRECISION
C     YI      (input) DOUBLE PRECISION
C             These scalars define the real and imaginary part of the
C             complex root.
C
C     REFERENCE
C
C     Adapted from EISPACK subroutine CSROOT.
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C
C     REVISIONS
C
C     1993, November 15.
C     1998, August 26.
C
C***********************************************************************
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION XI,XR,YI,YR
C     ..
C     .. Parameters ..
C
      DOUBLE PRECISION ZERO,HALF
      PARAMETER (ZERO=0.0D0,HALF=1.0D0/2.0D0)
C     ..
C     .. Local Scalars ..
C
      DOUBLE PRECISION S
C     ..
C     .. External Subroutines ..
C     . LAPACK .
C
C     ..
C     .. Intrinsic Functions ..
C
C***********************************************************************
      INTRINSIC ABS,SQRT
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAPY2
      EXTERNAL DLAPY2
C     ..
C
      S = SQRT(HALF* (DLAPY2(XR,XI)+ABS(XR)))
      IF (XR.GE.ZERO) YR = S
      IF (XI.LT.ZERO) S = -S
      IF (XR.LE.ZERO) YI = S
      IF (XR.LT.ZERO) YR = HALF* (XI/YI)
      IF (XR.GT.ZERO) YI = HALF* (XI/YR)
C
      RETURN

      END
      SUBROUTINE DHABL(JOBSCL,N,A,LDA,QG,LDQG,D,RWORK,IERR)
C
C     .. Scalar Arguments ..
      INTEGER IERR,LDA,LDQG,N
      CHARACTER JOBSCL
C     ..
C     .. Array Arguments ..
C
C     PURPOSE
C
C     Perform a symplectic scaling on the Hamiltonian matrix
C
C                 ( A    G  )
C     (*)     H = (       T ),
C                 ( Q   -A  )
C
C     i.e., perform either the symplectic scaling transformation
C
C                                     -1
C                    ( A'   G'  )   ( D   0 ) ( A   G ) ( D  0   )
C     (**)    H' <-- (        T ) = (       ) (      T) (     -1 )
C                    ( Q'  -A'  )   ( 0   D ) ( Q  -A ) ( 0  D   )
C
C     where D is a diagonal scaling matrix, or the symplectic norm
C     scaling transformation
C
C                     ( A''   G''  )    1  (  A   G/tau )
C     (***)   H'' <-- (          T ) = --- (          T )
C                     ( Q''  -A''  )   tau (tau Q   -A  )
C
C     where tau is a real scalar.  Note that if tau is not equal to 1
C     then (***) is NOT a similarity transformation.  The eigenvalues
C     of H are then tau times the eigenvalues of H''.
C
C     For symplectic scaling (**), D is chosen to give the rows and
C     columns of A' approximately equal 1-norms and to give Q' and G'
C     approximately equal norms.  (See METHOD below for details.) For
C     norm scaling, tau = MAX(1, ||A||, ||G||, ||Q||) where ||.||
C     denotes the 1-norm (column sum norm).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBSCL  CHARACTER*1
C             Indicates which scaling strategy is used.
C             = 'A' or 'a':  do the symplectic scaling (**);
C             = 'B' or 'b':  do the norm scaling (***);
C             = 'N' or 'n': do nothing; set IERR = 0 and return.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  Normally, N.GE.1.
C             If N .EQ. 0, then IERR is set to zero and DHABL returns
C             without referencing any other argument.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input the leading N-by-N part of this array contains
C             the upper left block A of the Hamiltonian matrix H in
C             (*).
C             On output A has been overwritten by the leading N-by-N
C             part of the scaled Hamiltonian matrix H' in (**) or H'' in
C             (***) depending on the setting of JOBSCL.
C
C     LDA     (input) INTEGER
C             The leading dimension of array A as declared in the
C             calling program.  LDA .GE. N.
C
C     QG      (input/output) DOUBLE PRECISION array, 
C             dimension (LDQG,N+1)
C             On input QG contains the upper right symmetric block G and
C             the lower left symmetric block Q of the Hamiltonian matrix
C             in (*).  On output these are overwritten by the
C             corresponding blocks of the scaled Hamiltonian matrix H'
C             in (**) or H'' in (***) depending on the setting of
C             JOBSCL.  The diagonal and lower triangle of the submatrix
C             in columns 1 to N contain the lower triangular part of Q.
C             The diagonal and upper triangle of the submatrix in
C             columns  2 to N+1 contain the upper triangle of G.
C             So, if I. GE. J, then Q(i,j) is stored in QG(i,j)
C             and G(i,j) is in QG(i,j+1).
C
C     LDQG    (input) INTEGER
C             The leading dimension of array QG just as declared in the
C             calling procedure.  LDQG .GE. N.
C
C     D       (output) DOUBLE PRECISION array, dimension nd
C             If JOBSCL = 'A' or 'a', then nd = N and on output, D
C             contains the diagonal elements of the diagonal scaling
C             matrix in (**).
C             If JOBSCL = 'B' or 'b', then nd = 1 and D(1) is set to tau
C             from (***). In this case, no other elements of D are
C             referenced.
C
C     Workspace
C
C     RWORK   DOUBLE PRECISION array, dimension at least N
C
C     Error Indicator
C
C     IERR    INTEGER
C             = 0:  successful exit;
C             < 0:  if IERR = -j, then the j-th argument had an illegal
C                   value;
C
C     METHOD
C
C     1. Symplectic scaling (JOBSCL = 'A' or 'a'):
C       First, LAPACK subroutine DGEBAL is used to equilibrate the
C       1-norms of the rows and columns of A using a diagonal scaling
C       matrix D_A.  Then, H is similarily transformed by the symplectic
C       diagonal matrix D1 = diag(D_A,D_A**(-1)).  Next, the
C       off-diagonal blocks of the resulting Hamiltonian matrix are
C       equilibrated in the 1-norm using the symplectic diagonal matrix
C       D2 of the form
C
C                   ( I/rho    0   )
C              D2 = (              )
C                   (   I    rho*I )
C
C       where rho is a real scalar. Thus, in (**), D = D1*D2.
C
C     2. Norm scaling (JOBSCL = 'B' or 'b'):
C       The norm of the matrices A and G of (*) is reduced by setting
C       A := A/tau  and  G := G/(tau**2) where tau is the power of the
C       base of the arithmetic closest to MAX(1, ||A||, ||G||, ||Q||)
C       and ||.|| denotes the 1-norm.
C
C     NUMERICAL ASPECTS
C
C     For symplectic scaling, the complexity of the used algorithms is
C     hard to estimate and depends upon how well the rows and columns of
C     A in (*) are equilibrated.  In one sweep, each row/column of A is
C     scaled once, i.e., the cost of one sweep is N**2 multiplications.
C     Usually, 3-6 sweeps are enough to equilibrate the norms of the
C     rows and columns of a matrix.  Roundoff errors are possible as
C     LAPACK routine DGEBAL does NOT use powers of the machine base for
C     scaling. The second stage (equilibrating ||G|| and ||Q||) requires
C     N**2 multiplications.
C     For norm scaling, 3*N**2 + O(N) multiplications are required and
C     NO rounding errors occur as all multiplications are performed with
C     powers of the machine base.
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C
C     REVISIONS
C
C     1997, December 5.
C     1998, August 25.
C
C     KEYWORDS
C
C     Hamiltonian matrix, symplectic similarity transformation,
C     scaling.
C
C    *************************************************************
*
      DOUBLE PRECISION A(LDA,N),D(*),QG(LDQG,N+1),RWORK(N)
C     ..
C     .. Parameters ..
*
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Local Scalars ..
*
      DOUBLE PRECISION ANRM,BASE,EPS,GNRM,OFL,QNRM,RHO,SFMAX,SFMIN,TAU,
     +                 UFL,Y
      INTEGER I,IHI,ILO,INFO,J
      LOGICAL NONE,NORM,SYMP
C     ..
C     .. External Functions ..
*     . BLAS, LAPACK .
*
      DOUBLE PRECISION DLAMCH,DLANGE,DLANSY
      LOGICAL LSAME
      EXTERNAL DLAMCH,DLANGE,DLANSY,LSAME
C     ..
C     .. External Subroutines ..
*     . BLAS, LAPACK .
*     . other .
*
      EXTERNAL DGEBAL,DLABAD,DLASCL,DRSCL
C     ..
C     .. Intrinsic Functions ..
*
      INTRINSIC ABS,MAX,SQRT
C     ..
*
      IERR = 0
      IF (N.NE.0) THEN
          SYMP = LSAME(JOBSCL,'A')
          NORM = LSAME(JOBSCL,'B')
          NONE = LSAME(JOBSCL,'N')
*
          IF (.NOT.SYMP .AND. .NOT.NORM .AND. .NOT.NONE) IERR = -1
          IF (N.LT.0) IERR = -2
          IF (LDA.LT.N) IERR = -4
          IF (LDQG.LT.N) IERR = -6
*
          IF (IERR.EQ.0) THEN
*
*            .. some machine dependant constants ..
              BASE = DLAMCH('B')
              EPS = DLAMCH('P')
              UFL = DLAMCH('S')
              OFL = ONE/UFL
              CALL DLABAD(UFL,OFL)
              SFMAX = (EPS/BASE)/UFL
              SFMIN = ONE/SFMAX

*            .. set tolerance TOL ..
              IF (NORM) THEN
                  ANRM = DLANGE('1',N,N,A,LDA,RWORK)
                  GNRM = DLANSY('1','U',N,QG(1,2),LDQG,RWORK)
                  QNRM = DLANSY('1','L',N,QG,LDQG,RWORK)
                  Y = MAX(ONE,ANRM,GNRM,QNRM)
                  TAU = ONE
*
   10             CONTINUE
*               .. while TAU < Y ..
                  IF ((TAU.LT.Y) .AND. (TAU.LT.SFMAX)) THEN
                      TAU = TAU*BASE
                      GO TO 10

                  END IF
*               .. end while ..
                  IF (TAU.GT.ONE) THEN
                      IF (ABS(TAU/BASE-Y).LT.ABS(TAU-Y)) TAU = TAU/BASE
                      CALL DLASCL('G',N,N,TAU,ONE,N,N,A,LDA,INFO)
                      CALL DLASCL('U',N,N,TAU,ONE,N,N,QG(1,2),LDQG,INFO)
                      CALL DLASCL('U',N,N,TAU,ONE,N,N,QG(1,2),LDQG,INFO)
                  END IF
*
                  D(1) = TAU
*
              ELSE IF (SYMP) THEN
                  CALL DGEBAL('S',N,A,LDA,ILO,IHI,D,INFO)
                  DO 30 J = 1,N
                      DO 20 I = J,N
                          QG(I,J) = QG(I,J)*D(J)*D(I)
                          QG(J,I+1) = QG(J,I+1)/D(J)/D(I)
   20                 CONTINUE
   30             CONTINUE
*
                  GNRM = DLANSY('1','U',N,QG(1,2),LDQG,RWORK)
                  QNRM = DLANSY('1','L',N,QG,LDQG,RWORK)
                  IF (GNRM.EQ.ZERO) THEN
                      IF (QNRM.EQ.ZERO) RHO = ONE
                      IF (QNRM.NE.ZERO) RHO = SFMAX

                  ELSE IF (QNRM.EQ.ZERO) THEN
                      IF (GNRM.EQ.ZERO) RHO = ONE
                      IF (GNRM.NE.ZERO) RHO = SFMIN

                  ELSE
                      RHO = SQRT(QNRM)/SQRT(GNRM)
                  END IF
*
                  CALL DLASCL('L',0,0,RHO,ONE,N,N,QG,LDQG,INFO)
                  CALL DLASCL('U',0,0,ONE,RHO,N,N,QG(1,2),LDQG,INFO)
                  CALL DRSCL(N,SQRT(RHO),D,1)
*
              END IF
*
          END IF
*
      END IF
*     *** Last Line of DHABL ***
      END
      SUBROUTINE DHAEVD(N,A,LDA,QG,LDQG,WR,WI,RWORK,IERR)
C
C     PURPOSE
C
C     This is an ``easy-to-use'' driver that uses the method of [1] to
C     calculate the eigenvalues of a Hamiltonian matrix.  A matrix H
C     is Hamiltonian if it has the form of
C
C                  ( A   G )                T        T
C     (*)    H  =  (      T),   where  G = G ,  Q = Q .
C                  ( Q  -A )
C
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q in (*).  Ordinarily,
C             N must be positive.  If N equals zero, then DHAEVD sets
C             IERR = 0 and returns immediately without accessing the
C             other arguments.
C
C      A      (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input, the leading N-by-N part of this array contains
C             the upper left block A of the Hamiltonian matrix H in (*).
C             On output A is overwritten by the leading N-by-N upper
C             left block of the square reduced Hamiltonian H' in (**).
C             (See Section METHOD below.)
C
C     LDA     (input) INTEGER
C             The leading dimension of array A as declared in the
C             calling program.  LDA .GE. N.
C
C     QG      (input/output) DOUBLE PRECISION array, 
C             dimension (LDQG,N+1)
C             On input QG contains the upper right symmetric block G and
C             the lower left symmetric block Q of the Hamiltonian H in
C             (*).
C             On output QG contains the upper right symmetric block G'
C             and lower left symmetric block Q' of the square reduced
C             Hamiltonian matrix H' (**) as produced by DHASRD; see
C             Section METHOD below.
C             The diagonal and lower triangle of the submatrix in
C             columns 1 to N contain the lower triangular part of Q.
C             The diagonal and upper triangle of the submatrix in
C             columns 2 to N+1 contain the upper triangle of G.
C             So, if I. GE. J, then Q(i,j) = Q(j,i) is stored in QG(i,j)
C             and G(i,j) = G(j,i) is in QG(j,i+1).
C
C     LDQG    (input) INTEGER
C             The leading dimension of array QG as declared in the
C             calling procedure.  LDQG .GE. N.
C
C     WR      (output) DOUBLE PRECISION array, dimension at least N
C     WI      (output) DOUBLE PRECISION array, dimension at least N.
C             WR and WI are overwritten by N eigenvalues of H with
C             non-negative real part.  The remaining N eigenvalues are
C             the negatives of these eigenvalues.  The real parts
C             overwrite WR.  The imaginary parts overwrite WI.
C             Eigenvalues are stored in WR and WI in decreasing order of
C             magnitude of the real parts, i.e., WR(I) .GE. WR(I+1).
C             (In particular, an eigenvalue closest to the imaginary
C              axis is WR(N)+WI(N)i.)
C             In addition, eigenvalues with zero real part are sorted in
C             decreasing order of magnitude of imaginary parts.  Note
C             that non-real eigenvalues with non-zero real part appear
C             in complex conjugate pairs, but eigenvalues with zero real
C             part do not, in general, appear in complex conjugate
C             pairs.
C
C     Workspace
C
C     RWORK   DOUBLE PRECISION array, dimension at least N*(N+1)
C
C     Error Indicator
C
C     IERR    INTEGER
C             = 0:  successful exit;
C             < 0:  if IERR = -j, then the j-th argument had an illegal
C                   value;
C             > 0:  if IERR =  j, then LAPACK subroutine DHSEQR failed
C                   to converge while computing the j-th eigenvalue.
C
C     METHOD
C
C     DHAEVD calls DHASRD to transform the Hamiltonian matrix H in (*)
C     into a square-reduced Hamiltonian matrix
C
C                       (A'  G' )
C     (**)         H' = (      T)
C                       (Q' -A' )
C                                                                 T
C     by an orthogonal symplectic similarity transformation H' = U H U
C     where
C                      (  U1   U2 )
C     (***)        U = (          ).
C                      ( -U2   U1 )
C                                                               T
C     The square reduced Hamiltonian matrix satisfies  Q'A' - A' Q' = 0,
C     and
C
C           2       T     2     ( A''   G'')
C         H'  :=  (U  H U)   =  (        T ).
C                               ( 0     A'')
C
C     In addition, A'' is upper Hessenberg and G'' is skew symmetric.
C     The square roots of the eigenvalues of  A'' = A'A' + G'Q'  are the
C     eigenvalues of H.
C     DHAEVD calls DHAEVS to form A'', to compute its eigenvalues, and
C     to recover the eigenvalues of H from those of A''.
C
C     REFERENCES
C
C     [1] Van Loan, C. F.
C         A symplectic method for approximating all the eigenvalues of
C         a Hamiltonian matrix.
C         Linear Algebra and its Applications 61 (1984), pp. 233-251.
C
C     [2] Byers, R.
C         Hamiltonian and symplectic algorithms for the algebraic
C         Riccati equation.
C         Ph. D. Thesis, Cornell University, Ithaca, NY, January 1983.
C
C     [3] Benner, P., Byers, R., and Barth, E.
C         Fortran 77 Subroutines for Computing the Eigenvalues of
C         Hamiltonian Matrices {I}: The Square-Reduced Method.
C         Submitted to ACM Trans. Math. Software, 1998.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires (86/3)*N**3 + O(N**2) floating point
C     operations.
C     Eigenvalues computed by this subroutine are exact eigenvalues
C     of a perturbed Hamiltonian matrix  H + E  where
C
C                 || E || <= c sqrt(eps) || H ||,
C
C     c is a modest constant depending on the dimension N and eps is the
C     machine precision. Moreover, if the norm of H and an eigenvalue
C     are of roughly the same magnitude, the computed eigenvalue is
C     essentially as accurate as the computed eigenvalue obtained by
C     traditional methods. See [1] or [2].
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C
C     REVISIONS
C
C     1997, December 5.
C     1998, August 25.
C
C     KEYWORDS
C
C     Eigenvalues, (square-reduced) Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER IERR,LDA,LDQG,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,N),QG(LDQG,N+1),RWORK(N* (N+1)),WI(N),WR(N)
C     ..
C     .. Local Scalars ..
C     ..
C     .. External Subroutines ..
      EXTERNAL DHAEVS,DHASRD
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(1)
C     ..
*
      IERR = 0
      IF (N.NE.0) THEN
*
          IF (N.LT.0) IERR = -1
          IF (LDA.LT.N) IERR = -3
          IF (LDQG.LT.N) IERR = -5
*
          IF (IERR.EQ.0) CALL DHASRD('N',N,A,LDA,QG,LDQG,DUMMY,1,RWORK,
     +                               IERR)
          IF (IERR.EQ.0) CALL DHAEVS('S',N,A,LDA,QG,LDQG,WR,WI,RWORK,
     +                               IERR)
*
      END IF
*
      RETURN
*     *** Last Line of DHAEVD ***
      END
      SUBROUTINE DHAEVS(JOBSCL,N,A,LDA,QG,LDQG,WR,WI,RWORK,IERR)
C
C     .. Scalar Arguments ..
      INTEGER IERR,LDA,LDQG,N
      CHARACTER JOBSCL
C     ..
C     .. Array Arguments ..
C
C     PURPOSE
C
C     This subroutine computes the eigenvalues of an N-by-N square-
C     reduced Hamiltonian matrix
C
C                  ( A   G  )
C     (*)    H  =  (      T )
C                  ( Q  -A  )
C
C     Here, A is an N-by-N matrix, and G and Q are symmetric N-by-N
C     matrices.  It is assumed without a check that H is square-reduced,
C     i.e., that
C
C             2    ( A'   G'  )
C     (**)   H  =  (        T )    with A' upper Hessenberg.
C                  ( 0    A'  )
C
C                            T                2
C     (Equivalently, Q'A'- A' Q' = 0, A'' = A' + G'Q', and for i > j+1,
C      A''(i,j) = 0.)  Ordinarily, H is the output from DHASRD.
C     The eigenvalues of H are computed as the square roots of the
C     eigenvalues of A''.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBSCL  CHARACTER*1
C             JOBSCL is the mode parameter passed to LAPACK subroutine
C             DGEBAL to optionally permute and balance the Hessenberg
C             matrix A' in (**).  See LAPACK subroutine DGEBAL and
C             Section METHOD below.
C             = 'N' or 'n':  do neither scaling nor balancing;
C             = 'P' or 'p':  attempt to permute the A' in (**) to block
C                            triangular form;
C             = 'S' or 's':  do scaling in order to equilibrate the rows
C                            and columns of A';
C             = 'B' or 'b':  perform both, 'P' and 'S'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  Ordinarily N
C             must be positive.  If N equals zero, then DHAEVS sets
C             IERR = 0 and returns immediately without accessing the
C             other arguments.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array contains the
C             upper left block A of the square-reduced Hamiltonian
C             matrix H in (*) as produced by DHASRD.
C
C     LDA     (input) INTEGER
C             The leading dimension of array A as declared in the
C             calling program.  LDA .GE. N.
C
C     QG      (input) DOUBLE PRECISION array, dimension (LDQG,N+1)
C             QG contains the upper right symmetric block G and the
C             lower left symmetric block Q of the square-reduced
C             Hamiltonian matrix H in (*) as produced by DHASRD.
C             The diagonal and lower triangle of the submatrix in
C             columns 1 to N contain the lower triangular part of Q.
C             The diagonal and upper triangle of the submatrix in
C             columns 2 to N+1 contain the upper triangle of G.
C             So, if I. GE. J, then Q(i,j) = Q(j,i) is stored in QG(i,j)
C             and G(i,j) = G(j,i) is in QG(j,i+1).
C
C     LDQG    (input) INTEGER
C             The leading dimension of array QG just as declared in the
C             calling procedure.  LDQG. GE. N.
C
C     WR      (output) DOUBLE PRECISION array, dimension at least N.
C     WI      (output) DOUBLE PRECISION array, dimension at least N.
C             WR and WI are overwritten by N eigenvalues of H with
C             non-negative real part.  The remaining N eigenvalues are
C             the negatives of these eigenvalues.  The real parts
C             overwrite WR.  The imaginary parts overwrite WI.
C             Eigenvalues are stored in WR and WI in decreasing order of
C             magnitude of the real parts, i.e., WR(I) .GE. WR(I+1).
C             (In particular, an eigenvalue closest to the imaginary
C              axis is WR(N)+WI(N)i.)
C             In addition, eigenvalues with zero real part are sorted in
C             decreasing order of magnitude of imaginary parts.  Note
C             that non-real eigenvalues with non-zero real part appear
C             in complex conjugate pairs, but eigenvalues with zero real
C             part do not, in general, appear in complex conjugate
C             pairs.
C
C     Workspace
C
C     RWORK   DOUBLE PRECISION array, dimension at least N*(N+1)
C
C     Error Indicator
C
C     IERR    INTEGER
C             = 0:  successful exit;
C             < 0:  if IERR = -j, then the j-th argument had an illegal
C                   value;
C             > 0:  if IERR =  j, then LAPACK subroutine DHSEQR failed
C                   to converge while computing the j-th eigenvalue.
C
C     METHOD
C
C     DHAEVS forms the upper Hessenberg matrix A' in (**) and calls
C     LAPACK subroutines to calculate its eigenvalues.  The eigenvalues
C     of H are the square roots of the eigenvalues of A'.
C
C     REFERENCES
C
C     [1] Van Loan, C. F.
C         A symplectic method for approximating all the eigenvalues of
C         a Hamiltonian matrix.
C         Linear Algebra and its Applications 61 (1984), pp. 233-251.
C
C     [2] Byers, R.
C         Hamiltonian and symplectic algorithms for the algebraic
C         Riccati equation.
C         Ph. D. Thesis, Cornell University, Ithaca, NY, January 1983.
C
C     [3] Benner, P., Byers, R., and Barth, E.
C         Fortran 77 Subroutines for Computing the Eigenvalues of
C         Hamiltonian Matrices {I}: The Square-Reduced Method.
C         Submitted to ACM Trans. Math. Software, 1998.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires (32/3)*N**3 + O(N**2) floating point
C     operations.
C     Eigenvalues computed by this subroutine are exact eigenvalues
C     of a perturbed Hamiltonian matrix  H + E  where
C
C                 || E || <= c sqrt(eps) || H ||,
C
C     c is a modest constant depending on the dimension N and eps is the
C     machine precision. Moreover, if the norm of H and an eigenvalue
C     are of roughly the same magnitude, the computed eigenvalue is
C     essentially as accurate as the computed eigenvalue obtained by
C     traditional methods. See [1] or [2].
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C
C     REVISIONS
C
C     1994, February 1.
C     1994, October 4.
C     1995, March 22.
C     1996, April 10.
C     1997, December 15.
C     1998, August 24.
C
C     KEYWORDS
C
C     Eigenvalues, (square-reduced) Hamiltonian matrix.
C
C     ******************************************************************
*
      DOUBLE PRECISION A(LDA,N),QG(LDQG,N+1),RWORK(N* (N+1)),WI(N),WR(N)
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SWAP,X,Y
      INTEGER I,IGNORE,IHI,ILO,J,M
      LOGICAL SORTED
C     ..
C     .. External Functions ..
*     . BLAS, LAPACK .
      LOGICAL LSAME
      EXTERNAL LSAME
C     ..
C     .. External Subroutines ..
*     . BLAS, LAPACK .
*     . others .
      EXTERNAL DCOPY,DCROOT,DGEBAL,DGEMM,DHSEQR,DLACPY,DSYMV
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(1)
C     ..
*
      IERR = 0
      IF (N.NE.0) THEN
          IF (.NOT. (LSAME(JOBSCL,'S').OR.LSAME(JOBSCL,
     +        'P').OR.LSAME(JOBSCL,'B').OR.LSAME(JOBSCL,'N'))) IERR = -1
          IF (N.LT.0) IERR = -2
          IF (LDA.LT.N) IERR = -4
          IF (LDQG.LT.N) IERR = -6
*
          IF (IERR.EQ.0) THEN
*             ..                               2
*             .. form the eigenvalues of A' = A  + GQ ..
              CALL DLACPY('L',N,N,QG,LDQG,RWORK,N)
              DO 10 I = 1,N
                  CALL DCOPY(N-I+1,RWORK(I+N* (I-1)),1,
     +                       RWORK(I+N* (I-1)),N)
   10         CONTINUE
              DO 20 J = 1,N
                  CALL DSYMV('U',N,ONE,QG(1,2),LDQG,RWORK(1+N* (J-1)),1,
     +                       ZERO,WR,1)
                  CALL DCOPY(N,WR,1,RWORK(1+N* (J-1)),1)
   20         CONTINUE
              CALL DGEMM('N','N',N,N,N,ONE,A,LDA,A,LDA,ONE,RWORK,N)
*             ..                          2
*             .. Find the eigenvalues of A  + GQ ..
              CALL DGEBAL(JOBSCL,N,RWORK,N,ILO,IHI,RWORK(1+N*N),IGNORE)
              CALL DHSEQR('E','N',N,ILO,IHI,RWORK,N,WR,WI,DUMMY,1,
     +                    RWORK(1+N*N),N,IERR)
              IF (IERR.EQ.0) THEN
*                 ..
*                 .. eigenvalues of H are the square roots ..
                  DO 30 I = 1,N
                      X = WR(I)
                      Y = WI(I)
                      CALL DCROOT(X,Y,WR(I),WI(I))
   30             CONTINUE
*                 ..
*                 .. Sort eigenvalues into decreasing order by real part
*                    and, for eigenvalues with zero real part only,
*                    decreasing order of imaginary part.  (This simple
*                    bubble sort preserves the relative order of
*                    eigenvalues with equal but nonzero real part.
*                    This ensures that complex conjugate pairs remain
*                    together.  ..
                  SORTED = .FALSE.
                  DO 50 M = N,1,-1
                      IF (SORTED) GO TO 60
                      SORTED = .TRUE.
                      DO 40 I = 1,M - 1
                          IF (((WR(I).LT.WR(I+1)).OR.
     +                        ((WR(I).EQ.ZERO).AND.
     +                        (WR(I+1).EQ.ZERO).AND.
     +                        (WI(I).LT.WI(I+1))))) THEN
                              SWAP = WR(I)
                              WR(I) = WR(I+1)
                              WR(I+1) = SWAP
                              SWAP = WI(I)
                              WI(I) = WI(I+1)
                              WI(I+1) = SWAP
*
                              SORTED = .FALSE.
*
                          END IF
*
   40                 CONTINUE
   50             CONTINUE
   60             CONTINUE
*
              END IF
*
          END IF
*
      END IF
*
      RETURN
*     *** Last Line of DHAEVS ***
      END
      SUBROUTINE DHASRD(COMPU,N,A,LDA,QG,LDQG,U,LDU,RWORK,IERR)
C
C     .. Scalar Arguments ..
      INTEGER IERR,LDA,LDQG,LDU,N
      CHARACTER COMPU
C     ..
C     .. Array Arguments ..
C
C
C     PURPOSE
C
C     To transform a Hamiltonian matrix
C
C                  ( A   G  )
C     (*)      H = (      T )
C                  ( Q  -A  )
C
C     into a square-reduced Hamiltonian matrix
C
C                   (A'  G' )
C     (**)     H' = (      T)
C                   (Q' -A' )
C                                                                 T
C     by an orthogonal symplectic similarity transformation H' = U H U
C     where
C                  (  U1   U2 )
C     (***)    U = (          ).
C                  ( -U2   U1 )
C                                                              T
C     The square-reduced Hamiltonian matrix satisfies Q'A' - A' Q' = 0,
C     and
C
C           2       T     2     ( A''   G'')
C         H'  :=  (U  H U)   =  (        T ).
C                               ( 0     A'')
C
C     In addition, A'' is upper Hessenberg and G'' is skew symmetric.
C     The square roots of the eigenvalues of A'' = A'*A'+G'*Q' are the
C     eigenvalues of H.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPU   CHARACTER*1
C             Indicates whether the orthogonal symplectic similarity
C             transformation matrix U in (***) is returned or
C             accumulated into an orthogonal symplectic matrix or if the
C             transformation matrix is not required.
C             = 'N' or 'n':             U is not referenced;
C             = 'V', 'v', 'A', or 'a':  the orthogonal symplectic
C                                       similarity transformations are
C                                       accumulated into U.  On input, U
C                                       must contain an orthogonal
C                                       symplectic matrix S; on exit , U
C                                       contains S*U with U from (**);
C             = 'I', 'i', 'F', or 'f':  on entry, U need not be set, and
C                                       on exit, U contains the
C                                       orthogonal symplectic matrix U
C                                       from (**).
C             See the description of U below for details.
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  Normally, N.GE.1.
C             If N = 0, then IERR is set to 0, DHASRD returns 
C             immediately without referencing any other argument.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input, the leading N-by-N part of this array contains
C             the upper left block A of the Hamiltonian matrix H in (*).
C             On output, the leading N-by-N part of this array is
C             overwritten by the upper right block A' of the
C             square-reduced Hamiltonian matrix H' in (**).
C
C     LDA     (input) INTEGER
C             The leading dimension of array A as declared in the
C             calling program. LDA .GE. N.
C
C     QG      (input/output) DOUBLE PRECISION array, 
C             dimension (LDQG,N+1)
C             On input QG contains the upper right symmetric block G and
C             the lower left symmetric block Q of the Hamiltonian matrix
C             in (*).  On output these are overwritten by the
C             corresponding blocks of the square reduced Hamiltonian
C             matrix H' in (**).
C             The diagonal and lower triangle of the submatrix in
C             columns 1 to N contain the lower triangular part of Q.
C             The diagonal and upper triangle of the submatrix in
C             columns 2 to N+1 contains the upper triangle of G.
C             So, if I. GE. J, then Q(i,j) = Q(j,i) is stored in QG(i,j)
C             and G(i,j) = G(j,i) is in QG(j,i+1).
C
C     LDQG    (input) INTEGER
C             The leading dimension of array QG just as declared in the
C             calling procedure.  LDQG. GE. N.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,2*N)
C             If COMPU = 'N' or 'n', then U is not referenced.
C             If COMPU = 'I', 'i', 'F' or 'f', then the input contents
C             of U are not specified.  On output the leading (N-by-2*N)
C             segment of array U is overwritten by the the first N rows
C             of the orthogonal symplectic matrix U in (**).
C             If COMPU = 'V', 'v', 'A' or 'a', then, on input, the
C             leading N-by-(2*N) segment of this array is assumed to
C             contain the first N rows of an orthogonal symplectic
C             matrix S. On output this array is overwritten by the first
C             N rows of the product S*U where U is the orthogonal
C             symplectic matrix from (**).
C             The storage scheme implied by (***) is used for orthogonal
C             symplectic matrices, i.e., only the first N rows are
C             stored as they contain all relevant information.
C
C     LDU     (input) INTEGER
C             The leading dimension of array U as declared in the
C             calling program.  If COMPU is not 'N' or 'n', then
C             LDU .GE. N.
C             If COMPU = 'N' or 'n', then the array U is not referenced
C             and LDU may be set to 1.
C
C     Workspace
C
C     RWORK   DOUBLE PRECISION array, dimension at least 2*N
C
C     Error Indicator
C
C     IERR    INTEGER
C             = 0:  successful exit;
C             < 0:  if IERR = -j, then the j-th argument had an illegal
C                   value.
C
C     METHOD
C
C     DHASRD transforms a Hamiltonian matrix H into a square-reduced
C     Hamiltonian matrix H' using the implicit version of Van Loan's
C     method as proposed in [1,2,3].
C
C     REFERENCES
C
C     [1] Van Loan, C. F.
C         A symplectic method for approximating all the eigenvalues of
C         a Hamiltonian matrix.
C         Linear Algebra and its Applications 61 (1984), pp. 233-251.
C
C     [2] Byers, R.
C         Hamiltonian and symplectic algorithms for the algebraic
C         Riccati equation.
C         Ph. D. Thesis, Cornell University, Ithaca, NY, January 1983.
C
C     [3] Benner, P., Byers, R., and Barth, E.
C         Fortran 77 Subroutines for Computing the Eigenvalues of
C         Hamiltonian Matrices {I}: The Square-Reduced Method.
C         Submitted to ACM Trans. Math. Software, 1998.
C
C     NUMERICAL ASPECTS
C
C     This algorithm requires approximately 20*N**3 flops for
C     transforming H into square-reduced form. If the transformations
C     are required, this adds another 8*N**3 flops. The method is
C     strongly backward stable in the sense that if H' and U are the
C     computed square-reduced Hamiltonian and computed orthogonal
C     symplectic similarity transformation, then there is an orthogonal
C     symplectic matrix T and a Hamiltonian matrix M such that
C
C                  H T  =  T M
C
C        || T - U ||   <=  c1 * eps
C
C        || H' - M ||  <=  c2 * eps * || H ||
C
C     where c1, c2 are modest constants depending on the dimension N and
C     eps is the machine precision.
C
C     Eigenvalues computed by explicitly forming the upper Hessenberg
C     matrix  W = A'A' + G'Q' with A', G', Q' as in (**) and applying
C     the Hessenberg QR iteration to W are exactly eigenvalues of a
C     perturbed Hamiltonian matrix H + E,  where
C
C        || E ||  <=  c3 * sqrt(eps) * || H ||,
C
C     and c3 is a modest constant depending on the dimension N and eps
C     is the machine precision.  Moreover, if the norm of H and an
C     eigenvalue lambda are of roughly the same magnitude, the computed
C     eigenvalue is essentially as accurate as the computed eigenvalue
C     from traditional methods.  See [1] or [2].
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C
C     REVISIONS
C
C     1982, July.
C     1992, June.
C     1995, March.
C     1996, April.
C     1997, December.
C     1998, April.
C     1998, August 26.
C
C     KEYWORDS
C
C     (Square-reduced) Hamiltonian matrix, symplectic orthogonal
C     similarity transformation.
C
C     ******************************************************************
C
      DOUBLE PRECISION A(LDA,N),QG(LDQG,N+1),RWORK(2*N),U(LDU,*)
C     ..
C     .. Parameters ..
C
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION COSINE,SINE,TAU,TEMP,X,Y
      INTEGER I,J
      LOGICAL ACCUM,FORGET,FORM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION T(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      LOGICAL LSAME
      EXTERNAL DDOT,LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMV,DLARFG,DLARFX,DLARTG,DROT,DSYMV,DSYR2
C     ..
C
      IERR = 0
      IF (N.NE.0) THEN
          ACCUM = LSAME(COMPU,'A') .OR. LSAME(COMPU,'V')
          FORM = LSAME(COMPU,'F') .OR. LSAME(COMPU,'I')
          FORGET = LSAME(COMPU,'N')

          IF (.NOT.ACCUM .AND. .NOT.FORM .AND. .NOT.FORGET) IERR = -1
          IF (N.LT.0) IERR = -2
          IF (LDA.LT.N) IERR = -4
          IF (LDQG.LT.N) IERR = -6
          IF (LDU.LT.N .AND. .NOT.FORGET) IERR = -8
          IF (IERR.EQ.0) THEN
C
C             .. TRANSFORM TO SQUARE REDUCED FORM ..
              DO 10 J = 1,N - 1
C                                     T
C                 .. RWORK <- (Q A - A Q)(J+1:N,J) ..
                  CALL DCOPY(J-1,QG(J,1),LDQG,RWORK(N+1),1)
                  CALL DCOPY(N-J+1,QG(J,J),1,RWORK(N+J),1)
                  CALL DGEMV('T',N,N-J,-ONE,A(1,J+1),LDA,RWORK(N+1),1,
     +                       ZERO,RWORK(J+1),1)
                  CALL DGEMV('N',N-J,J,ONE,QG(J+1,1),LDQG,A(1,J),1,ONE,
     +                       RWORK(J+1),1)
                  CALL DSYMV('L',N-J,ONE,QG(J+1,J+1),LDQG,A(J+1,J),1,
     +                       ONE,RWORK(J+1),1)
C
C                 .. SYMPLECTIC REFLECTION TO ZERO (H*H)((N+J+2):2N,J)
                  CALL DLARFG(N-J,RWORK(J+1),RWORK(J+1+1),1,TAU)
                  Y = RWORK(J+1)
                  RWORK(J+1) = ONE
C
                  CALL DLARFX('L',N-J,N,RWORK(J+1),TAU,A(J+1,1),LDA,
     +                        RWORK(N+1))
                  CALL DLARFX('R',N,N-J,RWORK(J+1),TAU,A(1,J+1),LDA,
     +                        RWORK(N+1))
C
                  CALL DLARFX('L',N-J,J,RWORK(J+1),TAU,QG(J+1,1),LDQG,
     +                        RWORK(N+1))
                  CALL DSYMV('L',N-J,TAU,QG(J+1,J+1),LDQG,RWORK(J+1),1,
     +                       ZERO,RWORK(N+J+1),1)
                  CALL DAXPY(N-J,-TAU*DDOT(N-J,RWORK(N+J+1),1,
     +                       RWORK(J+1),1)/TWO,RWORK(J+1),1,
     +                       RWORK(N+J+1),1)
                  CALL DSYR2('L',N-J,-ONE,RWORK(J+1),1,RWORK(N+J+1),1,
     +                       QG(J+1,J+1),LDQG)
C
                  CALL DLARFX('R',J,N-J,RWORK(J+1),TAU,QG(1,J+2),LDQG,
     +                        RWORK(N+1))
                  CALL DSYMV('U',N-J,TAU,QG(J+1,J+2),LDQG,RWORK(J+1),1,
     +                       ZERO,RWORK(N+J+1),1)
                  CALL DAXPY(N-J,-TAU*DDOT(N-J,RWORK(N+J+1),1,
     +                       RWORK(J+1),1)/TWO,RWORK(J+1),1,
     +                       RWORK(N+J+1),1)
                  CALL DSYR2('U',N-J,-ONE,RWORK(J+1),1,RWORK(N+J+1),1,
     +                       QG(J+1,J+2),LDQG)
C
                  IF (FORM) THEN
C                     .. SAVE REFLECTION ..
                      CALL DCOPY(N-J,RWORK(J+1),1,U(J+1,J),1)
                      U(J+1,J) = TAU
C
                  ELSE IF (ACCUM) THEN
C                     .. ACCUMULATE REFLECTION ..
                      CALL DLARFX('R',N,N-J,RWORK(J+1),TAU,U(1,J+1),LDU,
     +                            RWORK(N+1))
                      CALL DLARFX('R',N,N-J,RWORK(J+1),TAU,U(1,N+J+1),
     +                            LDU,RWORK(N+1))
                  END IF
C
C                 .. (X,Y) := ((J+1,J),(N+J+1,J)) COMPONENT OF H*H ..
                  X = DDOT(J,QG(1,J+2),1,QG(J,1),LDQG) +
     +                DDOT(N-J,QG(J+1,J+2),LDQG,QG(J+1,J),1) +
     +                DDOT(N,A(J+1,1),LDA,A(1,J),1)
C
C                 .. SYMPLECTIC ROTATION TO ZERO (H*H)(N+J+1,J) ..
                  CALL DLARTG(X,Y,COSINE,SINE,TEMP)
C
                  CALL DROT(J,A(J+1,1),LDA,QG(J+1,1),LDQG,COSINE,SINE)
                  CALL DROT(J,A(1,J+1),1,QG(1,J+2),1,COSINE,SINE)
                  IF (N-J-1.GT.0) THEN
                      CALL DROT(N-J-1,A(J+1,J+2),LDA,QG(J+2,J+1),1,
     +                          COSINE,SINE)
                      CALL DROT(N-J-1,A(J+2,J+1),1,QG(J+1,J+3),LDQG,
     +                          COSINE,SINE)
                  END IF
C
                  T(1,1) = A(J+1,J+1)
                  T(1,2) = QG(J+1,J+2)
                  T(2,1) = QG(J+1,J+1)
                  T(2,2) = -T(1,1)
                  CALL DROT(2,T(1,1),1,T(1,2),1,COSINE,SINE)
                  CALL DROT(2,T(1,1),2,T(2,1),2,COSINE,SINE)
                  A(J+1,J+1) = T(1,1)
                  QG(J+1,J+2) = T(1,2)
                  QG(J+1,J+1) = T(2,1)
C
                  IF (FORM) THEN
C                     .. SAVE ROTATION ..
                      U(J,J) = COSINE
                      U(J,N+J) = SINE
C
                  ELSE IF (ACCUM) THEN
C                     .. ACCUMULATE ROTATION ..
                      CALL DROT(N,U(1,J+1),1,U(1,N+J+1),1,COSINE,SINE)
                  END IF
C
C                 .. RWORK := (A  + GQ)(J+1:N,J) ..
                  CALL DGEMV('N',N-J,N,ONE,A(J+1,1),LDA,A(1,J),1,ZERO,
     +                       RWORK(J+1),1)
                  CALL DCOPY(J-1,QG(J,1),LDQG,RWORK(N+1),1)
                  CALL DCOPY(N-J+1,QG(J,J),1,RWORK(N+J),1)
                  CALL DGEMV('T',J,N-J,ONE,QG(1,J+2),LDQG,RWORK(N+1),1,
     +                       ONE,RWORK(J+1),1)
                  CALL DSYMV('U',N-J,ONE,QG(J+1,J+2),LDQG,RWORK(N+J+1),
     +                       1,ONE,RWORK(J+1),1)
C
C                 .. SYMPLECTIC REFLECTION TO ZERO (H*H)(J+2:N,J)
                  CALL DLARFG(N-J,RWORK(J+1),RWORK(J+1+1),1,TAU)
                  RWORK(J+1) = ONE
                  CALL DLARFX('L',N-J,N,RWORK(J+1),TAU,A(J+1,1),LDA,
     +                        RWORK(N+1))
                  CALL DLARFX('R',N,N-J,RWORK(J+1),TAU,A(1,J+1),LDA,
     +                        RWORK(N+1))
C
                  CALL DLARFX('L',N-J,J,RWORK(J+1),TAU,QG(J+1,1),LDQG,
     +                        RWORK(N+1))
                  CALL DSYMV('L',N-J,TAU,QG(J+1,J+1),LDQG,RWORK(J+1),1,
     +                       ZERO,RWORK(N+J+1),1)
                  CALL DAXPY(N-J,-TAU*DDOT(N-J,RWORK(N+J+1),1,
     +                       RWORK(J+1),1)/TWO,RWORK(J+1),1,
     +                       RWORK(N+J+1),1)
                  CALL DSYR2('L',N-J,-ONE,RWORK(J+1),1,RWORK(N+J+1),1,
     +                       QG(J+1,J+1),LDQG)
C
                  CALL DLARFX('R',J,N-J,RWORK(J+1),TAU,QG(1,J+2),LDQG,
     +                        RWORK(N+1))
                  CALL DSYMV('U',N-J,TAU,QG(J+1,J+2),LDQG,RWORK(J+1),1,
     +                       ZERO,RWORK(N+J+1),1)
                  CALL DAXPY(N-J,-TAU*DDOT(N-J,RWORK(N+J+1),1,
     +                       RWORK(J+1),1)/TWO,RWORK(J+1),1,
     +                       RWORK(N+J+1),1)
                  CALL DSYR2('U',N-J,-ONE,RWORK(J+1),1,RWORK(N+J+1),1,
     +                       QG(J+1,J+2),LDQG)
C
                  IF (FORM) THEN
C                     .. SAVE REFLECTION ..
                      CALL DCOPY(N-J,RWORK(J+1),1,U(J+1,N+J),1)
                      U(J+1,N+J) = TAU
C
                  ELSE IF (ACCUM) THEN
C                     .. ACCUMULATE REFLECTION ..
                      CALL DLARFX('R',N,N-J,RWORK(J+1),TAU,U(1,J+1),LDU,
     +                            RWORK(N+1))
                      CALL DLARFX('R',N,N-J,RWORK(J+1),TAU,U(1,N+J+1),
     +                            LDU,RWORK(N+1))
                  END IF
C
   10         CONTINUE
C
C
              IF (FORM) THEN
C                 .. FORM S BY ACCUMULATING TRANSFORMATIONS ..
                  DO 40 J = N - 1,1,-1
C                     .. INITIALIZE (J+1)ST COLUMN OF S ..
                      DO 20 I = 1,N
                          U(I,J+1) = ZERO
   20                 CONTINUE
                      U(J+1,J+1) = ONE
                      DO 30 I = 1,N
                          U(I,N+J+1) = ZERO
   30                 CONTINUE
C
C                     .. SECOND REFLECTION ..
                      TAU = U(J+1,N+J)
                      U(J+1,N+J) = ONE
                      CALL DLARFX('L',N-J,N-J,U(J+1,N+J),TAU,U(J+1,J+1),
     +                            LDU,RWORK(N+1))
                      CALL DLARFX('L',N-J,N-J,U(J+1,N+J),TAU,
     +                            U(J+1,N+J+1),LDU,RWORK(N+1))
C                     .. ROTATION ..
                      CALL DROT(N-J,U(J+1,J+1),LDU,U(J+1,N+J+1),LDU,
     +                          U(J,J),U(J,N+J))
C                     .. FIRST REFLECTION ..
                      TAU = U(J+1,J)
                      U(J+1,J) = ONE
                      CALL DLARFX('L',N-J,N-J,U(J+1,J),TAU,U(J+1,J+1),
     +                            LDU,RWORK(N+1))
                      CALL DLARFX('L',N-J,N-J,U(J+1,J),TAU,U(J+1,N+J+1),
     +                            LDU,RWORK(N+1))
   40             CONTINUE
C                 .. FIRST COLUMN IS FIRST COLUMN OF IDENTITY ..
                  DO 50 I = 1,N
                      U(I,1) = ZERO
                      U(I,N+1) = ZERO
   50             CONTINUE
                  U(1,1) = ONE
              END IF
C
          END IF
C
      END IF
C     .. LAST LINE OF DHASRD ..
      END
