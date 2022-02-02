*       test driver for DHASRD.f
*     *****************************************************************
*     This program sets up several ad-hoc Hamiltonian matrices and
*     transforms it to square-reduced form by Van Loan's method.
*
*     *****************************************************************
*     There are three consistency checks.
*     Test 1: Similarity Check
*             Check that the output Hamiltonian matrix is similar to the
*             input  Hamiltonian matrix ``modulo rounding errors.''
*     Test 2: Symplectic Orthogonality
*             Check that the similarity transformation is symplectic
*             orthogonal ``modulo rounding errors.''
*     Test 3: Square Reduced Check
*             Check that the output Hamiltonian matrix is square reduced
*             ``modulo rounding errors.''
*
*     ******************************************************************
C     .. Parameters ..
*        NOUT is the standard output, e.g., a terminal
      INTEGER NOUT
      PARAMETER (NOUT=6)
      INTEGER NMAX,NMIN
      PARAMETER (NMAX=15,NMIN=-1)
      INTEGER LD
      PARAMETER (LD=NMAX)
      INTEGER LDA,LDQG,LDU
      PARAMETER (LDA=LD,LDQG=LD,LDU=LD)
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      CHARACTER*(*) EQULS
      PARAMETER (EQULS=' =============================================')
C     ..
C     .. Local Scalars ..
*
      DOUBLE PRECISION ANRM,EPS,GNRM,HNRM,QNRM,SIMERR,SMPERR,SQRERR,Z11,
     +                 Z12,Z21
      INTEGER I,IERR,IFUDGE,ITER,J,N,NFAILS
      CHARACTER COMPU
C     ..
C     .. Local Arrays ..
*
      DOUBLE PRECISION A(LDA,LDA),AP(LDA,LDA),QG(LDQG,LDQG+1),
     +                 QPGP(LDQG,LDQG+1),RWORK(2*LD),U(LDU,2*LDU),
     +                 Y(LD,LD),Z(LD,LD)
      INTEGER ISEED(4)
C     ..
C     .. External Functions ..
*
      DOUBLE PRECISION DLAMCH,DLANGE,DLANSY,DLAPY2,DNRM2
      LOGICAL LSAME
      EXTERNAL DLAMCH,DLANGE,DLANSY,DLAPY2,DNRM2,LSAME
C     ..
C     .. External Subroutines ..
*
      EXTERNAL DCOPY,DGEMM,DHASRD,DLACPY,DLARNV,DLASET,DSYMM,DSYRK
C     ..
C     .. Intrinsic Functions ..
*
      INTRINSIC DBLE,SQRT
C     ..
C     .. Data statements ..
*     . (seed for random number generation) .
*
*     ******************************************************************
      DATA ISEED/86,1967,2001,1995/
C     ..
      EPS = DLAMCH('P')
*
      NFAILS = 0
      DO 50 ITER = 1,2
          IF (ITER.EQ.1) COMPU = 'F'
          IF (ITER.EQ.2) COMPU = 'A'
*         .. 2 test runs, 1 for generating transformation matrix, 1 for
*         accumulating transformations ..
          WRITE (NOUT,FMT=9020) EQULS
          WRITE (NOUT,FMT=9020) EQULS
          IF (LSAME(COMPU,'F')) THEN
              WRITE (NOUT,FMT=9020)
     +          'TEST RUN 1: GENERATE SYMPLECTIC ORTHOGONAL'
              WRITE (NOUT,FMT=9020)
     +          '            TRANSFORMATION MATRIX.'

          ELSE
              WRITE (NOUT,FMT=9020)
     +          'TEST RUN 2: ACCUMULATE SYMPLECTIC ORTHOGONAL '
              WRITE (NOUT,FMT=9020) '            TRANSFORMATIONS.'
          END IF

          DO 40 N = NMIN,NMAX
              WRITE (NOUT,FMT=9020)
              WRITE (NOUT,FMT=9020) EQULS
              WRITE (NOUT,FMT=9030) 'DHASRD TESTS: ORDER N = ',N
              IF (N.LT.0) THEN
                  WRITE (NOUT,FMT=9020)
     +              'OF COURSE, THIS IS NOT A VALID TEST CASE.'
                  WRITE (NOUT,FMT=9020)
     +             'IT CHECKS THAT DHASRD RETURNS WITH ERROR IERR = -2.'
              END IF
*
*             .. set up Hamiltonian ..
              DO 10 I = 1,N
                  CALL DLARNV(3,ISEED,N,A(1,I))
                  CALL DLARNV(3,ISEED,N,QG(1,I))
   10         CONTINUE
              IF(N.GE.0) THEN
                CALL DLARNV(3,ISEED,N,QG(1,N+1))
              ENDIF

              CALL DLACPY('A',N,N,A,LD,AP,LD)
              CALL DLACPY('A',N,N+1,QG,LDQG,QPGP,LDQG)

*
*             .. set up symplectic orthogonal matrix for accumulating
*             transformations ..
              IF (COMPU.EQ.'A') CALL DLASET('A',N,2*N,ZERO,ONE,U,LDU)
*
*             .. Frobenius norm of Hamiltonian Matrix ..
              ANRM = DLANGE('F',N,N,A,LD,RWORK)
              GNRM = DLANSY('F','L',N,QG,LDQG,RWORK)
              QNRM = DLANSY('F','U',N,QG(1,2),LDQG,RWORK)
              HNRM = DLAPY2(DLAPY2(ANRM,GNRM),DLAPY2(ANRM,QNRM))
*
*             .. square reduce ..
              CALL DHASRD(COMPU,N,A,LDA,QG,LDQG,U,LDU,RWORK,IERR)

*     ******************************************************************
*             .. Test 0: Has DHASRD something to tell? =====
              WRITE (NOUT,FMT=9010) 'IERR = ',IERR
              IF (IERR.EQ.-1) WRITE (NOUT,FMT=9020) 'COMPU = '//COMPU//
     +            ' IS NOT ONE OF ''A'', ''F'', OR ''N''.'
              IF (IERR.EQ.-2) WRITE (NOUT,FMT=9040) 'N = ',N,
     +            ' IS NEGATIVE.'
              IF (IERR.EQ.-4) WRITE (NOUT,FMT=9050) 'LDA = ',LDA,
     +            ' IS LESS THAN N = ',N,'.'
              IF (IERR.EQ.-6) WRITE (NOUT,FMT=9050) 'LDQG = ',LDQG,
     +            ' IS LESS THAN N = ',N,'.'
              IF (IERR.EQ.-8) WRITE (NOUT,FMT=9050) 'LDU = ',LDU,
     +            ' IS LESS THAN N = ',N,'.'
              IF (IERR.EQ.0) THEN
                  WRITE (NOUT,FMT=9020) 'NO ERRORS WERE DETECTED.'
*     ******************************************************************
*                 .. Test 1:  similarity check (modulo rounding) ..
*
*
*                 .. Z11 <- || (AP*S1 - GG*S2) - (S1*A + S2*Q)|| ..
                  CALL DGEMM('N','N',N,N,N,ONE,AP,LDA,U,LDU,ZERO,Z,LD)
                  CALL DSYMM('L','U',N,N,-ONE,QPGP(1,2),LDQG,U(1,N+1),
     +                       LDU,ONE,Z,LD)
                  CALL DGEMM('N','N',N,N,N,-ONE,U,LDU,A,LDA,ONE,Z,LD)
                  CALL DSYMM('R','L',N,N,-ONE,QG,LDQG,U(1,N+1),LDU,ONE,
     +                       Z,LD)
                  Z11 = DLANGE('F',N,N,Z,LD,RWORK)
*
*                 .. Z12 <- ||(AP*S2 + GG*S1) - (S1*G - S2*A') || ..
                  CALL DGEMM('N','N',N,N,N,ONE,AP,LDA,U(1,N+1),LDU,ZERO,
     +                       Z,LD)
                  CALL DSYMM('L','U',N,N,ONE,QPGP(1,2),LDQG,U,LDU,ONE,Z,
     +                       LD)
                  CALL DSYMM('R','U',N,N,-ONE,QG(1,2),LDQG,U,LDU,ONE,Z,
     +                       LD)
                  CALL DGEMM('N','T',N,N,N,ONE,U(1,N+1),LDU,A,LDA,ONE,Z,
     +                       LD)
                  Z12 = DLANGE('F',N,N,Z,LD,RWORK)
*
*                 .. Z21 <- || (QQ*S1 + AP'*S2) - (-S2*A + S1*Q) || ..
                  CALL DSYMM('L','L',N,N,ONE,QPGP,LDQG,U,LDU,ZERO,Z,LD)
                  CALL DGEMM('T','N',N,N,N,ONE,AP,LDA,U(1,N+1),LDU,ONE,
     +                       Z,LD)
                  CALL DGEMM('N','N',N,N,N,ONE,U(1,N+1),LDU,A,LDA,ONE,Z,
     +                       LD)
                  CALL DSYMM('R','L',N,N,-ONE,QG,LDQG,U,LDU,ONE,Z,LD)
                  Z21 = DLANGE('F',N,N,Z,LD,RWORK)
*
*                 .. || HH*S - S*H || / || H || ..
                  SIMERR = DLAPY2(DLAPY2(Z11,Z12),DLAPY2(Z11,Z21))
                  IF (HNRM.NE.ZERO) SIMERR = SIMERR/HNRM
                  IFUDGE = 10*N
*
                  WRITE (NOUT,FMT=9020)
                  WRITE (NOUT,FMT=9020) 'SIMILARITY CHECK:'
                  WRITE (NOUT,FMT=9000) 'RELATIVE ERROR  = ',SIMERR
*
                  IF (SIMERR.GT.DBLE(IFUDGE)*EPS) THEN
                      NFAILS = NFAILS + 1
                      WRITE (NOUT,FMT=9020)
     +                  'SUSPICIOUS: THE RELATIVE ERROR OF SIMILARITY'
                      WRITE (NOUT,FMT=9020)
     +                  '            IS LARGER THAN EXPECTED. '

                  ELSE
                      WRITE (NOUT,FMT=9020)
     +                  'OK: THE OUTPUT APPEARS TO BE SIMILAR TO THE'
                      WRITE (NOUT,FMT=9020)
     +                  '    INPUT ''MODULO ROUNDING ERRORS'' AS IT '
                      WRITE (NOUT,FMT=9020) '    SHOULD BE.'
                  END IF
*     ******************************************************************
*                 .. Test 2: symplectic orthogonal check ..
*
*                 .. Z11 <- || S2'*S1 - S1'*S2 || ..
                  CALL DGEMM('T','N',N,N,N,ONE,U(1,N+1),LDU,U,LDU,ZERO,
     +                       Z,LD)
                  CALL DGEMM('T','N',N,N,N,-ONE,U,LDU,U(1,N+1),LDU,ONE,
     +                       Z,LD)
                  Z11 = DLANGE('F',N,N,Z,LD,RWORK)
*
*                 .. Z12 <- || S2'S2 + S1'S1 - I || ===
                  CALL DLASET('U',N,N,ZERO,-ONE,Z,LD)
                  CALL DSYRK('U','T',N,N,ONE,U(1,N+1),LDU,ONE,Z,LD)
                  CALL DSYRK('U','T',N,N,ONE,U,LDU,ONE,Z,LD)
                  Z12 = DLANSY('F','U',N,Z,LD,RWORK)
*
                  SMPERR = SQRT(TWO)*DLAPY2(Z11,Z12)
                  IF (N.GT.0) SMPERR = SMPERR/SQRT(DBLE(N))
                  IFUDGE = 10*N
*
                  WRITE (NOUT,FMT=9020)
                  WRITE (NOUT,FMT=9020)
     +              'SYMPLECTIC ORTHOGONALITY CHECK:'
                  WRITE (NOUT,FMT=9000) 'RELATIVE ERROR  = ',SMPERR
                  IF (SMPERR.GT.DBLE(IFUDGE)*EPS) THEN
                      NFAILS = NFAILS + 1
                      WRITE (NOUT,FMT=9020)
     +                  'SUSPICIOUS: THE RELATIVE DISTANCE TO '
                      WRITE (NOUT,FMT=9020)
     +                  '            THE SYMPLECTIC ORTHOGONAL CLASS IS'
                      WRITE (NOUT,FMT=9020)
     +                  '            LARGER THAN EXPECTED.'

                  ELSE
                      WRITE (NOUT,FMT=9020)
     +                  'OK: THE OUTPUT SIMILARITY TRANSFORMATION IS '
                      WRITE (NOUT,FMT=9020)
     +                  '    SYMPLECTIC ORTHOGONAL ''MODULO ROUNDING '
                      WRITE (NOUT,FMT=9020)
     +                  '    ERRORS'' AS IT SHOULD BE.'
                  END IF
*     ******************************************************************
*                 .. Test 3: square reduced check ..
*
*                 .. Z <-  A*A + G*Q  ..
                  CALL DLACPY('L',N,N,QG,LDQG,Y,LD)
                  DO 20 I = 1,N - 1
                      CALL DCOPY(N-I,QG(I+1,I),1,Y(I,I+1),LD)
   20             CONTINUE
                  CALL DSYMM('L','U',N,N,ONE,QG(1,2),LDQG,Y,LD,ZERO,Z,
     +                       LD)
                  CALL DGEMM('N','N',N,N,N,ONE,A,LDA,A,LDA,ONE,Z,LD)
*
*                 .. Z11 <- FROB NORM BELOW SUB DIAGONAL ..
                  Z11 = ZERO
                  DO 30 J = 1,N - 2
                      Z11 = DLAPY2(Z11,DNRM2(N-J-1,Z(J+2,J),1))
   30             CONTINUE
*
*                 .. R21 <- || Q*A - A'*Q || ..
                  CALL DGEMM('T','N',N,N,N,-ONE,A,LDA,Y,LD,ZERO,Z,LD)
                  CALL DSYMM('L','L',N,N,ONE,QG,LDQG,A,LDA,ONE,Z,LD)
                  Z21 = DLANGE('F',N,N,Z,LD,RWORK)
*
*                 .. square reduced check ..
                  SQRERR = DLAPY2(Z11,Z21)
                  IF (HNRM.NE.ZERO) SQRERR = SQRERR/HNRM/HNRM
                  IFUDGE = 10*N
*
                  WRITE (NOUT,FMT=9020)
                  WRITE (NOUT,FMT=9020) 'SQUARE REDUCED CHECK:'
                  WRITE (NOUT,FMT=9000) 'RELATIVE ERROR  = ',SQRERR
                  IF (SQRERR.GT.DBLE(IFUDGE)*EPS) THEN
                      NFAILS = NFAILS + 1
                      WRITE (NOUT,FMT=9020)
     +                  'SUSPICIOUS: THE RELATIVE DIFFERENCE FROM'
                      WRITE (NOUT,FMT=9020)
     +                  '            BEING SQUARE REDUCED IS '
                      WRITE (NOUT,FMT=9020)
     +                  '            LARGER THAN EXPECTED. '

                  ELSE
                      WRITE (NOUT,FMT=9020)
     +                  'OK: THE OUTPUT IS SQUARE REDUCED ''MODULO '
                      WRITE (NOUT,FMT=9020)
     +                  '    ROUNDING ERRORS'' AS IT SHOULD BE.'
                  END IF

              END IF

              WRITE (NOUT,FMT=9020)
   40     CONTINUE
*         .. repeat test for accumulated transformations ..
   50 CONTINUE
*     ******************************************************************
      WRITE (NOUT,FMT=9020)
      WRITE (NOUT,FMT=9020)
      WRITE (NOUT,FMT=9020) EQULS
      IF (NFAILS.EQ.0) THEN
          WRITE (NOUT,FMT=9020) 'NO SUSPICIOUS TEST RESULTS.'

      ELSE
          WRITE (NOUT,FMT=9040) 'SIGH..., THERE WERE ',NFAILS,
     +      ' SUSPICIOUS TEST RESULTS.'
      END IF

      STOP
*

 9000 FORMAT (A,E10.4)
 9010 FORMAT (A,I5)
 9020 FORMAT (A)
 9030 FORMAT (A,I3)
 9040 FORMAT (A,I3,A)
 9050 FORMAT (A,I3,A,I3,A)
      END
