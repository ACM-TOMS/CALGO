***********************************************************************
*                                                                     *
*     Test program for ScaLAPACK-style routine PDHSEQR                *
*     (random matrix)                                                 *
*                                                                     *
*     Authors: Robert Granat                                          *
*              Bo Kagstrom                                            *
*              Daniel Kressner                                        *
*              Meiyue Shao                                            *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*                                                                     *
***********************************************************************
*
      PROGRAM TESTRANDOM
*
*     Declarations
*
      IMPLICIT NONE
*     ...Parameters...
      LOGICAL           COMPRESI, COMPORTH, TEST_CHKRESI
      LOGICAL           DEBUG, PRN, TIMESTEPS, UNI_LAPACK
      PARAMETER         ( DEBUG = .FALSE.,
     $                    PRN = .FALSE.,
     $                    TIMESTEPS = .TRUE.,
     $                    COMPRESI = .TRUE.,
     $                    COMPORTH = .TRUE.,
     $                    TEST_CHKRESI = .TRUE.,
     $                    UNI_LAPACK = .TRUE. )
      INTEGER           COMPHESS, BALANCE, N, NB, ARSRC, ACSRC
*     What processor should hold the first element in A?
      PARAMETER         ( ARSRC = 0, ACSRC = 0 )
      INTEGER           BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                  LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER         ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                    CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                    RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER*8         DPALLOC, INTALLC
      INTEGER           DPSIZ, INTSZ, NOUT, IZERO
      PARAMETER         ( DPSIZ = 8, DPALLOC = 400 000 000,
     $                    INTSZ = 8, INTALLC = 1 000 000,
     $                    NOUT = 6, IZERO = 0 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+00, ONE = 1.0D+00 )
*
*     ...Local Scalars...
      INTEGER           ICTXT, IAM, NPROW, NPCOL, MYROW, MYCOL,
     $                  SYS_NPROCS, NPROCS, AROWS, ACOLS
      INTEGER           THREADS
      INTEGER           INFO, KTOP, KBOT, ILO, IHI, SOLVER, I
      INTEGER           IPA, IPACPY, IPQ, WR1, WI1, WR2, WI2, IPW1,
     $                  IPW2, IPIW, LWORK, LIWORK
      INTEGER           HESS
      DOUBLE PRECISION  EPS, TOL1, TOL2
      DOUBLE PRECISION  STAMP, TOTTIME, T_BA, T_GEN, T_HS, T_SCH, T_QR,
     $                  T_RES
      DOUBLE PRECISION  RNORM, ANORM, R1, ORTH, O1, O2, DPDUM, ELEM1,
     $                  ELEM2, ELEM3
      INTEGER           TOTAL, CURR, PASS
*
*     ...Local Arrays...
      INTEGER           DESCA( DLEN_ ), DESCQ( DLEN_ ), DESCVEC( DLEN_ )
      DOUBLE PRECISION  SCALE( 1 000 000 )
      DOUBLE PRECISION, ALLOCATABLE :: MEM(:)
      INTEGER, ALLOCATABLE :: IMEM(:)
      INTEGER           TEMP( 7 )
      CHARACTER*( 512 ) LINE
*
*     ...Intrinsic Functions...
      INTRINSIC         INT, DBLE, SQRT, MAX, MIN
*
*     ...External Functions...
      INTEGER           NUMROC
      DOUBLE PRECISION  PDLAMCH, PDLANGE, MPI_WTIME, PCHKRESI
      EXTERNAL          BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT,
     $                  BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT,
     $                  IGEBS2D, IGEBR2D
      EXTERNAL          NUMROC, PDLAMCH, PDLASET, PDGEHRD, PDLANGE
      EXTERNAL          DGEBAL, DGEHRD
      EXTERNAL          MPI_WTIME
      EXTERNAL          PDGEBAL, PDLAHQR, PDHSEQR
      EXTERNAL          PDMATGEN2, PCHKRESI
#ifdef USE_OMP
      INTEGER           OMP_GET_NUM_THREADS
      EXTERNAL          OMP_GET_NUM_THREADS
#endif
*
*     ...Executable statements...
*
      CALL BLACS_PINFO( IAM, SYS_NPROCS )
      NPROW = INT( SQRT( DBLE(SYS_NPROCS) ) )
      NPCOL = SYS_NPROCS / NPROW
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, '2D', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF ( ICTXT.LT.0 ) GO TO 777
*
*     Read out the number of underlying threads and set stack size in
*     kilobytes.
*
#ifdef USE_OMP
*$OMP PARALLEL SHARED(THREADS)
*$OMP MASTER
      THREADS = OMP_GET_NUM_THREADS()
      IF (DEBUG) WRITE(*,*), '% OMP THREADS =', THREADS
*$OMP END MASTER
*$OMP END PARALLEL
#else
      THREADS = 1
#endif
*
*     Allocate and initialize memory with zeros.
*
      INFO = 0
      ALLOCATE ( MEM( DPALLOC ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate MEM. INFO = ', INFO
         GO TO 777
      END IF
      ALLOCATE ( IMEM( INTALLC ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate IMEM. INFO = ', INFO
         GO TO 777
      END IF
      MEM( 1:DPALLOC ) = ZERO
      IMEM( 1:INTALLC ) = IZERO
*
*     Print welcoming message.
*
      IF( IAM.EQ.0 ) THEN
         WRITE(*,*)
         WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         WRITE(*,*) '%%         TESTPROGRAM FOR PDHSEQR          %%'
         WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         WRITE(*,*)
         WRITE(*,*) '%', NPROW, 'x', NPCOL, 'processor grid'
         WRITE(*,*)
      END IF
*
*     Get machine epsilon.
*
      EPS = PDLAMCH( ICTXT, 'Epsilon' )
      TOL1 = EPS*50
      TOL2 = EPS*50
*
      IF( MYROW+MYCOL.EQ.0 ) THEN
 1       CONTINUE
*           Skip comments
            READ( *, '(A)' ) LINE
            IF( INDEX(LINE, '#').NE.0 )
     $   GO TO 1
         READ( LINE, FMT = * ) TOTAL
         IF( NPROW*NPCOL.GT.1 )
     $      CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, TOTAL, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, TOTAL, 1, 0, 0 )
      END IF
      PASS = 0
*
*     Loop over problem parameters.
*
      DO CURR = 1, TOTAL
         IF( MYROW+MYCOL.EQ.0 ) THEN
 2          CONTINUE
*              Skip comments
               READ( *, '(A)' ) LINE
               IF( INDEX(LINE, '#').NE.0 )
     $      GO TO 2
            READ( LINE, FMT = * ) TEMP(1), TEMP(2), TEMP(3), TEMP(4),
     $           TEMP(5), TEMP(6), TEMP(7)
            IF( NPROW*NPCOL.GT.1 )
     $         CALL IGEBS2D( ICTXT, 'All', ' ', 7, 1, TEMP, 7 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 7, 1, TEMP, 7, 0, 0 )
         END IF
         SOLVER = TEMP(1)
         N = TEMP(2)
         NB = TEMP(3)
         COMPHESS = TEMP(4)
         BALANCE = TEMP(5)
         KTOP = TEMP(6)
         KBOT = TEMP(7)
         KTOP = MAX(KTOP, 1)
         KTOP = MIN(KTOP, N)
         KBOT = MIN(KBOT, N)
         KBOT = MAX(KBOT, 1)
         NB = MAX(NB, 1)
         NB = MIN(NB, N)
*
         TOTTIME = MPI_WTIME()
         T_GEN = 0.0D+00
         T_RES = 0.0D+00
         T_SCH = 0.0D+00
*
*        Set INFO to zero for this run.
*
         INFO = 0
         NPROCS = NPROW*NPCOL
*
*        Count the number of rows and columns of current problem
*        for the current block sizes and grid properties.
*
         STAMP = MPI_WTIME()
         AROWS = NUMROC( N, NB, MYROW, 0, NPROW )
         ACOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
*        Set up matrix descriptors.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Set up descriptors...'
         CALL DESCINIT( DESCA, N, N, NB, NB, MIN(ARSRC,NPROW-1),
     $        MIN(NPCOL-1,ACSRC), ICTXT, MAX(1, AROWS), INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% DESCINIT DESCA failed, INFO =", INFO
            GO TO 999
         END IF
         CALL DESCINIT( DESCQ, N, N, NB, NB, MIN(ARSRC,NPROW-1),
     $        MIN(NPCOL-1,ACSRC), ICTXT, MAX(1, AROWS), INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% DESCINIT DESCQ failed, INFO =", INFO
            GO TO 999
         END IF
         CALL DESCINIT( DESCVEC, N, 1, N, 1, MIN(ARSRC,NPROW-1),
     $        MIN(NPCOL-1,ACSRC), ICTXT, N, INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% DESCINIT DESCVEC failed, INFO =", INFO
            GO TO 999
         END IF
*
*        Assign pointer for ScaLAPACK arrays - first set DP memory.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Assign pointers...'
         IPA    = 1
         IPACPY = IPA + DESCA( LLD_ ) * ACOLS
         IPQ    = IPACPY + DESCA( LLD_ ) * ACOLS
         WR1    = IPQ + DESCQ( LLD_ ) * ACOLS
         WI1    = WR1 + N
         WR2    = WI1 + N
         WI2    = WR2 + N
         IPW1   = WI2 + N
         IPW2   = IPW1 + DESCA( LLD_ ) * ACOLS
         IF( DEBUG ) WRITE(*,*) '% (IPW2,DPALLOC):', IPW2, DPALLOC
         IF( IPW2+DESCA(LLD_)*ACOLS .GT. DPALLOC+1 ) THEN
            WRITE(*,*) '% Not enough DP memory!'
            GO TO 999
         END IF
*
*        Then set integer memory pointers.
*
         IPIW = 1
*
*        Generate test problem.
*
         CALL PDLASET( 'All over', N, N, ZERO, ONE, MEM(IPQ), 1, 1,
     $        DESCQ )
         CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant',
     $        N, N, NB, NB, MEM(IPA), DESCA( LLD_ ), 0, 0, 7, 0,
     $        AROWS, 0, ACOLS, MYROW, MYCOL, NPROW, NPCOL )
         IF( COMPHESS.EQ.0 ) THEN
            CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO,
     $           MEM(IPA), 3, 1, DESCA )
            CALL PDLASET( 'All over', N, N, ZERO, ONE, MEM(IPQ),
     $           1, 1, DESCQ )
            IF( KTOP.GT.1 )
     $         CALL PDLASET( 'Lower triangular', KTOP-1, KTOP-1,
     $              ZERO, ZERO, MEM(IPA), 2, 1, DESCQ )
            IF( KBOT.LT.N )
     $         CALL PDLASET( 'Lower triangular', N-KBOT, N-KBOT,
     $              ZERO, ZERO, MEM(IPA), KBOT+1, KBOT, DESCQ )
         END IF
*
*        Do balancing if general matrix.
*
         T_BA = MPI_WTIME()
         IF( COMPHESS.NE.0 .AND. BALANCE.NE.0 ) THEN
            IF( NPROCS.EQ.1 .AND. SOLVER.NE.2 .AND. UNI_LAPACK ) THEN
               IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == dgebal =='
               CALL DGEBAL( 'Both', N, MEM(IPA), DESCA(LLD_), ILO,
     $              IHI, SCALE, INFO )
               IF ( INFO.NE.0 ) THEN
                  WRITE(*,*) "% DGEBAL failed, INFO =", INFO
                  GO TO 999
               END IF
            ELSE
               IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdgebal =='
               CALL PDGEBAL( 'Both', N, MEM(IPA), DESCA, ILO, IHI,
     $              SCALE, INFO )
               IF ( INFO.NE.0 ) THEN
                  WRITE(*,*) "% PDGEBAL failed, INFO =", INFO
                  GO TO 999
               END IF
            END IF
         ELSEIF( COMPHESS.NE.0 ) THEN
            ILO = 1
            IHI = N
         ELSE
            ILO = KTOP
            IHI = KBOT
         END IF
         T_BA = MPI_WTIME() - T_BA
         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
     $      '%%% Balancing took in seconds:',T_BA
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': ILO,IHI=',ILO,IHI
*
*        Make a copy of A.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Copy matrix A'
         CALL PDLACPY( 'All', N, N, MEM(IPA), 1, 1, DESCA, MEM(IPACPY),
     $                 1, 1, DESCA )
*
*        Print matrices to screen in debugging mode.
*
         IF( PRN )
     $      CALL PDLAPRNT( N, N, MEM(IPACPY), 1, 1, DESCA, 0, 0,
     $           'A', NOUT, MEM(IPW1) )
         T_GEN = T_GEN + MPI_WTIME() - STAMP - T_BA
         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
     $      '%%% Matrix generation took in seconds:',T_GEN
*
*        Only compute the Hessenberg form if necessary.
*
         T_HS = MPI_WTIME()
         IF( COMPHESS.EQ.0 ) GO TO 30
*
*        Reduce A to Hessenberg form.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM,
     $      ': Reduce to Hessenberg form...N=',N, ILO,IHI
         IF( NPROCS.EQ.1 .AND. SOLVER.NE.2 .AND. UNI_LAPACK ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == dgehrd =='
            CALL DGEHRD( N, ILO, IHI, MEM(IPA), DESCA(LLD_),
     $           MEM(IPW1), MEM(IPW2), -1, INFO )
            IF (DPALLOC-IPW2.LT.MEM(IPW2)) THEN
               WRITE(*,*) "% Not enough memory for DGEHRD"
               GO TO 999
            END IF
            CALL DGEHRD( N, ILO, IHI, MEM(IPA), DESCA(LLD_),
     $           MEM(IPW1), MEM(IPW2), DPALLOC-IPW2, INFO )
            IF ( INFO.NE.0 ) THEN
               WRITE(*,*) "% DGEHRD failed, INFO =", INFO
               GO TO 999
            END IF
         ELSE
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdgehrd =='
            CALL PDGEHRD( N, ILO, IHI, MEM(IPA), 1, 1, DESCA, MEM(IPW1),
     $           MEM(IPW2), -1, INFO )
            IF (DPALLOC-IPW2.LT.MEM(IPW2)) THEN
               WRITE(*,*) "% Not enough memory for PDGEHRD"
               GO TO 999
            END IF
            CALL PDGEHRD( N, ILO, IHI, MEM(IPA), 1, 1, DESCA, MEM(IPW1),
     $           MEM(IPW2), DPALLOC-IPW2, INFO )
            IF ( INFO.NE.0 ) THEN
               WRITE(*,*) "% PDGEHRD failed, INFO =", INFO
               GO TO 999
            END IF
         END IF
*
*        Form Q explicitly.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ':Form Q explicitly'
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdormhr =='
         CALL PDORMHR( 'L', 'N', N, N, ILO, IHI, MEM(IPA), 1, 1,
     $        DESCA, MEM(IPW1), MEM(IPQ), 1, 1, DESCQ, MEM(IPW2),
     $        -1, INFO )
         IF (DPALLOC-IPW2.LT.MEM(IPW2)) THEN
            WRITE(*,*) "% Not enough memory for PDORMHR"
            GO TO 999
         END IF
         CALL PDORMHR( 'L', 'N', N, N, ILO, IHI, MEM(IPA), 1, 1,
     $        DESCA, MEM(IPW1), MEM(IPQ), 1, 1, DESCQ, MEM(IPW2),
     $        DPALLOC-IPW2, INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% PDORMHR failed, INFO =", INFO
            GO TO 999
         END IF
*
*        Extract the upper Hessenberg part of A.
*
         CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO,
     $        MEM(IPA), 3, 1, DESCA )
*
*        Print reduced matrix A in debugging mode.
*
         IF( PRN ) THEN
            CALL PDLAPRNT( N, N, MEM(IPA), 1, 1, DESCA, 0, 0, 'H', NOUT,
     $           MEM(IPW1) )
            CALL PDLAPRNT( N, N, MEM(IPQ), 1, 1, DESCQ, 0, 0, 'Q', NOUT,
     $           MEM(IPW1) )
         END IF
*
 30      CONTINUE
         T_HS = MPI_WTIME() - T_HS
         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
     $      '%%% Hessenberg took in seconds:',T_HS
*
*        Compute the real Schur form of the Hessenberg matrix A.
*
         T_QR = MPI_WTIME()
         IF( SOLVER.EQ.1 ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdlahqr =='
            CALL PDLAHQR( .TRUE., .TRUE., N, ILO, IHI, MEM(IPA), DESCA,
     $           MEM(WR1), MEM(WI1), ILO, IHI, MEM(IPQ), DESCQ,
     $           MEM(IPW1), DPALLOC-IPW1+1, IMEM, INTALLC, INFO )
            IF (INFO.NE.0) THEN
               WRITE(*,*) "% PDLAHQR: INFO =", INFO
            END IF
         ELSEIF( SOLVER.EQ.2 ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdhseqr =='
            CALL PDHSEQR( 'Schur', 'Vectors', N, ILO, IHI, MEM(IPA),
     $           DESCA, MEM(WR2), MEM(WI2), MEM(IPQ), DESCQ, MEM(IPW1),
     $           -1, IMEM, -1, INFO )
            LWORK = DPALLOC-IPW1
            LIWORK = INTALLC
            IF (LWORK.LT.MEM(IPW1)) THEN
               WRITE(*,*) "% Not enough DP memory for PDHSEQR"
               GO TO 999
            END IF
            IF (LIWORK.LT.IMEM(1)) THEN
               WRITE(*,*) "% Not enough INT memory for PDHSEQR"
               GO TO 999
            END IF
            CALL PDHSEQR( 'Schur', 'Vectors', N, ILO, IHI, MEM(IPA),
     $           DESCA, MEM(WR2), MEM(WI2), MEM(IPQ), DESCQ, MEM(IPW1),
     $           LWORK, IMEM, LIWORK, INFO )
            IF (INFO.NE.0) THEN
               WRITE(*,*) "% PDHSEQR: INFO =", INFO
            END IF
         ELSE
             WRITE(*,*) '% ERROR: Illegal SOLVER number!'
             GO TO 999
         END IF
         T_QR = MPI_WTIME() - T_QR
         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
     $      '%%% QR algorithm took in seconds:',T_QR
         T_SCH = T_SCH + T_QR + T_HS + T_BA
*
*        Print reduced matrix A in debugging mode.
*
         IF( PRN ) THEN
            CALL PDLAPRNT( N, N, MEM(IPA), 1, 1, DESCA, 0, 0, 'T',
     $           NOUT, MEM(IPW1) )
            CALL PDLAPRNT( N, N, MEM(IPQ), 1, 1, DESCQ, 0, 0, 'Z',
     $           NOUT, MEM(IPW1) )
         END IF
*
*        Check that returned Schur form is really a quasi-triangular
*        matrix.
*
         HESS = 0
         DO I = 1, N-1
            IF( I.GT.1 ) THEN
               CALL PDELGET( 'All', '1-Tree', ELEM1, MEM(IPA), I, I-1,
     $              DESCA )
            ELSE
               ELEM1 = ZERO
            END IF
            CALL PDELGET( 'All', '1-Tree', ELEM2, MEM(IPA), I+1, I,
     $           DESCA )
            IF( I.LT.N-1 ) THEN
               CALL PDELGET( 'All', '1-Tree', ELEM3, MEM(IPA), I+2, I+1,
     $              DESCA )
            ELSE
               ELEM3 = ZERO
            END IF
            IF( ELEM2.NE.ZERO .AND. ABS(ELEM1)+ABS(ELEM2)+ABS(ELEM3).GT.
     $         ABS(ELEM2) ) HESS = HESS + 1
         END DO
*
*        Compute residual norms and other results:
*
*           1) RNORM = || T - Q'*A*Q ||_F / (SQRT(N) * ||A||_F)
*           2) ORTH  = MAX( || I - Q'*Q ||_F, || I - Q*Q' ||_F ) / N
*
         STAMP = MPI_WTIME()
         IF( COMPRESI .AND. .NOT. TEST_CHKRESI ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Compute residuals 1'
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': pdgemm 3'
            CALL PDGEMM( 'N', 'N', N, N, N, ONE, MEM(IPACPY), 1, 1,
     $           DESCA, MEM(IPQ), 1, 1, DESCQ, ZERO, MEM(IPW1), 1, 1,
     $           DESCA )
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': pdgemm 4'
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': N=',N
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': DESCA=',DESCA(1:DLEN_)
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': DESCQ=',DESCQ(1:DLEN_)
            CALL PDGEMM( 'T', 'N', N, N, N, -ONE, MEM(IPQ), 1, 1,
     $           DESCQ, MEM(IPW1), 1, 1, DESCA, ONE, MEM(IPA), 1, 1,
     $           DESCA )
            R1 = PDLANGE( 'Frobenius', N, N, MEM(IPA), 1, 1, DESCA,
     $           DPDUM )
            ANORM = PDLANGE( 'Frobenius', N, N, MEM(IPACPY), 1, 1,
     $           DESCA, DPDUM )
            IF( ANORM.GT.ZERO )THEN
               RNORM = R1 / (ANORM + EPS) / DBLE(N)
            ELSE
               RNORM = R1 / DBLE(N)
            END IF
         ELSEIF( COMPRESI .AND. TEST_CHKRESI ) THEN
            RNORM = PCHKRESI( N, MEM(IPACPY), 1, 1, DESCA, MEM(IPA),
     $           1, 1, DESCA, MEM(IPQ), 1, 1, DESCQ, MEM(IPW1),
     $           DPALLOC-IPW1+1 ) / DBLE(N)
         ELSE
            RNORM = 0.0D0
         END IF
*
         IF( COMPORTH ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Compute residuals 2'
            CALL PDLASET( 'All', N, N, ZERO, ONE, MEM(IPW1), 1, 1,
     $           DESCQ )
            CALL PDLACPY( 'All', N, N, MEM(IPQ), 1, 1, DESCQ, MEM(IPW2),
     $           1, 1, DESCQ )
            CALL PDGEMM( 'T', 'N', N, N, N, -ONE, MEM(IPQ), 1, 1, DESCQ,
     $           MEM(IPW2), 1, 1, DESCQ, ONE, MEM(IPW1), 1, 1, DESCQ )
            O1 = PDLANGE( 'Frobenius', N, N, MEM(IPW1), 1, 1, DESCQ,
     $           DPDUM )
            CALL PDLASET( 'All', N, N, ZERO, ONE, MEM(IPW1), 1, 1,
     $           DESCQ )
            CALL PDGEMM( 'N', 'T', N, N, N, -ONE, MEM(IPQ), 1, 1, DESCQ,
     $           MEM(IPW2), 1, 1, DESCQ, ONE, MEM(IPW1), 1, 1, DESCQ )
            O2 = PDLANGE( 'Frobenius', N, N, MEM(IPW1), 1, 1, DESCQ,
     $           DPDUM )
            ORTH = MAX(O1,O2) / DBLE(N)
         ELSE
            ORTH = 0.0D0
         END IF
*
         T_RES = T_RES + MPI_WTIME() - STAMP
         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
     $      '%%% Residual check took in seconds:',T_RES
         TOTTIME = MPI_WTIME() - TOTTIME
         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
     $      '%%% Total execution time in seconds:', TOTTIME
*
*        Print residuals.
*
         IF( DEBUG .AND. IAM.EQ.0 ) THEN
            WRITE(*,*) '%%% || T - Q^T*A*Q ||_F',
     $         '/ (N * ||A||_F) =', RNORM
            WRITE(*,*)
     $         '%%% MAX( || I - Q^T*Q ||_F, || I - Q*Q^T ||_F ) / N =',
     $         ORTH
            WRITE(*,*)
         END IF
         IF( HESS.EQ.0 .AND. RNORM.LE.TOL1 .AND. ORTH.LE.TOL2 ) THEN
            PASS = PASS+1
            IF( IAM.EQ.0 )
     $         WRITE(*,*) '% Test', CURR, 'passed!'
         ELSE
            IF( IAM.EQ.0 )
     $         WRITE(*,*) '% Test', CURR, 'failed! (',
     $              HESS, RNORM/TOL1, ORTH/TOL2, ')'
         END IF
         GO TO 888
 999     CONTINUE
            IF( IAM.EQ.0 )
     $         WRITE(*,*) '% Test', CURR, 'failed!'
 888     CONTINUE
      END DO
*
      IF( IAM.EQ.0 )
     $   WRITE(*,*) '%', PASS, 'out of', TOTAL, 'tests passed!'
*
*     Deallocate MEM and IMEM.
*
      DEALLOCATE( MEM, IMEM )
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
 777  CONTINUE
*
      CALL BLACS_EXIT( 0 )
*
      END
