CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PDBCKRD( M, N, C, IC, JC, DESCC, EXTINFA, EXTINFB, 
     $                    EXC, LEXC, DWORK, LDWORK, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 21, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar arguments ..
      INTEGER              M, N, IC, JC, LEXC, LDWORK, INFO
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION     C( * ), EXC( * ), DWORK( * )
      INTEGER              DESCC( * ), EXTINFA( * ), EXTINFB( * ) 
C     ..   
C  
C  Purpose and description
C  =======================
C  This subroutine does an implicit backredistribution of the elements
C  in the matrix sub(C) which has been implicitely redistributed by the
C  subroutine PDIMPRED.
C
C  Notes
C  =====
C
C  Each global data object is described by an associated description
C  vector called DESC_.  This vector stores the information required to 
C  establish the mapping between an object element and its corresponding 
C  process and memory location. 
C
C  Let A be a generic term for any 2D block cyclicly distributed array.
C  Such a global array has an associated description vector DESCA.  In
C  the following comments, the character _ should be read as "of the
C  global array".
C
C  NOTATION        STORED IN      EXPLANATION
C  --------------- -------------- --------------------------------------
C  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
C                                 DTYPE_A = 1.
C  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
C                                 the BLACS process grid A is distribu-
C                                 ted over. The context itself is glo-
C                                 bal, but the handle (the integer
C                                 value) may vary.
C  M_A    (global) DESCA( M_ )    The number of rows in the global
C                                 array A.
C  N_A    (global) DESCA( N_ )    The number of columns in the global
C                                 array A.
C
C
C  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
C                                 the rows of the array.
C  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
C                                 the columns of the array.
C  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
C                                 row of the array A is distributed.
C  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
C                                 first column of the array A is
C                                 distributed.
C  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
C                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
C
C  Let K be the number of rows or columns of a distributed matrix, and
C  assume that its process grid has dimension p x q.
C  LOCr( K ) denotes the number of elements of K that a process would
C  receive if K were distributed over the p processes of its process
C  column.
C  Similarly, LOCc( K ) denotes the number of elements of K that a
C  process would receive if K were distributed over the q processes of
C  its process row.
C  The values of LOCr() and LOCc() may be determined via a call to the
C  ScaLAPACK tool function, NUMROC:
C          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ). 
C  
C  An upper bound for these quantities may be computed by:
C          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C
C  Arguments
C  =========
C
C  Input/Output parameters
C
C  M         (global input) INTEGER
C            The number of rows of the global distributed matrix C. 
C            M >= 0.
C
C  N         (global input) INTEGER
C            The number of columns of the global distributed matrix C. 
C            N >= 0.
C
C  C         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_C,LOCc(N)). 
C            On entry C contains the local pieces of the global 
C            distributed matrix C. 
C
C  IC        (global input) INTEGER
C            Row start index for sub(C), i.e., the submatrix to operate 
C            on. MOD(IC,MB_A) = MOD(IA,MB_A) must hold. 
C
C  JC        (global input) INTEGER
C            Column start index for sub(C), i.e., the submatrix to 
C            operate on. MOD(JC,MB_B) = MOD(JB,MB_B) must hold.
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix C.
C
C  EXTINFA   (global input) INTEGER array of dimension
C            ICEIL( M + MOD( IA - 1, MB_A ), MB_A ). Carries the 
C            extension information for the matrix A computed by
C            PDEXTCHK.
C  
C  EXTINFB   (global input) INTEGER array of dimension
C            ICEIL( N + MOD( IJ - 1, MB_B ), MB_B ). Carries the 
C            extension information for the matrix B computed by
C            PDEXTCHK.
C
C  EXC       (local input) DOUBLE PRECISION array, dimension LEXC.
C            Stores the extension elements associated with the matrix
C            C for usage in calls to the extended block building 
C            routines DBEXMAT and DUBEXMA.
C
C  LEXC      (local input) INTEGER, the length of the array EXC.
C            LEXC >= ICEIL( LOCr( M + MOD( IC - 1, MB_C ), MB_C ) ) * 
C                    ICEIL( LOCc( N + MOD( JC - 1, NB_C ), NB_C ) ) * 
C                    ( MB_C + NB_C + 1 )
C  
C  Workspace
C
C  DWORK     (local workspace) DOUBLE PRECISION array, dimension
C            LDWORK. 
C
C  LDWORK    (local or global input) INTEGER
C            The dimension of the array DWORK.
C            LDWORK >= MAX( MB_C, NB_C )
C 
C            If LDWORK = -1, LDWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in DWORK(1) and return
C            immediately. No error will then be signaled by PXERBLA.
C
C  Error indicator
C
C  INFO      (global output) INTEGER
C             = 0:  successful exit
C             < 0:  unsuccessful exit. 
C             If the i-th argument is an array and the j-entry had
C             an illegal value, then INFO = -(i*100+j), if the i-th
C             argument is a scalar and had an illegal value, then
C             INFO = -i. If INFO = 1, we had no valid BLACS context.
C
C  Method
C  ======
C  See the references.
C
C  Additional requirements
C  =======================
C
C  A and B must be distributed using the same blocking factor in 
C  each direction, i.e., MB_A=NB_A, MB_B=NB_B. Moreover, for C the 
C  blocksize in the row direction must agree with A's, i.e. MB_C=MB_A
C  must hold, and the blocksize in the column direction must agree with
C  B's, i.e. NB_C=NB_B must hold.
C
C  Limitations
C  ===========
C  None.
C
C  References
C  ==========
C
C  [1] Robert Granat and Bo Kågström, Parallel ScaLAPACK-style Algorithms 
C      for Standard and Generalized Sylvester-Type Matrix Equations, in 
C      preparation, Department of Computing Science and HPC2N, Umeå 
C      University, 2006.
C 
C  [2] Robert Granat and Bo Kågström, SCASY - A ScaLAPACK-style High 
C      Performance Library for Sylvester-Type Matrix Equations, in 
C      preparation, Department of Computing Science and HPC2N, Umeå 
C      University, 2006. 
C
C  See also: http://www.cs.umu.se/research/parallel/scasy
C
C  Parallel execution recommendations
C  ==================================
C  None.
C
C  Revisions
C  =========
C  Please report bugs to <granat@cs.umu.se>.
C
C  Keywords
C  ========
C  Implicit redistribution, real Schur form.
C
C  ======================================================================
C     
C     .. Local Parameters ..
      INTEGER          CTXT_, MB_, NB_, RSRC_, CSRC_, LLD_, DLEN_
      PARAMETER        ( CTXT_ = 2, MB_ = 5, NB_ = 6, RSRC_ = 7, 
     $                   CSRC_ = 8, LLD_ = 9, DLEN_ = 9 )
C     ..
C     .. Local Scalars ..
      LOGICAL          LQUERY
      INTEGER          LLDC, CROWS, CCOLS, RSRC, CSRC, ERSRC, ECSRC, 
     $                 SRSRC, SCSRC, SERSRC, SECSRC, IX, JX, IX2, JX2, 
     $                 IX3, JX3, IX4, JX4, I, J, K, DBA, DBB, NPROCS, 
     $                 POS, ISTART, JSTART, LBI, LBJ, MB, NB, NBC, MYI, 
     $                 MYJ, LENG, LEN, ST, MYROW, MYCOL, NPROW, NPCOL,
     $                 ICTXT, IROFFC, ICOFFC, LIC, LJC, CRSRC, CCSRC, 
     $                 SCROWS, SCCOLS, JSX, ISX, EXMEMC, RSRC1, CSRC1
      DOUBLE PRECISION CORNER
C     ..
C     .. Local Arrays ..
      INTEGER          DESCSC( DLEN_ )
C     ..
C     .. External Subroutines ..
      EXTERNAL         INFOG2L, DGESD2D, DGERV2D, DLACPY, BLACS_GRIDINFO
C     ..
C     .. External functions ..
      INTEGER ICEIL, NUMROC
      EXTERNAL ICEIL, NUMROC
C     ..
C     .. Executable Statements ..
C 
      INFO = 0
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C
C     Check if the context is valid
C
      IF( NPROW.EQ.-1 ) THEN
         INFO = 1
      END IF
C
C     Check input arguments
C 
      LQUERY = LDWORK.EQ.-1
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( M, 1, N, 2, IC, JC, DESCC, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( IC.LT.1 ) THEN
               INFO = -4
            ELSEIF( JC.LT.1 ) THEN
               INFO = -5
            END IF
         END IF
      END IF
C      
C     Compute the number of diagonal blocks in A and B and some other
C     values
C     
      IF( INFO.EQ.0 ) THEN
         MB = DESCC( MB_ )
         NB = DESCC( NB_ )
         IROFFC = MOD( IC - 1, MB ) 
         DBA = ICEIL( M + IROFFC, MB )
         ICOFFC = MOD( JC - 1, NB ) 
         DBB = ICEIL( N + ICOFFC, NB )
         LLDC = DESCC( LLD_ )
         NPROCS = NPROW * NPCOL
C     
C     Compute the number of column of sub(A) and sub(B) that 
C     belong to me
C     
         CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, LIC, 
     $                 LJC, CRSRC, CCSRC )
         SCROWS = NUMROC( M + IROFFC, DESCC(MB_), MYROW, CRSRC, NPROW )
         SCCOLS = NUMROC( N + ICOFFC, DESCC(NB_), MYCOL, CCSRC, NPCOL )
         EXMEMC = ICEIL(SCROWS,MB)*ICEIL(SCCOLS,NB)*(MB+NB+1)
C     
C     Check length of output arrays
C     
         IF( LEXC.LT.EXMEMC ) INFO = -10
      END IF
C     
C     Check workspace
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT. LQUERY .AND. LDWORK.LT. MAX( MB, NB ) ) THEN
            INFO = -18
         ELSEIF( LQUERY ) THEN
            INFO = 0
            DWORK( 1 ) = MAX( MB, NB )
            RETURN
         END IF
      END IF
C     
C     Error messages
C
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDBCKRD', -INFO )
      END IF              
C     
C     Init matrix descriptor for smallest submatrix of C including 
C     sub(C) which conforms with ScaLAPACK conventions. We need this
C     descriptor for simplifying the computations of indicies in EXC 
C     later on.
C
      CALL DESCINIT ( DESCSC, M + IROFFC, N + ICOFFC, DESCC(MB_), 
     $                DESCC(NB_), CRSRC, CCSRC, ICTXT, LLDC, INFO )
C
C     Compute the number of block columns of sub(C) that belongs to me
C
      NBC = ICEIL( SCCOLS, NB )
C
C     Now do the 4 sweeps over C
C
      DO 10 K = 1, 4
C
C     Select which MBxNB block in each 2x2 block to operate on
C     in this sweep
C
        IF( K.EQ.1.OR.K.EQ.2 ) ISTART = 1
        IF( K.EQ.3.OR.K.EQ.4 ) ISTART = 2
        IF( K.EQ.1.OR.K.EQ.3 ) JSTART = 1
        IF( K.EQ.2.OR.K.EQ.4 ) JSTART = 2
C
C     Loop through all those blocks in C
C
         DO 20 I = ISTART, DBA, 2
            DO 30 J = JSTART, DBB, 2
C
C     Check if this block was extended somehow
C
               IF( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 .OR.
     $              EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 ) THEN     
C     
C     Check what processes that owns Cij
C     
                  CALL INFOG2L( (I-1) * MB + IC - IROFFC, 
     $                          (J-1) * NB + JC - ICOFFC,
     $                          DESCC, NPROW, NPCOL, MYROW, MYCOL, IX, 
     $                          JX, RSRC, CSRC )
                  CALL INFOG2L( (I-1) * MB + 1, 
     $                          (J-1) * NB + 1,
     $                          DESCSC, NPROW, NPCOL, MYROW, MYCOL, ISX, 
     $                          JSX, RSRC1, CSRC1 )
C
C     Check how Cij was extended and choose branch from this.
C     Here it is possible that Cij might have been extended such
C     that not all extension elements were addressed during the
C     solves and updates. This means that when the extension
C     row or column has been reveived by its rightful owner, that
C     guy has to check if he shall address the whole array or skip
C     the first element. This is checked by looking at the extension
C     of the to Cmyimyj corresponding Bmyjmyj and Amyimyi. If Bmyjmyj
C     is diminished we skip the first element of the incoming row.
C     Similarly, if Amyimyi is diminished we skip the first element
C     of the incoming column.
C
C     Was Cij extended with an extra row? Then check who shall have that
C     row back and send it back.
C     
                  IF( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ) THEN
                     CALL INFOG2L( I * MB + IC - IROFFC, 
     $                             (J-1) * NB + JC - ICOFFC, 
     $                             DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IX2, JX2, SRSRC, SCSRC )
                     IF( ( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ).OR.
     $                    ( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) ) THEN
C
C     Set the message length and start exchangning data
C
                        LEN = MIN( NB, N + ICOFFC - (J-1)*NB )
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           LBI = ICEIL( ISX, MB )
                           LBJ = ICEIL( JSX, NB )
                           POS = ( (LBI-1)*NBC + (LBJ-1) )*
     $                           ( MB + NB + 1) + 1  
                           CALL DLACPY( 'All', LEN, 1, EXC(POS), NB,
     $                                  DWORK, NB )
                           IF( NPROW.GT.1 ) 
     $                          CALL DGESD2D( ICTXT, LEN, 1, DWORK, NB,
     $                                        SRSRC, SCSRC )
                        END IF
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           IF( NPROW.GT.1 ) 
     $                          CALL DGERV2D( ICTXT, LEN, 1, DWORK, NB,
     $                                        RSRC, CSRC )
                           MYJ = J
                           IF( EXTINFB( MYJ ).EQ.2 .OR. 
     $                          EXTINFB( MYJ ).EQ.3 ) THEN
                              LENG = LEN - 1
                              JX2 = JX2 + 1
                              ST = 2
                           ELSE
                              LENG = LEN
                              ST = 1
                           END IF
                           CALL DLACPY( 'All', 1, LENG, DWORK( ST ), 1, 
     $                                  C((JX2-1)*LLDC + IX2), LLDC )
                        END IF
                     END IF
                  END IF
C
C     Was Cij extended with an extra column? Then send it back.
C
                  IF( EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 ) THEN
                     CALL INFOG2L( (I-1) * MB + IC - IROFFC, 
     $                             J * NB + JC - ICOFFC, 
     $                             DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IX3, JX3, ERSRC, ECSRC )
                     IF( ( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) .OR.
     $                    ( MYROW.EQ.ERSRC.AND.MYCOL.EQ.ECSRC )) THEN
C
C     Set the message length and start exchangning data
C
                        LEN = MIN( MB, M + IROFFC - (I-1)*MB )
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           LBI = ICEIL( ISX, MB )
                           LBJ = ICEIL( JSX, NB )
                           POS = ( (LBI-1)*NBC + (LBJ-1) )*
     $                           ( MB + NB + 1) + 1 + NB + 1
                           CALL DLACPY( 'All', LEN, 1, EXC(POS), MB,
     $                                  DWORK, MB )
                           IF( NPCOL.GT.1 ) 
     $                          CALL DGESD2D( ICTXT, LEN, 1, DWORK, MB,
     $                                        ERSRC, ECSRC )
                        END IF
                        IF( MYROW.EQ.ERSRC.AND.MYCOL.EQ.ECSRC ) THEN
                           IF( NPCOL.GT.1 ) 
     $                          CALL DGERV2D( ICTXT, LEN, 1, DWORK, MB,
     $                                        RSRC, CSRC )
                           MYI = I
                           IF( EXTINFA( MYI ).EQ.2 .OR.
     $                          EXTINFA( MYI ).EQ.3 ) THEN
                              LENG = LEN - 1
                              IX3 = IX3 + 1
                              ST = 2
                           ELSE
                              LENG = LEN
                              ST = 1
                           END IF
                           CALL DLACPY( 'All', LENG, 1, DWORK( ST ), MB, 
     $                          C((JX3-1)*LLDC + IX3), LLDC)
                        END IF
                     END IF
                  END IF 
C
C     Was Cij extended with an extra corner element? Then send it back.
C
                  IF( ( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ).AND.
     $                ( EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 ) )THEN   
                     CALL INFOG2L( I * MB + IC - IROFFC, 
     $                             J * NB + JC - ICOFFC,
     $                             DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IX4, JX4, SERSRC, SECSRC)
                     IF( ( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ).OR.
     $                    ( MYROW.EQ.SERSRC.AND.MYCOL.EQ.SECSRC ) ) THEN
                         IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           LBI = ICEIL( ISX, MB )
                           LBJ = ICEIL( JSX, NB )
                           POS = ( (LBI-1)*NBC + (LBJ-1) )*
     $                           ( MB + NB + 1) + 1 + NB 
                           CORNER = EXC( POS )
                           IF( NPROCS.GT.1 ) 
     $                          CALL DGESD2D( ICTXT, 1, 1, CORNER, 1,
     $                                        SERSRC, SECSRC )
                        END IF
                        IF( MYROW.EQ.SERSRC.AND.MYCOL.EQ.SECSRC ) THEN
                           IF( NPROCS.GT.1 ) 
     $                          CALL DGERV2D( ICTXT, 1, 1, CORNER, 1,
     $                                        RSRC, CSRC )
                           C((JX4-1)*LLDC + IX4) = CORNER
                        END IF
                     END IF
                  END IF
               END IF
C     
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C     
      END
C
C     End of PDBCKRD
C
C *** Last line of PDBCKRD ***
