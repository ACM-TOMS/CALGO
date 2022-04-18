CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PDEXTCHK( M, A, IA, JA, DESCA, EXTINFO, LEXTINFO, 
     $                     DOEXT, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 20, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar arguments ..
      INTEGER              M, IA, JA, INFO, LEXTINFO
      LOGICAL              DOEXT
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION     A( * )
      INTEGER              DESCA( * ), EXTINFO( * )
C     ..
C     
C  Purpose
C  =======
C  This subroutine takes as input the real Schurform of a matrix
C  and checks what diagonal blocks that have to be extended one row
C  and one column to not lose elements that belong to any 2-by-2 
C  diagonal block and are split between blocks (processors). 
C
C  The output is in EXTINFO where its i:th element corresponds
C  to the i:th diagonal block of A and has the following meaning:
C
C  If EXTINFO(i) = 0 - Aii shall be unchanged
C  If EXTINFO(i) = 1 - Aii shall be extended one row and one 
C                      column to the South and the East
C  If EXTINFO(i) = 2 - Aii shall be diminished one row and one
C                      column to the North and to the West
C  If EXTINFO(i) = 3 - Aii shall be extended and diminished
C                      as in the cases above.
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
C  Input arguments
C
C  M         (global input) INTEGER
C            The number of rows and columns of the global distributed 
C            matrix sub(A). M >= 0.
C
C  A         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A. On output,
C            the local parts of the distributed matrix A in upper
C            Hessenberg form.
C
C  IA        (global input) INTEGER
C            Row start index for sub(A), i.e., the submatrix to operate
C            on. IA >= 1.
C
C  JA        (global input) INTEGER
C            Column start index for sub(A), i.e., the submatrix to 
C            operate on. JA = IA must hold.
C
C  DESCA     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix A.
C
C  Output information
C
C  EXTINFO   (global output) INTEGER array of dimension LEXTINFO
C            
C
C  LEXTINFO  (global input) INTEGER
C            The length of the array EXTINFO. 
C            LEXTINFO >= ICEIL( M + MOD( IA - 1, MB_A ), MB_A ).
C
C  DOEXT     (global output) LOGICAL
C            On output, this scalar is .TRUE. if an implicit 
C            redistribution is necessary, otherwise DOEXT is
C            returned as .FALSE..
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
C  This subroutine checks the need of an implicit redistribution of
C  a real Schur form. See the references for details.
C
C  Additional requirements
C  =======================
C  None.
C
C  Limitations
C  ===========
C  No limitations.
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
C  =====================================================================
C
C     .. Local Parameters ..
      INTEGER          CTXT_, MB_, LLD_
      DOUBLE PRECISION ZERO
      PARAMETER        ( CTXT_ = 2, MB_ = 5, LLD_ = 9, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER          MB, LLDA, I, DBA, IX, JX, RSRC, CSRC, NPROC,
     $                 ICTXT, NPROW, NPCOL, MYROW, MYCOL, IROFFA,
     $                 ISTART, II
C     ..
C     .. External Subroutines ..
      EXTERNAL         INFOG2L, IGEBS2D, IGEBR2D, BLACS_GRIDINFO
C     ..
C     .. External Functions ..
      INTEGER          ICEIL
      EXTERNAL         ICEIL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        MOD
C     ..
C     .. Executable Statements ..
C
      ICTXT = DESCA( CTXT_ )
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
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( M, 1, M, 1, IA, JA, DESCA, 5, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( IA.LT.1 ) THEN
               INFO = -3
            ELSEIF( JA.NE.IA ) THEN
               INFO = -4
            END IF
         END IF
      END IF
C 
C     Get the blocksize and the local leading dimension of A.
C     Compute the number of diagonal blocks of A and set the logical
C     variable to false.
C     
      IF( INFO.EQ.0 ) THEN
         MB = DESCA( MB_ )
         LLDA = DESCA( LLD_ )
         IROFFA = MOD( IA - 1, MB ) 
         DBA = ICEIL( M + IROFFA, MB )
         NPROC = NPROW*NPCOL
         DOEXT = .FALSE.
      END IF
C
C     Check LEXTINFO
C
      IF( INFO.EQ.0 ) THEN
         IF( LEXTINFO.LT.DBA ) INFO = -7
      END IF
C     
C     Error messages
C
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDEXTCHK', -INFO )
      END IF
C     
C     Only do extension check if we have more than one block
C     
      IF( ICEIL(IA,MB).NE.ICEIL(IA+M-1,MB) ) THEN
C     
C     Handle the first block by itself since it has to predecessor
C     
         ISTART = ICEIL( IA, MB )
         I = 1
         CALL INFOG2L( ISTART * MB + 1 , ISTART * MB, DESCA, NPROW, 
     $        NPCOL, MYROW, MYCOL, IX, JX, RSRC, CSRC )
         IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
            IF( A((JX - 1) * LLDA + IX).NE.ZERO ) THEN
               EXTINFO( I ) = 1
            ELSE
               EXTINFO( I ) = 0
            END IF
         END IF   
C     
C     Since the information shall be global we must broadcast it
C     
         IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC .AND. NPROC.GT.1 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, EXTINFO( I ), 1 )
         ELSEIF( NPROC.GT.1 ) THEN
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, EXTINFO( I ), 1, 
     $                    RSRC, CSRC )
         ENDIF
C     
C     Handle rest of the diagonal blocks, except for the last one: 
C     if I shall extend but my predecessor won't I will be declared 
C     as extended (1), if both I and my predecessor shall extend I will
C     be declared as both extended and diminished (3), if I shall not 
C     extend and my predecessor won't extend either I will be delared as
C     unchanged (0) and, finally, if I shall not extend and my 
C     predecessor shall extend I will declared as diminished (2). 
C     The information shall be global so we must broadcast each entry 
C     when it has been set.
C     
         DO 10 II = ISTART + 1, ISTART + DBA - 1
            I = II - ISTART + 1
            CALL INFOG2L( II * MB + 1 , II * MB , DESCA, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RSRC, CSRC )
            IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
               IF( A((JX - 1) * LLDA + IX).NE.ZERO ) THEN
                  IF( EXTINFO( I-1 ).EQ.0 .OR. EXTINFO( I-1 ).EQ.2 )
     $                 EXTINFO( I ) = 1
                  IF( EXTINFO( I-1 ).EQ.1 .OR. EXTINFO( I-1 ).EQ.3 )
     $                 EXTINFO( I ) = 3
               ELSE
                  IF( EXTINFO( I-1 ).EQ.0 .OR. EXTINFO( I-1 ).EQ.2 ) 
     $                 EXTINFO( I ) = 0
                  IF( EXTINFO( I-1 ).EQ.1 .OR. EXTINFO( I-1 ).EQ.3 )
     $                 EXTINFO( I ) = 2
               END IF
               IF( NPROC.GT.1 ) 
     $              CALL IGEBS2D( ICTXT,'All', ' ', 1, 1, EXTINFO( I ), 
     $                            1 )
            ELSE
               IF( NPROC.GT.1 ) THEN
                  CALL IGEBR2D( ICTXT,'All', ' ', 1, 1, EXTINFO( I ), 
     $                            1, RSRC, CSRC )
               END IF
            END IF
 10      CONTINUE
C     
C     Now handle the last diagonal block by itself, since it cannot
C     be extended, just diminished or unchanged. All the needed 
C     information is now global - no need for BC.
C     
         I = DBA
         IF( EXTINFO( I-1 ).EQ.1 .OR. EXTINFO( I-1 ).EQ.3 ) THEN
            EXTINFO( I ) = 2
         ELSE
            EXTINFO ( I ) = 0
         END IF
C
      ELSE
         EXTINFO( 1 ) = 0
      END IF
C     
C     Go through the information and set the logical variable
C
      DO 20 I = 1, DBA
         IF( EXTINFO(I).NE.0 .AND. .NOT.DOEXT ) THEN
            DOEXT = .TRUE.
         END IF
 20   CONTINUE
C     
      END
C     
C     End of PDEXTCHK
C     
C *** Last line of PDEXTCHK ***
