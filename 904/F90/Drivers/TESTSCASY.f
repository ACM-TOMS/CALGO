CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 1.0beta, July 31, 2009.                        C
C         Copyright 2009, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM TESTSCASY
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SCASY version 1.0beta, July 31, 2009.                            C
C     Written by Robert Granat.                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Purpose
C     =======
C
C     Testprogram for ScaLAPACK-style software package SCASY.
C     The matrix equations considered in SCASY are the following:
C
C     Standard Sylvester (CT): op(A)*X +/- X*op(B)  = C            (SYCT)
C     Standard Lyapunov  (CT): op(A)*X + X*op(A')   = C            (LYCT)
C     Standard Sylvester (DT): op(A)*X*op(B) +/- X  = C            (SYDT)
C     Standard Lyapunov  (DT): op(A)*X*op(A') +/- X = C            (LYDT)
C                              |op(A)*X +/- Y*op(B) = C 
C     Gen. Coupled Sylvester : |                                   (GCSY)
C                              |op(D)*X +/- Y*op(E) = F
C     Gen. Sylvester         : op(A)*X*op(B) +/- op(C)*X*op(D) = E (GSYL)
C     Gen. Lyapunov      (CT): op(A)*X*op(E') - op(E)*X*op(A') = C (GLYCT)
C     Gen. Lyapunov      (DT): op(A)*X*op(A') + op(E)*X*op(E') = C (GLYCT)
C
C     This testprogram runs through a user-given set of parameters and 
C     tests all solver routines and condition estimators in SCASY using
C     the included matrix and matrix pair generators. 
C
C     References
C     ==========
C
C     [1] Robert Granat and Bo Kågström, Parallel ScaLAPACK-style Algorithms 
C     for Standard and Generalized Sylvester-Type Matrix Equations, in 
C     preparation, Department of Computing Science and HPC2N, Umeå 
C     University, 2008.
C 
C     [2] Robert Granat and Bo Kågström, SCASY - A ScaLAPACK-style High 
C     Performance Library for Sylvester-Type Matrix Equations, in 
C     preparation, Department of Computing Science and HPC2N, Umeå 
C     University, 2008. 
C
C     See also: http://www.cs.umu.se/research/parallel/scasy
C
C     Revision history
C     ================
C
C     None. Please send comments and bug-reports to granat@cs.umu.se.
C
C     User instructions
C     =================
C     
C     Modify the parameters given below for each matrix equation
C     specifying the range of problem size, submatrix starting 
C     indicies, grid parameters etc. The program will then run
C     through each instance of the parameter space and write the 
C     results to the specified output. The following parameters
C     may be modified by the user:
C
C     NAME               TYPE         DESCRIPTION   
C     ------------------------------------------------------------------
C     SYCT                 LOGICAL    Decides if SYCT is considered or
C                                     not (.TRUE. or .FALSE.). 
C     LYCT                 LOGICAL    Decides if LYCT is considered or
C                                     not (.TRUE. or .FALSE.). 
C     SYDT                 LOGICAL    Decides if SYDT is considered or
C                                     not (.TRUE. or .FALSE.).
C     LYDT                 LOGICAL    Decides if LYDT is considered or
C                                     not (.TRUE. or .FALSE.).
C     GCSY                 LOGICAL    Decides if GCSY is considered or
C                                     not (.TRUE. or .FALSE.).
C     GSYL                 LOGICAL    Decides if GSYL is considered or
C                                     not (.TRUE. or .FALSE.).
C     GLYCT                LOGICAL    Decides if GLYCT is considered or
C                                     not (.TRUE. or .FALSE.).
C     GLYDT                LOGICAL    Decides if GLYDT is considered or
C                                     not (.TRUE. or .FALSE.).
C     NODEMEM              INTEGER    The number of bytes in each node's
C                                     memory of the current machine to 
C                                     use for double precision storage.
C     INODEMEM             INTEGER    The number of bytes in each node's
C                                     memory of the current machine to 
C                                     use for integer number storage.
C     PRNTSIZ              INTEGER    All matrices will be printed to
C                                     output NOUT if the dimensions are
C                                     smaller than or equal to PRNTSIZ.
C     MAX_CPUS             INTEGER    The maximum number of available
C                                     CPUs on the current machine.
C     NOUT                 INTEGER    The prefered output unit.
C     [ACRO]_CONDEST       LOGICAL    Perform condition estimation of
C                                     matrix equation [ACRO] or not
C                                     (.TRUE. or .FALSE.).
C     [ACRO]_REDUCE_ONLY   LOGICAL    Do not solve the equation - only
C                                     reduce it to triangular form 
C                                     or .not. (.TRUE. or .FALSE.).  
C     [ACRO]_SOLVE         LOGICAL    Do solve the equation or not 
C                                     (.TRUE. or .FALSE.).
C     [ACRO]_PIPELINE      LOGICAL    Turn pipelining on and off
C                                     (.TRUE. or .FALSE.). In case of
C                                     .FALSE., the loops over the 
C                                     pipelining block factors are 
C                                     ignored.
C     [ACRO]_SQUARE_DIM    LOGICAL    Generate only square right hand 
C                                     sides for [ACRO].
C     [ACRO]_SQUARE_BLK    LOGICAL    Use only square block in right
C                                     hand side of [ACRO].
C     [ACRO]_COMM          CHARACTER  Decides which communication pattern
C                                     to utilize in triangular solver. 
C     [ACRO]_[X]FORM       CHARACTER  The generated matrix or matrix 
C                                     pair from P1SQMATGD will be 
C                                     returned in upper triangular form 
C                                     ('U') or not ('_').   
C     [ACRO]_[X]DIAG       CHARACTER  Supply explicit diagonal element
C                                     to P1SQMATGD ('D') or not ('_').
C     [ACRO]_UPLOSIGN_MIN  INTEGER    The sign and transpose variants 
C     [ACRO]_UPLOSIGN_MAX             going from min to max to be 
C                                     considered as described in RECSY.
C     [ACRO]_[X]SCHUR      CHARACTER  The matrix or matrix pair will be
C                                     passed on to the solver routine
C                                     in Schur form ('S') or not ('N').
C     [ACRO]_[X]_SYMM      CHARACTER  The right hand side matrix is 
C                                     assumed to be symmetric ('S') or
C                                     ('N').
C     [ACRO]_TRANSZ        CHARACTER  The transposed Kronecker product
C                                     representation will be used in the
C                                     solution process (only used for 
C                                     GCSY).
C     [ACRO]_[X]_NQTRBL    INTEGER    Fraction of (pairs of) diagonal 
C                                     blocks of matrix (pair) being
C                                     (pairs of) 2-by-2 blocks 
C                                     corresponding to complex 
C                                     conjugate pairs of eigenvalues 
C                                     in percentage (0-100).
C     [ACRO]_M_MIN         INTEGER    The minumum number of rows of
C                                     right hand side matrix (pair).
C     [ACRO]_M_MAX         INTEGER    The maximum number of rows of
C                                     right hand side matrix (pair).
C     [ACRO]_M_STEP        INTEGER    Loop increment going from
C                                     [ACRO]_M_MIN to [ACRO]_M_MAX.
C     [ACRO]_N_MIN         INTEGER    The minimum number of columns of
C                                     right hand side matrix (pair). 
C     [ACRO]_N_MAX         INTEGER    The minimum number of columns of
C                                     right hand side matrix (pair).
C     [ACRO]_N_STEP        INTEGER    Loop increment going from
C                                     [ACRO]_N_MIN to [ACRO]_N_MAX.
C     [ACRO]_MB_MIN        INTEGER    The minimum row blocking factor
C                                     utilized for all matrices.
C     [ACRO]_MB_MAX        INTEGER    The maximum row blocking factor
C                                     utilized for all matrices.
C     [ACRO]_MB_STEP       INTEGER    Loop increment going from
C                                     [ACRO]_MB_MIN to [ACRO]_MB_MAX.
C     [ACRO]_NB_MIN        INTEGER    The minumum column blocking factor
C                                     utlized for all matrices.
C     [ACRO]_NB_MAX        INTEGER    The maxiumum column blocking factor
C                                     utlized for all matrices.
C     [ACRO]_NB_STEP       INTEGER    Loop increment going from
C                                     [ACRO]_NB_MIN to [ACRO]_NB_MAX.
C     [ACRO]_MB2_MIN       INTEGER    The minimum row blocking factor
C                                     utilized for all matrices in 
C                                     deep pipelining.
C     [ACRO]_MB2_MAX       INTEGER    The maximum row blocking factor
C                                     utilized for all matrices in
C                                     deep pipelining.
C     [ACRO]_MB2_STEP      INTEGER    Loop increment going from
C                                     [ACRO]_MB2_MIN to [ACRO]_MB2_MAX.
C     [ACRO]_NB2_MIN       INTEGER    The minumum column blocking factor
C                                     utlized for all matrices in
C                                     deep pipelining.
C     [ACRO]_NB2_MAX       INTEGER    The maxiumum column blocking factor
C                                     utlized for all matrices in
C                                     deep pipelining.
C     [ACRO]_NB2_STEP      INTEGER    Loop increment going from
C                                     [ACRO]_NB2_MIN to [ACRO]_NB2_MAX.
C     [ACRO]_I[X]_MIN      INTEGER    The minimum row and column starting 
C                                     index for sub([X]).
C     [ACRO]_I[X]_MAX      INTEGER    The maximum row and column starting 
C                                     index for sub([X]).
C     [ACRO]_I[X]_STEP     INTEGER    Loop increment going from
C                                     [ACRO]_I[X]_MIN to [ACRO]_I[X]_MAX.
C     [ACRO]_NPROW_MIN     INTEGER    The minimum number of process rows
C                                     utilized.
C     [ACRO]_NPROW_MAX     INTEGER    The maximum number of process rows
C                                     utilized.
C     [ACRO]_NPROW_STEP    INTEGER    Loop increment going from
C                                     [ACRO]_NPROW_MIN to [ACRO]_NPROW_MAX.
C     [ACRO]_NPCOL_MIN     INTEGER    The minimum number of process 
C                                     columns utilized.
C     [ACRO]_NPCOL_MAX     INTEGER    The maximum number of process
C                                     columns utilized.
C     [ACRO]_NPCOL_STEP    INTEGER    Loop increment going from
C                                     [ACRO]_NPCOL_MIN to [ACRO]_NCOL_MAX.
C     [ACRO]_REPEAT        INTEGER    Repeat each instance a given 
C                                     number of times
C
C     where [ACRO] is replaced with the acronym for the matrix equation 
C     in question, [X] is replaced with the matrix or matrix pair 
C     name(s) in question and ('_') denotes any character. The given 
C     set of parameters are checked for consistency. If any 
C     inconsistency is detected, the program is aborted.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === THE FOLLOWING PART MAY BE MODIFIED BY THE USER ===
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     General parameters: what equations to solve, node memory in bytes 
C     (for both double precision and integers), maximum number of CPUs, 
C     maximum number of threads for OpenMP versions of RECSY node 
C     solvers and SMP-BLAS, matrix sizes for printing, output channel. 
C       
C
      LOGICAL          SYCT, LYCT, SYDT, LYDT, GCSY, GSYL, GLYCT, GLYDT
      INTEGER          PRNTSIZ, MAX_CPUS, NOUT
#ifdef USE_INTEGER8
      INTEGER*8        NODEMEM, INODEMEM
#else
      INTEGER          NODEMEM, INODEMEM
#endif
      INTEGER          INTKIND
      PARAMETER        ( SYCT =  .TRUE., 
     $                   LYCT =  .TRUE.,
     $                   SYDT =  .TRUE.,  
     $                   LYDT =  .TRUE., 
     $                   GCSY =  .TRUE., 
     $                   GSYL =  .TRUE., 
     $                   GLYCT = .TRUE., 
     $                   GLYDT = .TRUE.,
#ifdef USE_INTEGER8
     $                   INTKIND  = SELECTED_INT_KIND(18),
     $                   NODEMEM  =  2500000000_INTKIND,
     $                   INODEMEM =   500000000_INTKIND,
#else
     $                   INTKIND  = SELECTED_INT_KIND(9),
     $                   NODEMEM  =  250000000_INTKIND,
     $                   INODEMEM =   50000000_INTKIND,
#endif
     $                   MAX_CPUS = 384,
     $                   PRNTSIZ  = 10, 
     $                   NOUT     = 6 )
C
C     Parameters for SYCT
C
      LOGICAL         SYCT_CONDEST, SYCT_REDUCE_ONLY, 
     $                SYCT_SOLVE, SYCT_SQUARE_DIM, SYCT_SQUARE_BLK,
     $                SYCT_SQUARE_GRID, SYCT_PIPELINE
      CHARACTER       SYCT_AFORM, SYCT_BFORM, SYCT_CFORM, SYCT_ADIAG, 
     $                SYCT_BDIAG, SYCT_CDIAG, SYCT_TRANSA, SYCT_TRANSB, 
     $                SYCT_ASCHUR, SYCT_BSCHUR, SYCT_COMM
      INTEGER         SYCT_A_NQTRBL, SYCT_B_NQTRBL, SYCT_ISGN,
     $                SYCT_M_MIN, SYCT_M_MAX, SYCT_N_MIN, SYCT_N_MAX, 
     $                SYCT_MB_MIN, SYCT_MB_MAX, SYCT_NB_MIN, 
     $                SYCT_NB_MAX, SYCT_IA_MIN, SYCT_IA_MAX, 
     $                SYCT_IB_MIN, SYCT_IB_MAX, SYCT_NPROW_MIN, 
     $                SYCT_NPROW_MAX, SYCT_NPCOL_MIN, SYCT_NPCOL_MAX, 
     $                SYCT_M_STEP, SYCT_N_STEP, SYCT_MB_STEP, 
     $                SYCT_NB_STEP, SYCT_IA_STEP, SYCT_IB_STEP,
     $                SYCT_NPROW_STEP, SYCT_NPCOL_STEP,
     $                SYCT_REPEAT, SYCT_UPLOSIGN_MIN, SYCT_MB2_MIN, 
     $                SYCT_MB2_MAX, SYCT_NB2_MIN, SYCT_NB2_MAX,
     $                SYCT_MB2_STEP, SYCT_NB2_STEP,
     $                SYCT_UPLOSIGN_MAX
      PARAMETER       ( SYCT_CONDEST = .TRUE., 
     $                  SYCT_REDUCE_ONLY = .FALSE., 
     $                  SYCT_SOLVE = .TRUE.,
     $                  SYCT_PIPELINE = .TRUE.,
     $                  SYCT_SQUARE_DIM = .TRUE., 
     $                  SYCT_SQUARE_BLK = .TRUE.,
     $                  SYCT_SQUARE_GRID = .TRUE.,
     $                  SYCT_A_NQTRBL = 77,
     $                  SYCT_B_NQTRBL = 77,
     $                  SYCT_COMM = 'Demand',
     $                  SYCT_AFORM = 'Not Upper', 
     $                  SYCT_BFORM = 'Not Upper', 
     $                  SYCT_CFORM = 'Random',
     $                  SYCT_ADIAG = 'Diagonal attached', 
     $                  SYCT_BDIAG = 'Diagonal attached', 
     $                  SYCT_CDIAG = 'No diagonal attached', 
     $                  SYCT_UPLOSIGN_MIN = 0, 
     $                  SYCT_UPLOSIGN_MAX = 0, 
     $                  SYCT_ASCHUR = 'Not Schur form', 
     $                  SYCT_BSCHUR = 'Not Schur form',
     $                  SYCT_M_MIN = 2048,
     $                  SYCT_M_MAX = 2048,
     $                  SYCT_M_STEP = 1,
     $                  SYCT_N_MIN = 2048,
     $                  SYCT_N_MAX = 2048,
     $                  SYCT_N_STEP = 1,
     $                  SYCT_MB_MIN = 64, 
     $                  SYCT_MB_MAX = 64,
     $                  SYCT_MB_STEP = 1,
     $                  SYCT_NB_MIN = 64, 
     $                  SYCT_NB_MAX = 64,
     $                  SYCT_NB_STEP = 1,
     $                  SYCT_MB2_MIN = 64, 
     $                  SYCT_MB2_MAX = 64,
     $                  SYCT_MB2_STEP = 1,
     $                  SYCT_NB2_MIN = 6, 
     $                  SYCT_NB2_MAX = 6,
     $                  SYCT_NB2_STEP = 1,
     $                  SYCT_IA_MIN = 1, 
     $                  SYCT_IA_MAX = 1,
     $                  SYCT_IA_STEP = 1,
     $                  SYCT_IB_MIN = 1, 
     $                  SYCT_IB_MAX = 1,
     $                  SYCT_IB_STEP = 1,
     $                  SYCT_NPROW_MIN = 1, 
     $                  SYCT_NPROW_MAX = 2,
     $                  SYCT_NPROW_STEP = 1,
     $                  SYCT_NPCOL_MIN = 1,
     $                  SYCT_NPCOL_MAX = 2,
     $                  SYCT_NPCOL_STEP = 1,
     $                  SYCT_REPEAT = 1 )
C
C     Parameters for LYCT
C
      LOGICAL         LYCT_CONDEST, LYCT_REDUCE_ONLY, LYCT_PIPELINE,
     $                LYCT_SOLVE, LYCT_SQUARE_GRID, LYCT_SQUARE_BLK
      CHARACTER       LYCT_AFORM, LYCT_ADIAG, LYCT_CDIAG, 
     $                LYCT_TRANSA, LYCT_ASCHUR,
     $                LYCT_C_SYMM, LYCT_COMM
      INTEGER         LYCT_A_NQTRBL, 
     $                LYCT_N_MIN, LYCT_N_MAX, 
     $                LYCT_NB_MIN, LYCT_NB_MAX, 
     $                LYCT_IA_MIN, LYCT_IA_MAX, 
     $                LYCT_NPROW_MIN, LYCT_NPROW_MAX, 
     $                LYCT_NPCOL_MIN, LYCT_NPCOL_MAX,
     $                LYCT_N_STEP, LYCT_NB_STEP,
     $                LYCT_IA_STEP,
     $                LYCT_NPROW_STEP, LYCT_NPCOL_STEP,
     $                LYCT_REPEAT, LYCT_UPLOSIGN_MIN,
     $                LYCT_UPLOSIGN_MAX, LYCT_NB2_MIN, 
     $                LYCT_NB2_MAX, LYCT_NB2_STEP
      PARAMETER       ( LYCT_CONDEST = .TRUE., 
     $                  LYCT_REDUCE_ONLY = .FALSE., 
     $                  LYCT_SOLVE = .TRUE.,
     $                  LYCT_PIPELINE = .TRUE.,
     $                  LYCT_SQUARE_GRID = .TRUE.,
     $                  LYCT_SQUARE_BLK = .TRUE.,
     $                  LYCT_AFORM = 'Not Upper',  
     $                  LYCT_C_SYMM = 'Symmetric',
     $                  LYCT_ADIAG = 'Diagonal attached',  
     $                  LYCT_CDIAG = 'No diagonal attached', 
     $                  LYCT_UPLOSIGN_MIN = 0,
     $                  LYCT_UPLOSIGN_MAX = 0,
     $                  LYCT_ASCHUR = 'Not Schur form', 
     $                  LYCT_A_NQTRBL = 77,
     $                  LYCT_N_MIN = 2048, 
     $                  LYCT_N_MAX = 2048,
     $                  LYCT_N_STEP = 1,
     $                  LYCT_NB_MIN = 64, 
     $                  LYCT_NB_MAX = 64,
     $                  LYCT_NB_STEP = 1,
     $                  LYCT_NB2_MIN = 32, 
     $                  LYCT_NB2_MAX = 32,
     $                  LYCT_NB2_STEP = 1,
     $                  LYCT_IA_MIN = 1, 
     $                  LYCT_IA_MAX = 1, 
     $                  LYCT_IA_STEP = 1,
     $                  LYCT_NPROW_MIN = 1, 
     $                  LYCT_NPROW_MAX = 2,
     $                  LYCT_NPROW_STEP = 1,
     $                  LYCT_NPCOL_MIN = 1,
     $                  LYCT_NPCOL_MAX = 2,
     $                  LYCT_NPCOL_STEP = 1,
     $                  LYCT_REPEAT = 1 )
C
C     Parameters for SYDT
C
      LOGICAL         SYDT_CONDEST, SYDT_SHIFTS, 
     $                SYDT_ON_DEMAND, SYDT_REDUCE_ONLY, 
     $                SYDT_SOLVE, SYDT_SQUARE_DIM, SYDT_SQUARE_BLK,
     $                SYDT_SQUARE_GRID, SYDT_PIPELINE
      CHARACTER       SYDT_AFORM, SYDT_BFORM, SYDT_CFORM, SYDT_ADIAG, 
     $                SYDT_BDIAG, SYDT_CDIAG, SYDT_TRANSA, SYDT_TRANSB, 
     $                SYDT_ASCHUR, SYDT_BSCHUR, SYDT_COMM
      INTEGER         SYDT_A_NQTRBL, SYDT_B_NQTRBL, SYDT_ISGN,
     $                SYDT_M_MIN, SYDT_M_MAX, SYDT_N_MIN, SYDT_N_MAX, 
     $                SYDT_MB_MIN, SYDT_MB_MAX, SYDT_NB_MIN, 
     $                SYDT_NB_MAX, SYDT_IA_MIN, SYDT_IA_MAX, 
     $                SYDT_IB_MIN, SYDT_IB_MAX, SYDT_NPROW_MIN, 
     $                SYDT_NPROW_MAX, SYDT_NPCOL_MIN, SYDT_NPCOL_MAX, 
     $                SYDT_M_STEP, SYDT_N_STEP, SYDT_MB_STEP, 
     $                SYDT_NB_STEP, SYDT_IA_STEP, SYDT_IB_STEP,
     $                SYDT_NPROW_STEP, SYDT_NPCOL_STEP,
     $                SYDT_REPEAT, SYDT_UPLOSIGN_MIN, SYDT_UPLOSIGN_MAX,
     $                SYDT_MB2_MIN, SYDT_MB2_MAX, SYDT_MB2_STEP
      PARAMETER       ( SYDT_CONDEST = .TRUE., 
     $                  SYDT_REDUCE_ONLY = .FALSE., 
     $                  SYDT_SOLVE = .TRUE.,
     $                  SYDT_PIPELINE = .FALSE.,
     $                  SYDT_SQUARE_DIM = .TRUE., 
     $                  SYDT_SQUARE_BLK = .TRUE.,
     $                  SYDT_SQUARE_GRID = .TRUE.,
     $                  SYDT_A_NQTRBL = 77,
     $                  SYDT_B_NQTRBL = 77,
     $                  SYDT_COMM = 'Demand',
     $                  SYDT_AFORM = 'Not Upper', 
     $                  SYDT_BFORM = 'Not Upper', 
     $                  SYDT_CFORM = 'Random',
     $                  SYDT_ADIAG = 'Diagonal attached', 
     $                  SYDT_BDIAG = 'Diagonal attached', 
     $                  SYDT_CDIAG = 'No diagonal attached', 
     $                  SYDT_UPLOSIGN_MIN = 3,
     $                  SYDT_UPLOSIGN_MAX = 3,
     $                  SYDT_ASCHUR = 'Not Schur form', 
     $                  SYDT_BSCHUR = 'Not Schur form',
     $                  SYDT_M_MIN = 2048,
     $                  SYDT_M_MAX = 2048,
     $                  SYDT_M_STEP = 1,
     $                  SYDT_N_MIN = 2048,
     $                  SYDT_N_MAX = 2048,
     $                  SYDT_N_STEP = 1,
     $                  SYDT_MB_MIN = 64, 
     $                  SYDT_MB_MAX = 64,
     $                  SYDT_MB_STEP = 1,
     $                  SYDT_NB_MIN = 64, 
     $                  SYDT_NB_MAX = 64,
     $                  SYDT_NB_STEP = 1,
     $                  SYDT_MB2_MIN = 64, 
     $                  SYDT_MB2_MAX = 64,
     $                  SYDT_MB2_STEP = 1,
     $                  SYDT_IA_MIN = 1, 
     $                  SYDT_IA_MAX = 1,
     $                  SYDT_IA_STEP = 1,
     $                  SYDT_IB_MIN = 1, 
     $                  SYDT_IB_MAX = 1,
     $                  SYDT_IB_STEP = 1,
     $                  SYDT_NPROW_MIN = 1, 
     $                  SYDT_NPROW_MAX = 2,
     $                  SYDT_NPROW_STEP = 1,
     $                  SYDT_NPCOL_MIN = 1,
     $                  SYDT_NPCOL_MAX = 2,
     $                  SYDT_NPCOL_STEP = 1,
     $                  SYDT_REPEAT = 1 )
C
C     Parameters for LYDT
C
      LOGICAL         LYDT_CONDEST, LYDT_REDUCE_ONLY, 
     $                LYDT_SOLVE, LYDT_SQUARE_GRID
      CHARACTER       LYDT_AFORM, LYDT_ADIAG, LYDT_CDIAG, 
     $                LYDT_TRANSA, LYDT_ASCHUR,
     $                LYDT_C_SYMM, LYDT_COMM
      INTEGER         LYDT_A_NQTRBL, 
     $                LYDT_N_MIN, LYDT_N_MAX, 
     $                LYDT_NB_MIN, LYDT_NB_MAX, 
     $                LYDT_IA_MIN, LYDT_IA_MAX, 
     $                LYDT_NPROW_MIN, LYDT_NPROW_MAX, 
     $                LYDT_NPCOL_MIN, LYDT_NPCOL_MAX,
     $                LYDT_N_STEP, LYDT_NB_STEP,
     $                LYDT_IA_STEP,
     $                LYDT_NPROW_STEP, LYDT_NPCOL_STEP,
     $                LYDT_REPEAT, LYDT_UPLOSIGN_MIN,
     $                LYDT_UPLOSIGN_MAX, LYDT_NB2_MIN, 
     $                LYDT_NB2_MAX, LYDT_NB2_STEP
      PARAMETER       ( LYDT_CONDEST = .TRUE.,
     $                  LYDT_REDUCE_ONLY = .FALSE., 
     $                  LYDT_SOLVE = .TRUE.,
     $                  LYDT_SQUARE_GRID = .TRUE.,
     $                  LYDT_AFORM = 'Not Upper',  
     $                  LYDT_C_SYMM = 'Symmetric',
     $                  LYDT_ADIAG = 'Diagonal attached',  
     $                  LYDT_CDIAG = 'No diagonal attached', 
     $                  LYDT_UPLOSIGN_MIN = 0,
     $                  LYDT_UPLOSIGN_MAX = 0,
     $                  LYDT_ASCHUR = 'Not Schur form', 
     $                  LYDT_A_NQTRBL = 77,
     $                  LYDT_N_MIN = 2048, 
     $                  LYDT_N_MAX = 2048,
     $                  LYDT_N_STEP = 1,
     $                  LYDT_NB_MIN = 64, 
     $                  LYDT_NB_MAX = 64,
     $                  LYDT_NB_STEP = 1, 
     $                  LYDT_NB2_MIN = 64, 
     $                  LYDT_NB2_MAX = 64,
     $                  LYDT_NB2_STEP = 1,
     $                  LYDT_IA_MIN = 1, 
     $                  LYDT_IA_MAX = 1, 
     $                  LYDT_IA_STEP = 1,
     $                  LYDT_NPROW_MIN = 1, 
     $                  LYDT_NPROW_MAX = 2,
     $                  LYDT_NPROW_STEP = 1,
     $                  LYDT_NPCOL_MIN = 1,
     $                  LYDT_NPCOL_MAX = 2,
     $                  LYDT_NPCOL_STEP = 1,
     $                  LYDT_REPEAT = 1 )
C
C     Parameters for GCSY
C
      LOGICAL         GCSY_CONDEST, GCSY_SHIFTS, 
     $                GCSY_ON_DEMAND, GCSY_REDUCE_ONLY, 
     $                GCSY_SOLVE, GCSY_SQUARE_DIM, GCSY_SQUARE_BLK,
     $                GCSY_SQUARE_GRID, GCSY_PIPELINE
      CHARACTER       GCSY_ADFORM, GCSY_BEFORM, GCSY_CFFORM, 
     $                GCSY_ADDIAG, GCSY_BEDIAG, GCSY_CFDIAG, 
     $                GCSY_TRANSAD, GCSY_TRANSBE, GCSY_ADSCHUR, 
     $                GCSY_BESCHUR, GCSY_TRANSZ, GCSY_COMM
      INTEGER         GCSY_AD_NQTRBL, GCSY_BE_NQTRBL, GCSY_ISGN,
     $                GCSY_M_MIN, GCSY_M_MAX, GCSY_N_MIN, GCSY_N_MAX, 
     $                GCSY_MB_MIN, GCSY_MB_MAX, GCSY_NB_MIN, 
     $                GCSY_NB_MAX, GCSY_IAD_MIN, GCSY_IAD_MAX, 
     $                GCSY_IBE_MIN, GCSY_IBE_MAX,  GCSY_NPROW_MIN, 
     $                GCSY_NPROW_MAX, GCSY_NPCOL_MIN, GCSY_NPCOL_MAX, 
     $                GCSY_M_STEP, GCSY_N_STEP, GCSY_MB_STEP, 
     $                GCSY_NB_STEP, GCSY_IAD_STEP, GCSY_IBE_STEP,
     $                GCSY_NPROW_STEP, GCSY_NPCOL_STEP,
     $                GCSY_REPEAT, GCSY_UPLOSIGN_MIN, GCSY_UPLOSIGN_MAX,
     $                GCSY_MB2_MIN, GCSY_MB2_MAX, GCSY_MB2_STEP,
     $                GCSY_NB2_MIN, GCSY_NB2_MAX, GCSY_NB2_STEP
      PARAMETER       ( GCSY_CONDEST = .TRUE., 
     $                  GCSY_REDUCE_ONLY = .FALSE., 
     $                  GCSY_SOLVE = .TRUE.,
     $                  GCSY_PIPELINE = .TRUE.,
     $                  GCSY_SQUARE_GRID = .TRUE.,
     $                  GCSY_SQUARE_DIM = .TRUE.,
     $                  GCSY_SQUARE_BLK = .TRUE.,
     $                  GCSY_ADFORM = 'Upper', 
     $                  GCSY_BEFORM = 'Upper', 
     $                  GCSY_CFFORM = 'Random',
     $                  GCSY_ADDIAG = 'Diagonal attached', 
     $                  GCSY_BEDIAG = 'Diagonal attached', 
     $                  GCSY_CFDIAG = 'No diagonal attached', 
     $                  GCSY_UPLOSIGN_MIN = 0, 
     $                  GCSY_UPLOSIGN_MAX = 0, 
     $                  GCSY_ADSCHUR = 'Schur form', 
     $                  GCSY_BESCHUR = 'Schur form',
     $                  GCSY_TRANSZ = 'Not Transposed',
     $                  GCSY_COMM = 'Demand',
     $                  GCSY_AD_NQTRBL = 77,
     $                  GCSY_BE_NQTRBL = 77,
     $                  GCSY_M_MIN = 2048, 
     $                  GCSY_M_MAX = 2048,
     $                  GCSY_M_STEP = 1,
     $                  GCSY_N_MIN = 2048, 
     $                  GCSY_N_MAX = 2048,
     $                  GCSY_N_STEP = 1,
     $                  GCSY_MB_MIN = 128, 
     $                  GCSY_MB_MAX = 128,
     $                  GCSY_MB_STEP = 1,
     $                  GCSY_NB_MIN = 128, 
     $                  GCSY_NB_MAX = 128,
     $                  GCSY_NB_STEP = 1,
     $                  GCSY_MB2_MIN = 64, 
     $                  GCSY_MB2_MAX = 64,
     $                  GCSY_MB2_STEP = 1,
     $                  GCSY_NB2_MIN = 64, 
     $                  GCSY_NB2_MAX = 64,
     $                  GCSY_NB2_STEP = 1,
     $                  GCSY_IAD_MIN = 1, 
     $                  GCSY_IAD_MAX = 1,
     $                  GCSY_IAD_STEP = 1,
     $                  GCSY_IBE_MIN = 1, 
     $                  GCSY_IBE_MAX = 1,
     $                  GCSY_IBE_STEP = 1,
     $                  GCSY_NPROW_MIN = 1, 
     $                  GCSY_NPROW_MAX = 2,
     $                  GCSY_NPROW_STEP = 1,
     $                  GCSY_NPCOL_MIN = 1,
     $                  GCSY_NPCOL_MAX = 2,
     $                  GCSY_NPCOL_STEP = 1,
     $                  GCSY_REPEAT = 1 )
C
C     Parameters for GSYL
C
      LOGICAL         GSYL_CONDEST, GSYL_SHIFTS, 
     $                GSYL_ON_DEMAND, GSYL_REDUCE_ONLY, 
     $                GSYL_SOLVE, GSYL_SQUARE_DIM, GSYL_SQUARE_BLK,
     $                GSYL_SQUARE_GRID, GSYL_PIPELINE
      CHARACTER       GSYL_ACFORM, GSYL_BDFORM, GSYL_EFORM, GSYL_ACDIAG, 
     $                GSYL_BDDIAG, GSYL_EDIAG, GSYL_TRANSAC, 
     $                GSYL_TRANSBD, GSYL_ACSCHUR, GSYL_BDSCHUR,
     $                GSYL_COMM 
      INTEGER         GSYL_AC_NQTRBL, GSYL_BD_NQTRBL, GSYL_ISGN,
     $                GSYL_M_MIN, GSYL_M_MAX, GSYL_N_MIN, GSYL_N_MAX, 
     $                GSYL_MB_MIN, GSYL_MB_MAX, GSYL_NB_MIN, 
     $                GSYL_NB_MAX, GSYL_IAC_MIN, GSYL_IAC_MAX, 
     $                GSYL_IBD_MIN, GSYL_IBD_MAX, GSYL_NPROW_MIN, 
     $                GSYL_NPROW_MAX, GSYL_NPCOL_MIN, GSYL_NPCOL_MAX, 
     $                GSYL_M_STEP, GSYL_N_STEP, GSYL_MB_STEP, 
     $                GSYL_NB_STEP, GSYL_IAC_STEP, GSYL_IBD_STEP,
     $                GSYL_NPROW_STEP, GSYL_NPCOL_STEP,
     $                GSYL_REPEAT, GSYL_UPLOSIGN_MIN, GSYL_UPLOSIGN_MAX,
     $                GSYL_MB2_MIN, GSYL_MB2_MAX, GSYL_MB2_STEP
      PARAMETER       ( GSYL_CONDEST = .TRUE., 
     $                  GSYL_REDUCE_ONLY = .FALSE., 
     $                  GSYL_SOLVE = .TRUE.,
     $                  GSYL_PIPELINE = .TRUE.,
     $                  GSYL_SQUARE_DIM = .TRUE., 
     $                  GSYL_SQUARE_BLK = .TRUE.,
     $                  GSYL_SQUARE_GRID = .TRUE.,
     $                  GSYL_ACFORM = 'Upper', 
     $                  GSYL_BDFORM = 'Upper', 
     $                  GSYL_EFORM = 'Random',
     $                  GSYL_ACDIAG = 'Diagonal attached', 
     $                  GSYL_BDDIAG = 'Diagonal attached', 
     $                  GSYL_EDIAG = 'No diagonal attached', 
     $                  GSYL_UPLOSIGN_MIN = 3, 
     $                  GSYL_UPLOSIGN_MAX = 3, 
     $                  GSYL_ACSCHUR = 'Schur form', 
     $                  GSYL_BDSCHUR = 'Schur form',
     $                  GSYL_COMM = 'Demand',
     $                  GSYL_AC_NQTRBL = 77,
     $                  GSYL_BD_NQTRBL = 77,
     $                  GSYL_M_MIN = 2048, 
     $                  GSYL_M_MAX = 2048,
     $                  GSYL_M_STEP = 1,
     $                  GSYL_N_MIN = 2048, 
     $                  GSYL_N_MAX = 2048,
     $                  GSYL_N_STEP = 1,
     $                  GSYL_MB_MIN = 128, 
     $                  GSYL_MB_MAX = 128,
     $                  GSYL_MB_STEP = 1,
     $                  GSYL_NB_MIN = 128, 
     $                  GSYL_NB_MAX = 128,
     $                  GSYL_NB_STEP = 1, 
     $                  GSYL_MB2_MIN = 128, 
     $                  GSYL_MB2_MAX = 128,
     $                  GSYL_MB2_STEP = 1,
     $                  GSYL_IAC_MIN = 1, 
     $                  GSYL_IAC_MAX = 1,
     $                  GSYL_IAC_STEP = 1,
     $                  GSYL_IBD_MIN = 1, 
     $                  GSYL_IBD_MAX = 1,
     $                  GSYL_IBD_STEP = 1,
     $                  GSYL_NPROW_MIN = 1, 
     $                  GSYL_NPROW_MAX = 2,
     $                  GSYL_NPROW_STEP = 1,
     $                  GSYL_NPCOL_MIN = 1,
     $                  GSYL_NPCOL_MAX = 2,
     $                  GSYL_NPCOL_STEP = 1,
     $                  GSYL_REPEAT = 1 )
C
C     Parameters for GLYCT
C
      LOGICAL         GLYCT_CONDEST, GLYCT_REDUCE_ONLY, 
     $                GLYCT_SOLVE, GLYCT_SQUARE_GRID
      CHARACTER       GLYCT_AEFORM, GLYCT_AEDIAG, 
     $                GLYCT_CDIAG, GLYCT_TRANSAE, GLYCT_AESCHUR,
     $                GLYCT_C_SYMM, GLYCT_COMM
      INTEGER         GLYCT_AE_NQTRBL, GLYCT_N_MIN, GLYCT_N_MAX, 
     $                GLYCT_NB_MIN, GLYCT_NB_MAX, 
     $                GLYCT_IAE_MIN, GLYCT_IAE_MAX, 
     $                GLYCT_NPROW_MIN, GLYCT_NPROW_MAX, 
     $                GLYCT_NPCOL_MIN, GLYCT_NPCOL_MAX,
     $                GLYCT_N_STEP, GLYCT_NB_STEP,
     $                GLYCT_IAE_STEP, GLYCT_NPROW_STEP, 
     $                GLYCT_NPCOL_STEP, GLYCT_REPEAT, 
     $                GLYCT_UPLOSIGN_MIN, GLYCT_UPLOSIGN_MAX,
     $                GLYCT_NB2_MIN, GLYCT_NB2_MAX, GLYCT_NB2_STEP
      PARAMETER       ( GLYCT_CONDEST = .TRUE., 
     $                  GLYCT_REDUCE_ONLY = .FALSE., 
     $                  GLYCT_SOLVE = .TRUE.,
     $                  GLYCT_SQUARE_GRID = .TRUE.,
     $                  GLYCT_AEFORM = 'Upper',  
     $                  GLYCT_AEDIAG = 'Diagonal attached',  
     $                  GLYCT_CDIAG = 'No diagonal attached', 
     $                  GLYCT_UPLOSIGN_MIN = 0,
     $                  GLYCT_UPLOSIGN_MAX = 0,
     $                  GLYCT_AESCHUR = 'Schur form',
     $                  GLYCT_C_SYMM = 'Symmetric',
     $                  GLYCT_AE_NQTRBL = 77,
     $                  GLYCT_N_MIN = 2048, 
     $                  GLYCT_N_MAX = 2048,
     $                  GLYCT_N_STEP = 1,
     $                  GLYCT_NB_MIN = 128, 
     $                  GLYCT_NB_MAX = 128,
     $                  GLYCT_NB_STEP = 1,
     $                  GLYCT_NB2_MIN = 128, 
     $                  GLYCT_NB2_MAX = 128,
     $                  GLYCT_NB2_STEP = 1,
     $                  GLYCT_IAE_MIN = 1, 
     $                  GLYCT_IAE_MAX = 1,
     $                  GLYCT_IAE_STEP = 1,
     $                  GLYCT_NPROW_MIN = 1, 
     $                  GLYCT_NPROW_MAX = 2,
     $                  GLYCT_NPROW_STEP = 1,
     $                  GLYCT_NPCOL_MIN = 1,
     $                  GLYCT_NPCOL_MAX = 2, 
     $                  GLYCT_NPCOL_STEP = 1,
     $                  GLYCT_REPEAT = 1 )
C
C     Parameters for GLYDT
C
      LOGICAL         GLYDT_CONDEST, GLYDT_REDUCE_ONLY, 
     $                GLYDT_SOLVE, GLYDT_SQUARE_GRID
      CHARACTER       GLYDT_AEFORM, GLYDT_AEDIAG, 
     $                GLYDT_CDIAG, GLYDT_TRANSAE, GLYDT_AESCHUR,
     $                GLYDT_C_SYMM, GLYDT_COMM
      INTEGER         GLYDT_AE_NQTRBL, GLYDT_N_MIN, GLYDT_N_MAX, 
     $                GLYDT_NB_MIN, GLYDT_NB_MAX, 
     $                GLYDT_IAE_MIN, GLYDT_IAE_MAX, 
     $                GLYDT_NPROW_MIN, GLYDT_NPROW_MAX, 
     $                GLYDT_NPCOL_MIN, GLYDT_NPCOL_MAX,
     $                GLYDT_N_STEP, GLYDT_NB_STEP,
     $                GLYDT_IAE_STEP, GLYDT_NPROW_STEP, 
     $                GLYDT_NPCOL_STEP, GLYDT_REPEAT,
     $                GLYDT_UPLOSIGN_MIN, GLYDT_UPLOSIGN_MAX,
     $                GLYDT_NB2_MIN, GLYDT_NB2_MAX, GLYDT_NB2_STEP
      PARAMETER       ( GLYDT_CONDEST = .TRUE., 
     $                  GLYDT_REDUCE_ONLY = .FALSE., 
     $                  GLYDT_SOLVE = .TRUE.,
     $                  GLYDT_SQUARE_GRID = .TRUE.,
     $                  GLYDT_AEFORM = 'Upper',  
     $                  GLYDT_AEDIAG = 'Diagonal attached',  
     $                  GLYDT_CDIAG = 'No diagonal attached', 
     $                  GLYDT_UPLOSIGN_MIN = 0,
     $                  GLYDT_UPLOSIGN_MAX = 0,  
     $                  GLYDT_AESCHUR = 'Schur form', 
     $                  GLYDT_C_SYMM = 'Symmetric',
     $                  GLYDT_AE_NQTRBL = 77,
     $                  GLYDT_N_MIN = 2048, 
     $                  GLYDT_N_MAX = 2048, 
     $                  GLYDT_N_STEP = 1,
     $                  GLYDT_NB_MIN = 128, 
     $                  GLYDT_NB_MAX = 128,
     $                  GLYDT_NB_STEP = 1, 
     $                  GLYDT_NB2_MIN = 128, 
     $                  GLYDT_NB2_MAX = 128,
     $                  GLYDT_NB2_STEP = 1,
     $                  GLYDT_IAE_MIN = 1, 
     $                  GLYDT_IAE_MAX = 1,
     $                  GLYDT_IAE_STEP = 1,
     $                  GLYDT_NPROW_MIN = 1, 
     $                  GLYDT_NPROW_MAX = 2,
     $                  GLYDT_NPROW_STEP = 1,
     $                  GLYDT_NPCOL_MIN = 1,
     $                  GLYDT_NPCOL_MAX = 2, 
     $                  GLYDT_NPCOL_STEP = 1,
     $                  GLYDT_REPEAT = 1 )
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === THE FOLLOWING PART MAY NOT BE MODIFIED BY THE USER! ===
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Internal declarations
C
C     .. Parameters ..
      INTEGER           BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                  LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER         ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                    CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                    RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
#ifdef USE_INTEGER8
      INTEGER*8         DBLESZ, INTGSZ, MEMSIZ, IMEMSIZ
#else
      INTEGER           DBLESZ, INTGSZ, MEMSIZ, IMEMSIZ
#endif
      INTEGER           ERRORNUM, IZERO
      DOUBLE PRECISION  ZERO, ONE, TEN
      PARAMETER         ( DBLESZ = 8, INTGSZ = 4, 
     $                    MEMSIZ = NODEMEM / DBLESZ,
     $                    IMEMSIZ = INODEMEM / INTGSZ,
     $                    ERRORNUM = 1, IZERO = 0, ZERO = 0.0D+00,
     $                    ONE = 1.0D+00, TEN = 10.0D+00 )
C    ..
C    .. Local Scalars ..
      CHARACTER         LYCT_C_SYMM2, LYDT_C_SYMM2, GLYCT_C_SYMM2, 
     $                  GLYDT_C_SYMM2, NOP, GLYCT_ITRANSAE,
     $                  GLYDT_ITRANSAE
      INTEGER           IAM, NPROCS, ICTXT, SYS_NPROW, SYS_NPCOL, M, N, 
     $                  MB, NB, IA, JA, IB, JB, IC, JC, NPROW, NPCOL, 
     $                  MYROW, MYCOL, TEMP_ICTXT, AROWS, ACOLS, BROWS, 
     $                  BCOLS, INFO, IA_MIN2, JA_MIN2, IA_MAX2, JA_MAX2, 
     $                  IB_MIN2, JB_MIN2, IB_MAX2, JB_MAX2, IC_MIN2, 
     $                  JC_MIN2, IC_MAX2, JC_MAX2, IPA, IPACPY, IPB, 
     $                  IPBCPY, IPC, IPCCPY, IPX, DDIAG, SDIAG, IPW, 
     $                  WORKSIZ, NOITER, NOEXTSYS, I, ID, JD, IE, JE, 
     $                  IF, JF, REPEAT, ADROWS, ADCOLS, BEROWS, BECOLS, 
     $                  IAD_MIN2, JAD_MIN2, IBE_MIN2, JBE_MIN2, JAD, 
     $                  IAD, IBE, JBE, ICF_MIN2, JCF_MIN2, ICF, JCF, 
     $                  IAD_MAX2, JAD_MAX2, IBE_MAX2, JBE_MAX2,
     $                  ICF_MAX2, JCF_MAX2, DIAG1, DIAG2, TEMP_NPROW,
     $                  IPD, IPDCPY, IPE, IPECPY, IPF, IPFCPY,
     $                  IPY, ACROWS, ACCOLS, BDROWS, BDCOLS,
     $                  IAC_MIN2, IAC, JAC_MIN2, JAC, IBD, IBD_MIN2,
     $                  JBD, JDB_MIN2, IE_MIN2, JE_MIN2, TEMP_NPCOL,
     $                  IAC_MAX2, JAC_MAX2, IBD_MAX2, JBD_MAX2,
     $                  IE_MAX2, JE_MAX2, JBD_MIN2, AEROWS, AECOLS,
     $                  IAE_MIN2, IAE, JAE_MIN2, JAE, IAE_MAX2,
     $                  JAE_MAX2, GLYCT_AENQTRBL, GLYDT_AENQTRBL,
     $                  BLKS2B2, UPLOSIGN, MB2, NB2, THREADS, K
      DOUBLE PRECISION  TOTTIME, T1, T2, MYTIME, EPS, CNORM, XNORM,
     $                  SCALE, ANORM, BNORM, ABSRESID, RESID, ABSERROR,
     $                  RELERROR, ESTTIME, EST, ADNORM, ENORM, DNORM,
     $                  CFNORM, XYNORM, YNORM, FNORM, BENORM, R1NORM,
     $                  R2NORM, RNORM, FORWARD, RELFORW, ACNORM, 
     $                  BDNORM, AENORM, TPROF1, TPROF2, TPROF3,
     $                  CPROF1, CPROF2, CPROF3, ALPHA, GFLOPS
C    ..
C    .. Local Arrays ..
      INTEGER           DESCA( DLEN_ ), DESCB( DLEN_ ), DESCC( DLEN_ ),
     $                  DESCD( DLEN_ ), DESCE( DLEN_ ), DESCF( DLEN_ ),
     $                  MBNB2( 2 )
#ifdef USE_DYNAMIC
      INTEGER, ALLOCATABLE :: IMEM(:)
      DOUBLE PRECISION, ALLOCATABLE :: MEM(:)
#else
      INTEGER IMEM( IMEMSIZ )
      DOUBLE PRECISION MEM( MEMSIZ )
#endif
C    ..
C    .. External Subroutines ..
      EXTERNAL          BLACS_PINFO, BLACS_GET, BLACS_ABORT,
     $                  BLACS_GRIDEXIT, BLACS_EXIT
C    ..
C    .. External Functions ..
      LOGICAL           LSAME
      INTEGER           NUMROC, OMP_GET_NUM_THREADS
      DOUBLE PRECISION  MPI_WTIME, PDLANGE, PDLAMCH, DLAPY2
      EXTERNAL          LSAME, MPI_WTIME, NUMROC, OMP_GET_NUM_THREADS
C    ..
C    .. Intrinsic Functions ..
C    ..
C    .. Executable Statements ..
C
C     Get process id and number of processes - set up underlying
C     mesh in a horizontal ring. The processor configurations to
C     use during the tests are set up later.
C
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL BLACS_GET( 0, 0, ICTXT )
#ifdef LOOPGRID
      SYS_NPROW = 1
      SYS_NPCOL = NPROCS
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', SYS_NPROW, SYS_NPCOL )
#else
      NPROW = INT( SQRT( REAL(NPROCS) ) )
      NPCOL = NPROCS / NPROW
      IF( NPROW*NPCOL.NE.NPROCS ) THEN
         NPROW = NPROCS
         NPCOL = 1
      END IF
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      TEMP_ICTXT = ICTXT
#endif
C
C     Start global timer for testprogram
C
      TOTTIME = MPI_WTIME()
      CALL CORE_UNLIMITED()
C
C     Print welcome message 
C
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT , FMT = *)
         WRITE( NOUT , FMT = *) 
     $        '========================================================'
         WRITE( NOUT , FMT = *)
     $        '   TESTSCASY PROGRAM: SCASY version 1.0, July 31 2009   '
         WRITE( NOUT , FMT = *)
     $        '========================================================'  
         WRITE( NOUT, FMT = *)        
         WRITE( NOUT, FMT = *)'The outputs from the tests are:'
         WRITE( NOUT, FMT = *)'ULS  - The UPLOSIGN variant (see RECSY)'
         WRITE( NOUT, FMT = *)'M    - The M matrix dimension'
         WRITE( NOUT, FMT = *)'N    - The N matrix dimension'
         WRITE( NOUT, FMT = *)'MB   - The blocking factor for M'  
         WRITE( NOUT, FMT = *)'NB   - The blocking factor for N'
         WRITE( NOUT, FMT = *)'MB2  - The pipelining block size for M'
         WRITE( NOUT, FMT = *)'NB2  - The pipelining block size for N'
         WRITE( NOUT, FMT = *)'I[X] - Row start index for matrix sub(X)'
         WRITE( NOUT, FMT = *)'I[X] - Col start index for matrix sub(X)'
         WRITE( NOUT, FMT = *)'P_r  - # processor (node) rows'
         WRITE( NOUT, FMT = *)'P_c  - # processor (node) columns'
         WRITE( NOUT, FMT = *)'Thr  - # OMP threads on each node'
         WRITE( NOUT, FMT = *)'Tme  - Total runtime in seconds'
         WRITE( NOUT, FMT = *)'Tp1  - Time for B-S''s method: step 1'
         WRITE( NOUT, FMT = *)'Tp2  - Time for B-S''s method: steps 2,4'
         WRITE( NOUT, FMT = *)'Tp3  - Time for B-S''s method: step 3'
         WRITE( NOUT, FMT = *)'Gfl  - Glops/sec. for triangular solver'
         WRITE( NOUT, FMT = *)'Sc   - Global scaling factor'
         WRITE( NOUT, FMT = *)'E_a  - Absolute error norm'
         WRITE( NOUT, FMT = *)'E_r  - Relative error norm'
         WRITE( NOUT, FMT = *)'R_a  - Absolute residual norm'
         WRITE( NOUT, FMT = *)'R_r  - Relative residual norm'
         WRITE( NOUT, FMT = *)'est  - Estimate of sep[ACRO]^{-1}'
         WRITE( NOUT, FMT = *)'etm  - Est: Time for computing estimate'
         WRITE( NOUT, FMT = *)'Cp1  - Est: Time in all-2-all broadcasts'
         WRITE( NOUT, FMT = *)'Cp2  - Est: Time in ScaLAPACK''s PDLACON'
         WRITE( NOUT, FMT = *)'Cp3  - Est: Time for triangular solver'
         WRITE( NOUT, FMT = *)'itr  - Est: # iterations in estimation'
         WRITE( NOUT, FMT = *)'INF  - INFO = 0 <=> Test was O.K.'
      END IF
C     
C     Checking parameters for inconstistency - checking if performed
C     only for the chosen matrix equations
C     
      IF( .NOT. (SYCT .OR. LYCT .OR. SYDT .OR. LYDT .OR.
     $     GCSY .OR. GSYL .OR. GLYCT .OR. GLYDT ) ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = * ) ' *** ERROR ***'
            WRITE( NOUT , FMT = *) 'No equation is considered!'
            WRITE( NOUT , FMT = *)
         END IF
         GO TO 9998
      END IF
      IF( NODEMEM .LT. 1 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = * ) ' *** ERROR ***'
            WRITE( NOUT , FMT = *) 'Negative NODEMEM!'
            WRITE( NOUT , FMT = *)
         END IF
         GO TO 9998
      END IF  
      IF( PRNTSIZ .LT. 0 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = * ) ' *** ERROR ***'
            WRITE( NOUT , FMT = *) 'Negative PRNTSIZ!'
            WRITE( NOUT , FMT = *)
         END IF
         GO TO 9998
      END IF
      IF( MAX_CPUS .LT. 1 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = * ) ' *** ERROR ***'
            WRITE( NOUT , FMT = *) 'Negative MAX_CPUS!'
            WRITE( NOUT , FMT = *)
         END IF
         GO TO 9998
      END IF
      IF( NPROCS .GT. MAX_CPUS ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
            WRITE( NOUT , FMT = *) 'MAX_CPUS is to small!'
            WRITE( NOUT , FMT = *)
         END IF
         GO TO 9998
      END IF
      IF( NOUT.NE.6 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *) ' *** Warning! ***'
            WRITE( NOUT , FMT = *) 'NOUT not standard output!'
            WRITE( NOUT , FMT = *)
         END IF
      END IF
C
C     Checking SYCT
C
      IF( SYCT ) THEN
         IF( SYCT_CONDEST .AND. .NOT. ( SYCT_REDUCE_ONLY .OR. 
     $        SYCT_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of SYCT not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. ( LSAME( SYCT_COMM, 'D' ) 
     $        .OR. LSAME( SYCT_COMM, 'S' ) ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Communication pattern in triangular solver must'//
     $              ' be specified for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_UPLOSIGN_MIN.LT.0 .OR. SYCT_UPLOSIGN_MIN.GT.7 .OR.
     $        SYCT_UPLOSIGN_MAX.LT.0 .OR. SYCT_UPLOSIGN_MAX.GT.7 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', SYCT_ASCHUR ) .AND. 
     $        .NOT. LSAME( 'S', SYCT_ASCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of A must be specified for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', SYCT_BSCHUR ) .AND. 
     $        .NOT. LSAME( 'S', SYCT_BSCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of B must be specified for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_A_NQTRBL.LT.0 .OR. SYCT_A_NQTRBL.GT.100 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= A_NQTRBL <= 100 must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_B_NQTRBL.LT.0 .OR. SYCT_B_NQTRBL.GT.100 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= B_NQTRBL <= 100 must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_M_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN >= 1 must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_M_MIN .GT. SYCT_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN <= M_MAX must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_N_MIN .GT. SYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_M_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_MB_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MIN >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_MB_MAX.GT.SYCT_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MAX <= M_MAX  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_NB_MAX.GT.SYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_MB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_MB2_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MIN >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_MB2_MAX.GT.SYCT_MB_MIN .AND. SYCT_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MAX <= MB_MIN  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_NB2_MAX.GT.SYCT_NB_MIN. AND. SYCT_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_MAX <= NB_MIN  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_MB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_NB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_IA_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MIN >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_IA_MAX.LT.SYCT_IA_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX >= IA_MIN  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_IA_MAX.GT.SYCT_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX <= M_MAX  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
          IF( SYCT_IB_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_MIN >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_IB_MAX.LT.SYCT_IB_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_MAX >= IB_MIN  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYCT_IB_MAX.GT.SYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_MAX <= N_MAX  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_IA_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_IB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_NPROW_MAX.LT.SYCT_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_NPCOL_MAX.LT.SYCT_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( SYCT_NPROW_MAX*SYCT_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYCT_NPROW_MAX*SYCT_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for SYCT !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( SYCT_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYCT_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for SYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Checking LYCT
C     
      IF( LYCT ) THEN
         IF( LYCT_CONDEST .AND. .NOT. ( LYCT_REDUCE_ONLY .OR. 
     $        LYCT_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of LYCT not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_UPLOSIGN_MIN.LT.0 .OR. LYCT_UPLOSIGN_MIN.GT.1 .OR.
     $        LYCT_UPLOSIGN_MAX.LT.0 .OR. LYCT_UPLOSIGN_MAX.GT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', LYCT_ASCHUR ) .AND. 
     $        .NOT. LSAME( 'S', LYCT_ASCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of A must be specified for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', LYCT_C_SYMM ) .AND. 
     $        .NOT. LSAME( 'S', LYCT_C_SYMM ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Symmetric mode of C must be specified for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_A_NQTRBL.LT.0 .OR. LYCT_A_NQTRBL.GT.100 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= A_NQTRBL <= 100 must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_N_MIN .GT. LYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_NB_MAX.GT.LYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_NB2_MAX.GT.LYCT_NB_MIN .AND. LYCT_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_MAX <= NB_MIN  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_NB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_STEP >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_IA_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MIN >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_IA_MAX.LT.LYCT_IA_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX >= IA_MIN  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYCT_IA_MAX.GT.LYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX <= N_MAX  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_IA_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_STEP >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_NPROW_MAX.LT.LYCT_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_NPCOL_MAX.LT.LYCT_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( LYCT_NPROW_MAX*LYCT_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'LYCT_NPROW_MAX*LYCT_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for LYCT !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( LYCT_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYCT_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for LYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Checking SYDT
C
      IF( SYDT ) THEN
         IF( SYDT_CONDEST .AND. .NOT. ( SYDT_REDUCE_ONLY .OR. 
     $        SYDT_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of SYDT not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. ( LSAME( SYDT_COMM, 'D' ) 
     $        .OR. LSAME( SYDT_COMM, 'S' ) ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Communication pattern in triangular solver must'//
     $              ' be specified for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_UPLOSIGN_MIN.LT.0 .OR. SYDT_UPLOSIGN_MIN.GT.7 .OR.
     $        SYDT_UPLOSIGN_MAX.LT.0 .OR. SYDT_UPLOSIGN_MAX.GT.7 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', SYDT_ASCHUR ) .AND. 
     $        .NOT. LSAME( 'S', SYDT_ASCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of A must be specified for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', SYDT_BSCHUR ) .AND. 
     $        .NOT. LSAME( 'S', SYDT_BSCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of B must be specified for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_A_NQTRBL.LT.0 .OR. SYDT_A_NQTRBL.GT.100 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= A_NQTRBL <= 100 must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_B_NQTRBL.LT.0 .OR. SYDT_B_NQTRBL.GT.100 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= B_NQTRBL <= 100 must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_M_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN >= 1 must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_M_MIN .GT. SYDT_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN <= M_MAX must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_N_MIN .GT. SYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_M_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_MB_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MIN >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_MB_MAX.GT.SYDT_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MAX <= M_MAX  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_NB_MAX.GT.SYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_MB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_MB2_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MIN >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_MB2_MAX.GT.SYDT_MB_MIN .AND. SYDT_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MAX <= MB_MIN  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_IA_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MIN >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_IA_MAX.LT.SYDT_IA_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX >= IA_MIN  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_IA_MAX.GT.SYDT_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX <= M_MAX  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
          IF( SYDT_IB_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_MIN >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_IB_MAX.LT.SYDT_IB_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_MAX >= IB_MIN  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( SYDT_IB_MAX.GT.SYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_MAX <= N_MAX  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_IA_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_IB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IB_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_NPROW_MAX.LT.SYDT_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_NPCOL_MAX.LT.SYDT_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( SYDT_NPROW_MAX*SYDT_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYDT_NPROW_MAX*SYDT_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for SYDT !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( SYDT_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( SYDT_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for SYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Checking LYDT
C
      IF( LYDT ) THEN
         IF( LYDT_CONDEST .AND. .NOT. ( LYDT_REDUCE_ONLY .OR. 
     $        LYDT_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of LYDT not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_UPLOSIGN_MIN.LT.0 .OR. LYDT_UPLOSIGN_MIN.GT.1 .OR.
     $        LYDT_UPLOSIGN_MAX.LT.0 .OR. LYDT_UPLOSIGN_MAX.GT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', LYDT_ASCHUR ) .AND. 
     $        .NOT. LSAME( 'S', LYDT_ASCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of A must be specified for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', LYDT_C_SYMM ) .AND. 
     $        .NOT. LSAME( 'S', LYDT_C_SYMM ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Symmetric mode of C must be specified for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_N_MIN .GT. LYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_NB_MAX.GT.LYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_NB2_MAX.GT.LYDT_NB_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_MAX <= NB_MIN  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_NB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_STEP >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_IA_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MIN >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_IA_MAX.LT.LYDT_IA_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX >= IA_MIN  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( LYDT_IA_MAX.GT.LYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_MAX <= N_MAX  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_IA_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IA_STEP >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_NPROW_MAX.LT.LYDT_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_NPCOL_MAX.LT.LYDT_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( LYDT_NPROW_MAX*LYDT_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'LYDT_NPROW_MAX*LYDT_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for LYDT !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( LYDT_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( LYDT_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for LYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Checking GCSY
C
      IF( GCSY ) THEN
         IF( GCSY_CONDEST .AND. .NOT. ( GCSY_REDUCE_ONLY .OR. 
     $        GCSY_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of GCSY not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. ( LSAME( GCSY_COMM, 'D' ) 
     $        .OR. LSAME( GCSY_COMM, 'S' ) ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Communication pattern in triangular solver must'//
     $              ' be specified for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_UPLOSIGN_MIN.LT.0 .OR. GCSY_UPLOSIGN_MIN.GT.7 .OR.
     $        GCSY_UPLOSIGN_MAX.LT.0 .OR. GCSY_UPLOSIGN_MAX.GT.7 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GCSY_ADSCHUR ) .AND. 
     $        .NOT. LSAME( 'S', GCSY_ADSCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of (A,D) must be'//
     $              ' specified for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GCSY_BESCHUR ) .AND. 
     $        .NOT. LSAME( 'S', GCSY_BESCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of (B,E) must be'//
     $              ' specified for GCSY !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GCSY_TRANSZ ) .AND. 
     $        .NOT. LSAME( 'T', GCSY_TRANSZ ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Transpose mode of Z must be'//
     $              ' specified for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_AD_NQTRBL.LT.0 .OR. GCSY_AD_NQTRBL.GT.100 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= AD_NQTRBL <= 100 must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_BE_NQTRBL.LT.0 .OR. GCSY_BE_NQTRBL.GT.100 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= BE_NQTRBL <= 100 must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_M_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN >= 1 must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_M_MIN .GT. GCSY_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN <= M_MAX must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_N_MIN .GT. GCSY_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_M_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_MB_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MIN >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_MB_MAX.GT.GCSY_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MAX <= M_MAX  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_NB_MAX.GT.GCSY_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_MB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_MB2_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MIN >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_MB2_MAX.GT.GCSY_MB_MIN .AND. GCSY_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MAX <= MB_MIN  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_NB2_MAX.GT.GCSY_NB_MIN .AND. GCSY_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_MAX <= NB_MIN  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_MB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_NB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_IAD_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAD_MIN >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_IAD_MAX.LT.GCSY_IAD_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAD_MAX >= IAD_MIN  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_IAD_MAX.GT.GCSY_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAD_MAX <= M_MAX  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
          IF( GCSY_IBE_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBE_MIN >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_IBE_MAX.LT.GCSY_IBE_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBE_MAX >= IBE_MIN  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GCSY_IBE_MAX.GT.GCSY_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBE_MAX <= N_MAX  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_IAD_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAD_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_IBE_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBE_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_NPROW_MAX.LT.GCSY_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_NPCOL_MAX.LT.GCSY_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( GCSY_NPROW_MAX*GCSY_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY_NPROW_MAX*GCSY_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for GCSY !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( GCSY_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GCSY_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for GCSY!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Checking GSYL
C
      IF( GSYL ) THEN
         IF( GSYL_CONDEST .AND. .NOT. ( GSYL_REDUCE_ONLY .OR. 
     $        GSYL_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of GSYL not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. ( LSAME( GSYL_COMM, 'D' ) 
     $        .OR. LSAME( GSYL_COMM, 'S' ) ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Communication pattern in triangular solver must'//
     $              ' be specified for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_UPLOSIGN_MIN.LT.0 .OR. GSYL_UPLOSIGN_MIN.GT.7 .OR.
     $        GSYL_UPLOSIGN_MAX.LT.0 .OR. GSYL_UPLOSIGN_MAX.GT.7 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GSYL_ACSCHUR ) .AND. 
     $        .NOT. LSAME( 'S', GSYL_ACSCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of (A,C) must'//
     $              ' be specified for GSYL !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GSYL_BDSCHUR ) .AND. 
     $        .NOT. LSAME( 'S', GSYL_BDSCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of (B,D) must'//
     $              ' be specified for GSYL !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_AC_NQTRBL.LT.0 .OR. GSYL_AC_NQTRBL.GT.100 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= AC_NQTRBL <= 100 must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_BD_NQTRBL.LT.0 .OR. GSYL_BD_NQTRBL.GT.100 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= BD_NQTRBL <= 100 must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_M_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN >= 1 must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_M_MIN .GT. GSYL_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_MIN <= M_MAX must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_N_MIN .GT. GSYL_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_M_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'M_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_MB_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MIN >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_MB_MAX.GT.GSYL_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_MAX <= M_MAX  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_NB_MAX.GT.GSYL_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_MB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_MB2_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MIN >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_MB2_MAX.GT.GSYL_MB_MIN .AND. GSYL_PIPELINE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_MAX <= MB_MIN  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_MB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'MB2_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_IAC_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAC_MIN >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_IAC_MAX.LT.GSYL_IAC_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAC_MAX >= IAC_MIN  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_IAC_MAX.GT.GSYL_M_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAC_MAX <= M_MAX  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
          IF( GSYL_IBD_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBD_MIN >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_IBD_MAX.LT.GSYL_IBD_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBD_MAX >= IBD_MIN  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GSYL_IBD_MAX.GT.GSYL_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBD_MAX <= N_MAX  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_IAC_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAC_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_IBD_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IBD_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_NPROW_MAX.LT.GSYL_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_NPCOL_MAX.LT.GSYL_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( GSYL_NPROW_MAX*GSYL_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GSYL_NPROW_MAX*GSYL_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for GSYL !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( GSYL_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GSYL_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for GSYL!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Checking GLYCT
C     
      IF( GLYCT ) THEN
         IF( GLYCT_CONDEST .AND. .NOT. ( GLYCT_REDUCE_ONLY .OR. 
     $        GLYCT_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of GLYCT not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_UPLOSIGN_MIN.LT.0 .OR. GLYCT_UPLOSIGN_MIN.GT.1 .OR.
     $        GLYCT_UPLOSIGN_MAX.LT.0 .OR. GLYCT_UPLOSIGN_MAX.GT.1 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GLYCT_AESCHUR ) .AND. 
     $        .NOT. LSAME( 'S', GLYCT_AESCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of (A,E) must be specified for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GLYCT_C_SYMM ) .AND. 
     $        .NOT. LSAME( 'S', GLYCT_C_SYMM ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Symmetric mode of C must be specified for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_AE_NQTRBL.LT.0 .OR. GLYCT_AE_NQTRBL.GT.100 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= AE_NQTRBL <= 100 must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_N_MIN .GT. GLYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_NB_MAX.GT.GLYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_NB2_MAX.GT.GLYCT_NB_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_MAX <= NB_MIN  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_NB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_STEP >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_IAE_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_MIN >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_IAE_MAX.LT.GLYCT_IAE_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_MAX >= IAE_MIN  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYCT_IAE_MAX.GT.GLYCT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_MAX <= N_MAX  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_IAE_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_STEP >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_NPROW_MAX.LT.GLYCT_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_NPCOL_MAX.LT.GLYCT_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( GLYCT_NPROW_MAX*GLYCT_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYCT_NPROW_MAX*GLYCT_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for GLYCT !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( GLYCT_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYCT_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for GLYCT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C      
C     Checking GLYDT
C
      IF( GLYDT ) THEN
         IF( GLYDT_CONDEST .AND. .NOT. ( GLYDT_REDUCE_ONLY .OR. 
     $        GLYDT_SOLVE ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Condition estimation of GLYDT not possible'//
     $              ' without reduction to triangular form!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_UPLOSIGN_MIN.LT.0 .OR. GLYDT_UPLOSIGN_MIN.GT.1 .OR.
     $        GLYDT_UPLOSIGN_MAX.LT.0 .OR. GLYDT_UPLOSIGN_MAX.GT.1 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'UPLOSIGN specifications are not valid for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GLYDT_AESCHUR ) .AND. 
     $        .NOT. LSAME( 'S', GLYDT_AESCHUR ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Schur mode of (A,E) must be specified for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( .NOT. LSAME( 'N', GLYDT_C_SYMM ) .AND. 
     $        .NOT. LSAME( 'S', GLYDT_C_SYMM ) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'Symmetric mode of C must be specified for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_AE_NQTRBL.LT.0 .OR. GLYDT_AE_NQTRBL.GT.100 ) 
     $        THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              '0 <= AE_NQTRBL <= 100 must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_N_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN >= 1 must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_N_MIN .GT. GLYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_MIN <= N_MAX must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_N_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'N_STEP >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_NB_MAX.GT.GLYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_MAX <= N_MAX  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_NB_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB_STEP >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_NB2_MAX.GT.GLYDT_NB_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_MAX <= NB_MIN  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_NB2_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NB2_STEP >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_IAE_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_MIN >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_IAE_MAX.LT.GLYDT_IAE_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_MAX >= IAE_MIN  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998 
         END IF
         IF( GLYDT_IAE_MAX.GT.GLYDT_N_MAX ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_MAX <= N_MAX  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_IAE_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'IAE_STEP >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_NPROW_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MIN >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_NPCOL_MIN.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MIN >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_NPROW_MAX.LT.GLYDT_NPROW_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_MAX >= NPROW_MIN  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_NPCOL_MAX.LT.GLYDT_NPCOL_MIN ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_MAX >= NPCOL_MIN  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#ifdef LOOPGRID
         IF( GLYDT_NPROW_MAX*GLYDT_NPCOL_MAX .GT. 
     $        MIN(NPROCS,MAX_CPUS) ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYDT_NPROW_MAX*GLYDT_NPCOL_MAX <= '//
     $              'MIN(NPROCS,MAX_CPUS)'//
     $              ' must hold for GLYDT !'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
#endif
         IF( GLYDT_NPROW_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPROW_STEP >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_NPCOL_STEP.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'NPCOL_STEP >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
         IF( GLYDT_REPEAT.LT.1 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'REPEAT >= 1  must hold for GLYDT!'
               WRITE( NOUT , FMT = *)
            END IF
            GO TO 9998  
         END IF
      END IF
C
C     Parameter checking done, allocate memory
C
#ifdef USE_DYNAMIC
      ALLOCATE( MEM(MEMSIZ), IMEM(IMEMSIZ) )
#endif
C
C     Read out the default number of OpenMP threads
C
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
      THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
#else
      THREADS = 1
#endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON SYCT ===
C
      IF( SYCT ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for SYCT ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6669 )
     $           'ULS','M','N','MB','NB','MB2','NB2','IA','IB','P_r',
     $           'P_c','Thr','Tme','Tp1','Tp2','Tp3','Gfl','Scl','E_a',
     $           'E_r','R_a','R_r','est','etm','Cp1','Cp2','Cp3','itr',
     $           'INF'
         END IF
C
C     Loop though the space of all combinations of parameters,
C     M, N, MB, NB, NPROW, NPCOL
C
         DO UPLOSIGN = SYCT_UPLOSIGN_MIN, SYCT_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 .OR. UPLOSIGN.EQ.1 ) THEN
            SYCT_TRANSA = 'N'
            SYCT_TRANSB = 'N'
         ELSEIF( UPLOSIGN.EQ.2 .OR. UPLOSIGN.EQ.3 ) THEN
            SYCT_TRANSA = 'N'
            SYCT_TRANSB = 'T'
         ELSEIF( UPLOSIGN.EQ.4 .OR. UPLOSIGN.EQ.5 ) THEN
            SYCT_TRANSA = 'T'
            SYCT_TRANSB = 'N'
         ELSEIF( UPLOSIGN.EQ.6 .OR. UPLOSIGN.EQ.7 ) THEN
            SYCT_TRANSA = 'T'
            SYCT_TRANSB = 'T'
         END IF
         SYCT_ISGN = 1.0D+00-(MOD(UPLOSIGN,2))*2.0D+00
         DO M = SYCT_M_MIN, SYCT_M_MAX, SYCT_M_STEP
         DO N = SYCT_N_MIN, SYCT_N_MAX, SYCT_N_STEP
         DO MB = SYCT_MB_MIN, SYCT_MB_MAX, SYCT_MB_STEP
         DO NB = SYCT_NB_MIN, SYCT_NB_MAX, SYCT_NB_STEP
         DO MB2 = SYCT_MB2_MIN, SYCT_MB2_MAX, SYCT_MB2_STEP
         DO NB2 = SYCT_NB2_MIN, SYCT_NB2_MAX, SYCT_NB2_STEP
#ifdef LOOPGRID
         DO TEMP_NPROW = SYCT_NPROW_MIN, SYCT_NPROW_MAX, SYCT_NPROW_STEP
         DO TEMP_NPCOL = SYCT_NPCOL_MIN, SYCT_NPCOL_MAX, SYCT_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0
C
C     For the current set of parameters, set up temporary grid
C     
#ifdef LOOPGRID      
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL 
#endif   
C
C     Exclude processes not belonging to the temporary grid OR skip
C     the current set of parameters if not valid
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $           MYCOL.LT.0 .OR. MYCOL.GE.NPCOL .OR.
     $           (SYCT_SQUARE_DIM .AND. M.NE.N) .OR.
     $           (SYCT_SQUARE_BLK .AND. MB.NE.NB)  .OR.
     $           (SYCT_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) .OR.
     $           (.NOT.SYCT_PIPELINE .AND. (MB.NE.MB2) ) ) GO TO 20
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            AROWS = NUMROC( SYCT_M_MAX, MB, MYROW, 0, NPROW )
            ACOLS = NUMROC( SYCT_M_MAX, MB, MYCOL, 0, NPCOL )
            BROWS = NUMROC( SYCT_N_MAX, NB, MYROW, 0, NPROW )
            BCOLS = NUMROC( SYCT_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, SYCT_M_MAX, SYCT_M_MAX, MB, MB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYCT: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 9998
            END IF 
            CALL DESCINIT( DESCB, SYCT_N_MAX, SYCT_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, BROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYCT: DESCINIT(B), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 9998
            END IF 
            IF( SYCT_SOLVE ) 
     $           CALL DESCINIT( DESCC, SYCT_M_MAX, SYCT_N_MAX, MB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYCT: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 9998
            END IF
C
C     Set internal blocking factors
C
            MBNB2(1) = MB2
            MBNB2(2) = NB2
C
C     Loop through the space of parameter IA, JA, IB, JB, IC, JC and
C     ignore invalid cases
C
            DO IA = SYCT_IA_MIN, SYCT_IA_MAX, SYCT_IA_STEP
            JA = IA
            IC = IA
            IF( .NOT. ( IA + M - 1 ).LE.SYCT_M_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IA = ',IA
               GO TO 17
            END IF
            DO IB = SYCT_IB_MIN, SYCT_IB_MAX, SYCT_IB_STEP
            JB = IB
            JC = JB
            IF( .NOT. ( IB + N - 1 ).LE.SYCT_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) 'IB,N,SYCT_N_MAX=',IB,N,SYCT_N_MAX
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IB = ',IB
               GO TO 15
            END IF
            DO REPEAT = 1, SYCT_REPEAT, 1
C
C     Assign pointer for ScaLAPACK arrays
C
               IPA = 1
               IF( .NOT. LSAME( SYCT_ASCHUR, 'S') ) THEN
                  IPACPY = IPA + DESCA( LLD_ ) * ACOLS
               ELSE
                  IPACPY = IPA
               END IF
               IPB = IPACPY + DESCA( LLD_ ) * ACOLS
               IF( .NOT. LSAME( SYCT_BSCHUR, 'S') ) THEN
                  IPBCPY = IPB + DESCB( LLD_ ) * BCOLS
               ELSE
                  IPBCPY = IPB
               END IF
               IF( SYCT_SOLVE ) THEN
                  IPC = IPBCPY + DESCB( LLD_ ) * BCOLS
                  IPCCPY = IPC + DESCC( LLD_ ) * BCOLS
                  IPX = IPCCPY + DESCC( LLD_ ) * BCOLS
                  DDIAG = IPX + DESCC( LLD_ ) * BCOLS 
                  SDIAG = DDIAG + MAX( SYCT_M_MAX, SYCT_N_MAX )
                  IPW = SDIAG + MAX( SYCT_M_MAX-1, SYCT_N_MAX-1 )
               ELSE
                  IPC = 1
                  IPCCPY = 1
                  IPX = 1
                  DDIAG = IPBCPY + DESCB( LLD_ ) * BCOLS 
                  SDIAG = DDIAG + MAX( SYCT_M_MAX, SYCT_N_MAX )
                  IPW = SDIAG + MAX( SYCT_M_MAX-1, SYCT_N_MAX-1 )
               END IF
C
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of A
C
               DO I = 1, SYCT_M_MAX
                  MEM(DDIAG+(I-1)) = DBLE(I)
               END DO
C
C     Init A
C
               BLKS2B2 = INT(DBLE(SYCT_M_MAX/2) * 
     $                   DBLE(SYCT_A_NQTRBL)/100.0D+00)
               CALL P1SQMATGD( SYCT_ADIAG, 'No subdiagonal', SYCT_AFORM,
     $                         SYCT_M_MAX, MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DDIAG), MEM(SDIAG), BLKS2B2, 17,
     $                         MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'SYCT: P1SQMATGD(A), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 20
               END IF 
C
C     Copy A
C
               IF( .NOT. LSAME( SYCT_ASCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', SYCT_M_MAX, SYCT_M_MAX, MEM(IPA),
     $                          1, 1, DESCA, MEM(IPACPY), 1, 1, DESCA )
               END IF
C
C     Set eigenvalues of B
C
               DO I = 1, SYCT_N_MAX
                  MEM(DDIAG+(I-1)) = DBLE(SYCT_ISGN)*(DBLE(I))
               END DO      
C     
C     Initialize B
C
               BLKS2B2 = INT(DBLE(SYCT_N_MAX/2) * 
     $                   DBLE(SYCT_B_NQTRBL)/100.0D+00)
               CALL P1SQMATGD( SYCT_BDIAG, 'No subdiagonal', SYCT_BFORM, 
     $                         SYCT_N_MAX, MEM(IPB), DESCB, ONE, ONE, 
     $                         MEM(DDIAG), MEM(SDIAG), BLKS2B2, 19,
     $                         MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'SYCT: P1SQMATGD(B), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 20
               END IF 
C
C     Copy B
C     
               IF( .NOT. LSAME( SYCT_BSCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', SYCT_N_MAX, SYCT_N_MAX, MEM(IPB),
     $                          1, 1, DESCB, MEM(IPBCPY), 1, 1, DESCB )
               END IF
C     
C     Init a matrix X as the given "solution"
C     
               IF( SYCT_SOLVE )
     $              CALL PDMATGEN2( TEMP_ICTXT, SYCT_CFORM, SYCT_CDIAG, 
     $                              SYCT_M_MAX, SYCT_N_MAX, MB, NB, 
     $                              MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, 
     $                              AROWS, 0, BCOLS, MYROW, MYCOL, 
     $                              NPROW, NPCOL )
C
C     Compute right hand side C as op(A) * X +/- X * op(B)
C
               IF( SYCT_SOLVE ) THEN
                  CALL PDGEMM( SYCT_TRANSA, 'N' , M, N, M, ONE, 
     $                         MEM(IPA), IA, JA, DESCA, MEM(IPX), IC, 
     $                         JC, DESCC, ZERO, MEM(IPC), IC, JC, 
     $                         DESCC )
                  CALL PDGEMM( 'N', SYCT_TRANSB, M, N, N, 
     $                         DBLE(SYCT_ISGN), MEM(IPX), IC, JC, DESCC, 
     $                         MEM(IPB), IB, JB, DESCB, ONE, MEM(IPC), 
     $                         IC, JC, DESCC )
C     
C     Copy C for future reference
C
                  CALL PDLACPY( 'All', M, N, MEM(IPC), IC, JC, DESCC, 
     $                          MEM(IPCCPY), IC, JC, DESCC )
               END IF
C 
C     If the matrices are small enough, print them to the screen
C     
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(A):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, M, MEM( IPA ), IA, JA, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF
C     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(B):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPB ), IB, JB, DESCB, 
     $                           0, 0, 'B', NOUT, MEM( IPW ) )
               END IF     
C     
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $              .AND. SYCT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(C):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. SYCT_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
C               
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $              .AND. SYCT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X0):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. SYCT_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPX ), IC, JC, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
C
C     Call SYCT solver routine
C
               INFO = 0; SCALE = ZERO
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( SYCT_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGESYCTD( 'Reduce', SYCT_ASCHUR, SYCT_BSCHUR, 
     $                           SYCT_TRANSA, SYCT_TRANSB, SYCT_ISGN, 
     $                           SYCT_COMM, M, N, MEM(IPA), IA, JA,
     $                           DESCA, MEM(IPB), IB, JB, DESCB, 
     $                           MEM(IPC), IC, JC,  DESCC, MBNB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( SYCT_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGESYCTD( 'Solve', SYCT_ASCHUR, SYCT_BSCHUR, 
     $                           SYCT_TRANSA, SYCT_TRANSB, SYCT_ISGN, 
     $                           SYCT_COMM, M, N, MEM(IPA), IA, JA,
     $                           DESCA, MEM(IPB), IB, JB, DESCB, 
     $                           MEM(IPC), IC, JC,  DESCC, MBNB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'SYCT: PGESYCTD, INFO =',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 20
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'SYCT: PGESYCTD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. SYCT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. 
     $             SYCT_SOLVE ) THEN
                  CALL PDLAPRNT( M, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Gflops/sec.
C
               GFLOPS = (DBLE(M)*DBLE(N)**2+DBLE(M)**2*DBLE(N)) /
     $              (10.0D+00 ** 9 * TPROF3)
C     
C     In reduction mode only, move on to condition estimation
C
               IF( SYCT_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 10 
               END IF
C
C     Compute results of computation
C     
C                          ||op(A) * X  +/- X * op(B) - SCALE * C||
C     Compute relres =  ---------------------------------------------
C                       ((||X||*(||A||+||B||) + SCALE * ||C||) * eps )
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
               CNORM = PDLANGE( 'I', M, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', M, N, MEM( IPC ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDGEMM( SYCT_TRANSA, 'N', M, N, M, -ONE, 
     $                      MEM( IPACPY ), IA, JA, DESCA, MEM( IPC ), 
     $                      IC, JC, DESCC, SCALE, MEM( IPCCPY ), IC, JC,
     $                      DESCC )
               CALL PDGEMM( 'N', SYCT_TRANSB, M, N, N, 
     $                      DBLE(-1*SYCT_ISGN), MEM( IPC ), IC, JC, 
     $                      DESCC, MEM( IPBCPY ), IB, JB, DESCB, ONE, 
     $                      MEM( IPCCPY ), IC, JC, DESCC )
               IF( LSAME( SYCT_TRANSA, 'T' ) ) THEN
                  ANORM = PDLANGE( 'I', M, M, MEM( IPACPY ), IA, JA, 
     $                            DESCA, MEM( IPW ) )
               ELSE
                  ANORM = PDLANGE( '1', M, M, MEM( IPACPY ), IA, JA, 
     $                             DESCA, MEM( IPW ) )
               END IF
               IF( LSAME( SYCT_TRANSB, 'T' ) ) THEN
                  BNORM = PDLANGE( 'I', N, N, MEM( IPBCPY ), IB, JB, 
     $                             DESCB, MEM( IPW ) )
               ELSE
                  BNORM = PDLANGE( '1', N, N, MEM( IPBCPY ), IB, JB, 
     $                             DESCB, MEM( IPW ) )
               END IF
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = PDLANGE( 'F', M, N, MEM( IPCCPY ), IC, JC, 
     $                             DESCC, MEM(IPW) ) 
               RESID =  PDLANGE( 'I', M, N, MEM( IPCCPY ), IC, JC, 
     $                           DESCC, MEM(IPW)) / (((ANORM + BNORM) * 
     $                           XNORM + SCALE * CNORM)*EPS)
C     
C     Compute real error norms
C
               XNORM = PDLANGE( 'F', M, N, MEM( IPX ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDMATADD( M, N, SCALE, MEM(IPX), IC, JC, DESCC, 
     $                        -ONE, MEM(IPC), IC, JC, DESCC )
               ABSERROR = PDLANGE( 'F', M, N, MEM( IPC ), IC, JC, DESCC, 
     $                             MEM(IPW) )
               RELERROR = ABSERROR / XNORM
C
C     From reduction mode we arrive here
C
 10            CONTINUE
C
C     Condition estimation
C
               IF( SYCT_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                  CALL PSYCTCON( SYCT_TRANSA, SYCT_TRANSB, SYCT_ISGN,
     $                           SYCT_COMM, M, N, MEM( IPA ), IA, JA, 
     $                           DESCA, MEM( IPB ), IB, JB, DESCB, 
     $                           MBNB2, MEM( IPW ), WORKSIZ, IMEM, 
     $                           IMEMSIZ, EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'SYCT: PSYCTCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 20
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'SYCT: PSYCTCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7779 )
     $                 UPLOSIGN,M,N,MB,NB,MB2,NB2,IA,IB,NPROW,NPCOL,
     $                 THREADS,MYTIME,TPROF1,TPROF2,TPROF3,GFLOPS,SCALE,
     $                 ABSERROR,RELERROR,ABSRESID,RESID,EST,ESTTIME,
     $                 CPROF1,CPROF2,CPROF3,NOITER,INFO
               END IF
C
            END DO
 15         CONTINUE
            END DO         
 17         CONTINUE
            END DO
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 20         CONTINUE
C
C     Synchronize all processes.
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO
         END DO
#endif
         END DO               
         END DO
         END DO
         END DO
         END DO
         END DO
         END DO
C         
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON LYCT === 
C
      IF( LYCT ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for LYCT ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6670 )
     $           'ULS','N','NB','NB2','IA','P_r','P_c','Thr','Tme',
     $           'Tp1','Tp2','Tp3','Gfl','Scl','E_a','E_r','R_a','R_r',
     $           'est','etm','Cp1','Cp2','Cp3','itr','INF'
         END IF
C
C     Loop though the space of all combinations of parameters,
C     N, NB, NPROW, NPCOL
C
         DO UPLOSIGN = LYCT_UPLOSIGN_MIN, LYCT_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 ) THEN
            LYCT_TRANSA = 'N'
         ELSEIF( UPLOSIGN.EQ.1 ) THEN
            LYCT_TRANSA = 'T'
         END IF
         DO N = LYCT_N_MIN, LYCT_N_MAX, LYCT_N_STEP
         DO NB = LYCT_NB_MIN, LYCT_NB_MAX, LYCT_NB_STEP
         DO NB2 = LYCT_NB2_MIN, LYCT_NB2_MAX, LYCT_NB2_STEP
#ifdef LOOPGRID
         DO TEMP_NPROW = LYCT_NPROW_MIN, LYCT_NPROW_MAX, LYCT_NPROW_STEP
         DO TEMP_NPCOL = LYCT_NPCOL_MIN, LYCT_NPCOL_MAX, LYCT_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C         
#ifdef LOOPGRID           
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid
C     
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $           MYCOL.LT.0 .OR. MYCOL.GE.NPCOL  .OR.
     $           (LYCT_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) .OR.
     $           (.NOT.LYCT_PIPELINE .AND. NB.NE.NB2) ) GO TO 40
C     
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            AROWS = NUMROC( LYCT_N_MAX, NB, MYROW, 0, NPROW )
            ACOLS = NUMROC( LYCT_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, LYCT_N_MAX, LYCT_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'LYCT: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 40
            END IF 
            IF( LYCT_SOLVE )
     $           CALL DESCINIT( DESCC, LYCT_N_MAX, LYCT_N_MAX, NB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'LYCT: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 40
            END IF 
C
C     Loop through the space of parameter IA, JA and
C     ignore invalid cases
C
            DO IA = LYCT_IA_MIN, LYCT_IA_MAX, LYCT_IA_STEP
            JA = IA
            IC = IA
            JC = JA
            IF( .NOT. ( IA + N - 1 ).LE.LYCT_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IA = ',IA
               GO TO 35
            END IF
            DO REPEAT = 1, LYCT_REPEAT, 1
C
C     Assign pointer for ScaLAPACK arrays
C
               IPA = 1
               IF( .NOT. LSAME( LYCT_ASCHUR, 'S') ) THEN
                  IPACPY = IPA + DESCA( LLD_ ) * ACOLS
               ELSE
                  IPACPY = IPA
               END IF
               IF( LYCT_SOLVE ) THEN
                  IPC = IPACPY + DESCA( LLD_ ) * ACOLS
                  IPCCPY = IPC + DESCC( LLD_ ) * ACOLS
                  IPX = IPCCPY + DESCC( LLD_ ) * ACOLS
                  DDIAG = IPX + DESCC( LLD_ ) * ACOLS 
                  SDIAG = DDIAG + LYCT_N_MAX 
                  IPW = SDIAG + LYCT_N_MAX-1
               ELSE
                  IPC = 1
                  IPCCPY = 1
                  IPX = 1
                  DDIAG = IPACPY + DESCA( LLD_ ) * ACOLS 
                  SDIAG = DDIAG + LYCT_N_MAX 
                  IPW = SDIAG + LYCT_N_MAX-1
               END IF
C     
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of A
C
               DO I = 1, LYCT_N_MAX
                  MEM(DDIAG+(I-1)) = DBLE(I)+1000.0D+0
               END DO
C
C     Init A
C
               BLKS2B2 = INT(DBLE(LYCT_N_MAX/2) * 
     $                   DBLE(LYCT_A_NQTRBL)/100.0D+00)
               CALL P1SQMATGD( LYCT_ADIAG, 'No subdiagonal', LYCT_AFORM,
     $                         LYCT_N_MAX, MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DDIAG), MEM(SDIAG), BLKS2B2, 17,
     $                         MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'LYCT: P1SQMATGD(A), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 40
               END IF 
C
C     Copy A
C
               IF( .NOT. LSAME( LYCT_ASCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', LYCT_N_MAX, LYCT_N_MAX, MEM(IPA),
     $                          1, 1, DESCA, MEM(IPACPY), 1, 1, DESCA )
               END IF
C     
C     Init a matrix X as the given solution
C     
               IF( LYCT_SOLVE ) THEN
                  CALL PDMATGEN2( TEMP_ICTXT, 'Symmetric', LYCT_CDIAG, 
     $                            LYCT_N_MAX, LYCT_N_MAX, NB, NB, 
     $                            MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, 
     $                            AROWS, 0, ACOLS, MYROW, MYCOL, NPROW, 
     $                            NPCOL )
C
C     Compute right hand side C as op(A) * X + X * op(A')
C
                  IF( LSAME( LYCT_TRANSA, 'N' ) ) THEN
                     CALL PDGEMM( 'N', 'N' , N, N, N, ONE, MEM(IPA), 
     $                    IA, JA, DESCA, MEM(IPX), IC, JC,
     $                    DESCC, ZERO, MEM(IPC), IC, JC, DESCC )
                     CALL PDGEMM( 'N', 'Transpose', N, N, N, 
     $                    ONE, MEM(IPX), IC, JC, DESCC,
     $                    MEM(IPA), IA, JA, DESCA, ONE, 
     $                    MEM(IPC), IC, JC, DESCC )
                  ELSEIF( LSAME( LYCT_TRANSA, 'T' ) ) THEN
                     CALL PDGEMM( 'Transpose', 'N' , N, N, N, ONE, 
     $                    MEM(IPA), IA, JA, DESCA, MEM(IPX), IC, 
     $                    JC, DESCC, ZERO, MEM(IPC), IC, JC, 
     $                    DESCC )
                     CALL PDGEMM( 'N', 'N', N, N, N, ONE, 
     $                    MEM(IPX), IC, JC, DESCC, MEM(IPA), IA, 
     $                    JA, DESCA, ONE, MEM(IPC), IC, JC, 
     $                    DESCC )
                  END IF
C     
C     Copy C for future reference
C
                  CALL PDLACPY( 'All', N, N, MEM(IPC), IC, JC, DESCC, 
     $                          MEM(IPCCPY), IC, JC, DESCC )
               END IF
C
C     Check if sub(C) is symmetric for this particular instance
C
               IF( LSAME( LYCT_C_SYMM, 'S' ) .AND. IC.EQ.JC ) THEN
                  LYCT_C_SYMM2 = 'S'
               ELSE
                  LYCT_C_SYMM2 = 'N'
               END IF
C 
C     If the matrices are small enough, print them to the screen
C     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(A):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPA ), IA, JA, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF
C     
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. LYCT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(C):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. LYCT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
               IF( N.LE.PRNTSIZ .AND. LYCT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPX ), IC, JC, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
C
C     Call LYCT solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( LYCT_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGELYCTD( 'Reduce', LYCT_C_SYMM2, LYCT_TRANSA,
     $                           LYCT_ASCHUR, N, MEM(IPA), IA, JA, 
     $                           DESCA, MEM(IPC), IC, JC, DESCC, NB2, 
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( LYCT_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGELYCTD( 'Solve', LYCT_C_SYMM2, LYCT_TRANSA,
     $                           LYCT_ASCHUR, N, MEM(IPA), IA, JA, 
     $                           DESCA, MEM(IPC), IC, JC, DESCC, NB2, 
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'LYCT: PGELYCTD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 40
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'LYCT: PGELYCTD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. LYCT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. LYCT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Gflops/sec.
C
               IF( LSAME( LYCT_C_SYMM, 'S' ) ) THEN
                  GFLOPS = (1.0D+00*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               ELSE
                  GFLOPS = (2.0D+00*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               END IF
C
C     In reduction mode only, move on to condition estimation
C
               IF( LYCT_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 30 
               END IF
C
C     Compute results of computation
C     
C                          ||op(A) * X + X * op(A') - SCALE * C||
C     Compute relres =  ---------------------------------------------
C                       ((||X||*(2*||A||) + SCALE * ||C||) * eps )
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
               CNORM = PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDGEMM( LYCT_TRANSA, 'N', N, N, N, -ONE, 
     $                      MEM( IPACPY ), IA, JA, DESCA, MEM( IPC ), 
     $                      IC, JC, DESCC, SCALE, MEM( IPCCPY ), IC, JC,
     $                      DESCC )
               IF( LSAME( LYCT_TRANSA, 'N') ) THEN
                  CALL PDGEMM( 'N', 'Transpose', N, N, N, 
     $                         -ONE, MEM( IPC ), IC, JC, 
     $                         DESCC, MEM( IPACPY ), IA, JA, DESCA, ONE, 
     $                         MEM( IPCCPY ), IC, JC, DESCC )
               ELSE
                  CALL PDGEMM( 'N', 'No Transpose', N, N, N, 
     $                         -ONE, MEM( IPC ), IC, JC, 
     $                         DESCC, MEM( IPACPY ), IA, JA, DESCA, ONE, 
     $                         MEM( IPCCPY ), IC, JC, DESCC )
               END IF
               IF( LSAME( LYCT_TRANSA, 'T' ) ) THEN
                  ANORM = PDLANGE( 'I', N, N, MEM( IPACPY ), IA, JA, 
     $                            DESCA, MEM( IPW ) )
               ELSE
                  ANORM = PDLANGE( '1', N, N, MEM( IPACPY ), IA, JA, 
     $                             DESCA, MEM( IPW ) )
               END IF
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = PDLANGE( 'F', N, N, MEM( IPCCPY ), IC, JC, 
     $                             DESCC, MEM(IPW) ) 
               RESID =  PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, 
     $                           DESCC, MEM(IPW)) / (((ANORM + ANORM) * 
     $                           XNORM + SCALE * CNORM)*EPS)
C     
C     Compute real error norms
C
               XNORM = PDLANGE( 'F', N, N, MEM( IPX ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDMATADD( N, N, SCALE, MEM(IPX), IC, JC, DESCC, 
     $                        -ONE, MEM(IPC), IC, JC, DESCC )
               ABSERROR = PDLANGE( 'F', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                             MEM(IPW) )
               RELERROR = ABSERROR / XNORM
C
C     From reduction mode we arrive here
C
 30            CONTINUE
C
C     Condition estimation
C
               IF( LYCT_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                  CALL PLYCTCON( LYCT_TRANSA, N, MEM( IPA ), IA, JA,
     $                           DESCA, NB2, MEM( IPW ), WORKSIZ, IMEM, 
     $                           IMEMSIZ, EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'LYCT: PLYCTCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 40
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'LYCT: PLYCTCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
              
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7780 )
     $                 UPLOSIGN,N,NB,NB2,IA,NPROW,NPCOL,THREADS,MYTIME,
     $                 TPROF1,TPROF2,TPROF3,GFLOPS,SCALE,ABSERROR,
     $                 RELERROR,ABSRESID,RESID,EST,ESTTIME,CPROF1,
     $                 CPROF2,CPROF3,NOITER,INFO
               END IF
C
            END DO 
 35         CONTINUE
            END DO
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 40         CONTINUE
C
C     Synchronize all processes. 
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO 
         END DO 
#endif 
         END DO
         END DO
         END DO
         END DO
C
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON SYDT ===
C
      IF( SYDT ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for SYDT ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6669 )
     $           'ULS','M','N','MB','NB','MB2','NB2','IA','IB','P_r',
     $           'P_c','Thr','Tme','Tp1','Tp2','Tp3','Gfl','Scl','E_a',
     $           'E_r','R_a','R_r','est','etm','Cp1','Cp2','Cp3','itr',
     $           'INF'
         END IF
C
C     Loop though the space of all combinations of parameters,
C     M, N, MB, NB, NPROW, NPCOL
C
         DO UPLOSIGN = SYDT_UPLOSIGN_MIN, SYDT_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 .OR. UPLOSIGN.EQ.1 ) THEN
            SYDT_TRANSA = 'N'
            SYDT_TRANSB = 'N'
         ELSEIF( UPLOSIGN.EQ.2 .OR. UPLOSIGN.EQ.3 ) THEN
            SYDT_TRANSA = 'N'
            SYDT_TRANSB = 'T'
         ELSEIF( UPLOSIGN.EQ.4 .OR. UPLOSIGN.EQ.5 ) THEN
            SYDT_TRANSA = 'T'
            SYDT_TRANSB = 'N'
         ELSEIF( UPLOSIGN.EQ.6 .OR. UPLOSIGN.EQ.7 ) THEN
            SYDT_TRANSA = 'T'
            SYDT_TRANSB = 'T'
         END IF
         SYDT_ISGN = INT(1.0D+00-(MOD(UPLOSIGN,2))*2.0D+00)
         DO M = SYDT_M_MIN, SYDT_M_MAX, SYDT_M_STEP
         DO N = SYDT_N_MIN, SYDT_N_MAX, SYDT_N_STEP
         DO MB = SYDT_MB_MIN, SYDT_MB_MAX, SYDT_MB_STEP
         DO NB = SYDT_NB_MIN, SYDT_NB_MAX, SYDT_NB_STEP
         DO MB2 = SYDT_MB2_MIN, SYDT_MB2_MAX, SYDT_MB2_STEP
#ifdef LOOPGRID
         DO TEMP_NPROW = SYDT_NPROW_MIN, SYDT_NPROW_MAX, SYDT_NPROW_STEP
         DO TEMP_NPCOL = SYDT_NPCOL_MIN, SYDT_NPCOL_MAX, SYDT_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C        
#ifdef LOOPGRID            
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid OR
C     skip current set of parameters if not valid
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $          MYCOL.LT.0 .OR. MYCOL.GE.NPCOL .OR.
     $          (SYDT_SQUARE_DIM .AND. M.NE.N) .OR.
     $          (SYDT_SQUARE_BLK .AND. MB.NE.NB)  .OR.
     $          (SYDT_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) .OR.
     $          (.NOT.SYDT_PIPELINE .AND. (MB.NE.MB2) ) ) GO TO 60
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            AROWS = NUMROC( SYDT_M_MAX, MB, MYROW, 0, NPROW )
            ACOLS = NUMROC( SYDT_M_MAX, MB, MYCOL, 0, NPCOL )
            BROWS = NUMROC( SYDT_N_MAX, NB, MYROW, 0, NPROW )
            BCOLS = NUMROC( SYDT_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, SYDT_M_MAX, SYDT_M_MAX, MB, MB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYDT: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 60
            END IF 
            CALL DESCINIT( DESCB, SYDT_N_MAX, SYDT_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, BROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYDT: DESCINIT(B), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 60
            END IF
            IF( SYDT_SOLVE ) 
     $           CALL DESCINIT( DESCC, SYDT_M_MAX, SYDT_N_MAX, MB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'SYDT: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 60
            END IF 
C
C     Loop through the space of parameter IA, JA, IB, JB, IC, JC and
C     ignore invalid cases
C
            DO IA = SYDT_IA_MIN, SYDT_IA_MAX, SYDT_IA_STEP
            JA = IA
            IC = IA
            IF( .NOT. ( IA + M - 1 ).LE.SYDT_M_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IA = ',IA
               GO TO 57
            END IF
            DO IB = SYDT_IB_MIN, SYDT_IB_MAX, SYDT_IB_STEP
            JB = IB
            JC = JB
            IF( .NOT. ( IB + N - 1 ).LE.SYDT_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IB = ',IB
               GO TO 55
            END IF
            DO REPEAT = 1, SYDT_REPEAT, 1
C
C     Assign pointer for ScaLAPACK arrays
C
               IPA = 1
               IF( .NOT. LSAME( SYDT_ASCHUR, 'S') ) THEN
                  IPACPY = IPA + DESCA( LLD_ ) * ACOLS
               ELSE
                  IPACPY = IPA
               END IF
               IPB = IPACPY + DESCA( LLD_ ) * ACOLS
               IF( .NOT. LSAME( SYDT_BSCHUR, 'S') ) THEN
                  IPBCPY = IPB + DESCB( LLD_ ) * BCOLS
               ELSE
                  IPBCPY = IPB
               END IF
               IF( SYDT_SOLVE ) THEN
                  IPC = IPBCPY + DESCB( LLD_ ) * BCOLS
                  IPCCPY = IPC + DESCC( LLD_ ) * BCOLS
                  IPX = IPCCPY + DESCC( LLD_ ) * BCOLS
                  DDIAG = IPX + DESCC( LLD_ ) * BCOLS 
                  SDIAG = DDIAG + MAX( SYDT_M_MAX, SYDT_N_MAX )
                  IPW = SDIAG + MAX( SYDT_M_MAX-1, SYDT_N_MAX-1 )
               ELSE
                  IPC = 1
                  IPCCPY = 1
                  IPX = 1
                  DDIAG = IPBCPY + DESCB( LLD_ ) * BCOLS 
                  SDIAG = DDIAG + MAX( SYDT_M_MAX, SYDT_N_MAX )
                  IPW = SDIAG + MAX( SYDT_M_MAX-1, SYDT_N_MAX-1 )
               END IF
C
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of A
C
               DO I = 1, SYCT_M_MAX
                  MEM(DDIAG+(I-1)) = DBLE(I)+100.0D+0
               END DO
C
C     Init A
C
               BLKS2B2 = INT(DBLE(SYDT_M_MAX/2) * 
     $                   DBLE(SYDT_A_NQTRBL)/100.0D+00)
               CALL P1SQMATGD( SYDT_ADIAG, 'No subdiagonal', SYDT_AFORM,
     $                         SYDT_M_MAX, MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DDIAG), MEM(SDIAG), BLKS2B2, 17, 
     $                         MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'SYDT: P1SQMATGD(A), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 60
               END IF
C
C     Copy A
C
               IF( .NOT. LSAME( SYDT_ASCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', SYDT_M_MAX, SYDT_M_MAX, MEM(IPA),
     $                          1, 1, DESCA, MEM(IPACPY), 1, 1, DESCA )
               END IF
C
C     Set eigenvalues of B
C
               DO I = 1, SYDT_N_MAX
                  MEM(DDIAG+(I-1)) = DBLE(I) + 100.0D+0
               END DO      
C     
C     Initialize B
C
               BLKS2B2 = INT(DBLE(SYDT_N_MAX/2) * 
     $                   DBLE(SYDT_B_NQTRBL)/100.0D+00)
               CALL P1SQMATGD( SYDT_BDIAG, 'No subdiagonal', SYDT_BFORM, 
     $                         SYDT_N_MAX, MEM(IPB), DESCB, ONE, ONE, 
     $                         MEM(DDIAG), MEM(SDIAG), BLKS2B2, 19,
     $                         MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'SYDT: P1SQMATGD(B), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 60
               END IF 
C
C     Copy B
C     
               IF( .NOT. LSAME( SYDT_BSCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', SYDT_N_MAX, SYDT_N_MAX, MEM(IPB),
     $                          1, 1, DESCB, MEM(IPBCPY), 1, 1, DESCB )
               END IF
C     
C     Init a matrix X as the given "solution"
C     
               IF( SYDT_SOLVE ) THEN
                  CALL PDMATGEN2( TEMP_ICTXT, SYDT_CFORM, SYDT_CDIAG, 
     $                            SYDT_M_MAX, SYDT_N_MAX, MB, NB, 
     $                            MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, 
     $                            AROWS, 0, BCOLS, MYROW, MYCOL, NPROW, 
     $                            NPCOL )
C     
C     Compute right hand side C as op(A) * X * op(B) +/- X
C     
                  CALL PDGEMM( SYDT_TRANSA, 'N', M, N, M, ONE, 
     $                         MEM(IPA), IA, JA, DESCA, MEM(IPX), IC, 
     $                         JC, DESCC, ZERO, MEM(IPCCPY), IC, JC, 
     $                         DESCC )
                  CALL PDGEMM( 'N', SYDT_TRANSB, M, N, N, ONE, 
     $                         MEM(IPCCPY), IC, JC, DESCC, MEM(IPB), IB, 
     $                         JB, DESCB, ZERO, MEM(IPC), IC, JC, 
     $                         DESCC )
                  CALL PDMATADD( M, N, DBLE(SYDT_ISGN), MEM(IPX), IC, 
     $                           JC, DESCC, ONE, MEM(IPC), IC, JC, 
     $                           DESCC )
C     
C     Copy C for future reference
C     
                  CALL PDLACPY( 'All', M, N, MEM(IPC), IC, JC, DESCC, 
     $                          MEM(IPCCPY), IC, JC, DESCC )
               END IF
C 
C     If the matrices are small enough, print them to the screen
C     
               IF(  IAM.EQ.0 .AND. M.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(A):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, M, MEM( IPA ), IA, JA, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF
C     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(B):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPB ), IB, JB, DESCB, 
     $                           0, 0, 'B', NOUT, MEM( IPW ) )
               END IF     
C
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. SYDT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X0):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. SYDT_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPX ), IC, JC, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
C
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. SYDT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(C):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. SYDT_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
C
C     Call SYDT solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( SYDT_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGESYDTD( 'Reduce', SYDT_ASCHUR, SYDT_BSCHUR, 
     $                           SYDT_TRANSA, SYDT_TRANSB, SYDT_ISGN, 
     $                           SYDT_COMM, M, N, MEM(IPA), IA, JA,
     $                           DESCA, MEM(IPB), IB, JB, DESCB, 
     $                           MEM(IPC), IC, JC,  DESCC, MB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( SYDT_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGESYDTD( 'Solve', SYDT_ASCHUR, SYDT_BSCHUR, 
     $                           SYDT_TRANSA, SYDT_TRANSB, SYDT_ISGN, 
     $                           SYDT_COMM, M, N, MEM(IPA), IA, JA,
     $                           DESCA, MEM(IPB), IB, JB, DESCB, 
     $                           MEM(IPC), IC, JC,  DESCC, MB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'SYDT: PGESYDTD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 60
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'SYDT: PGESYDTD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. SYDT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. SYDT_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Gflops/sec.
C
               GFLOPS = (DBLE(M)*DBLE(N)**2+DBLE(M)**2*DBLE(N)) /
     $              (10.0D+00 ** 9 * TPROF3)
C
C     In reduction mode only, move on to condition estimation
C
               IF( SYDT_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 50 
               END IF
C
C     Compute results of computation
C     
C                          ||op(A) * X * op(B) +/- X - SCALE * C||
C     Compute relres =  ---------------------------------------------
C                       ((||X||*(||A||*||B||) + SCALE * ||C||) * eps )
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
               CNORM = PDLANGE( 'I', M, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', M, N, MEM( IPC ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDGEMM( SYDT_TRANSA, 'N', M, N, M, ONE, 
     $                      MEM(IPACPY), IA, JA, DESCA, MEM(IPC), IC, 
     $                      JC, DESCC, ZERO, MEM(IPW), IC, JC, DESCC )
               CALL PDGEMM( 'N', SYDT_TRANSB, M, N, N, ONE, 
     $                      MEM(IPW), IC, JC, DESCC, MEM(IPBCPY), IB, 
     $                      JB, DESCB, -SCALE, MEM(IPCCPY), IC, JC, 
     $                      DESCC )
               CALL PDMATADD( M, N, DBLE(SYDT_ISGN), MEM(IPC), IC, 
     $                        JC, DESCC, ONE, MEM(IPCCPY), IC, JC, 
     $                        DESCC )
               IF( LSAME( SYDT_TRANSA, 'T' ) ) THEN
                  ANORM = PDLANGE( 'I', M, M, MEM( IPACPY ), IA, JA, 
     $                            DESCA, MEM( IPW ) )
               ELSE
                  ANORM = PDLANGE( '1', M, M, MEM( IPACPY ), IA, JA, 
     $                             DESCA, MEM( IPW ) )
               END IF
               IF( LSAME( SYDT_TRANSB, 'T' ) ) THEN
                  BNORM = PDLANGE( 'I', N, N, MEM( IPBCPY ), IB, JB, 
     $                             DESCB, MEM( IPW ) )
               ELSE
                  BNORM = PDLANGE( '1', N, N, MEM( IPBCPY ), IB, JB, 
     $                             DESCB, MEM( IPW ) )
               END IF
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = PDLANGE( 'F', M, N, MEM( IPCCPY ), IC, JC, 
     $                             DESCC, MEM(IPW) ) 
               RESID =  PDLANGE( 'I', M, N, MEM( IPCCPY ), IC, JC, 
     $                           DESCC, MEM(IPW)) / (((ANORM * BNORM) * 
     $                           XNORM + SCALE * CNORM)*EPS)
C     
C     Compute real error norms
C
               XNORM = PDLANGE( 'F', M, N, MEM( IPX ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDMATADD( M, N, SCALE, MEM(IPX), IC, JC, DESCC, 
     $                        -ONE, MEM(IPC), IC, JC, DESCC )
               ABSERROR = PDLANGE( 'F', M, N, MEM( IPC ), IC, JC, DESCC, 
     $                             MEM(IPW) )
               RELERROR = ABSERROR / XNORM
C
C     From reduction mode we arrive here
C
 50            CONTINUE
C
C     Condition estimation
C
               IF( SYDT_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                  CALL PSYDTCON( SYDT_TRANSA, SYDT_TRANSB, SYDT_ISGN,
     $                           SYDT_COMM, M, N, MEM( IPA ), IA, JA, 
     $                           DESCA, MEM( IPB ), IB, JB, DESCB, 
     $                           MB2, MEM( IPW ), WORKSIZ, IMEM, 
     $                           IMEMSIZ, EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'SYDT: PSYDTCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 60
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'SYDT: PSYDTCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7779 )
     $                 UPLOSIGN,M,N,MB,NB,MB2,NB,IA,IB,NPROW,NPCOL,
     $                 THREADS,MYTIME,TPROF1,TPROF2,TPROF3,GFLOPS,SCALE,
     $                 ABSERROR,RELERROR,ABSRESID,RESID,EST,ESTTIME,
     $                 CPROF1,CPROF2,CPROF3,NOITER,INFO
               END IF
C
            END DO
 55         CONTINUE
            END DO
 57         CONTINUE
            END DO
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 60         CONTINUE
C
C     Synchronize all processes.
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO
         END DO
#endif
         END DO
         END DO 
         END DO               
         END DO
         END DO
         END DO
C
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON LYDT ===
C
      IF( LYDT ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for LYDT ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6670 )
     $           'ULS','N','NB','NB2','IA','P_r','P_c','Thr','Tme',
     $           'Tp1','Tp2','Tp3','Gfl','Scl','E_a','E_r','R_a','R_r',
     $           'est','etm','Cp1','Cp2','Cp3','itr','INF'
         END IF
C
C     Loop though the space of all combinations of parameters,
C     N, NB, NPROW, NPCOL
C
         DO UPLOSIGN = LYDT_UPLOSIGN_MIN, LYDT_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 ) THEN
            LYDT_TRANSA = 'N'
         ELSEIF( UPLOSIGN.EQ.1 ) THEN
            LYDT_TRANSA = 'T'
         END IF
         DO N = LYDT_N_MIN, LYDT_N_MAX, LYDT_N_STEP
         DO NB = LYDT_NB_MIN, LYDT_NB_MAX, LYDT_NB_STEP
         DO NB2 = LYDT_NB2_MIN, LYDT_NB2_MAX, LYDT_NB2_STEP 
#ifdef LOOPGRID
         DO TEMP_NPROW = LYDT_NPROW_MIN, LYDT_NPROW_MAX, LYDT_NPROW_STEP
         DO TEMP_NPCOL = LYDT_NPCOL_MIN, LYDT_NPCOL_MAX, LYDT_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C                  
#ifdef LOOPGRID  
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $          MYCOL.LT.0 .OR. MYCOL.GE.NPCOL  .OR.
     $          (LYDT_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) ) GO TO 80
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            AROWS = NUMROC( LYDT_N_MAX, NB, MYROW, 0, NPROW )
            ACOLS = NUMROC( LYDT_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, LYDT_N_MAX, LYDT_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'LYDT: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 80
            END IF 
            IF( LYDT_SOLVE )
     $           CALL DESCINIT( DESCC, LYDT_N_MAX, LYDT_N_MAX, NB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, AROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'LYDT: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 80
            END IF 
C
C     Loop through the space of parameter IA, JA and
C     ignore invalid cases
C
            DO IA = LYDT_IA_MIN, LYDT_IA_MAX, LYDT_IA_STEP
            JA = IA
            IC = IA
            JC = JA
            IF( .NOT. ( IA + N - 1 ).LE.LYDT_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IA = ',IA
               GO TO 75
            END IF
            DO REPEAT = 1, LYDT_REPEAT, 1
C
C     Assign pointer for ScaLAPACK arrays
C
               IPA = 1
               IF( .NOT. LSAME( LYDT_ASCHUR, 'S') ) THEN
                  IPACPY = IPA + DESCA( LLD_ ) * ACOLS
               ELSE
                  IPACPY = IPA
               END IF
               IF( LYDT_SOLVE ) THEN
                  IPC = IPACPY + DESCA( LLD_ ) * ACOLS
                  IPCCPY = IPC + DESCC( LLD_ ) * ACOLS
                  IPX = IPCCPY + DESCC( LLD_ ) * ACOLS
                  DDIAG = IPX + DESCC( LLD_ ) * ACOLS 
                  SDIAG = DDIAG + LYDT_N_MAX 
                  IPW = SDIAG + LYDT_N_MAX-1
               ELSE
                  IPC = 1
                  IPCCPY = 1
                  IPX = 1
                  DDIAG = IPACPY + DESCA( LLD_ ) * ACOLS 
                  SDIAG = DDIAG + LYDT_N_MAX 
                  IPW = SDIAG + LYDT_N_MAX-1
               END IF
C     
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of A
C
               DO I = 1, LYDT_N_MAX
                  MEM(DDIAG+(I-1)) = 1000.0D+00 + DBLE(I)
               END DO
C
C     Init A
C
               BLKS2B2 = INT(DBLE(LYDT_N_MAX/2) * 
     $                   DBLE(LYDT_A_NQTRBL)/100.0D+00)
               CALL P1SQMATGD( LYDT_ADIAG, 'No subdiagonal', LYDT_AFORM,
     $                         LYDT_N_MAX, MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DDIAG), MEM(SDIAG), BLKS2B2, 17,
     $                         MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'LYDT: P1SQMATGD(A), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 80
               END IF 
C
C     Copy A
C
               IF( .NOT. LSAME( LYDT_ASCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', LYDT_N_MAX, LYDT_N_MAX, MEM(IPA),
     $                          1, 1, DESCA, MEM(IPACPY), 1, 1, DESCA )
               END IF
C     
C     Init a matrix X as the given solution
C     
               IF( LYDT_SOLVE ) THEN
                  IF( LSAME( LYDT_C_SYMM, 'N' ) ) THEN
                     CALL PDMATGEN2( TEMP_ICTXT, 'Symmetric', 
     $                    LYDT_CDIAG, LYDT_N_MAX, LYDT_N_MAX, NB, NB, 
     $                    MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, AROWS, 0, 
     $                    ACOLS, MYROW, MYCOL, NPROW, NPCOL )
                  ELSE
                     CALL PDLASET( 'All', LYDT_N_MAX, LYDT_N_MAX, 
     $                    ZERO, ONE, MEM(IPX), IC, JC, DESCC )
                  END IF
C
C     Compute right hand side C as op(A) * X * op(A') - X
C
                  IF( LSAME( LYDT_TRANSA, 'N' ) ) THEN 
                     CALL PDGEMM( 'N', 'N', N, N, N, ONE, MEM(IPA), 
     $                    IA, JA, DESCA, MEM(IPX), IC, JC, DESCC, 
     $                    ZERO, MEM(IPCCPY), IC, JC, DESCC )
                     CALL PDGEMM( 'N', 'Transpose', N, N, N, ONE, 
     $                    MEM(IPCCPY), IC, JC, DESCC, MEM(IPA), 
     $                    IA, JA, DESCA, ZERO, MEM(IPC), IC, JC, 
     $                    DESCC )
                     CALL PDMATADD( N, N, -ONE, MEM(IPX), IC, JC, DESCC, 
     $                    ONE, MEM(IPC), IC, JC, DESCC )
                  ELSEIF( LSAME( LYDT_TRANSA, 'T' ) ) THEN 
                     CALL PDGEMM( 'Transpose', 'N', N, N, N, ONE, 
     $                    MEM(IPA), IA, JA, DESCA, MEM(IPX), IC,
     $                    JC, DESCC, ZERO, MEM(IPCCPY), IC, JC, 
     $                    DESCC )
                     CALL PDGEMM( 'N', 'N', N, N, N, ONE, 
     $                    MEM(IPCCPY), IC, JC, DESCC, MEM(IPA), IA, 
     $                    JA, DESCA, ZERO, MEM(IPC), IC, JC, DESCC )
                     CALL PDMATADD( N, N, -ONE, MEM(IPX), IC, JC, DESCC, 
     $                    ONE, MEM(IPC), IC, JC, DESCC )
                  END IF
C     
C     Copy C for future reference
C
                  CALL PDLACPY( 'All', N, N, MEM(IPC), IC, JC, DESCC, 
     $                          MEM(IPCCPY), IC, JC, DESCC )
               END IF
C
C     Check if sub(C) is symmetric for this particular instance
C
               IF( LSAME( LYDT_C_SYMM, 'S' ) .AND. IC.EQ.JC ) THEN
                  LYDT_C_SYMM2 = 'S'
               ELSE
                  LYDT_C_SYMM2 = 'N'
               END IF
C 
C     If the matrices are small enough, print them to the screen
C     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(A):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPA ), IA, JA, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF
C     
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. LYDT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(C):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. LYDT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
                IF( N.LE.PRNTSIZ .AND. LYDT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPX ), IC, JC, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
C
C     Call LYDT solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( LYDT_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGELYDTD( 'Reduce', LYDT_C_SYMM2, LYDT_TRANSA,
     $                           LYDT_ASCHUR, N, MEM(IPA), IA, JA, 
     $                           DESCA, MEM(IPC), IC, JC, DESCC, NB2, 
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( LYDT_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGELYDTD( 'Solve', LYDT_C_SYMM2, LYDT_TRANSA,
     $                           LYDT_ASCHUR, N, MEM(IPA), IA, JA, 
     $                           DESCA, MEM(IPC), IC, JC, DESCC, NB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'LYDT: PGELYDTD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 80
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'LYDT: PGELYDTD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. LYDT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. LYDT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Gflops/sec.
C
               IF( LSAME( LYDT_C_SYMM, 'S' ) ) THEN
                  GFLOPS = (1.0D+00*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               ELSE
                  GFLOPS = (2.0D+00*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               END IF
C
C     In reduction mode only, move on to condition estimation
C
               IF( LYDT_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 70 
               END IF
C
C     Compute results of computation
C     
C                          ||op(A) * X * op(A') - X - SCALE * C||
C     Compute relres =  ---------------------------------------------
C                       ((||X||*(||A||**2) + SCALE * ||C||) * eps )
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
               CNORM = PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               IF( LSAME( LYDT_TRANSA, 'N' ) ) THEN
                  CALL PDGEMM( 'N', 'N', N, N, N, ONE, MEM(IPACPY), IA, 
     $                         JA, DESCA, MEM(IPC), IC, JC, DESCC, 
     $                         ZERO, MEM(IPW), IC, JC, DESCC )
                  CALL PDGEMM( 'N', 'Transpose', N, N, N, ONE, 
     $                         MEM(IPW), IC, JC, DESCC, MEM(IPACPY), 
     $                         IA, JA, DESCA, -SCALE, MEM(IPCCPY), IC, 
     $                         JC, DESCC )
                  CALL PDMATADD( N, N, -ONE, MEM(IPC), IC, JC, DESCC, 
     $                           ONE, MEM(IPCCPY), IC, JC, DESCC )
               ELSEIF( LSAME( LYDT_TRANSA, 'T' ) ) THEN
                  CALL PDGEMM( 'Transpose', 'N', N, N, N, ONE, 
     $                         MEM(IPACPY), IA, JA, DESCA, MEM(IPC), IC,
     $                         JC, DESCC, ZERO, MEM(IPW), IC, JC, 
     $                         DESCC )
                  CALL PDGEMM( 'N', 'N', N, N, N, ONE, MEM(IPW), IC, JC, 
     $                         DESCC, MEM(IPACPY), IA, JA, DESCA, 
     $                         -SCALE, MEM(IPCCPY), IC, JC, DESCC )
                  CALL PDMATADD( N, N, -ONE, MEM(IPC), IC, JC, DESCC, 
     $                           ONE, MEM(IPCCPY), IC, JC, DESCC )
               END IF
               IF( LSAME( LYDT_TRANSA, 'T' ) ) THEN
                  ANORM = PDLANGE( 'I', N, N, MEM( IPACPY ), IA, JA, 
     $                            DESCA, MEM( IPW ) )
               ELSE
                  ANORM = PDLANGE( '1', N, N, MEM( IPACPY ), IA, JA, 
     $                             DESCA, MEM( IPW ) )
               END IF
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = PDLANGE( 'F', N, N, MEM( IPCCPY ), IC, JC, 
     $                             DESCC, MEM(IPW) ) 
               RESID =  PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, 
     $                           DESCC, MEM(IPW)) / (((ANORM * ANORM) * 
     $                           XNORM + SCALE * CNORM)*EPS)
C     
C     Compute real error norms
C
               XNORM = PDLANGE( 'F', N, N, MEM( IPX ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDMATADD( N, N, SCALE, MEM(IPX), IC, JC, DESCC, 
     $                        -ONE, MEM(IPC), IC, JC, DESCC )
               ABSERROR = PDLANGE( 'F', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                             MEM(IPW) )
               RELERROR = ABSERROR / XNORM
C
C     From reduction mode we arrive here
C
 70            CONTINUE
C
C     Condition estimation
C
               IF( LYDT_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                  CALL PLYDTCON( LYDT_TRANSA, N, MEM( IPA ), IA, JA,
     $                           DESCA, NB2, MEM( IPW ), WORKSIZ, IMEM, 
     $                           IMEMSIZ, EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'LYDT: PLYDTCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 80
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'LYDT: PLYDTCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
              
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7780 )
     $                 UPLOSIGN,N,NB,NB2,IA,NPROW,NPCOL,THREADS,MYTIME,
     $                 TPROF1,TPROF2,TPROF3,GFLOPS,SCALE,ABSERROR,
     $                 RELERROR,ABSRESID,RESID,EST,ESTTIME,CPROF1,
     $                 CPROF2,CPROF3,NOITER,INFO
               END IF
C
            END DO 
 75         CONTINUE
            END DO
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 80         CONTINUE
C
C     Synchronize all processes. 
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO 
         END DO
#endif
         END DO               
         END DO
         END DO
         END DO
C
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON GCSY ===
C
      IF( GCSY ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for GCSY ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6669 )
     $           'ULS','M','N','MB','NB','MB2','NB2','IA','IB','P_r',
     $           'P_c','Thr','Tme','Tp1','Tp2','Tp3','Gfl','Scl','E_a',
     $           'E_r','R_a','R_r','est','etm','Cp1','Cp2','Cp3','itr',
     $           'INF'
         END IF
C
C     Loop though the space of all combinations of parameters,
C     M, N, MB, NB, NPROW, NPCOL
C
         DO UPLOSIGN = GCSY_UPLOSIGN_MIN, GCSY_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 .OR. UPLOSIGN.EQ.1 ) THEN
            GCSY_TRANSAD = 'N'
            GCSY_TRANSBE = 'N'
         ELSEIF( UPLOSIGN.EQ.2 .OR. UPLOSIGN.EQ.3 ) THEN
            GCSY_TRANSAD = 'N'
            GCSY_TRANSBE = 'T'
         ELSEIF( UPLOSIGN.EQ.4 .OR. UPLOSIGN.EQ.5 ) THEN
            GCSY_TRANSAD = 'T'
            GCSY_TRANSBE = 'N'
         ELSEIF( UPLOSIGN.EQ.6 .OR. UPLOSIGN.EQ.7 ) THEN
            GCSY_TRANSAD = 'T'
            GCSY_TRANSBE = 'T'
         END IF
         GCSY_ISGN = INT(1.0D+00-(MOD(UPLOSIGN,2))*2.0D+00)
         DO M = GCSY_M_MIN, GCSY_M_MAX, GCSY_M_STEP
         DO N = GCSY_N_MIN, GCSY_N_MAX, GCSY_N_STEP
         DO MB = GCSY_MB_MIN, GCSY_MB_MAX, GCSY_MB_STEP
         DO NB = GCSY_NB_MIN, GCSY_NB_MAX, GCSY_NB_STEP
         DO MB2 = GCSY_MB2_MIN, GCSY_MB2_MAX, GCSY_MB2_STEP
         DO NB2 = GCSY_NB2_MIN, GCSY_NB2_MAX, GCSY_NB2_STEP
C          MB2 = MB
C          NB2 = NB
#ifdef LOOPGRID
         DO TEMP_NPROW = GCSY_NPROW_MIN, GCSY_NPROW_MAX, GCSY_NPROW_STEP
         DO TEMP_NPCOL = GCSY_NPCOL_MIN, GCSY_NPCOL_MAX, GCSY_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C           
#ifdef LOOPGRID
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid OR
C     skip this run for the current parameter set
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $           MYCOL.LT.0 .OR. MYCOL.GE.NPCOL .OR.
     $           (GCSY_SQUARE_DIM .AND. M.NE.N) .OR.
     $           (GCSY_SQUARE_BLK .AND. MB.NE.NB)  .OR.
     $           (GCSY_SQUARE_BLK .AND. MB2.NE.NB2) .OR.
     $           (GCSY_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2)) .OR.
     $           (.NOT.(UPLOSIGN.LE.1 .OR. UPLOSIGN.GE.6) .AND. 
     $           LSAME(GCSY_TRANSZ,'T')) .OR.
     $           (.NOT.GCSY_PIPELINE .AND. 
     $           (MB.NE.MB2.OR.NB.NE.NB2) ) ) GO TO 100
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            ADROWS = NUMROC( GCSY_M_MAX, MB, MYROW, 0, NPROW )
            ADCOLS = NUMROC( GCSY_M_MAX, MB, MYCOL, 0, NPCOL )
            BEROWS = NUMROC( GCSY_N_MAX, NB, MYROW, 0, NPROW )
            BECOLS = NUMROC( GCSY_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, GCSY_M_MAX, GCSY_M_MAX, MB, MB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, ADROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 100
            END IF
            CALL DESCINIT( DESCD, GCSY_M_MAX, GCSY_M_MAX, MB, MB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, ADROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY: DESCINIT(D), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 100
            END IF
            CALL DESCINIT( DESCB, GCSY_N_MAX, GCSY_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, BEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY: DESCINIT(B), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 100
            END IF 
            CALL DESCINIT( DESCE, GCSY_N_MAX, GCSY_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, BEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY: DESCINIT(E), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 100
            END IF
            IF( GCSY_SOLVE )
     $           CALL DESCINIT( DESCC, GCSY_M_MAX, GCSY_N_MAX, MB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, ADROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 100
            END IF 
            IF( GCSY_SOLVE )
     $           CALL DESCINIT( DESCF, GCSY_M_MAX, GCSY_N_MAX, MB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, ADROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GCSY: DESCINIT(F), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 100
            END IF
C
C     Set internal blocking factors
C
            MBNB2(1) = MB2
            MBNB2(2) = NB2 
C
C     Loop through the space of parameter IAD, JAD, IBE, JBE, ICF, JCF 
C     and ignore invalid cases
C
            DO IAD = GCSY_IAD_MIN, GCSY_IAD_MAX, GCSY_IAD_STEP
            JAD = IAD
            ICF = IAD
            IF( .NOT. ( IAD + M - 1 ).LE.GCSY_M_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IAD = ',IAD
               GO TO 99
            END IF
            DO IBE = GCSY_IBE_MIN, GCSY_IBE_MAX, GCSY_IBE_STEP
            JBE = IBE
            JCF = JBE
            IF( .NOT. ( IBE + N - 1 ).LE.GCSY_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IBE = ',IBE
               GO TO 95
            END IF
            DO REPEAT = 1, GCSY_REPEAT, 1
C
C     Assign pointers for ScaLAPACK arrays
C
                IPA = 1
                IF( LSAME( GCSY_ADSCHUR, 'S' ) ) THEN
                   IPACPY = IPA
                ELSE
                   IPACPY = IPA + DESCA( LLD_) * ADCOLS
                END IF
                IPB = IPACPY + DESCA( LLD_) * ADCOLS
                IF( LSAME( GCSY_BESCHUR, 'S' ) ) THEN
                   IPBCPY = IPB
                ELSE
                   IPBCPY = IPB + DESCB( LLD_ ) * BECOLS
                END IF
                IPD = IPBCPY +  DESCB( LLD_ ) * BECOLS
                IF( LSAME( GCSY_ADSCHUR, 'S' ) ) THEN
                   IPDCPY = IPD
                ELSE
                   IPDCPY = IPD + DESCD( LLD_) * ADCOLS
                END IF
                IPE = IPDCPY + DESCD( LLD_ ) * ADCOLS
                IF( LSAME( GCSY_BESCHUR, 'S' ) ) THEN
                   IPECPY = IPE
                ELSE
                   IPECPY = IPE + DESCE( LLD_ ) * BECOLS
                END IF
                IF( GCSY_SOLVE ) THEN
                   IPC = IPECPY + DESCE( LLD_ ) * BECOLS
                   IPCCPY = IPC + DESCC( LLD_ ) * BECOLS
                   IPF = IPCCPY + DESCC( LLD_ ) * BECOLS
                   IPFCPY = IPF + DESCF( LLD_ ) * BECOLS
                   IPX = IPFCPY + DESCF( LLD_ ) * BECOLS
                   IPY = IPX + DESCC( LLD_ ) * BECOLS
                   DIAG1 = IPY + DESCF( LLD_ ) * BECOLS 
                   SDIAG = DIAG1 + MAX( GCSY_M_MAX-1, GCSY_N_MAX-1 )
                   DIAG2 = SDIAG + MAX( GCSY_M_MAX, GCSY_N_MAX )
                   IPW = DIAG2 + MAX( GCSY_M_MAX, GCSY_N_MAX )
                ELSE
                   IPC = 1
                   IPCCPY = 1
                   IPF = 1
                   IPFCPY = 1
                   IPX = 1
                   IPY = 1
                   DIAG1 = IPECPY + DESCE( LLD_ ) * BECOLS 
                   SDIAG = DIAG1 + MAX( GCSY_M_MAX-1, GCSY_N_MAX-1 )
                   DIAG2 = SDIAG + MAX( GCSY_M_MAX, GCSY_N_MAX )
                   IPW = DIAG2 + MAX( GCSY_M_MAX, GCSY_N_MAX )
                END IF
C     
C     Set some variables
C
               WORKSIZ = MEMSIZ - IPW + 1
               NOEXTSYS = 0
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of (A,D)
C
               DO I = 1, GCSY_M_MAX
                  MEM(DIAG1+(I-1)) = DBLE(I) + 1000.0D+00
                  MEM(DIAG2+(I-1)) = ONE
               END DO
C
C     Init (A,D)
C
               BLKS2B2 = INT(DBLE(GCSY_M_MAX/2) * 
     $                   DBLE(GCSY_AD_NQTRBL)/100.0D+00)
               CALL P2SQMATGD( GCSY_ADDIAG, 'No subdiagonal', 
     $                         GCSY_ADDIAG, GCSY_ADFORM, GCSY_M_MAX, 
     $                         MEM(IPA), DESCA, ONE, ONE / DBLE(M), 
     $                         MEM(DIAG1), MEM(SDIAG), BLKS2B2, 17,
     $                         MEM(IPD), DESCD, ONE, ONE / DBLE(M), 
     $                         MEM(DIAG2), 19, MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GCSY: P2SQMATGD(A,D), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 100
               END IF 
C
C     Copy (A,D)
C
               IF( .NOT. LSAME( GCSY_ADSCHUR, 'S' ) ) THEN
                   CALL PDLACPY( 'All', GCSY_M_MAX, GCSY_M_MAX, 
     $                           MEM(IPA), 1, 1, DESCA, MEM(IPACPY),
     $                           1, 1, DESCA )
                   CALL PDLACPY( 'All', GCSY_M_MAX, GCSY_M_MAX, 
     $                           MEM(IPD), 1, 1, DESCD, MEM(IPDCPY),
     $                           1, 1, DESCD )
                END IF
C
C     Set eigenvalues of (B,E)
C
               DO I = 1, GCSY_N_MAX
                  MEM(DIAG1+(I-1)) = DBLE(I) + 1000.0D+00
                  MEM(DIAG2+(I-1)) = -DBLE(GCSY_ISGN)
               END DO      
C     
C     Initialize (B,E)
C
               BLKS2B2 = INT(DBLE(GCSY_N_MAX/2) * 
     $                   DBLE(GCSY_BE_NQTRBL)/100.0D+00)
               CALL P2SQMATGD( GCSY_BEDIAG, 'No subdiagonal', 
     $                         GCSY_BEDIAG, GCSY_BEFORM, GCSY_N_MAX, 
     $                         MEM(IPB), DESCB, ONE, ONE / DBLE(N), 
     $                         MEM(DIAG1), MEM(SDIAG), BLKS2B2, 23, 
     $                         MEM(IPE), DESCE, ONE, ONE / DBLE(N), 
     $                         MEM(DIAG2), 31, MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GCSY: P2SQMATGD(B,E), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 100
               END IF 
C
C     Copy B
C     
               IF( .NOT. LSAME( GCSY_BESCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', GCSY_N_MAX, GCSY_N_MAX, MEM(IPB),
     $                          1, 1, DESCB, MEM(IPBCPY), 1, 1, DESCB )
                  CALL PDLACPY( 'All', GCSY_N_MAX, GCSY_N_MAX, MEM(IPE),
     $                          1, 1, DESCE, MEM(IPECPY), 1, 1, DESCE )
               END IF
C     
C     Init a matrix pair (X,Y) as the given "solution"
C     
               IF( GCSY_SOLVE ) THEN
                  CALL PDMATGEN2( TEMP_ICTXT, GCSY_CFFORM, GCSY_CFDIAG, 
     $                            GCSY_M_MAX, GCSY_N_MAX, MB, NB, 
     $                            MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, 
     $                            ADROWS, 0, BECOLS, MYROW, MYCOL, 
     $                            NPROW, NPCOL )
                  CALL PDMATGEN2( TEMP_ICTXT, GCSY_CFFORM, GCSY_CFDIAG, 
     $                            GCSY_M_MAX, GCSY_N_MAX, MB, NB, 
     $                            MEM(IPY), DESCC(LLD_), 0, 0, 13, 0, 
     $                            ADROWS, 0, BECOLS, MYROW, MYCOL, 
     $                            NPROW, NPCOL )
C     
C     Compute C as op(A) * X +/- Y * op(B) or op(A^T) * X +/- op(D^T) * Y 
C     Compute F as op(D) * X +/- Y * op(E) or X * op(B^T) +/- Y * op(E^T)
C     
                  IF( .NOT. LSAME( GCSY_TRANSZ, 'T') ) THEN
                     CALL PDGEMM( GCSY_TRANSAD, 'N' , M, N, M, ONE, 
     $                            MEM(IPA), IAD, JAD, DESCA, MEM(IPX), 
     $                            ICF, JCF, DESCC, ZERO, MEM(IPC), 
     $                            ICF, JCF, DESCC )
                     CALL PDGEMM( 'N', GCSY_TRANSBE, M, N, N, 
     $                            DBLE(GCSY_ISGN), MEM(IPY), ICF, JCF,
     $                            DESCC, MEM(IPB), IBE, JBE, DESCB, ONE, 
     $                            MEM(IPC), ICF, JCF, DESCC )
                     CALL PDGEMM( GCSY_TRANSAD, 'N' , M, N, M, ONE, 
     $                            MEM(IPD), IAD, JAD, DESCD, MEM(IPX), 
     $                            ICF, JCF, DESCC, ZERO, MEM(IPF), 
     $                            ICF, ICF, DESCF )
                     CALL PDGEMM( 'N', GCSY_TRANSBE, M, N, N, 
     $                            DBLE(GCSY_ISGN), MEM(IPY), ICF, JCF,
     $                            DESCC, MEM(IPE), IBE, JBE, DESCE, ONE, 
     $                            MEM(IPF), ICF, JCF, DESCF )
                  ELSE
                     NOP = 'N'
                     IF( LSAME(GCSY_TRANSAD,'N') ) NOP = 'T'
                     CALL PDGEMM( NOP, 'N', M, N, M, ONE, MEM(IPA), IAD, 
     $                           JAD, DESCA, MEM(IPX), ICF, JCF, DESCC, 
     $                           ZERO, MEM(IPC), ICF, JCF, DESCC )
                     CALL PDGEMM( NOP, 'N', M, N, M, ONE, MEM(IPD), IAD, 
     $                            JAD, DESCC, MEM(IPY), ICF, JCF, DESCF, 
     $                            ONE, MEM(IPC), ICF, JCF, DESCC )
                     NOP = 'N'
                     IF( LSAME(GCSY_TRANSBE,'N') ) NOP = 'T'
                     CALL PDGEMM( 'N', NOP, M, N, N, DBLE(GCSY_ISGN), 
     $                            MEM(IPX), ICF, JCF, DESCC, MEM(IPB), 
     $                            IBE, JBE, DESCB, ZERO, MEM(IPF), 
     $                            ICF, JCF, DESCF )
                     CALL PDGEMM( 'N', NOP, M, N, N, DBLE(GCSY_ISGN), 
     $                            MEM(IPY), ICF, JCF, DESCF, MEM(IPE), 
     $                            IBE, JBE, DESCE, ONE, MEM(IPF),
     $                            ICF, JCF, DESCF )
                  END IF
C     
C     Copy (C,F) for future reference
C     
                  CALL PDLACPY( 'All', M, N, MEM(IPC), ICF, JCF, DESCC, 
     $                          MEM(IPCCPY), ICF, JCF, DESCC )
                  CALL PDLACPY( 'All', M, N, MEM(IPF), ICF, JCF, DESCF, 
     $                          MEM(IPFCPY), ICF, JCF, DESCF )
               END IF
C     
C     If the matrices are small enough, print them to output
C     
               IF(  IAM.EQ.0 .AND. M.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix A:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, M, MEM( IPA ), IAD, JAD, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix B:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPB ), IBE, JBE, DESCB, 
     $                           0, 0, 'B', NOUT, MEM( IPW ) )
               END IF          
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GCSY_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix C:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GCSY_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPC ), ICF, JCF, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
               IF(  IAM.EQ.0 .AND. M.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix D:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, M, MEM( IPD ), IAD, JAD, DESCD, 
     $                           0, 0, 'D', NOUT, MEM( IPW ) )
               END IF     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix E:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPE ), IBE, JBE, DESCE, 
     $                           0, 0, 'E', NOUT, MEM( IPW ) )
               END IF          
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GCSY_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix F:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GCSY_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPF ), ICF, JCF, DESCF, 
     $                           0, 0, 'F', NOUT, MEM( IPW ) )
               END IF
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GCSY_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X0):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GCSY_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPX ), ICF, JCF, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GCSY_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(Y0):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GCSY_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPY ), ICF, JCF, DESCF, 
     $                           0, 0, 'Y0', NOUT, MEM( IPW ) )
               END IF
C     
C     Call GCSY solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( GCSY_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGCSYD( 'Reduce', GCSY_TRANSZ, GCSY_ADSCHUR, 
     $                           GCSY_BESCHUR, GCSY_TRANSAD, 
     $                           GCSY_TRANSBE, GCSY_ISGN, GCSY_COMM, 
     $                           M, N, MEM(IPA), IAD, JAD, DESCA, 
     $                           MEM(IPB), IBE, JBE, DESCB, MEM(IPC), 
     $                           ICF, JCF, DESCC, MEM(IPD), IAD, JAD, 
     $                           DESCD, MEM(IPE), IBE, JBE, DESCE, 
     $                           MEM(IPF), ICF, JCF, DESCF, MBNB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( GCSY_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGCSYD( 'Solve', GCSY_TRANSZ, GCSY_ADSCHUR, 
     $                           GCSY_BESCHUR, GCSY_TRANSAD, 
     $                           GCSY_TRANSBE, GCSY_ISGN, GCSY_COMM, 
     $                           M, N, MEM(IPA), IAD, JAD, DESCA, 
     $                           MEM(IPB), IBE, JBE, DESCB, MEM(IPC), 
     $                           ICF, JCF, DESCC, MEM(IPD), IAD, JAD, 
     $                           DESCD, MEM(IPE), IBE, JBE, DESCE, 
     $                           MEM(IPF), ICF, JCF, DESCF, MBNB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GCSY: PGEGCSYD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 100
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'GCSY: PGEGCSYD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrices if possible
C
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GCSY_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GCSY_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPC ), ICF, JCF, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GCSY_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(Y):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GCSY_SOLVE ) 
     $             THEN
                  CALL PDLAPRNT( M, N, MEM( IPF ), ICF, JCF, DESCF, 
     $                           0, 0, 'Y', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Gflops/sec.
C
               GFLOPS = 2.0D+00*(DBLE(M)*DBLE(N)**2+DBLE(M)**2*DBLE(N))/
     $              (10.0D+00 ** 9 * TPROF3)
C
C     In reduction mode only, move on to condition estimation
C
               IF( GCSY_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 90 
               END IF
C
C     Compute results of computation
C     
C     Compute the relative residual for GCSY as
C
C  ||(op(A)*X  +/- Y*op(B) - SCALE*C,op(D)*X  +/- Y*op(E) - SCALE*F)||_F
C  ---------------------------------------------------------------------
C  ((||(A,D)||_F + ||(B,E)||_F)*||(X,Y)||_F + SCALE*||(C,F)||_F ) * eps 
C          
C               
C     Find machine epsilon
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
C
C     Compute the norms of (C,F), (X,Y), (A,D) and (B,E)
C
               CNORM = PDLANGE( 'F', M, N, MEM( IPCCPY ), ICF, JCF, 
     $                          DESCC, MEM(IPW) )
               XNORM = PDLANGE( 'F', M, N, MEM( IPC ), ICF, JCF, DESCC, 
     $                          MEM( IPW ) )
               FNORM = PDLANGE( 'F', M, N, MEM( IPFCPY ), ICF, JCF, 
     $                          DESCF, MEM(IPW) )
               YNORM = PDLANGE( 'F', M, N, MEM( IPF ), ICF, JCF, DESCF, 
     $                          MEM( IPW ) )
               XYNORM = DLAPY2( XNORM, YNORM )
               CFNORM = DLAPY2( CNORM, FNORM )
C
               ANORM = PDLANGE( 'F', M, M, MEM( IPACPY ), IAD, JAD, 
     $                          DESCA, MEM( IPW ) )
               BNORM = PDLANGE( 'F', N, N, MEM( IPBCPY ), IBE, JBE, 
     $                          DESCB, MEM( IPW ) )
               DNORM = PDLANGE( 'F', M, M, MEM( IPDCPY ), IAD, JAD, 
     $                          DESCD, MEM( IPW ) )
               ENORM = PDLANGE( 'F', N, N, MEM( IPECPY ), IBE, JBE, 
     $                          DESCE, MEM( IPW ) )
               ADNORM = DLAPY2( ANORM, DNORM )
               BENORM = DLAPY2( BNORM, ENORM )
C     
C     Do the matrix multiplications from the relative residual and 
C     compute the associated norms
C     
               IF( LSAME( GCSY_TRANSZ, 'T') ) GO TO 97
               CALL PDGEMM( GCSY_TRANSAD, 'N', M, N, M, -ONE, 
     $                      MEM( IPACPY), IAD, JAD, DESCA, MEM( IPC ), 
     $                      ICF, JCF, DESCC, SCALE, MEM( IPCCPY ), 
     $                      ICF, JCF, DESCC )
               CALL PDGEMM( 'N', GCSY_TRANSBE, M, N, N, 
     $                      DBLE(-1*GCSY_ISGN), MEM(IPF), ICF, JCF, 
     $                      DESCF, MEM( IPBCPY ), IBE, JBE, DESCB, ONE,
     $                      MEM( IPCCPY ), ICF, JCF, DESCC )
               CALL PDGEMM( GCSY_TRANSAD, 'N', M, N, M, -ONE, 
     $                      MEM( IPDCPY), IAD, JAD, DESCD, MEM( IPC ), 
     $                      ICF, JCF, DESCC, SCALE, MEM( IPFCPY ), 
     $                      ICF, JCF, DESCF )
               CALL PDGEMM( 'N', GCSY_TRANSBE, M, N, N, 
     $                      DBLE(-1*GCSY_ISGN), MEM(IPF), ICF, JCF, 
     $                      DESCF, MEM( IPECPY ), IBE, JBE, DESCE, ONE,
     $                      MEM( IPFCPY ), ICF, JCF, DESCF )
 96            GO TO 98
 97            CONTINUE
               NOP = 'N'
               IF( LSAME( GCSY_TRANSAD, 'N' ) ) NOP = 'T'
               CALL PDGEMM( NOP, 'N', M, N, M, -ONE, MEM( IPACPY), IAD,
     $                      JAD, DESCA,MEM( IPC ), ICF, JCF, DESCC, 
     $                      SCALE, MEM( IPCCPY ), ICF, JCF, DESCC )
               CALL PDGEMM( NOP, 'N', M, N, M, -ONE, MEM( IPDCPY), IAD,
     $                      JAD, DESCD, MEM( IPF ), ICF, JCF, DESCF, 
     $                      ONE, MEM( IPCCPY ), ICF, JCF, DESCC )
               NOP = 'N'
               IF( LSAME( GCSY_TRANSBE,'N' ) ) NOP = 'T'
               CALL PDGEMM( 'N', NOP, M, N, N, DBLE(-1*GCSY_ISGN), 
     $                      MEM(IPC), ICF, JCF, DESCC, MEM( IPBCPY ), 
     $                      IBE, JBE, DESCB, SCALE, MEM( IPFCPY ), 
     $                      ICF, JCF, DESCC )
               CALL PDGEMM( 'N', NOP, M, N, N, DBLE(-1*GCSY_ISGN), 
     $                      MEM(IPF), ICF, JCF, DESCF, MEM( IPECPY ), 
     $                      IBE, JBE, DESCE, ONE, MEM( IPFCPY ), 
     $                      ICF, JCF, DESCF )
 98            CONTINUE
C
C     Compute the residual norms ||(R1,R2)||_F
C
               R1NORM = PDLANGE( 'F', M, N, MEM( IPCCPY ), ICF, JCF, 
     $                           DESCC, MEM( IPW ) )
               R2NORM = PDLANGE( 'F', M, N, MEM( IPFCPY ), ICF, JCF, 
     $                           DESCF, MEM( IPW ) )
               RNORM = DLAPY2( R1NORM, R2NORM )
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = RNORM
               RESID =  RNORM / 
     $              (((ADNORM+BENORM)*XYNORM + SCALE*CFNORM)*EPS)
C     
C     Compute forward error norms
C     
C     forward = ||(X-XX,Y-YY)||_F
C     relforw = ||(X-XX,Y-YY)||_F / ||(X,Y)||_F
C     
C     where X and Y is the exact solution and XX and YY is the computed
C     ones
C     
               XNORM = PDLANGE( 'F', M, N, MEM( IPX ), ICF, JCF, DESCC, 
     $                          MEM( IPW ) )
               YNORM = PDLANGE( 'F', M, N, MEM( IPY ), ICF, JCF, DESCF, 
     $                          MEM( IPW ) )
               XYNORM = DLAPY2( XNORM, YNORM )
               CALL PDMATADD( M, N, SCALE, MEM(IPX), ICF, JCF, DESCC, 
     $                        -ONE, MEM(IPC), ICF, JCF, DESCC )
               CALL PDMATADD( M, N, SCALE, MEM(IPY), ICF, JCF, DESCC, 
     $                        -ONE, MEM(IPF), ICF, JCF, DESCC )
               XNORM = PDLANGE( 'F', M, N, MEM( IPC ), ICF, JCF, DESCC, 
     $                          MEM(IPW) )
               YNORM = PDLANGE( 'F', M, N, MEM( IPF ), ICF, JCF, DESCC, 
     $                          MEM(IPW) )
               FORWARD = DLAPY2( XNORM, YNORM )
               RELFORW = FORWARD / XYNORM
C     
C     From reduction mode we arrive here
C     
 90            CONTINUE
C     
C     Condition estimation
C
               IF( GCSY_CONDEST .AND. (UPLOSIGN.LT.2 .OR.
     $              UPLOSIGN.GT.5) ) THEN
                  ESTTIME = MPI_WTIME()
                  CALL PGCSYCON( GCSY_TRANSAD, GCSY_TRANSBE, 
     $                           GCSY_ISGN, GCSY_COMM, M, N, MEM(IPA), 
     $                           IAD, JAD, DESCA, MEM(IPB), IBE, JBE, 
     $                           DESCB, MEM(IPD), IAD, JAD, DESCD, 
     $                           MEM(IPE), IBE, JBE, DESCE, MBNB2,
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, EST, 
     $                           NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'GCSY: PGCSYCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 100
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'GCSY: PGCSYCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7779 )
     $                 UPLOSIGN,M,N,MB,NB,MB2,NB2,IAD,IBE,NPROW,NPCOL,
     $                 THREADS,MYTIME,TPROF1,TPROF2,TPROF3,GFLOPS,
     $                 SCALE,FORWARD,RELFORW,ABSRESID,RESID,EST,
     $                 ESTTIME,CPROF1,CPROF2,CPROF3,NOITER,INFO
               END IF
C
            END DO
 95         CONTINUE
            END DO
 99         CONTINUE
            END DO
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 100        CONTINUE
C
C     Synchronize all processes. 
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO
         END DO
#endif
         END DO 
         END DO
         END DO
         END DO               
         END DO
         END DO
         END DO
C         
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON GSYL ===
C 
      IF( GSYL ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for GSYL ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6669 )
     $           'ULS','M','N','MB','NB','MB2','NB2','IA','IB',
     $           'P_r','P_c','Thr','Tme','Tp1','Tp2','Tp3','Gfl','Scl',
     $           'E_a','E_r','R_a','R_r','est','etm','Cp1','Cp2','Cp3',
     $           'itr','INF'
         END IF
C     
C     Loop though the space of all combinations of parameters,
C     M, N, MB, NB, NPROW, NPCOL
C
         DO UPLOSIGN = GSYL_UPLOSIGN_MIN, GSYL_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 .OR. UPLOSIGN.EQ.1 ) THEN
            GSYL_TRANSAC = 'N'
            GSYL_TRANSBD = 'N'
         ELSEIF( UPLOSIGN.EQ.2 .OR. UPLOSIGN.EQ.3 ) THEN
            GSYL_TRANSAC = 'N'
            GSYL_TRANSBD = 'T'
         ELSEIF( UPLOSIGN.EQ.4 .OR. UPLOSIGN.EQ.5 ) THEN
            GSYL_TRANSAC = 'T'
            GSYL_TRANSBD = 'N'
         ELSEIF( UPLOSIGN.EQ.6 .OR. UPLOSIGN.EQ.7 ) THEN
            GSYL_TRANSAC = 'T'
            GSYL_TRANSBD = 'T'
         END IF
         GSYL_ISGN = INT(1.0D+00-(MOD(UPLOSIGN,2))*2.0D+00)
         DO M = GSYL_M_MIN, GSYL_M_MAX, GSYL_M_STEP
         DO N = GSYL_N_MIN, GSYL_N_MAX, GSYL_N_STEP
         DO MB = GSYL_MB_MIN, GSYL_MB_MAX, GSYL_MB_STEP
         DO NB = GSYL_NB_MIN, GSYL_NB_MAX, GSYL_NB_STEP
         DO MB2 = GSYL_MB2_MIN, GSYL_MB2_MAX, GSYL_MB2_STEP 
#ifdef LOOPGRID
         DO TEMP_NPROW = GSYL_NPROW_MIN, GSYL_NPROW_MAX, GSYL_NPROW_STEP
         DO TEMP_NPCOL = GSYL_NPCOL_MIN, GSYL_NPCOL_MAX, GSYL_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C           
#ifdef LOOPGRID
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $          MYCOL.LT.0 .OR. MYCOL.GE.NPCOL .OR.
     $          (GSYL_SQUARE_DIM .AND. M.NE.N) .OR.
     $          (GSYL_SQUARE_BLK .AND. MB.NE.NB)  .OR.
     $          (GSYL_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) .OR.
     $           (.NOT.GCSY_PIPELINE .AND. (MB.NE.MB2) ) ) GO TO 120
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            ACROWS = NUMROC( GSYL_M_MAX, MB, MYROW, 0, NPROW )
            ACCOLS = NUMROC( GSYL_M_MAX, MB, MYCOL, 0, NPCOL )
            BDROWS = NUMROC( GSYL_N_MAX, NB, MYROW, 0, NPROW )
            BDCOLS = NUMROC( GSYL_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, GSYL_M_MAX, GSYL_M_MAX, MB, MB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, ACROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GSYL: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 120
            END IF
            CALL DESCINIT( DESCC, GSYL_M_MAX, GSYL_M_MAX, MB, MB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, ACROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GSYL: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 120
            END IF
            CALL DESCINIT( DESCB, GSYL_N_MAX, GSYL_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, BDROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GSYL: DESCINIT(B), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 120
            END IF 
            CALL DESCINIT( DESCD, GSYL_N_MAX, GSYL_N_MAX, NB, NB, 0, 0, 
     $                     TEMP_ICTXT, MAX(1, BDROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GSYL: DESCINIT(D), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 120
            END IF 
            IF( GSYL_SOLVE )
     $           CALL DESCINIT( DESCE, GSYL_M_MAX, GSYL_N_MAX, MB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, ACROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GSYL: DESCINIT(E), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 120
            END IF 
C
C     Loop through the space of parameters IAC, JAC, IBD, JBD, IE, JE,
C     and ignore invalid cases
C
            DO IAC = GSYL_IAC_MIN, GSYL_IAC_MAX, GSYL_IAC_STEP
            JAC = IAC
            IE = IAC
            IF( .NOT. ( IAC + M - 1 ).LE.GSYL_M_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IAC = ',IAC
               GO TO 117
            END IF
            DO IBD = GSYL_IBD_MIN, GSYL_IBD_MAX, GSYL_IBD_STEP
            JBD = IBD
            JE = JBD
            IF( .NOT. ( IBD + N - 1 ).LE.GSYL_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IBD = ',IBD
               GO TO 115
            END IF
            DO REPEAT = 1, GSYL_REPEAT, 1
C
C     Assign pointers for ScaLAPACK arrays
C
                IPA = 1
                IF( LSAME( GSYL_ACSCHUR, 'S' ) ) THEN
                   IPACPY = IPA
                ELSE
                   IPACPY = IPA + DESCA( LLD_) * ACCOLS
                END IF
                IPB = IPACPY + DESCA( LLD_) * ACCOLS
                IF( LSAME( GSYL_BDSCHUR, 'S' ) ) THEN
                   IPBCPY = IPB
                ELSE
                   IPBCPY = IPB + DESCB( LLD_ ) * BDCOLS
                END IF
                IPC = IPBCPY +  DESCB( LLD_ ) * BDCOLS
                IF( LSAME( GSYL_ACSCHUR, 'S' ) ) THEN
                   IPCCPY = IPC
                ELSE
                   IPCCPY = IPC + DESCC( LLD_) * ACCOLS
                END IF
                IPD = IPCCPY + DESCC( LLD_ ) * ACCOLS
                IF( LSAME( GSYL_BDSCHUR, 'S' ) ) THEN
                   IPDCPY = IPD
                ELSE
                   IPDCPY = IPD + DESCD( LLD_ ) * BDCOLS
                END IF
                IF( GSYL_SOLVE ) THEN
                   IPE = IPDCPY + DESCD( LLD_ ) * BDCOLS
                   IPECPY = IPE + DESCE( LLD_ ) * BDCOLS
                   IPX = IPECPY + DESCE( LLD_ ) * BDCOLS
                   DIAG1 = IPX + DESCE( LLD_ ) * BDCOLS 
                   SDIAG = DIAG1 + MAX( GSYL_M_MAX-1, GSYL_N_MAX-1 )
                   DIAG2 = SDIAG + MAX( GSYL_M_MAX, GSYL_N_MAX )
                   IPW = DIAG2 + MAX( GSYL_M_MAX, GSYL_N_MAX )
                ELSE
                   IPE = 1
                   IPECPY = 1
                   IPX = 1
                   DIAG1 = IPDCPY + DESCD( LLD_ ) * BDCOLS 
                   SDIAG = DIAG1 + MAX( GSYL_M_MAX-1, GSYL_N_MAX-1 )
                   DIAG2 = SDIAG + MAX( GSYL_M_MAX, GSYL_N_MAX )
                   IPW = DIAG2 + MAX( GSYL_M_MAX, GSYL_N_MAX )
                END IF
C     
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of (A,C)
C
               DO I = 1, GSYL_M_MAX
                  MEM(DIAG1+(I-1)) = DBLE(I)
                  MEM(DIAG2+(I-1)) = ONE
               END DO
C
C     Init (A,C)
C
               BLKS2B2 = INT(DBLE(GSYL_M_MAX/2) * 
     $                   DBLE(GSYL_AC_NQTRBL)/100.0D+00)
               CALL P2SQMATGD( GSYL_ACDIAG, 'No subdiagonal', 
     $                         GSYL_ACDIAG, GSYL_ACFORM, GSYL_M_MAX, 
     $                         MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DIAG1), MEM(SDIAG), BLKS2B2, 17, 
     $                         MEM(IPC), DESCC, ONE, ONE, 
     $                         MEM(DIAG2), 19, MEM(IPW), WORKSIZ, INFO )
C
C     Copy (A,C)
C
               IF( .NOT. LSAME( GSYL_ACSCHUR, 'S' ) ) THEN
                   CALL PDLACPY( 'All', GSYL_M_MAX, GSYL_M_MAX, 
     $                           MEM(IPA), 1, 1, DESCA, MEM(IPACPY),
     $                           1, 1, DESCA )
                   CALL PDLACPY( 'All', GSYL_M_MAX, GSYL_M_MAX, 
     $                           MEM(IPC), 1, 1, DESCC, MEM(IPCCPY),
     $                           1, 1, DESCC )
                END IF
C
C     Set eigenvalues of (B,D)
C
               DO I = 1, GSYL_M_MAX
                  MEM(DIAG1+(I-1)) = DBLE(GSYL_ISGN)*DBLE(I)
                  MEM(DIAG2+(I-1)) = ONE
               END DO      
C     
C     Initialize (B,D)
C
               BLKS2B2 = INT(DBLE(GSYL_N_MAX/2) * 
     $                   DBLE(GSYL_BD_NQTRBL)/100.0D+00)
               CALL P2SQMATGD( GSYL_BDDIAG, 'No subdiagonal', 
     $                         GSYL_BDDIAG, GSYL_BDFORM, GSYL_N_MAX, 
     $                         MEM(IPB), DESCB, ONE, ONE, 
     $                         MEM(DIAG1), MEM(SDIAG), BLKS2B2, 23, 
     $                         MEM(IPD), DESCD, ONE, ONE, 
     $                         MEM(DIAG2), 31, MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GSYL: P2SQMATGD(B,D), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 120
               END IF 
C
C     Copy (B,D)
C     
               IF( .NOT. LSAME( GSYL_BDSCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', GSYL_N_MAX, GSYL_N_MAX, MEM(IPB),
     $                          1, 1, DESCB, MEM(IPBCPY), 1, 1, DESCB )
                  CALL PDLACPY( 'All', GSYL_N_MAX, GSYL_N_MAX, MEM(IPD),
     $                          1, 1, DESCD, MEM(IPDCPY), 1, 1, DESCD )
               END IF
C     
C     Init a matrix X as the given "solution"
C     
               IF( GSYL_SOLVE ) THEN
                  CALL PDMATGEN2( TEMP_ICTXT, GSYL_EFORM, GSYL_EDIAG, 
     $                            GSYL_M_MAX, GSYL_N_MAX, MB, NB, 
     $                            MEM(IPX), DESCE(LLD_), 0, 0, 7, 0, 
     $                            ACROWS, 0, BDCOLS, MYROW, MYCOL, 
     $                            NPROW, NPCOL )
C     
C     Compute right hand side E as op(A)*X*op(B) +/- op(C)X*op*(D)
C     
                  CALL PDGEMM( GSYL_TRANSAC, 'N' , M, N, M, 1.0D+0, 
     $                         MEM(IPA), IAC, JAC, DESCA, MEM(IPX), IE, 
     $                         JE, DESCE, 0.0D+0, MEM(IPW), IE, JE, 
     $                         DESCE )
                  CALL PDGEMM( 'N', GSYL_TRANSBD, M, N, N, ONE, 
     $                         MEM(IPW), IE, JE, DESCE, MEM(IPB), IBD, 
     $                         JBD, DESCB, 0.0D+0, MEM(IPE), IE, JE, 
     $                         DESCE )
                  CALL PDGEMM( GSYL_TRANSAC, 'N' , M, N, M, 
     $                         DBLE(GSYL_ISGN), MEM(IPC), IAC, JAC, 
     $                         DESCC, MEM(IPX), IE, JE, DESCE, 0.0D+0, 
     $                         MEM(IPW), IE, JE, DESCE )
                  CALL PDGEMM( 'N', GSYL_TRANSBD, M, N, N, ONE, 
     $                         MEM(IPW), IE, JE, DESCE, MEM(IPD), IBD, 
     $                         JBD, DESCD, 1.0D+0, MEM(IPE), IE, JE, 
     $                         DESCE )
C     
C     Copy E for future reference
C     
                  CALL PDLACPY( 'All', M, N, MEM(IPE), IE, JE, DESCE, 
     $                          MEM(IPECPY), IE, JE, DESCE )
               END IF
C     
C     If the matrices are small enough, print them to output
C     
               IF(  IAM.EQ.0 .AND. M.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix A:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, M, MEM( IPA ), IAC, JAC, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix B:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPB ), IBD, JBD, DESCB, 
     $                           0, 0, 'B', NOUT, MEM( IPW ) )
               END IF          
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix C:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, M, MEM( IPC ), IAC, JAC, DESCE, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix D:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPD ), IBD, JBD, DESCD, 
     $                           0, 0, 'D', NOUT, MEM( IPW ) )
               END IF     
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GSYL_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix E:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GSYL_SOLVE ) 
     $              THEN
                  CALL PDLAPRNT( M, N, MEM( IPE ), IE, JE, DESCE, 
     $                           0, 0, 'E', NOUT, MEM( IPW ) )
               END IF
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GSYL_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix X0:'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ .AND. GSYL_SOLVE ) 
     $              THEN
                  CALL PDLAPRNT( M, N, MEM( IPX ), IE, JE, DESCE, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF         
C     
C     Call GSYL solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( GSYL_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGSYLD( 'Reduce', GSYL_ACSCHUR, 
     $                           GSYL_BDSCHUR, GSYL_TRANSAC, 
     $                           GSYL_TRANSBD, GSYL_ISGN, GSYL_COMM, 
     $                           M, N, MEM(IPA), IAC, JAC, DESCA, 
     $                           MEM(IPB), IBD, JBD, DESCB, MEM(IPC), 
     $                           IAC, JAC, DESCC, MEM(IPD), IBD, JBD, 
     $                           DESCD, MEM(IPE), IE, JE, DESCE, 
     $                           MB2, MEM(IPW), WORKSIZ, IMEM, 
     $                           IMEMSIZ, NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( GSYL_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGSYLD( 'Solve', GSYL_ACSCHUR, 
     $                           GSYL_BDSCHUR, GSYL_TRANSAC, 
     $                           GSYL_TRANSBD, GSYL_ISGN, GSYL_COMM, 
     $                           M, N, MEM(IPA), IAC, JAC, DESCA, 
     $                           MEM(IPB), IBD, JBD, DESCB, MEM(IPC), 
     $                           IAC, JAC, DESCC, MEM(IPD), IBD, JBD, 
     $                           DESCD, MEM(IPE), IE, JE, DESCE, 
     $                           MB2, MEM(IPW), WORKSIZ, IMEM, 
     $                           IMEMSIZ, NOEXTSYS, SCALE, INFO ) 
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GSYL: PGEGSYLD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 120
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'GSYL: PGEGSYLD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ 
     $             .AND. GSYL_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( M.LE.PRNTSIZ .AND. N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( M, N, MEM( IPE ), IE, JE, DESCE, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Gflops/sec.
C
               IF( M.GT.N ) THEN
                  GFLOPS = (4.0D+00*DBLE(M)*DBLE(N)**2+
     $                 2.0D+00*DBLE(M)**2*DBLE(N)) /
     $                 (10.0D+00 ** 9 * TPROF3)
               ELSE
                  GFLOPS = (2.0D+00*DBLE(M)*DBLE(N)**2+
     $                 4.0D+00*DBLE(M)**2*DBLE(N)) /
     $                 (10.0D+00 ** 9 * TPROF3)
               END IF
C
C     In reduction mode only, move on to condition estimation
C
               IF( GSYL_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 110 
               END IF
C
C     Compute results of computation
C     
C     Compute the relative residual for GSYL as
C
C         ||(op(A)*X*op(B) +/- op(C)*X*op(D) - SCALE*E)||
C  ------------------------------------------------------------
C  ((||A||*||B||+||C||*||D||)*||X|| + SCALE*||E|| ) * eps 
C            
C               
C     Find machine epsilon
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
C
C     Compute the needed norms
C
               ENORM = PDLANGE( 'I', M, N, MEM( IPECPY ), IE, JE, DESCE, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', M, N, MEM( IPE ), IE, JE, DESCE, 
     $                          MEM( IPW ) )  
               ANORM = PDLANGE( 'I', M, M, MEM( IPACPY ), IAC, JAC, 
     $                          DESCA, MEM( IPW ) )
               BNORM = PDLANGE( 'I', N, N, MEM( IPBCPY ), IBD, JBD, 
     $                          DESCB, MEM( IPW ) )
               CNORM = PDLANGE( 'I', N, N, MEM( IPCCPY ), IAC, JAC, 
     $                          DESCC, MEM( IPW ) )
               DNORM = PDLANGE( 'I', M, M, MEM( IPDCPY ), IBD, JBD, 
     $                          DESCD, MEM( IPW ) )
C
C     Do the matrix multiplications from the relative residual
C              
               CALL PDGEMM( GSYL_TRANSAC, 'N', M, N, M, ONE, 
     $                      MEM( IPACPY), IAC, JAC, DESCA, MEM( IPE ), 
     $                      IE, JE, DESCE, ZERO, MEM( IPW ), IE, JE, 
     $                      DESCE )
               CALL PDGEMM( 'N', GSYL_TRANSBD, M, N, N, ONE, MEM(IPW), 
     $                      IE, JE, DESCE, MEM( IPBCPY ), IBD, JBD,
     $                      DESCB, -SCALE, MEM( IPECPY ), IE, JE, 
     $                      DESCE )
               CALL PDGEMM( GSYL_TRANSAC, 'N', M, N, M, DBLE(GSYL_ISGN), 
     $                      MEM( IPCCPY), IAC, JAC, DESCC, MEM( IPE ), 
     $                      IE, JE, DESCE, ZERO, MEM( IPW ), IE, JE,
     $                      DESCE )
               CALL PDGEMM( 'N', GSYL_TRANSBD, M, N, N, ONE, MEM(IPW), 
     $                      IE, JE, DESCE, MEM( IPDCPY ), IBD, JBD, 
     $                      DESCD, ONE, MEM( IPECPY ), IE, JE, DESCE )
C     
               RNORM = PDLANGE( 'F', M, N, MEM( IPECPY ), IE, JE, DESCE, 
     $                          MEM( IPW ) )
C     
C     Compute residuals 
C     
               ABSRESID = RNORM
               RESID =  PDLANGE( 'I', M, N, MEM( IPECPY ), IE, JE, 
     $              DESCE, MEM( IPW ) ) / 
     $              (((ANORM*BNORM+CNORM*DNORM)*XNORM+SCALE*ENORM)*EPS)
C     
C     Compute forward error norms
C     
C     forward = ||X-XX||_F
C     relforw = ||X-XX||_F / ||X||_F
C     
C     where X is the exact solution and XX is the computed
C     one
C     
               XNORM = PDLANGE( 'F', M, N, MEM( IPX ), IE, JE, DESCE, 
     $                          MEM( IPW ) )
               CALL PDMATADD( M, N, SCALE, MEM(IPX), IE, JE, DESCE, 
     $                        -ONE, MEM(IPE), IE, JE, DESCE )
               FORWARD = PDLANGE( 'F', M, N, MEM( IPE ), IE, JE, DESCE, 
     $                            MEM(IPW) )
               RELFORW = FORWARD / XNORM
C     
C     From reduction mode we arrive here
C     
 110           CONTINUE
C     
C     Condition estimation
C     
               IF( GSYL_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                  CALL PGSYLCON( GSYL_TRANSAC, GSYL_TRANSBD, 
     $                           GSYL_ISGN, GSYL_COMM, M, N, MEM(IPA), 
     $                           IAC, JAC, DESCA, MEM(IPB), IBD, JBD, 
     $                           DESCB, MEM(IPC), IAC, JAC, DESCC, 
     $                           MEM(IPD), IBD, JBD, DESCD, MB2, 
     $                           MEM(IPW), WORKSIZ, IMEM, IMEMSIZ, 
     $                           EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'GSYL: PGSYLCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 120
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'GSYL: PGSYLCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7779 )
     $                 UPLOSIGN,M,N,MB,NB,MB2,NB,IAC,IBD,NPROW,NPCOL,
     $                 THREADS,MYTIME,TPROF1,TPROF2,TPROF3,GFLOPS,
     $                 SCALE,FORWARD,RELFORW,ABSRESID,RESID,EST,ESTTIME,
     $                 CPROF1,CPROF2,CPROF3,NOITER,INFO
               END IF
C
            END DO
 115        CONTINUE
            END DO 
 117        CONTINUE
            END DO
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 120        CONTINUE
C
C     Synchronize all processes.
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO
         END DO 
#endif
         END DO
         END DO
         END DO               
         END DO
         END DO
         END DO
C         
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON GLYCT ===
C
      IF( GLYCT ) THEN
C
C     Print message
C
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for GLYCT ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6670 )
     $           'ULS','N','NB','NB2','IAE','P_r','P_c','Thr','Tme',
     $           'Tp1','Tp2','Tp3','Gfl','Scl','E_a','E_r','R_a','R_r',
     $           'est','etm','Cp1','Cp2','Cp3','itr','INF'
         END IF
C     
C     Loop though the space of all combinations of parameters,
C     N, NB, NPROW, NPCOL
C
         DO UPLOSIGN = GLYCT_UPLOSIGN_MIN, GLYCT_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 ) THEN
            GLYCT_TRANSAE = 'N'
         ELSEIF( UPLOSIGN.EQ.1 ) THEN
            GLYCT_TRANSAE = 'T'
         END IF
         DO N = GLYCT_N_MIN, GLYCT_N_MAX, GLYCT_N_STEP
         DO NB = GLYCT_NB_MIN, GLYCT_NB_MAX, GLYCT_NB_STEP
         DO NB2 = GLYCT_NB2_MIN, GLYCT_NB2_MAX, GLYCT_NB2_STEP
#ifdef LOOPGRID
         DO TEMP_NPROW=GLYCT_NPROW_MIN,GLYCT_NPROW_MAX,GLYCT_NPROW_STEP
         DO TEMP_NPCOL=GLYCT_NPCOL_MIN,GLYCT_NPCOL_MAX,GLYCT_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C           
#ifdef LOOPGRID
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $          MYCOL.LT.0 .OR. MYCOL.GE.NPCOL  .OR.
     $          (GLYCT_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) ) GO TO 140
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            AEROWS = NUMROC( GLYCT_N_MAX, NB, MYROW, 0, NPROW )
            AECOLS = NUMROC( GLYCT_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, GLYCT_N_MAX, GLYCT_N_MAX, NB, NB, 
     $                     0, 0, TEMP_ICTXT, MAX(1, AEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYCT: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 140
            END IF 
            CALL DESCINIT( DESCE, GLYCT_N_MAX, GLYCT_N_MAX, NB, NB, 
     $                     0, 0, TEMP_ICTXT, MAX(1, AEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYCT: DESCINIT(E), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 140
            END IF
            IF( GLYCT_SOLVE )
     $           CALL DESCINIT( DESCC, GLYCT_N_MAX, GLYCT_N_MAX, NB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, AEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYCT: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 140
            END IF 
C
C     Loop through the space of parameter IAE, JAE and
C     ignore invalid cases
C
            DO IAE = GLYCT_IAE_MIN, GLYCT_IAE_MAX, GLYCT_IAE_STEP
            JAE = IAE
            IC = IAE
            JC = JAE
            IF( .NOT. ( IAE + N - 1 ).LE.GLYCT_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IAE = ',IAE
               GO TO 135
            END IF
            DO REPEAT = 1, GLYCT_REPEAT, 1
C
C     Assign pointer for ScaLAPACK arrays
C
               IPA = 1
               IF( .NOT. LSAME( GLYCT_AESCHUR, 'S' ) ) THEN
                  IPACPY = IPA + DESCA( LLD_) * AECOLS
               ELSE
                  IPACPY = IPA
               END IF
               IPE = IPACPY + DESCA( LLD_ ) * AECOLS
               IF( .NOT. LSAME( GLYCT_AESCHUR, 'S' ) ) THEN
                  IPECPY = IPE + DESCE( LLD_ ) * AECOLS
               ELSE
                  IPECPY = IPE
               END IF
               IF( GLYCT_SOLVE ) THEN
                  IPC = IPECPY + DESCE( LLD_ ) * AECOLS
                  IPCCPY = IPC + DESCC( LLD_ ) * AECOLS
                  IPX = IPCCPY + DESCC( LLD_ ) * AECOLS
                  DIAG1 = IPX + DESCC( LLD_ ) * AECOLS 
                  SDIAG = DIAG1 + GLYCT_N_MAX
                  DIAG2 = SDIAG + GLYCT_N_MAX
                  IPW = DIAG2 + GLYCT_N_MAX
               ELSE
                  IPC = 1
                  IPCCPY = 1
                  IPX = 1
                  DIAG1 = IPECPY + DESCE( LLD_ ) * AECOLS 
                  SDIAG = DIAG1 + GLYCT_N_MAX
                  DIAG2 = SDIAG + GLYCT_N_MAX
                  IPW = DIAG2 + GLYCT_N_MAX
               END IF
C     
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
               IF( LSAME( GLYCT_TRANSAE, 'T' ) ) THEN
                  GLYCT_ITRANSAE = 'N'
               ELSE
                  GLYCT_ITRANSAE = 'T'
               END IF
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of (A,E)
C
               DO I = 1, GLYCT_N_MAX
                  MEM(DIAG1+(I-1)) = DBLE(I)
                  MEM(DIAG2+(I-1)) = DBLE(I) + TEN ** 4
               END DO
C
C     Init (A,E)
C
               BLKS2B2 = INT(DBLE(GLYCT_N_MAX/2) * 
     $                   DBLE(GLYCT_AE_NQTRBL)/100.0D+00)
               CALL P2SQMATGD( GLYCT_AEDIAG, 'No subdiagonal', 
     $                         GLYCT_AEDIAG, GLYCT_AEFORM, GLYCT_N_MAX, 
     $                         MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DIAG1), MEM(SDIAG), BLKS2B2, 17,
     $                         MEM(IPE), DESCE, ONE, ONE, 
     $                         MEM(DIAG2), 19, MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GLYCT: P2SQMATGD(A,E), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 140
               END IF 
C
C     Copy A and E
C
               IF( .NOT. LSAME( GLYCT_AESCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', GLYCT_N_MAX, GLYCT_N_MAX, 
     $                          MEM(IPA), IAE, JAE, DESCA, MEM(IPACPY), 
     $                          IAE, JAE, DESCA )
                  CALL PDLACPY( 'All', GLYCT_N_MAX, GLYCT_N_MAX, 
     $                          MEM(IPE), IAE, JAE, DESCE, MEM(IPECPY), 
     $                          IAE, JAE, DESCE )
               END IF
C     
C     Init a matrix X as the given "solution"
C     
               IF( GLYCT_SOLVE ) THEN
                  IF( LSAME( GLYCT_C_SYMM, 'N' ) ) THEN
                     CALL PDMATGEN2( TEMP_ICTXT, 'Symmetric', 
     $                    GLYCT_CDIAG, GLYCT_N_MAX, GLYCT_N_MAX, NB, NB, 
     $                    MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, AEROWS, 0, 
     $                    AECOLS, MYROW, MYCOL, NPROW, NPCOL )
                  ELSE
                     CALL PDLASET( 'All', GLYCT_N_MAX, GLYCT_N_MAX, 
     $                    ZERO, ONE, MEM(IPX), IC, JC, DESCC )
                  END IF
C     
C     Compute the right hand side C
C     
                  CALL PDGEMM( GLYCT_TRANSAE, 'N', N, N, N, ONE, 
     $                         MEM(IPA), IAE, JAE, DESCA, MEM(IPX), 
     $                         IC, JC, DESCC, ZERO, MEM(IPW), 
     $                         IC, JC, DESCC )
                  CALL PDGEMM( 'N', GLYCT_ITRANSAE, N, N, N, ONE, 
     $                         MEM(IPW), IC, JC, DESCC, MEM(IPE), IAE, 
     $                         JAE, DESCE, ZERO, MEM(IPC), IC, JC, 
     $                         DESCC )
                  CALL PDGEMM( GLYCT_TRANSAE, 'N', N, N, N, ONE, 
     $                         MEM(IPE), IAE, JAE, DESCE, MEM(IPX), IC,
     $                         JC, DESCC, ZERO, MEM(IPW), IC, JC, 
     $                         DESCC )
                  CALL PDGEMM( 'N', GLYCT_ITRANSAE, N, N, N, ONE, 
     $                         MEM(IPW), IC, JC, DESCC, MEM(IPA), IAE, 
     $                         JAE, DESCA, ONE, MEM(IPC), IC, JC, 
     $                         DESCC )
C     
C     Copy C for future reference
C     
                  CALL PDLACPY( 'All', N, N, MEM(IPC), IC, JC, DESCC, 
     $                 MEM(IPCCPY), IC, JC, DESCC )
               END IF
C
C     Check if sub(C) is symmetric for this particular instance
C
               IF( LSAME( GLYCT_C_SYMM, 'S' ) .AND. IC.EQ.JC ) THEN
                  GLYCT_C_SYMM2 = 'S'
               ELSE
                  GLYCT_C_SYMM2 = 'N'
               END IF
C 
C     If the matrices are small enough, print them to the screen
C     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(A):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPA ), IAE, JAE, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF
C     
                IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(E):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPE ), IAE, JAE, DESCE, 
     $                           0, 0, 'E', NOUT, MEM( IPW ) )
               END IF
C
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. GLYCT_SOLVE ) 
     $              THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(C):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. GLYCT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. GLYCT_SOLVE ) 
     $              THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X0):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. GLYCT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPX ), IC, JC, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
C
C     Call GLYCT solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( GLYCT_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGLYCTD( 'Reduce', GLYCT_C_SYMM2, 
     $                            GLYCT_TRANSAE, GLYCT_AESCHUR, N, 
     $                            MEM(IPA), IAE, JAE, DESCA, MEM(IPE), 
     $                            IAE, JAE, DESCE, MEM(IPC), IC, JC, 
     $                            DESCC, NB2, MEM(IPW), WORKSIZ, IMEM, 
     $                            IMEMSIZ, NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( GLYCT_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGLYCTD( 'Solve', GLYCT_C_SYMM2, 
     $                            GLYCT_TRANSAE,  GLYCT_AESCHUR, N, 
     $                            MEM(IPA), IAE, JAE, DESCA, MEM(IPE), 
     $                            IAE, JAE, DESCE, MEM(IPC), IC, JC, 
     $                            DESCC, NB2, MEM(IPW), WORKSIZ, IMEM, 
     $                            IMEMSIZ, NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GLYCT: PGEGLYCTD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 140
               END IF
               IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'GLYCT: PGEGLYCTD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. GLYCT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Glops/sec.
C
               IF( LSAME( GLYCT_C_SYMM, 'S' ) ) THEN
                  GFLOPS = ((16.0D+00/3.0D+00)*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               ELSE
                  GFLOPS = ((8.0D+00/3.0D+00)*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               END IF
C
C     In reduction mode only, move on to condition estimation
C
               IF( GLYCT_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 130 
               END IF
C
C     Compute results of computation
C     
C     Compute the relative residual for GLYCT as
C
C            ||(op(A)*X*op(E^T) + op(E)*X*op(A^T) - SCALE*C)||
C  -------------------------------------------------------------------
C             ((2*||A||*||E||)*||X|| + SCALE*||C|| ) * eps 
C     
C     Find machine epsilon
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
C
C     Compute the norms of C and X for usage in the residual 
C     computations
C
               CNORM = PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
C
C     Do the matrix multiplications from the relative residual
C
               CALL PDGEMM( GLYCT_TRANSAE, 'N', N, N, N, ONE, 
     $                      MEM(IPACPY), IAE, JAE, DESCA, MEM(IPC), 
     $                      IC, JC, DESCC, ZERO, MEM(IPW), 
     $                      IC, JC, DESCC )
               CALL PDGEMM( 'N', GLYCT_ITRANSAE, N, N, N, ONE, 
     $                      MEM(IPW), IC, JC, DESCC, MEM(IPECPY), IAE, 
     $                      JAE, DESCE, -SCALE, MEM(IPCCPY), IC, JC, 
     $                      DESCC )
               CALL PDGEMM( GLYCT_TRANSAE, 'N', N, N, N, ONE, 
     $                      MEM(IPECPY), IAE, JAE, DESCE, MEM(IPC), IC,
     $                      JC, DESCC, ZERO, MEM(IPW), IC, JC, 
     $                      DESCC )
               CALL PDGEMM( 'N', GLYCT_ITRANSAE, N, N, N, ONE, 
     $                      MEM(IPW), IC, JC, DESCC, MEM(IPACPY), IAE, 
     $                      JAE, DESCA, ONE, MEM(IPCCPY), IC, JC, 
     $                      DESCC )
C     
C     Compute the rest of the norms
C     
               RNORM = PDLANGE( 'F', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               ANORM = PDLANGE( 'I', N, N, MEM( IPACPY ), IAE, JAE, 
     $                          DESCA, MEM( IPW ) )
               ENORM = PDLANGE( 'I', N, N, MEM( IPECPY ), IAE, JAE, 
     $                          DESCE, MEM( IPW ) )
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = RNORM
               RESID =  PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, 
     $              DESCC, MEM( IPW ) ) / 
     $              (((2*ANORM*ENORM)*XNORM + SCALE*CNORM)*EPS)
C     
C     Compute forward error norms
C     
C     forward = ||X-XX||_F
C     relforw = ||X-XX||_F / ||X||_F
C     
C     where X is the exact solution and XX is the computed
C     one
C     
               XNORM = PDLANGE( 'F', N, N, MEM( IPX ), IC, JC, DESCC, 
     $                         MEM( IPW ) )
               CALL PDMATADD( N, N, SCALE, MEM(IPX), IC, JC, DESCC, 
     $                        -ONE, MEM(IPC), IC, JC, DESCC )
               FORWARD = PDLANGE( 'F', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                            MEM(IPW) )
               RELFORW = FORWARD / XNORM
C     
C     From reduction mode we arrive here
C     
 130           CONTINUE
C     
C     Condition estimation
C     
               IF( GLYCT_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                   CALL PGLYCTCON( GLYCT_TRANSAE, N, MEM(IPA), IAE, 
     $                             JAE, DESCA, MEM(IPE), IAE, JAE, 
     $                             DESCE, NB2, MEM(IPW), WORKSIZ, IMEM, 
     $                             IMEMSIZ, EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'GLYCT: PGLYCTCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 140
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'GLYCT: PGLYCTCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7780 )
     $                 UPLOSIGN,N,NB,NB2,IAE,NPROW,NPCOL,THREADS,MYTIME,
     $                 TPROF1,TPROF2,TPROF3,GFLOPS,SCALE,FORWARD,
     $                 RELFORW,ABSRESID,RESID,EST,ESTTIME,CPROF1,CPROF2,
     $                 CPROF3,NOITER,INFO
               END IF
C
            END DO
 135        CONTINUE
            END DO  
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration
C
 140        CONTINUE
C
C     Synchronize all processes.
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO 
         END DO
#endif
         END DO
         END DO               
         END DO
         END DO
C
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     === GENERATING, SOLVING AND MAKING CONDITIONESTIMATION ON GLYDT ===
C
      IF( GLYDT ) THEN
C     
C     Print message
C     
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT , FMT = *)
            WRITE( NOUT , FMT = *)  ' *** Going for GLYDT ***'
            WRITE( NOUT , FMT = *)
C     
C     Print results table header
C     
            WRITE( NOUT, FMT = 6670 )
     $           'ULS','N','NB','NB2','IAE','P_r','P_c','Thr','Tme',
     $           'Tp1','Tp2','Tp3','Gfl','Scl','E_a','E_r','R_a','R_r',
     $           'est','etm','Cp1','Cp2','Cp3','itr','INF'
         END IF
C     
C     Loop though the space of all combinations of parameters,
C     N, NB, NPROW, NPCOL
C     
         DO UPLOSIGN = GLYDT_UPLOSIGN_MIN, GLYDT_UPLOSIGN_MAX, 1
         IF( UPLOSIGN.EQ.0 ) THEN
            GLYDT_TRANSAE = 'N'
         ELSEIF( UPLOSIGN.EQ.1 ) THEN
            GLYDT_TRANSAE = 'T'
         END IF
         DO N = GLYDT_N_MIN, GLYDT_N_MAX, GLYDT_N_STEP
         DO NB = GLYDT_NB_MIN, GLYDT_NB_MAX, GLYDT_NB_STEP
         DO NB2 = GLYDT_NB2_MIN, GLYDT_NB2_MAX, GLYDT_NB2_STEP
#ifdef LOOPGRID
         DO TEMP_NPROW=GLYDT_NPROW_MIN,GLYDT_NPROW_MAX,GLYDT_NPROW_STEP
         DO TEMP_NPCOL=GLYDT_NPCOL_MIN,GLYDT_NPCOL_MAX,GLYDT_NPCOL_STEP
#endif
C
C     Set INFO to zero for this run
C
            INFO = 0; SCALE = ZERO
C
C     For the current set of parameters, set up temporary grid
C           
#ifdef LOOPGRID
            CALL BLACS_GET( ICTXT, 10, TEMP_ICTXT )
            CALL BLACS_GRIDINIT( TEMP_ICTXT, 'Row-major', TEMP_NPROW, 
     $                           TEMP_NPCOL )
            CALL BLACS_GRIDINFO( TEMP_ICTXT, NPROW, NPCOL, MYROW, 
     $                           MYCOL )
            NPROCS = NPROW*NPCOL
#endif
C
C     Exclude processes not belonging to the temporary grid
C
            IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. 
     $          MYCOL.LT.0 .OR. MYCOL.GE.NPCOL  .OR.
     $          (GLYDT_SQUARE_GRID.AND.
     $           (NPROW.NE.NPCOL.AND.NPROW*NPCOL.NE.2) ) ) GO TO 160
C
C     Count the number of rows and columns of maximum sized problem
C     for the current block sizes and grid properties
C
            AEROWS = NUMROC( GLYDT_N_MAX, NB, MYROW, 0, NPROW )
            AECOLS = NUMROC( GLYDT_N_MAX, NB, MYCOL, 0, NPCOL )
C
C     Set up matrix descriptors for maximum sized problem
C
            CALL DESCINIT( DESCA, GLYDT_N_MAX, GLYDT_N_MAX, NB, NB, 
     $                     0, 0, TEMP_ICTXT, MAX(1, AEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYDT: DESCINIT(A), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 160
            END IF 
            CALL DESCINIT( DESCE, GLYDT_N_MAX, GLYDT_N_MAX, NB, NB, 
     $                     0, 0, TEMP_ICTXT, MAX(1, AEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYDT: DESCINIT(E), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 160
            END IF 
            IF( GLYDT_SOLVE ) 
     $           CALL DESCINIT( DESCC, GLYDT_N_MAX, GLYDT_N_MAX, NB, NB, 
     $                          0, 0, TEMP_ICTXT, MAX(1, AEROWS), INFO )
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
               WRITE( NOUT , FMT = *) 
     $              'GLYDT: DESCINIT(C), INFO=',INFO
               WRITE( NOUT , FMT = *)
               GO TO 160
            END IF 
C
C     Loop through the space of parameter IAE, JAE and
C     ignore invalid cases
C
            DO IAE = GLYDT_IAE_MIN, GLYDT_IAE_MAX, GLYDT_IAE_STEP
            JAE = IAE
            IC = IAE
            JC = JAE
            IF( .NOT. ( IAE + N - 1 ).LE.GLYDT_N_MAX ) THEN
               WRITE( NOUT, FMT = * ) ' *** WARNING ***, PID =',IAM
               WRITE( NOUT, FMT = * ) 'Ignoring IAE = ',IAE
               GO TO 155
            END IF
            DO REPEAT = 1, GLYDT_REPEAT, 1
C
C     Assign pointer for ScaLAPACK arrays
C
               IPA = 1
               IF( .NOT. LSAME( GLYDT_AESCHUR, 'S' ) ) THEN
                  IPACPY = IPA + DESCA( LLD_) * AECOLS
               ELSE
                  IPACPY = IPA
               END IF
               IPE = IPACPY + DESCA( LLD_ ) * AECOLS
               IF( .NOT. LSAME( GLYDT_AESCHUR, 'S' ) ) THEN
                  IPECPY = IPE + DESCE( LLD_ ) * AECOLS
               ELSE
                  IPECPY = IPE
               END IF
               IF( GLYDT_SOLVE ) THEN
                  IPC = IPECPY + DESCE( LLD_ ) * AECOLS
                  IPCCPY = IPC + DESCC( LLD_ ) * AECOLS
                  IPX = IPCCPY + DESCC( LLD_ ) * AECOLS
                  DIAG1 = IPX + DESCC( LLD_ ) * AECOLS 
                  SDIAG = DIAG1 + GLYDT_N_MAX
                  DIAG2 = SDIAG + GLYDT_N_MAX
                  IPW = DIAG2 + GLYDT_N_MAX
               ELSE
                  IPC = 1
                  IPCCPY = 1
                  IPX = 1
                  DIAG1 = IPECPY + DESCE( LLD_ ) * AECOLS 
                  SDIAG = DIAG1 + GLYDT_N_MAX
                  DIAG2 = SDIAG + GLYDT_N_MAX
                  IPW = DIAG2 + GLYDT_N_MAX
               END IF
C     
C     Set some variables
C
               WORKSIZ = MEMSIZ - (IPW+1)
               NOEXTSYS = 0
               IF( LSAME( GLYDT_TRANSAE, 'T' ) ) THEN
                  GLYDT_ITRANSAE = 'N'
               ELSE
                  GLYDT_ITRANSAE = 'T'
               END IF
C
C     Init memory with zeros
C
               DO I = 1, MEMSIZ
                  MEM(I) = ZERO
               END DO
               DO I = 1, IMEMSIZ
                  IMEM( I ) = IZERO
               END DO
C
C     Set eigenvalues of (A,E)
C
               DO I = 1, GLYDT_N_MAX
                  MEM(DIAG1+(I-1)) = DBLE(I) + DBLE(N)
                  MEM(DIAG2+(I-1)) = ONE
               END DO
C
C     Init (A,E)
C
               BLKS2B2 = INT(DBLE(GLYDT_N_MAX/2) * 
     $                   DBLE(GLYDT_AE_NQTRBL)/100.0D+00)
               CALL P2SQMATGD( GLYDT_AEDIAG, 'No subdiagonal', 
     $                         GLYDT_AEDIAG, GLYDT_AEFORM, GLYDT_N_MAX, 
     $                         MEM(IPA), DESCA, ONE, ONE, 
     $                         MEM(DIAG1), MEM(SDIAG), BLKS2B2, 17,
     $                         MEM(IPE), DESCE, ONE, ONE, 
     $                         MEM(DIAG2), 19, MEM(IPW), WORKSIZ, INFO )
               IF( INFO.NE.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GLYDT: P2SQMATGD(A,E), INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 160
               END IF 
C
C     Copy A and E
C
               IF( .NOT. LSAME( GLYDT_AESCHUR, 'S' ) ) THEN
                  CALL PDLACPY( 'All', GLYDT_N_MAX, GLYDT_N_MAX, 
     $                          MEM(IPA), IAE, JAE, DESCA, MEM(IPACPY), 
     $                          IAE, JAE, DESCA )
                  CALL PDLACPY( 'All', GLYDT_N_MAX, GLYDT_N_MAX, 
     $                          MEM(IPE), IAE, JAE, DESCE, MEM(IPECPY), 
     $                          IAE, JAE, DESCE )
               END IF
C     
C     Init a matrix X as the given "solution"
C     
               IF( GLYDT_SOLVE ) THEN
                  IF( LSAME( GLYDT_C_SYMM, 'N' ) ) THEN
                     CALL PDMATGEN2( TEMP_ICTXT, 'Symmetric', 
     $                    GLYDT_CDIAG, GLYDT_N_MAX, GLYDT_N_MAX, NB, NB, 
     $                    MEM(IPX), DESCC(LLD_), 0, 0, 7, 0, AEROWS, 0, 
     $                    AECOLS, MYROW, MYCOL, NPROW, NPCOL )
                  ELSE
                     CALL PDLASET( 'All', GLYDT_N_MAX, GLYDT_N_MAX, 
     $                    ZERO, ONE, MEM(IPX), IC, JC, DESCC )
                  END IF
C     
C     Compute the right hand side C
C     
                  CALL PDGEMM( GLYDT_TRANSAE, 'N', N, N, N, ONE, 
     $                         MEM(IPA), IAE, JAE, DESCA, MEM(IPX), 
     $                         IC, JC, DESCC, ZERO, MEM(IPW), 
     $                         IC, JC, DESCC )
                  CALL PDGEMM( 'N', GLYDT_ITRANSAE, N, N, N, ONE, 
     $                         MEM(IPW), IC, JC, DESCC, MEM(IPA), IAE, 
     $                         JAE, DESCA, ZERO, MEM(IPC), IC, JC, 
     $                         DESCC )
                  CALL PDGEMM( GLYDT_TRANSAE, 'N', N, N, N, -ONE, 
     $                        MEM(IPE), IAE, JAE, DESCE, MEM(IPX), IC, 
     $                        JC, DESCC, ZERO, MEM(IPW), IC, JC, DESCC )
                  CALL PDGEMM( 'N', GLYDT_ITRANSAE, N, N, N, ONE, 
     $                         MEM(IPW), IC, JC, DESCC, MEM(IPE), IAE, 
     $                         JAE, DESCE, ONE, MEM(IPC), IC, JC, 
     $                         DESCC )
C     
C     Copy C for future reference
C     
                  CALL PDLACPY( 'All', N, N, MEM(IPC), IC, JC, DESCC, 
     $                          MEM(IPCCPY), IC, JC, DESCC )
               END IF
C
C     Check if sub(C) is symmetric for this particular instance
C
               IF( LSAME( GLYDT_C_SYMM, 'S' ) .AND. IC.EQ.JC ) THEN
                  GLYDT_C_SYMM2 = 'S'
               ELSE
                  GLYDT_C_SYMM2 = 'N'
               END IF
C 
C     If the matrices are small enough, print them to the screen
C     
               IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(A):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPA ), IAE, JAE, DESCA, 
     $                           0, 0, 'A', NOUT, MEM( IPW ) )
               END IF
C     
                IF(  IAM.EQ.0 .AND. N.LE.PRNTSIZ ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(E):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPE ), IAE, JAE, DESCE, 
     $                           0, 0, 'E', NOUT, MEM( IPW ) )
               END IF
C
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. GLYDT_SOLVE ) 
     $              THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(C):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. GLYDT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'C', NOUT, MEM( IPW ) )
               END IF
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. GLYDT_SOLVE ) 
     $              THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X0):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. GLYDT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPX ), IC, JC, DESCC, 
     $                           0, 0, 'X0', NOUT, MEM( IPW ) )
               END IF
C
C     Call GLYDT solver routine
C
               CALL BLACS_BARRIER( TEMP_ICTXT, 'All' )
               IF( GLYDT_REDUCE_ONLY ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGLYDTD( 'Reduce', GLYDT_C_SYMM2, 
     $                            GLYDT_TRANSAE, GLYDT_AESCHUR, N, 
     $                            MEM(IPA), IAE, JAE, DESCA, MEM(IPE), 
     $                            IAE, JAE, DESCE, MEM(IPC), IC, JC, 
     $                            DESCC, NB2, MEM(IPW), WORKSIZ, IMEM, 
     $                            IMEMSIZ, NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME()
               ELSEIF( GLYDT_SOLVE ) THEN
                  T1 = MPI_WTIME()
                  CALL PGEGLYDTD( 'Solve', GLYDT_C_SYMM2, 
     $                            GLYDT_TRANSAE, GLYDT_AESCHUR, N, 
     $                            MEM(IPA), IAE, JAE, DESCA, MEM(IPE), 
     $                            IAE, JAE, DESCE, MEM(IPC), IC, JC, 
     $                            DESCC, NB2, MEM(IPW), WORKSIZ, IMEM, 
     $                            IMEMSIZ, NOEXTSYS, SCALE, INFO )
                  T2 = MPI_WTIME() 
               END IF
               IF( INFO.LT.0 ) THEN
                  WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                  WRITE( NOUT , FMT = *) 
     $                 'GLYDT: PGEGLYDTD, INFO=',INFO
                  WRITE( NOUT , FMT = *)
                  GO TO 160
               END IF
                IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN 
                  WRITE( NOUT , FMT = *)  
     $                 'GLYDT: PGEGLYDTD, INFO=',INFO
               END IF
               MYTIME = T2 - T1
               TPROF1 = MEM(IPW)
               TPROF2 = MEM(IPW+1)
               TPROF3 = MEM(IPW+2)
C
C     Print out the solution matrix if possible
C
               IF( IAM.EQ.0 .AND. N.LE.PRNTSIZ .AND. GLYDT_SOLVE ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = * ) 'Matrix sub(X):'
                  WRITE( NOUT, FMT = * )
               END IF
               IF( N.LE.PRNTSIZ .AND. GLYDT_SOLVE ) THEN
                  CALL PDLAPRNT( N, N, MEM( IPC ), IC, JC, DESCC, 
     $                           0, 0, 'X', NOUT, MEM( IPW ) )
               END IF
C
C     Global maximum reduction on timings
C
               IF( NPROCS .GT. 1 ) THEN
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, MYTIME, 1,
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF1, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF2, 1, 
     $                          -1, -1, -1, -1, -1 )
                  CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, TPROF3, 1, 
     $                          -1, -1, -1, -1, -1 )
               END IF
C
C     Compute Glops/sec.
C
               IF( LSAME( GLYDT_C_SYMM, 'S' ) ) THEN
                  GFLOPS = ((16.0D+00/3.0D+00)*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               ELSE
                  GFLOPS = ((8.0D+00/3.0D+00)*DBLE(N)**3) / 
     $                 (10.0D+00 ** 9 * TPROF3)
               END IF
C
C     In reduction mode only, move on to condition estimation
C
               IF( GLYDT_REDUCE_ONLY ) THEN
                  ABSERROR = ZERO
                  RELERROR = ZERO
                  ABSRESID = ZERO
                  RESID = ZERO
                  GO TO 150 
               END IF
C
C     Compute results of computation
C     
C     Compute the relative residual for GLYDT as
C
C            ||(op(A)*X*op(A^T) + op(E)*X*op(E^T) - SCALE*C)||
C  ---------------------------------------------------------------------
C             ((||A||**2+||E||**2)*||X|| + SCALE*||C|| ) * eps 
C     
C     Find machine epsilon
C     
               EPS = PDLAMCH( TEMP_ICTXT, 'Epsilon' )
C
C     Compute the norms of C and X for usage in the residual 
C     computations
C
               CNORM = PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM(IPW) )
               XNORM = PDLANGE( 'I', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
C
C     Do the matrix multiplications from the relative residual
C
               CALL PDGEMM( GLYDT_TRANSAE, 'N', N, N, N, ONE, 
     $                      MEM(IPACPY), IAE, JAE, DESCA, MEM(IPC), 
     $                      IC, JC, DESCC, ZERO, MEM(IPW), 
     $                      IC, JC, DESCC )
               CALL PDGEMM( 'N', GLYDT_ITRANSAE, N, N, N, ONE, 
     $                      MEM(IPW), IC, JC, DESCC, MEM(IPACPY), IAE, 
     $                      JAE, DESCA, -SCALE, MEM(IPCCPY), IC, JC, 
     $                      DESCC )
               CALL PDGEMM( GLYDT_TRANSAE, 'N', N, N, N, -ONE, 
     $                      MEM(IPECPY), IAE, JAE, DESCE, MEM(IPC), IC, 
     $                      JC, DESCC, ZERO, MEM(IPW), IC, JC, DESCC )
               CALL PDGEMM( 'N', GLYDT_ITRANSAE, N, N, N, ONE, 
     $                      MEM(IPW), IC, JC, DESCC, MEM(IPECPY), IAE, 
     $                      JAE, DESCE, ONE, MEM(IPCCPY), IC, JC, 
     $                      DESCC )
C     
C     Compute the rest of the norms
C     
               RNORM = PDLANGE( 'F', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               ANORM = PDLANGE( 'I', N, N, MEM( IPACPY ), IAE, JAE, 
     $                          DESCA, MEM( IPW ) )
               ENORM = PDLANGE( 'I', N, N, MEM( IPECPY ), IAE, JAE, 
     $                          DESCE, MEM( IPW ) )
C     
C     Compute residuals and draw conclusions
C     
               ABSRESID = RNORM
               RESID = PDLANGE( 'I', N, N, MEM( IPCCPY ), IC, JC, DESCC, 
     $              MEM( IPW ) ) / 
     $              (((ANORM**2+ENORM**2)*XNORM + SCALE*CNORM)*EPS)
C     
C     Compute forward error norms
C     
C     forward = ||X-XX||_F
C     relforw = ||X-XX||_F / ||X||_F
C     
C     where X is the exact solution and XX is the computed
C     one
C
               XNORM = PDLANGE( 'F', N, N, MEM( IPX ), IC, JC, DESCC, 
     $                          MEM( IPW ) )
               CALL PDMATADD( N, N, SCALE, MEM(IPX), IC, JC, DESCC, 
     $                       -ONE, MEM(IPC), IC, JC, DESCC )
               FORWARD = PDLANGE( 'F', N, N, MEM( IPC ), IC, JC, DESCC, 
     $                            MEM(IPW) )
               RELFORW = FORWARD / XNORM
C
C     From reduction mode we arrive here
C
 150           CONTINUE
C
C     Condition estimation
C
               IF( GLYDT_CONDEST ) THEN
                  ESTTIME = MPI_WTIME()
                   CALL PGLYDTCON( GLYDT_TRANSAE, N, MEM(IPA), IAE, 
     $                             JAE, DESCA, MEM(IPE), IAE, JAE, 
     $                             DESCE, NB2, MEM(IPW), WORKSIZ, IMEM, 
     $                             IMEMSIZ, EST, NOITER, INFO )
                  ESTTIME = MPI_WTIME() - ESTTIME
                  IF( INFO.LT.0 ) THEN
                     WRITE( NOUT , FMT = *) ' *** ERROR ***, PID = ',IAM
                     WRITE( NOUT , FMT = *) 
     $                    'GLYDT: PGLYDTCON, INFO=',INFO
                     WRITE( NOUT , FMT = *)
                     GO TO 160
                  END IF 
                  IF( IAM.EQ.0 .AND. INFO.GT.0 .AND. .FALSE. ) THEN
                     WRITE( NOUT , FMT = *) 
     $                    'GLYDT: PGLYDTCON, INFO=',INFO
                  END IF 
                  CPROF1 = MEM(IPW)
                  CPROF2 = MEM(IPW+1)
                  CPROF3 = MEM(IPW+2)
                  IF( NPROCS .GT. 1 ) THEN
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             ESTTIME, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF1, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF2, 1, -1, -1, -1, -1, -1 )
                     CALL DGAMX2D( TEMP_ICTXT, 'All', ' ', 1, 1, 
     $                             CPROF3, 1, -1, -1, -1, -1, -1 )
                  END IF
               ELSE
                  NOITER = IZERO
                  EST = -ONE
                  ESTTIME = ZERO
                  CPROF1 = ZERO
                  CPROF2 = ZERO
                  CPROF3 = ZERO
               END IF
C
C     Write results to output
C
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 7780 )
     $                 UPLOSIGN,N,NB,NB2,IAE,NPROW,NPCOL,THREADS,MYTIME,
     $                 TPROF1,TPROF2,TPROF3,GFLOPS,SCALE,FORWARD,
     $                 RELFORW,ABSRESID,RESID,EST,ESTTIME,CPROF1,CPROF2,
     $                 CPROF3,NOITER,INFO
               END IF
C
            END DO
 155        CONTINUE
            END DO  
C
C     Release temporary grid
C
#ifdef LOOPGRID
            CALL BLACS_GRIDEXIT( TEMP_ICTXT )
#endif
C
C     Leftover processes go here and wait for next iteration.
C
 160        CONTINUE            
C
C     Synchronize all processes. 
C
            CALL BLACS_BARRIER( ICTXT, 'All' )
C 
#ifdef LOOPGRID
         END DO 
         END DO
#endif
         END DO
         END DO
         END DO
         END DO
C
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Clean exit
C
      GO TO 9999
C
C     If an error occured, leave catastrophically
C
 9998 CONTINUE
#ifdef USE_DYNAMIC
      IF( ALLOCATED(MEM) ) DEALLOCATE( MEM )
      IF( ALLOCATED(IMEM) ) DEALLOCATE( IMEM )
#endif
      TOTTIME = MPI_WTIME() - TOTTIME
      CALL DGAMX2D(ICTXT,'All',' ',1,1, TOTTIME,1,-1,-1,-1,-1,-1)
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * ) 
         WRITE( NOUT, FMT = * ) 
     $        'TESTSCASY executed in seconds:',TOTTIME
         WRITE( NOUT, FMT = * ) 
      END IF
      CALL BLACS_ABORT( ICTXT, ERRORNUM )
C
C     If everything went on normally, leave nicely
C
 9999 CONTINUE
C
#ifdef USE_DYNAMIC
      IF( ALLOCATED(MEM) ) DEALLOCATE( MEM )
      IF( ALLOCATED(IMEM) ) DEALLOCATE( IMEM )
#endif
      TOTTIME = MPI_WTIME() - TOTTIME
      CALL DGAMX2D(ICTXT,'All',' ',1,1, TOTTIME,1,-1,-1,-1,-1,-1)
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * ) 
         WRITE( NOUT, FMT = * ) 
     $        'TESTSCASY executed in seconds:',TOTTIME
          WRITE( NOUT, FMT = * ) 
      END IF
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_EXIT( 0 )
C
C     Formats
C
 6666 FORMAT (A4,A5,A5,A4,A4,A5,A5,A5,A5,A5,A7,A7,A7,A7,A7,A7,A9,A9,A9,
     $        A9,A7,A7,A7,A7,A4,A4)
 6667 FORMAT (A4,A5,A4,A5,A5,A5,A5,A7,A7,A7,A7,A7,A7,A9,A9,A9,A9,A7,A7,
     $        A7,A7,A4,A4)
 6669 FORMAT (A4,A5,A5,A4,A4,A4,A4,A5,A5,A5,A5,A5,A7,A7,A7,A7,A7,A9,A9,
     $        A9,A9,A9,A9,A7,A7,A7,A7,A4,A4)
 6670 FORMAT (A4,A5,A4,A5,A5,A5,A5,A5,A7,A7,A7,A7,A7,A9,A9,A9,A9,A9,A9,
     $        A7,A7,A7,A7,A4,A4)
 7777 FORMAT (I4,I5,I5,I4,I4,I5,I5,I5,I5,I5,F7.1,F7.1,F7.1,F7.1,F7.1,
     $        F7.3,E9.1,E9.1,E9.1,E9.1,F7.1,F7.1,F7.1,F7.1,I4,I4) 
 7778 FORMAT (I4,I5,I4,I5,I5,I5,I5,F7.1,F7.1,F7.1,F7.1,F7.1,F7.3,E9.1,
     $        E9.1,E9.1,E9.1,F7.1,F7.1,F7.1,F7.1,I4,I4)  
 7779 FORMAT (I4,I5,I5,I4,I4,I4,I4,I5,I5,I5,I5,I5,F7.1,F7.1,F7.1,F7.1,
     $        F7.3,E9.1,E9.1,E9.1,E9.1,E9.1,E9.1,F7.1,F7.1,F7.1,F7.1,I4,
     $        I4) 
 7780 FORMAT (I4,I5,I4,I5,I5,I5,I5,I5,F7.1,F7.1,F7.1,F7.1,F7.3,E9.1,
     $        E9.1,E9.1,E9.1,E9.1,E9.1,F7.1,F7.1,F7.1,F7.1,I4,I4) 
C
C     End of program TESTSCASY
C
      END
C
C *** Last line of TESTSCASY ***
