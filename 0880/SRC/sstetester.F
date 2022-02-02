PROGRAM SSTETESTER
!
USE GSTEDEFINITIONS
USE SSTEDEFINITIONS
USE SSTEINTERFACES3
! 
!==============================================================================!
!                                                                              !
! Version 08/08/2007, by Osni Marques (LBNL)                                   !
!                                                                              !
! Purpose:                                                                     !
! =======                                                                      !
!                                                                              !
! SSTETESTER is a program to test the symmetric tridiagonal eigenvalue solvers !
! implemented in LAPACK.                                                       !
!                                                                              !
! The test matrices used by SSTETESTER can be generated automatically (as      !
! described below) or generated from data read from files. In this case, the   !
! the data input a symmetric tridiagonal matrix or an eigenvalue distribution, !
! from which a symmetric tridiagonal matrix will then be generated.            !
!                                                                              !
! Input data for SSTETESTER:                                                   !
! =========================                                                    !
!                                                                              !
! Input data for SSTETESTER are specified by means of macro commands (or key   !
! key words). These macros are groups of characters that uniquely define a     !
! specific subset of the input data, such as matrix types and dimensions,      !
! matrices to be read from files, etc. The input data are usually read from    !
! a file.                                                                      !
!                                                                              !
! Macros can be specified in any order, either in lower or upper case, with    !
! the corresponding subset of data they define. Data can be separated by       !
! blanks or commas, and the character % is interpreted as the beginning        !
! of a comment. Therefore, a line beginning with % is ignored.                 !
!                                                                              !
! In the directory INPUT the user should find three examples of input files:   !
!                                                                              !
!    input.test.easy                                                           !
!    input.test.mild                                                           !
!    input.test.hard                                                           !
!                                                                              !
! which exercise various features of SSTETESTER. After a successful            !
! compilation, the user should do, for example,                                !
!                                                                              !
!    dstetester.x < input.easy                                                 !
!                                                                              !
! The macros currently supported by SSTETESTER are described below. In the     !
! sequence, macros and other character entries are written between quotes      !
! to distiguish them from numerical data.                                      !
!                                                                              !
! "CALLST"  "CALLST_LIST"                                                      !
!     Sets the subroutines to be tested, where "CALLST_LIST" can be            !
!         "STEQRV"  (calls SSTEQR with COMPZ='V'),                             !
!         "STEVXA"  (calls SSTEVX with RANGE='A'),                             !
!         "STEVXI"  (calls SSTEVX with RANGE='I'),                             !
!         "STEVXV"  (calls SSTEVX with RANGE='V'),                             !
!         "STEDCI"  (calls SSTEDC with COMPZ='I'),                             !
!         "STEGRA"  (calls SSTEGR with RANGE='A'),                             !
!         "STEGRI"  (calls SSTEGR with RANGE='I'),                             !
!         "STEGRV"  (calls SSTEGR with RANGE='V'),                             !
!         "ALL"     (all the tests above).                                     !
! "DUMP"  "LIST"                                                               !
!     Defines data to be written into files, where "LIST" can be one of        !
!     the following                                                            !
!         "T"    file "_stetester.out.T", tridiagonal matrix (i,d_i,e_i)       !
!         "W"    file "_stetester.out.W", eigenvalues                          !
!         "Z"    file "_stetester.out.Z", eigenvectors                         !
!         "LOG"  file "_stetester.out.log", timing, residuals, orthogonality   !
!         "T.M"  file "_stetester.out.m", T in Matlab format                   !
!         "W.M"  file "_stetester.out.m", W in Matlab format                   !
!         "Z.M"  file "_stetester.out.m", Z in Matlab format                   !
! "ECOND"  K                                                                   !
!     Rets COND, which is used in the generation of built-in eigenvalue        !
!     distributions 1 to 5. K (integer) can be                                 !
!         = 1, COND = 1 / SQRT(ULP), default,                                  !
!         = 2, COND = 1 / (N*SQRT(ULP)),                                       !
!         = 3, COND = 1 / (10*N*SQRT(ULP)),                                    !
!         = 4, COND = 1 / ULP,                                                 !
!         = 5, COND = 1 / (N*ULP),                                             !
!         = 6, COND = 1 / (10*N*ULP),                                          !
!     where                                                                    !
!         ULP = (relative machine precision)*(base of the machine)             !
!         N   = dimension of the matrix to be generated                        !
! "EDIST"  K                                                                   !
!     Resets the type of the distribution to be used to generate random        !
!     eigenvalues. K (integer) can be                                          !
!         = 1, uniform distribution (-1,1), default,                           !
!         = 2, uniform distribution ( 0,1),                                    !
!         = 3, normal distribution  ( 0,1)                                     !
! "ESIGN"  K                                                                   !
!     Assigns signs to the eigenvalues. K (integer) can be                     !
!         = 0, the eigenvalues will not be negative, default,                  !
!         = 1, the eigenvalues may be positive, negative, or zero              !
! "EIGVAL"                                                                     !
!     Defines built-in eigenvalue distributions to be used in the generation   !
!     of test matrices. The next two lines must define                         !
!          ETYPE( 1 )  ETYPE( 2 )  ETYPE( 3 ) ...                              !
!          ESIZE( 1 )  ESIZE( 2 )  ESIZE( 3 ) ...                              !
!     where ETYPE (integer) is a list of built-in eigenvalue distributions     !
!     (the implemented distributions are described below), and ESIZE (integer) !
!     is a list of dimensions of the matrices to be generated. The dimensions  !
!     can also be defined by NMIN[:NINC]:NMAX, NMIN>0, NINC>0, NMAX>=NMIN.     !
! "EIGVALF"  "EIGVALF_FILE"                                                    !
!     Defines a file containing an eigenvalue distribution to be used in the   !
!     generation of a tridiagonal matrix, where "EIGVF_FILE" is a file name.   !
!     This file should contain only one entry per line as follows:             !
!          NEIG                                                                !
!          EIG( 1 )                                                            !
!          EIG( 2 )                                                            !
!          ...                                                                 !
!     which will then be used to generate a symmetric tridiagonal matrix with  !
!     eigenvalues ( EIG( 1 ) EIG( 2 ) ... EIG( NEIG ) ).                       !
! "EIGVI"                                                                      !
!     Defines indices of the smallest and largest eigenvalues to be computed.  !
!     The next two lines must define                                           !
!          IL( 1 )  IL( 2 )  ...  IL( NI )                                     !
!          IU( 1 )  IU( 2 )  ...  IU( NI )                                     !
!     with 1 <= IL( I ) <= IU( I ), I=1,2,...NI. These indices are used only   !
!     in the tests where RANGE='I'.                                            !
! "EIGVV"                                                                      !
!     Defines lower and upper bounds of intervals to be searched for           !
!     eigenvalues. The next two lines must define                              !
!          VL( 1 )  VL( 2 )  ...  VL( NV )                                     !
!          VU( 1 )  VU( 2 )  ...  VU( NV )                                     !
!     with VL( I ) < VU( I ), I=1,2,...NV. These intervals are used only in    !
!     the tests where RANGE='V'.                                               !
! "GLUED"                                                                      !
!     Defines glued tridiagonal matrices, for example                          !
!              | T_1  g_1      |                                               !
!          T = | g_1  T_2  g_2 |                                               !
!              |      g_2  T_3 |                                               !
!     The next four lines must define                                          !
!          GFORM( 1 )  GFORM( 2 ) ... GFORM( K-1 )  GFORM( K )                 !
!          GTYPE( 1 )  GTYPE( 2 ) ... GTYPE( K-1 )  GTYPE( K )                 !
!          GSIZE( 1 )  GSIZE( 2 ) ... GSIZE( K-1 )  GSIZE( K )                 !
!          GAMMA( 1 )  GAMMA( 2 ) ... GAMMA( K-1 )                             !
!     where GFORM (integer) defines how the matrix is generated                !
!          = 1, built-in eigenvalue distribution,                              !
!          = 2, built-in tridiagonal matrix,                                   !
!     GTYPE(integer) defines the built-in eigenvalue distribution or matrix    !
!     (accordingly to GFORM), GSIZE (integer) defines the dimension of the     !
!     matrix, and GAMMA (real) defines the glue factor.                        !
! "HBANDA"  K                                                                  !
!     Sets the halfbandwidth of the symmetric matrix to be generated and       !
!     then tridiagonalized. K (integer) should be a number between 1 and 100,  !
!     then a matrix with max(1,N*(K/100)) subdiagonals is generated. By        !
!     default, K=100.                                                          !
! "HBANDR"  K                                                                  !
!     Sets the halfbandwidth of the matrices Z'*T*Z and Z'*Z used in the       !
!     tests (1) and (2) defined below. K (integer) should be a number between  !
!     1 and 100, then max(1,N*(K/100)) subdiagonals of those matrices are      !
!     computed. If K=0 the tests are not performed and are simply set to       !
!     0. By default, K=100.                                                    !
! "ISEED"  ISEED( 1 )  ISEED( 2 )  ISEED( 3 )  ISEED( 4 )                      !
!     Resets the seed of the random number generator. Each entry of            !
!     ISEED (integer) should lie between 0 and 4095 inclusive and ISEED(4)     !
!     should be odd. By default the initial seed is 4, 3, 2, 1.                !
! "MATRIX"                                                                     !
!     Defines built-in tridiagonal matrices to be used in the tests. The next  !
!     two lines must define                                                    !
!          MTYPE( 1 )  MTYPE( 2 )  MTYPE( 3 ) ...                              !
!          MSIZE( 1 )  MSIZE( 2 )  MSIZE( 3 ) ...                              !
!     where MTYPE(integer) is a list of built-in tridiagonal matrices (the     !
!     implemented types are described below), and MSIZE (integer) is a list    !
!     of dimensions of the matrices to be generated. The dimensions can        !
!     also be defined by NMIN[:NINC]:NMAX, NMIN>0, NINC>0, NMAX>=NMIN.         !
! "MATRIXF"  "MATRIXF_FILE"                                                    !
!     Defines a file containing a tridiagonal matrix, where "MATRIXF_FILE"     !
!     is a file name. This file should contain:                                !
!          NDE                                                                 !
!          row_index_1  D( row_index_1 )  E( row_index_1 )                     !
!          row_index_2  D( row_index_2 )  E( row_index_2 )                     !
!          ...                                                                 !
!     which will then be used to generate a symmetric tridiagonal matrix of    !
!     the form                                                                 !
!                      E( 1 )  E( 2 )  ... E( NDE-1 )                          !
!          tridiag( D( 1 )  D( 2 )  D( 3 ) ... D( NDE ) ).                     !
!                      E( 1 )  E( 2 )  ... E( NDE-1 )                          !
! "MINTT"  MINTT                                                               !
!     Sets the minimum accumulated timing for each test specified by CALLST    !
!     such that time_per_call * number_of_calls > MINTT (this option allows    !
!     for more accurate timings for relatively small matrices in particular    !
!     on fast processors). The default is 0.05s.                               !
! "NRILIU"  K                                                                  !
!     Defines the number of K random indices of the smallest and largest       !
!     eigenvalues to be computed. These indices are used only in the tests     !
!     where RANGE='I'.                                                         !
! "NRVLVU"  K                                                                  !
!     Defines the number of K random lower and upper bounds of intervals       !
!     to be searched for eigenvalues. These intervals are used only in the     !
!     tests where RANGE='V'.                                                   !
! "TOLEBZ"  TOLEBZ                                                             !
!     Sets TOLEBZ (real) the absolute error tolerance for the eigenvalues      !
!     to  be computed by SSTEBZ. The default is twice the underflow threshold  !
!     (see the leading comments of SSTEBZ for details).                        !
! "TOLEGR"  TOLEGR ** This option has been deprecated **                       !
!     Sets TOLEGR (real) the absolute error tolerance for the eigenvalues to   !
!     be computed by SSTEGR. The default is twice the underflow threshold      !
!     (see the leading comments of SSTEGR for details).                        !
! "END"                                                                        !
!     End of data (subsequent lines are ignored).                              !
!                                                                              !
! Built-in eigenvalue distributions:                                           !
! ---------------------------------                                            !
!                                                                              !
! Type  Description                                                            !
!                                                                              !
! (1)   EIG( 1 ) = 1 and EIG( 2:N ) = 1.0/COND;                                !
! (2)   EIG( 1:N-1 ) = 1 and EIG( N ) = 1.0/COND;                              !
! (3)   EIG( I ) = COND**(-(I-1)/(N-1));                                       !
! (4)   EIG( I ) = 1-(I-1)/(N-1)*(1-1/COND);                                   !
! (5)   Eigenvalue distribution set to random numbers in the range (1/COND,1), !
!       their logarithms are uniformly distributed;                            !
! (6)   Eigenvalue distribution set to random numbers from the same            !
!       distribution as of the rest of the matrix;                             !
! (7)   EIG( I ) = ULP*I, I = 1,2,...N-1, and EIG( N ) = 1.0.                  !
! (8)   EIG( 1 ) = ULP, EIG( I ) = 1+SQRT(ULP)*I, I=2,3,...N-1, and EIG( N )=2 !
! (9)   EIG( 1 ) = 1, EIG( I ) = EIG( I-1 )+100*ULP, I=1,2,...N-1.             !
!                                                                              !
! TYPE < 0 has the same meaning as ABS(TYPE), except that the elements of EIG  !
! are reversed. For example, TYPE = -1 implies EIG( 1:N-1 ) = 1.0/COND and     !
! EIG( N ) = 1. If ESIGN = 1, then the elements of EIG are multiplied by a     !
! random sign (i.e., +1 or -1).                                                !
!                                                                              !
! Built-in matrices:                                                           !
! -----------------                                                            !
!                                                                              !
! Type  Description                                                            !
!                                                                              !
! (0)   The zero matrix.                                                       !
! (1)   The identity matrix.                                                   !
! (2)   (1,2,1) tridiagonal matrix.                                            !
! (3)   Wilkinson-type tridiagonal matrix.                                     !
! (4)   Clement-type tridiagonal matrix.                                       !
! (5)   Legendre orthogonal polynomials                                        !
! (6)   Laguerre orthogonal polynomials                                        !
! (7)   Hermit orthogonal polynomials                                          !
!                                                                              !
! Tests performed:                                                             !
! ---------------                                                              !
!                                                                              !
! For each subroutine defined by the macro "CALLST", the following values are  !
! calculated:                                                                  !
!                                                                              !
! (1)   | T - Z*W*Z' | / ( |T|*N*ULP ) if M=N, otherwise                       !
!       | Z'*T*Z - W | / ( |T|*N*ULP )                                         !
! (2)   | I - Z*Z' | / ( N*ULP ), if M=N, otherwise                            !
!       | I - Z'*Z | / ( N*ULP )                                               !
!                                                                              !
! for the eigenvalue problem T*Z = Z*W.                                        !
!                                                                              !
! Subsets of eigenvalues:                                                      !
! ----------------------                                                       !
!                                                                              !
! The tests where RANGE='I' or RANGE='V' require the values of IL, IU, VL      !
! and VU defined by the macros "EIGVI" (and/or "NRILIU") and "EIGVV" (and/or   !
! "NRVLVU"), respectively. If VL and VU are not set, then their values are     !
! matched to the values set to IL and IU. The number of eigenvalues in         !
! RANGE are indicated by M in the output, together with GAPMIN, which          !
! is defined as the minimum gap between                                        !
!    |W(IL-1)-W(IL)| and |W(IU)-W(IU+1)|                                       !
! or                                                                           !
!    |W(right_of_VL)-W(left_of_VL)| and |W(right_of_VU)-W(left_of_VU)|         !
!                                                                              !
! Data stored in derived data type T_DATA, module SSTEDEFINITIONS:             !
! ---------------------------------------------------------------              !
!                                                                              !
! form    type    size      en   id                                            !
!                                                                              !
!   1    ETYPE   ESIZE   GAMMA   built-in eigenvalue distribution              !
!   2    MTYPE   MSIZE   GAMMA   built-in tridiagonal matrix                   !
!   3        1    NEIG       0   eigenvalue distribution read from file        !
!   4        1     NDE       0   tridiagonal matrix read from file             !
!                                                                              !
!==============================================================================!
!                                                                              !
! The user is welcome to incorporate additional test matrices in the routines  !
!    SSTEDEFEIGV (built-in eigenvalue distributions) and                       !
!    SSTEDEFMTRX (built-in tridiagonal matrices),                              !
! in which case the user is encouraged to send a note to OAMarques@lbl.gov so  !
! we can consider incorporating other test matrices in our distribution.       !
!                                                                              !
!==============================================================================!
!
!.. Scalars ..
LOGICAL :: IEEE
INTEGER :: ILAENV_10, ILAENV_11
INTEGER :: HBANDA = 100, HBANDR = 100, MAXN   = 0, NCASE = 0, &
           NILIU  =   0, NRILIU =   0, NRVLVU = 0, NVLVU = 0 
REAL( KIND=PREC ) :: MINTT = 0.05, TOLEBZ, TOLEGR
!
!.. Static Arrays ..
LOGICAL, DIMENSION( 8 ) :: CALLST = .FALSE. 
LOGICAL, DIMENSION( 7 ) :: DUMP = .FALSE.
INTEGER, DIMENSION( 4 ) :: ISEED = 0
!
!.. Derived Data Types ..
TYPE( EIGV_LIST ), POINTER :: EIGV
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( I_LIST ), POINTER :: ILIU
TYPE( V_LIST ), POINTER :: VLVU
TYPE( T_LIST ), POINTER :: T
!
!.. External Functions ..
INTEGER, EXTERNAL :: ILAENV
REAL( KIND=PREC ), EXTERNAL :: SLAMCH
!
!.. Intrinsic Function ..
INTRINSIC SQRT
!
!.. Executable Statements ......................................................
!
ISEED = ISEED_INIT
!
! Set TOLEBZ, TOLEGR and check for IEEE arithmetic support.
!
TOLEBZ = 2*SLAMCH('Safe minimum')
TOLEGR = 2*SLAMCH('Safe minimum')
ILAENV_10 = ILAENV( 10, 'SSTETESTER', 'N', 1, 2, 3, 4 )
ILAENV_11 = ILAENV( 11, 'SSTETESTER', 'N', 1, 2, 3, 4 )
IEEE = ILAENV_10.EQ.1 .AND. ILAENV_11.EQ.1      
IF ( .NOT.IEEE ) THEN
   WRITE( UNIT=FUOUT, FMT='(A)' ) '* SSTETESTER warning: IEEE = .FALSE.'
   WRITE( UNIT=FUOUT, FMT='(A,I2)' ) &
   '* IEEE NaN arithmetic can be trusted not to trap =', ILAENV_10, &
   '* Infinity arithmetic can be trusted not to trap =', ILAENV_11
END IF
!
! Allocate derived data types.
!
ALLOCATE( EIGV ); NULLIFY( EIGV%NEXT )
ALLOCATE( MTRX ); NULLIFY( MTRX%NEXT )
ALLOCATE( ILIU ); NULLIFY( ILIU%NEXT )
ALLOCATE( VLVU ); NULLIFY( VLVU%NEXT )
ALLOCATE( T ); NULLIFY( T%NEXT )
!
! Parse the data file.
!
CALL SSTEGETDATA( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                  TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, &
                  VLVU, NVLVU, EIGV, MTRX, T )
!
! Perform the tests.
!
CALL SSTERUNTESTS( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT,  &
                   TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, &
                   VLVU, NVLVU, EIGV, MTRX, T )
!
! Deallocate derived data types.
!
DEALLOCATE( EIGV )
DEALLOCATE( MTRX )
DEALLOCATE( ILIU )
DEALLOCATE( VLVU )
DEALLOCATE( T )
!
END PROGRAM SSTETESTER
