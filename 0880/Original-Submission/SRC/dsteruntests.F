SUBROUTINE DSTERUNTESTS( HBANDA, HBANDR, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                         TOLEBZ, TOLEGR, CALLST, DUMP, ILIU, NILIU, VLVU, &
                         NVLVU, EIGV, MTRX, T )
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
USE DSTEINTERFACES2, ONLY : DSTEBNDGAP, DSTEMATGEN, DSTEPRNRSLT 
!
!.. Scalar Arguments ..
INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU
REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
!
!.. Array Arguments ..
LOGICAL :: DUMP( 7 ), CALLST( 8 )
!
!.. Derived Data Type Arguments ..
TYPE( EIGV_LIST ), POINTER :: EIGV
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( I_LIST ), POINTER :: ILIU
TYPE( V_LIST ), POINTER :: VLVU
TYPE( T_LIST ), POINTER :: T
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTERUNTESTS runs the LAPACK symmetric tridiagonal eigensolver subroutines  !
!  and checks the results.                                                     !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  HBANDA   (input) INTEGER                                                    !
!           Sets the halfbandwidth of the symmetric matrix to be generated and !
!           then tridiagonalized, i.e. a matrix with max(1,N*(HBANDA/100))     !
!           subdiagonals is generated.                                         !
!                                                                              !
!  HBANDR   (input) INTEGER                                                    !
!           Sets the halfbandwidth of the matrices Z'*T*Z and Z'*Z used in the !
!           tests, i.e. max(1,N*(HBANDR/100)) subdiagonals are computed. If    !
!           0 then the tests are not performed.                                !
!                                                                              !
!  MAXN     (input) INTEGER                                                    !
!           Maximum dimension of the matrix to be tested.                      !
!                                                                              !
!  NCASE    (input) INTEGER                                                    !
!           Number of tridiagonal matrices defined in T.                       !
!                                                                              !
!  NRILIU   (input) INTEGER                                                    !
!           Number of pairs of random indices of the smallest and largest      !
!           eigenvalues to be returned, used only when RANGE='I'.              !
!                                                                              !
!  NRVLVU   (input) INTEGER                                                    !
!           Number of random intervals to be searched for eigenvalues,         !
!           used only when RANGE='V'.                                          !
!                                                                              !
!  MINTT    (input) REAL( KIND=PREC )                                          !
!           Minimum total timing for each subroutine call, i.e.                !
!           time_per_call * number_of_calls > MINTT                            !
!                                                                              !
!  TOLEBZ   (input) REAL( KIND=PREC )                                          !
!           Absolute error tolerance for the eigenvalues to be computed by     !
!           DSTEVX (see the leading comments of DSTEVX for details).           !                                         !
!                                                                              !
!  TOLEGR   (input) REAL( KIND=PREC )                                          !
!           Absolute error tolerance for the eigenvalues to be computed by     !
!           DSTEGR (see the leading comments of DSTEGR for details).           !
!                                                                              !
!  CALLST   (input) LOGICAL array, dimension ( 8 )                             !
!           Flags the subroutines to be tested,                                !
!           CALLST( 1 ) : call DSTEQR with COMPZ='V'                           !
!           CALLST( 2 ) : call DSTEVX with RANGE='A'                           !
!           CALLST( 3 ) : call DSTEVX with RANGE='I'                           !
!           CALLST( 4 ) : call DSTEVX with RANGE='V'                           !
!           CALLST( 5 ) : call DSTEDC with COMPZ='I'                           !
!           CALLST( 6 ) : call DSTEGR with RANGE='A'                           !
!           CALLST( 7 ) : call DSTEGR with RANGE='I'                           !
!           CALLST( 8 ) : call DSTEGR with RANGE='V'                           !
!                                                                              !
!  DUMP     (input) LOGICAL array, dimension ( 5 )                             !
!           Defines data to be written into files,                             !
!           DUMP( 1 ) : tridiagonal matrix (i,d_i,e_i)                         !
!           DUMP( 2 ) : eigenvalues                                            !
!           DUMP( 3 ) : eigenvectors                                           !
!           DUMP( 4 ) : timing, residuals, orthogonality                       !
!           DUMP( 5 ) : tridiagonal matrix (i,d_i,e_i) in Matlab format        !
!           DUMP( 6 ) : eigenvalues in Matlab format                           !
!           DUMP( 7 ) : eigenvectors in Matlab format                          !
!                                                                              !
!  ILIU     (input) I_LIST (derived data type)                                 !
!           Indices (in ascending order) of the smallest and largest           !
!           eigenvalues to be returned, used only when RANGE='I'.              !
!                                                                              !
!  NILIU    (input) INTEGER                                                    !
!           Number of pairs of indices in ILIU.                                !
!                                                                              !
!  VLVU     (input) V_LIST (derived data type)                                 !
!           Lower and upper bounds of intervals to be searched for             !
!           eigenvalues, used only when RANGE='V'.                             !
!                                                                              !
!  NVLVU    (input) INTEGER                                                    !
!           Number of intervals in VLVU.                                       !
!                                                                              !
!  EIGV     (input) EIGV_LIST (derived data type)                              !
!           Eigenvalue distributions read from files.                          !
!                                                                              !
!  MTRX     (input) MTRX_LIST (derived data type)                              !
!           Tridiagonal matrices read from files.                              !
!                                                                              !
!  T        (input) T_LIST (derived data type)                                 !
!           Properties of the tridiagonal matrices to be used in the tests.    !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
LOGICAL :: REPEAT, TESTILIU, TESTVLVU
INTEGER :: IBLOCK, ICASE, ID, IE, IFAIL, IINT, IIWORK, IL, INFO, &
           IRWORK, ISPLIT, ISUPPZ, IU, IW, IWEBZ, LDZ, LIWORK, &
           LRWORK, M, N, NRUN, NSPLIT
REAL( KIND=PREC ) :: FLOPS, GAPMIN, TIME0, TIME1, TIME2, TOL, TTIME, VL, VU
!
!.. Static Arrays ..
INTEGER :: ISEED( 4 )
REAL( KIND=PREC ) :: RESULT( 3 )
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: IWORK( : ), ILIU_A( :,: )
REAL( KIND=PREC ), ALLOCATABLE :: D( : ), E( : ), RWORK( : ), &
                                  VLVU_A( :,: ), Z( :,: )
!
!.. Local Derived Data Types ..
TYPE( EIGV_LIST ), POINTER :: FIRST_EIGV
TYPE( MTRX_LIST ), POINTER :: FIRST_MTRX
TYPE( I_LIST ), POINTER :: FIRST_ILIU 
TYPE( V_LIST ), POINTER :: FIRST_VLVU
TYPE( T_LIST ), POINTER :: FIRST_T
!
!.. External Functions ..
REAL( KIND=PREC ), EXTERNAL :: DLAMCH, DSTEFLOPS, DSTETIMER
!
!.. External Subroutines ..
EXTERNAL DSTEDC, DSTEGR, DSTEQR, DSTEVX
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
WRITE( UNIT=FUOUT, FMT='(A,A)' ) '========== Running test cases ============'
!
IF ( DUMP( 1 ) ) OPEN( UNIT=FUDUMP( 1 ), FILE=FDUMP( 1 ) )
IF ( DUMP( 2 ) ) OPEN( UNIT=FUDUMP( 2 ), FILE=FDUMP( 2 ) )
IF ( DUMP( 3 ) ) OPEN( UNIT=FUDUMP( 3 ), FILE=FDUMP( 3 ) )
IF ( DUMP( 4 ) ) OPEN( UNIT=FUDUMP( 4 ), FILE=FDUMP( 4 ) )
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) ) THEN
   OPEN( UNIT=FUDUMP( 5 ), FILE=FDUMP( 5 ) )
END IF
!
ALLOCATE( FIRST_EIGV ); NULLIFY( FIRST_EIGV%NEXT ); FIRST_EIGV => EIGV
ALLOCATE( FIRST_MTRX ); NULLIFY( FIRST_MTRX%NEXT ); FIRST_MTRX => MTRX
ALLOCATE( FIRST_ILIU ); NULLIFY( FIRST_ILIU%NEXT ); FIRST_ILIU => ILIU
ALLOCATE( FIRST_VLVU ); NULLIFY( FIRST_VLVU%NEXT ); FIRST_VLVU => VLVU
ALLOCATE( FIRST_T ); NULLIFY( FIRST_T%NEXT ); FIRST_T => T
ALLOCATE( ILIU_A( 2,MAX(1,NILIU+NRILIU+NRVLVU) ) )
ALLOCATE( VLVU_A( 2,MAX(1,NVLVU+NRILIU+NRVLVU) ) )
!
N = 0
IL = 0; IU = 0
VL = ZERO; VU = ZERO
TESTILIU = CALLST( 3 ) .OR. CALLST( 7 )
TESTVLVU = CALLST( 4 ) .OR. CALLST( 8 )
ISEED = ISEED_INIT
!
DO IINT = 1, NILIU
   ILIU_A( 1,IINT ) = ILIU%IL
   ILIU_A( 2,IINT ) = ILIU%IU
   ILIU => ILIU%NEXT
END DO
DO IINT = 1, NVLVU
   VLVU_A( 1,IINT ) = VLVU%VL
   VLVU_A( 2,IINT ) = VLVU%VU
   VLVU => VLVU%NEXT
END DO
!
LIWORK = MAX( 3 + 5*MAXN, & !* Work space required by DSTEDC *!
              10*MAXN     ) !* Work space required by DSTEGR *!
LRWORK = MAX( 1 + 4*MAXN + MAXN**2, & !* Work space required by DSTEDC *!
              18*MAXN,              & !* Work space required by DSTEGR *!
              2*MAXN-2,             & !* Work space required by DSTEQR *!
              5*MAXN,               & !* Work space required by DSTEVX *!
              (MAXN+1)*MAXN         ) !* Work space required by DSTECHKRSLT *!
!
LIWORK = LIWORK + MAXN*3
LRWORK = LRWORK + MAXN*4
!
ALLOCATE( D( MAXN ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'D (in DSTERUNTESTS)' )
ALLOCATE( E( MAXN ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'E (in DSTERUNTESTS)' ) 
ALLOCATE( IWORK( LIWORK ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'IWORK (in DSTERUNTESTS)' ) 
ALLOCATE( RWORK( LRWORK ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'RWORK (in DSTERUNTESTS)' ) 
ALLOCATE( Z( MAXN,MAXN ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'Z (in DSTERUNTESTS)' ) 
!
! If PAPI is used, initialize it.
!
TIME0 = DSTETIMER( 0 )
CALL DSTEPRNFLOP( 0, ZERO, 0, 0 )
!
! Run test cases.
!
DO ICASE = 1, NCASE
!
   INFO = 0
   RESULT = 0
!
!  Generate the test matrix.
!
   CALL DSTEMATGEN( N, D, E, HBANDA, ISEED, RWORK, LRWORK, &
                    EIGV, MTRX, T, ICASE, DUMP )
!
!  Set pointers for the real and integer work arrays.
!
   ID = 1 
   IE = ID + N 
   IW = IE + N 
   IWEBZ = IW + N 
   IRWORK = IWEBZ + N 
!
   IFAIL = 1
   IBLOCK = 1
   ISPLIT = IBLOCK + N
   ISUPPZ = 1
   IIWORK = ISUPPZ + N*2
!
   LDZ = N
!
!  Compute pairs IL,IU and VL,VU for cases RANGE="V"/"I".
!
   IF ( TESTILIU .OR. TESTVLVU ) THEN
!
!     Compute the eigenvalues of the test matrix using bisection. This 
!     information will be used later to compare the results.
!
      TOL = 2*DLAMCH('Safe minimum')
      CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
      CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
      CALL DSTEBZ( 'A', 'E', N, VL, VU, IL, IU, TOLEBZ, RWORK( ID ), &
                   RWORK( IE ), M, NSPLIT, RWORK( IWEBZ ), IWORK( IBLOCK ), &
                   IWORK( ISPLIT ), RWORK( IRWORK ), IWORK( IIWORK ), INFO )
!
!     Choose random values for IL and IU.
!
      CALL DSTERNDILIU( N, ISEED, NRILIU+NRVLVU, ILIU_A( 1,NILIU+1 ) )
!
!     Choose values for VL and VU. If NVLVU=NRVLVU=0 then the values of VL 
!     and VU are computed based on the values of IL and IU, otherwise the
!     values of VL and VU are chosen randomly. 
!
      IF ( MAX( NVLVU, NRVLVU ) == 0 ) THEN
         CALL DSTERNDVLVU( N, RWORK( IWEBZ ), NRILIU, VLVU_A( 1,1 ), &
                           ILIU_A( 1,1 ) )
         NRVLVU = NRILIU
      ELSE
         CALL DSTERNDVLVU( N, RWORK( IWEBZ ), NRVLVU, VLVU_A( 1,NVLVU+1 ), &
                           ILIU_A( 1,NILIU+NRILIU+1 ) )
      END IF 
!
   END IF
!
!  Test DSTEQR( COMPZ='I' ) ....................................................
!
   IF ( CALLST( 1 ) ) THEN
      NRUN = 0
      TTIME = ZERO
      REPEAT = .TRUE.
      TIME0 = DSTETIMER( 1 )
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )
         TIME1 = DSTETIMER( 1 )
         CALL DSTEQR( 'I', N, RWORK( ID ), RWORK( IE ), &
                      Z, LDZ, RWORK( IRWORK ), INFO )
         TIME2 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME2 < TIME1 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 IF ( NRUN == 1 ) THEN
                    TTIME = TTIME + ( TIME2 - TIME1 )
                 ELSE
                    TTIME = TTIME + ( TIME2 - TIME0 )
                 END IF
                 REPEAT = ( TTIME < MINTT )
         END IF
      END DO
      RESULT( 1 ) = TTIME / MAX( 1, NRUN )
      CALL DSTECHKRSLT( N, D, E, N, RWORK( ID ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), INFO )
      CALL DSTEPRNRSLT( N, N, RWORK( ID: ), Z, LDZ, RESULT, &
                        INFO, 1, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 1, ICASE )
   END IF
!
!  Test DSTEVX( RANGE='A' ) ....................................................
!
   IF ( CALLST( 2 ) ) THEN
      NRUN = 0
      TTIME = ZERO
      REPEAT = .TRUE.
      TIME0 = DSTETIMER( 1 )
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )
         TIME1 = DSTETIMER( 1 )
         CALL DSTEVX( 'V', 'A', N, RWORK( ID ), RWORK( IE ), VL, VU, IL, IU, &
                      TOLEBZ, M, RWORK( IW ), Z, LDZ, RWORK( IRWORK ), &
                      IWORK( IIWORK ), IWORK( IFAIL ), INFO)
         TIME2 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME2 < TIME1 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 IF ( NRUN == 1 ) THEN
                    TTIME = TTIME + ( TIME2 - TIME1 )
                 ELSE
                    TTIME = TTIME + ( TIME2 - TIME0 )
                 END IF
                 REPEAT = ( TTIME < MINTT )
         END IF
      END DO
      RESULT( 1 ) = TTIME / MAX( 1, NRUN )
      CALL DSTECHKRSLT( N, D, E, M, RWORK( IW ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), INFO )
      CALL DSTEPRNRSLT( N, N, RWORK( IW: ), Z, LDZ, RESULT, &
                        INFO, 2, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 2, ICASE )
   END IF
!
!  Test DSTEVX( RANGE='I' ) ....................................................
!
   IF ( CALLST( 3 ) ) THEN
      DO IINT = 1, NILIU+NRILIU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         IL = ILIU_A( 1,IINT )
         IU = ILIU_A( 2,IINT )
         TIME0 = DSTETIMER( 1 )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME1 = DSTETIMER( 1 )
            CALL DSTEVX( 'V', 'I', N, RWORK( ID ), RWORK( IE ), VL, VU, IL, IU, &
                         TOLEBZ, M, RWORK( IW ), Z, LDZ, RWORK( IRWORK ), &
                         IWORK( IIWORK ), IWORK( IFAIL ), INFO)
            TIME2 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME2 < TIME1 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    IF ( NRUN == 1 ) THEN
                       TTIME = TTIME + ( TIME2 - TIME1 )
                    ELSE
                       TTIME = TTIME + ( TIME2 - TIME0 )
                    END IF
                    REPEAT = ( TTIME < MINTT )
            END IF
         END DO
         RESULT( 1 ) = TTIME / MAX( 1, NRUN )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, ILO=IL, IUO=IU )
         CALL DSTECHKRSLT( N, D, E, M, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), INFO )
         CALL DSTEPRNRSLT( N, M, RWORK( IW: ), Z, LDZ, RESULT, INFO, 3, &
                           ICASE, DUMP, IINTO=IINT, ILO=IL, IUO=IU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 3, ICASE )
      END DO
   END IF
!
!  Test DSTEVX( RANGE='V' ) ....................................................
!
   IF ( CALLST( 4 ) ) THEN
      DO IINT = 1, NVLVU+NRVLVU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         VL = VLVU_A( 1,IINT )
         VU = VLVU_A( 2,IINT )
         TIME0 = DSTETIMER( 1 )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME1 = DSTETIMER( 1 )
            CALL DSTEVX( 'V', 'V', N, RWORK( ID ), RWORK( IE ), VL, VU, IL, IU, &
                         TOLEBZ, M, RWORK( IW ), Z, LDZ, RWORK( IRWORK ), &
                         IWORK( IIWORK ), IWORK( IFAIL ), INFO)
            TIME2 = DSTETIMER( 1 )  
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME2 < TIME1 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    IF ( NRUN == 1 ) THEN
                       TTIME = TTIME + ( TIME2 - TIME1 )
                    ELSE
                       TTIME = TTIME + ( TIME2 - TIME0 )
                    END IF
                    REPEAT = ( TTIME < MINTT )
            END IF
         END DO
         RESULT( 1 ) = TTIME / MAX( 1, NRUN )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, VLO=VL, VUO=VU )
         CALL DSTECHKRSLT( N, D, E, M, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), INFO )
         CALL DSTEPRNRSLT( N, M, RWORK( IW: ), Z, LDZ, RESULT, INFO, 4, &
                           ICASE, DUMP, IINTO=IINT, VLO=VL, VUO=VU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 4, ICASE )
      END DO
   END IF
!
!  Test DSTEDC( COMPZ='I' ) ....................................................
!
   IF ( CALLST( 5 ) ) THEN
      NRUN = 0
      TTIME = ZERO
      REPEAT = .TRUE.
      TIME0 = DSTETIMER( 1 )
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )
         TIME1 = DSTETIMER( 1 )
         CALL DSTEDC( 'I', N, RWORK( ID ), RWORK( IE ), Z, LDZ, RWORK( IRWORK ), &
                      LRWORK-N*2, IWORK, LIWORK, INFO )
         TIME2 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME2 < TIME1 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 IF ( NRUN == 1 ) THEN
                    TTIME = TTIME + ( TIME2 - TIME1 )
                 ELSE
                    TTIME = TTIME + ( TIME2 - TIME0 )
                 END IF
                 REPEAT = ( TTIME < MINTT )
         END IF
      END DO
      RESULT( 1 ) = TTIME / MAX( 1, NRUN )
      CALL DSTECHKRSLT( N, D, E, N, RWORK( ID ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), INFO )
      CALL DSTEPRNRSLT( N, N, RWORK( ID: ), Z, LDZ, RESULT, INFO, &
                        5, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 5, ICASE )
   END IF
!    
!  Test DSTEGR( RANGE='A' ) ....................................................
!
   IF ( CALLST( 6 ) ) THEN
      NRUN = 0       
      TTIME = ZERO
      REPEAT = .TRUE.
      TIME0 = DSTETIMER( 1 )
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )
         TIME1 = DSTETIMER( 1 )
         CALL DSTEGR( 'V', 'A', N, RWORK( ID ), RWORK( IE ), VL, VU, &
                      IL, IU, TOLEGR, M, RWORK( IW ), Z, LDZ, IWORK( ISUPPZ ), &
                      RWORK( IRWORK ), LRWORK-N*3, IWORK( IIWORK ), &
                      LIWORK-N*2, INFO )
         TIME2 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME2 < TIME1 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 IF ( NRUN == 1 ) THEN
                    TTIME = TTIME + ( TIME2 - TIME1 )
                 ELSE
                    TTIME = TTIME + ( TIME2 - TIME0 )
                 END IF
                 REPEAT = ( TTIME < MINTT )
         END IF
      END DO
      RESULT( 1 ) = TTIME / MAX( 1, NRUN )
      CALL DSTECHKRSLT( N, D, E, M, RWORK( IW ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), INFO )
      CALL DSTEPRNRSLT( N, M, RWORK( IW: ), Z, LDZ, RESULT, INFO, &
                        6, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 6, ICASE )
   END IF
!
!  Test DSTEGR( RANGE='I' ) ....................................................
!
   IF ( CALLST( 7 ) ) THEN
      DO IINT = 1, NILIU+NRILIU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         IL = ILIU_A( 1,IINT )
         IU = ILIU_A( 2,IINT )
         TIME0 = DSTETIMER( 1 )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME1 = DSTETIMER( 1 )
            CALL DSTEGR( 'V', 'I', N, RWORK( ID ), RWORK( IE ), VL, VU, &
                         IL, IU, TOLEGR, M, RWORK( IW ), Z, LDZ, IWORK( ISUPPZ ), &
                         RWORK( IRWORK ), LRWORK-N*3, IWORK( IIWORK ), &
                         LIWORK-N*2, INFO )
            TIME2 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME2 < TIME1 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    IF ( NRUN == 1 ) THEN
                       TTIME = TTIME + ( TIME2 - TIME1 )
                    ELSE
                       TTIME = TTIME + ( TIME2 - TIME0 )
                    END IF
                    REPEAT = ( TTIME < MINTT )
            END IF
         END DO
         RESULT( 1 ) = TTIME / MAX( 1, NRUN )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, ILO=IL, IUO=IU )
         CALL DSTECHKRSLT( N, D, E, M, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), INFO )
         CALL DSTEPRNRSLT( N, M, RWORK( IW: ), Z, LDZ, RESULT, INFO, 7, &
                           ICASE, DUMP, IINTO=IINT, ILO=IL, IUO=IU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 7, ICASE )
      END DO
   END IF
!
!  Test DSTEGR( RANGE='V' ) ....................................................
!
   IF ( CALLST( 8 ) ) THEN
      DO IINT = 1, NVLVU+NRVLVU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         VL = VLVU_A( 1,IINT )
         VU = VLVU_A( 2,IINT )
         TIME0 = DSTETIMER( 1 )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME1 = DSTETIMER( 1 )
            CALL DSTEGR( 'V', 'V', N, RWORK( ID ), RWORK( IE ), VL, VU, &
                         IL, IU, TOLEGR, M, RWORK( IW ), Z, LDZ, IWORK( ISUPPZ ), &
                         RWORK( IRWORK ), LRWORK-N*3, IWORK( IIWORK ), &
                         LIWORK-N*2, INFO )
            TIME2 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME2 < TIME1 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    IF ( NRUN == 1 ) THEN
                       TTIME = TTIME + ( TIME2 - TIME1 )
                    ELSE
                       TTIME = TTIME + ( TIME2 - TIME0 )
                    END IF
                    REPEAT = ( TTIME < MINTT )
            END IF
         END DO
         RESULT( 1 ) = TTIME / MAX( 1, NRUN )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, VLO=VL, VUO=VU )
         CALL DSTECHKRSLT( N, D, E, M, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), INFO )
         CALL DSTEPRNRSLT( N, M, RWORK( IW: ), Z, LDZ, RESULT, INFO, 8, &
                           ICASE, DUMP, IINTO=IINT, VLO=VL, VUO=VU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 8, ICASE )
      END DO
   END IF
!
END DO
! 
IF ( DUMP( 1 ) ) CLOSE( UNIT=FUDUMP( 1 ) )
IF ( DUMP( 2 ) ) CLOSE( UNIT=FUDUMP( 2 ) )
IF ( DUMP( 3 ) ) CLOSE( UNIT=FUDUMP( 3 ) )
IF ( DUMP( 4 ) ) CLOSE( UNIT=FUDUMP( 4 ) )
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) ) THEN
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'clear N;'
   CLOSE( UNIT=FUDUMP( 5 ) )
END IF
CALL DSTEPRNFLOP( 0, ZERO, 0, NCASE )
!
EIGV => FIRST_EIGV
MTRX => FIRST_MTRX
ILIU => FIRST_ILIU
VLVU => FIRST_VLVU
T => FIRST_T
!
END SUBROUTINE DSTERUNTESTS
