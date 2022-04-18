SUBROUTINE SSTEMATGEN( N, D, E, HBANDA, ISEED, WORK, LWORK, &
                       EIGV, MTRX, T, ICASE, DUMP )
!
USE GSTEDEFINITIONS
USE SSTEDEFINITIONS
!
!.. Scalar Arguments ..
INTEGER :: HBANDA, ICASE, LWORK, N
!
!.. Array Arguments ..
LOGICAL :: DUMP( 7 )
INTEGER :: ISEED( 4 )
REAL( KIND=PREC ) :: D( * ), E( * ), WORK( * )
!
!.. Derived Data Type Argument ..
TYPE( EIGV_LIST ), POINTER :: EIGV
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( T_LIST ), POINTER :: T
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMATGEN generates a symmetric matrix A with specified eigenvalues and    !
!  halfbandwith HBANDA and then tridiagonalizes A.                             !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N        (output) INTEGER                                                   !
!           The dimension of the matrix.                                       !
!                                                                              !
!  D        (output) REAL( KIND=PREC ) array, dimension ( N )                  !
!           The N diagonal elements of the tridiagonal matrix.                 !
!                                                                              !
!  E        (output) REAL( KIND=PREC ) array, dimension ( N )                  !
!           The N-1 off-diagonal elements of the tridiagonal matrix in         !
!           elements 1 to N-1, E(N) is set to zero.                            !
!                                                                              !
!  HBANDA   (input) INTEGER                                                    !
!           Sets the halfbandwidth of the symmetric matrix to be generated and !
!           then tridiagonalized, i.e. a matrix with max(1,N*(HBANDA/100))     !
!           subdiagonals is generated.                                         !
!                                                                              !
!  ISEED    (input/output) INTEGER array, dimension ( 4 )                      !
!           Seed for the random number generator. Each entry of ISEED should   !
!           lie between 0 and 4095 inclusive and ISEED(4) should be odd.       !
!                                                                              !
!  WORK     (workspace) REAL( KIND=PREC ) array, dimension ( LWORK )           !
!           Workspace.                                                         !
!                                                                              !
!  LWORK    (input) INTEGER                                                    !
!           Dimension of WORK.                                                 !
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
!  ICASE    (input) INTEGER                                                    !
!           Case number.                                                       !
!                                                                              !
!  DUMP     (input) LOGICAL, dimension ( 5 )                                   !
!           Defines data to be written into files,                             !
!           DUMP( 1 ) : tridiagonal matrix (i,d_i,e_i)                         !
!           DUMP( 2 ) : eigenvalues                                            !
!           DUMP( 3 ) : eigenvectors                                           !
!           DUMP( 4 ) : timing, residuals, orthogonality                       !
!           DUMP( 5 ) : tridiagonal matrix (i,d_i,e_i) in Matlab format        !
!           DUMP( 6 ) : eigenvalues in Matlab format                           !
!           DUMP( 7 ) : eigenvectors in Matlab format                          !
!                                                                              !
!==============================================================================!
!  
!.. Local Scalars ..
CHARACTER( LEN=1 ) :: PACK
INTEGER :: ECOND, EDIST, ESIGN, IA, INFO, IW, IWORK, &
           J, MFORM, MSIZE, MTYPE, NDIAG
REAL( KIND=PREC ) :: GAMMA
!
!.. External Subroutines ..
EXTERNAL HANDLER, SLATMS, SSBTRD, SSYTRD
!
!.. Executable Statements ......................................................
!
WRITE( UNIT=FUOUT, FMT='(A5,I4)' ) 'Case:', ICASE
!
IF ( DUMP( 1 ) ) WRITE( UNIT=FUDUMP( 1 ), FMT='(A,I4)' ) 'Case:', ICASE
IF ( DUMP( 2 ) ) WRITE( UNIT=FUDUMP( 2 ), FMT='(A,I4)' ) 'Case:', ICASE
IF ( DUMP( 3 ) ) WRITE( UNIT=FUDUMP( 3 ), FMT='(A,I4)' ) 'Case:', ICASE
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) ) &
   WRITE( UNIT=FUDUMP( 5 ), FMT='(''% '',A,I5,1X,64(''#''))' ) 'Case:', ICASE
!
N = 0
!
DO
!
!  Get data from T while GAMMA != 0
!
   MFORM = T%DATA%FORM
   MTYPE = T%DATA%TYPE
   MSIZE = T%DATA%SIZE
   ECOND = T%DATA%COND
   EDIST = T%DATA%DIST
   ESIGN = T%DATA%SIGN
   ISEED = T%DATA%SEED
   GAMMA = T%DATA%EN  
!
   IF ( GAMMA == ZERO ) THEN
      WRITE( UNIT=FUOUT, FMT='(4X,A,I2,A,I2,A,I6)' ) &
             'Matrix data: form =', MFORM, ', type =', MTYPE, &
             ', size =', MSIZE
   ELSE
      WRITE( UNIT=FUOUT, FMT='(4X,A,I2,A,I2,A,I6,A,1P,E11.4)' ) &
             'Matrix data: form =', MFORM, ', type =', MTYPE, &
             ', size =', MSIZE, ', glue =', GAMMA 
   END IF
   IF ( DUMP( 1 ) ) THEN
      IF ( GAMMA == ZERO ) THEN
         WRITE( UNIT=FUDUMP(1), &
                FMT='(4X,A,I2,A,I2,A,I6,/,4X,A,4I5)' ) &
                'Matrix data:  form =', MFORM, ', type =', MTYPE, &
                ', size =', MSIZE, 'SEED =', ISEED
      ELSE
         WRITE( UNIT=FUDUMP(1), &
                FMT='(4X,A,I2,A,I2,A,I6,A,1P,E11.4,/,4X,A,4I5)' ) &
                'Matrix data:  form =', MFORM, ', type =', MTYPE, &
                ', size =', MSIZE, ', glue =', GAMMA, &
                'SEED =', ISEED
      END IF
   END IF
!
   SELECT CASE ( MFORM )
!
   CASE ( 1 )
!
!       Built-in eigenvalue distribution.
!
        CALL SSTEDEFEIGV( ECOND, EDIST, ESIGN, MTYPE, ISEED, MSIZE, WORK( 1 ) )
!
   CASE ( 2 )
!
!       Built-in tridiagonal matrix.
!
        CALL SSTEDEFMTRX( MTYPE, MSIZE, D( N+1 ), E( N+1 ) )
!
   CASE ( 3 )
!
!       Eigenvalue distribution read from file.
!
        WORK( 1:MSIZE ) = EIGV%S( 1:MSIZE )
        EIGV => EIGV%NEXT
!
   CASE ( 4 )
!
!       Tridiagonal matrix read from file.
!
        D( 1:MSIZE ) = MTRX%D( 1:MSIZE )
        E( 1:MSIZE ) = MTRX%E( 1:MSIZE )
        MTRX => MTRX%NEXT
!
   END SELECT
!
   IF ( MFORM == 1 .OR. MFORM == 3 ) THEN
!
      IW = 1
      IA = IW + MSIZE
      IWORK = IA + MSIZE*MSIZE
!
!     Generate symmetric matrix, halfbandwidth max(1,MSIZE*(HBANDA/100)).
!
      IF ( HBANDA > 50 ) THEN
         PACK = 'N'  ! No packing.
      ELSE
         PACK = 'B'  ! Store the lower triangle.
      END IF
      NDIAG = MAX( 1, INT( MSIZE*(HBANDA/HNDRD) ) )
      CALL SLATMS( MSIZE, MSIZE, DIST( EDIST ), ISEED, SYMM( ESIGN ), &
                   WORK( IW ), 0, ZERO, ONE, NDIAG, NDIAG, PACK, &
                   WORK( IA ), MSIZE, WORK( IWORK ), INFO )  
      IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'SLATMS' )
!
!     Tridiagonalize matrix.
!
      IF ( PACK == 'N' ) THEN
         CALL SSYTRD( 'L', MSIZE, WORK( IA ), MSIZE, D( N+1 ), E( N+1 ), &
                      WORK( IW ), WORK( IWORK ), LWORK-IWORK+1, INFO )
         IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'SSYTRD' )
      ELSE
         CALL SSBTRD( 'N', 'L', MSIZE, NDIAG, WORK( IA ), MSIZE, D( N+1 ), &
                      E( N+1 ), WORK( IW ), 1, WORK( IW ), INFO )
         IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'SSBTRD' )
      END IF
!
   END IF
!
   N = N + MSIZE
   E( N ) = GAMMA
   IF ( GAMMA == ZERO ) EXIT
   T => T%NEXT
!
END DO
!
IF ( DUMP( 1 ) ) THEN
   WRITE( UNIT=FUDUMP(1), FMT='(5X,A,5X,A,11X,A,/,(I6,1P,2E17.7))' ) &
          'J', 'D( J )', 'E( J )', ( J, D( J ), E( J ), J = 1,N )
END IF
IF ( DUMP( 4 ) ) THEN
   IF ( N == MSIZE ) THEN
      WRITE( UNIT=FUDUMP( 4 ), FMT='(A,I4,A,I6,2X,5(A,I2),A)' ) &
             'Case:', ICASE, ', N = ', N, '(form=', MFORM, ', type=', MTYPE, &
             ', cond=', ECOND, ', dist=', EDIST, ', sign=', ESIGN, ')'
   ELSE
      WRITE( UNIT=FUDUMP( 4 ), FMT='(A,I4,A,I6,2X,5(A,I2),A)' ) &
             'Case:', ICASE, ', N = ', N, '(form=', 9, ', type=', 9, &
             ', cond=', ECOND, ', dist=', EDIST, ', sign=', ESIGN, &
             ', glued matrix)'
   END IF
END IF
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) ) THEN
   WRITE( UNIT=FUDUMP(5), FMT='(A,I5,A)' ) 'N =', N, ';'
   WRITE( UNIT=FUDUMP(5), FMT='(A,I3.3,A)' ) 'N_', ICASE, ' = N;'
END IF
IF ( DUMP( 5 ) ) THEN
   WRITE( UNIT=FUDUMP(5), FMT='(A)' ) 'D = zeros(N,1); E = zeros(N,1);'
   WRITE( UNIT=FUDUMP(5), FMT='(2(A,I5,A,1P,E15.7E3,A))' ) &
          ( 'D(', J, ')=', D( J ), '; ', &
            'E(', J, ')=', E( J ), '; ', &
            J = 1, N )
   WRITE( UNIT=FUDUMP(5), FMT='(2(A,I3.3,A))' ) &
          'D_', ICASE, ' = D; ', 'E_', ICASE, ' = E; clear D E;'
END IF
!
T => T%NEXT
!
END SUBROUTINE SSTEMATGEN
