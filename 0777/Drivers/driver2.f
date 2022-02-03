C MAIN ROUTINE TO TEST POLSYS1H.
C
C THIS ROUTINE REQUIRES ONE INPUT FILE, DATA2 .
C
C A SAMPLE INPUT FILE AND ASSOCIATED OUTPUT ARE GIVEN
C IN THE COMMENTS THAT FOLLOW.  THIS SAMPLE PROBLEM IS
C CITED IN THE HOMPACK REPORT.
C
C***** SAMPLE INPUT DATA (READ FROM THE FILE 'DATA2'):
C ' TWO QUADRICS, NO SOLUTIONS AT INFINITY, TWO REAL SOLUTIONS.'
C &PROBLEM  
C       IFLGHM = 1
C       IFLGSC = 1
C       TOTDG  = 4
C       MAXT = 6
C       EPSBIG = 1.D-04
C       EPSSML = 1.D-14    
C       SSPAR(5) = 1.D+00    
C       NUMRR = 10
C       N = 2 /
C 00006                     NUMTRM(1)
C 00002                     DEG(1,1,1)
C 00000                     DEG(1,2,1)
C            -.00098D+00
C 00000                     DEG(1,1,2)
C 00002                     DEG(1,2,2)
C            978000.D+00
C 00001                     DEG(1,1,3)
C 00001                     DEG(1,2,3)
C               -9.8D+00
C 00001                     DEG(1,1,4)
C 00000                     DEG(1,2,4)
C             -235.0D+00
C 00000                     DEG(1,1,5)
C 00001                     DEG(1,2,5)
C            88900.0D+00
C 00000                     DEG(1,1,6)
C 00000                     DEG(1,2,6)
C             -1.000D+00
C 00006                     NUMTRM(2)
C 00002                     DEG(2,1,1)
C 00000                     DEG(2,2,1)
C             -.0100D+00
C 00000                     DEG(2,1,2)
C 00002                     DEG(2,2,2)
C             -.9840D+00
C 00001                     DEG(2,1,3)
C 00001                     DEG(2,2,3)
C             -29.70D+00
C 00001                     DEG(2,1,4)
C 00000                     DEG(2,2,4)
C             .00987D+00
C 00000                     DEG(2,1,5)
C 00001                     DEG(2,2,5)
C             -.1240D+00
C 00000                     DEG(2,1,6)
C 00000                     DEG(2,2,6)
C             -.2500D+00
C***** END OF SAMPLE INPUT DATA.
C
C***** ASSOCIATED SAMPLE OUTPUT (WRITTEN TO THE FILE 'OUTHP.DAT'):
C
!     POLSYS1H TEST ROUTINE 7/7/95
!
! TWO QUADRICS, NO SOLUTIONS AT INFINITY, TWO REAL SOLUTIONS.            
!
! IF IFLGHM=1, HOMOGENEOUS; IF IFLGHM=0, INHOMOGENEOUS; IFLGHM= 1
!
! IF IFLGSC=1, SCLGNP USED; IF IFLGSC=0, NO SCALING; IFLGSC=    1
!
! TOTDG=    4          MAXT=    6
!
! EPSBIG, EPSSML =  1.000000000000000E-04  1.000000000000000E-14
!
! SSPAR(5) =  1.000000000000000E+00
!
! NUMBER OF EQUATIONS =    2
!
! NUMBER OF RECALLS WHEN IFLAG=3:   10
!
!
!  ****** COEFFICIENT TABLEAU ******
!
!  NUMT( 1) =    6
!  KDEG( 1, 1, 1) =    2
!  KDEG( 1, 2, 1) =    0
!  COEF( 1, 1) =  -9.800000000000000E-04
!  KDEG( 1, 1, 2) =    0
!  KDEG( 1, 2, 2) =    2
!  COEF( 1, 2) =   9.780000000000000E+05
!  KDEG( 1, 1, 3) =    1
!  KDEG( 1, 2, 3) =    1
!  COEF( 1, 3) =  -9.800000000000001E+00
!  KDEG( 1, 1, 4) =    1
!  KDEG( 1, 2, 4) =    0
!  COEF( 1, 4) =  -2.350000000000000E+02
!  KDEG( 1, 1, 5) =    0
!  KDEG( 1, 2, 5) =    1
!  COEF( 1, 5) =   8.890000000000000E+04
!  KDEG( 1, 1, 6) =    0
!  KDEG( 1, 2, 6) =    0
!  COEF( 1, 6) =  -1.000000000000000E+00
!
!  NUMT( 2) =    6
!  KDEG( 2, 1, 1) =    2
!  KDEG( 2, 2, 1) =    0
!  COEF( 2, 1) =  -1.000000000000000E-02
!  KDEG( 2, 1, 2) =    0
!  KDEG( 2, 2, 2) =    2
!  COEF( 2, 2) =  -9.840000000000000E-01
!  KDEG( 2, 1, 3) =    1
!  KDEG( 2, 2, 3) =    1
!  COEF( 2, 3) =  -2.970000000000000E+01
!  KDEG( 2, 1, 4) =    1
!  KDEG( 2, 2, 4) =    0
!  COEF( 2, 4) =   9.870000000000000E-03
!  KDEG( 2, 1, 5) =    0
!  KDEG( 2, 2, 5) =    1
!  COEF( 2, 5) =  -1.240000000000000E-01
!  KDEG( 2, 1, 6) =    0
!  KDEG( 2, 2, 6) =    0
!  COEF( 2, 6) =  -2.500000000000000E-01
!
!
!
!  IFLG1 =   11
!
!  PATH NUMBER =    1
!
!  FINAL VALUES FOR PATH
!
!  ARCLEN =   1.005533190562901E+01
!  NFE =   53
!  IFLG2 =  1
! REAL, FINITE SOLUTION
!  LAMBDA =  1.000000000000003E+00
! X( 1) =  2.342338519591276E+03  8.841149143431121E-13
! X( 2) = -7.883448240941412E-01 -9.356862757018485E-16
!
! X( 3) = -9.493594594086552E-03 -1.064475509002627E-03
!
!
!  PATH NUMBER =    2
!
!  FINAL VALUES FOR PATH
!
!  ARCLEN =   1.721129286057142E+00
!  NFE =   37
!  IFLG2 =  1
! COMPLEX, FINITE SOLUTION
!  LAMBDA =  1.000000000000006E+00
! X( 1) =  1.614785792344189E-02  1.684969554988811E+00
! X( 2) =  2.679947396144760E-04  4.428029939736605E-03
!
! X( 3) = -3.819489729424030E-01  3.720689434572830E-01
!
!
!  PATH NUMBER =    3
!
!  FINAL VALUES FOR PATH
!
!  ARCLEN =   2.023295279367267E+00
!  NFE =   35
!  IFLG2 =  1
! COMPLEX, FINITE SOLUTION
!  LAMBDA =  1.000000000000000E+00
! X( 1) =  1.614785792343521E-02 -1.684969554988812E+00
! X( 2) =  2.679947396144598E-04 -4.428029939736611E-03
!
! X( 3) = -3.293704938476598E-01  5.566197755230126E-01
!
!
!  PATH NUMBER =    4
!
!  FINAL VALUES FOR PATH
!
!  ARCLEN =   4.163266156958467E+00
!  NFE =   46
!  IFLG2 =  1
! REAL, FINITE SOLUTION
!  LAMBDA =  9.999999999999998E-01
! X( 1) =  9.089212296153869E-02  1.153793567884107E-16
! X( 2) = -9.114970981974997E-02  1.887399041592030E-17
!
! X( 3) = -5.736733957279616E-02  1.362436637092185E-01
!
!
! TOTAL NFE OVER ALL PATHS =        171
C
C***** END OF ASSOCIATED SAMPLE OUTPUT.
C
C *************************************************************
C
C  PROGRAM DESCRIPTION:  1. READS IN DATA.
C                        2. GENERATES POLSYS1H INPUT.
C                        3. CALLS POLSYS1H.
C                        4. WRITES POLSYS1H OUTPUT.
C
C DIMENSIONS SHOULD BE SET AS FOLLOWS:
C
C     DIMENSION NUMT(NN),COEF(NN,MMAXT),KDEG(NN,NN+1,MMAXT)
C     DIMENSION IFLG2(TTOTDG)
C     DIMENSION LAMBDA(TTOTDG),ROOTS(2,NN+1,TTOTDG),ARCLEN(TTOTDG),
C    & NFE(TTOTDG)
C
C WHERE:
C    N   IS THE NUMBER OF EQUATIONS.  NN .GE. N.
C    MAXT  IS THE MAXIMUM NUMBER OF TERMS IN ANY ONE EQUATION.
C       MMAXT  .GE.  MAXT.
C    TOTDG  IS THE TOTAL DEGREE OF THE SYSTEM.  TTOTDG .GE. TOTDG.
C
C THIS TEST CODE HAS DIMENSIONS SET AS FOLLOWS:
C
C NN=10, MMAXT=30, TTOTDG=1024
C
      PROGRAM TESTP
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90, ONLY : POLSYS1H
      INTEGER, PARAMETER:: NN=10,MMAXT=30,TTOTDG=1024
      INTEGER:: IFLG1,IFLG2(TTOTDG),IFLGHM,IFLGSC,ITOTIT,J,K,
     &  KDEG(NN,NN+1,MMAXT),L,M,MAXT,N,NFE(TTOTDG),NP1,NT,
     &  NUMRR,NUMT(NN),TOTDG
      REAL (KIND=R8):: ARCLEN(TTOTDG),COEF(NN,MMAXT),EPSBIG,EPSSML,
     &  LAMBDA(TTOTDG),ROOTS(2,NN+1,TTOTDG),SSPAR(8)
      CHARACTER (LEN=72):: TITLE
! If using a subroutine library of the HOMPACK90 subroutines rather than
! the MODULE HOMPACK90 (as above), then the following INTERFACE
! statements are necessary.
!     INTERFACE
!       SUBROUTINE POLSYS1H(N,NUMT,COEF,KDEG,IFLG1,IFLG2,EPSBIG,EPSSML,
!    &     SSPAR,NUMRR,LAMBDA,ROOTS,ARCLEN,NFE)
!       USE HOMOTOPY
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN):: N,NUMT(:),KDEG(:,:,:),NUMRR
!       REAL (KIND=R8), INTENT(IN):: COEF(:,:),EPSBIG,EPSSML
!       INTEGER, INTENT(IN OUT):: IFLG1,IFLG2(:)
!       REAL (KIND=R8), INTENT(IN OUT):: SSPAR(8)
!       REAL (KIND=R8), INTENT(OUT):: LAMBDA(:),ROOTS(:,:,:),ARCLEN(:)
!       INTEGER, INTENT(OUT):: NFE(:)
!       END SUBROUTINE POLSYS1H
!     END INTERFACE
C
      NAMELIST /PROBLEM/ IFLGHM,IFLGSC,TOTDG,MAXT,N,
     &  EPSBIG,EPSSML,SSPAR,NUMRR
C
      OPEN (UNIT=7,FILE='DATA2',ACTION='READ',POSITION='REWIND',
     &  DELIM='APOSTROPHE',STATUS='OLD')
      OPEN (UNIT=6,FILE='RES2.OUT',ACTION='WRITE',
     &  DELIM='APOSTROPHE',STATUS='REPLACE')
C
      SSPAR(1:8) = 0.0
      READ (7,*) TITLE
      WRITE (6,10) TITLE
  10  FORMAT(5X,'POLSYS1H TEST ROUTINE 7/7/95',//,A72)
C
      READ (7, NML=PROBLEM)
C
      WRITE (6,100) IFLGHM
 100  FORMAT(/
     &' IF IFLGHM=1, HOMOGENEOUS; IF IFLGHM=0, INHOMOGENEOUS; IFLGHM='
     & ,I2)
      WRITE (6,102) IFLGSC
 102  FORMAT(/
     &' IF IFLGSC=1, SCLGNP USED; IF IFLGSC=0, NO SCALING; IFLGSC=',I5)
      WRITE (6,104) TOTDG,MAXT
 104  FORMAT(/,' TOTDG=',I5,10X,'MAXT=',I5)
      WRITE (6,106) EPSBIG,EPSSML,SSPAR(5),N,NUMRR
 106  FORMAT(/,' EPSBIG, EPSSML =',2ES22.14,
     &       //,' SSPAR(5) =',ES22.14,
     &       //,' NUMBER OF EQUATIONS =',I5,
     &       //,' NUMBER OF RECALLS WHEN IFLAG=3:',I5)
C
      NP1=N+1
C
C NOTE THAT THE DEGREES OF VARIABLES IN EACH TERM OF EACH EQUATION
C ARE DEFINED BY THE FOLLOWING INDEXING SCHEME:
C
C     KDEG(J,  L,  K)
C
C          ^   ^   ^
C
C          E   V   T
C          Q   A   E
C          U   R   R
C          A   I   M
C          T   A
C          I   B
C          O   L
C          N   E
C
      WRITE(6,200)
 200  FORMAT(//,'  ****** COEFFICIENT TABLEAU ******')
      KDEG = 0  !SET UNUSED DEGREES TO ZERO
      EQN: DO J=1,N
        READ (7,1000) NUMT(J)
        WRITE (6,210) J,NUMT(J)
 210    FORMAT(/,'  NUMT(',I2,') =',I5)
        NT=NUMT(J)
        TERMS: DO K=1,NT
          VARS: DO L=1,N
            READ (7,1000) KDEG(J,L,K)
            WRITE (6,220) J,L,K,KDEG(J,L,K)
 220        FORMAT('  KDEG(',I2,',',I2,',',I2,') =',I5)
          END DO VARS
          READ (7,2000) COEF(J,K)
          WRITE (6,230) J,K,COEF(J,K)
 230      FORMAT('  COEF(',I2,',',I2,') =',ES22.14)
        END DO TERMS
      END DO EQN
      WRITE (6,FMT="(//)")
C
      IFLG1=10*IFLGHM+IFLGSC
      DO M=1,TOTDG
        IFLG2(M)=-2
      END DO
      CALL POLSYS1H(N,NUMT(1:N),COEF(1:N,1:MAXT),
     &  KDEG(1:N,1:N+1,1:MAXT),IFLG1,IFLG2(1:TOTDG),EPSBIG,EPSSML,
     &  SSPAR,NUMRR,LAMBDA(1:TOTDG),ROOTS(1:2,1:N+1,1:TOTDG),
     &  ARCLEN(1:TOTDG),NFE(1:TOTDG))
C
      WRITE (6,240) IFLG1
 240  FORMAT('  IFLG1 =',I5,/)
      ITOTIT = SUM(NFE(1:TOTDG))
      DO M=1,TOTDG
        WRITE (6,260) M
 260    FORMAT('  PATH NUMBER =',I5,//'  FINAL VALUES FOR PATH'/)
        WRITE (6,280) ARCLEN(M)
 280    FORMAT('  ARCLEN =',ES22.14)
        WRITE (6,290) NFE(M)
 290    FORMAT('  NFE =',I5)
        WRITE (6,300) IFLG2(M)
 300    FORMAT('  IFLG2 =',I3)
C
C   DESIGNATE SOLUTIONS "REAL" OR "COMPLEX"
C
        IF (ANY(ABS(ROOTS(2,1:N,M)) .GE. 1.0E-4)) THEN
          WRITE (6,779,ADVANCE='NO')
 779      FORMAT(' COMPLEX, ')
        ELSE
          WRITE (6,780,ADVANCE='NO')
 780      FORMAT(' REAL, ')
        END IF
C
C   DESIGNATE SOLUTION "FINITE" OR "INFINITE"
C
        IF (SUM(ABS(ROOTS(1:2,NP1,M))) .LT. 1.0E-6) THEN
          WRITE (6,781)
 781      FORMAT('INFINITE SOLUTION')
        ELSE
          WRITE (6,782)
 782      FORMAT('FINITE SOLUTION')
        END IF
C
        WRITE (6,320) LAMBDA(M),(J,(ROOTS(L,J,M),L=1,2),J=1,N)
 320    FORMAT('  LAMBDA =',ES22.14,/,(' X(',I2,') =',2ES22.14))
        WRITE (6,330) NP1,ROOTS(1:2,NP1,M)
 330    FORMAT(/,' X(',I2,') =',2ES22.14,//)
      END DO
      WRITE (6,400) ITOTIT
 400  FORMAT(' TOTAL NFE OVER ALL PATHS = ',I10)
      STOP
 1000 FORMAT(I5)
 2000 FORMAT(ES22.14)
      END PROGRAM TESTP
!
! HOMOTOPY subroutines for the polynomial system driver POLSYS1H.
! These subroutines should be used verbatim with POLSYS1H for solving
! polynomial systems of equations.  The polynomial coefficients, defined
! as input to POLSYS1H, are accessed by the routines here via the global
! arrays in HOMPACK90_GLOBAL.
!
C ###################################################################
C ONLY THE SUBROUTINES RHO AND RHOJAC ARE USED BY THE POLYNOMIAL
C SYSTEM DRIVER POLSYS1H.  ALL THE OTHER ROUTINES HERE ARE PROVIDED
C SIMPLY AS TEMPLATES.
C ###################################################################
!
      SUBROUTINE F(X,V)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(IN):: X(:)
      REAL (KIND=R8), INTENT(OUT):: V(:)
C
C EVALUATE  F(X)  AND RETURN IN THE VECTOR  V .
C
      V(1)=X(1) ! INTENT(OUT) VARIABLE MUST BE DEFINED.
      RETURN
      END SUBROUTINE F

      SUBROUTINE FJAC(X,V,K)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(IN):: X(:)
      REAL (KIND=R8), INTENT(OUT):: V(:)
      INTEGER, INTENT(IN):: K
C
C RETURN IN  V  THE KTH COLUMN OF THE JACOBIAN MATRIX OF
C F(X) EVALUATED AT  X .
C
      V(1)=X(1) ! INTENT(OUT) VARIABLE MUST BE DEFINED.
      RETURN
      END SUBROUTINE FJAC

      SUBROUTINE RHO(A,LAMBDA,X,V)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR
      REAL (KIND=R8), INTENT(IN):: A(:),X(:)
      REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
      REAL (KIND=R8), INTENT(OUT):: V(:)
C
C EVALUATE  RHO(A,LAMBDA,X)  AND RETURN IN THE VECTOR  V .
C
C THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
C  POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER IS
C CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
C REPLACEMENT CODE HERE.
      INTERFACE
        SUBROUTINE HFUNP(N,A,LAMBDA,X)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: N
        REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
        END SUBROUTINE HFUNP
      END INTERFACE
      INTEGER:: J,NPOL
C FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
      IF (LAMBDA .LT. 0.0) LAMBDA=0.0
      NPOL=IPAR(1)
      CALL HFUNP(NPOL,A,LAMBDA,X)
      DO J=1,2*NPOL
        V(J)=PAR(IPAR(3 + (4-1)) + (J-1))
      END DO
C
      RETURN
      END SUBROUTINE RHO

      SUBROUTINE RHOA(A,LAMBDA,X)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(OUT):: A(:)
      REAL (KIND=R8), INTENT(IN):: LAMBDA,X(:)
C
C CALCULATE AND RETURN IN  A  THE VECTOR Z SUCH THAT
C  RHO(Z,LAMBDA,X) = 0 .
C
      A(1)=LAMBDA ! INTENT(OUT) VARIABLE MUST BE DEFINED.
      RETURN
      END SUBROUTINE RHOA

      SUBROUTINE RHOJAC(A,LAMBDA,X,V,K)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR
      REAL (KIND=R8), INTENT(IN):: A(:),X(:)
      REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
      REAL (KIND=R8), INTENT(OUT):: V(:)
      INTEGER, INTENT(IN):: K
C
C RETURN IN THE VECTOR  V  THE KTH COLUMN OF THE JACOBIAN
C MATRIX [D RHO/D LAMBDA, D RHO/DX] EVALUATED AT THE POINT
C (A, LAMBDA, X).
C
C THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
C  POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER IS
C CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
C REPLACEMENT CODE HERE.
      INTERFACE
        SUBROUTINE HFUNP(N,A,LAMBDA,X)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: N
        REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
        END SUBROUTINE HFUNP
      END INTERFACE
      INTEGER:: J,NPOL,N2
      NPOL=IPAR(1)
      N2=2*NPOL
      IF (K .EQ. 1) THEN
C FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
        IF (LAMBDA .LT. 0.0) LAMBDA=0.0
        CALL HFUNP(NPOL,A,LAMBDA,X)
        DO J=1,N2
          V(J)=PAR(IPAR(3 + (6-1)) + (J-1))
        END DO
        RETURN
      ELSE
        DO J=1,N2
          V(J)=PAR(IPAR(3 + (5-1)) + (J-1) + N2*(K-2))
        END DO
      ENDIF
C
      RETURN
      END SUBROUTINE RHOJAC

      SUBROUTINE FJACS(X)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
      REAL (KIND=R8), INTENT(IN):: X(:)
C
C If MODE = 1,
C evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
C the result in packed skyline storage format in QRSPARSE.  LENQR is the
C length of QRSPARSE, and ROWPOS contains the indices of the diagonal
C elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
C ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPOS
C is not used by this storage format.
C
C If MODE = 2,
C evaluate the N x N Jacobian matrix of F(X) at X, and return the result
C in sparse row storage format in QRSPARSE.  LENQR is the length of
C QRSPARSE, ROWPOS contains the indices of where each row begins within
C QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
C the corresponding elements in QRSPARSE.  Even if zero, the diagonal
C elements of the Jacobian matrix must be stored in QRSPARSE.
C
      RETURN
      END SUBROUTINE FJACS

      SUBROUTINE RHOJS(A,LAMBDA,X)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
      REAL (KIND=R8), INTENT(IN):: A(:),LAMBDA,X(:)
C
C If MODE = 1,
C evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
C the result in packed skyline storage format in QRSPARSE.  LENQR is the
C length of QRSPARSE, and ROWPOS contains the indices of the diagonal
C elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
C ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPOS
C is not used by this storage format.
C
C If MODE = 2,
C evaluate the N x N Jacobian matrix of F(X) at X, and return the result
C in sparse row storage format in QRSPARSE.  LENQR is the length of
C QRSPARSE, ROWPOS contains the indices of where each row begins within
C QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
C the corresponding elements in QRSPARSE.  Even if zero, the diagonal
C elements of the Jacobian matrix must be stored in QRSPARSE.
C
      RETURN
      END SUBROUTINE RHOJS
