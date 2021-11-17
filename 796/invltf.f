      SUBROUTINE INVLTF(TOL,VALT,FZ,SIGMA0,SSBAR,NMAX,FZINV,ERROR,
     +                  IFZEVAL,WORK,IFAIL)
C
C *********************************************************************
C
C    PURPOSE:
C    =======
C    THIS ROUTINE COMPUTES AN APPROXIMATE VALUE OF AN INVERSE LAPLACE
C         TRANSFORM, BY MEANS OF AN EXPANSION IN SINE AND COSINE
C         FOURIER SERIES.
C         THE METHOD WHICH THIS SOFTWARE IS BASED ON UTILIZES THE Q-D
C         ALGORITHM TO ACCELERATE THE CONVERGENCE OF THE SERIES.
c         THE SUMMATION OF SUCCESSIVE TERMS IS TRUNCATED WHEN THE
C         ESTIMATED TRUNCATION ERROR IS LESS OR EQUAL THAN AN INPUT
C         PROVIDED TOLERANCE. THE DISCRETIZATION ERROR IS MADE LESS OR
C         EQUAL THAN THIS TOLERANCE BY SUITABLY CHOOSING SOME METHOD'S
C         PARAMETERS.
C
C
C   THE CALLING SEQUENCE IS
C   =======================
C
C     CALL  INVLTF(TOL,VALT,FZ,SIGMA0,SSBAR,NMAX,FZINV,ERROR,
C    +                  IFZEVAL,WORK,IFAIL)
C
C
C
C
C
C   INPUT PARAMETERS
C   ================
C
C
C   TOL: DOUBLE PRECISION.   
C                         ON ENTRY TOL CONTAINS THE RELATIVE 
C                         REQUIRED ACCURACY.
C
C   VALT:   DOUBLE PRECISION.
C                         ON ENTRY, VALT CONTAINS A POSITIVE VALUE OF T
C                         FOR WHICH THE INVERSE LAPLACE TRANSFORM IS REQUIRED. 
C                         VALT  HAS TO  BE GREATER THAN ZERO.
C
C   FZ:  COMPLEX*16 (DOUBLE PRECISION COMPLEX).
C                         NAME OF THE FUNCTION SUBPROGRAMS FOR THE COMPLEX 
C                         VALUED LAPLACE TRANSFORM TO BE INVERTED.
C                         ITS SPECIFICATION IS:
C                         COMPLEX*16 FUNCTION FZ(Z)
C                         COMPLEX*16 Z
C
C                         Z: COMPLEX*16. ON ENTRY, Z MUST SPECIFY THE POINT AT 
C                         WHICH THE LAPLACE TRANSFORM FUNCTION VALUE IS
C                         REQUIRED.
C
C                         FZ MUST BE DECLARED AS EXTERNAL IN THE PROGRAM
C                         FROM WHICH INVLTF IS
C                         CALLED.
C
C   SIGMA0: DOUBLE PRECISION.
C                         ON ENTRY, SIGMA0 CONTAINS THE VALUE OF THE
C                         ABSCISSA OF CONVERGENCE OF 
C                         THE LAPLACE TRANSFORM FUNCTION TO
C                         BE INVERTED OR AN UPPER BOUND OF THIS.
C                         IT IS RECOMMENDED THAT A CORRECT VALUE TO SIGMA0
C                         OR A CORRECT UPPER BOUND TO SIGMA0 IS PROVIDED.
C                         IF AN INCORRECT VALUE IS USED THE ROUTINE APPEARS TO 
C                         WORK WELL BUT CONVERGES TO COMPLETELY WRONG RESULTS.
C                         THERE IS NO WAY IN WHICH THE ROUNTINE CAN DETECT THIS. 
C
C   SSBAR:  DOUBLE PRECISION. 
C                         ON ENTRY, IT SPECIFIES THE VALUE OF THE
C                         PARAMETER SS (GREATER THAN 2) TO BE USED IN 
C                         CALCULATING THE PARAMETER D*T.
C                         TO OBTAIN DEFAULT OPTION (SS=4.1d0) ONE
C                         MAY SET   SSBAR = 0.
C
C   NMAX: INTEGER. 
C                         ON ENTRY, NMAX SPECIFIES THE MAXIMUM 
C                         NUMBER OF FZ EVALUATIONS ALLOWED.
C
C
C   OUTPUT PARAMETERS
C   =================
C
C   FZINV: DOUBLE PRECISION.
C                         ON EXIT, FZINV CONTAINS  THE APPROXIMATION
C                         OF THE INVERSE LAPLACE TRANSFORM  AT THE POINT VALT.
C
C   ERROR: DOUBLE PRECISION ARRAY OF DIMENSION 3.
C                        ON EXIT, ERROR(1) CONTAINS AN ESTIMATE OF THE RELATIVE ERROR
C                        WHICH SHOULD BE AN  UPPER BOUND FOR
C
C                        ABS(TRUEVALUE AT VALT - FZINV)/ABS(TRUEVALUE AT VALT)
C
C                        ON EXIT, ERROR(2) CONTAINS AN ESTIMATE OF THE ABSOLUTE ERROR
C                        WHICH SHOULD BE AN UPPER BOUND FOR
C
C                        ABS(TRUEVALUE AT VALT - FZINV)
C
C                        ON EXIT, ERROR(3) CONTAINS AN ESTIMATE OF THE TRUNCATION ERROR 
C                        MADE IN  CALCULATING FZINV.
C
C   IFZEVAL: INTEGER.
C                        ON EXIT, IFZEVAL CONTAINS THE
C                        NUMBER OF EVALUATIONS OF FZ USED TO OBTAIN FZINV.
C
C   IFAIL: INTEGER. 
C                        ON EXIT, A DIAGNOSTIC.
C                        = 0   ALL INPUTS WITHIN LIMITS. 
C                              ALGORITHM APPERENTLY WORKED
C                              PROPERLY.
C                        =  1  TOL IS EQUAL OR GREATER THAN 1.
C                        =  2  VALT NOT POSITIVE. 
C                        = -1  ACCURACY NOT REACHED AFTER NMAX FUNCTION EVALUATIONS.
C                        = -2  ALGORITHM DID NOT APPEAR TO BE CONVERGIING, 
C                              POSSIBLY DUE TO AN 
C                              UNSUITABLE VALUE OF PARAMETER SSBAR.
C
C    WORK : DOUBLE COMPLEX  ARRAY OF DIMENSION(2,0:2*NMAX). WORKSPACE AREA.
C
C
C   SUBROUTINES OR FUNCTIONS NEEDED
C
C     FZ      : USER PROVIDED FUNCTION
C     QDACC   : Q-D ALGORITHM
C     BACKCF  : COMPUTATION OF THE CONTINUED FRACTION
C     D1MACH  : PROVIDES MACHINE EPSILON
C
C ***************************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION SIGMA0,SS,TOL,SSBAR,VALT,FZINV
      INTEGER IFAIL,NMAX,IFZEVAL
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX  WORK(2,0:NMAX)
      DOUBLE PRECISION ERROR(3)
C     ..
C     .. Function Arguments ..
      DOUBLE COMPLEX FZ
      EXTERNAL FZ
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,GAMMA,ALPHA,DT,ERR1,ERR2,FA3,PART2,
     +                 PART0,PART1,RELP,T,TOL1,TOLL,D1MACH
      INTEGER MA,MP
      LOGICAL FIRST
C     ..
C     .. External Subroutines ..
      EXTERNAL QDACC,BACKCF,D1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DATAN,DLOG
C     ..
C     .. Common blocks ..
      COMMON /PARAM/FIRST,MA,MP
      COMMON /PARAM1/GAMMA,PI
C     ..
C
      ERROR(1) = 0.0D0
      ERROR(2) = 0.0D0
      ERROR(3) = 0.0D0
C
      SS= SSBAR + 4.1D0
C ****************************************************************************
C
C  INPUT DATA CHECK
C
C ****************************************************************************
      IFAIL = 0

          IF (TOL .GE. 1) THEN
              IFAIL = 1
              RETURN

          END IF
          IF (VALT.LT.0) THEN
              IFAIL = 2
              RETURN

          END IF

      
C
C
C ************************************************************************
C  SETTING THE PARAMETERS : ALPHA, TOL1, RELP
C  THE PRODUCT D*T IS COMPUTED IN SUCH A WAY THAT
C  THE ERROR IS LESS THAN TOL
C ************************************************************************
C

      RELP = D1MACH(4)
C     RELP = X02AJF()
C     RELP  =0.222044d-15 
C     FOR AN IBM 6000 RISC STATION MOD. 32H
      IF (TOL.LT.RELP) TOL = RELP
      TOL1 = TOL/10.D0
      ALPHA = 1.4D0
c
      DT = -DLOG(TOL1*2*DATAN(1.D0))/ (SS-DLOG(10.D0))
c
      DT = DT*ALPHA
          GAMMA = DT/VALT + SIGMA0
          T = SS/2.D0*VALT
          IFZEVAL = 0
          PI = 4.D0*DATAN(1.D0)/T
      FIRST = .TRUE.
      MA = 2
      MP = MA - 1
C ***********************************************************************
C
C    CALL TO THE Q-D ALGORITHM
C
C ***********************************************************************
10     CONTINUE

       IF (IFAIL .NE. 0) RETURN 
       CALL QDACC(FZ,IFZEVAL,IFAIL,WORK,NMAX,PART0,PART1,PART2)
       IF (IFAIL .NE. 0) RETURN 
       CALL BACKCF(VALT,WORK,NMAX,PART2)
C
C
C
C**********************************************************************
C
C     FZINV IS THE APPROXIMATED VALUE OF THE INVERSE LAPLACE TRANSFORM
C
C *********************************************************************
      FA3 = DEXP(DT)/T
      FZINV = DEXP(SIGMA0*VALT)*FA3*PART2
C
C
C
C *********************************************************************
C
C    COMPUTATION OF THE TRUNCATION ERROR
C
C ERR1 IS THE DIFFERENCE BETWEEN THE LAST APPROXIMATION AND THE PREVIOUS
C ONE,
C WHILE ERR2 IS THE DIFFERENCE BETWEEN THE LAST APPROXIMATION
C AND THE LAST-THIRD
C
C ***************************************************************************
C
      ERR1 = FA3*ABS(PART2-PART1)
      ERR2 = FA3*ABS(PART2-PART0)
      ERROR(2) = TOL1*ABS(FA3*PART2)
      TOLL = TOL1*ABS(FA3*PART2) + RELP*ABS(FA3*PART2)
      IF (TOLL.EQ.0.D0) TOLL = RELP
C
C ****************************************************************************
C
C               STOPPING CRITERION
C
C     THE ALGORITHM STOPS WHEN BOTH ERR1 AND ERR2 ARE LESS THAN RELP
C                         OR
C     THE NUMBER OF FUNCTION EVALUATIONS, OF COURSE, IS GREATER OR EQUAL TO NMAX
C
C *****************************************************************************
C

      IF (IFZEVAL.LT.NMAX) THEN
          IF (ERR1.NE.0.D0 .AND. ERR2.NE.0.D0) THEN
              IF (ERR1.GT.RELP .AND. ERR2.GT.RELP) THEN
                  IF (ERR1.GT.TOLL .OR. ERR2.GT.TOLL) GO TO 10
              END IF

          ELSE IF (ERR1.EQ.0.D0 .AND. ERR2.NE.0.D0) THEN
              IF (ERR2.GT.RELP .OR. ERR2.GT.TOLL) THEN
                  GO TO 10

              ELSE IF (ERR2.EQ.0.D0 .AND. ERR1.NE.0.D0) THEN
                  IF (ERR1.GT.RELP .OR. ERR1.GT.TOLL) GO TO 10
              END IF

          END IF

      ELSE
          IF (ERR1.NE.0.D0 .AND. ERR2.NE.0.D0) THEN
              IF (ERR1.GT.RELP .AND. ERR2.GT.RELP) THEN
                  IF (ERR1.GT.TOLL .OR. ERR2.GT.TOLL) THEN
                      IFAIL = - 1
                  END IF

              END IF

          ELSE IF (ERR1.EQ.0.D0 .AND. ERR2.NE.0.D0) THEN
                IF (ERR2.GT.RELP .OR. ERR2.GT.TOLL) THEN
                       IFAIL = -1

              ELSE IF (ERR2.EQ.0.D0 .AND. ERR1.NE.0.D0) THEN
                  IF (ERR1.GT.RELP .OR. ERR1.GT.TOLL) THEN
                      IFAIL = -1
                  END IF

              END IF

          END IF

      END IF
C ****************************************************************************
C
C      ESTIMATION OF THE ABSOLUTE ERROR, TRUNCATION ERROR AND THE SUM OF THEM
C
C ****************************************************************************

          ERROR(2) = ERROR(2) + (ERR1+ERR2)
          IF (PART2.NE.0.D0) THEN
              ERROR(3) = (ERR1+ERR2)/ABS(FA3*PART2)

          ELSE
              ERROR(3) = ERR1 + ERR2
          END IF

          ERROR(1) = TOL1 + ERROR(3)
      RETURN

      END
C
C ****************************************************************
C
      SUBROUTINE QDACC(FZ,IFV,IFAIL,WORK,NMAX,PART0,PART1,PART2)
C
C ****************************************************************
C    PURPOSE:
C    =======
C    THIS SUBROUTINE IMPLEMENTS A COLUMN-VERSION OF THE Q-D ALGORITHM.
C    THE QD SCHEME IS BUILT UP PROCEEDING FROM LEFT TO RIGHT ALONG THE
C    DIAGONAL DIRECTION. THE ELEMENTS OF EACH DIAGONALS ARE COMPUTED 
C    PROCEEDING BOTTOM-UP IN THE TABLE. THE LAST ELEMENT COMPUTED IN
C    EACH DIAGONAL IS THE COEFFICIENT OF THE CONTINUED FRACTION.
C                                                                       
C    DESCRIPTION OF PARAMETERS:
C    =========================
C
C    INPUT PARAMETERS:
C    
C   FZ:  COMPLEX*16 (DOUBLE PRECISION COMPLEX) FUNCTION, SUPPLIED BY
C        THE USER.
C        FZ MUST RETURN THE VALUE OF THE LAPLACE TRANSFORM FUNCTION TO
C        BE INVERTED, AT A GIVEN POINT. ITS SPECIFICATION IS:
C                  COMPLEX*16 FUNCTION FZ(Z)
C                  COMPLEX*16 Z
C
C                  Z: COMPLEX*16. ON ENTRY, Z MUST SPECIFY THE POINT AT
C                     WHICH THE LAPLACE TRANSFORM FUNCTION VALUE IS
C                     REQUIRED.
C
C        FZ MUST BE DECLARED AS EXTERNAL IN THE PROGRAM FROM WHICH INVLTF IS
C        CALLED.
C 
C   IFV  : INTEGER. ON EXIT, IFV CONTAINS THE
C          NUMBER OF EVALUATIONS OF FZ USED TO OBTAIN PART2.
C 
C   IFAIL: INTEGER. ON EXIT, IFAIL CONTAINS POSSIBLE ERRORS DETECTED
C   
C   WORK : DOUBLE COMPLEX ARRAY OF DIMENSION(2,0:2*NMAX). WORKSPACE AREAS.
C  
C   PART0: DOUBLE PRECISION. ON EXIT PART0 CONTAINS THE APPROXIMATION
C          OF THE LAST-THIRD ACCELERATED TERM.
C 
C   PART1: DOUBLE PRECISION.  ON EXIT PART1 CONTAINS THE APPROXIMATION
C          OF THE LAST-SECOND ACCELERATED TERM.
C
C   PART2: DOUBLE PRECISION. ON EXIT PART2 CONTAINS THE APPROXIMATION 
C          OF THE LAST ACCELERATED.
C ****************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION GAMMA,PART2,PI
      INTEGER IFAIL,IFV,MA,NMAX
      LOGICAL FIRST
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX WORK(2,0:2*NMAX)
C     ..
C     .. Function Arguments ..
      DOUBLE COMPLEX FZ
      EXTERNAL FZ
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX AUX1,AUX2,AUX3,H,MX
      DOUBLE PRECISION PART0,PART1,TOLL
      INTEGER IBEGIN,I,J,K,MM,MP,UL
      LOGICAL RIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DCMPLX,DEXP,DREAL,CDEXP,CDSQRT
C     ..
C     .. Common blocks ..
      COMMON /PARAM/FIRST,MA,MP
      COMMON /PARAM1/GAMMA,PI
C **********************************************************
C   INITIALIZATION STEP
C **********************************************************
      IBEGIN = 2*MP
      MM = 2*MA
      IF (FIRST) THEN
          H = FZ(DCMPLX(GAMMA,0.D0))/ (2.D0,0.D0)
          IFV = IFV + 1
          MX = H
          WORK(2,0) = H
C **************************************************************
C   INITIALIZE THE WORKSPACE AREA: THE ARRAYS Q and D OF 
C   THE COLUMN VERSION OF THE QD ALGORITHM
C   THE VECTOR Q HAS BEEN STORED IN THE FIRST ROW OF THE ARRAY WORK
C   THE VECTOR D HAS BEEN STORED IN THE SECOND ROW OF THE ARRAY WORK
C ***************************************************************
          DO 20 K = 0,1
              WORK(1,K) = MX
              MX = FZ(DCMPLX(GAMMA,PI* (K+1)))
              IFV = IFV + 1
              WORK(1,K) = MX/ WORK(1,K)
              WORK(2,K+1) =  WORK(1,K)
   20     CONTINUE
          AUX1 =  WORK(1,1)
          WORK(1,1) =  WORK(1,1) -  WORK(1,0)
          WORK(2,2) = - WORK(1,1)
          WORK(1,0) = AUX1
          WORK(2,1) = -WORK(2,1)
          WORK(1,2) = 0.D0
      END IF

      UL = MM - 1
C     DO 15 K = 3, UL
C       WORK(1,K) = 0.0D0
C15    CONTINUE
C
C *****************************************************************
C COMPUTATION OF THE COEFFICIENTS NEEDED IN THE QD ALGORITHM
C *****************************************************************
          IFAIL = 0
          I= IBEGIN
25        CONTINUE
          RIV = .TRUE.
          AUX1 = WORK(1,0)
          WORK(1,0) = MX
          MX = FZ(DCMPLX(GAMMA,PI* (I+1)))
          IFV = IFV + 1
          WORK(1,0) = MX/WORK(1,0)
          AUX2 = WORK(1,1)
          WORK(1,1) = WORK(1,0) - AUX1
          J=2
30        CONTINUE
              IF (J .EQ. I)  WORK(1,J) = 0.D0
              AUX3 = WORK(1,J)
              IF (RIV) THEN
                  IF (AUX2.NE. (0.D0,0.D0)) THEN
                      WORK(1,J) = AUX1*WORK(1,J-1)/AUX2
                  ELSE
                      IFAIL = -2
                  END IF

              ELSE
                  WORK(1,J) = WORK(1,J-1) - AUX2 + AUX1
              END IF

              AUX1 = AUX2
              AUX2 = AUX3
              RIV = .NOT. RIV
          J=J+1
          IF(J .LE.I .AND. IFAIL .EQ. 0) GOTO 30
          WORK(2,I+1) = -WORK(1,I)
          I = I+1
          IF (I .LE. UL .AND. IFAIL .EQ. 0 ) GOTO 25
      IF (.NOT.FIRST) THEN
         IBEGIN = IBEGIN + 1
      END IF
      IF (FIRST) THEN
          PART0 = 0.D0
          PART1 = 0.D0

      ELSE
          PART0 = PART1
          PART1 = PART2
      END IF
C ***************************************************************************
C ************* END OF THE QD ALGORITHM  *************************************
C ***************************************************************************
      RETURN 
      END

C *****************************************************************************

      SUBROUTINE BACKCF(TV,WORK,NMAX,PART2)

C *****************************************************************************
C  COMPUTATION OF THE  CONTINUED FRACTION BY THE  BACKWARD FORMULA
C *****************************************************************************
      INTEGER MA,MM,MP,NMAX,I
      DOUBLE COMPLEX  WORK(2,0:2*NMAX),CC,Z,H2M,R2M

      DOUBLE PRECISION PART2,GAMMA,PI,TV
      LOGICAL FIRST
      COMMON /PARAM/FIRST,MA,MP
      COMMON /PARAM1/ GAMMA,PI

      
      MM=2*MA
      Z = CDEXP(DCMPLX(0.D0,PI*TV))
      H2M = ((WORK(2,MM-1)-WORK(2,MM))*Z+ (1.D0,0.D0))/ (2.D0,0.D0)
      R2M = CDSQRT((WORK(2,MM)*Z/ (H2M*H2M))+ (1.D0,0.D0))
      R2M = (1.D0,0.D0) - R2M
      R2M = -H2M*R2M
      CC = (R2M*Z)
      DO 50 I = MM - 1,1,-1
          CC = (WORK(2,I)*Z)/ ((1.D0,0.D0)+CC)
   50 CONTINUE
      CC = WORK(2,0)/ ((1.D0,0.D0)+CC)
      PART2 = DREAL(CC)
      MP = MA
      MA = MA + 1
      FIRST = .FALSE.
      RETURN
      END

