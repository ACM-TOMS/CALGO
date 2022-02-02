! This program tests many cases that arise with the module
! INTERVAL_ARITHMETIC. Each of the possible combinations of mixed-mode
! operations is tested.
!
! This version assumes that small integers are stored exactly into
! floating point numbers.
!
PROGRAM TEST_INTERVAL_SYSTEM

   USE INTERVAL_ARITHMETIC
   IMPLICIT NONE

      INTEGER :: OUTPUT_UNIT = 6
      CHARACTER(LEN=12) :: OUTPUT_FILE_NAME = 'INTARITH.OUT'

      CHARACTER(LEN=10) CURRENT_DATE, CURRENT_TIME

      TYPE(INTERVAL) :: X, Y, Z, W, XX
      DOUBLE PRECISION R, RR
      INTEGER I
      LOGICAL ALL_TESTS_OK
      DOUBLE PRECISION TOL

      DOUBLE PRECISION TINY, TEST
      COMMON /MACH2/ TINY, TEST

      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN

      CALL SIMINI
      ISEVER = 4
      IERPUN = OUTPUT_UNIT
      IPRTCL = 0

      CALL DATE_AND_TIME(DATE=CURRENT_DATE)
      CALL DATE_AND_TIME(TIME=CURRENT_TIME)

      OPEN(UNIT=OUTPUT_UNIT, FILE=OUTPUT_FILE_NAME)

      WRITE(OUTPUT_UNIT,'(4(1x,a))') &
                           'Output from TEST_INTERVAL_SYSTEM on ', &
                           CURRENT_DATE(5:6)// &
                           "/"//CURRENT_DATE(7:8)// &
                           "/"//CURRENT_DATE(1:4), ' at ', &
                           CURRENT_TIME(1:2)// &
                           ":"//CURRENT_TIME(3:4)// &
                           ":"//CURRENT_TIME(5:6)//"."
      WRITE(OUTPUT_UNIT,'(I1)')

      X = INTERVAL(1,3)
      Y = INTERVAL(-4,-2)
      Z = INTERVAL(-3,2)
      W = INTERVAL(0,.5D0)
      R = 2D0
      I = 3
      ALL_TESTS_OK = .TRUE.
      TOL = 100D0*EPSILON(1D0)

      WRITE(OUTPUT_UNIT,'(1X,A)') 'X, Y, Z, W:'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  X
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  Y
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  Z
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  W
      WRITE(OUTPUT_UNIT,'(1X,A)') 'R:'
      WRITE(OUTPUT_UNIT, '(1X,ES28.20E2)')  R
      WRITE(OUTPUT_UNIT,'(1X,A)') 'I:'
      WRITE(OUTPUT_UNIT, '(1X,I10)')  I
      WRITE(OUTPUT_UNIT,'(1X)')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing the four elementary operations and exponentiation ---'
      CALL CHECK_RESULT(X+Y,INTERVAL(-3,-1),'X+Y')
      CALL CHECK_RESULT(X-Y,INTERVAL(3,7),'X-Y')
      CALL CHECK_RESULT(X*Y,INTERVAL(-12,-2),'X*Y')
      CALL CHECK_RESULT(X/Y,INTERVAL(-1.5D0,-.25D0),'X/Y')

      CALL CHECK_RESULT(X+R,INTERVAL(3,5),'X+R')
      CALL CHECK_RESULT(R+X,INTERVAL(3,5),'R+X')
      CALL CHECK_RESULT(X-R,INTERVAL(-1,1),'X-R')
      CALL CHECK_RESULT(R-X,INTERVAL(-1,1),'R-X')
      CALL CHECK_RESULT(X*R,INTERVAL(2,6),'X*R')
      CALL CHECK_RESULT(R*X,INTERVAL(2,6),'R*X')
      CALL CHECK_RESULT(X/R,INTERVAL(0.5D0,1.5D0),'X/R')
      CALL CHECK_RESULT(R/X,&
       INTERVAL(0.666666666666666666666666666667D0,2),'R/X')

      CALL CHECK_RESULT(X+I,INTERVAL(4,6),'X+I')
      CALL CHECK_RESULT(I+X,INTERVAL(4,6),'I+X')
      CALL CHECK_RESULT(X-I,INTERVAL(-2,0),'X-I')
      CALL CHECK_RESULT(I-X,INTERVAL(0,2),'I-X')
      CALL CHECK_RESULT(X*I,INTERVAL(3,9),'X*I')
      CALL CHECK_RESULT(I*X,INTERVAL(3,9),'I*X')
      CALL CHECK_RESULT(X/I,&
       INTERVAL(0.333333333333333333333333333333D0,1),'X/I')
      CALL CHECK_RESULT(I/X,INTERVAL(1,3),'I/X')

      CALL CHECK_RESULT(X**Y,&
       INTERVAL(0.123456790123456790123456790123D0,1),'X**Y')
      CALL CHECK_RESULT(X**I,INTERVAL(1,27),'X**I')
      CALL CHECK_RESULT(I**X,INTERVAL(3,27),'I**X')
      CALL CHECK_RESULT(X**R,INTERVAL(1,9),'X**R')
      CALL CHECK_RESULT(R**X,INTERVAL(2,8),'R**X')


      WRITE(OUTPUT_UNIT,'(1X,A)') 'Testing X <-- X op X ---'
      XX = X;  XX = XX + XX
      CALL CHECK_RESULT(XX,INTERVAL(2,6),'X=X+X')
      XX = X;  XX = XX - XX
      CALL CHECK_RESULT(XX,INTERVAL(-2,2),'X=X-X')
      XX = X;  XX = XX * XX
      CALL CHECK_RESULT(XX,INTERVAL(1,9),'X=X*X')
      XX = X;  XX = XX / XX
      CALL CHECK_RESULT(XX,&
       INTERVAL(0.333333333333333333333333333333D0,3D0),'X=X/X')

      WRITE(OUTPUT_UNIT,'(1x,a)') &
          'Testing all cases of multiplication ---'
      CALL CHECK_RESULT(INTERVAL(1,3)*INTERVAL(2,4),&
       INTERVAL(2,12), 'Multiplication case 1--- [1,3]*[2,4]')
      CALL CHECK_RESULT(INTERVAL(-1,3)*INTERVAL(2,4),&
       INTERVAL(-4,12), 'Multiplication case 2--- [-1,3]*[2,4]')
      CALL CHECK_RESULT(INTERVAL(-2,-1)*INTERVAL(3,4),&
       INTERVAL(-8,-3), 'Multiplication case 3--- [-2,-1]*[3,4]')
      CALL CHECK_RESULT(INTERVAL(1,2)*INTERVAL(-1,1),&
       INTERVAL(-2,2), 'Multiplication case 4--- [1,2]*[-1,1]')
      CALL CHECK_RESULT(INTERVAL(-2,2)*INTERVAL(-1,1),&
       INTERVAL(-2,2), 'Multiplication case 5--- [-2,2]*[-1,1]')
      CALL CHECK_RESULT(INTERVAL(-2,-1)*INTERVAL(1,2),&
       INTERVAL(-4,-1), 'Multiplication case 6--- [-2,-1]*[1,2]')
      CALL CHECK_RESULT(INTERVAL(-1,1)*INTERVAL(-2,-1),&
       INTERVAL(-2,2), 'Multiplication case 7--- [-1,1]*[-2,-1]')
      CALL CHECK_RESULT(INTERVAL(-2,-1)*INTERVAL(-4,-3),&
       INTERVAL(3,8), 'Multiplication case 8--- [-2,-1]*[-4,-3]')
      CALL CHECK_RESULT(INTERVAL(-3,1)*INTERVAL(-2,3),&
       INTERVAL(-9,6), 'Multiplication case 9--- [-3,1]*[-2,3]')

      CALL CHECK_RESULT(INTERVAL(-1,5)*INTERVAL(-1,1),&
       INTERVAL(-5,5), '[-1,5]*[-1,1]')

      WRITE(OUTPUT_UNIT,'(1X,A)') 'Additional testing of X<--X*X ---'
      XX = INTERVAL(-3,-1); XX = XX*XX
      CALL CHECK_RESULT(XX,INTERVAL(1,9),&
       'Case 8 of multiplication, U=U*U with U=[-3,-1]')
      XX = INTERVAL(-3,1); XX = XX*XX
      CALL CHECK_RESULT(XX,INTERVAL(-3,9),&
       'Case 9 of multiplication, U=U*U with U=[-3,1]')

      WRITE(OUTPUT_UNIT,'(1X,A)') 'Testing double*interval ---'
      CALL CHECK_RESULT(2.0D0*INTERVAL(0,3),INTERVAL(0,6),&
       '2*[0,3]')
      CALL CHECK_RESULT(-2.0D0*INTERVAL(0,3),INTERVAL(-6,0), &
       '-2*[0,3]')
      CALL CHECK_RESULT(2.0D0*INTERVAL(-3,0),INTERVAL(-6,0), &
       '2*[-3,0]')
      CALL CHECK_RESULT(-2.0D0*INTERVAL(-3,0),INTERVAL(0,6), &
       '-2*[-3,0]')
      CALL CHECK_RESULT(2.0D0*INTERVAL(1,3),INTERVAL(2,6), &
       '2*[1,3]')
      CALL CHECK_RESULT(-2.0D0*INTERVAL(1,3),INTERVAL(-6,-2), &
       '-2*[1,3]')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
         'Testing additional cases of division ---'
        CALL CHECK_RESULT(Y/X, &
       INTERVAL(-4,-0.666666666666666666666666666667D0),'Y/X')
      XX = Y; XX = XX/XX
      CALL CHECK_RESULT(XX,INTERVAL(0.5D0,2D0),'Y<--Y/Y')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing MAG, WID, MID, MIG, ABS ---'
      CALL CHECK_REAL_RESULT(MAG(X),3D0,'MAG(X)')
      CALL CHECK_REAL_RESULT(MAG(Y),4D0,'MAG(Y)')
      CALL CHECK_REAL_RESULT(MAG(Z),3D0,'MAG(Z)')
      CALL CHECK_REAL_RESULT(WID(X),2D0,'WID(X)')
      CALL CHECK_REAL_RESULT(MID(X),2D0,'MID(X)')
      CALL CHECK_REAL_RESULT(MIG(X),1D0,'MIG(X)')
      CALL CHECK_REAL_RESULT(MIG(Y),2D0,'MIG(Y)')
      CALL CHECK_REAL_RESULT(MIG(Z),0D0,'MIG(Z)')
      CALL CHECK_RESULT(ABS(X),INTERVAL(1,3),'ABS(X)')
      CALL CHECK_RESULT(ABS(Y),INTERVAL(2,4),'ABS(Y)')
      CALL CHECK_RESULT(ABS(Z),INTERVAL(0,3),'ABS(Z)')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing MAX, MIN ---'
      CALL CHECK_RESULT(MAX(X,Z),INTERVAL(1,3),'MAX(X,Z)')
      CALL CHECK_RESULT(MAX(X,R),INTERVAL(2,3),'MAX(X,R)')

      CALL CHECK_RESULT(MAX(R,X),INTERVAL(2,3),'MAX(R,X)')
      CALL CHECK_RESULT(MAX(X,I),INTERVAL(3,3),'MAX(X,I)')
      CALL CHECK_RESULT(MAX(I,X),INTERVAL(3,3),'MAX(I,X)')

      CALL CHECK_RESULT(MIN(X,Z),INTERVAL(-3,2),'MIN(X,Z)')
      CALL CHECK_RESULT(MIN(X,R),INTERVAL(1,2),'MIN(X,R)')

      CALL CHECK_RESULT(MIN(R,X),INTERVAL(1,2),'MIN(R,X)')
      CALL CHECK_RESULT(MIN(X,I),INTERVAL(1,3),'MIN(X,I)')
      CALL CHECK_RESULT(MIN(I,X),INTERVAL(1,3),'MIN(I,X)')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing ACOS, ACOT, ASIN, ATAN ---'
      CALL CHECK_RESULT(ACOS(W), &
       INTERVAL(1.04719755119659774615421446109D0, &
                1.57079632679489661923132169164D0),'ACOS(W)')
      CALL CHECK_RESULT(ACOT(W), &
       INTERVAL(1.10714871779409050301706546018D0, &
                1.57079632679489661923132169164D0),'ACOT(W)')
      CALL CHECK_RESULT(ASIN(W), &
       INTERVAL(0, &
                0.523598775598298873077107230547D0),'ASIN(W)')
      CALL CHECK_RESULT(ATAN(W), &
       INTERVAL(0, &
                0.463647609000806116214256231461D0),'ATAN(W)')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing COS, COT, EXP, LOG, SIN, SQRT, TAN, SINH ---'
      CALL CHECK_RESULT(COS(ACOS(W)),W,'COS(ACOS(W))')
      CALL CHECK_RESULT(COT(ACOT(W)),W,'COT(ACOT(W))')
      CALL CHECK_RESULT(EXP(X), &
       INTERVAL(2.71828182845904523536028747135D0, &
                20.0855369231876677409285296546D0),'EXP(X)')
      CALL CHECK_RESULT(LOG(EXP(X)),X,'LOG(EXP(X))')
      CALL CHECK_RESULT(SIN(ASIN(W)),W,'SIN(ASIN(W))')
      CALL CHECK_RESULT(SQRT(X**2),X,'SQRT(X**2)')
      CALL CHECK_RESULT(TAN(ATAN(W)),W,'TAN(ATAN(W))')
      CALL CHECK_RESULT(SINH(X), &
       INTERVAL(1.17520119364380145688238185060D0, &
                10.0178749274099018989745936195D0),'SINH(X)')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing .IS. and .CH. ---'
      CALL CHECK_RESULT(X.IS.Z,INTERVAL(1,2),'X.IS.Z')
      CALL CHECK_RESULT(X.CH.Y,INTERVAL(-4,3),'X.CH.Y')
      CALL CHECK_RESULT(Y.CH.R,INTERVAL(-4,2),'Y.CH.R')
      CALL CHECK_RESULT(R.CH.Y,INTERVAL(-4,2),'R.CH.Y')
      CALL CHECK_RESULT(Y.CH.I,INTERVAL(-4,3),'Y.CH.I')
      CALL CHECK_RESULT(I.CH.Y,INTERVAL(-4,3),'I.CH.Y')
      CALL CHECK_RESULT(1D0.CH.2D0,INTERVAL(1,2),'1D0.CH.2D0')
      CALL CHECK_RESULT(3.CH.4,INTERVAL(3,4),'3.CH.4')
      CALL CHECK_RESULT(1.CH.2D0,INTERVAL(1,2),'1.CH.2D0')
      CALL CHECK_RESULT(3D0.CH.4,INTERVAL(3,4),'3D0,.CH.4')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing .SB., .SP., .DJ., and .IN. ---'
      CALL CHECK_LOGICAL_RESULT(W.SB.Z,.TRUE.,'W.SB.Z')
      CALL CHECK_LOGICAL_RESULT(X.SB.Z,.FALSE.,'X.SB.Z')
      CALL CHECK_LOGICAL_RESULT(Z.SP.W,.TRUE.,'Z.SP.W')
      CALL CHECK_LOGICAL_RESULT(X.DJ.Y,.TRUE.,'X.DJ.Y')
      CALL CHECK_LOGICAL_RESULT(X.DJ.Z,.FALSE.,'X.DJ.Z')
      CALL CHECK_LOGICAL_RESULT(R.IN.X,.TRUE.,'R.IN.X')
      CALL CHECK_LOGICAL_RESULT(I.IN.X,.TRUE.,'I.IN.X')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing .LT., .GT., .LE., .GE., .EQ., and .NE. ---'
      CALL CHECK_LOGICAL_RESULT(Y.LT.X,.TRUE.,'Y.LT.X')
      CALL CHECK_LOGICAL_RESULT(Y.LT.R,.TRUE.,'Y.LT.R')
      CALL CHECK_LOGICAL_RESULT(R.LT.Y,.FALSE.,'R.LT.Y')
      CALL CHECK_LOGICAL_RESULT(I.LT.X,.FALSE.,'I.LT.X')
      CALL CHECK_LOGICAL_RESULT(X.LT.I,.FALSE.,'X.LT.I')

      CALL CHECK_LOGICAL_RESULT(Z.GT.X,.FALSE.,'Z.GT.X')
      CALL CHECK_LOGICAL_RESULT(Y.GT.R,.FALSE.,'Y.GT.R')
      CALL CHECK_LOGICAL_RESULT(R.GT.Y,.TRUE.,'R.GT.Y')
      CALL CHECK_LOGICAL_RESULT(I.GT.Y,.TRUE.,'I.GT.Y')
      CALL CHECK_LOGICAL_RESULT(Y.GT.I,.FALSE.,'Y.GT.I')

      CALL CHECK_LOGICAL_RESULT(Y.LE.X,.TRUE.,'Y.LE.X')
      CALL CHECK_LOGICAL_RESULT(Y.LE.R,.TRUE.,'Y.LE.R')
      CALL CHECK_LOGICAL_RESULT(R.LE.Y,.FALSE.,'R.LE.Y')
      CALL CHECK_LOGICAL_RESULT(I.LE.X,.FALSE.,'I.LE.X')
      CALL CHECK_LOGICAL_RESULT(X.LE.I,.TRUE.,'X.LE.I')

      CALL CHECK_LOGICAL_RESULT(Z.GE.X,.FALSE.,'Z.GE.X')
      CALL CHECK_LOGICAL_RESULT(Y.GE.R,.FALSE.,'Y.GE.R')
      CALL CHECK_LOGICAL_RESULT(R.GE.Y,.TRUE.,'R.GE.Y')
      CALL CHECK_LOGICAL_RESULT(I.GE.Y,.TRUE.,'I.GE.Y')
      CALL CHECK_LOGICAL_RESULT(Y.GE.I,.FALSE.,'Y.GE.I')

      CALL CHECK_LOGICAL_RESULT(X.EQ.X,.TRUE.,'X.EQ.X')
      CALL CHECK_LOGICAL_RESULT(Y.EQ.X,.FALSE.,'Y.EQ.X')
      CALL CHECK_LOGICAL_RESULT(Y.EQ.R,.FALSE.,'Y.EQ.R')
      CALL CHECK_LOGICAL_RESULT(R.EQ.Y,.FALSE.,'R.EQ.Y')
      CALL CHECK_LOGICAL_RESULT(I.EQ.X,.FALSE.,'I.EQ.X')
      CALL CHECK_LOGICAL_RESULT(X.EQ.I,.FALSE.,'X.EQ.I')

      CALL CHECK_LOGICAL_RESULT(X.NE.X,.FALSE.,'X.NE.X')
      CALL CHECK_LOGICAL_RESULT(Y.NE.X,.TRUE.,'Y.NE.X')
      CALL CHECK_LOGICAL_RESULT(Y.NE.R,.TRUE.,'Y.NE.R')
      CALL CHECK_LOGICAL_RESULT(R.NE.Y,.TRUE.,'R.NE.Y')
      CALL CHECK_LOGICAL_RESULT(I.NE.X,.TRUE.,'I.NE.X')
      CALL CHECK_LOGICAL_RESULT(X.NE.I,.TRUE.,'X.NE.I')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing INF and SUP ---'
      CALL CHECK_REAL_RESULT(INF(X),1D0,'INF(X)')
      CALL CHECK_REAL_RESULT(SUP(X),3D0,'SUP(X)')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Additional testing of .CH. dealing with empty intervals ---'
      CALL CHECK_RESULT(INTERVAL(2,1).CH.INTERVAL(1,-1), &
       INTERVAL(1,-1),'[2,1].CH.[1,-1]')
      CALL CHECK_RESULT(INTERVAL(2,1).CH.INTERVAL(-1,1), &
       INTERVAL(-1,1),'[2,1].CH.[-1,1]')
      CALL CHECK_RESULT(INTERVAL(1,2).CH.INTERVAL(1,-1), &
       INTERVAL(1,2),'[1,2].CH.[1,-1]')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Testing some implicit conversions ---'
      XX = 1
      CALL CHECK_RESULT(XX,INTERVAL(1,1),'[1,1]<-- 1')
      XX = 1D0
      CALL CHECK_RESULT(XX,INTERVAL(1,1), '[1,1]<-- 1D0')

      WRITE(OUTPUT_UNIT,'(1X,A)') &
         'Additional testing of mixed mode multiplication ---'
      XX = INTERVAL(-1,1)
      RR = 0D0
      CALL CHECK_RESULT(RR*XX,INTERVAL(0,0),'0*[-1,1]')
      XX = INTERVAL(0,1)
      RR = -1D0
      CALL CHECK_RESULT(RR*XX,INTERVAL(-1,0),'-1*[0,1]')
      RR = 1D0
      CALL CHECK_RESULT(RR*XX,INTERVAL(0,1),'1*[0,1]')
      XX = INTERVAL(-1,0)
      RR = -1D0
      CALL CHECK_RESULT(RR*XX,INTERVAL(0,1),'-1*[-1,0]')
      RR = 1D0
      CALL CHECK_RESULT(RR*XX,INTERVAL(-1,0),'1*[-1,0]')
      XX = INTERVAL(2,3); RR = -1D0
      CALL CHECK_RESULT(RR*XX,INTERVAL(-3,-2),'-1*[2,3]')

      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Check the conversion function:'
      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(.3D0):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(.3D0)
      CALL CHECK_RESULT(IVL(.3D0),INTERVAL(.3D0,.3D0),'IVL(.3D0)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(.3D0,.6D0):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(.3D0,.6D0)
      CALL CHECK_RESULT(IVL(.3D0,.6D0),INTERVAL(.3D0,.6D0), &
       'IVL(.3D0,.6D0)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(1,2):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(1,2)
      CALL CHECK_RESULT(IVL(1,2),INTERVAL(1,2),'IVL(1,2)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(1):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(1)
      CALL CHECK_RESULT(IVL(1),INTERVAL(1,1),'IVL(1)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(.3D0,1):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(.3D0,1)
      CALL CHECK_RESULT(IVL(.3D0,1),INTERVAL(.3D0,1),'IVL(.3D0,1)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(1,3.1D0):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(1,3.1D0)
      CALL CHECK_RESULT(IVL(1,3.1D0),INTERVAL(1,3.1D0),&
       'IVL(1,3.1D0)')
      WRITE(OUTPUT_UNIT,'(1X,A)') 'IVL(X):'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))')  IVL(X)
      CALL CHECK_RESULT(IVL(X),X,'IVL(X)')
      WRITE(OUTPUT_UNIT,'(I1)')

      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1x,A)') &
         'Check RNDOUT near underflow and overflow thresholds:'
      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1x,A)') '[TINY,1]/[5,6]:'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))') &
             INTERVAL(TINY,1)/INTERVAL(5,6)
      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1x,A)') '[-1,-TINY]/[5,6]:'
      WRITE(OUTPUT_UNIT, '(2(1X,ES28.20E2))') &
             INTERVAL(-1,-TINY)/INTERVAL(5,6)

      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Precipitate a division error by trying [1,2]/[-1,1]:'
      XX = INTERVAL(1,2)/INTERVAL(-1,1)

      WRITE(OUTPUT_UNIT,'(I1)')
      WRITE(OUTPUT_UNIT,'(1X,A)') &
       'Precipitate an error by trying [1,2].IS.[3,4]:'
      XX = INTERVAL(1,2).IS.INTERVAL(3,4)

      WRITE(OUTPUT_UNIT,'(I1)')

      IF(ALL_TESTS_OK) THEN
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         'All tests completed satisfactorily.  Module'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         'INTERVAL_ARITHMETIC appears to be installed correctly.'
      ELSE
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         'Some of the tests failed.  This could be due to one'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         'or more of the following:'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '1. The machine constant routine D1MACH is not installed'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   properly in the Fortran 77 package INTLIB.  See the'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   Fortran 90 version with this package.'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '2. The arithmetic is not IEEE arithmetic, and the value'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   MAXERR in INTLIB routine SIMINI should be some number'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   other than 1.'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '3. The "exact result" in the testing routine is'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   represented as as two double precision numbers.'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   Conversion to the internal representation may not'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   be as exact as the interval operations.'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '4. Some "exact results" are irrational numbers, and are'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   represented to only 30 digits.  If double precision'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   contains more than 30 digits, then the "exact result"'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   may not be sufficiently precise.  In this case, the'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   constants in the testing routine "testsys.f90" should'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   be made more precise.  (In this case, the INTLIB'
         WRITE(OUTPUT_UNIT,'(1X,A)') &
         '   constants in SIMINI should also be made more precise.)'
      END IF

      CLOSE(OUTPUT_UNIT)

CONTAINS

   SUBROUTINE CHECK_RESULT(A,B,STRING)
      TYPE(INTERVAL) :: A, B
      CHARACTER(LEN=*) STRING
      LOGICAL CHECKS_OUT

      CHECKS_OUT = A%LOWER.LE.B%LOWER .AND. A%UPPER.GE.B%UPPER

      IF(.NOT.CHECKS_OUT) THEN
         ALL_TESTS_OK = .FALSE.
         WRITE(OUTPUT_UNIT,'(2(1X,A))') STRING, &
            'did not contain the exact result.'
         WRITE(OUTPUT_UNIT,'(5X,A)') 'COMPUTED RESULT:'
         WRITE(OUTPUT_UNIT, '(5X,2(1X,ES28.20E2))') A
         WRITE(OUTPUT_UNIT,'(5X,A)') 'EXACT_RESULT:'
         WRITE(OUTPUT_UNIT, '(5X,2(1X,ES28.20E2))') B
      END IF

   END SUBROUTINE CHECK_RESULT

   SUBROUTINE CHECK_REAL_RESULT(A,B,STRING)
      DOUBLE PRECISION :: A, B
      CHARACTER(LEN=*) STRING
      LOGICAL CHECKS_OUT

      CHECKS_OUT = ABS(B-A).LT.TOL*MAX(ABS(A),1D0)

      IF(.NOT.CHECKS_OUT) THEN
         ALL_TESTS_OK = .FALSE.
         WRITE(OUTPUT_UNIT,'(2(1X,A))') STRING, &
            'was not approximately equal to the exact result.'
         WRITE(OUTPUT_UNIT,'(5X,A,1X,ES28.20E2)') 'COMPUTED RESULT:', A
         WRITE(OUTPUT_UNIT,'(5X,A,1X,ES28.20E2)') '   EXACT_RESULT:', B
      END IF

   END SUBROUTINE CHECK_REAL_RESULT

   SUBROUTINE CHECK_LOGICAL_RESULT(A,B,STRING)
      LOGICAL :: A, B
      CHARACTER(LEN=*) STRING
      LOGICAL CHECKS_OUT

      CHECKS_OUT = A .EQV. B

      IF(.NOT.CHECKS_OUT) THEN
         ALL_TESTS_OK = .FALSE.
         WRITE(OUTPUT_UNIT,'(2(1X,A))') STRING, &
            'was not correct.'
         WRITE(OUTPUT_UNIT,'(5X,A,1X,L1)') 'COMPUTED RESULT:', A
         WRITE(OUTPUT_UNIT,'(5X,A,1X,L1)') '   EXACT_RESULT:', B
      END IF

   END SUBROUTINE CHECK_LOGICAL_RESULT


END PROGRAM TEST_INTERVAL_SYSTEM
