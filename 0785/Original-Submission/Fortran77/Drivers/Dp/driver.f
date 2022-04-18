C
C   THIS IS THE FINAL VERSION SUBMITTED TO ACM TOMS.
      PROGRAM DRIVER
C   THIS PROGRAM DRIVES DSCPACK. IT READS FROM THE FILE DSCDATA. ON
C   OUTPUT,IT GIVES THE COMPUTED ACCESSORY PARAMETERS.THE ACCURACY
C   OF THE PARAMETERS IS ALSO TESTED.
C
C
C     .. Scalars in Common ..
      DOUBLE PRECISION DLAM
      INTEGER IU
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION UARY(8),VARY(3)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX C,WW,ZZ,ZZ0
      DOUBLE PRECISION EPS,PI,TOL,U
      INTEGER I,IGUESS,INV,IPOLY,ISHAPE,ISOLV,ITRY,K,LINEARC,M,N,NPTQ
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX W0(30),W1(30),Z0(30),Z1(30)
      DOUBLE PRECISION ALFA0(30),ALFA1(30),PHI0(30),PHI1(30),QWORK(1660)
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX WDSC,ZDSC
      EXTERNAL WDSC,ZDSC
C     ..
C     .. External Subroutines ..
      EXTERNAL ANGLES,CHECK,DSCDATA,DSCSOLV,DSCTEST,QINIT,THDATA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ACOS
C     ..
C     .. Common blocks ..
      COMMON /PARAM4/UARY,VARY,DLAM,IU
C     ..
      PI = ACOS(-1.D0)
C
C  SPECIFY THE GOEMETRY OF THE POLYGON REGION(READ FROM DATA FILE):
   10 PRINT '(/)'
      WRITE (6,FMT=*) 'TYPE 1,2,3,4,5,6 OR 7 TO SELECT AN EXAMPLE'
      WRITE (6,FMT=*) 'AVAILABLE SELECTION:'
      WRITE (6,FMT=*) '1  A SQUARE-SYMMETRIC REGION'
      WRITE (6,FMT=*) '2  A MILDLY CROWDED INFINITE REGION'
      WRITE (6,FMT=*) '3  A HEAVILY CROWDED REGION'
      WRITE (6,FMT=*) '4  A CHINESE-CHARACTER-STRUCTURED REGION'
      WRITE (6,FMT=*) '5  A 4-DIRECTION-INFINITE REGION'
      WRITE (6,FMT=*) '6  AN EX. GIVEN BY REFEREE FOR CHECKING INV. MAP'
      WRITE (6,FMT=*) '7  THE UPPER HALF PLANE WITH A HORIZONTAL SLIT'
C
C  INITIALISING ARRAYS:
      DO 20 K = 1,30
          ALFA0(K) = 0.D0
          ALFA1(K) = 0.D0
   20 CONTINUE
      READ (5,FMT=*) IPOLY
      CALL DSCDATA(IPOLY,M,N,Z0,Z1,ALFA0,ALFA1)
      CALL ANGLES(N,Z1,ALFA1,1)
      IF (IPOLY.EQ.2 .OR. IPOLY.EQ.5 .OR. IPOLY.EQ.7) GO TO 30
      CALL ANGLES(M,Z0,ALFA0,0)
   30 CONTINUE
C
C  GENERATE THE GAUSS-JACOBI WEIGHTS & NODES AND CHECK THE INPUT:
      WRITE (6,FMT=*) 'INPUT NPTQ(THE NUMBER OF G-J POINTS)'
      WRITE (6,FMT=*)
     +  '(RECOMMENDED VALUES FOR NPTQ ARE: 2,3,4,5,6,7,OR 8)'
      READ (5,FMT=*) NPTQ
      CALL QINIT(M,N,ALFA0,ALFA1,NPTQ,QWORK)
      ISHAPE = 0
      IF (IPOLY.EQ.2 .OR. IPOLY.EQ.5 .OR. IPOLY.EQ.7) ISHAPE = 1
      CALL CHECK(ALFA0,ALFA1,M,N,ISHAPE)
C
C  SPECIFY SOME PARAMETERS OF THE CALLING SEQUENCE OF DSCSOLV:
      IGUESS = 1
      LINEARC = 1
      TOL = 1.D-10
C
C   SOLVE THE ACCESSORY PARAMETER PROBLEM:
      PRINT '(/)'
      WRITE (6,FMT=*)
     +  'SOLVE THE NONLINEAR SYSTEM? (1 FOR "YES", 2 FOR "NO")'
      READ (5,FMT=*) ISOLV
      IF (ISOLV.EQ.1) CALL DSCSOLV(TOL,IGUESS,M,N,U,C,W0,W1,PHI0,PHI1,
     +                             Z0,Z1,ALFA0,ALFA1,NPTQ,QWORK,ISHAPE,
     +                             LINEARC)
C
C   OUTPUT WILL BE ARRANGED IN DSCSOLV.
C
C   COMPUTE DATA FOR THETA-FUNCTION AND TEST THE ACCURACY:
      CALL THDATA(U)
      CALL DSCTEST(M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,NPTQ,QWORK)
C
C  FOLLOWING PARAGRAPH IS FOR TESTING INVERSE EVALUATIONS:
      WRITE (6,FMT=*)
     +  'WANT TO INVERSE THE MAP? (1 FOR "YES",2 FOR "NO")'
      READ (5,FMT=*) INV
      IF (INV.EQ.1) THEN
          DO 50 I = 1,10
              WRITE (6,FMT=*)
     +          'INPUT THE POINT INSIDE THE DOMAIN IN THE FORM OF'
              WRITE (6,FMT=*) '          ( X-COORD., Y-COORD. )'
              READ (5,FMT=*) ZZ
              EPS = 1.D-6
              WW = WDSC(ZZ,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,PHI1,
     +             NPTQ,QWORK,EPS,1)
              WRITE (6,FMT=*) 'THE PREIMAGE OF ZZ=',WW
              IF (ABS(WW).LE.1.D-12) GO TO 40
              WRITE (6,FMT=*) 'CHECK BY MAPPING THE PREIMAGE BACK'
              ZZ0 = ZDSC(WW,0,2,M,N,U,C,W0,W1,Z0,Z1,ALFA0,ALFA1,PHI0,
     +              PHI1,NPTQ,QWORK,1)
              WRITE (6,FMT=*) 'THE POINT ENTERED=',ZZ0
   40         WRITE (6,FMT=*)
     +          'WANT TO TRY ANOTHER PT ? (1 FOR "YES",2 FOR "NO")'
              READ (5,FMT=*) INV
              IF (INV.EQ.2) GO TO 60
   50     CONTINUE
      END IF
C
      PRINT '(/)'
   60 WRITE (6,FMT=*)
     +  'WANT TO TRY ANOTHER EXAMPLE? (1 FOR "YES",2 FOR "NO")'
      READ (5,FMT=*) ITRY
      IF (ITRY.EQ.1) GO TO 10
      STOP

      END
C.....................................................................
C THE FOLLOWING SUBROUTINE GENERATES DATA.
C
      SUBROUTINE DSCDATA(IPOLY,M,N,Z0,Z1,ALFA0,ALFA1)
C
C   SPECIFY INPUT DATA:
C
C      CALL KILLUN
C
C     .. Scalar Arguments ..
      INTEGER IPOLY,M,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX Z0(30),Z1(30)
      DOUBLE PRECISION ALFA0(30),ALFA1(30)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,Q,S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ACOS,DCMPLX,SQRT
C     ..
      PI = ACOS(-1.D0)
C
      IF (IPOLY.EQ.1) THEN
          M = 4
          N = 4
          Q = SQRT(2.D0)
          Z0(1) = DCMPLX(1.D0+Q,1.D0+Q)
          Z0(2) = DCMPLX(-1.D0-Q,1.D0+Q)
          Z0(3) = DCMPLX(-1.D0-Q,-1.D0-Q)
          Z0(4) = DCMPLX(1.D0+Q,-1.D0-Q)
          Z1(1) = DCMPLX(Q,0.D0)
          Z1(2) = DCMPLX(0.D0,Q)
          Z1(3) = DCMPLX(-Q,0.D0)
          Z1(4) = DCMPLX(0.D0,-Q)

      ELSE IF (IPOLY.EQ.2) THEN
          M = 12
          N = 6
          Z0(1) = (0.6875D0,0.875D0)
          Z0(2) = (0.D0,0.D0)
          Z0(3) = (1.D0,0.D0)
          Z0(4) = (0.875D0,0.875D0)
          Z0(5) = (1.125D0,1.375D0)
          Z0(6) = (2,1.375D0)
          Z0(7) = (1.25D0,2.D0)
          Z0(8) = (2.25D0,2.D0)
          Z0(9) = (2.375D0,2.75D0)
          Z0(10) = (1.625D0,2.25D0)
          Z0(11) = (1.125D0,2.625D0)
          Z0(12) = (-0.5D0,2.75D0)
          Z1(1) = (0.375D0,1.875D0)
          Z1(2) = (0.5D0,2.D0)
          Z1(3) = (1.D0,1.5D0)
          Z1(4) = (0.5D0,2.1875D0)
          Z1(5) = (0.5D0,2.5D0)
          Z1(6) = Z1(4)
C
          ALFA0(1) = 1.39169261159339475D0
          ALFA0(2) = 0.28801540784794967D0
          ALFA0(3) = 0.454832764699133488D0
          ALFA0(4) = 1.19275085295129979D0
          ALFA0(5) = 1.35241638234956651D0
          ALFA0(6) = 0.D0
          ALFA0(7) = 2.D0
          ALFA0(8) = 0.552568456711253445D0
          ALFA0(9) = 0.260264501477747753D0
          ALFA0(10) = 1.39199980651013222D0
          ALFA0(11) = 0.819604487273064009D0
          ALFA0(12) = 0.295854728586457991D0

      ELSE IF (IPOLY.EQ.3) THEN
          M = 11
          N = 6
          Z0(1) = (0.5D0,2.5D0)
          Z0(2) = (0.5D0,0.5D0)
          Z0(3) = (1.D0,0.5D0)
          Z0(4) = (1.D0,1.D0)
          Z0(5) = (1.D0,0.5D0)
          Z0(6) = (0.5D0,0.5D0)
          Z0(7) = (0.5D0,2.5D0)
          Z0(8) = (0.D0,2.5D0)
          Z0(9) = (0.D0,0.D0)
          Z0(10) = (2.D0,0.D0)
          Z0(11) = (2.D0,2.5D0)
          Z1(1) = (1.D0,2.D0)
          Z1(2) = (1.D0,1.5D0)
          Z1(3) = (1.D0,2.D0)
          Z1(4) = (1.5D0,2.D0)
          Z1(5) = (1.5D0,0.5D0)
          Z1(6) = (1.5D0,2.D0)
          ALFA0(1) = 1.D0/2.D0
          ALFA0(2) = 1.D0/2.D0
          ALFA0(3) = 1.D0/2.D0
          ALFA0(4) = 2.D0
          ALFA0(5) = 3.D0/2.D0
          ALFA0(6) = 3.D0/2.D0
          ALFA0(7) = 1.D0/2.D0
          ALFA0(8) = 1.D0/2.D0
          ALFA0(9) = 1.D0/2.D0
          ALFA0(10) = 1.D0/2.D0
          ALFA0(11) = 1.D0/2.D0
          ALFA1(1) = 3.D0/2.D0
          ALFA1(2) = 2.D0
          ALFA1(3) = 1.D0/2.D0
          ALFA1(4) = 1.D0/2.D0
          ALFA1(5) = 2.D0
          ALFA1(6) = 3.D0/2.D0

      ELSE IF (IPOLY.EQ.4) THEN
          M = 4
          N = 17
          Z0(1) = (-1.D0,-1.D0)
          Z0(2) = (1.D0,-1.D0)
          Z0(3) = (1.D0,1.D0)
          Z0(4) = (-1.D0,1.D0)
          Z1(1) = (0.D0,0.5D0)
          Z1(2) = (0.D0,0.D0)
          Z1(3) = (-0.5D0,0.D0)
          Z1(4) = Z1(2)
          Z1(5) = (0.D0,-0.5D0)
          Z1(6) = (-0.5D0,-0.5D0)
          Z1(7) = (0.5D0,-0.5D0)
          Z1(8) = (0.25D0,-0.5D0)
          Z1(9) = (0.25D0,-0.25D0)
          Z1(10) = Z1(8)
          Z1(11) = Z1(5)
          Z1(12) = Z1(2)
          Z1(13) = (0.5D0,0.D0)
          Z1(14) = Z1(2)
          Z1(15) = Z1(1)
          Z1(16) = (0.5D0,0.5D0)
          Z1(17) = (-0.5D0,0.5D0)

      ELSE IF (IPOLY.EQ.5) THEN
          M = 12
          N = 4
          S = 1.D0/6.D0
          Z0(1) = DCMPLX(0.D0,-8.D0*S)
          Z0(2) = (0.D0,0.D0)
          Z0(3) = (1.D0,-1.5D0)
          Z0(4) = (1.D0,-0.5D0)
          Z0(5) = (0.D0,0.D0)
          Z0(6) = DCMPLX(0.D0,2.D0*S)
          Z0(7) = (0.D0,0.D0)
          Z0(8) = (-1.D0,0.D0)
          Z0(9) = (0.D0,0.D0)
          Z0(10) = DCMPLX(-7.D0*S,-1.D0)
          Z0(11) = (0.D0,0.D0)
          Z0(12) = DCMPLX(-2.D0*S,-10.D0*S)
          ALFA0(1) = 1.75D0
          ALFA0(2) = 0.D0
          ALFA0(3) = 1.D0
          ALFA0(4) = 1.D0
          ALFA0(5) = (PI-3.5D0)/PI
          ALFA0(6) = 3.5D0/PI
          ALFA0(7) = 1.5D0
          ALFA0(8) = 1.D00
          ALFA0(9) = 0.D0
          ALFA0(10) = 1.75D0
          ALFA0(11) = 0.D0
          ALFA0(12) = 1.D0
          Z1(1) = DCMPLX(2.D0*S,-0.5D0)
          Z1(2) = DCMPLX(-S,-2.D0*S)
          Z1(3) = DCMPLX(-4.D0*S,-5.D0*S)
          Z1(4) = (0.D0,-1.D0)

      ELSE IF (IPOLY.EQ.6) THEN
          M = 7
          N = 2
          Z0(1) = (-2.D0,-1.D0)
          Z0(2) = (2.D0,-1.D0)
          Z0(3) = (2.D0,2.D0)
          Z0(4) = (-0.8D0,2.D0)
          Z0(5) = (1.D0,0.5D0)
          Z0(6) = (-1.D0,2.D0)
          Z0(7) = (-2.D0,2.D0)
          Z1(1) = (0.D0,0.D0)
          Z1(2) = (-1.D0,0.D0)

      ELSE
          M = 3
          N = 2
          Z0(1) = (1.01D0,0.D0)
          Z0(2) = (100.D0,100.D0)
          Z0(3) = (-1.01D0,0.D0)
          Z1(1) = (0.D0,2.D0)
          Z1(2) = (0.D0,1.D0)
          ALFA0(1) = 1.D0
          ALFA0(2) = -1.D0
          ALFA0(3) = 1.D0
      END IF

      RETURN

      END
