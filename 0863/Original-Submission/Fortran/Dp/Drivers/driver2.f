      PROGRAM DRIVER2
C
C     Double Precision
C
C     SAMPLE PROGRAM THAT GENERATES PSEUDORANDOMLY DATA FROM SIN(PI*X)
C      IN ORDER TO TEST THE SUBROUTINE L2WPMA.
C
C     CALLS FUNCTION RND.
C.......................................................................
C
C.... P A R A M E T E R S  (SET BY THE USER) ....
C     I1        INTEGER, LOWER DATA INDX (USUALLY I1 = 1).
C     NX        INTEGER, UPPER DATA INDX.
C     XA        REAL, LEFT LIMIT OF ABSCISSAE X(.), X(I1) = XA.
C     XB        REAL, RIGHT LIMIT OF ABSCISSAE X(.), X(N) = XB.
C
C.... I N P U T  (BY THE USER ON PROGRAM REQUEST) ....
C     N         INTEGER, NUMBER OF DATA POINTS TO BE GENERATED.
C     SIZE      REAL, RELATIVE MAGNITUDE OF NOISE.
C                RECOMMENTED VALUES FOR SIZE THAT COVER A WIDE RANGE
C                OF DATA BEHAVIOUR ARE
C                = 0
C                = 4/N
C                = 10/N
C                = 50/N
C                = 100/N
C                = 250/N
C
C.... O U T P U T  (DIRECTED TO THE FILE "XFWDAT") ....
C     X(I1:N)   DATA POINTS (ABSCISSAE).
C     F(I1:N)   FUNCTION MEASUREMENTS (INCLUDING NOISE).
C     WF(I1:N)  WEIGHTS ASSOCIATED WITH FUNCTION MEASUREMENTS.
C
C.... M E T H O D  (USER INTERFACE) ....
C     THIS PROGRAM GENERATES DATA AS FOLLOWS:
C      A CONTINOUS FUNCTION F(X), DEFINED ON [XA,XB], IS DEFINED
C      BY THE USER. THEN F(X) IS EVALUATED ON THE GRID
C      XA = X(I1) < X(I1 + 1) <... < X(N) = XB AND THEN RANDOM NUMBERS
C      FROM THE UNIFORM DISTRIBUTION OVER THE INTERVAL (-SIZE,SIZE)
C      ARE ADDED TO THE FUNCTION VALUES. THE TYPE OF THE FUNCTION,
C      THE LIMITS OF THE DATA AND THE SIZE OF THE MAGNITUDE OF THE
C      RANDOM NUMBERS ARE UPON THE USER'S DECISION. THE GRID NEED
C      NOT BE EQUALLY SPACED.
C
C     THE METHOD THAT PRODUCES THE PSEUDO-RANDOM NUMBERS MAY BE
C      FOUND IN FUNCTION RND(.).
C
C     AFTER THE DATA ARE GENERATED, THEY ARE DIRECTED TO FILE "XFWDAT",
C      WHICH IS INPUT TO PROGRAM DRIVER1.
C.......................................................................
C
C     .. Parameters ..
      INTEGER I1,NX
      DOUBLE PRECISION XA,XB
      PARAMETER (I1 = 1,NX = 2000,XA = 0.0,XB = 4.0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,SIZE,XSTEP,XVAL
      INTEGER I,ISEED,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION F(I1:NX),X(I1:NX),WF(I1:NX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION RND
      EXTERNAL RND
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DSIN
C
C     DISPLAY MESSAGES AND INPUT VALUES FOR N AND SIZE.
C
      PRINT 9000,XA,XB
      PRINT 9010,NX
      READ *,N
      PRINT 9020
      READ *,SIZE
C
C.... SET DATA POINTS..................................................
C
C      EVALUATE F(X)=SIN(PI*X), X in [XA,XB] WITH UNIFORM MESH SIZE.
C       THE MEASUREMENTS ARE KEPT IN THE ARRAYS X(I1:N) AND F(I1:N),
C       WHERE I1 <= N. THE WEIGHTS ARE SET EQUAL TO UNITY.
C
      PI = 4*DATAN(1.0D0)
      XSTEP = (XB - XA)/ (N - 1)
      XVAL = PI*XA
      X(I1) = XVAL
      F(I1) = DSIN(XVAL)
      WF(I1) = 1.0D0
      DO 10 I = I1 + 1,N - 1
         XVAL = XVAL + PI*XSTEP
         X(I) = XVAL
         F(I) = DSIN(XVAL)
         WF(I) = 1.0D0
   10 CONTINUE
      X(N) = PI*XB
      F(N) = DSIN(X(N))
      WF(N) = 1.0D0
C
C     GENERATE DATA BY ADDING NOISE TO THE FUNCTION MEASUREMENTS F(.).
C     ISEED IS SET BY THE USER. IT CAN BE ANY INTEGER SUCH THAT
C       1 .LE. ISEED .LE. 65535 .
C
      ISEED = 1
      DO 20 I = I1,N
         F(I) = F(I) + SIZE*RND(ISEED)
   20 CONTINUE
C
C     SEND (X,F,W) DATA TRIADS TO A DATAFILE.
C
      OPEN (2,FILE = 'XFWDAT')
      WRITE (2,FMT=9030) I1,N
      DO 30 I = I1,N
         WRITE (2,FMT=9040) X(I),F(I),WF(I)
   30 CONTINUE
      CLOSE (2)
      PRINT 9050
C
      STOP
C
 9000 FORMAT (//5X,'Data is going to be generated by adding random ',
     +       'noise ',/5X,'to measurements of the function SIN(pi*x), ',
     +       F5.1,' < = x < = ',F5.1)
 9010 FORMAT (//5X,'Just below, input number of data points N, ','wher',
     +       'e 0< N < = ',I5)
 9020 FORMAT (/5X,'Just below, input magnitude of noise ','(eg. 0.25, ',
     +       '0.5, 1.0 etc)')
 9030 FORMAT (2I5)
 9040 FORMAT (3E20.10)
 9050 FORMAT (/5X,'Data has been generated and kept in file "XFWDAT".')
      END
